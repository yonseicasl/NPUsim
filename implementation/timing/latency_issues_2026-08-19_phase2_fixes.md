# 2026-08-19 latency 재평가 Phase 2 구현 기록

기준 문서: [../../assessment/TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md](../../assessment/TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md)
선행 기록: [latency_issues_2026-08-19_phase1_fixes.md](latency_issues_2026-08-19_phase1_fixes.md)

Phase 2(L5/L6/L7/L10)를 구현했다. Phase 1이 analytical 모델 안의 수치 오류를 제거한 것이었다면
Phase 2는 **평균 per-tile 폐형식을 실제 event timeline으로 교체**하고, 그 timeline 위에서
bypass·NoP 계약을 표현하는 작업이다.

## L5. Per-tile producer/consumer timeline

### 이전 상태

`temporal_pipeline_run_cycles()`는 stage 총량을 repetition 수로 나눠 평균 tile 비용을 만든 뒤
`fill(sum) + (tiles-1)*max` 를 계산했다. 여기에는 다음이 없었다.

- consumer/resource `available_cycle`
- buffer depth에 따른 허용 overlap
- consumer back-pressure와 producer stall
- tile별 비용 차이 / tail tile
- rate가 다른 stage의 참여 (DRAM/multi-chip을 포함하는 run은 flat `max()`로 후퇴했다)

### 구현

`pipeline_timeline_cycles()` — stage chain 위의 per-tile event recurrence.

```text
start(s,k)  = max( finish(s,k-1),              // stage는 하나의 직렬 resource
                   finish(s-1,k),              // 입력 tile이 도착해야 한다
                   finish(s+1,k-depth[s]) )    // downstream slot이 비어야 한다
finish(s,k) = start(s,k) + cost(s,k)
```

세 번째 항이 평균 폐형식으로는 표현할 수 없던 **back-pressure**다. 느린 consumer 앞의 빠른
producer는 queue가 차면 멈추며, 병렬로 끝까지 달리지 않는다.

**buffer depth가 기존 boolean overlap 플래그의 일반화다.**

| depth | 의미 | 기존 대응 |
|---|---|---|
| 1 | single buffer. producer는 consumer가 tile k를 가져가기 전에 k+1을 시작할 수 없다 | "이 경계는 overlap하지 않음" |
| 2 | double buffer. 한쪽을 채우며 다른 쪽을 비운다 | "이 경계는 overlap함" |
| N | N-deep queue | 표현 불가였음 |

균일 per-tile 비용에서 depth 2는 기존 `fill+bottleneck` 값을 그대로 재현하고(unit test로 고정),
2-stage depth 1은 stage 총량의 합 — 즉 기존 non-overlap 분기의 값 — 을 그대로 만든다. 따라서
두 기존 case는 같은 engine의 특수해이며, **새로 추가된 것**은 다음이다.

1. **rate가 다른 stage의 참여.** off-chip stage(DRAM/multi-chip)는 datatype별로 refetch하므로
   on-chip stage가 소비하는 tile 수보다 적게 서비스한다. 공유 tile clock 위에서 자기가 서비스하는
   tile에만 비용을 싣고 나머지는 0을 실어(`offchip_repetition_tiles`) timeline에 참여한다.
   이전에는 이런 stage가 run에 들어오면 전체가 flat `max()`로 후퇴했다.
2. **tile별 비용 차이 / tail tile.** engine이 비용 vector를 받는다.
3. **stall attribution.** stage별로 "downstream buffer가 차서 막힌 시간"을 보고한다.

메모리는 stage당 ring buffer(`depth+2`)만 쓰므로 tile 수(GLB repetition)에 무관하다.

구현: [utils/interconnect_timing.h](../../utils/interconnect_timing.h) /
[utils/interconnect_timing.cc](../../utils/interconnect_timing.cc) `pipeline_timeline_cycles()`,
[scheduler/stats.cc](../../scheduler/stats.cc) `finalize_layer_timeline()`.

### 출력 추가

```text
Buffer depth (tiles in flight across each boundary)
 * DRAM -> Multi-chip :          1 tiles
 * Multi-chip -> GLB  :          1 tiles
 * GLB -> PE array    :          2 tiles
 * PE array -> PE     :          2 tiles
Back-pressure stall (blocked by a full downstream buffer)
 * DRAM               :     1152.0 cycles
 ...
```

depth는 overlap 계약 그 자체이므로 config flag 안에 묻어두지 않고 출력한다. network scope에서는
depth의 min을 보고한다(경계는 가장 덜 decouple된 layer만큼만 decouple된다).

### live 영향

Gemmini `gemm_64x64x64` critical path 87,852 → 59,010. 이전 모델은 non-overlap 경계에서 **segment
전체를 합산**했는데, 이는 3개 이상 stage chain에서 과도한 직렬화다. single buffer라도 producer
tile k+1은 consumer tile k만 기다리며 그 아래 stage 전체를 기다리지 않으므로, 실제로는 staggered
pipeline이 된다. Compute-schedule metric과 traffic counter는 불변이고 golden gate도 불변이다.

## L6. Bypass direct stream을 boundary 계약으로 연결

bypass **accounting** 계약은 Phase 1(P1-B)에서 정리됐다. 남은 문제는 global timeline의
multi-chip↔GLB overlap이 datatype별 bypass 여부가 아니라 GLB 전체의 `double_buffer` flag 하나로
결정된다는 것이었다.

### 확정 계약

bypass된 datatype은 GLB에 **buffer가 없다.** 따라서 그 stream은 tile을 붙들고 다음 tile을 받을 수
없다. 하나라도 bypass된 datatype이 있으면 GLB의 `double_buffer` flag와 무관하게 multi-chip→GLB
경계의 tile-level decoupling이 사라진다(depth 1).

구현: [scheduler/stats.cc](../../scheduler/stats.cc) `finalize_layer_timeline()`의 depth 결정.

### 검증 — 왜 latency 비교만으로는 안 되는가

GLB `double_buffer`는 **accounting도** 바꾼다(GB3: double-buffered destination은 fill access를
숨긴다). 그래서 latency만 비교하면 "bypass가 decoupling을 없앴다"와 "flag가 fill을 숨겼다"를
구분할 수 없다. 그래서 timeline이 실제로 사용한 **depth를 출력**하고 그것을 직접 검사한다.

[validation/bypass/check.py](../../validation/bypass/check.py)는 두 설정을 교차한 4개 config
(bypass × GLB double_buffer)로 BP1~BP6을 검사한다.

| config | multi-chip→GLB depth | weight GLB access |
|---|---:|---:|
| bypass + double_buffer=1 | **1** | 0 |
| bypass + double_buffer=0 | 1 | 0 |
| resident + double_buffer=1 | **2** | 256 |
| resident + double_buffer=0 | 1 | 512 |

`bypass_db`의 depth 1과 `resident_db`의 depth 2가 같은 flag에서 갈리는 것이 L6 계약이며,
`resident_db`가 2를 보이므로 flag를 무시해서 통과하는 것도 아니다.

- BP1: bypass된 weight의 GLB access = 0, resident는 > 0
- BP2: GLB↔PE-array weight link transaction이 4개 run 모두 동일 (bypass는 SRAM만 우회)
- BP3: NoP weight transaction이 4개 run 모두 동일 (tile은 어차피 package를 건넌다)
- BP4: 위 depth 표
- BP5: resident+sb > resident+db > 0 (double-buffered destination이 fill을 숨긴다)
- BP6: SRAM 작업을 없애는 것이 layer를 늘리지 않는다

## L7. NoP multicast와 shared-link 계약

### 이전 상태

broadcast datatype(split 차원에 의존하지 않는 tile)의 source read는 한 번만 과금했지만
**link transaction과 hop fill은 chip마다** 합산했다. bus에서 이는 단순히 틀렸다 — 한 번의 전송이
물리적으로 모든 receiver에게 도달한다.

### 확정 계약

package ingress가 공유 bottleneck이며, 이것이 두 delivery 계약을 나눈다.

| 계약 | ingress를 건너는 copy 수 | 총 link traversal (energy) |
|---|---|---|
| routed unicast (chip마다 다른 chunk) | chip 수 | 각 route의 hop 합 |
| multicast (모든 chip이 같은 tile) | **1** | multicast tree의 link 수 |

route depth(fill)는 **delivery마다 한 번** 과금한다. 동시 route의 hop은 pipeline되며 서로 뒤에
줄서지 않는다(SP1). 이전에는 chip마다 fill을 더했다.

`nop_multicast` config 설정을 추가했다(기본 1). multicast router가 없는 package는 0으로 모델링한다.
single-chip config는 chip 수 1이므로 어느 쪽이든 동일 — golden baseline은 불변이다.

구현: [utils/interconnect_timing.cc](../../utils/interconnect_timing.cc) `nop_delivery_cost()`,
[components/multi_chip.cc](../../components/multi_chip.cc) `account_descriptor_dense_distribution()`.

### 검증

수기 unit test: BUS unicast/multicast, MESH 1×4 line (unicast 7 traversal / multicast 4,
fill 2), MESH 2×2 (unicast 5 / multicast 4).

live A/B — `gemmini_2chip` vs `gemmini_2chip_unicast`(같은 하드웨어·같은 매핑, `nop_multicast`만 다름):

| datatype | multicast | unicast |
|---|---:|---:|
| broadcast weight | 256 | **512** (= 256 × 2 chips) |
| partitioned input | 256 | 256 (불변) |

MC7을 이 A/B로 교체했다. 기존 MC7(link mirror)은 multicast 적용 후 값이 같아져 판별력을 잃는다 —
2배 큰 broadcast tile을 한 번 보내는 비용이 1배 tile을 2 chip에 보내는 비용과 같아지기 때문이다.
mirror 자체는 MC9로 남겼고, 판별력이 없다는 점을 checker docstring에 명시했다.

### 남은 범위

router arbitration, virtual channel, store-and-forward buffer는 구현하지 않았다. destination-buffer
pressure와 producer/consumer stall propagation은 L5 engine의 boundary depth가 담당한다.

## L10. 비대칭 multi-chip fixture

### 왜 config만으로는 불가능한가

`global_buffer_t::update_tile_size()`가 모든 chip에 동일한 scheduler tile을 배정하므로 live per-chip
값은 구조적으로 같고, 올바른 entity/datatype 결합 순서와 잘못된 순서를 구분할 수 없다. 이를 바꾸려면
per-chip mapping 지원이 필요하며 그것은 test fixture가 아니라 모델 기능이다.

### 구현한 것 — synthetic stats injection

[unittest/stats_timeline_test.cc](../../unittest/stats_timeline_test.cc)는 `stats_t`에 비대칭
per-chip GLB 값을 주입하고 **production 경로**
(`scale_serial_repetitions()` → `merge_global_buffer_fill()` → `finalize_layer_timeline()`)를 그대로
실행한다. 순수 reduction unit test가 검증하지 못하는 다음을 덮는다.

- 두 GLB side가 **서로 다른 repetition factor**로 scaling되는 동안 entity/datatype 차원이 유지되는가
- live 코드가 scaling 종료 시점에 실제로 적용하는 결합 순서가 맞는가
- 양쪽을 합친 뒤에야 bottleneck이 되는 chip이 stage axis와 critical path에 도달하는가

수기 예제 (2 chip, base는 균일 ×2, fill은 datatype별 ×(1,2,1)):

```text
chip 0 base (10, 0,0) fill (0,20,0) -> scaling 후 (20,0,0) + (0,40,0) = (20,40,0)
chip 1 base (0,30,0)  fill (40,0,0) -> scaling 후 (0,60,0) + (40,0,0) = (40,60,0)

SEPARATE: max_chip(max_type(base+fill)) = max(40,60) = 60
          잘못된 순서 max_chip(max_type(base)+max_type(fill)) = max(60,100) = 100
SHARED  : 두 형태 모두 100 (separate에서만 드러나는 이유)
```

chip 순서 교환에 대한 불변성도 검사한다.

검출력 확인: helper를 L1 수정 전 형태로 되돌리면 `separate GLB: GLB access axis after scaling
(measured 100, expected 60)`으로 실패한다. 원복 후 통과.

여전히 검증하지 못하는 것: **live scheduler가 per-chip 값을 올바르게 생성하는가**. 이는 per-chip
mapping이 생긴 뒤에만 닫을 수 있으며 checker docstring에 명시했다.

## 검증 결과

| 검증 | 결과 |
|---|---|
| Unit/config validation | 통과, 219개 config |
| Gemmini RTL timing | 통과, MAPE 4.40%, max 7.86% (baseline 불변) |
| Eyeriss silicon latency | 통과, MAPE 4.26%, max 6.39% (baseline 불변) |
| Eyeriss GLB/DRAM traffic | max 133.94% (불변, informational) |
| Traffic T1~T10 | 8개 workload 통과 |
| Multi-chip MC1~MC9 | 통과 (MC7 multicast A/B 신규) |
| Buffering LB-A~LB-E | 통과 |
| Bypass BP1~BP6 | 통과 (신규) |
| Asymmetric stats/timeline (L10) | 통과 (신규) |
| Energy E1~E5 | 통과 |
| Network rollup (L3) | 통과 |

## 남은 항목

### Phase 3 — Detailed component latency
- **L8** DRAM: address 기반 channel/rank/bank/row mapping, bank별 open-row state, bank conflict,
  read/write turnaround, burst scheduling. 또는 analytical 범위를 sequential/conflict-free로 명시
  제한하거나 DRAMSim3를 공식 gate로 연결.
- **L9** systolic operand skew / per-PE pipeline register / partial-sum propagation / array drain /
  tile-boundary bubble, active-lane 기반 reduction depth, non-zero format 단가 micro-test.

### Phase 4 — 외부 완료 판정
- **L4** memory-inclusive `Critical-path latency`의 외부 reference 확보와 calibration.
  현재 외부 golden은 compute schedule만 검증한다. **reference 데이터 확보가 선행되어야 하므로
  코드 수정만으로 닫을 수 없다.**
- compute-bound / memory-bound / multi-chip acceptance gate 추가
- external traffic 오차 축소 또는 해당 범위를 명시적 비지원으로 제한
- **L11** 지원 layer 확장 후 전체 network end-to-end latency 검증

### Phase 2에서 의도적으로 남긴 것
- NoP router arbitration, virtual channel, store-and-forward buffer (L7 일부)
- per-chip mapping — L10의 "live scheduler가 per-chip 값을 생성하는가"를 닫는 전제
- tail tile은 engine이 표현하지만 현재 mapping은 dimension을 padding하므로 실제로는 균일하다
  (padding warning으로 보고됨). 비균일 tile은 mapper가 tail tile을 만들 때 활성화된다.
