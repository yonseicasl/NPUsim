# 2026-08-19 재평가 P1 수정 기록

기준 문서: [../../assessment/TIMING_SIMULATION_REEVALUATION_2026-08-19.md](../../assessment/TIMING_SIMULATION_REEVALUATION_2026-08-19.md)

재평가에서 "critical-path cycle을 직접 왜곡할 수 있는 문제"로 분류된 P1 4건과,
P1 검증에 필요한 P2/P4 항목을 구현했다. 아래는 각 항목에서 **확정한 계약**과
그 계약을 실행 가능하게 검증하는 회귀다.

## P1-A. Single-buffer LB cycle 중복 계상

### 확정 계약

PE stage의 세 transfer axis는 **하나의 전송을 세 가지로 본 값**이며 서로 더할 수
있는 별개 작업량이 아니다.

```text
stage_axis_access[4]  : LB port 총 access cycle
stage_axis_link[4]    : LB<->MAC bus 총 link cycle
stage_axis_overlap[4] : 위 두 작업 + MAC register access의 pipelined MAKESPAN
```

`stage_axis_overlap[4]`(`cycle_mac_lb`)는 구조적으로 개별 axis를 모두 dominate한다.
Single buffering이 바꾸는 것은 이 전송의 **비용**이 아니라 이 전송을 **compute
뒤에 숨길 수 있는지**이므로 직렬화 식은 다음과 같다.

```text
double_buffer : busy_pe = max(compute, access, link, overlap, format)
single_buffer : busy_pe = compute + max(access, link, overlap) + format
```

구현: [scheduler/stats.cc](../../scheduler/stats.cc) `finalize_layer_timeline()`.
기존 식은 `compute + access + link + overlap + format`이었으므로 동일 transaction을
2~3회 과금했다.

### Capacity 계약 (재평가 항목 4)

PE의 `double_buffer`는 **timing 속성**이며 LB tile 사본을 두 개 요구하지 않는다.
근거는 다음과 같다.

- 이 flag가 지배하는 overlap은 `LB -> MAC register fill` 대 `MAC compute`이므로
  두 사본이 필요한 저장소는 MAC register file이다.
- 이 모델의 MAC register는 물리 lane 수(`number_of_macs*mac_width`)로 tile 하나를
  정확히 담도록 표현되며 2-deep register 모델이 없다.
- 두 사본을 요구하면 검증된 config가 전부 거부된다. Eyeriss conv3는 이미 input LB
  partition을 100% 사용하고, 모든 config의 MAC tile은 lane capacity를 정확히 채운다.
  즉 검증된 실리콘을 "모델링 불가"로 선언하게 된다.

따라서 `check_tile_size()`는 두 mode 모두 사본 1개를 요구하며, GLB의
`double_buffer`(실제로 두 tile 반쪽을 할당)와는 의미가 다르다. 이 차이는
[components/pe.cc](../../components/pe.cc) `pe_t::init()`에 명문화했고, 활성 mode를
layer timeline 결과에 출력해 run마다 보이게 했다.

### Config hook과 회귀

- accelerator config 11개 전부의 PE section에 `double_buffer`를 명시했다(기본값 1).
- A/B fixture: [configs/accelerators/gemmini_single_buffer.cfg](../../configs/accelerators/gemmini_single_buffer.cfg)는
  `gemmini.cfg`와 `double_buffer` 한 줄만 다르다.
- [validation/buffering/check.py](../../validation/buffering/check.py)가 LB-A~LB-E를 검사한다.
  LB-D는 수정 전 식(`sum of all axes`)보다 작고 double-buffer 값보다 큼을 동시에 요구하므로
  "직렬화가 보이는가"와 "중복 계상이 사라졌는가"를 함께 고정한다.

측정값: `compute 3518`, `overlap 9264` -> single `12782`, 수정 전 식은 `18942`.

## P1-B. GLB bypass end-to-end 계약

### 확정 계약 (재평가 선택지 2)

bypass는 **GLB storage를 우회**하며 chip의 물리 interconnect는 우회하지 않는다.

| 자원 | bypass 시 |
|---|---|
| multi-chip -> GLB fill(write) access | 과금하지 않음 |
| GLB read access | 과금하지 않음 |
| GLB <-> PE-array link transaction | **과금함** |
| PE-array temporal buffer write | **과금함** |
| NoP source read / link | 과금함 (원래 위치에서) |

수정 전 의미론은 정확히 반대였다. GLB fill을 유지하고 GLB->PE-array link/destination을
제거했으므로, 우회한다고 선언한 storage를 과금하면서 실제로 일어나는 전달은 무료로
만들었고, bypass된 stream에는 T4/T5가 검사할 traffic이 하나도 남지 않았다.

또한 `account_output_writeback_link()`의 early return은 bypass된 OUTPUT의 off-chip
store 전체를 삭제하고 있었다. 이제 GLB read만 제외하고 multi-chip write와 NoP link는
유지한다.

구현: [components/global_buffer.cc](../../components/global_buffer.cc),
[components/multi_chip.cc](../../components/multi_chip.cc).

### 보고와 검증

bypass된 datatype의 fabric traffic은 물리적으로 실재하지만 SRAM access는 아니다.
두 값을 혼동하지 않도록 결과 파일에 bypass된 datatype을 명시한다.

```text
GLB-bypassed (direct stream) :       weight
```

- [validation/check_timing.py](../../validation/check_timing.py)의 Eyeriss GLB traffic
  비교(JSSC 2017은 SRAM access를 센다)는 이 목록을 읽어 해당 datatype을 제외한다.
  따라서 외부 기준 대비 수치는 수정 전과 완전히 동일하다.
- [validation/traffic/check.py](../../validation/traffic/check.py)는 T4/T5의 bypass
  면제를 삭제했고(이제 bypass된 weight도 T5 = 1.000으로 검증됨), T9를 추가해
  "bypass = SRAM access 0 + fabric traffic 비영" 계약을 직접 검사한다.

## P1-C. Packet-group width conversion

### 확정 계약

downstream transaction 수가 줄어드는 경계는 barrier지만 그 dependency는 **stream 전체가
아니라 packet 단위**다. stream을 `packets = min(비영 stage count)` 개의 독립 group으로
나눈다. 가장 좁은 stage가 group마다 정확히 1개 transaction을 담당하며, 이것이 곧
assemble/gather 단위다.

```text
makespan = group0의 barrier-aware 내부 makespan
         + sum(나머지 group의 bottleneck stage period)
```

group별 transaction은 가능한 한 균등 분배하고 큰 group을 앞에 둔다(tail group이 저렴한
쪽이 된다). period는 group index에 대해 조각별 상수이므로 `O(packets)` loop 없이
`O(1)` 구간 합으로 닫는다.

구현: [utils/interconnect_timing.cc](../../utils/interconnect_timing.cc)
`pipelined_transfer_cycles()` / `packet_group_transfer_cycles()`.
등폭 overload는 이제 이 경로의 퇴화 case로 forward한다.

### 수기 회귀 ([unittest/validation_test.cc](../../unittest/validation_test.cc))

| case (src/link/dst, cost 5/1/2) | 값 | 근거 |
|---|---:|---|
| 32 / 1 / 32 | 225 | single packet, 32*5 + (1 + 32*2) |
| 13 / 1 / 13 | 92 | tail 포함 single packet |
| **64 / 2 / 64** | **385** | 2 packet, 225 + 160. 전역 barrier는 449로 과대평가 |
| **65 / 2 / 65** | **392** | 33+32 분배, 232 + 160 |
| 10 / 5 / 2 | 55 | 2 packet, 30 + 25 (전역 barrier는 59) |
| 8 / 8 / 0 | 41 | destination stage 제거(pass-through) |

## P1-D. CE4 GLB fill 집계 순서

### 확정 계약

GLB access axis는 다음 순서로만 결합한다.

```text
max_chip( combine_type(base[chip][t]) + combine_type(fill[chip][t]) )
```

chip을 먼저 collapse하면 `combine_type(max_chip(...))`가 되어 존재하지 않은 elapsed
time을 만든다. shared GLB에서 chip 0이 input fill 100, chip 1이 weight fill 100이면
잘못된 순서는 200을 보고하지만 두 chip은 병렬로 100에 끝난다.

base(GLB read side)는 균일 배수로, fill(multi-chip -> GLB write side)은 datatype별
배수로 scaling되므로 **entity 차원을 scaling 이후까지 보존**해야 한다. 이를 위해
`stats_t`에 per-chip 값을 유지한다.

구현: [scheduler/stats.h](../../scheduler/stats.h)의
`chip_combined_access_global_buffer`, `chip_fill_access_cycle_global_buffer`와
[scheduler/stats.cc](../../scheduler/stats.cc) `merge_global_buffer_fill()`.

### 검증

집계 순서 자체를 순수 함수 `entity_combined_cycles()`로 분리해 재평가 문서의 100/100
예제를 포함한 수기 unit test로 고정했다.

**비대칭 config fixture는 현재 모델에서 불가능하다.**
`global_buffer_t::update_tile_size()`가 모든 chip에 동일한 scheduler tile을 배정하므로
per-chip GLB 값은 구조적으로 동일하고 두 결합 순서는 항상 일치한다. 즉 이 수정은
per-chip mapping이 생겼을 때를 대비한 hardening이며, 순서 검증은 unit test가 담당한다.
이 제약은 [validation/multi_chip/check.py](../../validation/multi_chip/check.py) docstring에
기록했다.

## P2. Multi-chip 회귀 검출력 (항목 1~2)

대칭 fixture로는 broadcast/partitioned 차이를 볼 수 없으므로 **mirror 매핑 A/B**를
추가했다. 동일 하드웨어, 동일 총 연산량, split 차원만 반대다.

- `gemmini_2chip` : B split -> weight가 broadcast
- `gemmini_2chip_ksplit` : K split -> weight가 partitioned, input이 broadcast

source read가 공유되면 broadcast tensor의 processor tile은 chip tile과 같고 한 번
읽히며, partitioned tensor의 processor tile은 `chip 수 x chip tile`이고 chip마다
읽힌다. 따라서 **두 run의 datatype별 multi-chip access cycle이 정확히 같아야 한다.**

- MC6: 위 대칭성 (공유가 깨지면 B-split weight가 K-split보다 커진다)
- MC7: NoP **link**는 공유되지 않으므로 두 run의 link transaction이 mirror되어야 한다
  (MC6가 "전부 공유"로 통과하는 것을 막는다)
- MC8: 두 매핑의 총 연산량 동일

검출력 확인: `fresh_chunk`를 항상 true로 바꿔 공유를 제거하면 MC6가
`B-split weight 768 vs K-split 512`로 실패한다. 원복 후 통과.

또한 multi-chip temporal buffer의 datatype별 access cycle이 stage axis에만 들어가고
보고되지 않던 gap을 메웠다(`Access cycle (DRAM fill + NoP source reads)`).

## P4-1. Format axis를 network/report/checker에 연결

`stage_axis_format`은 `busy_cycle_pe`에 참여하는데도 network 합산과 출력에서 빠져
있었다. 현재 config의 format 단가가 0이라 회귀가 통과했을 뿐이다.

- `update_network_stats()`가 `stage_axis_format`을 합산하고 `pe_double_buffer`를
  AND-reduce한다.
- 결과 파일에 `PE format-IP axis`와 `PE local buffer` (double/single) 를 출력한다.
- traffic checker의 T7이 stage 4를 출력된 규칙 그대로 재구성한다. format 단가를
  켜더라도 계약이 유지된다.

## 검증 결과

| 검증 | 결과 |
|---|---|
| Unit/config validation | 통과, 209개 config |
| Gemmini RTL timing | 통과, MAPE 4.40%, max 7.86% (baseline 불변) |
| Eyeriss silicon latency | 통과, MAPE 4.26%, max 6.39% (baseline 불변) |
| Eyeriss GLB/DRAM traffic | 수정 전과 동일 (max 133.94%, informational) |
| Traffic T1~T9 | 8개 workload 통과, T5가 bypass weight에서 `nan` -> `1.000` |
| Multi-chip MC1~MC8 | 통과 |
| Buffering LB-A~LB-E | 통과 |
| ASan/UBSan smoke | 통과 |

## 남은 항목

재평가 문서의 나머지는 fix가 아니라 모델 확장이므로 착수하지 않았다.

- P2 항목 3~5: multicast link 계약, shared-link contention/arbitration,
  destination pressure와 stall propagation
- P3: per-tile `ready_cycle`/`available_cycle` timeline, buffer depth 기반 overlap 제한,
  tail tile, leakage/MAC availability 재연결
- P4 항목 2~5: active accumulator/lane 기반 reduction energy, address 기반 DRAM
  bank/row state 또는 DRAMSim3 경로, systolic skew/drain/bubble, 단가 calibration
- Deferred: sparse timing, 외부 데이터 의존 항목

CE4 비대칭 검증을 config fixture로 끝내려면 per-chip mapping 지원이 선행되어야 한다.
