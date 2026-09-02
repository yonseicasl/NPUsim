# 2026-08-19 latency 재평가 Phase 1 수정 기록

기준 문서: [../../assessment/TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md](../../assessment/TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md)
선행 기록: [reevaluation_2026-08-19_p1_fixes.md](reevaluation_2026-08-19_p1_fixes.md)

Phase 1로 분류된 세 항목(L1/L2/L3)을 수정했다. 세 항목 모두 topology 확장이 아니라
현재 analytical 모델 안의 **수치 오류 또는 결과 계약 모순**이다.

## L1. Separate GLB base+fill datatype 대응관계 소실

### 문제

이전 수정(P1-D)은 chip 차원은 보존했지만 **datatype 차원을 base side에서 먼저 축약**했다.

```text
이전: max_chip( combine_type(base[chip][t]) + combine_type(fill[chip][t]) )
```

Shared GLB에서는 두 항이 모두 sum이라 `sum_type(base) + sum_type(fill) == sum_type(base+fill)`로
성립한다. 그래서 오류가 드러나지 않았다. 그러나 separate GLB에서는

```text
max_type(base) + max_type(fill) != max_type(base[t] + fill[t])
```

이고, 서로 **다른 partition**의 peak를 더하게 된다. Partition은 병렬로 동작하므로 존재하지
않은 elapsed time이 만들어진다.

### 수정

base와 fill을 모두 `[chip][type]`로 최종 scaling 시점까지 보존하고, datatype별로 먼저 더한다.

```text
현재: max_chip( combine_type( base[chip][t] + fill[chip][t] ) )
```

- [scheduler/stats.h](../../scheduler/stats.h): `chip_access_cycle_global_buffer`(신규, `[chip][type]`),
  `chip_fill_access_cycle_global_buffer`
- [scheduler/stats.cc](../../scheduler/stats.cc): `update_stats()`가 per-chip/per-datatype base를
  저장, `scale_serial_repetitions()`가 base는 균일·fill은 datatype별로 scaling,
  `merge_global_buffer_fill()`이 위 순서로 결합
- [utils/interconnect_timing.cc](../../utils/interconnect_timing.cc): 두 side를 datatype별로 더한 뒤
  결합하는 `entity_combined_cycles(a, b, serialized_types)` 추가

### 검증

Gemmini `gemm_64x64x64`(separate GLB) live 결과:

| 항목 | 수정 전 | 수정 후 |
|---|---:|---:|
| GLB per-type 최종 access (in/wt/out) | 512 / 512 / 1024 | 512 / 512 / 1024 |
| GLB access axis | **1,280** | **1,024** |
| Critical-path latency | 87,916 | 87,852 |
| Bottleneck level | PE array | PE array |

회귀:

- 수기 micro-test: shared/separate 각각, 그리고 위 live 값(base `(256,256,1024)` +
  fill `(256,256,0)` → 1,024, 이전 식은 1,280) — [unittest/validation_test.cc](../../unittest/validation_test.cc)
- live invariant **T10**: 결과 파일에 `GLB datatype rule` (sum/max)을 출력하고, GLB access axis가
  보고된 per-datatype access cycle의 datatype 결합값을 넘지 않는지 검사한다. 이 bound는 chip 수와
  무관하게 성립하며 single-chip에서는 등식이다. 수정 전 값 1,280은 bound 1,024를 넘으므로 실패한다.
  — [validation/traffic/check.py](../../validation/traffic/check.py)
- T10 검출력 확인: Eyeriss conv1 결과의 GLB axis를 23,633,280 → 23,700,000으로 임시 변경하면
  `T10 ... exceeds its datatype combination 23633280.0`로 실패한다. 원복 후 전체 통과.

## L2. Packet grouping과 recurrence

### 문제

이전 구현은 **aggregate transaction count를 균등 분배**해 group을 만들고,
`첫 group latency + 이후 group들의 period 합`으로 makespan을 닫았다. 두 가지가 틀렸다.

- **L2-A**: 균등 분배는 실제 line width가 만들지 않는 packet 경계를 발명한다.
  8b line / 256b link로 200 B를 옮기면 실제 grouping은 `32,32,32,32,32,32,8`인데
  균등 분배는 `29,29,29,29,28,28,28`을 만든다.
- **L2-B**: group 크기가 다르면 `first + period 합`이 닫히지 않는다. 값싼 tail packet에
  full packet의 period를 물리거나, 아직 밀려 있는 downstream stage에 gate되는 것을 놓친다.

### 수정

count-only API를 **physical packet API**로 교체했다.

```text
packet_bits = max(source line, link, destination line)
```

가장 넓은 stage의 transaction 하나가 모든 좁은 stage의 transaction을 정수 개 포함하므로,
이 지점이 곧 pack/unpack 단위이고 가장 넓은 stage는 packet마다 transaction 1개를 갖는다.
마지막 packet만 tail이다.

makespan은 **stage resource availability recurrence**로 계산한다.

```text
packet마다, 존재하는 stage를 순서대로:
  start  = max(자기 resource의 available, dependency)
    dependency = gather 경계(하위 count < 상위 count)  -> 상위 packet 전체 완료
                 fan-out 경계                          -> 상위의 첫 transaction
  finish = start + count*cost,  단 finish >= 상위 finish + cost
  available = finish
```

동일한 full packet 위에서 이 recurrence는 (max,+) 시스템이므로 `period = max(count*cost)`로
주기화된다. 전이 구간만 명시적으로 계산하고 나머지는 외삽하므로 packet 수가 payload에
비례해 커져도 비용이 늘지 않는다.

구현: [utils/interconnect_timing.h](../../utils/interconnect_timing.h) /
[utils/interconnect_timing.cc](../../utils/interconnect_timing.cc)의
`transfer_packet_groups_t`, `transfer_packet_groups()`, `pipelined_transfer_cycles(groups, ...)`.
`datatype_transfer_timing_t`가 `groups`를 함께 반환하므로 6개 call site
(pe, pe_array x2, global_buffer, multi_chip, dram)는 width를 이미 아는 지점에서 그대로 넘긴다.

물리적으로 없는 stage는 zero-cost가 아니라 **제거**한다(`without_source()`,
`without_destination()`): pass-through PE-array의 temporal buffer write, GLB-bypass된 GLB read,
공유 chunk를 다시 읽지 않는 NoP source read. 자기 몫만 옮기는 destination
(PE-array 분배: 각 PE는 array tile의 slice만 자기 LB에 씀)은 `with_destination_total()`로
공유 packet 경계 위에 총량을 분배한다(문서화된 근사).

또한 group별 count는 caller가 cycle/energy를 과금하는 aggregate count에 **정확히 합산되도록
reconcile**한다. makespan이 과금된 작업량과 다른 일을 서술하는 것을 구조적으로 막는다.

### 수기 회귀

모든 case를 `(payload bits, source line, link, destination line)`로 기술한다. cost는 `(5,1,2)`.

| case | packets | 값 | 근거 |
|---|---|---:|---|
| 32 B, 8b/256b/8b | 1 × (32,1,32) | 225 | 양방향 barrier, 32*5 + 1 + 32*2 |
| 64 B, 8b/256b/8b | 2 × (32,1,32) | 385 | 225 + period 160 |
| **200 B, 8b/256b/8b** | 6 × (32,1,32) + (8,1,8) | **1041** | 균등분배 식은 1059 |
| 13 B, 8b/256b/8b | 1 × (13,1,13) | 92 | packet보다 작은 payload |
| 32 B, 256b/32b/32b | 1 × (1,8,8) | 22 | fan-out, barrier 없음 |
| 32 B, 256b/32b/256b | 1 × (1,8,1) | 15 | fan-out 후 gather |
| 32 B, 64b 등폭 | 4 × (1,1,1) | 23 | CE5 등폭 퇴화 |
| 200 B, destination 제거 | — | 1001 | 마지막 link 착지 시점 |
| 200 B, source 제거 | — | 401 | tail packet drain |

`(65,2,65)`처럼 물리적으로 만들 수 없는 count 조합은 회귀에서 제거했다. 8b line과 256b link로
65 B를 옮기면 link transaction은 2개가 아니라 3개다. 대신 `packet_counts_consistent()`가
분해 결과가 aggregate count에 합산되는지 검사한다.

동일한 full packet 구간을 외삽하는 최적화는 **모든 packet을 그대로 시뮬레이션하는
brute-force reference와의 등가성**으로 고정했다(`packet_extrapolation_exact()`).
bottleneck stage 3종, packet 정수배가 아닌 payload, packet 65,536개 규모, stage 제거 변형까지
포함해 오차 없이 일치한다. 비어 있는 stage는 전진하지 않으므로 주기 판정에서 제외해야 하며
(그렇지 않으면 외삽이 영구히 막혀 loop가 payload에 비례한다), 이 경우도 회귀에 포함된다.

### live 영향

Compute schedule과 traffic counter는 불변이고, 과대 직렬화가 있던 overlap axis가 줄어든다.
Eyeriss conv1 DRAM overlap axis 19,560,387 → 19,539,267, critical path 86,233,532 → 86,212,412.

## L3. Network busy와 axes의 계약

### 문제

Layer scope에서는 `busy = max(자기 axes)`로 재구성되지만, network scope에서는

```text
network busy[stage] = sum_layer(layer busy[stage])
network axis[x]     = sum_layer(layer axis[x])
```

두 값이 **독립적인 합**이고 `sum_layer(max(...)) != max(sum_layer(...))`이다.
Eyeriss AlexNet에서 network GLB busy 640,561,336 vs `max(network GLB axes)` 640,421,248
(차이 140,088). 그런데 출력과 주석은 axes로 busy를 재구성할 수 있다고 설명하고 있었다.

### 수정 (평가 문서 권장안 1)

의미를 계약으로 확정하고 출력에 명시한다. `max(network axes)`로 network busy를 바꾸지는
않는다 — layer는 직렬 실행이므로 그렇게 하면 모든 stage의 실제 busy를 과소계상한다.

```text
Busy cycles (ratio of critical path)  [network: summed over layers]
Busy-cycle axes (access / link / overlap)  [network: per-axis work summed over layers; NOT reducible to busy]
Busy-cycle axes (access / link / overlap)  [layer: busy = max of these axes]
```

- [scheduler/stats.h](../../scheduler/stats.h): `network_rollup` flag와 두 scope의 계약 주석
- [scheduler/stats.cc](../../scheduler/stats.cc): `update_network_stats()`가 flag를 세우고
  `print_results()`가 scope를 명시. 재구성 가능하다고 설명하던 주석 정정
- [validation/check_timing.py](../../validation/check_timing.py): `check_network_rollup()`가
  `network busy == sum(layer busy)`, `network axis == sum(layer axis)`, format axis 합,
  그리고 두 scope의 계약 표기 자체를 검사한다
- [validation/traffic/check.py](../../validation/traffic/check.py): T7이 layer 계약 표기를 요구하도록
  하고, network에는 적용되지 않는 이유를 문서화

## 검증 결과

| 검증 | 결과 |
|---|---|
| Unit/config validation | 통과, 209개 config |
| Gemmini RTL timing | 통과, MAPE 4.40%, max 7.86% (baseline 불변) |
| Eyeriss silicon latency | 통과, MAPE 4.26%, max 6.39% (baseline 불변) |
| Eyeriss GLB/DRAM traffic | max 133.94% (불변, informational) |
| Traffic T1~T10 | 8개 workload 통과 |
| Multi-chip MC1~MC8 | 통과 |
| Buffering LB-A~LB-E | 통과 |
| Network rollup (L3) | 통과 |
| ASan/UBSan smoke | 통과 |

## 남은 항목

평가 문서의 Phase 2~4는 수치 오류가 아니라 모델 확장이므로 착수하지 않았다.

- **Phase 2**: L5 per-tile ready/available timeline과 depth 제한 buffer,
  L6 bypass direct stream의 datatype별 boundary, L7 NoP multicast/contention/back-pressure,
  L10 비대칭 multi-chip fixture(per-chip mapping 또는 synthetic stats injection 필요)
- **Phase 3**: L8 DRAM address/bank/row 또는 DRAMSim3 gate 확정,
  L9 systolic skew/drain/bubble·active-lane reduction·non-zero format micro-test
- **Phase 4**: L4 memory-inclusive `Critical-path latency`의 외부 reference 확보와
  compute-bound/memory-bound acceptance gate, external traffic 오차 축소, L11 지원 layer 확장

L4는 평가 문서가 "완료 선언 blocker"로 분류했지만 외부 reference 데이터 확보가 선행되어야
하므로 코드 수정만으로 닫을 수 없다.
