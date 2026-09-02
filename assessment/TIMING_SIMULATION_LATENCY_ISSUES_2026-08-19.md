# NPUsim timing simulation latency 재평가 및 잔여 이슈

- 평가일: 2026-08-19
- 평가 기준: 현재 working tree(커밋되지 않은 수정 및 미추적 validation 파일 포함)
- 선행 문서:
  - [TIMING_SIMULATION_REEVALUATION_2026-08-19.md](TIMING_SIMULATION_REEVALUATION_2026-08-19.md)
  - [../implementation/timing/reevaluation_2026-08-19_p1_fixes.md](../implementation/timing/reevaluation_2026-08-19_p1_fixes.md)
- 포함: dense timing의 compute schedule, memory-inclusive critical path, stage busy/axis,
  buffer overlap, width conversion, NoC/NoP/DRAM 및 latency 검증 계약
- 제외: functional 값 정확성, dynamic-energy 단가/총계, sparse timing 구현
- 참고: leakage energy는 여기서 평가하지 않지만 최종 critical-path latency를 시간창으로
  사용하므로 latency 수정 후 반드시 재검증해야 한다.

## 최종 판정

현재 timing simulation에는 서로 다른 두 latency가 있다.

```text
Compute-schedule latency
  = Computation cycle + Fold fill cycle
  = RTL/silicon golden과 비교하는 공식 검증 metric

Critical-path latency
  = DRAM -> multi-chip -> GLB -> PE-array -> PE를 결합한 분석 wall-clock
  = 현재 외부 golden이 없는 informational metric
```

`Compute-schedule latency`는 Gemmini RTL과 Eyeriss silicon 회귀에서 안정화됐다.
그러나 사용자가 일반적으로 기대하는 end-to-end latency에 가까운
`Critical-path latency`에는 아직 수치 결함과 모델 공백이 남아 있다.

현재 가장 중요한 결론은 다음과 같다.

1. **Separate GLB의 base+fill 결합이 실제로 cycle을 과대계상한다.**
2. **Multi-packet width conversion은 full packet만 개선됐고 tail/group recurrence는 틀리다.**
3. **Network busy와 network axes는 각각 layer 합으로는 맞지만 서로 재구성되지 않는다.**
4. **GLB bypass, NoP, buffer overlap은 datatype별 실제 event timeline이 아니라 stage 총량으로
   결합되므로 direct stream, contention, back-pressure를 표현하지 못한다.**
5. **외부 검증은 compute schedule만 대상으로 하며 memory-inclusive critical path 자체는
   검증되지 않았다.**

따라서 현재 판정은 다음과 같다.

| 영역 | 현재 판정 |
|---|---|
| Dense single-chip compute schedule | 신뢰도 높음 / 외부 회귀 통과 |
| Layer 결과 계약과 공식 metric | 해결 |
| Dense memory hierarchy의 개별 transaction/accounting | 대체로 개선, 잔여 결함 있음 |
| Memory-inclusive critical-path latency | 미완료 |
| Multi-chip latency | 제한적 분석 모델 / 미완료 |
| Cycle-accurate 또는 event-accurate latency | 미지원 |

## 이번 재평가에서 새로 확인한 실제 cycle 결함

### L1. Separate GLB base+fill datatype 대응관계 소실

- 심각도: **높음 / P1**
- 영향: GLB busy와 critical path 과대계상 가능
- 현재 fixture에서 재현: **예**

`stats_t::update_stats()`는 chip별 GLB base access를 datatype vector로 보존하지 않고
먼저 하나의 scalar로 축약한다.

```text
shared GLB   : sum_type(base[type])
separate GLB : max_type(base[type])
```

그 뒤 `merge_global_buffer_fill()`은 이 scalar에 fill의 datatype 결합값을 더한다.

```text
현재 구현
  max_chip(
    combine_type(base[chip][type])
    + combine_type(fill[chip][type])
  )
```

Shared GLB에서는 두 항이 모두 합이므로 이 식이 맞다.

```text
sum_type(base) + sum_type(fill)
== sum_type(base + fill)
```

하지만 separate GLB에서는 일반적으로 성립하지 않는다.

```text
현재: max_type(base) + max_type(fill)
필요: max_type(base[type] + fill[type])
```

예를 들어 input partition의 base가 100이고 weight partition의 fill이 100이면 두 partition은
병렬로 100 cycle에 끝난다. 현재 식은 서로 다른 partition의 최대값을 더해 200을 만든다.

현재 Gemmini `gemm_64x64x64` 결과에도 이 문제가 보인다.

| 항목 | cycle |
|---|---:|
| GLB input 최종 access | 512 |
| GLB weight 최종 access | 512 |
| GLB output 최종 access | 1,024 |
| separate GLB의 올바른 access axis | 1,024 |
| 현재 출력된 GLB access axis | 1,280 |
| 과대계상 | **256** |

관련 구현:

- [../scheduler/stats.cc](../scheduler/stats.cc) `chip_combined_access_global_buffer`
- [../scheduler/stats.cc](../scheduler/stats.cc) `merge_global_buffer_fill()`
- [../utils/interconnect_timing.cc](../utils/interconnect_timing.cc) `entity_combined_cycles()`

현재 unit test는 `base`를 이미 datatype 결합이 끝난 entity scalar로 받으므로 이 오류를
검출할 수 없다. Per-chip/per-datatype base matrix를 repetition scaling 이후까지 보존해야 한다.

필요한 식은 다음과 같다.

```text
max_chip(
  combine_type(
    base[chip][type] + fill[chip][type]
  )
)
```

완료 조건:

1. base와 fill을 모두 `[chip][type]` 형태로 최종 scaling 시점까지 보존한다.
2. shared/separate 각각 수기 micro-test를 추가한다.
3. 현재 Gemmini case의 GLB access axis가 1,280이 아닌 1,024가 되는지 검증한다.
4. 수정 전후 critical-path와 bottleneck 변화를 기록한다.

### L2. Multi-packet tail의 grouping과 recurrence 오류

- 심각도: **높음 / P1**
- 영향: width conversion이 있는 전 계층의 transfer makespan 과대/과소 가능
- 현재 unit test에서 잘못된 expected value가 고정됨

현재 `pipelined_transfer_cycles()`는 다음 방식으로 group을 만든다.

```text
packets = min(nonzero stage transaction counts)
각 stage transaction을 packets개 group에 가능한 한 균등 분배
group 0의 전체 latency + 이후 group들의 자체 bottleneck period 합
```

이 방식은 `(64 source, 2 link, 64 destination)`처럼 모든 packet 크기가 같은 경우에는
전역 barrier를 제거하고 올바른 overlap을 만든다.

```text
32/1/32 packet 두 개
첫 packet latency 225 + 다음 packet source period 160 = 385
```

그러나 tail이 있는 경우 두 문제가 생긴다.

#### L2-A. Aggregate count를 균등 분배하면 실제 line-width packet 경계와 달라짐

예를 들어 source/destination line이 1 byte, link packet이 32 byte이고 총 200 byte라면
실제 grouping은 다음과 같다.

```text
32, 32, 32, 32, 32, 32, 8
```

현재 count-only helper는 `(200, 7, 200)`을 다음처럼 균등 분배한다.

```text
29, 29, 29, 29, 28, 28, 28
```

따라서 실제 full-packet/tail dependency를 복원할 수 없다. 현재
`datatype_transfer_timing()`은 endpoint/link width를 알고 있지만 makespan helper에는
transaction count만 전달한다.

#### L2-B. Variable-size group은 `first latency + period sum`으로 닫히지 않음

현재 unit test는 `(65,2,65)`, cost `(5,1,2)`를 `33+32`로 나누고 392 cycle을 기대한다.
하지만 test가 주장하는 바로 그 grouping을 resource availability recurrence로 배치하면
390 cycle이다.

```text
packet 1 source       0 -> 165
packet 1 link       165 -> 166
packet 1 destination166 -> 232

packet 2 source     165 -> 325
packet 2 link       325 -> 326
packet 2 destination326 -> 390
```

현실적인 `(200,7,200)`의 `32 x 6 + 8 tail`을 같은 recurrence로 계산하면 1,041 cycle이고,
현재 균등분배 수식은 1,059 cycle을 만든다.

관련 구현/회귀:

- [../utils/interconnect_timing.cc](../utils/interconnect_timing.cc)
  `pipelined_transfer_cycles()` / `packet_group_transfer_cycles()`
- [../unittest/validation_test.cc](../unittest/validation_test.cc) `(65,2,65) == 392`

필요한 수정:

1. element bit count와 source/link/destination width를 helper에 전달한다.
2. full packet과 마지막 tail packet의 정확한 transaction group을 만든다.
3. 각 stage의 `available_cycle`로 packet별 recurrence를 수행한다.
4. `(32,1,32)`, `(64,2,64)`, 실제 width 기반 tail, 양쪽 경계가 축소되는 case를 검증한다.
5. 인위적인 count 조합은 물리적으로 가능한 width/element 조합인지 함께 검사한다.

## 결과/검증 계약의 잔여 문제

### L3. Network busy와 network axes의 비가역 집계

- 심각도: **중간 / P1~P2**
- 영향: `network.txt`의 busy-axis 자기 일관성 및 downstream 분석
- layer latency 자체보다 network report 계약 문제에 가까움

Layer에서는 각 stage busy를 해당 layer axes로 재구성할 수 있다.

```text
memory stage busy = max(access, link, overlap)
PE busy           = buffer mode에 따른 compute/access/link/overlap/format 결합
```

Network에서는 다음 두 값을 독립적으로 합산한다.

```text
network busy[stage] = sum_layer(layer busy[stage])
network axis[x]     = sum_layer(layer axis[x])
```

일반적으로 다음 식은 성립하지 않는다.

```text
sum_layer(max(axis[layer])) != max(sum_layer(axis[layer]))
```

현재 Eyeriss AlexNet network 결과에서 실제 차이가 있다.

| 항목 | cycle |
|---|---:|
| Network GLB busy | 640,561,336 |
| max(network GLB access/link/overlap axes) | 640,421,248 |
| 차이 | **140,088** |

현재 `check_network_timeline()`은 network axes가 전부 0이 아닌지만 검사한다.
Traffic T7의 `busy == combined axes`는 layer 파일에만 적용된다.

이는 둘 중 하나로 계약을 확정해야 한다.

1. **권장:** network busy는 layer busy의 합, network axes는 단순 axis work 합으로 정의하고
   network axes만으로 busy를 재구성할 수 없음을 출력/문서에 명시한다.
2. Network에서도 재구성 가능한 새 통계를 원한다면 layer별 dominant-axis contribution 등
   정보를 별도 보존한다. 단순히 `max(network axes)`로 network busy를 바꾸면 serial layer의
   실제 busy 합을 과소계상하므로 허용하면 안 된다.

완료 조건:

- network busy와 axes의 의미가 출력에 명시된다.
- checker가 `network busy == sum(layer busy)`와 `network axis == sum(layer axis)`를 각각 검사한다.
- reconstructible하지 않은 값을 reconstructible하다고 설명하는 주석과 T7 문구를 수정한다.
- format axis와 single/double-buffer mode도 같은 network 계약을 따른다.

### L4. Memory-inclusive Critical-path latency는 외부 검증 대상이 아님

- 심각도: **높음 / 완료 선언 blocker**
- 영향: end-to-end latency의 절대 정확도 판단 불가

현재 외부 golden 비교는 `Critical-path latency`가 아니라 `Compute-schedule latency`를 사용한다.

현재 통과 수치:

| 검증 | 대상 metric | 결과 |
|---|---|---|
| Gemmini RTL 6개 GEMM | Compute-schedule | MAPE 4.40%, max 7.86% |
| Eyeriss silicon 5개 CONV | Compute-schedule | MAPE 4.26%, max 6.39% |

따라서 위 결과는 computation과 calibrated fold/setup schedule을 검증하지만, 다음 항목의
절대 정확도는 검증하지 않는다.

- DRAM access/link latency
- NoP 및 GLB latency
- PE-array/PE data movement latency
- stage overlap과 전체 critical path
- bottleneck level

특히 Eyeriss external traffic strict gate는 현재 `133.94% > 15%`로 실패한다. Memory
latency가 traffic count에 비례하는 현재 analytical 모델에서 외부 traffic 오차가 큰 상태로
memory-inclusive critical latency의 절대 정확도를 주장할 수 없다.

완료 조건:

1. `Compute-schedule latency`와 `Critical-path latency`의 지원/검증 범위를 계속 분리한다.
2. 최소 한 accelerator에 대해 memory-inclusive end-to-end latency reference를 확보한다.
3. component unit cycle을 reference에 맞춰 calibration하되 하나의 workload에만 과적합하지 않는다.
4. compute-bound와 memory-bound fixture를 모두 acceptance gate에 포함한다.

## 구조적 latency 모델 공백

### L5. Per-tile producer/consumer timeline 부재

- 심각도: **중간~높음 / P2**
- 영향: fill/drain, tail, buffer depth, stall propagation

`temporal_pipeline_run_cycles()`는 stage 총량을 repetition 수로 나눠 평균 tile 비용을 만든다.

```text
fill(sum of average stage costs)
+ (tiles - 1) * max(average stage cost)
```

이전 flat `max()`보다 one-time fill을 더 잘 표현하지만 실제 event timeline은 아니다.

표현하지 못하는 항목:

- tile별 비용 변화와 마지막 tail tile
- producer `ready_cycle`
- consumer/resource `available_cycle`
- buffer depth 1/2/N에 따른 허용 overlap
- queue occupancy와 producer stall
- consumer back-pressure
- DRAM부터 PE까지 이어지는 end-to-end interval

DRAM/multi-chip처럼 datatype별 repetition rate가 다른 stage를 포함하는 overlap run은 현재도
flat `max()`로 제한된다. 이것은 잘못된 비율을 섞지 않기 위한 방어이지만 정확한 timeline을
제공하지는 않는다.

필요한 최소 엔진:

```text
for each tile/packet:
  producer_start = producer_available
  producer_end   = producer_start + producer_cost(tile)
  consumer_start = max(producer_end, consumer_available, buffer_slot_available)
  consumer_end   = consumer_start + consumer_cost(tile)
```

모든 stage를 한 번에 cycle-accurate simulator로 바꾸기보다, 먼저 descriptor가 이미 알고 있는
tile/packet 크기를 event로 승격하고 각 boundary에 깊이 제한 queue를 두는 방식이 적합하다.

### L6. GLB bypass direct stream의 stage 경계가 datatype별로 표현되지 않음

- 심각도: **중간 / P2**
- 영향: bypass stream의 overlap 및 critical path

현재 bypass 계약 자체는 정리됐다.

```text
GLB SRAM fill/write : 없음
GLB SRAM read       : 없음
NoP link            : 있음
GLB<->PE-array link : 있음
PE-array destination: 있음
```

즉 “GLB storage만 우회하고 물리 fabric은 유지”한다. T9도 이 계약을 검증한다.

하지만 global timeline의 multi-chip↔GLB boundary overlap은 datatype별 bypass 여부가 아니라
GLB 전체의 `double_buffer` flag 하나로 결정된다. Bypassed stream에는 GLB storage가 없으므로
NoP delivery와 on-chip fabric delivery가 하나의 direct-stream pipeline이어야 하지만 현재 stage
총량 모델은 이를 별도로 표현하지 못한다.

따라서 P1-B accounting 계약은 해결됐지만 bypass의 end-to-end latency는 L5 timeline 과제로
남아 있다. Bypass on/off A/B가 다음을 수기 expected로 검증해야 한다.

- SRAM access 제거량
- 남는 NoP/on-chip link latency
- direct stream fill/drain
- 다른 resident datatype과의 shared-port contention
- 최종 critical-path delta

### L7. NoP multicast, contention, arbitration, back-pressure 부재

- 심각도: **중간 / P2**
- 영향: multi-chip latency의 절대 정확도

현재 구현된 범위:

- broadcast datatype의 shared source read를 distinct chunk 수에 따라 한 번만 과금
- chip별 link/destination은 routed-unicast로 과금
- mesh route는 ingress `(0,0)` 기준 Manhattan hop fill
- mirrored B/K split fixture에서 source sharing과 link mirror 검증

남은 범위:

- multicast tree/shared-link fanout
- 여러 route가 같은 physical link를 공유할 때 contention
- router arbitration과 queue
- virtual channel 또는 store-and-forward buffer
- destination-buffer pressure
- producer/consumer stall propagation

현재 MC1~MC8 통과는 source-read sharing과 aggregate transaction 계약의 검증이지 NoP
performance model 검증이 아니다.

### L8. DRAM bank/row/address timing은 제한적 analytical hook만 존재

- 심각도: **중간 / P2~P3**
- 영향: memory-bound workload latency

현재 analytical DRAM에는 다음 hook이 있다.

- row-buffer miss cost
- `t_ras_cycle`, `t_rp_cycle`
- `num_banks`

그러나 row miss를 bank에 round-robin으로 균등 분산한다고 가정한다. 다음은 없다.

- 실제 address 기반 channel/rank/bank/row mapping
- bank별 open-row state
- bank conflict와 queue
- read/write turnaround
- burst scheduling
- random/strided access

현재 accelerator config는 새 tRAS/tRP/bank hook을 활성화하지 않으므로 주요 regression은
기본 flat model만 검증한다. DRAMSim3 compile path가 존재하지만 기본 timing validation은
analytical path를 사용한다.

완료 방향은 둘 중 하나를 명확히 선택해야 한다.

1. Analytical 모델의 지원 범위를 sequential/conflict-free stream으로 제한하고 이를 출력한다.
2. Address trace와 DRAMSim3를 공식 memory-bound latency gate에 연결한다.

### L9. Systolic/PE reduction/format의 상세 pipeline 부족

- 심각도: **낮음~중간 / P3**

현재 systolic timing은 active array shape의 diameter fill과 평균 hop energy를 반영한다.
Gemmini에는 별도 fold-fill/setup calibration도 있다. 그러나 다음 event는 없다.

- operand skew
- per-PE pipeline register
- partial-sum propagation
- array drain
- tile-boundary bubble
- output accumulation과 다음 tile 주입의 충돌

PE lane reduction latency는 구조적 `mac_width`의 `ceil(log2(mac_width))`를 issue마다 더하지만
partial active lane의 실제 tree depth는 사용하지 않는다. Format-IP cycle은 layer/network axis에
연결됐지만 모든 현재 config의 단가가 0이므로 non-zero latency 회귀가 없다.

따라서 다음 설정에 대한 수기 latency fixture가 필요하다.

- `mac_width > 1`이며 마지막 accumulator가 partial인 case
- systolic active height/width와 tail tile을 함께 바꾸는 case
- non-zero format payload/metadata/spill cycle
- format stage와 LB transfer/compute의 overlap 규칙

## 검증 범위의 공백

### L10. Multi-chip 비대칭 live fixture 부재

현재 two-chip fixture의 computation cycle은 여전히 대칭이다.

```text
MIN = AVG = MAX = 512
```

Mirrored B/K split은 source-sharing 계약을 검증하지만 chip 0과 chip 1이 서로 다른
datatype/compute latency를 갖는 live 실행은 아니다. 현재 scheduler가 모든 chip에 동일한 GLB
tile을 배정하므로 config만으로 비대칭 GLB cycle을 만드는 것도 불가능하다.

순수 reduction unit test는 수식 hardening에는 유효하지만 다음을 검증하지 못한다.

- live scheduler에서 per-chip 값이 올바르게 생성되는가
- scaling 이후 entity/type 차원이 유지되는가
- 서로 다른 chip bottleneck이 critical path에 반영되는가

Per-chip mapping 또는 synthetic stats injection fixture가 필요하다.

### L11. 지원 workload와 layer 범위가 제한됨

현재 latency external validation 범위:

- Gemmini: dense GEMM, weight-stationary, 6개 point
- Eyeriss: dense AlexNet CONV 5개, row-stationary

AlexNet 실행은 convolution/connected 이외 5개 layer를 timing에서 제외하고 partial result로
표시한다. 따라서 출력된 network latency는 pooling, activation, normalization 등을 포함한
실제 end-to-end AlexNet latency가 아니다.

또한 다음은 external latency regression이 없다.

- Eyeriss FC의 latency
- TPU/TPUv3/Simba/MAERI/FSD/EyerissV2
- input/output/undefined stationary 조합
- multi-chip 실측 또는 cycle reference
- sparse/compressed execution

지원하지 않는 layer는 계속 명시적으로 제외할 수 있지만, “network end-to-end latency”라는
표현은 모든 layer가 지원될 때만 사용해야 한다.

## 해결된 항목 — 다시 열지 않아도 되는 현재 계약

다음 항목은 현재 범위와 회귀 안에서 해결된 것으로 유지한다.

### R1. 공식 compute metric과 layer/network identity

```text
Compute-schedule latency == Computation cycle + Fold fill cycle
```

Layer와 network 모두 runtime/checker identity gate가 있고 golden checker가 공식 metric을 직접
사용한다.

### R2. Single-buffer LB 중복 계상

```text
double buffer:
  busy_pe = max(compute, access, link, overlap, format)

single buffer:
  busy_pe = compute + max(access, link, overlap) + format
```

`overlap`은 동일 LB↔MAC transaction의 pipelined makespan이므로 access/link와 다시 더하지 않는다.
A/B fixture 결과:

| 항목 | cycle |
|---|---:|
| Compute schedule | 3,518 |
| Double-buffer PE busy | 9,264 |
| Single-buffer PE busy | 12,782 |
| 수정 전 중복 식 | 18,942 |

PE `double_buffer`는 현재 **timing-only overlap 계약**이며 두 LB tile 사본을 요구하지 않는다.
물리적인 2-deep MAC register/LB 모델은 L5 범위다.

### R3. GLB storage-bypass accounting

Bypassed datatype은 GLB SRAM read/write energy/cycle을 내지 않지만 NoP와 on-chip fabric을
통과한다. T9가 “GLB access 0 + fabric traffic non-zero”를 검사한다. 상세 direct-stream
latency만 L6으로 남는다.

### R4. Pass-through PE-array destination

Temporal buffer가 없는 PE-array는 존재하지 않는 destination buffer write cycle을 overlap
pipeline에 넣지 않는다.

### R5. Shared GLB의 chip/type fill 순서

Shared GLB에 대해서는 per-chip fill datatype을 결합한 뒤 chip 간 max를 선택하는 순서가
복원됐다. L1은 이 수정이 separate GLB의 base+fill type 대응까지 보존하지 못한 별도 문제다.

### R6. Full-size multi-packet overlap

`(64,2,64)`처럼 모든 packet group 크기가 같은 case는 전역 barrier 대신 packet 간 overlap을
허용해 449가 아닌 385 cycle을 만든다. Tail/variable group만 L2에 남는다.

### R7. Format axis의 기본 wiring

Format axis는 layer PE busy, network 합산, 결과 출력, layer T7에 연결됐다. Non-zero fixture와
network 계약은 L3/L9에 남는다.

## 최신 검증 실행 결과

이번 재평가에서 현재 바이너리로 다음을 다시 실행했다.

```bash
bash unittest/run_validation.sh
python3 validation/buffering/check.py
bash unittest/run_multi_chip_validation.sh
bash unittest/run_timing_validation.sh
python3 validation/traffic/check.py
python3 validation/check_timing.py --check-baseline --check-traffic
bash unittest/run_full_sanitizers.sh
```

| 검증 | 결과 | 해석 |
|---|---|---|
| Unit/config validation | 209개 config 통과 | helper/config 기본 계약 통과 |
| Buffering LB-A~LB-E | 통과 | single-buffer 중복 제거 확인 |
| Multi-chip MC1~MC8 | 통과 | source sharing과 mirrored link 계약 확인 |
| Gemmini RTL | MAPE 4.40%, max 7.86% | compute schedule 통과 |
| Eyeriss silicon latency | MAPE 4.26%, max 6.39% | compute schedule 통과 |
| Traffic T1~T9 | 8개 workload 통과 | 내부 transaction identity 통과 |
| ASan/UBSan smoke | 통과 | Gemmini/Eyeriss smoke에서 memory/UB 오류 없음 |
| Eyeriss traffic strict gate | **실패, max 133.94% > 15%** | 외부 memory traffic 정확도 미달 |

모든 회귀가 통과해도 L1/L2/L3가 남는 이유는 현재 checker가 해당 수식을 검증하지 않거나,
잘못된 tail expected value를 정답으로 고정하고 있기 때문이다.

## 권장 수정 순서

### Phase 1 — 현재 analytical 모델 안의 수치 오류 제거

1. **L1 separate GLB base+fill**
   - per-chip/per-type base 보존
   - `combine_type(base[type] + fill[type])`
   - shared/separate micro-test 및 Gemmini live regression
2. **L2 packet/tail recurrence**
   - count-only API를 width/element-aware packet API로 교체
   - packet별 stage availability recurrence
   - full/tail/both-boundary fixture
3. **L3 network rollup 계약**
   - busy와 axes를 각각 layer 합으로 검증
   - reconstructibility에 대한 출력/문서 정정

이 세 항목은 topology 확장 전에 끝낼 수 있고, 현재 존재하는 숫자 자체의 오류 또는 결과
계약 모순이므로 가장 먼저 처리해야 한다.

### Phase 2 — End-to-end analytical latency 계약 완성

4. **L5 최소 per-tile timeline**
   - ready/available cycle
   - depth-limited buffer slot
   - tail tile
5. **L6 bypass direct stream을 datatype별 boundary로 연결**
6. **L7 NoP multicast/shared-link contention 및 back-pressure**
7. **L10 비대칭 multi-chip fixture**

### Phase 3 — Detailed component latency

8. **L8 DRAM analytical 계약 확정 또는 DRAMSim3 gate**
9. **L9 systolic skew/drain/bubble, active-lane reduction, format pipeline**
10. component knob별 non-zero micro-test와 calibration

### Phase 4 — 외부 완료 판정

11. Memory-inclusive `Critical-path latency` reference 확보
12. compute-bound/memory-bound/multi-chip acceptance gate 추가
13. external traffic 오차 축소 또는 해당 traffic 범위를 명시적 비지원으로 제한
14. 지원 layer 확장 후 전체 network end-to-end latency 검증

## Timing latency 완료 판정 기준

다음 조건이 모두 충족되기 전에는 “timing simulation latency 완료”를 선언하지 않는다.

1. Shared/separate buffer에서 entity/type/base/fill 집계 순서가 수기식과 일치한다.
2. Equal-width, width conversion, full packet, tail packet makespan이 packet event 계산과 일치한다.
3. Layer와 network의 metric/busy/axis 의미가 모순 없이 checker로 고정된다.
4. Single/double/bypass/pass-through 설정이 수기 예측 가능한 critical-path delta를 만든다.
5. Per-tile ready/available과 buffer depth가 fill/drain/back-pressure를 결정한다.
6. Multi-chip route 공유와 physical-link contention이 latency에 반영된다.
7. DRAM latency의 지원 계약이 address/bank/row 수준 또는 명시적 제한으로 확정된다.
8. Systolic/adder/reduction/format의 non-zero knob이 수기 micro-test와 일치한다.
9. `Critical-path latency`가 최소 하나의 외부 end-to-end reference로 검증된다.
10. 제외 layer와 sparse 경로가 partial result 또는 unsupported로 명확히 보고된다.
11. Latency 변경 후 leakage 시간창과 MAC availability가 같은 final wall-clock으로 재검증된다.

## 종합 의견

현재 구현은 “검증된 compute schedule + 개선된 analytical memory timeline”으로 보는 것이
정확하다. 이전 P1 수정으로 single-buffer 중복과 GLB bypass accounting은 좋아졌지만,
separate GLB와 packet tail에서 아직 직접적인 cycle 오류가 확인된다. 또한 현재 외부 latency
수치는 memory-inclusive critical path가 아니라 compute schedule을 검증한다.

따라서 latency 작업은 다음 순서로 진행하는 것이 가장 안전하다.

```text
separate GLB base+fill
-> physical packet/tail recurrence
-> network busy/axis 계약
-> per-tile ready/available timeline
-> bypass/NoP contention/back-pressure
-> DRAM/systolic/format 상세화
-> memory-inclusive external latency calibration
```

이 순서는 먼저 현재 analytical contract의 수치 오류를 제거한 뒤, 같은 비용식을 event
timeline에 올리고, 마지막에 외부 end-to-end latency로 완료를 판정하도록 구성한다.
