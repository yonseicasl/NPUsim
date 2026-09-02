# NPUsim timing simulation 재평가

- 평가일: 2026-08-19
- 기준 문서: [TIMING_SIMULATION_REAUDIT_2026-08-18.md](TIMING_SIMULATION_REAUDIT_2026-08-18.md)
- 평가 대상: 현재 working tree의 timing simulation 구현과 최신 실행 결과
- 포함: dense timing, cycle/axis 집계, overlap, buffer/NoC/NoP/DRAM, runtime format, timing 관련 energy hook, 검증 계약
- 제외: functional simulation 정확성
- 주의: 이 평가는 커밋된 `HEAD`만이 아니라 현재 수정·미추적 파일을 포함한 working tree를 기준으로 한다.

## 최종 판정

2026-08-18 재감사에서 가장 높은 우선순위였던 P0 결과 계약은 해결됐다. Network의 compute schedule과 busy-cycle axis가 레이어 결과의 합으로 채워지고, golden checker도 공식 `Compute-schedule latency`를 직접 검증한다. Gemmini RTL과 Eyeriss silicon compute-schedule 회귀 및 sanitizer도 통과한다.

그러나 timing simulation 전체 완료를 선언할 수는 없다. Dense single-chip compute-schedule 경로는 신뢰도가 높지만, 다음 영역에는 cycle 정합성 또는 검증 범위의 공백이 남아 있다.

1. single-buffer LB cycle의 중복 계상
2. end-to-end가 아닌 GLB bypass 구현
3. multi-packet width conversion의 과도한 전역 barrier
4. GLB fill에 남아 있는 entity/datatype 집계 순서 문제
5. 비대칭 multi-chip 회귀와 NoP contention/back-pressure 부재
6. 실제 per-tile event timeline 부재
7. format axis의 network 출력/검증 누락
8. 상세 DRAM, PE reduction, systolic 모델 및 calibration 부족

우선순위별 현재 판정은 다음과 같다.

| 우선순위 | 판정 | 요약 |
|---|---|---|
| P0 | 해결 | network axis, 공식 metric 직접 검증, identity gate가 동작함 |
| P1 | 미완료 | pass-through는 해결됐지만 CE4 fill, GLB bypass, single-buffer, width conversion 문제 잔존 |
| P2 | 부분 해결 | 2-chip fixture와 source-read sharing 추가. contention/back-pressure 없음 |
| P3 | 부분 해결 | 평균 tile 비용 기반 fill+bottleneck만 존재. per-tile event timeline 없음 |
| P4 | 부분 해결 | 여러 설정 hook이 추가됐지만 물리 모델과 calibration은 부족함 |
| Deferred | 변화 없음 | sparse timing 및 외부 데이터 의존 항목 |

## 해결된 항목

### 1. Network timing 총계 axis

`stats_t::update_network_stats()`가 다음 값을 레이어별로 합산한다.

- `stage_axis_compute`
- `stage_axis_access[5]`
- `stage_axis_link[5]`
- `stage_axis_overlap[5]`

구현 위치: [scheduler/stats.cc](../scheduler/stats.cc), `update_network_stats()`.

Eyeriss AlexNet 결과를 독립적으로 다시 합산한 결과는 다음과 같다.

| 검산 항목 | 레이어 합 | network 출력 | 판정 |
|---|---:|---:|---|
| Compute-schedule latency | 23,999,585.4 | 23,999,585.3 | 출력 반올림 범위 내 일치 |
| DRAM busy | 568,047,803.0 | 568,047,803.0 | 일치 |
| Multi-chip busy | 281,715,200.0 | 281,715,200.0 | 일치 |
| Global buffer busy | 762,330,560.0 | 762,330,560.0 | 일치 |
| PE-array busy | 291,982,800.0 | 291,982,800.0 | 일치 |
| PE busy | 414,594,560.0 | 414,594,560.0 | 일치 |

각 stage의 access/link/overlap axis도 network 출력과 레이어 합이 일치했다. 따라서 기존 CE2/CE7의 network zero-axis 문제는 해결된 것으로 판정한다.

### 2. 공식 metric과 golden checker 연결

[validation/check_timing.py](../validation/check_timing.py)는 golden 비교 시 `Computation cycle + Fold fill cycle`을 자체적으로 재구성하지 않고 `Compute-schedule latency`를 직접 사용한다. 동시에 다음 identity도 layer와 network에서 검사한다.

```text
Compute-schedule latency == Computation cycle + Fold fill cycle
```

`stats_t::update_network_stats()`에도 동일한 identity에 대한 런타임 검사가 추가됐다. 기존 CE3/P0 검증 계약 문제는 해결된 것으로 판정한다.

### 3. Pass-through PE-array destination stage

[components/global_buffer.cc](../components/global_buffer.cc)의 `account_descriptor_dense_transfer()`는 `pe_array->exist_temporal_buffer == false`일 때 다음 항목을 overlap pipeline에서 제거한다.

- destination transaction count
- PE-array write-cycle stage

존재하지 않는 temporal-buffer write stage가 `cycle_pe_array_global_buffer`에 남던 문제는 해결됐다.

### 4. Entity-local 집계의 기본 경로

[scheduler/stats.cc](../scheduler/stats.cc)는 PE 또는 chip 내부에서 datatype을 먼저 결합한 뒤 entity 간 최대값을 선택하도록 별도 `entity_combined_*` 값을 계산한다.

```text
max_entity(combine_type(cycle[entity][type]))
```

PE local-buffer/link, PE-array fabric, GLB base access/link/overlap에는 올바른 집계 순서가 적용됐다. 다만 GLB fill은 아래의 미해결 항목과 같이 예외로 남아 있으므로 CE4 전체는 부분 해결이다.

### 5. 2-chip 실행 fixture

다음 파일이 추가됐다.

- [configs/accelerators/gemmini_2chip.cfg](../configs/accelerators/gemmini_2chip.cfg)
- [unittest/run_multi_chip_validation.sh](../unittest/run_multi_chip_validation.sh)
- [validation/multi_chip/check.py](../validation/multi_chip/check.py)

현재 fixture에서 다음 항목은 통과한다.

- `MIN <= AVG <= MAX`
- 전체 computation count
- Multi-chip access/link axis non-zero
- layer/network compute-schedule identity
- network busy-cycle axis populated

### 6. 기타 개선

- NoP shared tile의 source read를 distinct chunk 수에 따라 공유함
- DRAM tRAS/tRP 및 bank 수 설정 hook 추가
- PE lane-reduction energy hook 추가
- systolic output writeback에 MESH topology 적용
- format-IP cycle을 layer PE busy 계산에 포함
- temporal stage run에 fill+bottleneck helper 추가

이 항목들은 구현 진전은 있으나 아래 한계 때문에 각각 P2~P4 부분 해결로 분류한다.

## 아직 해결되지 않은 핵심 문제

### 1. Single-buffer LB cycle 중복 계상

심각도: **높음 / P1**

[components/pe.cc](../components/pe.cc)의 `account_descriptor_dense_mac_transfer()`는 동일 전송에 대해 다음 값을 모두 누적한다.

- `access_cycle_lb`: source 또는 destination LB access
- `transfer_cycle`: link serialization
- `cycle_mac_lb`: access/link/destination을 이미 결합한 pipeline makespan

그런데 [scheduler/stats.cc](../scheduler/stats.cc)의 single-buffer PE busy 계산은 다음처럼 세 값을 다시 모두 합한다.

```text
compute + access + link + overlap + format
```

`overlap`에 access/link 작업이 이미 들어 있으므로 동일 작업이 중복 계상된다. Single-buffer 직렬화는 일반적으로 `compute + transfer_makespan`처럼 중복되지 않는 축으로 재구성해야 한다.

추가 문제:

- `pe_t::double_buffer` 기본값이 `true`이다.
- 기존 accelerator config에는 PE-level `double_buffer` 설정이 없다.
- local-buffer capacity 검사는 double buffer일 때 두 복사본을 요구하지 않는다.
- single/double A/B cycle delta를 검증하는 실행 회귀가 없다.

따라서 LB7은 해결되지 않은 것으로 판정한다.

### 2. GLB bypass가 end-to-end direct stream이 아님

심각도: **높음 / P1**

[components/global_buffer.cc](../components/global_buffer.cc)는 bypass datatype에 대해 GLB→PE-array descriptor accounting을 즉시 반환한다. 그러나 upstream [components/multi_chip.cc](../components/multi_chip.cc)는 여전히 해당 datatype을 GLB destination에 기록하고 `fill_access_cycle/energy`를 과금한다.

Eyeriss conv1 weight의 실제 결과:

```text
GLB -> PE-array weight transactions : 0
GLB weight access cycle             : 69,696
```

현재 의미론은 다음과 같다.

```text
multi-chip -> GLB fill은 유지
GLB -> PE-array link/access는 제거
```

이는 “GLB를 거치지 않고 source에서 PE-array로 직접 stream”하는 end-to-end bypass와 다르다. 다음 중 하나를 명시적으로 구현해야 한다.

1. 실제 direct stream 경로와 source/link/destination timing을 별도로 계상
2. GLB storage만 bypass하고 physical link/interface는 유지하는 별도 계약 정의
3. 지원하지 않는 경우 config의 bypass를 거부

현재 traffic checker는 bypass weight의 T4/T5를 면제하기 때문에 이 불일치를 검출하지 못한다.

### 3. Multi-packet width conversion dependency

심각도: **높음 / P1**

[utils/interconnect_timing.cc](../utils/interconnect_timing.cc)의 count-aware helper는 downstream transaction 수가 감소하면 해당 경계를 전체 stream에 대한 barrier로 처리한다.

이 방식은 `(32 source, 1 link, 32 destination)` 같은 single-packet 예제와 tail micro-test는 해결한다. 그러나 `(64, 2, 64)`처럼 여러 packet이 있는 경우 다음 overlap을 허용하지 않는다.

```text
packet 1 link/destination 처리
    ||
packet 2 source assembly
```

즉 모든 source transaction이 끝난 후에야 모든 downstream 처리를 시작하는 과대평가가 발생한다. Aggregate count만 받는 현재 API로는 packet별 dependency와 tail group을 정확히 표현하기 어렵다.

필요한 수정:

- source/link/destination line width 또는 group ratio를 API에 전달
- packet별 assemble → transfer → disassemble dependency 구성
- 1 packet, 2 packet, tail packet 수기 makespan 테스트 추가

CE5/SP3/LB5/GB7/MC5는 부분 해결로 유지한다.

### 4. CE4의 GLB fill 집계 순서

심각도: **높음 / P1**

GLB base access는 entity-local datatype 결합 후 chip 간 max를 선택한다. 그러나 GLB fill은 아직 다음 순서로 집계된다.

```text
fill_access_cycle_global_buffer[type] = max_chip(fill[chip][type])
stage_fill = combine_type(fill_access_cycle_global_buffer[type])
```

Shared GLB에서 필요한 순서는 다음과 같다.

```text
max_chip(sum_type(fill[chip][type]))
```

예를 들어 chip 0의 input fill이 100이고 chip 1의 weight fill이 100이면 현재 식은 200을 만들 수 있지만 실제 병렬 elapsed는 100이다. CE4는 base 경로만 해결됐으며 fill 경로는 부분 해결이다.

### 5. Multi-chip 회귀의 검출력 부족

심각도: **중간 / P2**

현재 2-chip fixture는 모든 PE에서 computation cycle이 동일하다.

```text
MIN = AVG = MAX = 512
```

이 대칭 fixture에서는 다음 회귀를 놓칠 수 있다.

- chip loop마다 `first_active_pe`를 잘못 초기화하는 CE6 오류
- chip 0은 input-bound, chip 1은 weight-bound인 CE4 비대칭 집계 오류
- broadcast와 partitioned datatype의 source/link 차이

최소한 두 chip의 datatype별 cycle 분포가 서로 다른 비대칭 fixture와 수기 expected value가 필요하다.

### 6. NoP multicast, contention, back-pressure

심각도: **중간 / P2**

[components/multi_chip.cc](../components/multi_chip.cc)는 shared tile의 source read만 한 번 과금한다. Link transaction과 destination write는 chip별로 계속 합산하는 routed-unicast 모델이다.

남은 기능:

- multicast tree 또는 shared-link 전달 계약
- 동일 physical link를 공유하는 route 간 contention
- arbitration/queue
- destination-buffer pressure
- producer/consumer stall propagation

따라서 “source-read sharing”은 구현됐지만 multi-chip multicast timing 전체가 해결된 것은 아니다.

### 7. Per-tile global timeline 부재

심각도: **중간 / P3**

새 `temporal_pipeline_run_cycles()`는 stage 전체 비용을 repetition 수로 나눠 평균 per-tile 비용을 만든 후 다음 식을 사용한다.

```text
fill(sum of average per-tile stage costs)
+ (tiles - 1) * bottleneck average per-tile cost
```

이전의 단순 flat `max()`보다 fill 비용을 더 잘 반영하지만 실제 per-tile timeline은 아니다.

표현하지 못하는 항목:

- tile별 비용 차이와 tail tile
- `ready_cycle` / `available_cycle`
- buffer depth에 따른 overlap 허용량
- producer stall과 consumer back-pressure
- DRAM/NoP부터 PE까지 이어지는 end-to-end interval

특히 DRAM/multi-chip을 포함하는 overlap run은 repetition rate 차이 때문에 여전히 기존 flat `max()`를 사용한다. P3는 부분 해결로 유지한다.

### 8. Format axis의 network 계약 누락

심각도: **중간 / P4**

`stage_axis_format`은 layer의 `busy_cycle_pe`에 반영되지만 다음 위치에는 포함되지 않는다.

- `update_network_stats()` 합산
- `Busy-cycle axes (access / link / overlap)` 출력
- traffic T7의 busy-axis invariant

현재 accelerator config에서 format cycle 단가가 0이므로 기존 회귀는 통과한다. Format 단가를 활성화하면 network의 printed axes만으로 PE busy를 재구성할 수 없고 T7 계약도 불완전해진다.

### 9. 상세 DRAM 모델

심각도: **낮음~중간 / P4**

다음 설정은 추가됐다.

- `t_ras_cycle`
- `t_rp_cycle`
- `num_banks`

그러나 현재 analytical 모델은 row miss를 bank에 round-robin으로 균등 분배한다고 가정한다. Address 기반 bank mapping, row-buffer state, bank conflict, random/strided request는 없다. 또한 현재 accelerator config는 새 설정을 사용하지 않으므로 기본 실행 결과에는 변화가 없다.

DR2는 hook 확장으로는 진전됐지만 상세 timing 해결로 볼 수 없다.

### 10. PE lane reduction energy

심각도: **중간 / P4**

`lane_reduction_energy(mac_width, unit_energy)`는 한 computation issue마다 `mac_width - 1`개의 addition energy를 과금한다. 하지만 다음 상태를 사용하지 않는다.

- `active_accumulator_units`
- 마지막 accumulator의 `active_mac_width`
- 여러 independent accumulator unit

예를 들어 `number_of_macs=8`, `mac_width=8`인 구조에서 여러 accumulator가 동시에 활성화돼도 현재는 한 개 8-lane tree의 reduction energy만 과금한다. 반대로 8개보다 적은 lane만 활성화된 경우에도 full-width tree 비용을 과금한다.

현재 config에는 `mac_reduction_energy` calibration도 없으므로 단가는 기본 0이다. PE1/PE2는 부분 해결이다.

### 11. Systolic 세부 pipeline

심각도: **낮음~중간 / P4**

Output writeback topology가 MESH로 일관되게 적용됐지만 다음은 여전히 모델링되지 않는다.

- operand skew
- per-PE pipeline register
- partial-sum propagation
- array drain
- tile-boundary bubble

SY2의 analytical diameter-fill 범위는 개선됐지만 detailed systolic pipeline은 미완료다.

## 검증 실행 결과

### 실행 명령

```bash
bash unittest/run_validation.sh
bash unittest/run_timing_validation.sh
python3 validation/traffic/check.py
bash unittest/run_multi_chip_validation.sh
bash unittest/run_full_sanitizers.sh
python3 validation/check_timing.py --check-baseline --check-traffic
```

### 결과

| 검증 | 결과 |
|---|---|
| Unit/config validation | 통과, 205개 config |
| Gemmini RTL timing | 통과, MAPE 4.40%, max 7.86% |
| Eyeriss silicon latency | 통과, MAPE 4.26%, max 6.39% |
| Traffic T1~T8 | 8개 workload 통과 |
| 2-chip MC1~MC5 | 통과 |
| ASan/UBSan smoke | 통과 |
| Network axis 독립 합산 | 통과 |
| Eyeriss silicon traffic strict gate | 실패, max 133.94% > 15% |
| AlexNet 전체 layer timing | 부분 결과, 5개 non-convolution/connected layer 제외 |

기본 `run_timing_validation.sh`에서 Eyeriss traffic 오차는 informational이므로 timing gate 자체는 성공한다. `--check-traffic`을 활성화하면 exit code 1로 실패한다.

Traffic T1~T8 통과는 simulator 내부 transaction identity가 일관되다는 의미다. Eyeriss silicon traffic strict gate 실패는 외부 traffic reference와의 절대 오차가 아직 크다는 뜻이므로 두 결과는 모순되지 않는다.

## 수정 우선순위

### P1-A — Single-buffer cycle 재구성

1. `cycle_mac_lb`와 access/link axes의 포함 관계를 명문화한다.
2. single-buffer busy가 같은 transaction을 한 번만 과금하도록 수정한다.
3. `double_buffer=0/1` A/B fixture를 추가한다.
4. double-buffer LB capacity에 두 복사본을 요구할지 계약을 확정한다.

### P1-B — GLB bypass end-to-end 계약

1. bypass datatype의 실제 producer와 consumer를 명시한다.
2. multi-chip→GLB fill과 GLB→PE-array link 중 어떤 물리 자원이 남는지 결정한다.
3. direct-stream descriptor 경로를 별도로 구현한다.
4. bypass on/off의 access/link/critical-path delta를 수기 expected 값으로 검증한다.

### P1-C — Packet-group width conversion

1. aggregate count API를 packet/group-aware API로 교체한다.
2. assemble, transfer, disassemble dependency를 packet별로 계산한다.
3. `(32,1,32)`, `(64,2,64)`, tail packet을 모두 회귀에 포함한다.

### P1-D — CE4 GLB fill

1. chip별 fill datatype 결합값을 update 시점에 보존한다.
2. repetition scaling 후에도 `max_chip(combine_type(...))` 순서를 유지한다.
3. 비대칭 2-chip 수기 fixture로 검증한다.

### P2 — Multi-chip 검증과 contention

1. 비대칭 2-chip CE4/CE6 fixture 추가
2. broadcast/partitioned datatype별 source/link expected value 추가
3. multicast link 계약 정의
4. shared-link contention/arbitration 구현
5. destination pressure와 stall propagation 연결

### P3 — 실제 per-tile timeline

1. tile별 `ready_cycle`과 `available_cycle` 도입
2. producer/consumer interval 기록
3. buffer depth에 따른 overlap 제한
4. tail tile 및 per-tile 비용 차이 반영
5. 새 wall-clock에 leakage와 MAC availability 재연결

### P4 — Detailed model 및 calibration

1. format axis를 network/report/checker에 연결
2. active accumulator/lane 기반 reduction energy 수정
3. address 기반 DRAM bank/row state 또는 DRAMSim3 경로 확정
4. systolic skew/drain/bubble 구현
5. component 단가를 측정 또는 문헌 기반으로 calibration

## 완료 판정 기준

Timing simulation 완료를 선언하려면 최소한 다음 조건이 모두 충족돼야 한다.

1. layer/network 공식 metric과 모든 busy axis가 자기 일관적이다.
2. single/double/bypass/pass-through 설정이 수기 예측 가능한 cycle delta를 만든다.
3. 비대칭 multi-chip에서 entity/datatype 집계 순서를 검증한다.
4. equal-width 및 multi-packet width conversion이 packet-level 수기 계산과 일치한다.
5. NoP shared-link contention과 buffer pressure가 timeline에 반영된다.
6. per-tile timeline 또는 명시적으로 제한된 analytical 계약이 fill/drain/back-pressure를 일관되게 정의한다.
7. format, DRAM, PE reduction, systolic 상세 모델이 출력과 검증 계약에 연결된다.
8. 지원하지 않는 layer는 partial result로 명시되고, strict gate 범위가 문서화된다.

## 종합 의견

이번 변경으로 P0 결과 계약과 dense single-chip compute-schedule 검증은 안정화됐다. 이는 중요한 진전이다. 다만 현재 남은 P1 문제들은 단순한 모델 정밀도 문제가 아니라 동일 transfer의 중복 계상, 불완전한 bypass 경로, multi-packet 직렬화처럼 critical-path cycle을 직접 왜곡할 수 있는 문제다.

따라서 다음 구현 순서는 아래와 같이 유지한다.

```text
single-buffer 중복 계상
→ GLB bypass end-to-end 경로
→ packet-group width conversion
→ CE4 GLB fill 집계
→ 비대칭 multi-chip 회귀
→ contention/back-pressure 및 per-tile timeline
→ 상세 모델/calibration
```
