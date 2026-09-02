# NPUsim energy·power 정합성 재평가 및 잔여 이슈

- 평가일: 2026-08-19
- 평가 기준: 현재 working tree의 timing simulation dense 경로
- 대상: compute, register/local buffer, PE-array NoC/temporal buffer, GLB, NoP,
  DRAM, format/reduction 및 static energy
- 제외: functional 값 정확성, 아직 명시적으로 거부되는 sparse 실행, 실제 RTL의
  transistor-level power 분석
- 선행 문서:
  - [ENERGY_TRAFFIC_VALIDATION.md](ENERGY_TRAFFIC_VALIDATION.md)
  - [TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md](TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md)
  - [../issue/timing/static_energy.md](../issue/timing/static_energy.md)

## 최종 판정

현재 energy 모델은 다음 범위에서는 상당히 개선됐다.

- dense runtime datatype 기반의 storage/access/link transaction 계산
- MAC 수에 비례하는 기본 computation energy
- GLB base read와 off-chip fill write의 repetition scaling 분리
- GLB storage bypass 시 SRAM energy 제거 및 물리 link traffic 유지
- NoP routed-unicast/multicast의 aggregate link traversal energy
- 최종 critical-path latency를 사용하는 PE/PE-array/GLB/multi-chip leakage 시간창
- layer energy의 network 합산

현재 `validation/energy/check.py`의 E1~E5b도 통과한다.

```text
E1~E5b: ALL ENERGY CHECKS PASSED
```

그러나 이 결과는 전체 energy 정합성 완료를 의미하지 않는다. Checker가 검증하는 것은
Eyeriss의 지원 layer에서 ALU, RF, GIN, GLB, DRAM access의 일부 상대 항등식과 정성적 분해다.
현재 모델에는 다음 세 종류의 문제가 남아 있다.

1. **계산되지만 출력되지 않는 energy가 있어 결과 파일만으로 total energy를 만들 수 없다.**
2. **Mesh multicast, computation precision, accumulator format에는 실제 energy 산식 공백이 있다.**
3. **절대 pJ reference와 power 계약이 없으므로 absolute energy/power 정확도를 주장할 수 없다.**

따라서 현재 지원 수준은 다음처럼 구분해야 한다.

| 사용 목적 | 현재 판정 |
|---|---|
| Dense event/transaction energy accounting | 대체로 개선됨 |
| Eyeriss 상대 energy breakdown | 제한적으로 사용 가능 |
| GLB fill/bypass energy 비교 | 신뢰도 비교적 높음 |
| Mesh/systolic NoC energy | 미완료 |
| Precision별 전체 energy 비교 | 미완료 |
| Multi-chip absolute energy | 외부 검증 없음 |
| Absolute pJ | 검증되지 않음 |
| Average/peak power | 미지원 |
| Sparse/compressed energy | 미지원 |

## 실제 계산 및 결과 계약 이슈

### E1. 수집되지만 결과에 출력되지 않는 energy

- 심각도: **높음 / P1**
- 영향: layer/network total energy 재구성 불가, 활성화된 knob의 효과가 결과에서 소실

다음 energy는 component와 `stats_t`에 존재하고 network에도 합산되지만 결과 파일에
출력되지 않는다.

| 누락 항목 | 수집/합산 | 현재 출력 |
|---|---|---|
| PE-array temporal-buffer access energy | `access_energy_pe_array` | interconnection/static만 출력 |
| Multi-chip temporal-buffer access energy | `access_energy_multi_chip` | interconnection/static만 출력 |
| DRAM link 및 row-activation energy | `transfer_energy_dram` | DRAM access energy만 출력 |

특히 `dram_t::account_row_activations()`는 다음처럼 row activation energy를
`transfer_energy`에 누적한다.

```text
transfer_energy[type] += row_activations * row_miss_energy
```

`gemmini_dram_detail.cfg`는 실제로 `row_miss_energy = 20`을 활성화하지만 DRAM 결과에는
access energy만 출력된다. 따라서 knob는 내부 값은 바꾸면서 사용자에게는 보이지 않는다.

현재 출력에는 다음도 없다.

- layer dynamic energy 합계
- layer static energy 합계
- layer total energy
- network dynamic/static/total energy
- component subtotal identity

완료 조건:

1. 모든 수집 energy를 component별로 출력한다.
2. `component subtotal = access + link + compute/format/reduction + static`을 정의한다.
3. `layer total = sum(component subtotal)`을 출력하고 checker로 재계산한다.
4. `network total = sum(layer total)`을 checker로 고정한다.
5. 동일한 물리 event가 두 component subtotal에 중복 포함되지 않는지 source/destination/link
   경계를 문서화한다.

### E2. Mesh multicast NoC energy의 physical-link traversal 과소계상

- 심각도: **높음 / P1**
- 영향: systolic/spatial mesh의 operand distribution energy 과소계상
- 현재 fixture에서 검출되지 않는 이유: 주요 config의 `noc_energy = 0`

PE-array dense distribution은 공유 operand tile을 temporal buffer에서 한 번 읽고 fabric에서
multicast하는 계약을 사용한다. Source read와 serialized stream을 PE 수만큼 복제하지 않는 것은
올바른 방향이다.

그러나 mesh energy는 한 multicast stream에 `average Manhattan hops`만 곱한다.

```text
transfer energy
  = serialized transactions
  * noc_energy
  * average Manhattan hops
```

이 식은 한 목적지로 가는 평균 unicast에는 적합하지만, 모든 active PE에 fanout하는 multicast
tree가 사용한 전체 physical link 수는 아니다. 예를 들어 16x16 mesh의 현재 multiplier는 15지만,
256개 endpoint를 잇는 tree는 최소 255개의 내부 grid link를 사용한다. `noc_energy`가 코드 주석의
계약대로 per-hop/per-link energy라면 현 식은 multicast energy를 크게 과소계상한다.

반면 multi-chip NoP는 `bottleneck_link_tiles`와 `total_link_traversals`를 분리하여 latency와
energy를 다르게 계산한다. PE-array NoC에도 같은 구분이 필요하다.

완료 조건:

1. NoC cost에 latency용 route depth와 energy용 total physical-link traversal을 분리한다.
2. BUS, crossbar, routed unicast, mesh multicast의 energy 계약을 각각 정의한다.
3. 1x1, 1xN, Nx1, NxM active shape의 수기 tree-link fixture를 추가한다.
4. `noc_energy > 0`인 live A/B config에서 active shape 변화가 예상 energy delta를 만든다.
5. output writeback의 per-source unicast와 input/weight multicast를 같은 식으로 취급하지 않는다.

### E3. Computation energy가 runtime precision과 무관함

- 심각도: **높음 / P1**
- 영향: FP16/INT8/INT4/MXFP 비교에서 compute energy가 동일하게 계산됨

Runtime datatype은 hierarchy의 bit traffic과 capacity에는 연결됐다. 하지만 MAC energy는
format과 무관한 scalar `computation_energy` 하나를 active operation 수에 곱한다.

```text
computation_energy += active_scalar_MACs * u_computation_energy
```

따라서 같은 config에서 input/weight format만 바꾸면 memory energy는 변하지만 MAC energy는
변하지 않는다. INT4와 FP16의 연산 energy가 같아도 모델은 이를 이상으로 판단하지 않는다.

필요한 계약은 다음 중 하나다.

1. `(input format, weight format, accumulator format)`별 MAC energy table
2. 기준 precision energy와 bit-width/operation-kind scaling rule
3. 외부 estimator가 계산한 per-operation energy를 config에 명시하고 provenance를 출력

완료 조건:

- 최소 FP16, INT8, INT4에 대해 서로 다른 수기 expected energy fixture를 둔다.
- MAC count가 같은 precision A/B에서 memory transaction 차이와 compute unit-cost 차이를 각각
  독립적으로 검증한다.
- 지원하지 않는 format 조합은 임의 scalar를 재사용하지 않고 명시적으로 거부하거나
  `uncalibrated`로 출력한다.

### E4. `accumulator_format`이 energy/storage 경로에 연결되지 않음

- 심각도: **높음 / P1~P2**
- 영향: accumulator RF/SRAM, spill/reload, reduction energy의 precision 정합성 부재

`runtime_datatypes()`는 `accumulator_format`을 파싱하고 report에 이름을 출력한다. 그러나 현재
component 계산에서 이 format을 소비하는 경로가 없다. Accumulator spill event도 accumulator
format이 아니라 output datatype의 storage transaction 수를 사용한다.

그 결과 다음 설정은 출력 문자열만 바꾸고 energy를 바꾸지 않는다.

```text
accumulator_format = fp32
accumulator_format = fp16
```

완료 조건:

1. accumulator resident capacity와 read/write transaction을 accumulator format으로 계산한다.
2. 최종 output cast/pack과 accumulator spill/reload를 서로 다른 event로 분리한다.
3. reduction tree energy와 accumulator write energy를 별도 출력한다.
4. FP32 accumulator + INT8 output 같은 mixed-precision 수기 fixture를 추가한다.

### E5. Fold/setup/drain latency에 대응하는 dynamic energy 부재

- 심각도: **중간 / P2**
- 영향: calibrated compute schedule을 사용해도 해당 동작의 dynamic energy는 0

`weight_fold_fill_cycle`, `layer_setup_cycle`과 analytical systolic drain은 compute schedule에
실제 cycle을 추가한다. 이 cycle은 최종 latency를 늘리므로 leakage에는 간접 반영된다. 하지만
weight reload, control/config, pipeline register, drain/flush 자체의 dynamic energy 항목은 없다.

따라서 latency calibration이 정밀해질수록 energy 쪽에서는 동일 동작이 leakage 외에는 무료로
남는 비대칭이 생긴다.

완료 조건:

- `weight_fold_fill_energy`, `layer_setup_energy` 또는 event-count 기반 equivalent를 정의한다.
- calibrated latency bubble과 analytical bubble을 동시에 과금하지 않는 것처럼 energy도 동일한
  선택 규칙을 사용한다.
- fold 수와 setup 횟수에 대한 수기 identity를 추가한다.

## 검증 및 calibration 이슈

### E6. 외부 traffic 오차가 memory energy로 직접 전파됨

- 심각도: **높음 / P1 완료 blocker**
- 영향: DRAM/GLB absolute energy 정확도 주장 불가

현재 Eyeriss external traffic 비교 결과는 다음과 같다.

| 항목 | 현재 결과 |
|---|---:|
| DRAM traffic MAPE | 22.95% |
| DRAM 최대 오차 | 50.00% |
| GLB reference band 내부 | 5개 CONV 중 1개 |

Energy의 기본식은 다음과 같다.

```text
dynamic energy = event/transaction count * unit energy
```

따라서 event count가 실제 silicon traffic과 다르면 unit-energy 산식이 정확해도 memory energy
오차가 거의 선형으로 전달된다. 현재 GLB 차이는 batch filter reuse와 psum spill/reload 같은
누락 traffic source/retention 모델과 연결되어 있다.

완료 조건:

1. 내부 traffic identity와 외부 traffic accuracy를 별도 gate로 유지한다.
2. Eyeriss DRAM/GLB traffic 오차를 명시한 목표 범위 안으로 줄인다.
3. energy checker가 traffic gate 실패 상태에서 absolute memory energy를 `validated`로 표시하지
   않도록 한다.
4. 최소 한 accelerator에서 component별 measured energy 또는 calibrated estimator reference를
   확보한다.

### E7. Eyeriss energy fixture의 normalized unit과 `pJ` 표기 혼용

- 심각도: **중간~높음 / P1~P2**
- 영향: 상대 energy 결과가 absolute pJ로 오해될 수 있음

`eyeriss_energy.cfg`는 공개 자료의 정규화 접근 비용을 사용한다.

```text
MAC = 1
RF  = 1x
GIN = 2x
GLB = 6x
DRAM = 200x per 16b word
```

이는 상대 분석 fixture로는 유효하다. 그러나 stats 출력은 모든 값을 `pJ`로 표시한다. 실제
MAC=1 pJ calibration을 수행한 것이 아니라면 출력 단위와 fixture 의미가 다르다.

완료 방향:

- config에 `energy_unit = pJ|normalized`와 reference/provenance를 추가하거나,
- 모든 unit cost를 실제 pJ로 calibration한 뒤 pJ 출력만 허용한다.

Checker도 normalized fixture에서는 absolute total pJ를 검증했다고 표현하지 않아야 한다.

### E8. Energy unit-cost 설정이 누락/오타/음수를 안전하게 거부하지 않음

- 심각도: **중간 / P2**
- 영향: 잘못된 config가 0 또는 음의 energy를 정상 결과처럼 출력 가능

대부분의 dynamic energy key는 미설정 시 0으로 남는다. Config validation은 energy key의 존재,
단위, provenance, 음수 여부를 일반적으로 검사하지 않는다. Static energy와 DRAM row energy 등
일부만 음수 guard가 있다.

이 때문에 다음 문제가 가능하다.

- key 오타가 silent zero-cost component를 만듦
- energy 분석용 config인데도 NoC/NoP/static energy가 전부 0
- 음수 dynamic unit cost가 음의 total energy를 생성
- `0`이 실제 zero인지 아직 calibration되지 않은 값인지 구분 불가

완료 조건:

1. energy 분석 모드에서 필수 component unit cost를 선언한다.
2. 모든 dynamic/static unit cost에 finite/non-negative 검사를 적용한다.
3. `0`, `unknown`, `not-modeled`을 구분한다.
4. report에 config provenance와 energy calibration scope를 출력한다.

### E9. 주요 energy knob의 non-zero 회귀 부재

- 심각도: **중간 / P2**
- 영향: 구현은 존재하지만 scaling/aggregation/output regression이 없음

현재 non-zero live regression이 부족한 항목:

- `noc_energy` / `nop_energy`
- `mac_reduction_energy`
- `format_payload_energy` / `format_metadata_energy`
- `accumulator_spill_energy`
- PE/LB/PE-array/GLB/multi-chip static energy
- DRAM row-activation energy의 출력 및 total 반영
- multi-chip multicast/unicast energy delta

`validation/format/check.py`는 non-zero format **cycle**을 검증하지만 energy는 활성화하지 않는다.
`validation/dram/check.py`도 row/bank **timing** 수기식은 검사하지만 row energy identity는 검사하지
않는다. Static energy unit test는 `unit * elapsed` helper 한 건뿐이며 production scaling,
final-latency rescale, network 합산을 검증하지 않는다.

완료 조건:

- 각 knob에 대해 `unit cost x event count` 수기 identity를 하나 이상 둔다.
- unit cost만 N배 한 A/B에서 해당 energy만 N배가 되는지 검사한다.
- latency만 바꾼 A/B에서는 dynamic energy가 유지되고 leakage만 최종 latency에 비례하는지
  검사한다.

## Power 모델 공백

### E10. Total energy와 wall-clock을 잇는 power 계약 부재

- 심각도: **높음 / P1~P2**
- 현재 상태: **power 미지원**

현재 출력에는 다음 metric이 없다.

- average power
- peak power
- dynamic power와 leakage power
- component별 active/idle power
- clock-tree/control/DMA power
- EDP/ED2P

Static energy는 `pJ/cycle * critical-path cycles`로 계산되지만, frequency는 주로 bitwidth 유도와
specification 출력에 사용된다. Stats는 cycle을 실제 시간으로 변환하지 않는다. 현재 각 config는
component frequency를 동일하게 두므로 문제가 숨겨지지만, 서로 다른 clock domain을 허용하면
같은 elapsed cycle 수를 모든 component에 적용할 수 없다.

Power 완료에 필요한 최소 계약:

```text
layer_time_seconds = authoritative_cycles / authoritative_frequency
average_power      = total_energy / layer_time_seconds
```

그 전에 다음을 결정해야 한다.

1. 모든 stage가 같은 clock domain인지, 아니면 local cycle을 시간으로 정규화할지
2. `static_energy`가 pJ/cycle인지 leakage power인지
3. stall 중 clock-tree/idle dynamic energy를 어느 domain에 과금할지
4. peak power를 event concurrency에서 계산할지 단순 unsupported로 둘지

### E11. DRAM background/refresh 및 system power 부재

- 심각도: **중간 / P2~P3**

현재 DRAM energy는 read/write access, link transaction, optional row activation 정도다. 다음은
모델링되지 않는다.

- background/standby power
- refresh energy
- I/O termination 및 PHY energy
- command/address bus energy
- controller/scheduler energy
- read/write turnaround의 energy

또한 accelerator 쪽에도 clock network, controller, DMA, interconnect idle dynamic power가 없다.
따라서 현재 dynamic+leakage 합계를 추가하더라도 chip 또는 system total power와 바로 비교할 수
없다. 지원 범위를 core datapath energy로 제한해 출력하거나, 별도 background power model을
추가해야 한다.

## 지원 범위 이슈

### E12. Sparse 및 unsupported layer의 energy 부재

- 심각도: **중간 / P2**

Timing entry는 sparse/compressed 실행을 명시적으로 거부하므로 dense 비용으로 잘못 계산하는
오류는 막고 있다. 그러나 다음 energy는 지원되지 않는다.

- sparse value/index/pointer traffic
- zero-detection 및 gating
- encoder/decoder energy
- metadata cache/buffer energy
- irregular NoC/DRAM traffic

AlexNet의 pooling, activation, normalization 등 unsupported layer도 network energy에서 제외된다.
따라서 partial timing report와 동일하게 network energy도 실제 end-to-end network energy가 아니다.
향후 total energy를 출력할 때 반드시 `Timing/Energy scope: N of M layers`를 함께 표시해야 한다.

## 현재 해결된 것으로 유지할 항목

다음은 현재 dense analytical 범위에서 다시 열 필요가 없는 계약이다.

### R1. GLB fill energy repetition scaling

GLB read side는 uniform GLB repetition으로, off-chip fill write side는 datatype repetition으로
각각 scaling한 뒤 합쳐진다. Retained tensor fill을 uniform 반복해 출처 없는 write energy를
만들던 문제는 해결됐다.

### R2. GLB storage bypass energy

Bypassed datatype은 GLB SRAM read/write energy를 내지 않지만 NoP와 GLB-PE-array fabric energy는
계속 낸다. Bypass는 storage만 제거한다.

### R3. Dynamic energy는 overlap에 의해 사라지지 않음

Double buffering은 latency overlap을 바꾸지만 수행된 read/write/link event의 dynamic energy는
항상 과금한다. Latency와 energy의 max/sum 규칙을 혼용하지 않는다.

### R4. Leakage final-latency rescale

PE, PE-array temporal buffer, GLB, multi-chip leakage는 초기 시간창에서 계산된 뒤 최종
per-tile critical path가 달라지면 같은 비율로 재조정된다. Network에서는 layer leakage 합을
사용한다.

### R5. NoP aggregate energy 방향

Multi-chip NoP는 busiest physical link가 latency를 정하고, total link traversal이 energy를
정하는 계약을 사용한다. Routed-unicast와 multicast A/B의 latency/traffic 계약은 존재한다.
남은 것은 이를 non-zero energy fixture로 고정하는 일이다.

## 이번 검증 실행 결과

현재 바이너리로 다음을 실행했다.

```bash
./npusim.sh run eyeriss_energy alexnet silicon
python3 validation/energy/check.py
```

결과:

| 검사 | 결과 | 한계 |
|---|---|---|
| E1 ALU energy/MAC identity | 통과 | scalar unit cost identity |
| E2 RF access/MAC band | 통과 | 넓은 정성 범위 |
| E3 CONV RF:rest band | 통과 | absolute energy 검증 아님 |
| E4 FC DRAM dominance | 통과 | 정성적 dominance |
| E5a DRAM access energy/cycle | 통과 | logical access 단가비 |
| E5b GLB interconnection energy | 통과 | serialized transaction identity |

현재 Eyeriss energy breakdown은 CONV에서 RF:rest `3.46~4.21:1`, FC에서 DRAM share 약
`74%`를 보인다. 이는 상대 breakdown sanity check로는 유용하지만, E1~E5b에 포함되지 않은
component와 absolute pJ 정확도를 검증하지 않는다.

## 권장 해결 순서

### Phase 1 — Energy 결과 계약 완성

1. **E1 누락 energy 출력**
   - PE-array access, multi-chip access, DRAM transfer/row energy 출력
2. **Component subtotal과 layer/network total 추가**
   - dynamic/static/total 분리
   - layer total 및 network layer-sum identity
3. **Energy scope와 unit provenance 출력**
   - `pJ`와 `normalized` 구분
   - partial layer 범위 명시
4. **Energy checker 확장**
   - 모든 출력 항목을 다시 합산해 total identity 검증

이 단계가 먼저 필요한 이유는 이후 산식을 수정해도 결과에서 관측되지 않으면 회귀를 만들 수
없기 때문이다.

### Phase 2 — 현재 산식의 직접 결함 제거

5. **E2 Mesh multicast physical-link energy**
   - route depth와 total link traversal 분리
   - 1xN/NxM non-zero NoC fixture
6. **E3 precision-dependent MAC energy**
   - FP16/INT8/INT4 unit cost 또는 scaling rule
7. **E4 accumulator format 연결**
   - accumulator storage/spill/cast/reduction 분리
8. **E5 fold/setup/drain dynamic energy**
   - event count와 unit energy 추가

### Phase 3 — Non-zero regression과 config 안전성

9. **E9 knob별 micro/A-B gate**
   - NoC, NoP, reduction, format, spill, leakage, DRAM row
10. **E8 config validation**
    - finite/non-negative, required/unknown/zero 구분
11. **Leakage production-path gate**
    - final latency 비례, dynamic 불변, network layer 합
12. **Multi-chip energy gate**
    - multicast/unicast total traversal 수기식

### Phase 4 — 외부 energy calibration

13. **E6 traffic 정확도 개선**
    - batch filter reuse
    - psum spill/reload 및 GLB source 보강
14. **최소 한 accelerator의 absolute component energy reference 확보**
15. **precision별 compute/memory energy calibration**
16. **지원 범위를 벗어난 accelerator/config는 unvalidated로 표시**

### Phase 5 — Power 모델

17. **공통 시간 단위 및 clock-domain 계약 확정**
18. **average dynamic/leakage/total power 출력**
19. **DRAM background/refresh 및 clock/control power 범위 결정**
20. **EDP/ED2P와 power regression 추가**

Peak power는 실제 event concurrency가 필요하므로 per-tile/per-packet event timeline이 충분히
정밀해진 뒤 진행하는 것이 안전하다. 초기 power 단계에서는 average power만 지원하고 peak는
명시적 unsupported로 두는 편이 낫다.

## Energy 완료 판정 기준

다음 조건이 모두 충족되기 전에는 “energy simulation 완료”를 선언하지 않는다.

1. 모든 내부 energy accumulator가 결과에 출력되거나 명시적으로 diagnostic-only로 분류된다.
2. Component subtotal, layer total, network total이 checker로 재구성된다.
3. Dynamic energy는 event count와 unit cost의 수기식에 일치한다.
4. Static energy는 최종 wall-clock 및 physical always-on component 수와 일치한다.
5. Mesh multicast와 routed-unicast energy가 physical-link traversal에 일치한다.
6. Runtime operand/accumulator precision이 compute와 storage energy를 실제로 변화시킨다.
7. Non-zero format/reduction/spill/row/static knob이 모두 regression에 포함된다.
8. Energy 단위와 calibration provenance가 출력에 명시된다.
9. 최소 한 accelerator의 memory traffic과 component energy가 외부 reference로 검증된다.
10. Partial layer와 sparse 미지원 범위가 total energy 옆에 명시된다.

Power 완료는 여기에 추가로 다음을 요구한다.

11. Cycle을 seconds로 변환하는 authoritative clock 계약이 있다.
12. Average dynamic, leakage, total power가 energy/time identity와 일치한다.
13. 포함되지 않은 background/clock/PHY power가 명시된다.

## 종합 의견

Latency 작업으로 leakage 시간창을 신뢰할 수 있는 방향으로 정리한 것은 energy 단계의 좋은
기반이다. 또한 dense traffic과 GLB fill/bypass energy의 내부 항등식도 이미 상당 부분 확보됐다.

그러나 현재 가장 먼저 해결해야 할 문제는 더 정교한 energy 공식을 추가하는 것이 아니라,
**이미 계산되는 모든 energy를 빠짐없이 출력하고 total identity를 만드는 것**이다. 이 기반을
고정한 뒤 Mesh multicast, precision-dependent MAC/accumulator energy를 수정하고, 마지막으로
external traffic/energy calibration과 power 모델로 진행하는 순서가 가장 안전하다.

```text
energy observability + total identity
-> mesh multicast physical-link energy
-> precision/accumulator energy
-> non-zero knob regressions
-> external traffic/energy calibration
-> average power
-> detailed background/peak power
```
