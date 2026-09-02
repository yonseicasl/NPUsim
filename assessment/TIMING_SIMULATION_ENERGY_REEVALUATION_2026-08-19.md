# Timing simulation energy 재평가 및 잔여 이슈

- 평가일: 2026-08-19
- 평가 기준: 현재 working tree의 dense timing simulation 경로
- 기준 문서: [ENERGY_POWER_ISSUES_2026-08-19.md](ENERGY_POWER_ISSUES_2026-08-19.md)
- 구현 기록: [../implementation/energy/energy_power_issues_2026-08-19_fixes.md](../implementation/energy/energy_power_issues_2026-08-19_fixes.md)
- 평가 범위: 기존 assessment의 Phase 1, 2, 3, 5 구현과 Phase 4 제외 상태 재검토

## 최종 판정

현재 버전은 이전 평가보다 크게 개선됐다. Component별 dynamic/static subtotal, layer/network
합산, NoC/NoP energy, precision별 MAC 단가, reduction/fold/setup/leakage 및 average power 산식에
대한 내부 회귀가 추가됐고 관련 검사는 모두 통과한다.

그러나 **Phase 4를 제외한 모든 energy 문제가 완료된 상태는 아니다.** 특히 다음 두 항목은
결과의 수치 자체를 잘못 해석하게 만들 수 있어 우선 수정해야 한다.

1. Accumulator spill과 output cast가 실제 event가 아니라 output 요청마다 함께 과금된다.
2. 보정되지 않은 기본 config도 absolute pJ로 인정되어 average power가 출력된다.

따라서 현재 지원 수준은 다음과 같다.

| 사용 목적 | 판정 |
|---|---|
| Dense transaction/event energy 내부 비교 | 사용 가능 |
| 기존 component별 상대 energy 비교 | 조건부 사용 가능 |
| Accumulator/format energy 비교 | 미완료 |
| Precision별 MAC energy 문법 및 scaling | 부분 완료 |
| Absolute pJ 정확도 | 검증되지 않음 |
| Average power 산술식 | 구현 완료 |
| Average power 결과의 calibration 자격 | 미완료 |
| Peak/chip/system power | 미지원 |
| Sparse/compressed energy | 미지원 |

## 검증 실행 결과

현재 바이너리로 다음 검증을 다시 실행했다.

- Unit/config validation: 267개 config 통과
- Energy, NoC, precision, accumulator, fold/setup, reduction, NoP, leakage, power 검사 통과
- DRAM 및 traffic 내부 항등식 통과
- Gemmini RTL latency: MAPE 4.40%, 최대 오차 7.86%
- Eyeriss silicon latency: MAPE 4.26%, 최대 오차 6.39%
- SCALE-Sim cross-simulator gate 통과

이 통과 결과는 구현이 정의한 내부 항등식이 유지된다는 뜻이다. 실제 event 의미가 올바른지와
absolute unit cost가 보정됐는지는 별도 문제이며, 아래 잔여 이슈가 바로 그 공백이다.

## 잔여 이슈

### RE1. Accumulator spill/output cast의 event 의미가 잘못됨

- 우선순위: **P1 / 가장 먼저 수정**
- 기존 issue: E4
- 상태: **부분 해결**

`accumulator_format`은 이제 spill byte 계산에 연결되어 FP32와 FP16 결과가 달라진다. 하지만
`account_format_events(OUTPUT, ...)`가 실제 output reload 여부를 확인하기 전에 호출되고, 이 함수가
매 호출마다 다음 두 event를 동시에 과금한다.

```text
accumulator spill/reload bytes
final output cast/pack bytes
```

관련 경로:

- `components/pe.cc:804-815`
- `components/pe.cc:1810-1943`
- `components/pe.cc:2723-2752`

현재 GEMM 64x64x64 fixture의 수기 대조 결과는 다음과 같다.

| 항목 | 현재 결과 | 실제 의미상 기준 |
|---|---:|---:|
| MAC 연산 | 262,144 | 64x64x64 |
| 최종 output element | 4,096 | 64x64 |
| INT8 output cast bytes | 262,144 | 최종 cast가 1회라면 4,096 bytes |
| FP32 accumulator spill bytes | 1,048,576 | 실제 eviction/reload boundary에서 산출해야 함 |

현재 output cast는 최종 output 수가 아니라 사실상 MAC 수를 따라 64배 과금된다. Accumulator
spill도 실제 partial-sum eviction/reload event가 아닌 output request 횟수를 따라간다. 특히
`edge_accumulation=1`에서는 accumulator가 array edge에 존재한다고 정의하면서 energy event는 PE의
output 요청 경로에서 생성되므로 component ownership도 일치하지 않는다.

현재 accumulator checker는 다음 비율만 검증한다.

```text
fp32 spill bytes / fp16 spill bytes = 2
spill bytes / cast bytes = accumulator width / output width
```

따라서 두 잘못된 event count가 같은 호출에서 만들어져도 검사가 통과한다.

완료 조건:

1. `output accumulator create`, `partial-sum spill`, `partial-sum reload`, `final cast/pack` event를
   서로 분리한다.
2. 새로운 output을 0으로 초기화하는 경로에는 spill/reload energy를 과금하지 않는다.
3. Retained output에는 reload를 과금하지 않고 실제 eviction boundary에서만 spill을 과금한다.
4. Final cast는 최종 output writeback 시점에 output element당 한 번만 과금한다.
5. `edge_accumulation=1`이면 accumulator access/spill energy를 array-edge component에 귀속한다.
6. GEMM 64x64x64에서 INT8 final cast가 4,096 element를 따르는 수기 fixture를 추가한다.
7. C/K tiling을 바꾼 A/B에서 최종 cast는 불변이고 spill/reload만 retention boundary에 따라
   변하는지 검사한다.

### RE2. Uncalibrated config가 absolute pJ 및 power로 인정됨

- 우선순위: **P1**
- 기존 issue: E7, E8, E10
- 상태: **부분 해결**

`energy_unit=normalized`인 Eyeriss fixture에서는 power가 올바르게 차단된다. 그러나
`energy_units_t`의 기본 단위가 pJ이고 `is_absolute()`는 단위 enum만 확인한다.

관련 경로:

- `utils/energy_units.cc:5-40`
- `scheduler/stats.cc:1660-1713`

대부분의 accelerator config에는 `energy_unit`과 `energy_reference` 선언이 없다. 그 결과 동일한
리포트가 다음 내용을 동시에 출력한다.

```text
MAC energy basis : computation_energy (UNCALIBRATED ...)
Energy unit      : pJ (absolute; totals and power are meaningful)
provenance       : not declared (uncalibrated)
Average power    : <numeric mW>
```

이는 산술식의 오류가 아니라 **출력 자격 조건의 오류**다. `energy/time` 계산은 맞지만 입력 energy가
absolute pJ라는 근거가 없으므로 mW로 표시하면 안 된다. 현재 power regression도 명시적인 absolute
fixture가 아니라 기본값 pJ인 synthetic Gemmini fixture를 사용해 이 문제를 고정하지 못한다.

완료 조건:

1. `energy_unit` 미선언 상태를 pJ가 아닌 `unspecified`로 둔다.
2. Absolute power 출력에는 최소한 다음 조건을 모두 요구한다.
   - `energy_unit=pJ` 명시
   - non-empty `energy_reference`
   - power scope에 필요한 component cost 선언 완료
   - 현재 operand/accumulator precision의 compute cost가 calibrated 상태
3. 조건이 부족하면 energy total은 `uncalibrated/estimated`로 표시하고 power는 unsupported로 둔다.
4. `pJ + provenance 없음`, `pJ + precision fallback`, `normalized`, `fully calibrated pJ` 네 가지
   fixture를 추가한다.
5. Power checker의 absolute fixture가 `energy_unit=pJ`와 provenance를 명시하도록 변경한다.

### RE3. Layer setup event가 unit energy 값에 종속됨

- 우선순위: **P2**
- 기존 issue: E5
- 상태: **부분 해결**

현재 setup cycle은 존재하지만 setup event는 `u_layer_setup_energy > 0`일 때만 1로 증가한다.

관련 경로:

- `scheduler/stats.cc:792-802`

기본 Gemmini는 `layer_setup_cycle=2270`을 선언하지만 `layer_setup_energy`가 0이므로 다음처럼 출력된다.

```text
Layer-setup energy : 0.00 pJ over 0 setup event(s)
```

실제 의미는 `1 setup event, unit energy uncalibrated`다. 반대로 setup cycle이 0인데 energy만 양수이면
현재 구현은 setup event와 energy를 만들 수 있어 latency와 energy event source가 분리된다.

완료 조건:

1. Setup event count는 unit energy가 아니라 실제 layer setup 수행 여부에서 생성한다.
2. `layer_setup_cycle > 0, energy = 0`이면 `1 event / uncalibrated`로 출력한다.
3. Energy만 양수이고 setup event source가 없으면 config error 또는 명시적 event 모델을 요구한다.
4. Fold/setup checker가 uncalibrated setup의 event count도 검사하도록 확장한다.

### RE4. MAC precision key가 accumulator precision을 포함하지 않음

- 우선순위: **P2**
- 기존 issue: E3, E4
- 상태: **부분 해결**

현재 compute energy key는 다음 두 operand만 포함한다.

```text
mac_energy_<input format>_<weight format>
```

관련 경로:

- `components/pe.cc:401-420`

따라서 input/weight가 같으면 accumulator FP16과 FP32도 같은 MAC energy를 사용한다. 현재 precision
fixture는 INT4/INT8/FP16 operand width scaling은 검증하지만 accumulator datapath 차이는 검증하지
않는다.

완료 방향:

1. `mac_energy_<input>_<weight>_<accumulator>` key를 사용하거나,
2. multiplier와 accumulation energy를 분리해 accumulator precision cost를 별도 과금한다.

완료 조건:

- INT8 x INT8에서 FP16/FP32 accumulator A/B를 추가한다.
- Operand traffic은 불변이고 compute/reduction/accumulator energy 중 정의된 축만 변해야 한다.
- 지원하지 않는 조합은 scalar fallback을 사용하더라도 power-calibrated로 인정하지 않는다.

### RE5. Energy config 검증이 schema 수준이 아님

- 우선순위: **P2**
- 기존 issue: E8
- 상태: **부분 해결**

현재 validator는 key 이름에 `energy` 문자열이 포함된 기존 설정의 값만 검사한다.

관련 경로:

- `utils/energy_units.cc:46-84`

따라서 negative, NaN, non-numeric 값은 막지만 다음은 검출하지 못한다.

- energy key 오타
- component별 필수 key 누락
- datatype vector 길이 부족 또는 빈 중간 field
- 일부 cost만 선언된 partial calibration
- `0`, `not-modeled`, `uncalibrated`의 의미 차이

`Uncosted components` 역시 설정 선언 여부가 아니라 계산 결과가 0인지로 판단한다. 선언된 cost가
있지만 해당 layer에 activity가 없는 component도 `no energy cost declared`로 표시될 수 있고, 반대로
한 개 cost만 선언된 component는 나머지가 누락돼도 costed로 보일 수 있다.

완료 조건:

1. Component 종류와 활성 feature별 energy schema를 정의한다.
2. Scalar/vector key의 타입, 길이, 필수 여부를 파싱 단계에서 검증한다.
3. Unknown key를 거부하거나 warning으로 출력해 오타를 검출한다.
4. `modeled zero`, `uncalibrated`, `not modeled`, `no activity` 상태를 분리한다.
5. `Uncosted components`를 결과값이 아닌 schema/calibration 상태에서 계산한다.

### RE6. Mesh multicast link traversal의 topology 계약이 불완전함

- 우선순위: **P2**
- 기존 issue: E2
- 상태: **대폭 개선됐으나 부분 해결**

Multicast distribution과 output writeback-unicast를 서로 다른 multiplier로 분리한 방향은 올바르다.
그러나 multicast traversal을 항상 다음처럼 계산한다.

```text
link_traversals = active_height * active_width
```

관련 경로:

- `utils/interconnect_timing.cc:104-127`

이 값은 endpoint ingress link까지 한 개씩 포함한다는 계약에서는 가능하다. 반면 mesh 내부 router
link만 세는 spanning tree는 N개 endpoint에 대해 N-1개 edge가 필요하다. 현재 모델에는 source/root
위치, multicast route/tree, endpoint injection link 포함 여부가 없어 `N`이 generic physical-link
traversal인지 확인할 수 없다.

완료 조건:

1. `noc_energy`가 router-to-router link, endpoint injection/ejection 또는 둘 모두 중 무엇의 단가인지
   정의한다.
2. Multicast source 위치와 route/tree 정책을 정의한다.
3. 1x1, 1xN, Nx1, NxM에서 실제 edge set을 생성하거나 수기 edge 목록과 대조한다.
4. `N`과 `N-1` 중 선택한 계약을 report에 출력하고 fixture 설명이 같은 계약을 사용하게 한다.

### RE7. 결과 label과 specification unit 표기가 완전히 일관되지 않음

- 우선순위: **P3**
- 기존 issue: E1, E7
- 상태: **부분 해결**

Network rollup 값은 layer 합과 일치하지만 `network.txt`에서도 다음 label을 사용한다.

```text
Layer dynamic energy
Layer static energy
Layer total energy
```

관련 경로:

- `scheduler/stats.cc:1624-1629`

또한 result file의 energy unit은 config 선언을 따르지만 accelerator specification의 여러 unit-cost
출력은 여전히 `pJ`가 하드코딩되어 있다.

예시 경로:

- `components/pe.cc:2628-2643`
- `components/global_buffer.cc:1946-1950`
- `components/spatial_arch.cc:1549`
- `components/systolic_array.cc:914`
- `components/multi_chip.cc:1888`

완료 조건:

1. Network rollup에서는 `Network dynamic/static/total energy`를 출력한다.
2. Checker가 layer와 network label을 구분해 검사한다.
3. Specification과 result의 모든 energy unit label을 `energy_units()`에서 가져온다.
4. Normalized fixture의 stdout/specification에도 `pJ`가 남지 않는 회귀를 추가한다.

### RE8. Non-zero checker가 event 의미와 total containment를 모두 검증하지 않음

- 우선순위: **P2~P3**
- 기존 issue: E9
- 상태: **부분 해결**

현재 신규 checker는 unit-cost scaling과 일부 total delta를 잘 검증한다. 하지만 다음 공백이 있다.

- Accumulator checker는 precision 비율만 보고 최종 output 수와 event 수를 대조하지 않음
- Reduction checker는 reduction energy의 layer/network total containment를 확인하지 않음
- Energy summary checker의 주 fixture에서는 reduction/fold/setup cost가 0이라 해당 subtotal 누락을
  검출할 수 없음
- Power checker는 명시적 calibrated-pJ fixture가 아니라 default-pJ fixture를 absolute로 사용함
- Fold/setup checker는 uncalibrated setup event count를 확인하지 않음

완료 조건:

1. 각 non-zero knob에 `event count x unit cost = component subtotal delta`를 적용한다.
2. 같은 delta가 layer total과 network total에 정확히 한 번 포함되는지 검사한다.
3. Event count 자체를 workload/mapping에서 수기로 계산해 검증한다.
4. Ratio-only fixture에는 최소 한 개의 absolute-count identity를 함께 둔다.

## 의도적으로 남아 있는 지원 범위

### Phase 4 — 외부 traffic 및 absolute energy calibration

사용자의 이번 구현 범위에서 제외한 항목이다. 기존 결과는 그대로 유지된다.

| 항목 | 현재 상태 |
|---|---:|
| Eyeriss DRAM traffic MAPE | 22.95% |
| Eyeriss DRAM 최대 오차 | 50.00% |
| GLB reference band 내부 | 5개 CONV 중 1개 |
| Component별 measured/calibrated absolute energy reference | 없음 |

따라서 RE1~RE8을 해결해도 Phase 4 전에는 absolute pJ 정확도를 주장할 수 없다.

### Peak 및 system power

현재 구현은 average core-datapath power만 지원하고 다음을 명시적으로 제외한다.

- peak power
- DRAM background/refresh 및 I/O termination
- clock network
- controller/DMA 및 PHY

이는 현재 report에 명시되어 있으므로 silent 오류는 아니다. 다만 chip/system power 완료로 간주할
수는 없다.

### Sparse 및 unsupported layer

Sparse timing entry는 명시적으로 거부되며 pooling, activation, normalization 등 unsupported layer는
energy scope에서 제외된다. Dense 결과로 잘못 계산하지는 않지만 end-to-end network energy는 아니다.

## 해결된 것으로 유지할 항목

다음 항목은 현재 dense analytical 범위에서 다시 열 필요가 없다.

- PE-array/Multi-chip temporal-buffer access energy 출력
- DRAM link/row-activation energy 출력
- Component subtotal 및 layer/network 합산 산식
- Energy scope와 partial network 표시
- GLB fill/bypass energy 계약
- Dynamic energy와 overlap 분리
- NoP latency/aggregate energy 분리
- Reduction energy의 별도 축
- Final critical-path leakage rescale
- Normalized fixture의 power 차단
- Average power/EDP/ED2P 산술 identity

## 권장 해결 순서

### 1단계 — 실제 energy event 복구

1. RE1 accumulator spill/reload/final-cast event 분리
2. Absolute event-count fixture 추가
3. Component/layer/network total containment 검증

이 단계가 끝나기 전에는 accumulator 및 format energy 비교를 사용하지 않는 것이 안전하다.

### 2단계 — Absolute energy와 power 자격 조건 강화

1. RE2 `unspecified/normalized/pJ-calibrated` 상태 분리
2. Power 출력 gate 강화
3. Calibrated-pJ 전용 power fixture 추가

이 단계가 끝나기 전에는 현재 출력되는 mW 값을 검증된 absolute power로 사용하면 안 된다.

### 3단계 — Event/provenance 안전성

1. RE3 setup event source 수정
2. RE5 component energy schema 도입
3. RE8 non-zero total/event checker 강화

### 4단계 — Precision 및 topology 계약

1. RE4 accumulator precision을 compute-energy 계약에 포함
2. RE6 multicast source/tree/link 범위 확정

### 5단계 — Report 계약 정리

1. RE7 network label 분리
2. Specification의 hard-coded pJ 제거
3. Calibration/scope 상태를 machine-readable하게 출력

### 6단계 — Phase 4

Traffic source/retention 모델과 외부 absolute reference를 확보한 뒤에 수행한다.

## 완료 판정 기준

Phase 4를 제외한 energy 구현 완료를 선언하려면 최소한 다음 조건을 만족해야 한다.

1. Final cast count가 최종 output element count와 일치한다.
2. Spill/reload count가 실제 retention/eviction boundary에서 생성된다.
3. Accumulator precision이 storage와 compute energy 양쪽의 명시적 계약에 포함된다.
4. 모든 dynamic axis가 component, layer, network total에 정확히 한 번 포함된다.
5. Energy config가 component별 calibration 상태를 판별한다.
6. Uncalibrated/default config는 absolute pJ 또는 mW를 출력하지 않는다.
7. Fold/setup event count가 unit cost와 독립적이다.
8. Mesh multicast의 counted-link 범위와 route가 명시되고 수기 edge fixture로 검증된다.
9. Layer/network label과 모든 unit 표기가 일관된다.
10. 위 조건을 검증하는 absolute-count 및 negative regression이 존재한다.

## 종합 의견

현재 구현은 기존의 “energy가 일부 계산되지만 결과 계약과 검증이 부족한 상태”에서
“dense relative energy 분석을 수행할 수 있는 상태”로 발전했다. 그러나 accumulator event와
absolute-power gating은 단순한 calibration 문제가 아니라 현재 코드의 event 및 validity 계약
문제다. 이 두 항목을 먼저 해결한 뒤 config schema와 검증 공백을 닫는 것이 가장 안전하다.
