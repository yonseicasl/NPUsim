# Special Function Unit(SFU) 추가 계획

## 목적

NPUSim에 activation을 처리하는 명시적인 Special Function Unit(SFU)을 추가한다. 목표는 activation의 functional 결과만 계산하는 것이 아니라, 최종 output 생성 시점에 발생하는 SFU 연산의 cycle, internal traffic, dynamic/static energy를 동일한 hardware event에서 일관되게 계산하는 것이다.

초기 구현은 convolution/connected layer에 결합된 element-wise activation을 대상으로 한다. standalone softmax와 shortcut/vector 연산은 서로 다른 실행·traffic 특성을 가지므로 후속 단계에서 별도로 확장한다.

우선 지원 순서는 다음과 같다.

1. LINEAR bypass
2. ReLU
3. Leaky ReLU
4. Sigmoid, Tanh, Hard Sigmoid, Hard Swish 등 element-wise nonlinear activation
5. Softmax와 같은 reduction 기반 복합 연산

현재 제공되는 network config에서는 ReLU, Leaky ReLU, Linear가 대부분을 차지하므로 1차 구현만으로 기존 convolution/connected workload의 activation 범위를 상당 부분 처리할 수 있다.

## 현재 상태와 문제

현재 [`components/pe.cc`](../components/pe.cc)의 `pe_t::activation()`은 모든 MAC register에 ReLU를 적용하도록 하드코딩되어 있다. 그러나 이 함수는 실제 accelerator 실행 경로에서 호출되지 않으므로 activation의 timing, traffic, energy는 계산되지 않는다.

반면 `FUNCTIONAL` 빌드에서는 [`scheduler/npu.cc`](../scheduler/npu.cc)가 각 Nebula layer의 `forward()`를 호출하고, Nebula가 layer에 지정된 activation을 수행한다. 따라서 새 SFU가 동일한 output에 activation을 다시 적용하면 functional 결과가 중복 변환된다.

현재 구조의 주요 문제는 다음과 같다.

- ReLU 이외의 activation 종류를 PE activation 함수가 표현하지 못한다.
- activation 연산 수와 invocation 수가 통계에 없다.
- activation throughput, pipeline latency, stall이 layer latency에 반영되지 않는다.
- SFU 내부 ingress/egress traffic과 energy가 없다.
- nonlinear activation이 실행되더라도 비용이 없다는 사실이 energy completeness에 표시되지 않는다.
- partial-sum spill/reload와 최종 output commit을 구분하지 않으면 activation이 reduction tile마다 중복 청구될 수 있다.
- mapping padding으로 계산된 가상 output과 실제 network output을 구분해야 한다.
- standalone softmax layer는 현재 convolution/connected timing 경로에서 제외된다.

## 설계 원칙과 범위

### 핵심 원칙

1. **최종 output당 한 번 실행**: activation은 partial sum마다 실행하지 않고 모든 reduction이 끝난 유효 output element마다 정확히 한 번 실행한다.
2. **functional owner 단일화**: 초기 구현에서는 Nebula가 functional activation의 유일한 실행 주체가 되고, SFU는 동일 이벤트의 hardware cost를 모델링한다.
3. **operation별 비용 분리**: ReLU와 Sigmoid/Tanh/Softmax에 하나의 공통 `activation_cycle` 또는 `activation_energy`를 적용하지 않는다.
4. **fused traffic 보존**: fused activation은 중간 tensor를 DRAM에 materialize하지 않으므로 기존 external DRAM traffic을 증가시키지 않는다.
5. **내부 traffic 명시**: SFU ingress/egress element·byte·transaction은 기존 hierarchy traffic과 분리해 보고한다.
6. **cycle과 energy의 동일 event source**: operation count, cycle, dynamic energy는 동일한 SFU invocation에서 생성한다.
7. **미지원·미가격 상태 표시**: 활성 nonlinear operation의 SFU 또는 energy cost가 없으면 silent zero가 아니라 `NOT_MODELED`, `PARTIAL`, `UNPRICED_ACTIVE` 상태를 보고한다.
8. **기존 정합성 보존**: `[sfu]`가 없는 기존 accelerator config의 수치 결과는 유지하되 activation이 모델링 범위 밖이라는 사실을 명시한다.

### 초기 지원 범위

| 구분 | 1차 지원 | 후속 지원 |
|---|---|---|
| layer | convolution, connected의 fused activation | shortcut activation, standalone layer |
| operation | Linear, ReLU, Leaky ReLU | ELU, Sigmoid, Tanh, Hard Sigmoid, Hard Swish 등 |
| 복합 연산 | 지원하지 않음 | Softmax multi-pass microprogram |
| execution | valid output element 기준 | masked/padded lane의 물리 energy option |
| traffic | fused internal traffic | materialized GLB/DRAM round trip |
| overlap | 직렬·보수적 모델 | queue/backpressure 기반 streaming overlap |

Pooling은 activation SFU 범위에 포함하지 않는다. Shortcut addition도 activation과 함께 실행될 수 있지만 본질적으로 vector ALU 연산이므로 SFU activation과 별도의 operation class로 다룬다.

## 목표 아키텍처

초기 SFU는 chip당 하나의 post-processing component로 두고, 최종 accumulation과 output cast 사이에 배치한다.

```text
MAC / reduction
        │
        ▼
final accumulator value
        │ final_output_tile event
        ▼
SFU activation
        │
        ▼
output format cast
        │
        ▼
GLB / NoP / DRAM output writeback
```

초기 의미론은 `activation-before-output-cast`로 고정한다. 향후 NVDLA SDP처럼 conversion, bias, scale, activation의 순서가 다른 구조를 지원할 때는 post-processing pipeline의 stage order를 config에 명시한다.

### 새 SFU component

다음 파일을 추가한다.

- `components/sfu.h`
- `components/sfu.cc`

`npu_t`는 기본적으로 chip당 SFU 하나를 소유한다.

```cpp
std::vector<sfu_t*> sfus;
```

`sfu_t`는 최소한 다음 상태와 API를 제공한다.

- 지원 activation operation 조회 및 validation
- layer/tile별 `final_output_tile` event 수신
- valid element 수, invocation 수, active chip 수 기록
- operation별 pipeline latency와 initiation interval 적용
- input/output byte 및 internal transaction 계산
- dynamic operation/read/write/setup energy 계산
- layer critical path를 위한 busy cycle과 stall 정보 제공
- layer duration에 따른 static energy 계산
- reset, specification 출력, stats export

### Activation operation descriptor

Nebula의 `activation_type_t`를 SFU가 직접 하드코딩된 switch로 반복하지 않도록 adapter 또는 descriptor를 둔다.

- operation name
- element-wise 여부
- reduction 필요 여부
- bypass 여부
- 지원 precision/profile
- pipeline latency
- initiation interval
- per-element dynamic energy
- setup cycle/energy

`LINEAR`는 bypass operation으로 선언하고 element, cycle, energy를 청구하지 않는다. Unknown operation은 ReLU로 fallback하지 않고 config validation 또는 layer 시작 시 실패시킨다.

## 최종 output 이벤트 계약

SFU 정확성에서 가장 중요한 부분은 activation 실행 시점이다. SFU를 기존 `pe_t::activation()`에 단순히 연결하면 PE tile과 reduction 순서에 따라 중복 실행될 수 있다.

Scheduler 또는 output commit 경로에 다음 의미의 event를 추가한다.

```text
final_output_tile {
    layer_index
    chip_index
    output_tile_identity
    valid_element_count
    activation_type
    input_format
    output_format
}
```

이 event는 다음 조건을 모두 만족할 때 한 번만 발생해야 한다.

- 해당 output tile의 C/R/S reduction이 모두 완료됨
- accumulator가 더 이상 partial sum으로 재사용되지 않음
- accumulator spill/reload 여부와 무관하게 최종 값이 확정됨
- mapping padding을 제외한 valid output 영역이 계산됨
- output cast 및 최종 writeback 이전임

`flush_psum_writeback()`만 event source로 사용하면 마지막 resident tile만 보이거나 중간 spill과 최종 commit을 구분하지 못할 수 있다. 따라서 scheduler가 output tile identity와 reduction completion을 추적하고, final commit 경로에서 SFU event를 발행해야 한다.

기본 activation element 수는 network의 유효한 `B × K × P × Q`로 계산한다. Mapping이 K/P/Q를 padding하더라도 padding output이 외부로 commit되지 않는다면 SFU에 청구하지 않는다. 향후 실제 hardware가 masked lane에도 dynamic energy를 소비한다는 근거가 있으면 `charge_masked_lanes`와 같은 별도 profile로 추가한다.

## 설정 형식

Accelerator config에 `[sfu]` section을 추가한다.

```ini
[sfu]
num_units_per_chip = 1
lanes = 16
queue_depth = 1
placement = post_accumulator
fusion = fused

relu_pipeline_latency = 1
relu_initiation_interval = 1

leaky_pipeline_latency = 2
leaky_initiation_interval = 1

sfu_op_energy_relu = 0.10
sfu_op_energy_leaky = 0.18
sfu_read_energy = 0.02
sfu_write_energy = 0.02
sfu_setup_energy = 0.00
sfu_static_energy = 0.001
```

초기 설정 계약은 다음과 같다.

- `num_units_per_chip`: chip당 병렬 SFU pipeline 수
- `lanes`: SFU pipeline 하나가 동시에 처리하는 element 수
- `queue_depth`: producer와 SFU 사이의 output tile queue depth
- `placement`: 초기에는 `post_accumulator`만 허용
- `fusion`: 초기에는 `fused`만 허용
- `<op>_pipeline_latency`: 첫 결과가 나올 때까지의 cycle
- `<op>_initiation_interval`: 다음 chunk를 투입할 수 있는 cycle 간격
- `sfu_op_energy_<op>`: operation별 유효 element dynamic energy
- `sfu_read_energy`, `sfu_write_energy`: SFU 내부 input/output element access energy
- `sfu_setup_energy`: invocation당 setup energy
- `sfu_static_energy`: 물리 SFU 하나의 pJ/cycle leakage

0 lane, 0 initiation interval, 음수 cycle/energy, 지원되지 않는 placement/fusion은 fail-fast한다. 활성 operation에 필요한 timing profile이 없으면 strict mode에서는 실패시키고 compatibility mode에서는 미모델링 상태를 명시한다.

## Cycle 모델

하나의 연속적인 SFU invocation에 대해 다음 기본식을 사용한다.

```text
chunks = ceil(valid_elements / lanes)

SFU cycles =
    setup_cycles
  + pipeline_latency
  + (chunks - 1) × initiation_interval
```

`valid_elements = 0` 또는 `LINEAR`이면 SFU cycle은 0이다. 여러 SFU unit이 있는 경우 output chunk를 균등 분배하고, chip 내부 latency는 unit별 busy cycle의 최댓값으로 계산한다.

여러 chip이 병렬로 실행되는 경우:

- SFU latency는 active chip별 SFU busy cycle의 최댓값
- SFU energy와 traffic은 active chip 전체의 합

으로 결합한다.

### Timeline 통합

초기 단계에서는 `queue_depth=1`인 직렬 SFU 모델을 사용해 과도한 overlap을 방지한다. 그러나 SFU cycle을 `finalize_layer_timeline()` 결과에 사후 덧셈만 해서는 안 된다. 그렇게 하면 다음 값이 서로 불일치한다.

- layer critical-path latency
- MAC available cycle과 utilization
- SFU 및 기존 component static energy window
- output producer stall

따라서 [`scheduler/stats.cc`](../scheduler/stats.cc)의 timeline에 SFU를 post-processing/output-path resource로 추가한다. 최종 모델에서는 각 output tile의 final-reduction 완료 시각을 SFU release time으로 사용하고, queue depth에 따라 producer backpressure를 계산한다.

1차 구현에서는 직렬 실행을 명시적으로 보고하고, 후속 단계에서 다음 streaming 모델을 추가한다.

```text
output producer
    → bounded queue
    → SFU pipeline
    → output cast/writeback consumer
```

SFU가 느리면 queue가 채워져 PE-array output drain이 stall하고, SFU가 충분히 빠르면 pipeline fill/drain을 제외한 대부분의 시간이 기존 계산과 겹쳐져야 한다.

## Traffic 모델

Fused mode에서는 activation 전후 tensor가 DRAM에 별도로 materialize되지 않는다. 따라서 다음 불변조건을 유지한다.

- input/weight DRAM traffic 변화 없음
- 최종 output DRAM payload 변화 없음
- activation 때문에 추가 DRAM read/write가 발생하지 않음
- 기존 final output hierarchy transaction을 중복 청구하지 않음

SFU는 다음 internal traffic을 별도로 기록한다.

- ingress elements/bytes/transactions
- egress elements/bytes/transactions
- input/output bitwidth
- tail lane utilization
- optional internal buffer access

SFU placement가 array edge인지 GLB output port인지에 따라 기존 PE-array↔GLB link ownership이 달라질 수 있다. 초기 `post_accumulator` placement에서는 기존 final output payload를 그대로 전달하면서 SFU 내부 read/write만 추가하고, 동일 payload의 hierarchy transfer를 한 번 더 청구하지 않는다.

향후 `fusion=materialized`를 지원할 때만 중간 output의 GLB 또는 DRAM write/read traffic을 명시적으로 추가한다.

## Energy 모델

SFU dynamic energy는 다음 식을 기본으로 한다.

```text
dynamic energy =
    valid_elements × operation_energy
  + input_elements × sfu_read_energy
  + output_elements × sfu_write_energy
  + invocations × sfu_setup_energy
```

Static energy는 물리 SFU 개수와 최종 layer critical-path latency를 기준으로 계산한다.

```text
static energy =
    physical_sfu_count × layer_latency × sfu_static_energy
```

Power gating을 모델링하기 전에는 config에 존재하는 모든 물리 SFU를 layer 동안 always-on으로 취급한다. 향후 clock/power gating을 추가할 때는 idle/static energy 정책을 별도 profile로 분리한다.

[`utils/energy_units.cc`](../utils/energy_units.cc)의 energy schema에 다음 key family와 SFU component row를 추가한다.

- `sfu_op_energy_*`
- `sfu_read_energy`
- `sfu_write_energy`
- `sfu_setup_energy`
- `sfu_static_energy`

활성 operation에 operation/read/write energy가 선언되지 않았다면 SFU subtotal을 0으로 보이지 않고 unpriced active event로 표시한다. 이러한 run은 absolute energy, average power, EDP qualification을 통과하지 못해야 한다.

## Functional 모델 정책

초기 구현에서는 Nebula layer의 `forward()`가 functional value 계산의 유일한 owner이다. SFU는 activation type과 element count를 읽어 cost event만 생성한다.

이 정책을 사용하는 이유는 현재 NPUSim PE 내부 결과가 Nebula network output의 유일한 source of truth가 아니기 때문이다. SFU가 PE buffer에 activation을 적용한 뒤 Nebula `forward()`가 다시 실행되면 activation이 중복 적용될 수 있다.

SFU component에는 unit test용 pure activation evaluator를 둘 수 있지만, main functional path에서는 호출하지 않는다. 향후 NPUSim component data가 network output을 직접 생산하도록 변경할 경우 다음 중 하나만 선택한다.

1. Nebula `forward()`를 convolution/activation으로 분리하고 activation을 SFU에 위임
2. Nebula를 reference oracle로만 실행하고 NPUSim output과 비교

어느 경우에도 동일 output을 두 경로가 함께 수정하지 않는다.

## 단계별 실행 계획

### Phase 0 — 계약 고정과 baseline 확보

1. activation별 지원 matrix와 SFU placement/cast order를 문서화한다.
2. `[sfu]`가 없는 Eyeriss/Gemmini/NVDLA 관련 baseline의 cycle, traffic, energy 결과를 보존한다.
3. convolution/connected layer별 activation type과 valid output element 수를 출력하는 진단 경로를 만든다.
4. 기존 `pe_t::activation()`이 호출되지 않는다는 regression 검사를 추가하고 제거/deprecated 여부를 결정한다.
5. nonlinear activation이 존재하지만 SFU가 없는 run의 report 정책을 고정한다.

완료 조건: 기존 수치 baseline이 고정되고, 모든 실행 layer의 activation type과 expected element 수가 추적 가능하다.

### Phase 1 — SFU component와 config parser

1. `components/sfu.{h,cc}`를 추가한다.
2. `npu_t`에 chip별 SFU 생성, 소멸, connect, reset 경로를 추가한다.
3. `[sfu]` section과 operation별 latency/II parser를 구현한다.
4. Linear/ReLU/Leaky ReLU 지원 및 unknown operation fail-fast를 구현한다.
5. SFU specification 출력에 unit, lane, throughput, supported operations를 표시한다.

완료 조건: SFU config의 유효·무효 조합을 unit test로 검증하고, `[sfu]`가 없는 기존 config가 실행된다.

### Phase 2 — 최종 output event와 정확한 연산 수

1. scheduler에 output tile identity와 reduction completion 상태를 추가한다.
2. final accumulator commit 시 `final_output_tile` event를 발생시킨다.
3. mapping padding을 제외한 tile별 valid output element를 계산한다.
4. SFU invocation, valid element, chunk, tail lane utilization을 기록한다.
5. accumulator spill/reload와 reduction tiling이 activation count를 증가시키지 않는지 검증한다.

완료 조건: 모든 layer에서 SFU valid element 합이 실제 output element 수와 일치하고, 각 output identity가 정확히 한 번 처리된다.

### Phase 3 — Traffic와 energy 통합

1. SFU ingress/egress byte와 transaction을 추가한다.
2. fused mode에서 기존 hierarchy/DRAM transaction이 증가하지 않도록 ownership을 정리한다.
3. operation/read/write/setup/static energy를 동일 invocation event에서 계산한다.
4. `stats_t`와 network rollup에 SFU component row를 추가한다.
5. energy schema와 unpriced active event 검사를 추가한다.

완료 조건: energy perturbation이 해당 SFU subtotal에만 선형적으로 반영되고, fused SFU enable 전후 DRAM traffic이 동일하다.

### Phase 4 — Latency와 critical path 통합

1. operation별 latency/II 공식을 구현한다.
2. chip 내부 parallel SFU unit과 chip 간 병렬 실행을 반영한다.
3. 우선 `queue_depth=1` 직렬 모델을 timeline에 추가한다.
4. layer latency 변경에 맞춰 모든 component leakage window와 MAC availability를 재계산한다.
5. SFU busy/stall/overlap 정보를 layer report에 출력한다.

완료 조건: hand calculation과 SFU cycle이 일치하며 SFU latency 변화가 critical path, static energy, utilization에 일관되게 반영된다.

### Phase 5 — Streaming overlap과 backpressure

1. tile별 final-output release time을 timeline에 전달한다.
2. bounded SFU queue와 output writeback consumer를 모델링한다.
3. queue depth 1/2/N에 따른 producer stall을 계산한다.
4. SFU가 bottleneck인 경우와 숨겨지는 경우를 분리해 보고한다.
5. output cast와 SFU의 순서 및 resource overlap을 명시적으로 처리한다.

완료 조건: queue depth와 SFU throughput perturbation이 예상된 bottleneck 전환을 만들고, stall attribution 합계가 layer latency와 일치한다.

### Phase 6 — Operation 확장

1. ELU, Sigmoid, Tanh, Loggy, Hard Sigmoid, Hard Swish 등 Nebula element-wise activation을 추가한다.
2. operation별 precision, approximation mode, LUT/polynomial profile을 정의한다.
3. functional evaluator를 Nebula reference와 비교한다.
4. GELU 등 Nebula frontend에 없는 operation은 frontend enum/semantics가 추가된 뒤 지원한다.

완료 조건: operation별 latency/energy profile이 독립적이고, 미지원 approximation이나 precision 조합은 fail-fast한다.

### Phase 7 — Standalone Softmax

Softmax는 element-wise activation으로 계산하지 않고 다음 micro-operation으로 분해한다.

```text
max reduction
→ subtract
→ exp/LUT
→ sum reduction
→ reciprocal
→ normalize
```

1. standalone softmax layer를 accelerator timing 대상으로 등록한다.
2. reduction pass, element-wise pass, scratchpad read/write를 분리한다.
3. tensor가 on-chip에 유지되는 경우와 DRAM materialization되는 경우를 구분한다.
4. batch와 softmax axis에 따른 invocation 수를 계산한다.
5. mapping section이 없는 vector/reduction layer를 위한 별도 scheduler contract를 만든다.

완료 조건: softmax의 multi-pass cycle·traffic·energy가 hand calculation과 일치하고 convolution/connected mapping index를 깨뜨리지 않는다.

### Phase 8 — RTL 및 silicon calibration

1. NVDLA SDP의 activation/conversion datapath를 우선 calibration 대상으로 사용한다.
2. RTL performance counter 또는 trace에서 SFU input/output element, stall, busy cycle을 추출한다.
3. synthesis/power activity 또는 검증 가능한 외부 자료로 operation energy를 보정한다.
4. Gemmini activation 지원 경로가 존재하는 경우 별도의 microbenchmark로 검증한다.
5. Eyeriss/Gemmini 기존 baseline의 latency residual에 SFU 비용을 임의로 끼워 맞추지 않는다.

완료 조건: 각 SFU profile에 frequency, datatype, operation, source와 tool/version provenance가 기록되고, uncalibrated profile은 absolute fidelity 근거로 사용되지 않는다.

## 수정 대상 파일

| 파일 | 주요 변경 |
|---|---|
| `components/sfu.h`, `components/sfu.cc` | 새 SFU component, timing/traffic/energy event |
| `scheduler/npu.h`, `scheduler/npu.cc` | SFU ownership, config parse, connect/reset, layer activation 전달 |
| `scheduler/scheduler.h`, `scheduler/scheduler.cc` | final output tile identity와 reduction completion event |
| `scheduler/stats.h`, `scheduler/stats.cc` | SFU stats, critical path, network rollup, report |
| `utils/energy_units.h`, `utils/energy_units.cc` | SFU energy component와 cost schema |
| `components/pe.h`, `components/pe.cc` | 미사용 하드코딩 ReLU 제거 또는 deprecated 처리 |
| `configs/accelerators/*` | opt-in SFU profile과 calibration provenance |
| `validation/sfu/*` | formula, event-count, traffic, energy, timeline regression |

현재 build는 `components/*.cc` wildcard를 사용하므로 새 SFU source는 별도의 source list 수정 없이 포함될 수 있다. 실제 구현 시 build 결과로 이를 확인한다.

## 테스트 전략

### Formula와 boundary test

| 항목 | 검증 내용 |
|---|---|
| bypass | Linear 및 0 element에서 cycle/energy/event가 0 |
| lane tail | `N < lanes`, `N = lanes`, `N = lanes + 1`의 chunk/cycle 정확성 |
| operation | ReLU와 Leaky의 latency/energy profile 독립성 |
| parallelism | unit/chip 병렬 시 latency=max, energy=sum |
| invalid config | 0 lane/II, 음수 energy, unknown operation, unsupported placement 거부 |

### Event 정확성 test

- `SFU valid elements = B × K × P × Q`
- output tile identity당 invocation 정확히 한 번
- C/R/S reduction tiling 변화에도 activation count 불변
- accumulator spill/reload 변화에도 activation count 불변
- mapping padding은 valid element에 포함하지 않음
- active chip partition의 element 합이 layer output과 일치

### Traffic·energy test

- fused SFU enable 전후 input/weight/output DRAM traffic 동일
- internal ingress/egress byte가 runtime datatype bitwidth와 일치
- `sfu_op_energy_relu` 2배 perturbation 시 SFU operation energy만 2배
- read/write energy perturbation이 해당 internal event에만 반영
- missing SFU energy key가 unpriced active event로 검출됨
- SFU 추가 후 dynamic/static/network total의 합계 identity 유지

### Timing test

- latency/II 공식이 hand calculation과 일치
- lane 또는 unit 수 증가 시 SFU busy cycle이 예상대로 감소
- queue depth 1은 직렬 contract 재현
- queue depth 증가 시 허용된 범위에서만 overlap 증가
- 느린 SFU가 producer stall과 layer latency를 증가시킴
- 빠른 SFU는 fill/drain 외 cycle이 기존 경로에 대부분 숨겨짐
- SFU latency 변화 후 MAC utilization과 leakage energy window 일치

### Functional test

- SFU pure evaluator와 Nebula ReLU/Leaky/Sigmoid/Tanh 결과 비교
- `FUNCTIONAL` main path에서 activation 중복 적용 없음
- unknown activation이 ReLU로 silent fallback하지 않음
- datatype별 tolerance와 saturation/rounding 정책 확인

### Regression test

- `[sfu]`가 없는 기존 Eyeriss/Gemmini/NVDLA config의 numeric baseline 유지
- SFU가 없는 nonlinear layer는 scope warning 출력
- 기존 timing, traffic, energy, power validation suite 통과
- standalone softmax는 Phase 7 이전까지 명시적으로 timing scope에서 제외

## 호환성 및 migration 정책

- 기존 accelerator config에 `[sfu]`가 없으면 SFU cycle/traffic/energy를 추가하지 않는다.
- 다만 nonlinear activation이 존재하면 report에 activation이 모델링되지 않았음을 표시한다.
- `energy_unit=pJ`와 provenance가 있더라도 활성 SFU event가 unpriced이면 absolute energy/power를 qualified result로 출력하지 않는다.
- 기존 Eyeriss/Gemmini config에 임의의 default SFU cost를 추가해 기존 latency 정합성을 변경하지 않는다.
- calibration이 끝난 SFU-enabled accelerator config는 기존 baseline config와 분리한다.
- `pe_t::activation()` 제거 전까지 호출을 금지하고, 새 SFU event 경로와 동시에 사용되지 않도록 assertion을 둔다.
- Softmax 지원 전에는 단순 per-element activation 비용으로 대체하지 않는다.

## 위험과 의사결정 게이트

1. **최종 output 판별 실패**: reduction completion event가 정확하지 않으면 SFU 구현을 timing에 연결하지 않는다.
2. **functional 중복 적용**: Nebula와 SFU의 output mutation owner가 하나로 고정되지 않으면 functional integration을 진행하지 않는다.
3. **output cast 순서 불명확**: activation과 quantization 순서를 profile에서 확정하기 전에는 low-precision functional fidelity를 주장하지 않는다.
4. **traffic double-counting**: SFU 내부 link와 기존 PE-array↔GLB output link의 ownership이 확정되지 않으면 energy subtotal을 활성화하지 않는다.
5. **낙관적 overlap**: output release time과 queue depth가 없으면 streaming overlap을 추정하지 않고 직렬 모델을 사용한다.
6. **미보정 energy**: 기존 cycle 오차의 residual을 SFU energy/cycle calibration 값으로 사용하지 않는다.
7. **softmax 범위 확대**: element-wise SFU regression이 끝나기 전에 reduction-capable softmax 구현을 시작하지 않는다.

## 권장 구현 순서

```text
activation contract + 기존 baseline
→ SFU component/config parser
→ final_output_tile event와 정확한 element count
→ fused internal traffic + dynamic/static energy
→ 직렬 SFU critical-path latency
→ queue/backpressure 기반 streaming overlap
→ element-wise nonlinear operation 확장
→ standalone softmax multi-pass 모델
→ NVDLA SDP RTL cycle/energy calibration
```

가장 먼저 달성해야 할 milestone은 다음과 같다.

> ReLU/Leaky ReLU가 최종 유효 output element마다 정확히 한 번 실행된 것으로 계수되고, fused DRAM traffic을 변경하지 않으면서 SFU cycle과 energy가 독립적으로 검증되는 상태.

## 최종 완료 조건

- convolution/connected의 모든 fused activation이 명시적인 SFU event로 기록된다.
- activation event는 partial sum이나 mapping repetition에 의해 중복되지 않는다.
- SFU cycle, stall, internal traffic, dynamic/static energy가 layer와 network report에 포함된다.
- SFU를 포함한 critical path, utilization, leakage window가 하나의 timeline contract에서 일치한다.
- fused mode는 activation 때문에 추가 DRAM traffic을 만들지 않는다.
- 활성 SFU event의 비용이 없으면 absolute energy/power qualification이 차단된다.
- `FUNCTIONAL` 결과에서 activation이 중복 적용되지 않는다.
- 기존 SFU-disabled baseline은 numeric regression을 유지한다.
- Softmax는 reduction·exp·normalization·scratchpad traffic을 포함한 별도 multi-pass 모델로 검증된다.
- RTL 또는 외부 근거가 있는 SFU profile만 absolute cycle/energy fidelity의 근거로 사용된다.
