# Runtime Datatype 지원 진단

> **재판정 (2026-09-02):** 이 문서는 functional/timing 분리([README.md](README.md))
> 이전의 원본이다. **timing/storage/accounting 절반은 해결 완료**
> ([timing/runtime_datatype.md](timing/runtime_datatype.md) — dense descriptor 전 계층
> 연결, sparse는 명시적 거부); **functional 절반(P0: 값이 format을 따르지 않음)은
> 미해결로 유지** ([functional/runtime_datatype.md](functional/runtime_datatype.md)).
> 현황은 분할 문서와 [assessment/ISSUE_REAUDIT_2026-09-02.md](../assessment/ISSUE_REAUDIT_2026-09-02.md)를 보라.


## 판정

현재 codebase는 accelerator configuration의 `[accelerator]` section에서 input, weight, output, accumulator format을 **runtime에 선택**할 수 있다. 지원되는 format 이름은 `fp32`, `fp16`, `bf16`, `int8`, `int4`, `int2`, `uint8`, `mxfp8`, `mxfp4`이며, PE local-buffer capacity/utilization은 선택한 payload bit 수와 MXFP scale metadata bit 수를 반영한다.

그러나 이는 **runtime datatype의 storage/accounting 지원**이다. 실제 FUNCTIONAL tensor representation, quantize/dequantize, rounding/saturation, packed memory, format-aware MAC 및 conversion/scale IP는 아직 구현되지 않았다. 따라서 현재 결과를 low-precision functional simulation 또는 datatype-specific performance/energy model로 해석하면 안 된다.

- 심각도: **P0 / FUNCTIONAL datatype 정확성**, **P1 / traffic·cycle·energy 정확성**
- 영향 범위: FUNCTIONAL build, INT/FP/MX 비교 실험, low-precision accelerator IP 연구

## 현재 구현된 범위

| 항목 | 상태 | 근거 |
|---|---|---|
| Runtime format 선택 | 완료 | `[accelerator]`의 `data_format`, `input_format`, `weight_format`, `output_format`, `accumulator_format`을 읽는다. [`datatype.cc`](../utils/datatype.cc#L67-L82), [`npu.cc`](../scheduler/npu.cc#L59-L63) |
| Scalar format descriptor | 완료 | FP32/FP16/BF16/INT2/4/8/UINT8 parser가 있다. [`datatype.cc`](../utils/datatype.cc#L45-L63) |
| MXFP descriptor/accounting | 부분 완료 | MXFP4/8, 32-element block, 8-bit scale metadata를 storage bit 계산에 포함한다. [`datatype.cc`](../utils/datatype.cc#L54-L60), [`datatype.cc`](../utils/datatype.cc#L104-L119) |
| PE local-buffer capacity/utilization | 부분 완료 | 선택 format의 storage byte 수를 capacity check와 utilization에 사용한다. [`pe.cc`](../components/pe.cc#L499-L508), [`pe.cc`](../components/pe.cc#L623-L629) |
| Config/accounting regression | 완료 | BF16, INT4, MXFP metadata byte accounting을 validation test가 확인한다. [`unittest/validation_test.cc`](../tests/validation_test.cc#L47-L68) |

## 미해결 문제

| 우선순위 | 문제 | 결과 |
|---|---|---|
| P0 | FUNCTIONAL value가 runtime format으로 표현·변환되지 않음 | `fp16`/`bf16`/`int4`/MXFP config를 선택해도 실제 layer value와 MAC 수학은 해당 format이 아님 |
| P0 | `USER_INTEGER`/`USER_FLOAT`와 전역 `data_t`가 남아 있음 | 실제 C++ object type은 여전히 compile time에 정해짐 |
| P0 | MAC IP와 accumulator narrowing이 format-aware하지 않음 | INT overflow/saturation, FP rounding, BF16/FP16 widening, MX scale apply를 모델링하지 못함 |
| P1 | packed storage가 없음 | INT4/INT2 sub-byte packing, FP bit layout, MX payload/scale allocation이 실제 buffer layout과 다름 |
| P1 | memory hierarchy traffic이 전면 전환되지 않음 | 다수의 `sizeof(data_t)` 기반 transfer/DRAM/sparse 계산이 runtime descriptor를 사용하지 않음 |
| P1 | MXFP profile 의미가 고정·검증되지 않음 | exact payload encoding, scale encoding, scale layout, block-axis, special-value 정책이 없음 |
| P1 | conversion/scale IP 비용이 없음 | conversion, packing, decoding, scale load/reuse의 cycle/energy가 누락됨 |
| P2 | datatype별 report와 config compatibility 정책이 부족함 | profile provenance, payload/metadata 분리 통계, legacy config의 의미가 불명확 |

## 상세 근거

### 1. Runtime config는 functional datatype을 바꾸지 않는다

runtime registry는 format descriptor와 storage bit 수만 제공한다. [`datatype.h`](../utils/datatype.h#L21-L50) 반면 [`data.h`](../utils/data.h#L6-L60)는 여전히 compile-time `USER_INTEGER`/`USER_FLOAT` macro와 전역 `data_t`에 의존한다. 기본 build에서는 `uint8_t data_t`를 사용한다. [`user-def.h`](../utils/user-def.h#L4-L14)

따라서 다음 config는 buffer capacity accounting에는 영향을 주지만, FP16 rounding이나 INT4 quantization을 실행하지 않는다.

```ini
[accelerator]
input_format = bf16
weight_format = int4
accumulator_format = fp32
output_format = bf16
```

해결 방향:

1. host FP32 reference tensor와 simulator packed storage를 분리한다.
2. load/store/MAC boundary에서 descriptor 기반 `quantize`, `dequantize`, `convert`를 호출한다.
3. `USER_INTEGER`/`USER_FLOAT`를 deprecated 처리한 뒤 제거한다.

### 2. MAC과 accumulator IP가 없다

현재 MAC은 `data_t`에 대해 직접 곱셈·누산한다. [`data.h`](../utils/data.h#L47-L52), [`pe.cc`](../components/pe.cc#L2306-L2332) operand format, product width, accumulator width, overflow/NaN/Inf, final narrowing을 runtime descriptor로 결정하지 않는다.

필요한 IP/module:

- quantize/dequantize 및 cast IP
- INT signed/unsigned MAC, zero-point 보정, INT32 accumulation, saturation IP
- FP16/BF16 operand decode 및 FP32 accumulation IP
- final output narrowing/rounding IP
- MX payload decode, block-scale load/decode/apply, scale cache/reuse IP

각 IP는 지원 format pair, throughput, pipeline latency, dynamic/static energy를 profile로 가져야 한다.

### 3. physical buffer/traffic은 아직 descriptor와 일관되지 않다

PE local-buffer의 capacity check는 runtime format byte 수를 사용하지만, 실제 allocation과 다른 hierarchy의 transaction 계산에는 여전히 `sizeof(data_t)`가 광범위하게 사용된다. 예를 들어 PE array allocation은 `sizeof(data_t)`로 entry 수를 정한다. [`pe_array.cc`](../components/pe_array.cc#L89-L95) DRAM·global buffer·NoC 및 sparse metadata 수식도 같은 가정을 사용한다.

필요한 변경:

1. logical element count, packed payload bits, metadata bits, aligned allocation bits를 별도 API로 둔다.
2. PE, PE array, global buffer, multi-chip, DRAM의 capacity와 transfer를 `ceil(storage_bits / link_width_bits)` 공통 helper로 전환한다.
3. INT4/INT2 packing tail, line alignment, MX scale stream/interleaving을 명시적으로 모델링한다.
4. sparse metadata와 MX metadata를 별도 stream으로 계산해 double-counting을 방지한다.

### 4. MXFP는 이름과 storage metadata만 있다

현재 MXFP4/8 parser는 32-element block과 8-bit scale metadata라는 accounting profile만 정의한다. [`datatype.cc`](../utils/datatype.cc#L54-L60) 다음 정책은 아직 없다.

- payload encoding (예: exponent/mantissa profile)
- scale encoding과 rounding
- block axis 및 tensor edge padding
- scale metadata layout (separate/interleaved)
- zero/NaN/Inf/subnormal 처리
- scale load/decode/apply cycle 및 energy

MXFP는 scalar alias가 아니라 block format이므로, 위 항목과 reference vector를 profile별로 고정하기 전에는 functional 지원으로 표시하면 안 된다.

## 수정 우선순위

```text
runtime descriptor (완료)
→ FP16/BF16/INT8/INT4 quantize·convert·MAC/accumulator IP
→ packed storage와 memory hierarchy bit accounting 전환
→ MX block-scale encode/decode·traffic·MAC IP
→ calibration, report, CI regression, legacy macro 제거
```

## 완료 조건

- format이 input/weight/output/accumulator별로 runtime config에서 선택되고, unsupported pair는 fail-fast한다.
- FP16/BF16/INT8/INT4의 pack/unpack, rounding, saturation, MAC 결과가 reference vector와 일치한다.
- functional convolution/FC가 selected operand·accumulator·output format을 실제로 사용한다.
- 모든 component의 capacity, line/NoC/DRAM transaction, sparse metadata accounting이 `sizeof(data_t)`가 아닌 descriptor bit 수를 사용한다.
- MXFP profile별 payload/scale encoding, block boundary, padding, scale traffic, decode/MAC 비용이 테스트와 report에 반영된다.
- conversion, scale, MAC, accumulator spill event의 cycle/energy가 별도 보고된다.
- `USER_INTEGER`/`USER_FLOAT` compile-time selection이 제거되고, legacy config의 기본 의미가 문서화된다.

## 최종 결론

현재는 **runtime-configurable datatype accounting의 기반**만 구현됐다. 다양한 datatype의 정확한 functional 동작과 datatype-specific hardware IP를 지원했다고 볼 수 없다. 다음 구현은 scalar FP16/BF16/INT8/INT4의 conversion·MAC·packed storage를 먼저 완성한 뒤 MXFP block-scale IP로 확장해야 한다.
