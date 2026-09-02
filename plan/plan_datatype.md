# 데이터 정밀도·데이터타입 확장 계획

## 목적

NPUSim이 단일 compile-time `data_t`와 고정 `sizeof(data_t)` 가정에서 벗어나, layer·tensor·hardware component별로 서로 다른 데이터 정밀도를 모델링하도록 확장한다. 지원 대상은 다음 순서로 한다.

1. 기준 형식: FP32 (reference / accumulation)
2. 표준 scalar 형식: FP16, BF16, INT8, INT4
3. 확장 scalar 형식: INT2, UINT8 및 필요 시 FP8 profile
4. block-scaled 형식: MXFP 계열(MXFP8, MXFP4 등)

목표는 단지 functional value를 변환하는 것이 아니다. 저장 용량, packing, scale metadata, conversion/rounding/saturation IP, MAC datapath, accumulator, memory/NoC traffic, latency, energy까지 같은 datatype contract에서 계산하는 것이다.

## 현재 상태와 문제

현재 [`utils/data.h`](../utils/data.h)는 `USER_INTEGER` 또는 `USER_FLOAT` compile-time macro로 하나의 전역 `data_t`만 선택한다. 기본 타입은 `uint8_t`이고, `DATA_BIT=3` 또는 `16`도 실제 C++ object size와 storage traffic을 일치시키지 못한다. [`utils/user-def.h`](../utils/user-def.h)

또한 component/scheduler/DRAM 전반에서 `sizeof(data_t)`가 buffer entry 수, line transaction, sparse metadata cost 계산에 직접 사용된다. 이 의존이 약 552개이므로, `data_t`만 교체하면 다음 문제가 생긴다.

- tensor의 logical element 수와 physical storage byte 수가 분리되지 않는다.
- INT4/INT2 packing, FP8/BF16 bit layout, MXFP block scale metadata를 표현할 수 없다.
- input/weight/output/accumulator가 서로 다른 format일 수 없다.
- Nebula layer의 float pointer를 `data_t*`로 reinterpret cast하는 현재 DRAM 경로는 type-safe하지 않다.
- conversion, scale, rounding, saturation, widening MAC의 cycle/energy가 없다.

## 설계 원칙과 범위

### 핵심 원칙

1. **runtime descriptor 우선**: format 선택을 build macro가 아니라 accelerator/network config의 descriptor로 표현한다. FP32 reference mode만 예외적으로 고정한다.
2. **logical value와 physical storage 분리**: functional reference value, packed payload, transfer bit 수를 하나의 C++ object size로 표현하지 않는다.
3. **operand·accumulator 분리**: input/weight format, partial-sum accumulator format, output format 및 conversion 시점을 독립적으로 지정한다.
4. **명시적 rounding 정책**: round-to-nearest-even, stochastic rounding, truncate 및 saturation/overflow/NaN/Inf 처리를 datatype profile의 일부로 만든다.
5. **MX는 block format으로 취급**: MXFP는 scalar format의 alias가 아니다. payload element format, block shape/axis, shared scale format, scale layout, scale traffic을 함께 가져야 한다.
6. **functional/cost 동치**: 한 transfer/MAC/conversion event가 functional path와 traffic/cycle/energy path에 동시에 기록돼야 한다.
7. **미지원 조합은 fail-fast**: 모델링되지 않은 format, kernel, scale layout, hardware IP 조합은 silent fallback 없이 config validation에서 실패시킨다.

### 초기 지원 계약

| 역할 | 1차 지원 | 비고 |
|---|---|---|
| input / weight storage | FP32, FP16, BF16, INT8, INT4 | per-tensor 또는 per-channel scale 지원 |
| MAC operand | FP16, BF16, INT8, INT4 | 실제 PE IP profile이 지원하는 조합만 허용 |
| accumulator | FP32, INT32 | 저정밀 operand의 widened accumulation 기본값 |
| output storage | FP32, FP16, BF16, INT8 | quantize 위치를 output commit으로 고정 |
| MX | MXFP8/MXFP4 profile | block scale metadata와 payload/scale traffic을 함께 모델링 |

FP8/MXFP의 exact encoding, scale type, block size는 vendor 이름을 코드에 박아 넣지 않는다. named profile과 명시적 필드로 정의하고, 각 profile의 사양과 reference vector를 테스트에 고정한다.

## 목표 아키텍처

```text
Network host tensor (FP32 reference)
        │ import / quantize
        ▼
Tensor storage descriptor ── packed payload + optional scale/zero-point metadata
        │ transfer (payload bits + metadata bits)
        ▼
Component buffer / register file
        │ decode / dequantize / cast IP
        ▼
MAC operand format ── multiply / widen / accumulate ── accumulator format
        │ final conversion / activation
        ▼
Output storage descriptor
```

### 새 핵심 타입

`utils/datatype.{h,cc}` 또는 동등한 독립 module을 만든다.

- `data_format_kind_t`: FP32, FP16, BF16, INT, UINT, FP8 profile, MX
- `scalar_format_t`: sign, exponent, mantissa 또는 integer bitwidth/signedness
- `quantization_spec_t`: scale/zero-point granularity, rounding, saturation, calibration source
- `block_format_spec_t`: payload format, block shape/axis, scale format, scale layout, scale count 규칙
- `tensor_format_t`: payload format + quantization/block metadata + alignment/packing policy
- `mac_format_t`: input format, weight format, product format, accumulator format, output conversion format
- `format_profile_t`: config에서 참조 가능한 named profile

storage API는 다음 역할을 제공한다.

- `logical_element_count()`
- `payload_bits(element_count)` / `metadata_bits(element_count, shape)` / `total_storage_bits()`
- `pack_from_reference()` / `unpack_to_reference()`
- `quantize()` / `dequantize()` / `convert()`
- `multiply_accumulate()` 및 overflow/NaN/Inf 상태 기록

C++ `data_t`는 최종적으로 제거하거나 reference-simulation 전용 `float` alias로 축소한다. 모든 component buffer는 raw packed byte storage와 `tensor_view_t`/descriptor를 통해 접근하도록 한다.

## 단계별 실행 계획

### Phase 0 — 요구사항 고정과 baseline 확보

1. datatype 지원 matrix를 작성한다: storage, register, MAC operand, accumulator, output 각각에 허용되는 조합을 명시한다.
2. MX profile의 block size, scale type, scale granularity, rounding 및 special-value 정책을 결정한다.
3. 현재 FP32/기본 integer 결과, tile size, traffic, cycle, energy를 작은 convolution/FC workload로 golden baseline화한다.
4. `sizeof(data_t)`와 `USER_INTEGER`/`USER_FLOAT` 의존 목록을 자동 검사로 고정한다.
5. 기존 PE issue의 activation/reduction 정확성 문제가 해결되기 전에는 datatype 정확도 baseline을 확정하지 않는다.

산출물: `docs/datatype-contract.md`, support matrix, reference test vector, migration inventory.

### Phase 1 — descriptor와 reference conversion library

1. `datatype` module에 scalar/quantized/block format descriptor와 strict parser를 구현한다.
2. FP32↔FP16, FP32↔BF16, FP32↔INT8/INT4 변환을 구현한다.
3. rounding, saturation, signedness, zero-point, per-tensor/per-channel scale을 테스트 가능하게 분리한다.
4. host tensor는 FP32 reference representation을 유지하고, simulator boundary에서 명시적으로 import/export한다. reinterpret cast를 제거한다.
5. config에 named profile 및 per-data-type format을 추가한다.

예시:

```ini
[datatype_profile.int8_sym]
kind = int
bits = 8
signed = true
scale_granularity = per_tensor
rounding = nearest_even
saturation = clamp

[spatial_arch]
input_format = bf16
weight_format = mxfp4_e2m1_b32
accumulator_format = fp32
output_format = bf16
mac_profile = bf16_mxfp4_to_fp32
```

완료 조건: conversion unit test가 bit-exact reference vector와 일치하고, invalid profile/조합이 config parse에서 실패한다.

### Phase 2 — storage와 transfer bit accounting 분리

1. buffer capacity를 element count/byte allocation이 아닌 `storage_bits`와 alignment로 관리한다.
2. `sizeof(data_t)` 기반 local buffer, PE array, global buffer, multi-chip, DRAM의 capacity 및 transaction 계산을 descriptor API로 치환한다.
3. line size와 bus bitwidth를 bit 단위로 통일하고 `ceil(total_bits / link_width_bits)`를 공통 helper로 만든다.
4. INT4/INT2의 sub-byte packing, alignment, line tail padding을 구현한다.
5. MX payload와 scale metadata를 별도 stream 또는 interleaved layout으로 모델링하고, 두 방식 중 지원 layout을 config에 명시한다.
6. sparse metadata와 datatype metadata를 별도 accounting해 double-counting을 방지한다.

완료 조건: 동일 logical tile에 대해 FP32/FP16/BF16/INT8/INT4의 payload bit 수와 transaction 수가 hand calculation과 일치한다. MX tile은 payload와 scale traffic이 각각 보고된다.

### Phase 3 — functional PE datapath와 IP module

1. PE에 format-aware load/decode, operand conversion, MAC, accumulator narrowing, output encode 단계를 추가한다.
2. `mac_ip_t`를 도입해 지원 operand pair, accumulator type, lane count, pipeline latency, throughput, dynamic energy를 profile로 정의한다.
3. INT MAC은 signedness, zero-point 보정, INT32 accumulation 및 saturation을 명시한다.
4. FP16/BF16 MAC은 operand rounding과 FP32 accumulation을 우선 지원한다.
5. MX MAC은 block scale load/decode, payload decode, scale apply 위치(pre-multiply 또는 equivalent), scale reuse/cache를 모델링한다.
6. conversion/scale decode/pack IP도 MAC과 별도의 cycle/energy event로 기록한다.
7. PE stationary mode가 format metadata와 scale block 경계를 보존하도록 scheduler tile validation을 추가한다.

완료 조건: small GEMM/convolution의 functional result가 선택한 reference backend의 허용 오차 또는 bit-exact integer 기준을 만족한다. MAC/conversion/scale event 수가 trace로 검증된다.

### Phase 4 — scheduler, mapping, memory hierarchy 통합

1. scheduler가 tensor shape뿐 아니라 format descriptor와 MX block partition을 전달하도록 확장한다.
2. tile을 MX block 경계에 맞추지 못하는 mapping은 padding/fragmentation 비용을 계산하거나 fail-fast한다.
3. output partial sum은 storage format과 분리된 accumulator format으로 유지하고 최종 commit에서만 narrow/quantize한다.
4. PE array, global buffer, multi-chip, DRAM의 transfer API를 typed tensor view로 전환한다.
5. NoC/DRAM bandwidth, local-buffer port, scale metadata port contention을 component profile로 모델링한다.
6. report에 payload bytes, metadata bytes, conversion cycles/energy, MAC cycles/energy, accumulator spills를 분리해 출력한다.

완료 조건: dataflow/stationary 변경에도 동일 logical format과 MX scale semantics가 유지되고, report 합계가 component event 합계와 일치한다.

### Phase 5 — performance/energy calibration과 regression

1. datatype/IP profile별 unit cycle·energy·area parameter schema를 정의한다.
2. published RTL, synthesis, vendor documentation 또는 사용자가 제공한 characterization 중 출처를 profile metadata에 기록한다.
3. FP32 baseline 대비 FP16/BF16/INT8/INT4/MX의 traffic, throughput, conversion overhead를 비교하는 microbenchmark를 추가한다.
4. functional accuracy, bit accounting, latency/energy accounting, config rejection을 CI regression에 포함한다.
5. optional stochastic rounding은 deterministic seed와 statistical test를 사용한다.

완료 조건: profile별 provenance가 없는 성능/에너지 수치는 report에서 "uncalibrated"로 표시되고, CI가 datatype regression을 실행한다.

## 테스트 전략

| 계층 | 필수 테스트 |
|---|---|
| encoding | FP16/BF16/INT4/INT8 pack-unpack, rounding tie, saturation, NaN/Inf/subnormal 정책 |
| quantization | per-tensor/per-channel scale, zero-point, signed/unsigned, calibration failure |
| MX | block 경계, scale count/layout, zero block, mixed scale, padding tile |
| MAC | INT4×INT4→INT32, BF16/FP16→FP32, overflow, final narrowing |
| transfer | payload/metadata bits, line alignment, link-width ceiling, sparse+MX accounting |
| integration | convolution/FC across three stationary modes, heterogeneous input/weight/output/accumulator formats |
| negative | unsupported format pair, invalid exponent/mantissa, block misalignment, unavailable IP profile |

reference 비교는 integer는 bit-exact, floating/MX는 profile별 ULP 또는 absolute/relative tolerance를 사용한다. PyTorch 등 외부 reference backend는 optional oracle로 두되, core conversion vector는 repository 내부에 고정한다.

## 호환성 및 migration 정책

- 기존 config에 datatype field가 없으면 legacy profile을 명시적으로 선택한다. 기본값은 문서화된 FP32 reference 또는 기존 byte behavior 중 하나로 결정하되, 두 의미를 혼용하지 않는다.
- `USER_INTEGER`/`USER_FLOAT`는 Phase 1 이후 deprecated warning, Phase 4 이후 제거한다.
- 새 format 도입 중에는 old/new traffic report를 동시에 출력하는 compatibility mode를 제공해 결과 차이를 검토한다.
- packed storage migration은 component 단위로 feature flag를 두지 않는다. 혼합 상태에서 capacity는 old path, traffic은 new path로 계산되는 것을 금지한다.

## 위험과 의사결정 게이트

1. **MX 사양 선택**: block size와 scale encoding이 정해지지 않으면 MX 구현을 시작하지 않는다.
2. **functional 정확도와 hardware cost 분리 실패**: reference value만 FP32로 유지하고 low-precision event를 생략하는 shortcut을 금지한다.
3. **metadata 누락**: MX scale 및 sparse metadata의 storage/traffic/energy ownership을 명확히 지정한다.
4. **accumulator overflow**: operand precision보다 accumulator 정책을 먼저 확정한다.
5. **과도한 범위 확대**: FP8/MX는 Phase 3의 scalar FP16/BF16/INT 검증이 끝난 뒤 추가한다.

## 권장 구현 순서

```text
format contract + reference vectors
→ descriptor/conversion library
→ bit-based storage/transfer accounting
→ INT8/INT4 + FP16/BF16 PE MAC IP
→ scheduler/memory hierarchy typed views
→ MX block-scale storage + decode/MAC IP
→ calibration, CI regression, legacy macro 제거
```

## 최종 완료 조건

- format이 tensor/component/MAC accumulator 단위로 명시되고 unsupported 조합은 fail-fast한다.
- `sizeof(data_t)`와 compile-time global datatype 선택이 capacity/traffic/MAC 의미를 결정하지 않는다.
- FP16, BF16, INT8, INT4 및 선택한 MX profile이 functional, storage, traffic, cycle, energy report에 일관되게 반영된다.
- MX의 payload와 scale metadata가 allocation, transfer, decode, energy report에 모두 포함된다.
- heterogeneous precision layer의 output이 reference 정책과 일치하고, 각 datatype IP event가 traceable하다.
