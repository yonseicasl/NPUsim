# Runtime datatype: functional simulation 이슈

## 판정

현재 runtime format 설정은 storage capacity/accounting만 바꾼다. `fp16`, `bf16`, `int4`, MXFP를 선택해도 functional tensor의 표현·변환·MAC 수학은 선택한 format을 따르지 않는다. 따라서 low-precision functional 결과는 신뢰할 수 없다.

- 심각도: **P0**
- 영향: INT/FP/MX 결과 비교, low-precision convolution/FC functional build

## 확인된 문제

| 우선순위 | 문제 | 영향 |
|---|---|---|
| P0 | tensor value가 runtime format으로 quantize/dequantize되지 않음 | config가 바뀌어도 값과 rounding이 바뀌지 않음 |
| P0 | `USER_INTEGER`/`USER_FLOAT` 및 전역 `data_t`가 compile-time type으로 남아 있음 | 실제 object representation이 runtime descriptor와 다름 |
| P0 | MAC/accumulator가 format-aware하지 않음 | INT overflow/saturation, FP rounding, narrowing을 검증할 수 없음 |
| P0 | MXFP payload·scale profile이 정의되지 않음 | block format의 functional semantics가 없음 |

## 근거와 필요한 변경

- Format parser는 `fp32`, `fp16`, `bf16`, `int8`, `int4`, `int2`, `uint8`, `mxfp8`, `mxfp4` descriptor를 만들지만, [`data.h`](../../utils/data.h#L6-L60)의 functional value는 여전히 compile-time `data_t`에 의존한다. 기본 build는 [`user-def.h`](../../utils/user-def.h#L4-L14)의 `uint8_t`를 사용한다.
- [`pe.cc`](../../components/pe.cc#L2306-L2332)의 MAC은 operand/product/accumulator descriptor를 사용하지 않는다. 따라서 signedness, zero-point, accumulator widening, overflow, NaN/Inf와 output narrowing 정책이 없다.
- MXFP parser는 block 32개와 8-bit scale metadata라는 accounting profile만 제공한다. payload encoding, scale encoding·rounding, block axis/edge padding, special-value 정책과 scale 적용 동작이 정의돼 있지 않다.

다음 경계에서 descriptor 기반 동작을 도입해야 한다.

1. Host FP32 reference tensor와 simulator의 packed/encoded storage를 분리한다.
2. Load/store와 operand/accumulator 경계에서 `quantize`, `dequantize`, `convert`를 수행한다.
3. INT MAC(부호·zero-point·saturation), FP16/BF16 decode와 accumulation, final narrowing을 format pair별로 정의한다.
4. MXFP는 payload/scale encoding과 block layout을 profile로 고정한 후 decode·scale apply를 구현한다.

## 완료 조건

- input, weight, accumulator, output format별 지원 조합과 미지원 조합의 fail-fast 정책이 있다.
- FP16/BF16/INT8/INT4 pack/unpack, rounding, saturation, MAC 결과가 reference vector와 일치한다.
- functional convolution/FC가 선택한 operand·accumulator·output format을 실제로 사용한다.
- MXFP profile별 payload/scale encoding, block boundary, padding, special value가 regression으로 고정된다.
- `USER_INTEGER`/`USER_FLOAT` compile-time selection을 제거하거나 runtime functional 경로에서 사용하지 않는다.
