# Data precision timing implementation

## 상태: hierarchy-wide dense datatype 완료

Runtime format descriptor를 PE/MAC부터 local buffer, PE array, global buffer, multi-chip/NoP, DRAM까지 같은 bit transaction 계약으로 연결했다. 지원 dense format의 payload, MXFP scale metadata, line access, link serialization, capacity 및 통계를 모든 modeled component boundary에서 일관되게 계산한다.

## 구현된 변경

### Runtime format 및 MXFP layout

- input, weight, output, accumulator별 runtime format을 지원한다: BF16, FP16/FP32, INT8/INT4/INT2, UINT8, MXFP8, MXFP4.
- `mxfp_metadata_layout = separate|interleaved`를 `[accelerator]` profile에서 설정한다. 기본값은 `separate`이다.
- separate layout은 payload와 scale metadata의 byte/link tail을 각각 정렬한다.
- interleaved layout은 payload와 scale을 합친 bit stream을 정렬한다.
- `payload_bits()`, `metadata_bits()`, `storage_bits()`, `payload_bytes()`, `metadata_bytes()`, `storage_bytes()`를 제공한다.
- `payload_transactions()`, `metadata_transactions()`, `storage_transactions()`를 제공한다.

예를 들어 MXFP4 block-32의 33 elements는 32-bit link에서 separate layout이면 payload 5 + scale 1 = 6 transaction, interleaved layout이면 5 transaction이다.

### 적용 경로

- PE MAC/register/local-buffer capacity, utilization, access 및 MAC↔local-buffer link
- Spatial, systolic-array, adder-tree PE-array의 dense INPUT/WEIGHT/OUTPUT distribution과 write-back
- GLB↔PE-array의 source/destination line access 및 link serialization
- global-buffer separate/shared capacity와 utilization
- multi-chip temporal buffer capacity와 GLB↔multi-chip NoP distribution
- DRAM↔multi-chip INPUT/WEIGHT load 및 OUTPUT write-back
- MAC↔local-buffer, PE-array NoC, GLB, NoP, DRAM 경계별 payload/metadata/serialized transaction report

각 component의 `line_size`와 interconnect `bitwidth`는 bits로 해석한다. memory line은 8-bit 이상 power-of-two인지 검사하며, source access, destination access, link transaction 중 가장 긴 stage count로 안전한 pipeline cycle을 계산한다.

### PE format-IP accounting

PE profile은 다음 optional per-transaction 비용을 지원한다. 미설정 시 모두 0이다.

- `format_payload_cycle`, `format_payload_energy`: payload conversion/pack/unpack/decode 경로
- `format_metadata_cycle`, `format_metadata_energy`: MX scale metadata load/reuse/apply 경로
- `accumulator_spill_cycle`, `accumulator_spill_energy`: output accumulator spill 경로

각 PE request에서 payload/metadata/storage transaction 수로 cycle과 energy를 누적하고, layer 결과 파일에 input/weight/output별 `Format-IP cycle`과 `Format-IP energy`를 기록한다.

## 검증

`unittest/validation_test.cc`는 FP32, FP16, BF16, INT8, INT4, INT2, UINT8 전 형식과 MXFP8/MXFP4의 payload/metadata bit 수, packing tail, separate/interleaved transaction boundary를 확인한다. INT2/MXFP8/FP32가 서로 다른 source line, destination line, link width를 통과하는 hierarchy micro-test도 포함한다.

```bash
./unittest/run_validation.sh
```

실제 Eyeriss AlexNet timing 실행에서도 INT4와 MXFP4를 각각 검증했으며, MXFP4 결과에서 PE-array NoC, NoP, DRAM 모두 scale metadata transaction이 별도로 보고된다.

## 이 작업과 분리된 남은 범위

- CSC/COO 같은 sparse encoding의 value/index/pointer stream과 decoder 비용. timing entry는 현재 sparse를 명시적으로 거부하므로 dense 비용으로 잘못 계산하지 않는다.
- format-IP의 상세 throughput 및 다른 pipeline과의 overlap. 현재 configurable per-transaction cycle/energy 비용은 제공한다.
- bank contention, routed NoC/NoP, multicast/back-pressure처럼 datatype과 독립적인 architecture timing 정밀도
