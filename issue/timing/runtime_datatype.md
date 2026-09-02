# Runtime datatype: timing simulation 이슈

## 판정

지원 dense runtime datatype의 hierarchy-wide timing 연결은 완료됐다. PE/MAC, local buffer, PE array, global buffer, multi-chip/NoP, DRAM은 동일한 descriptor bit transaction을 사용하며 payload, MXFP metadata, serialized traffic을 경계별로 보고한다.

- 완료: **dense datatype storage/capacity/transfer/reporting**
- 별도 미지원: **sparse encoding timing, 상세 format-IP pipeline overlap**

## 확인된 문제

| 우선순위 | 문제 | 영향 |
|---|---|---|
| 완료 | packed storage와 hierarchy-wide transfer를 공통 descriptor로 연결 | INT4/INT2/MXFP line·NoC·NoP·DRAM traffic을 bit transaction으로 계산 |
| 완료 | dense timing 경로의 host `sizeof(data_t)` 의존 제거 | runtime format과 endpoint access 수가 일치 |
| 완료 | MXFP payload와 scale metadata stream 분리 | 경계별 payload/metadata/serialized 통계 제공 |
| 부분 완료 | configurable per-transaction format-IP 비용 제공 | 상세 throughput/overlap은 남음 |
| 완료 | format provenance 및 payload/metadata report 추가 | 계층별 결과 해석 가능 |

## 구현 결과

`datatype_transfer_timing()`은 logical element 수와 runtime format으로 payload/metadata bit 수를 만든 뒤 source line, destination line, link width에 각각 정렬한다. PE-array, GLB, multi-chip, DRAM은 이 결과를 access cycle·energy와 transfer cycle·energy에 공통 사용한다. MXFP `separate` layout은 payload와 scale stream tail을 별도로 정렬하고, `interleaved`는 결합된 physical stream을 직렬화한다.

Sparse compression은 datatype과 다른 encoding 계층이다. 아직 index/pointer traffic이 구현되지 않았으므로 timing 진입점에서 명시적으로 거부하여 dense 비용으로 묵시 계산하는 오류를 막는다.

> **A7 sparse timing 종결 판단 (2026-08-15):** sparse timing(CSC/COO value/index/pointer stream + decoder 비용)은 **functional sparse 실행이 선행 조건**이다 — 현재 `npu.cc`가 `compression_type != dense`를 최상위 init에서 거부하므로("sparse PE execution is not implemented") sparse 데이터 흐름 자체가 시뮬레이터에 존재하지 않는다. functional sparse 이전에 timing만 만들면 검증 불가능한 코드가 된다. 전 계층(DRAM OUTPUT 포함, DR1)의 명시적 거부 guard가 완비되어 silent-dense 오류는 차단된 상태이며, functional sparse 구현 시 timing 모델을 후속한다.

## 완료 조건

- [x] PE, PE array, global buffer, multi-chip, DRAM의 dense capacity와 transfer가 descriptor bit 수만 사용한다.
- [x] INT4/INT2 packing tail와 alignment를 포함한 transaction 수가 micro-case와 일치한다.
- [x] MXFP payload와 scale traffic이 분리되어 보고된다. Sparse는 구현 전까지 명시적으로 거부한다.
- [~] format conversion·scale·spill의 configurable cycle/energy가 보고된다. 상세 throughput과 pipeline overlap은 남는다.
- [x] report에 input/weight/accumulator/output format과 hierarchy별 payload/metadata breakdown이 포함된다.
