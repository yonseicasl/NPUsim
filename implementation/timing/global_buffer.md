# Global buffer timing implementation

## 상태: dense **load** 경계 완료, **write-back 미완료**

GLB↔PE-array와 GLB↔multi-chip/DRAM의 dense INPUT/WEIGHT/OUTPUT **load** 경계는 runtime datatype descriptor의 bit transaction으로 계산한다. Sparse compression과 bank/contention 모델은 별도 범위다.

> ⚠️ **정정 (2026-08-15):** GLB→multi-chip **OUTPUT write-back(`flush_data`)** 의 link transfer 비용은 실제로는 주석 처리되어 있어(`global_buffer.cc:1497-1500` 등) `write_back_cycle`이 0으로 남고, write-back access는 여전히 host-pointer `mask_bits` 주소 walk를 쓴다. 따라서 "OUTPUT 경계까지 descriptor로 완료"는 **load 방향에만 해당**한다. 이슈: [issue/timing/global_buffer.md](../../issue/timing/global_buffer.md) GB1/GB2.

## 구현 계약

- `line_size`와 link `bitwidth`의 단위는 bits다. GLB line size는 8 이상 power-of-two bits, link width는 0이 아닌 값으로 초기화에서 검증한다.
- `datatype_transfer_timing()`은 source line access, destination line access, serialized link transaction, pipeline transaction 수를 한 번에 만든다.
- GLB read/write access cycle·energy는 endpoint transaction 수를, transfer cycle·energy는 serialized link transaction 수를 사용한다.
- pipeline latency는 공통 `pipelined_transfer_cycles()`로 0/1/다수 transaction을 안전하게 처리한다.
- MXFP separate layout은 payload와 scale metadata stream의 tail을 독립 정렬한다. interleaved layout은 합친 physical stream을 직렬화한다.
- 결과 파일은 `Link transactions (payload/metadata/serialized)`를 GLB 섹션에 기록한다.
- timing build에서 sparse CSC/COO는 dense 비용으로 묵시적으로 계산하지 않고 오류로 중단한다.

## 검증

`unittest/validation_test.cc`의 micro-test는 MXFP4 33 elements에 대해 source 16-bit, destination 32-bit, link 64-bit일 때 `10/6/4` transaction과 `payload=3`, `metadata=1`을 검증한다. INT4 tail case도 별도로 검증한다.

```bash
./unittest/run_validation.sh
./npusim.sh build lib
./npusim.sh run eyeriss alexnet energy
```

## 남은 범위

- sparse payload/index/pointer의 hierarchy-wide traffic model
- GLB bank/arbitration/read-write conflict/multicast/back-pressure
