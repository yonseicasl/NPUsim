# PE timing implementation

## 상태: 부분 완료

PE lane 계약, overlap cycle 경계, static-energy 시간축, runtime datatype transaction 및 선택적 format-IP event accounting을 구현했다. sparse execution 자체는 simulator 초기화 단계에서 아직 허용하지 않으므로 PE timing issue의 전체 완료 조건은 충족하지 않았다.

## 구현된 변경

### MAC lane contract

- `number_of_macs`: 독립 accumulator/MAC unit 수
- `mac_width`: accumulator당 scalar FMA lane 수
- mapping의 active MAC factor: 활성 scalar FMA lane 수
- scalar register/lane capacity: `number_of_macs × mac_width`

`update_tile_size()`는 active lane capacity를 검사하고 active accumulator unit, 마지막 active unit의 lane 수, `active scalar lanes / total scalar lanes` utilization을 계산한다. Computation energy는 활성 scalar FMA 수에 비례하고 computation cycle은 병렬 lane issue step마다 누적한다.

### Cycle 및 static energy

- local-buffer↔MAC overlap은 0, 1, 2 이상 transaction을 분기해 unsigned underflow 없이 `double` cycle로 계산한다.
- PE elapsed cycle은 compute, access, link transfer, overlap/write-back, format-IP cycle 중 최댓값이다.
- `[pe] static_energy`는 pJ/cycle이며 이 **PE-array elapsed cycle**에 한 번만 곱한다. inactive PE leakage도 포함한다.

> ⚠️ **정정 (2026-08-15):** 위 static energy 시간창은 "layer elapsed cycle"이 아니라 **PE-array 내부 elapsed cycle**이다(`stats.cc:400-421`). DRAM/NoP/GLB stall 시간이 빠지므로 memory-bound 레이어에서 leakage가 크게 과소 집계된다. 또한 아래 MAC lane grouping(`active_mac_width`/`utilization_mac`)은 계산되지만 실제 cycle/utilization 경로에서 읽히지 않는다(utilization은 `stats.cc:559-560`에서 별도 재계산). 이슈: [issue/timing/static_energy.md](../../issue/timing/static_energy.md), [issue/timing/pe.md](../../issue/timing/pe.md).

### Datatype 및 format-IP event

PE input/weight/output request는 runtime descriptor의 payload/metadata/storage transaction 수로 optional format-IP cost를 누적한다. Dense timing의 local-buffer↔MAC load, MAC→local-buffer accumulator write-back, local-buffer→PE-array flush도 `datatype_transfer_timing()`으로 source/destination line access와 serialized link transaction을 계산한다.

- payload event: conversion, pack/unpack, decode를 profile unit cost로 모델링
- metadata event: MX scale load/reuse/apply를 profile unit cost로 모델링
- output event: accumulator spill을 profile unit cost로 모델링

Layer 결과는 input/weight/output별 `Format-IP cycle`과 `Format-IP energy`, `MAC <-> local-buffer transactions (payload/metadata/serialized)`를 출력한다.

## 검증

```bash
./unittest/run_validation.sh
git diff --check
```

또한 수정한 `datatype`, `pe`, `global_buffer`, `multi_chip`, `stats` translation unit을 `-Wall -Wextra -Werror`로 컴파일한다.

INT4/MXFP4 Eyeriss AlexNet 대조 실험은 MAC↔local-buffer부터 DRAM까지 INT4 metadata가 0이고 MXFP scale metadata가 별도 집계됨을 확인한다.

`unittest/validation_test.cc`의 PE lane micro-test는 8×8 full lane, 20/64 partial lane(3 accumulator, 마지막 4 lane, 31.25% utilization), 0/65 lane capacity boundary, 0.5-cycle issue의 fractional 누적을 검증한다.

## 남은 완료 조건

- dense와 sparse가 동일한 functional/cost event graph를 공유하도록 정리
