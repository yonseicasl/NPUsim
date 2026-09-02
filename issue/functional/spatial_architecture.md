# Spatial architecture: functional simulation 이슈

## 판정

`spatial_arch_t`의 dense BUS scatter는 동작하지만, PE별 tensor 좌표와 reduction이 보장되지 않고 COO는 데이터를 복사하지 않아도 완료 상태가 된다. 따라서 reduction dimension을 공간에 배치한 convolution 및 sparse COO functional 결과는 신뢰할 수 없다.

- 심각도: **P0 reduction/좌표**, **P1 COO**

## 확인된 문제

| 우선순위 | 문제 | 영향 |
|---|---|---|
| P0 | PE 간 partial-sum reduction이 없음 | 여러 PE의 같은 output 기여가 STORE 순서에 따라 덮어써짐 |
| P0 | input/weight/output offset을 독립 modulo로 선택 | 한 PE의 operand와 output 좌표 결합이 보장되지 않음 |
| P0 | output STORE가 대입 | 같은 output의 partial sum을 보존하지 못함 |
| P1 | `SPARSE_COO` transfer 분기가 비어 있음 | 데이터를 옮기지 않고 data-ready 상태로 전환됨 |
| P1 | output load batch/group index 식 오류 | multi-batch/grouped convolution 초기 partial sum이 틀릴 수 있음 |

## 근거와 수정 방향

Eyeriss mapping은 `PE_Y`에 `R`을 배치할 수 있으므로 C/R/S reduction을 수행해야 한다. 그러나 [`scheduler.cc`](../../scheduler/scheduler.cc#L1611-L1644)의 output store는 누산이 아닌 대입이다.

Scheduler는 세 tensor offset vector를 독립적으로 만든다([`scheduler.cc`](../../scheduler/scheduler.cc#L753-L821)). Spatial transfer도 각 vector 크기로 PE index를 따로 modulo한다([`spatial_arch.cc`](../../components/spatial_arch.cc#L226-L230), [`spatial_arch.cc`](../../components/spatial_arch.cc#L776-L780), [`spatial_arch.cc`](../../components/spatial_arch.cc#L1327-L1332)). PE별 `(B,K,P,Q,C,R,S)` 좌표를 먼저 만들고 세 주소를 함께 유도해야 한다. Weight offset의 filter stride와 group 순회도 함께 수정한다.

`SPARSE_COO` input/weight branch는 빈 블록인데도 이후 상태가 해제된다([`spatial_arch.cc`](../../components/spatial_arch.cc#L333-L335), [`spatial_arch.cc`](../../components/spatial_arch.cc#L881-L883)). 실제 COO value/index를 전달하거나 초기화 단계에서 unsupported로 거부해야 한다.

## 완료 조건

- active PE마다 검증 가능한 공동 tensor coordinate와 세 data offset이 있다.
- C/R/S를 PE에 분산한 작은 convolution이 reference와 일치하며 PE 실행 순서에 독립적이다.
- multi-batch, grouped, depthwise convolution의 output load/store가 reference와 일치한다.
- COO는 실제 compressed functional transfer를 수행하거나 사전 오류를 낸다.
