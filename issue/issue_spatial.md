# Spatial architecture 논리 구현 진단

> **재판정 (2026-09-02):** 이 문서는 functional/timing 분리 이전의 원본이다.
> timing 절반은 해결 완료([timing/spatial_architecture.md](timing/spatial_architecture.md)
> — SP1–3 + NVDLA nv_small holdout 0.55% 외부 검증); functional 절반(reduction/좌표
> 결합 P0)은 [functional/spatial_architecture.md](functional/spatial_architecture.md)에
> 잔존. 현황: [assessment/ISSUE_REAUDIT_2026-09-02.md](../assessment/ISSUE_REAUDIT_2026-09-02.md).


## 판정

현재 `spatial_arch_t`는 세 배열 구현 중 초기화와 BUS 기반 데이터 scatter가 가장 완성되어 있다. 그러나 **지원 범위는 사실상 BUS 기반 generic PE array에 한정**된다. Mesh는 빈 분기이고 crossbar/store-and-forward 실행 경로는 없으며, PE 간 reduction과 topology-aware 좌표·hop 모델도 없다.

따라서 다음과 같이 구분해야 한다.

- BUS 기반 PE-array memory/access 통계: **부분 구현**
- 일반적인 spatial dataflow 표현: **불완전**
- Mesh/crossbar/store-and-forward NoC: **미구현**
- PE 간 partial-sum reduction: **미구현**

- 심각도: **P0 / 제공 설정 실행 불가 및 functional 결과 신뢰성**
- 분석 방법: 실행 결과를 사용하지 않은 정적 코드·상태 전이·수식 분석
- 영향 범위: Eyeriss/EyerissV2/Simba/FSD 설정과 mapping

## 긍정적으로 구현된 부분

Spatial 구현에는 다음 기본 요소가 존재한다.

- height×width에 따른 PE 생성([`spatial_arch.cc`](../components/spatial_arch.cc#L26-L29))
- `pe_stationary` 및 `mac_stationary` 설정 처리
- 공통 temporal buffer 초기화
- `is_waiting_data`, transfer/access 통계 vector 초기화([`spatial_arch.cc`](../components/spatial_arch.cc#L122-L178))
- BUS 경로의 dense input/weight/output scatter
- 일부 sparse format의 metadata traffic 추정
- active PE 기반 utilization과 per-PE 실행

따라서 완전히 비어 있는 클래스는 아니다. 다만 위 구현은 spatial topology 자체보다 PE-array와 bus memory traffic을 모델링한다.

## 확인된 문제

| 우선순위 | 문제 | 결과 |
|---|---|---|
| P0 | Mesh 분기가 비어 있음 | mesh 설정 시 요청과 데이터 상태가 진행되지 않음 |
| P0 | Crossbar/store-and-forward 실행 경로 없음 | 제공 FSD 설정이 첫 전송에서 실패 |
| P0 | PE 간 partial-sum reduction 부재 | reduction 차원을 spatial mapping에 배치하면 결과 덮어쓰기 |
| P0 | output STORE가 대입 연산 | 동일 output에 대한 여러 PE 기여를 보존하지 못함 |
| P0 | independent modulo offset | PE의 input/weight/output 공동 좌표가 보장되지 않음 |
| P1 | physical `(x,y)` 위치가 데이터 경로에 사용되지 않음 | array height/width가 hop·routing·latency에 영향 없음 |
| P1 | BUS contention/broadcast 모델이 단순 access count | arbitration, fanout, multicast group을 표현하지 못함 |
| P1 | sparse COO가 빈 분기 | FUNCTIONAL에서 데이터를 옮기지 않고 존재한다고 표시 |
| P1 | tile-read flag가 layer 수명 동안 reset되지 않음 | 반복 전송 access/cycle 통계 과소 집계 가능 |
| P1 | active shape 검증 부재 | active_x/active_y가 물리 width/height와 각각 불일치 가능 |

## 상세 근거

### 1. 실제 지원되는 NoC는 BUS뿐이다

`data_transfer()`는 NoC type으로 분기하지만 mesh 분기는 비어 있다([`spatial_arch.cc`](../components/spatial_arch.cc#L198-L203)). 실제 input/weight/output 요청 처리 전체가 BUS 분기 내부에만 존재한다.

함수 마지막에는 BUS가 아닌 모든 타입에 대해 `Undefined NoC Topology`로 종료한다([`spatial_arch.cc`](../components/spatial_arch.cc#L1450-L1456)). 따라서 사양 출력에서 crossbar를 인식하는 코드([`spatial_arch.cc`](../components/spatial_arch.cc#L1495-L1507))와 실제 실행 지원 범위가 다르다.

특히 제공 FSD 설정은 [`noc = store_and_forward`](../configs/accelerators/fsd.cfg#L59-L62)를 사용한다. enum에는 유효한 타입이지만 spatial 실행 경로에는 해당 분기가 없으므로 첫 PE-array 전송 시 종료한다.

### 2. BUS 경로는 topology가 아니라 direct scatter다

Dense input은 모든 active PE를 순회하면서 PE-array buffer에서 PE local buffer로 직접 복사한다([`spatial_arch.cc`](../components/spatial_arch.cc#L218-L299)). Weight와 output도 같은 구조다.

배열의 `height`, `width`, PE의 row/column, source와 destination 사이 거리, multicast route는 데이터 전달에 사용되지 않는다. 따라서 1×168, 12×14, 14×12처럼 PE 수가 같은 다른 shape를 구분하는 topology 비용이 없다.

BUS 모델로 제한하더라도 다음 항이 명시되지 않는다.

- 한 cycle에 허용되는 transaction 수
- arbitration
- multicast/broadcast group
- fanout에 따른 driver/link energy
- 동시 요청 serialization
- PE별 수신 거리

현재 cycle은 주로 unique buffer access 수와 `noc_cycle`의 곱이고, energy는 PE local access 수와 `noc_energy`의 곱이다([`spatial_arch.cc`](../components/spatial_arch.cc#L301-L329)). 이는 coarse traffic estimator로는 사용할 수 있지만 topology-accurate NoC 모델은 아니다.

### 3. Eyeriss mapping도 reduction을 요구하지만 reduction 경로가 없다

제공 Eyeriss AlexNet 첫 convolution은 `PE_Y`에 filter-height `R=11`을 배치한다([`energy.map`](../configs/mappings/eyeriss/alexnet/energy.map#L7-L10)). 즉 여러 PE가 같은 output의 서로 다른 R partial product를 계산할 수 있으므로 PE-array 또는 상위 계층에서 reduction해야 한다.

그러나 PE 결과를 PE-array output buffer에 저장하는 함수는 덧셈이 아니라 대입을 수행한다([`scheduler.cc`](../scheduler/scheduler.cc#L1611-L1644)). Spatial 클래스에도 row reduction, cluster accumulator, reduction network가 없다.

따라서 mapping factor의 곱으로 computation count를 만드는 것과 결과 tensor를 수학적으로 합치는 것이 분리되어 있다.

### 4. PE별 tensor 좌표가 하나의 좌표계에서 만들어지지 않는다

Scheduler는 input, weight, output offset vector를 각각 독립적으로 생성한다([`scheduler.cc`](../scheduler/scheduler.cc#L753-L821)). Spatial data transfer는 같은 PE index를 각 vector 크기로 별도 modulo한다.

- input: [`spatial_arch.cc`](../components/spatial_arch.cc#L226-L230)
- weight: [`spatial_arch.cc`](../components/spatial_arch.cc#L776-L780)
- output: [`spatial_arch.cc`](../components/spatial_arch.cc#L1327-L1332)

이 방식은 하나의 PE에 대응하는 `(B,K,P,Q,C,R,S)` 좌표를 먼저 만든 뒤 세 tensor 주소를 유도하지 않는다. 각 data type의 unique offset 수가 다르면 modulo 주기가 달라져, PE의 weight가 의도한 input/output과 결합된다는 보장이 없다.

가중치 offset 자체도 C/R 차원의 stride로 filter height/width가 아니라 input height/width를 사용한다([`scheduler.cc`](../scheduler/scheduler.cc#L789-L795)). Group loop 역시 source group 전체를 step으로 사용하여 group별 weight offset을 충분히 생성하지 못할 수 있다.

### 5. FUNCTIONAL output load의 batch/group index 식이 잘못되어 있다

공통 scheduler의 `output_data_load()` 목적지 index 식은 batch 항 다음에 group 항을 더하지 않고 곱한다([`scheduler.cc`](../scheduler/scheduler.cc#L1596-L1601)). `g=0`이면 batch offset까지 0으로 소거될 수 있다.

이 문제는 spatial 고유 코드는 아니지만 spatial의 output 초기 partial sum 전달이 해당 함수를 사용하므로 multi-batch/grouped convolution의 FUNCTIONAL 결과를 신뢰할 수 없게 한다.

### 6. sparse COO 분기는 비어 있지만 상태는 완료 처리된다

FUNCTIONAL input의 `SPARSE_COO` 분기는 빈 블록이다([`spatial_arch.cc`](../components/spatial_arch.cc#L333-L335)). Weight COO 분기도 동일하다([`spatial_arch.cc`](../components/spatial_arch.cc#L881-L883)).

하지만 compression 분기 이후에는 모든 PE에 대해 해당 데이터가 존재한다고 표시하고 요청을 해제한다. 즉 데이터 복사와 metadata 처리 없이 계산 가능 상태로 전환된다. COO는 명시적 unsupported 오류로 막거나 실제 format-aware 전달을 구현해야 한다.

### 7. tile-granular read 상태의 수명 주기가 불명확하다

`read_tile_granular_pe_input/weight/output`은 scheduler 생성 시 한 번 false로 초기화된다([`scheduler.cc`](../scheduler/scheduler.cc#L88-L98)). 이후 각 array 전송에서 true로 바뀌지만 layer 내 새로운 tile 또는 NPU reset 때 false로 되돌리는 경로가 없다.

따라서 vector slot을 재사용하는 후속 tile 전송은 PE-array buffer read/access가 이미 계산된 것으로 간주되어 cycle/energy를 과소 집계할 수 있다. 이 flag가 tile identity를 나타내려면 tile generation 또는 epoch과 함께 관리해야 한다.

### 8. active PE 검증은 shape를 확인하지 않는다

현재 검증은 `active_x × active_y <= num_pes`만 검사한다([`npu.cc`](../scheduler/npu.cc#L254-L268)). 다음 조건은 검사하지 않는다.

- `active_x <= physical width`
- `active_y <= physical height`

BUS 모델에서는 PE를 일렬 vector로 취급하므로 즉시 범위 밖 접근이 발생하지 않을 수 있지만, 출력되는 2차원 utilization과 향후 mesh/routing 구현의 전제가 깨진다.

## 수정 방향

### 1단계: 지원 범위와 설정 계약 정리

1. 당장은 `noc=bus`만 허용하고 다른 타입은 초기화 시 명확한 unsupported 오류로 거부한다.
2. FSD 설정을 실제 지원 topology로 변경하거나 store-and-forward 구현을 추가한다.
3. 사양 출력은 실제 선택된 실행 분기와 동일한 값을 표시한다.
4. active_x/active_y를 물리 width/height와 각각 검증한다.

### 2단계: PE 좌표와 reduction 의미 복구

1. PE마다 명시적인 `(x,y)`와 tensor coordinate를 생성한다.
2. input/weight/output offset을 동일 coordinate에서 함께 유도한다.
3. K/B/P/Q는 output identity, C/R/S는 reduction identity로 구분한다.
4. reduction 차원이 여러 PE에 배치되면 accumulator 또는 reduction network를 선택하도록 한다.
5. output STORE는 reduction 완료 후 output owner가 한 번만 수행하도록 한다.

### 3단계: BUS analytical model 구체화

BUS만 지원하더라도 다음 계약을 정의해야 한다.

- transaction width와 line size
- arbitration 및 동시 요청 처리
- broadcast/multicast latency
- fanout energy
- source buffer read 공유 조건
- local buffer write 수
- overlap 가능한 stage

현재 `read_tile_granular` boolean 대신 `(tile_id, data_offset)` 기반 중복 제거가 필요하다.

### 4단계: 추가 topology 구현

Mesh/crossbar/store-and-forward를 지원하려면 topology별 전략 클래스로 분리하는 것이 안전하다.

- route 계산
- hop 수
- link/router cycle 및 energy
- contention
- multicast replication 위치
- deadlock-free progress 조건

빈 분기를 유지하면 지원된 enum과 실제 기능을 구분할 수 없으므로 제거하거나 명시적으로 실패시켜야 한다.

## 완료 조건

- 모든 제공 accelerator config의 NoC가 초기화 단계에서 지원 여부를 명확히 판정한다.
- FSD가 실행 중간의 `Undefined NoC Topology`가 아니라 지원 topology로 동작하거나 사전 오류를 낸다.
- 각 active PE의 tensor coordinate와 세 data offset을 독립 계산기로 대조할 수 있다.
- C/R/S를 PE에 분산한 작은 convolution이 reference와 일치한다.
- multi-batch 및 grouped/depthwise convolution의 output load/store가 reference와 일치한다.
- 반복 tile에서 buffer access/cycle이 tile identity에 맞게 다시 집계된다.
- BUS fanout, transaction 수, cycle, energy가 수작업 micro-case와 일치한다.
- mesh/crossbar를 지원한다고 표시할 경우 hop 및 contention에 따라 결과가 변화한다.

## 최종 결론

현재 spatial architecture는 **BUS 기반 PE scatter와 memory traffic을 근사하는 제한적 모델**로는 의미가 있다. 하지만 일반 spatial architecture 또는 Eyeriss/Simba/FSD의 실제 데이터플로·NoC·reduction을 논리적으로 재현한다고 보기는 어렵다. 우선 BUS-only 지원 범위를 명시하고, PE 공동 좌표와 reduction 경로를 복구한 뒤 topology 확장을 진행해야 한다.
