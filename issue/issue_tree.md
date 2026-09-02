# Adder tree 논리 구현 진단

> **재판정 (2026-09-02):** 이 문서는 functional/timing 분리 이전의 원본이다.
> timing 절반은 해결 완료([timing/adder_tree.md](timing/adder_tree.md) — AT1–7,
> 잔여는 AT5 MAERI 단가 calibration 사용자 과제); functional 절반(tree reduction
> 계약 P0)은 [functional/adder_tree.md](functional/adder_tree.md)에 잔존. 현황:
> [assessment/ISSUE_REAUDIT_2026-09-02.md](../assessment/ISSUE_REAUDIT_2026-09-02.md).


## 판정

현재 `adder_tree_t`는 **reduction tree로서 핵심 동작이 구현되어 있지 않다.** 별도 클래스 이름과 사양 출력은 존재하지만 실제 실행은 입력·가중치·출력을 각 PE로 scatter하는 일반 PE array 경로이며, tree node, level, fan-in, partial-sum reduction, depth 기반 cycle/energy가 없다.

- 심각도: **P0 / functional reduction 및 timing/energy 신뢰성**
- 분석 방법: 실행 결과를 사용하지 않은 정적 코드·상태 전이·수식 분석
- 영향 범위: `maeri.cfg`, MAERI 계열 mapping, `adder_tree_t` 통계 및 FUNCTIONAL 경로

## 정상 adder tree가 만족해야 하는 조건

활성 leaf가 `N`개이고 2-input adder를 사용한다면 최소한 다음 관계가 모델에 존재해야 한다.

- 총 덧셈 수: `N - 1`
- balanced tree depth: `ceil(log2(N))`
- 각 leaf가 어느 output/reduction group에 속하는지 나타내는 좌표
- non-power-of-two leaf 처리
- level별 valid/partial sum 또는 이에 동등한 analytical latency
- adder operation과 tree link의 cycle/energy

현재 클래스에는 위 정보를 표현하는 필드나 연산이 없다.

## 확인된 문제

| 우선순위 | 문제 | 결과 |
|---|---|---|
| P0 | tree node/level/reduction 연산 부재 | partial product를 최종 output으로 합칠 수 없음 |
| P0 | output STORE가 누산이 아닌 덮어쓰기 | 동일 output의 여러 partial sum 중 마지막 값만 남을 수 있음 |
| P0 | PE별 input/weight/output 좌표 결합 부재 | reduction group 및 output owner가 보장되지 않음 |
| P0 | `is_waiting_data` 미초기화 및 해제 누락 | 범위 밖 접근과 요청 상태 고착 가능 |
| P1 | tree depth와 무관한 cycle 식 | 활성 leaf 수가 늘어도 reduction latency가 모델링되지 않음 |
| P1 | adder 및 tree link energy 부재 | 계산/통신 energy가 과소 또는 잘못 집계됨 |
| P1 | 설정은 bus지만 출력은 adder tree | NoC 설정과 실제 동작의 의미가 불일치 |
| P1 | 일반 spatial/broadcast 경로의 대량 복제 | topology별 invariant를 강제하기 어려움 |

## 상세 근거

### 1. 클래스에 tree 구조가 없다

[`adder_tree_t` 선언](../components/adder_tree.h#L7-L24)은 init, tile update, data transfer, print 함수만 가진다. 다음 상태는 존재하지 않는다.

- leaf-to-reduction-group mapping
- internal node
- tree level/depth
- fan-in
- node pipeline register
- partial-sum valid state
- root output

구현 전체에도 `log2(active_leaf)`, `N-1` adder count 또는 level 순회가 없다.

### 2. 데이터 경로는 tree가 아니라 PE별 직접 scatter다

Dense input은 모든 active PE를 순회하면서 PE-array buffer에서 각 local buffer로 직접 복사한다([`adder_tree.cc`](../components/adder_tree.cc#L202-L282)). Weight와 output도 동일한 구조를 사용한다.

전송 cycle은 memory line 접근 횟수와 `noc_cycle`의 곱으로 계산된다([`adder_tree.cc`](../components/adder_tree.cc#L284-L311)). 이 식에는 다음 요소가 없다.

- reduction group의 leaf 수
- tree depth
- level별 pipeline
- root contention
- non-power-of-two balancing

따라서 이 경로는 scatter/broadcast 비용 모델이지 reduction tree 모델이 아니다.

### 3. 제공 MAERI mapping은 실제 reduction을 요구한다

MAERI BERT mapping의 첫 layer는 `PE_X = K 8 × C 16`으로 128개 PE를 활성화한다([`cycle.map`](../configs/mappings/maeri/bert/cycle.map#L7-L10)). 각 output channel에 대해 C 방향 16개 partial product를 합쳐야 한다.

하지만 PE 결과는 PE-array output buffer로 반환될 때 덧셈 없이 대입된다([`scheduler.cc`](../scheduler/scheduler.cc#L1611-L1644)). 따라서 16-way reduction이 tree에서 수행되지도 않고 software accumulator에서 수행되지도 않는다.

### 4. PE index와 tensor 좌표가 공동으로 생성되지 않는다

Scheduler는 input, weight, output offset을 각각 별도 순회로 생성한다([`scheduler.cc`](../scheduler/scheduler.cc#L753-L821)). PE는 각 vector에 `index % vector.size()`를 독립 적용한다. 예를 들어 PE 결과 저장은 [`output_offset_pe_array[index % size]`](../components/pe.cc#L1601-L1605)을 사용한다.

MAERI 예에서 weight offset은 K×C에 대응하는 128개가 될 수 있지만 output offset은 K에 대응하는 8개다. 세 vector를 하나의 PE 좌표로 결합하지 않으므로 `(K,C)` leaf가 어느 K reduction root로 가야 하는지 보장되지 않는다.

가중치 offset의 C/R stride도 filter dimensions가 아닌 input dimensions를 사용한다([`scheduler.cc`](../scheduler/scheduler.cc#L789-L795)).

### 5. tree cycle과 energy가 없다

PE의 computation energy는 각 PE 내부 MAC 비용이며([`pe.cc`](../components/pe.cc#L2564-L2602)), tree internal adder 비용과 동일하지 않다. `adder_tree_t`에는 다음 비용 파라미터가 없다.

- adder cycle/energy
- tree node read/write
- level 간 link cycle/energy
- root accumulator cycle/energy

현재 `noc_cycle`과 `noc_energy`는 데이터 scatter에만 곱해진다. `N=1`, `N=16`, `N=168`의 reduction depth 차이를 표현할 수 없다.

### 6. 요청 상태 machine이 완결되지 않는다

`adder_tree_t::init()`은 `request_to_global_buffer`까지만 초기화하고 `is_waiting_data`를 초기화하지 않는다([`adder_tree.cc`](../components/adder_tree.cc#L116-L123)). 공통 `wait_data()`는 이 vector를 즉시 인덱싱한다([`pe_array.cc`](../components/pe_array.cc#L181-L188)).

또한 adder-tree의 `data_transfer()`는 spatial/systolic 구현과 달리 전송 후 `is_waiting_data[type]=false`로 되돌리는 처리가 없다. 초기화만 추가하더라도 첫 요청 이후 waiting 상태가 고착될 수 있다.

### 7. NoC 설정과 표시가 모순된다

제공 MAERI 설정은 [`noc = bus`](../configs/accelerators/maeri.cfg#L54-L57)다. 구현은 `noc_type`을 읽지만 `data_transfer()`에서 topology 선택에 사용하지 않는다. 반면 사양 출력은 설정과 관계없이 [`NoC type = Adder Tree`](../components/adder_tree.cc#L1466-L1471)라고 표시한다.

이는 다음 세 개념을 구분하지 못한 상태다.

1. operand distribution network
2. multiplier/PE 배열
3. partial-sum reduction network

MAERI류 모델이라면 distribution network와 reduction tree를 별도 구성요소와 비용으로 모델링해야 한다.

### 8. FUNCTIONAL output load/store도 reduction을 대체하지 못한다

Output 초기값은 PE-array buffer에서 각 PE로 직접 load된다([`adder_tree.cc`](../components/adder_tree.cc#L1309-L1316)). 여러 PE가 같은 초기 partial sum을 독립적으로 읽은 뒤 계산하고 STORE하면, STORE가 대입이므로 서로의 결과가 합쳐지지 않는다.

PE 실행 순서에 의존해 우연히 순차 누산되는 구조도 아니다. 모든 active PE는 같은 outer simulation step에서 독립적으로 실행된다([`npu.cc`](../scheduler/npu.cc#L390-L398)).

## 수정 방향

### 1단계: 컴포넌트 의미 분리

다음 세 계층을 명시적으로 나눈다.

1. operand distribution network: bus, tree, crossbar 등
2. multiplier/PE leaf
3. reduction network: adder tree와 root accumulator

`noc = bus`와 `reduction = adder_tree`처럼 설정 계약도 분리하는 것이 적절하다.

### 2단계: reduction group 생성

1. mapping의 B/K/P/Q를 output identity로 정의한다.
2. C/R/S 및 분할된 MAC lane을 reduction identity로 정의한다.
3. 각 PE에 `(output_id, reduction_id)`를 직접 부여한다.
4. 같은 `output_id`에 속한 leaf만 하나의 tree에 연결한다.
5. group별 실제 leaf 수와 non-power-of-two padding을 계산한다.

### 3단계: tree 기능 및 비용 구현

Cycle-accurate 구현이라면 level별 partial-sum state를 둔다. Analytical 구현이라면 최소한 다음 식을 명시적으로 사용한다.

- `num_additions = max(0, leaves - 1)`
- `depth = ceil(log2(leaves))`
- pipelined/unpipelined latency 구분
- level별 link traversal 수
- adder operation energy
- root write energy

FUNCTIONAL STORE는 reduction 결과가 root에 도달했을 때 한 번만 수행되어야 한다.

### 4단계: 상태 및 공통 scheduler 수정

1. `is_waiting_data`를 base class에서 공통 초기화·reset한다.
2. 전송 완료 시 해당 waiting flag를 해제한다.
3. 독립 modulo offset을 제거하고 PE 좌표에서 세 data offset을 함께 생성한다.
4. weight offset의 filter stride와 group 순회를 수정한다.

## 완료 조건

- 1, 2, 3, 5, 16, 168 leaf에 대해 addition count와 depth가 수작업 계산과 일치한다.
- K가 여러 개인 mapping에서 서로 다른 output group의 partial sum이 섞이지 않는다.
- C/R/S를 PE에 분산한 작은 conv/GEMM 결과가 reference와 일치한다.
- FUNCTIONAL 결과가 PE 실행 순서와 무관하다.
- reduction cycle이 tree depth 및 pipeline 설정에 따라 변한다.
- reduction energy가 `N-1` adder operation과 level별 link 수에 따라 변한다.
- 출력된 distribution/reduction topology가 실제 설정 및 실행 경로와 일치한다.

## 최종 결론

현재 `adder_tree_t`는 **이름이 다른 generic PE-array scatter 구현**이다. MAERI와 같은 reduction-centric architecture를 표현하려면 partial-sum group, internal node, tree depth, root output을 새로 모델링해야 한다.
