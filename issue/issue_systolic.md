# Systolic array 논리 구현 진단

> **재판정 (2026-09-02):** 이 문서는 functional/timing 분리 이전의 원본이다.
> timing 절반은 해결 완료([timing/systolic_array.md](timing/systolic_array.md) —
> SY1–5·V1–3 + Gemmini RTL MAPE 4.4%); functional 절반(wavefront propagation P0)은
> [functional/systolic_array.md](functional/systolic_array.md)에 잔존(단, FUNCTINOAL
> 오타·tpu.cfg `pe_stationary` 키는 수정 확인됨). 현황:
> [assessment/ISSUE_REAUDIT_2026-09-02.md](../assessment/ISSUE_REAUDIT_2026-09-02.md).


## 판정

현재 `systolic_array_t`는 별도 컴포넌트와 일부 행·열 방향 전송 통계를 갖고 있지만, **논리적으로 유효한 systolic array 구현으로 볼 수 없다.** 제공 TPU 설정에서는 stationary type이 정상적으로 설정되지 않으며, 이를 수정하더라도 실제 데이터 경로와 cycle 모델에는 PE 간 wavefront, skew, fill/drain, partial-sum 전달이 없다.

- 심각도: **P0 / 결과 신뢰성 및 실행 안전성**
- 분석 방법: 실행 결과를 사용하지 않은 정적 코드·상태 전이·인덱스 분석
- 영향 범위: `tpu.cfg`, `tpuv3.cfg`, systolic timing/energy 통계, FUNCTIONAL 경로

## 정상 systolic array가 만족해야 하는 조건

1. 입력과 가중치는 지정된 배열 경계에 주입되어야 한다.
2. 데이터는 매 cycle 인접 PE로 한 hop씩 이동해야 한다.
3. PE의 `(row, column)`과 매핑된 tensor 좌표가 일관되어야 한다.
4. partial sum은 PE 내부 또는 지정된 배열 방향으로 전달·누산되어야 한다.
5. latency에는 데이터 skew, pipeline fill, steady state, drain이 반영되어야 한다.
6. energy에는 실제 hop 수와 활성 링크 수가 반영되어야 한다.

현재 구현은 위 조건 중 PE 생성과 일부 방향별 통계 코드만 부분적으로 갖고 있다.

## 확인된 문제

| 우선순위 | 문제 | 결과 |
|---|---|---|
| P0 | 설정 키 불일치로 stationary type이 undefined | 제공 TPU 설정에서 핵심 전송 분기가 실행되지 않음 |
| P0 | `is_waiting_data` 미초기화 | 공통 상태 확인에서 범위 밖 접근 |
| P0 | physical width와 active width 혼용 | partial-width 매핑에서 unsigned underflow와 잘못된 PE 접근 |
| P0 | PE 간 partial-sum 누산 부재 | reduction 차원을 PE에 분산하면 결과를 덮어씀 |
| P0 | FUNCTIONAL 전처리기 오타 | input-stationary weight가 복사되지 않음 |
| P1 | 실제 이동은 PE-array에서 각 PE로 직접 복사 | systolic wavefront를 재현하지 못함 |
| P1 | fill/drain 및 hop latency 부재 | 배열 크기에 따른 cycle 변화가 구조적으로 누락됨 |
| P1 | NoC type이 동작을 선택하지 않음 | 설정과 모델 사이의 계약이 없음 |
| P1 | 통계 및 skip flag 오기록 | 데이터별 cycle/재사용 통계 왜곡 |

## 상세 근거

### 1. 제공 설정과 stationary key가 일치하지 않는다

구현은 [`pe_array_stationary`](../components/systolic_array.cc#L53-L57)를 읽는다. 그러나 제공 TPU 설정은 [`pe_stationary`](../configs/accelerators/tpu.cfg#L43-L45), TPUv3도 동일한 키를 사용한다.

기본 클래스의 초기값은 [`UNDEFINED_STATIONARY`](../components/pe_array.cc#L18-L21)이다. 따라서 제공 설정을 사용하면 다음 세 분기 중 어느 것도 선택되지 않는다.

- input-stationary 입력 전송: [`systolic_array.cc`](../components/systolic_array.cc#L206)
- weight/output-stationary 입력 전송: [`systolic_array.cc`](../components/systolic_array.cc#L324)
- weight 및 output의 대응 분기: [`systolic_array.cc`](../components/systolic_array.cc#L392-L602)

그럼에도 함수는 PE에 실제 데이터를 전달하지 않은 상태에서 모든 PE를 `exist_data_lb=true`로 표시하고 요청을 해제한다([`systolic_array.cc`](../components/systolic_array.cc#L374-L377)). undefined 값은 사양 출력에서도 별도 오류가 아니라 output-stationary로 표시된다([`systolic_array.cc`](../components/systolic_array.cc#L826-L837)).

### 2. 요청 상태 벡터가 초기화되지 않는다

`systolic_array_t::init()`은 `exist_data`와 `request_to_global_buffer`는 초기화하지만 `is_waiting_data`를 초기화하지 않는다([`systolic_array.cc`](../components/systolic_array.cc#L124-L132)). 공통 `wait_data()`는 크기 확인 없이 세 원소를 인덱싱한다([`pe_array.cc`](../components/pe_array.cc#L181-L188)).

따라서 정상적인 첫 데이터 요청 경로에서도 정의되지 않은 동작이 발생할 수 있다. `cycle_temporal_pe` 역시 systolic 클래스에서는 초기화되지 않지만 여러 전송 분기에서 인덱싱된다.

### 3. 2차원 PE 인덱싱이 일관되지 않는다

세로 전달 경계는 `i / num_active_pe_x == 0`으로 판별하지만 이전 행 PE는 `pes[i-width]`로 접근한다([`systolic_array.cc`](../components/systolic_array.cc#L412-L423)). 활성 PE는 물리 배열의 row stride를 보존하지 않고 vector 앞부분에 조밀하게 배치된다.

예를 들어 physical width가 256이고 active width가 192이면 두 번째 active row의 첫 PE는 `i=192`이다. 경계 판정상 두 번째 행으로 처리되지만 `i-width`는 unsigned underflow가 된다. 제공 TPU 매핑에는 `PE_X` active count가 192인 항목이 존재하므로([`cycle.map`](../configs/mappings/tpu/vit/cycle.map#L35-L36)) 실제 제공 구성과도 충돌한다.

같은 문제가 output 전달의 [`pes[i-width]`](../components/systolic_array.cc#L642-L653)에도 존재한다.

### 4. FUNCTIONAL 데이터는 인접 PE를 통과하지 않는다

입력과 가중치는 각 PE에 대해 scheduler의 `transfer_data()`를 직접 호출하여 PE-array buffer에서 local buffer로 복사한다. 예시는 다음과 같다.

- 입력 직접 복사: [`systolic_array.cc`](../components/systolic_array.cc#L332-L336)
- weight-stationary 가중치 직접 복사: [`systolic_array.cc`](../components/systolic_array.cc#L447-L451)

이후 통계만 `pes[i-1]` 또는 `pes[i-width]`가 읽은 것처럼 기록한다. 실제 데이터 상태에는 다음 정보가 없다.

- 현재 cycle
- PE별 input/weight/partial-sum pipeline register
- 유효 데이터의 row/column 위치
- injection skew
- backpressure 또는 link contention

따라서 FUNCTIONAL 결과가 맞더라도 그것은 systolic 이동을 실행한 결과가 아니라 최종 PE 데이터를 직접 배치한 결과다.

### 5. partial sum 전달 및 reduction이 없다

각 PE는 독립적으로 [`mac_operation()`](../components/pe.cc#L2118-L2127)을 수행한다. PE가 local output을 PE-array output buffer에 반환할 때 scheduler는 목적지 값에 더하지 않고 단순 대입한다([`scheduler.cc`](../scheduler/scheduler.cc#L1611-L1644)).

즉 C/R/S 같은 reduction 차원이 여러 PE에 분산되면 동일 output에 대한 partial sum을 전달하거나 합치는 경로가 없다. 마지막 STORE가 앞선 PE의 값을 덮어쓸 수 있다.

### 6. cycle/energy 모델이 systolic geometry를 반영하지 않는다

row/column 전달 분기의 transfer cycle은 대체로 한 번의 `noc_cycle × line/bitwidth`로 계산된다. 예를 들어 입력의 수평 전달은 [`systolic_array.cc`](../components/systolic_array.cc#L361-L370)에 배열 폭이나 목적 PE까지의 거리 항이 없다.

최종 computation cycle도 모든 PE의 cycle 중 최댓값으로 집계된다([`stats.cc`](../scheduler/stats.cc#L356-L370)). 별도의 systolic fill/drain 항이 없다. 따라서 active height/width가 바뀌어도 예상되는 wavefront latency가 보장되지 않는다.

### 7. 명확한 구현 오타와 통계 오류가 있다

- [`#ifdef FUNCTINOAL`](../components/systolic_array.cc#L399-L410): `FUNCTIONAL` 빌드에서도 input-stationary weight 복사가 제외된다.
- [`skip_transfer[OUTPUT] = false`](../components/systolic_array.cc#L577-L591): weight 전달 후 weight가 아닌 output flag를 변경한다.
- [`cycle_temporal_pe[INPUT]`](../components/systolic_array.cc#L543-L554): weight-stationary weight pipeline cycle 일부를 input 통계에 더한다.
- `noc_type`은 파싱하지만 동작 선택에는 사용하지 않으며 사양 출력은 항상 Store-and-Forward다([`systolic_array.cc`](../components/systolic_array.cc#L841-L846)).

### 8. 공통 scheduler의 offset 생성도 좌표 보존을 보장하지 않는다

입력, 가중치, 출력 offset은 서로 독립된 vector로 생성된다([`scheduler.cc`](../scheduler/scheduler.cc#L753-L821)). 각 PE는 자신의 index를 각 vector 크기로 별도 modulo하여 사용한다. 하나의 공통 `(B,K,P,Q,C,R,S)` 좌표에서 세 offset을 동시에 유도하지 않으므로 vector 길이가 다르면 동일 PE의 input/weight/output 좌표 결합이 보장되지 않는다.

가중치 offset 계산에는 `c`와 `r`의 stride로 filter shape가 아니라 input height/width를 사용하는 오류도 있다([`scheduler.cc`](../scheduler/scheduler.cc#L789-L795)).

## 수정 방향

### 1단계: 실행 안전성과 설정 계약 복구

1. `pe_stationary`를 공통 키로 사용하거나 legacy alias를 명시적으로 지원한다.
2. 필수 stationary/NoC 값을 검증하고 undefined를 허용하지 않는다.
3. `is_waiting_data`와 `cycle_temporal_pe`를 base helper에서 공통 초기화·reset한다.
4. `active_x <= width`, `active_y <= height`를 각각 검증한다.
5. PE 좌표 변환 함수를 만들고 physical index와 compact logical index 중 하나만 사용한다.
6. 전처리기 오타, skip flag, 잘못된 통계 data type을 수정한다.

### 2단계: systolic 데이터 경로 구현

1. PE마다 north/south/east/west pipeline register와 valid bit를 둔다.
2. stationary mode별 injection edge와 propagation direction을 명시한다.
3. cycle마다 한 hop만 이동하도록 상태 전이를 구현한다.
4. partial sum 전달 방향과 accumulation owner를 정의한다.
5. non-square 및 partially active array의 boundary를 처리한다.

### 3단계: 분석 모델 구현

cycle-accurate 모델이 목적이 아니라면 최소한 다음 항을 갖는 검증 가능한 analytical model이 필요하다.

- injection bandwidth
- operand skew
- active row/column 수
- pipeline fill/drain
- per-hop link latency/energy
- partial-sum reduction/propagation latency
- 타일 경계의 bubble 및 재주입 비용

## 완료 조건

- TPU/TPUv3 설정을 읽은 직후 stationary type이 weight-stationary임을 검증한다.
- 2×2, 2×3, partial-width array에서 모든 predecessor index가 유효하다.
- 작은 GEMM에서 cycle별 operand 위치와 partial sum의 예상 trace가 일치한다.
- 결과가 독립 reference GEMM/conv와 일치하며 STORE 순서에 영향을 받지 않는다.
- active dimension을 바꾸면 이론적인 fill/drain 및 hop latency만큼 cycle이 변한다.
- 각 link traversal 수로 계산한 energy가 수작업 계산과 일치한다.

## 최종 결론

현재 코드는 systolic 전용 데이터 경로라기보다 **PE-array 직접 배치와 일부 인접 전달 통계를 결합한 초기 모델**이다. 설정 및 메모리 안전 문제를 먼저 해결한 뒤, PE 좌표·cycle 상태·partial-sum 경로를 새로 정의해야 한다.
