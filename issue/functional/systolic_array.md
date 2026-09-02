# Systolic array: functional simulation 이슈

## 판정

현재 systolic 구현은 PE-array buffer에서 각 PE local buffer로 직접 데이터를 복사한다. boundary injection, per-cycle neighbor propagation, partial-sum 전달이 없으며 제공 TPU config의 stationary key도 읽지 못한다. 따라서 systolic functional 결과는 실제 systolic 데이터 경로를 실행한 결과가 아니다.

- 심각도: **P0**

## 확인된 문제

| 문제 | 영향 |
|---|---|
| `pe_array_stationary`와 제공 config의 `pe_stationary`가 불일치 | 핵심 transfer 분기가 실행되지 않음 |
| `is_waiting_data`가 초기화되지 않음 | 첫 요청에서도 범위 밖 접근 가능 |
| physical width와 active width를 혼용 | partial-width mapping에서 predecessor index underflow |
| input/weight를 PE에 직접 복사 | operand가 인접 PE를 한 hop씩 통과하지 않음 |
| PE 간 partial sum 전달·누산 없음 | reduction을 PE에 분산하면 결과가 덮어써짐 |
| `#ifdef FUNCTINOAL` 오타 | input-stationary weight copy가 functional build에서 제외됨 |

## 수정 방향

1. `pe_stationary`를 공통 key로 읽고 undefined stationary/NoC를 fail-fast한다.
2. `is_waiting_data`와 공통 vector를 base helper에서 init/reset한다.
3. physical 또는 compact logical index 중 하나를 선택하고 좌표 변환 API로 predecessor를 계산한다.
4. PE마다 operand/partial-sum pipeline register와 valid bit를 두고, stationary mode별 boundary injection과 한-cycle one-hop propagation을 구현한다.
5. output identity와 reduction identity를 분리해 accumulation owner와 final STORE 시점을 고정한다.

## 완료 조건

- TPU/TPUv3 config가 weight-stationary로 파싱되고 undefined를 허용하지 않는다.
- 2×2, 2×3, partial-width array의 모든 predecessor index가 유효하다.
- 작은 GEMM/conv에서 cycle별 operand 위치와 partial sum trace가 reference와 일치한다.
