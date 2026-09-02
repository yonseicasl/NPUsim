# Adder tree: functional simulation 이슈

## 판정

`adder_tree_t`는 실제 reduction tree가 아니라 PE별 input/weight/output direct scatter를 수행한다. tree node, reduction group, partial-sum state, root output이 없으므로 MAERI처럼 C/R/S를 leaf에 분산하는 mapping의 functional 결과는 신뢰할 수 없다.

- 심각도: **P0**

## 확인된 문제

| 문제 | 영향 |
|---|---|
| tree node/level/fan-in/reduction 연산이 없음 | leaf partial product를 final output으로 합칠 수 없음 |
| output STORE가 대입 | 동일 output의 마지막 partial sum만 남을 수 있음 |
| PE별 input/weight/output 좌표 결합이 없음 | leaf가 어느 output reduction group에 속하는지 보장되지 않음 |
| `is_waiting_data` 미초기화 및 해제 누락 | 범위 밖 접근 또는 요청 상태 고착 가능 |

## 근거와 수정 방향

MAERI BERT mapping은 `K 8 × C 16`의 128 PE를 활성화해 K별 C reduction을 요구한다. 그러나 `adder_tree_t`에는 leaf-to-output mapping이나 internal node 상태가 없고, [`scheduler.cc`](../../scheduler/scheduler.cc#L1611-L1644)는 output을 누산하지 않는다.

다음 순서로 functional 경로를 만든다.

1. B/K/P/Q를 output identity, C/R/S 및 분할 MAC lane을 reduction identity로 정의한다.
2. 각 PE에 `(output_id, reduction_id)`를 직접 부여하고, 같은 output leaf만 tree에 연결한다.
3. level별 valid partial sum 또는 동등한 reduction state를 두고 non-power-of-two leaf를 처리한다.
4. root에 도달한 reduction 결과만 한 번 STORE하며, waiting state는 base class에서 init/reset한다.

## 완료 조건

- K가 여러 개인 mapping에서 서로 다른 output group의 partial sum이 섞이지 않는다.
- 1, 2, 3, 5, 16, 168 leaf의 small conv/GEMM 결과가 reference와 일치한다.
- functional 결과가 PE 실행 순서와 무관하다.
- 전송 완료 뒤 waiting flag가 해제되고 다음 tile 요청이 진행된다.
