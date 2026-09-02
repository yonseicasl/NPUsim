# PE: functional simulation 이슈

> **재판정 (2026-09-02):** P0 activation은 **해결됨** — functional activation의 유일
> owner는 Nebula `forward()`로 확정(plan/plan_sfu.md functional 정책; FUNCTIONAL
> 빌드에서 layer마다 호출되어 `activation(sum(products))` reference를 충족), PE의
> 하드코딩 ReLU는 선언째 제거 + grep 가드. `linear`에 activation이 적용되지 않는
> 것과 activation 중복 적용 없음도 같은 계약으로 보장. SFU evaluator ↔ 실제 Nebula
> 구현의 직접 비교가 validation/sfu에 상주. **잔존:** sparse(P1 — dense-copy+cost
> 분기, 전 계층 명시적 거부로 오계산 차단), PE-local shared/bypass(P2, 미지원 유지).

## 판정

기본 dense tile reduction의 배열 범위와 reset은 개선됐지만, ReLU activation과 sparse execution의 functional 계약이 빠져 있다. 따라서 ReLU layer 및 supported라고 표시되는 sparse functional 결과를 신뢰할 수 없다.

- 심각도: **P0 activation**, **P1 sparse execution**

## 확인된 문제

| 우선순위 | 문제 | 영향 |
|---|---|---|
| P0 | `activation()`이 구현만 되고 호출되지 않음 | `activation=relu`의 convolution/FC가 linear 결과를 냄 |
| P1 | sparse functional path가 dense copy 후 cost 분기만 수행 | compressed execution 결과가 아님 |
| P2 | PE-local `shared`/`bypass`가 미구현 | 해당 구성을 실행할 수 없음 |

## 근거와 수정 방향

- [`pe.cc`](../../components/pe.cc#L2334-L2338)는 accumulator에 ReLU를 적용하지만 stationary compute/store 경로에서 호출되지 않는다. partial sum마다 적용하지 않고, 하나의 output identity의 C/R/S reduction이 끝난 뒤 최종 commit 직전에 layer policy를 한 번 적용해야 한다.
- FUNCTIONAL input transfer는 dense tensor를 MAC register로 복사한 뒤 compression type에 따라 비용만 분기한다([`pe.cc`](../../components/pe.cc#L633-L639), [`pe.cc`](../../components/pe.cc#L743-L770)). COO는 명시적으로 미지원이지만 다른 sparse 결과도 compressed value/index/pointer, decoder, nonzero matching, irregular issue를 실행하지 않는다.
- `shared`/`bypass`는 현재 초기화에서 명시적으로 거부한다([`pe.cc`](../../components/pe.cc#L242-L246), [`pe.cc`](../../components/pe.cc#L314-L325)). 이 상태는 유지하되 설정 문서에서 unsupported로 표시하거나, shared allocation/live-range 및 bypass forwarding path를 구현한다.

## 완료 조건

- 음수·양수 product가 섞인 convolution이 `activation(sum(products))` reference와 일치한다.
- `linear` layer에는 activation이 적용되지 않는다.
- input/weight/output-stationary가 동일 layer에 대해 같은 최종 tensor를 만든다.
- 지원 sparse format은 compressed data path로 functional event를 실행하며, 미지원 format은 초기화에서 거부한다.
- shared/bypass의 지원 범위가 설정 검증과 실제 functional path에서 일치한다.
