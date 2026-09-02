# Timing validation contract

## 상태

Dense Conv/FC latency의 외부 기준과 현재 NPUsim counter 기준을 CSV로 분리해
고정한다. RTL 또는 실리콘 측정치는 변경하지 않으며, 의도적인 simulator 모델
변경만 NPUsim baseline 갱신을 허용한다.

## Acceptance gate

| 경로 | 기준 | 완료 조건 |
|---|---|---|
| Gemmini systolic/WS | Verilator RTL GEMM 6점 | MAPE ≤ 8%, 각 점 ≤ 8% |
| Eyeriss spatial/RS | JSSC 2017 silicon Conv 5층 | latency MAPE ≤ 5%, 각 층 ≤ 8% |
| Eyeriss memory traffic | JSSC 2017 silicon Conv 5층 | 각 GLB·DRAM 축 ≤ 15% (traffic milestone부터 gate) |

`validation/check_timing.py`는 외부 정확도와 현재 counter 회귀를 함께 검사한다.
`unittest/run_timing_validation.sh`는 모든 대상 workload를 다시 실행한 후 이
검사를 호출한다. RTL 시뮬레이션 자체는 장시간 실행이므로 매 CI에서 반복하지
않고 저장된 cycle-exact 측정치를 비교 기준으로 사용한다.

현재 traffic 오차는 알려진 미완료 항목이므로 출력만 하고 실패시키지 않는다.
`--check-traffic`은 partial-sum residency와 DRAM reuse 구현이 끝난 뒤 CI gate로
활성화한다.
