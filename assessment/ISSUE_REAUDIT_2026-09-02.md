# issue/ 문서 재감사 (2026-09-02)

`issue/` 트리의 판정을 **현재 코드 상태**(HEAD, SFU Phase 0–8 + 2026-08-15~28의 timing/
fidelity 작업 반영)와 대조한 결과. 각 문서에는 동일 날짜의 "재판정" 섹션을 추가했다.
여기서는 마스터 표만 유지한다. 검증 방법: 각 주장에 대해 코드 grep + 해당 validation
스위트 실행 결과(2026-08-28~09-02 전부 통과 상태)를 대조했다.

## 요약

- **stale(문서가 뒤처짐) → 재판정으로 해소**: 8개 문서의 판정 헤더가 이후 구현을
  반영하지 못하고 있었다. 대표: global_cycle_overlap("전역 clock 없음" → per-tile
  timeline까지 구현·검증됨), pe("mac_width>1 신뢰 불가" → lane 계약 구현·검증됨),
  systolic("flat scatter로 해석" → fill/hop + RTL 4.4% 검증).
- **실제 미해결로 남는 것**: sparse 전 계층(functional 선행 조건, 명시적 거부로 차단),
  functional low-precision 수학(runtime_datatype P0), functional
  adder-tree/spatial/systolic의 reduction·propagation 계약(P0, timing과 무관),
  LB3/GB6(FUNCTIONAL 빌드 dead-code 경로), bank-conflict/datatype-arbitration
  (선언된 모델 경계), AT5(MAERI 단가 calibration — 사용자 과제).

## 마스터 표

| 문서 | 이전 판정 | 재판정 (2026-09-02) |
|---|---|---|
| timing/global_cycle_overlap | stage-granular timeline까지; per-tile은 후속 | **종결** — per-tile event timeline(L5 `pipeline_timeline_cycles`) + back-pressure stall + boundary depth(+MC2 outstanding) + SFU 6번째 stage까지 구현, unittest hand-calc 고정 |
| timing/cycle_engine | CE1–7 해결 | **현행 유지** (문서 정확) |
| timing/static_energy | SE1–4 해결 | **현행 유지**; SFU leakage(always-on, 최종 window)도 동일 계약으로 합류 |
| timing/runtime_datatype | dense 완료, sparse 거부 | **현행 유지** |
| timing/pe | PE2 부분·판정 "mac_width>1 신뢰 불가" | **판정 갱신** — L9 active-lane 과금(`lane_state` 기반 fill/energy)으로 PE2 잔여 해소, `mac_energy_<in>_<wt>` precision family(E3)·reduction energy 축(E4) 추가. mac_width>1은 이제 계약·검증 존재 |
| timing/pe_array | PA1–9 해결(체크박스 일부 stale) | **체크박스 정정** + stripe-transition 보정(2026-08-26, nv_small RTL) 참조 추가 |
| timing/local_buffer | LB7(단일/이중 구분) 미해결 | **LB7 해결** — LB7/P1-A: 단일 버퍼는 `compute + transfer-makespan + format` 직렬, 이중은 max; `pe_double_buffer`가 timeline 계약. LB4 체크박스 정정. LB3(FUNCTIONAL dead-code)만 잔존 |
| timing/global_buffer | bypass P4 미적용, contention 남음 | **bypass 해결** — P1-B: bypass 스트림은 GLB SRAM 소스 생략(NoP ingress 직결, packet `without_source`), boundary-depth L6까지 반영. contention 잔여는 bank-conflict/datatype-arbitration으로 축소(CE4 shared-sum + num_banks + per-tile back-pressure 반영됨) |
| timing/multi_chip | multicast fanout 미구현, back-pressure 남음 | **multicast 해결** — L7 `nop_multicast`(broadcast 1회 스트림) + RE6 link contract 보고. back-pressure는 per-tile timeline + MC2 outstanding으로 부분 반영; chip-간 arbitration 세분화만 잔존 |
| timing/dram | DR2 "tRC/tRAS/tRP 남음" | **DR2 갱신** — tRAS+tRP(tRC) + `num_banks` 구현(P4-1), SFU softmax streaming도 동일 모델 소비. 잔여: per-address bank conflict/random·strided(선언된 경계, `dram_timing_limits`로 보고) |
| timing/spatial_architecture | "잔여 이슈" 헤더 | **헤더 정정** — SP1–3 해결 완료, 잔여는 모델 경계(router queue/VC)뿐. NVDLA nv_small(spatial 경로) holdout 0.55% 외부 검증 추가 |
| timing/systolic_array | 판정 "flat scatter", SY1/SY2 미해결 표기 | **판정 갱신** — SY1–5·V1–3 해결(자체 표와 일치), Gemmini RTL MAPE 4.4%로 외부 검증. '우려' 노트의 "SY1/SY2 여전히 미해결" 문구는 수정 이전 기록 |
| timing/adder_tree | AT5 이관, PA 공통 참조 | **참조 정정** — PA2–4 해결 반영. AT5(MAERI 단가)는 사용자 calibration 과제 유지 |
| issue_pe (top) | P0 activation 미호출, P1 mac_width 등 | **P0/P1 대부분 해소** — activation ownership 계약 확정(functional=Nebula `forward()` 단독, cost=SFU, `pe_t::activation` 선언째 제거+grep 가드; plan_sfu.md), mac_width=PE1/L9, 절삭=PE5/B12, static=SE1–4. 잔존: sparse·PE-local shared/bypass(명시적 미지원 유지) |
| functional/pe | P0 activation, P1 sparse | **P0 해소**(위와 동일 근거; `activation(sum(products))` reference는 Nebula가 충족, SFU evaluator와 교차 비교됨). sparse는 미지원 정책 유지 |
| functional/runtime_datatype | P0 (값이 format을 안 따름) | **미해결 유지** — 저장/회계만 format-aware; functional quantize/dequantize/MAC 수학 미구현 (유일하게 남은 P0급) |
| functional/adder_tree · spatial · systolic | P0 reduction/propagation | **대체로 잔존** — timing과 분리된 functional 계약 이슈. 단, `FUNCTINOAL` 오타는 수정 확인, tpu.cfg `pe_stationary` 키 정정 확인. 상세 재검증은 FUNCTIONAL 빌드 회귀와 함께 후속 |
| frontend/pytorch_vllm | — | **본 감사 범위 외** — executable-IR/workload-graph 작업이 병행 세션에서 진행 중 |

## 검증 기준선 (재감사 시점)

Gemmini RTL 4.40% / Eyeriss silicon 4.26% / NVDLA holdout 0.55%·traffic EXACT /
SCALE-Sim ρ 0.929 / SFU: nv_small X1 fill 151·II 1, nv_small_256_full LUT fill 156·
II 1(전부 spec 교차검증) / 전체 unittest·validation 게이트 통과.
