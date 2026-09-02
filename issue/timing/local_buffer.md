# Local buffer: timing simulation 이슈

> **재판정 (2026-09-02):** LB7(단일/이중 버퍼 미구분)은 **해결됨** — LB7/P1-A 계약:
> `pe.cc:238`이 `double_buffer`를 파싱하고, `finalize_layer_timeline()`이 단일 버퍼면
> `compute + transfer-makespan + format-IP` 직렬, 이중이면 max로 결합하며 결과 파일에
> "PE local buffer: single/double-buffered"로 계약을 출력한다. 완료 조건의 LB4
> 체크박스도 자체 표와 어긋난 stale — `lb_static_energy`는 구현·검증됨(표 참조).
> 잔존: LB2 per-tile 세분화(현재 benign한 max-peak), LB3(FUNCTIONAL 빌드 dead-code
> 경로).

## 판정

PE 내부 local buffer(LB)의 dense access/transfer는 datatype descriptor 기반(`storage_transactions(type, elements, line_size_lb)`, packing-aware)으로 계산되어 정상이다. 그러나 **LB utilization 통계가 깨져 있고**(평균값이 항상 0), tile마다 덮어써 마지막 tile만 반영되며, single-buffer인데 full-overlap을 가정하고, LB 전용 leakage 파라미터가 없다. 그동안 전용 문서가 없었다.

- 심각도: **P2**
- 관련: [pe.md](pe.md), [runtime_datatype.md](runtime_datatype.md), [static_energy.md](static_energy.md)

## 확인된 문제

| # | 위치 | 심각도 | 문제 | 영향 |
|---|---|---|---|---|
| LB1 | `stats.cc` (print_results) | **해결됨 (2026-08-15)** | ~~`total_utilization_local_buffer`가 대입 안 됨(항상 0.0%)~~ → SEPARATE는 per-type 평균, SHARED는 합으로 계산해 출력. 검증: eyeriss LB Average 0.0%→39.5%. GLB SEPARATE의 동일 결함도 함께 수정 |
| LB2 | `pe.cc` | **해결됨 (2026-08-15, hardening)** | (정정: `update_tile_size`가 레이어당 1회 호출되어 tile 크기가 레이어 내 고정 → overwrite는 **현재 benign**) `std::max` peak-누적으로 hardening + `(float)` 절삭 제거(double). PE3(`utilization_mac`)도 동일 hardening. 향후 per-tile 변동 도입 시 안전 | eyeriss LB util 45.8/39.3/33.3/39.5 무회귀 |
| LB3 | `pe.cc` FUNCTIONAL blocks (`:756-1120`, `:2045+`, flush `:1956+`) | **P3** | FUNCTIONAL 빌드의 LB access가 host `sizeof(data_t)` 주소 walk(`mask_bits_lb`)로 packing 무시. write-back(`:2124` 등)은 LB **store**에 `u_read_cycle_lb`(방향 오류) 부과 | FUNCTIONAL 빌드에서 INT4 LB access가 4-byte stride 기준으로 집계되어 descriptor 경로와 발산. timing 빌드에선 dead code |
| LB4 | `pe.{cc,h}` | **해결됨 (2026-08-15)** | ~~LB 전용 leakage 파라미터 부재~~ → `lb_static_energy`(pJ/cycle) 추가, MAC-side `static_energy`와 분리 스윕 가능·동일 elapsed 창에 합산(리포팅 호환 유지). per-datatype 벡터 형식은 기존 키와의 일관성 위해 유지 | 검증: `lb_static_energy=0.5` 단독 → 56,678,529,600 pJ = 0.5×674,744,400(latency)×168(물리 PE) 정확 일치, 기본 무회귀 |
| LB5 | `interconnect_timing.cc` | **해결됨 (2026-08-15)** | ~~1-transaction 과다~~ → `pipelined_transfer_cycles`를 raw (source,link,dest) 인자로 리팩터, 1-transaction=source+link+dest. 공유 헬퍼 1곳 수정으로 GB7/MC5/SP3 동시 해결. unittest가 `==8`(이전 12) 직접 assert |
| LB6 | `stats.cc` | **해결됨 (2026-08-15)** | ~~`avg_access_cycle_lb[i] /= num_active_pe` 0 미방어~~ → `if(num_active_pe > 0)` guard 추가(avg_computation/mac/lb 전부, PE6과 함께) |
| LB7 | `pe.cc` (단일 `*_lb` region), `modeled_elapsed_cycles` `:2663-2675` | **P3** | LB는 **single-buffered**이나 `modeled_elapsed_cycles()`가 compute·LB access·LB↔MAC transfer를 `max`(full overlap)로 결합. utilization 분모도 전체 버퍼를 한 tile에 할당 가정 | one buffer에선 tile-N load와 tile-N compute가 겹칠 수 없는데 serialization이 숨겨짐. double-buffer 설계와 cycle·utilization상 구분 불가 |

## 정상 확인 항목

- `line_size_lb`는 bits, init에서 power-of-two ≥8 검증(div-by-zero 없음). non-FUNCTIONAL LB access는 `storage_transactions(type, elements, line_size_lb)`로 packing-aware.
- `access_cycle_lb`/`transfer_cycle`/`cycle_mac_lb`/`write_back_cycle_lb`는 `std::max`로만 결합되어 elapsed window에서 이중 계산 없음.

## 완료 조건

- [x] `total_utilization_local_buffer`가 실제 평균으로 계산된다(LB1). (2026-08-15)
- [ ] LB utilization이 tile 누적 기준으로 레이어 전체를 반영한다(LB2).
- [ ] (선택) LB 전용 leakage 파라미터를 도입한다(LB4).
- [ ] single/double buffer가 cycle·utilization에서 구분된다(LB7).
- [x] 공유 `pipelined_transfer_cycles` 1-transaction 과다를 수정한다(LB5, 횡단). (2026-08-15)
