# DRAM: timing simulation 이슈

> **재판정 (2026-09-02):** DR2의 "남음: tRC/tRAS/tRP" 목록이 stale — **JEDEC
> tRC(=tRAS+tRP) 타이밍과 `num_banks` bank 병렬화가 구현됨**(P4-1/DR2,
> `dram.h:88-103`; energy는 activation 전액, latency는 busiest-bank). SFU standalone
> softmax의 operand streaming도 동일 row 모델을 소비한다. 잔존은 선언된 모델 경계 —
> per-address bank conflict, random/strided 패턴 — 이며 `dram_timing_limits`로
> 리포트에 명시된다(idealized-lower-bound임을 출력).

## 판정

Dense INPUT/WEIGHT/OUTPUT load는 descriptor bit transaction으로 계산되지만, **OUTPUT 경로에만 sparse guard가 빠져** sparse config에서 dense 비용으로 조용히 계산될 수 있다. 또한 비-DRAMSIM3 DRAM은 row-buffer/bank/timing이 없는 flat fixed-cost 모델이다.

- 심각도: **P2**

## 확인된 문제

| # | 위치 | 문제 | 영향 |
|---|---|---|---|
| DR1 (=T8) | `components/dram.cc` (OUTPUT load) | **해결됨 (2026-08-15)** | ~~OUTPUT load 경로에만 sparse 거부 guard 없음~~ → INPUT/WEIGHT와 동일한 guard를 OUTPUT request 블록에 추가. 참고: `npu.cc:88`이 이미 최상위 init에서 non-DENSE를 거부하므로 이 guard는 INPUT/WEIGHT와 대칭인 심층 방어. 검증: dense 무회귀(36300), sparse config는 init에서 거부 |
| DR2 (=T21) | `components/dram.{cc,h}` | **부분 해결 (2026-08-15)** | ~~flat fixed-cost~~ → **open-page row-buffer 모델** 추가: dense stream이 `⌈bytes/row_buffer_size⌉` row activate, activation당 `row_miss_cycle/energy`(기본 0=off). load·write-back 공용. **남음**: random/strided 패턴, bank conflict, tRC/tRAS/tRP — cycle-accurate는 DRAMSIM3 권장 | 검증: row_buffer=1KB·miss=10 → INPUT transfer Δ=217,800=10×21,780 events(1 row/event) 수작업 일치, 기본값 무회귀 |
| DR3 (=T6/T13) | `components/dram.cc` | **해결됨 (2026-08-15)** | `frequency>0` guard + `bitwidth==0` clean error + 공용 `derived_link_bitwidth()`의 fractional 절삭 경고(B12) | frequency=0 → clean error; 327.68→327 경고 확인 |
| DR4 | `components/memory_controller.cc` (DRAMSIM3 빌드만) | **해결됨 (2026-08-15)** | ~~빈 deque `back()` UB~~ → empty guard 추가 + read/write 구분. (정정: "비인접 중복 재삽입"은 결함이 아님 — 뒤늦게 재방문하는 line은 정당한 별도 DRAM transaction이므로 인접-coalescing 의미를 유지) | 검증: DRAMSIM3 미설치로 본 빌드 검증 불가 → ASan/UBSan stub 테스트로 empty-insert·coalescing·write-구분·revisit-distinct 4ケース assert 통과 |
| DR4 (Phase-3, 신규→해결) | `stats.cc` `scale_serial_repetitions`, `mapping_table.cc` | **해결됨 (2026-08-17)** — ~~GLB 반복 스케일(×R)이 off-chip 트래픽에 일괄 적용되어, tensor가 의존하지 않는 차원의 반복(weight에 대한 Q/B 등)까지 DRAM 재인출로 계상~~ → `datatype_repetitions()`(weight K·C·R·S / input B·C·P·Q / output K·B·P·Q)로 DRAM·multi-chip 트래픽·에너지·cycle 벡터를 per-type 스케일(요청 수·leakage는 uniform 유지). busy/timeline은 보수적으로 uniform 유지(문서화) | 검증: Eyeriss conv3 DRAM weight 92MB→1.77MB(=volume fetch-once 정확), DRAM 축 33×→0.8~3.7×; latency 축 불변; phase1/2 compare·unittest 회귀 통과 |
| DR5 (Phase-3) | `mapping_table.cc` legacy GLB row | **해결됨 (2026-08-17)** — grouped conv(g>1)의 legacy행 K·GROUP 상호작용을 실험으로 확정: live 머신은 DRAM행 K 루프를 누적 GROUP으로 나누므로 **/G는 정확히 한 행에만** 적용해야 함(K가 DRAM행에 남으면 legacy g=1, K 전부 legacy면 legacy g=G). `datatype_repetitions()`에 K-보유 타입(weight/output)의 legacy-G 나눗셈 + 나눗셈 가능성 검사 추가. silicon.map conv2/4/5를 legacy K × DRAM K2(출력 writeback 트리거 유지) 패턴으로 재매핑 | 검증: conv2/4/5 comp·fold **불변**(latency ≤6.4% 유지), weight DRAM 76,800/165,888/110,592 txn = **볼륨 fetch-once 정확 일치**; DRAM 축 0.8~3.7× → **0.8/1.3/1.0/1.6/2.2×**; phase1/2·unittest 회귀 통과 |
| DR6 (Phase-3 발견) | output writeback 경계 | **해결됨 (2026-08-17)** — ~~레이어 최종 output tile flush 미계상(0.5~0.97×)~~ → `multi_chip::flush_output_writeback()`을 live 루프 종료 시 1회 호출(repetition 스케일 전이라 반복별 최종 tile도 커버; 계정 블록은 `account_output_writeback_to_dram()`으로 공용화) | 검증: traffic 하네스 T1 — eyeriss 5 conv + gemm 3점 **전부 output DRAM = 볼륨 1.000×** ([validation/traffic](../../validation/traffic/README.md)) |

## 완료 조건

- [x] OUTPUT load 경로가 INPUT/WEIGHT와 동일한 sparse 거부 guard를 갖는다(DR1). (2026-08-15)
- [~] 0 `bitwidth`와 `frequency==0`은 init에서 검증한다(완료, 2026-08-15). non-power-of-two 검증은 남음.
- [ ] (선택) row-buffer/bank/timing을 반영하는 DRAM latency 모델 또는 DRAMSIM3 경로를 명시적으로 요구한다.
- [ ] `insert_pending_queue`의 empty-deque 및 중복 삽입을 수정한다.
