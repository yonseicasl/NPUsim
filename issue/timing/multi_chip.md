# Multi-chip / NoP: timing simulation 이슈

> **재판정 (2026-09-02):** 남은 P2 중 **multicast fanout은 해결됨** — L7
> `nop_multicast`(broadcast tile을 공유 매체에 1회 스트림; `multi_chip.cc:507`)와
> RE6 `nop_link_contract` 보고, `gemmini_nop*_mc/uc` fixture로 검증. back-pressure는
> per-tile timeline(L5) + MC2 `max_outstanding_requests`(DRAM 경계 depth)로 부분
> 반영됐다. 잔존: chip-간 링크 arbitration/contention 세분화, sparse stream.

## 판정

Multi-chip의 dense INPUT/WEIGHT/OUTPUT **load** 경계는 runtime datatype descriptor로 전환됐다. 그러나 GLB→multi-chip **OUTPUT write-back의 link 비용이 주석 처리되어 `multi_chip->write_back_cycle`이 0**으로 남고, chip grid/route/contention은 반영되지 않으며, 패키지 temporal buffer의 static energy도 미집계다. 이전 "dense 완료" 표기는 write-back 방향에서 사실과 다르다.

- datatype 상태: **dense load·OUTPUT write-back 완료** (MC1/MC3 해결, 2026-08-15)
- 남은 심각도: **P2 topology·contention 부재** (route/hop/fanout/back-pressure)

## ⚠️ 완료 오표기 정정

| # | 위치 | 상태 | 문제 | 영향 |
|---|---|---|---|---|
| MC1 (=T1) | `global_buffer.cc` `flush_data`, `multi_chip` 대상 | **해결됨 (2026-08-15)** | ~~GLB→multi-chip output store의 NoP 전송 비용이 0~~ → `multi_chip->transfer_cycle/energy[OUTPUT]` + link-transaction 카운터가 채워짐. 상세·검증은 [global_buffer.md](global_buffer.md) GB1 | OUTPUT write-back NoP traffic 집계됨 |
| MC2 (=SE3) | `multi_chip.{cc,h}`, `stats.{cc,h}` | **해결됨 (2026-08-15)** | ~~패키지 temporal buffer의 static energy 미집계~~ → `u_static_energy` 파싱 + `update_static_energy()` + stats set/accumulate/scale/print 추가 | `[multi_chip] static_energy` 반영. [static_energy.md](static_energy.md) SE3 |

## 확인된 문제

| 우선순위 | 위치 | 문제 | 영향 |
|---|---|---|---|
| ~~P2~~ 부분 해결 (=T9) | route/hop | **routed-unicast mesh 구현 (2026-08-15)**: `nop=mesh` 수용, chip별 ingress(0,0)로부터 Manhattan hop — energy는 per-hop×per-txn, latency는 SP1 pipelined-hop fill. 검증: unittest가 1×4 vs 2×2에서 chip3 hop 3 vs 2를 assert(1×N≠M×N), 1×1 grid는 bus와 동일(130,680), store_and_forward 거부. **남음**: arbitration/link contention/back-pressure(P2, A1 per-tile timeline 필요) |
| P2 (=T9) | `multi_chip.cc` distribution | broadcast/multicast fanout, shared source-read 없음(PE-array PA1과 동형의 per-chip 복제) — chip 간 공유 operand의 multicast 모델 미구현 | 여러 chip에 같은 input/weight 공급 시 source read N배, NoP traffic N배로 과다 |
| MC6 (신규) | `multi_chip.cc` init | **해결됨 (2026-08-15)** — NoP 단가 키 불일치: 코드는 `nop_cycle/nop_energy`만 읽는데 **배포된 7개 config 전부 `noc_cycle/noc_energy`** 사용 → 모든 run에서 NoP 단가 silent 0(Interconnection cycle 0.0의 원인). `noc_*` fallback 추가 + **미초기화 `u_transfer_cycle/energy`**(GB1 writeback이 사용, UB)를 ctor 초기화 후 nop 단가로 설정 | 검증: Interconnection cycle 0.0→INPUT 130,680/WEIGHT 2,635,380/OUTPUT 36,300 = link-txn 카운터와 정확 일치 |
| ~~P3~~ 완료 (=SE2) | `stats.cc` | ~~leakage 도메인 불일치~~ — always-on으로 통일 해결(2026-08-15) | [static_energy.md](static_energy.md) SE2 |

## 추가 확인된 문제 (2026-08-15 재감사, 미문서화)

| # | 위치 | 심각도 | 문제 | 영향 |
|---|---|---|---|---|
| MC3 | `multi_chip.cc` `request_data` | **해결됨 (2026-08-15)** | ~~DRAM write-back 비용이 `!initial && !equal_output_tile` guard 밖에서 무조건 부과~~ → num_request + timing + 모든 cost 부과를 guard 안으로 이동(실제 write-back일 때만 과금). 검증: eyeriss layer_0 DRAM Output access cycle 3,484,800→3,405,600(spurious pass 제거) |
| MC4 | write-back double_buffer gating | **해결됨 (2026-08-15, 관례 명문화)** | 코드 관례는 "S→D 전송의 access cycle을 **destination**의 double_buffer로 gate(energy 항상)"이며 load(DRAM→MC, MC→GLB)에 이미 적용돼 있었음. GB3 수정으로 writeback 중 destination이 flag를 가진 두 경로(PE-array→**GLB**, GLB→**MC**)를 gate. MC→**DRAM** writeback은 destination(DRAM)에 double-buffer 개념이 없으므로 무조건 부과가 관례상 일관(비대칭 아님) |
| MC5 | `interconnect_timing.cc` | **해결됨 (2026-08-15)** | `cycle_temporal_chips`의 1-transaction 과다 — 공유 헬퍼 `pipelined_transfer_cycles` 수정으로 해결([local_buffer.md](local_buffer.md) LB5) |

> **참고:** `~600-1660`의 hand-rolled `nop_cycle*ceil(line_size/bitwidth)` 블록과 `num_access_multi_chip` counter reuse, SPARSE div-by-zero는 전부 `#ifdef FUNCTIONAL`/`#elif NPUSIM_LEGACY_ADDRESS_TIMING` 안이라 timing 빌드에서 dead code다. `utilization`은 `active_chips/num_chips`(공간 비율)만이며 traffic/contention 무관(설계상, topology P2에 인접).

## 정상 확인 항목 (dense load)

- INPUT/WEIGHT/OUTPUT **load**의 GLB↔chip↔DRAM traffic이 descriptor를 사용한다.
- NoP transfer cycle/energy가 모두 `num_access_multi_chip`을 사용한다(counter 일치).
- line/link width가 bits로 통일되고 각 active chip endpoint를 사용한다.
- `nop=bus`만 허용하고 mesh 등 미구현 NoP를 init에서 거부한다.
- overlap 식의 0 transaction guard, grid/shape/frequency/capacity/overflow 검증.

## 구현 결과

Dense **load** distribution과 DRAM load는 `datatype_transfer_timing()`으로 계산된다. Grid route/hop/contention, multicast fanout, OUTPUT write-back 비용은 미구현으로 남는다. Contention/route는 [global_cycle_overlap.md](global_cycle_overlap.md)의 timeline 위에서 정의된다.

## 완료 조건

- [x] INPUT/WEIGHT/OUTPUT **load**의 dense 방향이 descriptor byte·transaction model을 공유한다.
- [x] **OUTPUT write-back의 NoP 전송 비용을 복원한다(GB1과 연동).** (2026-08-15, MC1·MC3)
- [x] 패키지 temporal buffer static energy를 집계한다. (2026-08-15, SE3)
- [ ] sparse value/index/pointer stream을 모델링한다.
- [~] route/hop, router/link latency·energy를 구현한다(mesh routed-unicast, 2026-08-15). arbitration/back-pressure·multicast/fanout은 남음. 미구현 NoP는 명시적으로 거부(유지).
