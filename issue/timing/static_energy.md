# Static energy: timing simulation 이슈 (횡단)

## 판정

**SE1~SE4 모두 해결됨 (2026-08-15).** leakage는 이제 memory 포함 `layer_elapsed_cycles`(SE1)에 비례하고, always-on 도메인으로 PE·GLB가 전체 물리 chip에 통일됐으며(SE2), multi-chip temporal buffer leakage도 집계된다(SE3). PE·GLB·multi-chip static energy가 결과 파일에 출력된다(SE4). 정확한 wall-clock(component별 직렬/overlap)은 [global_cycle_overlap.md](global_cycle_overlap.md)의 timeline 도입 시 max-of-components 근사를 대체한다.

- 성격: 여러 컴포넌트(PE, GLB, multi-chip)에 걸친 횡단 이슈.
- 신규 관찰 **해결됨 (2026-08-15)**: `update_network_stats()`가 이제 PE/PE-array/GLB/multi-chip static energy를 layer 간 합산한다(network.txt total에 반영). 검증: network PE static = layer 합(이전 0). PE-array leakage도 포함([pe_array.md](pe_array.md) PA2).

## 확인된 문제

| # | 위치 | 상태 | 문제 | 영향 |
|---|---|---|---|---|
| SE1 | `scheduler/stats.cc:395-430` | **해결됨 (2026-08-15)** | ~~`pe_elapsed_cycles`(PE 내부 축 max)만 PE·GLB leakage에 사용~~ → `layer_elapsed_cycles = max(PE, GLB, multi-chip, DRAM)`로 교체 | memory-bound 레이어에서 leakage 시간창이 대폭 커짐(측정: eyeriss FC에서 pe_elapsed=1,926 → layer_elapsed=2.49M, 약 1300× 회복) |
| SE2 | `scheduler/stats.cc` | **해결됨 (2026-08-15, always-on)** | ~~PE leakage는 active chip만, GLB는 전체 물리 chip~~ → PE leakage도 `get_number_of_chips()`(전체 물리 chip)로 set·accumulate하도록 통일(GLB와 일치). PE static 누적을 active loop 밖의 전용 all-physical loop로 분리 | 비활성 chip도 leak(always-on). 검증: 2-chip(active=1) 구성에서 PE leakage all-physical=5.15e8 vs active-only baseline=2.58e8 = 정확히 2× (=num_chips/active) |
| SE3 | `components/multi_chip.{cc,h}`, `stats.{cc,h}` | **해결됨 (2026-08-15)** | ~~패키지 temporal buffer의 static energy 미집계(config `static_energy` 키가 파싱조차 안 됨)~~ → `multi_chip`에 `u_static_energy` 파싱 + `static_energy` + `update_static_energy()` 추가, stats에서 layer_elapsed로 set·accumulate·scale·print | `[multi_chip] static_energy`가 반영됨. 검증: 2.0pJ/cycle → 506,167,200 pJ(= 2.0×layer_elapsed×repetitions), 이전엔 config 무시로 항상 0 |
| SE4 | `scheduler/stats.cc` (print_results) | **해결됨 (2026-08-15)** | ~~`static_energy_pe`/`static_energy_global_buffer`가 결과 파일에 출력되지 않음~~ → PE·GLB 섹션에 "Static energy (leakage over layer elapsed cycles)" 블록 추가 | leakage가 결과에 출력됨. 검증: static_energy=1.0pJ/cycle에서 PE static=42.5e9 pJ(= layer_elapsed 1.53M 기준, PE-only 139k 대비 정확히 11× → SE1 window 반영 확인) |

집계 자체(한 번만 누적, `scale_serial_repetitions`가 static을 이중 스케일하지 않음)는 정상이다.

## SE1 해결 내용 (2026-08-15)

`global_buffer_t`, `multi_chip_t`, `dram_t`에 `modeled_elapsed_cycles()`(각자 cost 축의 max)를 추가하고, `stats.cc`에서 `layer_elapsed_cycles = max(pe_elapsed, 모든 GLB, multi-chip, DRAM)`를 구해 PE와 GLB leakage에 공통 사용한다. 전역 timeline이 아직 없으므로([global_cycle_overlap.md](global_cycle_overlap.md)) 이는 "component 최댓값" 하한 근사다 — component가 직렬화되는 경우(합) 여전히 과소일 수 있으나, memory 시간을 완전히 누락하던 PE-only 창보다 엄격히 정확하다. 정확한 wall-clock은 timeline 도입 시 대체한다.

## 남은 근거

SE1~SE4 모두 해결됐고, network-level 합산과 PE-array temporal buffer leakage([pe_array.md](pe_array.md) PA2)도 완료됐다. 남은 것은 (a) LB 전용 leakage 파라미터 분리([local_buffer.md](local_buffer.md) LB4), (b) DRAM background power 미모델([dram.md](dram.md) DR2 인접), (c) max-of-components 근사를 정확한 critical-path로 대체하는 [global_cycle_overlap.md](global_cycle_overlap.md) 뿐이다.

## 문서 정정 대상

- `implementation/timing/pe.md` — "static_energy ... **layer elapsed cycle**에 한 번만 곱한다"는 서술은 실제로 **PE-array elapsed**이므로 정정 필요.
- `implementation/timing/multi_chip.md` — "GLB static energy는 ... **modeled layer elapsed cycle** 기준"도 같은 정정 필요.

## 완료 조건

- [x] leakage 시간창이 memory 시간을 포함한다. **(2026-08-15 갱신: A1 timeline 도입으로 max-of-components 근사가 double-buffer-aware critical-path latency로 대체됨)** (SE1→A1)
- [x] PE와 GLB의 leakage 대상 chip 집합이 일치한다(always-on = 전체 물리 chip). (SE2)
- [x] multi-chip temporal buffer의 static energy가 elapsed에 비례해 집계된다. (SE3)
- [x] static energy가 결과 파일(layer)에 출력된다. (SE4)
- [x] network-level total에 static energy를 합산한다. (2026-08-15)
