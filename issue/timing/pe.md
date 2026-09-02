# PE: timing simulation 이슈

> **재판정 (2026-09-02): 판정 헤더가 stale — mac_width 계약은 구현·검증 완료.**
> PE1(lane→accumulator reduction fill)에 더해 L9/PE1-PE2로 **active lane state 기반
> 과금**이 구현됐다(`calculate_mac_lane_state()` → `pe.cc:575` 소비; latency=가장 깊은
> live tree, energy=live tree 전체 작업량 — `utils/pe_lane.h` L9 계약·unittest 고정).
> PE2의 "informational" 잔여가 이로써 해소됐고, `mac_energy_<input>_<weight>`
> precision family(E3)와 reduction energy 축(E4, `mac_reduction_energy`)도 추가됐다.
> PE5는 B12 경고로 설계 확정. 남는 것: PE3의 per-tile utilization 세분화(현재
> max-peak hardening, tile 크기가 레이어 내 고정이라 benign)뿐. 아래 원문 판정의
> "`mac_width > 1` 신뢰 불가"는 더 이상 사실이 아니다.

## 판정

`num_macs`와 `mac_width`의 마이크로아키텍처 의미가 cycle/energy에 반영되지 않는다(계산만 하고 안 읽는 dead field 존재). `mac_width > 1` 구성의 cycle·energy는 신뢰할 수 없다. static energy 시간축(PE4)과 utilization 분모(PE8)는 해결됐다(2026-08-15).

- 심각도: ~~P1~~ → **재판정으로 해소 (2026-09-02, 상단 참조)**

## 확인된 문제

| # | 위치 | 문제 | 영향 |
|---|---|---|---|
| PE1 (=T5) | `pe.cc`, `utils/pe_lane.{h,cc}` | **해결됨 (2026-08-15)** | ~~mac_width가 latency에 무영향~~ → 문서화된 계약(num_macs=accumulator 수, mac_width=lane 수/accumulator)대로 **lane→accumulator adder-tree fill `⌈log₂(mac_width)⌉`단**을 issue당 latency에 가산(`lane_reduction_fill_cycles`, pipelined이라 throughput·per-FMA energy 불변). 13개 computation site 공통 적용 | 검증: 64×1 cycle 3,833,280 vs 8×8 **15,333,120=정확히 4×**(1+3 fill), energy 양쪽 동일 236,130,048 pJ. unittest가 fill(1→0, 2→1단, 5/8→3단) 직접 assert |
| PE2 (=T5) | `pe.cc:496-498` | **부분 해결** — `mac_width`(하드웨어 구조)가 이제 cycle 경로를 구동. `active_mac_width`/`active_mac_units`는 여전히 informational(트리 깊이는 구조적이라 active lane 수와 무관 — 의도된 설계) | — |
| PE3 | `pe.cc` | **해결됨 (2026-08-15, hardening)** | (정정: tile 크기가 레이어 내 고정이라 overwrite는 현재 benign) `std::max` peak-누적으로 hardening. 필드가 소비되는지 여부는 PE2/A3에서 결정 | — |
| PE4 (=SE1) | static energy 시간축 | **해결됨 (2026-08-15)** — leakage가 이제 memory 포함 `layer_elapsed_cycles`로 스케일됨 | 횡단 이슈 → [static_energy.md](static_energy.md) |
| PE5 (=T13/T18) | `pe.cc`, `interconnect_timing` | **부분 해결 (2026-08-15)** — 공용 `derived_link_bitwidth()`가 fractional 절삭 시 **경고** 출력(명시 `bitwidth`는 무음 우선). 절삭 자체는 기존 결과 보존 위해 유지, non-power-of-two 허용도 유지(경고로 가시화) | unittest: 2.0/1.0→16, 0.9/1.0→7(+warn), freq 0→0 |
| PE6 (=T19) | `stats.cc` | **해결됨 (2026-08-15)** | ~~`avg_computation_cycle /= num_active_pe` 0 미방어~~ → `if(num_active_pe > 0)` guard (avg_computation/mac/lb 전부) |
| PE7 | `pe.cc` (`undefined_stationary_t::computation`) | **해결됨 (2026-08-15, 방어적)** — ~~루프 없이 1회만 증가하여 energy/utilization 64× 과소~~ → `num_active_macs` 루프 추가(다른 3개 stationary와 동일). **단, 이 경로는 현재 도달 불가**: `pe_t::init`(`pe.cc:293`)이 `mac_stationary=undefined_stationary`를 `exit(1)`로 거부하므로 실행 전 abort됨(검증: "Error: PE requires a valid mac_stationary"). 수정은 향후 거부가 풀릴 경우를 위한 정합성 보강 |
| PE8 | `stats.cc` | **해결됨 (2026-08-15)** | ~~`mac_available_cycle`이 `pe_elapsed`(memory 제외) 사용~~ → `layer_elapsed_cycles`(memory 포함, SE1 window)로 교체. 검증: eyeriss(memory-bound) MAC utilization 1.0%로 정정(이전 near-full) |

## 해소 확인 (이전 issue 목록에서)

- `double` latency의 `unsigned` 절삭 — 해소됨. cycle 통계는 `double`, `pipeline_transactions`는 cast 전 guard.
- 0/1/2 transaction overlap 경계(`pipelined_transfer_cycles`, 수기 등가물) — off-by-one 일관.
- `bitwidth`/`line_size` divide-by-zero(init 경로) — guard 있음.

## 근거와 수정 방향

`mac_register_capacity = num_macs*mac_width`로 storage는 할당하지만 `active_mac_width`는 timing에 쓰이지 않는다. Simba 8×8을 64 scalar FMA로 볼지, 8 accumulator에 8 multiplier lane으로 볼지 먼저 정하고, reduction latency·register port·active component·energy를 그 정의에 맞춰 **실제 cycle/energy 경로에 연결**해야 한다(현재는 계산 후 폐기). static energy는 [static_energy.md](static_energy.md)의 layer critical-path 시간창을 사용한다.

## 완료 조건

- [ ] `num_macs`/`mac_width`의 lane·accumulator mapping이 cycle·utilization·energy를 **실제로 구동**한다(PE1/PE2, dead field 제거 또는 연결).
- [ ] `mac_width>1`이 reduction latency/energy를 변화시킨다(PE1).
- [ ] utilization이 tile 누적 기준으로 레이어 전체를 반영한다(PE3).
- [x] static energy가 elapsed simulation time에 비례한다(PE4=SE1). (2026-08-15)
