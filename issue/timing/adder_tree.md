# Adder tree: timing simulation 이슈

> **재판정 (2026-09-02):** "남은 것" 중 PE-array 공통 이슈(PA1 등)는 전부 해결됐다
> ([pe_array.md](pe_array.md) — PA1 shared-stream multicast, PA2–4 leakage/window/
> dead-stat 배선 완료; 하단 공통 이슈 인용문도 stale). 실질 잔존은 **AT5 하나**:
> MAERI adder 단가(`adder_cycle/adder_energy`)의 실측 calibration — 모델이 단가에
> 정확히 반응함은 검증됐으므로 사용자 calibration 과제로 유지.

## 판정

**AT1~AT4, AT6 해결됨 (2026-08-15).** delta-splice를 제거하고 `account_descriptor_dense_writeback`을 virtual override하여 reduction 비용을 실제 발생 지점(OUTPUT write-back)에 부과한다. leaf 주입은 base의 line-granular 회계(data-volume 스케일), internal tree는 mapping에서 유도한 fan-in `F = N·e/E` 기준 element-granular adds(`e·(F−1)/F`, 합산 시 정확히 `E·(F−1)`) + `⌈log₂F⌉` depth fill. 남은 것: AT5/AT7 잔여(config calibration·bitwidth 절삭), PE-array 공통 이슈(PA1 등).

**검증 (eyeriss-기반 adder_tree swap, alexnet layer_0):** fan-in 유도 n=110·e=16·E=160 → F=11, adds=10, depth=4 (수작업 일치). `adder_cycle=2, adder_energy=3` 시 OUTPUT transfer energy 0→104,544,000 pJ = 3×(16·10/11)×2,395,800 events, cycle Δ=19,166,400 = 동일 event 수×8 — 두 독립 식이 동일 event 수로 교차 검증됨. in-place 누적에서도 매 write-back마다 부과(AT1 해소).

## 확인된 문제

| # | 위치 | 문제 | 영향 |
|---|---|---|---|
| AT1 (=T3) | writeback override | **해결됨 (2026-08-15)** — reduction이 OUTPUT write-back 이벤트에서 부과됨(`transfer_output`과 무관, in-place 누적 포함) | — |
| AT2 (=T3) | writeback override | **해결됨 (2026-08-15)** — leaf 주입은 base line-granular 회계, internal adds는 element 수 `e`에 비례(`e·(F−1)/F`) | — |
| AT3 (=T6) | writeback override | **해결됨 (2026-08-15)** — fan-in을 mapping에서 유도(`F = N·e/E`), disjoint output이면 F=1→비용 0 | 검증: F=11 (n=110,e=16,E=160) |
| AT4 (=T11) | writeback override | **해결됨 (2026-08-15)** — splice 제거로 OUTPUT의 access/link-transaction 통계가 base와 단일 경로로 일관됨 | — |
| AT5 (=T17) | `adder_tree.cc:157-166`, `231-233` | 기본 `u_adder_cycle=u_adder_energy=0`이고 어떤 config도 calibrate하지 않음 → reduction이 `depth·noc_cycle`, `(N-1)·noc_energy`로 붕괴 | 기본 config에서 OUTPUT latency가 이전 scatter 대비 silent 과소(예: N=168 → 8·noc_cycle) |
| AT6 | writeback override | **해결됨 (2026-08-15)** — 덮어쓰기 제거; base의 `pipelined_transfer_cycles` overlap 위에 depth fill을 **가산** | — |
| AT7 | `adder_tree.cc` | **해결됨 (2026-08-15)** | div-by-zero/UB guard + `derived_link_bitwidth()` 공용화로 fractional 절삭 시 경고(B12). 절삭 값 자체는 의도적 유지 |
| AT5 | config | **종결 (config 과제로 이관)** — `u_adder_cycle/energy`의 기본 0은 코드 결함이 아니라 **calibration 미실시**. 모델은 단가에 정확히 반응함이 검증됨(2/3 pJ 실험). MAERI 대상 실측/문헌 기반 단가는 사용자 calibration 항목 — `maeri.cfg`에 `adder_cycle/adder_energy` 키로 설정 |

> **PE-array 공통 이슈:** adder_tree는 `pe_array_t` base이므로 [pe_array.md](pe_array.md)의 PA2(temporal buffer leakage 미모델링)·PA3(reduction cycle이 leakage window 미포함)·PA4(`cycle_temporal_pe`가 dead stat → AT6 splice가 실질 no-op)가 그대로 적용된다. 특히 AT6의 `cycle_temporal_pe[OUTPUT]` splice는 PA4로 인해 어떤 통계에도 반영되지 않는다.

## 근거

균형 2-input tree는 활성 leaf `N`에 대해 `N-1` additions, `ceil(log2 N)` depth를 갖는다 — 헬퍼는 정확하다. 문제는 이 비용을 **어디에, 무엇에 비례해, 어떤 통계와 함께** 부과하느냐다. reduction은 partial-sum writeback 경로(`pe.cc`의 output write 이벤트)에서 tile data-volume에 비례해, reduction 차원의 실제 leaf 수로 계산되어야 한다.

## 완료 조건

- reduction 비용이 partial-sum accumulation/writeback 이벤트에서 부과되며, `transfer_output` 설정에 좌우되지 않는다.
- fan-in이 reduction 차원의 active PE 수와 일치한다(1/2/3/5/16/168 leaf에서 addition count·depth가 수작업과 일치).
- reduction cycle/energy가 tile data-volume(요소 수/line 수)에 따라 스케일한다.
- OUTPUT의 `access_*`·link-transaction·`num_data_transfer` 통계가 트리 모델과 일관된다.
- temporal buffer overlap 의미가 보존된다.
