# PE array: timing simulation 이슈

## 판정

`pe_array_t`(spatial_arch / systolic_array / adder_tree의 공통 base)의 timing 경로는 `account_descriptor_dense_distribution()`과 `account_descriptor_dense_writeback()`이다. 이 경로는 **broadcast/multicast reuse를 표현하지 못하고**(operand를 목적 PE 수만큼 복제 집계), **PE-array temporal buffer의 leakage를 모델링하지 않으며**, PE-array 자체의 overlap/write-back cycle 통계가 **소비되지 않는 dead stat**이다. 이 컴포넌트는 그동안 전용 문서가 없었다.

- 심각도: **P1**
- 관련: [spatial_architecture.md](spatial_architecture.md), [systolic_array.md](systolic_array.md), [adder_tree.md](adder_tree.md), [static_energy.md](static_energy.md)

## 확인된 문제

| # | 위치 | 심각도 | 문제 | 영향 |
|---|---|---|---|---|
| PA1 | `pe_array.cc` `account_descriptor_dense_distribution` | **해결됨 (2026-08-15)** | ~~source read/link를 목적 PE마다 복제~~ → **shared-stream multicast 계약**: distinct data(배열-레벨 tile)를 temporal buffer에서 1회 읽어 fabric에 1회 스트리밍, per-PE는 LB write만. overlap은 max(src,link,dest) 단일 pipeline 항 | 검증(eyeriss layer_0): INPUT NoC txn 2,395,800→217,800(**정확히 11×**, 공유 factor와 일치); WEIGHT 14,374,800→**2,635,380 = GLB↔MC weight 전송량과 정확히 일치**(각 weight가 배열에 1회 진입); OUTPUT −0.7%(disjoint, tail-rounding만). SY3도 동일 base라 함께 해소 |
| PA2 | `pe_array.{cc,h}`, `stats.{cc,h}` | **해결됨 (2026-08-15)** | ~~PE-array temporal buffer에 static/leakage energy 부재~~ → `u_static_energy`(`pe_array_static_energy` 파싱, `initialize_temporal_buffer`) + `static_energy` + `update_static_energy()` 추가, stats에서 set·accumulate·scale·print. 검증: `pe_array_static_energy=3.0` → 1,770,975,360 pJ(이전 항상 0) |
| PA3 | `stats.cc`, `pe_array.cc` | **해결됨 (2026-08-15)** | ~~PE-array cycle이 leakage window 미포함~~ → `pe_array_t::modeled_elapsed_cycles()` 추가, `layer_elapsed_cycles` max에 포함. 검증: eyeriss layer_0 window 1.53M→3.58M(pe_array가 실제 병목) |
| PA4 | `pe_array.h` + `stats.{cc,h}` | **해결됨 (2026-08-15)** | PA3의 `modeled_elapsed_cycles()`가 leakage window에서 소비 + `cycle_temporal_pe`를 stats로 배선해 "Overlapped cycle (temporal buffer)"로 출력(수집/스케일/네트워크 합산 포함). 검증: 기본 0.0, `exist_temporal_buffer=1` → 6.5M/84.4M/43.2M 출력 |
| PA5 | `pe_array.cc` `account_descriptor_dense_writeback` | **해결됨 (2026-08-15)** | ~~output writeback이 raw noc_cycle만 사용~~ → `spatial_noc_cost()` topology multiplier 적용(non-mesh는 multiplier 1이라 무영향). 검증: noc=mesh에서 OUTPUT writeback transfer 2,740,650→45,865,050 (≈16.7× hop scaling), noc=bus 불변 |
| PA6 | `pe_array.{cc,h}`, `stats.{cc,h}` | **해결됨 (2026-08-15)** | ~~temporal buffer 점유 utilization 부재~~ → `buffer_utilization[type]`(peak across layer)을 distribution 경로에서 계산, stats 집계·"Temporal-buffer occupancy"로 출력 | 검증: eyeriss layer_0 → 7.4/5.1/2.0% 신규 출력 |
| PA7 (=V2) | `pe_array.{cc,h}`, `stats.{cc,h}` | **해결됨 (2026-08-17)** | RTL-캘리브레이션 weight-residency fold fill: WEIGHT 배포 이벤트당 `u_weight_fold_fill_cycle × per-PE weight tile 원소수` 과금 + per-layer `layer_setup_cycle`(repetition 스케일 후 1회 가산), "Fold fill cycle"로 출력. 기본 off(기존 결과 byte-identical) | 검증: Gemmini RTL 6점 `comp+fold fill` ≤8%(MAPE 4.4%); fill 값 수기 산술 전점 일치. 상세: [systolic_array.md](systolic_array.md) V2 |
| PA8 (=V3) | `pe_array.{cc,h}`, `pe.{cc,h}` | **해결됨 (2026-08-17)** | `edge_accumulation` knob: OUTPUT을 per-PE LB·array 집계 용량 검사에서 면제(array-edge accumulator 모델, 출력은 GLB output partition으로 스트림) → `b_pe=M` weight-residency 매핑 표현 가능. 기본 off | 검증: [systolic_array.md](systolic_array.md) V3, [validation/phase2](../../validation/phase2/README.md) |
| PA9 (최종 판정 2026-08-17) | GLB↔array 트래픽 모델 | **판정 완료 (경계 검증으로 종결; point-match는 외부 데이터 부재로 불가)** — GLB 축의 point 오차(2.0~7.7×)는 칩의 per-layer mapper 파라미터(m,n,e,p,q,r,t — ISCA/JSSC/공개 thesis 모두 미게재)가 있어야만 해소 가능함을 확인. 임의 knob 적합 대신 **경계 검증**으로 종결: (i) 칩 GLB 실측이 전 5층에서 [완전 on-chip 보존 하한, literal 재스트리밍 상한] 내 (compare.py band verdict ALL ok); (ii) DRAM 축은 칩의 공개 zeros%로 RLC 압축을 되돌린 dense-등가 기준 **0.95~1.50×** (conv3 −5%); (iii) latency 축 무영향. 원인 분해(batch-interleave 필터 공유 + psum 왕복)는 유효하나 크기 분리는 동일 데이터 부재로 불가 — 15% traffic gate는 informational 유지 | [validation/phase3](../../validation/phase3/README.md). 재개 조건: 칩 mapper 파라미터 확보(예: MIT dspace thesis 접근) 시 pass-구조 재구성으로 point-match 시도 |

## 근거

두 concrete driver(spatial_arch, systolic_array)와 adder_tree는 `#ifndef FUNCTIONAL`에서 곧바로 공통 base 메서드로 early-return하므로, timing 숫자는 전부 `pe_array_t`의 이 경로에서 나온다. 따라서 위 결함은 세 variant 모두에 공통 적용된다.

> **재판정 (2026-09-02):** 아래 첫 체크박스는 stale — PA1은 자체 표대로 해결됨
> (shared-stream multicast, 2026-08-15 검증 포함). 이후 추가: **stripe-transition
> 보정**(2026-08-26, nv_small RTL 실측 — `stripe_transition_cycle`, NVDLA holdout
> 0.55%의 구성 요소)과 SFU output-path 통합(plan_sfu.md).

## 완료 조건

- [x] shared operand의 broadcast/multicast reuse가 목적 PE 수에 선형 복제되지 않는다(PA1). (2026-08-15)
- [x] PE-array temporal buffer가 `u_static_energy` 파싱 + `update_static_energy()` + stats 집계·출력을 갖는다(PA2). (2026-08-15)
- [x] `pe_array_t::modeled_elapsed_cycles()`가 `layer_elapsed_cycles`에 반영된다(PA3). (2026-08-15)
- [x] overlap/write-back 카운터가 소비·보고된다(PA4). (2026-08-15)
- [x] output writeback도 `noc_type` geometry를 적용한다(PA5). (2026-08-15)
