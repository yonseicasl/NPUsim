# Systolic array: timing simulation 이슈

> **재판정 (2026-09-02): 아래 판정 헤더와 '⚠️ 최근 uncommitted 변경' 노트의
> "SY1/SY2 여전히 미해결" 문구는 수정 이전 기록이다.** 자체 표대로 SY1–SY5,
> V1–V3 전부 해결 상태이며(diameter fill + pipelined-hop, fold-fill/edge-accumulation
> RTL 캘리브레이션), Gemmini RTL 6점 MAPE 4.40%가 외부 검증으로 게이트에 상주한다
> (`unittest/run_timing_validation.sh`). event-level wavefront 정밀화(operand skew,
> per-hop register state)는 여전히 A1 per-tile timeline 위의 후속 범위 — timeline
> 자체는 구현됐으나 systolic wavefront를 tile보다 세밀하게 배치하는 작업은 남아 있다.

## 판정

Timing 빌드의 systolic 경로는 전적으로 `pe_array_t::account_descriptor_dense_distribution`이다(`systolic_array.cc:214-217`이 `#ifndef FUNCTIONAL`에서 여기로 early-return). 실제 systolic geometry(wavefront·skew·fill/drain·hop 수)가 cycle/energy에 반영되지 않고, **`noc_type` geometry(mesh hop)조차 적용되지 않는다.** 현 timing 결과는 flat store-and-forward scatter로 해석해야 한다.

- 심각도: **P1**

## 확인된 문제

| # | 위치 | 문제 | 영향 |
|---|---|---|---|
| SY1 (=T4) | `systolic_array.cc` | **해결됨 (2026-08-15; SP1로 정밀화)** | ~~mesh geometry 미적용~~ → 1차 수정(multiplicative 19×) 후 SP1의 pipelined-hop 계약으로 정밀화: 현재는 additive fill(SY2 행 참조). 초기 19× 검증치는 superseded |
| SY2 | `systolic_array.cc` + `pe_array.cc` | **해결됨 (2026-08-15, analytical)** | ~~wavefront fill/drain 부재~~ → systolic은 구조적으로 2D grid이므로 **noc 라벨과 무관하게** 매 stream에 diameter fill `((ax−1)+(ay−1))`·noc_cycle + per-hop(avg Manhattan) energy 적용. active shape가 latency를 바꿈. per-hop pipeline register state·operand skew·tile-boundary bubble의 event-level 정밀화는 A1 per-tile timeline 범위로 남음 | 검증: layer_0(11×10) INPUT 609,840 = T+18×events(mesh 기하와 동일), layer_4(다른 shape) 559,104 — shape 의존성 성립 |
| **V1 (Phase-1 validation)** | 비교 축 | **재판정: 해소됨 (2026-08-16, 코드 변경 0)** — 최초 판정("fold-level fill 모델 부재, FC −93%")은 **비교 축 오류**였음. fill/drain은 SY2가 이미 operand stream당 정확히 과금하며(`IC_w = txn + events×fill`, fc2 산술 검증), FC의 array-schedule 시간은 comp가 아닌 stream 축에 기록됨. 올바른 비교량 = **max(comp, IC_in, IC_w)** → FC 오차 **+6.3/+6.1/−5.6%**, 전체 MAPE 74.5→**21.9%**, ρ 0.024→**0.929** | [validation/phase1](../../validation/phase1/README.md). 남은 문서화된 차이: OUTPUT stream(write-through partial sum) vs WS in-array accumulation. (선택 개선: array-sched 파생값 직접 리포팅) |
| **V2 (Phase-2 RTL validation)** | compute fold overhead | **해결됨 (2026-08-17, RTL-캘리브레이션 구현)** — ~~fold-fill 항 미구현으로 소형/저-T GEMM −22~−73% 과소~~ → `pe_array_t` 배포 경로에 **weight-residency fold fill** 구현: WEIGHT 이벤트당 `u_weight_fold_fill_cycle × per-PE weight tile 원소수`(원소 각각이 실제 WS의 순차 residency; fold 수 = `(N/16)(K/16)` = RTL과 일치) + per-layer `layer_setup_cycle`(stats에서 repetition 스케일 후 1회). config knob(`weight_fold_fill_cycle=14`, `layer_setup_cycle=2270` — RTL 6점 그리드 적합의 유일 만족역), **기본 off로 기존 결과 byte-identical** | 검증: 6점 fill 값(224/896/3,584/14,336/1,792/14,336)이 `folds×14` 수기 산술과 전점 정확 일치; `comp+fold fill` vs RTL **전 6점 ≤8%, MAPE 4.4%**; alexnet matched/energy 회귀 byte-identical; unittest 5/5. 주의: per-원소 fold 해석은 GEMM/WS류 매핑(per-PE 시간축 weight=출력채널)에 유효 — conv류(per-PE R·S set)는 A1 범위 |
| **V3 (Phase-2 발견)** | 출력 누적 모델 | **해결됨 (2026-08-17)** — ~~edge-accumulator 부재로 Gemmini류 WS 매핑 표현 불가~~ → `edge_accumulation` knob: OUTPUT을 per-PE LB(`pe.cc check_tile_size` SEPARATE/SHARED)·array 집계(`pe_array.cc`) 용량 검사에서 면제(출력은 GLB output partition=accumulator로 스트림). `b_pe=M` 매핑으로 weight가 M 전체 상주 → 재스트리밍 3~4× 과대 해소. 기본 off | 검증: gemmini.cfg+`b_pe=M` 매핑으로 6점 전부 용량 검사 통과·실행(512³은 64KB accumulator 슬랩 제약이 K 루프를 DRAM으로 스필 — RTL의 C mvout과 동일 구조); weight 이벤트 수 4/8/16/64/32/128로 fold 수 RTL 일치 |
| SY3 | 공통 base | **해결됨 (2026-08-15)** — PA1의 shared-stream multicast 계약으로 해소([pe_array.md](pe_array.md)) | — |
| SY4 (=T4) | `systolic_array.cc` | **해결됨 (2026-08-15)** | ~~`noc_type` parse만 하고 미검증~~ → `is_supported_spatial_noc` init 거부 추가(spatial_arch 패턴). 검증: `noc=adder_tree` → "undefined NoC timing model" 오류 |
| SY5 (=T13/T18) | `systolic_array.cc` | **해결됨 (2026-08-15)** | div-by-zero/UB guard + 공용 `derived_link_bitwidth()`로 fractional 절삭 시 **경고**(B12; explicit `bitwidth`는 무음 우선). 절삭 값 자체는 기존 결과 보존 위해 의도적 유지 |

## ⚠️ 최근 uncommitted 변경에 대한 주의

`systolic_array.cc`의 최근 diff(data_type 라벨 `[INPUT]→[WEIGHT]`, `skip_transfer` 수정, `pes[i-width]→pes[i-num_active_pe_x]`, NoC 라벨 출력)는 모두 `#ifndef FUNCTIONAL … return;`(214-217) **뒤의 `#ifdef FUNCTIONAL` 블록** 안에 있다. 이는 **timing 빌드에서 실행되지 않는 dead code**이므로 timing cycle/energy를 바꾸지 않는다(diff 전후 byte-identical). FUNCTIONAL 정확성엔 유효한 수정이나, timing 관점의 SY1/SY2는 여전히 미해결이다. 오히려 `print_specification`이 이제 "Mesh"를 출력하므로 SY1(모델이 mesh를 무시)의 표기 불일치가 더 두드러진다.

## 수정 방향

Analytical model이라도 injection bandwidth, operand skew, active row/column, per-hop link latency/energy, fill/drain, partial-sum propagation, tile boundary bubble을 명시적으로 계산해야 한다(전역 timeline은 [global_cycle_overlap.md](global_cycle_overlap.md) 참조). 최소한 `spatial_arch.cc:220`처럼 `spatial_noc_cost()`를 timing 경로에 연결하고, `noc_type`을 init에서 검증한다.

## 완료 조건

- active dimension 변화가 이론적 fill/drain·hop latency만큼 total cycle을 바꾼다.
- link traversal 수 기반 energy가 수작업 micro-case와 일치한다.
- timing 경로가 `noc_type` geometry(mesh 등)를 적용하고, 미지원 topology를 init에서 거부한다.
- 출력 NoC 표기와 실제 timing path가 같은 topology를 나타낸다.
