# Global buffer: timing simulation 이슈

> **재판정 (2026-09-02):** GB7(신규, bypass P4)은 **해결됨** — P1-B로 bypass 스트림이
> descriptor 회계 경로에서 GLB SRAM 소스를 실제로 생략한다(`global_buffer.cc:190-225`:
> NoP ingress 직결 소스, packet `without_source`, fill/read 미과금 + L6로 boundary
> depth에서도 decoupling 제거). contention 잔여는 축소됨: CE4 shared-port 합산 규칙 +
> `num_banks` ideal interleave + per-tile timeline back-pressure가 반영됐고, 남는 것은
> per-address bank conflict와 datatype-level arbitration(선언된 모델 경계). sparse는
> A7 정책(functional 선행) 유지.

## 판정

GLB↔PE-array의 dense load와 **GLB→multi-chip OUTPUT write-back(`flush_data`)** 이 모두 datatype-aware descriptor bit transaction으로 계산된다. GB1(write-back link transfer 비용 주석 처리)·GB2(write-back access가 host-pointer 주소 walk)는 **해결됐다(2026-08-15)**. sparse encoding과 contention은 남는다.

- 완료 범위: GLB↔PE-array dense load, GLB→multi-chip/DRAM dense **load 및 OUTPUT write-back**
- 미완료: sparse encoding, contention(bank/arbitration/multicast)

## ⚠️ 완료 오표기 정정 (해결됨)

이전 문서는 "GLB↔multi-chip/DRAM dense 완료"라 했으나 write-back 경로에 아래 결함이 있었다. 모두 2026-08-15에 수정·검증됐다.

| # | 위치 | 상태 | 문제 | 영향 |
|---|---|---|---|---|
| GB1 (=T1) | `global_buffer.cc` `flush_data` (5개 블록) | **해결됨 (2026-08-15)** | ~~write-back link 비용이 주석 처리~~ → `account_output_writeback_link()` 헬퍼가 serialized NoP link transfer(`multi_chip->transfer_cycle/energy[OUTPUT]`, `overlapped_transfer_cycle`, payload/metadata/storage link transactions)와 write-back 오버랩 카운터를 부과 | off-chip output-store traffic이 집계됨. 검증: eyeriss AlexNet layer_0의 GLB↔multi-chip OUTPUT NoP transactions `0/0/0`→`36300/0/36300`, input/weight 무변화 |
| GB2 (=T7) | `global_buffer.cc` `flush_data` | **해결됨 (2026-08-15)** | ~~write-back access가 host-pointer `mask_bits` 주소 walk로 계산~~ → 5개 블록의 address walk를 제거하고 `account_output_writeback_link()`가 `datatype_transfer_timing()`의 `source_accesses`/`destination_accesses`로 GLB read·multi-chip write access를 계산(multi-chip→DRAM write-back 패턴과 동일) | packing-aware. 검증: 동일 160-element output tile에서 write-back access lines가 uint8=160 → int4=80 (2 elements/byte)로 스케일, 이전 walk는 host `sizeof(data_t)` 기준이라 format 무관하게 160 고정이었을 것 |

> **참고:** GB1 수정 시점엔 `write_back_cycle`이 dead stat이었으나, **SE1 이후 `global_buffer_t::modeled_elapsed_cycles()`가 이를 읽어** leakage window에 반영한다(효과는 benign — `access_cycle[OUTPUT]`의 부분집합이라 max를 지배하지 않음). 관측 가능한 GB1 효과는 `multi_chip->transfer_cycle/energy[OUTPUT]`와 link-transaction 카운터다.

## 추가 확인된 문제 (2026-08-15 재감사, 미문서화)

| # | 위치 | 심각도 | 문제 | 영향 |
|---|---|---|---|---|
| GB3 | `global_buffer.cc`, `pe_array.cc` | **해결됨 (2026-08-15)** | (정정: "미소비" 주장은 부정확 — MC→GLB load는 `multi_chip.cc:405`에서 이미 GLB flag로 gate됨. 실제 gap은 **writeback 경로**) PE-array→GLB writeback과 GLB→MC writeback(`account_output_writeback_link`)의 access cycle을 destination flag로 gate(energy 항상 부과 — load 관례와 통일). `check_tile_size`가 double buffer 시 2-copy 용량 요구 | 검증: `double_buffer=1` → GLB OUTPUT access cycle 4,646,400→26,400·energy 불변 14,010,480; 용량 1.5KB에서 1×(1012B) 통과, 2×(2024B) 거부 |
| GB4 | `global_buffer.{cc,h}` | **해결됨 (2026-08-15)** | ~~separate buffer utilization 분모가 합산 size~~ → base에 `capacity_per_type` 도입(separate=partition, shared=total), 0-capacity guard 포함 | shared(eyeriss) 0.3% 무회귀 확인 |
| GB5 | `global_buffer.cc` | **해결됨 (2026-08-15)** | ~~`bitwidth = 8*bandwidth/frequency`의 frequency 0 미검증(inf→unsigned UB)~~ → `frequency>0` ternary로 divide 보호 + 최종 `bitwidth==0` 거부. 검증: frequency=0 config가 UB 대신 clean error. (dram/spatial/adder/systolic/multi_chip 동일 수정) |
| GB6 | `global_buffer.cc` (FUNCTIONAL load 경로) | **P3 (미해결)** | FUNCTIONAL 빌드의 load 경로가 여전히 host-pointer `mask_bits` walk. timing 빌드에선 dead code | FUNCTIONAL 빌드에서 INT4 load access가 descriptor 경로와 발산 |
| GB7 | `interconnect_timing.cc` | **해결됨 (2026-08-15)** | `cycle_pe_array_global_buffer`의 1-transaction 과다 — 공유 헬퍼 수정으로 해결([local_buffer.md](local_buffer.md) LB5) |

## 정상 확인 항목 (dense load)

| 항목 | 구현 계약 | 검증 |
|---|---|---|
| Line/link 단위 | `line_size`, `bitwidth`를 bits로 단일 해석 | 16/32/64-bit MXFP micro-test |
| Dense load | INPUT/WEIGHT/OUTPUT load가 `account_descriptor_dense_transfer` 사용 | GLB source, PE-array destination, link transaction 일치 |
| MXFP | payload/scale metadata를 separate/interleaved에 맞춰 계산 | 33-element MXFP4 separate boundary |
| Sparse 거부 (load) | `global_buffer.cc:181-185`이 non-DENSE `compression_type`을 top에서 거부 | timing entry guard (단, DRAM OUTPUT 경로엔 hole — [dram.md](dram.md) DR1) |

## 남은 문제

| 우선순위 | 문제 | 영향 |
|---|---|---|
| ~~P1~~ 완료 | GB1 write-back transfer 비용 복원 및 descriptor화 | off-chip output-store timing/energy 복원됨 |
| ~~P1~~ 완료 | GB2 write-back access의 descriptor 전환 | packing-aware access 카운트 |
| P1 | CSC/COO payload·index/pointer stream 및 sparse decoder 비용 | sparse GLB timing 미지원 |
| ~~P2~~ 부분 해결 (=T9) | bank/contention — **`num_banks` knob 추가(2026-08-15)**: ideal line-interleaving으로 access cycle 1/B(energy 불변), conflict-free 상한 명시. 검증: num_banks=4 → Access cycle 10,062,360→2,515,590(정확히 /4)·energy 불변 | **남음**: bank conflict, read/write port, arbitration, back-pressure — A1 per-tile timeline 위에서 정의 |
| ~~P1~~ 완료 (=SE1) | GLB static energy 시간창 — layer critical-path 창으로 해결(2026-08-15) | [static_energy.md](static_energy.md) |
| GB5 (traffic 검증 적발, 신규→해결) | `global_buffer.cc` data_transfer WEIGHT | **해결됨 (2026-08-17)** — ~~legacy host-address timing 경로가 `#else`(timing 빌드)로 활성이어서 descriptor 계정과 **이중 계상**(num_data_transfer 2×, GLB weight access cycle/energy ~2×)~~ → INPUT과 동일하게 `#if 0` 처리 | 검증: traffic 하네스 T5(GLB weight 바이트=이벤트×array tile)가 0.500→**1.000** 전 8케이스 |
| GB6 (energy 검증 적발, 신규→해결) | `multi_chip.cc`/`global_buffer.{h,cc}`/`stats.{h,cc}` | **해결됨 (2026-08-17)** — ~~mc→GLB fill(GLB write)이 uniform ×R 스케일되어 DR4의 datatype별 off-chip credit과 자기모순(공급 없는 데이터가 52× 쓰임)~~ → GLB `fill_access_{cycle,energy}` 분리 계정, stats에서 per-type 스케일 후 access에 합산 | 검증: conv3 GLB weight 184MB→94MB(read 92 + fill 1.77 ✓); energy E3 비율 1.78~2.57→2.20~3.04 ([validation/energy](../../validation/energy/README.md)) |
| GB7 (신규) | GLB `bypass` | **P4 — bypass 플래그가 descriptor 계정 경로에서 미적용**(용량 검사에서만 소비). eyeriss.cfg `bypass=0:1:0`인데 weight GLB 접근이 과금됨 — 칩은 실제로 filter를 GLB에 저장하므로 결과값은 오히려 현실적이나, 플래그 의미론 불일치 | 수정 방향: descriptor 경로에서 bypass 시 GLB access 생략(직결 스트림) 또는 cfg에서 bypass 제거 |
| GB8 (energy 검증 적발, 신규→해결) | `stats.cc` GLB 섹션 출력 | **해결됨 (2026-08-17)** — GLB Interconnection(GIN) energy가 수집만 되고 미출력 → "Interconnection energy" 3행 추가 | 검증: eyeriss_energy run에서 GIN 2.3~6.2% 가시화 |
| GB9 (외부 assessment 적발, 신규→해결) | `global_buffer.cc` `modeled_elapsed_cycles` | **해결됨 (2026-08-17)** — ~~GB6의 fill 분리 후 fill cycle이 busy/leakage timeline에서 누락되어 `busy_cycle_glb ≥ max(access_cycle_glb)` 불변식 위반(전 8층)~~ → busy 축에 `access + fill` 합 반영(read/write 포트 직렬화의 보수적 상한; 스케일 전 합산이라 per-type ≤ uniform 관계로 스케일 후 불변식 보장) | 검증: 전 8층 불변식 성립; latency 축(comp+fold)·게이트(check_timing base/baseline) 불변 통과 |

## 완료 조건

- [x] single-chip dense GLB↔PE-array load가 descriptor bit transaction을 사용한다.
- [x] line/link width가 bits로 고정되고 validation·micro-test가 있다.
- [x] **OUTPUT write-back의 transfer 비용을 복원(descriptor 기반 link transfer)한다.** (2026-08-15)
- [x] write-back access 카운트를 descriptor bit transaction으로 전환한다(GB2). (2026-08-15)
- [ ] sparse value/index/pointer traffic을 구현한다.
- [ ] contention model 또는 명시적 serialized-requester 제약을 도입한다.
