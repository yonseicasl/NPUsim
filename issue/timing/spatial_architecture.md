# Spatial architecture: timing simulation 이슈

> **재판정 (2026-09-02):** 판정 헤더의 "잔여 이슈" 표현은 stale — SP1–SP3 전부 자체
> 표에서 해결 상태다. 이후 외부 검증도 추가됨: **NVDLA nv_small(spatial 경로)이
> holdout MAPE 0.55%·DRAM traffic EXACT**로 게이트에 상주한다
> (`validation/nvdla/compare_full.py`). 잔존은 모델 경계(router queue/VC/adaptive
> route)로, per-tile timeline 위의 별도 topology model 범위다.

## 판정

**명시된 timing abstraction 범위 내에서는 동작한다.** BUS/store-and-forward/crossbar는 serialized-unicast, MESH는 ingress `(0,0)`에서 active PE까지 Manhattan route를 따르는 routed-unicast로 모델링하고, 미지원 topology는 init에서 거부한다(`is_supported_spatial_noc`). systolic 경로와 달리 `spatial_noc_cost()`를 timing에 실제로 연결한다(`spatial_arch.cc:220`). 다만 hop 비용 적용 방식과 init 검증에 잔여 이슈가 있다.

## 해결된 항목

| 항목 | 구현 계약 | 검증 |
|---|---|---|
| Topology | BUS, store-and-forward, crossbar, MESH 지원 | topology support micro-test |
| Route/hop | MESH latency=최대 Manhattan hop, energy=평균 Manhattan hop | 2×3 MESH: latency 3, energy 1.5 |
| Multicast | ~~destination별 replication~~ → **shared-stream multicast**(2026-08-15, PA1): distinct 배열 tile을 1회 read·1회 stream, per-PE는 LB write만 | [pe_array.md](pe_array.md) PA1 검증 참조 |
| Tile state | tile update마다 read-granular flag reset | 이전 tile read state 재사용 방지 |

## Phase-3 실리콘 검증 (2026-08-17)

Eyeriss 실칩 실측(JSSC'17 Table V, AlexNet b4, 200MHz) 대비 spatial_arch(RS)
경로의 Processing-Latency가 **전 5개 conv layer ≤6.4% (MAPE 4.3%)**로 검증됨
([validation/phase3](../../validation/phase3/README.md)).

- 매핑: 칩식 spatial padding(conv1 P 55→7×8, conv2 P 27→28 + C의 PE_Y 이동)으로
  active PE를 칩과 일치(154/140/156/156/156). 매핑 곱>차원의 padded 표현이 기존
  tile 기계에서 일관 동작(`silicon.map`); `npu.cc`에 padding/부분-커버리지
  **경고 가드** 추가(`mapping_table::calculate_total_parameter_size`).
- 캘리브레이션: pass ramp-up(논문 명시 미달 원인)을 V2 fold-fill knob으로 과금 —
  `eyeriss.cfg weight_fold_fill_cycle = 0.11`(5점 적합 허용역 [0.089, 0.127]).
- 메모리축(2026-08-17, DR4+DR5 해결 후): **DRAM 0.8/1.3/1.0/1.6/2.2×** —
  전 5층 weight fetch-once 정확([dram.md](dram.md) DR4·DR5 해결, 잔여 DR6).
  **GLB↔array 2.0~7.7×**는 batch-인터리브 필터 공유 + psum 왕복의 상쇄 쌍으로
  분해 확정([pe_array.md](pe_array.md) PA9, P4). latency 축은 전 수정에서 불변.

## 확인된 문제 (잔여)

| # | 위치 | 문제 | 영향 |
|---|---|---|---|
| SP1 (=T10) | `interconnect_timing.{h,cc}` 등 | **해결됨 (2026-08-15)** | ~~MESH latency = T×max-hop~~ → **pipelined-hop 계약**: stream latency = `(T + max_hops−1)`·noc_cycle (route depth는 1회성 fill), energy는 per-txn×avg-hop 유지(이제 metric 불일치가 아니라 의도된 분리: latency=critical path, energy=총 traversal). distribution·writeback 공통 적용 | 검증: eyeriss mesh INPUT 609,840 = 217,800(T) + **18(fill)×21,780(events)** 수작업 일치(이전 19×T=4.1M 과다 제거). 2×3 mesh unittest: fill=2, energy=1.5 |
| SP2 (=T12) | `spatial_arch.cc` | **해결됨 (2026-08-15)** | ~~`frequency>0`/`bitwidth!=0` init 미검증(inf→unsigned UB)~~ → `frequency>0` ternary + `bitwidth==0` clean error(T1-4 공통 batch, dram/glb/adder/systolic/multi_chip 동일) |
| SP3 (=T16) | `interconnect_timing.cc` | **해결됨 (2026-08-15)** | ~~`pipelined_transfer_cycles`의 1-transaction 과다~~ → 함수를 raw (source, link, dest) 인자로 리팩터해 1-transaction=source+link+dest로 수정. N≥2는 기존과 동일. 6개 call site 모두 갱신, unittest가 `==8`(이전 12) 직접 assert |

> **횡단(shared helper) 이슈 — 해결됨:** SP3는 `interconnect_timing.cc` 공유 함수 버그였고, **한 곳 수정으로 모든 호출 경로가 동시 해결**됐다 — GLB `cycle_pe_array_global_buffer`([global_buffer.md](global_buffer.md) GB7), PE local buffer `cycle_mac_lb`([local_buffer.md](local_buffer.md) LB5), multi-chip `cycle_temporal_chips`([multi_chip.md](multi_chip.md) MC5), PE-array `cycle_temporal_pe`, DRAM `cycle_chip_dram`.

## 모델 경계

router queue, virtual channel, adaptive route, cycle-accurate back-pressure는 제공하지 않는다. 이를 요구하는 비교에는 [global_cycle_overlap.md](global_cycle_overlap.md)의 timeline 위에서 별도 topology model이 필요하다.

## 완료 조건

- MESH latency의 hop 적용이 물리적으로 정당화된다(전체 traffic×max-hop이 아니라 파이프라인 hop 지연으로), latency/energy의 hop metric이 일관된다.
- `frequency`/`bitwidth` 유효성이 init에서 검증된다.
- `pipelined_transfer_cycles`의 1-transaction 케이스가 `read+link+write`를 반환한다.
