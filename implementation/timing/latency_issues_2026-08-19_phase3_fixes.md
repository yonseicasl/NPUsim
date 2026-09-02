# 2026-08-19 latency 재평가 Phase 3 구현 기록

기준 문서: [../../assessment/TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md](../../assessment/TIMING_SIMULATION_LATENCY_ISSUES_2026-08-19.md)
선행 기록: [latency_issues_2026-08-19_phase1_fixes.md](latency_issues_2026-08-19_phase1_fixes.md),
[latency_issues_2026-08-19_phase2_fixes.md](latency_issues_2026-08-19_phase2_fixes.md)

Phase 3(L8/L9 + 항목 10 "component knob별 non-zero micro-test")를 구현했다.

## L8. DRAM analytical 계약 — 선택지 1을 택함

평가 문서는 둘 중 하나를 **명확히 선택**하라고 요구했다.

1. analytical 모델의 지원 범위를 sequential/conflict-free stream으로 제한하고 이를 출력한다.
2. address trace와 DRAMSim3를 공식 memory-bound latency gate에 연결한다.

**선택지 1을 택했다.** 근거: timing 경로에는 per-request address가 없다
(`account_descriptor_dense_load(type, elements)`는 datatype과 element 수만 받는다). address를
만들어내서 bank conflict를 계산하면 모델이 갖고 있지 않은 정보를 근거로 정밀도를 주장하는 것이 된다.
선택지 2는 address trace 인프라가 선행되어야 하며 DRAMSim3 compile path는 이미 존재한다.

### 1) 계약을 출력에 명시

```text
DRAM timing model     : analytical open-page, JEDEC tRC=tRAS+tRP, ideal bank interleaving, bus turnaround
  not modeled here    : request queueing, address-based bank conflicts, burst scheduling, strided/random
                        access (row misses are spread evenly across banks, not mapped by address)
```

지원 범위를 코드에서 추론하게 두지 않고 결과 파일에 적는다. 특히 **idealized bank interleaving을
bank conflict 모델로 읽지 않도록** 명시한다. 문자열은 knob에 따라 변하므로 계약이 knob-sensitive함을
검증할 수 있다.

구현: [components/dram.cc](../../components/dram.cc)
`describe_timing_model()` / `describe_timing_limits()`,
[scheduler/stats.cc](../../scheduler/stats.cc) 출력 및 network rollup.

### 2) address가 필요 없는 물리적 요소 추가 — bus turnaround

data bus 방향이 read↔write로 뒤집히면 tWTR/tRTW가 든다. 이것은 **address가 전혀 필요 없는** 유일한
detailed DRAM timing 요소다 — 모델은 자기가 load를 서비스하는지 output write-back을 서비스하는지
이미 알고 있다. 따라서 근사 없이 정확히 과금한다.

- config: `t_wtr_cycle`, `t_rtw_cycle` (기본 0 = 비활성, 기존 config 불변)
- load는 read 방향, output write-back은 write 방향으로 기록하고 전환 시에만 과금

### 3) row activation 산술을 순수 helper로 분리

`dram_row_activations(stream_bytes, row_buffer_bytes, num_banks)`가
`{activations, busiest_bank}`를 돌려준다. **energy는 activations**(모든 activation은 겹쳐도 에너지를
쓴다), **latency는 busiest_bank**(bank 병렬성은 지연만 숨긴다)에 비례한다. 수기 unit test로 고정했다.

### 4) knob을 켠 fixture와 수기 항등식

기존에는 `row_buffer_size` / `t_ras_cycle` / `t_rp_cycle` / `num_banks`를 켠 config이 하나도 없어서
**knob이 임의로 틀려도 아무 회귀도 잡지 못했다.** knob만 다른 3개 config을 추가했다.

| config | knob |
|---|---|
| `gemmini` | off (flat per-access) |
| `gemmini_dram_detail` | 256 B row, tRC=10, 8 banks, tWTR=tRTW=2 |
| `gemmini_dram_serial` | 같은 설정에 num_banks = 1 |

256 B row는 하나의 DRAM tile이 여러 row에 걸치게 하려는 **의도적으로 작은 값**이다. 이 fixture는
knob 산술을 검증하는 것이고 특정 DDR3 part를 모델링하는 것이 아니다(문서화함).

[validation/dram/check.py](../../validation/dram/check.py) DR1~DR6:

- DR1: transaction count 3개 run 모두 동일 (row activation은 traffic이 아니라 timing)
- DR2: compute-schedule latency 3개 run 모두 동일 (DRAM timing은 검증 metric을 건드리지 않는다)
- DR3: off < 8 banks < 1 bank (activation은 비용이고 bank 병렬성이 일부를 숨긴다)
- **DR4 수기 항등식**: bank 병렬성 절감 = `sum_type(streams × (rows − ceil(rows/8)) × tRC)`
  → 측정 380 == 수기 380
- **DR5 수기 항등식**: knob 활성화 초과분 = `sum_type(streams × ceil(rows/8) × tRC) + turnaround 1회`
  → 측정 102 == 수기 102
- DR6: 출력된 계약이 knob을 따라간다 (flat run은 row model을 주장하지 않고, 1-bank run은 bank
  interleaving을 주장하지 않는다)

`streams`와 `rows`는 보고된 transaction 수와 processor tile 크기에서 역산하므로 항등식은 리포트만으로
닫힌다.

## L9-a. Active-lane 기반 reduction

### 문제

모든 issue가 **구조적** `mac_width`로 `lane_reduction_fill_cycles` / `lane_reduction_energy`를
과금했다. 이는 양방향으로 틀렸다.

- 여러 accumulator가 동시에 활성일 때 **tree 하나의 work만** 과금 (과소)
- lane이 일부만 활성일 때 **full-width tree 비용**을 과금 (과대)

### 확정 계약

```text
LATENCY = 가장 깊은 live tree의 depth   (accumulator들은 병렬 reduce)
        = ceil(log2(mac_width))          단, 2개 이상 unit이 활성일 때
        = ceil(log2(final_accumulator_lanes))  단일 unit일 때
ENERGY  = 모든 live tree의 work 합
        = (units-1)*(mac_width-1) + (final_accumulator_lanes-1)
```

`mac_lane_state_t`에 `lane_width`를 추가하고, pe_t가 resolve된 lane state를 보관해 10개 call site를
전부 active-lane 형태로 옮겼다. 구조적 overload는 남겨 회귀에서 두 형태를 대조한다.

수기 회귀 (8 accumulator × 8 lane):

| 활성 lane | depth | work | 이전 (구조적) work |
|---|---:|---:|---:|
| 64 (8 unit) | 3 | **56** | 7 |
| 8 (1 unit) | 3 | 7 | 7 |
| 3 (partial) | **2** | **2** | 7 |
| 9 (1 full + 1 lane) | 3 | 7 | 7 |

기존 config은 `mac_reduction_energy`가 0이고 `mac_width`가 1이라 결과 불변이다.

## L9-b. Systolic drain / tile-boundary bubble

### 이전 상태

injection wavefront(operand skew)만 과금했다. 마지막 operand가 들어간 뒤 마지막 partial sum이 edge
accumulator까지 전파되는 **drain**이 없었고, 그 동안 array는 다른 weight residency로 넘어갈 수 없다.

### 확정 계약

```text
skew_fill_hops = (active_h - 1) + (active_w - 1)   // 가장 먼 PE까지의 operand skew (기존)
drain_hops     = active_h - 1                      // partial sum이 column을 따라 edge까지
```

drain이 **height**를 따르는 이유: partial sum은 column을 따라 흐른다. Gemmini의 RTL 측정치
(DIM=16에서 residency당 ~14 cycle)가 바로 이 column depth이며, 이것이 diameter가 아님을 뒷받침한다.

### 두 경로가 서로 다른 배수를 쓰는 이유 (중요)

- **calibrated `weight_fold_fill_cycle`**: per-PE weight **element**당 과금. Gemmini는 element마다
  array를 reload하므로 각 element가 독립적인 residency다. RTL 측정치이므로 우선하며, analytical drain을
  더하면 같은 stall을 이중 계상한다.
- **analytical drain**: array-wide weight **tile**당 한 번(= tile 경계) 과금.

처음에 analytical drain도 element당 과금해봤더니 scalesim(12×14 array, PE당 weight 352개)에서
fold bubble이 computation의 **10배**가 됐다. 이는 같은 tile 안의 연속된 weight 사이에서도 pipeline이
완전히 drain한다는 가정이고, 이 모델이 다른 곳에서 주장하는 pipelining과 모순된다. tile 경계당으로
바꾸면 99,000 cycle(computation의 2.8%)로, 방어 가능한 tile-boundary bubble이 된다.
따라서 analytical drain은 **하한**이며, 실제로 element마다 drain하는 설계는
`weight_fold_fill_cycle`을 calibration해야 한다 — 코드에 명문화했다.

Gemmini(calibrated)와 Eyeriss(spatial_arch, drain 없음)는 불변이다.

## L9-c. Non-zero format-IP 단가 fixture

format axis는 Phase 1(P4-1)에서 PE busy·network 합산·출력에 연결했지만 **모든 config의 단가가 0**이라
그 배선을 검증하는 회귀가 없었다. 단가와 buffering만 다른 3개 config을 추가했다.

[validation/format/check.py](../../validation/format/check.py) FM1~FM6:

- FM1: 3개 run 모두 axis > 0 (knob이 accounting에 도달)
- FM2: 단가에 **선형** — 8,256 → 82,560 (10배)
- FM3: double buffer에서 `busy = max(compute, axes, format)`. 단가 4에서는 transfer makespan
  (9,264) 뒤에 숨고, 단가 40에서는 axis가 **지배**해 busy == 82,560. 느린 format IP가 보이지 않거나
  통째로 숨겨지지 않음을 두 case로 구분한다.
- **FM4 수기 항등식**: single buffer에서 `busy = compute + max(axes) + format`
  = 3,518 + 9,264 + 8,256 = **21,038** (측정 일치)
- FM5: network format axis가 layer 합을 담는다 — **non-zero 값**으로 검증(0인 config로는 불가능했던 것)
- FM6: 느린 format IP가 critical path를 늘린다 (59,010 → 104,136)

## 검증 결과

| 검증 | 결과 |
|---|---|
| Unit/config validation | 통과, 229개 config |
| Gemmini RTL timing | 통과, MAPE 4.40%, max 7.86% (baseline 불변) |
| Eyeriss silicon latency | 통과, MAPE 4.26%, max 6.39% (baseline 불변) |
| Eyeriss GLB/DRAM traffic | max 133.94% (불변, informational) |
| Traffic T1~T10 | 8개 workload 통과 |
| Multi-chip MC1~MC9 | 통과 |
| Buffering LB-A~LB-E | 통과 |
| Bypass BP1~BP6 | 통과 |
| Asymmetric stats/timeline (L10) | 통과 |
| **DRAM DR1~DR6** | **통과 (신규)** |
| **Format FM1~FM6** | **통과 (신규)** |
| Energy E1~E5 | 통과 |
| Network rollup (L3) | 통과 |

## 남은 항목 — Phase 4

- **L4** memory-inclusive `Critical-path latency`의 외부 reference 확보와 calibration.
  현재 외부 golden은 compute schedule만 검증한다. **측정 데이터 확보가 선행되어야 하므로 코드
  수정만으로 닫을 수 없다.** 평가 문서가 완료 blocker로 분류한 항목이다.
- compute-bound / memory-bound / multi-chip acceptance gate 추가
- external traffic 오차(현재 max 133.94%) 축소 또는 해당 범위를 명시적 비지원으로 제한
- **L11** 지원 layer 확장(pooling/activation/normalization) 후 전체 network end-to-end latency 검증

## Phase 3에서 의도적으로 남긴 것

- **DRAM address 기반 bank/row mapping과 queueing**: timing 경로에 per-request address가 없으므로
  선택지 2(DRAMSim3 gate)로만 닫을 수 있다. 현재 계약은 출력에 명시했다.
- **systolic per-PE pipeline register의 개별 모델링**: skew fill과 drain으로 그 효과를 집계 수준에서
  표현하며, register 단위 event 모델은 cycle-accurate 영역이다.
- **element별 drain**: 위에 적은 이유로 하한을 택했고, 필요한 설계는 calibration으로 표현한다.
