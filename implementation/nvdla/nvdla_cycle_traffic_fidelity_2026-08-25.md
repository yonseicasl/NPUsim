# 2026-08-25 NVDLA cycle/traffic 정합성 — CRC 근본원인 해결과 첫 전수 비교

기준 문서: `assessment/NPUSIM_EYERISS_GEMMINI_NVDLA_FIDELITY_ASSESSMENT_2026-08-25.md` §5, §8 P0.

## 한 줄 요약

P0의 blocker였던 multi-Atomic-C CRC mismatch를 **계측으로 근본원인을 찾아 수정**했고
(TSMC E2 RAM 모델의 `#1` 펄스 로직이 Verilator에서 쓰기를 유실 — RTL 무변경,
`+define+EMULATION`으로 해결), 공식 nv_small dc trace **13/13 전부 CRC PASS**로
golden set을 완성했으며, 현재 HEAD의 NPUsim과 전수 비교해
**compute-bound conv 7종에서 fitted parameter 0개로 compute-schedule MAPE 6.8%
(max −15.3%)**, weight/output DRAM traffic **정확 일치**를 얻었다.

## P0 항목별 상태

| P0 항목 | 상태 |
| --- | --- |
| 1. worktree를 현재 HEAD로 | **완료** — config/mapping/validation 이식, 재현 확인 |
| 2. multi-Atomic-C CRC 해결 | **완료** — 근본원인 증명 + 수정 (아래) |
| 3. current-head 재현 | **완료** — 이식 직후 +14.2%/traffic 일치 재현 |
| 4. workload 수집 | **완료** — 14 trace golden (기존 1 → 14) |
| 5-6. calibration/holdout | 미착수 (이제 가능해짐 — 후속) |

## CRC 근본원인 — 가설이 아니라 계측으로

전 단계 기록은 "C>8이면 실패, 원인 미상"이었다. 실제로 판정한 과정:

1. **실패 재현 후 로그 분해**: DBB 읽기 8,064 beat 전부 실데이터 도착(weight
   6,912 + input 1,152 = 전량 정확히 1회), 출력은 SDP가 프로그래밍된 cube
   기하(라인/서피스 stride, beat 수)를 **정확히** 지키며 **값만 전부 0**.
2. **VCD 계측 빌드**(`--trace`, 별도 Mdir): CSC→CMAC operand 버스에서
   `sc2mac_wt_a_data = 0` (8/8 샘플), `sc2mac_dat_a_data`는 실데이터,
   handshake 864회 정상. CDMA weight 로더 내부에서 `dma_rd_rsp_vld` 6,912회
   발화하는데 64-bit 데이터 스테이지는 전 구간 0.
3. **모델 대조**: 죽는 경로(MCIF read-response FIFO)는 `RAMDP_*_E2` 모델 —
   `assign #1` 지연 체인 11개로 write pulse를 성형하고 transparent latch로
   주소/데이터를 잡는다. Verilator가 정확히 이 array들에 MULTIDRIVEN 경고.
   동작하는 경로(CBUF)는 `RAMPDP_*_D2` — `#1`이 **0개**, 순수
   `always @(posedge WECLK)`.
4. **작은 trace가 통과했던 이유**: FIFO의 near-empty flop bypass. 경계는
   "C>8"이 아니라 "**RAM 경로로 넘어갈 만큼의 in-flight beat 수**"였다.
5. **수정**: `+define+EMULATION`. 이 매크로는 RAM 모델 파일 93개에서만
   소비되고 NVDLA RTL 본체에서는 0회 — NVIDIA가 같은 파일에 심어둔 단순 동기
   모델로의 스위치라 RTL 의미 불변.
6. **수정으로 증명**: 실패했던 run들의 tick 수를 **비트 단위로 재현**하며
   CRC만 PASS로 바뀐다(dc_6x8x192: 53,724 tick 동일). 즉 결함은 값 전용,
   타이밍 무영향 — 기존 169-cycle baseline도 불변.

## Golden set

14 trace(dc 13종 + sdp 1종) 전부 PASS. `validation/nvdla/golden_rtl_full.csv`.
물리 DBB traffic이 **전 trace에서 닫힌형 법칙과 정확히 일치**:

    input  = W·H·ceil(C/8)·8      (Atomic-C 패딩)
    weight = ceil(K·C·R·S/128)·128 (128B fetch granule)
    output = Wo·Ho·ceil(K/8)·8    (Atomic-K 패딩)

## NPUsim 비교 (현재 HEAD, fitted parameter 0)

trace 치수는 이름이 아니라 **register 프로그래밍에서 해독**했다(이름과 실제가
다르고 — dc_13x15x64_5x3은 실제 R=3,S=5 — don't-care 필드에 junk가 흔하다:
통과하던 1x1x8조차 dilation 4x14가 프로그래밍돼 있으나 1x1 kernel이라 무의미).

표현 불가 4종은 사유와 함께 제외(dilation 2, 비대칭 stride 1, R>H 1). nebula가
`filter_height/filter_width/padding_height/width`를 이미 지원해 비정사각
kernel은 패치 없이 표현했고, 비대칭 pad는 총량 분할 + 홀수분 입력 1행 확장으로
인코딩했다(왜곡은 enc/vol 열로 상시 표시).

| 축 | 결과 |
| --- | --- |
| cycle (compute-bound 7) | compute-schedule **MAPE 6.8% / max −15.3%**, 전부 음수(미모델 per-stripe overhead), K=270 대형에서 **−0.5%** |
| cycle (memory-bound 2) | critical path −76.9%/+32.4% — 고정 ~33-cycle DBB latency·engine fill 미모델 (informational) |
| DRAM weight/output | 선언 인코딩과 **정확 일치** (7/8; K=1 trace의 출력은 8B beat 최소단위) |
| DRAM input | **KNOWN-OPEN** — legacy row R/S는 halo 계약 밖(`input_halo_reuse()` early-out)이라 0.61~0.83×로 미달. E20-3급 모델 확장 항목 |

## config 정정 (결과 fitting 아님)

`[dram] read_cycle=64`는 `PRIMARY_MEMIF_LATENCY_64`의 **오독**이었다(FIFO 크기
산정용 파라미터). golden 환경의 실제 메모리(`nvdla.cpp` AXIResponder)는
**고정 ~33 cycle latency + 1 beat/cycle 스트림**(r0_fifo 32-깊이 파이프,
매 cycle 요청 수락)이다. beat당 점유 1로 정정하고
`max_outstanding_requests=128`을 선언했다. 미모델로 남는 것은 1회성 ~33 cycle
fill — micro trace에서만 유의미하고 config 주석에 남겼다.

## 남은 것 (README known-open과 동일)

1. input 축 halo 계약 확장 (legacy R/S)
2. memory-bound latency/fill 모델
3. per-stripe overhead — **calibration/holdout 분리로만** 닫는다 (plan §6.2/6.3)
4. CACC-as-temporal-buffer (critical-path 축의 psum round-trip 제거)
5. energy 축 전체 (plan Phase 5 — 미착수 유지)

## 산출물

- `validation/nvdla/{golden_rtl_full.csv, extract_golden.py, gen_nvdla_workload.py, compare_full.py, README.md, PROVENANCE.md}`
- `configs/accelerators/nvdla_small.cfg` (정정 주석 포함), `configs/networks/nvdla_dc_*.cfg` 9종, `configs/mappings/nvdla_small/*` 9종, `models/weights/nvdla_dc_*.wgh`
- RTL: `ext/nvdla/outdir/nv_small/verilator_emu/` (EMULATION binary + 14 run 로그), `verilator_trace/` (VCD 계측 빌드)
- 구 스크립트(compare_cycle/compare_traffic)와 n=1 golden CSV는 compare_full/golden_rtl_full로 대체 (이력은 git과 worktree에 보존)

## 커밋 시 주의

`.gitignore`가 `configs/`, `models/weights`를 통째로 무시한다(기존 config들은
tracked라 override). 이번에 새로 만든 `configs/accelerators/nvdla_small.cfg`,
`configs/networks/nvdla_*.cfg`, `configs/mappings/nvdla_small/`,
`models/weights/nvdla_*.wgh`는 커밋할 때 `git add -f`가 필요하다.
`validation/nvdla/`와 이 문서는 무시 대상이 아니다.

---

# 2026-08-26 후속 — input halo 계약 확장 (known-open 1번 종결)

## 무엇을 바꿨나

`mapping_table_t::input_halo_reuse()`의 legacy-row FILTER early-out을 제거하고 계약을 확장했다.

* **union extent에 legacy R/S를 접는다**: `full_h = (P_tot−1)·stride + R_tot`
  (R_tot = base_R × legacy_R). legacy filter가 1이면 기존 식으로 정확히 환원된다.
* **양방향 적용, 비대칭 guard** (`stats_t::scale_serial_repetitions`):
  - 축소(coalesce down)는 **retention 주장**이므로 기존대로 ring working set이 GLB에
    들어갈 때만.
  - 확대(union까지 올림)는 **무조건**. union은 loop 순서·버퍼 크기와 무관한 물리적
    하한이다 — window가 만지는 모든 byte는 최소 1회 인터페이스를 건넌다.
* legacy filter가 있으면 working set = union 전체(레벨 내 stripe 순서를 모델링하지
  않으므로 partial-ring retention을 주장할 수 없다 — nv_small CDMA/CBUF의 실제
  동작과도 일치).
* `active`가 `unique < replicated`에서 `unique != replicated`로 — 미달 방향도 계약
  안으로 들어온다.

## 왜 old early-out이 보수적이지 않고 틀렸었나

R/S가 legacy row에 있으면 below-legacy tile에 kernel extent가 없어서, repetition
fallback은 OH×OW×C tile만 옮긴다 — **union보다 0.61~0.83배 미달**. RTL의 CDMA는
전체 cube를 정확히 1회 fetch한다. "적용하지 않으면 안전"이 아니라 미달 방향으로
틀리고 있었다.

## 검증

* **unittest 손계산 항등식 추가**: legacy-R/S 케이스에서
  replicated 10,752 < unique 16,384 — dc_13x15x64 **실측과 정확히 일치**.
* **NVDLA 9 trace: input 축 전부 exact** (0.61~0.83× → 1.00×). cycle 오차 불변
  (traffic 전용 수정임을 확인 — compute-bound MAPE 6.8% 유지).
* **Eyeriss/Gemmini 외부 baseline 비트 단위 불변**
  (`check_timing --check-baseline --check-traffic` exit 0) — 그들의 legacy row에는
  filter loop가 없어 계약 확장이 아무 것도 바꾸지 않는다. SCALE-Sim PASSED.
* 24 gate 중 23 통과. `knobs` KN9 실패 1건은 **이 변경과 무관** — 같은 tree에서
  진행 중인 별도 작업(`components/sfu.cc`, untracked)의 `strict_profiles` knob를
  undeclared-knob 스캔이 잡은 것. KN9는 정적 스캔(코드 read 키 − config 선언 키)
  이라 halo 변경과 교차하지 않으며, 해당 작업 쪽에서 분류해야 한다.
* T11의 "not needed" 문구를 확장 의미에 맞게 갱신
  ("per-repetition streaming already moves exactly the union").

## 남는 report 형식

up-방향 적용 시: `applied; 10752 -> 16384 input elements, working set N B union
lower bound, DRAM serialized A -> B` (down-방향은 기존 "fits GLB" 유지 — T11 regex
불변).

---

# 2026-08-26 후속 2 — memory-bound latency/fill 모델 (known-open 2번 종결)

## 측정 (validation/nvdla/measure_pipeline_constants.py)

end-to-end 창을 로그의 네 사건으로 분해했다: enable→첫 read 요청(launch), 첫
요청→첫 데이터(mem_latency), 데이터 스트림, 마지막 데이터→마지막 write(drain).

| 상수 | 값 | 출처 |
| --- | --- | --- |
| launch | **28** | conv trace 13/13에서 **정확히 상수** (구조 상수) |
| mem_latency | **32** | 14/14에서 정확히 상수 **+ harness 소스(AXI_R_LATENCY=32)에서 독립 유도** — 고정 latency 1 beat/cycle 스트림이므로 layer당 1회 지불 |
| drain | **77** | compute tail이 0인 유일한 trace(dc_1x1x8)에서 **n=1 calibration** — conv 파이프 깊이(CBUF→CSC→CMAC→CACC→SDP→WDMA) |

drain은 상수가 아니었다(33~292,497) — compute-bound trace의 drain gap은 compute
tail이 섞인 값이라 파이프 깊이가 아니다. 이를 스크립트 docstring에 명시했다.

## 선언과 규율

`[spatial_arch] layer_setup_cycle = 137` (+ `layer_setup_energy = 0` modeled-zero).
기존 knob 재사용 — **모델 코드 무변경, config 전용**. Gemmini의
`layer_setup_cycle=2270`과 같은 기계로, compute-schedule과 critical path 양쪽에
1회 가산된다.

규율: **dc_1x1x8은 holdout에서 나와 calibration 점으로 명시 전환**(비교표에 `cal`
표기, 집계에서 제외). 구조 상수 2개(launch/mem_latency)는 fitting이 아니라 측정이며
전 trace에서 정확히 상수임을 확인했다.

## 결과

| | 전 | 후 |
| --- | --- | --- |
| dc_1x1x8 critical path | −76.9% | **−18.3%** (cal 점) |
| compute-bound holdout 7종 | MAPE 6.8% / max 15.3% | **MAPE 6.7% / max 15.0%** |
| dc_4x1x8192 | +32.4% | +32.4% (per-stripe overhead bound — open item 3, latency 문제 아님) |

잔여 −18.3%는 크기가 아니라 **합성의 문제**다: RTL은 launch+latency가 스트림과
직렬인데 NPUsim의 setup은 PE stage 비용이라 DRAM 스트림과 겹친다. 상수를 줄여
이를 흡수하면 curve-fitting이므로 측정값 그대로 둔다.

traffic 3축·물리 법칙은 전부 불변(exact 유지). 이 변경은 nvdla_small.cfg와 비교
스크립트만 건드려 24개 gate·외부 baseline과 교차하지 않는다.

---

# 2026-08-26 후속 3 — per-stripe overhead calibration/holdout (open item 3 종결)

## 형태 규명 — 데이터가 모델을 고른다

setup 137 제거 후 잔차를 후보 count들로 나눠본 표에서:

* **32x26x76의 K16/K270 쌍이 결정적**: weight block 수가 17배 차이인데 잔차는
  22,210 vs 22,898로 거의 동일 → **잔차는 Kg 무관**. per-weight-block 계열 모델 전부 기각.
* stride-1/무dilation conv 7종에서 잔차/(Cg·R·S·Ho) = **3.40~4.82로 밀집**.
* 비단위 stride/dilation trace는 같은 축에서 9.4~29.4 — 다른 stripe 기계.
  degenerate shape(1×1 kernel C-heavy, H=1)는 또 다른 영역.

**모델**: overhead = B × ceil(C/8)·R·S·Ho — output row마다 각 reduction slice를
스트리밍할 때의 고정 bubble.

## 규율

* **형태 출처를 명시**: 형태 선택에 K16/K270 쌍을 썼으므로 **두 쌍 모두 calibration
  set에 넣었다**. holdout이 형태 선택에 오염되지 않도록.
* CAL = {13x15x64, 32x26x76 K16, 32x26x76 K270, 8x16x128},
  HOLD = {14x7x49, 24x44x14, 35x22x54}. B는 CAL에서만 원점 통과 최소제곱:
  **B = 4.356**.
* scope 선언: stride-1/dilation-1 direct conv. 밖의 trace는 informational.

## NPUsim 구현 (일반 knob으로)

* `[spatial_arch] stripe_transition_cycle` / `stripe_transition_energy` —
  `pe_array_t`가 파싱, count는 `mapping_table_t::stripe_transition_count()`
  (legacy C·R·S·P; layer마다 `set_psum_retention_scope`에서 기록).
* stats: layer setup과 같은 방식으로 **repetition scaling 후 1회** fold-fill에
  가산 → compute-schedule과 critical path 양쪽 반영. 미선언 config(u=0)는
  비트 단위로 아무 것도 안 바뀐다.
* house 규약 완비: energy schema 키 등록(`optional[4]`→`[5]` 확장), RE3 교차
  규칙(energy>0 && cycle==0 거부), unpriced-event 항목, modeled-zero 주석,
  report 라인("0.00 uncalibrated over 1680 transition(s) [modeled zero]").

## 결과 (시뮬레이터 자체로 재현)

| 구분 | trace | err |
| --- | --- | --- |
| **HOLDOUT** | 14x7x49 / 24x44x14 / 35x22x54 | **−0.8 / −0.3 / +0.6%** → **MAPE 0.55%, max 0.75%** |
| cal | 13x15x64 / K16 / K270 / 8x16x128 | +0.4 / −0.1 / −0.0 / +3.0% |
| cal(drain) | 1x1x8 | −15.7% (합성 잔차, 이전 기록 참조) |
| scope 밖 | 4x1x8192 | −64.2% informational |

holdout 기준 **6.7% → 0.55%** (max 15.0% → 0.75%).

## 회귀

23/24 gate(knobs 1건은 계속 무관한 `sfu.cc` WIP의 `strict_profiles`), unit test
OK, `check_timing --check-baseline --check-traffic` exit 0 — **Eyeriss/Gemmini
비트 단위 불변** (그들의 config는 knob 미선언 → u=0 경로), SCALE-Sim PASSED.

---

# 2026-08-26 후속 4 — CACC-as-temporal-buffer (psum residency 정정)

## 두 가지 수정

**1. `psum_must_leave_array()`의 같은-레벨 규칙을 순서 인지형으로.** 종전 규칙은
"레벨 내 순서를 모델링하지 않으므로 reduction·output 공존 시 보수적으로 spill"
이었다. 그런데 그 순서는 미지가 아니다 — **scheduler가 legacy-GLB offset 순서를
`[spatial_arch] pe_stationary`로 실제 생성한다**(scheduler.cc의 pe_array 분기).
OUTPUT_STATIONARY면 각 출력 tile의 reduction이 완결된 후 전진 → psum은 array
가장자리(nv_small: CACC)에 머물고 GLB를 건너지 않는다. 다른-레벨 중첩(Eyeriss
conv3: DRAM C=64 안에 GLB output loop)은 순서가 아니라 **중첩**이므로 무조건
spill — 불변. 두 경우 모두 unittest 손계산 항등식으로 고정했다.

**2. nvdla_small.cfg의 `pe_stationary`를 weight→output stationary로 정정.**
nv_small CSC는 output stripe마다 모든 (r,s,c-group) slice를 완결하고 전진한다 —
array↔CBUF 수준의 순서는 output-stationary다. stripe **내부**의 lane별 weight
유지는 `mac_stationary`(다른 레벨)가 담당한다. v1의 weight_stationary는 실기가
하지 않는 psum GLB 왕복을 만들어냈다.

## 검증

* **psum 왕복 소멸**: 8x16x128에서 psum legs **0 loads / 476 stores** — 476 =
  출력 tile 수와 정확히 일치(T12 항등식: 마지막 write-out만 tile당 1회).
* **compute-schedule·DRAM traffic·calibration 모두 비트 불변** — holdout
  MAPE 0.55% 유지, traffic 3축 exact 유지.
* Eyeriss/Gemmini **비트 단위 불변**(`--check-baseline` exit 0) — 그들의 map은
  같은-레벨 공존이 없고(across-level 또는 reduction-only), 규칙 변경이 닿지
  않는다. unit test·23 gate·SCALE-Sim 통과.
* knobs 실패가 1→3건이 됐으나 셋 다 동시 진행 중인 `components/sfu.{cc,h}`
  작업의 knob(`strict_profiles`, `profile_reference`,
  `softmax_operand_residency`) — 이 변경과 무관.

## 후속 항목으로 넘긴 것 — array streaming 세분성

critical path의 잔여 5.0M 축(68,544 pass × 73 cycle)은 psum이 아니라
**position 단위 OS 순서의 정직한 비용**이다: pass마다 weight tile이 바뀌므로
array fabric이 출력 위치당 weight 64 beat를 재스트리밍한다. RTL은 이것도, WS
대안의 psum spill도 내지 않는다 — **stripe 단위 hybrid**(stripe 간 OS + stripe
내 weight 유지)이기 때문이고, 레벨당 stationary 하나인 NPUsim mapping으로는
표현할 수 없다. stripe 세분성(또는 array streaming의 파이프라인화)이 필요한
별도 모델 작업으로 README에 등재했다. compute-schedule 비교량은 영향이 없어
holdout 지표는 그대로 선다.

---

# 2026-08-26 후속 5 — array streaming 세분성 (critical path 사다리 완주)

## 원인은 하나의 부류, 다섯 겹

critical path의 ~70배 인플레이션은 "스트림을 op 단위로, 실기에 없는 직렬화 공유
매체 위에서 과금"하는 인공물이 다섯 레벨에 겹친 것이었다. 각각을 **구조 선언**으로
풀었다 (모든 knob 기본값 = 기존 동작; Eyeriss/Gemmini **비트 불변** 확인):

| 겹 | 사실 (출처) | 선언 | 8x16x128 crit |
| --- | --- | --- | ---: |
| (시작) | | | 5,003,743 |
| output 수집 | lane별 **전용 wire**로 tree 주입 (VCD 동시 발화) | `writeback_injection_parallel` — 시간만 중첩, 에너지/txn 유지 | |
| wt 분배 폭 | CSC wt 버스 = **64B/cycle** (RTL 포트) | intra-array `bitwidth 512`, `line_size 64:512:64` | 822,548 |
| array 패브릭 | dat/wt/accu **별도 버스** | `array_fabric_separate` — link 축 sum→max | 221,282 |
| PE operand | register가 compute와 lockstep (1 op/lane/cycle 지속) | `operand_streams_pipelined` — **단일 억제 지점**(`suppress_streaming_cycles`)이 cycle 뷰만 0 | 152,738* |
| GLB 포트/스트림 | **D_BANK가 dat/wt를 분리 뱅크로**, 동시 서비스; 연속 전송은 1 beat/cycle 파이프라인 | `[shared]`→`[separate]`(68/52/8 KB 정적 근사) + `fabric_separate` + `streams_pipelined` | **93,947** |

(* PE 억제와 GLB 첫 단계가 같은 측정 구간에 걸쳐 있어 사다리 수치는 적용 순서 기준)

곁들여 **상속 오류 정정**: `[shared] bypass=0:1:0`(weight GLB 우회)은 Eyeriss
템플릿의 잔재 — nv_small의 weight는 CBUF 상주(측정된 fetch-once + stripe당 재독)라
`0:0:0`으로 바로잡았다.

구현 세부: `fabric_separate`류는 **config-static**(per-layer reset이 지우지 않도록
in-class 기본값) — 처음에 reset 블록에 초기화를 뒀다가 parse 직후 지워지는 것을
probe로 잡아 고쳤다. PE 억제는 흩어진 10여 accrual을 개별 gating하지 않고 stats
수집 직전 단일 메서드로 — dense/sparse/descriptor 경로가 이 의미에서 어긋날 수
없다.

## 결과

* **critical path, 전 9 trace: −11.7% ~ +22.2%, MAPE 11.7%** (informational).
  잔차가 일관되게 양수 = 남은 미모델 중첩(보수 방향).
  1x1x8 −6.3%(최초 −76.9%), 4x1x8192 −11.7%(stripe-bound로 밀어뒀던 trace가
  crit 축에서는 유의미해짐), K270 +3.4%.
* **compute-schedule·DRAM traffic·calibration/holdout 전부 비트 불변** —
  holdout MAPE 0.55%, traffic 3축 exact, 물리 법칙 14/14 유지.
* 회귀: 23/24 gate(knobs 3건은 동일한 SFU WIP 항목), unit test OK,
  `--check-baseline` exit 0, SCALE-Sim PASSED.

## 남는 것

crit 잔차(+4~22%)의 후보: stage 경계의 잔여 fill 합성, stripe 세분성의 잔재
(GLB access 축 71,280 = in+out fill 성분), DRAM/NoP 축의 고정 latency 합성.
전부 informational crit 축의 개선 항목이며, 검증 지표(compute-schedule holdout)와
분리되어 있다.

---

# 2026-08-26 후속 6 — energy 축 (P2): component-level 외부 calibration

합성 툴체인이 없는 환경에서 정직하게 가능한 범위 = gemmini_cacti22와 같은
**component-level 외부 calibration + 나머지의 명시적 거부**.

## Calibrated (전 항목 재현 가능)

| 축 | 출처 | 값 (22nm) |
| --- | --- | --- |
| CBUF SRAM | CACTI 7 (pinned commit 재빌드) × **RTL macro 기하**(nv_ram_rws_256x64: 2KB/64b) | read 0.6911 / write 1.0221 pJ; 512b wt line = 8 bank 병렬 → ×8 |
| CACC assembly | CACTI(nv_ram_rws_16x256: 512B; 256b 버스는 CACTI가 거부 → 128b 행을 **per-byte 소비**) | accumulator_spill_energy = 0.0530 pJ/B (read/write 평균; reload/spill byte가 근사 동일해 평균이 합에 충실) |
| CBUF leakage | 파티션당 macro 수 × 0.8038 mW @1GHz | 27.33/20.90/3.22 pJ/cycle |
| DBBIF | DRAMsim3 DDR3-1600 (기존 reference 행) | 452.25 pJ/8B beat |

reference CSV 행은 `gen_nvdla_sram_rows.py`로 재생성 가능(도구/커밋/오버라이드
전부 PROVENANCE). **의도적 미선언**(unpriced 기계가 발화): MAC compute, operand
register, CMAC adder tree, stripe/setup 제어, NoC wire, CACC delivery(544B 기하
CACTI 거부), MCIF FIFO, clock network. RTL은 TSMC16FFLR, 22nm는 최근접 CACTI 노드
(gemmini_cacti22와 동일한 명시).

## 곁들여 잡은 결함 2건

1. **schema 분류 로직**: optional 키만 교정 선언된 컴포넌트가 자기 주석("declared
   optional DOES make it costed")과 달리 NOT MODELED로 찍혔고, optional 스캔
   루프 경계가 넓힌 배열(5)이 아닌 4로 남아 있었다. 수정 후 MAC 행이
   "PARTIAL (UNDERCOUNT): missing computation_energy, ..."로 정직하게 읽힌다.
   기존 config 주석 전부 비트 불변(회귀로 확인).
2. **자체 config 버그**: `[separate]` GLB static을 넣을 자리를
   `memory_type = separate` 검색으로 찾다가 PE LB의 같은 키에 먼저 걸려 GLB
   leakage가 **PE당** 과금(×64)됐다 — 64×51.44 pJ/cycle 항등으로 잡아 정정.

## 검증

* **P4-10 신설**: 선언 단가 == reference 유도값(CBUF ×1/×8, CACC per-byte 평균,
  DRAM 행), wattage **거부**, ESTIMATED 표기, per-component UNDERCOUNT 가시성.
* breakdown이 물리적으로 타당: DMA-heavy 1×1 trace DRAM 72.6% ↔ K270
  compute-heavy 8.4%. DRAM 에너지는 beat 수 항등으로 정확(7,820×452.25).
* 25 스위트 중 24 통과(knobs 3건은 동일 SFU WIP), unit test OK,
  `--check-baseline` exit 0(스키마 수정이 기존 주석을 하나도 안 움직임),
  NVDLA holdout 0.55% 불변.

## 남는 것 (P2의 상한)

절대 total energy golden(RTL 전체 에너지)은 이 환경에서 만들 수 없다 — 합성
netlist power가 필요하며, 이는 평가서가 gemmini에도 동일하게 남긴 한계다. 모델은
그 부재를 숨기지 않고 wattage 거부 + UNDERCOUNT 표기로 말한다.
