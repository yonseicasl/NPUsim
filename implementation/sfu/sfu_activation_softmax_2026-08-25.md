# SFU 구현 기록 — element-wise activation + standalone softmax (2026-08-25)

[plan/plan_sfu.md](../../plan/plan_sfu.md)의 1차 구현. 사용자 우선순위에 따라 element-wise
경로(Phase 1–4의 직렬 모델)와 함께 Phase 7의 standalone softmax multi-pass 모델을 먼저
구현했다. DNN/LLM에서 자주 쓰는 activation(ReLU, Leaky, ELU, Sigmoid, Tanh, HardSigmoid,
HardSwish, GELU, SiLU)을 operation descriptor로 지원한다.

## 구현 범위

### components/sfu.{h,cc} — 새 컴포넌트

- chip당 SFU. `[sfu]` section은 opt-in — 없으면 어떤 수치도 변하지 않는다.
- Operation descriptor 테이블: bypass 여부, micro-op 여부, `<op>_pipeline_latency`,
  `<op>_initiation_interval`, `sfu_op_energy_<op>`, declared 플래그.
- Cycle 모델 (plan 공식): `chunks = ceil(valid/lanes)`,
  `cycles = setup + latency + (chunks-1) x II`. Unit 분할 시 latency = 최대 unit,
  event/energy = 합.
- Traffic 모델: ingress/egress element/byte/transaction을 accumulator precision으로 기록
  (SFU는 output cast 이전에 위치). Fused invariant: DRAM traffic 불변(검증됨).
- Energy 모델: `ops x op_energy + reads x read + writes x write + invocations x setup`,
  static은 `물리 SFU 수 x 최종 layer latency x sfu_static_energy` (always-on).
- Softmax microprogram: `vmax → vadd(sub) → exp → vadd(sum) → recip → vmul` 6-pass 직렬,
  row 단위 reduction은 lane-tree fold tail(`ceil(log2(min(lanes,N))) x II`/row) 포함.
- 순수 functional evaluator는 unit test 전용. FUNCTIONAL 경로의 owner는 Nebula
  `forward()` 하나(중복 적용 없음).
- Fail-fast: lanes/units 0, `queue_depth != 1`, placement/fusion contract 위반, II 0,
  linear에 profile 선언, `strict_profiles=1`에서 미선언 profile 사용.
- GELU/SiLU는 Nebula frontend enum이 없어 config/SFU API 레벨에서만 지원(계획서 Phase 6
  정책대로 frontend 추가 후 network에서 선택 가능).

### 이벤트 계약 (Phase 2의 1차 구현)

- Fused activation은 layer의 **최종(반복 스케일링 후)** 시점에 정확히 한 번 발생:
  `npu_t::apply_fused_sfu_activation()`이 `scale_serial_repetitions()` **이후** 실행되어
  reduction tiling/GLB repetition/psum spill 횟수와 무관하게 element 수 =
  network의 `B x K x P x Q` (mapping padding 제외). Connected는 `B x output_size`.
- Chip 분할: active chip에 균등 분배, 합 = layer valid output (unit test로 고정).
- Tile-identity 단위 event(부분 commit별)는 후속 단계 — 현재는 layer-granular 직렬
  invocation이며 보고서에 serial contract로 명시.

### Timeline 통합 (Phase 4의 직렬 모델)

`stats_t::apply_sfu_activation()`:

- `layer_latency += SFU busy` (queue_depth=1 직렬 output-path 자원).
- 모든 component leakage window를 새 latency로 rescale(`finalize_layer_timeline()`과
  동일한 규칙) — SFU가 느리면 다른 component의 static energy도 함께 커진다.
- `mac_available_cycle`/`utilization_mac` 재계산 — 느린 SFU가 MAC utilization을 낮춘다.
- Producer drain stall = SFU busy로 보고(직렬 contract의 명시).
- Streaming overlap/backpressure(Phase 5)는 미구현 — 보고서에 serial로 명시.

### Standalone softmax (Phase 7 1차)

- `npu_t::run_standalone_softmax()`: mapping section 불필요. rows = batch,
  length = output_size. chip 0의 SFU에서 실행(보고서에 명시).
- SFU-only `stats_t`를 즉석 생성(`record_sfu_only_layer`) — clock contract는 5개
  component frequency 검사로 채움, network rollup에 포함(timing scope의 지원 layer 수
  증가), rollup의 config-contract 필드/boundary depth min은 오염시키지 않음
  (`sfu_only_layer`, `network_rollup_mapped_seeded`).
- Scope 명시: operand streaming(DRAM/GLB↔SFU) traffic은 미포함 — 보고서에 바이트 수와
  함께 "NOT included" 문구 출력. materialized 회계는 후속 단계.

### Stats/report/energy 통합 (Phase 3)

- `stats_t`에 SFU counter/traffic/energy/calibration 필드 + network rollup 합산.
- Layer/network report에 `SFU (activation)` 블록: operation, units x lanes, valid
  elements, ops, invocations/chunks, busy/stall, tail lane utilization,
  ingress/egress(elem/bytes/txn), timing profile calibration, energy 5축.
- Energy summary에 opt-in 7번째 row `SFU (activation)`; `[sfu]` 없는 config는 6-row와
  `N of 6` undercount 그대로 (`energy_cost_schema_t::components_in_scope()`).
- Energy schema에 `sfu_read/write/setup/static_energy` + `sfu_op_energy_` prefix family
  추가; `[sfu]` section role. 필수 키 read/write — 누락 시 PARTIAL.
- Unpriced active event: 활성 op의 `sfu_op_energy_<op>` 미선언, read/write/setup 미선언
  모두 개별 이벤트로 보고되어 absolute energy/power qualification을 차단.
- `[sfu]` 없이 nonlinear activation이 실행되면 SFU 블록/energy summary에
  `Activation scope : NOT MODELED -- ...` 명시(수치는 불변; legacy power qualification도
  계획서의 호환성 정책대로 불변).
- Defaulted(1/1) timing profile 사용 시 `UNCALIBRATED` 표시; `strict_profiles=1`이면
  실패.

### Deprecated

- `pe_t::activation()` (하드코딩 ReLU): 호출 시 abort하도록 변경, call-site 재출현을
  `unittest/run_validation.sh`의 grep 회귀 검사로 차단.

## 검증

- `validation/sfu/run_sfu_validation.sh` (신규): formula/boundary, lane tail, profile
  독립성, energy 선형성, unit/chip 병렬성, element-count identity, softmax hand
  calculation(22,572 cycles 사례), serial merge, unpriced/uncalibrated 검출, activation
  name mapping(silent ReLU fallback 금지), functional evaluator vs reference, invalid
  config fail-fast(8종 subprocess).
- E2E: `gemmini_sfu x gemm256_relu` — relu 65,536 element → SFU 4,096 cycle이 critical
  path에 정확히 가산(32,748,190.6 → 32,752,286.6), DRAM traffic은 `gemmini` 실행과
  bit-identical(fused invariant). `gemmini_sfu x gemm256` — linear bypass 0-cost +
  standalone softmax(hand calc 일치). `gemmini_sfu_unpriced` — unpriced SFU event 보고.
- 기존 스위트 전부 통과: run_validation.sh, run_stats_timeline_validation.sh,
  run_timing_validation.sh(baseline 유지), run_multi_chip_validation.sh,
  validation/{unpriced,power,energy_schema}/check.py, run_full_sanitizers.sh(SFU 실행
  2종 추가).

## 코드리뷰 반영 (2026-08-26, /code-review high — 10건 검증됨, 전부 반영)

1. **Reduction-partition chip 분배** — CHIPS 행이 C/R/S를 분할하면 chip들은 같은 output의
   partial sum 복제본이므로 element를 나누면 SFU window가 과소평가됨. →
   `mapping_table_t::calculate_output_partition_chips()` 신설(K/B/P/Q chip factor 곱),
   `apply_fused_sfu_activation()`이 output-partition chip에만 분배. validation_test에
   C-split(1 partition) vs K×B-split(4 partition) 검사 추가.
2. **Partial mapping coverage** — mapping이 layer를 부분만 시뮬레이션하면 SFU도 그만큼만
   과금해야 함. → valid elements = 차원별 `min(mapped, layer)` (padding 제외 + coverage
   clamp 동시 처리). 리포트 label도 갱신.
3. **미지원 activation의 mid-run abort** — layer 시뮬레이션이 끝난 뒤 exit하면 완료 결과가
   버려짐. → `run()` 시작 시 전체 layer의 activation을 사전 검증(fail-fast를 시뮬레이션
   이전으로 이동).
4. **SFU-only rollup의 boundary/contract 오염** — softmax-only rollup이 init값 depth 1을
   측정치처럼 출력, softmax-first 순서에서 placeholder contract가 seed됨. → contract 복사
   guard를 `first_mapped_layer` 기준으로 통일, mapped layer가 없는 scope는
   `Buffer depth : n/a` 출력.
5. **Substring dedup** — "elu"가 "relu"에 포함되어 operation 목록에서 탈락. → 정확한
   segment 단위 dedup(`append_unique_segment`/`merge_unique_segments`).
6. **Multi-unit tail utilization** — 전체 share의 나머지가 아니라 가장 바쁜 unit의 최종
   chunk 기준으로 계산 (units=2, N=17 → 9/16). 테스트 추가.
7. **길이 1 softmax** — reduction pass(max/sum)는 op 0인데 cycle을 과금하고 미보정
   profile 플래그도 누락. → op이 0인 pass는 cycle/chunk/traffic/setup 전부 생략 (모든
   cycle에 event source 보장). busy=17(4-pass) 테스트 추가.
8. **Reduction-tree fold 비용** — fold step은 직렬 데이터 의존이므로 II가 아닌
   **pipeline latency**로 과금. latency≠II 케이스(29 cycle) 테스트 추가.
9. **Linear bypass의 sfu_active** — bypass layer는 dynamic event가 없으므로
   `sfu_active`를 세우지 않음(leakage는 always-on 모델대로 유지, sfu.h 문구 명확화).
10. **pe_t::activation() 잔재** — runtime abort stub과 약한 grep 대신 선언 자체를 제거
    (재도입 = 컴파일 에러), grep 패턴을 무인자 호출/선언 전반으로 확장.

## Phase 5 + Phase 7 확장 구현 (2026-08-26)

### Phase 5 — Streaming overlap / backpressure

- `[sfu] queue_depth`가 실제 계약이 됨: 1 = 직렬(기존과 bit-identical), ≥2 = streaming.
- Streaming 시 SFU는 `finalize_layer_timeline()`의 per-tile 타임라인에 **6번째 stage**로
  참여한다(`pipeline_timeline_cycles`, boundary depth = queue_depth). Fused 이벤트는
  `stats_t::set_sfu_activation()`으로 `scale_serial_repetitions()` **이전에** staging되고,
  타이밍 통합은 finalize가 소유한다(leakage window/MAC availability가 한 번에 정합).
- Final-output release time: output tile은 reduction group의 **마지막 tile**에서 commit
  된다 — SFU per-tile 비용은 `datatype_repetitions[OUTPUT]`개의 output group 마감 tile에
  배치(`output_repetition_tiles`). Reduction repetition만 있는 layer는 streaming이어도
  overlap이 거의 없다(활성화는 reduction 완료 후에만 가능).
- Attribution: `sfu_serial_cycle`(critical-path 노출) + `sfu_hidden_cycle`(producer 뒤에
  숨은 양) = busy 항등식; `sfu_stall_cycle` = 6-stage 실행의 PE-stage stall(= SFU queue
  backpressure). 리포트에 exposed/hidden/stall과 SFU bottleneck 여부 명시.
- Hand-calc 검증(`unittest/stats_timeline_test.cc::run_sfu_streaming_cases`):
  직렬 60 = 40+20 재현, fast-SFU 42 = T·p + s(2 노출/6 은닉), slow-SFU 90 = p + T·s
  (50 노출/30 은닉/stall 20), reduction-bunched 48 = T·p + S(은닉 0).
- E2E: `gemmini_sfu_stream × gemm256_relu` — 4,096 SFU cycle 중 16 cycle(drain tail)만
  critical path에 노출, 4,080 hidden, stall 0. 직렬 twin은 32,752,286.6으로 불변.

### Phase 7 확장 — softmax operand streaming + multi-chip 분배

- `[sfu] softmax_operand_residency = dram | glb` (기본 dram).
  - `dram`: DRAM→GLB staging→SFU, 결과는 역방향. `datatype_transfer_timing` +
    `pipelined_transfer_cycles`로 hop별 3-자원(소스 access/link/목적지 access) 회계.
    DRAM device·off-chip link 비용은 DRAM row로, GLB staging/feed/result/drain port
    비용은 GLB row로 energy summary에 귀속. row-activation·NoP hop은 미모델 명시.
  - `glb`: fused-schedule 시나리오. GLB feed/result port만 과금, `2×tensor bytes >
    GLB size`이면 fail-fast.
- Softmax row를 모든 물리 chip의 SFU에 분배(latency = busiest chip, work/energy 합).
- Layer critical path = ingress makespan + multi-pass window + egress makespan(직렬 계약
  명시); SFU leakage는 확장된 window 기준.
- E2E hand-check (`gemmini_sfu × gemm256`, 256×256 int8): DRAM access 393,216 cycle
  (= 131,072 access × 3), link 16,384 txn × 3 = 49,152, GLB 16,384 access(4 pass × 4,096),
  GLB energy 49,643.52 = 8,192×3.00 + 8,192×3.06, critical path 423,991 =
  200,711 + 22,572 + 200,708 — 전부 unit cost와 일치. softmax layer가 이제
  DRAM-bound(92.7%)로 정직하게 보고됨.

### 검증 (2026-08-26 추가분 포함 전부 통과)

validation/sfu(+invalid-residency, queue_depth=0 거부), stats_timeline(+streaming 4케이스),
run_validation, timing baseline, multi-chip, unpriced/power/energy_schema 게이트,
ASan/UBSan(+`gemmini_sfu_stream` 실행 추가).

## Phase 2 완전판 + Phase 8 RTL calibration (2026-08-26)

### Phase 2 — final_output_tile event

- Event source: `multi_chip_t::account_output_writeback_to_dram()` — output이 OUTPUT
  precision으로 element당 정확히 한 번 commit되는 경계(RE1/DR6/T1)이자 final cast가
  과금되는 지점. 이 호출 = `final_output_tile` event 1회
  (`final_output_tile_events/elements` 카운터, layer마다 reset).
- 총 event 수 = live pass의 commit 수 × output-datatype repetition. **Identity gate**:
  committed elements(padding 포함, mapped 기준) == mapped output volume일 때만 event
  단위 모델을 사용, 불일치 시 layer-granular 단일 invocation으로 폴백하고 리포트에
  MISMATCH 명시(plan 게이트 1: 불확실한 reduction-completion event는 timing에 연결
  금지).
- `sfu_t::elementwise_invocation(op, valid, commit_events)`: element를 event에 균등
  분배(최대 2가지 event 크기), event마다 자체 setup + pipeline fill, 직렬 drain.
  chunks/transaction은 event별 tail 반올림 반영, tail utilization은 최악 event 기준.
- Spill/reload 불변은 구조적(psum 트래픽은 commit 경계를 지나지 않음).
- E2E: gemm256_relu — 1 commit × 256 output reps = 256 events, identity OK, relu
  (L=II=1)는 event fill이 0이라 busy 4,096/critical path 불변(예측대로), invocations
  256. Leaky처럼 L>II인 op은 event당 fill이 정직하게 누적(단위 테스트 8 vs 5 케이스).

### Phase 8 — NVDLA nv_small SDP RTL calibration

- 기존 NVDLA fidelity 인프라(nvdla/hw@1a65f1f + Verilator 5.022 EMULATION 하네스,
  validation/nvdla/PROVENANCE.md)를 재사용해 공식 nv_small `sdp_*` trace 12개를
  CRC-PASS로 실행 (`validation/sfu/run_sdp_rtl_traces.sh`).
- 추출/분석 `validation/sfu/calibrate_sdp.py` — cycle 경계: 마지막 SDP-block
  D_OP_ENABLE(0x8008/0x9038) accept → 마지막 DBB output write. Regime 인식 분석
  (혼합 회귀 금지):
  - **코어 처리량 1 elem/cycle 정확히**: (32,851−85)/32,767 = 1.0000
    (sdp_4x1x8192 pass-through).
  - **Streaming skeleton fill 85 cycle**: N=1/N=8 동일.
  - **X1-ALU 경로 fill 151 cycle**: 448 − 297×1 (bs reg-operand). ReLU는 같은 X1
    ALU stage → `relu_pipeline_latency = 151`.
  - W=1/H-major surface(8.01 cyc/elem)와 mem-operand trace는 DMA-bound 참고치로만
    기록(residual 끼워맞춤 금지).
- `configs/accelerators/nvdla_small_sfu.cfg`: 위 constant + `profile_reference`
  provenance 선언. LUT-class op은 **의도적으로 미선언**(공식 ew_* trace 4개가 CRC
  FAIL — open item) → 사용 시 UNCALIBRATED 표시. Energy는 측정치 없음 → 키 미선언
  → active event가 UNPRICED로 보고되고 absolute qualification 차단.
- Provenance 기계화: `[sfu] profile_reference` → spec/리포트에
  `Timing provenance : <ref | NOT DECLARED ...>` 출력. 상세는
  `validation/sfu/SDP_CALIBRATION.md`.

## Phase 6 — Operation 확장 (2026-08-26)

### Nebula frontend: GELU/SiLU (item 4)

- `ext/nebula` 로컬 커밋 `de5d7ee` (branch `npusim`): `GELU_ACTIVATION`/`SILU_ACTIVATION`
  enum + `activation_type_str` + `gelu/silu_activation()` (GELU는 tanh 근사, SiLU는
  x·σ(x)) + `layer_t::activate()/gradient()` dispatch. **Forward/inference 전용** —
  `gradient()` 인터페이스는 post-activation output만 받아 두 함수의 도함수를 복원할 수
  없으므로 gradient stub은 명확한 메시지와 함께 abort(조용한 학습 오염 방지).
- 이제 network config에서 `activation=gelu`/`activation=silu`가 선택 가능 — E2E:
  `gemmini_sfu × gemm256_gelu` → SFU gelu invocation, 256 commit × (8+15) = 5,888
  cycle, op energy 65,536×0.9 = 58,982.4 (hand-check 일치).
- 주의: `./npusim.sh build nebula|ext|all`은 `origin/npusim`으로 re-pin하므로 로컬
  커밋(ae24ba1, de5d7ee)을 upstream(yonsei-icsl/nebula npusim branch)에 push하기 전에는
  해당 빌드 타깃을 피할 것.

### Loggy (item 1 마무리)

- `SFU_OP_LOGGY` (2σ(x)−1, transcendental class) — Nebula의 loggy가 SFU로 매핑되는
  마지막 element-wise activation. 계획서 미지원 목록(hardtan/lhtan/plse/ramp/relie/
  stair)은 여전히 fail-fast.

### Approximation / precision profile (item 2 + 완료 조건)

- `<op>_approximation = exact|piecewise` (ALU-class: relu/leaky/hsigmoid/hswish/
  vmax/vadd/vmul) 또는 `lut|polynomial|piecewise` (transcendental: elu/sigmoid/tanh/
  loggy/gelu/silu/exp/recip). 클래스 밖 조합·linear 선언은 **fail-fast**;
  `strict_profiles=1`이면 활성 transcendental op의 approximation 미선언도 실패.
  스펙 출력에 op별 `[lut]` 등 표시.
- `[sfu] profile_precision = <format>`: profile이 특성화된 피연산자 정밀도.
  `parse_data_format`으로 정규화·검증(오타 fail-fast), 런타임 **accumulator format**과
  비교 — strict면 불일치 시 실패, 아니면 모든 invocation `timing_calibrated=false` +
  리포트 note("characterized for int8 but runtime accumulator is fp32 — UNCALIBRATED").

### Nebula reference 비교 (item 3)

- `validation/sfu/sfu_test.cc`가 이제 `ext/nebula/common/activations.cc`를 **직접
  링크**해 재유도 공식이 아닌 실제 Nebula 구현과 element 단위 비교: relu, leaky,
  sigmoid(=logistic), tanh, hsigmoid, hswish, loggy, gelu, silu. Nebula elu는 구현
  버그(x>0 → x·x)로 명시적 제외.

### 검증 (Phase-6 추가분)

- 신규 테스트: Nebula 직접 비교, loggy/gelu 매핑, precision mismatch(비-strict 
  UNCALIBRATED / 일치 시 calibrated 유지), fail-fast 4종 추가(invalid-approximation,
  invalid-linear-approx, strict-precision-mismatch, strict-missing-approximation).
- validation_test.cc: 모든 shipped config의 `*_approximation` 값 검사.
- 전체 스위트 재통과 (nebula 헤더 변경으로 npusim clean 재빌드 후): run_validation,
  stats_timeline, timing baseline, multi-chip, unpriced/power/energy_schema,
  NVDLA compare(HOLDOUT 0.55% 유지), ASan.

## 잔여 세부 개선 반영 (2026-08-26)

### LUT trace CRC FAIL — 근본원인 해결

- **nv_small에는 LUT/EW 엔진이 존재하지 않는다** (`spec/defs/nv_small.spec`:
  `SDP_LUT_DISABLE`, `SDP_EW_DISABLE`). 증거: 하네스 FAIL 메시지의 첫 값은 *메모리에서
  계산된* CRC인데, 서로 다른 입력 텐서·형상의 4개 trace가 전부 동일한 계산 CRC
  (`c43621c5`)를 냈다 — 입력-무관 출력 = 프로그래밍된 EW/LUT stage가 없는 datapath의
  converter fallout. trace.bin의 golden은 바이트 수준으로 올바름을 확인(패킹/비교기
  무죄). stock 바이너리 대조는 기존 RAMDP 병리로 hang이라 판정 불가, spec 정의만으로
  확정. 상세: `validation/sfu/SDP_CALIBRATION.md`.
- 귀결: sigmoid/tanh/exp calibration은 nv_small 재실행이 아니라 **EW/LUT 활성 spec
  (nv_medium/large) 빌드**가 필요 — 실제 후속 작업으로 기록.

### `supported_ops` 아키텍처 allowlist (신규 계약)

- `[sfu] supported_ops = linear,relu,leaky`: 목록 밖 op은 **존재하지 않는 실행
  유닛**으로 취급되어 up-front fail-fast (activation 사전 검증 + softmax micro-op 5종
  검사 + invocation 방어 검사). linear는 항상 암묵 포함, 미선언 토큰은 config 오류.
  `nvdla_small_sfu.cfg`에 적용 — sigmoid 네트워크는 시뮬레이션 전에
  "execution unit does not exist"로 거부됨(E2E 확인).

### Softmax groups (axis 세분화)

- Nebula `d0a8967`: `[softmax] groups`가 실제로 동작 (batch×group span별 독립 정규화,
  비약분 groups는 fail-fast; 이전엔 파싱조차 안 되어 groups=4가 조용히 1로 동작).
- NPUsim: rows = batch×groups, length = size/groups. E2E `gemm256_g4`:
  1024×64 → multi-pass 28,812 cycle (hand-calc 정확 일치).

### Softmax streaming row-activation 비용

- dram이 row model을 선언한 경우(`row_buffer_size` 등) operand stream 양방향에
  open-page activation을 과금 — dram_t와 동일 해상도(tRC 보정 시 tRC, 아니면 flat
  row_miss; bank 병렬성은 latency만 완화, energy는 전액). transfer axis에 귀속,
  event 수는 `row_activation_events`로 흘러 `row_miss_energy` unpriced 검사와 연동.
- E2E hand-check: row 256B×8bank, 65,536B/방향 → 512 activations, 방향당
  32×(7+3)=320 cycle → critical path 423,991→424,631 정확 일치. Scope note에서
  "row-activation 미모델" 문구 제거(NoP hop 제외는 유지).

### SDP energy — 실행 불가 판정과 프로토콜 문서화

- `ext/nvdla/syn/`은 Synopsys DC 스크립트를 제공하지만 이 머신에 dc_shell/liberty가
  없어 실행 불가. 향후 프로토콜(합성 → VM_TRACE=1 하네스로 PASS trace SAIF/VCD →
  switching power에서 per-element energy)을 SDP_CALIBRATION.md에 명시. 그때까지 energy
  키 미선언 → UNPRICED 기제가 absolute 주장 차단.

### Tile-identity timeline 배치 — 해상도 한계 명시

- Commit 이벤트의 repetition 내부(sub-tile) 배치는 공유 tile clock의 해상도 아래에
  있어 현 timeline 구조로는 표현 불가 — 허위 정밀도를 만들지 않고 리포트의 commit
  라인과 본 문서로 한계를 명시하는 것으로 종결. (Streaming 릴리스는 계속
  output-repetition group 마감 tile 기준.)

## LUT-class RTL calibration — nv_small_256_full bring-up (2026-08-28)

- 가장 작은 LUT/EW 활성 spec인 `nv_small_256_full`을 기존 하네스 위에 bring-up:
  tmake 대체 드라이버(`nvdla_tmake_driver.py`, perl YAML 부재 우회) + `.f` 생성기
  (`gen_verilator_filelist.py`, 모든 RAM/FIFO를 -v 선등록) + 로컬 패치 3건
  (config-불변 `NVDLA_CFGROM.vh` 복사, 미컴파일 CVSRAM 블록의 `AXI_WDATA_MKPTR` 수정,
  트레이스 패커 SEC_MEM 허용). 전 과정 `SDP_CALIBRATION.md`에 기록.
- **nv_small에서 실패했던 공식 ew_* 트레이스 4종 + 대형 LUT 트레이스가 전부 golden
  CRC PASS** — 근본원인(유닛 부재) 판정의 폐루프 확인.
- 측정 상수 (spec 교차검증 정확 일치):
  - X1(BS/BN/ReLU) 경로 **2 elem/cycle** (0.5001 측정 = SDP_BS/BN_THROUGHPUT_2),
    stream fill 81, X1 fill 119.
  - **EW/LUT 경로 II = 1 elem/cycle 정확** ((453−156)/297 = 1.0000 =
    SDP_EW_THROUGHPUT_1), **fill = 156 cycle, LUT 테이블 내용 무관**
    (exp/le_lin/lo_lin 전부 453). LUT lookup 자체는 fill +6 cycle, throughput 영향 0.
- 산출물: sigmoid/tanh/exp/recip용 ready-to-declare `[sfu]` profile
  (SDP_CALIBRATION.md; GELU/SiLU는 다중 프리미티브라 단일 측정으로 미커버 명시).
  `nvdla_small_sfu.cfg`는 불변(nv_small엔 여전히 유닛 없음), 주석으로 참조 연결.
- 부수 해결: nv_small `sdp_8x8x32_bypass`의 "short read"도 같은 패커 assert 계열로
  근본원인 판정.

## 잔여 항목 (plan 기준)

- ~~Phase 2 완전판~~: 완료 (2026-08-26; commit-boundary `final_output_tile` event +
  identity gate. per-tile 릴리스의 timeline 배치는 여전히 repetition-group 단위 근사).
- ~~Phase 5~~: 완료 (2026-08-26; queue_depth 기반 streaming overlap/backpressure).
- ~~Phase 6~~: 완료 (2026-08-26; Nebula frontend GELU/SiLU + loggy + approximation/
  precision profile + Nebula 직접 비교. GELU/SiLU backward는 인터페이스 제약으로
  forward 전용 — abort로 명시).
- ~~Phase 7 확장~~: operand streaming(dram/glb residency)·multi-chip row 분배 완료
  (2026-08-26). 남은 세분화: softmax axis/groups(현재 Nebula 의미론 = batch당 전체
  output_size 1회), softmax streaming의 SFU-pass와의 overlap(현재 직렬 계약),
  NoP hop·row-activation 비용.
- ~~Phase 8 (timing)~~: nv_small(X1) + **nv_small_256_full(LUT/EW) calibration 완료**
  (2026-08-28). 남은 것: SDP energy는 DC/liberty 툴링 확보 시 문서화된 프로토콜로,
  Gemmini activation microbenchmark. `gemmini_sfu*.cfg`는 여전히 PLACEHOLDER로 명시.
