# NVDLA RTL 기반 NPUsim fidelity 검증 계획

**작성일:** 2026-08-20  
**대상 metric:** cycle, traffic, energy  
**대상 구현:** NPUsim timing/traffic/energy model과 Open NVDLA RTL

**2026-08-21 개정 (코드 감사 + 실제 소스코드 수정 1건):** 7.2/8.4/9절의 "예상 NPUsim 보완
항목"을 실제 `components/`·`scheduler/`·`utils/` 코드 감사로 재검증했다. 16개 항목 중 절반
가량은 이미 architecture-neutral 메커니즘으로 존재해 NVDLA config/mapping만 있으면 되고
(7.2-A, 8.4-A), 1개(memory request concurrency / outstanding-request)는 구조적으로 불가능한
상태였던 것을 이번에 실제로 구현·검증했으며(8.4-B, `validation/memory_concurrency/check.py`),
나머지는 genuine gap으로 재확인하거나(7.2-B, 8.4-C) 애초에 "매핑 작성 책임"으로 재분류했다
(7.2-C). 4.2절은 이 문서 작성 이후 같은 세션에서 실제로 진행된 첫 RTL 데이터 포인트(n=1,
stale baseline)를 반영해 정정했다. 상세 내역은 각 절 참고.

## 1. 목적

이 계획의 목적은 동일한 convolution workload를 NPUsim과 NVDLA RTL에서 실행하고 다음 세 축의
정합성을 독립적으로 평가하는 것이다.

1. **Cycle fidelity** — NPUsim의 layer 실행 cycle과 cycle-exact RTL 실행 cycle 비교
2. **Traffic fidelity** — NPUsim의 계층별 traffic과 RTL interface/SRAM handshake에서 관측한
   실제 전송량 비교
3. **Energy fidelity** — NPUsim의 절대 energy와 합성 netlist, switching activity, SRAM 및
   DRAM energy model로 구성한 RTL implementation energy 비교

Energy의 최종 주장은 실측 silicon이 아니라 다음 범위로 제한한다.

> 명시된 technology library, voltage, frequency, process corner 및 power-analysis flow에서의
> NVDLA RTL implementation fidelity

NVDLA silicon 실측값이 없는 상태에서 `silicon absolute fidelity`라고 표현하지 않는다.

## 2. 초기 검증 범위

### 2.1 기준 RTL

초기 기준은 `nvdla/hw`의 `master` branch에서 특정 commit을 고정한 `nv_small` configuration으로
한다. 실험을 시작할 때 commit hash를 `PROVENANCE.md`와 모든 golden CSV에 기록한다.

`nv_small.spec`의 주요 특성은 다음과 같다.

| 항목 | 설정 |
| --- | --- |
| Precision | INT8 |
| MAC atomic geometry | Atomic-C=8, Atomic-K=8, 총 64 INT8 MAC |
| Convolution buffer | 32 banks × 8 bytes × 512 entries = 128 KiB |
| Primary memory interface | 64 bit |
| Weight compression | Disabled |
| Winograd | Disabled |
| Batch | Disabled |
| Secondary memory interface | Disabled |

근거:

- [공식 `nv_small.spec`](https://raw.githubusercontent.com/nvdla/hw/master/spec/defs/nv_small.spec)
- [NVDLA scalability parameters](https://nvdla.org/hw/v2/scalability.html)
- [NVDLA Hardware Manual](https://nvdla.org/hw/contents.html)

### 2.2 NPUsim 지원 범위

초기 검증은 NPUsim이 현재 timing cost를 부여하는 범위와 맞추어 다음으로 제한한다.

- INT8 dense direct convolution
- batch=1
- compression disabled
- Winograd disabled
- single primary memory interface
- convolution과 output write-back에 필요한 최소 SDP pass-through
- pooling, activation, normalization 및 다른 post-processing은 metric에서 제외

`nv_large`, FP16/INT16, Winograd, compression, batch 및 secondary memory interface는 `nv_small`
holdout gate를 통과한 뒤 별도 확장 단계에서 다룬다.

## 3. 검증 원칙

### 3.1 동일 실행 조건

다음 조건이 양쪽에서 같아야 한다.

- tensor shape, padding, stride, dilation, group
- datatype과 physical alignment
- dense/compressed storage mode
- CBUF 및 외부 memory capacity
- memory bandwidth, latency 및 outstanding request limit
- layer 시작·종료 경계
- output write-back 포함 여부

### 3.2 Calibration과 holdout 분리

모델 parameter를 정하는 calibration workload와 최종 fidelity를 평가하는 holdout workload를
분리한다. Holdout 결과를 본 뒤 unit cost나 timing constant를 다시 fitting하면 해당 결과는
검증값이 아니라 재보정값이므로 별도로 표시한다.

### 3.3 비교 경계의 명시

Cycle, traffic, energy 각각에 대해 포함 경계를 명시한다. 특히 다음을 혼합하지 않는다.

- CSB configuration 시간과 accelerator execution 시간
- testbench의 zero-time memory preload와 NVDLA DMA traffic
- NVDLA core energy와 외부 DRAM device energy
- logical tensor volume과 physical bus traffic
- CBUF와 CACC traffic

## 4. Phase 0 — 재현 환경 구축과 provenance 고정

### 4.1 고정할 항목

- NVDLA RTL repository URL과 commit
- `nv_small.spec` 원본 및 hash
- local RTL/testbench patch
- Verilator, Clang, Java 및 build tool version
- memory-model plusarg
- workload/trace/input file hash
- 입력 데이터 seed
- core/CSB clock frequency
- synthesis 및 power tool version
- standard-cell/SRAM library와 PVT corner

### 4.2 현재 환경 점검 결과

2026-08-20 현재 머신은 다음 상태다.

- Memory: 157 GiB
- Workspace filesystem 여유 공간: 약 145 GiB
- Clang: 18.1.3
- Java: OpenJDK 21
- Verilator: 설치 확인되지 않음
- Yosys: 설치 확인되지 않음
- 저장소에 `ext/nvdla`는 아직 없음

NVDLA 공식 build 문서는 오래된 Java/SystemC/Verilator/Clang 조합을 기준으로 하므로, 먼저
현재 toolchain에서 build 가능한지 확인하고 필요한 compatibility patch를 모두 기록해야 한다.
공식 문서는 Verilator build가 최대 20 GiB 이상의 메모리를 요구할 수 있다고 명시한다.

근거: [NVDLA Integration Guide](https://github.com/nvdla/doc/blob/master/doc/hw/v1/integration_guide.rst)

### 4.2-후속 (2026-08-21 업데이트, 이 문서 평가 시점)

위 4.2절은 작성 시점(2026-08-20 23:01) 스냅샷이며, **같은 세션 내에서 8분 뒤부터 바로 낡은
정보가 됐다.** 그대로 남겨두면 "아직 아무것도 못 했다"는 잘못된 인상을 주므로 사실 관계만
정정한다.

- `ext/nvdla`는 2026-08-20 23:09에 실제로 채워졌다 (`nvdla/hw` master 브랜치,
  commit `1a65f1f5b48268accaa47c95f95c2601918be095`). `nv_small.spec`의 실측 내용은 위 2.1절
  표와 정확히 일치한다 (Atomic-C=8, Atomic-K=8, CBUF 32×8B×512=128KiB, primary memif 64-bit
  등).
- Verilator는 시스템 `PATH`에는 없지만 `ext/chipyard/.conda-env/bin/verilator`
  (5.022, conda-forge)로 실제 동작이 확인됐다. Yosys는 여전히 미설치로 확인됐다 —
  R4(energy synthesis 리스크)는 그대로 유효하다.
- `nv_small` RTL을 이 Verilator로 빌드해 공식 INT8 trace(`dc_1x1x8_1x1x8x1_int8_0`)를
  실행했고, NVIDIA가 trace에 내장한 golden CRC와 **정확히 일치**하는 결과를 얻었다
  (`0x8f68a2ae`) — 4.3절의 "RTL output memory comparison PASS" 조건 중 하나를 실제로
  통과한 사례다.
- 이 결과로 `end_to_end_cycles`(RTL 169 vs NPUsim 첫 시도 193, +14.2%)와 세 stream(input/
  weight/output) 모두의 logical DBBIF traffic 완전 일치를 얻었다. Physical byte가 logical
  byte보다 큰 두 경우(weight 8B→128B, output 1B→8B)는 7.2절이 미리 지목한 atomic/64-bit
  alignment 항목과 정확히 일치한다.
- **단, 이 작업은 branch `worktree-nvdla-fidelity`(base commit `5e3becb`, 즉 `fc7b3e2`
  "New branch for NPU/LPU support" 이전 — 이 문서가 있는 `xPU`의 최근 timing/energy fidelity
  수정 커밋들을 전혀 포함하지 않은 상태)에서 수행됐고, 아직 `xPU`에 병합되지 않았다.** 그
  +14.2%라는 수치는 지금의 `xPU` NPUsim으로 다시 낸 값이 아니므로, calibration/holdout
  데이터로 그대로 인용해서는 안 된다 — 재현하려면 같은 trace/config를 현재 `xPU` HEAD로 다시
  실행해야 한다.
- 더 큰 convolution(C=128, C=192, 즉 atomic-C 그룹이 2개 이상)에서는 CRC mismatch(전체 출력
  0)가 재현됐고 원인은 아직 못 찾았다 — R2/R3 리스크가 실제로 발현된 사례이며, C-tail/K-tail
  calibration(6.2절)을 시작하기 전에 반드시 먼저 풀어야 한다.
- `core_cycles`(CACC 완료 시점) 추출, internal CBUF/CACC traffic 계측, sensitivity suite(8.2절),
  energy(Phase 5)는 전부 미착수 상태 그대로다.

즉 Phase 0은 "완료"는 아니지만 4.3절 완료 조건 5개 중 처음 3~4개(RTL 생성, Verilator build,
공식 trace 실행, output comparison PASS)는 이미 최소 한 번 통과된 적이 있다 — "동일 trace
반복 실행 시 bit-exact 재현"만 아직 확인되지 않았다. 이 문서의 나머지 절(5절 이후)이 서술하는
계획 자체는 이 결과로 인해 달라지지 않는다: n=1, 다른 base에서 낸 데이터 포인트 하나는 여전히
holdout 증거가 아니다.

### 4.3 완료 조건

- `nv_small` RTL 생성 성공
- Verilator testbench build 성공
- 공식 convolution trace 실행 성공
- RTL output memory comparison PASS
- 동일 trace 반복 실행 시 cycle과 traffic이 bit-exact하게 재현

## 5. Phase 1 — RTL golden collector

### 5.1 Cycle collector

`dla_core_clk`의 positive edge를 세는 64-bit counter를 testbench에 추가한다. 한 layer에 대해
두 cycle window를 수집한다.

| Metric | 시작 | 종료 | 용도 |
| --- | --- | --- | --- |
| `core_cycles` | CDMA/CSC launch | CACC completion | convolution core 원인 분석 |
| `end_to_end_cycles` | layer engine의 마지막 `D_OP_ENABLE` accept | SDP write-back completion interrupt | NPUsim 최종 비교 |

CSB register programming과 trace의 `load_mem`/`dump_mem`은 execution cycle에서 제외한다.
NVDLA trace player의 직접 memory load/dump는 simulation time을 사용하지 않으므로 NVDLA DMA
traffic에도 포함하지 않는다.

RTL이 제공하는 performance counter도 함께 저장한다.

- CDMA data/weight read stall
- CDMA data/weight read latency
- SDP read/write DMA stall
- 지원되는 PDP/CDP/RUBIK counter는 해당 engine을 검증할 때만 사용

근거:

- [NVDLA Programming Guide — performance debug](https://nvdla.org/hw/v1/ias/programming_guide.html)
- [NVDLA trace-player testbench](https://github.com/nvdla/doc/blob/master/doc/hw/v1/integration_guide.rst)

### 5.2 External traffic collector

DBBIF slave-agent 또는 DUT memory interface에서 **실제로 handshake가 완료된 transaction만** 센다.

- read request 수와 요청 byte
- returned read beat 수와 byte
- write request 수
- accepted write-data beat 수
- write strobe/byte-enable가 활성화한 실제 byte
- outstanding request와 stall cycle

메모리 주소 영역을 사전에 분리하여 다음 stream별로 집계한다.

- input activation
- weight
- output
- intermediate/metadata

Testbench가 수행하는 초기 `load_mem`, 결과 확인용 `read_mem` 및 `dump_mem`은 모두 제외한다.

### 5.3 Internal traffic collector

초기에는 CBUF SRAM port의 read/write enable과 byte width를 직접 계수한다. 이후 필요하면 다음
boundary를 추가한다.

- CSC → CMAC activation delivery
- CSC → CMAC weight delivery
- CMAC → CACC partial sum/result delivery
- CACC access
- CACC → SDP output delivery

NPUsim과 NVDLA의 component 대응은 다음과 같이 정의한다.

| NPUsim component/boundary | NVDLA RTL boundary |
| --- | --- |
| DRAM | DBBIF external memory traffic |
| Global buffer input/weight | CBUF data/weight banks |
| PE-array/NoC | CSC → CMAC delivery |
| Accumulator/output storage | CACC |
| Output write-back | SDP WDMA |
| Multi-chip/NoP | 대상 없음 |

NVDLA의 CBUF는 input/weight 저장소이고 CACC는 별도 accumulator이므로, NPUsim의 GLB traffic
하나를 CBUF와 직접 동일시해서는 안 된다.

## 6. Phase 2 — Workload suite

### 6.1 Bring-up set

1. 공식 `nv_small` sanity convolution
2. 공식 INT8 direct-convolution trace 한 개
3. 8×8 atomic geometry에 정확히 맞는 최소 convolution

먼저 공식 trace를 재사용한다. 새로운 shape가 필요할 때는 NVDLA verification suite의 trace
generator를 사용하고, 생성 명령과 결과 hash를 보존한다.

근거: [NVDLA Verification Suite User Guide](https://nvdla.org/hw/v2/verif_guide.html)

### 6.2 Calibration set

| 분류 | 검증 목적 |
| --- | --- |
| 단일 atomic tile | 기본 MAC/pipeline latency |
| C가 8의 배수 | Atomic-C 정상 이용률 |
| C tail | channel padding과 idle MAC |
| K tail | kernel padding과 idle MAC |
| 1×1 convolution | filter reuse와 weight-heavy 경로 |
| 3×3 convolution | spatial/filter reuse |
| stride=2 | input window 이동과 traffic |
| explicit padding | input halo 및 zero-padding |
| CBUF-fit | resident reuse |
| CBUF-overflow | refill/refetch와 memory stall |

### 6.3 Holdout set

- C tail과 K tail이 동시에 존재하는 layer
- stride와 padding이 동시에 존재하는 layer
- activation-heavy layer
- weight-heavy layer
- FC-equivalent 1×1 layer
- calibration에 사용하지 않은 AlexNet/ResNet representative convolution
- unseen memory latency/bandwidth 조합

### 6.4 Energy 입력 패턴

Switching activity의 데이터 의존성을 보기 위해 각 energy workload를 최소 다음 패턴으로 실행한다.

- all-zero
- 고정 seed random 3개 이상
- 가능한 경우 representative activation/weight distribution

Cycle과 dense traffic은 seed에 불변이어야 한다. Energy는 seed별 결과와 평균, 표준편차를 모두
기록한다.

## 7. Phase 3 — Traffic fidelity

Traffic을 먼저 해결한다. Cycle의 memory stall과 energy event 수가 traffic에 종속되므로 external
traffic이 맞지 않는 상태에서 cycle/energy constant를 조정하지 않는다.

### 7.1 Metric과 gate

| Metric | 권장 통과 기준 |
| --- | ---: |
| Input logical bytes | RTL과 exact |
| Weight logical bytes | RTL과 exact |
| Output logical bytes | RTL과 exact |
| DBBIF physical read/write bytes | exact; 허용 시 마지막 partial beat 이내 |
| DBBIF accepted read/write beat 수 | exact |
| CBUF read/write bytes | holdout MAPE ≤ 5%, max ≤ 10% |
| Stream/source conservation | 모든 layer PASS |

Logical volume과 physical traffic을 별도 열로 저장한다. Atomic alignment, tail padding 또는 bus
beat padding으로 인해 physical byte가 logical byte보다 클 수 있다.

### 7.2 예상 NPUsim 보완 항목 (2026-08-21 코드 감사로 재평가)

작성 시점의 목록은 "NVDLA 고유 구조라서 NPUsim에 없을 것"이라는 추정으로 9개 항목을 균등하게
나열했다. 실제 `components/`·`scheduler/`·`utils/` 코드를 감사한 결과 그 추정은 절반 정도만
맞았다 — 나머지는 이미 architecture-neutral 메커니즘으로 존재하고, NVDLA config/mapping이
없어서 아직 안 써봤을 뿐이다. 항목별로 구분한다.

**A. 이미 일반 메커니즘으로 커버됨 — 코드가 아니라 NVDLA config/mapping이 필요**

- **Atomic-C/K tail padding** — `scheduler/npu.cc`의 "spatial-padding guard"가 mapping factor
  곱이 실제 layer 차원을 넘으면 그 초과분을 compute로 과금한다. Eyeriss silicon 검증
  (`validation/phase3`)이 이미 이 경로로 55→7×8=56 padding을 재현했다. NVDLA의 Atomic-C=8/
  Atomic-K=8을 PE-array-level 매핑 factor로 선언하면 그대로 적용된다.
- **64-bit memory atomic alignment** — `utils/datatype.cc`의 `storage_transactions()`가
  임의의 `transaction_bits`에 대해 `ceil_div`로 atomic 정렬을 계산하고, `components/dram.cc`는
  `bitwidth`를 config에서 직접 읽는다. `bitwidth=64`로 설정하면 nv_small의 실제 primary
  memory interface 폭과 일치한다. (실측 확인: 아래 4.2-후속절의 첫 RTL 데이터 포인트에서
  weight 8B→128B, output 1B→8B 물리 byte 팽창이 정확히 이 메커니즘으로 재현됐다.)
- **partial weight reuse** — `mapping_table_t::datatype_repetitions()`가 (K,C,R,S) loop
  factor만으로 weight의 off-chip 재fetch 배수를 스칼라로 계산하므로, B/P/Q만 도는 revisit은
  추가 DRAM traffic 없이 재사용된다. `validation/traffic`의 T2(=1.00, "fetch-once")로 8개
  workload에서 검증됨.
- **CACC partial-sum residency/spill** — `mapping_table_t::reduction_tiled_above_array()` +
  `pe_array_t::set_psum_retention_scope()`가 reduction 차원이 array 위에서 tiling되는지 여부로
  psum 보존/spill을 결정하고, `components/pe.cc`가 이를 accumulator spill 이벤트로 과금한다.
  둘 다 architecture-neutral.
- **SDP output write-back traffic** — `global_buffer_t::account_output_writeback_link()` →
  `multi_chip_t::account_output_writeback_to_dram()` → `multi_chip_t::flush_output_writeback()`
  로 이어지는 전용 write-back 경로가 이미 존재하고, layer 마지막 resident output tile을 놓치는
  버그(DR6)도 이미 고쳐져 있다.
- **request burst와 final partial beat** — `utils/interconnect_timing.h`의
  `transfer_packet_groups_t`/`transfer_packet_groups()`가 모든 transfer를 "마지막 한 개만
  나머지 크기의 tail packet"으로 이미 분해한다. 특정 아키텍처에 종속되지 않은 순수 bit-width
  산술이다.

**B. 실제로 코드 보완이 필요함**

- **CBUF bank allocation과 data/weight partition** — `separate_buffer_t`는 input/weight/output
  byte 단위 partition(`capacity_per_type`)은 갖고 있지만, `num_banks`는 "conflict-free 병렬성
  divisor"로만 쓰이고 bank 단위 용량이나 layer별 재할당 개념이 없다 (`global_buffer.cc`
  주석: "Bank conflicts, arbitration, and read/write port contention are NOT modeled"). 실제
  NVDLA CBUF(32 bank)는 firmware가 layer마다 weight/activation bank 배분을 바꿀 수 있는데,
  NPUsim의 `capacity_per_type`는 실행 전체에 걸쳐 고정된 정적 config 값이다. 이건 genuine gap
  이지만, 구현 우선순위는 아래 8.4절 개정본의 "실제 구현" 항목들보다 낮다 — 단일 layer 단위의
  Phase 1~3 holdout 실험(6.1~6.3절)에는 영향이 없고, 여러 layer가 서로 다른 bank 분할을
  요구하는 네트워크 단위 실험에서만 드러난다.
- **convolution stripe 및 channel reuse (CBUF-stripe/atomic-C 스트리밍 세부 단위)** —
  input halo union(`mapping_table_t::input_halo_reuse()`)과 datatype repetition 억제는 GLB
  tile/pass 단위에서는 이미 정확하지만, CBUF 내부에서 한 stripe를 CMAC이 아직 소비하는 동안
  다음 stripe를 DMA로 겹쳐 채우는 sub-tile 단위 겹침은 명시적으로 모델링돼 있지 않다. 다만 이
  겹침은 이미 일반적으로 존재하는 double-buffer 기반 파이프라인 오버랩(8.4절 항목 1 참고)이
  tile을 stripe 크기로 잡으면 대체로 흡수한다 — traffic *총량*에는 영향이 없고 그 총량이
  compute와 얼마나 겹치는지(즉 cycle)에만 영향을 준다는 뜻이므로, 이 항목은 7절(traffic)보다
  8절(cycle)의 문제로 재분류하는 것이 맞다.

**C. 목록에서 재분류 — 보완 항목이 아니라 매핑 작성 책임**

- **CBUF capacity에 따른 refill/refetch** — 처음에는 "누락된 메커니즘"으로 분류했으나, 코드
  감사 후 재검토한 결과 이는 NPUsim의 기존 설계와 일치한다: `separate_buffer_t::check_tile_size()`
  /`shared_buffer_t::check_tile_size()`는 매핑이 지정한 tile이 CBUF partition에 안 맞으면
  `exit(1)`로 즉시 실패한다. 이것은 버그가 아니라 **NVDLA 컴파일러의 실제 동작과 같은 계약**이다
  — NVDLA 컴파일러도 하나의 layer descriptor가 CBUF에 안 맞으면 하드웨어가 실행 중에 알아서
  나누는 게 아니라 컴파일 시점에 여러 개의 작은 descriptor로 미리 쪼갠다. 즉 "CBUF보다 큰
  weight/input을 어떻게 나눠 realize할지"는 이미 존재하는 tile-반복(datatype_repetitions) +
  용량 assert 메커니즘으로 표현되고, 필요한 것은 새 코드가 아니라 CBUF 128KiB에 맞게 tile
  factor를 고른 NVDLA `.map`이다. (Eyeriss/Gemmini 검증도 정확히 같은 방식으로 매핑을
  작성했다.) 6.2절의 "CBUF-overflow" calibration 항목은 이 계약이 실제로 지켜지는지 — 즉 CBUF에
  맞지 않는 매핑을 의도적으로 시도했을 때 NPUsim이 실패하고 RTL/컴파일러도 그 layer를 애초에
  받아들이지 않는지 — 를 확인하는 것으로 범위를 좁힌다.

## 8. Phase 4 — Cycle fidelity

### 8.1 Baseline memory condition

실험값은 반드시 plusarg로 명시하고 golden에 저장한다. 초기 예시는 다음과 같다.

- memory bandwidth utilization: 100%
- read latency: 64 cycles
- write response latency: 30 cycles
- max outstanding request: 128

이는 예시 baseline이며 실제 testbench가 해당 값을 지원하는지 Phase 0에서 확인한 뒤 고정한다.
NVDLA testbench는 bandwidth, read/write latency 및 outstanding transaction 수를 runtime option으로
조절할 수 있다.

근거: [NVDLA Integration Guide — memory model options](https://github.com/nvdla/doc/blob/master/doc/hw/v1/integration_guide.rst)

### 8.2 Sensitivity suite

- bandwidth utilization: 25%, 50%, 100%
- read latency: 32, 64, 128, 180 cycles
- CBUF-fit과 CBUF-overflow workload
- compute-bound와 memory-bound workload

### 8.3 Metric과 gate

Layer별 signed error는 다음과 같이 계산한다.

```text
error[%] = 100 × (NPUsim cycles - RTL cycles) / RTL cycles
```

Holdout 기준 권장 gate:

- MAPE ≤ 10%
- max absolute error ≤ 15%
- signed mean error 절댓값 ≤ 5%
- 모든 sensitivity test에서 병목 변화 방향 일치

Calibration 결과와 holdout 결과를 분리해 보고한다. Tail, CBUF miss, memory stall 분류별 오차도
따로 출력한다.

### 8.4 예상 NPUsim 보완 항목 (2026-08-21 코드 감사로 재평가, 1건 구현 완료)

7.2절과 같은 기준으로 재분류한다.

**A. 이미 일반 메커니즘으로 커버됨 — NVDLA config만 필요**

- **CDMA load와 CMAC compute overlap** — 모든 계층 경계(DRAM/multi-chip/GLB/PE-array/PE)가
  `double_buffer` flag를 갖고, `scheduler/stats.cc`가 이를 4개 경계의 overlap 배열로 엮어
  `utils/interconnect_timing.h`의 `pipeline_timeline_cycles()`(N-deep staging buffer의 per-tile
  타임라인, 아래 B항목이 바로 이 함수를 확장한 결과다)에 넘긴다. 완전히 architecture-neutral.
- **CBUF refill stall** — 같은 5-stage per-tile 타임라인이 GLB(=CBUF 대응) 단계의 refill-vs-drain
  stall을 N-deep buffer와 함께 모델링하고, `Back-pressure stall`로 리포트에 출력한다.
  `validation/buffering`(LB-A~E), `validation/bottleneck`(BN1-BN6, DRAM-bound 재분류 포함)로
  검증됨.
- **layer setup 고정 비용** — `components/pe_array.h`의 `u_layer_setup_cycle`이 layer당 1회,
  config로 지정하는 고정 비용을 이미 charge한다 (Gemmini RTL 보정값 2270 cycle 예시 존재).
  다만 "interrupt" 왕복 지연을 setup 비용과 별도로 분리하지는 않는다 — 하나의 scalar
  knob으로 합쳐 충당해야 한다.

**B. 이번 검토에서 구현 완료 — memory request concurrency와 backpressure**

착수 전 감사에서 `pipeline_timeline_cycles()`가 이미 임의의 N-deep staging buffer를 지원함에도
(`boundary_depths[b]`가 "N: an N-deep queue"라고 문서화돼 있음), DRAM↔multi-chip 경계의 실제
depth는 `multi_chip_t::double_buffer` bool 하나에서 `? 2 : 1`로만 도출되고 있어 8.1절이 요구하는
"max outstanding request: 128" 축을 **구조적으로 표현할 수 없었다** — `components/dram.cc`도
"request queueing... not modeled"라고 명시하고 있었다.

`components/multi_chip.h`/`.cc`에 `max_outstanding_requests` config 키(기본값 0 = 미설정)를,
`scheduler/stats.{h,cc}`에 해당 경계 depth override를 추가해 이 축을 열었다. 0(미설정)이면 기존
`double_buffer` 기반 1/2 depth를 정확히 재현하므로 기존 config는 전부 이전과 bit-identical하다
(`validation/bottleneck`, `validation/traffic`, `validation/buffering` 전부 수정 전후 숫자까지
동일 — 회귀 없음 확인). 새 knob 자체는 `validation/memory_concurrency/check.py`(MC1-MC5)로
검증했다: depth를 1→128로 올리면 DRAM back-pressure stall이 정확히 그만큼 줄고 critical path도
그 감소분만큼만 줄어드는 것을 end-to-end로 확인했다. 상세 기록:
[implementation/timing/memory_concurrency_2026-08-21_nvdla_prep.md](../implementation/timing/memory_concurrency_2026-08-21_nvdla_prep.md).

이 fix는 여전히 **tile 단위** depth cap이다 — RTL의 실제 per-request address 기반 큐잉/bank
conflict는 여전히 모델링하지 않으며, `dram_t::describe_timing_limits()`의 해당 disclaimer는
그대로 유효하다. 8.2절 sensitivity suite를 *실행 가능*하게 만들었을 뿐, 실제 NVDLA RTL 대비
이 축의 정합성(8.3절 gate)은 Phase 0/1 golden data 없이는 아직 검증되지 않았다.

**병합 상태:** 이 fix는 branch `worktree-nvdla-fidelity-review`(commit `9f1897a`, `xPU` HEAD
`8c9102d`에서 분기)에 커밋돼 있고, 아직 `xPU`에 병합되지 않았다. 즉 지금 이 문서를 보고 있는
`xPU` 작업 디렉터리의 `components/multi_chip.h`·`.cc`, `scheduler/stats.h`·`.cc`에는 아직
`max_outstanding_requests`가 없다 — 코드를 실제로 쓰려면 먼저 이 branch를 review/merge해야
한다.

**C. 실제로 코드 보완이 필요함 (미착수)**

- **CSC/CMAC/CACC pipeline fill/drain** — 현재 `utils/interconnect_timing.cc`의
  `systolic_pipeline_cost()`는 weight-residency 경계당 fill+drain을 **하나로 뭉쳐** 계산한다
  (Gemmini RTL 보정, `weight_fold_fill_cycle=14`). NVDLA의 CSC(주소 생성)→CMAC(MAC array)→
  CACC(accumulator)는 서로 다른 fill/drain latency를 갖는 독립된 3단 파이프라인이다. nv_small의
  8×8 atomic array처럼 작은 array에서는 이 뭉침이 실측 오차(첫 데이터 포인트의 +14.2%, 4.2-후속절
  참고)의 상당 부분을 설명할 수 있는 유력한 후보다 — 다만 n=1 데이터로는 원인을 단정할 수 없다.
- **Atomic-C/K tail utilization** — `utils/pe_lane.h`/`.cc`의 lane 모델은 1차원이다
  (`active_scalar_lanes`가 `mac_width`의 배수가 아닐 때만 tail을 계산). NVDLA는 Atomic-C(reduction
  방향)와 Atomic-K(병렬 복제 방향)가 독립적인 2차원 8×8 array다. 대략적인 tail padding은 이미
  7.2절 A항목의 spatial-padding guard(mapping factor 단위, 이미 2차원)가 흡수하지만, PE 내부
  lane 단위의 세밀한 2차원 반영은 없다 — 필요성은 C 항목 원인 규명(첫 항목) 이후 실측으로
  재평가한다.
- **SDP pass-through/write-back latency** — write-back *링크* 자체(GLB→multi-chip→DRAM)는 이미
  전용 경로로 모델링돼 있다(7.2절 A항목의 SDP output write-back traffic 참고). 그러나 SDP의
  *기능* 단계(bias/activation/pooling 자체의 pipeline latency)는 timing path에 전혀 없다 — 다만
  이 계획서 2.2절이 이미 pooling/activation/normalization을 초기 검증 범위에서 명시적으로
  제외했으므로, 이 항목은 nv_small holdout 1차 통과 전까지는 낮은 우선순위다.

## 9. Phase 5 — Energy fidelity

### 9.1 두 단계의 energy 검증

Energy는 다음 두 단계로 나눈다.

1. **Activity/relative energy fidelity**
   - RTL toggle 및 event 변화 방향
   - workload와 memory setting에 따른 component energy 순위
   - 아직 absolute hardware energy라고 주장하지 않음
2. **Absolute implementation energy fidelity**
   - 합성 netlist와 power library
   - measurement-window switching activity
   - SRAM macro/CACTI energy
   - DRAM device energy
   - leakage 포함

Plain RTL 또는 Verilator toggle count만으로는 absolute energy golden을 만들 수 없다.

### 9.2 Activity dump

- reset과 CSB programming 구간 제외
- layer launch부터 completion까지만 VCD/FST/SAIF 생성
- 정상 clock gating 활성화
- full-run VCD 대신 measurement window와 관심 hierarchy만 dump
- 입력 seed별 activity file 및 hash 저장

NVDLA는 내부 clock gating을 사용하므로 이를 강제로 해제한 activity를 기본 golden으로 사용하지
않는다.

근거: [NVDLA Integration Guide — clock and power control](https://github.com/nvdla/doc/blob/master/doc/hw/v1/integration_guide.rst)

### 9.3 Synthesis와 power analysis

필요한 입력은 다음과 같다.

- standard-cell Liberty
- target voltage/temperature/process corner
- synthesis constraint와 target frequency
- mapped gate-level netlist
- activity annotation
- clock-tree/clock-power 처리 정책
- SRAM macro power model 또는 독립 CACTI reference

NVDLA repository는 Design Compiler용 reference synthesis flow와 partition별 constraint를 제공한다.
그러나 RAM은 behavioral model이므로 실제 SRAM macro/wrapper와 power model을 integrator가 제공해야
한다.

근거: [NVDLA Integration Guide — synthesis and memories](https://github.com/nvdla/doc/blob/master/doc/hw/v1/integration_guide.rst)

### 9.4 Energy boundary

두 종류의 total을 별도 보고한다.

```text
NVDLA core energy
  = logic dynamic
  + clock dynamic
  + CBUF/CACC SRAM dynamic
  + core/SRAM leakage over execution window

System-inclusive modeled energy
  = NVDLA core energy
  + external DRAM device energy
  + explicitly modeled off-chip link/PHY energy
```

DRAM energy는 RTL DBBIF transaction trace를 동일한 DRAMsim3 device configuration에 replay하여
구한다. Off-chip PHY/link가 reference에 없으면 0으로 간주하지 말고 `UNPRICED/OUT OF SCOPE`로
표시한다.

### 9.5 NPUsim unit-cost 설정 원칙

- MAC, logic, SRAM 단가는 독립 characterization에서 얻는다.
- layer total을 맞추기 위해 component unit cost를 역산하지 않는다.
- 동일 component cost를 calibration과 holdout에 그대로 사용한다.
- technology node, voltage, frequency scaling을 묵시적으로 섞지 않는다.
- power tool이 보고한 hierarchy와 NPUsim component boundary를 mapping table로 남긴다.

### 9.6 Metric과 gate

Holdout 기준 초기 gate는 다음과 같다.

| Metric | 권장 통과 기준 |
| --- | ---: |
| Total dynamic energy MAPE | ≤ 15% |
| Total energy max absolute error | ≤ 25% |
| Logic/CBUF/DRAM component MAPE | ≤ 20% |
| Workload별 energy rank | 100% 일치 |
| Energy provenance와 measurement boundary | 모든 결과에 존재 |

Reference tool/library uncertainty가 gate보다 큰 경우에는 단일 pass/fail 대신 uncertainty interval과
함께 보고한다.

### 9.7 기존 NPUsim energy 인프라와의 관계 (2026-08-21 코드 감사)

9절은 "무엇을 검증할지"를 새로 정의하지만, 이를 뒷받침할 provenance/schema 인프라는 상당 부분
이미 존재한다 — 처음부터 만들 필요는 없고, 아래 세 가지 차이만 메우면 된다.

**이미 있는 것 (재사용 대상)**

- `utils/energy_units.h`/`.cc`의 `energy_unit_t`(UNSPECIFIED/PICOJOULE/NORMALIZED) + 자유
  텍스트 `energy_reference` provenance 문자열. `is_absolute()`는 `pJ` 단위와 비어 있지 않은
  provenance 둘 다 있어야 absolute로 인정한다 — bare `energy_unit=pJ`는 거부된다.
- `energy_cost_schema_t`가 컴포넌트별(MAC/PE/PE_ARRAY/GLOBAL_BUFFER/MULTI_CHIP/DRAM) 상태를
  NOT_MODELED/PARTIAL/MODELED_ZERO/CALIBRATED로 분류하고, `validation/energy_schema/check.py`
  (ES1-ES7)로 검증된다.
- `validation/unpriced/check.py`(UP1-UP5) + `unpriced_event_t`가 "이벤트는 발생했는데 pricing
  key가 선언 안 됨"을 이름까지 지목해 flag하고 wattage 출력을 막는다. 명시적 `key = 0`("modeled
  zero")과 key 부재("unpriced")를 구분한다 — 이것이 9.4절이 요구하는 "UNPRICED/OUT OF SCOPE
  명시" 원칙의 컴포넌트-cost-key 단위 구현체다.
- `validation/knobs/check.py`의 KN2(energy scale invariance: 모든 unit cost ×10 → 모든 리포트
  energy ×10, cycle은 불변)와 KN3(zero collapse)는 9.5절이 요구하는 "component unit cost를
  역산하지 않는다" 원칙을, 총량에 맞춰 몰래 상수를 끼워 넣으면 깨지는 구조적 불변식으로 이미
  강제하고 있다.

**메워야 할 차이**

- **9.4절의 core-vs-system-inclusive 이분법이 지금은 없다.** 현재는 리포트가 "core datapath만
  포함, DRAM background/refresh·clock network·controller/DMA·PHY는 제외"라는 단일 경계 하나만
  갖고, DRAM device access energy는 이미 그 안에 포함돼 있다(외부 시스템 항목으로 별도 표시되지
  않음). 9.4절 형식을 그대로 쓰려면 "NVDLA core energy" 행과 "system-inclusive" 행을 별도로
  뽑는 리포트 항목이 새로 필요하다 — provenance 스키마 확장이지 재작성이 아니다.
- **DRAM device별 energy 차등화가 안 돼 있다.** `validation/phase4/check.py`의 P4-8이 이미
  발견해 놓은 사실: 현재 배포된 6개 config 전부가 선언된 DRAM 종류(HBM2 vs DDR3 등)와 무관하게
  동일한 8.0 pJ/byte를 쓴다 — 실측 스펙 상 대역폭에 따라 최대 4.5배 차이가 나는데도. 9.4절의
  "DRAM energy는 RTL DBBIF trace를 동일 DRAMsim3 device configuration에 replay"를 실제로
  하려면 이 calibration 부채부터 갚아야 한다. 코드 로직 문제가 아니라 config 데이터 문제이므로
  이번 소스코드 수정 범위에는 포함하지 않았다 — 별도로 처리할 항목으로 남긴다.
- **calibration/holdout 분리를 강제하는 장치가 없다.** 9.5절의 "동일 component cost를
  calibration과 holdout에 그대로 사용" 원칙은 현재 CACTI/DRAMsim3라는 calibration 소스가 하나뿐
  이라 아직 테스트할 대상 자체가 없다 (KN2/KN3는 이 원칙이 깨지면 잡아내지만, calibration↔holdout
  분리 자체를 강제하지는 않는다).
- **9.6절이 요구하는 절대 energy MAPE gate는 아직 없다.** 이는 계획서가 처음 지적한 대로 새로
  만들어야 하는 게 맞다 — 다만 그 gate는 위 provenance/unpriced 인프라 위에 얹으면 되고, 처음부터
  설계할 필요는 없다.

## 10. 자동화 산출물

권장 repository 구조는 다음과 같다.

```text
validation/nvdla/
├── README.md
├── PROVENANCE.md
├── workloads.csv
├── run_rtl.sh
├── extract_rtl_metrics.py
├── generate_npusim_inputs.py
├── compare_cycle.py
├── compare_traffic.py
├── compare_energy.py
├── golden_rtl_cycle.csv
├── golden_rtl_traffic.csv
└── golden_rtl_energy.csv

configs/accelerators/
└── nvdla_small.cfg

models/mappings/nvdla_small/
└── *.map
```

각 golden row에는 최소 다음 필드를 포함한다.

- RTL commit/config hash/patch hash
- trace와 input file hash
- complete layer shape
- datatype/storage/compression mode
- memory latency/bandwidth/outstanding 설정
- cycle window와 RTL cycle
- stream별 logical/physical traffic
- CBUF/CACC traffic
- input seed
- synthesis/power tool와 library/PVT
- frequency와 activity window
- component별 dynamic/static/total energy

`compare_*.py`는 golden 파일이 없거나 provenance가 불완전하면 skip이 아니라 실패해야 한다.

**2026-08-21 참고:** `validation/nvdla/`는 아직 `xPU`에 없다 (4.2-후속절 참고 — 다른 base에서
한 번 만들어졌으나 미병합). 대신 이번 개정에서 `validation/memory_concurrency/check.py`가
새로 생겼다 — NVDLA용은 아니지만, "config 한 줄 차이로 실행 → 리포트 파싱 → 방향성/등식
assert" 패턴이 `compare_*.py`가 따라야 할 것과 같은 스타일이므로 참고 템플릿으로 쓸 수 있다.

## 11. 권장 해결 순서

1. NVDLA source와 `nv_small` commit 고정
2. RTL/Verilator build 및 공식 convolution trace PASS
3. cycle/DBBIF collector 구현
4. 단일 trace의 golden cycle/traffic 생성
5. `nvdla_small.cfg`와 대응 mapping 작성
6. external traffic exactness 해결
7. CBUF/CACC traffic 계측 및 계층 정합성 해결
8. baseline cycle 오차 해결
9. bandwidth/latency sensitivity 해결
10. activity dump 및 relative energy 검증
11. synthesis/SRAM power reference 구축
12. absolute component/total energy 검증
13. calibration과 독립 holdout 평가
14. regression gate와 최종 assessment 작성

Traffic을 cycle보다 먼저 고치는 이유는 memory stall cycle이 traffic/request sequence에 종속되기
때문이다. Energy는 traffic event와 실행 시간이 모두 확정된 뒤 진행한다.

## 12. Phase별 종료 조건

| Phase | 종료 조건 |
| --- | --- |
| Phase 0 | 공식 trace가 반복 가능하게 PASS하고 모든 dependency/provenance가 기록됨 |
| Phase 1 | cycle, DBBIF, CBUF collector가 수기 micro-test 산술과 일치 |
| Phase 2 | calibration/holdout workload 및 input hash 고정 |
| Phase 3 | external traffic gate 통과, internal CBUF traffic 오차 범위 충족 |
| Phase 4 | holdout cycle MAPE/max/bias gate 통과 |
| Phase 5A | data/memory 변화에 대한 RTL·NPUsim energy 방향과 rank 일치 |
| Phase 5B | 절대 component/total energy gate와 provenance 계약 통과 |

## 13. 주요 리스크

### R1. Toolchain 노후화

NVDLA build infrastructure가 오래된 tool version을 전제로 한다. 현재 toolchain patch가 필요할 수
있으며, patch가 RTL 의미론을 바꾸지 않는지 공식 trace로 확인해야 한다.

### R2. Trace 생성 난이도

임의 layer의 valid register/weight/data trace 생성이 build보다 어려울 수 있다. 따라서 공식
`nv_small` trace를 먼저 사용하고 trace generator를 재현한 뒤 workload를 확장한다.

### R3. NPUsim abstraction gap

CBUF bank scheduling, CSC/CMAC/CACC pipeline, SDP write-back 및 memory concurrency가 현재 NPUsim
component model에 직접 표현되지 않을 수 있다. 한 개의 fitted setup constant로 여러 구조적 차이를
숨기지 말고 traffic/stall/pipeline counter를 분리해 구현한다.

### R4. Energy reference의 불완전성

RTL에는 공정과 SRAM macro power가 없다. Standard-cell library 또는 SRAM reference가 없으면
absolute energy 단계는 완료할 수 없으며, 그 상태의 결과는 normalized/relative로만 표시한다.

### R5. Waveform 용량

전체 hierarchy/full-run VCD는 디스크를 빠르게 소모할 수 있다. measurement window와 관심 hierarchy로
제한하고 가능한 경우 FST/SAIF를 사용한다.

## 14. 최종 주장 형식

모든 holdout gate를 통과했을 때 다음 수준으로 결과를 기술한다.

> NPUsim was validated against cycle-exact Open NVDLA `nv_small` RTL for the documented INT8 dense
> direct-convolution holdout suite. It achieved the reported cycle error bounds and exact external
> traffic accounting under matched memory conditions. Energy agreement is relative to the stated
> synthesized NVDLA implementation, activity traces, SRAM characterization, and DRAM model; it is
> not a claim of measured NVDLA silicon energy fidelity.

`nv_small` convolution 이외의 architecture, precision, compression, Winograd, post-processing 및
silicon operating point로 이 주장을 일반화하지 않는다.

## 15. 바로 시작할 최소 작업

첫 구현 단위는 Phase 0과 Phase 1의 일부로 제한한다.

1. `ext/nvdla` dependency를 reproducible commit으로 준비
2. `nv_small` Verilator build
3. 공식 INT8 convolution trace 한 개 PASS
4. layer 실행 cycle과 DBBIF input/weight/output byte 추출
5. 한 줄짜리 `golden_rtl_cycle.csv`와 `golden_rtl_traffic.csv` 생성

이 결과가 확보된 뒤에야 NPUsim의 `nvdla_small.cfg`와 mapping을 확정한다. 먼저 config를 추정해
만들면 RTL scheduling과 다른 가정을 고정할 위험이 있다.

**2026-08-21 참고:** 이 다섯 단계는 branch `worktree-nvdla-fidelity`에서 한 번 실제로 시도됐다
(4.2-후속절 참고) — `ext/nvdla` 준비, Verilator build, 공식 trace CRC PASS, cycle/traffic golden
추출까지는 성공했지만, base가 `xPU`의 최근 timing/energy 수정 커밋들을 포함하지 않고 아직
병합도 안 됐다. 그대로 재사용하지 말고, RTL golden(3번까지의 산출물)은 참고하되 **NPUsim
쪽(4~5번)은 현재 `xPU` HEAD로 다시 뽑아야** 한다. 또한 그 시도에서 atomic-C 그룹이 2개 이상인
convolution(C=128/192)의 CRC mismatch가 재현됐다 — 이 원인을 먼저 규명하지 않으면 6.2절의
C-tail/K-tail calibration set을 만들 수 없으므로, 이번 최소 작업의 실질적인 최우선 순위는 이
버그 규명이다.
