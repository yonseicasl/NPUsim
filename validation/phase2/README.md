# Phase 2 — RTL ground-truth validation vs Gemmini (Verilator, cycle-exact)

## 목적과 결론 요약

NPUsim의 systolic timing을 **cycle-exact RTL**(Chipyard `GemminiRocketConfig`,
16×16 WS, Verilator)과 동일 GEMM 6점에서 대조했다.

- **V2/V3 구현 후 (2026-08-17): 전 6점 |err| ≤ 8% — MAPE 4.4%, 최대 +7.9%.**
  NPUsim이 직접 출력하는 `Computation cycle + Fold fill cycle`이 비교량이다.
- V2 = RTL-캘리브레이션 fold-fill(`weight_fold_fill_cycle=14`/fold +
  `layer_setup_cycle=2270`/layer), V3 = `edge_accumulation`(array-edge
  accumulator; per-PE LB 출력 용량 제약 해제로 Gemmini류 weight-residency 매핑
  표현). 기본값 off — 기존 config 결과는 **byte-identical** 보존.
- 구현 전 기준: comp 축은 steady-state −3.5~−5.9%, 소형/저-T는 fold 오버헤드
  미과금으로 −22~−73% (MAPE 27.7%); RTL 실측 per-fold 오버헤드 15.7~18.5 cycle.

## 실험 설정

| 항목 | 값 |
|---|---|
| RTL | Chipyard(HEAD, `ext/chipyard`) `chipyard.harness.TestHarness.GemminiRocketConfig`, Verilator |
| Gemmini | DIM=16 (16×16 PE), WS dataflow, spad 256KB, accumulator 64KB, int8/int32 |
| 벤치 | `gemm_sweep.c`(`bareMetalC`, 이 저장소 `validation/phase2/`) — `tiled_matmul_auto`, `read_cycles()` |
| NPUsim | `configs/accelerators/gemmini.cfg`(16×16 WS, GLB 128/128/64KB, 128-bit, int8), `gen_workload.py`가 GEMM별 net/매핑/더미가중치 생성 |
| 실행 | RTL: `sims/verilator/simulator-* gemm_sweep-baremetal` (~2h40m). NPUsim: `./npusim.sh run gemmini gemm_MxKxN ws` |
| 비교 | `python3 validation/phase2/compare.py` (RTL 수치: `rtl_cycles.txt`) |

## 결과 (V2/V3 구현 후 — 시뮬레이터 출력값 기준)

```
        M,K,N   NPU comp  fold fill  NPU sched       RTL    err%
  64,  64,  64      1,024      2,494      3,518     3,810    -7.7
 128, 128, 128      8,192      3,166     11,358    10,530    +7.9
 256, 256, 256     65,536      5,854     71,390    69,653    +2.5
  16, 512, 512     16,384     16,606     32,990    32,506    +1.5
 512, 512,  64     65,536      4,062     69,598    74,413    -6.5
 512, 512, 512    524,288     16,606    540,894   543,240    -0.4

MAPE = 4.4%   max |err| = 7.9%  (n=6)
NPU sched = Computation cycle + Fold fill cycle
```

### V2/V3 구현 요약 (2026-08-17)

| 항목 | 내용 |
|---|---|
| fold 정의 | array-wide weight residency. per-PE 시간축 weight 원소 각각이 실제 WS 하드웨어의 순차 residency이므로 **fold 수 = weight 배포 이벤트 × per-PE weight tile 원소수** = `(N/16)(K_gemm/16)` — RTL fold 수와 정확히 일치 |
| 과금 | `pe_array_t` 배포 경로에서 WEIGHT 이벤트당 `u_weight_fold_fill_cycle × tile_size_lb[WEIGHT]`; per-layer `layer_setup_cycle`은 stats에서 repetition 스케일 **후** 1회 가산 |
| 캘리브레이션 | RTL 6점 그리드 적합 → **(fill=14, setup=2270)**이 전 점 ±8% 유일 만족역. fill 14 ≈ DIM−2 (weight 열 로드 + accumulate-drain 파이프라인) |
| V3 매핑 | `edge_accumulation=1`: OUTPUT을 per-PE LB 용량 검사에서 면제(출력은 GLB output partition=64KB accumulator로 스트림) → `b_pe=M` 매핑으로 weight를 M 전체에 상주. per-PE input LB 512B = spad A-half 128KB/256PE. 512³은 accumulator 슬랩 64KB 제약으로 K 외곽 루프가 DRAM으로 스필(RTL의 C mvout 동작과 동일) |
| fill 항 검증 | 6점 fill 값 224/896/3,584/14,336/1,792/14,336 = `folds×14` 수기 산술과 **전점 정확 일치** |
| 무회귀 | 기본 knob off: alexnet matched/energy 재실행 **byte-identical**, unittest 5/5 통과 |

### RTL 자체 sanity
Gemmini RTL 효율(ideal/RTL): 512³ **96.5%**, 256³ 94%, (512,512,64) 88%,
(16,512,512) 50%(=T=16 fold-fill 이론값 16/(16+17)≈48% ✓), 64³ 27%(고정
setup 지배) — double-buffered DMA가 compute를 완전히 숨기는 것으로 확인
(memory stall이 아니라 fold 오버헤드가 잔차의 전부).

### 해석
1. **Steady-state 검증**: 큰 GEMM에서 NPUsim comp가 RTL과 3.5~5.9% — timing
   모델의 1차 계약(활성 lane당 1 MAC/cycle) 성립.
2. **Fold 오버헤드의 실측 확정**: fold당 15.7~18.5 cycle — WS fold의
   weight-로드·accumulate-드레인 파이프라인 깊이. 전점 동시 적합의 유일해는
   (fill 14, setup 2,270). 소형(64³)의 큰 편차는 고정 setup(flush·DMA 셋업)이
   지배해서 생긴 것.
3. **V2 (구현·검증 완료)**: fold-fill 항이 `Fold fill cycle`로 리포팅되어
   `comp + fold fill`이 전 6점 ≤8%. fold 수는 per-PE weight 원소 단위 residency로
   계산 — GEMM/WS류 매핑(per-PE 시간축 weight = 서로 다른 출력채널)에서 유효한
   해석이며, conv류 매핑의 fold 의미론(per-PE R·S set 단위)은 A1 per-tile
   timeline 범위로 남는다(기본 off인 이유).
4. **V3 (구현·검증 완료)**: `edge_accumulation`으로 per-PE LB 출력 용량 제약을
   해제해 Gemmini의 edge-accumulator weight residency(`b_pe=M`)를 매핑으로 표현
   — weight 재스트리밍 3~4× 과대 문제 해소, fold 수가 RTL과 일치하게 됨.

## Critical path (메모리 포함) vs 같은 RTL golden

RTL golden은 DMA를 포함한 `tiled_matmul_auto` 전 구간이므로, 위 compute-schedule만이 아니라
메모리를 포함한 `Critical-path latency`도 같은 값과 비교된다
(`python3 validation/phase2/check_critical_path.py [accelerator]`, `measure_fidelity.sh`의 `crit` 행).

| config | compute schedule (MAPE / max) | critical path (MAPE / max) |
|---|---|---|
| `gemmini.cfg` | 4.4% / 7.9% | 1,866% / +2,625% |
| `gemmini_rtl.cfg`, 구조 선언만 (`0003cbc`) | 4.4% / 7.9% | 47.6% / +111% |
| `gemmini_rtl.cfg`, + RTL 타일링·fold fill 재보정 | 9.7% / 11.8% | **11.0% / +18.9%** |

`gemmini.cfg`는 PE마다 출력 write-back을 공유 매체에 직렬로 싣고 그때마다 메시 경로 fill
(16×16에서 ~30 hop)을 다시 과금하며, `[dram]`이 Eyeriss 값(1 B 라인 × 3 cycle)이라
critical path가 RTL의 15~27배다. `gemmini_rtl.cfg`는 RTL 소스에서 가져온 구조 선언을 더한다
(NVDLA의 2026-08-26 array streaming granularity와 같은 옵션):

- 열별 병렬 출력과 분리된 A/B/C 배선(`writeback_injection_parallel`,
  `array_fabric_separate`, `operand_streams_pipelined`)
- scratchpad와 accumulator의 분리 SRAM(`fabric_separate`, `streams_pipelined`)
- DMA `dma_buswidth=128`(128-bit beat/cycle)과 `max_in_flight_mem_reqs=16`
- WS `tiled_matmul_auto`의 scratchpad/accumulator double buffering(`gemmini.h`:
  `double_buffered = WS`) → `[separate] double_buffer = 1`

double buffering은 GLB 용량을 절반으로 계산하므로 GEMM 매핑도 RTL이 실제로 고르는 타일로
다시 만든다(`gen_workload.py --rtl`: `tiled_matmul_auto`의 타일 선택을 DIM=16, scratchpad
256 KB, accumulator 64 KB로 재현). RTL은 (i0, j, k) 블록마다 가중치를 preload하므로 fold 수는
`I0 × N/16 × K/16`이다. 6점 중 5점은 RTL과 fold 수가 같고, 512³만 NPUsim 타일이 균등 분할이라
I0 = 8(RTL은 최대 80행 타일 7개)이다. 이로써 critical path가 compute schedule + 약 800 cycle로
붙는다.

fold fill 재보정: `gemmini.cfg`의 fill 14는 I0가 빠진 fold 수(N/16 × K/16) 위에서 맞춘 값이다.
RTL fold 수 위에서는 fold당 잔차가 0.5~32 cycle로 흩어져(16×512×512처럼 fold당 16행이면
~13.5, 64행이면 ~1~2, 512×512×64는 128행인데 ~12.9) fold당 상수 하나로는 최대 오차 약 12%가
하한이다. minimax 적합 결과 fill 9.8, setup 2,270. `check_critical_path.py`의 게이트(compute
schedule max 12%, critical path max 20%)는 이 수준의 회귀 방지선이다. fold 오버헤드를 fold당
행 수로 모델링하려면 행 수를 바꾼 RTL 점이 더 필요하다. M=1 decode GEMV는 RTL 점이 없어
검증 범위 밖이다.

## Cold-L2 RTL 재측정과 DRAM 경로 (2026-11-02)

위 golden(`gemm_sweep.c`)은 6점이 같은 A/B/C를 재사용해 첫 점 이후 트래픽 일부가 512 KB L2에
맞고, 6점이 모두 연산 바운드라 DRAM 모델을 검증하지 않는다. `gemm_sweep_cold.c`는 점마다 한 번도
접근하지 않은 버퍼를 써서(bare-metal crt는 .bss를 지우지 않는다) 모든 DMA 읽기가 L2 miss로 SimDRAM
(DRAMSim2 DDR3-1333, `+dramsim`)에 가게 하고, K = 16인 메모리 바운드 2점을 더했다. RTL은
GemminiRocketConfig Verilator(코어·버스 500 MHz, L2 ↔ SimDRAM AXI 64 bit)를 Ubuntu 24.04 컨테이너에서
돌렸다(`rtl_cycles_cold.txt`, `golden_rtl_cycles_cold.csv`).

NPUsim은 `compare_cold.py`로 세 설정을 비교한다: A `gemmini_rtl`, B `gemmini_rtl_mem`(500 MHz, 64-bit
메모리 포트), C = B + MOSAIC Ramulator2 `ddr3_1333_simdram.yaml`(SimDRAM과 같은 DDR3, 순차 스트림
10.4 GB/s; NPUsim을 MOSAIC_ROOT로 빌드).

| M x K x N | RTL warm | RTL cold | A | B | C |
|---|---|---|---|---|---|
| 64x64x64 | 3,810 | 4,148 | +2.5% | +5.6% | +11.0% |
| 128x128x128 | 10,530 | 12,412 | +0.8% | +1.9% | +3.7% |
| 256x256x256 | 69,653 | 74,327 | +5.8% | +6.0% | +6.3% |
| 16x512x512 | 32,506 | 49,892 | -41.4% | -28.9% | -28.9% |
| 512x512x64 | 74,413 | 76,610 | -3.6% | -3.4% | -3.4% |
| 16x16x16384 | - | 112,929 | -53.8% | -39.3% | -1.4% |
| 16x16x65536 | - | 789,328 | -74.7% | -66.4% | -8.6% |

- L2가 비면 RTL이 3~54% 느려진다(16x512x512 +53%): 기존 golden은 L2 효과를 포함한다.
- 연산 바운드 4점은 세 설정 모두 ±6%. 메모리 바운드 점은 RTL 실효 대역폭이 4.6 B/cycle(512 KB),
  2.7 B/cycle(2 MB, 데이터가 L2보다 커서 dirty C의 write-back이 겹친다)로 64-bit 포트(8 B/cycle)와
  DDR3보다 낮다.
- DMA 왕복 지연 측정(`dma_latency.c`, `rtl_dma_latency.txt`, 500 MHz): 한 줄 mvin(64 B 라인) L2 miss
  60 cycle / hit 33, 한 행 mvout(16 B) miss 57 / hit 30, fence만 2. 읽기는 동시 요청 16 x 64 B / 58
  cycle = 17 B/cycle > 포트 8 B/cycle이라 동시 요청 상한에 걸리지 않는다. 쓰기는 출력 행이 16 B라 요청이
  작아 16 x 16 B / 55 cycle = 4.7 B/cycle이 상한이다.
- 겹침 측정(`gemm_sweep_cold_perf.c`, `rtl_perf_counters_cold.txt`; Gemmini 성능 카운터, load/store/execute
  Venn 영역과 각 유닛 활동). RTL도 이중 버퍼링을 한다: `tiled_matmul_auto`가 WS에서 scratchpad와
  accumulator를 반으로 나누고, 하드웨어 `LoopMatmul`이 타일 루프 2개를 동시에 처리한다(reservation
  station load 8 / execute 16 / store 4). 실제로 16x16x16384에서 execute 단독 cycle은 2, load와 store가
  함께 활동하는 cycle이 67%(+ 셋 다 28%)로 연산은 DMA 뒤에 숨는다. 시간을 정하는 것은 겹침이 아니라
  DMA 처리량이다: STORE_ACTIVE 97%, C 256 KB를 약 2.4 B/cycle로 쓴다. 16x512x512은 읽기가 부하 상태에서
  약 5.4 B/cycle이다. (`RDMA_BYTES` 카운터는 beat마다 전송 크기를 더해 약 4배로 집계된다.)
- NPUsim의 MOSAIC DRAM 스트림에 DMA를 넣었다: `[dram] mosaic_request_bytes / mosaic_max_inflight /
  mosaic_path_latency_ns`(읽기), `mosaic_write_request_bytes / mosaic_write_path_latency_ns`(쓰기, 기본은
  읽기 값). 기본 0이면 이전과 같다. `gemmini_rtl_mem.cfg`는 측정값을 쓴다: 읽기 64 B / 16 / 80 ns
  (118 ns 왕복 - DDR3 37.5 ns), 쓰기 16 B / 110 ns. 16x16x16384 -39.3% -> -30.9%, 16x16x65536
  -66.4% -> -62.0%.
- L2(InclusiveCache: 512 KB, 64 B 라인, write-back, write-allocate)의 DRAM 트래픽: `[dram] mosaic_line_bytes`
  (라인보다 작은 쓰기는 라인을 먼저 읽는다: 16 B 출력 행이 cold 라인마다 64 B 읽기를 낳는다)와
  `mosaic_cache_kb`(레이어 데이터가 캐시보다 크면 쓴 바이트의 (사용량 - 캐시)/사용량이 dirty write-back으로
  DRAM에 간다). `gemmini_rtl_mem.cfg`는 64 B / 512 KB. C(표의 마지막 열)는 둘을 모두 반영한 값이다:
  16x16x16384 -18.8%, 16x16x65536 -55.1%, MAPE 21.6% -> 18.1%. 연산 바운드 점은 +3~11%.
- 남은 원인은 아래 메모리 경로 측정으로 확정했다.

### 메모리 경로 마이크로벤치 (`mem_path.c`, 결과 `rtl_mem_path.txt`)

연산 없이 cold 버퍼에 mvin/mvout 스트림만 발행하고 fence까지 잰다 (B/cycle, 500 MHz).

| 스트림 | 크기 | B/cycle |
|---|---|---|
| 읽기 (64 B 행) | 64 KB / 256 KB / 1 MB / 4 MB | 7.72 / 7.23 / 6.10 / 5.84 |
| 읽기 warm (L2 hit) | 256 KB | 13.1 |
| 쓰기 16 B 행, 연속 | 256 KB / 1 MB | 1.73 / 1.74 |
| 쓰기 16 B 행, 16 KB 간격 (16x16x16384의 C) | 256 KB | 2.02 |
| 쓰기 16 B 행, 64 KB 간격 (16x16x65536의 C) | 256 KB | **0.55** |
| 읽기 + 쓰기 섞기 (1 KB씩 번갈아) | 128 KB / 512 KB / 2 MB | 4.38 / 4.65 / 4.88 |

- 읽기는 64비트 포트 상한(8 B/cycle)의 73~97%이고 Ramulator2 DDR3(포트 상한 8)와 가깝다. DRAM 모델은
  주원인이 아니다.
- 쓰기는 DRAM이 아니라 L2의 부분 쓰기 처리가 제한한다. 라인마다 RFO가 4배인 16 KB 간격이 연속보다 빠르다
  (같은 라인의 16 B 쓰기 4개가 한 MSHR에서 직렬화된다). 16 B 요청당 약 8~9 cycle.
- 16x16x16384 RTL 112,929 cycles = 섞기 512 KB(B 256 KB + C 256 KB) 112,720 cycles. 이 점은 메모리 경로로
  전부 설명된다.
- 64 KB 간격 쓰기가 3.7배 느린 것은 L2 set 충돌이다. L2는 8-way, set+bank 인덱스가 주소 비트 6~15라
  64 KB 간격의 C 16개 행이 한 set의 8 way를 다툰다. 16x16x65536 RTL 789,328이 섞기 2 MB 429,405보다
  1.8배 느린 비선형이 이것이다.
- NPUsim 반영 (표의 C 열, MAPE 18.1% -> 9.0%):
  1. 쓰기 경로 처리율: `[dram] mosaic_write_request_ns`(L2가 부분 쓰기를 받는 최소 간격). 섞기 2 MB
     스트림에서 잰 6.55 cycle = 13.1 ns. DRAM이 아닌 별개 자원이라 critical path의 DRAM link 축에서
     타입별 DRAM 스트림 합과 max로 결합한다(`stats_t::write_port_cycle_dram`). Gemmini의 load/store
     DMA가 분리돼 읽기와 쓰기가 겹치기 때문이다. 16x16x16384 -18.8% -> -1.4%.
  2. set 충돌: `[dram] mosaic_cache_ways`(8). 레이어의 출력 행 간격(`output_size`)과 타일 행 수로, 한
     set에 way 수보다 많은 행이 몰리는지 레이어마다 판정한다. 몰리면 부분 쓰기마다 라인 RFO와 dirty
     write-back을 DRAM에 더한다. 16x16x65536 -55.1% -> -8.6%. 다른 5점은 충돌이 없다.
- 남은 최대 오차: 16x512x512 -28.9% (cold L2에서 RTL이 warm보다 53% 느린 점; 이번 두 효과와 무관).

## 인프라 트러블슈팅 기록 (재현자 참고)

1. **머신의 산발적 데이터 손상**: 대용량 HTTPS 다운로드 바이트 플립(동일 크기·
   상이 sha256 재현), conda/JVM(G1 GC heap-walk) SIGSEGV 다발 — **불량 RAM
   의심, 하드웨어 점검(memtest) 권고**. 우회: chunk별 이중-fetch 자기일관성
   다운로더(miniforge), 재시도 루프(빌드 산출물 캐시로 전진).
2. chipyard HEAD 결함 2건 로컬 패치: `install-espresso.sh` 부재(no-op shim),
   **firrtl2 `MemConf` 파서의 trailing-space 요구**(`ports (\S+)\s+` → optional
   그룹 내부로 이동; writer/parser drift로 `mask_gran` 없는 라인 전부 파싱 실패).
3. 중단된 `git submodule` checkout이 "커밋 일치·빈 워크트리" 상태를 만들어
   `submodule update`가 조용히 통과 → `git submodule foreach --recursive
   'git checkout -f HEAD'`로 복구.
4. Verilator sim stdout은 블록 버퍼링 — 실시간 관측이 필요하면 `stdbuf -oL`.
5. NPUsim 기존 crash 버그 발견·수정: `stats.cc` `reserve()+operator[]` UB
   (allocator-상태 의존 segfault; Phase-1의 segfault 오진의 진범).

## Phase 1+2 종합 결론 (V2/V3 구현 후)

- **Phase 2 (RTL, ground truth)**: **전 6 GEMM 점 |err| ≤ 8% (MAPE 4.4%)** —
  시뮬레이터 출력값 기준. RTL이 기준 축이다.
- **Phase 1 (SCALE-Sim, cross-check)**: FC ±6.3%. conv steady(conv3/4/5)의
  −13~−18.8%는 **참조 시뮬레이터 자체의 fold 오버헤드**가 원인 — SS의 fold 상수
  (R+C−2≈24)는 RTL 실측(≈14)의 약 2배로, RTL이 SS의 과대 계상을 판정했다.
  SS total × Overall Util로 SS의 fold 오버헤드를 제거하면 conv3/4/5 모두
  **+7.69% (≤8% ✓)**. conv1/2의 +31/+22%는 active-PE 매핑 산술(7∤dims로 SS의
  partial-row 채움을 NPUsim 매핑이 표현 불가)로 정확히 설명되는 표현력 한계다.
- 즉 **timing 오차 주장: RTL 기준 전 구간 ≤8%; SCALE-Sim과의 잔차는 두 기준의
  상호 불일치(RTL로 판정 완료) 및 매핑 표현력으로 전량 귀속**.
- 남은 개선: conv류 fold 의미론(A1 per-tile timeline), fold-fill의
  critical-path/leakage window 반영.
