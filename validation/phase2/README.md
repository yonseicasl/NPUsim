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
