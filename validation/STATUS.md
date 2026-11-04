# NPUsim 정합성 현황 (2026-11-03)

NPUsim을 MOSAIC에 붙인 뒤 진행한 정합성 작업의 현재 상태. 각 항목의 자세한 근거와 수치는 링크한 디렉터리의
README에 있다. 브랜치: npusim `mosaic`, MOSAIC `npusim-ddr3` -> `mosaic`.

## 한눈에

| 대상 | 기준 (ground truth) | 현재 오차 | 상태 |
|---|---|---|---|
| Gemmini GEMM, L2 warm (6점) | GemminiRocketConfig RTL (Verilator) | MAPE 4.4%, 최대 7.9% | 완료 (`phase2/`) |
| Gemmini GEMM, L2 cold + DDR3 (7점) | 같은 RTL + SimDRAM/DRAMSim2 DDR3-1333 | MAPE 9.0% (메모리 바운드 -1.4%, -8.6%) | 16x512x512 -28.9% 남음 |
| TPU v4 PaLM decode | Pope et al. 2023 공개 지연 | NPUsim 41-60% 짧음 | attention 미포함 시점의 결과, 재평가 필요 |
| TPU v5e decode GEMM (M 1..128) | Colab v5e-1 실측, weight 스트리밍 | MAPE 약 5% (-0.3..-8%) | 완료 (`tpu_v5e/`) |
| TPU v5e compute-bound GEMM | 같은 실측 | -1..+12% | 완료 |
| TPU v5e decode attention | 같은 실측, KV 스트리밍 | B >= 8 MAPE 8.3%, B = 1 17.6% | **이번에 모델 연결** |
| TPU v5e 작은/비정렬 GEMM | 같은 실측 | -10..-36% | 원인 후보만, 미반영 |

## 1. Gemmini RTL (`validation/phase2/`)

- warm L2 golden에 대해 compute-schedule MAPE 4.4% (V2/V3: RTL 타일링, GLB 이중 버퍼, fold fill 재보정).
- cold L2 + DRAM 재측정 벤치(`gemm_sweep_cold*.c`)와 MOSAIC Ramulator2 DDR3 비교(`compare_cold.py`, 설정 C):
  MAPE 26.1% -> 9.0%. 반영한 구조 모델 (모두 `[dram]` 키, 기본값 꺼짐, `gemmini_rtl_mem.cfg`에서 켬):
  - DMA 요청 크기, 동시 요청 수, 경로 지연 (`mosaic_request_bytes`, `mosaic_max_inflight`, `mosaic_path_latency_ns`, 쓰기용 키)
  - L2 write-allocate RFO와 dirty write-back (`mosaic_line_bytes`, `mosaic_cache_kb`)
  - L2 쓰기 포트 처리율 (`mosaic_write_request_ns`, DRAM과 별개 자원으로 critical path에서 max 결합)
  - 출력 행의 L2 set 충돌 (`mosaic_cache_ways`)
- 근거 측정: RTL 성능 카운터(겹침), DMA 지연, 메모리 경로 마이크로벤치(`mem_path.c`, `rtl_mem_path.txt`).
- 남은 것: 16x512x512 -28.9%. golden DRAMSim2의 `scheme2` 매핑이 512 B 간격 B 읽기를 한 bank에 몰고
  `TOTAL_ROW_ACCESSES=4`가 row를 닫는다. 해결에는 (a) MOSAIC Ramulator2에 scheme2 같은 RoCoBaRaCh 매퍼와
  `ddr3_1333_simdram.yaml` 정렬, (b) NPUsim `mosaic_stream`의 행 단위 strided 주소 스트림이 둘 다 필요.

## 2. MOSAIC DRAM 연동

- `NPUSIM_MOSAIC_DRAM=<yaml>`이면 dense 타일 스트림의 DRAM 시간을 MOSAIC Ramulator2(CosimFrontend)로 계산.
- `configs/architecture/dram/ddr3_1333_simdram.yaml`: Chipyard SimDRAM의 DDR3-1333 (MOSAIC 쪽).
- MOSAIC 변환기 버그 수정: `ClosedRowPolicy` -> Ramulator 2.1 `ClosedCAP(cap)`.
- MOSAIC의 기존 2.0 형식 HBM2 설정은 HBM2_R20 모델로 가는데, 16 B 트랜잭션에 nBL = 4 tCK라 이론의 1/4만 낸다
  (`hbm2_cosim_test` 23%). 다른 엔진이 쓰므로 고치지 않았다.
- TPU v5e HBM2e 설정 `configs/architecture/dram/hbm2e_tpu_v5e.{py,yaml}` (Ramulator 2.1 원본 HBM2, 공개 사양만):
  순차 읽기 698 GB/s (85.2%, 전체 뱅크 refresh), refresh 없이 748 (91.3%), 실측 732-757 (89-92%). TPU는 refresh를
  대부분 숨긴다. NPUsim decode GEMM에 연결하면 MAPE 5.4% (819 GB/s) -> 4.4% (`tpu_v5e/compare_hbm.py`).

## 3. TPU v5e 단일 칩 (`validation/tpu_v5e/`, 종합 보고서 `tpu_v5e/REPORT.md`)

측정: Colab v5e-1, `tpu_v5e1_stream_bench.py` 3회 (결과 `results/2026-11-03_colab_v5e1_followup*`, 원본 zip
포함). 병렬 세션의 `measure_tpu.py` warm/rot 결과(브랜치 `tpu-measure`, zip은 `results/`)와 대조했다.

확정된 사실:
- 실효 HBM 읽기 732-757 GB/s (사양 819의 89-92%). weight 스트리밍(690-750)과 KV 읽기(733-797)도 같은 속도.
- decode GEMM은 HBM 바운드, NPUsim -0.3..-8% (819 GB/s 가정만큼 짧다).
- decode 레이어 = GEMM만 레이어 + attention (12개 (B, L) 모두 2.4% 이내). attention = KV 바이트 / 약 745 GB/s.
- int8 weight는 decode 크기(M <= 32)에서 bf16의 절반 시간.
- 측정 재현성: 실행 간 중앙값 0.3-1.3% (작은 GEMM 2.8-6.3%).

attention 모델 (이번 작업): NPUsim에 이미 있던 attention 소비자 모드(`[kvcache] attention = 1`: K/V 공급을
DRAM 모델 위 전용 스트림으로, QK/softmax/AV를 선언된 처리율로)에 `batch`(시퀀스 수, 기본 1)를 추가했다. 읽기,
cache append, MAC, softmax, 점유량이 batch배가 된다. `compare_attention.py`가 측정과 같은 Llama-3-8B 구조로
15개 (B, L)를 돌린다: B >= 8은 -6..-8% (B 8 / L 512만 -17.8%), B = 1은 -12..-30%. 처리율은 구조값
(MXU 4 x 128^2 MAC/cycle, softmax 0.005 cycle/score)이고, 단계가 KV 읽기에 묶여 있어 결과에 거의 영향이 없다.
일관된 -7%는 GEMM과 같은 HBM 효율 차이다.

측정 방법에서 얻은 교훈 (다음 TPU 측정에 그대로 적용):
- 같은 입력 버퍼를 반복하면 XLA cross-program prefetch가 copy를 생략해 weight가 칩 위에 남는다 (warm: 3-4배 빠름).
- 프로그램 하나에 GEMM 하나면 HBM -> VMEM copy와 연산이 직렬이 된다 (rot: HBM 바운드 GEMM이 25-50% 느림).
  실제 decode처럼 한 프로그램 안에서 여러 weight를 스트리밍해야 한다.
- 루프 의존성을 결과의 원소 하나로만 이으면 XLA가 slice를 연산 안으로 밀어 원소 하나만 계산할 수 있다.
- 큰 버퍼를 쌓아 동적 slice로 꺼내면 HBM 복사본이 생긴다. KV를 (B, L, kv, d)로 두면 einsum용 relayout copy가
  생긴다 (attention 3.6배로 측정됨). 서로 다른 버퍼를 반복 본문에 펼치고 KV는 (B, kv, L, d)로 둔다.
- 루프 결과와 무관한 읽기는 루프 밖으로 빠진다. 읽기마다 carry를 섞는다. 컴파일된 HLO로 확인할 것.
- Colab에서 JAX가 TPU를 못 찾으면 경고만 내고 CPU로 돈다. 벤치는 TPU가 없으면 첫 셀에서 멈춘다.

## 4. 남은 과제 (우선순위 순)

1. HBM 효율: HBM2e 설정이 refresh 없이는 실측과 맞고(91%), 전체 뱅크 refresh로는 6%p 낮다(85%). TPU처럼
   refresh를 스트림에 숨기는 refresh 관리자(또는 근거 있는 REFSB 스케줄)가 남은 일. compute-bound GEMM이 MOSAIC
   경로에서 4-12% 빨라지는 원인도 확인할 것.
2. PaLM/TPU v4 비교에 attention 모델을 넣어 41-60% 공백 재평가 (tpu_v4 + `[kvcache]` attention).
3. 작은/비정렬 GEMM: GEMM당 고정 비용(약 1-2 us)과 4 MXU 열 분할 단위 (1024 -> 1100: NPUsim 1.6배, 실측 1.27배).
4. decode 레이어의 원소 단위 연산(잔차, gate/up 분리; GEMM만 레이어가 GEMM 4개 합보다 1-12% 김).
5. int8 weight의 tpu_v5e 설정 비교 (실측은 있음).
6. 8칩 collective: Kaggle v5e-8 (`kaggle_tpu_bench.ipynb`) 실행 대기.
7. Gemmini 16x512x512: 1절의 (a), (b).

## 5. 재현

```
# Gemmini cold + MOSAIC DDR3 (설정 A/B/C)
LD_LIBRARY_PATH=<npusim libs> python3 validation/phase2/compare_cold.py
# TPU v5e attention: NPUsim vs 실측
LD_LIBRARY_PATH=<npusim libs> python3 validation/tpu_v5e/compare_attention.py
# TPU v5e GEMM: NPUsim vs 실측 (compare_v5e.py, 브랜치 tpu-measure)
python3 validation/tpu_v5e/colab_summary.py validation/tpu_v5e/results/2026-11-03_colab_v5e1_followup3
# TPU 측정: Colab v5e-1에서 tpu_v5e1_stream_bench.ipynb Run All (단계 선택: %env TPU_BENCH_SECTIONS=...)
```

브랜치 관계: 병렬 세션의 `tpu-measure` (measure_tpu.py, compare_v5e.py, tpu_v5e.cfg)는 아직 `mosaic`에 합치지 않았다.
`tpu_v5e.cfg`는 같은 내용으로 `mosaic`에도 있다. 합칠 때 `validation/tpu_v5e/README.md`는 손으로 합쳐야 한다.
