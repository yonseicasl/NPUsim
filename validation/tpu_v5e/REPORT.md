# NPUsim vs Cloud TPU v5e 정합성 비교 보고서 (단일 칩, 2026-11-03/04)

## 요약

NPUsim의 `tpu_v5e.cfg`(공개 사양만으로 만든 구조 설정, 측정값에 맞춘 파라미터 없음)를 Colab TPU v5e-1에서 잰
실측과 비교했다.

| 항목 | 실측 | NPUsim 오차 | 비고 |
|---|---|---|---|
| decode GEMM (Llama-3-8B projection, M 1..128, 16개) | HBM 바운드, weight 690-750 GB/s | **MAPE 5.4%** (-2.3..-8.3%) | MOSAIC HBM2e 사용 시 MAPE 4.4% (+0.8..+13%) |
| compute-bound GEMM (M 512, 2048, 1024^3 이상) | 187-192 TFLOPS (사양 197) | -6..+12% | |
| decode attention (B >= 8) | KV 바이트 / 약 745 GB/s | **MAPE 8.3%** (-6..-8%, 1개 -17.8%) | B = 1은 MAPE 17.6% |
| decode 레이어 | = GEMM만 레이어 + attention (2.4% 이내) | projection GEMM 4개가 GEMM만 레이어의 -4..-11% | 원소 단위 연산 미모델 |
| 작은/비정렬 GEMM (M 100..512, K=N 1000..2048) | | -35..+15% (평균 -0.5..-25%) | 고정 오버헤드, MXU 분할 |
| HBM 순차 읽기 | 732-757 GB/s (사양 819의 89-92%) | 분석 모델은 819 가정 | MOSAIC HBM2e 698 GB/s (85%) |

남은 차이의 주원인은 하나다: NPUsim의 분석 DRAM 모델이 HBM 사양 대역폭(819 GB/s)을 그대로 쓰는데 실측 실효
대역폭은 89-92%다. decode GEMM과 attention의 일관된 -5..-8%가 이것이다.

## 1. 비교 대상

**하드웨어**: Google Colab v5e-1 TPU (`TPU v5 lite`, 1칩), JAX 0.11.2. 측정 3회(2026-11-03 08:29, 12:27, 13:16 UTC).

**NPUsim 설정** (`configs/accelerators/tpu_v5e.cfg`, tpu_v4.cfg에서 MXU 수, 클록, HBM만 바꿈):
- TensorCore 1개 x MXU 4개(128 x 128 bf16), 1500 MHz (197 TFLOPS / (4 x 128^2 x 2))
- HBM2e 16 GiB, 819 GB/s (`[dram]`: 1024 B 라인당 1.875 cycle)
- GEMM 매핑: validation/palm의 방식(출력 열을 128열 weight strip으로 4개 MXU에 분배, 행은 512행 타일)
- attention: `[kvcache] attention = 1`(K/V를 DRAM 모델 위 전용 스트림으로 읽고 QK/softmax/AV 계산), 이번에
  `batch` 추가. 처리율은 구조값(MXU 4 x 128^2 MAC/cycle, softmax 0.005 cycle/score).
- 선택: `NPUSIM_MOSAIC_DRAM=hbm2e_tpu_v5e.yaml`이면 GEMM 타일 스트림의 DRAM 시간을 MOSAIC Ramulator2 HBM2e로 계산.

**워크로드**: Llama-3-8B decode (D 4096, 32 q-head, 8 kv-head, d 128, FFN 14336), bf16.

## 2. 측정 방법

`tpu_v5e1_stream_bench.py`(Colab 노트북): 측정마다 한 XLA 프로그램 안의 `fori_loop`을 n번, 2n번 돌려 반복당
시간 = (t(2n) - t(n)) / n, 5회 중앙값. 매 반복이 서로 다른 피연산자 사본을 읽어(총 512 MB 이상) 칩 위 재사용이
없다. 반복은 이전 결과 전체에 의존한다.

측정 방식이 결과를 크게 좌우했다. 같은 decode GEMM 16개에 대한 NPUsim 오차:

| 측정 방식 | 8x4096x4096 실측 | NPUsim MAPE (M <= 128) | 문제 |
|---|---|---|---|
| warm (같은 버퍼 반복, `measure_tpu.py`) | 11.7 us | 121.9% | XLA cross-program prefetch가 copy를 생략, weight가 칩 위에 남음 |
| rot (버퍼 순환, 프로그램당 GEMM 1개) | 56.4 us (copy 44.8 + 연산 11.6) | 16.9% | 프로그램에 연산이 하나라 HBM -> VMEM copy와 연산이 직렬 |
| **stream (한 프로그램 안 다중 weight)** | **47.4 us** | **5.4%** | 실제 decode와 같은 구조 |
| NPUsim | 43.4 us | | |

이 밖에 벤치를 고치며 확인한 함정(컴파일된 HLO로 확인): 결과의 원소 하나만 이어 쓰면 XLA가 slice를 연산 안으로
밀어 넣을 수 있다; 쌓아 둔 큰 버퍼를 동적 slice로 꺼내면 HBM 복사본이 생긴다; KV를 (B, L, kv, d)로 두면 einsum용
relayout copy가 생긴다(attention이 3.6배로 측정됨); 루프 결과와 무관한 읽기는 루프 밖으로 빠진다.

재현성(2차 대비 3차, 중앙값 / 최대 차이): stream GEMM 0.9% / 5.8%, decode 레이어 0.7% / 5.1%, attention 0.3% / 28.8%,
작은 GEMM 6.3% / 70.9%. 큰 최대값은 수 us짜리 경우(B 1 L 512 attention, M 100 K=N 1000-1024)다.

## 3. 결과

### 3.1 HBM 순차 읽기

| 크기 | 64 MB | 256 MB | 1 GB | 4 GB |
|---|---|---|---|---|
| GB/s (사양 819 대비) | 732 (89.3%) | 749 (91.4%) | 755 (92.1%) | 757 (92.4%) |

weight 스트리밍(690-750 GB/s), attention의 KV 읽기(733-797 GB/s)도 같은 속도다.

### 3.2 decode GEMM (weight 스트리밍, 3차 실행)

`compare_hbm.py` (`results/2026-11-03_colab_v5e1_followup3/compare_hbm.log`), us:

| M | K | N | 실측 | NPUsim 819 GB/s | 오차 | NPUsim + MOSAIC HBM2e | 오차 |
|---|---|---|---|---|---|---|---|
| 1 | 4096 | 6144 | 69.6 | 64.9 | -6.8% | 71.3 | +2.4% |
| 1 | 4096 | 4096 | 47.0 | 43.3 | -8.0% | 47.5 | +1.1% |
| 1 | 4096 | 28672 | 313.8 | 302.3 | -3.7% | 332.7 | +6.0% |
| 1 | 14336 | 4096 | 157.5 | 147.2 | -6.6% | 162.1 | +2.9% |
| 8 | 4096 | 6144 | 69.3 | 65.1 | -6.1% | 71.5 | +3.3% |
| 8 | 4096 | 4096 | 47.4 | 43.4 | -8.3% | 47.7 | +0.8% |
| 8 | 4096 | 28672 | 313.3 | 303.1 | -3.3% | 333.6 | +6.5% |
| 8 | 14336 | 4096 | 158.2 | 147.4 | -6.8% | 162.5 | +2.7% |
| 32 | 4096 | 6144 | 69.0 | 65.9 | -4.5% | 72.4 | +5.0% |
| 32 | 4096 | 4096 | 46.5 | 44.0 | -5.3% | 48.4 | +4.1% |
| 32 | 4096 | 28672 | 313.1 | 306.0 | -2.3% | 336.6 | +7.5% |
| 32 | 14336 | 4096 | 158.6 | 148.6 | -6.3% | 163.8 | +3.3% |
| 128 | 4096 | 6144 | 72.4 | 69.6 | -3.9% | 76.4 | +5.4% |
| 128 | 4096 | 4096 | 49.8 | 46.9 | -5.8% | 51.4 | +3.3% |
| 128 | 4096 | 28672 | 310.6 | 319.8 | +2.9% | 351.1 | +13.0% |
| 128 | 14336 | 4096 | 164.3 | 153.8 | -6.4% | 169.5 | +3.1% |
| | | | | **MAPE** | **5.4%** | | **4.4%** |

- 한 프로그램 안에서는 warm(같은 weight 반복)과 cold가 같다(비율 0.99-1.02). on-chip 재사용은 프로그램 경계를
  넘는 prefetch에서만 생긴다.
- int8 weight: bf16 대비 M 1-32에서 0.50-0.54배(HBM 바운드, 바이트가 절반), M 128에서 0.58-0.61배, M 512에서
  0.96-1.02배. dequant 비용은 decode 크기에서 보이지 않는다. (NPUsim int8 설정과의 비교는 아직)

### 3.3 compute-bound GEMM

| 경우 | NPUsim 오차 |
|---|---|
| M 512, Llama projection 4개 (stream) | -2.1, -5.6, +5.0, -2.6% (MOSAIC HBM2e: -4..-12%) |
| M 2048, 같은 4개 (`measure_tpu.py`, 연산 시간) | +1.0..+7.0% |
| 정사각 1024^3 .. 8192^3 (같음) | -0.8, +12.0, +6.8, -0.0% |

큰 GEMM은 187-192 TFLOPS(사양의 95-97%)에 이른다.

### 3.4 decode attention (KV 스트리밍)

`compare_attention.py` (`results/2026-11-03_colab_v5e1_followup3/compare_attention.log`), us:

| B | L | KV MB | 실측 | NPUsim | roofline (819) | 오차 |
|---|---|---|---|---|---|---|
| 1 | 512 | 2 | 3.1 | 2.7 | 2.6 | -13.9% |
| 1 | 2048 | 8 | 15.0 | 10.5 | 10.2 | -30.0% |
| 1 | 8192 | 32 | 48.8 | 41.9 | 41.0 | -14.3% |
| 1 | 32768 | 128 | 190.7 | 167.4 | 163.9 | -12.2% |
| 8 | 512 | 16 | 25.6 | 21.0 | 20.5 | -17.8% |
| 8 | 2048 | 64 | 89.0 | 83.8 | 81.9 | -5.9% |
| 8 | 8192 | 256 | 365.6 | 334.8 | 327.8 | -8.4% |
| 8 | 32768 | 1024 | 1440.7 | 1338.8 | 1311.0 | -7.1% |
| 32 | 512 | 64 | 91.3 | 83.9 | 81.9 | -8.2% |
| 32 | 2048 | 256 | 361.1 | 334.9 | 327.8 | -7.3% |
| 32 | 8192 | 1024 | 1442.5 | 1338.9 | 1311.0 | -7.2% |
| 32 | 32768 | 4096 | 5761.2 | 5354.9 | 5244.2 | -7.1% |
| 128 | 512 | 256 | 364.4 | 335.4 | 327.8 | -8.0% |
| 128 | 2048 | 1024 | 1449.7 | 1339.4 | 1311.0 | -7.6% |
| 128 | 8192 | 4096 | 5764.6 | 5355.4 | 5244.2 | -7.1% |

MAPE: B >= 8 8.3% (n=11), B = 1 17.6% (n=4). 실측은 B >= 8에서 KV roofline의 1.03-1.12배, 즉 attention은 실효
HBM 속도의 KV 읽기다. NPUsim의 연산 처리율은 결과에 거의 영향이 없다(KV 읽기가 지배).

### 3.5 decode 레이어 (Llama-3-8B 1 레이어, 3차 실행), us

| B | L | 레이어 | GEMM만 레이어 | attention | (GEMM만 + attn) / 레이어 | GEMM 4개 실측 합 | NPUsim GEMM 4개 | NPUsim / GEMM만 |
|---|---|---|---|---|---|---|---|---|
| 1 | 512 | 625.3 | 620.3 | 3.1 | -0.3% | 588.0 | 557.6 | -10.1% |
| 1 | 2048 | 643.3 | 620.7 | 15.0 | -1.2% | 588.0 | 557.6 | -10.2% |
| 1 | 8192 | 675.9 | 620.2 | 48.8 | -1.0% | 588.0 | 557.6 | -10.1% |
| 8 | 512 | 621.0 | 581.4 | 25.6 | -2.3% | 588.1 | 559.0 | -3.8% |
| 8 | 2048 | 676.8 | 595.2 | 89.0 | +1.1% | 588.1 | 559.0 | -6.1% |
| 8 | 8192 | 945.2 | 582.5 | 365.6 | +0.3% | 588.1 | 559.0 | -4.0% |
| 32 | 512 | 684.4 | 609.7 | 91.3 | +2.4% | 587.2 | 564.5 | -7.4% |
| 32 | 2048 | 966.4 | 605.5 | 361.1 | +0.0% | 587.2 | 564.5 | -6.8% |
| 32 | 8192 | 2039.8 | 603.0 | 1442.5 | +0.3% | 587.2 | 564.5 | -6.4% |
| 128 | 512 | 1031.2 | 654.1 | 364.4 | -1.2% | 597.1 | 590.0 | -9.8% |
| 128 | 2048 | 2113.2 | 657.7 | 1449.7 | -0.3% | 597.1 | 590.0 | -10.3% |
| 128 | 8192 | 6447.5 | 659.2 | 5764.6 | -0.4% | 597.1 | 590.0 | -10.5% |

- 레이어 = GEMM만 레이어 + attention이 12개 조건 모두 2.4% 안에서 성립한다. norm, softmax 외 원소 단위 연산은 작다.
- GEMM만 레이어는 GEMM 4개를 따로 잰 합보다 -1..+10% 길다(잔차 덧셈, gate/up 분리). NPUsim은 이 부분이 없다.
- 따라서 NPUsim으로 레이어 시간을 낼 때: projection GEMM 4개 + attention 단계 (+ 원소 단위 연산, 미모델).

### 3.6 작은/비정렬 GEMM (M 100, 128, 200, 256, 512 x K=N)

| K=N | 1000 | 1024 | 1100 | 1152 | 1280 | 2000 | 2048 |
|---|---|---|---|---|---|---|---|
| NPUsim 오차 평균 | -25.3% | -23.2% | -0.5% | -2.4% | -11.5% | -14.3% | -13.7% |
| 범위 | -35..-12% | -32..-13% | -12..+15% | -19..+11% | -20..-3% | -18..-9% | -18..-8% |

1000과 1024의 실측이 같아 128 단위 padding은 맞다. 원인 후보: GEMM당 고정 비용(약 1-2 us)이 NPUsim에 없고, 출력
열을 4 MXU에 128열 단위로 나누는 방식이 XLA보다 거칠다(1024 -> 1100에서 NPUsim 1.6배, 실측 1.27배). 수 us 규모라
실행 간 잡음(중앙값 6%)도 크다.

## 4. MOSAIC HBM2e

`configs/architecture/dram/hbm2e_tpu_v5e.{py,yaml}`(MOSAIC): HBM2e 2 stack = 16채널 x 2 pseudo channel, 3.2 Gbps,
16 GiB, Ramulator 2.1 원본 HBM2, JEDEC 타이밍, 전체 뱅크 refresh. 공개 사양만 썼다.

| 순차 읽기 | GB/s | 사양 대비 |
|---|---|---|
| HBM2e 설정 (전체 뱅크 refresh) | 698 | 85.2% |
| 같은 설정, refresh 없음 | 748 | 91.3% |
| 실측 | 732-757 | 89-92% |

refresh 이외의 손실은 실측과 맞고, TPU는 refresh를 대부분 숨긴다. Ramulator의 뱅크별 refresh 관리자는 연속
스트림에서 오히려 느리다(655, 575 GB/s). 이 설정으로 decode GEMM MAPE가 5.4% -> 4.4%가 되지만 부호가 바뀐다
(+1..+7%). compute-bound M 512는 MOSAIC 경로에서 -4..-12%로 더 빨라진다(원인 조사 필요).

## 5. 결론

1. decode의 두 주성분(projection GEMM, attention)은 NPUsim이 5-8% 안에서 맞춘다. 두 경우 모두 HBM 바운드이고,
   오차는 NPUsim이 사양 대역폭을 쓰는 데서 온다.
2. 레이어 시간은 GEMM + attention으로 분해되므로, NPUsim의 PaLM/decode 흐름에 attention 단계를 넣으면 레이어 단위
   비교가 가능하다.
3. TPU 측정에서는 측정 방식이 오차보다 크다(warm 122%, rot 17%, stream 5%). 한 프로그램 안에서 서로 다른 피연산자를
   스트리밍해야 실제 decode를 잰다.

## 6. 남은 과제

1. HBM 효율: refresh를 스트림에 숨기는 모델(또는 근거 있는 REFSB 스케줄)로 HBM2e를 89-92%에 맞추기. 분석 모델의
   819를 피팅으로 낮추지는 않는다.
2. compute-bound GEMM이 MOSAIC 경로에서 빨라지는 원인.
3. 작은/비정렬 GEMM: GEMM당 고정 비용, MXU 열 분할 단위.
4. 레이어의 원소 단위 연산(잔차, gate/up 분리).
5. int8 weight: NPUsim tpu_v5e int8 설정과 비교.
6. 8칩: tensor-parallel decode와 collective (Kaggle v5e-8, `kaggle_tpu_bench.ipynb`).
7. TPU v4 PaLM 비교(validation/palm, NPUsim 41-60% 짧음)에 attention 단계 넣어 재평가.

## 7. 재현

```
# 측정 (Colab v5e-1 TPU 런타임): tpu_v5e1_stream_bench.ipynb Run All -> zip을 results/<date>_.../에 풀기
# NPUsim 비교 (MOSAIC 빌드, LD_LIBRARY_PATH=<npusim libs>):
python3 validation/tpu_v5e/compare_hbm.py         # decode/compute GEMM: 819 GB/s vs MOSAIC HBM2e vs 실측
python3 validation/tpu_v5e/compare_attention.py   # attention 단계 vs 실측
python3 validation/tpu_v5e/colab_summary.py validation/tpu_v5e/results/2026-11-03_colab_v5e1_followup3
# HBM2e 순차 읽기: hbm_stream_probe.cc (빌드 방법은 파일 머리말)
```

| 파일 | 내용 |
|---|---|
| `results/2026-11-03_colab_v5e1_followup{,2,3}/` (+ `.zip`) | stream 벤치 1-3차 (1차는 attention/HBM 읽기 무효, 2차는 HBM 읽기 무효) |
| `results/2026-11-03_colab_v5e1_{warm,rot}.zip` | 병렬 세션의 `measure_tpu.py` 원본 (풀린 결과와 비교는 브랜치 `tpu-measure`) |
| `compare_hbm.py`, `compare_attention.py`, `colab_summary.py` | 비교 스크립트 |
| `tpu_v5e1_stream_bench.py`, `kaggle_tpu_bench.py` | 측정 벤치 (Colab 1칩, Kaggle 8칩) |
| `README.md` | 실행별 상세 기록 |
