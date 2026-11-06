# TPU 정합성 참고 자료

NPUsim을 Google TPU에 맞추기 위한 공개 자료 모음 (조사일 2026-11-02). arXiv 논문은
`bash refs/tpu/fetch.sh`로 `refs/tpu/pdf/`에 받는다(git에는 넣지 않는다).

요약: 최신 세대(v5e/v5p/v6e/Ironwood)는 아키텍처 파라미터만 공개되어 있고, 커널·레이어 단위
실측 golden은 공개된 것이 없다. 최근 TPU 검증 논문들도 Cloud TPU에서 직접 측정했고 데이터는
공개하지 않았다. 아키텍처와 실측이 함께 공개된 것은 TPUv1(카운터 기반 사이클 분해)과 TPU v4의
LLM 추론(PaLM)이다.

## 1. 아키텍처 파라미터

| 자료 | 세대 | 내용 | PDF |
|---|---|---|---|
| Jouppi, Lakshmanamurthy, Young, Patterson, *Google's Training Supercomputers from TPU v2 to Ironwood*, IEEE Micro 2026 ([arXiv 2606.15870](https://arxiv.org/abs/2606.15870)) | v2, v3, v4, v5p, Ironwood | Table 1: 세대별 MXU·VMEM·HBM·TensorCore/SparseCore·ICI | `jouppi2026_tpuv2_to_ironwood` |
| Jouppi 외, *In-Datacenter Performance Analysis of a Tensor Processing Unit*, ISCA 2017 ([arXiv 1704.04760](https://arxiv.org/abs/1704.04760)) | v1 | 전체 마이크로아키텍처 + 6개 앱 실측 (아래 3절) | `jouppi2017_tpuv1` |
| Jouppi 외, *TPU v4: An Optically Reconfigurable Supercomputer ...*, ISCA 2023 ([arXiv 2304.01433](https://arxiv.org/abs/2304.01433)) | v4 | 칩·OCS·SparseCore, MLPerf | `jouppi2023_tpuv4` |
| Jouppi 외, *Ten Lessons From Three Generations Shaped Google's TPUv4i*, ISCA 2021 ([DOI](https://doi.org/10.1109/ISCA52012.2021.00010)) | v4i | 추론 칩, CMEM 128 MiB, 앱별 성능 | (arXiv 없음) |
| Norrie 외, *The Design Process for Google's Training Chips: TPUv2 and TPUv3*, IEEE Micro 2021 ([DOI](https://doi.org/10.1109/MM.2021.3058217)) | v2, v3 | 코어 구조, MXU, HBM | (arXiv 없음) |
| [Cloud TPU 시스템 아키텍처](https://docs.cloud.google.com/tpu/docs/system-architecture-tpu-vm), [TPU7x 문서](https://docs.cloud.google.com/tpu/docs/tpu7x) | v4 ~ Ironwood | MXU 크기(v5p까지 128x128, v6e·Ironwood 256x256), SparseCore 수 | - |
| [How to Think About TPUs (JAX scaling book)](https://jax-ml.github.io/scaling-book/tpus/) | v5e, v6e 등 | 피크 FLOPS, HBM 대역폭, 실행 모델 | - |

### Jouppi 2026 Table 1 (그대로 옮김)

| | TPU v2 | TPU v3 | TPU v4 | TPU v5p | Ironwood |
|---|---|---|---|---|---|
| Year | 2017 | 2018 | 2021 | 2023 | 2025 |
| Peak BF16 TFLOPS | 46 | 123 | 275 | 459 | 2307 |
| Peak FP8 TFLOPS | - | - | - | 459 | 4614 |
| MXU per TPU (number, size) | 2, 128x128 | 4, 128x128 | 8, 128x128 | 8, 128x128 | 4, 256x256 bf16 / 4, 512x512 fp8 |
| VMEM per TPU (MiB) | 32 | 32 | 32 | 128 | 128 |
| HBM | HBM2 x2 | HBM2 x4 | HBM2 x4 | HBM2E x6 | HBM3E x8 |
| HBM capacity (GiB) | 16 | 32 | 32 | 96 | 192 |
| HBM BW (GB/s) | 700 | 900 | 1200 | 2765 | 7300 |
| TensorCores / SparseCores per TPU | 2 / 2 | 2 / 2 | 2 / 4 | 2 / 4 | 2 / 4 |
| ICI links x BW per link (GB/s) | 4 x 62 | 4 x 70 | 6 x 50 | 6 x 100 | 6 x 100 |
| Pod size (chips) | 256 | 1024 | 4096 | 8960 | 9216 |

본문: 벡터 유닛은 128 lane x 8 sublane, 레인당 ALU 2개(제한형)에서 4개(범용)로 늘었고, 벡터
레지스터는 v2~v5p의 8x128에서 Ironwood의 16x256으로 커졌다. 클럭은 공개하지 않는다.

주의: "MXU per TPU"를 그대로 쓰면 v4는 8 x 128² x 2 FLOP x 1.05 GHz = 275 TFLOPS로 맞지만,
Ironwood는 4 x 256² x 2 FLOP로 2307 TFLOPS가 되려면 약 4.4 GHz가 필요해 맞지 않는다. Cloud 문서도
256x256 MXU를 "128K MAC/cycle"로 적어(원소 65,536개) 단위가 다르다. Ironwood 모델은 MXU 수와
클럭을 피크 FLOPS로 교차 확인해야 한다.

## 2. TPU 실측과 비교한 시뮬레이터 (검증 방법론)

| 자료 | 측정 대상 | 측정 내용 | 보고 오차 | 데이터 공개 | PDF |
|---|---|---|---|---|---|
| Choi, Oh, *EONSim*, 2025 ([arXiv 2511.06679](https://arxiv.org/abs/2511.06679)) | TPU v6e | 추론 시간, on-chip 메모리 접근 수 | 1.4% / 2.2% | 없음 | `choi2025_eonsim` |
| Dang, Raj, Man, Tong, Krishna, *SCALE-Sim TPU*, 2026 ([arXiv 2603.22535](https://arxiv.org/abs/2603.22535)) | TPU v4 | GEMM 스윕(M·K·N 32~4096, 커널 시간, HBM 전송 제외), elementwise | GEMM R² 0.89 / MAPE 32% | 없음 | `dang2026_scalesim_tpu` |
| Yang 외, *PyTorchSim*, MICRO 2025 ([ACM](https://dl.acm.org/doi/10.1145/3725843.3756045)) | TPUv3 (2코어, 940 MHz, 코어당 128x128 x2) | GEMM 등 | MAE 11.5% | 코드 공개 | (ACM) |
| Zhang 외, *LLMCompass*, ISCA 2024 ([arXiv 2312.03134](https://arxiv.org/abs/2312.03134)) | TPUv3 8코어 (A100, MI210도) | 연산자, LLM prefill/decode | 연산자 10.9%, LLM 4.1% | 코드 공개 | `zhang2024_llmcompass` |
| Zhou 외, implicit convolution / TPUSim, IISWC 2021 ([arXiv 2110.03901](https://arxiv.org/abs/2110.03901)) | TPUv2 | GEMM, CNN 레이어 | 4.42% / 5.8% | 일부 | `zhou2021_implicit_conv_tpusim` |
| Park 외, *LENS*, 2026 ([arXiv 2606.18042](https://arxiv.org/abs/2606.18042)) | NPU 시스템 LLM 추론 | 전체 지연 | 2.15% | 미확인 | `park2026_lens` |
| Wang, Wei, Brooks, *ParaDnn*, MLSys 2020 ([arXiv 1907.10701](https://arxiv.org/abs/1907.10701)) | TPUv2/v3 | FC/CNN/RNN 파라미터 스윕 | (벤치마크) | 코드 공개 | `wang2020_paradnn_tpu_gpu_cpu` |

## 3. 공개 실측 데이터 (golden 후보)

### 3.1 LLM 추론: Pope 외, *Efficiently Scaling Transformer Inference*, MLSys 2023 ([arXiv 2211.05102](https://arxiv.org/abs/2211.05102), `pope2023_scaling_transformer_inference`)

TPU v4 (275 TFLOPS bf16, HBM 32 GiB / 1200 GB/s, ICI 270 GB/s, 3D torus), context 2048.
Prefill 지연은 2048 토큰 처리, decode 지연은 64 토큰 생성 기준. WS 2D = weight-stationary 2D,
WG XYZ = weight-gathered XYZ.

PaLM 540B (Table 2)

| Scenario | Phase | Chips | Batch | FFN | Attention | Weights | MFU | Latency |
|---|---|---|---|---|---|---|---|---|
| Low-latency | Prefill | 64 | 1 | WS 2D | Head | int8 | 43% | 0.29 s |
| Low-latency | Decode | 64 | 64 | WS 2D | Batch | int8 | 14% | 1.82 s |
| High-throughput | Prefill | 64 | 512 | WG XYZ | Batch | bf16 | 76% | 85.2 s |
| High-throughput | Decode | 64 | 512 | WS 2D | Batch | bf16 | 33% | 6.0 s |

PaLM 62B (Table 3)

| Scenario | Phase | Chips | Batch | FFN | Attention | Weights | MFU | Latency |
|---|---|---|---|---|---|---|---|---|
| Low-latency | Prefill | 16 | 1 | WS 2D | Head | int8 | 36% | 0.16 s |
| Low-latency | Decode | 16 | 32 | WS 2D | Batch | int8 | 8% | 0.73 s |
| High-throughput | Prefill | 32 | 512 | WG XYZ | Batch | bf16 | 73% | 20.2 s |
| High-throughput | Decode | 8 | 512 | WS 2D | Batch | bf16 | 37% | 5.1 s |

NPUsim 비교: `validation/palm/` (decode 4점, 결과와 해석은 그 README).

본문: 540B, 64칩, batch 64 decode는 int8 28.5 ms/token, bf16 36.9 ms/token. 모델 하이퍼파라미터는
Table D.1. 다중 칩 분할(partitioning)과 ICI 통신이 포함된 값이므로 단일 칩 모델과 비교하려면
분할 방식과 통신을 함께 모델링해야 한다.

### 3.2 TPUv1 하드웨어 카운터: Jouppi 외, ISCA 2017

256x256 int8 MAC, 700 MHz(피크 92 TOPS), unified buffer 24 MiB, accumulator 4 MiB, DDR3 8 GiB.

Table 1 (앱 사양)

| App | FC | Conv | Vector | Pool | Total layers | Nonlinear | Weights | Ops/weight byte | Batch |
|---|---|---|---|---|---|---|---|---|---|
| MLP0 | 5 | | | | 5 | ReLU | 20M | 200 | 200 |
| MLP1 | 4 | | | | 4 | ReLU | 5M | 168 | 168 |
| LSTM0 | 24 | | 34 | | 58 | sigmoid, tanh | 52M | 64 | 64 |
| LSTM1 | 37 | | 19 | | 56 | sigmoid, tanh | 34M | 96 | 96 |
| CNN0 | | 16 | | | 16 | ReLU | 8M | 2888 | 8 |
| CNN1 | 4 | 72 | | 13 | 89 | ReLU | 100M | 1750 | 32 |

Table 3 (하드웨어 카운터 기반 사이클 분해, 1·4·5·6행 합 100%)

| | MLP0 | MLP1 | LSTM0 | LSTM1 | CNN0 | CNN1 |
|---|---|---|---|---|---|---|
| Array active cycles | 12.7% | 10.6% | 8.2% | 10.5% | 78.2% | 46.2% |
| Useful MACs (% peak) | 12.5% | 9.4% | 8.2% | 6.3% | 78.2% | 22.5% |
| Unused MACs | 0.3% | 1.2% | 0.0% | 4.2% | 0.0% | 23.7% |
| Weight stall cycles | 53.9% | 44.2% | 58.1% | 62.1% | 0.0% | 28.1% |
| Weight shift cycles | 15.9% | 13.4% | 15.8% | 17.1% | 0.0% | 7.0% |
| Non-matrix cycles | 17.5% | 31.9% | 17.9% | 10.3% | 21.8% | 18.7% |
| RAW stalls | 3.3% | 8.4% | 14.6% | 10.6% | 3.5% | 22.8% |
| Input data stalls | 6.1% | 8.8% | 5.1% | 2.4% | 3.4% | 0.6% |
| TeraOps/s (92 peak) | 12.3 | 9.7 | 3.7 | 2.8 | 86.0 | 14.1 |

레이어 단위 shape은 공개되지 않아(앱 이름과 레이어 수·가중치 총량만) 앱 단위 비교만 가능하다.

### 3.3 기타

- **TpuGraphs** (Phothilimthana 외, NeurIPS 2023 Datasets, [arXiv 2308.13490](https://arxiv.org/abs/2308.13490), [GitHub](https://github.com/google-research-datasets/tpu_graphs)): XLA 그래프 x 컴파일 설정별 TPU 실측 실행 시간(ns). 그래프 단위라 사이클 시뮬레이터 검증에는 간접적. 관련: Kaufman 외, *A Learned Performance Model for TPUs*, MLSys 2021 ([arXiv 2008.01040](https://arxiv.org/abs/2008.01040)).
- **Google 공식 처리량**: Trillium Llama 2 70B 처리량이 v5e 대비 2.9배 ([Cloud 블로그](https://cloud.google.com/blog/products/compute/trillium-tpu-is-ga?hl=en)). v6e·Ironwood의 공식 MLPerf Llama 2 70B 결과는 찾지 못했다.
- **측정 도구**: [TPUMS](https://github.com/bytetwin/accelerator-microbenchmarks) (JAX, v6e·Ironwood의 GEMM 처리량·지연, 메모리, 인터커넥트).

## 4. Cloud TPU 직접 측정 비용 (v5e, v6e)

가격은 미국 리전 칩-시간 기준, 2026년 제3자 집계 자료에서 옮김. 공식 [가격 페이지](https://cloud.google.com/tpu/pricing)는
조사 시 표가 받아지지 않았으므로 사용 전 공식 페이지(또는 Cloud Billing Catalog API)로 확인할 것.

| | v5e | v6e (Trillium) |
|---|---|---|
| On-demand | 약 $1.20 | 약 $2.70 (europe-west4 $2.97, asia-northeast1 $3.24) |
| Spot | 약 $0.58 (미확인, 변동) | 약 $1.40 (미확인, 변동) |
| 1년 / 3년 약정 | - | 약 $1.89 / $1.22 |
| 최소 단위 | 1칩 (v5litepod-1) | 1칩 (v6e-1) |
| 칩당 HBM | 16 GiB | 32 GiB |

TPU VM 가격에는 호스트가 포함된다. 디스크·네트워크 비용은 작다. 새 프로젝트는 TPU quota를 먼저
요청해야 할 수 있다.

측정 계획별 추정 (개발·디버깅 포함, on-demand 기준. spot이면 약 절반):

| 계획 | 슬라이스 | 칩-시간 | v5e | v6e |
|---|---|---|---|---|
| GEMM·elementwise 스윕 (형태 수백 개, XLA profiler 커널 시간) | 1칩 | 10~20 | $12~24 | $27~54 |
| LLM 레이어·단일 칩 추론 (예: 8B 모델, v6e는 bf16, v5e는 int8) | 1칩 | 5~10 | $6~12 | $14~27 |
| 70B급 다중 칩 추론 (TP 8) | 8칩 x 2~3시간 | 16~24 | $19~29 | $43~65 |

무료 경로:

- [TPU Research Cloud](https://sites.research.google/trc/about/): 연구자 대상 무료 Cloud TPU 지원(신청제). 현재 v5e/v6e 포함 여부와 조건은 신청 페이지에서 확인할 것. 스토리지·VM 비용은 별도라는 보고가 있다.
- Kaggle Notebooks: [Kaggle TPU 문서](https://www.kaggle.com/docs/tpu) "up to 20 hours per week of TPUs and up to 9h at a time in a single session"(오래된 문서). 2026-10 사용자 보고([GitHub issue](https://github.com/TranBaDat2607/claude-code-kaggle-grandmaster/issues/18))에서 `kaggle quota` 출력이 "TPU 0.00h used, 20.00h remaining", 머신은 `TpuV5E8`(v5e 8칩), batch 세션 8.8시간 상한, TPU batch 세션 동시 1개. `TpuV5E1`/`TpuV6E1` 요청은 `TpuV5E8`로 바뀌어 v6e는 쓸 수 없다. profiler 사용 가능 여부는 확인 필요.

측정 키트: `validation/tpu_v5e/` (1칩은 Colab·TRC VM에서 `measure_tpu_v5e1.py`, 8칩은 Kaggle `TpuV5E8`에서 `measure_tpu_v5e8.py` 실행 -> `compare_v5e.py`로 NPUsim `tpu_v5e.cfg`와 비교).
