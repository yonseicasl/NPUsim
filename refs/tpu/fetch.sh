#!/usr/bin/env bash
# Download the arXiv papers indexed in README.md into refs/tpu/pdf/ (not tracked by git).
# Papers that are not on arXiv (IEEE/ACM only) are linked from README.md instead.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p pdf
while read -r id name; do
    [[ -f "pdf/$name.pdf" ]] || curl -sfL "https://arxiv.org/pdf/$id" -o "pdf/$name.pdf"
done <<'EOF'
2606.15870 jouppi2026_tpuv2_to_ironwood
1704.04760 jouppi2017_tpuv1
2304.01433 jouppi2023_tpuv4
2211.05102 pope2023_scaling_transformer_inference
2603.22535 dang2026_scalesim_tpu
2511.06679 choi2025_eonsim
2312.03134 zhang2024_llmcompass
2110.03901 zhou2021_implicit_conv_tpusim
2606.18042 park2026_lens
2308.13490 phothilimthana2023_tpugraphs
2008.01040 kaufman2021_learned_tpu_perf_model
1907.10701 wang2020_paradnn_tpu_gpu_cpu
EOF
ls -la pdf
