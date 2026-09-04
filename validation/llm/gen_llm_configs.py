#!/usr/bin/env python3
"""Generate NPUsim network + Gemmini mapping configs for LLM projection GEMMs.

Scope (agreed): projection-level only -- Q/K/V/O projections and FFN matmuls, which are
weight matmuls (input x weight) and map 1:1 onto NPUsim's [connected] layer. The
attention score/context matmuls (QK^T, S*V) are activation x activation and are OUT OF
SCOPE for this evaluation.

Each unique (K, N) GEMM shape becomes one single-layer network config so the shape is
independent (a [connected] layer takes K = network input_size; chaining layers would
couple shapes). Two M sets are emitted: decode (M=1) and prefill (M=128).

GEMM (M, K, N) -> NPUsim:
  [net] batch = M, input_size (height*width*channels) = K
  [connected] output = N
Gemmini 16x16 mapping -- BASELINE, execution-guaranteed (NOT performance-optimal):
  N_out = PE_X(16 spatial) x DRAM fold;  K_in = PE_Y(16 spatial) x DRAM fold;  M, folds
  all sit at DRAM so the on-chip (GLB/LB) tile is just the 16x16 spatial slab and always
  fits every shape's buffers. Performance-oriented mappings are the mapper's job
  (NeuroSpector etc., per evaluation.md); this baseline exists so every config RUNS and
  a roofline can be produced before mapper integration. Shapes are all multiples of 16,
  so the array fills exactly.

Run from repo root:  python3 validation/llm/gen_llm_configs.py
Writes configs/networks/llm_*.cfg and configs/mappings/gemmini/llm_*/ws.map
(configs/ is gitignored; force-add when committing, like the SFU configs).
"""
import os

ARRAY = 16  # Gemmini 16x16 systolic array

# Unique projection GEMM shapes per model, as (label, K, N). Deduplicated:
# Q==O and K==V share shapes where dims match; Gate==Up; etc. Source: model_config.md.
MODELS = {
    "bert": [          # BERT-base, 12 layers, hidden 768
        ("qkvo", 768, 768),
        ("fc1", 768, 3072),
        ("fc2", 3072, 768),
    ],
    "tinyllama": [     # TinyLlama-1.1B, 22 layers, hidden 2048, GQA kv 256
        ("qo", 2048, 2048),
        ("kv", 2048, 256),
        ("ffn_up", 2048, 5632),   # gate == up
        ("ffn_down", 5632, 2048),
    ],
    "qwen3": [         # Qwen3-0.6B, 28 layers
        ("q", 1024, 2048),
        ("kv", 1024, 1024),
        ("o", 2048, 1024),
        ("ffn_up", 1024, 3072),   # gate == up
        ("ffn_down", 3072, 1024),
    ],
    "gemma_slide": [   # Gemma E2B sliding-attention layer (28), hidden 1536
        ("q", 1536, 2048),
        ("kv", 1536, 256),
        ("o", 2048, 1536),
        ("ffn_up", 1536, 6144),   # gate == up
        ("ffn_down", 6144, 1536),
    ],
    "gemma_full": [    # Gemma E2B full-attention layer (7); FFN shared with sliding
        ("q", 1536, 4096),
        ("kv", 1536, 512),
        ("o", 4096, 1536),
    ],
}

M_SETS = {"decode": 1, "prefill": 128}

NET_TEMPLATE = """# {model} {label} projection GEMM  M(batch)={m} K(input)={k} N(output)={n}
# projection-level LLM workload (attention score/context are out of scope)
[net]
height={k}
width=1
channels=1
batch={m}

[data]
dataset=imagenet
test=datasets/gemm/test.lst
labels=datasets/gemm/labels.lst
top=1

[connected]
output={n}
activation=linear

[softmax]
groups=1

[cost]
type=l2
"""

MAP_HEADER = """# K, B, P, Q, C, R, S, H, W, GROUP, STRIDE
# {model} {label}: N(out)={n}=PE_X({px})xDRAM({nf}); K(in)={k}=PE_Y({py})xDRAM({kf}); M(batch)={m}
# BASELINE mapping (execution-guaranteed, not performance-optimal) -- see gen_llm_configs.py
[connected]
"""


def row(name, k, b, c):
    return f"{name:<7} = {k}, {b}, 1, 1, {c}, 1, 1, 0, 0, 1, 1,\n"


def make_mapping(model, label, k, n, m):
    px = min(ARRAY, n)
    nf = n // px         # N fold -> DRAM
    py = min(ARRAY, k)
    kf = k // py         # K fold -> DRAM
    assert px * nf == n and py * kf == k, f"{model} {label}: {k}x{n} not array-aligned"
    body = MAP_HEADER.format(model=model, label=label, n=n, px=px, nf=nf, k=k, py=py,
                             kf=kf, m=m)
    body += row("MAC", 1, 1, 1)
    body += row("PE", 1, 1, 1)          # no temporal fold in PE (keeps LB tile minimal)
    body += row("PE_X", px, 1, 1)       # N spatial (array width)
    body += row("PE_Y", 1, 1, py)       # K spatial (array height)
    body += row("GLB", 1, m, 1)         # batch staged on-chip (reuse; matches gemm256)
    body += row("CHIPS_X", 1, 1, 1)
    body += row("CHIPS_Y", 1, 1, 1)
    body += row("DRAM", nf, 1, kf)      # N fold + K fold off-chip (weight reload)
    return body


def main():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    net_dir = os.path.join(root, "configs", "networks")
    map_root = os.path.join(root, "configs", "mappings", "gemmini")
    count = 0
    names = []
    for model, shapes in MODELS.items():
        for label, k, n in shapes:
            for phase, m in M_SETS.items():
                name = f"llm_{model}_{label}_{phase}"
                with open(os.path.join(net_dir, name + ".cfg"), "w") as f:
                    f.write(NET_TEMPLATE.format(model=model, label=label, m=m, k=k, n=n))
                mdir = os.path.join(map_root, name)
                os.makedirs(mdir, exist_ok=True)
                with open(os.path.join(mdir, "ws.map"), "w") as f:
                    f.write(make_mapping(model, label, k, n, m))
                names.append(name)
                count += 1
    print(f"generated {count} network+mapping configs:")
    for nm in names:
        print(f"  {nm}")


if __name__ == "__main__":
    main()
