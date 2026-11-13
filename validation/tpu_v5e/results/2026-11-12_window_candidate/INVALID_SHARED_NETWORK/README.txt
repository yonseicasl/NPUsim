INVALID -- do not cite. Everything here was produced by balance_npusim.py while it gave every case the same
network name ("g"): gemm_files() writes configs/networks/<name>.cfg and weights/<name>.wgh, shared by all accels,
so each case overwrote the previous shape and every run used the LAST shape's network (e.g. batch 512, 8192 x 8192)
with its own mapping. Found by Codex in the live files of the running batch (2026-11-12). Kept only as a record of
the bug; the latencies, error tables and any explanation derived from these reports are withdrawn.
