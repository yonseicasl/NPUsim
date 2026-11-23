// MOSAIC DRAM stream helper (components/dram.cc; NPUsim built with MOSAIC_ROOT). Modes:
//   one  YAML BYTES          DRAM ns of one sequential read stream, default DMA (64-transaction requests)
//   pair YAML_A YAML_B BYTES the same for two configurations IN ONE PROCESS, A then B then A again
//   diag YAML                samples for the attention consumer's read streams (see run_mosaic_kv_stream_test.sh)
// The process-static memo of the helper must be keyed by the configuration: `pair` equals two `one` runs.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

// components/dram.h
double mosaic_dram_stream_ns(const std::string &yaml, size_t bytes, bool is_write, size_t request_bytes,
                             unsigned max_inflight, double path_latency_ns);

static double ns(const std::string &yaml, size_t bytes) {
    return mosaic_dram_stream_ns(yaml, bytes, false, 0, 0, 0.0);
}

int main(int argc, char **argv) {
    const std::string mode = argc > 1 ? argv[1] : "";
    const size_t MiB = static_cast<size_t>(1) << 20;
    if(mode == "one" && argc == 4) {
        std::printf("%.6f\n", ns(argv[2], std::strtoull(argv[3], 0, 10)));
    } else if(mode == "pair" && argc == 5) {
        const size_t bytes = std::strtoull(argv[4], 0, 10);
        const double a = ns(argv[2], bytes), b = ns(argv[3], bytes), a2 = ns(argv[2], bytes);
        std::printf("%.6f %.6f %.6f\n", a, b, a2);
    } else if(mode == "diag" && argc == 3) {
        // 0 bytes cost nothing; a non-transaction-aligned stream is accepted and costs at least its aligned part
        if(ns(argv[2], 0) != 0.0) { std::fprintf(stderr, "0 bytes must cost 0\n"); return 1; }
        std::printf("unaligned,%zu,%.6f\n", 4*MiB + 13, ns(argv[2], 4*MiB + 13));
        // one stream of both halves against two fresh streams of one half each, same helper and settings
        for(size_t half = MiB; half <= 16*MiB; half *= 2) {
            std::printf("two_streams,%zu,%.6f,%.6f\n", half, ns(argv[2], half), ns(argv[2], 2*half));
        }
        // the long-stream rule: T(32 MiB) + (bytes - 32 MiB) x slope(16 -> 32 MiB), against an exact 64 MiB stream
        const double t16 = ns(argv[2], 16*MiB), t32 = ns(argv[2], 32*MiB), t64 = ns(argv[2], 64*MiB);
        std::printf("extrapolation,%.6f,%.6f,%.6f,%.6f\n", t16, t32, t32 + 32.0*MiB*(t32 - t16)/(16.0*MiB), t64);
    } else {
        std::fprintf(stderr, "usage: %s one YAML BYTES | pair YAML_A YAML_B BYTES | diag YAML\n", argv[0]);
        return 2;
    }
    return 0;
}
