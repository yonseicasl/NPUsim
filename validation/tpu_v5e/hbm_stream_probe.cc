// HBM/DRAM probe for MOSAIC configs (validation/tpu_v5e/README.md, "MOSAIC HBM2e"). Sequential read stream through MOSAIC's CosimFrontend with a chosen request size (split into DRAM
// transactions by the frontend), at most `cap` requests outstanding (0 = as many as the frontend takes).
// build (from the MOSAIC root, after building NPUsim with MOSAIC_ROOT):
//   g++ -std=c++17 -O2 -Isimulator/src simulator/3rd/npusim/validation/tpu_v5e/hbm_stream_probe.cc \
//     simulator/3rd/npusim/library/obj/mosaic/{cosim_frontend,dram_config}.o simulator/build/lib/libramulator.so \
//     simulator/build/lib/libfmt.a -Wl,-rpath,$PWD/simulator/build/lib -o hbm_stream_probe
// usage: hbm_stream_probe <yaml> <bytes> <request bytes> <cap>
#include <cstdio>
#include <cstdlib>
#include "dram/cosim_frontend.h"
int main(int argc, char **argv) {
    const size_t bytes = std::strtoull(argv[2], 0, 10), req = std::strtoull(argv[3], 0, 10);
    const size_t cap = std::strtoull(argv[4], 0, 10);
    mosaicsim::CosimFrontend fe(argv[1]);
    size_t pos = 0, pushed = 0, done = 0;
    while (pos < bytes || done < pushed) {
        while (pos < bytes && (cap == 0 || pushed - done < cap) && !fe.full(false)) {
            fe.push(pos, false, req, &fe);
            pos += req; ++pushed;
        }
        fe.tick();
        while (fe.return_queue_pop()) ++done;
    }
    const double ns = fe.get_clk() * fe.get_tCK();
    std::printf("req=%zu cap=%zu: %.2f GB/s\n", req, cap, bytes / ns);
}
