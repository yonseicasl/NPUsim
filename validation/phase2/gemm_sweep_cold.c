// Gemmini RTL sweep with a cold L2: every point runs on its own A/B/C, which nothing has touched
// before (the bare-metal crt does not clear .bss), so each DMA read misses the 512 KB L2 and goes to
// SimDRAM (DRAMSim2 DDR3, run with +dramsim). gemm_sweep.c reuses one A/B/C for all its points, so
// after the first point part of its traffic hits the L2.
//
// Points: five of gemm_sweep.c's (cold counterparts of golden_rtl_cycles.csv) and two memory-bound
// ones -- K = 16 leaves 16 MACs per byte of B and C, below what the DDR3 delivers per compute cycle.
// Build against the DIM=16 gemmini_params.h (GemminiRocketConfig) with the bare-metal flags of
// gemmini-rocc-tests/bareMetalC/Makefile; run on simulator-chipyard.harness-GemminiRocketConfig with
// +dramsim.
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "include/gemmini_testutils.h"

#define POINT(name, M, K, N)                                       \
    static elem_t name##_A[M][K] row_align(1);                     \
    static elem_t name##_B[K][N] row_align(1);                     \
    static elem_t name##_C[M][N] row_align(1);                     \
    static void name(void) { bench(M, K, N, (elem_t*)name##_A, (elem_t*)name##_B, (elem_t*)name##_C); }

static void bench(size_t M, size_t K, size_t N, elem_t *A, elem_t *B, elem_t *C) {
    gemmini_flush(0);
    uint64_t start = read_cycles();
    tiled_matmul_auto(M, N, K,
                      A, B, NULL, C,
                      K, N, N, N,
                      MVIN_SCALE_IDENTITY, MVIN_SCALE_IDENTITY, MVIN_SCALE_IDENTITY,
                      NO_ACTIVATION, ACC_SCALE_IDENTITY, 0, false,
                      false, false,
                      false, false,
                      0,
                      WS);
    uint64_t end = read_cycles();
    printf("GEMMPOINT M=%lu K=%lu N=%lu cycles=%lu cold=1\n",
           (unsigned long)M, (unsigned long)K, (unsigned long)N, (unsigned long)(end - start));
}

POINT(p64, 64, 64, 64)
POINT(p128, 128, 128, 128)
POINT(p256, 256, 256, 256)
POINT(p16x512, 16, 512, 512)
POINT(p512x64, 512, 512, 64)
POINT(mem16k, 16, 16, 16384)
POINT(mem64k, 16, 16, 65536)

int main(void) {
    p64();
    p128();
    p256();
    p16x512();
    p512x64();
    mem16k();
    mem64k();
    printf("SWEEP DONE\n");
    exit(0);
}
