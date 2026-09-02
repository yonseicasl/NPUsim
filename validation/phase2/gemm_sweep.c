// Gemmini RTL validation sweep: prints WS-dataflow cycle counts for the same
// GEMM points NPUsim runs. Copy into gemmini-rocc-tests/bareMetalC/ and add
// "gemm_sweep" to that directory's Makefile tests list.
// NOTE: tiled_matmul_auto's argument list differs across gemmini versions --
// adjust to match the checked-out include/gemmini.h if compilation fails.
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef BAREMETAL
#include <sys/mman.h>
#endif
#include "include/gemmini_testutils.h"

#define MAXD 512
static elem_t A[MAXD][MAXD] row_align(1);
static elem_t B[MAXD][MAXD] row_align(1);
static elem_t C[MAXD][MAXD] row_align(1);

static void bench(size_t M, size_t K, size_t N) {
    gemmini_flush(0);
    uint64_t start = read_cycles();
    tiled_matmul_auto(M, N, K,
                      (elem_t*)A, (elem_t*)B, NULL, (elem_t*)C,
                      MAXD, MAXD, MAXD, MAXD,
                      MVIN_SCALE_IDENTITY, MVIN_SCALE_IDENTITY, MVIN_SCALE_IDENTITY,
                      NO_ACTIVATION, ACC_SCALE_IDENTITY, 0,
                      false, false,
                      false, false,
                      0,
                      WS);
    uint64_t end = read_cycles();
    printf("GEMMPOINT M=%lu K=%lu N=%lu cycles=%lu\n",
           (unsigned long)M, (unsigned long)K, (unsigned long)N,
           (unsigned long)(end - start));
}

int main(void) {
#ifndef BAREMETAL
    if (mlockall(MCL_CURRENT | MCL_FUTURE) != 0) { perror("mlockall failed"); exit(1); }
#endif
    bench(64, 64, 64);
    bench(128, 128, 128);
    bench(256, 256, 256);
    bench(16, 512, 512);
    bench(512, 512, 64);
    bench(512, 512, 512);
    printf("SWEEP DONE\n");
    exit(0);
}
