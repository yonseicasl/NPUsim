// gemm_sweep_cold.c with Gemmini's performance counters: how much the load, execute and store
// controllers overlap on a cold L2. Two builds (8 counters each):
//   PASS=1  the 7 disjoint load/store/execute Venn regions (Controller.scala incr_*_cycles) + RDMA bytes
//   PASS=2  each controller's active cycles, reservation-station-full cycles, DMA waits, WDMA bytes
// Every point runs on its own never-touched buffers, as in gemm_sweep_cold.c. Build with -DPASS=1 or 2.
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "include/gemmini_testutils.h"

#ifndef PASS
#define PASS 1
#endif

#if PASS == 1
static const int events[8] = {MAIN_LD_CYCLES, MAIN_ST_CYCLES, MAIN_EX_CYCLES, MAIN_LD_ST_CYCLES,
                              MAIN_LD_EX_CYCLES, MAIN_ST_EX_CYCLES, MAIN_LD_ST_EX_CYCLES, RDMA_BYTES_REC};
static const char *names[8] = {"LD", "ST", "EX", "LD_ST", "LD_EX", "ST_EX", "LD_ST_EX", "RDMA_BYTES"};
#else
static const int events[8] = {LOAD_ACTIVE_CYCLE, EXE_ACTIVE_CYCLE, STORE_ACTIVE_CYCLE,
                              RESERVATION_STATION_FULL_CYCLES, LOAD_DMA_WAIT_CYCLE, STORE_DMA_WAIT_CYCLE,
                              LOOP_MATMUL_ACTIVE_CYCLES, WDMA_BYTES_SENT};
static const char *names[8] = {"LOAD_ACTIVE", "EXE_ACTIVE", "STORE_ACTIVE", "RS_FULL", "LOAD_DMA_WAIT",
                               "STORE_DMA_WAIT", "LOOP_MATMUL_ACTIVE", "WDMA_BYTES"};
#endif

#define POINT(name, M, K, N)                                       \
    static elem_t name##_A[M][K] row_align(1);                     \
    static elem_t name##_B[K][N] row_align(1);                     \
    static elem_t name##_C[M][N] row_align(1);                     \
    static void name(void) { bench(M, K, N, (elem_t*)name##_A, (elem_t*)name##_B, (elem_t*)name##_C); }

static void bench(size_t M, size_t K, size_t N, elem_t *A, elem_t *B, elem_t *C) {
    gemmini_flush(0);
    for (int i = 0; i < 8; i++) counter_configure(i, events[i]);
    counter_reset();
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
    printf("PERF%d M=%lu K=%lu N=%lu cycles=%lu", PASS, (unsigned long)M, (unsigned long)K, (unsigned long)N,
           (unsigned long)(end - start));
    for (int i = 0; i < 8; i++) printf(" %s=%u", names[i], counter_read(i));
    printf("\n");
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
