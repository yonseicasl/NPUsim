// Gemmini DMA round trip on the RTL memory path (GemminiRocketConfig, +dramsim): cycles from issuing a
// one-row mvin to the fence that waits for it, for a never-touched 64 B line (L2 miss -> SimDRAM), the
// same line again (L2 hit), and an empty fence (overhead). Lines are consecutive, as in a stream.
// The L2-miss round trip, minus the DRAM access itself, is the path latency NPUsim's MOSAIC DRAM
// stream charges (gemmini_rtl_mem.cfg mosaic_path_latency_ns). Build like gemm_sweep_cold.c.
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "include/gemmini_testutils.h"

#define TRIALS 32
static elem_t lines[TRIALS][64] __attribute__((aligned(64)));   // one 64 B line per trial, untouched

static uint64_t median(uint64_t *v, int n) {
    for (int i = 1; i < n; i++)
        for (int j = i; j > 0 && v[j - 1] > v[j]; j--) { uint64_t t = v[j]; v[j] = v[j - 1]; v[j - 1] = t; }
    return v[n / 2];
}

int main(void) {
    static uint64_t fence[TRIALS], cold[TRIALS], warm[TRIALS];
    gemmini_flush(0);
    gemmini_config_ld(64);
    for (int i = 0; i < TRIALS; i++) {
        uint64_t t0 = read_cycles();
        gemmini_fence();
        fence[i] = read_cycles() - t0;

        t0 = read_cycles();
        gemmini_extended_mvin(&lines[i][0], 0, DIM, 1);
        gemmini_fence();
        cold[i] = read_cycles() - t0;

        t0 = read_cycles();
        gemmini_extended_mvin(&lines[i][0], 1, DIM, 1);
        gemmini_fence();
        warm[i] = read_cycles() - t0;
    }
    printf("DMALAT fence=%lu cold=%lu warm=%lu (median cycles of %d)\n", (unsigned long)median(fence, TRIALS),
           (unsigned long)median(cold, TRIALS), (unsigned long)median(warm, TRIALS), TRIALS);
    for (int i = 0; i < TRIALS; i++)
        printf("DMALAT_TRIAL %d cold=%lu warm=%lu\n", i, (unsigned long)cold[i], (unsigned long)warm[i]);
    exit(0);
}
