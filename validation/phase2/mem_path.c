// Sustained bandwidth of the Gemmini RTL memory path (DMA -> 512 KB L2 -> 64-bit AXI -> SimDRAM DDR3,
// +dramsim) under load, without compute: streams of mvin (reads) and mvout (writes) on never-touched
// buffers, timed from the first command to the fence that waits for the last one.
//   read  N      mvin of N bytes, 1 KB per command (16 rows of 64 B), N = 64 KB .. 4 MB (L2 is 512 KB)
//   read warm    the same 256 KB again (L2 hits)
//   write N      mvout of N bytes in 16 B rows (DIM int8), contiguous rows
//   write strided N  16 B rows 16 KB apart, as the C rows of a K = 16, N = 16384 GEMM tile
//   write strided64k N  rows 64 KB apart (N = 65536): all 16 rows fall in one 8-way L2 set
//   mixed N      read N and write N, commands interleaved (a memory-bound GEMM's traffic), N = 64 KB .. 1 MB
// Scratchpad addresses rotate so no two commands in flight touch the same rows.
// Build like gemm_sweep_cold.c; run on simulator-chipyard.harness-GemminiRocketConfig with +dramsim.
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "include/gemmini_testutils.h"

#define KB 1024
#define SPAD_ROWS (BANK_NUM * BANK_ROWS)

static elem_t rd64k[64 * KB] __attribute__((aligned(64)));
static elem_t rd256k[256 * KB] __attribute__((aligned(64)));
static elem_t rd1m[1024 * KB] __attribute__((aligned(64)));
static elem_t rd4m[4096 * KB] __attribute__((aligned(64)));
static elem_t wr256k[256 * KB] __attribute__((aligned(64)));
static elem_t wr1m[1024 * KB] __attribute__((aligned(64)));
static elem_t wrs256k[256 * KB] __attribute__((aligned(64)));      // strided rows: 16 B every 16 KB
static elem_t wrs1m[1024 * KB] __attribute__((aligned(64)));       // strided rows: 16 B every 64 KB
static elem_t mxr[256 * KB]__attribute__((aligned(64)));
static elem_t mxw[256 * KB] __attribute__((aligned(64)));
static elem_t mxr64k[64 * KB] __attribute__((aligned(64)));
static elem_t mxw64k[64 * KB] __attribute__((aligned(64)));
static elem_t mxr1m[1024 * KB] __attribute__((aligned(64)));
static elem_t mxw1m[1024 * KB] __attribute__((aligned(64)));

static size_t spad = 0;
static size_t next_rows(size_t rows) {         // rotating scratchpad rows
    if (spad + rows > SPAD_ROWS) spad = 0;
    size_t r = spad;
    spad += rows;
    return r;
}

// mvin `bytes` from buf: 16 rows x 64 B (4 x DIM columns) per command, rows contiguous
static void reads(elem_t *buf, size_t bytes) {
    for (size_t off = 0; off < bytes; off += 16 * 64)
        gemmini_extended_mvin(buf + off, next_rows(16 * 4), 64, 16);
}

// mvout `bytes` as 16 B rows, 16 rows per command, rows `stride` bytes apart
static void writes(elem_t *buf, size_t bytes, size_t stride, size_t span) {
    gemmini_config_st(stride);
    for (size_t row = 0; row < bytes / DIM; row += 16) {
        size_t r = row % (span / stride);           // strided: wrap to the next column block
        size_t col = (row / (span / stride)) * DIM;
        gemmini_extended_mvout(buf + r * stride + col, next_rows(16), DIM, 16);
    }
}

// read `bytes` and write the same amount as 16 B rows, commands interleaved per 1 KB
static void mixed(elem_t *rd, elem_t *wr, size_t bytes) {
    for (size_t off = 0; off < bytes; off += 16 * 64) {
        gemmini_extended_mvin(rd + off, next_rows(16 * 4), 64, 16);
        for (size_t w = 0; w < 4; w++)
            gemmini_extended_mvout(wr + off + w * 16 * DIM, next_rows(16), DIM, 16);
    }
}

static void report(const char *what, size_t bytes, uint64_t cycles) {
    printf("MEMPATH %s bytes=%lu cycles=%lu B_per_cycle_x100=%lu\n", what, (unsigned long)bytes,
           (unsigned long)cycles, (unsigned long)(100 * bytes / cycles));
}

#define TIME(what, bytes, body) do { gemmini_fence(); uint64_t t0 = read_cycles(); body; gemmini_fence(); \
                                     report(what, bytes, read_cycles() - t0); } while (0)

int main(void) {
    gemmini_flush(0);
    gemmini_config_ld(64);
    TIME("read", 64 * KB, reads(rd64k, 64 * KB));
    TIME("read", 256 * KB, reads(rd256k, 256 * KB));
    TIME("read_warm", 256 * KB, reads(rd256k, 256 * KB));
    TIME("read", 1024 * KB, reads(rd1m, 1024 * KB));
    TIME("read", 4096 * KB, reads(rd4m, 4096 * KB));
    TIME("write", 256 * KB, writes(wr256k, 256 * KB, DIM, 256 * KB));
    TIME("write", 1024 * KB, writes(wr1m, 1024 * KB, DIM, 1024 * KB));
    TIME("write_strided", 256 * KB, writes(wrs256k, 256 * KB, 16 * KB, 256 * KB));
    TIME("write_strided64k", 256 * KB, writes(wrs1m, 256 * KB, 64 * KB, 1024 * KB));
    gemmini_config_st(DIM);
    TIME("mixed", 128 * KB, mixed(mxr64k, mxw64k, 64 * KB));
    TIME("mixed", 512 * KB, mixed(mxr, mxw, 256 * KB));
    TIME("mixed", 2048 * KB, mixed(mxr1m, mxw1m, 1024 * KB));
    printf("MEMPATH DONE\n");
    exit(0);
}
