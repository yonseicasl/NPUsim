#ifndef __FUNCTIONAL_ARTIFACT_H__
#define __FUNCTIONAL_ARTIFACT_H__

// G1 (functional-sim gaps plan Step 4): the npusim.tensor.v1 artifact -- the value
// companion of an npusim.exec.v1 executable. The manifest binds real tensor payloads
// (graph inputs, parameters, buffers, constants) and per-operation CPU goldens to one
// specific executable by hash. Everything is validated BEFORE simulation: schema and
// executable hash, per-tensor dtype/shape/byte-size, payload sha256, NaN/Inf policy,
// and exact coverage (every input and parameter provided exactly once, no extras).

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

class workload_graph_t;

// Streaming SHA-256 (FIPS 180-4), used to pin payloads to their manifest digests.
std::string sha256_hex(const void *m_data, size_t m_bytes);

struct functional_artifact_t {
    std::string schema_version;
    std::string executable_sha256;
    // Value payloads keyed by executable tensor id (float32, little-endian).
    std::map<std::string, std::vector<float>> tensors;
    std::map<std::string, std::string> tensor_roles;   // id -> input|parameter|buffer|constant
    // Per-operation goldens keyed by operation id (stage post_activation).
    std::map<std::string, std::vector<float>> golden;
    // Arithmetic semantics of the artifact (manifest "semantics" block; absent = fp32).
    // int8: integer operands stored in float, optional requant shift/clamp and
    // asymmetric zero points, restricted to all-mapped (linear/conv) executables.
    // fp16/bf16: operands pre-rounded to the grid, outputs rounded by the simulator.
    std::string profile;                                // fp32 | int8 | fp16 | bf16
    int requant_shift;
    int requant_min;
    int requant_max;
    int input_zero_point;
    int weight_zero_point;

    functional_artifact_t() : profile("fp32"), requant_shift(0), requant_min(-127),
                              requant_max(127), input_zero_point(0), weight_zero_point(0) {}

    // Parse + validate a manifest against the loaded executable. Any violation prints a
    // diagnostic and exits non-zero -- an invalid artifact must never reach compute.
    void load(const std::string &m_manifest_path, const workload_graph_t &m_graph);
};

#endif
