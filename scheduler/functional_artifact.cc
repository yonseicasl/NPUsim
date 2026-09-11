#include "functional_artifact.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>

#include "workload_graph.h"

namespace {

// Compact SHA-256 (FIPS 180-4). Self-contained so artifact hashing needs no new
// external dependency; toy-fixture payloads are small, so throughput is irrelevant.
struct sha256_state_t {
    uint32_t h[8];
    uint64_t bits;
    uint8_t block[64];
    size_t fill;
};

const uint32_t sha256_k[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

inline uint32_t rotr(uint32_t x, unsigned n) { return (x >> n) | (x << (32 - n)); }

void sha256_compress(sha256_state_t *s, const uint8_t *p) {
    uint32_t w[64];
    for(unsigned i = 0; i < 16; ++i) {
        w[i] = (uint32_t(p[i*4]) << 24) | (uint32_t(p[i*4+1]) << 16) |
               (uint32_t(p[i*4+2]) << 8) | uint32_t(p[i*4+3]);
    }
    for(unsigned i = 16; i < 64; ++i) {
        const uint32_t s0 = rotr(w[i-15], 7) ^ rotr(w[i-15], 18) ^ (w[i-15] >> 3);
        const uint32_t s1 = rotr(w[i-2], 17) ^ rotr(w[i-2], 19) ^ (w[i-2] >> 10);
        w[i] = w[i-16] + s0 + w[i-7] + s1;
    }
    uint32_t a = s->h[0], b = s->h[1], c = s->h[2], d = s->h[3];
    uint32_t e = s->h[4], f = s->h[5], g = s->h[6], h = s->h[7];
    for(unsigned i = 0; i < 64; ++i) {
        const uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        const uint32_t ch = (e & f) ^ (~e & g);
        const uint32_t t1 = h + S1 + ch + sha256_k[i] + w[i];
        const uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        const uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
        const uint32_t t2 = S0 + mj;
        h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    s->h[0] += a; s->h[1] += b; s->h[2] += c; s->h[3] += d;
    s->h[4] += e; s->h[5] += f; s->h[6] += g; s->h[7] += h;
}

} // namespace

std::string sha256_hex(const void *m_data, size_t m_bytes) {
    sha256_state_t s;
    const uint32_t init[8] = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
                              0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    std::memcpy(s.h, init, sizeof(init));
    s.bits = uint64_t(m_bytes)*8; s.fill = 0;
    const uint8_t *p = static_cast<const uint8_t*>(m_data);
    size_t remaining = m_bytes;
    while(remaining >= 64) { sha256_compress(&s, p); p += 64; remaining -= 64; }
    uint8_t tail[128];
    std::memcpy(tail, p, remaining);
    tail[remaining] = 0x80;
    size_t total = remaining + 1;
    while(total % 64 != 56) tail[total++] = 0;
    for(int i = 7; i >= 0; --i) tail[total++] = uint8_t(s.bits >> (i*8));
    for(size_t off = 0; off < total; off += 64) sha256_compress(&s, tail + off);
    char hex[65];
    for(unsigned i = 0; i < 8; ++i) std::snprintf(hex + i*8, 9, "%08x", s.h[i]);
    return std::string(hex, 64);
}

namespace {

void artifact_fail(const std::string &m_message) {
    std::cerr << "Error: tensor artifact: " << m_message << std::endl;
    exit(1);
}

std::string manifest_directory(const std::string &m_path) {
    const std::string::size_type slash = m_path.find_last_of("/\\");
    return slash == std::string::npos ? std::string(".") : m_path.substr(0, slash);
}

// Read one payload file (raw little-endian float32), verify byte size, digest and the
// NaN/Inf policy, and return the values.
std::vector<float> read_payload(const std::string &m_path, size_t m_expected_elements,
                                const std::string &m_expected_sha, const std::string &m_what) {
    std::ifstream in(m_path.c_str(), std::ios::binary | std::ios::ate);
    if(!in.is_open()) artifact_fail(m_what + ": cannot open payload " + m_path);
    const std::streamsize bytes = in.tellg();
    if(static_cast<size_t>(bytes) != m_expected_elements*sizeof(float)) {
        artifact_fail(m_what + ": payload " + m_path + " holds " +
                      std::to_string(bytes/sizeof(float)) + " floats, manifest shape needs " +
                      std::to_string(m_expected_elements));
    }
    in.seekg(0);
    std::vector<float> values(m_expected_elements);
    in.read(reinterpret_cast<char*>(values.data()), bytes);
    if(!in.good()) artifact_fail(m_what + ": short read on " + m_path);
    const std::string digest = sha256_hex(values.data(), static_cast<size_t>(bytes));
    if(digest != m_expected_sha) {
        artifact_fail(m_what + ": payload " + m_path + " sha256 " + digest +
                      " does not match the manifest's " + m_expected_sha);
    }
    for(size_t i = 0; i < values.size(); ++i) {
        if(!std::isfinite(values[i])) {
            artifact_fail(m_what + ": payload " + m_path + " has NaN/Inf at element " +
                          std::to_string(i));
        }
    }
    return values;
}

} // namespace

void functional_artifact_t::load(const std::string &m_manifest_path,
                                 const workload_graph_t &m_graph) {
    boost::property_tree::ptree root;
    try {
        boost::property_tree::read_json(m_manifest_path, root);
    } catch(const boost::property_tree::json_parser::json_parser_error &error) {
        artifact_fail("cannot parse " + m_manifest_path + ": " + error.message());
    }
    schema_version = root.get<std::string>("schema_version", "");
    if(schema_version != "npusim.tensor.v1") {
        artifact_fail("schema_version '" + schema_version + "' is not npusim.tensor.v1");
    }
    executable_sha256 = root.get<std::string>("executable_sha256", "");
    if(executable_sha256.empty() || executable_sha256 != m_graph.executable_sha256) {
        artifact_fail("executable_sha256 '" + executable_sha256 +
                      "' does not match the loaded executable ('" +
                      m_graph.executable_sha256 + "'); the values belong to a different"
                      " executable");
    }
    const std::string base = manifest_directory(m_manifest_path);

    // The coverage contract: every graph input and every parameter/buffer/constant of the
    // executable must be supplied exactly once; activations must not be.
    std::set<std::string> required;
    for(const workload_tensor_t &tensor : m_graph.tensors) {
        if(tensor.kind == "parameter" || tensor.kind == "buffer" || tensor.kind == "constant")
            required.insert(tensor.id);
    }
    for(const std::string &input : m_graph.inputs) required.insert(input);

    if(root.find("tensors") == root.not_found()) artifact_fail("manifest lists no tensors");
    for(const auto &entry : root.get_child("tensors")) {
        const boost::property_tree::ptree &node = entry.second;
        const std::string id = node.get<std::string>("id", "");
        if(id.empty()) artifact_fail("tensor entry without id");
        if(tensors.count(id)) artifact_fail("tensor " + id + " supplied twice");
        const std::string role = node.get<std::string>("role", "");
        const std::string dtype = node.get<std::string>("dtype", "");
        if(dtype != "float32") {
            artifact_fail("tensor " + id + " dtype '" + dtype +
                          "' unsupported (npusim.tensor.v1 carries float32)");
        }
        if(!required.count(id)) {
            artifact_fail("tensor " + id + " is not a graph input/parameter of this"
                          " executable (activations are computed, not supplied)");
        }
        const workload_tensor_t &declared = m_graph.tensor(id);
        const std::string expected_role =
            declared.kind == "activation" ? "input" : declared.kind;
        if(role != expected_role) {
            artifact_fail("tensor " + id + " role '" + role + "' disagrees with the"
                          " executable's kind '" + expected_role + "'");
        }
        size_t elements = 1;
        if(node.find("shape") == node.not_found()) artifact_fail("tensor " + id + " has no shape");
        std::vector<size_t> shape;
        for(const auto &dim : node.get_child("shape"))
            shape.push_back(dim.second.get_value<size_t>());
        for(size_t d : shape) elements *= d;
        if(shape != declared.shape) {
            artifact_fail("tensor " + id + " manifest shape disagrees with the executable");
        }
        const std::string payload = node.get<std::string>("payload", "");
        const std::string sha = node.get<std::string>("sha256", "");
        if(payload.empty() || sha.empty())
            artifact_fail("tensor " + id + " needs payload and sha256");
        tensors[id] = read_payload(base + "/" + payload, elements, sha, "tensor " + id);
        tensor_roles[id] = role;
    }
    for(const std::string &id : required) {
        if(!tensors.count(id)) artifact_fail("required tensor " + id + " has no payload");
    }

    // Goldens: at least the graph output's producing operation must be covered so the
    // acceptance gate always compares something.
    std::set<std::string> output_producers;
    for(const workload_operation_t &operation : m_graph.operations) {
        for(const std::string &out : operation.outputs) {
            for(const std::string &graph_out : m_graph.outputs) {
                const workload_tensor_t &storage = m_graph.storage_tensor(graph_out);
                if(out == graph_out || out == storage.id) output_producers.insert(operation.id);
            }
        }
    }
    if(root.find("golden") != root.not_found()) {
        for(const auto &entry : root.get_child("golden")) {
            const boost::property_tree::ptree &node = entry.second;
            const std::string operation_id = node.get<std::string>("operation_id", "");
            const std::string stage = node.get<std::string>("stage", "post_activation");
            if(stage != "post_activation")
                artifact_fail("golden for " + operation_id + " uses unsupported stage " + stage);
            const workload_operation_t *operation = NULL;
            for(const workload_operation_t &candidate : m_graph.operations)
                if(candidate.id == operation_id) operation = &candidate;
            if(operation == NULL)
                artifact_fail("golden references unknown operation " + operation_id);
            if(golden.count(operation_id))
                artifact_fail("operation " + operation_id + " has two goldens");
            const workload_tensor_t &out = m_graph.tensor(operation->outputs.front());
            const std::string payload = node.get<std::string>("payload", "");
            const std::string sha = node.get<std::string>("sha256", "");
            if(payload.empty() || sha.empty())
                artifact_fail("golden " + operation_id + " needs payload and sha256");
            golden[operation_id] = read_payload(base + "/" + payload, out.elements(), sha,
                                                "golden " + operation_id);
        }
    }
    for(const std::string &producer : output_producers) {
        if(!golden.count(producer)) {
            artifact_fail("graph-output operation " + producer + " has no golden; the"
                          " acceptance gate requires every graph output to be compared");
        }
    }
    std::cout << "# Tensor artifact: " << tensors.size() << " payload(s), "
              << golden.size() << " golden(s), executable " << executable_sha256.substr(0, 16)
              << "..." << std::endl;
}
