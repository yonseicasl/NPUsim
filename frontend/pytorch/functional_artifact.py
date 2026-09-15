"""npusim.tensor.v1 -- the value artifact that accompanies an npusim.exec.v1 executable.

The manifest binds real tensor payloads (graph inputs, parameters, buffers, constants)
and per-operation CPU goldens to ONE executable by sha256. Goldens come from a pure-
Python reference interpreter of the executable IR itself (deterministic scalar CPU,
fp32-rounded like the simulator datapath); when PyTorch is available, parameter/input
values come from the exported program and the interpreter's final outputs are cross-
checked against the eager model's forward pass.

No numpy/torch is required for `synthesize` (seeded deterministic values), so the C++
acceptance suite can exercise the executable-IR functional path on any machine.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
import struct
from pathlib import Path
from typing import Any, Mapping

ARTIFACT_SCHEMA_VERSION = "npusim.tensor.v1"


def _f32(value: float) -> float:
    return struct.unpack("<f", struct.pack("<f", value))[0]


def _round_bf16(value: float) -> float:
    x = struct.unpack("<I", struct.pack("<f", value))[0]
    if (x & 0x7fffffff) > 0x7f800000:
        return value
    x = (x + 0x7FFF + ((x >> 16) & 1)) & 0xFFFFFFFF
    return struct.unpack("<f", struct.pack("<I", x & 0xFFFF0000))[0]


def _round_fp16(value: float) -> float:
    return struct.unpack("<e", struct.pack("<e", value))[0]


def _round_lowp(profile: str, value: float) -> float:
    if profile == "bf16":
        return _round_bf16(value)
    if profile == "fp16":
        return _round_fp16(value)
    return value


def _pack(values: list[float]) -> bytes:
    return struct.pack("<%df" % len(values), *values)


def _numel(shape: list[int]) -> int:
    total = 1
    for dimension in shape:
        total *= int(dimension)
    return total


class FunctionalArtifactError(ValueError):
    """Raised when an artifact cannot be produced consistently."""


# --------------------------------------------------------------------------------------
# Reference interpreter of npusim.exec.v1 (the canonical scalar CPU oracle, plan §5.1).
# --------------------------------------------------------------------------------------

def _apply_activation(values: list[float], activation: str) -> list[float]:
    if activation == "relu":
        return [v if v > 0.0 else 0.0 for v in values]
    if activation == "leaky":
        return [v if v > 0.0 else _f32(0.1 * v) for v in values]
    return values


def _linear(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    x = fetch(op["inputs"][0])
    w = fetch(op["inputs"][1])
    bias = fetch(op["inputs"][2]) if len(op["inputs"]) >= 3 else None
    g = op["geometry"]
    m_rows, k, n = g["batch"], g["input_features"], g["output_features"]
    if len(x) != m_rows * k or len(w) != n * k:
        raise FunctionalArtifactError(f"linear {op['id']} operand sizes disagree")
    out = [0.0] * (m_rows * n)
    for m in range(m_rows):
        for j in range(n):
            acc = 0.0
            base_x, base_w = m * k, j * k
            for c in range(k):
                acc = _f32(acc + _f32(x[base_x + c] * w[base_w + c]))
            if bias is not None:
                acc = _f32(acc + bias[j])
            out[m * n + j] = acc
    return out


def _conv2d(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    x = fetch(op["inputs"][0])
    w = fetch(op["inputs"][1])
    bias = fetch(op["inputs"][2]) if len(op["inputs"]) >= 3 else None
    g = op["geometry"]
    B, C, H, W = g["batch"], g["input_channels"], g["input_height"], g["input_width"]
    N, P, Q = g["output_channels"], g["output_height"], g["output_width"]
    R, S = g["filter_height"], g["filter_width"]
    sh, sw = g["stride_height"], g["stride_width"]
    ph, pw = g["padding_height"], g["padding_width"]
    dh, dw = g["dilation_height"], g["dilation_width"]
    groups = g["groups"]
    Cg, Ng = C // groups, N // groups
    out = [0.0] * (B * N * P * Q)
    for b in range(B):
        for n in range(N):
            grp = n // Ng
            for p in range(P):
                for q in range(Q):
                    acc = 0.0
                    for c in range(Cg):
                        for r in range(R):
                            ih = p * sh + r * dh - ph
                            if ih < 0 or ih >= H:
                                continue
                            for s in range(S):
                                iw = q * sw + s * dw - pw
                                if iw < 0 or iw >= W:
                                    continue
                                acc = _f32(acc + _f32(
                                    x[((b * C + grp * Cg + c) * H + ih) * W + iw]
                                    * w[((n * Cg + c) * R + r) * S + s]))
                    if bias is not None:
                        acc = _f32(acc + bias[n])
                    out[((b * N + n) * P + p) * Q + q] = acc
    return out


def _softmax(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    x = fetch(op["inputs"][0])
    g = op["geometry"]
    rows, length = g["rows"], g["row_length"]
    # B-4 causal mask: within each [Tq][Tk] block, query r attends to keys 0..(r%span).
    span = (g.get("causal_span") or length) if g.get("causal") else 0
    out = [0.0] * (rows * length)
    for r in range(rows):
        valid = min(length, (r % span) + 1) if span else length
        row = x[r * length:r * length + valid]
        peak = max(row)
        exps = [_f32(math.exp(_f32(v - peak))) for v in row]
        total = _f32(sum(exps))
        for i, e in enumerate(exps):
            out[r * length + i] = _f32(e / total)
    return out


def _pool2d(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    x = fetch(op["inputs"][0])
    g = op["geometry"]
    B, C = g["batch"], g["channels"]
    H, W = g["input_height"], g["input_width"]
    P, Q = g["output_height"], g["output_width"]
    max_mode = g["mode"] == "max"
    out = [0.0] * (B * C * P * Q)
    for b in range(B):
        for c in range(C):
            for p in range(P):
                for q in range(Q):
                    best, total, valid = None, 0.0, 0
                    for kh in range(g["kernel_height"]):
                        ih = p * g["stride_height"] - g["padding_height"] + kh * g["dilation_height"]
                        if ih < 0 or ih >= H:
                            continue
                        for kw in range(g["kernel_width"]):
                            iw = q * g["stride_width"] - g["padding_width"] + kw * g["dilation_width"]
                            if iw < 0 or iw >= W:
                                continue
                            v = x[((b * C + c) * H + ih) * W + iw]
                            best = v if best is None or v > best else best
                            total += v
                            valid += 1
                    if max_mode:
                        out[((b * C + c) * P + p) * Q + q] = best if best is not None else 0.0
                    else:
                        window = g["kernel_height"] * g["kernel_width"]
                        samples = window if g.get("count_include_pad", True) else valid
                        out[((b * C + c) * P + p) * Q + q] = _f32(total / samples) if samples else 0.0
    return out


def _elementwise(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    a, b = fetch(op["inputs"][0]), fetch(op["inputs"][1])
    if op["geometry"]["operator"] == "multiply":
        return [_f32(x * y) for x, y in zip(a, b)]
    return [_f32(x + y) for x, y in zip(a, b)]


def _concat(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    axis = op["geometry"]["axis"]
    shapes = [shape_of(t) for t in op["inputs"]]
    outer = 1
    for d in shapes[0][:axis]:
        outer *= d
    inners = []
    for shape in shapes:
        inner = 1
        for d in shape[axis:]:
            inner *= d
        inners.append(inner)
    out: list[float] = []
    blocks = [fetch(t) for t in op["inputs"]]
    for o in range(outer):
        for block, inner in zip(blocks, inners):
            out.extend(block[o * inner:(o + 1) * inner])
    return out


def _batch_norm(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    if len(op["inputs"]) != 5:
        raise FunctionalArtifactError(
            f"batch_norm {op['id']} needs the full affine form (x, weight, bias, mean, var)")
    x = fetch(op["inputs"][0])
    scale, shift = fetch(op["inputs"][1]), fetch(op["inputs"][2])
    mean, var = fetch(op["inputs"][3]), fetch(op["inputs"][4])
    g = op["geometry"]
    channels = g["channels"]
    batch = shape_of(op["inputs"][0])[0]
    spatial = len(x) // batch // channels
    eps = g["epsilon"]
    out = [0.0] * len(x)
    for i, v in enumerate(x):
        c = (i // spatial) % channels
        out[i] = _f32(scale[c] * _f32(v - mean[c]) / _f32(math.sqrt(_f32(var[c] + eps))) + shift[c])
    return out


def _layer_norm(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    x = fetch(op["inputs"][0])
    gamma, beta = fetch(op["inputs"][1]), fetch(op["inputs"][2])
    L = op["geometry"]["normalized_size"]
    eps = op["geometry"]["epsilon"]
    if len(gamma) != L or len(beta) != L:
        raise FunctionalArtifactError(f"layer_norm {op['id']} weight/bias length != normalized size")
    out = [0.0] * len(x)
    for r in range(len(x) // L):
        row = x[r * L:(r + 1) * L]
        mean = sum(row) / L
        var = sum((v - mean) ** 2 for v in row) / L
        inv = 1.0 / math.sqrt(var + eps)
        for i, v in enumerate(row):
            out[r * L + i] = _f32(_f32((v - mean) * inv) * gamma[i] + beta[i])
    return out


def _matmul(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    g = op["geometry"]
    Bt, M, K, N = g["matmul_batch"], g["matmul_m"], g["matmul_k"], g["matmul_n"]
    tb = bool(g.get("matmul_transpose_b", False))
    a, b = fetch(op["inputs"][0]), fetch(op["inputs"][1])
    out = [0.0] * (Bt * M * N)
    for bt in range(Bt):
        ab, bb, cb = bt * M * K, bt * K * N, bt * M * N
        for m in range(M):
            for n in range(N):
                acc = 0.0
                for k in range(K):
                    bv = b[bb + n * K + k] if tb else b[bb + k * N + n]
                    acc = _f32(acc + _f32(a[ab + m * K + k] * bv))
                out[cb + m * N + n] = acc
    return out


def _transpose(op: Mapping[str, Any], fetch, shape_of) -> list[float]:
    # B-4: swap two axes (multi-head split/merge). Row-major reorder.
    x = fetch(op["inputs"][0])
    shape = shape_of(op["inputs"][0])
    a0, a1 = op["geometry"]["axis0"], op["geometry"]["axis1"]
    in_stride = [1] * len(shape)
    for d in range(len(shape) - 2, -1, -1):
        in_stride[d] = in_stride[d + 1] * shape[d + 1]
    out_shape = list(shape); out_shape[a0], out_shape[a1] = out_shape[a1], out_shape[a0]
    out_stride = [1] * len(shape)
    for d in range(len(shape) - 2, -1, -1):
        out_stride[d] = out_stride[d + 1] * out_shape[d + 1]
    out = [0.0] * len(x)
    for flat in range(len(x)):
        rem, idx = flat, [0] * len(shape)
        for d in range(len(shape)):
            idx[d] = rem // in_stride[d]; rem %= in_stride[d]
        idx[a0], idx[a1] = idx[a1], idx[a0]
        out[sum(i * s for i, s in zip(idx, out_stride))] = x[flat]
    return out


_KERNELS = {
    "npusim.linear": _linear,
    "npusim.conv2d": _conv2d,
    "npusim.softmax": _softmax,
    "npusim.pool2d": _pool2d,
    "npusim.elementwise": _elementwise,
    "npusim.concat": _concat,
    "npusim.batch_norm": _batch_norm,
    "npusim.layer_norm": _layer_norm,
    "npusim.matmul": _matmul,
    "npusim.transpose": _transpose,
}


def interpret_executable(
    executable: Mapping[str, Any],
    values: Mapping[str, list[float]],
    semantics: Mapping[str, Any] | None = None,
) -> dict[str, list[float]]:
    """Run every executable operation on the reference kernels; returns {op_id: output}.

    `semantics` mirrors the manifest's semantics block (profile int8/fp16/bf16). int8:
    the zero points are subtracted from the operand copies before the MAC kernels (the
    payloads carry the RAW offset operands), bias adds in, and requant shift/clamp runs
    after the fused activation -- the same order the simulator's finalize uses. fp16/
    bf16: every operation's output is rounded to the reduced-mantissa grid.
    """
    semantics = dict(semantics or {})
    profile = semantics.get("profile", "fp32")
    shift = int(semantics.get("requant_shift", 0))
    requant_min = int(semantics.get("requant_min", -127))
    requant_max = int(semantics.get("requant_max", 127))
    zp_i = int(semantics.get("input_zero_point", 0))
    zp_w = int(semantics.get("weight_zero_point", 0))
    tensors = {t["id"]: t for t in executable["tensors"]}
    weight_ids: set[str] = set()
    data_input_ids: set[str] = set()
    for op in executable["operations"]:
        if op["kind"] in ("npusim.linear", "npusim.conv2d"):
            data_input_ids.add(op["inputs"][0])
            if len(op["inputs"]) >= 2:
                weight_ids.add(op["inputs"][1])

    def storage(tensor_id: str) -> str:
        return tensors[tensor_id].get("alias_of") or tensor_id

    store: dict[str, list[float]] = {storage(k): list(v) for k, v in values.items()}

    def fetch(tensor_id: str) -> list[float]:
        key = storage(tensor_id)
        if key not in store:
            raise FunctionalArtifactError(f"tensor {tensor_id} has no value yet")
        raw = store[key]
        if profile == "int8":
            if zp_i and tensor_id in data_input_ids:
                return [v - zp_i for v in raw]
            if zp_w and tensor_id in weight_ids:
                return [v - zp_w for v in raw]
        return raw

    def shape_of(tensor_id: str) -> list[int]:
        return [int(d) for d in tensors[tensor_id]["shape"]]

    goldens: dict[str, list[float]] = {}
    for op in executable["operations"]:
        kernel = _KERNELS.get(op["kind"])
        if kernel is None:
            raise FunctionalArtifactError(f"no reference kernel for {op['kind']}")
        if profile == "int8" and not op["mapping_required"]:
            raise FunctionalArtifactError(
                f"int8 semantics covers linear/conv operations only; {op['id']} is not mapped")
        out = kernel(op, fetch, shape_of)
        out = _apply_activation(out, op.get("activation", "linear"))
        if profile == "int8" and shift > 0:
            requantized = []
            for v in out:
                a = (int(round(v))*1 + (1 << (shift - 1))) >> shift
                requantized.append(float(max(requant_min, min(requant_max, a))))
            out = requantized
        elif profile in ("fp16", "bf16"):
            out = [_round_lowp(profile, v) for v in out]
        expected = _numel(shape_of(op["outputs"][0]))
        if len(out) != expected:
            raise FunctionalArtifactError(
                f"operation {op['id']} produced {len(out)} elements, expected {expected}")
        store[storage(op["outputs"][0])] = out
        goldens[op["id"]] = out
    return goldens


# --------------------------------------------------------------------------------------
# Value sources
# --------------------------------------------------------------------------------------

def synthesize_values(
    executable: Mapping[str, Any], seed: int,
    semantics: Mapping[str, Any] | None = None,
) -> dict[str, list[float]]:
    """Deterministic seeded values for every graph input/parameter (torch-free path).

    int8 semantics: integer operands (uint8-range when the matching zero point is set,
    signed int8 otherwise), small integer biases. fp16/bf16: values pre-rounded to the
    grid (mixed-precision contract: exact operands, fp32 accumulate, rounded output).
    """
    semantics = dict(semantics or {})
    profile = semantics.get("profile", "fp32")
    zp_i = int(semantics.get("input_zero_point", 0))
    zp_w = int(semantics.get("weight_zero_point", 0))
    rng = random.Random(seed)
    values: dict[str, list[float]] = {}
    inputs = set(executable["inputs"])
    weight_ids: set[str] = set()
    bias_ids: set[str] = set()
    for op in executable["operations"]:
        if op["kind"] in ("npusim.linear", "npusim.conv2d"):
            if len(op["inputs"]) >= 2:
                weight_ids.add(op["inputs"][1])
            if len(op["inputs"]) >= 3:
                bias_ids.add(op["inputs"][2])
    for tensor in executable["tensors"]:
        tensor_id, kind, count = tensor["id"], tensor["kind"], _numel(tensor["shape"])
        if tensor_id in inputs:
            if profile == "int8":
                low, high = (0, 255) if zp_i else (-127, 127)
                values[tensor_id] = [float(rng.randint(low, high)) for _ in range(count)]
            else:
                values[tensor_id] = [_round_lowp(profile, _f32(rng.uniform(-1.0, 1.0)))
                                     for _ in range(count)]
        elif kind in {"parameter", "buffer", "constant"}:
            if profile == "int8":
                if tensor_id in bias_ids:
                    values[tensor_id] = [float(rng.randint(-1000, 1000)) for _ in range(count)]
                else:
                    low, high = (0, 255) if zp_w else (-127, 127)
                    values[tensor_id] = [float(rng.randint(low, high)) for _ in range(count)]
            else:
                # Small weights keep deep-model activations in a well-conditioned range;
                # running_var-style buffers must stay positive (sqrt in BN).
                positive = "var" in tensor_id
                low, high = (0.5, 1.5) if positive else (-0.5, 0.5)
                values[tensor_id] = [_round_lowp(profile, _f32(rng.uniform(low, high)))
                                     for _ in range(count)]
    return values


def collect_torch_values(exported_program: Any, example_args: tuple) -> dict[str, list[float]]:
    """Pull real parameter/buffer/input values from a torch.export ExportedProgram."""
    signature = exported_program.graph_signature
    state = exported_program.state_dict
    qualified: dict[str, str] = {}
    qualified.update(getattr(signature, "inputs_to_parameters", {}))
    qualified.update(getattr(signature, "inputs_to_buffers", {}))
    constants = getattr(exported_program, "constants", {})
    lifted = dict(getattr(signature, "inputs_to_lifted_tensor_constants", {}))

    values: dict[str, list[float]] = {}
    for placeholder, name in qualified.items():
        values[placeholder] = [float(v) for v in state[name].detach().reshape(-1).tolist()]
    for placeholder, name in lifted.items():
        values[placeholder] = [float(v) for v in constants[name].detach().reshape(-1).tolist()]
    user_inputs = [
        spec.arg.name
        for spec in signature.input_specs
        if getattr(getattr(spec, "kind", None), "name", "") == "USER_INPUT"
    ]
    if len(user_inputs) != len(example_args):
        raise FunctionalArtifactError("example args do not match the exported graph inputs")
    for name, tensor in zip(user_inputs, example_args):
        values[name] = [float(v) for v in tensor.detach().reshape(-1).tolist()]
    return values


# --------------------------------------------------------------------------------------
# Artifact writer
# --------------------------------------------------------------------------------------

def _payload_name(identifier: str) -> str:
    return "".join(c if c.isalnum() or c in "_-" else "_" for c in identifier) + ".bin"


def write_artifact(
    executable: Mapping[str, Any],
    values: Mapping[str, list[float]],
    goldens: Mapping[str, list[float]],
    output_dir: str | Path,
    generator: Mapping[str, Any],
    semantics: Mapping[str, Any] | None = None,
) -> str:
    """Write payloads + goldens + manifest; returns the manifest path."""
    out = Path(output_dir)
    (out / "payload").mkdir(parents=True, exist_ok=True)
    (out / "golden").mkdir(parents=True, exist_ok=True)
    tensors_by_id = {t["id"]: t for t in executable["tensors"]}
    inputs = set(executable["inputs"])

    manifest_tensors = []
    for tensor_id, tensor_values in sorted(values.items()):
        declared = tensors_by_id[tensor_id]
        expected = _numel(declared["shape"])
        if len(tensor_values) != expected:
            raise FunctionalArtifactError(
                f"tensor {tensor_id} has {len(tensor_values)} values, shape needs {expected}")
        payload = _pack([_f32(v) for v in tensor_values])
        relative = "payload/" + _payload_name(tensor_id)
        (out / relative).write_bytes(payload)
        role = "input" if tensor_id in inputs else declared["kind"]
        manifest_tensors.append({
            "id": tensor_id,
            "role": role,
            "dtype": "float32",
            "shape": [int(d) for d in declared["shape"]],
            "layout": "contiguous",
            "payload": relative,
            "sha256": hashlib.sha256(payload).hexdigest(),
        })

    manifest_goldens = []
    for op in executable["operations"]:
        if op["id"] not in goldens:
            continue
        payload = _pack([_f32(v) for v in goldens[op["id"]]])
        relative = "golden/" + _payload_name(op["id"])
        (out / relative).write_bytes(payload)
        manifest_goldens.append({
            "operation_id": op["id"],
            "stage": "post_activation",
            "tensor_id": op["outputs"][0],
            "payload": relative,
            "sha256": hashlib.sha256(payload).hexdigest(),
        })

    manifest = {
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "executable_sha256": executable["executable_sha256"],
        "generator": dict(generator),
        "tensors": manifest_tensors,
        "golden": manifest_goldens,
    }
    if semantics and semantics.get("profile", "fp32") != "fp32":
        manifest["semantics"] = dict(semantics)
    manifest_path = out / "tensors.json"
    with manifest_path.open("w", encoding="utf-8") as target:
        json.dump(manifest, target, indent=2, sort_keys=True)
        target.write("\n")
    return str(manifest_path)


def synthesize_artifact(
    executable_path: str, output_dir: str, seed: int,
    semantics: Mapping[str, Any] | None = None,
) -> str:
    """Torch-free artifact: seeded deterministic values + reference-interpreter goldens."""
    from .executable_ir import load_executable_ir

    executable = load_executable_ir(executable_path)
    if "executable_sha256" not in executable:
        raise FunctionalArtifactError(
            "executable has no executable_sha256; regenerate it with the compile command")
    values = synthesize_values(executable, seed, semantics)
    goldens = interpret_executable(executable, values, semantics)
    return write_artifact(executable, values, goldens, output_dir, {
        "framework": "reference-interpreter",
        "seed": seed,
        "device": "cpu",
    }, semantics)


def export_functional_artifact(
    factory_spec: str, executable_output: str, artifact_output: str,
    model_name: str | None = None, tolerance: float = 1e-4,
) -> dict[str, str]:
    """Torch path: export -> lower -> executable + artifact with REAL model values.

    The reference interpreter's final outputs are cross-checked against the eager
    model's forward pass, so an interpreter/lowering drift fails here, not in CI.
    """
    import importlib

    torch = importlib.import_module("torch")
    from .executable_ir import dump_executable_ir, load_executable_ir
    from .export import export_model, load_callable
    from .lowering import lower_graph

    produced = load_callable(factory_spec)()
    model, example_args = produced[0], produced[1]
    graph = export_model(model, example_args, model_name=model_name)
    dump_executable_ir(lower_graph(graph), executable_output)
    executable = load_executable_ir(executable_output)

    exported = torch.export.export(model, args=example_args, strict=True)
    values = collect_torch_values(exported, example_args)
    goldens = interpret_executable(executable, values)

    # Cross-check: interpreter vs eager torch on the graph outputs.
    model.eval()
    with torch.no_grad():
        eager = model(*example_args)
    eager_outputs = eager if isinstance(eager, (tuple, list)) else (eager,)
    graph_outputs = executable["outputs"]
    tensors = {t["id"]: t for t in executable["tensors"]}

    def producer_of(tensor_id: str) -> str:
        storage = tensors[tensor_id].get("alias_of") or tensor_id
        for op in executable["operations"]:
            if op["outputs"][0] in (tensor_id, storage):
                return op["id"]
        raise FunctionalArtifactError(f"graph output {tensor_id} has no producer")

    for tensor_id, eager_tensor in zip(graph_outputs, eager_outputs):
        reference = goldens[producer_of(tensor_id)]
        flat = [float(v) for v in eager_tensor.detach().reshape(-1).tolist()]
        worst = max(abs(a - b) for a, b in zip(reference, flat))
        if worst > tolerance:
            raise FunctionalArtifactError(
                f"interpreter disagrees with eager torch on {tensor_id}: max |diff| {worst}")

    manifest = write_artifact(executable, values, goldens, artifact_output, {
        "framework": "pytorch",
        "framework_version": torch.__version__,
        "device": "cpu",
    })
    return {"executable": executable_output, "manifest": manifest}
