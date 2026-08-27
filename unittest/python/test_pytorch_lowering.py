"""Executable IR and strict PyTorch lowering regression tests."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))

from frontend.pytorch.executable_ir import (  # noqa: E402
    ExecutableIRError,
    dump_executable_ir,
    load_executable_ir,
    validate_executable_ir,
)
from frontend.pytorch.graph_ir import load_graph_ir  # noqa: E402
from frontend.pytorch.lowering import LoweringError, lower_graph  # noqa: E402


class PyTorchLoweringTest(unittest.TestCase):
    def setUp(self) -> None:
        self.graph = load_graph_ir(REPOSITORY / "frontend/fixtures/linear_relu_graph.json")

    def test_linear_relu_fuses_with_exact_geometry(self) -> None:
        executable = lower_graph(self.graph)
        self.assertEqual(len(executable["operations"]), 1)
        operation = executable["operations"][0]
        self.assertEqual(operation["id"], "linear")
        self.assertEqual(operation["kind"], "npusim.linear")
        self.assertEqual(operation["activation"], "relu")
        self.assertEqual(operation["outputs"], ["relu_output"])
        self.assertEqual(operation["source_nodes"], ["linear", "relu"])
        self.assertEqual(
            operation["geometry"],
            {"batch": 64, "input_features": 64, "output_features": 64},
        )
        self.assertEqual(executable["coverage"]["captured_nodes"], 2)
        self.assertEqual(executable["coverage"]["lowered_source_nodes"], 2)

    def test_executable_round_trip_checks_hash(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "linear_relu.exec.json"
            dump_executable_ir(lower_graph(self.graph), path)
            loaded = load_executable_ir(path)
            self.assertEqual(loaded["operations"][0]["activation"], "relu")
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["operations"][0]["geometry"]["batch"] = 32
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ExecutableIRError, "executable_sha256"):
                load_executable_ir(path)

    def test_standalone_activation_is_not_silently_skipped(self) -> None:
        graph = copy.deepcopy(self.graph)
        graph["nodes"] = [graph["nodes"][1]]
        graph["inputs"] = ["linear_output"]
        graph["tensors"] = [
            tensor for tensor in graph["tensors"] if tensor["id"] in {"linear_output", "relu_output"}
        ]
        with self.assertRaisesRegex(LoweringError, "standalone activation"):
            lower_graph(graph)

    def test_unknown_operation_fails(self) -> None:
        graph = copy.deepcopy(self.graph)
        graph["nodes"][0]["op"] = "aten.magic.default"
        with self.assertRaisesRegex(LoweringError, "unsupported PyTorch operation"):
            lower_graph(graph)

    def test_conv2d_attributes_lower_to_geometry(self) -> None:
        graph = {
            "schema_version": "npusim.graph.v1",
            "producer": {"name": "test", "version": "1"},
            "model": {"name": "conv", "structure_sha256": "conv"},
            "inputs": ["x"],
            "outputs": ["y"],
            "tensors": [
                {"id": "x", "shape": [1, 4, 8, 8], "dtype": "float32", "layout": "contiguous", "kind": "input", "logical_bytes": 1024},
                {"id": "w", "shape": [6, 2, 3, 3], "dtype": "float32", "layout": "contiguous", "kind": "parameter", "logical_bytes": 432},
                {"id": "y", "shape": [1, 6, 4, 4], "dtype": "float32", "layout": "contiguous", "kind": "activation", "logical_bytes": 384},
            ],
            "nodes": [
                {
                    "id": "conv",
                    "op": "aten.conv2d.default",
                    "inputs": ["x", "w"],
                    "outputs": ["y"],
                    "attributes": {"stride": [2, 2], "padding": [1, 1], "dilation": [1, 1], "groups": 2},
                }
            ],
        }
        operation = lower_graph(graph)["operations"][0]
        self.assertEqual(operation["geometry"]["groups"], 2)
        self.assertEqual(operation["geometry"]["output_height"], 4)
        self.assertEqual(operation["geometry"]["output_width"], 4)

        # aten.convolution has transposed/output_padding before groups, unlike
        # aten.conv2d. Preserve that schema distinction during positional lowering.
        convolution_graph = copy.deepcopy(graph)
        convolution_node = convolution_graph["nodes"][0]
        convolution_node["op"] = "aten.convolution.default"
        convolution_node.pop("attributes")
        convolution_node["arguments"] = [
            {"tensor": "x"},
            {"tensor": "w"},
            None,
            [2, 2],
            [1, 1],
            [1, 1],
            False,
            [0, 0],
            2,
        ]
        convolution = lower_graph(convolution_graph)["operations"][0]
        self.assertEqual(convolution["geometry"]["groups"], 2)
        convolution_node["arguments"][6] = True
        with self.assertRaisesRegex(LoweringError, "transposed convolution"):
            lower_graph(convolution_graph)

    def test_validator_rejects_symbolic_executable_shape(self) -> None:
        executable = lower_graph(self.graph)
        executable["tensors"][0]["shape"][0] = "batch"
        with self.assertRaisesRegex(ExecutableIRError, "positive integer"):
            validate_executable_ir(executable)


class PyTorchCompileCommandTest(unittest.TestCase):
    def test_cli_compiles_without_torch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "linear_relu.exec.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "frontend.pytorch.cli",
                    "compile",
                    "--graph",
                    "frontend/fixtures/linear_relu_graph.json",
                    "--output",
                    str(output),
                ],
                cwd=REPOSITORY,
                text=True,
                capture_output=True,
                check=True,
            )
            result = json.loads(completed.stdout)
            self.assertEqual(result["output"], str(output))
            self.assertTrue(result["executable_sha256"])
            self.assertEqual(load_executable_ir(output)["operations"][0]["id"], "linear")


if __name__ == "__main__":
    unittest.main()
