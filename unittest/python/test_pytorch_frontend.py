"""Independent validation for the optional PyTorch graph frontend."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))

from frontend.pytorch.graph_ir import (  # noqa: E402
    GraphIRError,
    dump_graph_ir,
    graph_sha256,
    load_graph_ir,
    validate_graph_ir,
)


def minimal_graph() -> dict[str, object]:
    return {
        "schema_version": "npusim.graph.v1",
        "producer": {"name": "test", "version": "1"},
        "model": {"name": "linear", "structure_sha256": "test-model"},
        "inputs": ["input"],
        "outputs": ["output"],
        "tensors": [
            {"id": "input", "shape": [1, 4], "dtype": "float32", "layout": "contiguous", "kind": "input", "logical_bytes": 16},
            {"id": "weight", "shape": [4, 2], "dtype": "float32", "layout": "contiguous", "kind": "parameter", "logical_bytes": 32},
            {"id": "output", "shape": [1, 2], "dtype": "float32", "layout": "contiguous", "kind": "activation", "logical_bytes": 8},
        ],
        "nodes": [
            {"id": "linear", "op": "aten.linear.default", "inputs": ["input", "weight"], "outputs": ["output"], "attributes": {"fx_op": "call_function"}},
        ],
    }


class GraphIRValidationTest(unittest.TestCase):
    def test_valid_graph_has_stable_hash_and_round_trips(self) -> None:
        graph = minimal_graph()
        validate_graph_ir(graph)
        expected_hash = graph_sha256(graph)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "linear.graph.json"
            self.assertEqual(dump_graph_ir(graph, path), expected_hash)
            persisted = load_graph_ir(path)
        self.assertEqual(persisted["graph_sha256"], expected_hash)
        self.assertEqual(graph_sha256(persisted), expected_hash)

    def test_rejects_wrong_logical_bytes(self) -> None:
        graph = minimal_graph()
        graph["tensors"][0]["logical_bytes"] = 15  # type: ignore[index]
        with self.assertRaisesRegex(GraphIRError, "logical_bytes"):
            validate_graph_ir(graph)

    def test_rejects_non_topological_node(self) -> None:
        graph = minimal_graph()
        graph["nodes"][0]["inputs"] = ["output"]  # type: ignore[index]
        with self.assertRaisesRegex(GraphIRError, "topologically"):
            validate_graph_ir(graph)


class GraphIRCommandTest(unittest.TestCase):
    def test_cli_validates_fixture_without_torch(self) -> None:
        command = [
            sys.executable,
            "-m",
            "frontend.pytorch.cli",
            "validate",
            "frontend/fixtures/minimal_graph.json",
        ]
        completed = subprocess.run(
            command, cwd=REPOSITORY, text=True, capture_output=True, check=True
        )
        self.assertEqual(json.loads(completed.stdout)["status"], "valid")

    def test_cli_rejects_invalid_graph(self) -> None:
        graph = minimal_graph()
        graph["tensors"][0]["dtype"] = "float7"  # type: ignore[index]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.graph.json"
            path.write_text(json.dumps(graph), encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, "-m", "frontend.pytorch.cli", "validate", str(path)],
                cwd=REPOSITORY,
                text=True,
                capture_output=True,
            )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("unsupported IR dtype", completed.stderr)


if __name__ == "__main__":
    unittest.main()
