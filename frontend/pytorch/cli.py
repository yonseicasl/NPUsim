"""Command line interface for the NPUsim PyTorch frontend."""

from __future__ import annotations

import argparse
import json
from typing import Any

from .export import export_to_file, load_callable
from .graph_ir import GraphIRError, graph_sha256, load_graph_ir


def _build_export(factory_spec: str) -> tuple[Any, tuple[Any, ...], dict[str, Any]]:
    produced = load_callable(factory_spec)()
    if not isinstance(produced, tuple) or len(produced) not in {2, 3}:
        raise TypeError(
            "factory must return (torch.nn.Module, example_args) or "
            "(torch.nn.Module, example_args, example_kwargs)"
        )
    model, example_args = produced[0], produced[1]
    example_kwargs = produced[2] if len(produced) == 3 else {}
    if not isinstance(example_args, tuple) or not isinstance(example_kwargs, dict):
        raise TypeError("factory example_args must be a tuple and example_kwargs must be a dict")
    return model, example_args, example_kwargs


def _validate_command(path: str) -> int:
    graph = load_graph_ir(path)
    print(json.dumps({"graph_sha256": graph_sha256(graph), "status": "valid"}, sort_keys=True))
    return 0


def _export_command(arguments: argparse.Namespace) -> int:
    model, example_args, example_kwargs = _build_export(arguments.factory)
    graph_hash = export_to_file(
        model,
        example_args,
        arguments.output,
        example_kwargs,
        model_name=arguments.model_name,
    )
    print(json.dumps({"graph_sha256": graph_hash, "output": arguments.output}, sort_keys=True))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="NPUsim PyTorch graph frontend")
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate", help="validate an NPUsim graph IR artifact")
    validate.add_argument("graph", help="path to graph JSON")
    export = commands.add_parser("export", help="export a torch module factory to graph JSON")
    export.add_argument("--factory", required=True, help="module:callable returning model and example inputs")
    export.add_argument("--output", required=True, help="destination graph JSON")
    export.add_argument("--model-name", default=None, help="stable model name for result provenance")
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "validate":
            return _validate_command(arguments.graph)
        return _export_command(arguments)
    except (GraphIRError, RuntimeError, TypeError, ValueError, OSError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
