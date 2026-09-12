"""Command line interface for the NPUsim PyTorch frontend."""

from __future__ import annotations

import argparse
import json
from typing import Any

from .executable_ir import (
    ExecutableIRError,
    dump_executable_ir,
    executable_sha256,
    load_executable_ir,
)
from .export import export_to_file, load_callable
from .functional_artifact import (
    FunctionalArtifactError,
    export_functional_artifact,
    synthesize_artifact,
)
from .graph_ir import GraphIRError, graph_sha256, load_graph_ir
from .lowering import LoweringError, lower_graph


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


def _compile_command(arguments: argparse.Namespace) -> int:
    executable_hash = dump_executable_ir(lower_graph(load_graph_ir(arguments.graph)), arguments.output)
    print(json.dumps({"executable_sha256": executable_hash, "output": arguments.output}, sort_keys=True))
    return 0


def _validate_executable_command(path: str) -> int:
    executable = load_executable_ir(path)
    print(json.dumps({"executable_sha256": executable_sha256(executable), "status": "valid"}, sort_keys=True))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="NPUsim PyTorch graph frontend")
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate", help="validate an NPUsim capture graph artifact")
    validate.add_argument("graph", help="path to capture graph JSON")
    export = commands.add_parser("export", help="export a torch module factory to capture graph JSON")
    export.add_argument("--factory", required=True, help="module:callable returning model and example inputs")
    export.add_argument("--output", required=True, help="destination graph JSON")
    export.add_argument("--model-name", default=None, help="stable model name for result provenance")
    compile_command = commands.add_parser(
        "compile", help="lower a captured graph into framework-neutral executable IR"
    )
    compile_command.add_argument("--graph", required=True, help="validated capture graph JSON")
    compile_command.add_argument("--output", required=True, help="destination executable JSON")
    validate_executable = commands.add_parser(
        "validate-executable", help="validate an NPUsim executable IR artifact"
    )
    validate_executable.add_argument("executable", help="path to executable JSON")
    synthesize = commands.add_parser(
        "synthesize-functional",
        help="write an npusim.tensor.v1 artifact with seeded deterministic values and"
             " reference-interpreter goldens (no torch required)",
    )
    synthesize.add_argument("--executable", required=True, help="npusim.exec.v1 JSON")
    synthesize.add_argument("--output-dir", required=True, help="artifact directory")
    synthesize.add_argument("--seed", type=int, default=7, help="deterministic value seed")
    synthesize.add_argument("--profile", default="fp32",
                            choices=["fp32", "int8", "fp16", "bf16"],
                            help="arithmetic semantics of the artifact (int8: integer"
                                 " operands + optional requant/zero points; fp16/bf16:"
                                 " grid-rounded operands and outputs)")
    synthesize.add_argument("--requant-shift", type=int, default=0,
                            help="int8: rounding right-shift applied after activation")
    synthesize.add_argument("--input-zero-point", type=int, default=0,
                            help="int8: asymmetric input zero point (uint8 operands)")
    synthesize.add_argument("--weight-zero-point", type=int, default=0,
                            help="int8: asymmetric weight zero point (uint8 weights)")
    export_functional = commands.add_parser(
        "export-functional",
        help="export a torch module factory to an executable + npusim.tensor.v1 artifact"
             " with the model's REAL parameter/input values (requires torch)",
    )
    export_functional.add_argument("--factory", required=True,
                                   help="module:callable returning model and example inputs")
    export_functional.add_argument("--executable-output", required=True,
                                   help="destination executable JSON")
    export_functional.add_argument("--output-dir", required=True, help="artifact directory")
    export_functional.add_argument("--model-name", default=None,
                                   help="stable model name for result provenance")
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "validate":
            return _validate_command(arguments.graph)
        if arguments.command == "compile":
            return _compile_command(arguments)
        if arguments.command == "validate-executable":
            return _validate_executable_command(arguments.executable)
        if arguments.command == "synthesize-functional":
            semantics = {"profile": arguments.profile}
            if arguments.profile == "int8":
                semantics.update({
                    "requant_shift": arguments.requant_shift,
                    "requant_min": -127,
                    "requant_max": 127,
                    "input_zero_point": arguments.input_zero_point,
                    "weight_zero_point": arguments.weight_zero_point,
                })
            manifest = synthesize_artifact(arguments.executable, arguments.output_dir,
                                           arguments.seed, semantics)
            print(json.dumps({"manifest": manifest, "profile": arguments.profile,
                              "seed": arguments.seed}, sort_keys=True))
            return 0
        if arguments.command == "export-functional":
            result = export_functional_artifact(
                arguments.factory, arguments.executable_output, arguments.output_dir,
                model_name=arguments.model_name)
            print(json.dumps(result, sort_keys=True))
            return 0
        return _export_command(arguments)
    except (
        ExecutableIRError,
        FunctionalArtifactError,
        GraphIRError,
        LoweringError,
        RuntimeError,
        TypeError,
        ValueError,
        OSError,
    ) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
