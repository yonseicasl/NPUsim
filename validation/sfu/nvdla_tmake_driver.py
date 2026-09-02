#!/usr/bin/env python3
"""Minimal faithful replacement for ext/nvdla/tools/bin/tmake (whose perl deps -- YAML,
Capture::Tiny -- are not installed on this machine). Reproduces exactly what tmake does:
topologically walk tools/etc/build.config dependencies for the requested targets and run
`cd <sandbox>; make PROJECT=<project>` per sandbox, honoring per-target `command` and
`skip` entries. Nothing here deletes outdir (tmake's clean logic is intentionally NOT
reproduced -- the nv_small golden artifacts live there).

Usage: nvdla_tmake_driver.py <project> <target> [<target>...]
Run from ext/nvdla.
"""
import subprocess
import sys


def parse_build_config(path):
    """Dependency-free parser for build.config's flat 2-level YAML shape:
    top-level `target:` blocks containing `sandbox:`/`dependencies:`/`skip:` lists
    (`- item` entries) and scalar `command:`/`desc:` values."""
    tree = {}
    target = None
    field = None
    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            indent = len(line) - len(line.lstrip(" "))
            if indent == 0 and stripped.endswith(":"):
                target = stripped[:-1]
                tree[target] = {}
                field = None
            elif target is not None and stripped.startswith("- "):
                if field is None:
                    raise SystemExit(f"[driver] list item outside a field: {line}")
                tree[target].setdefault(field, []).append(stripped[2:].strip())
            elif target is not None and stripped.endswith(":"):
                field = stripped[:-1]
            elif target is not None and ":" in stripped:
                key, value = stripped.split(":", 1)
                tree[target][key.strip()] = value.strip()
                field = None
    return tree


def main():
    argv = list(sys.argv[1:])
    assume_done = set()
    if "--assume-done" in argv:
        i = argv.index("--assume-done")
        assume_done = set(argv[i + 1].split(","))
        del argv[i:i + 2]
    if len(argv) < 2:
        print(__doc__)
        return 2
    project = argv[0]
    targets = argv[1:]
    tree = parse_build_config("tools/etc/build.config")

    done = set(assume_done)
    if assume_done:
        print(f"[driver] assuming already built: {', '.join(sorted(assume_done))}")

    def build(key):
        if key in done:
            return
        node = tree.get(key)
        if node is None:
            raise SystemExit(f"[driver] unknown build target '{key}'")
        for dep in node.get("dependencies", []) or []:
            build(dep)
        if project in (node.get("skip") or []):
            print(f"[driver] SKIP {key} for {project} (build.config skip list)")
            done.add(key)
            return
        for sandbox in node.get("sandbox", []) or []:
            if "command" in node:
                cmd = node["command"].replace("<project>", project)
            else:
                cmd = f"make PROJECT={project}"
            print(f"[driver] {key}: (cd {sandbox}; {cmd})", flush=True)
            result = subprocess.run(cmd, shell=True, cwd=sandbox,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True)
            tail = "\n".join(result.stdout.splitlines()[-15:])
            if result.returncode != 0:
                print(tail)
                raise SystemExit(f"[driver] FAIL in {sandbox} for target {key}")
        done.add(key)

    for target in targets:
        build(target)
    print(f"[driver] PASS: {' '.join(targets)} for {project}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
