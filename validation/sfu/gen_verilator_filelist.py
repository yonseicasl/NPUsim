#!/usr/bin/env python3
"""Generate verif/verilator/verilator_<project>.f for an NVDLA project outdir.

Mirrors the structure of the shipped verilator_nv_small.f (plus the 2026-08-20 patch's
lesson): -I directories double as Verilator library search paths for modules whose file
name matches the module name, and every RAM model and generated FIFO library file is
registered with -v EXPLICITLY, because their module names do not match the file names
(-y/-I search cannot find them). Registering ALL of vmod/rams/model and vmod/fifos up
front avoids the 67-entry whack-a-mole the nv_small bring-up went through.

Usage: gen_verilator_filelist.py <project>       (run from ext/nvdla)
Writes verif/verilator/verilator_<project>.f
"""
import glob
import os
import sys


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    project = sys.argv[1]
    outdir = f"outdir/{project}/vmod"
    if not os.path.isdir(outdir):
        raise SystemExit(f"{outdir} does not exist; build vmod first")
    rel = f"../../outdir/{project}/vmod"
    lines = []
    for sub in sorted(os.listdir(f"{outdir}/nvdla")):
        if os.path.isdir(f"{outdir}/nvdla/{sub}"):
            lines.append(f"-I{rel}/nvdla/{sub}")
    for extra in ("rams/synth", "vlibs", "include"):
        lines.append(f"-I{rel}/{extra}")
    lines.append(f"-v {rel}/vlibs/RANDFUNC.vlib")
    lines.append(f"-v {rel}/vlibs/nv_assert_no_x.vlib")
    if os.path.exists(f"{outdir}/nvdla/nocif/NV_NVDLA_XXIF_libs.v"):
        lines.append(f"-v {rel}/nvdla/nocif/NV_NVDLA_XXIF_libs.v")
    for pattern, label in ((f"{outdir}/rams/model/*.v", "RAM models"),
                           (f"{outdir}/fifos/*.v", "generated FIFO libraries")):
        files = sorted(glob.glob(pattern))
        if not files:
            raise SystemExit(f"no files match {pattern}")
        lines.append(f"// --- all {label} ({len(files)}) ---")
        for path in files:
            lines.append(f"-v {rel}/{os.path.relpath(path, outdir)}")
    lines += [
        "",
        "-DNO_PLI_OR_EMU",
        "-DNO_PLI",
        "-DDESIGNWARE_NOEXIST",
        "-DSYNTHESIS",
        "-Wno-moddup",
        "-Wno-fatal",
        "--top-module NV_nvdla",
        f"{rel}/nvdla/top/NV_nvdla.v",
    ]
    target = f"verif/verilator/verilator_{project}.f"
    with open(target, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {target}: {len(lines)} lines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
