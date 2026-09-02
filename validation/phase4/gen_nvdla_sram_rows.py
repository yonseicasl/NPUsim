#!/usr/bin/env python3
"""Generate the nv_small SRAM rows of cacti_sram_reference.csv.

Same method as the original phase-4 rows (PROVENANCE.md): CACTI 7 commit 1ffd8dfb,
the repo's own cache.cfg with ONLY size/block/tech/bus/type/associativity overridden,
itrs-hp cells at 360 K as shipped. Run with the CACTI checkout path as argv[1];
appends nothing itself -- prints CSV rows for a human-reviewed paste.

Geometries (from the generated nv_small RTL, not from NPUsim's config):
  nvdla_cbuf_bank      2048 B, 8 B block, 64 b bus   -- nv_ram_rws_256x64, the CBUF
                       bank macro (32 banks x 2 macros = 128 KiB total)
  nvdla_cacc_assembly   512 B, 32 B block, 128 b bus -- nv_ram_rws_16x256 (16 x 256 b).
                       CACTI finds no valid organization for the physical 256 b bus on
                       an array this small; the 128 b-bus row is used and consumed
                       PER BYTE (accumulator_spill_energy is a per-byte key), which the
                       narrower bus leaves unchanged to first order. Stated, not hidden.
"""
import re, subprocess, os, sys

CACTI = sys.argv[1] if len(sys.argv) > 1 else "cacti"
GEOMS = [
    ("nvdla_cbuf_bank",     "separate", "read/write_energy",          2048, 8,  64),
    ("nvdla_cacc_assembly", "spatial_arch", "accumulator_spill_energy", 512, 32, 128),
]
base = open(os.path.join(CACTI, "cache.cfg")).read()

def override(txt, key, val):
    out, done = [], False
    for line in txt.split("\n"):
        if line.strip().startswith(f"-{key} "):
            if not done:
                out.append(f"-{key} {val}"); done = True
        else:
            out.append(line)
    assert done, key
    return "\n".join(out)

for name, section, key, size, block, bus in GEOMS:
    for tech in (0.045, 0.032, 0.022):
        cfg = base
        for k, v in [("size (bytes)", size), ("block size (bytes)", block),
                     ("technology (u)", tech), ("output/input bus width", bus),
                     ("cache type", '"ram"'), ("associativity", 1)]:
            cfg = override(cfg, k, v)
        path = os.path.join(CACTI, "gen_row.cfg")
        open(path, "w").write(cfg)
        r = subprocess.run(["./cacti", "-infile", "gen_row.cfg"], cwd=CACTI,
                           capture_output=True, text=True)
        def grab(pat):
            m = re.search(pat, r.stdout)
            return float(m.group(1)) if m else None
        rd = grab(r"Total dynamic read energy per access \(nJ\): ([\d.]+)")
        wr = grab(r"Total dynamic write energy per access \(nJ\): ([\d.]+)")
        lk = grab(r"Total leakage power of a bank.*\(mW\): ([\d.]+)")
        ac = grab(r"Access time \(ns\): ([\d.]+)")
        assert rd is not None, (name, tech)
        print(f"{name},{section},{key},{size},{block},{bus},{int(tech*1000)},"
              f"{rd*1000:.4f},{wr*1000:.4f},{lk:.4f},{ac:.4f}")
