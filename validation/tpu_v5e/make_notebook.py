#!/usr/bin/env python3
"""<script>.py -> <script>.ipynb for kaggle_tpu_bench.py and tpu_v5e1_stream_bench.py (or the scripts given); cells split
on '# %%', '# %% [markdown]' cells as markdown."""
import json
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
for script in sys.argv[1:] or ("kaggle_tpu_bench.py", "tpu_v5e1_stream_bench.py", "tpu_v5e1_hbm_probe.py"):
    src = open(os.path.join(here, script)).read()
    cells = []
    for chunk in re.split(r"(?m)^# %%", src)[1:]:
        head, _, body = chunk.partition("\n")
        body = body.strip("\n")
        if head.strip() == "[markdown]":
            text = "\n".join(l[2:] if l.startswith("# ") else l.lstrip("#") for l in body.splitlines())
            cells.append({"cell_type": "markdown", "metadata": {}, "source": text.splitlines(True)})
        else:
            cells.append({"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
                          "source": body.splitlines(True)})
    meta = {"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
            "language_info": {"name": "python"}}
    if script.startswith("tpu_v5e1"):          # Colab notebook
        meta["accelerator"] = "TPU"
    with open(os.path.join(here, script.replace(".py", ".ipynb")), "w") as f:
        json.dump({"cells": cells, "nbformat": 4, "nbformat_minor": 5, "metadata": meta}, f, indent=1)
    print(script, len(cells), "cells")
