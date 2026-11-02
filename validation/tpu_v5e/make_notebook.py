#!/usr/bin/env python3
"""kaggle_tpu_bench.py -> kaggle_tpu_bench.ipynb (cells split on '# %%', '# %% [markdown]' cells as markdown)."""
import json
import os
import re

here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "kaggle_tpu_bench.py")).read()
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
nb = {"cells": cells, "nbformat": 4, "nbformat_minor": 5,
      "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
                   "language_info": {"name": "python"}}}
with open(os.path.join(here, "kaggle_tpu_bench.ipynb"), "w") as f:
    json.dump(nb, f, indent=1)
print(len(cells), "cells")
