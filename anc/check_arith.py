#!/usr/bin/env python3
"""TeX-only R=8 arithmetic checker for Li2(1/3).

The manuscript is copied to a local temporary directory and parsed there.
The JSON files made by ``parse_tex.py`` and ``arith.py`` are transient
implementation outputs, never required inputs to the certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parent


def run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(cmd, cwd=cwd, check=True, text=True,
                          capture_output=True, env=env)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tex", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    started = time.monotonic()
    raw = args.tex.read_bytes()
    raw_sha = hashlib.sha256(raw).hexdigest()
    with tempfile.TemporaryDirectory(prefix="li3-arith-") as td:
        work = Path(td)
        tex = work / "input_main.tex"
        tex.write_bytes(raw)
        parsed = work / "parsed_tex.json"
        arithmetic = work / "arith_results.json"
        run([sys.executable, str(ROOT / "parse_tex.py"), str(tex), str(parsed)], work)
        run([sys.executable, str(ROOT / "arith.py"), str(parsed), str(arithmetic)], work)
        p = json.loads(parsed.read_text())
        a = json.loads(arithmetic.read_text())
    item4 = a["item4_arithmetic"]
    item5 = a["item5_alpha0"]
    if not (item4.get("pass") and item5.get("pass")):
        raise AssertionError("paper-only arithmetic replay did not pass")
    result = {
        "status": "PASS",
        "scope": "Exact R=8 arrangement, slab table, tau bounds and alpha assembled from the printed TeX; uniform valuation proof remains in the manuscript.",
        "raw_tex_sha256": raw_sha,
        "tex_sha256": p["tex_sha256"],
        "item4_arithmetic": item4,
        "item5_alpha0": item5,
        "seconds": round(time.monotonic() - started, 6),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
