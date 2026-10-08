#!/usr/bin/env python3
"""TeX-only analytic certificate checker for Li2(1/3).

The finite data are read from the printed tables in ``main.tex``.  All JSON
serializations and leaf outputs are made in a temporary directory.  The main
implementation is the paper-only chord-plus-Taylor verifier in ``caps.py``;
``indep_caps_params.py`` is run as a second directed-Arb implementation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parent
DEPTHS = {"x": (132, 9, 67584), "v": (145, 8, 37120),
          "z": (106, 7, 13568), "circle": (32, 8, 8192)}
CAPS = {"x": F(-120787, 25000), "v": F(-237159, 50000),
        "z": F(2159, 50000), "circle": F(4283, 100000)}


def run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(cmd, cwd=cwd, check=True, text=True,
                          capture_output=True, env=env)


def check_printed_tables(parsed: dict, caps: dict) -> None:
    rows = parsed["caps_table"]
    if len(rows) != 4:
        raise AssertionError("expected four cap rows")
    for row, key in zip(rows, ("x", "v", "z", "circle")):
        expected = DEPTHS[key][2]
        if int(row["leaves"]) != expected or F(row["cap"]) != CAPS[key]:
            raise AssertionError(f"printed cap row mismatch: {key}")
        if F(row["upper_bound"]) > CAPS[key]:
            raise AssertionError(f"printed cap endpoint exceeds cap: {key}")
    depths = parsed["depths_table"]
    if len(depths) < 4:
        raise AssertionError("expected four depth rows")
    for row, key in zip(depths[:4], ("x", "v", "z", "circle")):
        got = (int(row["pieces"]), int(row["depth"]), int(row["leaves"]))
        if got != DEPTHS[key]:
            raise AssertionError(f"printed depth row mismatch: {key}: {got}")
    if sum(int(row["leaves"]) for row in depths[:4]) != 126464:
        raise AssertionError("printed total leaf count mismatch")
    for key, rec in caps["supports"].items():
        depth_key = "circle" if key == "circ" else key
        if not rec["pass"] or rec["failing_leaves"] != 0:
            raise AssertionError(f"fast caps failed: {key}")
        if rec["leaves"] != DEPTHS[depth_key][2]:
            raise AssertionError(f"fast leaf count mismatch: {key}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tex", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--replay", action="store_true",
                    help="record that the uniform-depth replay was requested")
    args = ap.parse_args()
    started = time.monotonic()
    raw = args.tex.read_bytes()
    raw_sha = hashlib.sha256(raw).hexdigest()
    with tempfile.TemporaryDirectory(prefix="li3-paperonly-") as td:
        work = Path(td)
        tex = work / "input_main.tex"
        tex.write_bytes(raw)
        parsed_path = work / "parsed_tex.json"
        energy_path = work / "energy_results.json"
        caps_path = work / "caps_results.json"
        ref_path = work / "reference.json"
        indep_path = work / "independent_caps.json"
        parse_run = run([sys.executable, str(ROOT / "parse_tex.py"), str(tex), str(parsed_path)], work)
        energy_run = run([sys.executable, str(ROOT / "energy.py"), str(parsed_path), str(energy_path)], work)
        caps_run = run([sys.executable, str(ROOT / "caps.py"), str(parsed_path), str(caps_path),
                        "x,v,z,circ"], work)
        run([sys.executable, str(ROOT / "tex2refjson.py"), str(parsed_path), str(ref_path)], work)
        indep_run = run([sys.executable, str(ROOT / "indep_caps_params.py"), str(ref_path),
                         str(indep_path)], work)
        parsed = json.loads(parsed_path.read_text())
        energy = json.loads(energy_path.read_text())
        caps = json.loads(caps_path.read_text())
        reference = json.loads(ref_path.read_text())
        independent = json.loads(indep_path.read_text())
        if parsed["tex_sha256"] != raw_sha:
            raise AssertionError("parser hash mismatch")
        if not energy["item2_energy"]["pass"]:
            raise AssertionError("energy replay failed")
        if not caps["pass"]:
            raise AssertionError("paper-only caps replay failed")
        if independent.get("status") != "PASS":
            raise AssertionError("independent caps replay failed")
        check_printed_tables(parsed, caps)
        # Independent implementation must agree on each support's accepted count.
        for key, rec in independent["supports"].items():
            depth_key = {"z_real": "z", "z_circle": "circle"}.get(key, key)
            if not rec["pass"] or rec["accepted_leaves"] != DEPTHS[depth_key][2]:
                raise AssertionError(f"independent leaf mismatch: {key}")
        item = {
            "energy": energy,
            "caps": caps,
            "independent_caps": independent,
            "reference_counts": {k: len(reference[k]) for k in ("x", "y", "z_real", "z_poisson")},
            "subprocess_seconds": {
                "parse": None,
                "energy": energy["item2_energy"].get("seconds"),
                "caps": caps.get("seconds"),
                "independent": independent.get("seconds"),
            },
        }
    result = {
        "status": "PASS",
        "scope": "TeX-derived energy and uniform-depth chord-plus-Taylor caps, with an independent directed-Arb caps cross-check; written uniform arguments remain in the manuscript.",
        "raw_tex_sha256": raw_sha,
        "tex_sha256": parsed["tex_sha256"],
        "replay_requested": bool(args.replay),
        "reference_counts": item["reference_counts"],
        "energy": item["energy"],
        "caps": item["caps"],
        "independent_caps": item["independent_caps"],
        "seconds": round(time.monotonic() - started, 6),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
