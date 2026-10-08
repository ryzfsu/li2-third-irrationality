#!/usr/bin/env python3
"""Convert parsed_tex.json (from main.tex only) into the input format of the v6 reference
implementation indep_caps.py, for a cross-check run. Lower-half Poisson centres by conjugation."""
import sys, json
import refmeasure
d, chains, rows, poi = refmeasure.load(sys.argv[1])
def iv(key):
    T, W = chains[key]['T'], chains[key]['W']
    return [{'a': str(T[i]), 'b': str(T[i + 1]), 'w': str(W[i])} for i in range(len(W))]
out = {'radius': '2/5', 'x': iv('x'), 'y': iv('v'), 'z_real': iv('z_int'),
       'z_poisson': [{'re': str(a), 'im': str(b), 'w': str(w)} for a, b, w in poi]}
json.dump(out, open(sys.argv[2], 'w'), indent=0)
