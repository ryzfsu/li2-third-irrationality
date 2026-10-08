#!/usr/bin/env python3
"""Parse the printed data of the Li2(1/3) manuscript into exact rationals.

Sole data source: the TeX file given as argv[1]. No JSON/anc input is read.
Output: argv[2] (JSON), containing every parsed table plus the theorem-level
macros, all as exact fraction strings.
"""
import sys, re, json, hashlib
from fractions import Fraction as Fr

TEX = sys.argv[1]
OUT = sys.argv[2]
raw = open(TEX, 'rb').read()
src = raw.decode('utf-8')


def F(s):
    s = s.strip().strip('$').strip()
    return Fr(s)


def block(label, end=r'\end{longtable}'):
    i = src.index(r'\label{%s}' % label)
    j = src.index(end, i)
    return src[i:j]


def macro(name):
    m = re.search(r'\\newcommand\{\\%s\}\{(.*)\}\s*$' % name, src, re.M)
    return m.group(1)


def fracmacro(name):
    s = macro(name)
    m = re.fullmatch(r'\\frac\{(\d+)\}\{(\d+)\}', s)
    return Fr(int(m.group(1)), int(m.group(2))) if m else Fr(s)


out = {'tex_sha256': hashlib.sha256(raw).hexdigest()}
macros = {}
for nm in ['Rseam', 'nslabs', 'innervalue', 'Avalue', 'taubound', 'qminus', 'Kx', 'Kv', 'Kz',
           'Kzcirc', 'alphazero', 'alphabound', 'marginvalue', 'Narcs', 'nleaves', 'NE', 'NF',
           'precE', 'precF', 'vzero']:
    macros[nm] = str(fracmacro(nm))
macros['tauvalue_decimal'] = macro('tauvalue')
macros['capmargin'] = '2e-6'
out['macros'] = macros

# ---------------- interval chains ----------------
b = block('tab:reference-intervals')
chains = {}
cur = None
for line in b.splitlines():
    m = re.search(r'\\multicolumn\{7\}\{@\{\}l\}\{\$\\nu_(\w|z\^\{\\rm int\})\$: \$b_0=([^$]+)\$', line)
    if m:
        key = {'x': 'x', 'v': 'v'}.get(m.group(1), 'z_int')
        cur = key
        chains[cur] = {'b0': F(m.group(2)), 'rows': []}
        continue
    if cur and re.match(r'^\s*\d+\s*&', line):
        cells = [c.strip() for c in line.rstrip().rstrip('\\').split('&')]
        k = int(cells[0])
        vals = [c for c in cells[1:]]
        for t in range(3):
            bk, wk = vals[2 * t], vals[2 * t + 1]
            if bk.strip() == '':
                continue
            chains[cur]['rows'].append((k + t, F(bk), F(wk)))
out['chains'] = {}
for key, ch in chains.items():
    out['chains'][key] = {'b0': str(ch['b0']),
                          'rows': [[k, str(bk), str(wk)] for k, bk, wk in ch['rows']]}

# ---------------- Poisson table ----------------
b = block('tab:reference-poisson')
prow = []
body = b.split(r'\endhead', 1)[1]
for line in body.splitlines():
    if not re.match(r'^\s*\d+\s*&', line):
        continue
    cells = [c.strip() for c in line.rstrip().rstrip('\\').split('&')]
    for t in range(0, len(cells), 4):
        if t + 3 < len(cells) and cells[t] != '':
            prow.append((int(cells[t]), F(cells[t + 1]), F(cells[t + 2]), F(cells[t + 3])))
prow.sort()
out['poisson_rows'] = [[k, str(a), str(bb), str(w)] for k, a, bb, w in prow]

# ---------------- slab table ----------------
b = block('tab:slabs')
body = b.split(r'\endhead', 1)[1]
slabcols = [[], [], []]
for line in body.splitlines():
    if '&' not in line:
        continue
    cells = [c.strip() for c in line.rstrip().rstrip('\\').split('&')]
    for t in range(3):
        g = cells[4 * t:4 * t + 4]
        if len(g) == 4 and g[0].strip() != '':
            slabcols[t].append([str(F(x)) for x in g])
out['slabs'] = slabcols[0] + slabcols[1] + slabcols[2]

# ---------------- energy table ----------------
b = src[src.index(r'\label{tab:energy}'):src.index(r'\end{tabular}', src.index(r'\label{tab:energy}'))]
energy = {}
names = ['I_xx', 'I_vv', 'I_int', 'I_mix', 'I_Poi', 'I_vzint', 'I_vzPoi', 'D_xv', 'Q']
vals = re.findall(r'\\texttt\{\[([-0-9.]+),([-0-9.]+)\]\}', b)
assert len(vals) == 9
for nm, (lo, hi) in zip(names, vals):
    energy[nm] = [lo, hi]
out['energy_table'] = energy

# ---------------- caps table ----------------
b = src[src.index(r'\label{tab:caps}'):src.index(r'\end{tabular}', src.index(r'\label{tab:caps}'))]
caps = []
for line in b.splitlines():
    m = re.match(r'^(.*?)&\s*\$(\d+)\$\s*&\s*\$([-0-9.]+)\$\s*&\s*\$([-0-9.]+)\$', line)
    if m:
        caps.append({'support': m.group(1).strip(), 'leaves': int(m.group(2)),
                     'upper_bound': m.group(3), 'cap': m.group(4)})
out['caps_table'] = caps

# ---------------- depths table ----------------
b = src[src.index(r'\label{tab:depths}'):src.index(r'\end{tabular}', src.index(r'\label{tab:depths}'))]
dep = []
for line in b.splitlines():
    m = re.match(r'^(.*?)&\s*\$(\d+)\$\s*&\s*\$(\d+)\$\s*&\s*\$(\d+)\$', line)
    if m:
        dep.append({'support': m.group(1).strip(), 'pieces': int(m.group(2)), 'depth': int(m.group(3)),
                    'leaves': int(m.group(4))})
out['depths_table'] = dep

# ---------------- G4 coefficients (parsed only for completeness) ----------------
b = src[src.index(r'\label{tab:g4-coefficients}'):src.index(r'\end{tabular}', src.index(r'\label{tab:g4-coefficients}'))]
out['g4_q'] = re.findall(r'\$(-?\d{10,})\$', b)

json.dump(out, open(OUT, 'w'), indent=1)
print('parsed: chains', {k: len(v['rows']) for k, v in chains.items()}, 'poisson rows', len(prow),
      'slabs', len(out['slabs']), 'caps rows', len(caps), 'depth rows', len(dep))
