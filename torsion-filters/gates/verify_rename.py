# -*- coding: utf-8 -*-
"""Verifica el renombrado en el FUENTE: aplica el mapeo inverso al .tex actual y lo compara con la
version de git.  Si el renombrado fue exactamente lo previsto y nada mas, el diff sale VACIO."""
import difflib
import subprocess
import sys

cur = open(sys.argv[1], encoding='utf-8').read()
old = subprocess.run(['git', 'show', 'HEAD:research/orbit-pair/note_t2/t2_criterion.tex'],
                     capture_output=True, cwd=sys.argv[2]).stdout.decode('utf-8')

inv = cur
for a, b in [(r'\Phi_{t,r}', r'\Phi_t'), (r'\Phi_{2,r}', r'\Phi_2'), ('N_{t,r}', 'N_t'),
             (r'\Delta_{(1)}\ge\dots\ge\Delta_{(2r)}', r'd_1\ge\dots\ge d_{2r}'),
             (r'\Delta_{(r)}=\Delta_{(r+1)}', r'd_r=d_{r+1}'),
             (r'>\Delta_{(r)}', r'>d_r'), (r'=\Delta_{(r)}', r'=d_r'),
             (r'\Delta_{(r)}$', r'd_r$'),
             (r'\mathcal{E}', 'E'), (r'\mathcal{S}', 'S'), (r'\tau', 'V')]:
    inv = inv.replace(a, b)

a_lines = old.splitlines()
b_lines = inv.splitlines()
diff = [l for l in difflib.unified_diff(a_lines, b_lines, lineterm='', n=1)
        if l and l[0] in '+-' and not l.startswith(('+++', '---'))]

print('lineas git: %d   invertidas: %d' % (len(a_lines), len(b_lines)))
print('DIFERENCIAS tras deshacer el mapeo: %d' % len(diff))
for l in diff[:40]:
    print('   ', l[:160])
if not diff:
    print('>>> RENOMBRADO EXACTO: no se toco nada mas que los simbolos previstos.')
