# -*- coding: utf-8 -*-
"""Renombra las colisiones de notacion de note_t2 antes de integrarla en el paper.

Se escribe A FICHERO y se ejecuta: pasar esto por un heredoc de bash come las barras invertidas y
el regex se rompe en silencio (13303 coincidencias de V en 44 KB en el intento anterior).

  V  = valor del empate      ->  \\tau            (el paper usa V para la constante de Vandermonde)
  E  = clases de exceso      ->  \\mathcal{E}     (el paper usa E para las entradas PARES de beta)
  S  = valores de exceso     ->  \\mathcal{S}     (el paper usa S para un conjunto de t columnas)
  d_i = incrementos ordenados ->  \\Delta_{(i)}   (el paper usa d_1,d_2,d_3 para el invariante;
                                                  y \\delta ya esta cogido en la propia nota)
  Phi_t (r general)          ->  \\Phi_{t,r}      (el paper usa Phi_t para el objeto con r=1)
  N_t                        ->  N_{t,r}          (coherencia con lo anterior)

PROTEGIDOS, porque contienen las letras y NO se pueden tocar:
  \\Vdm  ·  \\S\\ref{...}  ·  \\S8  ·  D.~E.~Littlewood (inicial de autor)  ·  la bibliografia entera
"""
import re
import sys

p = sys.argv[1]
s = open(p, encoding='utf-8').read()
orig_len = len(s)

PROT = {r'\Vdm': '@@VDM@@', r'\S\ref': '@@SREF@@', r'\S8': '@@S8@@', r'D.~E.~': '@@DE@@'}
for k, v in PROT.items():
    s = s.replace(k, v)

i = s.index(r'\begin{thebibliography}')
body, bib = s[:i], s[i:]

# tokens sueltos: ni precedidos de letra o barra, ni seguidos de letra o digito
LETTERS = r'(?<![A-Za-z\\])%s(?![A-Za-z0-9])'
for sym, out in [('V', r'\tau'), ('E', r'\mathcal{E}'), ('S', r'\mathcal{S}')]:
    pat = LETTERS % sym
    n = len(re.findall(pat, body))
    if n > 200:
        raise SystemExit('ABORTA: %d coincidencias de %s, el patron esta roto' % (n, sym))
    body = re.sub(pat, lambda m, o=out: o, body)
    print('  %s -> %s : %d ocurrencias' % (sym, out, n))

for a, b in [(r'd_1\ge\dots\ge d_{2r}', r'\Delta_{(1)}\ge\dots\ge\Delta_{(2r)}'),
             (r'd_r=d_{r+1}', r'\Delta_{(r)}=\Delta_{(r+1)}'),
             (r'>d_r', r'>\Delta_{(r)}'),
             (r'=d_r', r'=\Delta_{(r)}'),
             (r'd_r$', r'\Delta_{(r)}$')]:
    k = body.count(a)
    body = body.replace(a, b)
    print('  %-24s -> %-28s : %d' % (a, b, k))

for a, b in [(r'\Phi_t', r'\Phi_{t,r}'), (r'\Phi_2', r'\Phi_{2,r}'), ('N_t', 'N_{t,r}')]:
    k = body.count(a)
    body = body.replace(a, b)
    print('  %-10s -> %-12s : %d' % (a, b, k))

s = body + bib
for k, v in PROT.items():
    s = s.replace(v, k)

if len(s) > orig_len * 1.15:
    raise SystemExit('ABORTA: el fichero crecio de %d a %d' % (orig_len, len(s)))

open(p, 'w', encoding='utf-8').write(s)
print('OK  %d -> %d bytes' % (orig_len, len(s)))
