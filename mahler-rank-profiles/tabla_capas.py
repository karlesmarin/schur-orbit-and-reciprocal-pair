# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: tablas de las capas del orden A = Z_p[c_1..c_m] por caracter (instrumento mahler_masas v3 desde el
#      2-oct: perfiles mod P por masas de clase, mismo primo P = (p, g(zeta_m)) que mahler_padico; antes
#      mahler_padico con Bernoulli exactos: 190/190 filas iguales. Forma
#      cerrada M0; control 31/31 frente al orden). Una fila por (n, p, piso a, orbita de Frobenius de caracteres impares
#      primitivos mod n), con etiquetas de Conrey (convencion LMFDB) para que cualquiera la cruce con otras bases.
#      Columnas:
#        n, p, a, conrey (etiquetas de la orbita, ';'), tam_orbita, orden, chi2 (exponente e: chi(2) = exp(2 pi i e/orden)),
#        chip_1 (chi(p) = 1: cero trivial), vueltas (P/A en J_k = p^{a-1}-1 + (k-1) p^{a-2}(p-1), k = 1..p-1),
#        lambda_tilde_min (primera vuelta presente - 1; '>=p-1' si ninguna), pares_ausentes (filas pares no pivote,
#        J <= Jmax, ';'), Jmax.
#      Salida: results/tabla_capas/capas.csv y capas.json (mismo contenido) y README_tabla.md (definiciones).
# FECHA: 1-oct-2026
# Uso (docker, /work): sage -python tabla_capas.py [n_max] [p...] [--piso1] [--dest=DIR] [--pisos=p:a:a,...]
import json
import os
import sys
import time
import warnings

sys.path.insert(0, '/work')
warnings.filterwarnings('ignore')
import mahler_padico as MP  # noqa: E402
import mahler_masas as MM  # noqa: E402
from sage.all import euler_phi  # noqa: E402


def filas_para(n, p, a):
    Jk = [p ** (a - 1) - 1 + (k - 1) * p ** (a - 2) * (p - 1) for k in range(1, p)]
    Jmax = Jk[-1]
    obs, _ = MP.orbitas(n, p)
    out = []
    for chi, orb, m, Km, chim in obs:
        ex = MP.exponentes(chim, Km, n)
        e2 = ex[2]
        piv = MM.matriz_B(ex, e2, m, n, p, a, Jmax).pivot_rows()  # v3, masas de clase; control 31/31
        v = ''.join('P' if J in piv else 'A' for J in Jk)
        primera = v.find('P')
        out.append({
            'n': n, 'p': p, 'a': a,
            'conrey': ';'.join(str(o.conrey_number()) for o in sorted(orb, key=lambda o: o.conrey_number())),
            'tam_orbita': len(orb), 'orden': int(m), 'chi2': int(e2), 'chip_1': bool(chim(p) == 1),
            'vueltas': v, 'lambda_tilde_min': (str(primera) if primera >= 0 else '>=%d' % (p - 1)),
            'pares_ausentes': ';'.join(str(J) for J in range(0, Jmax + 1, 2) if J not in piv), 'Jmax': Jmax,
        })
    return out


if __name__ == '__main__':
    # argumentos: nmax p1 p2 ... [--piso1] [--dest=DIR] [--pisos=p:a:a,...]
    args = [x for x in sys.argv[1:] if not x.startswith('--')]
    con_piso1 = '--piso1' in sys.argv
    dest = next((x.split('=', 1)[1] for x in sys.argv if x.startswith('--dest=')), '/work/results/tabla_capas')
    nmax = int(args[0]) if args else 50
    ps = [int(x) for x in args[1:]] or [5, 7]
    pisos = {5: [2, 3], 7: [2, 3], 11: [2], 13: [2]}
    for x in sys.argv:  # --pisos=11:2:3,13:2
        if x.startswith('--pisos='):
            for s in x.split('=', 1)[1].split(','):
                q = [int(v) for v in s.split(':')]
                pisos[q[0]] = q[1:]
    os.makedirs(dest, exist_ok=True)
    if con_piso1:
        from piso1 import piso1
    cols1 = ['d', 'e', 'vp_indice', 'Z_pn', 'Z_n', 'cert_K']
    filas, errores = [], []
    for p in ps:
        for n in range(3, nmax + 1, 2):
            if n % p == 0 or euler_phi(n) % p == 0:
                continue
            f1 = {}
            if con_piso1:
                try:
                    f1 = piso1(n, p)
                except Exception as e:  # noqa: BLE001
                    errores.append((n, p, 1, str(e)[:80]))
                    f1 = {c: 'ERROR' for c in cols1}
            for a in pisos.get(p, [2]):
                t0 = time.time()
                try:
                    fs = filas_para(n, p, a)
                except Exception as e:  # noqa: BLE001
                    errores.append((n, p, a, str(e)[:80]))
                    print("n=%d p=%d a=%d: ERROR %s" % (n, p, a, e), flush=True)
                    continue
                for f in fs:
                    f.update({('piso1_' + c): f1.get(c, '') for c in cols1} if con_piso1 else {})
                filas += fs
                print("n=%d p=%d a=%d: %d orbitas (%.1f s)" % (n, p, a, len(fs), time.time() - t0), flush=True)
    cols = ['n', 'p', 'a', 'conrey', 'tam_orbita', 'orden', 'chi2', 'chip_1', 'vueltas', 'lambda_tilde_min',
            'pares_ausentes', 'Jmax'] + (['piso1_' + c for c in cols1] if con_piso1 else [])
    print("ERRORES: %d %s" % (len(errores), errores[:10]), flush=True)
    with open(dest + '/capas.csv', 'w', encoding='utf-8') as fh:
        fh.write(','.join(cols) + '\n')
        for f in filas:
            fh.write(','.join(str(f[c]) for c in cols) + '\n')
    with open(dest + '/capas.json', 'w', encoding='utf-8') as fh:
        json.dump(filas, fh, indent=0)
    print("filas: %d -> %s" % (len(filas), dest), flush=True)
