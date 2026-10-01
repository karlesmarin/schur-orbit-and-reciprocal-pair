# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: la regla de Frobenius como PROPAGACION (sustituye a regla_frobenius_tabla.py).
#      S_a = filas PARES ausentes (no pivote) del piso a, J <= L. Los dos pisos a y a+1 se miden hasta el MISMO L
#      (rango completo; L = p^{a+1} - 1 por defecto) con mahler_masas (v3, control 31/31, P = (p, g(zeta_m))).
#      Se registran, por orbita de Frobenius:
#        hechos formales (prueba de regresion): (F1) S_{a+1} contenido en S_a; (F2) S_{a+1} cap [0,p^a) = S_a cap [0,p^a);
#        I_a     = {J in S_a : Phi(J) <= L, Phi(J) no en S_{a+1}}   (fallos de la conjetura (3) Phi(S_a) c S_{a+1});
#        B_{a+1} = S_{a+1} \ ((S_a cap [0,p^a)) cup Phi(S_a))       (ausencias sin explicar), cada una con su
#                  certificado de dependencia: fila J = sum_{K<J} c_K fila K en la matriz B del piso a+1 (se verifica).
#      Phi(J) = p(J+1) - 1. Las filas impares (caracter impar) se omiten.
#      Salida: results/propagacion/propagacion.csv y .json.
# FECHA: 2-oct-2026
# Uso (docker, /work): sage -python propagacion_frobenius.py [nmax] [p:a ...]      (defecto 50 5:2 7:2 5:3)
import json
import os
import sys
import time
import warnings

sys.path.insert(0, '/work')
warnings.filterwarnings('ignore')
import mahler_padico as MP  # noqa: E402
import mahler_masas as MM  # noqa: E402
from sage.all import euler_phi, kronecker_character  # noqa: E402


def medir(ex, e2, m, n, p, a, L):
    Ba = MM.matriz_B(ex, e2, m, n, p, a, L)
    Bb = MM.matriz_B(ex, e2, m, n, p, a + 1, L)
    pa, pb = set(Ba.pivot_rows()), set(Bb.pivot_rows())
    Sa = [J for J in range(0, L + 1, 2) if J not in pa]
    Sb = [J for J in range(0, L + 1, 2) if J not in pb]
    phi = lambda J: p * (J + 1) - 1  # noqa: E731
    sa, sb = set(Sa), set(Sb)
    F1 = sb <= sa
    F2 = {J for J in sb if J < p ** a} == {J for J in sa if J < p ** a}
    I = [J for J in Sa if phi(J) <= L and phi(J) not in sb]
    expl = {J for J in sa if J < p ** a} | {phi(J) for J in sa}
    Bn = [J for J in Sb if J not in expl]
    certs = {}
    for J in Bn:
        c = Bb[:J].solve_left(Bb[J])
        assert c * Bb[:J] == Bb[J]
        certs[J] = [(K, str(c[K])) for K in range(J) if c[K] != 0]
    return dict(L=L, S_a=Sa, S_b=Sb, F1=F1, F2=F2, I_a=I, B_b=Bn, cert=certs)


def fila(n, p, a, conrey, tam, m, e2, r):
    return dict(n=int(n), p=int(p), a=int(a), conrey=conrey, tam_orbita=int(tam), orden=int(m), chi2=int(e2),
                L=int(r['L']), n_S_a=len(r['S_a']), n_S_b=len(r['S_b']), F1=bool(r['F1']), F2=bool(r['F2']),
                I_a=';'.join(map(str, r['I_a'])), B_b=';'.join(map(str, r['B_b'])),
                cert={str(k): [(int(K), v) for K, v in c] for k, c in r['cert'].items()},
                S_a=';'.join(map(str, r['S_a'])), S_b=';'.join(map(str, r['S_b'])))


if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    pas = [tuple(int(v) for v in s.split(':')) for s in sys.argv[2:]] or [(5, 2), (7, 2), (5, 3)]
    filas = []
    for p, a in pas:
        L = p ** (a + 1) - 1
        for n in range(3, nmax + 1, 2):
            if n % p == 0 or euler_phi(n) % p == 0:
                continue
            t0 = time.time()
            obs, _ = MP.orbitas(n, p)
            for chi, orb, m, Km, chim in obs:
                ex = MP.exponentes(chim, Km, n)
                r = medir(ex, ex[2], m, n, p, a, L)
                conrey = ';'.join(str(o.conrey_number()) for o in sorted(orb, key=lambda o: o.conrey_number()))
                filas.append(fila(n, p, a, conrey, len(orb), m, ex[2], r))
            print("n=%d p=%d a=%d->%d L=%d: %d orbitas (%.1f s)" % (n, p, a, a + 1, L, len(obs), time.time() - t0),
                  flush=True)
    # lambda alta (5, 75619), caracter cuadratico
    n, p = 75619, 5
    chi = kronecker_character(-n)
    ex = [-1 if v == 0 else (0 if v == 1 else 1) for v in chi.values()]
    for a in (2, 3):
        t0 = time.time()
        r = medir(ex, ex[2], 2, n, p, a, p ** (a + 1) - 1)
        filas.append(fila(n, p, a, 'kronecker(-%d)' % n, 1, 2, ex[2], r))
        print("n=%d p=%d a=%d->%d: S_a %s | S_b %s | I_a %s | B_b %s (%.1f s)" % (
            n, p, a, a + 1, r['S_a'], r['S_b'], r['I_a'], r['B_b'], time.time() - t0), flush=True)
    dest = '/work/results/propagacion'
    os.makedirs(dest, exist_ok=True)
    cols = ['n', 'p', 'a', 'conrey', 'tam_orbita', 'orden', 'chi2', 'L', 'n_S_a', 'n_S_b', 'F1', 'F2', 'I_a', 'B_b',
            'S_a', 'S_b']
    with open(dest + '/propagacion.csv', 'w', encoding='utf-8') as fh:
        fh.write(','.join(cols) + '\n')
        for f in filas:
            fh.write(','.join(str(f[c]) for c in cols) + '\n')
    with open(dest + '/propagacion.json', 'w', encoding='utf-8') as fh:
        json.dump(filas, fh, indent=0)
    assert filas
    print("orbitas: %d; F1 falla: %d; F2 falla: %d; con I_a no vacio: %d; con B_{a+1} no vacio: %d" % (
        len(filas), sum(not f['F1'] for f in filas), sum(not f['F2'] for f in filas),
        sum(bool(f['I_a']) for f in filas), sum(bool(f['B_b']) for f in filas)), flush=True)
