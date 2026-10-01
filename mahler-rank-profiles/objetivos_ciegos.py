# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: test A CIEGAS: 16 caracteres NO cuadraticos (p = 7, 11) con lambda >= 2 no forzada (ni chi(2) = 1 ni
#      chi(p) = 1), con lambda calculada aparte por Bernoulli generalizados e interpolacion (objetivos_ciegos.csv;
#      las vueltas no se usaron para elegirlos). Por cada objetivo (n, p, Conrey):
#        - comprueba chi(2), chi(p) mod P frente a los suyos (convencion P = (p, g(zeta_m)), mahler_padico.primo_sobre_p);
#        - vueltas del piso 2 (mahler_masas, v3) -> lambda~ = primera presente - 1 (tope p-1);
#        - lambda por nuestra serie de Iwasawa (filas_negativas, rama <x>), independiente de la suya;
#        - y la compara con su lambda_value.
# FECHA: 2-oct-2026
# Uso (docker, /work, -e OPENBLAS_NUM_THREADS=4): sage -python objetivos_ciegos.py
import csv
import sys
import time

sys.path.insert(0, '/work')
import filas_negativas as FN  # noqa: E402
import mahler_masas as MM  # noqa: E402
import mahler_padico as MP  # noqa: E402
from sage.all import ZZ, CyclotomicField, DirichletGroup, GF, Integers  # noqa: E402

tot = ok_v = ok_s = ok_c = 0
for r in csv.DictReader(open('/work/objetivos_ciegos.csv', encoding='utf-8')):
    p, n, conrey, lam = int(r['p']), int(r['n']), int(r['conrey']), int(r['lambda_value'])
    t0 = time.time()
    # su definicion: n primo, generador gen de (Z/n)^x, chi(gen) = zeta_m^e (zeta_m = exp(2 pi i/m))
    assert ZZ(n).is_prime()
    m, gen, e = int(r['character_order']), int(r['generator']), int(r['generator_exponent'])
    dlog = [0] * n
    x = 1
    for t in range(n - 1):
        dlog[x] = t
        x = x * gen % n
    assert x == 1 and len(set(dlog[1:])) == n - 1, "gen no es raiz primitiva"
    ex = [-1] + [(e * dlog[u]) % m for u in range(1, n)]
    assert ex[n - 1] == m // 2, "chi no es impar"
    # etiqueta de Conrey por Sage (chi en el generador de Sage)
    Km = CyclotomicField(m)
    G = DirichletGroup(n, Km)
    u0 = int(Integers(n).unit_gens()[0])
    chi = G([Km.gen() ** ex[u0]])
    assert chi.order() == m and chi.conductor() == n and chi(-1) == -1
    conrey_sage = int(chi.conrey_number())
    g = MP.primo_sobre_p(m, p)
    raiz = [x for x in GF(p) if g(x) == 0] if g.degree() == 1 else []
    c2 = (int(raiz[0]) ** ex[2]) % p if raiz else None
    cp = (int(raiz[0]) ** ex[p % n]) % p if raiz else None
    conv = (c2 == int(r['chi2_mod_p']) and cp == int(r['chip_mod_p'])) if raiz else None
    Jk = [p - 1 + (k - 1) * (p - 1) for k in range(1, p)]
    piv = set(MM.matriz_B(ex, ex[2], m, n, p, 2, Jk[-1]).pivot_rows())
    v = ''.join('P' if J in piv else 'A' for J in Jk)
    lt = v.find('P') if 'P' in v else p - 1
    mu, gg, f = FN.medida(ex, m, n, p, 4)
    ls = FN.lam(FN.serie_iwasawa(mu, p, 4, 0, p + 3))
    tot += 1
    ok_v += (lt == min(lam, p - 1))
    ok_s += (ls == lam)
    ok_c += bool(conv)
    print("p=%d n=%d conrey=%d orden=%d: lambda fichero=%d | vueltas %s -> %d | nuestra serie %s | chi(2),chi(p) mod P = %s,%s "
          "(fichero %s,%s) %s | conrey Sage %d %s (%.1f s)" % (
              p, n, conrey, m, lam, v, lt, ls, c2, cp, r['chi2_mod_p'], r['chip_mod_p'], 'OK' if conv else 'DISTINTO',
              conrey_sage, '=' if conrey_sage == conrey else 'DISTINTA', time.time() - t0), flush=True)
assert tot > 0
print("OBJETIVOS CIEGOS: %d; vueltas == lambda del fichero: %d; serie de Iwasawa == lambda del fichero: %d; convencion P coincide: %d" % (
    tot, ok_v, ok_s, ok_c), flush=True)
