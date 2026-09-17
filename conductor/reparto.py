# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# EL REPARTO LOCAL.  Para cada primo p del conductor:  a = v_p([O+:A]) = delta_p,  b = v_p([A:c]),
# con a + b = v_p(N(c)) = sum_P f.e.  Gorenstein local <=> a = b.  La pregunta ya no es "cuanto vale
# el indice" sino "QUE decide el reparto", porque (3,21) y (6,28) tienen el MISMO c local -- un solo
# P con f=3, e=2 -- y reparten distinto: 3+3 contra 4+2.  Luego no lo decide c.
#
# Se agrupa por la forma local de c y se imprimen TODOS los repartos vistos con esa forma.  Una
# forma con un solo reparto no dice nada (no se ha contrastado); las interesantes son las que
# tienen dos o mas, porque ahi hay una variable escondida y se puede buscar cual es.
import io, re, sys
from math import gcd
from collections import defaultdict

def leer(rutas):
    filas = {}
    for ruta in rutas:
        try:
            txt = io.open(ruta, encoding='utf-8', errors='replace').read()
        except IOError:
            continue
        for ln in txt.splitlines():
            ln = re.sub(r'^((sage:|\.\.\.\.:)\s*)+', '', ln).rstrip()
            mt = re.match(r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(SI|no)\s+(.*)$', ln)
            if not mt:
                continue
            m, q, d, idx, nc, acx, gor, resto = mt.groups()
            tri = re.findall(r'\((\d+),(\d+),(\d+)\)', resto)
            if not tri:
                continue
            filas[(int(m), int(q))] = dict(m=int(m), q=int(q), d=int(d), idx=int(idx),
                                           nc=int(nc), acx=int(acx), gor=(gor == 'SI'),
                                           fac=[(int(a), int(b), int(c)) for (a, b, c) in tri])
    return [filas[k] for k in sorted(filas)]

def vp(n, p):
    v = 0
    while n and n % p == 0:
        n //= p
        v += 1
    return v

def orden(a, n):
    if n <= 1:
        return 1
    if gcd(a, n) != 1:
        return 0
    k, x = 1, a % n
    while x != 1:
        x = (x * a) % n
        k += 1
        if k > 4 * n:
            return -1
    return k

filas = leer(sys.argv[1:])
ng = [r for r in filas if not r['gor']]
print("filas: %d   (no-Gorenstein %d)" % (len(filas), len(ng)))
print("")
print("TODAS las no-Gorenstein, con su reparto primo a primo:")
print("   m    q     p   forma de c en p      f.e    a=v_p([O+:A])  b=v_p([A:c])   a-b")
for r in ng:
    porp = defaultdict(list)
    for (p, f, e) in r['fac']:
        porp[p].append((f, e))
    for p in sorted(porp):
        fe = sum(f * e for (f, e) in porp[p])
        a, b = vp(r['idx'], p), vp(r['acx'], p)
        print("  %3d %4d %5d   %-18s %4d %10d %13d %8d"
              % (r['m'], r['q'], p, " ".join("(%d,%d)" % t for t in porp[p]), fe, a, b, a - b))

print("")
print("AGRUPADO por forma local de c.  Una forma con UN solo reparto no ha contrastado nada:")
grupos = defaultdict(list)
for r in filas:
    porp = defaultdict(list)
    for (p, f, e) in r['fac']:
        porp[p].append((f, e))
    for p in sorted(porp):
        forma = tuple(sorted(porp[p]))
        grupos[forma].append((r['m'], r['q'], p, vp(r['idx'], p), vp(r['acx'], p), r['gor']))
for forma in sorted(grupos, key=lambda k: (-len(grupos[k]), str(k))):
    reps = set((a, b) for (_, _, _, a, b, _) in grupos[forma])
    marca = "*** DOS REPARTOS: hay variable escondida ***" if len(reps) > 1 else ("un solo reparto (%d casillas)" % len(grupos[forma]) if len(grupos[forma]) > 1 else "una sola casilla: no mide")
    print("  forma %-22s  repartos %-22s  %s" % (str(forma), str(sorted(reps)), marca))
    if len(reps) > 1:
        for (m, q, p, a, b, g) in grupos[forma]:
            print("       m=%d q=%d p=%d -> a=%d b=%d  %s   ord(p mod q/p^v)=%s  d=phi(q)/2"
                  % (m, q, p, a, b, "Gor" if g else "noGor",
                     orden(p, q // (p ** vp(q, p))) if q // (p ** vp(q, p)) > 1 else 1))
