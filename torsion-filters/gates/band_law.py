# -*- coding: utf-8 -*-
# LA BANDA +2s: de donde sale el valor de saturacion.  13 de agosto de 2026.
#
# QUE CIERRA.  saturation.py dejo abierto el valor en que se congela 'cancelan': los histogramas
# {4}, {2,18,24}, {17,18,24}, {2,32,36,38,42} no tenian explicacion, y CUATRO candidatas mias
# murieron (0 de 98 las tres primeras; la cuarta, sat==2 <-> ¬H*, necesaria y no suficiente).
#
# LAS CUATRO FALLABAN POR LO MISMO: contaban BANDAS MEZCLADAS.  Bajo eps_s = s*(1,0,..,0,-1) un bloque
# de Laplace (T, A) se desplaza segun donde caen los extremos de S, y eso los ordena en bandas que SE
# SEPARAN al iterar:
#
#     +2s   los dos extremos en T, beta_max en A y beta_min NO en A     <- contiene D
#      +s   UN extremo                                                  <- aqui esta 'first'
#      0, -s, -2s   el resto                                            <- por DEBAJO de first
#
# Yo contaba todos los de dos extremos, mezclando +2s con 0 y -2s.  Solo la primera esta por encima de
# first, y toda ella cancela.  De donde
#
#     SATURACION  =  numero de GRADOS DISTINTOS de la banda +2s          (alcance: sat > 2)
#
# Y POR QUE ESE NUMERO.  En un bloque de la banda +2s, beta_max y beta_min son ESPECTADORES FIJOS
# (alto y bajo).  Con T' = T sin los dos extremos (tamaño 2r-2) y A' = A sin beta_max (tamaño r-1):
#
#     deg = 2*sum_A(T) - sum(T) = (beta_max - beta_min) + ( 2*sum_A'(T') - sum(T') )
#
# IDENTIDAD, no ajuste: la banda +2s ES el espectro de grados del INTERIOR CONGELADO A RANGO r-1,
# trasladado.  Luego la saturacion es el tamaño de ese espectro -- finito y fijo PORQUE el interior
# esta congelado, que es exactamente por que satura.
#
# COLUMNAS
#   C0  fatal: probe() contra scan(), y los testigos publicados.
#   N1  la IDENTIDAD, verificada elemento a elemento (no por conteo): el multiconjunto de grados de la
#       banda +2s es el del interior a rango r-1 trasladado por beta_max - beta_min.
#   N2  LA LEY: sat == |banda +2s|, con su alcance sat > 2, y las tres candidatas muertas como
#       SEÑUELOS CON DENOMINADOR: contar bloques en vez de grados, la banda 0, y la banda +s.
#       Y se DICE que la banda -2s acierta lo mismo y NO es señuelo: es su espejo bajo sigma_C.
#   N3  la poblacion degenerada sat == 2, que la formula no describe.  Con su numero.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python band_law.py

import itertools
import json
import os
import sys
from collections import Counter

from second_stratum import setup, all_transversals
from survivors_wide import scan

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
assert "def probe(" in _head and "def shapes_of_width(" in _head, "k_vs_m.py cambio de forma"
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]
shapes_of_width = _ns["shapes_of_width"]

CFG = [(4, 2, 26, 30), (6, 3, 21, 20), (8, 3, 22, 16), (6, 4, 21, 16)]
SMIN = {4: 2, 6: 6, 8: 8}
OUT_JSON = "band_law_RESULT.json"


def extremos_S(beta, t):
    cl, E, Cd = setup(beta, t)
    S = sorted({v for k in E for v in Cd[k]})
    return S[-1], S[0]


def bandas(beta, t, r):
    """{desplazamiento: multiconjunto de grados}.  El desplazamiento se lee de donde caen los
    extremos de S en (T, A), que es lo unico que eps_s toca."""
    cl, E, Cd = setup(beta, t)
    hi, lo = extremos_S(beta, t)
    out = {}
    for (_, T, _, _) in all_transversals(beta, cl, r, t):
        v = set(T)
        nh, nl = hi in v, lo in v
        for A in itertools.combinations(T, r):
            sA = set(A)
            if nh and nl:
                sh = 2 if (hi in sA and lo not in sA) else (-2 if (lo in sA and hi not in sA) else 0)
            elif nh:
                sh = 1 if hi in sA else -1
            elif nl:
                sh = 1 if lo not in sA else -1
            else:
                sh = 0
            out.setdefault(sh, Counter())[2 * sum(A) - sum(T)] += 1
    return out


def interior_rango(beta, t, r, rango):
    """multiconjunto de grados del INTERIOR (T sin los dos extremos de S) al rango pedido."""
    cl, E, Cd = setup(beta, t)
    hi, lo = extremos_S(beta, t)
    acc = Counter()
    for (_, T, _, _) in all_transversals(beta, cl, r, t):
        v = set(T)
        if hi not in v or lo not in v:
            continue
        Tp = tuple(x for x in T if x != hi and x != lo)
        if len(Tp) < rango:
            continue
        for A in itertools.combinations(Tp, rango):
            acc[2 * sum(A) - sum(Tp)] += 1
    return acc


def traj(beta, t, r, s, J):
    hi, lo = extremos_S(beta, t)
    out = []
    for j in range(J + 1):
        b = tuple(sorted([(x + s * j if x == hi else (x - s * j if x == lo else x)) for x in beta],
                         reverse=True))
        rec = probe(b, t, r)
        if rec is None or not rec['surv'] or rec['prof'] is None:
            break
        out.append(rec['vac_cancelan'])
    return out


def sat(c):
    return c[-1] if len(c) >= 4 and c[-1] == c[-2] == c[-3] else None


def survivors(t, r, Wmax):
    N = t + 2 * r
    out = []
    for W in range(N - 1, Wmax + 1):
        for beta in shapes_of_width(W, N):
            rec = probe(beta, t, r)
            if rec and rec['surv'] and rec['prof'] is not None:
                out.append(beta)
    return out


# ===================================================================== C0 ========================
print("=" * 116)
print("C0  ACEPTACION -- fatal")
print("=" * 116)
print("")
bad = 0
for (t, r, M) in [(4, 2, 15), (6, 3, 18)]:
    n_ref, cont_ref, sv_ref = scan(t, r, M)
    mine, mb = Counter(), []
    for comb in itertools.combinations(range(M + 1), t + 2 * r):
        beta = tuple(sorted(comb, reverse=True))
        rec = probe(beta, t, r, deep=False)
        if rec is None:
            continue
        mine[(rec['e'] == t, rec['surv'])] += 1
        if rec['surv']:
            mb.append(beta)
    ok = (sum(mine.values()) == n_ref and mine == cont_ref
          and sorted(mb) == sorted(x['beta'] for x in sv_ref))
    bad += not ok
    print("  C0a  probe() == scan()  t=%d r=%d M=%d : %s" % (t, r, M, "ok" if ok else "*** FALLA ***"))
for (beta, p, c) in [((18, 17, 11, 8, 7, 6, 1, 0), 6, 2), ((38, 23, 21, 18, 17, 16, 15, 0), 16, 4)]:
    rec = probe(beta, 4, 2)
    ok = rec is not None and rec['prof'] == p and rec['vac_cancelan'] == c
    bad += not ok
    print("  C0b  testigo %-38s prof %s/%d cancelan %s/%d  %s"
          % (str(beta), rec['prof'] if rec else "-", p, rec['vac_cancelan'] if rec else "-", c,
             "ok" if ok else "*** FALLA ***"))
print("")
if bad:
    print("  C0 FALLA -- el resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA")
RES = {}

SV = {(t, r): survivors(t, r, Wmax) for (t, r, Wmax, _) in CFG}

# ===================================================================== N1 ========================
print("")
print("=" * 116)
print("N1  LA IDENTIDAD, elemento a elemento -- no por conteo, que seria mucho mas debil")
print("=" * 116)
print("")
print("     banda +2s  ==  espectro del interior a rango r-1  TRASLADADO por (beta_max - beta_min)")
print("")
for (t, r, Wmax, J) in CFG:
    n = ok = 0
    for beta in SV[(t, r)][:40]:
        hi, lo = extremos_S(beta, t)
        b2 = bandas(beta, t, r).get(2, Counter())
        it = interior_rango(beta, t, r, r - 1)
        tras = Counter({d + hi - lo: m for d, m in it.items()})
        n += 1
        ok += (b2 == tras)
    print("     t=%d r=%d : el MULTICONJUNTO coincide en %d de %d" % (t, r, ok, n))
    RES["N1_%d_%d" % (t, r)] = dict(n=n, ok=ok)
    sys.stdout.flush()

# ===================================================================== N2 ========================
print("")
print("=" * 116)
print("N2  LA LEY, con su alcance, y los tres SEÑUELOS con denominador")
print("=" * 116)
print("")
print("     t  r | sat>2 | LEY |banda+2s| | señ. #bloques+2s | señ. banda 0 | señ. banda +s | espejo -2s")
print("  " + "-" * 112)
TOT = Counter()
for (t, r, Wmax, J) in CFG:
    n = l = s1 = s2 = s3 = sp = 0
    hs = Counter()
    for beta in SV[(t, r)][:40]:
        v = sat(traj(beta, t, r, SMIN[t], J))
        if v is None:
            continue
        hs[v] += 1
        TOT[v > 2] += 1
        if v <= 2:
            continue
        n += 1
        B = bandas(beta, t, r)
        g2 = B.get(2, Counter())
        l += (v == len(g2))
        s1 += (v == sum(g2.values()))
        s2 += (v == len(B.get(0, Counter())))
        s3 += (v == len(B.get(1, Counter())))
        sp += (v == len(B.get(-2, Counter())))
    print("    %2d %2d | %5d | %-15s | %-16s | %-12s | %-13s | %s"
          % (t, r, n, "%d de %d" % (l, n), "%d de %d" % (s1, n), "%d de %d" % (s2, n),
             "%d de %d" % (s3, n), "%d de %d" % (sp, n)))
    print("           histograma de saturacion: %s" % dict(sorted(hs.items())))
    RES["N2_%d_%d" % (t, r)] = dict(n=n, ley=l, sen_bloques=s1, sen_b0=s2, sen_bs=s3, espejo=sp,
                                    hist={str(k): v for k, v in hs.items()})
    sys.stdout.flush()
print("")
print("     El espejo -2s acierta lo mismo que la ley y NO es señuelo: es la imagen de la banda +2s")
print("     bajo sigma_C.  Se dice para que nadie lo lea como una confirmacion independiente.")

# ===================================================================== N3 ========================
print("")
print("=" * 116)
print("N3  VEREDICTO")
print("=" * 116)
print("")
n = sum(RES["N2_%d_%d" % (t, r)]['n'] for (t, r, _, _) in CFG)
l = sum(RES["N2_%d_%d" % (t, r)]['ley'] for (t, r, _, _) in CFG)
i1 = sum(RES["N1_%d_%d" % (t, r)]['ok'] for (t, r, _, _) in CFG)
i0 = sum(RES["N1_%d_%d" % (t, r)]['n'] for (t, r, _, _) in CFG)
print("     N1  la identidad banda+2s == interior(r-1) trasladado : %d de %d" % (i1, i0))
print("     N2  sat == |banda +2s|, alcance sat > 2               : %d de %d" % (l, n))
print("     la poblacion degenerada sat == 2                      : %d formas, NO descritas" % TOT[False])
print("")
if l == n and i1 == i0:
    print("     LA SATURACION ES EL TAMAÑO DEL ESPECTRO DE GRADOS DEL INTERIOR CONGELADO, A RANGO r-1.")
    print("     Y satura PORQUE el interior esta congelado: su espectro es finito y fijo, y los")
    print("     extremos solo lo trasladan.  El eslabon MEDIDO es que toda la banda +2s cancela;")
    print("     el otro es una identidad de una linea.")
else:
    print("     NO CUADRA -- mirar N1 y N2, que es para lo que estan.")
print("")
json.dump(RES, open(OUT_JSON, "w"), indent=1)
print("DONE")
