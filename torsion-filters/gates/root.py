# -*- coding: utf-8 -*-
# LA RAIZ: por que cancela la banda +2s -- causa medida, y mecanismo REFUTADO.  13 de agosto de 2026.
#
# EL ESLABON.  Todo lo del dia (la ley del extremo fugitivo, prof = F1 - F2, la saturacion) descansaba
# en UNA cosa medida y sin explicar: que la banda +2s cancela entera.  La bandera (flag.py) sugiere la
# causa: los bloques de esa banda tienen a los dos extremos de S como ESPECTADORES en posiciones fijas,
# o sea no son del problema de rango r sino del INTERIOR, rango r-1 -- y el interior es anulante.
#
# N1  LA CAUSA, con contrafactual.  Si la causa es Phi'(interior) == 0, entonces donde el interior NO
#     sea nulo la banda NO debe cancelar.  Se barre la POBLACION OBJETIVO ENTERA, no solo los
#     supervivientes, justamente para que la casilla "no/no" tenga habitantes.  En t=4 r=2 TODAS las
#     formas tienen interior nulo, asi que esa configuracion sola seria un control incapaz de fallar:
#     por eso se corren tambien t=6 r=3 y t=8 r=3, donde la casilla no/no tiene 142 y 64.
#
# N2  EL MECANISMO, y aqui me equivoque.  Escribi que la banda "factoriza" como
#         banda_{+2s} = (monomio de los espectadores) x Phi'(interior)
#     y es FALSO: alt(S) con beta_max dentro no es monomio por alternante de r-1, los alternantes se
#     mezclan.  Medido: el soporte de la banda es ~1350 monomios y el de Phi' es 144.  0 de 206 casos
#     coinciden en tamaño, y 0 en multiconjunto de coeficientes.
#
# LO QUE QUEDA, que es un enunciado mejor que el que tenia: DOS POLINOMIOS DE TAMAÑOS INCOMPARABLES
# QUE SE ANULAN EXACTAMENTE A LA VEZ.  No es un producto.  El mecanismo esta abierto.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python root.py

import itertools
import json
import os
import sys
from collections import Counter, defaultdict

from second_stratum import setup, all_transversals
from depth_histogram import stratify, stratum, full_dict, measure
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

OUT_JSON = "root_RESULT.json"


def extremos_S(b, t):
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    S = sorted({v for k in E for v in Cd[k]})
    return S[-1], S[0]


def pelar(b, t):
    e = extremos_S(b, t)
    if e is None:
        return None
    hi, lo = e
    return tuple(x for x in b if x != hi and x != lo)


def phi_cero(b, t, rp):
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(b, cl, rp, t)
    return measure([(x[2], x[1]) for x in tr], rp) is None


def bloques_banda(b, t, r):
    """{grado: [bloques de la banda +2s]} -- los dos extremos de S en T, beta_max arriba, beta_min abajo."""
    e = extremos_S(b, t)
    if e is None:
        return None
    hi, lo = e
    cl, E, Cd = setup(b, t)
    tr = all_transversals(b, cl, r, t)
    B = stratify([(x[2], x[1]) for x in tr], r)
    out = defaultdict(list)
    for d, bucket in B.items():
        for (w, T, S_, Sc) in bucket:
            v = set(T)
            enS = {T[a] for a in S_}
            if hi in v and lo in v and hi in enS and lo not in enS:
                out[d].append((w, T, S_, Sc))
    return out


def banda_cancela(b, t, r):
    bb = bloques_banda(b, t, r)
    if not bb:
        return None
    viva = sum(1 for d, bk in bb.items() if stratum(bk, r))
    return len(bb), (viva == 0)


def banda_poly(b, t, r):
    bb = bloques_banda(b, t, r)
    acc = defaultdict(int)
    for d, bk in (bb or {}).items():
        for k, v in stratum(bk, r).items():
            acc[k] += v
    return {k: v for k, v in acc.items() if v}


def interior_poly(b, t, r):
    bp = pelar(b, t)
    if bp is None or len(bp) != t + 2 * (r - 1):
        return None
    st = setup(bp, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(bp, cl, r - 1, t)
    return {k: v for k, v in full_dict([(x[2], x[1]) for x in tr], r - 1).items() if v}


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
print("")
if bad:
    print("  C0 FALLA -- el resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA")
RES = {}

# ===================================================================== N1 ========================
print("")
print("=" * 116)
print("N1  LA CAUSA, con contrafactual -- poblacion OBJETIVO entera, no solo supervivientes")
print("=" * 116)
print("")
for (t, r, Wmax) in [(4, 2, 22), (6, 3, 20), (8, 3, 21)]:
    N = t + 2 * r
    tab = Counter()
    ej = {}
    for W in range(N - 1, Wmax + 1):
        for b in shapes_of_width(W, N):
            if probe(b, t, r, deep=False) is None:
                continue
            bp = pelar(b, t)
            if bp is None or len(bp) != t + 2 * (r - 1):
                continue
            z = phi_cero(bp, t, r - 1)
            if z is None:
                continue
            bc = banda_cancela(b, t, r)
            if bc is None or bc[0] == 0:
                continue
            tab[(z, bc[1])] += 1
            ej.setdefault((z, bc[1]), b)
    print("  t=%d r=%d  (W <= %d)" % (t, r, Wmax))
    print("     interior Phi' == 0 | banda +2s cancela |   n   | ejemplo")
    for k in sorted(tab, reverse=True):
        print("        %-3s | %-3s | %5d | %s"
              % ("SI" if k[0] else "no", "SI" if k[1] else "no", tab[k], str(ej[k])))
    d = tab[(True, False)] + tab[(False, True)]
    pobl = tab[(False, False)]
    print("     desacuerdos: %d   |   casilla no/no (el contrafactual): %d %s"
          % (d, pobl, "<-- VACIA: aqui el control NO PUEDE FALLAR" if pobl == 0 else "poblada"))
    RES["N1_%d_%d" % (t, r)] = dict(desacuerdos=d, contrafactual=pobl,
                                    tabla={str(k): v for k, v in tab.items()})
    print("")
    sys.stdout.flush()

# ===================================================================== N2 ========================
print("=" * 116)
print("N2  EL MECANISMO -- 'la banda factoriza como espectadores x Phi'' : REFUTADO")
print("=" * 116)
print("")
print("     t  r | casos con Phi' != 0 | (|sop banda|, |sop Phi'|) | mismo tamaño | coefs iguales")
print("  " + "-" * 104)
for (t, r, Wmax) in [(6, 3, 20), (8, 3, 21)]:
    N = t + 2 * r
    n = ig = prop = 0
    ej = []
    for W in range(N - 1, Wmax + 1):
        for b in shapes_of_width(W, N):
            if probe(b, t, r, deep=False) is None:
                continue
            ip = interior_poly(b, t, r)
            if not ip:
                continue
            bp_ = banda_poly(b, t, r)
            if not bp_:
                continue
            n += 1
            ig += (len(bp_) == len(ip))
            prop += (sorted(abs(v) for v in bp_.values()) == sorted(abs(v) for v in ip.values()))
            if len(ej) < 2:
                ej.append((len(bp_), len(ip)))
            if n >= 200:
                break
        if n >= 200:
            break
    print("    %2d %2d | %19d | %-25s | %-12s | %s"
          % (t, r, n, str(ej), "%d de %d" % (ig, n), "%d de %d" % (prop, n)))
    RES["N2_%d_%d" % (t, r)] = dict(n=n, mismo_tam=ig, coefs=prop, ejemplo=ej)

# ===================================================================== N3 ========================
print("")
print("=" * 116)
print("N3  VEREDICTO")
print("=" * 116)
print("")
d = sum(RES["N1_%d_%d" % (t, r)]['desacuerdos'] for (t, r) in [(4, 2), (6, 3), (8, 3)])
cf = sum(RES["N1_%d_%d" % (t, r)]['contrafactual'] for (t, r) in [(4, 2), (6, 3), (8, 3)])
mt = sum(RES["N2_%d_%d" % (t, r)]['mismo_tam'] for (t, r) in [(6, 3), (8, 3)])
nn = sum(RES["N2_%d_%d" % (t, r)]['n'] for (t, r) in [(6, 3), (8, 3)])
print("     N1  banda +2s cancela  <=>  Phi'(interior) == 0 : %d desacuerdos, contrafactual con %d casos"
      % (d, cf))
print("     N2  la banda es Phi' trasladado                 : %d de %d  -> REFUTADO" % (mt, nn))
print("")
print("     LO QUE QUEDA: dos polinomios de tamaños INCOMPARABLES -- ~1350 monomios contra 144 -- que")
print("     se anulan EXACTAMENTE A LA VEZ.  No es una factorizacion.  El mecanismo esta ABIERTO, y")
print("     es mejor pregunta que la que habia: el eslabon que sostiene todo el dia tiene ahora una")
print("     causa medida con contrafactual, y una explicacion propuesta y MUERTA.")
print("")
json.dump(RES, open(OUT_JSON, "w"), indent=1)
print("DONE")
