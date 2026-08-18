# -*- coding: utf-8 -*-
# EL CONSTRUCTOR: supervivientes de profundidad A LA CARTA, sin barrer.  13 de agosto de 2026.
#
# DE DONDE SALE.  La recursion medida hoy: quitando a un superviviente los DOS EXTREMOS DE S quedan
# t + 2(r-1) elementos -- una configuracion del mismo t un rango mas abajo -- y esa configuracion es
# ANULANTE: [Phi']_top = 0 en 96 de 96, y Phi' == 0 en 82 de 96 (las 14 restantes son la poblacion
# degenerada sat = 2).  O sea:
#
#     un superviviente profundo = una configuracion ANULANTE de rango r-1 + dos extremos que le
#     rompen la concentricidad
#
# y la ley del extremo fugitivo dice que esos dos extremos se pueden alejar sin fin, subiendo prof
# de s en s, PORQUE EL NUCLEO ANULANTE NO SE TOCA.
#
# PARA QUE SIRVE.  Hasta hoy un superviviente de profundidad prescrita costaba un barrido exhaustivo:
# C(W-1, N-2) formas, 1.1 millones para W=34 y 14 millones para W=50, 35 minutos de reloj para llegar
# a K = 22.  Con la recursion se FABRICA: nucleo anulante + dos extremos + iterar.  Coste O(1) por
# caso y profundidad a la carta.  Es lo que le falta a la suite de regresion del programa t >= 4.
#
# Y ALCANZA DONDE EL BARRIDO NO LLEGA: un caso fabricado es (13,8,6,3,2,1,0,-5), con una entrada
# NEGATIVA.  Ningun barrido sobre subconjuntos de {0..M} puede producirlo.
#
# LO QUE ESTE GUION TIENE QUE DEMOSTRAR, y por eso esta partido en dos:
#   N2  RECICLADO -- nucleos sacados de supervivientes ya barridos.  Es el control de que el mecanismo
#       es el que digo, pero NO demuestra que se pueda construir: los nucleos ya venian de casos buenos.
#   N3  DE CERO -- nucleos ANULANTES generados enumerando el rango r-1, que no han salido de ningun
#       barrido de supervivientes.  Si estos
#       tambien producen supervivientes, el constructor construye de verdad.  Es la columna que vale.
# Cada pieza fabricada se verifica con probe(), que es el criterio del barrido, no uno mio.
#
# Y SE DICE EL RENDIMIENTO, no solo los aciertos: cuantos nucleos generados dan superviviente y
# cuantos no.  Un constructor que acierta 3 de 500 tambien "construye", y no es lo mismo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python construct.py

import itertools
import json
import os
import sys
import time
from collections import Counter

from second_stratum import setup, all_transversals, inv_of
from depth_histogram import measure
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

SMIN = {4: 2, 6: 6, 8: 8, 10: 10}      # paso del regimen robusto: s == 0 (mod t) y par
OUT_JSON = "construct_RESULT.json"


def extremos_S(beta, t):
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    S = sorted({v for k in E for v in Cd[k]})
    return S[-1], S[0]


def nucleo(beta, t):
    """el interior: beta sin los dos extremos de S."""
    e = extremos_S(beta, t)
    if e is None:
        return None
    hi, lo = e
    return tuple(x for x in beta if x != hi and x != lo)


def es_anulante(bp, t, rp):
    """([Phi']_top == 0, Phi' == 0) para el nucleo a rango r-1.  None si no aplica."""
    st = setup(bp, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(bp, cl, rp, t)
    D = max(x[3] for x in tr)
    G = [x for x in tr if x[3] == D]
    top0 = (len(G) == 2 and inv_of(G[0][1], rp) == inv_of(G[1][1], rp) and G[0][2] == -G[1][2])
    return top0, (measure([(x[2], x[1]) for x in tr], rp) is None)


def vestir(bp, t, r, ancho=2):
    """añade DOS extremos al nucleo y devuelve los que salen supervivientes: [(beta, prof)].
    hi por encima del maximo y lo por debajo del minimo, barriendo 'ancho' periodos de t."""
    M, m = max(bp), min(bp)
    out = []
    for hi in range(M + 1, M + ancho * t + 1):
        for lo in range(m - ancho * t, m):
            beta = tuple(sorted((hi,) + bp + (lo,), reverse=True))
            rec = probe(beta, t, r)
            if rec and rec['surv'] and rec['prof'] is not None:
                out.append((beta, rec['prof']))
    return out


def subir(beta, t, r, K):
    """itera el extremo fugitivo hasta profundidad K exacta.  None si K no es alcanzable."""
    s = SMIN[t]
    rec = probe(beta, t, r)
    p0 = rec['prof']
    if K < p0 or (K - p0) % s:
        return None
    j = (K - p0) // s
    e = extremos_S(beta, t)
    hi, lo = e
    b = tuple(sorted([(x + s * j if x == hi else (x - s * j if x == lo else x)) for x in beta],
                     reverse=True))
    r2 = probe(b, t, r)
    if r2 is None or not r2['surv'] or r2['prof'] != K:
        return None
    return b


def concentricos(t, rp, Wmax, tope):
    """NUCLEOS GENERADOS DE CERO: se enumeran TODAS las configuraciones de rango rp (t + 2rp enteros,
    todas las clases ocupadas) y se filtran las ANULANTES.  No sale ni una de un barrido de
    supervivientes de rango rp+1: el espacio de partida es otro.

    PRIMERA VERSION EQUIVOCADA, y su fallo es instructivo: genere beta' SIMETRICO bajo x -> C - x, o
    sea concentrico COMO CONJUNTO, y salieron 0 supervivientes de 145 nucleos.  La concentricidad no
    se mide sobre beta' sino sobre su conjunto de EXCESO S' -- el mismo error de variable compuesta
    que ya me costo diez residuos falsos por la manana.  El nucleo real (17,11,8,7,6,1) NO es
    simetrico (11+6 = 17 != 18) y es anulante igual."""
    n = t + 2 * rp
    out = []
    for W in range(n - 1, Wmax + 1):
        for b in shapes_of_width(W, n):
            an = es_anulante(b, t, rp)
            if an and an[0]:
                out.append(b)
                if len(out) >= tope:
                    return out
    return out


# ===================================================================== C0 ========================
print("=" * 116)
print("C0  ACEPTACION -- fatal.  Y el control de ida y vuelta: reconstruir un superviviente conocido")
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
# C0b, IDA Y VUELTA.  Primera version mal escrita, y su fallo fue instructivo: exigia que 'vestir'
# encontrase el beta original con los extremos DONDE ESTAN, y en (38,23,21,18,17,16,15,0) estan a
# distancia 15 del nucleo -- fuera de cualquier ventana razonable, porque a ese ya lo empujo el
# generador.  Lo que hay que exigir es lo que el mecanismo dice: vestir CERCA y luego ITERAR.
def orbita(b, t, r, Wtope):
    hi, lo = extremos_S(b, t)
    s = SMIN[t]
    out = []
    for j in range(200):
        x = tuple(sorted([(y + s * j if y == hi else (y - s * j if y == lo else y)) for y in b],
                         reverse=True))
        if x[0] - x[-1] > Wtope:
            break
        out.append(x)
    return out

for (beta, t, r) in [((18, 17, 11, 8, 7, 6, 1, 0), 4, 2), ((38, 23, 21, 18, 17, 16, 15, 0), 4, 2),
                     ((26, 25, 15, 12, 11, 10, 1, 0), 4, 2)]:
    bp = nucleo(beta, t)
    an = es_anulante(bp, t, r - 1)
    hechos = vestir(bp, t, r)
    W = beta[0] - beta[-1]
    vuelve = any(beta in orbita(b, t, r, W) for (b, p) in hechos)
    bad += not vuelve
    print("  C0b  %-38s nucleo %-28s anulante %s; %d vestidos cerca, y ITERANDO vuelve a salir: %s"
          % (str(beta), str(bp), an, len(hechos), "si" if vuelve else "*** NO ***"))
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
print("N1  LA REGLA DE COLOCACION -- donde caen los dos extremos respecto al nucleo, medido")
print("=" * 116)
print("")
CFG = [(4, 2, 26), (6, 3, 21), (8, 3, 22)]
for (t, r, Wmax) in CFG:
    N = t + 2 * r
    cl_hi = Counter()
    cl_lo = Counter()
    n = 0
    for W in range(N - 1, Wmax + 1):
        for beta in shapes_of_width(W, N):
            rec = probe(beta, t, r)
            if not (rec and rec['surv'] and rec['prof'] is not None):
                continue
            e = extremos_S(beta, t)
            bp = nucleo(beta, t)
            if bp is None or len(bp) != t + 2 * (r - 1):
                continue
            hi, lo = e
            C = hi + lo
            n += 1
            cl_hi[(2 * (hi % t) - C) % t == 0] += 1
            cl_lo[(2 * (lo % t) - C) % t == 0] += 1
    print("     t=%d r=%d (%d formas) : beta_max en clase FIJA de sigma_C en %d, fuera en %d"
          % (t, r, n, cl_hi[True], cl_hi[False]))
    print("                              beta_min en clase FIJA en %d, fuera en %d"
          % (cl_lo[True], cl_lo[False]))
    RES["N1_%d_%d" % (t, r)] = dict(n=n, hi_fija=cl_hi[True], lo_fija=cl_lo[True])
    sys.stdout.flush()

# ===================================================================== N2 ========================
print("")
print("=" * 116)
print("N2  RECICLADO -- nucleos de supervivientes ya barridos.  Control del mecanismo, NO del constructor")
print("=" * 116)
print("")
for (t, r, Wmax) in CFG:
    N = t + 2 * r
    nucleos = []
    for W in range(N - 1, Wmax + 1):
        for beta in shapes_of_width(W, N):
            rec = probe(beta, t, r)
            if rec and rec['surv'] and rec['prof'] is not None:
                bp = nucleo(beta, t)
                if bp and len(bp) == t + 2 * (r - 1) and bp not in nucleos:
                    nucleos.append(bp)
        if len(nucleos) >= 8:
            break
    hechos = 0
    profs = set()
    for bp in nucleos[:8]:
        v = vestir(bp, t, r)
        hechos += len(v)
        profs |= {p for (_, p) in v}
    print("     t=%d r=%d : %d nucleos -> %d supervivientes vestidos, profundidades %s"
          % (t, r, len(nucleos[:8]), hechos, sorted(profs)))
    RES["N2_%d_%d" % (t, r)] = dict(nucleos=len(nucleos[:8]), hechos=hechos, profs=sorted(profs))
    sys.stdout.flush()

# ===================================================================== N3 ========================
print("")
print("=" * 116)
print("N3  DE CERO -- nucleos ANULANTES generados enumerando rango r-1.  LA COLUMNA QUE VALE")
print("=" * 116)
print("")
for (t, r, Wmax) in CFG:
    t0 = time.time()
    anul = concentricos(t, r - 1, 3 * t + 8, 300)
    gen = anul
    hechos = []
    for bp in anul[:60]:
        hechos += vestir(bp, t, r)
    print("     t=%d r=%d : %d nucleos ANULANTES generados de cero (enumerando rango r-1, no reciclando),"
          % (t, r, len(anul)))
    print("               %d vestidos dan SUPERVIVIENTE, profundidades %s   (%.0f s)"
          % (len(hechos), sorted({p for (_, p) in hechos}), time.time() - t0))
    if hechos:
        b, p = hechos[0]
        print("               ejemplo FABRICADO: %s  prof=%d" % (str(b), p))
    RES["N3_%d_%d" % (t, r)] = dict(generados=len(gen), anulantes=len(anul), supervivientes=len(hechos),
                                    profs=sorted({p for (_, p) in hechos}))
    sys.stdout.flush()

# ===================================================================== N4 ========================
print("")
print("=" * 116)
print("N4  PROFUNDIDAD A LA CARTA -- se pide K, se entrega beta, y se VERIFICA con probe()")
print("=" * 116)
print("")
PED = [(4, 2, 30), (4, 2, 60), (4, 2, 100), (6, 3, 40), (6, 3, 100), (8, 3, 60)]
print("     pedido           beta entregado                                     prof verificada  coste")
print("  " + "-" * 112)
entregas = []
for (t, r, K) in PED:
    t0 = time.time()
    N = t + 2 * r
    base = None
    for W in range(N - 1, 3 * t + 14):
        for beta in shapes_of_width(W, N):
            rec = probe(beta, t, r)
            if rec and rec['surv'] and rec['prof'] is not None:
                b = subir(beta, t, r, K)
                if b:
                    base = (beta, b)
                    break
        if base:
            break
    dt = time.time() - t0
    if base:
        beta0, b = base
        v = probe(b, t, r)
        ok = v and v['surv'] and v['prof'] == K
        entregas.append(ok)
        print("     t=%d r=%d K=%-4d %-50s %-16s %.1f s"
              % (t, r, K, str(b)[:50], "%d  %s" % (v['prof'], "ok" if ok else "*** MAL ***"), dt))
    else:
        entregas.append(False)
        print("     t=%d r=%d K=%-4d NO ENTREGADO (K no alcanzable con paso s=%d desde las semillas vistas)"
              % (t, r, K, SMIN[t]))
RES['N4'] = dict(pedidos=len(PED), entregados=sum(entregas))

# ===================================================================== N5 ========================
print("")
print("=" * 116)
print("N5  VEREDICTO")
print("=" * 116)
print("")
gen_ok = sum(1 for (t, r, _) in CFG if RES["N3_%d_%d" % (t, r)]['supervivientes'] > 0)
print("     N3  construccion DE CERO (nucleos generados, no reciclados) : funciona en %d de %d configuraciones"
      % (gen_ok, len(CFG)))
print("     N4  profundidad a la carta                                  : %d de %d entregados y verificados"
      % (RES['N4']['entregados'], RES['N4']['pedidos']))
print("")
print("     COSTE, que es el punto: alcanzar K = 22 por barrido exhaustivo costo 44 anchuras y 2123 s")
print("     (support_ladder.py, t=4 r=2 hasta W=50, 14 millones de formas en la ultima anchura).")
print("     Aqui la profundidad se pide y se entrega.")
print("")
print("     LO QUE ESTO NO ES: una biyeccion.  No se afirma que TODO superviviente se construya asi,")
print("     solo que estos se construyen y se verifican.  El reciproco -- que todo nucleo anulante")
print("     vestido de dos extremos da superviviente -- es FALSO y su rendimiento esta en N3.")
print("")
json.dump(RES, open(OUT_JSON, "w"), indent=1)
print("DONE")
