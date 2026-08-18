# -*- coding: utf-8 -*-
# EL DESCENSO Y LA FRONTERA: la ruta VERTICAL para t >= 4.  13 de agosto de 2026.
#
# POR QUE.  La ruta que hemos matado tres veces esta semana era HORIZONTAL: reducir "Phi == 0" a k
# estratos con k fijo.  Murio porque no hay k -- la profundidad no esta acotada.  Lo de hoy sugiere
# otra: la anulacion es HEREDITARIA HACIA ABAJO, y entonces el programa es una INDUCCION EN r.
#
# UN TEST MIO ANTERIOR FUE VACIO, y el fallo esta escrito aqui porque es instructivo: lo monte sobre
# la "poblacion objetivo" de probe(), y probe() devuelve None cuando (i) es cierta -- o sea su
# poblacion es EXACTAMENTE donde Phi != 0.  El antecedente A = "Phi_r == 0" no tenia un solo
# habitante y "A => B, 0 fallos" era un control incapaz de fallar.  Aqui se enumera SIN probe().
#
# N1  EL DESCENSO.  A = (Phi_r == 0) ; B = (Phi_{r-1} == 0 del interior, o sea de beta sin los dos
#     extremos de S).  Se cuentan los habitantes de A, y solo entonces se lee A => B.
#
# N2  LA FRONTERA, formulada SIN CIRCULARIDAD.  En una induccion se puede asumir el criterio a rango
#     r-1, luego B da que el INTERIOR es concentrico, con su propio centro C' = max S' + min S'.  Y
#     beta es concentrico si S = S' + {hi, lo} es simetrico respecto de C = hi + lo.  El interior ya
#     es simetrico -- pero respecto de C'.  Luego lo unico que falta es que los dos centros coincidan:
#
#         dado B :   Phi_r == 0   <=>   C' == C     ( + la condicion (ii), que se mide aparte )
#
#     Eso es aritmetica entre dos centros, NO el criterio: no es circular.  Si falla, se dice, y la
#     frontera es otra cosa.
#     SEÑUELOS con denominador: "C' == C mod t" (mas debil) y "el interior es concentrico" a secas
#     (que B ya implica, luego no puede discriminar y tiene que salir mal).
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python descent.py

import itertools
import json
import os
import sys
from collections import Counter, defaultdict

from second_stratum import setup, all_transversals
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

CFG = [(4, 2, 17), (6, 2, 17), (6, 3, 16)]
OUT_JSON = "descent_RESULT.json"


def datos(b, t, rp):
    """(Phi==0, concentrico, C, S, E) -- calculado SIN probe(), para no heredar su filtro.
    E son las clases de EXCESO, que la condicion (ii) completa necesita."""
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    S = sorted({v for k in E for v in Cd[k]})
    C = S[0] + S[-1]
    tr = all_transversals(b, cl, rp, t)
    z = measure([(x[2], x[1]) for x in tr], rp) is None
    return z, (set(C - v for v in S) == set(S)), C, S, set(E)


def pelar(b, t):
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    S = sorted({v for k in E for v in Cd[k]})
    hi, lo = S[-1], S[0]
    return tuple(x for x in b if x != hi and x != lo)


def cond_ii(C, t, E):
    """(ii) COMPLETA, verbatim de note_t2/t2_criterion.tex lineas 82-93:

        #{k in Z/t : 2k = C (mod t)} = 2   Y   ambas k estan en E (las clases de EXCESO)

    ARREGLADO el 13 de agosto tras una auditoria adversarial.  Yo habia escrito
    `any((2*k - C) % t == 0 for k in range(t))`: falta el "= 2" Y falta la pertenencia a E, y esa
    segunda cláusula es el 100 % del daño.  Con la version truncada, C0b fallaba en 118 de 5448
    (t=4 r=2) y 162 de 3033 (t=6 r=2) -- y para t IMPAR la version truncada es una TAUTOLOGIA
    (x -> 2x es biyectiva en Z/t), asi que predecia anulacion donde no se anula nada: 175 fallos de
    4585 en t=3 r=2, 305 de 4509 en t=5 r=2, con CERO ceros reales.
    Con la (ii) completa: 0 desacuerdos en ~334 000 formas, t = 1..10, r = 1..4.
    El criterio no tenia defecto; mi transcripcion estaba truncada."""
    sols = [k for k in range(t) if (2 * k - C) % t == 0]
    return len(sols) == 2 and all(k in E for k in sols)


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

# C0b: el criterio publicado, re-verificado por esta via independiente -- Phi == 0 <=> (i) y (ii).
# Si esto fallara, mi datos() no mide lo que digo y N1/N2 no valdrian.
print("")
for (t, r, Wmax) in [(4, 2, 15), (6, 2, 15)]:
    N = t + 2 * r
    m = 0
    n = 0
    for W in range(N - 1, Wmax + 1):
        for b in shapes_of_width(W, N):
            d = datos(b, t, r)
            if d is None:
                continue
            z, conc, C, S, E = d
            n += 1
            m += (z == (conc and cond_ii(C, t, E)))
    bad += (m != n)
    print("  C0b  Phi == 0  <=>  (i) y (ii)   t=%d r=%d : %d de %d   %s"
          % (t, r, m, n, "ok" if m == n else "*** FALLA ***"))
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
print("N1  EL DESCENSO -- con el numero de habitantes de A ANTES de leer la implicacion")
print("=" * 116)
print("")
DATA = {}
for (t, r, Wmax) in CFG:
    N = t + 2 * r
    tab = Counter()
    filas = []
    # F6 (auditoria): antes se hacia 'continue' en silencio cuando el interior perdia ocupacion, y
    # eso tiraba el 31-38 % de los habitantes de A sin decirlo.  Ahora se CUENTAN y se imprimen.
    caidas = Counter()
    for W in range(N - 1, Wmax + 1):
        for b in shapes_of_width(W, N):
            d = datos(b, t, r)
            if d is None:
                continue
            A, concA, C, S, EA = d
            bp = pelar(b, t)
            if bp is None or len(bp) != t + 2 * (r - 1):
                caidas[('tamaño', A)] += 1
                continue
            d2 = datos(bp, t, r - 1)
            if d2 is None:
                caidas[('el interior pierde ocupacion', A)] += 1
                continue
            B, concB, Cp, Sp, EB = d2
            tab[(A, B)] += 1
            filas.append((b, A, B, C, Cp, concA, concB, EA))
    nA = tab[(True, True)] + tab[(True, False)]
    print("  t=%d r=%d  (%d formas, W <= %d)" % (t, r, len(filas), Wmax))
    print("     A: Phi_r == 0 | B: Phi' == 0 |   n")
    for k in sorted(tab, reverse=True):
        print("        %-3s | %-3s | %6d" % ("SI" if k[0] else "no", "SI" if k[1] else "no", tab[k]))
    perdidos = sum(v for k, v in caidas.items() if k[1])
    print("     habitantes de A: %d %s" % (nA, "<-- VACIO, el test NO vale" if nA == 0 else ""))
    print("     descartados ANTES de la tabla: %s   (de ellos con A cierto: %d, o sea el %.0f%% de A)"
          % (dict(caidas) or "ninguno", perdidos,
             100.0 * perdidos / (nA + perdidos) if (nA + perdidos) else 0))
    if nA:
        print("     A => B : %s  (%d fallos de %d)"
              % ("SI" if tab[(True, False)] == 0 else "NO", tab[(True, False)], nA))
    DATA[(t, r)] = filas
    RES["N1_%d_%d" % (t, r)] = dict(n=len(filas), nA=nA, fallos=tab[(True, False)],
                                    tabla={str(k): v for k, v in tab.items()})
    print("")
    sys.stdout.flush()

# ===================================================================== N2 ========================
print("=" * 116)
print("N2  LA FRONTERA -- dado B, es  Phi_r == 0  <=>  C' == C ?")
print("=" * 116)
print("")
for (t, r, Wmax) in CFG:
    filas = [f for f in DATA[(t, r)] if f[2]]          # solo B cierto
    tab = Counter()
    s1 = Counter()
    s2 = Counter()
    for (b, A, B, C, Cp, concA, concB, EA) in filas:
        tab[(A, Cp == C)] += 1
        s1[(A, (Cp - C) % t == 0)] += 1                # señuelo: solo mod t
        s2[(A, concB)] += 1                            # señuelo: el interior concentrico a secas
    d = tab[(True, False)] + tab[(False, True)]
    print("  t=%d r=%d  (%d formas con B cierto)" % (t, r, len(filas)))
    print("     A: Phi_r == 0 | C' == C |   n")
    for k in sorted(tab, reverse=True):
        print("        %-3s | %-3s | %6d" % ("SI" if k[0] else "no", "SI" if k[1] else "no", tab[k]))
    print("     desacuerdos: %d  ->  %s" % (d, "LA FRONTERA ES C' == C" if d == 0 else "NO es C' == C sola"))
    print("     señuelo C' == C (mod t)         : %d desacuerdos"
          % (s1[(True, False)] + s1[(False, True)]))
    print("     señuelo interior concentrico    : %d desacuerdos  (B ya lo implica, tiene que fallar)"
          % (s2[(True, False)] + s2[(False, True)]))
    # y si no es C'==C sola, se prueba con (ii)
    if d:
        tab2 = Counter()
        for (b, A, B, C, Cp, concA, concB, EA) in filas:
            tab2[(A, (Cp == C) and cond_ii(C, t, EA))] += 1
        d2 = tab2[(True, False)] + tab2[(False, True)]
        print("     con (ii) añadida: C' == C y (ii)  : %d desacuerdos %s"
              % (d2, "<-- ESA es la frontera" if d2 == 0 else ""))
        RES["N2b_%d_%d" % (t, r)] = d2
    RES["N2_%d_%d" % (t, r)] = dict(n=len(filas), desacuerdos=d,
                                    sen_modt=s1[(True, False)] + s1[(False, True)],
                                    sen_conc=s2[(True, False)] + s2[(False, True)])
    print("")
    sys.stdout.flush()

# ===================================================================== N3 ========================
print("=" * 116)
print("N3  VEREDICTO")
print("=" * 116)
print("")
nA = sum(RES["N1_%d_%d" % (t, r)]['nA'] for (t, r, _) in CFG)
fa = sum(RES["N1_%d_%d" % (t, r)]['fallos'] for (t, r, _) in CFG)
print("     DESCENSO   Phi_r == 0 => Phi_{r-1} == 0 del interior : %d fallos sobre %d habitantes" % (fa, nA))
for (t, r, _) in CFG:
    k2 = RES.get("N2b_%d_%d" % (t, r))
    print("     FRONTERA   t=%d r=%d : C'==C -> %d desacuerdos%s"
          % (t, r, RES["N2_%d_%d" % (t, r)]['desacuerdos'],
             ("   |   C'==C y (ii) -> %d" % k2) if k2 is not None else ""))
print("")
print("     LO QUE ESTO SERIA SI AGUANTA: la ruta deja de ser HORIZONTAL (acotar estratos, muerta")
print("     tres veces esta semana) y pasa a ser una INDUCCION EN r, con t=2 cerrado y r=1 de base.")
print("     LO QUE NO ES: una prueba.  El descenso esta MEDIDO, no probado, y la frontera tambien.")
print("")
json.dump(RES, open(OUT_JSON, "w"), indent=1)
print("DONE")
