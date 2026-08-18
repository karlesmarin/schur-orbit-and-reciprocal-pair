# -*- coding: utf-8 -*-
# LA BANDERA: beta no es una configuracion con un numero, es una CADENA.  13 de agosto de 2026.
#
# DE DONDE SALE.  Siete lecturas mias han muerto hoy -- el cono, la ocupacion, s=2<=>t=2r, K_sop<=4,
# las cuatro candidatas a ley de saturacion, la multiplicidad de factores y "una sola pareja rota" --
# y las siete tenian la misma forma: buscaba UN ESCALAR.  Lo unico que no ha muerto es la recursion:
# quitando los dos extremos de S queda una configuracion ANULANTE de rango r-1.  La pregunta que no
# habia hecho es si eso se puede volver a hacer.
#
# SE PUEDE.  Pelando repetidamente -- quitar los dos extremos de S, bajar un rango -- sale una TORRE, y
# [Phi]_top = 0 EN TODOS LOS PISOS.  Lo que cambia de piso en piso es cuando entra la anulacion TOTAL
# Phi == 0.  Ejemplo (t=6 r=4):
#
#     rango 4 : |58| T      [Phi]_top = 0  pero Phi != 0     <- el superviviente
#     rango 3 : |40| T      [Phi]_top = 0  pero Phi != 0
#     rango 2 : |24| Z      Phi == 0
#     rango 1 : | 4| Z      Phi == 0
#
# y el numero de pisos hasta Z VARIA con beta: es un invariante entero, no una constante.  O sea el
# objeto es una BANDERA, y los escalares que se me morian eran sombras suyas.  Eso ademas re-justifica
# el vocabulario que el gate de literatura habia dado por prestado -- "flag valuations that record
# successive orders of vanishing along nested sequences of subvarieties" no era vocabulario: era la
# descripcion.
#
# Y LA BANDERA PAGA: cierra el agujero que saturation.py y band_law.py dejaron abierto.  La ley
# sat = |banda +2s| valia en 90 de 90 con alcance "sat > 2", y las 8 formas con sat = 2 quedaban SIN
# DESCRIBIR.  Ya no:
#
#     sat = |banda +2s|   <=>   EL PRIMER PELADO ES TOTALMENTE ANULANTE (Phi' == 0)
#
# 164 si / 8 no, CERO desacuerdos.  Y el señuelo -- pedir solo [Phi']_top == 0 -- FALLA EN 6, asi que
# lo portante es la anulacion total, no la del estrato de arriba.
#
# COLUMNAS
#   C0  fatal: probe() contra scan(), y los testigos publicados.
#   N1  LA TORRE, piso a piso, con el tamaño del espectro de cada uno.  Se enseña entera.
#   N2  cuantos pisos aguantan, por configuracion.  Si fuera constante seria una columna; sale un
#       histograma, y por eso es un invariante de beta.
#   N3  LA HIPOTESIS: sat = |banda+2s| <=> primer pelado Z.  Con su señuelo y su denominador.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python flag.py

import itertools
import json
import os
import sys
from collections import Counter, defaultdict

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

CFG = [(4, 2, 26, 30), (6, 3, 21, 20), (8, 3, 22, 16), (6, 4, 21, 16)]
SMIN = {4: 2, 6: 6, 8: 8}
OUT_JSON = "flag_RESULT.json"


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


def estado(b, t, rp):
    """([Phi]_top == 0, Phi == 0) a rango rp.  None si no aplica."""
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(b, cl, rp, t)
    D = max(x[3] for x in tr)
    G = [x for x in tr if x[3] == D]
    top0 = (len(G) == 2 and inv_of(G[0][1], rp) == inv_of(G[1][1], rp) and G[0][2] == -G[1][2])
    return top0, (measure([(x[2], x[1]) for x in tr], rp) is None)


def espectro(b, t, rp):
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    g = set()
    for (_, T, _, _) in all_transversals(b, cl, rp, t):
        for A in itertools.combinations(T, rp):
            g.add(2 * sum(A) - sum(T))
    return len(g)


def bandera(b, t, r):
    """[(rango, beta, |espectro|, estado)] pelando hasta que no se pueda."""
    out = []
    cur, rp = b, r
    while rp >= 1:
        out.append((rp, cur, espectro(cur, t, rp), estado(cur, t, rp)))
        nx = pelar(cur, t)
        if nx is None or len(nx) != t + 2 * (rp - 1):
            break
        cur, rp = nx, rp - 1
    return out


def banda2(b, t, r):
    e = extremos_S(b, t)
    if e is None:
        return None
    hi, lo = e
    cl, E, Cd = setup(b, t)
    g = set()
    for (_, T, _, _) in all_transversals(b, cl, r, t):
        v = set(T)
        if hi not in v or lo not in v:
            continue
        for A in itertools.combinations(T, r):
            if hi in A and lo not in A:
                g.add(2 * sum(A) - sum(T))
    return len(g)


def sat_de(b, t, r, s, J):
    hi, lo = extremos_S(b, t)
    o = []
    for j in range(J + 1):
        x = tuple(sorted([(y + s * j if y == hi else (y - s * j if y == lo else y)) for y in b],
                         reverse=True))
        rec = probe(x, t, r)
        if rec is None or not rec['surv'] or rec['prof'] is None:
            break
        o.append(rec['vac_cancelan'])
    return o[-1] if len(o) >= 4 and o[-1] == o[-2] == o[-3] else None


def survivors(t, r, Wmax):
    N = t + 2 * r
    out = []
    for W in range(N - 1, Wmax + 1):
        for b in shapes_of_width(W, N):
            rec = probe(b, t, r)
            if rec and rec['surv'] and rec['prof'] is not None:
                out.append(b)
    return out


# ===================================================================== C0 ========================
# EL GUARD __main__, anadido el 15 de agosto de 2026, y es el arreglo de un DEFECTO REAL.
# peel_zero.py hace  "from flag import pelar, extremos_S"  -- y sin guard eso EJECUTABA la corrida
# entera de este guion en cada import.  Consecuencia medida: toda salida de un guion que importe
# peel_zero (cancel_vs_vanishing.py entre ellos) llevaba por delante la tabla de LA TORRE, que no
# es suya.  Es el mismo arreglo que ya se le hizo a survivors_wide.py y por la misma razon: cambio
# de PURA INDENTACION, verificado por diff inverso contra git.
if __name__ == "__main__":
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
    for (beta, p) in [((18, 17, 11, 8, 7, 6, 1, 0), 6), ((38, 23, 21, 18, 17, 16, 15, 0), 16)]:
        rec = probe(beta, 4, 2)
        ok = rec is not None and rec['prof'] == p
        bad += not ok
        print("  C0b  testigo %-38s prof %s/%d  %s"
              % (str(beta), rec['prof'] if rec else "-", p, "ok" if ok else "*** FALLA ***"))
    print("")
    if bad:
        print("  C0 FALLA -- el resto NO vale.")
        print("DONE (veredicto suspendido)")
        raise SystemExit(1)
    print("  C0 PASA")
    RES = {}
    SV = {(t, r): survivors(t, r, W) for (t, r, W, _) in CFG}

    # ===================================================================== N1 ========================
    print("")
    print("=" * 116)
    print("N1  LA TORRE, entera.  Z = Phi == 0 ;  T = solo [Phi]_top == 0 ;  |n| = tamaño del espectro")
    print("=" * 116)
    for (t, r, Wmax, J) in CFG:
        print("")
        print("  t=%d r=%d" % (t, r))
        for b in SV[(t, r)][:4]:
            s = sat_de(b, t, r, SMIN[t], J)
            F = bandera(b, t, r)
            pisos = " ".join("r%d:|%s|%s" % (rp, sp, "Z" if (e and e[1]) else ("T" if (e and e[0]) else "-"))
                             for (rp, x, sp, e) in F)
            print("     %-46s sat=%-5s banda+2s=%-4s  %s" % (str(b), s, banda2(b, t, r), pisos))
        sys.stdout.flush()

    # ===================================================================== N2 ========================
    print("")
    print("=" * 116)
    print("N2  ALTURA DE LA BANDERA -- si fuera constante saldria UNA columna; sale un histograma")
    print("=" * 116)
    print("")
    for (t, r, Wmax, J) in CFG:
        c = Counter()
        for b in SV[(t, r)]:
            F = bandera(b, t, r)
            z = sum(1 for (_, _, _, e) in F if e and e[1])
            c[(len(F), z)] += 1
        print("     t=%d r=%d (%3d formas) : (pisos, pisos con Phi == 0) -> %s"
              % (t, r, len(SV[(t, r)]), dict(sorted(c.items()))))
        RES["N2_%d_%d" % (t, r)] = {str(k): v for k, v in c.items()}

    # ===================================================================== N3 ========================
    print("")
    print("=" * 116)
    print("N3  LA HIPOTESIS QUE FALTABA -- y cierra la poblacion degenerada de saturation.py")
    print("=" * 116)
    print("")
    tab = defaultdict(int)
    sen = defaultdict(int)
    for (t, r, Wmax, J) in CFG:
        for b in SV[(t, r)]:
            s = sat_de(b, t, r, SMIN[t], J)
            if s is None:
                continue
            bp = pelar(b, t)
            e = estado(bp, t, r - 1) if (bp is not None and len(bp) == t + 2 * (r - 1)) else None
            tab[(bool(e and e[1]), s == banda2(b, t, r))] += 1
            sen[(bool(e and e[0]), s == banda2(b, t, r))] += 1
    print("     primer pelado Z (Phi' == 0) | sat == |banda+2s| |  n")
    for k in sorted(tab, reverse=True):
        print("        %-3s | %-3s | %d" % ("SI" if k[0] else "no", "SI" if k[1] else "no", tab[k]))
    d = tab[(True, False)] + tab[(False, True)]
    print("     desacuerdos: %d  ->  %s" % (d, "LA HIPOTESIS ES EXACTA" if d == 0 else "no separa"))
    print("")
    print("     SEÑUELO: pedir solo [Phi']_top == 0")
    for k in sorted(sen, reverse=True):
        print("        %-3s | %-3s | %d" % ("SI" if k[0] else "no", "SI" if k[1] else "no", sen[k]))
    ds = sen[(True, False)] + sen[(False, True)]
    print("     desacuerdos del señuelo: %d  ->  %s"
          % (ds, "tambien separa, luego la Z NO es portante" if ds == 0
             else "NO separa: la anulacion TOTAL es lo portante, no la del estrato de arriba"))
    RES['N3'] = dict(desacuerdos=d, senuelo=ds, n=sum(tab.values()))

    # ===================================================================== N4 ========================
    print("")
    print("=" * 116)
    print("N4  VEREDICTO")
    print("=" * 116)
    print("")
    print("     sat = |banda +2s|  <=>  el primer pelado es TOTALMENTE anulante : %d desacuerdos en %d formas"
          % (RES['N3']['desacuerdos'], RES['N3']['n']))
    print("     el señuelo ([Phi']_top == 0 a secas) falla en %d, luego la Z es portante" % RES['N3']['senuelo'])
    print("")
    print("     Con esto, la ley de saturacion deja de tener poblacion sin describir: las formas con")
    print("     sat = 2 son exactamente aquellas cuya BANDERA tarda mas de un piso en llegar a Phi == 0.")
    print("")
    print("     LO QUE NO ES: una prueba de que la bandera gobierne prof.  Eso NO esta testado.  Y el")
    print("     alcance sigue siendo el de los barridos, con tope de reloj.")
    print("")
    json.dump(RES, open(OUT_JSON, "w"), indent=1)
    print("DONE")
