# -*- coding: utf-8 -*-
# C = tau: EL HUECO DEL 12, y su reduccion a H*.  13 de agosto de 2026.
#
# EL HUECO.  note_t2/t2_criterion.tex prueba (cor:ii) la condicion (ii) CON tau EN LUGAR DE C, y solo
# "si ademas e = 2 entonces C = tau".  Verbatim: "The two halves of the corollary must be kept apart.
# The statement about tau is unconditional; the statement about C --- which is the one that appears in
# (ii) as we defined it --- is not."  Cerrar C = tau vuelve INCONDICIONAL el enunciado
# Phi_t == 0 => (ii), y quita la correccion F4 que la auditoria encontro.
#
# LA REDUCCION, que es lo que este guion mide.  Con los dos maximizadores g_A, g_B:
#     S = K + g_com + {x1, x2, y1, y2},     K = T_A ∩ T_B,  g_com = la parte comun de los picks
# y ya esta PROBADO que sigma_V(K) = K (sigma_V es involucion y intercambia T_A y T_B, luego fija su
# interseccion) y que x1+x2 = y1+y2 = V (el empate).  El paso que faltaba se escribia como
# "sigma_V(g_com) = g_com, que nada obliga".  PERO NO HACE FALTA:
#
#     sigma_V(K) = K  =>  min K + max K = V
#     si ademas  min S, max S  estan EN K,  entonces  C = min S + max S = min K + max K = V = tau.
#
# "los extremos de S estan en K" es exactamente H*, la hipotesis aislada por otra via en
# closed_form_prof.py: los DOS maximizadores conservan los DOS extremos de S.  O sea:
#
#     H*  =>  C = tau        (algebra, dos lineas)
#
# y lo unico que queda medido es H* mismo sobre las formas con [Phi]_top = 0.
#
# COLUMNAS
#   C0  fatal: probe() contra scan(); y el criterio contra (i)&(ii) COMPLETA -- con las dos clausulas,
#       que es donde me equivoque esta tarde.
#   N1  sobre las formas con [Phi]_top = 0: se cuenta H*, se cuenta C == tau, y se cruzan.  Con el
#       DENOMINADOR de cuantas tienen e > 2, que es el caso que la prueba NO cubre.
#   N2  la cadena entera, paso a paso, con su tabla: sigma_V(K) = K ; extremos en K ; C = tau.
#       Si algun paso falla, se ve CUAL.
#   N3  SEÑUELOS: (a) sigma_V(g_com) = g_com -- lo que el .tex decia que hacia falta, y que segun esta
#       reduccion NO hace falta; si sale que tampoco se cumple siempre, mejor: confirma que el atajo
#       es necesario.  (b) e == 2, el caso ya cubierto, para ver cuanto aporta la parte nueva.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python c_eq_tau.py

import itertools
import json
import os
import sys
from collections import Counter

from second_stratum import setup, all_transversals, inv_of
from survivors_wide import scan
from _control import Contingencia, Implicacion

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
assert "def probe(" in _head and "def shapes_of_width(" in _head, "k_vs_m.py cambio de forma"
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]
shapes_of_width = _ns["shapes_of_width"]

CFG = [(4, 2, 20), (6, 2, 20), (4, 3, 18), (6, 3, 18), (8, 2, 20), (8, 3, 19), (10, 2, 20)]
OUT_JSON = "c_eq_tau_RESULT.json"


def anatomia(b, t, r):
    """Devuelve None si no aplica; si no, un dict con todo lo que la cadena necesita.
    'top0' = [Phi]_top == 0 con el mismo criterio que probe()."""
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    S = sorted({v for k in E for v in Cd[k]})
    C = S[0] + S[-1]
    tr = all_transversals(b, cl, r, t)
    D = max(x[3] for x in tr)
    G = [x for x in tr if x[3] == D]
    if len(G) != 2:
        return dict(top0=False, e=len(E), S=S, C=C)
    a, bb = G
    top0 = (inv_of(a[1], r) == inv_of(bb[1], r) and a[2] == -bb[2])
    TA, TB = set(a[1]), set(bb[1])
    K = TA & TB
    gA = set(S) - TA                     # los picks de A dentro de S
    gB = set(S) - TB
    gcom = gA & gB
    dif = sorted((gA | gB) - gcom)       # los cuatro del empate
    tau = (sum(dif) / 2.0) if len(dif) == 4 else None
    if tau is not None and tau == int(tau):
        tau = int(tau)
    return dict(top0=top0, e=len(E), S=S, C=C, K=sorted(K), gcom=sorted(gcom), dif=dif, tau=tau,
                TA=sorted(TA), TB=sorted(TB))


def simetrico(X, c):
    return set(c - v for v in X) == set(X)


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
print("N1  LA CADENA, paso a paso, sobre las formas con [Phi]_top == 0")
print("=" * 116)
print("")
print("     P1  sigma_V(K) = K            (PROBADO el 12; se re-mide por si acaso)")
print("     P2  min S y max S estan en K  (= H*, la hipotesis aislada hoy por otra via)")
print("     P3  C = tau                   (lo que se quiere)")
print("     y la reduccion dice:  P1 y P2  =>  P3, por algebra.")
print("")
print("     t   r  | [Phi]_top=0 | e>2  |  P1   |  P2   |  P3   | P1&P2 pero no P3")
print("  " + "-" * 104)
TOT = Counter()
for (t, r, Wmax) in CFG:
    N = t + 2 * r
    n = ne2 = p1 = p2 = p3 = mal = 0
    for W in range(N - 1, Wmax + 1):
        for b in shapes_of_width(W, N):
            a = anatomia(b, t, r)
            if a is None or not a['top0']:
                continue
            n += 1
            ne2 += (a['e'] > 2)
            V = a['tau']
            q1 = simetrico(a['K'], V) if V is not None else False
            q2 = (a['S'][0] in a['K']) and (a['S'][-1] in a['K'])
            q3 = (a['C'] == V)
            p1 += q1
            p2 += q2
            p3 += q3
            mal += (q1 and q2 and not q3)
            TOT[(q1, q2, q3)] += 1
    print("    %2d  %2d  | %11d | %4d | %5d | %5d | %5d | %d"
          % (t, r, n, ne2, p1, p2, p3, mal))
    RES["N1_%d_%d" % (t, r)] = dict(n=n, e_mayor_2=ne2, P1=p1, P2=p2, P3=p3, contraejemplos=mal)
    sys.stdout.flush()
print("")
print("     contingencia global (P1, P2, P3): %s" % dict(sorted(TOT.items(), reverse=True)))
malos = sum(v for k, v in TOT.items() if k[0] and k[1] and not k[2])
print("     formas con P1 y P2 y SIN P3 : %d   %s"
      % (malos, "<-- la reduccion FALLA" if malos else "<-- la reduccion aguanta"))
RES['N1_total'] = {str(k): v for k, v in TOT.items()}

# ===================================================================== N2 ========================
print("")
print("=" * 116)
print("N2  H*  =>  C = tau,  con la libreria de controles (avisa si la tabla no puede fallar)")
print("=" * 116)
print("")
im = Implicacion("H* (extremos de S en K)", "C == tau")
c = Contingencia("H* (extremos de S en K)", "C == tau")
for (t, r, Wmax) in CFG:
    N = t + 2 * r
    for W in range(N - 1, Wmax + 1):
        for b in shapes_of_width(W, N):
            a = anatomia(b, t, r)
            if a is None or not a['top0'] or a['tau'] is None:
                continue
            q2 = (a['S'][0] in a['K']) and (a['S'][-1] in a['K'])
            q3 = (a['C'] == a['tau'])
            im.add(q2, q3, b)
            c.add(q2, q3, b)
im.informe()
print("")
c.informe()
RES['N2'] = dict(tabla={str(k): v for k, v in c.t.items()}, avisos=c.avisos())

# ===================================================================== N3 ========================
print("")
print("=" * 116)
print("N3  SEÑUELOS -- y el primero mide si el atajo era NECESARIO")
print("=" * 116)
print("")
print("  (a) sigma_V(g_com) = g_com : es lo que el .tex decia que haria falta.  Si NO se cumple")
print("      siempre y aun asi C = tau, el atajo por los extremos no es un adorno: es lo que salva")
print("      el caso e > 2.")
ca = Contingencia("sigma_V(g_com) = g_com", "C == tau")
cb = Contingencia("e == 2 (el caso ya cubierto)", "C == tau")
for (t, r, Wmax) in CFG:
    N = t + 2 * r
    for W in range(N - 1, Wmax + 1):
        for b in shapes_of_width(W, N):
            a = anatomia(b, t, r)
            if a is None or not a['top0'] or a['tau'] is None:
                continue
            q3 = (a['C'] == a['tau'])
            ca.add(simetrico(a['gcom'], a['tau']) if a['gcom'] else True, q3, b)
            cb.add(a['e'] == 2, q3, b)
ca.informe()
print("")
print("  (b) e == 2 : el caso que cor:ii YA cubre.  El denominador dice cuanto aporta lo nuevo.")
cb.informe()
RES['N3'] = dict(gcom={str(k): v for k, v in ca.t.items()}, e2={str(k): v for k, v in cb.t.items()})

# ===================================================================== N4 ========================
print("")
print("=" * 116)
print("N4  VEREDICTO")
print("=" * 116)
print("")
n = sum(RES["N1_%d_%d" % (t, r)]['n'] for (t, r, _) in CFG)
ne = sum(RES["N1_%d_%d" % (t, r)]['e_mayor_2'] for (t, r, _) in CFG)
ce = sum(RES["N1_%d_%d" % (t, r)]['contraejemplos'] for (t, r, _) in CFG)
p2 = sum(RES["N1_%d_%d" % (t, r)]['P2'] for (t, r, _) in CFG)
p3 = sum(RES["N1_%d_%d" % (t, r)]['P3'] for (t, r, _) in CFG)
print("     formas con [Phi]_top == 0 : %d,  de ellas con e > 2 : %d  (el caso que cor:ii NO cubre)"
      % (n, ne))
print("     H* (extremos en K)        : %d de %d" % (p2, n))
print("     C == tau                  : %d de %d" % (p3, n))
print("     P1 y P2 sin P3            : %d   %s" % (ce, "<-- LA REDUCCION FALLA" if ce else ""))
print("")
if ce == 0 and p2 == n:
    print("     LA REDUCCION AGUANTA Y H* SE CUMPLE SIEMPRE.  Entonces C = tau esta probado MODULO H*,")
    print("     y probar H* -- 'los dos maximizadores conservan los dos extremos de S' -- cierra el")
    print("     hueco del 12 y vuelve INCONDICIONAL el enunciado Phi_t == 0 => (ii).")
elif ce == 0:
    print("     LA REDUCCION AGUANTA, pero H* FALLA en %d formas: hay que mirar QUE pasa alli, porque" % (n - p2))
    print("     C = tau se cumple igual (%d de %d): hay una segunda razon, y N5 la busca." % (p3, n))
else:
    print("     LA REDUCCION FALLA en %d formas.  El atajo por los extremos no basta." % ce)
print("")
# ===================================================================== N5 ========================
print("=" * 116)
print("N5  LA FORMULACION GENERAL -- los dos RECORTES, que contiene a H* como caso 0 = 0")
print("=" * 116)
print("")
print("     sigma_V(K) = K vale SIEMPRE (P1), luego min K + max K = V.  Con K subconjunto de S:")
print("         C = tau   <=>   (max S - max K)  ==  (min K - min S)")
print("     o sea: el recorte de arriba y el de abajo se COMPENSAN.  H* es el caso 0 = 0.")
print("")
rec = Counter()
cc = Contingencia("recorte arriba == recorte abajo", "C == tau")
for (t, r, Wmax) in CFG:
    N = t + 2 * r
    for W in range(N - 1, Wmax + 1):
        for b in shapes_of_width(W, N):
            a = anatomia(b, t, r)
            if a is None or not a['top0'] or a['tau'] is None or not a['K']:
                continue
            arr = a['S'][-1] - a['K'][-1]
            aba = a['K'][0] - a['S'][0]
            rec[(arr, aba)] += 1
            cc.add(arr == aba, a['C'] == a['tau'], b)
print("     histograma de (recorte arriba, recorte abajo), los 12 mas frecuentes:")
for k, v in sorted(rec.items(), key=lambda x: -x[1])[:12]:
    print("        %-12s : %5d%s" % (str(k), v, "   <-- H*" if k == (0, 0) else ""))
print("     pares con recorte arriba != recorte abajo : %d de %d"
      % (sum(v for k, v in rec.items() if k[0] != k[1]), sum(rec.values())))
print("")
cc.informe()
RES['N5'] = dict(hist={str(k): v for k, v in rec.items()})
print("")
print("     ALCANCE: barrido con tope de reloj (W <= %s).  Todo enunciado de aqui es MEDIDO."
      % ", ".join(str(w) for (_, _, w) in CFG))
print("")
json.dump(RES, open(OUT_JSON, "w"), indent=1, default=str)
print("DONE")
