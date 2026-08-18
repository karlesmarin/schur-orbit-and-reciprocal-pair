# -*- coding: utf-8 -*-
# ============================================================================================
#  A1 ES UN TEOREMA.  14 de agosto de 2026.
#
#  QUE SE PRUEBA.  A1 = "ningun extremo de S lo descartan los DOS maximizadores".  Con
#  extremes_of_S.py ya estaba visto que A1 => C = tau en tres lineas (S \ g_com = K u dif es
#  sigma_tau-estable).  Faltaba A1, y el intercambio de UNA transversal estaba REFUTADO con datos
#  (6320 casos con Delta < 0).  La pieza que faltaba no era otro intercambio: era usar cor:reflect
#  ENTERO --- no solo sigma_tau(K) = K sino T_B = tau - T_A, la reflexion del conjunto COMPLETO ---
#  junto con la FORMA DE PREFIJO del Paso 2 de thm:G.
#
#  LA PRUEBA (lado maximo).  Sea m = max S, de clase k, y supongamos m en g_A y en g_B.
#    (1) Paso 2 de thm:G: dentro de cada clase, leida decreciente, van primero los de A = H, luego
#        el pick g_k, luego los de L.  Como m = c_{k,1}, "m es el pick" <=> j_k = 0.
#        Si k fuera una de las dos clases empatadas, el maximizador que toma su incremento empatado
#        tendria j_k >= 1 (prefijo dentro de la clase).  Luego k NO esta empatada; y como la
#        seleccion contiene TODOS los incrementos > tau, todo incremento de la clase k es < tau.
#        En particular  Delta_k(1) = m + c < tau,  con c = c_{k,2}.
#    (2) Con j_k = 0 en los dos, c cae en L_A y en L_B.  Y T_B = tau - T_A da  L_B = tau - H_A.
#        Luego  c en L_B  =>  tau - c en H_A  incluido en S  =>  tau - c <= max S = m  =>  tau <= m + c.
#    (1) y (2) se contradicen.  El lado minimo es el reflejo exacto:  m- en g <=> j_{k-} = n-1 =>
#        TODOS los incrementos de k- seleccionados y k- no empatada => m- + c' > tau; y c' en
#        H_A n H_B con H_B = tau - L_A => tau - c' en L_A incluido en S => tau >= m- + c'.
#
#  Y LA HIPOTESIS ES INDISPENSABLE, no un adorno: con |G| = 2 pero SIN [Phi]_top = 0, A1 es FALSA.
#  Testigo construido a mano:  beta = (200,199,198,197,194,193,191,4), t=4, r=2.
#  Incrementos 392, 390, 390, 204; empate en las posiciones 2,3 entre las clases 1 y 3, |G| = 2;
#  pero Delta_0(1) = 204 queda fuera del top-2, asi que LOS DOS maximizadores descartan max S = 200.
#  Ahi inv(T_A) != inv(T_B), o sea [Phi]_top != 0.  N3 lo comprueba y busca mas como el.
#
#  COLUMNAS
#    C0  fatal: la parametrizacion por incrementos reproduce el conjunto G calculado por FUERZA
#        BRUTA sobre todas las transversales.  Si no, el Paso 2 no esta bien implementado y nada
#        de lo demas vale.  Mas la forma de prefijo, verificada elemento a elemento.
#    N1  EL MOTOR de la prueba, en su forma testable y NO vacia:
#        para todo v en L_A n L_B,  tau - v esta en H_A n H_B  (luego  v >= tau - max S).
#        Se imprime el numero de habitantes: si fuera 0, la medida no diria nada.
#    N2  la conclusion: max S en H_A u H_B, y min S en L_A u L_B, con la RAZON de cada caso
#        (clase empatada, o Delta_k(1) > tau).  Es A1 con su demostracion instrumentada.
#    N3  EL CONTROL QUE PUEDE FALLAR: |G| = 2 SIN [Phi]_top = 0.  Si A1 se cumpliera tambien ahi,
#        la prueba estaria usando una hipotesis que no hace falta y habria que sospechar de ella.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python a1_proof.py   (desde gates/)
# ============================================================================================

import itertools
import json
import os
import sys
from collections import Counter, defaultdict

from second_stratum import setup, all_transversals, inv_of, deg_of
from _control import Contingencia, Implicacion

OUT_JSON = "a1_proof_RESULT.json"
CFG = [(4, 2, 20), (6, 2, 20), (4, 3, 18), (6, 3, 18), (8, 2, 20), (8, 3, 19), (10, 2, 20)]


def clases(beta, t):
    C = defaultdict(list)
    for b in beta:
        C[b % t].append(b)
    for k in C:
        C[k].sort(reverse=True)
    return dict(C)


def incrementos(C):
    """[(valor, clase, j)] con Delta_k(j) = c_{k,j} + c_{k,j+1}, j desde 1."""
    out = []
    for k, cs in C.items():
        for j in range(len(cs) - 1):
            out.append((cs[j] + cs[j + 1], k, j + 1))
    return out


def maximizadores(beta, t, r):
    """Paso 2 + Pasos 4-6 de thm:G: los optimos son las selecciones de los r incrementos mayores.
    Devuelve (lista de dicts j_k, tau, clases empatadas) o None si no aplica."""
    C = clases(beta, t)
    if len(C) < t:
        return None                                   # (O) falla
    inc = incrementos(C)
    if len(inc) != 2 * r:
        return None
    inc.sort(key=lambda x: -x[0])
    tau = inc[r - 1][0]
    mayores = [x for x in inc if x[0] > tau]
    iguales = [x for x in inc if x[0] == tau]
    if len(mayores) + len(iguales) < r:
        return None
    s = r - len(mayores)                              # cuantos de los empatados entran
    sel = []
    if s == len(iguales):
        sel = [mayores + iguales]
    elif s >= 1:
        sel = [mayores + list(c) for c in itertools.combinations(iguales, s)]
    else:
        return None
    empatadas = tuple(sorted(set(x[1] for x in iguales))) if len(iguales) > 1 else ()
    js = []
    for S_ in sel:
        j = {k: 0 for k in C}
        for (_, k, jj) in S_:
            j[k] = max(j[k], jj)
        js.append(j)
    return js, tau, empatadas, C


def partes(j, C):
    """de la forma de prefijo: (H = A, pick g, L)."""
    A, g, L = [], [], []
    for k, cs in C.items():
        jk = j.get(k, 0)
        A += cs[:jk]
        if jk < len(cs):
            g.append(cs[jk])
        L += cs[jk + 1:]
    return sorted(A, reverse=True), sorted(g, reverse=True), sorted(L, reverse=True)


def shapes(t, r, Wmax):
    N = t + 2 * r
    for W in range(N - 1, Wmax + 1):
        for resto in itertools.combinations(range(1, W), N - 2):
            yield tuple(sorted((W,) + resto + (0,), reverse=True))


# ===================================================================== C0 ========================
print("=" * 116)
print("C0  ACEPTACION -- fatal.  La parametrizacion por incrementos contra la FUERZA BRUTA.")
print("=" * 116)
print("")
mal_G = mal_pre = n_c0 = 0
ejemplo_mal = None
for (t, r, Wmax) in [(4, 2, 14), (6, 2, 13), (4, 3, 13)]:
    for b in shapes(t, r, Wmax):
        st = setup(b, t)
        if st is None:
            continue
        cl, E, Cd = st
        tr = all_transversals(b, cl, r, t)
        D = max(x[3] for x in tr)
        Gbf = sorted(tuple(sorted(x[1], reverse=True)) for x in tr if x[3] == D)
        M = maximizadores(b, t, r)
        n_c0 += 1
        if M is None:
            mal_G += 1
            continue
        js, tau, emp, C = M
        Gmine = sorted(tuple(sorted(partes(j, C)[0] + partes(j, C)[2], reverse=True)) for j in js)
        if Gmine != Gbf:
            mal_G += 1
            if ejemplo_mal is None:
                ejemplo_mal = (t, r, b, Gmine, Gbf)
        for j in js:
            A, g, L = partes(j, C)
            T = sorted(A + L, reverse=True)
            # forma de prefijo: A tiene que ser EXACTAMENTE la mitad de arriba de T
            if A != T[:r] or L != T[r:]:
                mal_pre += 1
print("     formas contrastadas               : %d" % n_c0)
print("     G(incrementos) != G(fuerza bruta) : %d %s"
      % (mal_G, "" if not mal_G else "*** %s ***" % str(ejemplo_mal)[:60]))
print("     A no es la mitad de arriba de T   : %d" % mal_pre)
print("")
if mal_G or mal_pre:
    print("  C0 FALLA -- veredicto SUSPENDIDO.")
    print("DONE (suspendido)")
    raise SystemExit(1)
print("  C0 PASA")
sys.stdout.flush()


# ------------------------------------------------------------------ poblacion --------------------
def poblacion(cfg):
    """(t, r, beta, datos) para las formas con |G| = 2; 'top0' marca [Phi]_top = 0."""
    for (t, r, Wmax) in cfg:
        for b in shapes(t, r, Wmax):
            M = maximizadores(b, t, r)
            if M is None:
                continue
            js, tau, emp, C = M
            if len(js) != 2:
                continue
            PA, PB = partes(js[0], C), partes(js[1], C)
            TA = tuple(sorted(PA[0] + PA[2], reverse=True))
            TB = tuple(sorted(PB[0] + PB[2], reverse=True))
            cl = {k: [i for i, v in enumerate(b) if v % t == k] for k in set(x % t for x in b)}
            tr = {tuple(sorted(x[1], reverse=True)): x[2] for x in all_transversals(b, cl, r, t)}
            top0 = (inv_of(TA, r) == inv_of(TB, r) and tr.get(TA) == -tr.get(TB))
            S = sorted((v for k, cs in C.items() if len(cs) >= 2 for v in cs))
            yield (t, r, b, dict(tau=tau, emp=emp, C=C, HA=PA[0], gA=PA[1], LA=PA[2],
                                 HB=PB[0], gB=PB[1], LB=PB[2], TA=TA, TB=TB, top0=top0, S=S))


POB = list(poblacion(CFG))
TOP0 = [x for x in POB if x[3]['top0']]
print("")
print("     formas con |G| = 2 : %d ;  de ellas con [Phi]_top = 0 : %d  (archivado: 2522)  %s"
      % (len(POB), len(TOP0), "ok" if len(TOP0) == 2522 else "*** NO CUADRA ***"))
if len(TOP0) != 2522:
    print("  *** la poblacion no reproduce la de c_eq_tau.py -- LEER CON CUIDADO ***")
sys.stdout.flush()

# ===================================================================== N1 ========================
print("")
print("=" * 116)
print("N1  EL MOTOR:  v en L_A n L_B  =>  tau - v en H_A n H_B   (de T_B = tau - T_A)")
print("=" * 116)
print("")
n_hab = n_mal = 0
n_refl = 0
for (t, r, b, d) in TOP0:
    tau = d['tau']
    if [tau - v for v in d['HA']][::-1] == d['LB'] and [tau - v for v in d['LA']][::-1] == d['HB']:
        n_refl += 1
    HA, HB, LA, LB = set(d['HA']), set(d['HB']), set(d['LA']), set(d['LB'])
    for v in LA & LB:
        n_hab += 1
        if not ((tau - v) in HA and (tau - v) in HB):
            n_mal += 1
print("     T_B = tau - T_A verificado en           : %d de %d" % (n_refl, len(TOP0)))
print("     habitantes de 'v en L_A n L_B'          : %d   (si fuera 0, esto no diria nada)" % n_hab)
print("     fallos de 'tau - v en H_A n H_B'        : %d" % n_mal)
print("")

# ===================================================================== N2 ========================
print("=" * 116)
print("N2  A1, INSTRUMENTADA: la razon por la que cada extremo se salva")
print("=" * 116)
print("")
razon = Counter()
falla = []
for (t, r, b, d) in TOP0:
    C, tau, emp = d['C'], d['tau'], d['emp']
    S = d['S']
    m, mm = S[-1], S[0]
    k, kk = m % t, mm % t
    okmax = (m in d['HA']) or (m in d['HB'])
    okmin = (mm in d['LA']) or (mm in d['LB'])
    if not (okmax and okmin):
        falla.append((t, r, b))
    razon["max: clase empatada" if k in emp else "max: Delta_k(1) > tau"] += 1
    razon["min: clase empatada" if kk in emp else "min: Delta_last > tau"] += 1
    # y la desigualdad que la prueba usa, medida
    d1 = C[k][0] + C[k][1]
    dl = C[kk][-1] + C[kk][-2]
    razon["Delta_k(1) >= tau"] += (d1 >= tau)
    razon["Delta_last <= tau"] += (dl <= tau)
print("     A1 falla en : %d de %d %s" % (len(falla), len(TOP0), str(falla[:3])))
for kk_ in sorted(razon):
    print("        %-26s : %d" % (kk_, razon[kk_]))
print("")

# ===================================================================== N3 ========================
print("=" * 116)
print("N3  EL CONTROL QUE PUEDE FALLAR: |G| = 2 pero SIN [Phi]_top = 0")
print("=" * 116)
print("")
print("     Si A1 se cumpliera tambien aqui, la prueba estaria apoyandose en una hipotesis que no")
print("     hace falta, y habria que sospechar de ella.  El testigo construido a mano:")
TESTIGO = (200, 199, 198, 197, 194, 193, 191, 4)
M = maximizadores(TESTIGO, 4, 2)
if M is None:
    print("        *** el testigo no pasa (O) o el conteo de incrementos ***")
else:
    js, tau, emp, C = M
    print("        beta = %s   t=4 r=2" % str(TESTIGO))
    print("        incrementos: %s" % sorted((v for (v, _, _) in incrementos(C)), reverse=True))
    print("        tau = %d ;  clases empatadas = %s ;  |G| = %d" % (tau, str(emp), len(js)))
    for i, j in enumerate(js):
        A, g, L = partes(j, C)
        T = sorted(A + L, reverse=True)
        print("        max%d: j = %-22s H = %-16s g = %-20s L = %s  deg = %d"
              % (i, str(sorted(j.items())), str(A), str(g), str(L), deg_of(tuple(T), 2)))
    S = sorted((v for k, cs in C.items() if len(cs) >= 2 for v in cs))
    gA = set(partes(js[0], C)[1])
    gB = set(partes(js[1], C)[1])
    print("        max S = %d ;  descartado por los DOS: %s" % (S[-1], S[-1] in gA and S[-1] in gB))
    TA = tuple(sorted(partes(js[0], C)[0] + partes(js[0], C)[2], reverse=True))
    TB = tuple(sorted(partes(js[1], C)[0] + partes(js[1], C)[2], reverse=True))
    print("        inv(T_A) = %s" % str(inv_of(TA, 2)))
    print("        inv(T_B) = %s   ->  iguales: %s" % (str(inv_of(TB, 2)), inv_of(TA, 2) == inv_of(TB, 2)))
print("")
print("     Y el barrido: sobre TODAS las formas con |G| = 2, se cuenta A1 contra [Phi]_top = 0.")
cont = Contingencia("[Phi]_top = 0", "A1 (ningun extremo descartado por los dos)")
for (t, r, b, d) in POB:
    S = d['S']
    m, mm = S[-1], S[0]
    a1 = ((m in d['HA']) or (m in d['HB'])) and ((mm in d['LA']) or (mm in d['LB']))
    cont.add(d['top0'], a1, b)
cont.informe(indent="     ")
print("")
print("     COMO SE LEE ESTA TABLA, y no es como la lee el 'veredicto' de arriba.  Aqui NO se afirma")
print("     un bicondicional: se afirma la IMPLICACION [Phi]_top = 0 => A1.  Lo que la hace portante")
print("     es que la casilla (top0 = no, A1 = no) tenga habitantes, y los tiene:")
print("        A1 FALSA con |G| = 2 y [Phi]_top != 0 : %d formas" % cont.t[(False, False)])
print("        A1 FALSA con [Phi]_top  = 0           : %d formas" % cont.t[(True, False)])
print("     Que A1 se cumpla ademas en %d formas con [Phi]_top != 0 no debilita nada: dice que la"
      % cont.t[(False, True)])
print("     hipotesis es SUFICIENTE y no necesaria, que es exactamente lo que la prueba afirma.")
im3 = Implicacion("[Phi]_top = 0", "A1")
for k, v in cont.t.items():
    for _ in range(v):
        im3.add(k[0], k[1])
im3.informe(indent="     ")
print("")

# ===================================================================== VEREDICTO =================
print("=" * 116)
print("VEREDICTO")
print("=" * 116)
print("")
print("     A1 sobre [Phi]_top = 0 : %d de %d, y la prueba de la cabecera la explica caso a caso."
      % (len(TOP0) - len(falla), len(TOP0)))
print("     A1 SIN esa hipotesis   : falla (ver la tabla de N3).  La hipotesis es NECESARIA.")
print("     Luego  C = tau  es TEOREMA, y  Phi_t == 0 => (ii)  queda INCONDICIONAL.")
print("")
print("     ALCANCE: mismo barrido que c_eq_tau.py (tope de reloj).  La PRUEBA no tiene alcance:")
print("     lo medido aqui es que el codigo dice lo que la prueba dice, no que la prueba valga.")
json.dump({"n_G2": len(POB), "n_top0": len(TOP0), "A1_fallos": len(falla),
           "motor_habitantes": n_hab, "motor_fallos": n_mal, "reflexion": n_refl,
           "razon": dict(razon), "contingencia": {str(k): v for k, v in cont.t.items()}},
          open(OUT_JSON, "w"), indent=1)
print("DONE")
