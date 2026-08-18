# -*- coding: utf-8 -*-
# ============================================================================================
#  AUDITORIA DE LAS FORMULAS DE LA v2, una por una.  14 de agosto de 2026.
#
#  POR QUE.  Los bloques 1-3 de la v2 anaden una cadena de identidades con signos delicados
#  (constantes de Laplace, (-1)^{binom r 2}, el convenio de V, la normalizacion de alpha*).  Que el
#  .tex COMPILE no dice nada de si son ciertas.  Aqui cada formula DISPLAYED de los tres bloques se
#  comprueba con aritmetica entera exacta, y cada una lleva su propia linea de salida.
#
#  EL CONVENIO DE V, leido de lem:L1 y NO supuesto:  V = prod_{k<k'} (zeta^{k'} - zeta^k).
#  A t=2 eso es  zeta^1 - zeta^0 = -1 - 1 = -2.   (lem:L2 lo confirma por otra via:
#  prod_{k<k'}(zeta^k - zeta^{k'}) = (-1)^{binom t 2} V, y a t=2 da +2 = -V.)
#
#  QUE SE COMPRUEBA, y a que alcance:
#    F1  lem:L2  det(x_i^{N-j}) = (-1)^{binom t2} V (z^t-1)(z^{-t}-1)(z-z^{-1})      t=2, r=1
#    F2  lem:L1  M_S = (-1)^{inv(b_S)} V                                            t=2, todo r
#    F3  eq:laplacegen  det(x_i^{beta_j}) = (-1)^{binom{t+1}{2}} V sum_S w(S) A(S^c) t=2, todo r
#    F4  eq:deggen  deg(T) es el grado total maximo de A(T)                          todo r
#    F5  [A(T)]_max = eps * a_H(z) a_L(z^{-1}), eps dependiendo SOLO de r            todo r
#    F6  lem:reflgen  P(c-T) = P(T)                                                  todo r
#    F7  eq:dictgen  el diccionario, con su (-1)^{binom r2} y su (prod z)^{h_r-l_1}  todo r
#    F8  eq:incrgen + eq:congrgen  Delta_i(k) = 2i (mod t)                            todo t
#    F9  eq:topgen  [det]_{D_1} = (-1)^{binom{t+1}{2}} V sum_{g in G} w(g) P(T_g)     t=2, todo r
#    F10 la traduccion final del teorema en t=2:
#          beta_j + beta_{N+1-j} = C  para todo j   <=>   lambda_j + lambda_{N+1-j} = w
#          C par  <=>  w impar        (con N = 2r+2 par)
#    F11 el enunciado de rem:rankone contra la expansion simplectica directa:
#          r=1 -> +- caracter GENUINO ;  r>=2 -> propiamente VIRTUAL
#
#  NO SE COMPRUEBA AQUI, y se dice: nada con t>=3 que necesite Z[zeta] (F1-F3 y F9 van a t=2, que es
#  el caso del teorema); y thm:PvW, que es externo.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python v2_formulas.py   (desde gates/)
# ============================================================================================

import itertools
import json
import os
import sys
from collections import Counter
from itertools import combinations, permutations

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "folding_t2.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0 =")[0]
assert "def alt(" in _head and "def ldiv(" in _head, "folding_t2.py cambio de forma"
_ns = {"__name__": "folding_t2_preamble", "__file__": SRC}
exec(compile(_head, SRC, "exec"), _ns)
padd, pscale, psub, ldiv = _ns["padd"], _ns["pscale"], _ns["psub"], _ns["ldiv"]
alt, det_pares, perms_signed = _ns["alt"], _ns["det_pares"], _ns["perms_signed"]
restriccion, to_sp, signo_de = _ns["restriccion"], _ns["to_sp"], _ns["signo_de"]

OUT_JSON = "v2_formulas_RESULT.json"
RES = {}
FALLOS = 0


def binom2(n):
    return n * (n - 1) // 2


def a_alt(X, r, inverso=False):
    """a_X(z) = det(z_j^{X_i})_{r x r}; con inverso=True, a_X(z^{-1})."""
    out = {}
    s = -1 if inverso else 1
    for p, sg in perms_signed(r):
        e = [0] * r
        for i in range(r):
            e[p[i]] = s * X[i]
        k = tuple(e)
        v = out.get(k, 0) + sg
        if v:
            out[k] = v
        elif k in out:
            del out[k]
    return out


def pmulL(a, b):
    o = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            k = tuple(x + y for x, y in zip(ea, eb))
            v = o.get(k, 0) + ca * cb
            if v:
                o[k] = v
            elif k in o:
                del o[k]
    return o


def P_de(T, r):
    H, L = list(T[:r]), list(T[r:])
    return pmulL(a_alt(H, r), a_alt(L, r, inverso=True))


def schur(nu, r):
    """s_nu(z_1..z_r) por el bialternante, como Laurent (de hecho polinomio)."""
    delta = [r - 1 - i for i in range(r)]
    num = a_alt([nu[i] + delta[i] for i in range(r)], r)
    den = a_alt(delta, r)
    return ldiv(num, den)


def inv_word(w):
    return sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])


def betas(N, top, n=400, seed=12345):
    """N-subconjuntos decrecientes de {0..top}, deterministas y variados."""
    out, x = [], seed
    vistos = set()
    while len(out) < n:
        x = (1103515245 * x + 12345) % (2 ** 31)
        s = set()
        y = x
        while len(s) < N:
            y = (1103515245 * y + 12345) % (2 ** 31)
            s.add(y % (top + 1))
        b = tuple(sorted(s, reverse=True))
        if b not in vistos:
            vistos.add(b)
            out.append(b)
    return out


def linea(nombre, malos, n, extra=""):
    global FALLOS
    FALLOS += malos
    print("  %-6s %-58s %6d fallos de %-6d %s"
          % (nombre, "", malos, n, "ok" if malos == 0 else "*** FALLA ***"))
    RES[nombre] = {"malos": malos, "n": n}


print("=" * 116)
print("AUDITORIA DE FORMULAS -- bloques 1, 2 y 3 de la v2.  Aritmetica entera exacta.")
print("=" * 116)
print("")

# ---------------------------------------------------------------- F1  lem:L2 --------------------
print("  F1   lem:L2   det(x_i^{N-j}) = (-1)^{binom t2} V (z^t-1)(z^{-t}-1)(z-z^{-1})   [t=2, r=1]")
t, r = 2, 1
N = t + 2 * r
V = -2                                   # prod_{k<k'} (zeta^{k'} - zeta^k) a t=2
izq = alt([N - 1 - i for i in range(N)], r)
# (z^2-1)(z^{-2}-1)(z-z^{-1})
f1 = {(2,): 1, (0,): -1}
f2 = {(-2,): 1, (0,): -1}
f3 = {(1,): 1, (-1,): -1}
der = pscale(pmulL(pmulL(f1, f2), f3), ((-1) ** binom2(t)) * V)
malos = 0 if izq == der else 1
print("        izquierda = %s" % sorted(izq.items(), reverse=True))
print("        derecha   = %s" % sorted(der.items(), reverse=True))
linea("F1", malos, 1)
print("")

# ---------------------------------------------------------------- F2  lem:L1 --------------------
print("  F2   lem:L1   M_S = (-1)^{inv(b_S)} V   si los residuos son distintos, 0 si no   [t=2]")
malos = n = 0
for r in (1, 2, 3):
    N = 2 * r + 2
    for b in betas(N, 3 * N, 60, seed=7 + r):
        for S in combinations(range(N), 2):          # S = t columnas, t = 2
            bS = [b[j] % 2 for j in S]
            M = 0
            # M_S = det( zeta^{k beta_j} )_{k=0,1 ; j in S} = det [[1,1],[(-1)^b1,(-1)^b2]]
            M = (-1) ** (b[S[1]] % 2) - (-1) ** (b[S[0]] % 2)
            esp = 0 if bS[0] == bS[1] else ((-1) ** inv_word(bS)) * V
            n += 1
            malos += (M != esp)
linea("F2", malos, n)
print("")

# ---------------------------------------------------------------- F3  eq:laplacegen -------------
print("  F3   eq:laplacegen   det(x_i^{beta_j}) = (-1)^{binom{t+1}{2}} V sum_S w(S) A(S^c)   [t=2]")
print("        con w(S) = (-1)^{sum_{j in S} j + inv(b_S)} y j EMPEZANDO EN 1 (columnas del paper)")
malos = n = 0
for r in (1, 2, 3):
    N = 2 * r + 2
    for b in betas(N, 3 * N, 40, seed=101 + r):
        izq = alt(list(b), r)
        acc = {}
        for S in combinations(range(N), 2):
            bS = [b[j] % 2 for j in S]
            if bS[0] == bS[1]:
                continue
            w = (-1) ** (sum(j + 1 for j in S) + inv_word(bS))
            resto = [b[j] for j in range(N) if j not in S]
            acc = padd(acc, pscale(det_pares(resto, r), w))
        der = pscale(acc, ((-1) ** binom2(t + 1)) * V)
        n += 1
        malos += (izq != der)
linea("F3", malos, n)
print("")

# ---------------------------------------------------------------- F4, F5, F6, F7 ----------------
print("  F4   eq:deggen   deg(T) = sum_{a<=r} u_a - sum_{a>r} u_a es el grado total MAXIMO de A(T)")
print("  F5   [A(T)]_max = eps * a_H(z) a_L(z^{-1}), con eps dependiendo SOLO de r")
print("  F6   lem:reflgen   P(c-T) = P(T)")
print("  F7   eq:dictgen   P(T) = (-1)^{binom r2} (prod z)^{h_r-l_1} a_delta(z)^2 s_atil s_astar")
m4 = m5 = m6 = m7 = n45 = 0
eps_por_r = {}
for r in (1, 2, 3):
    for T in [tuple(sorted(s, reverse=True)) for s in
              (set(x) for x in (betas(2 * r, 6 * r + 6, 60, seed=31 + r)))]:
        if len(T) != 2 * r:
            continue
        n45 += 1
        A = det_pares(list(T), r)
        H, L = list(T[:r]), list(T[r:])
        d = sum(H) - sum(L)
        # F4: grado total maximo
        if A:
            gm = max(sum(e) for e in A)
            m4 += (gm != d)
        # F5: la componente de grado maximo
        Amax = {e: c for e, c in A.items() if sum(e) == d}
        P = P_de(T, r)
        eps = None
        if Amax and P:
            e0 = max(P)
            if e0 in Amax and P[e0] and Amax[e0] % P[e0] == 0:
                eps = Amax[e0] // P[e0]
        if eps is None or Amax != pscale(P, eps) or eps not in (1, -1):
            m5 += 1
        else:
            eps_por_r.setdefault(r, set()).add(eps)
        # F6: reflexion
        for c in (0, 7, 23):
            Tc = tuple(sorted((c - u for u in T), reverse=True))
            if P_de(Tc, r) != P:
                m6 += 1
        # F7: el diccionario
        delta = [r - 1 - i for i in range(r)]
        alpha = [H[i] - delta[i] for i in range(r)]
        atil = [a - alpha[-1] for a in alpha]
        Lst = [L[0] - L[r - 1 - i] for i in range(r)]
        astar = [Lst[i] - delta[i] for i in range(r)]
        if astar[-1] != 0:
            m7 += 1
            continue
        ad = a_alt(delta, r)
        mono = {tuple([H[-1] - L[0]] * r): 1}
        der = pscale(pmulL(pmulL(pmulL(ad, ad), mono),
                           pmulL(schur(atil, r), schur(astar, r))), (-1) ** binom2(r))
        m7 += (der != P)
linea("F4", m4, n45)
linea("F5", m5, n45)
print("        eps observado por r: %s   (tiene que ser UN solo valor por r)"
      % {k: sorted(v) for k, v in eps_por_r.items()})
FALLOS += sum(1 for v in eps_por_r.values() if len(v) != 1)
linea("F6", m6, n45 * 3)
linea("F7", m7, n45)
print("")

# ---------------------------------------------------------------- F8 ----------------------------
print("  F8   eq:incrgen + eq:congrgen   Delta_i(k) = c_{i,k} + c_{i,k+1} = 2i (mod t)")
malos = n = 0
for t2 in (2, 3, 4, 5, 6, 8):
    for r in (1, 2, 3):
        N = t2 + 2 * r
        for b in betas(N, 4 * N, 25, seed=t2 * 100 + r):
            cl = {}
            for v in b:
                cl.setdefault(v % t2, []).append(v)
            for i, cs in cl.items():
                cs.sort(reverse=True)
                for k in range(len(cs) - 1):
                    n += 1
                    malos += ((cs[k] + cs[k + 1] - 2 * i) % t2 != 0)
linea("F8", malos, n)
print("")

# ---------------------------------------------------------------- F9 ----------------------------
print("  F9   eq:topgen   [det]_{D_1} = (-1)^{binom{t+1}{2}} V sum_{g in G} w(g) P(T_g)   [t=2]")
malos = n = 0
for r in (1, 2, 3):
    N = 2 * r + 2
    for b in betas(N, 3 * N, 40, seed=555 + r):
        terminos = []
        for S in combinations(range(N), 2):
            bS = [b[j] % 2 for j in S]
            if bS[0] == bS[1]:
                continue
            w = (-1) ** (sum(j + 1 for j in S) + inv_word(bS))
            T = tuple(b[j] for j in range(N) if j not in S)
            terminos.append((sum(T[:r]) - sum(T[r:]), w, T))
        if not terminos:
            continue
        D1 = max(x[0] for x in terminos)
        izq = alt(list(b), r)
        izq = {e: c for e, c in izq.items() if sum(e) == D1}
        acc = {}
        for (d, w, T) in terminos:
            if d == D1:
                acc = padd(acc, pscale(P_de(T, r), w))
        der = pscale(acc, ((-1) ** (binom2(t + 1) + binom2(r))) * V)
        n += 1
        malos += (izq != der)
linea("F9", malos, n)
print("        [DEFECTO CAZADO AQUI EL 14-ago: eq:topgen se escribio SIN el (-1)^{binom r2} de")
print("         eq:Amaxgen, y fallaba 78 de 116 -- 0 en r=1, 38 de 40 en r=2, 40 de 40 en r=3.")
print("         eps_r = (-1)^{binom r2} verificado ademas en r=4 y r=5, que es donde podia romperse.")
print("         Corregido en v2_new_subsections.tex; sin el, la identidad es falsa para r>=2.]")
print("")

# ---------------------------------------------------------------- F10 ---------------------------
print("  F10  la traduccion del teorema en t=2:  beta simetrico <=> lambda auto-complementario,")
print("       y  C par <=> w impar  (N = 2r+2 par).  Aritmetica pura, sobre particiones reales.")
malos = n = 0
for r in (1, 2, 3):
    N = 2 * r + 2
    for lam in itertools.combinations_with_replacement(range(6, -1, -1), N):
        lam = list(lam)
        beta = [lam[j] + (N - 1 - j) for j in range(N)]
        C = beta[0] + beta[-1]
        sim = all(beta[j] + beta[N - 1 - j] == C for j in range(N))
        w = lam[0] + lam[-1]
        auto = all(lam[j] + lam[N - 1 - j] == w for j in range(N))
        n += 1
        malos += (sim != auto)
        if sim:
            malos += ((C % 2 == 0) != (w % 2 == 1))
            malos += (C != w + N - 1)
linea("F10", malos, n)
print("")

# ---------------------------------------------------------------- F11 ---------------------------
print("  F11  rem:rankone, PUBLICADO:  r=1 -> +- caracter GENUINO ;  r>=2 -> propiamente VIRTUAL.")
print("       Se contrasta contra la expansion en la base sp_mu, que es una via independiente.")
print("       Y los SIETE testigos que la propia observacion lista tienen que salir virtuales.")
malos = n = 0
mezcl = Counter()
for r, LM in ((1, 8), (2, 5), (3, 3)):
    N = 2 * r + 2
    for lam in itertools.combinations_with_replacement(range(LM, -1, -1), N):
        P, _ = restriccion(list(lam), r)
        if not P:
            continue
        n += 1
        s = signo_de(to_sp(P, r))
        mezcl[(r, s == 0)] += 1
        if r == 1 and s == 0:
            malos += 1                                  # r=1 NO puede ser virtual
TESTIGOS = [(2, (5, 4, 3, 1, 0, 0)), (2, (5, 5, 4, 2, 1, 0)), (2, (10, 5, 3, 1, 0, 0)),
            (2, (6, 5, 4, 2, 1, 1)), (3, (2, 2, 2, 1, 1, 0, 0, 0)), (3, (5, 5, 4, 1, 1, 0, 0, 0)),
            (3, (3, 3, 3, 2, 2, 1, 1, 1))]
mt = 0
for r, lam in TESTIGOS:
    P, _ = restriccion(list(lam), r)
    mt += (not P) or (signo_de(to_sp(P, r)) != 0)
print("        formas no nulas por r (r, virtual): %s" % dict(mezcl))
print("        los 7 testigos de rem:rankone que NO salen virtuales: %d" % mt)
linea("F11", malos + mt, n + 7)
print("")

# ---------------------------------------------------------------- VEREDICTO ---------------------
print("=" * 116)
print("VEREDICTO")
print("=" * 116)
print("")
print("     fallos totales: %d" % FALLOS)
print("     %s" % ("TODAS LAS FORMULAS DE LOS TRES BLOQUES PASAN."
                  if FALLOS == 0 else "*** HAY FORMULAS MAL -- no insertar hasta arreglarlas ***"))
print("")
print("     NO AUDITADO AQUI, y se dice: t >= 3 en F1, F2, F3 y F9 (harian falta enteros")
print("     ciclotomicos; F1-F3 y F9 van a t=2, que es el caso del teorema de la v2), y thm:PvW,")
print("     que es un input EXTERNO y no se comprueba: se cita.")
json.dump({"fallos": FALLOS, "detalle": RES}, open(OUT_JSON, "w"), indent=1)
print("DONE")
