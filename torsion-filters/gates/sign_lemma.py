# -*- coding: utf-8 -*-
"""EL LEMA DE SIGNO w(sigma P) = -w(P), y el ataque que lo pudo tumbar.

LA FORMULA CERRADA.  Con F la reflexion x -> C-x sobre S (identidad fuera), alpha la reflexion
a -> C-a sobre E (identidad fuera), q = N-t, nu = #{x en S : 2x = C} y eta = #{a en E : 2a = C mod t}:

    w(sigma P) / w(P) = sgn(F) sgn(alpha) (-1)^{binom(q,2)}
                      = (-1)^{binom(q,2) + (|S|-nu)/2 + (e-eta)/2}          (**)

INDEPENDIENTE DE P.  Verificada aqui contra el w medido: 438/438 en trece configuraciones.

EL ATAQUE, que es la razon de ser de este guion.  De (**) sale que el cociente es -1 exactamente
cuando eta = 2, o sea cuando LAS DOS soluciones de 2*kappa = C (mod t) estan en E.  Pero nuestro
criterio solo pide un incremento = C, que fuerza UNA clase auto-pareada de tamano par.  Si existiera
una beta con el criterio y eta = 1, el cociente seria +1, no habria cancelacion, y el criterio
tendria falsos positivos fuera del rango barrido.  Se busca: eta=1 sale 0 de 438.

Y NO SE DEJA COMO MEDIDA.  eta = 2 se demuestra: F es involucion sobre S con nu fijos, luego
|S| = nu (mod 2); alpha lo es sobre E con eta fijos, luego e = eta (mod 2); y en el regimen
N = t+2r se tiene |S| = 2r+e, de donde e = nu (mod 2).  Un incremento = C da una clase auto-pareada
j0 en E de tamano PAR, que por C-simetrica de tamano par no contiene C/2.  Si nu=1, el C/2 vive en S
con residuo j0 o j1, no puede ser j0, luego j1 en E.  Si nu=0, e es par, y E = {j0} + parejas daria
e impar salvo que j1 este en E.  En los dos casos eta = 2.

Con eta=2 sale e par, nu=0, y el exponente de (**) es e-1, IMPAR: w(sigma P) = -w(P).

Authors: Carles Marin, Claude (AI assistant).  La formula cerrada viene de una consulta externa
gateada; la prueba de eta=2 y la verificacion son nuestras.
"""
import sys
sys.path.insert(0, r"E:\proyectos\Curiosity\research\orbit-pair\gates")
from fractions import Fraction as F
from math import comb
from second_stratum import setup, all_transversals
from criterion_control import betas, lam_of, value

PTS = [[F(3, 2), F(5, 3), F(7, 4), F(11, 6)], [F(5, 2), F(7, 3), F(9, 4), F(13, 5)],
       [F(4, 3), F(9, 5), F(11, 6), F(15, 7)]]


def is_zero(b, t, r):
    lam = lam_of(list(b))
    for p in PTS:
        if value(lam, t, r, p[:r]) != 0:
            return False
    return True


def anat(b, t, r):
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    S = sorted(x for k in E for x in Cd[k])
    C = S[0] + S[-1]
    inc = []
    for k in E:
        c = sorted(Cd[k], reverse=True)
        inc += [c[i] + c[i + 1] for i in range(len(c) - 1)]
    return cl, E, Cd, S, C, inc


CFG = [(4, 2, 17), (6, 2, 17), (8, 2, 17), (10, 2, 18), (12, 2, 18),
       (4, 3, 15), (6, 3, 15), (8, 3, 16), (2, 2, 16), (2, 3, 15), (2, 4, 15),
       (4, 4, 15), (6, 4, 16)]
print("=" * 104)
print("EL ATAQUE:  nuestro criterio pide S C-simetrico + algun incremento = C.")
print("  eta = #{kappa en E : 2kappa = C mod t}.  La formula externa dice w(sigma P)/w(P) = -1  <=>  eta = 2.")
print("  Si existe beta con el criterio y eta = 1, o el criterio tiene falsos positivos o la formula falla.")
print("  Se verifica ademas la formula cerrada (**) contra el w medido.")
print("=" * 104)
print("  t  r     criterio   eta=2   eta=1   eta=0 |  con eta=1: Phi=0  | formula (**) ok")
TE = {0: 0, 1: 0, 2: 0}
Z1 = 0
FOK = FTOT = 0
for (t, r, W) in CFG:
    n = e2 = e1 = e0 = z1 = fok = ftot = 0
    for b in betas(t, r, W):
        a = anat(b, t, r)
        if a is None:
            continue
        cl, E, Cd, S, C, inc = a
        if sorted(C - x for x in S) != S or C not in inc:
            continue
        n += 1
        sols = [k for k in range(t) if (2 * k - C) % t == 0]
        eta = sum(1 for k in sols if k in E)
        if eta == 2:
            e2 += 1
        elif eta == 1:
            e1 += 1
            z1 += is_zero(b, t, r)
        else:
            e0 += 1
        # formula (**)
        nu = sum(1 for x in S if 2 * x == C)
        pred = (-1) ** (comb(len(b) - t, 2) + (len(S) - nu) // 2 + (len(E) - eta) // 2)
        tr = all_transversals(b, cl, r, t)
        idx = {tuple(sorted(x[0].items())): x for x in tr}
        ok = True
        for (sel, T, w, deg) in tr:
            ssel = {}
            for k, v in sel.items():
                ssel[(C - k) % t if k in E else k] = (C - v) if k in E else v
            y = idx[tuple(sorted(ssel.items()))]
            if y[2] != pred * w:
                ok = False
        ftot += 1
        fok += ok
    for k, v in ((2, e2), (1, e1), (0, e0)):
        TE[k] += v
    Z1 += z1; FOK += fok; FTOT += ftot
    print("  %2d %2d %10d %7d %7d %7d | %16d  | %8d/%d" % (t, r, n, e2, e1, e0, z1, fok, ftot))
print()
print("  eta=2: %d    eta=1: %d    eta=0: %d" % (TE[2], TE[1], TE[0]))
print("  de las de eta=1 (donde la formula predice +1, o sea NO cancelacion) se anulan: %d" % Z1)
print("  formula cerrada (**): %d/%d" % (FOK, FTOT))
