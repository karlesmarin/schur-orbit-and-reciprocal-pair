# -*- coding: utf-8 -*-
# LA IDENTIDAD DE HISTOGRAMAS QUE EXPLICA LA FACTORIZACION RESIDUAL.   9 de septiembre de 2026.
#
# POR QUE.  residual_factorizacion.py MIDE que  |chi| = kappa . D . chi_eta^H  en 1814 fichas.  Eso es
# una formula comprobada, no un teorema.  La pieza que dice POR QUE es la identidad de
# histogramas: al comparar ell con rho, las vueltas completas se cancelan y lo que
# queda es EL HISTOGRAMA DE RAICES DE OTRO GRUPO, mas pequeno y regular.
#
#     R_ell^{C_m} - R_rho^{C_m}  =  R_s^H - R_{rho_H}^H       en Z[X]/(X^q - 1),   s = eta + rho_H
#
# Y la Proposicion 8.5 del articulo ya da la forma de R:
#     R_x^{C_m} = (P_x^2 + P_x(X^2))/2 - m        con  P_x = sum_i (X^{x_i} + X^{-x_i})
#     R_x^{D_r} = (P_x^2 - P_x(X^2))/2 - r        <- el mismo objeto con el signo cambiado
# El signo del segundo sumando es LA UNICA diferencia entre C y D: el termino P(X^2) son las
# raices largas 2e_i, que D no tiene.  Por eso el caso v<=a sale de tipo C y el v>a de tipo D.
#
# QUE COMPRUEBA
#   H0  la forma de la Prop. 8.5: R_x reproduce el histograma contado a mano sobre las raices
#   H1  la identidad de residuos, sobre los supervivientes con perfil minimo, q cualquiera
#   H2  el reparto por tipo: v<=a da C, v>a da D --- y no al reves
#   H3  SENUELOS: tres maneras de romperla
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python -u residual_histogramas.py > residual_histogramas_OUT.txt 2>&1

import os
import sys
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from mpmath import mpf
from teorema_II import perfil, D_de
from galois import chi_lim, dominantes

fallos = []


def P_de(x, q):
    """P_x = sum_i (X^{x_i} + X^{-x_i})  en Z[X]/(X^q - 1), como vector de q coeficientes."""
    c = [0] * q
    for t in x:
        c[t % q] += 1
        c[(-t) % q] += 1
    return c


def cuad(c, q):
    """el cuadrado en Z[X]/(X^q-1)."""
    out = [0] * q
    for i, a in enumerate(c):
        if not a:
            continue
        for j, b in enumerate(c):
            if b:
                out[(i + j) % q] += a * b
    return out


def sub2(c, q):
    """P(X^2): duplica los exponentes."""
    out = [0] * q
    for i, a in enumerate(c):
        out[(2 * i) % q] += a
    return out


def R_de(x, q, tipo):
    """R_x segun la Prop. 8.5.  tipo 'C' lleva +P(X^2) (raices largas); 'D' lleva -P(X^2)."""
    n = len(x)
    P = P_de(x, q)
    A, B = cuad(P, q), sub2(P, q)
    s = 1 if tipo == "C" else -1
    out = [(A[i] + s * B[i]) // 2 for i in range(q)]
    assert all((A[i] + s * B[i]) % 2 == 0 for i in range(q)), "R no sale entero"
    out[0] -= n
    return out


def hist_a_mano(x, q, tipo):
    """el mismo histograma, contando <x,alpha> sobre TODAS las raices, sin formula."""
    n = len(x)
    c = [0] * q
    for i in range(n):
        for j in range(i + 1, n):
            for v in (x[i] - x[j], x[j] - x[i], x[i] + x[j], -x[i] - x[j]):
                c[v % q] += 1
        if tipo == "C":
            for v in (2 * x[i], -2 * x[i]):
                c[v % q] += 1
    return c


# ------------------------------------------------------------------- H0
print("=" * 96)
print("H0 -- LA FORMA DE LA PROP. 8.5, CONTRA EL HISTOGRAMA CONTADO A MANO")
mal = viv = 0
for q in range(2, 15):
    for n in range(1, 7):
        for x in dominantes(n, 6):
            for tipo in ("C", "D"):
                if tipo == "D" and n < 2:
                    continue
                viv += 1
                if R_de(list(x), q, tipo) != hist_a_mano(list(x), q, tipo):
                    mal += 1
print("   %d casos (q, rango, x, tipo)" % viv)
print(("   OK   " if not mal else "  FALLA ")
      + "H0: (P^2 +- P(X^2))/2 - n es el histograma de <x,alpha> sobre las raices")
if mal:
    fallos.append("H0")


# --------------------------------------------------- el perfil residual
def datos(lam, q):
    m = len(lam)
    u, v = divmod(m, q)
    a = (q - 1) // 2
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    P = perfil(ell, q)
    fijas = {0} | ({q // 2} if q % 2 == 0 else set())
    oc = {c: len(P.get(c, [])) for c in range(0, q // 2 + 1)}
    if v <= a:
        S = sorted([c for c in range(1, a + 1) if oc.get(c, 0) == 2 * u + 1], reverse=True)
        for c in range(0, q // 2 + 1):
            if oc.get(c, 0) != (u if c in fijas else 2 * u + (1 if c in S else 0)):
                return None
        if len(S) != v:
            return None
        return S, "C", [S[i] - (v - i) for i in range(v)], ell
    S = sorted([c for c in range(0, q // 2 + 1)
                if oc.get(c, 0) == (u if c in fijas else 2 * u + 1)], reverse=True)
    for c in range(0, q // 2 + 1):
        if oc.get(c, 0) != ((u + 1) if c in fijas else 2 * u + 2) - (1 if c in S else 0):
            return None
    r = q - v
    if len(S) != r:
        return None
    return S, "D", [S[i] - (r - 1 - i) for i in range(r)], ell


def rho_H(H, r):
    return [r - i for i in range(r)] if H == "C" else [r - 1 - i for i in range(r)]


def resta(a, b, q):
    return [a[i] - b[i] for i in range(q)]


# ------------------------------------------------------------------- H1/H2
print("=" * 96)
print("H1 -- LA IDENTIDAD:  R_ell^{C_m} - R_rho^{C_m}  =  R_s^H - R_{rho_H}^H")
print("H2 -- y el reparto por tipo: v<=a da C, v>a da D")
mal = viv = 0
tipos = defaultdict(int)
ejemplos = []
for m in range(2, 8):
    for q in range(2, 15):
        rho_C = [m - i for i in range(m)]
        for lam in dominantes(m, 5 if m <= 5 else 3):
            d = datos(lam, q)
            if d is None:
                continue
            S, H, eta, ell = d
            if abs(chi_lim(lam, q)) < mpf(10) ** (-12):
                continue
            r = len(eta)
            rH = rho_H(H, r)
            s = [eta[i] + rH[i] for i in range(r)]
            izq = resta(R_de(ell, q, "C"), R_de(rho_C, q, "C"), q)
            der = resta(R_de(s, q, H), R_de(rH, q, H), q) if r else [0] * q
            viv += 1
            a_lim = (q - 1) // 2
            tipos[("v<=a" if m % q <= a_lim else "v>a", H)] += 1
            if izq != der:
                mal += 1
                if len(ejemplos) < 4:
                    ejemplos.append((m, q, tuple(lam), S, H, tuple(eta)))
print("   %d supervivientes con perfil minimo" % viv)
print(("   OK   " if not mal else "  FALLA ") + "H1: la identidad de histogramas se cumple en los %d" % viv)
for e in ejemplos:
    print("        m=%d q=%d lam=%s S=%s H=%s eta=%s" % e)
if mal:
    fallos.append("H1")
print("   reparto (regimen, tipo residual): %s" % dict(tipos))
cruzados = [k for k in tipos if (k[0] == "v<=a") != (k[1] == "C")]
print(("   OK   " if not cruzados else "  FALLA ")
      + "H2: v<=a siempre da tipo C y v>a siempre da tipo D"
      + ("" if not cruzados else "   cruzados: %s" % cruzados))
if cruzados:
    fallos.append("H2")

# ------------------------------------------------------------------- H3
print("=" * 96)
print("H3 -- SENUELOS: tres maneras de romper la identidad")


def senuelo(nombre, trans):
    roto = tot = 0
    for m in range(2, 6):
        for q in range(2, 12):
            rho_C = [m - i for i in range(m)]
            for lam in dominantes(m, 4):
                d = datos(lam, q)
                if d is None:
                    continue
                S, H, eta, ell = d
                r = len(eta)
                if r == 0:
                    continue
                tot += 1
                izq = resta(R_de(ell, q, "C"), R_de(rho_C, q, "C"), q)
                der = trans(H, eta, r, q)
                if der is None:
                    continue
                if izq != der:
                    roto += 1
    print(("   OK   " if roto else "  FALLA ") + "H3: %-40s rompe %d de %d" % (nombre, roto, tot))
    if not roto:
        fallos.append("senuelo " + nombre)


def _tipo_cambiado(H, eta, r, q):
    if r < 2:
        return None
    K = "D" if H == "C" else "C"
    rK = rho_H(K, r)
    return resta(R_de([eta[i] + rK[i] for i in range(r)], q, K), R_de(rK, q, K), q)


def _eta_desplazada(H, eta, r, q):
    rH = rho_H(H, r)
    e = [eta[0] + 1] + list(eta[1:])
    return resta(R_de([e[i] + rH[i] for i in range(r)], q, H), R_de(rH, q, H), q)


def _sin_restar_rho(H, eta, r, q):
    rH = rho_H(H, r)
    return R_de([eta[i] + rH[i] for i in range(r)], q, H)


senuelo("usar el tipo equivocado (C por D)", _tipo_cambiado)
senuelo("desplazar eta en uno", _eta_desplazada)
senuelo("no restar el rho del grupo residual", _sin_restar_rho)

print("=" * 96)
if fallos:
    print("RESULTADO: %d FALLOS -> %s" % (len(fallos), fallos))
    sys.exit(1)
print("RESULTADO: la identidad de histogramas se sostiene, y es la que explica la factorizacion.")
