# -*- coding: utf-8 -*-
# LA FACTORIZACION RESIDUAL, MEDIDA CONTRA NUESTRAS FICHAS.   9 de septiembre de 2026.
#
# QUE PROPONE.  Que el corrector U que en §14.1 aparece como una constante medida en {1,2,3,4,6}
# es en realidad EL CARACTER DE UN GRUPO RESIDUAL PEQUENO Y REGULAR:
#
#     |chi_lambda^{C_m}(g_{p/q})|  =  kappa_q(lambda) . D(lambda,q) . chi_eta^H(h^H_{p/q})
#
# con H = Sp(2r) cuando v <= a  y  H = SO(2r) cuando v > a, y eta el perfil S menos su escalera.
#
# POR QUE HAY QUE MEDIRLO Y NO CREERLO.  Una formula que encaja demasiado bien es justo la que
# hay que medir antes de escribirla.  Aqui se mide contra
# el instrumento del propio articulo -- D_de() de teorema_II.py y chi_lim() de galois.py -- y con
# senuelos, que es lo que faltaba la vez anterior.
#
# EL INSTRUMENTO DEL RESIDUAL.  chi_lim evalua en la RECTA DE rho por la forma de senos,
#     chi = prod_{alpha>0} sin(pi.theta.<l,alpha>) / sin(pi.theta.<rho,alpha>),
# y h^H_{p/q} esta justo en la recta de rho_H.  Asi que basta cambiarle el sistema de raices y el
# rho: la misma formula sirve para Sp(2r) y para SO(2r).
#
# QUE COMPRUEBA
#   V0  el evaluador residual, contra chi_lim en tipo C, donde los dos deben coincidir
#   V1  la identidad (R_1) sobre TODOS los supervivientes, q cualquiera
#   V2  ...y en particular sobre las fichas con la condicion (I) -- las 828 -- donde
#       gamma = kappa . chi_eta^H  tiene que dar exactamente {1,2,3,4,6}
#   V3  su tabla de SO(4) para q=6, v=4, las seis filas
#   V4  el caso trabajado Sp(8), q=6, lambda=(2,1,0,0) -> -6
#   V5  SENUELOS: cuatro maneras de romperla que TIENEN que saltar
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python -u residual_factorizacion.py > residual_factorizacion_OUT.txt 2>&1

import os
import sys
from collections import defaultdict
from fractions import Fraction
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from mpmath import mp, mpf, sin as msin, pi as mpi
mp.dps = 40

from teorema_II import D_de, perfil, dim_sp, dim_so_odd, dim_gl
from galois import chi_lim, raices_C, dominantes, condicion_I

fallos = []


def raices_D(r):
    """raices positivas de D_r: e_i +- e_j, i<j.  (D_1 no tiene ninguna.)"""
    R = []
    for i in range(r):
        for j in range(i + 1, r):
            R.append(tuple(1 if x == i else (-1 if x == j else 0) for x in range(r)))
            R.append(tuple(1 if x in (i, j) else 0 for x in range(r)))
    return R


def chi_linea(lam, R, rho, q, p=1):
    """chi_lam en el elemento de la recta de rho con denominador q, por la forma de senos.
    Es chi_lim generalizado: mismo metodo, otro sistema de raices y otro rho."""
    n = len(lam)
    if n == 0 or not R:
        return mpf(1)
    ell = [lam[i] + rho[i] for i in range(n)]
    th = mpf(p) / q + mpf(10) ** (-25)
    v = mpf(1)
    for a in R:
        vl = sum(ell[i] * a[i] for i in range(n))
        vr = sum(rho[i] * a[i] for i in range(n))
        v *= msin(mpi * th * vl) / msin(mpi * th * vr)
    return v


# ------------------------------------------------------------------ V0
print("=" * 96)
print("V0 -- EL EVALUADOR RESIDUAL, CONTRA chi_lim EN TIPO C")
mal = viv = 0
for r in range(1, 5):
    for q in range(2, 13):
        for lam in dominantes(r, 3):
            a = chi_linea(lam, raices_C(r), [r - i for i in range(r)], q)
            b = chi_lim(lam, q)
            viv += 1
            if abs(a - b) > mpf(10) ** (-18):
                mal += 1
print("   %d evaluaciones" % viv)
print(("   OK   " if not mal else "  FALLA ") + "V0: la forma generalizada reproduce chi_lim en tipo C")
if mal:
    fallos.append("V0")


# ---------------------------------------------------- el perfil y el residual
def datos_residuales(lam, q):
    """(S, H, eta, kappa) segun su receta.  None si el perfil no es minimo."""
    m = len(lam)
    u, v = divmod(m, q)
    a = (q - 1) // 2
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    P = perfil(ell, q)
    fijas = {0} | ({q // 2} if q % 2 == 0 else set())
    ocup = {c: len(P.get(c, [])) for c in range(0, q // 2 + 1)}
    if v <= a:
        S = sorted([c for c in range(1, a + 1) if ocup.get(c, 0) == 2 * u + 1], reverse=True)
        # el perfil tiene que ser exactamente el minimo
        for c in range(0, q // 2 + 1):
            esperado = (u if c in fijas else 2 * u + (1 if c in S else 0))
            if ocup.get(c, 0) != esperado:
                return None
        if len(S) != v:
            return None
        r = v
        eta = [S[i] - (r - i) for i in range(r)]
        return S, "C", eta, 1
    else:
        S = sorted([c for c in range(0, q // 2 + 1)
                    if ocup.get(c, 0) == (u if c in fijas else 2 * u + 1)], reverse=True)
        for c in range(0, q // 2 + 1):
            esperado = ((u + 1) if c in fijas else 2 * u + 2) - (1 if c in S else 0)
            if ocup.get(c, 0) != esperado:
                return None
        if len(S) != q - v:
            return None
        r = q - v
        eta = [S[i] - (r - 1 - i) for i in range(r)]
        kappa = 1 if 0 in S else 2
        return S, "D", eta, kappa


def residual(H, eta, q, p=1):
    r = len(eta)
    if r == 0:
        return mpf(1)
    if H == "C":
        return chi_linea(eta, raices_C(r), [r - i for i in range(r)], q, p)
    return chi_linea(eta, raices_D(r), [r - 1 - i for i in range(r)], q, p)


# ------------------------------------------------------------------ V1/V2
print("=" * 96)
print("V1 -- (R_1) SOBRE TODOS LOS SUPERVIVIENTES,  q CUALQUIERA")
print("V2 -- ...y gamma = kappa . chi_eta^H sobre las fichas con la condicion (I)")
TOL = mpf(10) ** (-12)
viv = mal = sinperfil = 0
gammas = defaultdict(int)
gammas_I = defaultdict(int)
vivo_I = 0
ejemplos = []
for m in range(2, 7):
    for q in range(2, 13):
        rho = [m - i for i in range(m)]
        for lam in dominantes(m, 5 if m <= 4 else 3):
            c = chi_lim(lam, q)
            if abs(c) < TOL:
                continue
            ell = [lam[i] + rho[i] for i in range(m)]
            D = D_de(ell, q)
            if D is None:
                continue
            dat = datos_residuales(lam, q)
            if dat is None:
                sinperfil += 1
                continue
            S, H, eta, kap = dat
            res = residual(H, eta, q)
            pred = kap * D * res
            viv += 1
            if abs(abs(c) - abs(pred)) > mpf(10) ** (-10) * max(mpf(1), abs(c)):
                mal += 1
                if len(ejemplos) < 5:
                    ejemplos.append((m, q, tuple(lam), S, H, tuple(eta), kap, D,
                                     float(res), float(c), float(pred)))
            g = kap * res
            gammas[round(float(g), 6)] += 1
            if condicion_I(m, q):
                vivo_I += 1
                gammas_I[round(float(g), 6)] += 1

print("   supervivientes con perfil minimo: %d   (descartados por perfil no minimo: %d)"
      % (viv, sinperfil))
print(("   OK   " if not mal else "  FALLA ") + "V1: |chi| = kappa . D . chi_eta^H en los %d" % viv
      + ("" if not mal else "   %d fallos" % mal))
for e in ejemplos:
    print("        m=%d q=%d lam=%s S=%s H=%s eta=%s kap=%d D=%d res=%.6f  chi=%.6f pred=%.6f" % e)
if mal:
    fallos.append("V1")

print("   reparto de gamma = kappa . chi_eta^H, en TODOS: %s"
      % dict(sorted(gammas.items())[:12]))
print("   y en las fichas con la condicion (I) (%d): %s" % (vivo_I, dict(sorted(gammas_I.items()))))
esperados = {1.0, 2.0, 3.0, 4.0, 6.0}
fuera = [g for g in gammas_I if abs(g - round(g)) > 1e-6 or round(g) not in {1, 2, 3, 4, 6}]
print(("   OK   " if not fuera else "  FALLA ")
      + "V2: bajo (I) los gamma son exactamente {1,2,3,4,6}"
      + ("" if not fuera else "   fuera: %s" % fuera))
if fuera:
    fallos.append("V2")

# ------------------------------------------------------------------ V3
print("=" * 96)
print("V3 -- SU TABLA DE SO(4) PARA q=6, v=4:  las seis filas")
TABLA = {(0, 1): (1, 1, 1), (0, 2): (3, 1, 3), (0, 3): (4, 1, 4),
         (1, 2): (2, 2, 4), (1, 3): (3, 2, 6), (2, 3): (1, 2, 2)}
mal = 0
print("      S        eta        chi^{SO(4)}  kappa   gamma   su tabla")
for Sset, (res_esp, kap_esp, gam_esp) in sorted(TABLA.items()):
    S = sorted(Sset, reverse=True)
    r = 2
    eta = [S[i] - (r - 1 - i) for i in range(r)]
    res = residual("D", eta, 6)
    kap = 1 if 0 in S else 2
    gam = kap * res
    ok = abs(res - res_esp) < 1e-9 and kap == kap_esp and abs(gam - gam_esp) < 1e-9
    print("      %-8s %-10s %-12.6f %-7d %-7.3f %s"
          % (str(Sset), str(tuple(eta)), float(res), kap, float(gam),
             "SI" if ok else "NO  (su fila: %s)" % (TABLA[Sset],)))
    if not ok:
        mal += 1
print(("   OK   " if not mal else "  FALLA ") + "V3: las seis filas de su tabla")
if mal:
    fallos.append("V3")

# ------------------------------------------------------------------ V4
print("=" * 96)
print("V4 -- EL CASO TRABAJADO:  Sp(8), q=6, lambda=(2,1,0,0)")
lam = (2, 1, 0, 0)
ell = [lam[i] + 4 - i for i in range(4)]
dat = datos_residuales(lam, 6)
c = chi_lim(lam, 6)
if dat is None:
    print("  FALLA V4: el perfil no sale minimo")
    fallos.append("V4")
else:
    S, H, eta, kap = dat
    D = D_de(ell, 6)
    res = residual(H, eta, 6)
    print("   ell=%s  S=%s  H=%s  eta=%s  kappa=%d  D=%s  chi_eta=%.6f"
          % (ell, S, H, tuple(eta), kap, D, float(res)))
    print("   chi medido = %.6f     kappa.D.chi_eta = %.6f" % (float(c), float(kap * D * res)))
    ok = abs(abs(c) - kap * D * res) < 1e-9 and abs(c + 6) < 1e-9
    print(("   OK   " if ok else "  FALLA ") + "V4: da -6 y la factorizacion cuadra")
    if not ok:
        fallos.append("V4")

# ------------------------------------------------------------------ V5
print("=" * 96)
print("V5 -- SENUELOS: cuatro maneras de romperla, que TIENEN que saltar")


def rompe(nombre, transforma):
    roto = 0
    tot = 0
    for m in range(2, 6):
        for q in range(2, 11):
            rho = [m - i for i in range(m)]
            for lam in dominantes(m, 4):
                c = chi_lim(lam, q)
                if abs(c) < TOL:
                    continue
                ell = [lam[i] + rho[i] for i in range(m)]
                D = D_de(ell, q)
                dat = datos_residuales(lam, q)
                if D is None or dat is None:
                    continue
                tot += 1
                pred = transforma(dat, D, q)
                if pred is None:
                    continue
                if abs(abs(c) - abs(pred)) > mpf(10) ** (-10) * max(mpf(1), abs(c)):
                    roto += 1
    print(("   OK   " if roto else "  FALLA ") + "V5: %-42s rompe %d de %d" % (nombre, roto, tot))
    if not roto:
        fallos.append("senuelo " + nombre)


def eta_desplazada(d, D, q):
    if not d[2]:                       # eta vacio: no hay nada que desplazar, el caso no cuenta
        return None
    return d[3] * D * residual(d[1], [d[2][0] + 1] + list(d[2][1:]), q)


def tipo_cambiado(d, D, q):
    if len(d[2]) < 2:                  # con rango 0 o 1 los dos tipos coinciden; no discrimina
        return None
    return d[3] * D * residual("D" if d[1] == "C" else "C", d[2], q)


rompe("desplazar eta en uno", eta_desplazada)
rompe("usar el tipo equivocado (C por D y al reves)", tipo_cambiado)
rompe("kappa siempre 1",
      lambda d, D, q: 1 * D * residual(d[1], d[2], q))
rompe("residual siempre 1 --- que es el enunciado VIEJO",
      lambda d, D, q: d[3] * D * mpf(1))

print("=" * 96)
if fallos:
    print("RESULTADO: %d FALLOS -> %s" % (len(fallos), fallos))
    sys.exit(1)
print("RESULTADO: la factorizacion residual se sostiene en todo lo que aqui se ha medido.")
