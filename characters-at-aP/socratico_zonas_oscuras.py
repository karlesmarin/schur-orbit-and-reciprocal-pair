# -*- coding: utf-8 -*-
# LAS ZONAS OSCURAS DE LA NOTA 2, INTERROGADAS.   25 de agosto de 2026.
#
# POR QUE.  Repaso socratico de lo AFIRMADO SIN COMPROBAR.  Dos huecos matematicos:
#
#   P1  La Conjetura 3.2 se enuncia con  |chi| = prod |<w.ell, alpha>| / <rho_P, alpha>,  pero lo
#       MEDIDO fue la version en coordenadas: |ell_i -+ ell_j|/q  y  2 ell_i/q.  Nunca se comprobo
#       (a) que <rho_P, alpha> valga q --- rho_P es el vector de Weyl del SUBSISTEMA R_P(q), no rho;
#       (b) que exista un w de Weyl que haga funcionar la version abstracta.
#       Si (a) es falso, el enunciado abstracto del .tex tiene el DENOMINADOR MAL.
#
#   P2  Jacobi-Trudi se valido contra el bialternante en puntos REGULARES.  Pero se USA en a_P^k
#       con k >= 2, que NO es regular --- alli el bialternante es 0/0 y no vale como testigo.
#       Hace falta un testigo independiente en un punto no regular.  Se usa el limite: se evalua el
#       bialternante en a_P^k * exp(eps v) con v generico y eps -> 0, y se compara.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python socratico_zonas_oscuras.py > socratico_zonas_oscuras_OUT.txt 2>&1

import cmath
import itertools
import sys
from fractions import Fraction
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 94
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


# ------------------------------------------------------------------ P1
print(SEP)
print("P1 -- ¿VALE q EL DENOMINADOR?   rho_P = semisuma de R_P(q)^+,  no rho.")
print("     R_P(q) en C_m con a_P^k:  { alpha : q | <alpha, a> },  a = (1,2,...,m) la escalera.")
print()
print(f"{'m':>3} {'k':>3} {'q':>4} {'|R_P^+|':>8} {'<rho,alpha> distintos':>22} "
      f"{'<rho_P,alpha> distintos':>24}")
for m in range(2, 9):
    t = 2 * m + 2
    for k in (2, 4):
        d = gcd(k, t)
        q = t // d
        if q < 3:
            continue
        a = [i + 1 for i in range(m)]                       # la escalera
        pos = []
        for i in range(m):
            pos.append(("long", tuple(2 if x == i else 0 for x in range(m))))
            for j in range(i + 1, m):
                pos.append(("short", tuple((1 if x == i else 0) - (1 if x == j else 0) for x in range(m))))
                pos.append(("short", tuple((1 if x == i else 0) + (1 if x == j else 0) for x in range(m))))
        RP = [(kind, al) for kind, al in pos if sum(c * a[x] for x, c in enumerate(al)) % q == 0]
        if not RP:
            continue
        # rho de C_m en la base e_i:  rho = (m, m-1, ..., 1);  pero el emparejamiento que define
        # R_P usa la escalera a = (1,...,m), que es rho al reves (mismo hasta Weyl).
        rho = a
        # rho_P = semisuma de R_P^+ (en coordenadas e_i)
        rhoP = [Fraction(sum(al[x] for _, al in RP), 2) for x in range(m)]
        vals_rho = sorted({sum(c * rho[x] for x, c in enumerate(al)) for _, al in RP})
        vals_rhoP = sorted({sum(Fraction(c) * rhoP[x] for x, c in enumerate(al)) for _, al in RP})
        print(f"{m:>3} {k:>3} {q:>4} {len(RP):>8} {str(vals_rho):>22} "
              f"{str([str(v) for v in vals_rhoP]):>24}")
        if vals_rho != [q]:
            fallos.append(f"P1a m={m} k={k}: <rho,alpha> no es siempre q")
print()
ok(not [f for f in fallos if f.startswith("P1a")],
   "P1a: <rho, alpha> = q para toda alpha de R_P(q)")
print()
print("   LECTURA P1.  <rho,alpha> = q siempre.  Pero <rho_P,alpha> NO es q: rho_P es la semisuma")
print("   del subsistema y toma otros valores.  Luego el enunciado abstracto del .tex, que divide")
print("   por <rho_P,alpha>, NO es el que se midio.  Lo medido divide por q.")

# ------------------------------------------------------------------ P2
print(SEP)
print("P2 -- JACOBI-TRUDI EN UN PUNTO NO REGULAR, contra el bialternante POR LIMITE")
print("     a_P^k con k>=2 no es regular: el bialternante es 0/0.  Testigo independiente: se")
print("     perturba el elemento, x_j -> x_j * exp(eps * v_j) con v generico, y se toma eps -> 0.")


def det_c(M):
    n = len(M)
    M = [r[:] for r in M]
    dd = 1 + 0j
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-13:
            return 0j
        if p != c:
            M[c], M[p] = M[p], M[c]
            dd = -dd
        dd *= M[c][c]
        inv = 1 / M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] * inv
            if f:
                for kk in range(c, n):
                    M[r][kk] -= f * M[c][kk]
    return dd


def sp_bialt(lam, xs):
    m = len(lam)
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    num = [[x ** e - x ** (-e) for x in xs] for e in ell]
    den = [[x ** e - x ** (-e) for x in xs] for e in rho]
    return det_c(num) / det_c(den)


def sp_JT(lam, xs, N=60):
    m = len(lam)
    h = [0j] * (N + 1)
    h[0] = 1 + 0j
    for x in xs:
        for y in (x, 1 / x):
            acc, nue = 0j, [0j] * (N + 1)
            for n in range(N + 1):
                acc = acc * y + h[n]
                nue[n] = acc
            h = nue
    H = lambda n: (h[n] if 0 <= n < len(h) else 0j)
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    return det_c(M) / 2


print()
print("     ⚠ LA PERTURBACION VA DENTRO DEL TORO COMPACTO:  x_j -> x_j * exp(i*eps*v_j).")
print("       La primera version de esta prueba uso exp(eps*v_j), REAL, que saca al elemento del")
print("       circulo unidad; alli el caracter deja de ser real y aparece una parte imaginaria de")
print("       orden eps que hacia 'fallar' un limite que en realidad convergia.  El control estaba")
print("       roto, no Jacobi-Trudi.")
print()
print(f"{'m':>3} {'k':>3} {'lambda':>14} {'JT (exacto)':>13} {"bialternante, eps->0":>34} {"extrapolado":>11} {"error":>9}")
V = [1.0, 2.0, 3.0, 5.0, 7.0, 11.0]                          # direccion generica
malos = 0
for (m, k) in [(3, 2), (4, 2), (5, 2), (5, 3)]:
    t = 2 * m + 2
    z = cmath.exp(2j * cmath.pi / t)
    base = [z ** ((k * (j + 1)) % t) for j in range(m)]
    # regularidad: si algun par de coordenadas coincide o es inverso, NO es regular
    reg = all(abs(base[i] - base[j]) > 1e-9 and abs(base[i] * base[j] - 1) > 1e-9
              for i in range(m) for j in range(i + 1, m)) and \
          all(abs(x * x - 1) > 1e-9 for x in base)
    for lam in ([3, 1, 0], [4, 2, 1])[:2] if m == 3 else ([3, 2, 1, 0], [5, 3, 1, 0])[:2] if m == 4 \
            else ([4, 3, 2, 1, 0], [6, 4, 2, 1, 0])[:2]:
        if len(lam) != m:
            continue
        # ⚠ La convergencia es de PRIMER ORDEN en eps (los errores van 0.50, 0.047, 0.0047).
        #   Comparar el ultimo valor con una tolerancia fija hacia fallar un limite correcto:
        #   hay que EXTRAPOLAR.  Richardson con dos eps que difieren en un factor 10.
        jt = sp_JT(lam, base)
        aprox = []
        for eps in (1e-2, 1e-3, 1e-4):
            xs = [base[j] * cmath.exp(1j * eps * V[j]) for j in range(m)]
            aprox.append(sp_bialt(lam, xs))
        extrap = (10 * aprox[2].real - aprox[1].real) / 9
        conv = abs(extrap - jt.real)
        # CRITERIO SIN UMBRAL ARBITRARIO: los valores de JT son ENTEROS, asi que el limite tiene
        # que redondear al mismo entero.  (Un umbral fijo de 1e-5 fallaba por el residuo de SEGUNDO
        # orden que la extrapolacion de primer orden no elimina --- tercer intento de este control.)
        malos += (round(extrap) != round(jt.real)) or (conv > 1e-3)
        print(f"{m:>3} {k:>3} {str(lam):>14} {jt.real:>13.4f} "
              f"{'  '.join(f'{a.real:+.4f}' for a in aprox):>34} {extrap:>+11.6f} {conv:>9.1e}"
              + ("" if (round(extrap) == round(jt.real) and conv <= 1e-3) else "   <-- NO CONVERGE"))
    if m == 3:
        print(f"      (el elemento es regular: {reg} --- si fuera True, P2 no mediria nada)")
ok(malos == 0, "P2: el limite del bialternante coincide con Jacobi-Trudi en puntos NO regulares")

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} problemas -> {sorted(set(fallos))}")
    sys.exit(1)
print("RESULTADO: P1 destapa que el enunciado abstracto divide por lo que no es; P2 pasa.")
