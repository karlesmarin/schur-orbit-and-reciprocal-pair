# -*- coding: utf-8 -*-
# .POR QUE TRES SENOS?  EL COLAPSO DEL ALTERNANTE, Y LAS DOS RAMAS COMO DOS SISTEMAS DE RAICES.
# 19 de agosto de 2026, noche.  Repaso socratico de las formulas del Paper I.
#
# LO QUE SE AFIRMA, ESCRITO ANTES DE MEDIR.
#
# (P1)  EL DENOMINADOR NO ES UN MISTERIO.  El Vandermonde del alfabeto mu_t U {z,1/z} es
#           prod_{i<j}(x_i-x_j)  propto  (z^t-1)(z^{-t}-1)(z-1/z)  =  -(z^{t/2}-z^{-t/2})^2 (z-1/z),
#       o sea  sinh^2(t.theta/2).sinh(theta),  que es EXACTAMENTE el denominador del Teorema.  Luego
#       todo el contenido esta en que el alternante de ARRIBA factorice en tres senos.
#
# (P2)  Y factoriza porque solo sobreviven CUATRO terminos.  Desarrollando por Laplace a lo largo de
#       las dos filas del par, el menor t x t de las raices de la unidad se anula salvo que los
#       exponentes restantes toquen TODAS las clases.  Con t+2 exponentes en t clases el exceso es 2:
#       hay que quitar las dos cuentas sobrantes, una de A y una de B, y eso deja 4 elecciones.
#
# (P3)  LA IDENTIDAD.  Con u = a1-b1, v = a1-b2, w = a2-b1, x = a2-b2:
#           sinh(u.th) - sinh(v.th) - sinh(w.th) + sinh(x.th)
#             = -4 sinh(d1.th/2) sinh(d2.th/2) sinh(d3t.th/2),
#       con d1 = a1-a2, d2 = b1-b2, d3t = a1+a2-b1-b2.  Dos aplicaciones de suma-a-producto.
#       Los tres d NO son tres cosas sueltas: son las dos brechas DENTRO de cada clase y la brecha
#       ENTRE clases, que es lo unico que la identidad deja en pie.
#
# (P4)  Y LA OBSERVACION QUE VALE POR SI SOLA.  En la rama (ii) --- una sola clase con n_i = 3, con
#       A = {p,q} y B = {q,r} SOLAPADAS --- se tiene a2 = b1 = q, luego
#           d1 = p-q,  d2 = q-r,  d3t = p+q-q-r = p-r = d1 + d2.
#       Las tres formas dejan de ser independientes.  Las dos ramas del Teorema no son dos casos:
#       son DOS SISTEMAS DE RAICES.  Dos clases -> tres formas independientes.  Una clase -> alpha,
#       beta, alpha+beta: las tres raices positivas de A_2.
#
# QUE SE MIDE
#   C1  (P1) el Vandermonde del alfabeto contra la forma cerrada, como polinomio de Laurent.
#   C2  (P3) la identidad de los cuatro senos, sobre enteros al azar reproducibles.
#   C3  (P2)+(P3) el alternante entero contra -4 sinh sinh sinh, dividido por el Vandermonde, contra
#       el valor DIRECTO s_lambda(alfabeto).  Es el Teorema, reconstruido desde el colapso.
#   C4  (P4) en cuantas lambda la rama es la solapada, y en TODAS ellas .vale d3 = d1 + d2?
#   D1  SENUELO: la identidad con los signos +,+,+,+ en vez de +,-,-,+.  Tiene que fallar.
#   D2  SENUELO: d3 = d1 - d2 en la rama solapada.  Tiene que fallar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python pI_tres_senos_porque.py > pI_tres_senos_porque_OUT.txt 2>&1

import json
import sys
from collections import Counter

CASOS = [2, 3, 4, 5, 6]
LMAX = 10


def padd(a, b):
    o = dict(a)
    for k in b:
        o[k] = o.get(k, 0) + b[k]
        if o[k] == 0:
            del o[k]
    return o


def pneg(a):
    return {k: -v for k, v in a.items()}


def pmul(a, b):
    o = {}
    for k1 in a:
        for k2 in b:
            k = k1 + k2
            o[k] = o.get(k, 0) + a[k1] * b[k2]
    return {k: v for k, v in o.items() if v != 0}


def S(d):
    """2*sinh(d.theta/2) en la variable w = z^{1/2}:  w^d - w^{-d}.
       OJO con d = 0: el literal {0: 1, -0: -1} colapsa a {0: -1} porque -0 == 0, y daba
       2sinh(0) = -1 en vez de 0.  Es el caso que aparece cuando dos exponentes coinciden y
       cuando d3 = 0 --- que es justo la rama en que el Teorema dice Phi = 0."""
    if d == 0:
        return {}
    return {d: 1, -d: -1}


def det(M):
    n = len(M)
    if n == 0:
        return {0: 1}
    if n == 1:
        return dict(M[0][0])
    acc = {}
    for j in range(n):
        if not M[0][j]:
            continue
        men = [[M[i][k] for k in range(n) if k != j] for i in range(1, n)]
        tm = pmul(M[0][j], det(men))
        acc = padd(acc, tm if j % 2 == 0 else pneg(tm))
    return acc


print("=" * 100)
print(".POR QUE TRES SENOS?  El colapso del alternante, y las dos ramas como dos sistemas de raices")
print("=" * 100)
print("")

# ---- C2: la identidad de los cuatro senos, en w = z^{1/2} (exponentes DOBLADOS para no partir) ---
# 2sinh(d.th/2) -> w^d - w^{-d}.  El producto de tres va con un factor 4 que se lleva la cuenta.
c2_ok, c2_bad, d1_ok = 0, [], 0
paso = 0
for a1 in range(6, 0, -1):
    for a2 in range(a1 - 1, -1, -1):
        for b1 in range(6, 0, -1):
            for b2 in range(b1 - 1, -1, -1):
                paso += 1
                u, v, w_, x = a1 - b1, a1 - b2, a2 - b1, a2 - b2
                izq = padd(padd(S(2 * u), pneg(S(2 * v))), padd(pneg(S(2 * w_)), S(2 * x)))
                dd1, dd2, dd3 = a1 - a2, b1 - b2, a1 + a2 - b1 - b2
                der = pmul(pmul(S(dd1), S(dd2)), S(dd3))
                # -4 sinh sinh sinh, con 2sinh en cada factor: -(1/2)*prod(2sinh) ... se compara
                # 2*izq contra -prod, que es la forma sin fracciones.
                # Todo en 2sinh: LHS = (1/2)(U-V-W+X), RHS = -4 (A/2)(B/2)(C/2) = -ABC/2.
                # Luego la identidad sin fracciones es  U - V - W + X = -ABC.
                if padd(izq, der) == {}:
                    c2_ok += 1
                else:
                    c2_bad.append((a1, a2, b1, b2))
                mal = padd(padd(S(2 * u), S(2 * v)), padd(S(2 * w_), S(2 * x)))
                d1_ok += (padd(mal, der) == {})
print("  C2  FATAL la identidad de los cuatro senos : %d de %d   %s"
      % (c2_ok, paso, "PASA" if not c2_bad else "FALLA en %s" % c2_bad[:5]))
print("      2*(sinh u - sinh v - sinh w + sinh x) = -(2sinh d1)(2sinh d2)(2sinh d3), en w = z^(1/2)")
print("  D1  senuelo con los cuatro signos + : %d de %d   %s"
      % (d1_ok, paso, "NO DISCRIMINA" if d1_ok == paso else "falla en %d (bien)" % (paso - d1_ok)))
print("")

resumen = {"C2": [int(c2_ok), int(paso)], "D1": int(d1_ok)}

# ---- C1 y C4 --------------------------------------------------------------------------------
for T in CASOS:
    N = T + 2
    # C1: el Vandermonde del alfabeto.  Se compara en w = z^{1/2}: el alfabeto es
    # {w^{2k*?}}... se hace en z entero: producto explicito sobre las letras.
    # letras: las t raices (simbolicas, se usa la identidad prod_omega (z - omega) = z^t - 1)
    # (z^t - 1)(z^{-t} - 1)(z - z^{-1})  contra  -(z^{t/2}-z^{-t/2})^2 (z - z^{-1})
    # TODO en la variable w = z^{1/2}: un z^k es w^{2k}.  Mezclarlas fue el fallo de la primera
    # version, y daba FALLA en las cinco t sin que el algebra tuviera nada malo.
    zt = {2 * T: 1, 0: -1}                 # z^t - 1
    ztm = {-2 * T: 1, 0: -1}               # z^{-t} - 1
    zz = {2: 1, -2: -1}                    # z - 1/z
    izq = pmul(pmul(zt, ztm), zz)
    der = pneg(pmul(pmul(S(T), S(T)), zz))
    ok1 = (izq == der)
    # C4: el reparto de ramas y la dependencia d3 = d1 + d2 en la solapada
    ramas = Counter()
    c4_ok, c4_bad, d2_ok = 0, [], 0
    # enumeracion de beta-sets estrictamente decrecientes con maximo LMAX
    import itertools
    for comb in itertools.combinations(range(LMAX + 1), N):
        beta = tuple(sorted(comb, reverse=True))
        clases = {}
        for x in beta:
            clases.setdefault(x % T, []).append(x)
        if any(len(v) == 0 for v in clases.values()) or len(clases) < T:
            ramas["clase vacia -> cero"] += 1
            continue
        oc = sorted(len(v) for v in clases.values())
        if oc == [1] * (T - 2) + [2, 2]:
            ramas["dos clases con 2"] += 1
        elif oc == [1] * (T - 1) + [3]:
            ramas["una clase con 3"] += 1
            cl = [v for v in clases.values() if len(v) == 3][0]
            p, q, r = sorted(cl, reverse=True)
            dd1, dd2, dd3 = p - q, q - r, abs(p + q - q - r)
            if dd3 == dd1 + dd2:
                c4_ok += 1
            else:
                c4_bad.append((p, q, r))
            d2_ok += (dd3 == abs(dd1 - dd2))
        else:
            ramas["otra"] += 1
    n3 = ramas["una clase con 3"]
    print("  t=%d   C1  Vandermonde del alfabeto = -(2sinh(t))^2 . (z-1/z) : %s" % (T, "OK" if ok1 else "FALLA"))
    print("        ramas : %s" % dict(ramas))
    print("        C4  en la rama SOLAPADA, .d3 = d1 + d2?  %d de %d   %s"
          % (c4_ok, n3, "SIEMPRE" if (n3 > 0 and not c4_bad) else ("sin casos" if n3 == 0 else "FALLA")))
    print("        D2  senuelo d3 = |d1 - d2| : %d de %d   %s"
          % (d2_ok, n3, "NO DISCRIMINA" if (n3 > 0 and d2_ok == n3) else "falla (bien)"))
    print("")
    sys.stdout.flush()
    resumen["t=%d" % T] = dict(C1=bool(ok1), ramas={k: int(v) for k, v in ramas.items()},
                               C4=[int(c4_ok), int(n3)], D2=int(d2_ok))

json.dump(resumen, open("pI_tres_senos_porque_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
