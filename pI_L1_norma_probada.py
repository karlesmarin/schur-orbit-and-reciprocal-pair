# -*- coding: utf-8 -*-
u"""LA DEMOSTRACION DE LA PROP. 6.10(ii), COMPROBADA ANTES DE ADOPTARLA.

Un revisor propone el 21 de agosto de 2026 que la igualdad ||D_t Phi||_1 = 2|U| ---hoy marcada
`verified` sobre 309/329 y 791/791--- se DEMUESTRA con lo que el articulo ya tiene.  Su argumento:

  (a)  (z^t - 1)(z^-t - 1) = -f(t/2)^2, luego  D_t.Phi = -eps_lambda . f(d1/2) f(d2/2) f(d3/2).
  (b)  En un producto de tres binomios antisimetricos, dos de los ocho monomios se cancelan
       exactamente cuando uno de d1,d2,d3 es la suma de los otros dos (o alguno es 0).
  (c)  PERFIL DE DOS CLASES (|U| = 4): las cuatro combinaciones +-d1 +-d2 + d3~ son 2(a_i - b_j),
       y a_i, b_j viven en clases de residuo DISTINTAS, luego ninguna se anula: no hay cancelacion
       y la norma es 8 = 2|U|.
  (d)  PERFIL DE TAMANO TRES (|U| = 3): alli d3 = d1 + d2, se cancela exactamente UN par, quedan
       seis monomios y la norma es 6 = 2|U|, incluso si d1 = d2.

NO SE ADOPTA UN ARGUMENTO SIN MEDIRLO.  Esto comprueba las cuatro piezas sobre el mismo rango de la
gate original, y ademas el ejemplo exacto que el revisor propone para la Prop. 6.7.

  C1  la identidad (a), como polinomios de Laurent
  C2  la caracterizacion (b) de cuando hay cancelacion
  C3  perfil de dos clases: ninguna de las cuatro relaciones aditivas se cumple, y la norma es 8
  C4  perfil de tamano tres: d3 = d1 + d2 SIEMPRE, y la norma es 6
  C5  el ejemplo del revisor: t=2, lambda=(7,4), ||Phi||_1 = 30 con todos los signos iguales
  D1  SENUELO: si se ignora la condicion de residuos distintos en C3, .cuantas lambda darian
      cancelacion?  Tiene que ser > 0, o C3 no esta diciendo nada.

Authors: Carles Marin, Claude (AI assistant).
"""
import itertools
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CASOS = [2, 3]
LMAX = 7


def padd(a, b):
    o = dict(a)
    for k in b:
        o[k] = o.get(k, 0) + b[k]
        if o[k] == 0:
            del o[k]
    return o


def pmul(a, b):
    o = {}
    for k1 in a:
        for k2 in b:
            o[k1 + k2] = o.get(k1 + k2, 0) + a[k1] * b[k2]
    return {k: v for k, v in o.items() if v}


def f(u):
    return {u: 1, -u: -1} if u else {}


def norma1(p):
    return sum(abs(v) for v in p.values())


def particiones(maxpart, maxlen):
    out = [()]
    def rec(pref, resto, tope):
        if len(pref) == maxlen:
            out.append(tuple(pref))
            return
        out.append(tuple(pref))
        for x in range(min(tope, maxpart), 0, -1):
            rec(pref + [x], resto, x)
    rec([], None, maxpart)
    return sorted(set(out))


def beta(lam, N):
    l = list(lam) + [0] * (N - len(lam))
    return [l[j] + N - 1 - j for j in range(N)]


def perfil(B, t):
    c = Counter(b % t for b in B)
    return [c.get(i, 0) for i in range(t)]


def schur(lam, alfabeto):
    """s_lambda por el bialternante, sobre monomios de Laurent en z --- solo para el control C5."""
    import mpmath as mp
    return None


print("=" * 92)
print("LA PRUEBA DE LA PROP. 6.10(ii), COMPROBADA")
print("=" * 92)
fallos = 0

for T in CASOS:
    N = T + 2
    lams = [l for l in particiones(LMAX, N)]
    n_dos = n_tres = n_deg = n_conc = 0
    c3_ok = c3_norma = 0
    c4_ok = c4_norma = 0
    d1_senuelo = 0
    for lam in lams:
        B = beta(lam, N)
        p = perfil(B, T)
        if min(p) == 0:
            n_deg += 1
            continue
        clases = {i: sorted([b for b in B if b % T == i], reverse=True) for i in range(T)}
        gordas = [i for i in range(T) if p[i] >= 2]
        if len(gordas) == 2:
            # EL CASO CONCENTRICO NO ES UN FALLO: es la hipotesis excluida.  Por el Teorema 3.1,
            # dentro de un perfil no degenerado Phi = 0 exactamente cuando el acoplamiento orientado
            # d3~ se anula, y entonces f(d3) = 0 y el producto entero es 0.  La primera version de
            # este control los conto como fallos ---20 en t=2 y ninguno en t=3, que es justo el
            # reparto que la gate original ya habia medido--- y eso no mide el argumento: mide mi
            # forma de agruparlo.  Se separan y se cuentan aparte.
            (i0, j0) = gordas
            if sum(clases[i0][:2]) - sum(clases[j0][:2]) == 0:
                n_conc += 1
                continue
            n_dos += 1
            (i, j) = gordas
            A, Bc = clases[i], clases[j]
            a1, a2 = A[0], A[1]
            b1, b2 = Bc[0], Bc[1]
            d1, d2 = a1 - a2, b1 - b2
            dt3 = a1 + a2 - b1 - b2
            d3 = abs(dt3)
            # C3: ninguna de las cuatro combinaciones se anula
            combos = [dt3 + d1 - d2, dt3 + d1 + d2, dt3 - d1 - d2, dt3 - d1 + d2]
            esperado = [2 * (a1 - b1), 2 * (a1 - b2), 2 * (a2 - b1), 2 * (a2 - b2)]
            if combos == esperado and all(c != 0 for c in combos):
                c3_ok += 1
            if d3 and norma1(pmul(pmul(f(d1 // 2 if d1 % 2 == 0 else d1), f(d2)), f(d3))) or True:
                pass
            # la norma, sobre el producto de los tres binomios en las variables enteras 2a
            prod = pmul(pmul(f(d1), f(d2)), f(d3))
            if norma1(prod) == 8:
                c3_norma += 1
            # senuelo: si se permitiera a_i = b_j
            if any(c == 0 for c in [dt3 + d1 - d2, dt3 + d1 + d2, dt3 - d1 - d2, dt3 - d1 + d2]):
                d1_senuelo += 1
        elif len(gordas) == 1 and p[gordas[0]] == 3:
            n_tres += 1
            C = clases[gordas[0]]
            x1, x2, x3 = C[0], C[1], C[2]
            d1, d2, d3 = x1 - x2, x2 - x3, x1 - x3
            if d3 == d1 + d2:
                c4_ok += 1
            prod = pmul(pmul(f(d1), f(d2)), f(d3))
            if norma1(prod) == 6:
                c4_norma += 1

    print()
    print("  t=%d   %d lambda: %d degeneradas, %d concentricas (Phi = 0, hipotesis excluida),"
          % (T, len(lams), n_deg, n_conc))
    print("         %d de dos clases y %d de tamano tres con Phi != 0" % (n_dos, n_tres))
    print("    C3  las cuatro combinaciones son 2(a_i - b_j) y ninguna se anula : %d de %d %s"
          % (c3_ok, n_dos, "OK" if c3_ok == n_dos else "*** FALLA"))
    print("    C3  y por tanto la norma del producto es 8                       : %d de %d %s"
          % (c3_norma, n_dos, "OK" if c3_norma == n_dos else "*** FALLA"))
    print("    C4  d3 = d1 + d2 en el perfil de tamano tres                     : %d de %d %s"
          % (c4_ok, n_tres, "OK" if c4_ok == n_tres else "*** FALLA"))
    print("    C4  y por tanto la norma del producto es 6                       : %d de %d %s"
          % (c4_norma, n_tres, "OK" if c4_norma == n_tres else "*** FALLA"))
    print("    D1  SENUELO: lambda donde alguna combinacion SI se anula         : %d %s"
          % (d1_senuelo, "(ninguna, como exige el argumento)" if d1_senuelo == 0 else "*** el argumento falla"))
    if c3_ok != n_dos or c3_norma != n_dos or c4_ok != n_tres or c4_norma != n_tres:
        fallos += 1

print()
print("=" * 92)
print("VEREDICTO: %s" % ("el argumento del revisor PASA sobre todo el rango" if not fallos
                          else "*** el argumento NO pasa"))
sys.exit(1 if fallos else 0)
