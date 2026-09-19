# -*- coding: utf-8 -*-
# SOCRATICO SOBRE EL BARRIDO DE TIPO A: no cuanto vale el indice, sino POR QUE lo hay.
#
# P1  ¿El indice aparece SOLO cuando las cuentas del segmento de exponentes estan EQUILIBRADAS
#     modulo q'?  Esa es la causa comun que sale de leer NPP25 (su centralizador minimo es
#     sum n_i^2 minimo <=> |n_i - n_j| <= 1) y de nuestro parrafo de "cuenta por punto medio"
#     (cuenta constante => el primer momento colapsa en el diente de sierra).
#     Se cruza, casilla a casilla y primo a primo: equilibrado SI/NO contra p | indice SI/NO.
#     Las cuatro celdas se imprimen; una equivalencia solo se cree con las cuatro.
#
# P2  ¿Por que la SOMBRA?  A_SU (real) y A_GL (complejo) se diferencian en que cada generador
#     lleva un giro: e_k(x^{SU}) = zeta_{2q}^{-k(n-1)} e_k(x^{GL}).  Hipotesis H1: A_GL es
#     A_SU con una raiz de la unidad anadida, es decir A_GL = A_SU[e_n] con e_n = zeta_q^{n(n-1)/2}.
#     Se comprueba las dos inclusiones por reticulos, no por creencia.
#
# P3  ¿Que ES el generador de A_GL?  e_k(1, z, ..., z^{n-1}) = z^{k(k-1)/2} * [n sobre k]_z,
#     el binomio gaussiano en una raiz de la unidad.  Se verifica la identidad antes de usarla.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: docker ... sage socratico_tipoA.sage
import time


def reticulo(gens, K):
    """Z-span cerrado multiplicativamente de Z[gens] dentro de K, en la base de potencias."""
    L = matrix(ZZ, [K(1).list()]).hermite_form(include_zero_rows=False)
    for _ in range(40):
        filas = [list(f) for f in L.rows()]
        for f in L.rows():
            el = K(list(f))
            for g in gens:
                filas.append((el * g).list())
        L2 = matrix(ZZ, filas).hermite_form(include_zero_rows=False)
        if L2 == L:
            return L
        L = L2
    return L


def gens_de(exps, M, r):
    K = CyclotomicField(M)
    z = K.gen()
    xs = [z ** (e % M) for e in exps]
    S = PolynomialRing(K, 'T')
    T = S.gen()
    f = prod([T - x for x in xs])
    n = len(xs)
    g = [(-1) ** k * f.list()[n - k] for k in range(1, r + 1)]
    return K, [y for y in g if y != 0]


def indice(L):
    return prod([x for x in L.elementary_divisors() if x != 0])


def exps_su(n):
    if n % 2 == 1:
        return [(n + 1 - 2 * j) // 2 for j in range(1, n + 1)]
    return [(n + 1 - 2 * j) for j in range(1, n + 1)]


print("=" * 78)
print("P3  ¿Que es el generador de A_GL?  e_k(1,z,...,z^{n-1}) == z^{k(k-1)/2} [n,k]_z ?")
print("=" * 78)
fallos = 0
casos = 0
for n in range(2, 8):
    for q in (7, 9, 10, 12, 15):
        K = CyclotomicField(q)
        z = K.gen()
        xs = [z ** j for j in range(n)]
        S = PolynomialRing(K, 'T')
        T = S.gen()
        f = prod([T - x for x in xs])
        for k in range(1, n + 1):
            ek = (-1) ** k * f.list()[n - k]
            gauss = K(gaussian_binomial(n, k, z)) if hasattr(sage.all, "gaussian_binomial") else None
            if gauss is None:
                num = prod([1 - z ** (n - i) for i in range(k)])
                den = prod([1 - z ** (i + 1) for i in range(k)])
                gauss = num / den if den != 0 else None
            if gauss is None:
                continue
            casos += 1
            if ek != z ** (k * (k - 1) / 2) * gauss:
                fallos += 1
                if fallos <= 4:
                    print("   DISCREPA n=%d q=%d k=%d" % (n, q, k))
print("   comprobaciones: %d | discrepancias: %d" % (casos, fallos))
if fallos == 0:
    print("   la identidad se sostiene: los generadores de A_GL SON binomios gaussianos girados.")

print()
print("=" * 78)
print("P1  ¿El indice solo aparece con las cuentas EQUILIBRADAS modulo q'?")
print("=" * 78)
print("%-3s %-4s %-3s %-4s %-26s %-11s %s" % ("n", "q", "p", "q'", "cuentas mod q'", "equilibrado", "p|indice"))
tab = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
testigos = []
for n in range(3, 8):
    for q in range(3, 41):
        M = q if n % 2 == 1 else 2 * q
        if euler_phi(M) > 40:
            continue
        K, g = gens_de(exps_su(n), M, n - 1)
        idx = indice(reticulo(g, K))
        for p, _ in factor(q):
            v = valuation(q, p)
            qp = q // p ** v
            if qp <= 2:
                continue
            cuentas = [0] * qp
            for e in exps_su(n):
                cuentas[e % qp] += 1
            # EQUILIBRIO EN SU SOPORTE: para n par los exponentes son todos impares, asi que el
            # segmento solo puede tocar un coset y exigir |n_i - n_j| <= 1 sobre TODAS las clases
            # mide lo que no es.  La condicion correcta es dentro de las clases que toca.
            vivas = [c for c in cuentas if c > 0]
            equil = (max(vivas) - min(vivas)) <= 1
            divide = (idx % p == 0)
            tab[(equil, divide)] += 1
            if divide or (equil and qp <= 12):
                print("%-3d %-4d %-3d %-4d %-26s %-11s %s" % (
                    n, q, p, qp, str(cuentas)[:26], "SI" if equil else "no",
                    "SI" if divide else "no"))
            if divide and not equil:
                testigos.append((n, q, p, qp, cuentas))
print()
print("   LAS CUATRO CELDAS (equilibrado, p|indice):")
print("     SI/SI %4d     SI/no %4d" % (tab[(True, True)], tab[(True, False)]))
print("     no/SI %4d     no/no %4d" % (tab[(False, True)], tab[(False, False)]))
if tab[(False, True)] == 0:
    print("   => equilibrio es NECESARIO para que p divida al indice, en toda la muestra.")
else:
    print("   => HAY CONTRAEJEMPLOS al 'necesario':")
    for t in testigos[:6]:
        print("      n=%d q=%d p=%d q'=%d cuentas=%s" % t)
if tab[(True, False)] > 0:
    print("   => y no es SUFICIENTE: %d casos equilibrados sin que p divida." % tab[(True, False)])

print()
print("=" * 78)
print("P2  LA SOMBRA.  Primero lo que la ronda anterior enseno: para n PAR los dos anillos ni")
print("    siquiera viven en el mismo cuerpo -- el elemento SU tiene exponentes semienteros y")
print("    vive en Q(zeta_2q), el GL en Q(zeta_q).  Comparar sus indices es comparar ordenes de")
print("    cuerpos distintos.  Asi que la pregunta limpia es con n IMPAR, donde ambos estan en")
print("    Z[zeta_q], y el giro por generador es e_k^GL = zeta_q^{k(n-1)/2} e_k^SU.")
print("    Hipotesis: A_GL = A_SU[w] con w = zeta_q^{(n-1)/2}, el giro del PRIMER generador.")
print("%-3s %-4s %-8s %-8s %-9s %-11s %-11s %s" % (
    "n", "q", "idx SU", "idx GL", "SU<=GL?", "GL=SU[w1]?", "GL=SU[en]?", "orden w1"))
for n in range(3, 8, 2):
    for q in (8, 9, 10, 12, 15, 16, 20, 21, 24, 25):
        M = 2 * q if n % 2 == 0 else q
        Mbig = lcm(M, q)
        if euler_phi(Mbig) > 40:
            continue
        Kb = CyclotomicField(Mbig)
        zb = Kb.gen()
        def sube(exps, Mloc, r):
            fac = Mbig // Mloc
            xs = [zb ** ((e * fac) % Mbig) for e in exps]
            S = PolynomialRing(Kb, 'T')
            T = S.gen()
            f = prod([T - x for x in xs])
            nn = len(xs)
            return [(-1) ** k * f.list()[nn - k] for k in range(1, r + 1)]
        gsu = [y for y in sube(exps_su(n), M, n - 1) if y != 0]
        ggl = [y for y in sube([j for j in range(n)], q, n) if y != 0]
        Lsu = reticulo(gsu, Kb)
        Lgl = reticulo(ggl, Kb)
        fac = Mbig // q
        w1 = zb ** ((fac * ((n - 1) // 2)) % Mbig)
        en = zb ** ((fac * (n * (n - 1) // 2)) % Mbig)
        Lw1 = reticulo(gsu + [w1], Kb)
        Len = reticulo(gsu + [en], Kb)
        def contiene(A, B):
            return all(matrix(ZZ, list(A.rows()) + [list(r)]).hermite_form(
                include_zero_rows=False) == A for r in B.rows())
        print("%-3d %-4d %-8s %-8s %-9s %-11s %-11s %s" % (
            n, q, indice(Lsu), indice(Lgl),
            "si" if contiene(Lgl, Lsu) else "NO",
            "si" if Lw1 == Lgl else "NO",
            "si" if Len == Lgl else "NO",
            w1.multiplicative_order()))
