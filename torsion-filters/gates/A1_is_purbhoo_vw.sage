# -*- coding: utf-8 -*-
# A1 NO ES NUESTRA: es el Teorema 2.5 de Purbhoo-van Willigenburg (2008).  Aqui esta el diccionario,
# y esta comprobado a maquina, porque el sitio donde esto suele fallar es MI mapa de indices.
#
# ================================================================================================
# LO QUE DICEN ELLOS, verbatim, leido del PDF y no de un resumen.
#
# Kevin Purbhoo & Stephanie van Willigenburg, "On tensor products of polynomial representations",
# Canad. Math. Bull. 51 (2008) 584-592, arXiv:0706.3251.  Con  lambda^- := lambda / (lambda_n)^n
# (quitar todas las columnas de altura n):
#
#   TEOREMA 2.5.  "phi_lambda (x) phi_mu = phi_nu (x) phi_rho as representations of GL(n) if and
#   only if  lambda_n + mu_n = nu_n + rho_n  and  {lambda^-, mu^-} = {nu^-, rho^-} as multisets."
#
# y su prueba empieza: "We will show that (5) s_lambda(x_1,..,x_n) s_mu(x_1,..,x_n) =
# s_nu(x_1,..,x_n) s_rho(x_1,..,x_n) if and only if ..." -- o sea es literalmente una IGUALDAD DE
# PRODUCTOS DE DOS FUNCIONES DE SCHUR en n variables, sin teoria de representaciones por medio.
#
# C. S. Rajan, Annals of Math. 160 (2004) 683-704, arXiv:math/0305417, da la version de n factores:
#   TEOREMA 3.  "Let G = GL(r).  Suppose V = V_{lambda_1} (x) .. (x) V_{lambda_n} and
#   W = V_{mu_1} (x) .. (x) V_{mu_m} are tensor products of irreducible representations with NONZERO
#   highest weights .. Then n = m and there is a permutation tau .. such that
#   V_{lambda_i} = V_{mu_tau(i)} (x) det^{alpha_i}, for some integers alpha_i."
# *** ESTA FRASE LA TENIA AL REVES, y la debo a la auditoria socratica.  Escribi que la hipotesis
# "nonzero highest weights" de Rajan es la que cubre nuestro caso degenerado.  Es al reves: esa
# hipotesis EXCLUYE los factores triviales, asi que Rajan NO cubre el caso degenerado.  Quien lo
# cubre es PvW, que fija exactamente dos factores y lleva la cuenta de la columna del determinante.
# Y el caso degenerado es real y frecuente: atil es vacia siempre que H sea un bloque de enteros
# consecutivos, y astar siempre que lo sea L. ***
#
# Y UNA TRAMPA QUE HAY QUE DECIR: el atajo "el anillo de Laurent es un DFU, luego listo" es FALSO.
# Las funciones de Schur NO son irreducibles: s_{(3)}(x,y) = h_3 = (x+y)(x^2+y^2)... el contenido
# esta en PvW/Rajan, no en factorizacion unica.
#
# ------------------------------------------------------------------------------------------------
# EL DICCIONARIO, que es lo unico nuestro aqui, y es de dos lineas.
#
# Con H = h_1 > .. > h_r, L = l_1 > .. > l_r, delta = (r-1, .., 0):
#   a_H(z) = s_alpha(z) Vdm(z),   Vdm = Vandermonde (NO confundir con V, el valor de empate),                alpha = H - delta.
#   a_L(1/z): multiplicando la columna j por z_j^{l_1} y dando la vuelta a las r filas,
#             a_L(1/z) = (-1)^{r(r-1)/2} (prod_j z_j)^{-l_1} a_{L*}(z),   L* = reverse(l_1 - L),
#             o sea a_L(1/z) = +- (prod z)^{-l_1} s_{astar}(z) Vdm(z),    astar = L* - delta,
#             y astar ya esta NORMALIZADA: su ultima parte es 0 porque min(L*) = 0.
# Luego, sacando tambien las columnas llenas de alpha (alpha_r = h_r),
#
#       P(H,L) = a_H(z) a_L(1/z) = +- (prod z)^{h_r - l_1} Vdm(z)^2 s_{atil}(z) s_{astar}(z),
#       atil = alpha - (alpha_r)^r    (normalizada),      astar = L* - delta.
#
# Aplicando el Teorema 2.5 al producto s_{atil} s_{astar} (la potencia de (prod z) es justo su
# "lambda_n + mu_n"), sale que el invariante completo de T = H u L es el PAR
#
#       INV(T) = ( h_r - l_1 ,  {atil, astar} como MULTICONJUNTO ).
#
# Y las dos maneras de tener el mismo INV son exactamente las dos que habiamos medido:
#   * las dos partes iguales en el mismo orden  ->  T' = T + m   (traslacion),
#   * las dos partes intercambiadas             ->  T' = c - T   (reflexion),
# porque intercambiar atil con astar es intercambiar el papel de H y el de L, y eso es la reflexion.
# ------------------------------------------------------------------------------------------------
#
# VERIFICADO AQUI, cada cosa capaz de fallar:
#   D1  el diccionario como identidad de polinomios: P(H,L) = a_H(z) a_L(1/z) contra
#       +- (prod z)^{h_r-l_1} Vdm(z)^2 s_{atil}(z) s_{astar}(z), evaluado sobre GF(p) en puntos
#       aleatorios.  Si mi mapa de indices esta mal, muere aqui.
#   D2  EL CONTROL QUE DECIDE: INV(T) debe agrupar los T EXACTAMENTE igual que la clave real de
#       P(T).  Se comparan las dos particiones del conjunto de todos los T, en las dos direcciones.
#   D3  SENUELO: quitar el entero h_r - l_1 y quedarse solo con el multiconjunto.  Debe SOBRE-
#       disparar, o esa mitad del criterio de PvW no esta haciendo nada.
#   D4  SENUELO: usar el par ORDENADO (atil, astar) en vez del multiconjunto.  Debe INFRA-disparar,
#       perdiendo justo las reflexiones.
#   D5  no vacuidad: tiene que haber colisiones de los dos tipos, traslacion y reflexion.
#
# CONSECUENCIA, y hay que decirla en la primera linea de cualquier escrito: la cadena
#       Phi_t = 0  =>  ... =>  condicion (ii)
# queda CERRADA, y el eslabon que faltaba lo pone Purbhoo-van Willigenburg, no nosotros.
#
# UN DEFECTO MIO, cazado por D1 en la primera pasada: compare contra el determinante
# COMPLETO 2r x 2r A(T) en vez de contra P(H,L) = a_H(z) a_L(1/z).  A(T) tiene a P solo
# como su parte de grado maximo, asi que fallaba en los 495 puntos.  D1 estaba puesto
# exactamente para eso y disparo a la primera.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

p = next_prime(10 ** 9)
F = GF(p)
print("field GF(%d)" % p)


def perm_sign(q):
    n = len(q)
    seen = [False] * n
    s = 1
    for i in range(n):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = q[j]
            ln += 1
        if ln % 2 == 0:
            s = -s
    return s


def P_dict(T, r):
    """a_H(z) a_L(1/z) as {exponent tuple: coeff}, T decreasing, H = T[:r], L = T[r:]."""
    D = {}
    n = 2 * r
    for a in itertools.permutations(range(r)):
        for b in itertools.permutations(range(r)):
            q = [0] * n
            e = [0] * r
            for i in range(r):
                q[i] = 2 * a[i]
                e[a[i]] += T[i]
            for i in range(r):
                q[r + i] = 2 * b[i] + 1
                e[b[i]] -= T[r + i]
            k = tuple(e)
            D[k] = D.get(k, 0) + perm_sign(q)
    return dict((k, v) for k, v in D.items() if v != 0)


def key_of(T, r):
    D = P_dict(T, r)
    a = tuple(sorted(D.items()))
    b = tuple(sorted((k, -v) for k, v in D.items()))
    return min(a, b)


def inv_of(T, r):
    """the Purbhoo-van Willigenburg invariant of T, through the dictionary."""
    H, L = list(T[:r]), list(T[r:])
    delta = [r - 1 - i for i in range(r)]
    alpha = [H[i] - delta[i] for i in range(r)]
    atil = tuple(x - alpha[-1] for x in alpha)
    Lstar = sorted((L[0] - x for x in L), reverse=True)     # reverse(l_1 - L)
    astar = tuple(Lstar[i] - delta[i] for i in range(r))
    assert astar[-1] == 0, "astar must already be normalised"
    return (H[-1] - L[0], frozenset([atil, astar]) if atil != astar else frozenset([atil]),
            tuple(sorted([atil, astar])))


def schur(lam, xs):
    """s_lam(xs) by the bialternant, over F."""
    r = len(xs)
    lam = list(lam) + [0] * (r - len(lam))
    num = matrix(F, r, r, lambda i, j: xs[j] ** (lam[i] + r - 1 - i)).det()
    den = matrix(F, r, r, lambda i, j: xs[j] ** (r - 1 - i)).det()
    return num / den


# ---------------------------------------------------------------- D1 -----------------------------
print("")
print("=" * 104)
print("D1  the dictionary as a polynomial identity, over GF(p) at random points")
print("=" * 104)
print("")
print("     r   sets tested   points   D1 bad")
print("  " + "-" * 100)

set_random_seed(20260812)
d1tot = 0
for r, M in ((2, 12), (3, 10), (4, 9)):
    n = 2 * r
    nb = ns = 0
    for T in itertools.combinations(range(M + 1), n):
        Td = tuple(sorted(T, reverse=True))
        ns += 1
        if ns > 60:
            break
        H, L = list(Td[:r]), list(Td[r:])
        delta = [r - 1 - i for i in range(r)]
        alpha = [H[i] - delta[i] for i in range(r)]
        atil = [x - alpha[-1] for x in alpha]
        Lstar = sorted((L[0] - x for x in L), reverse=True)
        astar = [Lstar[i] - delta[i] for i in range(r)]
        for _ in range(3):
            xs = []
            while len(xs) < r:
                x = F.random_element()
                if x != 0 and all(x != y for y in xs) and all(x * y != 1 for y in xs) and x ** 2 != 1:
                    xs.append(x)
            # P(H,L) = a_H(z) a_L(1/z), NOT the full 2r x 2r alternant A(T) -- D1 caught me
            # comparing against A(T) on the first run, which is P only in its top-degree part.
            aH = matrix(F, r, r, lambda i, j: xs[j] ** H[i]).det()
            aLinv = matrix(F, r, r, lambda i, j: xs[j] ** (-L[i])).det()
            lhs = aH * aLinv
            V = matrix(F, r, r, lambda i, j: xs[j] ** (r - 1 - i)).det()
            rhs = prod(xs) ** (H[-1] - L[0]) * V ** 2 * schur(atil, xs) * schur(astar, xs)
            if lhs != rhs and lhs != -rhs:
                nb += 1
    d1tot += nb
    print("  %4d %13d %8d %8d" % (r, min(ns, 60), 3 * min(ns, 60), nb))

print("")
if d1tot:
    print("  D1 FAILED -- my index map is wrong, stop.")
    raise SystemExit(1)
print("  D1 PASS: P(H,L) = +- (prod z)^{h_r - l_1} Vdm(z)^2 s_atil(z) s_astar(z).")
print("  El +- : D1 acepta cualquiera de los dos signos en cada instancia, asi que NO establece")
print("  que el signo sea funcion de r sola.  Medido aparte por el auditor: es (-1)^{r(r-1)/2},")
print("  constante en T, r = 1..4.  Da igual para el resultado -- se cancela entre T_A y T_B --")
print("  pero el veredicto no debe decir mas de lo que el control prueba.")

# ---------------------------------------------------------------- D2, D3, D4, D5 -----------------
print("")
print("=" * 104)
print("D2  does the PvW invariant group the sets exactly as P does?     D3/D4 the two decoys")
print("=" * 104)
print("")
print("     r    M   sets   P-classes  INV-classes | D2 bad | D3 decoy (no det)  D4 decoy (ordered)")
print("  " + "-" * 100)

d2tot = 0
D3 = D4 = 0
NTR = NRE = 0
for r, M in ((2, 13), (3, 11), (4, 9)):
    n = 2 * r
    byP = {}
    byI = {}
    byI3 = {}
    byI4 = {}
    allT = []
    for T in itertools.combinations(range(M + 1), n):
        Td = tuple(sorted(T, reverse=True))
        allT.append(Td)
        kP = key_of(Td, r)
        iv = inv_of(Td, r)
        byP.setdefault(kP, []).append(Td)
        byI.setdefault((iv[0], iv[2] if iv[2][0] <= iv[2][1] else (iv[2][1], iv[2][0])),
                       []).append(Td)
        byI3.setdefault(iv[2], []).append(Td)                      # decoy: no det integer
        byI4.setdefault((iv[0], inv_of(Td, r)[2]), []).append(Td)  # placeholder, replaced below
    # rebuild the two decoys cleanly
    byI3 = {}
    byI4 = {}
    for Td in allT:
        H, L = list(Td[:r]), list(Td[r:])
        delta = [r - 1 - i for i in range(r)]
        alpha = [H[i] - delta[i] for i in range(r)]
        atil = tuple(x - alpha[-1] for x in alpha)
        Lstar = sorted((L[0] - x for x in L), reverse=True)
        astar = tuple(Lstar[i] - delta[i] for i in range(r))
        byI3.setdefault(tuple(sorted([atil, astar])), []).append(Td)        # drop the det integer
        byI4.setdefault((H[-1] - L[0], atil, astar), []).append(Td)         # ordered pair
    # the true invariant
    byI = {}
    for Td in allT:
        H, L = list(Td[:r]), list(Td[r:])
        delta = [r - 1 - i for i in range(r)]
        alpha = [H[i] - delta[i] for i in range(r)]
        atil = tuple(x - alpha[-1] for x in alpha)
        Lstar = sorted((L[0] - x for x in L), reverse=True)
        astar = tuple(Lstar[i] - delta[i] for i in range(r))
        byI.setdefault((H[-1] - L[0], tuple(sorted([atil, astar]))), []).append(Td)

    def parts(d):
        return set(frozenset(v) for v in d.values())

    pP, pI, p3, p4 = parts(byP), parts(byI), parts(byI3), parts(byI4)
    bad = 0 if pP == pI else len(pP.symmetric_difference(pI))
    d2tot += bad
    D3 += 0 if p3 == pP else len(p3.symmetric_difference(pP))
    D4 += 0 if p4 == pP else len(p4.symmetric_difference(pP))
    # non-vacuity: count collision types inside the true classes
    for cls in pP:
        v = sorted(cls)
        if len(v) < 2:
            continue
        A = v[0]
        for B in v[1:]:
            if any(tuple(sorted((x + m for x in A), reverse=True)) == B
                   for m in range(-2 * M, 2 * M + 1) if m != 0):
                NTR += 1
            if any(tuple(sorted((c - x for x in A), reverse=True)) == B
                   for c in range(0, 4 * M + 1)):
                NRE += 1
    print("  %4d %4d %6d %11d %12d | %6d | %17d %18d"
          % (r, M, len(allT), len(pP), len(pI), bad,
             0 if p3 == pP else len(p3.symmetric_difference(pP)),
             0 if p4 == pP else len(p4.symmetric_difference(pP))))

print("")
print("  D2 mismatched classes between the P-partition and the PvW-invariant partition: %d" % d2tot)
print("  D3 decoy (drop the det integer h_r - l_1): %d classes differ -- it must OVER-merge." % D3)
print("  D4 decoy (ordered pair instead of multiset): %d classes differ -- it must UNDER-merge." % D4)
print("  D5 non-vacuity: collisions realised by a translation %d, by a reflection %d." % (NTR, NRE))
print("")
if d2tot == 0 and D3 > 0 and D4 > 0 and NTR > 0 and NRE > 0:
    print("  A1 IS PURBHOO-VAN WILLIGENBURG THEOREM 2.5.  The dictionary is exact, both halves of")
    print("  their criterion are load-bearing, and the two collision types are exactly the")
    print("  translation (their det condition) and the reflection (their multiset).  A1 is NOT ours")
    print("  and must be cited, not claimed.  With it the chain closes:  Phi_t = 0  =>  (ii).")
else:
    print("  SOMETHING FAILED -- read the columns.")
print("")
print("DONE")

print("")
print("=" * 104)
print("F13  los dos pasos que estaban sin escribir entre  P(T_A) = +-P(T_B)  y  PvW")
print("=" * 104)
print("")
print("  (a) el +- se resuelve a +: los dos lados son un monomio por un polinomio de coeficientes no")
print("      negativos (producto de dos Schur), asi que un signo menos es imposible.")
print("  (b) (prod z)^d s_atil hay que reabsorberlo como s_{atil + (d^r)}, para que ese d SEA el")
print("      'lambda_n + mu_n' de PvW.  Requiere d = h_r - l_1 >= 1, cierto porque T es")
print("      estrictamente decreciente y h_r > l_1.  Ahi es donde aterriza la mitad del determinante")
print("      de su criterio, y sin decirlo el diccionario queda cojo.")
print("")
print("DONE")
