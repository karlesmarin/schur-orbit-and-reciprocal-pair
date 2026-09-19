# -*- coding: utf-8 -*-
# VERIFICACION DE LA DEMOSTRACION DEL HUECO (ii):  B^- = theta R.
#
# La demostracion, paso a paso, y cada paso se comprueba aqui en celdas de verdad.
#
#  (P1) REPRESENTANTE.  Por lem:div (multiplicar por un factor monico ENTERO no cambia el anillo):
#         n = Nr    : f_n(T) = (T-1) g(T) porque j=0 esta en el segmento => B = ring(g),
#                     y g es el segmento {1,...,Nr-1}.
#         n = Nr+1  : se deja {0,...,Nr} tal cual.
#       En los dos casos el representante S es SIMETRICO alrededor de Nr/2.
#
#  (P2) CENTRADO PURO p.  mu := zeta^{Nr/2} (2 invertible mod q porque q es impar) cumple
#       mu = rho^{Nr/2}: es una raiz de la unidad de orden potencia de p, SIN parte omega.
#       Centrando, e_k(S) = mu^k a_k con a_k REAL.
#
#  (P3) LA PRIMERA ES UNIDAD.  La reduccion de e_1(S) modulo el maximal de B es
#         -1  si n = Nr     (suma de todas las clases, menos el 1 que se quito)
#         +1  si n = Nr+1
#       En ambos casos no nula en F_p: e_1(S) es UNIDAD DE B.
#
#  (P4) mu ESTA EN B.  mu^2 = e_1(S) / iota(e_1(S)) pertenece a B por (P3); mu tiene orden impar
#       (divide q), luego mu = (mu^2)^{(ord+1)/2} pertenece a B.
#
#  (P5) LA EXTENSION.  B = A[mu] con A = Z_p[a_k] real; R := B^iota = A[tau], tau = mu + mu^{-1};
#       B = R (+) mu R = R (+) (mu - mu^{-1}) R, luego B^- = (mu - mu^{-1}) R.
#
#  (P6) theta GENERA.  theta := e_1 - iota(e_1) NO cambia al pasar de {0..Nr-1} a {1..Nr-1}
#       (la constante se va en la resta), y theta = a_1 (mu - mu^{-1}) con a_1 = e_1(S)/mu UNIDAD.
#       Luego theta R = (mu - mu^{-1}) R = B^-.   HUECO (ii) CERRADO.
#
# Se comprueba ademas si tau pertenece a A, es decir si R = A o R es estrictamente mayor.
#
# Authors: Carles Marin, Claude (AI assistant).
import sys


def vec(K, x, D):
    l = K(x).list()
    return l + [0] * (D - len(l))


def cierre(K, exps, q, D, ngen=None, vueltas=80):
    """reticulo entero del anillo generado por las e_k de los exponentes dados."""
    S = PolynomialRing(K, 'T')
    T = S.gen()
    z = K.gen()
    m = len(exps)
    f = prod([T - z ** (e % q) for e in exps])
    co = f.list()
    ng = ngen or m
    g = [(-1) ** k * co[m - k] for k in range(1, ng + 1)]
    g = [y for y in g if y != 0]
    L = matrix(ZZ, [vec(K, K(1), D)]).hermite_form(include_zero_rows=False)
    for _ in range(vueltas):
        filas = [list(fi) for fi in L.rows()]
        for fi in L.rows():
            el = K(list(fi))
            for gg in g:
                filas.append(vec(K, el * gg, D))
        L2 = matrix(ZZ, filas).hermite_form(include_zero_rows=False)
        if L2 == L:
            return L
        L = L2
    raise RuntimeError("el cierre no estabilizo")


def matriz_iota(K, D):
    z = K.gen()
    zi = z ** -1
    filas = []
    acc = K(1)
    for i in range(D):
        filas.append(vec(K, acc, D))
        acc = acc * zi
    return matrix(QQ, filas)


def p_ent(v, p):
    return all(c.denominator() % p != 0 for c in v)


def dentro(K, Linv, x, D, p):
    return p_ent(vector(QQ, vec(K, x, D)) * Linv, p)


def es_unidad_p(K, x, p):
    """x es unidad en O tensor Z_p."""
    if x == 0:
        return False
    return all(ex == 0 for (pr, ex) in K.ideal(x).factor() if pr.smallest_integer() == p)


def sub_fijo(L, Iota, D):
    A = L * (Iota - identity_matrix(QQ, D))
    Bk = A.kernel().basis_matrix()
    Bk = matrix(ZZ, Bk * Bk.denominator()).saturation()
    return matrix(ZZ, Bk * L).hermite_form(include_zero_rows=False)


def anti(L, Iota, D):
    A = L * (Iota + identity_matrix(QQ, D))
    Bk = A.kernel().basis_matrix()
    Bk = matrix(ZZ, Bk * Bk.denominator()).saturation()
    return matrix(ZZ, Bk * L).hermite_form(include_zero_rows=False)


def mult(K, L, x, D):
    return matrix(ZZ, [vec(K, K(list(fi)) * x, D) for fi in L.rows()]).hermite_form(
        include_zero_rows=False)


def igual_p(L1, L2, p):
    """L1 = L2 localmente en p."""
    try:
        S1 = matrix(QQ, L2).solve_left(matrix(QQ, L1))
        S2 = matrix(QQ, L1).solve_left(matrix(QQ, L2))
    except Exception:
        return False
    return all(p_ent(f, p) for f in S1.rows()) and all(p_ent(f, p) for f in S2.rows())


# Se incluyen A PROPOSITO los regimenes que la compuerta de capas excluia, para ver si la
# demostracion los necesita: E = phi(p^v) = 2 (v=1, p=3), r NO primo, y N > 1.
#
# Y AL FINAL, el regimen p | N.  Releyendo la demostracion paso a paso, p no divide N NO se usa
# en ninguno de los seis pasos: el cero de sum_{j<Nr} omega^j viene de sumar todas las clases, no
# de N.  Si eso es cierto, el teorema vale igual y lo unico que cambia es
#     s = v_p(mu - mu^{-1}) = v_p(zeta^{Nr} - 1) = p^{v_p(N)},
# que es exactamente el desplazamiento que `puente_general_s.sage` habia MEDIDO sin demostrar.
# Estas celdas lo deciden.
CELDAS = [(3, 2, 5, 1), (3, 2, 5, 2), (3, 2, 7, 1), (3, 2, 11, 1), (3, 2, 13, 1),
          (5, 1, 7, 1), (5, 1, 7, 2), (5, 1, 11, 1), (7, 1, 9, 1), (7, 1, 11, 1),
          (3, 3, 5, 1), (11, 1, 9, 1),
          (3, 1, 5, 1), (3, 1, 5, 2), (3, 1, 7, 1), (3, 1, 11, 1), (3, 1, 13, 1),
          (3, 1, 25, 1), (5, 1, 9, 1), (5, 1, 21, 1), (7, 1, 15, 1), (3, 1, 35, 1),
          (5, 1, 27, 1), (3, 2, 5, 4), (5, 1, 7, 3), (7, 1, 9, 2),
          # ---- p | N ----
          (3, 2, 5, 3), (3, 2, 7, 3), (3, 2, 11, 3), (3, 2, 13, 3), (3, 3, 5, 3),
          (3, 3, 5, 9), (5, 2, 7, 5), (3, 2, 5, 6), (7, 2, 5, 7)]

print("%-5s %-3s %-4s %-5s %-4s %-8s %-8s %-8s %-8s %-8s %-8s %-6s %s" % (
    "q", "p", "r", "n", "P", "e1(S)u", "mu en B", "a_1 u", "B=A[mu]", "B^-=thR", "R=A",
    "v(th)", "veredicto"))
tot = {"u": 0, "mu": 0, "a1": 0, "Amu": 0, "th": 0, "RA": 0, "n": 0, "vP": 0}
for (p, v, r, N) in CELDAS:
    q = p ** v * r
    K = CyclotomicField(q)
    z = K.gen()
    D = euler_phi(q)
    if D > 80:
        continue
    Iota = matriz_iota(K, D)
    inv2 = inverse_mod(2, q)
    mu = z ** ((N * r * inv2) % q)
    for off in (0, 1):
        n = N * r + off
        if not (2 <= n < q):
            continue
        tot["n"] += 1
        # (P1) el representante S, y el centrado
        S = list(range(1, N * r)) if off == 0 else list(range(0, N * r + 1))
        cent = [(e - N * r * inv2) % q for e in S]
        LB = cierre(K, list(range(n)), q, D)          # el anillo de verdad, segmento {0..n-1}
        LS = cierre(K, S, q, D)                       # el representante
        LA = cierre(K, cent, q, D)                    # el centrado (real)
        LBinv = LB.inverse()
        rep_ok = igual_p(LB, LS, p)                   # lem:div
        e1S = sum([z ** (e % q) for e in S])
        u_ok = es_unidad_p(K, e1S, p)                 # (P3)
        mu_ok = dentro(K, LBinv, mu, D, p)            # (P4)
        a1 = e1S * mu ** -1
        a1_ok = es_unidad_p(K, a1, p) and (a1 - K(list(vector(QQ, vec(K, a1, D)) * Iota)) == 0)
        LAmu = cierre(K, cent, q, D)
        LAmu = matrix(ZZ, LA.stack(mult(K, LA, mu, D))).hermite_form(include_zero_rows=False)
        Amu_ok = igual_p(LB, LAmu, p)                 # (P5) B = A[mu] (como A + mu A)
        LR = sub_fijo(LB, Iota, D)
        LBm = anti(LB, Iota, D)
        e1 = sum([z ** (j % q) for j in range(n)])
        th = e1 - K(list(vector(QQ, vec(K, e1, D)) * Iota))
        th_ok = igual_p(LBm, mult(K, LR, th, D), p)   # (P6)
        RA_ok = igual_p(LR, LA, p)                    # ?es tau real ya en A?
        # la valuacion de theta, normalizada con v(pi) = 1:  se compara con P = p^{v_p(N)}
        P = p ** valuation(N, p)
        vth = sum([ex for (pr, ex) in K.ideal(th).factor() if pr.smallest_integer() == p])
        npr = len([1 for (pr, ex) in K.ideal(p).factor()])
        vth1 = vth // npr if npr else -1
        vP_ok = (vth1 == P)
        for kk, vv in (("u", u_ok), ("mu", mu_ok), ("a1", a1_ok), ("Amu", Amu_ok),
                       ("th", th_ok), ("RA", RA_ok), ("vP", vP_ok)):
            tot[kk] += 1 if vv else 0
        ver = "ok" if (rep_ok and u_ok and mu_ok and a1_ok and Amu_ok and th_ok and vP_ok) \
            else "*** FALLA ***"
        if not rep_ok:
            ver += " (lem:div!)"
        print("%-5d %-3d %-4d %-5d %-4d %-8s %-8s %-8s %-8s %-8s %-8s %-6s %s" % (
            q, p, r, n, P, "SI" if u_ok else "NO", "SI" if mu_ok else "NO",
            "SI" if a1_ok else "NO", "SI" if Amu_ok else "NO",
            "SI" if th_ok else "NO", "SI" if RA_ok else "no",
            "%d%s" % (vth1, "" if vP_ok else "!"), ver))
        sys.stdout.flush()
print()
print("celdas: %d" % tot["n"])
for kk, txt in (("u", "(P3) e_1(S) unidad de O_p"), ("mu", "(P4) mu en B"),
                ("a1", "     a_1 = e_1(S)/mu unidad Y real"), ("Amu", "(P5) B = A + mu A"),
                ("th", "(P6) B^- = theta R"), ("RA", "     R = A (tau ya en A)"),
                ("vP", "     v(theta) = p^{v_p(N)}")):
    print("   %-34s %d de %d" % (txt, tot[kk], tot["n"]))
