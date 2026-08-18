# -*- coding: utf-8 -*-
# ============================================================================================
#  VERIFICACION EXTERNA DE UNA REFUTACION.  14 de agosto de 2026.
#
#  POR QUE ESTE GUION EXISTE.  folding_t2.py REFUTA "la restriccion es +- una combinacion no
#  negativa de caracteres de Sp(2r)": aguanta 2708 de 2708 en r = 1 y CAE en r = 2, con
#  s_(8,8,7,7,6,6)(1,-1,z1^{+-1},z2^{+-1}) = -sp_(2,2) + sp_(1,1).
#
#  Una refutacion merece el mismo control externo que una confirmacion -- en este proyecto TODOS
#  mis errores de rango del mes inflaban el resultado, y aqui el sesgo iria al reves.  El riesgo
#  concreto es que el bug este en to_sp() o en sp_char(), no en la restriccion (que ya esta
#  contrastada contra la DEFINICION por SSYT y contra evaluacion numerica exacta).
#
#  ESTE GUION NO COMPARTE UNA SOLA LINEA CON folding_t2.py:
#    - la restriccion se calcula en LaurentPolynomialRing de Sage por el cociente de alternantes,
#      con det() de Sage (no mi Laplace por columnas de paridad, no mi division greedy);
#    - los caracteres de Sp(2r) los da WeylCharacterRing, o sea la formula de Freudenthal de Sage,
#      no mi determinante de Weyl.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: MSYS_NO_PATHCONV=1 docker run --rm -v "$(pwd -W):/w" -w /w sagemath/sagemath:latest \
#           sage folding_sage_check.sage
# ============================================================================================

def restriccion_sage(lam, r):
    N = 2 * r + 2
    R = LaurentPolynomialRing(QQ, ['z%d' % (i + 1) for i in range(r)])
    zs = R.gens()
    u = [R(1), R(-1)]
    for z in zs:
        u += [z, z ** -1]
    beta = [lam[i] + (N - 1 - i) for i in range(N)]
    F = R.fraction_field()
    M = matrix(F, N, N, lambda i, j: u[j] ** beta[i])
    D = matrix(F, N, N, lambda i, j: u[j] ** (N - 1 - i))
    q = M.det() / D.det()
    return R(q), R, zs


def sp_sage(mu, r, R, zs):
    """caracter irreducible de Sp(2r) por WeylCharacterRing (Freudenthal), como Laurent."""
    A = WeylCharacterRing(['C', r])
    ch = A(tuple(mu))
    p = R(0)
    for wt, m in ch.weight_multiplicities().items():
        try:
            v = [int(x) for x in wt.to_vector()]
        except AttributeError:
            v = [int(x) for x in list(wt)]
        t = R(m)
        for a in range(r):
            t *= zs[a] ** v[a]
        p += t
    return p


def dic(p, r):
    d = {}
    for e, c in zip(p.exponents(), p.coefficients()):
        k = tuple(int(x) for x in (e if r > 1 else [e]))
        d[k] = int(c)
    return {k: v for k, v in d.items() if v}


def dom(e):
    return tuple(sorted((abs(x) for x in e), reverse=True))


def expandir(p, r, R, zs):
    rem = dic(p, r)
    coef = {}
    guard = 0
    while rem:
        guard += 1
        if guard > 5000:
            raise RuntimeError("no termina")
        best = max((dom(e) for e in rem), key=lambda d: (sum(d), d))
        c = rem.get(best)
        if not c:
            raise RuntimeError("no invariante de Weyl: falta z^%r" % (best,))
        coef[best] = coef.get(best, 0) + c
        sub = dic(sp_sage(list(best), r, R, zs), r)
        for e, cc in sub.items():
            v = rem.get(e, 0) - c * cc
            if v:
                rem[e] = v
            elif e in rem:
                del rem[e]
    return {m: c for m, c in coef.items() if c}


CASOS = [
    (1, [0, 0, 0, 0],           {(0,): 1},                     "control: s_vacia = sp_0"),
    (1, [2, 0, 0, 0],           {(2,): 1, (0,): 1},            "el contraejemplo MINIMO del dia 13"),
    (1, [1, 1, 1, 0],           {(1,): -1},                    "signo negativo global, r = 1"),
    (1, [1, 1, 0, 0],           {},                            "un CERO"),
    (2, [2, 1, 0, 0, 0, 0],     {(2, 1): 1, (1, 0): 1},        "r = 2, signo positivo limpio"),
    (2, [8, 8, 8, 7, 7, 7],     {},                            "r = 2, un CERO"),
    (2, [8, 8, 7, 7, 6, 6],     {(2, 2): -1, (1, 1): 1},       "*** LA REFUTACION: signos MEZCLADOS ***"),
    (2, [8, 8, 7, 7, 5, 5],     {(3, 3): -1, (3, 1): -1, (2, 2): 1, (1, 1): -1},
                                                               "*** segunda mezclada ***"),
    (2, [8, 8, 7, 7, 6, 4],     {(4, 2): -1, (2, 2): -1, (1, 1): 1},
                                                               "*** tercera mezclada ***"),
]

print("=" * 108)
print("VERIFICACION EXTERNA (Sage) de la refutacion de la 'combinacion no negativa'")
print("=" * 108)
print("")
malos = 0
for (r, lam, esp, nota) in CASOS:
    p, R, zs = restriccion_sage(lam, r)
    got = expandir(p, r, R, zs) if p != 0 else {}
    ok = (got == esp)
    malos += (not ok)
    print("  r=%d  lambda=%-22s %s" % (r, str(tuple(lam)), nota))
    print("       Sage    : %s" % (sorted(got.items(), reverse=True) if got else "0"))
    print("       esperado: %s   -> %s" % (sorted(esp.items(), reverse=True) if esp else "0",
                                           "COINCIDE" if ok else "*** DISCREPA ***"))
    print("")
print("=" * 108)
print("  discrepancias: %d de %d" % (malos, len(CASOS)))
if malos == 0:
    print("  -> Sage CONFIRMA la refutacion por una via que no comparte codigo con folding_t2.py:")
    print("     la restriccion NO es +- una combinacion no negativa de caracteres de Sp(2r).")
    print("     Aguanta en r = 1 y cae en r = 2: el contraejemplo estaba UNA TALLA ARRIBA.")
else:
    print("  -> HAY DISCREPANCIA: antes de escribir nada hay que averiguar cual de las dos")
    print("     implementaciones esta mal.")
print("=" * 108)
