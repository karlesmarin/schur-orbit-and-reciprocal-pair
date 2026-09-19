# -*- coding: utf-8 -*-
# COMPUERTA: ?es la DUALIDAD POR COMPLEMENTO nuestra un caso de LEVEL-RANK?
#
# prop:ops(b) dice  A_{SU(n)}(q) = A_{SU(q-n)}(q):  igualdad de ANILLOS.
# Level-rank empareja SU(N)_K con SU(K)_N, y con q = N+K eso es justo n <-> q-n.  En el elemento
# principal de orden q los valores de caracteres SON las dimensiones cuanticas, y la identidad
# clasica de q-dimensiones bajo transposicion de diagramas es level-rank.  La pregunta que decide
# si nuestra proposicion esta contenida en esa identidad es:
#
#     ?coinciden los CONJUNTOS DE VALORES de caracteres, o solo los anillos que generan?
#
# Porque [q-m]_u = -u^{-m}[m]_u: la identidad de q-dimensiones vale SALVO raiz de la unidad, y
# multiplicar por una raiz de la unidad puede sacarte del subanillo.  Si los valores coinciden
# exactamente (o salvo signo), nuestra (b) es level-rank y hay que citarla.  Si solo coinciden los
# anillos, la nuestra es un enunciado distinto y mas fino.
#
# Se calcula chi_lambda por el BIALTERNANTE en Q(zeta_{2q}), sin formulas de producto.
#
# Authors: Carles Marin, Claude (AI assistant).
import sys


def valores(n, q, K):
    """multiconjunto de chi_lambda(g_q) para SU(n) en el elemento principal de orden q."""
    t = K.gen()                       # zeta_{2q}
    xs = [t ** ((n - 1 - 2 * j) % (2 * q)) for j in range(n)]
    niv = q - n
    if niv < 0:
        return None
    out = []
    # lambda = (l_1 >= ... >= l_{n-1} >= 0), l_1 <= nivel   (alcoba de SU(n) a nivel q-n)
    def recorre(pref, resto, tope):
        if resto == 0:
            lam = pref + [0]
            num = matrix(K, n, n, lambda a, b: xs[a] ** (lam[b] + n - 1 - b))
            out.append(num.det())
            return
        for v in range(tope, -1, -1):
            recorre(pref + [v], resto - 1, v)
    recorre([], n - 1, niv)
    den = matrix(K, n, n, lambda a, b: xs[a] ** (n - 1 - b)).det()
    return [c / den for c in out]


def norm(S, modo):
    if modo == "exacto":
        return set(S)
    if modo == "signo":
        return set([x if str(x) <= str(-x) else -x for x in S])
    return None


print("%-4s %-4s %-6s %-7s %-7s %-12s %-12s %s" % (
    "q", "n", "q-n", "#SU(n)", "#SU(q-n)", "mismos vals", "salvo signo", "veredicto"))
ok_ex = ok_sg = no = 0
for q in range(4, 14):
    K = CyclotomicField(2 * q)
    for n in range(2, q):
        m = q - n
        if m < 2 or n > 6 or m > 6:
            continue
        A = valores(n, q, K)
        B = valores(m, q, K)
        if A is None or B is None:
            continue
        ex = (norm(A, "exacto") == norm(B, "exacto"))
        sg = (norm(A, "signo") == norm(B, "signo"))
        if ex:
            ok_ex += 1
        if sg:
            ok_sg += 1
        if not sg:
            no += 1
        print("%-4d %-4d %-6d %-7d %-7d %-12s %-12s %s" % (
            q, n, m, len(A), len(B), "SI" if ex else "no", "SI" if sg else "no",
            "los valores coinciden" if sg else "SOLO el anillo"))
        sys.stdout.flush()
print()
print("pares comparados: %d | valores identicos: %d | identicos salvo signo: %d | ni eso: %d"
      % (ok_sg + no, ok_ex, ok_sg, no))
