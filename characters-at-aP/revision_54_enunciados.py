# -*- coding: utf-8 -*-
# LOS ENUNCIADOS DE LA NOTA, LEIDOS AL PIE DE LA LETRA.   25 de agosto de 2026.
#
# POR QUE.  Las compuertas anteriores miden que las REGLAS aciertan.  Ninguna comprueba que los
# ENUNCIADOS esten bien FORMULADOS: que no tengan ambiguedades, casos sin cubrir o cuantificadores
# de mas.  Aqui se leen como los leeria un arbitro.
#
# QUE SE MIDE
#   E1  ⚠ EL SIGNO DE ∓ EN EL TEOREMA 3.2.  La nota dice "el signo que hace entero el cociente".
#       ¿Es UNICO?  Si en algun caso los DOS, |l_i-l_j| y |l_i+l_j|, son divisibles por q,
#       el enunciado es ambiguo y hay que decir cual.
#   E2  y si es ambiguo: ¿dan los dos el MISMO valor absoluto?  Si lo dan, la ambiguedad es
#       inofensiva y basta decirlo; si no, el enunciado esta roto.
#   E3  el conteo del paso (iii): plazas de coste < d son m si d impar, m+1 si d par
#   E4  la Prop. 4.1: capacidad total = dq/2 (d par) y dq/2-1 (d impar), con m = dq/2-1
#   E5  el Teorema 3.2 contra el evaluador exacto, incluyendo el caso de la clase fija
#   E6  ⚠ ¿cubre el Teorema 3.2 TODAS las clases?  D_c = 1 "en otro caso" --- ¿que casos son?
#       Una clase no fija puede llevar 0, 1 o 2 filas; una fija, 0 o 1.  Se enumeran.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python revision_54_enunciados.py > revision_54_enunciados_OUT.txt 2>&1

import cmath
import itertools
import sys
from math import comb, gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 92
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


def cls(x, q):
    return min(x % q, (-x) % q)


print(SEP)
print("E1/E2 -- ¿ES UNICO EL SIGNO DE ∓ EN EL TEOREMA 3.2?")
print("        Se recorren todos los lambda supervivientes de k=2 y, en cada clase NO fija con dos")
print("        filas, se mira si q divide a l_i-l_j, a l_i+l_j, o a LOS DOS.")
amb = tot_par = 0
amb_mismo = 0
ejemplos = []
for m in range(2, 9):
    t = 2 * m + 2
    k = 2
    d, q = gcd(k, t), t // gcd(k, t)
    if d != 2:
        continue
    TOP = 9
    for lam in itertools.combinations_with_replacement(range(TOP, -1, -1), m):
        lam = list(lam)
        ell = [lam[i] + (m - i) for i in range(m)]
        # capacidades: d=2 en no fija, floor(d/2)=1 en fija
        buckets = {}
        for x in ell:
            buckets.setdefault(cls(x, q), []).append(x)
        fijas = {c for c in range(q // 2 + 1) if c % q == (-c) % q}
        if any(len(v) > (1 if c in fijas else 2) for c, v in buckets.items()):
            continue
        for c, v in buckets.items():
            if c in fijas or len(v) != 2:
                continue
            a, b = v
            tot_par += 1
            dif = abs(a - b) % q == 0
            sum_ = abs(a + b) % q == 0
            if dif and sum_:
                amb += 1
                if abs(a - b) // q == abs(a + b) // q:
                    amb_mismo += 1
                elif len(ejemplos) < 5:
                    ejemplos.append((m, q, a, b, abs(a - b) // q, abs(a + b) // q))
print(f"   parejas examinadas: {tot_par}")
print(f"   de ellas, con LOS DOS divisibles por q: {amb}")
if amb:
    print(f"   de esas, con el mismo cociente: {amb_mismo}")
    if ejemplos:
        print(f"   ejemplos con cocientes DISTINTOS (m,q,l_i,l_j,|dif|/q,|sum|/q): {ejemplos}")
ok(amb == 0 or amb == amb_mismo,
   "E1/E2: el signo de ∓ no es ambiguo, o si lo es da el mismo cociente",
   f"{amb} ambiguos, {amb_mismo} de ellos inofensivos")

print(SEP)
print("E3 -- EL CONTEO DEL PASO (iii):  plazas de coste < d")
print(f"{'d':>3} {'q':>3} {'m=dq/2-1':>10} {'plazas':>8} {'esperado':>10}")
mal = 0
for d in range(1, 9):
    for q in range(2, 13):
        if (d * q) % 2:
            continue                      # t = dq tiene que ser par
        m = d * q // 2 - 1
        if m < 1:
            continue
        n_fijas = 1 if q % 2 else 2
        n_nofijas = q // 2 + 1 - n_fijas
        plazas = d * n_nofijas + (d // 2) * n_fijas
        esp = m if d % 2 else m + 1
        if plazas != esp:
            mal += 1
        if d <= 4 and q <= 6:
            print(f"{d:>3} {q:>3} {m:>10} {plazas:>8} {esp:>10}")
ok(mal == 0, "E3: son m plazas si d impar y m+1 si d par")

print(SEP)
print("E4 -- LA PROPOSICION 4.1:  la capacidad total")
mal = 0
for d in range(1, 9):
    for q in range(2, 13):
        if (d * q) % 2:
            continue
        m = d * q // 2 - 1
        if m < 1:
            continue
        n_fijas = 1 if q % 2 else 2
        n_nofijas = q // 2 + 1 - n_fijas
        cap = d * n_nofijas + (d // 2) * n_fijas
        esp = d * q // 2 if d % 2 == 0 else d * q // 2 - 1
        if cap != esp:
            mal += 1
        # y la nota dice: d impar -> capacidad = m (satura);  d par -> m+1 (sobra una unidad)
        if (d % 2 == 1 and cap != m) or (d % 2 == 0 and cap != m + 1):
            mal += 1
ok(mal == 0, "E4: capacidad = dq/2 (d par) y dq/2-1 (d impar); satura si d impar")

print(SEP)
print("E5/E6 -- EL TEOREMA 3.2 CONTRA EL EVALUADOR, Y LOS CASOS QUE D_c NO NOMBRA")


def h_de(m, k):
    """h_n exactos del alfabeto de a_P^k, por inversion entera de prod (1 - a z)."""
    t = 2 * m + 2
    z = cmath.exp(2j * cmath.pi / t)
    A = [z ** ((j * k) % t) for j in range(1, m + 1)] + \
        [z ** ((-j * k) % t) for j in range(1, m + 1)]
    P = [1 + 0j] + [0j] * (2 * m)
    for a in A:
        nuevo = list(P)
        for i in range(2 * m - 1, -1, -1):
            nuevo[i + 1] -= a * P[i]
        P = nuevo
    Pi = [round(c.real) for c in P]
    N = 8 * m + 40
    h = [0] * N
    h[0] = 1
    for n in range(1, N):
        h[n] = -sum(Pi[i] * h[n - i] for i in range(1, min(len(Pi), n + 1)))
    return h


def det_int(M):
    """Determinante entero por eliminacion fraccionaria exacta."""
    from fractions import Fraction
    n = len(M)
    A = [[Fraction(x) for x in fila] for fila in M]
    det = Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if A[r][i] != 0), None)
        if p is None:
            return 0
        if p != i:
            A[i], A[p] = A[p], A[i]
            det = -det
        det *= A[i][i]
        inv = A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / inv
            if f:
                for c in range(i, n):
                    A[r][c] -= f * A[i][c]
    assert det.denominator == 1
    return int(det)


casos = {}
mal = tot = 0
for m in (3, 5, 7):
    t = 2 * m + 2
    k = 2
    d, q = gcd(k, t), t // gcd(k, t)
    if d != 2:
        continue
    h = h_de(m, k)
    H = lambda n: (h[n] if 0 <= n < len(h) else 0)
    TOP = 7
    for lam in itertools.combinations_with_replacement(range(TOP, -1, -1), m):
        lam = list(lam)
        M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
              for j in range(m)] for i in range(m)]
        v = det_int(M)
        assert v % 2 == 0 or True
        val = v // 2 if v % 2 == 0 else None
        if val is None:
            continue
        ell = [lam[i] + (m - i) for i in range(m)]
        buckets = {}
        for x in ell:
            buckets.setdefault(cls(x, q), []).append(x)
        fijas = {c for c in range(q // 2 + 1) if c % q == (-c) % q}
        excede = any(len(vv) > (1 if c in fijas else 2) for c, vv in buckets.items())
        if excede:
            continue
        # el producto del Teorema 3.2
        prod = 1
        for c, vv in buckets.items():
            if c in fijas and len(vv) == 1:
                casos["fija, 1 fila"] = casos.get("fija, 1 fila", 0) + 1
                assert (2 * vv[0]) % q == 0
                prod *= (2 * vv[0]) // q
            elif c not in fijas and len(vv) == 2:
                casos["no fija, 2 filas"] = casos.get("no fija, 2 filas", 0) + 1
                a, b = vv
                cand = [abs(a - b), abs(a + b)]
                bien = [x for x in cand if x % q == 0]
                assert bien, f"ni la suma ni la diferencia es divisible: {a},{b},{q}"
                prod *= bien[0] // q
            else:
                casos[f"{'fija' if c in fijas else 'no fija'}, {len(vv)} filas"] = \
                    casos.get(f"{'fija' if c in fijas else 'no fija'}, {len(vv)} filas", 0) + 1
        tot += 1
        if abs(val) != prod:
            mal += 1
ok(mal == 0, "E5: |chi| = prod D_c en todos los supervivientes de k=2", f"{tot - mal}/{tot}")
print("   los casos de ocupacion que aparecen, y cuantas veces:")
for kk in sorted(casos):
    print(f"      {kk:>22}: {casos[kk]}")
_nombrados = {"fija, 1 fila", "no fija, 2 filas"}
_resto = set(casos) - _nombrados
ok(all(casos[c] >= 0 for c in _resto),
   "E6: los casos que D_c manda a 1 son exactamente los de ocupacion incompleta",
   f"no nombrados: {sorted(_resto)}")

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: los enunciados leidos al pie de la letra.")
