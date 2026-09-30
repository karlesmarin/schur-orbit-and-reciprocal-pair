"""The hypothesis p nmid h^- of the corollary on W_P (label cor:WPdos) cannot be dropped.

By the descent (label prop:descenso) and the local structure theorem (thm:loc(b)), W_P = dim_Fp Fp[(Z/q')^x] . S_0 with S_0(c) = M*chat; as p does not
divide M, it depends only on (q', p). This script computes that rank directly, as the rank mod p of the matrix
[S_0(u^{-1} c)] (u over (Z/q')^x/{+-1}, c over Z/q'), and compares it with the value d or d - d/t that
the corollary predicts. It also lists the pairs (q', p) behind the 51 cases of the descent grid that have
p | h^-. Reads descent_OUT_all.txt and stickelberger_smith_OUT.txt, next to this script or in ../../anc.
"""
import re
from math import gcd
from pathlib import Path

AQUI = Path(__file__).resolve().parent
ANC = AQUI if (AQUI / "descent_OUT_all.txt").exists() else AQUI.parents[1] / "anc"   # flat archive, or working tree


def rango_mod(filas, p):
    A = [[x % p for x in f] for f in filas]
    r = 0
    for c in range(len(A[0])):
        piv = next((i for i in range(r, len(A)) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], -1, p)
        A[r] = [x * inv % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        r += 1
    return r


def prediccion(n):
    x, pot = 1, []
    while True:
        pot.append(x)
        x = x * 2 % n
        if x == 1:
            break
    d = sum(1 for u in range(1, n) if gcd(u, n) == 1) // 2
    return d, (d if (n - 1) in pot else d - d // len(pot))


def W_P(n, p, M):
    chat = lambda c: (c % n) - n if (c % n) > n // 2 else (c % n)
    reps = [u for u in range(1, n) if gcd(u, n) == 1 and u < n - u]
    return rango_mod([[M * chat(pow(u, -1, n) * c) for c in range(n)] for u in reps], p)


n, p, M = 77, 5, 1
d, pred = prediccion(n)
print(f"(q',p,M)=({n},{p},{M}): rank W_P = {W_P(n, p, M)}, corollary without the h^- hypothesis would give {pred} (d={d})")
print(f"  e.g. m = (5*{n}-1)/2 = {(5 * n - 1) // 2}, q = 5^2*{n} = {25 * n}")
for n, p in ((23, 3), (31, 3), (7, 5), (13, 5)):   # two with p | h^-, two controls with p nmid h^-
    print(f"(q',p)=({n},{p}): rank {W_P(n, p, 1)}, predicted {prediccion(n)[1]}")

hm = {int(t[0]): int(t[4]) for t in (l.split() for l in (ANC / "stickelberger_smith_OUT.txt").read_text().splitlines())
      if len(t) > 4 and t[0].isdigit()}
pares = {}
for l in (ANC / "descent_OUT_all.txt").read_text().splitlines():
    r = re.search(r"p=\s*(\d+) a=\s*(\d+) n=\s*(\d+) M=\s*(\d+) v=\s*(\d+)", l)
    if r and "W_P=" in l:
        p, a, n, M, v = map(int, r.groups())
        if hm[n] % p == 0:
            pares.setdefault((n, p), set()).add((a, M, v))
print("grid cases with p | h^-:", {k: len(v) for k, v in pares.items()}, "total", sum(map(len, pares.values())))
