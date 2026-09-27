"""Tres cifras que la nota cita y que ninguna salida archivada imprimia (auditoria del 27-sep-2026).

  (1) pie de fig_stickel: los 21 q' no primos, 8 <= q' < 90, q' != 2 mod 4, en los que la matriz
      cuadrada sobre unidades [s(u^{-1}c)] (u, c unidades) es singular sobre Q.  Rango exacto (sympy).
  (2) fig_zonas y su parrafo: el reparto de los 168 perfiles por bandas de W_1, que los 17 de la banda
      baja tienen p = 2, y cuantos de ellos tienen W_2 calculado (todos = d) y con que E los demas.
      Lee perfil_W_OUT.txt con el mismo filtro que figuras.py.
  (3) Proposicion sinnott: la formula entera, potencias de 2 y de n incluidas, en los 66 n.
      Lee stickelberger_smith_OUT.txt y compara Delta con 2^a h^- (2n)^d / w_n.
Solo lee ficheros de este directorio.
"""
import re
from math import gcd
from pathlib import Path

from sympy import factorint, isprime, totient

AQUI = Path(__file__).resolve().parent


def rango_Q(filas):
    """Rango exacto sobre Q por eliminacion de Bareiss (sin fracciones, enteros de Python)."""
    A = [list(f) for f in filas]
    r, previo = 0, 1
    for c in range(len(A[0])):
        piv = next((i for i in range(r, len(A)) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(r + 1, len(A)):
            A[i] = [(A[r][c] * A[i][j] - A[i][c] * A[r][j]) // previo for j in range(len(A[0]))]
        previo, r = A[r][c], r + 1
    return r

# (1)
ns = [n for n in range(8, 90) if n % 4 != 2 and not isprime(n)]
singulares = []
for n in ns:
    reps = [u for u in range(1, n) if gcd(u, n) == 1 and u < n - u]
    unidades = [c for c in range(1, n) if gcd(c, n) == 1]
    s = lambda x: 0 if x % n == 0 else 2 * (x % n) - n
    U = [[s(pow(u, -1, n) * c) for c in unidades] for u in reps]
    if rango_Q(U) < len(reps):
        singulares.append(n)
print("(1) q' no primos considerados: %d; matriz sobre unidades singular sobre Q en %d: %s"
      % (len(ns), len(singulares), singulares))

# (2)
pat = re.compile(r"m=\s*(\d+) q=\s*(\d+) p=(\d+) v=(\d+) q'=\s*(\d+) E=\s*(\d+).*?a=(\d+) b=(\d+)\s+W=\[([0-9, ]*)\]")
perfiles = []
for linea in (AQUI / "perfil_W_OUT.txt").read_text(encoding="utf-8").splitlines():
    mm = pat.search(linea)
    if not mm:
        continue
    m, q, p, v, n, E = (int(x) for x in mm.groups()[:6])
    W = [int(x) for x in mm.group(9).split(",") if x.strip()]
    if n in (1, 2, 3, 4, 6) or len(W) < 2:
        continue
    perfiles.append((p, E, int(totient(n)) // 2, W))
bandas = {"W_1=0": [], "0<W_1<=d/2": [], "d/2<W_1<d": [], "W_1=d": []}
for p, E, d, W in perfiles:
    w = W[1]
    k = "W_1=0" if w == 0 else ("0<W_1<=d/2" if 2 * w <= d else ("d/2<W_1<d" if w < d else "W_1=d"))
    bandas[k].append((p, E, d, W))
print("(2) perfiles: %d; %s" % (len(perfiles), ", ".join("%s: %d" % (k, len(v)) for k, v in bandas.items())))
baja = bandas["0<W_1<=d/2"]
con_w2 = [t for t in baja if len(t[3]) > 2]
print("    banda baja: primos %s; con W_2 calculado %d (todos W_2 = d: %s); sin W_2 %d, con E en %s"
      % (sorted({t[0] for t in baja}), len(con_w2), all(t[3][2] == t[2] for t in con_w2),
         len(baja) - len(con_w2), sorted({t[1] for t in baja if len(t[3]) <= 2})))

# (3)
bien = mal = 0
for linea in (AQUI / "stickelberger_smith_OUT.txt").read_text(encoding="utf-8").splitlines():
    mm = re.match(r"^(\d+) (\d+) (\d+) (\d+) (\d+) (\d+) Delta/h\^- = (\d+)", linea)
    if not mm:
        continue
    n, _, d, _, h, delta, _ = map(int, mm.groups())
    t = len(factorint(n))
    a = 0 if t == 1 else 2 ** (t - 2) - 1
    w = n if n % 2 == 0 else 2 * n
    num = 2 ** a * (2 * n) ** d * h
    if num % w == 0 and num // w == delta:
        bien += 1
    else:
        mal += 1
        print("    NO CUMPLE n=%d: Delta=%d, formula=%s" % (n, delta, num / w))
print("(3) formula de Sinnott completa: cumple en %d de %d n; fallos = %d" % (bien, bien + mal, mal))
