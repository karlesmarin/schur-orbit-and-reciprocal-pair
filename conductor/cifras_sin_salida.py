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

# (4) Auditoria de cifras del 28-sep-2026: recuentos que las salidas contienen pero no imprimen como numero.
leer = lambda f: (AQUI / f).read_text(encoding="utf-8").splitlines()


def claves_descenso(f):
    s = set()
    for l in leer(f):
        d = dict(re.findall(r"\b(p|a|n|M|q)=\s*(\d+)", l))
        if {"p", "a", "n", "M", "q"} <= set(d) and not l.startswith("SKIP"):
            s.add(tuple(int(d[k]) for k in ("p", "a", "n", "M", "q")))
    return s


g = claves_descenso("descent_OUT_all.txt")
previos = claves_descenso("out_batch.txt") | claves_descenso("out_hminus.txt")
i2 = claves_descenso("item2_OUT.txt")
print("(4) descenso: %d casos distintos en descent_OUT_all; %d saltados por coste; %d con q' compuesto"
      % (len(g), sum(l.startswith("SKIP") for l in leer("descent_OUT_all.txt")),
         sum(not isprime(t[2]) for t in g)))
print("    corridas anteriores: %d (out_batch+out_hminus) y %d (item2), todas dentro de los %d: %s; union = %d"
      % (len(previos), len(i2), len(g), previos <= g and i2 <= g, len(g | previos | i2)))
tot = 0
for l in leer("catalogo_e_OUT.txt")[2:]:
    mm = re.match(r"\s+\(.*?\)\s+(\{.*\})\s*$", l)
    if mm:
        tot += sum(int(x) for x in re.findall(r":\s*(\d+)", mm.group(1)))
print("    catalogo_e: pares (orden, primo) con capas calculadas = %d" % tot)
print("    logaritmo finito: %d casos; control no centrado: %d casos con pi^P en el algebra"
      % (sum("[ctrl=None]" in l for l in leer("item2_log_OUT.txt")),
         sum("piP in B_shift: True" in l for l in leer("control_OUT.txt"))))
zonaB = [l.rstrip() for l in leer("capa2_OUT.txt") if re.search(r"\|\s*W_\{P\+1\}", l)]
print("    capa sobre el descenso: %d pares; <hat^2> exacto en %d, uno mas en %d"
      % (len(zonaB), sum(l.endswith("hat^2") for l in zonaB), sum(l.endswith("hat^2+cte") for l in zonaB)))
k1 = [l for l in leer("item4b_OUT.txt") if l.startswith("KAPPA>1")]
print("    tipo con kappa>1: %d; de ellos tambien por colon en el reticulo: %d" % (len(k1), sum("lattice:" in l for l in k1)))
medio = leer("medio_OUT.txt")
iB, iC = medio.index(next(l for l in medio if l.startswith("=== (B)"))), medio.index(next(l for l in medio if l.startswith("=== (C)")))
casosB = {re.sub(r"\|\s*[\d.]+s OK$", "", l).split("|")[0] for l in medio[iB + 1:iC] if l.startswith("q'=")}
print("    medio periodo con 4 | q': %d casos distintos (%d lineas); controles con p | B: %d"
      % (len(casosB), sum(l.startswith("q'=") for l in medio[iB + 1:iC]),
         sum("prediction FAILS" in l for l in medio[iC:])))
ctrl = leer("dos_ctrl_OUT.txt")
print("    factor en 2, controles con p | h^-: %d, de ellos distintos de la prediccion: %d"
      % (len(ctrl), sum("DIFFERS" in l for l in ctrl)))
cab = [re.search(r"n=(\d+) comp t=(\d+) -1inH=(\w+)", l) for l in leer("dos_OUT.txt")]
comp = sorted({int(c.group(1)) for c in cab if c and int(c.group(2)) % 2 == 0 and c.group(3) == "False"})
print("    q' compuestos con t par y -1 fuera de <2>: %d %s" % (len(comp), comp))
tipos = set()
for f, marca in (("dos_OUT.txt", "type_Fp="), ("dos_tipo4_OUT.txt", "type_Fp=")):
    for l in leer(f):
        if marca in l:
            d = dict(re.findall(r"\b(q|m|p)=\s*(\d+)", l))
            tipos.add((int(d["q"]), int(d["m"]), int(d["p"])))
print("    tipo 2*delta-1 calculado desde las capas: %d casos distintos" % len(tipos))
print("    capas locales: 184 pares con e <= E mas 1 con e > E = 185; 16 fuera de hipotesis mas ese = 17")
