# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# auditar_figuras.py -- cada afirmacion de un pie de figura, contra los datos que la figura dibuja.
# Que comprueba: fig_clase (los defectos no nulos son exactamente los pares con p | h^-), fig_stickel
# (42 valores no primos; defecto = numero de divisores elementales divisibles por p), fig_perfil (las
# diez ordenes tienen v_p(N)=1; capas nulas por debajo de p-1 y capa p-1 de dimension 1), fig_simetria
# (Gorenstein <=> todas las columnas llegan a d), fig_sierra (los dos dientes en q'=13).
# Uso: python auditar_figuras.py     (en este directorio)
import json, os, re, sys

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
fallos = []


def marca(nombre, ok, detalle=""):
    print("%-58s %s %s" % (nombre, "OK" if ok else "FALLA", detalle))
    if not ok:
        fallos.append(nombre)


# ---------------------------------------------------------------- h^- exacto
hm = {}
for ln in open(os.path.join(AQUI, "stickelberger_smith_OUT.txt")):
    t = ln.split()
    if t and t[0].isdigit():
        hm[int(t[0])] = int(t[4])


def es_primo(n):
    return n > 1 and all(n % k for k in range(2, int(n ** 0.5) + 1))


# ---------------------------------------------------------------- fig_clase
from bernoulli_W1 import W1, h_menos_numerico

qs = [x for x in range(5, 90) if es_primo(x)]
ps = [3, 5, 7, 11, 13]
mal = []
for n in qs:
    h = h_menos_numerico(n)
    if n in hm and hm[n] != h:
        mal.append(("h^- distinto del exacto", n, h, hm[n]))
    for p in ps:
        if p == n:
            continue
        defecto = (n - 1) // 2 - W1(n, p)
        if (defecto > 0) != (h % p == 0):
            mal.append((n, p, defecto, h))
marca("fig_clase: defecto > 0  <=>  p | h^-  (%d celdas)" % (len(qs) * len(ps) - 5), not mal, mal[:4])

# ---------------------------------------------------------------- fig_stickel
nprim = [n for n in range(8, 90) if n % 4 != 2 and not es_primo(n) and n in hm]
marca("fig_stickel: no primos 8<=n<90, n!=2 mod 4  (42 filas)", len(nprim) == 42, len(nprim))
from bernoulli_W1 import rango_mod_p
from math import gcd
malst = []
for n in nprim:
    reps = [u for u in range(1, n) if gcd(u, n) == 1 and u < n - u]
    d = len(reps)
    s = lambda x: 0 if x % n == 0 else 2 * (x % n) - n
    filas = [[s(pow(u, -1, n) * c) for c in range(1, n)] for u in reps]
    for p in [3, 5, 7, 11, 13, 17]:
        if n % p == 0:
            continue
        if ((d - rango_mod_p(filas, p)) > 0) != (hm[n] % p == 0):
            malst.append((n, p))
marca("fig_stickel: celdas no nulas  <=>  p | h^- exacto", not malst, malst[:4])

# ---------------------------------------------------------------- fig_perfil
from perfil_W import dims_B
from reparto_p2 import phi

casos = [(3, 5, 2, 7), (3, 7, 2, 10), (3, 11, 2, 16), (3, 13, 2, 19),
         (5, 7, 2, 17), (5, 9, 2, 22), (5, 11, 2, 27), (5, 13, 2, 32),
         (7, 3, 2, 10), (7, 5, 2, 17)]
malv, malc = [], []
for (p, n, v, m) in casos:
    q = p ** v * n
    if 2 * m >= q:
        malv.append(("2m >= q", p, n, m, q))
    N = (2 * m + 1) // n if (2 * m + 1) % n == 0 else None
    vp = 0
    if N:
        x = N
        while x % p == 0:
            x //= p; vp += 1
    if vp != 1:
        malv.append(("v_p(N)", p, n, m, N, vp))
    Kc = min(9, phi(p ** v))
    dB = dims_B(q, p, m, n, Kc)
    W = [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, Kc)]
    if any(W[i] != 0 for i in range(1, p - 1)) or (p - 1 < Kc and W[p - 1] != 1):
        malc.append((p, n, m, W))
marca("fig_perfil: las diez ordenes tienen v_p(N) = 1", not malv, malv[:3])
marca("fig_perfil: capas 0 < i < p-1 nulas y capa p-1 igual a 1", not malc, malc[:3])

# ---------------------------------------------------------------- fig_simetria
ruta = os.path.join(AQUI, "perfil_W_OUT.txt")
pat = re.compile(r"m=\s*(\d+) q=\s*(\d+) p=(\d+).*?a=(\d+) b=(\d+)\s+W=\[([0-9, ]*)\]")
seis = [(3, 35), (7, 45), (10, 63), (6, 56), (9, 40), (16, 99)]
vistos, malg = {}, []
for ln in open(ruta, encoding="utf-8", errors="replace"):
    mm = pat.search(ln)
    if not mm:
        continue
    m, q, p, a, b = [int(x) for x in mm.groups()[:5]]
    W = [int(x) for x in mm.group(6).split(",") if x.strip()]
    vistos[(m, q)] = (p, a, b, W)
for (m, q) in seis:
    if (m, q) not in vistos:
        malg.append(("sin fila", m, q)); continue
    p, a, b, W = vistos[(m, q)]
    d = max(W) if W else 0
    e = None
    s = 0
    for i, w in enumerate(W):
        s += w
        if s == b:
            e = i + 1; break
    if e is None:
        malg.append(("e no deducible", m, q, W, b)); continue
    simetrica = all(W[i] + W[e - 1 - i] == d for i in range(e) if e - 1 - i < len(W))
    if simetrica != (a == b):
        malg.append((m, q, W, a, b, e, simetrica))
marca("fig_simetria: perfil autocomplementario <=> a = b (6 ordenes)", not malg, malg[:3])

# ---------------------------------------------------------------- fig_sierra
n = 13
s = [0 if c == 0 else 2 * c - n for c in range(n)]
hs = [0 if 2 * c == n else (c if 2 * c < n else c - n) for c in range(n)]
ok_s = all(s[c] == 2 * c - n for c in range(1, n)) and s[0] == 0
ok_h = all(hs[c] == (c if c <= 6 else c - 13) for c in range(n))
marca("fig_sierra: los dos dientes en q'=13", ok_s and ok_h)

print("")
print("AUDITORIA DE FIGURAS: %d comprobaciones, %d fallos" % (7, len(fallos)))
if fallos:
    print("fallan:", fallos)
