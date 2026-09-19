# -*- coding: utf-8 -*-
r"""II_ruta_A_exacta.py -- la linea A de II comprobada EXACTAMENTE, fila a fila, con la ruta.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Para cada m = 0, 3 (mod 4), 3 <= m <= MMAX: signos base en los primos <= sqrt(m) (eps(2)=eps(3)=-1,
resto +1; si no cierra, se prueban variantes: ademas eps(5)=-1, eps(7)=-1, ...), equilibrado voraz de
los primos de (sqrt(m), m/2] con pesos p*T(floor(m/p)), y completacion EXACTA con los primos de
(m/2, m]: se busca un subconjunto A con sum A = (P + R)/2 (bitset).  Si la suma se alcanza, la eps
construida cumple S(m) = 0, y se RECOMPRUEBA S(m) = 0 recalculando la suma entera (control).
Uso: python II_ruta_A_exacta.py MMAX
"""
import sys, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
MMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 3


def criba(n):
    spf = list(range(n + 1))
    for i in range(2, int(n ** .5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


spf = criba(MMAX)
PR = [p for p in range(2, MMAX + 1) if spf[p] == p]
VARIANTES = [{2: -1, 3: -1}, {2: -1, 3: -1, 5: -1}, {2: -1, 3: -1, 7: -1}, {2: -1}, {3: -1},
             {2: -1, 3: -1, 5: -1, 7: -1}, {2: -1, 5: -1}, {3: -1, 5: -1}]


def intenta(m, base):
    raiz = math.isqrt(m)
    s = {p: base.get(p, 1) for p in PR if p <= raiz}
    eps = [0] * (m + 1); eps[1] = 1
    mayor = [1] * (m + 1)
    for c in range(2, m + 1):
        p = spf[c]
        eps[c] = (s[p] if p <= raiz else 1) * eps[c // p]
        mayor[c] = max(p, mayor[c // p])
    T = [0] * (m + 1)
    for k in range(1, m + 1):
        T[k] = T[k - 1] + k * eps[k]
    R = sum(c * eps[c] for c in range(1, m + 1) if mayor[c] <= raiz)
    medios = [p for p in PR if raiz < p <= m // 2]
    for p in sorted(medios, key=lambda p: -abs(p * T[m // p])):
        w = p * T[m // p]
        if abs(R - w) < abs(R + w):
            R -= w; s[p] = -1
        else:
            R += w; s[p] = 1
    L = [p for p in PR if m // 2 < p <= m]
    P = sum(L)
    if (P + R) % 2 or abs(R) > P:
        return None
    objetivo = (P + R) // 2
    # subconjunto de L con suma objetivo, con reconstruccion
    alcanz = [1]
    for p in L:
        alcanz.append(alcanz[-1] | (alcanz[-1] << p))
    if not (alcanz[-1] >> objetivo) & 1:
        return None
    t = objetivo
    for i in range(len(L) - 1, -1, -1):
        p = L[i]
        if (alcanz[i] >> t) & 1:
            s[p] = 1
        else:
            s[p] = -1; t -= p
    # control: recalcular S(m) con todos los signos
    e = [0] * (m + 1); e[1] = 1
    for c in range(2, m + 1):
        e[c] = s[spf[c]] * e[c // spf[c]]
    assert sum(c * e[c] for c in range(1, m + 1)) == 0, ("control falla", m)
    return True


filas = [m for m in range(max(3, MMIN), MMAX + 1) if m % 4 in (0, 3)]
sin, usos = [], {}
for m in filas:
    ok = None
    for i, b in enumerate(VARIANTES):
        if intenta(m, b):
            ok = i
            break
    if ok is None:
        sin.append(m)
    else:
        usos[ok] = usos.get(ok, 0) + 1
    if m % 1000 in (0, 999):
        print("hasta m=%d: filas %d, sin solucion por la ruta %d" % (m, sum(usos.values()) + len(sin), len(sin)))
        sys.stdout.flush()
print("filas admisibles 3..%d: %d; resueltas por la ruta con S(m)=0 recomprobado: %d; sin: %s"
      % (MMAX, len(filas), sum(usos.values()), sin[:30]))
print("variante usada -> filas:", {str(VARIANTES[k]): v for k, v in sorted(usos.items())})
