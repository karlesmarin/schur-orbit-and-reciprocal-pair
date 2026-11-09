# encaje: comprobación propia (8-nov) del párrafo «Prior art for the key step» de la v2: [2m+1]_q = (-1)^h [r]_q y que
# el único [r]_q racional con 1 < r < q−1 es [3]_6 = 2 (racional ⇔ fijo por todos sus conjugados de Galois)
import sys
from math import gcd

from mpmath import mp, mpf, sin, pi, fabs

sys.stdout.reconfigure(encoding='utf-8')
mp.dps = 40


def br(r, q, p=1):
    """[r]_q y su conjugado σ_p: sin(rpπ/q)/sin(pπ/q), p primo con 2q."""
    return sin(r * p * pi / q) / sin(p * pi / q)


malos = 0
for q in range(2, 61):
    for m in range(1, 80):
        h, r = divmod(2 * m + 1, q)
        lhs = br(2 * m + 1, q)
        rhs = (-1) ** h * (br(r, q) if r else mpf(0))
        if fabs(lhs - rhs) > mpf(10) ** -30:
            malos += 1
print('identidad [2m+1]_q = (-1)^h [r]_q: discrepancias', malos, '(q ≤ 60, m < 80)')
racionales = []
for q in range(4, 121):
    for r in range(2, q - 1):
        v = br(r, q)
        if all(fabs(br(r, q, p) - v) < mpf(10) ** -25 for p in range(1, 2 * q) if gcd(p, 2 * q) == 1):
            racionales.append((q, r, mp.nstr(v, 6)))
print('[r]_q racionales con 1 < r < q-1, q ≤ 120:', racionales)
# control: [1]_q y [q-1]_q valen 1 (racionales) y deben salir si se incluyen
ctrl = [(q, r) for q in (5, 7) for r in (1, q - 1)
        if all(fabs(br(r, q, p) - br(r, q)) < mpf(10) ** -25 for p in range(1, 2 * q) if gcd(p, 2 * q) == 1)]
print('control ([1]_q y [q-1]_q deben salir racionales):', ctrl)
