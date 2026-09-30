"""Pregunta 1 del articulo 1 en el marco de Prasad: rho de G como cocaracter del toro dual (plan del 27-sep).

Para G = Sp(2m), rho = (m,...,1) es un cocaracter del toro de G^v = SO(2m+1): rho(zeta_q) tiene valores propios
zeta^j, -m <= j <= m. Para G = SL(n+1), rho(zeta_q) en G^v = PGL(n+1) tiene valores propios zeta^j, 0 <= j <= n,
salvo escalar. Un elemento semisimple x de orden q es racional (conjugado de x^a para todo a primo con q) si:
  Sp(2m), SO(2m+1): el multiconjunto de valores propios queda fijo al multiplicar los exponentes por a;
  PGL(n+1): queda fijo salvo traslacion de los exponentes (escalar).
Comprueba, por fuerza bruta en q <= QMAX:
  (1) Sp(2m) en g_{1/q} y SO(2m+1) en rho(zeta_q) son racionales en los mismos q (lo fuerza
      Lambda^k(V+1) = Lambda^k V + Lambda^{k-1} V, con R(SO(2m+1)) = Z[Lambda^1..Lambda^m]);
  (2) tipo C (h = 2m): racional  <=>  q | h, h+1, h+2, o phi(q) <= 2;
  (3) tipo A (h = n+1): racional <=>  q | h-1, h, h+1, o phi(q) <= 2;
  (4) numeros regulares de W: tipo A, divisores de h-1 y de h; tipo B/C, divisores de h. Imprime los q
      racionales que no son regulares.
"""
from collections import Counter
from math import gcd

QMAX = 60
unidades = lambda q: [a for a in range(1, q) if gcd(a, q) == 1] or [1]
phi = lambda q: len(unidades(q))


def racional(exps, q, traslacion=False):
    base = Counter(e % q for e in exps)
    for a in unidades(q):
        img = Counter((a * e) % q for e in exps)
        if traslacion:
            if not any(Counter({(k + t) % q: v for k, v in img.items()}) == base for t in range(q)):
                return False
        elif img != base:
            return False
    return True


fallos = 0
print("tipo C (G^v = SO(2m+1)); regular = q | h")
for m in range(1, 13):
    h = 2 * m
    sp = [q for q in range(1, QMAX + 1) if racional([s * j for j in range(1, m + 1) for s in (1, -1)], q)]
    so = [q for q in range(1, QMAX + 1) if racional(range(-m, m + 1), q)]
    pred = [q for q in range(1, QMAX + 1) if h % q == 0 or (h + 1) % q == 0 or (h + 2) % q == 0 or phi(q) <= 2]
    if sp != so:
        fallos += 1
        print("  (1) FALLA m =", m, sp, so)
    if so != pred:
        fallos += 1
        print("  (2) FALLA m =", m, so, pred)
    print("  m=%2d h=%2d racionales %s | no regulares %s" % (m, h, so, [q for q in so if h % q]))

print("tipo A (G^v = PGL(n+1)); regular = q | h-1 o q | h")
for n in range(1, 16):
    h = n + 1
    pgl = [q for q in range(1, QMAX + 1) if racional(range(0, n + 1), q, traslacion=True)]
    pred = [q for q in range(1, QMAX + 1) if (h - 1) % q == 0 or h % q == 0 or (h + 1) % q == 0 or phi(q) <= 2]
    if pgl != pred:
        fallos += 1
        print("  (3) FALLA n =", n, pgl, pred)
    print("  n=%2d h=%2d racionales %s | no regulares %s"
          % (n, h, pgl, [q for q in pgl if (h - 1) % q and h % q]))

# controles: la prueba distingue (casos no racionales conocidos)
assert not racional(range(-3, 4), 11)
assert not racional(range(0, 4), 7, traslacion=True)
assert racional(range(0, 4), 5, traslacion=True)       # q = h+1 en A_3
print("fallos:", fallos)
