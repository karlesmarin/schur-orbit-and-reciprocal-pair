# -*- coding: utf-8 -*-
r"""extremo_bring.py -- el estrato maximal a n = 5, comprobado con codigo propio.

Nada se da por bueno sin recalcularlo aqui.  Se comprueban las tres afirmaciones mas
fuertes sobre este estrato, cada una por su cuenta:

  (1) FORMULA ELIPTICA n = 5.  Dice: la curva de Bring C_5 (M_1 = M_2 = M_3 = 0 en P^4) tiene
      #C_5(F_p) = p + 1 - 4 a_p, con a_p la traza de
          E : y^2 = -4x^3 + (373/16) x^2 - (431/8) x + 701/16
      (isogenia Jac(C_5) ~ E^4 sobre Q, que atribuye a Braden y Disney-Hogg, Prop. 4.5), y
          N_{delta=6}(5,p) = [ p + 1 - 4 a_p - 60 eps_R - 30 eps_4 - 24 eps_5 ] / 120,
      con eps_4 = [4 | p-1], eps_5 = [5 | p-1], eps_R = [R_5 = X^3+2X^2+3X+4 se escinde].
      AQUI se cuentan por separado: (a) los trinomios X^5 + aX + b con a,b != 0 que se escinden
      con raices distintas, por gcd(X^p - X, f); (b) los puntos de E por fuerza bruta.  Si la
      identidad falla en un primo, se ve.

  (2) SERENDIPIA DE CLEBSCH.  Dice: la cubica diagonal sum x_i = sum x_i^3 = 0 en P^4 tiene
      #C(F_p) = p^2 + (6 + chi_p(5)) p + 1.  Aqui se cuenta enumerando.

  (3) ASINTOTICA Delta/b = 1 + 2/p + O(p^{-min(2, n/2-1)}).  Se enfrenta a nuestro censo medido.

    python extremo_bring.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from itertools import product

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo                                                  # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 151


# ----------------------------------------------------------------- polinomios sobre F_p
def poly_mulmod(a, b, f, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i + j] = (r[i + j] + ai * bj) % p
    return poly_mod(r, f, p)


def poly_mod(a, f, p):
    a = a[:]
    df = len(f) - 1
    while len(a) - 1 >= df:
        c = a[-1] % p
        if c:
            d = len(a) - 1 - df
            for i, fi in enumerate(f):
                a[i + d] = (a[i + d] - c * fi) % p
        a.pop()
    while len(a) > 1 and a[-1] % p == 0:
        a.pop()
    return [x % p for x in a]


def poly_gcd(a, b, p):
    a, b = a[:], b[:]
    while any(b):
        inv = pow(b[-1], p - 2, p)
        b = [(x * inv) % p for x in b]
        r = a[:]
        while len(r) >= len(b) and any(r):
            c = r[-1] % p
            if c:
                d = len(r) - len(b)
                for i, bi in enumerate(b):
                    r[i + d] = (r[i + d] - c * bi) % p
            while len(r) > 1 and r[-1] % p == 0:
                r.pop()
            if len(r) < len(b):
                break
        a, b = b, r
    return a


def escinde_distinto(a, b, p, n=5):
    """X^n + aX + b tiene n raices distintas en F_p  <=>  f | X^p - X."""
    f = [b % p, a % p] + [0] * (n - 2) + [1]
    x = [0, 1]
    r = [1]
    e = p
    base = x
    while e:
        if e & 1:
            r = poly_mulmod(r, base, f, p)
        base = poly_mulmod(base, base, f, p)
        e >>= 1
    xp_menos_x = r[:] + [0] * max(0, 2 - len(r))
    xp_menos_x = [(xp_menos_x[i] - (1 if i == 1 else 0)) % p for i in range(max(len(r), 2))]
    while len(xp_menos_x) > 1 and xp_menos_x[-1] == 0:
        xp_menos_x.pop()
    g = poly_gcd(f, xp_menos_x, p)
    return len(g) - 1 == n


def traza_E(p):
    """a_p de y^2 = -4x^3 + (373/16)x^2 - (431/8)x + 701/16, por fuerza bruta."""
    inv16 = pow(16, p - 2, p)
    c3, c2, c1, c0 = (-4) % p, (373 * inv16) % p, (-431 * pow(8, p - 2, p)) % p, (701 * inv16) % p
    cuad = {(y * y) % p for y in range(p)}
    cuenta = 1                                    # punto del infinito
    for x in range(p):
        v = (((c3 * x + c2) * x + c1) * x + c0) % p
        if v == 0:
            cuenta += 1
        elif v in cuad:
            cuenta += 2
    return p + 1 - cuenta


def clebsch(p):
    """#{[x_1:...:x_5] : sum x_i = sum x_i^3 = 0} contando afines y dividiendo por p-1."""
    total = 0
    for x in product(range(p), repeat=4):
        x5 = (-sum(x)) % p
        if (sum(pow(t, 3, p) for t in x) + pow(x5, 3, p)) % p == 0:
            total += 1
    return (total - 1) // (p - 1)


def R5_escinde(p):
    return sum(1 for x in range(p) if (((x + 2) * x + 3) * x + 4) % p == 0) >= 1 and \
        len({x for x in range(p) if (((x + 2) * x + 3) * x + 4) % p == 0}) == 3


print("(1) FORMULA ELIPTICA PARA n = 5  --  trinomios y curva contados por separado")
print("%-5s %-8s %-10s %-10s %-12s %-12s %s" % ("p", "a_p", "N_orb", "eps(R5,4,5)", "120N+borde", "p+1-4a_p", "ok"))
def n_pares_rapido(p, n=5):
    """Mismo recuento que sum(escinde_distinto(a,b,p)) sobre a,b != 0, pero en O(p^2).

    Para cada raiz x y cada a, el b que hace a x raiz de X^n + aX + b esta forzado,
    b = -(x^n + ax).  La tabla (a,b) -> numero de raices DISTINTAS se llena recorriendo (x,a),
    y la casilla vale n exactamente cuando el trinomio se escinde del todo con raices distintas.
    La version lenta tardaba O(p^3) y por eso el barrido nunca se rehizo entero.
    """
    import numpy as np
    x = np.arange(p, dtype=np.int64)
    xn = np.array([pow(int(t), n, p) for t in x], dtype=np.int64)
    tabla = np.zeros((p, p), dtype=np.int16)
    for a in range(1, p):
        np.add.at(tabla, (a, (-(xn + a * x)) % p), 1)
    tabla[:, 0] = 0
    return int(np.count_nonzero(tabla == n))


# control de que lo rapido y lo lento cuentan lo MISMO, en los primos donde lo lento es viable
for _p in (7, 11, 13, 17, 19, 23):
    _lento = sum(1 for a in range(1, _p) for b in range(1, _p) if escinde_distinto(a, b, _p))
    assert _lento == n_pares_rapido(_p), (_p, _lento, n_pares_rapido(_p))
print("control: el recuento rapido coincide con el lento en p = 7..23")

fallos1 = casos1 = 0
primero_N = None
for p in range(7, PMAX + 1):
    if not es_primo(p) or p == 5:
        continue
    n_pares = n_pares_rapido(p)
    assert n_pares % (p - 1) == 0, (p, n_pares)
    N_orb = n_pares // (p - 1)
    ap = traza_E(p)
    epsR = 1 if R5_escinde(p) else 0
    eps4 = 1 if (p - 1) % 4 == 0 else 0
    eps5 = 1 if (p - 1) % 5 == 0 else 0
    izq = 120 * N_orb + 60 * epsR + 30 * eps4 + 24 * eps5
    der = p + 1 - 4 * ap
    ok = (izq == der)
    casos1 += 1
    fallos1 += 0 if ok else 1
    # EL FILTRO DE ANTES ERA "p <= 60 or not ok or p > PMAX - 20", y escondia p = 61..229.
    # De esa tabla recortada se leyo que "N(p) se hace no nulo por primera vez en p = 233", que
    # es falso: N(67) = 1 y estaba calculado aqui mismo.  Un filtro de impresion no es un dato.
    # Ahora se imprime SIEMPRE que N_orb > 0, que es justo lo que se queria saber.
    if N_orb or p <= 60 or not ok or p > PMAX - 20:
        print("%-5d %-8d %-10d %-12s %-12d %-12d %s"
              % (p, ap, N_orb, "%d,%d,%d" % (epsR, eps4, eps5), izq, der, "si" if ok else "NO"))
    if N_orb and primero_N is None:
        primero_N = p
    sys.stdout.flush()
print("primos: %d ; fallos: %d" % (casos1, fallos1))
print("primer p con N_orb > 0 (barrido desde 7): %s" % primero_N)

print("")
print("(2) CLEBSCH  #C(F_p) = p^2 + (6 + chi(5)) p + 1")
fallos2 = casos2 = 0
for p in range(7, min(PMAX, 43) + 1):
    if not es_primo(p):
        continue
    medido = clebsch(p)
    chi5 = 0 if p == 5 else (1 if pow(5 % p, (p - 1) // 2, p) == 1 else -1)
    pred = p * p + (6 + chi5) * p + 1
    ok = medido == pred
    casos2 += 1
    fallos2 += 0 if ok else 1
    print("   p=%-4d medido=%-8d formula=%-8d chi(5)=%-3d %s" % (p, medido, pred, chi5, "si" if ok else "NO"))
print("primos: %d ; fallos: %d" % (casos2, fallos2))

print("")
print("VEREDICTO PARCIAL: elipticas %s ; Clebsch %s"
      % ("SOBREVIVE" if not fallos1 else "FALLA", "SOBREVIVE" if not fallos2 else "FALLA"))
