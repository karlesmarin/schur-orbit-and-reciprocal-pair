# -*- coding: utf-8 -*-
# encaje: comprobaciones EXACTAS de las matemáticas que propone el arbitraje V253 antes de escribirlas en la v2
#   K1  gcd(q, r) = 1 => N_q por raíces = N_q por corraíces (tipo C, r = 2: todo q impar), en ρ y en λ+ρ; control q par
#   K2  Sp(4), q = 4, λ = (1,0): χ(g_{1/4}) exacto en Z[ζ_4] y Z_4(ℓ)/Z_4(ρ)
#   K3  G2, q = 3: valor del recuadro de §12 contra Fulton–Harris (Prop. 1 de [Kar24]) en el centro ωI del SL3 largo
#   K4  alfabeto: (2n+1)·μ_4 − {1} − {−1} = 2n·μ_4 + {i, −i} (Kumari, Thm 2.8 con t = 4, x_j = 1, y = i)
# Funciones copiadas de anc/certificar_exacto.py (lmul, ladd, pairings, prod_uno_menos, dominantes, ldiv_exacta, Ciclo).
# Ejecutar: PYTHONIOENCODING=utf-8 python comprobar_v253.py > comprobar_v253_OUT.txt
import sys
from collections import Counter
from fractions import Fraction

from sympy import Poly, cyclotomic_poly, symbols, Rational, sqrt, I, simplify, expand

sys.stdout.reconfigure(encoding='utf-8')
FALLOS = []


def ok(c, que, dato=''):
    print('   [%s] %s %s' % ('OK' if c else 'FALLA', que, dato))
    if not c:
        FALLOS.append(que)


def lmul(a, b):
    r = Counter()
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            r[e1 + e2] += c1 * c2
    return {e: c for e, c in r.items() if c}


def ladd(a, b, s=1):
    r = Counter(a)
    for e, c in b.items():
        r[e] += s * c
    return {e: c for e, c in r.items() if c}


def pairings(x, coraiz=False):
    m = len(x)
    out = []
    for i in range(m):
        out.append(x[i] if coraiz else 2 * x[i])
        for j in range(i + 1, m):
            out += [x[i] - x[j], x[i] + x[j]]
    return out


def prod_uno_menos(exps):
    r = {0: 1}
    for a in exps:
        r = lmul(r, {0: 1, a: -1})
    return r


def dominantes(m, top):
    if m == 0:
        yield ()
        return
    for v in range(top, -1, -1):
        for r in dominantes(m - 1, v):
            yield (v,) + r


def ldiv_exacta(num, den):
    num = dict(num)
    dmax = max(den)
    lc = den[dmax]
    q = {}
    while num:
        e = max(num)
        c = num[e]
        if c % lc:
            raise ValueError('no exacta')
        t = c // lc
        s = e - dmax
        q[s] = t
        num = ladd(num, lmul({s: t}, den), -1)
        if num and max(num) < min(den) + min(q) - 10 ** 4:
            raise ValueError('no exacta')
    return q


X = symbols('X')


class Ciclo:
    def __init__(self, q):
        self.q = q
        self.phi = [int(c) for c in Poly(cyclotomic_poly(q, X), X).all_coeffs()[::-1]]
        self.n = len(self.phi) - 1

    def red(self, coefs):
        c = list(coefs)
        n = self.n
        for k in range(len(c) - 1, n - 1, -1):
            a = c[k]
            if a:
                for i in range(n + 1):
                    c[k - n + i] -= a * self.phi[i]
        c = c[:n] + [0] * max(0, n - len(c))
        return tuple(c[:n])

    def de_laurent(self, d):
        c = [0] * self.q
        for e, v in d.items():
            c[e % self.q] += v
        return self.red(c)


def chi_laurent(lam):
    m = len(lam)
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    lr = sum(lam[i] * rho[i] for i in range(m))
    num = lmul({-lr: 1}, prod_uno_menos(pairings(ell)))
    return ldiv_exacta(num, prod_uno_menos(pairings(rho)))


def N(q, x, coraiz=False):
    return sum(1 for a in pairings(x, coraiz) if a % q == 0)


def Z(q, x):
    p = 1
    for a in pairings(x):
        if a % q == 0:
            p *= a
    return p


# ---------------------------------------------------------------- K1
print('=' * 100)
print('K1  tipo C_m (r = 2): N_q por raíces (x, α) frente a N_q por corraíces <x, α^∨>, en x = ρ y x = λ+ρ')
for paridad in ('impar', 'par'):
    casos = iguales = 0
    for m in range(2, 6):
        rho = [m - i for i in range(m)]
        for lam in dominantes(m, 6 if m < 4 else 4):
            ell = [lam[i] + rho[i] for i in range(m)]
            for q in range(2, 16):
                if (q % 2 == 1) != (paridad == 'impar'):
                    continue
                for x in (rho, ell):
                    # término a término: el mismo conjunto de raíces positivas
                    raiz = [a % q == 0 for a in pairings(x)]
                    cor = [a % q == 0 for a in pairings(x, True)]
                    casos += 1
                    iguales += (raiz == cor)
    if paridad == 'impar':
        ok(iguales == casos, 'K1: q impar (gcd(q,2)=1): las dos selecciones coinciden raíz a raíz',
           '(%d de %d pares (x, q), m = 2..5, q = 3..15)' % (iguales, casos))
    else:
        ok(iguales < casos, 'K1 control: q par: las selecciones deben diferir en algún caso',
           '(coinciden en %d de %d)' % (iguales, casos))
# y el criterio: con q impar, N_q(ℓ) = N_q(ρ) por raíces iff por corraíces
casos = acuerdo = 0
for m in range(2, 6):
    rho = [m - i for i in range(m)]
    for lam in dominantes(m, 6 if m < 4 else 4):
        ell = [lam[i] + rho[i] for i in range(m)]
        for q in range(3, 16, 2):
            casos += 1
            acuerdo += ((N(q, ell) == N(q, rho)) == (N(q, ell, True) == N(q, rho, True)))
ok(acuerdo == casos, 'K1: q impar: las dos condiciones de supervivencia coinciden', '(%d de %d)' % (acuerdo, casos))

# ---------------------------------------------------------------- K2
print('=' * 100)
print('K2  Sp(4), q = 4, λ = (1,0): valor exacto del carácter y cociente Z_4(ℓ)/Z_4(ρ)')
C4 = Ciclo(4)
chi = C4.de_laurent(chi_laurent((1, 0)))
ell, rho = [3, 1], [2, 1]
cociente = Fraction(Z(4, ell), Z(4, rho))
print('   χ_(1,0)(g_{1/4}) en Z[i] (coeficientes en 1, i):', chi, '; Z_4(ℓ) =', Z(4, ell), ', Z_4(ρ) =', Z(4, rho),
      ', N_4(ℓ) =', N(4, ell), ', N_4(ρ) =', N(4, rho))
ok(chi == (-2, 0), 'K2: χ_(1,0)(g_{1/4}) = −2, |χ| = 2', '(valor %s)' % (chi,))
ok(cociente == 1, 'K2: Z_4(ℓ)/Z_4(ρ) = 1', '(valor %s)' % cociente)
ok(4 * 1 % 4 == 0 and (2 * 2) % 4 == 0, 'K2: q = 4 divide 2m = 4 (está en la lista del Teorema 5.1)', '')
ok((2 * 2 + 2) % 4 != 0, 'K2: q = 4 no divide 2m+2 = 6 (no es potencia de a_P)', '')

# ---------------------------------------------------------------- K3
print('=' * 100)
print('K3  G2, q = 3: Fulton–Harris (Prop. 24.48; Prop. 1 de [Kar24]) en ωI del SL3 largo, contra el recuadro de §12')
w = Rational(-1, 2) + sqrt(3) * I / 2


def dimGL3(a, b, c):
    return Fraction((a - b + 1) * (b - c + 1) * (a - c + 2), 2)


def S_central(a, b, c):
    return w ** (a + b + c) * int(dimGL3(a, b, c))


def caja(a, b):
    v = Fraction((b + 1) * (a + b + 2) * (a + 2 * b + 3), 6)
    return v if a % 3 == 0 else (-v if a % 3 == 1 else 0)


def dimG2(a, b):
    # Weyl: con u = a+1, v = 3(b+1) los seis valores u, v, u+v, 2u+v, 3u+v, 3u+2v sobre 1,3,4,5,6,9
    u, v = a + 1, 3 * (b + 1)
    num = u * v * (u + v) * (2 * u + v) * (3 * u + v) * (3 * u + 2 * v)
    return Fraction(num, 1 * 3 * 4 * 5 * 6 * 9)


casos = buenos = ctrl = 0
for k in range(0, 7):
    for l in range(0, 7):
        val = (S_central(k + 2 * l + 1, k + l + 1, 0) - S_central(k + 2 * l + 1, l, 0)) / \
              (S_central(1, 1, 0) - S_central(1, 0, 0))
        val = simplify(expand(val))
        casos += 1
        buenos += (simplify(val - caja(k, l)) == 0)
        # control: el recuadro con los signos de las ramas a ≡ 0 y a ≡ 1 (mod 3) intercambiados
        ctrl += (simplify(val + caja(k, l)) == 0 and caja(k, l) != 0)
ok(buenos == casos, 'K3: Fulton–Harris en ωI = valor del recuadro de §12 (a = k, b = l)', '(%d de %d, a, b ≤ 6)' % (buenos, casos))
ok(ctrl == 0, 'K3 control: el recuadro con los signos intercambiados no coincide en ningún peso no nulo',
   '(coincide en %d)' % ctrl)

# ---------------------------------------------------------------- K4
print('=' * 100)
print('K4  (2n+1)·μ_4 − {1} − {−1} = 2n·μ_4 + {i, −i} como multiconjuntos de Z/4 (exponentes de i)')
bien = True
for n in range(1, 8):
    a = Counter({r: 2 * n + 1 for r in range(4)})
    a[0] -= 1
    a[2] -= 1
    b = Counter({r: 2 * n for r in range(4)})
    b[1] += 1
    b[3] += 1
    bien &= (a == b)
    # y es el alfabeto de a_P^k con d = 2n+1, q = 4: t = 4(2n+1), m = t/2 − 1, k = d
    t = 4 * (2 * n + 1)
    m = t // 2 - 1
    A = Counter()
    for j in range(1, m + 1):
        e = (j * (2 * n + 1)) % t       # exponente de ζ_t; ζ_t^e con e múltiplo de d = 2n+1 → i^{e/d}
        A[(e // (2 * n + 1)) % 4] += 1
        A[(-e // (2 * n + 1)) % 4] += 1
    bien &= (A == a)
ok(bien, 'K4: identidad de multiconjuntos y alfabeto de a_P^d con q = 4 (n = 1..7)', '')

print('=' * 100)
print('RESULTADO:', 'todas las comprobaciones pasan y sus controles fallan.' if not FALLOS else 'FALLAN: %s' % FALLOS)
