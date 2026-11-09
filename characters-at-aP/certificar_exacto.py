# -*- coding: utf-8 -*-
# encaje: certificado EXACTO (8-nov-2026, Carles: «formaliza todos los teoremas y las certificaciones») de las tres
# comprobaciones que la nota hacía en coma flotante; enteros y Z[X]/Φ_q, sin tolerancias. No es herramienta de Excovery:
# es el certificado de una nota (va al suplemento, anc/).
#
#   C1  la especialización principal (2): bialternante = z^{-(λ,ρ)} ∏ (1-z^{(ℓ,α)})/(1-z^{(ρ,α)}) como IDENTIDAD de
#       polinomios de Laurent (antes: revision_47_verificacion.py R2, a 1e-12 en dos valores reales de z).
#       Control: el estadístico de CORRAÍCES (2e_i ↦ ℓ_i) debe fallar.
#   C2  Teorema 5.1 en rangos pequeños: χ_λ(g_{1/q}) calculado EXACTO en Z[ζ_q]; «todos enteros» contra la estabilidad
#       de Galois del espectro y contra la condición (I) (antes: galois.py (2), a 40 dígitos con mpmath).
#       Control: la condición (I) sin la cláusula q = 6 debe discrepar.
#   C3  Proposición 8.2: M_s^2 = ∏_{w} (2 - w - w^{-1}) EXACTO en Z[ζ_q], igual para todo s (antes:
#       revision_49_prueba_Ms.py S7, a 1e-8/3e-15). Control: el producto complejo sin módulos ∏(1-w) SÍ depende de s.
#
# Ejecutar: PYTHONIOENCODING=utf-8 python certificar_exacto.py > certificar_exacto_OUT.txt
import sys
from collections import Counter
from itertools import permutations
from math import gcd

from sympy import cyclotomic_poly, Poly, symbols

sys.stdout.reconfigure(encoding='utf-8')
FALLOS = []


def ok(c, que, dato=''):
    print('   [%s] %s %s' % ('OK' if c else 'FALLA', que, dato))
    if not c:
        FALLOS.append(que)


# ---------------------------------------------------------------- polinomios de Laurent: dict exponente -> entero
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


def sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


def ldet(M):
    n = len(M)
    tot = {}
    for p in permutations(range(n)):
        t = {0: sign(p)}
        for i in range(n):
            t = lmul(t, M[i][p[i]])
        tot = ladd(tot, t)
    return tot


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


# ---------------------------------------------------------------- C1
print('=' * 100)
print('C1  (2) como identidad EXACTA de polinomios de Laurent:  N_ℓ(z)·∏(1-z^(ρ,α)) = N_ρ(z)·z^-(λ,ρ)·∏(1-z^(ℓ,α))')
print('    con N_x(z) = det( z^{j x_i} - z^{-j x_i} ),  la recta a_P(z) = (z, z^2, …, z^m)')
for coraiz in (False, True):
    casos = buenos = 0
    for m in (2, 3, 4):
        rho = [m - i for i in range(m)]
        Nr = ldet([[{(j + 1) * e: 1, -(j + 1) * e: -1} for j in range(m)] for e in rho])
        Pr = prod_uno_menos(pairings(rho, coraiz))
        for lam in dominantes(m, 4 if m < 4 else 3):
            ell = [lam[i] + rho[i] for i in range(m)]
            Nl = ldet([[{(j + 1) * e: 1, -(j + 1) * e: -1} for j in range(m)] for e in ell])
            lr = sum(lam[i] * rho[i] for i in range(m))
            izq = lmul(Nl, Pr)
            der = lmul(lmul(Nr, {-lr: 1}), prod_uno_menos(pairings(ell, coraiz)))
            casos += 1
            buenos += (izq == der)
    if not coraiz:
        ok(buenos == casos, 'C1: la identidad (2), estadístico de raíces', '(%d de %d pesos, m = 2, 3, 4)' % (buenos, casos))
    else:
        ok(buenos < casos, 'C1 control: con el estadístico de corraíces la identidad debe fallar',
           '(se cumple en %d de %d; λ = 0 la cumple siempre)' % (buenos, casos))


# ---------------------------------------------------------------- Z[ζ_q] = Z[X]/Φ_q
X = symbols('X')


class Ciclo:
    def __init__(self, q):
        self.q = q
        self.phi = [int(c) for c in Poly(cyclotomic_poly(q, X), X).all_coeffs()[::-1]]  # c_0 … c_n, mónico
        self.n = len(self.phi) - 1

    def red(self, coefs):
        """coefs: lista densa c_0…c_N (enteros) → resto mod Φ_q (Φ_q mónico: división exacta en Z)."""
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
        """Laurent en z → elemento de Z[ζ_q] (exponentes mod q)."""
        c = [0] * self.q
        for e, v in d.items():
            c[e % self.q] += v
        return self.red(c)

    def mul(self, a, b):
        c = [0] * (2 * self.n)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    c[i + j] += x * y
        return self.red(c)

    def uno(self):
        return tuple([1] + [0] * (self.n - 1))

    def potencia(self, r):
        c = [0] * self.q
        c[r % self.q] = 1
        return self.red(c)

    def es_entero(self, a):
        return all(x == 0 for x in a[1:])


def ldiv_exacta(num, den):
    """num/den de Laurent con cociente de Laurent (división exacta; error si no lo es)."""
    num = dict(num)
    dmax = max(den)
    dmin = min(den)
    lc = den[dmax]
    q = {}
    while num:
        e = max(num)
        c = num[e]
        if c % lc:
            raise ValueError('no exacta')
        t = c // lc
        s = e - dmax
        if e - dmax + dmin < min(num) - 10 ** 6:
            raise ValueError('no termina')
        q[s] = t
        num = ladd(num, lmul({s: t}, den), -1)
        if num and max(num) < dmax + min(q) - 2 * (dmax - dmin) - 10 ** 4:
            raise ValueError('no exacta')
    return q


def chi_laurent(lam):
    m = len(lam)
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    lr = sum(lam[i] * rho[i] for i in range(m))
    num = lmul({-lr: 1}, prod_uno_menos(pairings(ell)))
    return ldiv_exacta(num, prod_uno_menos(pairings(rho)))


def espectro_estable(m, q):
    A = Counter()
    for j in range(1, m + 1):
        A[j % q] += 1
        A[(-j) % q] += 1
    for p in range(2, q):
        if gcd(p, q) == 1 and Counter({(p * r) % q: n for r, n in A.items()}) != A:
            return False
    return True


def cond_I(m, q, con6=True):
    return (2 * m) % q == 0 or (2 * m + 1) % q == 0 or (2 * m + 2) % q == 0 or (con6 and q == 6 and m % 3 == 1)


print('=' * 100)
print('C2  Teorema 5.1: χ_λ(g_{1/q}) EXACTO en Z[ζ_q] (χ como polinomio de Laurent por división exacta, reducido mod Φ_q)')
CIC = {q: Ciclo(q) for q in range(2, 13)}
pares = acuerdo_g = acuerdo_I = acuerdo_I_sin6 = 0
chis = 0
for m in range(2, 6):
    leyes = [chi_laurent(lam) for lam in dominantes(m, 5)]
    chis += len(leyes)
    for q in range(2, 13):
        C = CIC[q]
        todos = all(C.es_entero(C.de_laurent(L)) for L in leyes)
        pares += 1
        acuerdo_g += (todos == espectro_estable(m, q))
        acuerdo_I += (todos == cond_I(m, q))
        acuerdo_I_sin6 += (todos == cond_I(m, q, con6=False))
print('   %d caracteres exactos (m = 2..5, λ_1 ≤ 5), %d pares (m, q), q = 2..12' % (chis, pares))
ok(acuerdo_g == pares, 'C2: «todos enteros» = espectro estable por Galois', '(%d de %d)' % (acuerdo_g, pares))
ok(acuerdo_I == pares, 'C2: «todos enteros» = condición (I) del Teorema 5.1', '(%d de %d)' % (acuerdo_I, pares))
ok(acuerdo_I_sin6 < pares, 'C2 control: sin la cláusula q = 6 la condición debe discrepar',
   '(acuerda en %d de %d)' % (acuerdo_I_sin6, pares))

# ---------------------------------------------------------------- C3
print('=' * 100)
print('C3  Proposición 8.2: M_s^2 = ∏_w (2 - w - w^-1) exacto en Z[ζ_q], w recorriendo los factores no nulos')
print('    {y/x : x<y en U_s} ∪ {xy : x<y, xy ≠ 1} ∪ {x^2 : x ∈ U_s, x^2 ≠ 1},  U_s = μ_q \\ {ζ^s}')


def multiconjunto(q, s):
    U = [r for r in range(q) if r != s % q]
    W = []
    for i in range(len(U)):
        for j in range(i + 1, len(U)):
            W.append((U[j] - U[i]) % q)
            if (U[i] + U[j]) % q:
                W.append((U[i] + U[j]) % q)
    for x in U:
        if (2 * x) % q:
            W.append((2 * x) % q)
    return W


todos_ok = True
control_dep = 0
valores = {}
for q in range(3, 16):
    C = Ciclo(q)
    M2, Z = [], []
    for s in range(q):
        W = multiconjunto(q, s)
        p = C.uno()
        z = C.uno()
        for w in W:
            f = C.red([2] + [0] * q)
            f = tuple(a - b - c for a, b, c in zip(f, C.potencia(w), C.potencia(-w)))
            p = C.mul(p, f)
            z = C.mul(z, tuple(a - b for a, b in zip(C.uno(), C.potencia(w))))
        M2.append(p)
        Z.append(z)
    iguales = all(v == M2[0] for v in M2) and C.es_entero(M2[0])
    todos_ok &= iguales
    control_dep += (len(set(Z)) > 1)
    valores[q] = M2[0][0] if C.es_entero(M2[0]) else None
print('   M_0^2 por q:', ', '.join('q=%d: %s' % (q, v) for q, v in valores.items()))
ok(todos_ok, 'C3: M_s^2 es un entero racional y no depende de s, q = 3..15 (exacto)')
ok(control_dep == len(valores), 'C3 control: el producto complejo sin módulos ∏(1-w) depende de s',
   '(en %d de %d valores de q)' % (control_dep, len(valores)))

# ---------------------------------------------------------------- C4
print('=' * 100)
print('C4  Lema de filas (Apéndice): D_{P(S)h}(b) = (∏_i P(T_i)) D_h(b), EXACTO con h racional aleatoria')
print('    D_h(b) = det( h_{b_i+j} + h_{b_i-j} )_{i,j=1..m};  P = (1-S)^2 y (1-S^2)  (antes: revision_47 R1, a 1e-9)')
import random
from fractions import Fraction


def fdet(M):
    M = [r[:] for r in M]
    n = len(M)
    d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return d


def Dh(h, b):
    H = lambda n: h[n] if 0 <= n < len(h) else Fraction(0)
    m = len(b)
    return fdet([[H(b[i] + j) + H(b[i] - j) for j in range(1, m + 1)] for i in range(m)])


def aplica(h, coef):
    return [sum(coef[r] * (h[n - r] if 0 <= n - r < len(h) else 0) for r in range(len(coef))) for n in range(len(h))]


def derecha(h, b, coef):
    from itertools import product as ip
    tot = Fraction(0)
    for eps in ip(range(len(coef)), repeat=len(b)):
        c = 1
        for e in eps:
            c *= coef[e]
        if c:
            tot += c * Dh(h, [b[i] - eps[i] for i in range(len(b))])
    return tot


random.seed(3)
casos = buenos = ctrl = 0
for coef in ([1, -2, 1], [1, 0, -1]):
    for _ in range(40):
        m = random.choice([2, 3, 4])
        h = [Fraction(random.randint(-9, 9), random.randint(1, 9)) for _ in range(60)]
        b = sorted(random.sample(range(6, 30), m), reverse=True)
        casos += 1
        buenos += Dh(aplica(h, coef), b) == derecha(h, b, coef)
for _ in range(20):
    m = random.choice([2, 3])
    h = [Fraction(random.randint(-9, 9), random.randint(1, 9)) for _ in range(60)]
    b = sorted(random.sample(range(6, 30), m), reverse=True)
    ctrl += Dh(aplica(h, [1, -2, 1]), b) != derecha(h, b, [1, 2, 1])
ok(buenos == casos, 'C4: el lema de filas, igualdad exacta', '(%d de %d)' % (buenos, casos))
ok(ctrl == 20, 'C4 control: con (1+S)^2 en la derecha debe fallar', '(falla en %d de 20)' % ctrl)

print('=' * 100)
print('RESULTADO:', 'las cuatro certificaciones exactas pasan y sus controles fallan.' if not FALLOS
      else 'FALLAN: ' + '; '.join(FALLOS))
if __name__ == '__main__':
    sys.exit(1 if FALLOS else 0)
