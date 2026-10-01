# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: instrumento de Mahler v2 (p-adico entero). Misma matematica que mahler_bernoulli (forma cerrada M0,
#      PRUEBA_VUELTAS.md): A = tau * Lambda diag(beta_chi) R. Cambios:
#      (1) Lambda y R son universales (solo p). Con D = diag(p^{floor(J/(p-1))}), Lambda' = D Lambda D^-1 y R' = D R D^-1
#          son p-ENTERAS (medido: _escala_test, p = 5, 7, J <= 400; y la guarda lo recomprueba): se guardan mod p^N,
#          N = floor(Jmax/(p-1)) + 2.
#      (2) beta_chi = 2 (1 - 2^j chi(2)) B_{j+1,chi}/(j+1) es p-entero (p no divide n); se calcula por CARACTER (no por
#          clase) y se lleva a Z/p^N[x]/(g), g = levantamiento de Hensel de un factor irreducible de Phi_m mod p
#          (m = orden de chi; el factor es el de primo_sobre_p: CONVENCION P = (p, g(zeta_m)), 2-oct). A' = sum_k x^k Lambda' diag(beta^(k)) R'  (f productos de matrices enteras, FLINT).
#      (3) A_{J,i} = A'_{J,i} / p^{floor(J/(p-1)) - floor(i/(p-1))}: se comprueba la divisibilidad (guarda) y se reduce
#          a F_{p^f} (x -> raiz de g mod p). Filas pivote = capas presentes (M1).
#      Caracteres: orbitas de Frobenius de DirichletGroup(n) (n moderado) o el cuadratico de Kronecker (n grande).
# FECHA: 1-oct-2026
# Uso (docker, /work): sage -python mahler_padico.py control
import sys
import time

sys.path.insert(0, '/work')
from sage.all import (ZZ, QQ, GF, PolynomialRing, PowerSeriesRing, matrix, diagonal_matrix, factorial, stirling_number2,
                      euler_phi, cyclotomic_polynomial, bernoulli, kronecker, Integers, DirichletGroup,
                      CyclotomicField)  # noqa: E402

_CACHE = {}


def universales(p, Jmax, imax):
    """Lambda', R' enteros (en ZZ, reducidos mod p^N) y N."""
    clave = (p, Jmax, imax)
    if clave in _CACHE:
        return _CACHE[clave]
    N = Jmax // (p - 1) + 2
    mod = ZZ(p) ** N
    prec = Jmax + 2
    R = PowerSeriesRing(QQ, 't', default_prec=prec)
    t = R.gen()
    h = sum(t ** k / factorial(k) for k in range(p))
    hm = sum((-t) ** k / factorial(k) for k in range(p))
    ell = ((h / hm).sqrt()).log()
    sh = [J // (p - 1) for J in range(Jmax + 1)]
    L = matrix(ZZ, Jmax + 1, Jmax + 1)
    pot = R(1).add_bigoh(prec)
    for j in range(Jmax + 1):
        cs = pot.list()
        for J in range(j, min(Jmax + 1, len(cs))):
            c = cs[J] * QQ(p) ** (sh[J] - sh[j])
            if c.denominator() % p == 0:
                raise ValueError("Lambda' no p-entera en (%d,%d)" % (J, j))
            L[J, j] = (c.numerator() * c.denominator().inverse_mod(mod)) % mod
        pot = pot * ell
    Rm = matrix(ZZ, Jmax + 1, imax)
    for j in range(Jmax + 1):
        for i in range(min(j + 1, imax)):
            c = factorial(i) * stirling_number2(j, i) / factorial(j) * QQ(p) ** (sh[j] - sh[i])
            if c.denominator() % p == 0:
                raise ValueError("R' no p-entera en (%d,%d)" % (j, i))
            Rm[j, i] = (c.numerator() * c.denominator().inverse_mod(mod)) % mod
    _CACHE[clave] = (L, Rm, N, sh)
    return _CACHE[clave]


def primo_sobre_p(m, p):
    """CONVENCION DE PRIMO (2-oct-2026). P = (p, g(zeta_m)) en Q(zeta_m), zeta_m = exp(2 pi i/m) (la de Sage y la de
    las etiquetas de Conrey), con g el factor monico irreducible de Phi_m mod p de MENOR tupla de coeficientes
    (c_0, c_1, ..., c_{f-1}) en [0, p), orden lexicografico empezando por el termino constante. No depende del orden
    interno de factor() ni de n: solo de (m, p). Cuerpo residual F_p[x]/(g), zeta_m -> x.
    Devuelve g en F_p[x]."""
    Fx = PolynomialRing(GF(p), 'x')
    fs = [h for h, _ in Fx(cyclotomic_polynomial(m)).factor()]
    return min(fs, key=lambda h: tuple(int(c) for c in h.list()))


def reducir(coefs, m, p):
    """elemento de Q(zeta_m) (coeficientes racionales en la base de potencias de zeta_m) -> F_p[x]/(g), g =
    primo_sobre_p(m, p). Exige p-integralidad de los coeficientes (si no, ValueError)."""
    g = primo_sobre_p(m, p)
    Fx = g.parent()
    acc = Fx(0)
    for k, c in enumerate(coefs):
        c = QQ(c)
        if c.denominator() % p == 0:
            raise ValueError("coeficiente no p-entero")
        acc += Fx.base_ring()(c.numerator()) / Fx.base_ring()(c.denominator()) * Fx.gen() ** k
    return acc % g


def hensel(m, p, N):
    """g en Z[x] monico, levantamiento mod p^N del factor primo_sobre_p(m, p) de Phi_m mod p; y su grado f."""
    Fx = PolynomialRing(GF(p), 'x')
    Zx = PolynomialRing(ZZ, 'x')
    Phi = Zx(cyclotomic_polynomial(m))
    g0 = primo_sobre_p(m, p)
    f = g0.degree()
    # levantamiento de Hensel lineal (Phi = g * h mod p^k)
    g = Zx(g0)
    hh = Zx(Fx(Phi) // g0)
    mod = p
    s0, t0 = [Zx(c) for c in (Fx(g0).xgcd(Fx(hh))[1:])]  # s g + t h = 1 mod p
    for _ in range(N):
        mod *= p
        e = (Phi - g * hh) % mod
        # correccion: g += t*e mod g, h += s*e mod h
        dg = Zx([c % mod for c in ((t0 * e) % g).list()])
        dh = Zx([c % mod for c in ((s0 * e) % hh).list()])
        g = Zx([c % mod for c in (g + dg).list()])
        hh = Zx([c % mod for c in (hh + dh).list()])
    assert all(c % (p ** N) == 0 for c in (Phi - g * hh).list()), "Hensel fallido"
    return g, f


def beta_coords(valores_B, chi2, p, N, m, g, f):
    """beta_j en Z/p^N[x]/(g) como lista de f coordenadas enteras, para j = 0..Jmax.
    valores_B[j] = B_{j+1,chi} como elemento de Q(zeta_m) (lista de coef. racionales en la base de potencias) o racional."""
    mod = ZZ(p) ** N
    Zx = PolynomialRing(ZZ, 'x')
    out = []
    for j, Bk in enumerate(valores_B):
        b = [QQ(c) for c in Bk]  # coeficientes en zeta_m
        # beta = 2/(j+1) * (B - 2^j chi(2) B) en Q(zeta_m): chi(2) = zeta^e (e = chi2)
        e = chi2
        Bpoly = PolynomialRing(QQ, 'z')(b)
        z = PolynomialRing(QQ, 'z').gen()
        Phi = PolynomialRing(QQ, 'z')(cyclotomic_polynomial(m))
        tot = (Bpoly - 2 ** j * (Bpoly * z ** e) % Phi) * QQ(2) / (j + 1)
        cs = tot.list()
        # reducir mod (g, p^N): tot(z) con z -> x, coeficientes p-enteros
        acc = Zx(0)
        for k, c in enumerate(cs):
            if c == 0:
                continue
            if c.denominator() % p == 0:
                raise ValueError("beta_%d no p-entera" % j)
            acc += Zx.gen() ** k * ((c.numerator() * c.denominator().inverse_mod(mod)) % mod)
        acc = acc % g
        coefs = [ZZ(c) % mod for c in acc.list()] + [0] * f
        out.append(coefs[:f])
    return out


def beta_masas(ex, m, n, p, N, K, g, f):
    """beta_j mod p^N, j = 0..K, SIN Bernoulli (serie de F_chi(1+T); deduccion abajo).
    ex[u] = e con chi(u) = zeta_m^e (u = 0..n-1), o -1 si chi(u) = 0; g = hensel(m, p, N)[0].
      C_r = sum_a chi(a) C(a, r)                                   (P_chi(1+T) = sum_r C_r T^r),  r = 0..K+1
      f_r = -n^{-1} (C_{r+1} + sum_{s=1}^{r} C(n, s+1) f_{r-s})   (F_chi(1+T) = sum_r f_r T^r)
      phi_j = sum_{r<=j} r! S(j, r) f_r,  T_{j+1,r} = r (T_{j,r} + T_{j,r-1})   (F_chi(e^s) = sum_j phi_j s^j/j!)
      beta_j = -2 (1 - 2^j chi(2)) phi_j.
    Todo entero mod p^N; la unica division es por n (unidad). Devuelve la lista de coordenadas como beta_coords."""
    mod = ZZ(p) ** N
    modi = int(mod)
    Zx = PolynomialRing(ZZ, 'x')
    X = Zx.gen()
    xp, c = [], Zx(1)
    for _ in range(m):
        xp.append(([int(v) % modi for v in c.list()] + [0] * f)[:f])
        c = (c * X) % g
    # sumas de binomiales por clase: S[e][r] = sum_{a: chi(a) = zeta^e} C(a, r) mod p^N (Pascal por filas, enteros)
    R = K + 2
    S = [[0] * R for _ in range(m)]
    fila = [1] + [0] * (R - 1)  # C(a, r), a = 0
    for a in range(n):
        e = ex[a]
        if e >= 0:
            Se = S[e]
            for r in range(min(a + 1, R)):
                Se[r] += fila[r]
        for r in range(min(a + 1, R - 1), 0, -1):  # fila de a+1
            fila[r] = (fila[r] + fila[r - 1]) % modi
    Cr = [[sum(S[e][r] * xp[e][k] for e in range(m)) % modi for k in range(f)] for r in range(R)]
    ninv = pow(int(n), -1, modi)
    bin_n = [1]
    for s in range(1, R + 1):
        bin_n.append(bin_n[-1] * (int(n) - s + 1) // s)
    fr = []
    for r in range(K + 1):
        acc = list(Cr[r + 1])
        for s in range(1, r + 1):
            b = bin_n[s + 1] % modi
            acc = [(u + b * v) % modi for u, v in zip(acc, fr[r - s])]
        fr.append([(-ninv * u) % modi for u in acc])
    e2 = ex[2]
    out, T = [], [1]  # T[r] = r! S(j, r)
    for j in range(K + 1):
        if j > 0:
            T = [0] + [r * (T[r] + T[r - 1]) % modi if r < len(T) else r * T[r - 1] % modi for r in range(1, j + 1)]
        phi = [sum(T[r] * fr[r][k] for r in range(1, j + 1)) % modi if j else fr[0][k] for k in range(f)]
        ph = Zx(phi)
        tot = (-2 * (ph - 2 ** j * ((ph * X ** e2) % g))) % g
        out.append(([int(v) % modi for v in tot.list()] + [0] * f)[:f])
    return out


def perfil_caracter(valores_B, e2, m, n, p, a, Jmax):
    """filas pivote para un caracter dado por sus B_{k,chi} en Q(zeta_m) (base de potencias) y chi(2) = zeta_m^e2."""
    imax = min(Jmax + 1, p ** a)
    L, Rm, N, sh = universales(p, Jmax, imax)
    g, f = hensel(m, p, N)
    return perfil_desde_beta(beta_coords(valores_B, e2, p, N, m, g, f), m, p, a, Jmax)


def perfil_desde_beta(bc, m, p, a, Jmax):
    """filas pivote a partir de beta_j mod (p^N, g) en coordenadas (beta_coords o beta_masas)."""
    imax = min(Jmax + 1, p ** a)
    L, Rm, N, sh = universales(p, Jmax, imax)
    mod = ZZ(p) ** N
    g, f = hensel(m, p, N)
    Fq = GF(p ** f, 'y', modulus=PolynomialRing(GF(p), 'y')([c % p for c in g.list()])) if f > 1 else GF(p)
    y = Fq.gen() if f > 1 else Fq(1)
    A = matrix(Fq, Jmax + 1, imax)
    for k in range(f):
        Ak = (L * diagonal_matrix(ZZ, [bc[j][k] for j in range(Jmax + 1)]) * Rm)
        for J in range(Jmax + 1):
            for i in range(min(J + 1, imax)):
                c = Ak[J, i] % mod
                s = sh[J] - sh[i]
                if c % (p ** s):
                    raise ValueError("A'(%d,%d) no divisible por p^%d: la escala no cierra" % (J, i, s))
                v = (c // p ** s) % p
                if v:
                    A[J, i] += Fq(v) * (y ** k if f > 1 else 1)
    return A.pivot_rows()


def bernoulli_caracter(chi, K, m):
    """[B_{k,chi} para k = 1..K] como listas de coeficientes en Q(zeta_m) (base de potencias)."""
    out = []
    for k in range(1, K + 1):
        B = chi.bernoulli(k)
        out.append(list(B.list()) if hasattr(B, 'list') else [B])
    return out


def exponentes(chim, Km, n):
    """ex[u] = e con chi(u) = zeta_m^e (zeta_m = Km.gen()), o -1 si chi(u) = 0; u = 0..n-1."""
    m = Km.gen().multiplicative_order()
    tab = {Km.gen() ** e: e for e in range(m)}
    vals = chim.values()
    return [tab[v] if v != 0 else -1 for v in vals]


def beta_caracter(chim, Km, m, n, p, Jmax, metodo='masas'):
    """(beta_j mod (p^N, g) en coordenadas, e2). metodo 'masas' (beta_masas, por defecto desde el 2-oct) o
    'bernoulli' (chi.bernoulli exactos: camino de control)."""
    N = Jmax // (p - 1) + 2
    g, f = hensel(m, p, N)
    ex = exponentes(chim, Km, n)
    if metodo == 'masas':
        return beta_masas(ex, m, n, p, N, Jmax, g, f), ex[2]
    vals = [list(Km(chim.bernoulli(k)).list()) for k in range(1, Jmax + 2)]
    return beta_coords(vals, ex[2], p, N, m, g, f), ex[2]


def orbitas(n, p):
    """orbitas de Frobenius (chi ~ chi^p) de caracteres impares primitivos mod n: [(chi, orbita, m, Km, chim)] y
    #primitivos. Orden fijo (el de DirichletGroup de Sage)."""
    lam = Integers(n).unit_group_exponent()
    K = CyclotomicField(lam)
    chis = [c for c in DirichletGroup(n, K) if c(-1) == -1 and c.conductor() == n]
    vistos, out = set(), []
    for chi in chis:
        if chi in vistos:
            continue
        orb, c = [], chi
        while c not in orb:
            orb.append(c)
            c = c ** p
        vistos.update(orb)
        m = chi.order()
        Km = CyclotomicField(m)
        chim = chi.change_ring(Km) if m != lam else chi
        out.append((chi, orb, m, Km, chim))
    return out, len(chis)


def perfil(n, p, a, Jmax, metodo='masas'):
    """como mahler_bernoulli.perfil (sin certificados): (d, [(chi(2), chi(p), tamano orbita, filas pivote)], #prim)."""
    obs, npr = orbitas(n, p)
    sal = []
    for chi, orb, m, Km, chim in obs:
        bc, _ = beta_caracter(chim, Km, m, n, p, Jmax, metodo)
        sal.append((str(chi(2)), str(chi(p)), len(orb), perfil_desde_beta(bc, m, p, a, Jmax)))
    return euler_phi(n) // 2, sal, npr


def perfil_cuadratico(n, p, a, Jmax, metodo='masas'):
    """caracter de Kronecker (./n), n = 3 mod 4 (impar): valores racionales, sin recorrer clases."""
    from sage.all import kronecker_character
    chi = kronecker_character(-n) if n % 4 == 3 else None
    if chi is None or chi(-1) != -1:
        raise ValueError("se espera n = 3 mod 4 (caracter impar)")
    e2 = 0 if chi(2) == 1 else 1
    if metodo == 'masas':
        N = Jmax // (p - 1) + 2
        g, f = hensel(2, p, N)
        ex = [-1 if v == 0 else (0 if v == 1 else 1) for v in chi.values()]
        return perfil_desde_beta(beta_masas(ex, 2, n, p, N, Jmax, g, f), 2, p, a, Jmax), chi(2), chi(p)
    vals = [[QQ(chi.bernoulli(k))] for k in range(1, Jmax + 2)]
    return perfil_caracter(vals, e2, 2, n, p, a, Jmax), chi(2), chi(p)


def snf_regresion(bc, p, Jmax):
    """Prueba de regresion (solo f = 1): con TODAS las columnas 0..Jmax, SNF(D (A/tau) D^-1) = SNF(diag beta)
    sobre Z_p, pues D Lambda D^-1 y D R D^-1 son unitriangulares enteras (v_p(Lambda_{J,i}), v_p(R_{J,i}) >= -floor((J-i)/(p-1)): l' = E'/E entera y
    v_p(m!) <= (m-1)/(p-1)).
    Compara los exponentes de los divisores elementales de [A' | p^N I] (topados en N) con v_p(beta_j) topados."""
    L, Rm, N, sh = universales(p, Jmax, Jmax + 1)
    mod = ZZ(p) ** N
    Ap = (L * diagonal_matrix(ZZ, [bc[j][0] for j in range(Jmax + 1)]) * Rm)
    Ap = Ap.apply_map(lambda c: c % mod)
    ed = Ap.augment(mod * matrix.identity(ZZ, Jmax + 1)).elementary_divisors()
    izq = sorted(min(ZZ(d).valuation(p), N) for d in ed[:Jmax + 1])
    der = sorted(N if bc[j][0] % mod == 0 else min(ZZ(bc[j][0]).valuation(p), N) for j in range(Jmax + 1))
    return izq == der, izq, der


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'control':
    # control 31/31 frente al orden; argv[2] = 'masas' (defecto) o 'bernoulli'; y la regresion SNF (17) en las orbitas
    # con f = 1
    from mahler_rapido import CONTROL
    metodo = sys.argv[2] if len(sys.argv) > 2 else 'masas'
    tot = ok = 0
    snf_n = snf_ok = 0
    for n, p, a, med in CONTROL:
        P = p ** a
        js = sorted(c - P for c in med)
        t0 = time.time()
        d, sal, npr = perfil(n, p, a, max(js), metodo)
        dt = time.time() - t0
        for chi, orb, m, Km, chim in orbitas(n, p)[0]:
            if hensel(m, p, 2)[1] == 1:
                bc, _ = beta_caracter(chim, Km, m, n, p, max(js), metodo)
                bueno, izq, der = snf_regresion(bc, p, max(js))
                snf_n += 1
                snf_ok += bueno
                if not bueno:
                    print("  SNF (17) FALLA n=%d p=%d: %s vs %s" % (n, p, izq, der), flush=True)
        for c, w in sorted(med.items()):
            W = d - npr + sum(tam for _, _, tam, piv in sal if (c - P) in piv)
            tot += 1
            ok += (W == w)
            if W != w:
                print("FALLA n=%d p=%d capa %d: %d vs %d" % (n, p, c, W, w), flush=True)
        print("n=%d p=%d a=%d: %d orbitas, %.1f s" % (n, p, a, len(sal), dt), flush=True)
    assert tot > 0 and snf_n > 0
    print("CONTROL padico (%s): %d de %d; regresion SNF (17): %d de %d orbitas con f = 1" % (
        metodo, ok, tot, snf_ok, snf_n), flush=True)
