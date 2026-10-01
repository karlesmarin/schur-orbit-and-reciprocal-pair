# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: filas negativas J = -k de la matriz de Mahler como MOMENTOS NEGATIVOS sobre las unidades:
#        fila -k (rama plena)    N_k = int_{Z_p^x} x^{-k} dmu_chi,
#        fila -k (rama <x>)      Nb_k = int_{Z_p^x} <x>^{-k} dmu_chi   (x = omega(x) <x>),
#      mu_chi = medida de H_chi = -2 (F_chi(Y) - chi(2) F_chi(Y^2)) (sum_i A_{J,i} C(x,i) = [t^J] H_chi(E^x), M0).
#   Medida exacta de los cosets: en la deduccion de (8)-(9) (mahler_masas), F_chi = sum_{a<M} Y^a R_a(Y^M) con
#      mu_F(a + M Z_p) = R_a(1) = w_a = -chi(M) (B_chi/n + C_chi(a/M mod n)) EXACTO en Z_(p)[zeta_m] (la reduccion mod p solo
#      entraba al pasar a T). mu_{F(Y^2)}(a + M Z_p) = w_{a/2 mod M}. Con M = p^N la suma de Riemann sobre (Z/p^N)^x de una
#      funcion constante mod p^N en cosets (x^{-k}, <x>^{-k}) es EXACTA mod p^N (mu entera).
#   Lectura: por Kummer/Kubota-Leopoldt, int_{Z_p^x} x^j dmu = -2 (1 - 2^j chi(2)) (1 - chi(p) p^j) L(-j, chi); la rama
#      <x>^{-1} (j = 0 mod p-1, j -> -1) da Nb_1 = -2 (1 - chi(2) <2>^{-1}) L_p(1, chi omega)  (Leopoldt, chi omega PAR).
#      Para chi impar, la rama plena x^{-1} es un cero trivial (paridad). Para zeta (G = Y/(1+Y)): N_1 = (1 - 1/p) log_p 2.
#      Mod p, Nb_1 = int_{Z_p^x} 1 dmu = (1 - chi(p)) beta_0: la informacion nueva esta en la VALORACION.
#   Primo: P = (p, g(zeta_m)), g = Hensel de mahler_padico.primo_sobre_p (convencion del 2-oct).
#   Modos: control (aditividad de mu; momentos positivos = (1 - chi(p) p^j) beta_j; zeta: N_1 = (1-1/p) log_p 2);
#          tabla  (Nb_1, N_2, Nb_2 por orbita de results/_reg/tabla_capas/capas.csv, frente a lambda~ y vueltas);
#          iwasawa (serie f(T) = int omega^r (1+T)^{l(x)} dmu mod P, a_i = int C(l(x), i) dmu por Lucas; lambda = ord_T;
#                   controles con las ramas de zeta; luego lambda de la rama <x> frente a lambda~ de las vueltas).
# FECHA: 2-oct-2026
# Uso (docker, /work): sage -python filas_negativas.py control | tabla | iwasawa
import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, '/work')
import mahler_padico as MP  # noqa: E402
from sage.all import Qp, ZZ, kronecker_character  # noqa: E402

LIM = 2 ** 29  # p^N < LIM: productos < 2^58, sumas de f <= 16 productos caben en int64


def precision(p):
    N = 1
    while p ** (N + 1) < LIM and p ** (N + 1) <= 2 * 10 ** 6:
        N += 1
    return N


def anillo(m, p, N):
    """g (Hensel mod p^N), f, tabla de coordenadas de x^e (e < m) mod p^N."""
    g, f = MP.hensel(m, p, N)
    mod = p ** N
    X = g.parent().gen()
    pot, c = [], g.parent()(1)
    for _ in range(m):
        pot.append(([int(v) % mod for v in c.list()] + [0] * f)[:f])
        c = (c * X) % g
    return g, f, np.array(pot, dtype=np.int64)


def por_potencia(v, e, g, f, m, mod):
    """v (filas de coordenadas) * x^e en Z/p^N[x]/(g)."""
    X = g.parent().gen()
    mm = [([int(u) % mod for u in ((X ** (k + e % m)) % g).list()] + [0] * f)[:f] for k in range(f)]
    out = np.zeros_like(v)
    for k in range(f):
        for i in range(f):
            if mm[k][i]:
                out[:, i] = (out[:, i] + v[:, k] * mm[k][i] % mod) % mod
    return out


def medida(ex, m, n, p, N):
    """mu_H(a + p^N Z_p), a = 0..p^N-1, coordenadas mod p^N; y (g, f)."""
    mod = M = p ** N
    g, f, pot = anillo(m, p, N)
    exa = np.array(ex, dtype=np.int64)
    ch = np.where((exa >= 0)[:, None], pot[np.maximum(exa, 0)], 0)
    u = np.arange(n, dtype=np.int64)
    Bchi = ((u[:, None] * ch) % mod).sum(axis=0) % mod
    C = np.zeros((n, f), dtype=np.int64)
    C[1:] = np.cumsum(ch[:-1], axis=0) % mod
    ninv = pow(int(n), -1, mod)
    Minv = pow(int(M), -1, int(n))
    c = (np.arange(M, dtype=np.int64) * Minv) % n
    v = (Bchi * ninv % mod + C[c]) % mod
    w = (-por_potencia(v, ex[M % n], g, f, m, mod)) % mod  # w_a = -chi(M) v_a
    inv2 = pow(2, -1, mod)
    w2 = w[(np.arange(M, dtype=np.int64) * inv2) % M]  # w_{a/2}
    mu = (-2 * (w - por_potencia(w2, ex[2], g, f, m, mod))) % mod
    return mu, g, f


def potmod(a, e, mod):
    """a^e mod mod, vectorial (a < 2^29)."""
    r = np.ones_like(a)
    b = a % mod
    while e:
        if e & 1:
            r = r * b % mod
        b = b * b % mod
        e >>= 1
    return r


def momento(mu, valores, mod):
    """sum_a valores[a] mu[a] mod p^N (valores = 0 fuera de las unidades)."""
    return [int(((valores * mu[:, k]) % mod).sum() % mod) for k in range(mu.shape[1])]


def funciones(p, N):
    mod = p ** N
    a = np.arange(mod, dtype=np.int64)
    uni = (a % p) != 0
    phi = (p - 1) * p ** (N - 1)
    om = potmod(a, p ** (N - 1), mod)  # Teichmuller mod p^N
    inv = potmod(a, phi - 1, mod)  # a^{-1}
    br = inv * om % mod  # <a>^{-1} = omega(a) / a
    fx = lambda k: np.where(uni, potmod(inv, k, mod), 0)  # noqa: E731
    fb = lambda k: np.where(uni, potmod(br, k, mod), 0)  # noqa: E731
    return fx, fb, uni, a


def val(coords, p, N):
    """v_P de un elemento en coordenadas (no ramificado): min v_p de las coordenadas; N si es 0 mod p^N."""
    vs = [ZZ(c).valuation(p) for c in coords if c % p ** N]
    return min(vs) if vs else N


def _control():
    # (i) aditividad: sum_{b = a mod p^N} mu^{(N+1)}(b) = mu^{(N)}(a)
    fallos = casos = 0
    for n, p in ((7, 5), (11, 7), (13, 7), (19, 5), (31, 11)):
        obs, _ = MP.orbitas(n, p)
        for chi, orb, m, Km, chim in obs:
            ex = MP.exponentes(chim, Km, n)
            mu3, _, _ = medida(ex, m, n, p, 3)
            mu2, _, _ = medida(ex, m, n, p, 2)
            M2 = p ** 2
            agg = np.zeros_like(mu2)
            for r in range(p):
                agg = (agg + mu3[r * M2:(r + 1) * M2]) % M2
            casos += 1
            fallos += not np.array_equal(agg % M2, mu2 % M2)
    print("(i) aditividad de la medida (N=2 vs N=3): %d orbitas, %d fallos" % (casos, fallos), flush=True)
    # (ii) momentos positivos: sum_{unidades} a^j mu = (1 - chi(p) p^j) beta_j mod p^N, j = 0..8
    fallos = casos = 0
    for n, p in ((7, 5), (11, 7), (13, 7), (19, 5), (31, 11), (23, 7)):
        N = 4
        mod = p ** N
        obs, _ = MP.orbitas(n, p)
        for chi, orb, m, Km, chim in obs:
            ex = MP.exponentes(chim, Km, n)
            mu, g, f = medida(ex, m, n, p, N)
            a = np.arange(mod, dtype=np.int64)
            uni = (a % p) != 0
            bc = MP.beta_masas(ex, m, n, p, N, 8, g, f)
            for j in range(9):
                lhs = momento(mu, np.where(uni, potmod(a, j, mod), 0), mod)
                b = np.array([bc[j]], dtype=np.int64)
                rhs = (b - por_potencia(b, ex[p % n], g, f, m, mod) * (p ** j % mod) % mod) % mod
                casos += 1
                if lhs != [int(x) for x in rhs[0]]:
                    fallos += 1
                    print("  (ii) FALLA n=%d p=%d j=%d: %s vs %s" % (n, p, j, lhs, rhs[0]), flush=True)
    print("(ii) momentos positivos frente a (1 - chi(p) p^j) beta_j: %d casos, %d fallos" % (casos, fallos), flush=True)
    # (iii) zeta: N_1 = int_{Z_p^x} x^{-1} dmu_G = (1 - 1/p) log_p 2 (mod p^{N-1})
    fallos = casos = 0
    for p in (5, 7, 11, 13, 17, 19, 23, 1093):
        N = precision(p) if p < 1000 else 2
        mod = p ** N
        fx, fb, uni, a = funciones(p, N)
        at =np.where(a == 0, mod, a)  # representante en [1, p^N]
        inv2 = pow(2, -1, mod)
        mu = np.where((at - 1) % 2 == 0, inv2, mod - inv2)[:, None] % mod
        N1 = momento(mu, fx(1), mod)[0]
        ref = (1 - 1 / Qp(p, N + 2)(p)) * Qp(p, N + 2)(2).log()
        refm = int(ZZ(ref.lift()) % p ** (N - 1)) if ref.valuation() >= 0 else None
        casos += 1
        ok = refm is not None and N1 % p ** (N - 1) == refm
        fallos += not ok
        q2 = (pow(2, p - 1, p * p) - 1) // p % p
        print("  p=%d N=%d: N_1 = %d, (1-1/p) log_p 2 = %s mod p^%d: %s; v_p(N_1) = %d; cociente de Fermat q_2 mod p = %d" % (
            p, N, N1 % p ** (N - 1), refm, N - 1, 'OK' if ok else 'FALLA', val([N1], p, N), q2), flush=True)
    print("(iii) zeta fila -1 frente a (1 - 1/p) log_p 2: %d primos, %d fallos" % (casos, fallos), flush=True)


def _capas_csv():
    """La tabla de capas: en el paquete publicado esta en datos/tabla_capas/; en el arbol de trabajo, en results/_reg/."""
    base = os.path.dirname(os.path.abspath(__file__))
    for rel in (('datos', 'tabla_capas', 'capas.csv'), ('results', '_reg', 'tabla_capas', 'capas.csv')):
        ruta = os.path.join(base, *rel)
        if os.path.exists(ruta):
            return ruta
    raise FileNotFoundError('capas.csv no esta ni en datos/tabla_capas/ ni en results/_reg/tabla_capas/')


def _tabla():
    filas = list(csv.DictReader(open(_capas_csv(), encoding='utf-8')))
    ref = {(int(r['n']), int(r['p']), r['conrey']): r for r in filas if r['a'] == '2'}
    ref3 = {(int(r['n']), int(r['p']), r['conrey']): r for r in filas if r['a'] == '3'}
    print("n p conrey chi(p)=1 | vueltas piso2 piso3 lambda~ | v(Nb_1) v(1-chi(2)<2>^-1) v(L_p(1,chi omega)) | "
          "v(N_2) v(Nb_2)", flush=True)
    pares = sorted({(k[0], k[1]) for k in ref})
    for n, p in pares:
        N = precision(p)
        mod = p ** N
        fx, fb, uni, a = funciones(p, N)
        obs, _ = MP.orbitas(n, p)
        for chi, orb, m, Km, chim in obs:
            conrey = ';'.join(str(o.conrey_number()) for o in sorted(orb, key=lambda o: o.conrey_number()))
            ex = MP.exponentes(chim, Km, n)
            mu, g, f = medida(ex, m, n, p, N)
            Nb1 = momento(mu, fb(1), mod)
            # factor (1 - chi(2) <2>^{-1}) en Z/p^N[x]/(g)
            br2 = int(fb(1)[2])
            uno = np.zeros((1, f), dtype=np.int64)
            uno[0, 0] = 1
            fac = (uno - por_potencia(uno, ex[2], g, f, m, mod) * br2 % mod) % mod
            vfac = val(fac[0], p, N)
            v1 = val(Nb1, p, N)
            r2, r3 = ref.get((n, p, conrey)), ref3.get((n, p, conrey))
            print("%d %d %s %s | %s %s %s | %d %d %s | %d %d" % (
                n, p, conrey, r2['chip_1'] if r2 else '?', r2['vueltas'] if r2 else '-', r3['vueltas'] if r3 else '-',
                r3['lambda_tilde_min'] if r3 else (r2['lambda_tilde_min'] if r2 else '-'), v1, vfac,
                (v1 - vfac) if v1 < N else '>=%d' % (N - vfac), val(momento(mu, fx(2), mod), p, N),
                val(momento(mu, fb(2), mod), p, N)), flush=True)


def ell_mod(p, N):
    """l(a) = log_p <a> / log_p(1+p) mod p^{N-1} para a unidad mod p^N (0 en las no unidades); y la mascara de unidades.
    log(1+z) = sum (-1)^{k+1} z^k/k con z = <a> - 1 (v >= 1): z^k/k es correcto mod p^N aunque z solo se conozca mod p^N
    (la variacion es multiplo de p^{N+k-1}); se trabaja mod p^{N+2} para dividir por la parte p de k."""
    mod, mod2 = p ** N, p ** (N + 2)
    a = np.arange(mod, dtype=np.int64)
    uni = (a % p) != 0
    om = potmod(a, p ** (N - 1), mod)
    br = potmod(a, (p - 1) * p ** (N - 1) - 1, mod) * om % mod  # <a>^{-1}
    ang = potmod(br, (p - 1) * p ** (N - 1) - 1, mod)  # <a>
    z = np.where(uni, (ang - 1) % mod, 0)
    lg = np.zeros_like(a)
    zk = np.ones_like(a)
    for k in range(1, 2 * N + 5):  # v(z^k/k) >= k - v_p(k) >= N mas alla; v_p(k) <= 2 para p >= 5 en este rango
        zk = zk * z % mod2
        e = ZZ(k).valuation(p)
        assert e <= 2
        t = (zk // p ** e) % mod * pow(int(k // p ** e), -1, mod) % mod
        lg = (lg + (t if k % 2 else mod - t)) % mod
    lu = Qp(p, N + 2)(1 + p).log()
    lu_p = int(ZZ((lu / p).lift()) % p ** (N - 1))  # log u / p, unidad
    ell = ((lg // p) % p ** (N - 1)) * pow(lu_p, -1, p ** (N - 1)) % p ** (N - 1)
    return np.where(uni, ell, 0), uni


def binom_lucas(ell, i, p, r):
    """C(ell, i) mod p (Lucas), ell < p^r vector."""
    tab = np.array([[ZZ(x).binomial(y) % p for y in range(p)] for x in range(p)], dtype=np.int64)
    out = np.ones_like(ell)
    e, ii = ell.copy(), i
    for _ in range(r):
        out = out * tab[e % p, ii % p] % p
        e //= p
        ii //= p
    return out


def serie_iwasawa(mu, p, N, rama, imax):
    """a_i mod P (coordenadas mod p), i = 0..imax, de f(T) = int_{Z_p^x} omega(x)^rama (1+T)^{l(x)} dmu."""
    mod = p ** N
    ell, uni = ell_mod(p, N)
    a = np.arange(mod, dtype=np.int64)
    om = potmod(a, p ** (N - 1), mod)
    w = np.where(uni, potmod(om, rama % (p - 1), mod), 0) % p
    mup = mu % p
    out = []
    for i in range(imax + 1):
        c = binom_lucas(ell, i, p, N - 1) * w % p
        out.append([int((c * mup[:, k] % p).sum() % p) for k in range(mu.shape[1])])
    return out


def lam(coefs):
    for i, c in enumerate(coefs):
        if any(c):
            return i
    return None


def mu_zeta(p, N):
    mod = p ** N
    a = np.arange(mod, dtype=np.int64)
    at = np.where(a == 0, mod, a)
    inv2 = pow(2, -1, mod)
    return np.where((at - 1) % 2 == 0, inv2, mod - inv2)[:, None] % mod


def _iwasawa():
    # controles con zeta (ramas omega^r, r impar: x^j con j = r mod p-1)
    print("CONTROL zeta: lambda de la rama omega^r (r impar) = ord_T(f mod p)", flush=True)
    for p, esperado in ((5, {}), (7, {}), (11, {}), (13, {}), (17, {7: 1}), (37, {31: 1}), (31, {9: 1, 19: 1})):
        N = precision(p)
        mu = mu_zeta(p, N)
        res = {}
        for r in range(1, p - 1, 2):
            if r == p - 2:
                continue  # rama del polo (j = -1 mod p-1): s = 1 es el polo de zeta_p
            res[r] = lam(serie_iwasawa(mu, p, N, r, 6))
        malos = {r: l for r, l in res.items() if l != esperado.get(r, 0)}
        print("  p=%d N=%d: ramas con lambda > 0: %s; esperado %s: %s" % (
            p, N, {r: l for r, l in res.items() if l}, esperado, 'OK' if not malos else 'FALLA %s' % malos), flush=True)
    # tabla: rama de las vueltas (<x>, rama 0) frente a lambda~ y v(Nb_1)
    filas = list(csv.DictReader(open(_capas_csv(), encoding='utf-8')))
    ref3 = {(int(r['n']), int(r['p']), r['conrey']): r for r in filas if r['a'] == '3'}
    tot = iguales = 0
    for n, p in sorted({(k[0], k[1]) for k in ref3}):
        N = precision(p)
        mod = p ** N
        fx, fb, uni, a = funciones(p, N)
        obs, _ = MP.orbitas(n, p)
        for chi, orb, m, Km, chim in obs:
            conrey = ';'.join(str(o.conrey_number()) for o in sorted(orb, key=lambda o: o.conrey_number()))
            ex = MP.exponentes(chim, Km, n)
            mu, g, f = medida(ex, m, n, p, N)
            l_s = lam(serie_iwasawa(mu, p, N, 0, 12))
            v1 = val(momento(mu, fb(1), mod), p, N)
            lt = ref3[(n, p, conrey)]['lambda_tilde_min']
            tot += 1
            iguales += str(l_s) == lt
            if l_s or str(l_s) != lt:
                print("  n=%d p=%d conrey=%s: lambda serie = %s | lambda~ vueltas = %s | v(fila -1) = %d%s" % (
                    n, p, conrey, l_s, lt, v1, '' if str(l_s) == lt else '  <-- DIFIERE'), flush=True)
    print("lambda (serie de Iwasawa mod p) == lambda~ (vueltas): %d de %d" % (iguales, tot), flush=True)


if __name__ == '__main__':
    t0 = time.time()
    {'control': _control, 'tabla': _tabla, 'iwasawa': _iwasawa}[sys.argv[1] if len(sys.argv) > 1 else 'control']()
    print("(%.1f s)" % (time.time() - t0), flush=True)
