# -*- coding: utf-8 -*-
# LAS DOS PROPUESTAS DE LA VUELTA 29, verificadas antes de escribirlas.   16 de agosto de 2026.
#
# (A)  EL CARACTER CUADRATICO.  Una ruta para cerrar gamma_t es el determinante seno
#
#          D_t = det( zeta^{ij} - zeta^{-ij} )_{1<=i,j<=n}
#
#      afirmando  sigma_k(D_t) = gamma_t(k) D_t   y   D_t^2 = (-t)^n,  de donde
#
#          gamma_t(k) = (k/t)^n ,      n = (t-1)/2,
#
#      o sea trivial si n es par (t = 1 mod 4) y el caracter cuadratico de Q(sqrt(-t)) si n es impar.
#      Nota: coincide con lo que ya teniamos, (k/t)^{(t+3)/2}, porque (t+3)/2 = n+2 y el simbolo es
#      +-1.  Su forma es mejor: dice QUE CUERPO es.
#
# (B)  LA INYECTIVIDAD DE  w -> mu_w.  Nuestra fila "245/245 verified" seria un TEOREMA si el mapa
#      es inyectivo.  Su argumento: W^1 para  B_{R'} > B_{m'} x D_r  se lee como
#
#          (un subconjunto S de m' coordenadas)  x  (una quiralidad de D_r),      |W^1| = 2 C(R', m'),
#
#      y  mu_w + rho_{D_r}  recuerda el complemento S^c mientras el signo de la ultima coordenada
#      recuerda la quiralidad.  Luego mu_w determina w.
#
# LO QUE SE MIDE
#   A1  D_t^2 == (-t)^n   exacto en el cuerpo ciclotomico.
#   A2  sigma_k(D_t) == gamma_t(k) D_t   para TODA unidad k.
#   A3  gamma_t(k) == (k/t)^n .
#   B1  |W^1| == 2 C(R', m').
#   B2  w -> mu_w  INYECTIVO sobre todo W^1  (no solo sobre los supervivientes).
#   B3  y la lectura estructural: w -> (S, quiralidad) es una biyeccion con W^1.
#
# CONTROLES
#   C0  A1 y B2 son fatales.
#   C1  SENUELO de B: se prueba tambien  w -> eta_w,  que NO debe ser inyectivo (si lo fuera, el
#       argumento no distinguiria nada).
#   C2  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage v29_two_lemmas.sage

import json
import sys
import itertools

OUT = {}
print("=" * 112)
print("(A)  EL DETERMINANTE SENO Y EL CARACTER CUADRATICO")
print("=" * 112)


def plegar(v, t):
    v = int(v) % t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def sgn_perm(perm):
    n = len(perm)
    s = 1
    visto = [False] * n
    for i in range(n):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            s = -s
    return s


def gamma_de(k, t, n):
    cl, ep = [], []
    for j in range(1, n + 1):
        c_, e_ = plegar(k * j, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, n + 1)):
        return None
    s = sgn_perm([c - 1 for c in cl])
    for e in ep:
        s *= e
    return int(s)


A = []
for t in [3, 5, 7, 9, 11, 13, 15, 17, 19, 21]:
    n = (t - 1) // 2
    K = CyclotomicField(t)
    z = K.gen()
    M = matrix(K, n, n, lambda i, j: z ** ((i + 1) * (j + 1)) - z ** (-(i + 1) * (j + 1)))
    D = M.determinant()
    a1 = (D ** 2 == K((-t) ** n))
    unidades = [k for k in range(1, t) if gcd(k, t) == 1]
    a2 = a3 = 0
    for k in unidades:
        g = gamma_de(k, t, n)
        # sigma_k : zeta -> zeta^k
        sig = K.hom([z ** k])
        a2 += 1 if sig(D) == g * D else 0
        a3 += 1 if g == (jacobi_symbol(k, t) ** n) else 0
    print("  t=%2d  n=%2d :  A1  D^2 == (-t)^n : %-5s |  A2  sigma_k(D) = gamma(k) D : %2d/%2d"
          " |  A3  gamma = (k/t)^n : %2d/%2d"
          % (t, n, str(a1), a2, len(unidades), a3, len(unidades)))
    sys.stdout.flush()
    A.append({"t": int(t), "n": int(n), "A1": bool(a1),
              "A2": [int(a2), int(len(unidades))], "A3": [int(a3), int(len(unidades))]})
OUT["A"] = A

print("")
print("=" * 112)
print("(B)  LA INYECTIVIDAD DE  w -> mu_w  SOBRE TODO W^1")
print("=" * 112)


def coset_reps(v, mp, rr):
    """[(u, sgn, S)] con u = w(v) estrictamente H-dominante; S = indices del bloque congelado."""
    Rp = mp + rr
    out = []
    for perm in itertools.permutations(range(Rp)):
        s = sgn_perm(list(perm))
        base = [v[perm[i]] for i in range(Rp)]
        for eps in itertools.product((1, -1), repeat=Rp):
            u = [base[i] * eps[i] for i in range(Rp)]
            if not (all(u[i] > u[i + 1] for i in range(mp - 1)) and (mp == 0 or u[mp - 1] > 0)):
                continue
            f = u[mp:]
            if rr >= 2 and not (all(f[i] > f[i + 1] for i in range(rr - 1))
                                and f[rr - 2] > abs(f[rr - 1])):
                continue
            sg = s
            for e in eps:
                sg *= e
            S = frozenset(perm[i] for i in range(mp))
            out.append((tuple(u), int(sg), S))
    return out


def dominantes(rk, cota):
    def rec(k, tope):
        if k == 0:
            yield ()
            return
        for a in range(tope, -1, -1):
            for resto in rec(k - 1, a):
                yield (a,) + resto
    return list(rec(rk, cota))


B = []
for (t, r, cota) in [(3, 2, 4), (5, 2, 3), (3, 3, 3), (7, 2, 2), (5, 3, 2)]:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    rho2 = [2 * (Rp - i) - 1 for i in range(Rp)]
    n = b1 = b2 = b3 = decoy = 0
    esperado = 2 * binomial(Rp, mp)
    for Lam in dominantes(Rp, cota):
        v = [2 * Lam[i] + rho2[i] for i in range(Rp)]
        reps = coset_reps(v, mp, r)
        n += 1
        b1 += 1 if len(reps) == esperado else 0
        mus = [tuple(u[mp:]) for u, sg, S in reps]          # = 2(mu_w + rho_{D_r})
        etas = [tuple(u[:mp]) for u, sg, S in reps]
        b2 += 1 if len(set(mus)) == len(mus) else 0
        # B3: (S, quiralidad) determina w, y mu_w determina (S, quiralidad)
        claves = [(S, 1 if u[Rp - 1] > 0 else -1) for u, sg, S in reps]
        b3 += 1 if len(set(claves)) == len(claves) else 0
        # senuelo: eta_w NO debe ser inyectivo
        decoy += 1 if len(set(etas)) == len(etas) else 0
    print("  t=%d r=%d (m'=%d R'=%d)  %3d pesos |  B1 |W^1| == 2C(R',m')=%3d : %3d/%3d"
          " |  B2  mu_w inyectivo : %3d/%3d |  B3 (S,quiralidad) inyectivo : %3d/%3d"
          " |  SENUELO eta_w inyectivo (debe ser bajo) : %3d/%3d"
          % (t, r, mp, Rp, n, esperado, b1, n, b2, n, b3, n, decoy, n))
    sys.stdout.flush()
    B.append({"t": int(t), "r": int(r), "n": int(n), "esperado": int(esperado),
              "B1": int(b1), "B2": int(b2), "B3": int(b3), "senuelo_eta": int(decoy)})
OUT["B"] = B

json.dump(OUT, open("v29_two_lemmas_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("  LECTURA, escrita ANTES de correr:")
print("   * si A1-A3 salen limpias, gamma_t es el caracter cuadratico de Q(sqrt(-t)) y se cierra.")
print("   * si B2 sale limpia, nu en {0,+-1} deja de ser medido y pasa a TEOREMA, y (L1) se queda")
print("     con una unica dificultad: la division por Delta_t.")
print("   * si el senuelo eta_w tambien sale inyectivo, el argumento de B no distingue nada.")
print("=" * 112)
print("DONE")
