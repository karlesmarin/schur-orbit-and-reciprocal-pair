# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: simplificacion de la matriz de Mahler A_{J,i} = pi_chi [t^J]((E_w - 1)^i U)(alpha^c)   (PRUEBA_VUELTAS.md, M1).
#      Conjetura de forma cerrada (Lema de Abel en forma generatriz):
#        sum_c conj(chi)(c) U(zeta^c e^s) = tau(conj chi) Phi_chi(s),
#        Phi_chi(s) = sum_j beta_j s^j / j!,   beta_j = 2 (1 - 2^j chi(2)) B_{j+1,chi} / (j+1),
#      y con E(t) = e^{l(t)}, E^x = e^{x l(t)}:  sum_i A_{J,i} C(x,i) = tau [t^J] Phi_chi(x l(t)), es decir
#        A_{J,i} = tau * sum_{j<=J} beta_j (i!/j!) S(j,i) [t^J] l(t)^j        (S = Stirling de segunda especie).
#      Sin alpha ni Q(zeta_n): solo Bernoulli de chi y la serie l = log E.
#      Contraste EXACTO (en el cuerpo ciclotomico, sin reducir) frente al calculo directo con zeta.
# FECHA: 1-oct-2026
# Uso (docker, /work): sage -python simplificacion_A.py
from sage.all import (CyclotomicField, DirichletGroup, PowerSeriesRing, QQ, binomial, carmichael_lambda, factorial,
                      lcm, stirling_number2)


def E_serie(p, prec):
    R = PowerSeriesRing(QQ, 't', default_prec=prec)
    t = R.gen()
    h = sum(t ** k / factorial(k) for k in range(p))
    hm = sum((-t) ** k / factorial(k) for k in range(p))
    return (h / hm).sqrt()


def compara(n, p, jmax):
    L = lcm(lcm(carmichael_lambda(n), n), 2)
    F = CyclotomicField(L)
    z = F.gen() ** (L // n)
    prec = jmax + 2
    E = E_serie(p, prec)
    ell = (E - 1).log() if False else E.log()
    R = PowerSeriesRing(F, 't', default_prec=prec)
    Et = R(E.list()).add_bigoh(prec)
    potl = [R(1).add_bigoh(prec)]
    for _ in range(jmax):
        potl.append(potl[-1] * R(ell.list()).add_bigoh(prec))
    res = []
    for chi in DirichletGroup(n, F):
        if chi(-1) != -1 or chi.conductor() != n:
            continue
        tau = sum(chi(c) ** -1 * z ** c for c in range(1, n) if chi(c) != 0)
        beta = [2 * (1 - 2 ** j * chi(2)) * chi.bernoulli(j + 1) / (j + 1) for j in range(jmax + 1)]
        # directo
        S = []
        for r in range(jmax + 1):
            w = Et ** r
            S.append(sum(chi(c) ** -1 * (2 * z ** c * w / ((z ** c * w) ** 2 - 1)) for c in range(1, n)
                         if chi(c) != 0))
        malos = 0
        tot = 0
        for i in range(jmax + 1):
            Ti = sum((-1) ** (i - r) * binomial(i, r) * S[r] for r in range(i + 1))
            for J in range(i, jmax + 1):
                directo = Ti[J]
                formula = tau * sum(beta[j] * factorial(i) / factorial(j) * stirling_number2(j, i) * potl[j][J]
                                    for j in range(i, J + 1))
                tot += 1
                malos += (directo != formula)
        res.append((chi.order(), tot, malos))
    return ell, res


for n, p in [(3, 5), (3, 7), (5, 11), (7, 5), (19, 7)]:
    jmax = 2 * (p - 1)
    ell, res = compara(n, p, jmax)
    print("n=%d p=%d: l(t) = %s" % (n, p, ell.add_bigoh(p + 3)), flush=True)
    for orden, tot, malos in res:
        print("   chi orden %d: %d entradas A_{J,i} (J <= %d), discrepancias exactas: %d" % (orden, tot, jmax, malos),
              flush=True)
