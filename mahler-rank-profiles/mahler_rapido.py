# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: version rapida del instrumento de Mahler (mahler_vueltas.py) calculada directamente en el cuerpo finito
#      F_{p^f} (f = orden de p mod lcm(n, exp de los caracteres)). Justificacion: todo es P-entero:
#        E(t) = raiz (termino 1) de h(t)/h(-t), h = exponencial truncada en grado p-1: coeficientes en Z_(p);
#        U(zeta^c w) = 2 zeta^c w / (zeta^{2c} w^2 - 1): el denominador en t = 0 es zeta^{2c} - 1, unidad sobre p (p no divide n);
#      asi que A_{j,i} mod P se obtiene reduciendo zeta y chi a F_{p^f}.
#      Salida: por caracter chi impar primitivo (p no divide phi(n)), presencia (P) / ausencia (A) en cada capa
#      P + j del sistema triangular de Mahler, para j en una lista; y W = suma de presencias (dimension de la capa)
#      para comparar con el orden medido.
#      Modo 'control': reproduce capas medidas del orden (vueltas_lambda_19_7_OUT, autosimilar_a3_OUT, cero_excepcional*
#      _OUT, vueltas_mazur_OUT, vueltas_masivo_OUT) y cuenta coincidencias.
# FECHA: 1-oct-2026
# Uso (docker, /work): sage -python mahler_rapido.py control
import sys

from sage.all import GF, PowerSeriesRing, QQ, binomial, factorial, lcm, matrix, euler_phi, primitive_root, \
    Integers, gcd


def E_mod_p(p, prec, F):
    R = PowerSeriesRing(QQ, 't', default_prec=prec)
    t = R.gen()
    h = sum(t ** k / factorial(k) for k in range(p))
    hm = sum((-t) ** k / factorial(k) for k in range(p))
    E = (h / hm).sqrt()
    return PowerSeriesRing(F, 't', default_prec=prec)([F(c) for c in E.list()]).add_bigoh(prec)


def caracteres_impares(n, F):
    """caracteres de (Z/n)^x con valores en F (ciclico o no: via DirichletGroup de Sage sobre F)."""
    from sage.all import DirichletGroup
    lam = Integers(n).unit_group_exponent()
    z = F.multiplicative_generator() ** ((F.order() - 1) // lam)
    G = DirichletGroup(n, F, zeta=z, zeta_order=lam)
    return [chi for chi in G if chi(-1) == -1]


def capas(n, p, js, a_max_i=None):
    """dict j -> (W, lista de (chi(2), chi(p), presencia)). a_max_i: si no es None, los funcionales i >= a_max_i no
    existen (piso a: i < p^a)."""
    lam = Integers(n).unit_group_exponent()
    N = lcm(n, lam)
    f = 1
    while (p ** f - 1) % N:
        f += 1
    F = GF(p ** f, 'g')
    z = F.multiplicative_generator() ** ((p ** f - 1) // n)  # raiz primitiva n-esima
    jmax = max(js)
    prec = jmax + 2
    Et = E_mod_p(p, prec, F)
    R = PowerSeriesRing(F, 't', default_prec=prec)
    chis = [c for c in caracteres_impares(n, F) if c.conductor() == n]
    d = euler_phi(n) // 2
    potE = [R(1).add_bigoh(prec)]
    for _ in range(jmax):
        potE.append(potE[-1] * Et)
    res = {j: [0, []] for j in js}
    # optimizacion (1-oct): las series U(zeta^c E^r) no dependen de chi; se calculan una vez
    unidades = [c for c in range(1, n) if gcd(c, n) == 1]
    Ucr = {(c, r): 2 * z ** c * potE[r] / ((z ** c * potE[r]) ** 2 - 1) for c in unidades for r in range(jmax + 1)}
    for chi in chis:
        S = []
        for r in range(jmax + 1):
            S.append(sum(chi(c) ** -1 * Ucr[(c, r)] for c in unidades))
        A = [[F(0)] * (jmax + 1) for _ in range(jmax + 1)]
        for i in range(jmax + 1):
            Ti = sum((-1) ** (i - r) * binomial(i, r) * S[r] for r in range(i + 1))
            for jj in range(i, jmax + 1):
                A[jj][i] = Ti[jj]
        for j in js:
            imax = j + 1 if a_max_i is None else min(j + 1, a_max_i)
            filas = [A[i0][:imax] for i0 in range(j)]
            r0 = matrix(F, filas).rank() if filas else 0
            r1 = matrix(F, filas + [A[j][:imax]]).rank()
            pres = r1 > r0
            res[j][0] += 1 if pres else 0
            res[j][1].append((str(chi(2)), str(chi(p)), 'P' if pres else 'A'))
    # caracteres no primitivos: se omiten (se contaran como presentes si su capa del piso 1 lo esta; aqui NO se
    # cuentan). W reportado = presencias de los primitivos + (d - #primitivos) (supuestos presentes).
    for j in js:
        res[j][0] += d - len(chis)
    return res


CONTROL = [  # (n, p, a, {capa: W medido})  -- capas medidas en el orden real
    (19, 7, 2, {55: 6, 61: 7, 67: 8, 73: 8, 79: 9}),
    (19, 7, 3, {349: 6, 355: 7, 361: 8, 367: 8, 373: 9, 379: 9, 385: 9, 391: 6, 397: 9, 403: 9, 409: 9, 415: 9}),
    (19, 11, 2, {131: 6, 141: 8}),
    (5, 11, 2, {131: 0, 141: 2}),
    (5, 19, 2, {379: 2}),
    (23, 5, 2, {29: 10, 33: 11}),
    (47, 5, 2, {29: 22, 33: 22, 37: 23}),
    (19, 5, 2, {29: 8, 33: 9}),
    (41, 11, 2, {131: 18, 141: 20}),
]

if __name__ == '__main__' and sys.argv[1] in ('control', 'resto'):
    if sys.argv[1] == 'resto':
        CONTROL = [c for c in CONTROL if (c[0], c[1]) in ((47, 5), (19, 5), (41, 11))]
    tot = ok = 0
    for n, p, a, med in CONTROL:
        P = p ** a
        js = sorted(c - P for c in med)
        res = capas(n, p, js, a_max_i=p ** a)
        for c, w in sorted(med.items()):
            Wp = res[c - P][0]
            tot += 1
            ok += (Wp == w)
            print("n=%d p=%d a=%d capa %d: Mahler W=%d, orden W=%d  %s" % (n, p, a, c, Wp, w,
                                                                       'OK' if Wp == w else 'FALLA'), flush=True)
    print("CONTROL: %d de %d" % (ok, tot))
