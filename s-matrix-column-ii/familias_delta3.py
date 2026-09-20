# -*- coding: utf-8 -*-
r"""familias_delta3.py -- buscar el fallo en las dos zonas senaladas.

Aviso recibido el 20-sep: hay dos problemas, en TRINOMIOS y en SEMIGRUPOS DE DEFECTO 3.  El
censo de la figura 2 y el recuento de semisistemas se reproducen exactos hasta 199, o sea que
lo que falla no son esos numeros.  Esta sonda ataca las afirmaciones que cargan peso en esas
dos zonas.

LA PRIMERA VERSION DE ESTA SONDA DIO DOS FALSOS POSITIVOS, y merece quedar escrito porque los
dos son el mismo error: PEDIR MENOS QUE LO QUE SE COMPRUEBA.

  * Llamaba "trinomio" a e_1 = ... = e_{n-2} = 0 y se dejaba el "ab != 0" que la Parte I SI
    pone.  Con esa definicion floja salian ocho primos con contraejemplo; todos eran cosets de
    un subgrupo (a = 0) o llevaban el cero (b = 0).  La Parte I tiene razon; el que pedia de
    menos era yo.
  * Llamaba "semisistema" solo a los T sin cero, cuando un espectro del censo puede llevar el
    cero adjunto y sus momentos pares se anulan igual.

Asi que aqui las dos nociones van con TODAS sus condiciones, y ademas se mide cada implicacion
por separado: una equivalencia que falla siempre en el mismo sentido dice donde esta la
hipotesis perdida, y una que falla en los dos dice que el enunciado esta mal.

  T1  "delta maximo  <=>  trinomio X^n + aX + b con ab != 0", para h = 1.  Las dos direcciones.
  T2  la misma equivalencia SIN el "ab != 0": es la forma en que la Nota II cita a la Parte I.

  D1  "se anula todo momento par  =>  semisistema (quiza con el cero adjunto)".
  D2  "las orbitas mas profundas SON los semisistemas": delta >= 3 que no lo sean.
  D3  "semisistema con h = 1 y M_3 M_5 M_7 != 0  =>  S_T = <3,5,7>, delta = 3".

    python familias_delta3.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from itertools import combinations
from math import comb

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import (centrados, es_primo, estabilizador, momento,           # noqa: E402
                 s_exacto, semigrupo)

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 19


def elementales(T, p):
    """e_1, ..., e_n de prod (X - t), con signo ya aplicado: coef de X^{n-i} es (-1)^i e_i."""
    e = [1] + [0] * len(T)
    for t in T:
        for i in range(len(T), 0, -1):
            e[i] = (e[i] + t * e[i - 1]) % p
    return e[1:]


def es_trinomio(T, p, con_ab=True):
    """prod (X - t) = X^n + aX + b.  Con con_ab, se exige ademas a != 0 y b != 0."""
    n = len(T)
    e = elementales(T, p)
    if any(e[i] % p for i in range(0, n - 2)):      # e_1 ... e_{n-2} han de anularse
        return False
    return not con_ab or (e[n - 2] % p != 0 and e[n - 1] % p != 0)


def semisistema(T, p, cero_permitido=True):
    """Uno de cada pareja {a,-a}; si cero_permitido, se admite el 0 adjunto."""
    S = set(T)
    if 0 in S:
        if not cero_permitido:
            return False
        S.discard(0)
    if len(S) != (p - 1) // 2:
        return False
    return all((a in S) != ((p - a) in S) for a in range(1, (p + 1) // 2))


def datos(T, p):
    h = len(estabilizador(T, p))
    return h, (p - 1) // h, semigrupo(s_exacto(T, p))


print("=" * 88)
print("T1/T2 -- el extremo: ¿defecto maximo <=> trinomio?  Las dos direcciones, por separado.")
print("=" * 88)
print("%-5s %-4s %-9s %-26s %-26s" % ("p", "n", "orb h=1", "CON ab!=0 (Parte I)",
                                      "SIN ab!=0 (como lo cita P3b)"))
t1_mal, t2_mal = [], []
for p in range(11, max(PMAX, 47) + 1):
    if not es_primo(p):
        continue
    for n in (4,):
        dmax = comb(n - 1, 2)
        libres = 0
        a1 = b1 = a2 = b2 = 0     # a: trinomio sin delta max ; b: delta max sin ser trinomio
        for T in centrados(p, n):
            try:
                h, e, S = datos(T, p)
            except ValueError:
                continue
            if h != 1:
                continue
            libres += 1
            top = S["genero"] == dmax
            for con_ab, (na, nb) in ((True, (0, 0)), (False, (0, 0))):
                pass
            tri_con = es_trinomio(T, p, True)
            tri_sin = es_trinomio(T, p, False)
            a1 += tri_con and not top
            b1 += top and not tri_con
            a2 += tri_sin and not top
            b2 += top and not tri_sin
        if a1 or b1:
            t1_mal.append((p, n, a1, b1))
        if a2 or b2:
            t2_mal.append((p, n, a2, b2))
        print("%-5d %-4d %-9d %-26s %-26s"
              % (p, n, libres,
                 "ok" if not (a1 or b1) else "tri sin max %d / max sin tri %d" % (a1, b1),
                 "ok" if not (a2 or b2) else "tri sin max %d / max sin tri %d" % (a2, b2)))
        sys.stdout.flush()

print("")
print("=" * 88)
print("D1/D2/D3 -- el centro")
print("=" * 88)
d1, d2, d3 = [], [], []
n_semis = n_prof = 0
for p in range(7, PMAX + 1):
    if not es_primo(p):
        continue
    for n in range(3, p - 1):
        for T in centrados(p, n):
            try:
                h, e, S = datos(T, p)
            except ValueError:
                continue
            delta = S["genero"]
            semi = semisistema(T, p)
            n_semis += semi
            # RANGO CORREGIDO.  La Proposicion 4 pide r par con 2 <= r < p-1, NO r < e.  Con el rango
            # malo, T = {1,3,9} en p=13 (donde e=4) pasaba el filtro mirando solo M_2, y se
            # publico como contraejemplo siendo falso: su M_6 = 3.
            if all(momento(T, r, p) == 0 for r in range(2, p - 1, 2)) and not semi:
                d1.append((p, n, T, delta))
            if delta >= 3:
                n_prof += 1
                if not semi:
                    d2.append((p, n, T, delta, tuple(S["huecos"])))
            if semi and h == 1 and all(momento(T, r, p) for r in (3, 5, 7)):
                if tuple(S["huecos"]) != (1, 2, 4):
                    d3.append((p, n, T, h, tuple(S["huecos"])))

print("  semisistemas centrados (con el cero admitido): %d ; orbitas con delta >= 3: %d"
      % (n_semis, n_prof))
print("")
print("  D1  momentos pares nulos y NO semisistema : %d" % len(d1))
for r in d1[:6]:
    print("        p=%-3d n=%-2d T=%-30s delta=%d" % r)
print("")
print("  D2  delta >= 3 y NO semisistema           : %d de %d" % (len(d2), n_prof))
for r in d2[:8]:
    print("        p=%-3d n=%-2d T=%-30s delta=%d huecos=%s" % r)
print("")
print("  D3  semisistema h=1, M3M5M7!=0, S_T != <3,5,7> : %d" % len(d3))
for r in d3[:6]:
    print("        p=%-3d n=%-2d T=%-30s h=%d huecos=%s" % r)

print("")
print("=" * 88)
print("VEREDICTO")
print("=" * 88)
print("  T1 con ab != 0, como lo enuncia la Parte I : %s"
      % ("FALLA en %d casos" % len(t1_mal) if t1_mal else "SE SOSTIENE"))
print("  T2 sin ab != 0, como lo cita la Nota II    : %s"
      % ("FALLA en %d (p, n) --- la Nota II perdio esa hipotesis al citar" % len(t2_mal)
         if t2_mal else "se sostiene"))
print("  D1 momentos pares nulos => semisistema     : %s"
      % ("FALLA en %d" % len(d1) if d1 else "se sostiene con el cero admitido"))
print("  D2 las mas profundas son los semisistemas  : %s"
      % ("FALSO: %d de %d orbitas de delta>=3 no lo son" % (len(d2), n_prof) if d2
         else "se sostiene en el rango"))
print("  D3 <3,5,7> con h=1 y M3M5M7 != 0           : %s"
      % ("FALLA en %d" % len(d3) if d3 else "SE SOSTIENE"))
