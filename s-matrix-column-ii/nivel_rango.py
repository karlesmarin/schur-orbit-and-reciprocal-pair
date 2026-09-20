# -*- coding: utf-8 -*-
r"""nivel_rango.py -- que dicen los datos de Delta(n, p), mirados en serio.

Cuatro preguntas a los datos, cada una con su control:

  (A) DESCOMPOSICION.  Delta = b - N0_orb + exceso, con b = orbitas, N0_orb = orbitas con delta = 0
      y exceso = sum (delta - 1) sobre las de delta >= 2.  Se mide cada pieza por separado, porque
      la primera tiene forma cerrada demostrada (P3a) y solo el exceso esta abierto.

  (B) NIVEL-RANGO.  El semigrupo se conserva bajo T -> -(F_p \ T), que lleva n a p-n.  Entonces el
      reparto de delta en n y en p-n debe coincidir ORBITA A ORBITA, y Delta(n,p)/b(n,p) tambien.
      Se comprueba en vez de suponerlo.

  (C) TERMINO DE ERROR.  A(n,p) = #{e_1 = e_2 = 0} frente a C1(n,p)/p, normalizado por sqrt(p):
      si el recuento es de tipo Weil, (A - C1/p)/(C1/p * sqrt(p) / sqrt(p)) queda acotado.  Se
      imprime la desviacion en unidades de sqrt(C1/p), que es la escala natural.

  (D) EL EXCESO, DE DONDE SALE.  Reparto de las orbitas con delta >= 2 segun QUE momentos se anulan
      (solo M_2, solo M_3, ambos, mas de dos), que es lo que decide el semigrupo.

    python nivel_rango.py [PMAX] [NMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from math import sqrt

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import centrados, es_primo, estabilizador, momento, semigrupo, s_exacto   # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 23
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 0          # 0 = todos los n hasta p-2


def datos(p, n):
    """por orbita: (delta, h, indices r <= n+2 con M_r = 0)."""
    out = []
    for T in centrados(p, n):
        d = semigrupo(s_exacto(T, p))["genero"]
        h = len(estabilizador(T, p))
        nulos = tuple(r for r in range(2, n + 3) if momento(T, r, p) % p == 0)
        out.append((d, h, nulos))
    return out


print("=" * 104)
print("(A) DESCOMPOSICION  Delta = b - N0 + exceso        (B) NIVEL-RANGO n <-> p-n")
print("=" * 104)
print("%-4s %-3s %-7s %-7s %-7s %-8s %-9s %s" % ("p", "n", "b", "N0", "exceso", "Delta", "Delta/b", "reparto de delta"))
tabla = {}
for p in range(7, PMAX + 1):
    if not es_primo(p):
        continue
    nmax = NMAX or (p - 2)
    for n in range(2, nmax + 1):
        if n > p - 2:
            break
        d = datos(p, n)
        b = len(d)
        N0 = sum(1 for x in d if x[0] == 0)
        Delta = sum(x[0] for x in d)
        exceso = sum(x[0] - 1 for x in d if x[0] >= 2)
        rep = {}
        for x in d:
            rep[x[0]] = rep.get(x[0], 0) + 1
        tabla[(p, n)] = (b, N0, exceso, Delta, rep, d)
        print("%-4d %-3d %-7d %-7d %-7d %-8d %-9.4f %s"
              % (p, n, b, N0, exceso, Delta, Delta / b if b else 0,
                 " ".join("%d:%d" % (k, rep[k]) for k in sorted(rep))))
    # control de nivel-rango dentro de este primo
    for n in range(2, (p - 2) // 2 + 1):
        if (p, n) in tabla and (p, p - n) in tabla:
            r1, r2 = tabla[(p, n)][4], tabla[(p, p - n)][4]
            print("     nivel-rango n=%d <-> n=%d : repartos %s" % (n, p - n, "IGUALES" if r1 == r2 else
                  "DISTINTOS  %s  vs  %s" % (r1, r2)))
    sys.stdout.flush()

print("")
print("=" * 104)
print("(D) DE DONDE SALE EL EXCESO: que momentos se anulan en las orbitas con delta >= 2")
print("=" * 104)
patrones = {}
for (p, n), (b, N0, exceso, Delta, rep, d) in sorted(tabla.items()):
    for delta, h, nulos in d:
        if delta >= 2:
            clave = (nulos, delta, h >= 2)
            patrones[clave] = patrones.get(clave, 0) + 1
for clave in sorted(patrones, key=lambda k: (-patrones[k], k)):
    nulos, delta, hmas = clave
    print("   M_r = 0 en r = %-14s delta = %-3d %-10s : %d orbitas"
          % (str(list(nulos)), delta, "h >= 2" if hmas else "h = 1", patrones[clave]))
