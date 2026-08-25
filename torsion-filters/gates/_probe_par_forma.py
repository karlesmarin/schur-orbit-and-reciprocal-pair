# -*- coding: utf-8 -*-
# SONDA: .QUE FORMA TIENE UNA COLUMNA LIBRE PAR?   19 de agosto de 2026.
#
# NO decide nada.  Inventaria la forma antes de enunciar, para no escribir un lema adivinado.  Lo que
# el nucleo par predice es "semirrecta(+m) menos semirrecta(-m)", o sea a lo sumo DOS bloques de
# clase, media recta dentro de cada uno, y signos opuestos.  Aqui se cuenta lo que hay de verdad.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_par_forma.py > _probe_par_forma_OUT.txt

import itertools
import sys
from collections import Counter

from divided_differences import plegar

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CASOS = [(4, 2, 4), (6, 2, 3), (8, 2, 3), (10, 2, 2), (6, 3, 2), (8, 3, 2)]


def col_libre_par(V, x, t, m):
    out = []
    for Vi in V:
        tot = 0
        for e in (1, -1):
            d = e * Vi - x
            if d >= m and (d - m) % t == 0:
                tot += e
            if d >= m + 2 and (d - m - 2) % t == 0:
                tot -= e
        out.append(tot)
    return out


valores = Counter()
nbloques = Counter()
forma = Counter()
signos_cruz = Counter()

for (t, r, cota) in CASOS:
    m = t // 2 - 1
    R = m + r
    for Lam in itertools.product(range(cota + 1), repeat=R):
        if any(Lam[i] < Lam[i + 1] for i in range(R - 1)):
            continue
        V = [Lam[i] + R - i for i in range(R)]
        pl = [plegar(v, t) for v in V]
        # orden por (clase plegada, w = sigma V)
        info = sorted([(pl[i][0], pl[i][1] * V[i], i) for i in range(R)])
        orden = [z[2] for z in info]
        clase_de = {k: info[k][0] for k in range(R)}
        for x in range(1, max(V) + t + 1):
            col = col_libre_par(V, x, t, m)
            if not any(col):
                continue
            for v in col:
                valores[v] += 1
            A = [pl[i][1] * col[i] for i in orden]          # firmando filas por sigma, como el impar
            nz = [k for k, v in enumerate(A) if v]
            bl = sorted({clase_de[k] for k in nz})
            nbloques[len(bl)] += 1
            # forma dentro de cada bloque: .es prefijo, sufijo, o ninguno?
            for c in bl:
                idx = [k for k in range(R) if clase_de[k] == c]
                dentro = [k for k in idx if A[k]]
                pref = dentro == idx[:len(dentro)]
                suf = dentro == idx[len(idx) - len(dentro):]
                cons = dentro == list(range(dentro[0], dentro[-1] + 1))
                forma[("prefijo" if pref else "sufijo" if suf else
                       "consecutivo" if cons else "ROTO")] += 1
                if len({A[k] for k in dentro}) != 1:
                    forma["signo_no_constante_en_bloque"] += 1
            if len(bl) == 2:
                s = []
                for c in bl:
                    dentro = [k for k in range(R) if clase_de[k] == c and A[k]]
                    s.append(A[dentro[0]])
                signos_cruz[("opuestos" if s[0] * s[1] < 0 else "iguales")] += 1

print("=" * 100)
print("SONDA: LA FORMA DE UNA COLUMNA LIBRE PAR")
print("=" * 100)
print("")
print("  valores que toma una entrada        : %s" % dict(sorted(valores.items())))
print("  bloques de clase que toca la columna: %s" % dict(sorted(nbloques.items())))
print("  forma dentro de cada bloque tocado  : %s" % dict(forma))
print("  cuando toca DOS bloques, los signos : %s" % dict(signos_cruz))
print("")
print("=" * 100)
print("DONE")
