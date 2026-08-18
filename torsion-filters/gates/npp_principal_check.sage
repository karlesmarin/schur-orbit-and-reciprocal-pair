# -*- coding: utf-8 -*-
# ¿ES NUESTRO ELEMENTO DE TORSION CONJUGADO A UN ELEMENTO PRINCIPAL?   16 de agosto de 2026.
#
# POR QUE HACE FALTA.  El Teorema 4.1 de NPP25 es sobre el elemento PRINCIPAL C_m = rho^v(e^{2pi i/m}),
# caracterizado por que TODAS las raices simples valen el mismo escalar en el.  Mirar solo la base
# estandar no basta: el grupo de Weyl mueve las bases, y g podria ser principal en OTRA camara.
#
# EL INVARIANTE CORRECTO.  Para un elemento del toro, W permuta las raices, luego el MULTICONJUNTO
#
#       M(g) = { alpha(g) : alpha en Phi+ }        (como numeros complejos exactos)
#
# es invariante por conjugacion en W salvo el paso Phi+ -> w(Phi+); tomando el multiconjunto sobre
# TODAS las raices (Phi, no solo Phi+) queda invariante de verdad.  Entonces:
#
#       g conjugado a C_m   =>   M(g) = M(C_m),     y  alpha(C_m) = zeta_m^{ht(alpha^v)}.
#
# Es condicion NECESARIA, que es lo que necesitamos: si los multiconjuntos difieren, g NO es
# principal, y el Teorema 4.1 no nos cubre.
#
# CONTROLES
#   C0  el caso de rango 1 (t=3 en B_1, t=4 en C_1) TIENE que salir principal: solo hay una raiz
#       simple.  Si no sale, el test esta mal escrito.
#   C1  se compara contra el principal del MISMO orden en el adjunto, y ademas contra todos los
#       m <= 2h, por si el orden que hay que mirar es otro.
#   C2  se imprime el multiconjunto cuando difiere, para poder mirarlo a mano.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage npp_principal_check.sage

import json
from collections import Counter

print("=" * 118)
print("¿ES PRINCIPAL NUESTRO ELEMENTO DE TORSION?   (necesario para saber si el Tma 4.1 de NPP25 aplica)")
print("=" * 118)
print("   t | grupo | h  | orden en el adjunto | principal en la base estandar | conjugado a ALGUN C_m")
print("   " + "-" * 112)

RES = []
for t in range(3, 12):
    if t % 2 == 0:
        typ, rk, expo = "C", (t - 2) // 2, None
    else:
        typ, rk, expo = "B", (t - 1) // 2, None
    if rk < 1:
        continue
    expo = [i + 1 for i in range(rk)]
    L = RootSystem("%s%d" % (typ, rk)).ambient_space()
    allroots = list(L.roots())
    simples = [L.simple_root(i) for i in L.index_set()]
    h = 2 * rk

    def val(a, e, mod):
        v = [int(u) for u in a.to_vector()]
        return sum(e[i] * v[i] for i in range(rk)) % mod

    # orden de g en el adjunto
    d_adj = None
    for k in range(1, 4 * t + 1):
        if all((k * val(a, expo, t)) % t == 0 for a in allroots):
            d_adj = k
            break

    # multiconjunto de valores de g sobre TODAS las raices, como fracciones k/t
    Mg = Counter(QQ(val(a, expo, t)) / t for a in allroots)
    principal_est = len(set(val(a, expo, t) for a in simples)) == 1

    # los principales C_m:  alpha(C_m) = zeta_m^{ht(alpha^v)}
    # ht(alpha^v) = <rho, alpha^v>
    rho = L.rho()
    rho_v = [QQ(u) for u in rho.to_vector()]

    def ht_coroot(a):
        av = [QQ(u) for u in a.to_vector()]
        num = sum(rho_v[i] * av[i] for i in range(rk))
        den = sum(av[i] * av[i] for i in range(rk))
        return 2 * num / den

    hits = []
    for m in range(2, 2 * h + 3):
        Mm = Counter((QQ(ht_coroot(a)) / m) - floor(QQ(ht_coroot(a)) / m) for a in allroots)
        if Mm == Mg:
            hits.append(int(m))
    print("   %2d | %s%-2d   | %2d | %19d | %-29s | %s"
          % (t, typ, rk, h, d_adj, "SI" if principal_est else "NO",
             str(hits) if hits else "NINGUNO  <-- fuera del Tma 4.1"))
    RES.append({"t": int(t), "tipo": typ, "rango": int(rk), "h": int(h), "d_adjunto": int(d_adj),
                "principal_base_estandar": bool(principal_est), "conjugado_a_C_m": hits})

print("")
print("=" * 118)
print("  LECTURA, escrita ANTES de correr: si la ultima columna sale NINGUNO para t >= 5, nuestro")
print("  elemento no es principal, el Teorema 4.1 no nos cubre, y estamos en la Question 8.1.")
json.dump(RES, open("npp_principal_check_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
