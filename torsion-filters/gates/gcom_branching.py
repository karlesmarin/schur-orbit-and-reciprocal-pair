# -*- coding: utf-8 -*-
# g_com EN EL LENGUAJE DEL BRANCHING.   16 de agosto de 2026.
#
# DE DONDE SALE.  Propuesta 4 de su reseña (vuelta 14), la unica de su hoja de ruta sin empezar: el
# Paper II no menciona g_com ni una vez, y para una continuacion eso es demasiado corte.  Pide
# traducir la Conjetura 8.43 del Paper I a un enunciado sobre los B_{eta,mu} tau_t(eta), y mostrar
# QUE tendria que ser cierto de ellos para que 8.43 siguiera.
#
# LA TRADUCCION, escrita antes de medir.
#
#   Paper I:   bajo ocupacion,   Phi_{t,r} = 0  =>  C - g_com = g_com.          (8.43)
#   Paper II:  Phi_{t,r} = 0    <=>  A_mu = 0 para todo mu   <=>   E^{(t)} a = 0,
#              o sea: el vector a del Paper I esta en el NUCLEO del operador de torsion.
#   Y con (H):  Phi = 0  <=>  A_{mu_max(beta)} = 0,   porque si Phi != 0 el peso superior sobrevive
#              con |A| = 1.  Luego la anulacion la decide UN SOLO coeficiente, el de arriba.
#
#   Traducida, 8.43 dice:  "si el compuesto mata el peso superior, entonces g_com es simetrica".
#
# LO QUE SE MIDE
#   G1  la tabla de contingencia sobre la poblacion ENTERA:  (Phi=0 o no) x (|G|) x (g_com vacia o
#       no) x (C - g_com = g_com o no).
#   G2  el CONTROL QUE PUEDE FALLAR, y es el que decide si 8.43 dice algo:  ¿es g_com simetrica
#       tambien en las formas que NO se anulan y tienen |G|=2?  Si lo fuera SIEMPRE, la conjetura no
#       seria sobre la anulacion sino una identidad, y habria que reformularla.
#   G3  el puente: entre las formas que se anulan, ¿se corresponde "g_com vacia" con la clase I
#       (heredada de t=2) y "g_com no vacia" con la clase III (por cancelacion)?  Esa es la
#       traduccion concreta que pide.
#   G4  el tamaño de g_com contra el numero de eta que sobreviven al filtro: si |g_com| midiera
#       algo del lado del branching, tendria que verse aqui.
#
# CONTROLES
#   C0  la anatomia se contrasta con phi_zero, que es la ruta de evaluacion independiente.
#   C1  n impreso en cada casilla; ninguna conclusion de una casilla con menos de 5 formas.
#   C2  g_com solo esta DEFINIDA cuando |G| = 2.  Las formas con |G| = 1 se cuentan aparte y no
#       entran en ninguna proporcion.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python gcom_branching.py

import json
import sys
from collections import Counter, defaultdict

from peel_gcom import anatomia, betas
from peel_zero import phi_zero

CFG = [(4, 2, 13), (6, 2, 13), (4, 2, 15)]

print("=" * 118)
print("g_com EN EL LENGUAJE DEL BRANCHING   --   la propuesta 4 de su reseña")
print("=" * 118)

RES = []
for (t, r, W) in CFG:
    tabla = Counter()
    tam_gcom = defaultdict(list)
    n_tot = n_ocup = 0
    ceros = 0
    ejemplos_asim = []
    for beta in betas(t, r, W):
        n_tot += 1
        a = anatomia(beta, t, r)
        if a is None:
            continue
        n_ocup += 1
        S, C, gcom, nG = a
        cero = phi_zero(beta, t, r) is True
        if cero:
            ceros += 1
        if nG != 2:
            tabla[("|G|=%d" % nG, "cero" if cero else "no cero", "g_com no definida")] += 1
            continue
        sim = (set(C - v for v in gcom) == set(gcom))
        clave = ("|G|=2",
                 "cero" if cero else "no cero",
                 "g_com vacia" if not gcom else ("simetrica" if sim else "ASIMETRICA"))
        tabla[clave] += 1
        tam_gcom[(cero, bool(gcom))].append(len(gcom))
        if cero and gcom and not sim:
            ejemplos_asim.append(beta)

    print("")
    print("-" * 118)
    print("  t=%d  r=%d  W=%d :  %d formas barridas, %d ocupadas, %d nulas" % (t, r, W, n_tot, n_ocup, ceros))
    print("-" * 118)
    print("  %-10s | %-8s | %-16s | n" % ("|G|", "Phi", "g_com"))
    print("  " + "-" * 60)
    for k in sorted(tabla):
        print("  %-10s | %-8s | %-16s | %d" % (k[0], k[1], k[2], tabla[k]))

    # G2: el control que puede fallar
    sim_nocero = tabla[("|G|=2", "no cero", "simetrica")]
    asim_nocero = tabla[("|G|=2", "no cero", "ASIMETRICA")]
    vac_nocero = tabla[("|G|=2", "no cero", "g_com vacia")]
    sim_cero = tabla[("|G|=2", "cero", "simetrica")]
    asim_cero = tabla[("|G|=2", "cero", "ASIMETRICA")]
    vac_cero = tabla[("|G|=2", "cero", "g_com vacia")]
    print("")
    print("  G2  EL CONTROL QUE PUEDE FALLAR -- g_com con |G|=2 y g_com NO vacia:")
    print("        entre las NULAS      : simetrica %d, ASIMETRICA %d" % (sim_cero, asim_cero))
    print("        entre las NO nulas   : simetrica %d, ASIMETRICA %d" % (sim_nocero, asim_nocero))
    if sim_nocero + asim_nocero == 0:
        print("        (no hay formas no nulas con g_com no vacia: el control esta VACIO aqui)")
    elif asim_nocero == 0:
        print("        -> g_com sale simetrica TAMBIEN sin anularse: 8.43 no estaria hablando de la")
        print("           anulacion.  Habria que reformularla.")
    else:
        print("        -> hay asimetricas entre las NO nulas: la hipotesis de anulacion SI hace")
        print("           trabajo, y 8.43 dice algo.  Proporcion: %d de %d asimetricas."
              % (asim_nocero, sim_nocero + asim_nocero))
    print("  G1  g_com VACIA: %d de %d nulas y %d de %d no nulas (con |G|=2)"
          % (vac_cero, vac_cero + sim_cero + asim_cero,
             vac_nocero, vac_nocero + sim_nocero + asim_nocero))
    if ejemplos_asim:
        print("  !!  CONTRAEJEMPLO A 8.43: %s" % str(ejemplos_asim[:3]))
    for k in sorted(tam_gcom):
        v = tam_gcom[k]
        if v:
            print("  G4  |g_com| en (%s, g_com no vacia=%s): min %d, max %d, media %.2f, n %d"
                  % ("nula" if k[0] else "no nula", k[1], min(v), max(v), sum(v) / len(v), len(v)))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "W": int(W), "n_barridas": int(n_tot),
                "n_ocupadas": int(n_ocup), "n_nulas": int(ceros),
                "tabla": {"|".join(k): int(v) for k, v in tabla.items()},
                "contraejemplos_843": [list(map(int, b)) for b in ejemplos_asim[:10]]})

print("")
print("=" * 118)
print("  LECTURA, escrita ANTES de correr:")
print("   * si aparece alguna ASIMETRICA entre las NULAS, la Conjetura 8.43 del Paper I es FALSA y")
print("     eso hay que mirarlo antes que ninguna otra cosa.")
print("   * si entre las NO nulas todas salen simetricas, 8.43 no usa la hipotesis y hay que")
print("     reformularla: seria una identidad sobre g_com, no un criterio de anulacion.")
print("   * si entre las NO nulas hay asimetricas, la hipotesis trabaja y la conjetura esta bien")
print("     planteada; entonces lo que falta es la traduccion al lado del branching.")
json.dump(RES, open("gcom_branching_RESULT.json", "w"), indent=1)
print("=" * 118)
print("DONE")
