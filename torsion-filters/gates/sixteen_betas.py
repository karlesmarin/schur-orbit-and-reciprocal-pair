# -*- coding: utf-8 -*-
# LAS 16 FORMAS DE  t=6, r=2  CON  Delta != 0.   15 de agosto de 2026.
#
# POR QUE EXISTE ESTE FICHERO.  La tabla de la vuelta 11
#
#     t=6 r=2   16 fallos   Delta: {(1,1):16}     omega_2 en 100 %
#
# se mando en una carta y su GUION NO SE GUARDO: quedo la salida, no el instrumento.  Es el fallo
# inverso al de [[save-the-outputs-not-just-the-scripts]] y cuesta lo mismo.  Aqui se reconstruye,
# se deja escrito, y se vuelca la lista de betas para que la tabla de paredes (wall_table.sage)
# trabaje sobre EL MISMO conjunto y no sobre uno parecido.
#
# QUE ES UN "FALLO".  Poblacion critica de dominant_vector.py N4:  C = tau, S\g_com simetrico,
# g_com ASIMETRICO, Phi != 0.  Un fallo es una forma donde TODO v maximal en dominancia entre las
# transversales tiene coeficiente 0 -- o sea, el candidato de arriba de NUESTRA presentacion de
# Laplace se cancela entero y hay que bajar.  Delta = v_maximal - v_primero_no_nulo.
#
# OJO, Y VA EN LA SALIDA.  coef(v) suma sobre las r! permutaciones de v: es el numero contaminado
# de la vuelta 09/11 ([[a-number-i-hand-over-carries-my-instrument]]).  Para DETECTAR si la orbita
# entera se cancela es lo correcto -- una orbita se anula o no como orbita.  Para dar un valor no
# lo es.  Aqui solo se usa para detectar, y no se exporta ningun coeficiente.
#
# CONTROLES
#   C0  se imprime n de la poblacion critica ANTES de filtrar, y cuantas se caen en cada filtro.
#   C1  ACEPTACION: el recuento tiene que dar  16 fallos  y  Delta = (1,1) en los 16, que es lo que
#       se mando en la vuelta 11.  Si no coincide, la reconstruccion NO es la de la carta y hay que
#       decirlo antes de usarla.
#   C2  se vuelca tambien t=4,r=2 (16 fallos, Delta variado) como CONTROL: si la reconstruccion
#       reprodujera los 16 de t=6 pero no los de t=4, seria casualidad y no reconstruccion.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python sixteen_betas.py

import itertools, json
from collections import defaultdict

from second_stratum import setup, all_transversals
from collision_graph import atoms, fibras, S_of
from peel_zero import betas, occupied, phi_zero
from dominant_vector import v_of, domina, gcom_de


def fallos_de(t, r, W):
    """Devuelve (n_critica, [(beta, v_max, v_primero, Delta, profundidad)])."""
    n = 0
    out = []
    for b in betas(t, r, W):
        if not occupied(b, t):
            continue
        st = setup(b, t)
        if st is None:
            continue
        cl, E, Cd = st
        if not E:
            continue
        C, tau, gc = gcom_de(b, t, r)
        if C is None:
            continue
        S = S_of(b, t)
        resto = [x for x in S if x not in gc]
        if C != tau or set(C - x for x in resto) != set(resto):
            continue
        if set(C - x for x in gc) == set(gc):        # g_com simetrico: fuera
            continue
        if phi_zero(b, t, r):                        # Phi = 0: fuera
            continue
        n += 1
        tr = all_transversals(b, cl, r, t)
        vs = defaultdict(int)
        for (_, T, w, _) in tr:
            vs[v_of(T, r)] += w
        at = atoms(b, t, r)
        f = fibras(at)

        def coef(v):
            return sum(sum(c for c, _ in f[p]) for p in itertools.permutations(v) if p in f)

        maxi = [v for v in vs if not any(u != v and domina(u, v) for u in vs)]
        if any(coef(v) != 0 for v in maxi):
            continue                                  # el certificado de arriba SI vale
        niveles = sorted(vs, key=lambda v: (-sum(v), tuple(-x for x in v)))
        vmax = niveles[0]
        for k, v in enumerate(niveles):
            if coef(v):
                out.append((tuple(b), vmax, v, tuple(vmax[i] - v[i] for i in range(r)), k))
                break
        else:
            out.append((tuple(b), vmax, None, None, -1))
    return n, out


print("=" * 108)
print("LAS FORMAS CON Delta != 0  --  reconstruccion del recuento de la vuelta 11")
print("=" * 108)
print("")

DUMP = {}
for (t, r, W) in [(6, 2, 15), (4, 2, 15)]:
    n, fal = fallos_de(t, r, W)
    dist = defaultdict(int)
    prof = defaultdict(int)
    for (_, _, _, D, k) in fal:
        dist[D] += 1
        prof[k] += 1
    print("  t=%d r=%d W=%d :  poblacion critica n = %d   |   fallos = %d   |   Delta: %s   |   profundidad: %s"
          % (t, r, W, n, len(fal), dict(dist), dict(prof)))
    DUMP["t%d_r%d" % (t, r)] = [{"beta": list(b), "v_max": list(vm), "v_1": list(v1) if v1 else None,
                                 "Delta": list(D) if D else None, "prof": k}
                                for (b, vm, v1, D, k) in fal]

print("")
print("  C1  ACEPTACION contra lo que se mando en la vuelta 11:")
c1 = DUMP["t6_r2"]
d6 = defaultdict(int)
for x in c1:
    d6[tuple(x["Delta"]) if x["Delta"] else None] += 1
ok = (len(c1) == 16 and d6.get((1, 1)) == 16)
print("      t=6 r=2 -> se esperaban 16 fallos con Delta=(1,1) en los 16.  Salen %d, y Delta %s.  %s"
      % (len(c1), dict(d6), "COINCIDE" if ok else "*** NO COINCIDE: no es la misma poblacion ***"))
c2 = DUMP["t4_r2"]
print("      t=4 r=2 -> se esperaban 16 fallos con Delta variado (control).  Salen %d." % len(c2))

json.dump(DUMP, open("sixteen_betas.json", "w"), indent=1)
print("")
print("  volcado en sixteen_betas.json")
print("")
print("=" * 108)
print("DONE" if ok else "DONE (C1 FALLA -- no usar la lista sin resolver esto)")
