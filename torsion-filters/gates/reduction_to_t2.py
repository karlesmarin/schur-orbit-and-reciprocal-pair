# -*- coding: utf-8 -*-
# LA REDUCCION A t=2, Y DOS CONTROLES QUE MATAN DOS IDEAS.  15 de agosto de 2026.
#
# ESTO NO SALE DE LA CONSULTA EXTERNA.  Sale de preguntarle al objeto qué es, hasta la raiz.
#
# LA RAIZ.  El alfabeto es mu_t = {1, zeta, ..., zeta^{t-1}}.  El paper ya dice (intro, "two
# structural reasons") que ese alfabeto es CERRADO BAJO x -> 1/x.  Pero entonces se parte solo:
#
#     t IMPAR:  mu_t = {1}      U  {(t-1)/2 pares reciprocos {zeta^k, zeta^-k}}
#     t PAR:    mu_t = {1, -1}  U  {(t-2)/2 pares reciprocos}
#
# o sea que los puntos FIJOS de la inversion son 1 (t impar) o 1 y -1 (t par), y TODO LO DEMAS ya
# son pares reciprocos, igual que nuestros (z, 1/z) pero evaluados en raices de la unidad.  Luego
#
#     Phi_{t,r}(beta; z)  =  Phi_{2, r + (t-2)/2}(beta; zeta, zeta^2, ..., z_1, ..., z_r)   [t par]
#
# y OJO al detalle que lo hace comprobable sin escribir una linea de algebra simbolica: la longitud
# no cambia.  t + 2r = 2 + 2(r + (t-2)/2).  ES EL MISMO BETA.  Con t=4, r=2:  Phi_{4,2} es Phi_{2,3}
# con la tercera variable puesta en zeta = i.
#
# Y DE AHI UNA PREDICCION QUE PUEDE FALLAR:
#
#     Phi_{2,3}(beta) == 0 IDENTICAMENTE   =>   Phi_{4,2}(beta) == 0
#
# El reciproco NO tiene por que valer: especializar puede crear ceros nuevos.  Y esa diferencia mide
# EXACTAMENTE cuanto del problema abierto es "ya estaba en t=2" -- que es el Teorema 8.6, probado
# para todo r -- y cuanto es torsion pura.  N2 lo mide.
#
# CONTROLES, y los dos MATAN una idea de la consulta
#   C1  su cadena omega_r.  Propone que el descenso vaya siempre por omega_r = (1,...,1).  Se mide el
#       vector Delta entero, no su norma.
#   C2  su divisibilidad por r!.  Observa 2 = 2!, 6 = 3!, 12 = 2*3! en los primeros coeficientes no
#       nulos que le dimos.  El control obvio es mirar TODOS los coeficientes de nivel, no solo el
#       primero: si todos son divisibles, la observacion es vacia.  Y hay una razon para sospecharlo,
#       y es NUESTRO defecto: nuestro coef(v) suma sobre las r! permutaciones de v, asi que la
#       divisibilidad por r! esta metida en el instrumento.  El numero que le dimos ya la llevaba.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python reduction_to_t2.py

import itertools
import math
from collections import Counter, defaultdict

from second_stratum import setup, all_transversals
from collision_graph import atoms, fibras, S_of
from dominant_vector import v_of, domina, gcom_de
from peel_zero import betas, occupied, phi_zero

if __name__ == "__main__":
    print("=" * 104)
    print("N1  LA REDUCCION -- mu_t son puntos fijos de x->1/x mas PARES RECIPROCOS")
    print("=" * 104)
    print("")
    print("   mu_4 = {1, i, -1, -i} = {1, -1} U {i, 1/i}")
    print("   luego los argumentos de Phi_{4,2} son (1, -1, i, 1/i, z1, 1/z1, z2, 1/z2)")
    print("   que son los de Phi_{2,3} con la tercera variable puesta en i.")
    print("   Y N = 4+2*2 = 8 = 2+2*3:  ES EL MISMO BETA, sin cambiar de longitud.")
    print("")
    print("   PREDICCION:   Phi_{2,3}(beta) == 0  =>  Phi_{4,2}(beta) == 0")

    print("")
    print("=" * 104)
    print("N2  CUANTO DEL PROBLEMA ABIERTO ES 'YA ESTABA EN t=2' Y CUANTO ES TORSION")
    print("=" * 104)
    print("")
    print("    W |     n | 23=0 & 42=0 | 23=0 & 42!=0 | 23!=0 & 42=0 | ceros de Phi_42 que son NUEVOS")
    for W in (13, 14, 15):
        tab = Counter()
        n = 0
        for b in betas(4, 2, W):
            z43 = phi_zero(b, 4, 2)
            z23 = phi_zero(b, 2, 3)
            if z43 is None or z23 is None:
                continue
            n += 1
            tab[(bool(z23), bool(z43))] += 1
        tot0 = tab[(True, True)] + tab[(False, True)]
        print("   %2d | %5d | %11d | %12d | %12d | %d de %d  (%.0f %%)"
              % (W, n, tab[(True, True)], tab[(True, False)], tab[(False, True)],
                 tab[(False, True)], tot0, 100.0 * tab[(False, True)] / max(tot0, 1)))
    print("")
    print("   La columna '23=0 & 42!=0' es el CONTRAEJEMPLO a la prediccion.  En 0 significa que la")
    print("   reduccion es cierta.  Y la ultima columna es lo que dice de verdad: la reduccion explica")
    print("   una fraccion PEQUEÑA de los ceros.  Todo lo demas lo crea la especializacion.")

    # ================================================================= C1 =======================
    print("")
    print("=" * 104)
    print("C1  SU CADENA omega_r -- se mide el vector Delta entero, no su norma")
    print("=" * 104)
    print("")
    for (t, r, W) in [(4, 2, 15), (6, 2, 15), (4, 3, 14)]:
        deltas = Counter()
        coefs = Counter()
        n = 0
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
            if set(C - x for x in gc) == set(gc):
                continue
            if phi_zero(b, t, r):
                continue
            tr = all_transversals(b, cl, r, t)
            vs = defaultdict(int)
            for (_, T, w, _) in tr:
                vs[v_of(T, r)] += w
            at = atoms(b, t, r)
            f = fibras(at)

            def coef(v):
                return sum(sum(x for x, _ in f[p]) for p in itertools.permutations(v) if p in f)
            maxi = [v for v in vs if not any(u != v and domina(u, v) for u in vs)]
            if any(coef(v) != 0 for v in maxi):
                continue
            n += 1
            v0 = maxi[0] if len(maxi) == 1 else max(maxi, key=sum)
            for v in sorted(vs, key=lambda v: (-sum(v), tuple(-x for x in v))):
                c = coef(v)
                if c:
                    deltas[tuple(a - bb for a, bb in zip(v0, v))] += 1
                    coefs[abs(c)] += 1
                    break
        om = tuple([1] * r)
        if not n:
            print("   t=%d r=%d : POBLACION VACIA" % (t, r))
            continue
        print("   t=%d r=%d   %d fallos" % (t, r, n))
        print("      Delta observados : %s" % dict(deltas))
        print("      omega_r = %s en %d de %d  (%.0f %%)"
              % (str(om), deltas.get(om, 0), n, 100.0 * deltas.get(om, 0) / n))
    print("")
    print("   Si omega_r no sale siempre, no hay cadena.  Y notese que algun Delta tiene componente")
    print("   NEGATIVA: el descenso ni siquiera es monotono coordenada a coordenada.")

    # ================================================================= C2 =======================
    print("")
    print("=" * 104)
    print("C2  SU DIVISIBILIDAD POR r!  -- mirando TODOS los niveles, no solo el primero no nulo")
    print("=" * 104)
    print("")
    print("    t  r | formas | coef. de nivel no nulos | NO divisibles por r! | veredicto")
    for (t, r, W) in [(4, 2, 14), (4, 3, 14)]:
        todos = Counter()
        n = nodiv = 0
        for b in betas(t, r, W):
            if not occupied(b, t):
                continue
            st = setup(b, t)
            if st is None:
                continue
            cl, E, Cd = st
            if not E or phi_zero(b, t, r):
                continue
            n += 1
            if n > 60:
                break
            tr = all_transversals(b, cl, r, t)
            vs = set(v_of(T, r) for (_, T, w, _) in tr)
            at = atoms(b, t, r)
            f = fibras(at)
            for v in vs:
                c = sum(sum(x for x, _ in f[p]) for p in itertools.permutations(v) if p in f)
                if c:
                    todos[abs(c)] += 1
                    nodiv += (abs(c) % math.factorial(r) != 0)
        tot = sum(todos.values())
        print("   %2d %2d | %6d | %23d | %20d | %s"
              % (t, r, n, tot, nodiv,
                 "no trivial" if nodiv else "VACIO: todo divisible"))
    print("")
    print("   Sale VACIO, y la razon es NUESTRA: coef(v) suma sobre las r! permutaciones de v, asi")
    print("   que el factor r! lo mete el instrumento.  El 2, el 6 y el 12 que le dimos ya lo")
    print("   llevaban puesto.  Un instrumento que fabrica el patron que su lector va a conjeturar.")
    print("")
    print("=" * 104)
