# -*- coding: utf-8 -*-
# EL ESTRATO DE ABAJO ES JUSTO LO QUE VE g_com.  Hipotesis del 12 de agosto de 2026.
#
# EL HUECO, tal como quedo escrito.  Con [Phi]_top = 0 se prueba (teorema del 12):
#       |G| = 2,  T_B = C - T_A,  sigma_C(K) = K,  x1+x2 = y1+y2 = C,
#       S  =  K  |_|  g_com  |_|  {x1,x2,y1,y2}
# y por tanto  sigma_C(S) = S  <=>  sigma_C(g_com) = g_com.  En t = 2 se tiene g_com = 0 y por eso
# el caso se cierra; en t >= 4 NO, y g_com es exactamente la parte de S que los dos maximizadores
# retiran EN COMUN, o sea la parte que el estrato de arriba no mira.
#
# LA HIPOTESIS QUE SE TESTA: el estrato de ABAJO si la mira.  Tres lecturas, de la mas debil a la
# mas fuerte, y cada una con su denominador:
#
#   H-A  (poblacion)   top = 0 y sigma_C(g_com) != g_com  =>  [Phi]_bot != 0.
#                      Es literalmente el enunciado que falta: con el, [Phi]_top = 0 y [Phi]_bot = 0
#                      fuerzan sigma_C(S) = S.
#   H-B  (mecanismo)   la VARIACION de G_bot toca las clases de g_com donde falla la simetria.
#                      Precisando: moved_bot = clases en que difieren dos minimizadores;
#                      asym_cls   = clases de los v in g_com con C-v no en g_com.
#                      Se afirma  moved_bot & asym_cls != vacio.
#   H-C  (localizacion) la clase de diferencias D que SOBREVIVE (suma con signo != 0) tiene a sus
#                      miembros difiriendo en una clase de g_com, no en la pareja empatada.
#
# CONTROLES QUE PUEDEN FALLAR, y por que estan puestos
#   A0  ACEPTACION, fatal.  [Phi]_bot calculado por clases de D (recuento con signo) tiene que
#       coincidir con el desarrollo monomio a monomio de sum_T w(T) [A(T)]_bot.  Sin esto, el
#       "recuento con signo" no es el estrato de abajo y todo lo demas es ruido.
#   A1  ACEPTACION, fatal.  [Phi]_top se calcula aqui por EXPANSION, no por el atajo INV+signo de
#       second_stratum.py.  Se comparan los dos: el atajo no resuelve el +- de P(T_A) = +- P(T_B).
#   A2  sigma_C conserva degmin y el multiconjunto D exactamente (una linea de prueba: las
#       diferencias de consecutivos se invierten en bloque).  Debe dar 0 fallos.
#   A3  la descomposicion  S = K |_| g_com |_| {x1,x2,y1,y2}  con |g_com| = e-2.  0 fallos.
#   C1  SENUELO, la hipotesis al reves:  moved_bot contenido en tied  (o sea, el de abajo mira lo MISMO
#       que el de arriba).  Tiene que fallar mucho, o H-B no dice nada.
#   C2  SENUELO conocido: sumar todos los w de G_bot de golpe, sin repartir por D.  Debe discrepar.
#   C3  NULO CONTROLADO en t = 2: alli g_com es vacio y el teorema publicado dice que el estrato de
#       arriba basta.  La poblacion objetivo debe salir VACIA.  Si sale algo, el teorema de t = 2
#       esta mal y eso es lo primero que hay que saber.
#   C4  no vacuidad de la poblacion objetivo: n se imprime SIEMPRE, tambien si es 0.
#   C5  la reflexion NO es el mecanismo de abajo (Gbot_anatomy ya lo midio: 83 formas con bot = 0
#       sin clausura).  Se re-mide aqui para que este al lado de H-B y no se lea de mas.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python bottom_sees_gcom.py     (aritmetica entera pura, no hace falta Sage)

import itertools
import sys
from collections import defaultdict

# beta como N-subconjuntos de {0..M}: |lambda| sin acotar y parte mayor acotada.  Es el muestreo
# que uso el auditor y produce beta mucho mas dispersos que barrer por |lambda|.
CONFIGS = [(2, 2, 15), (2, 3, 15),
           (4, 2, 15), (4, 3, 15),
           (6, 2, 15), (6, 3, 17),
           (8, 2, 17), (8, 3, 19)]


def perm_sign(seq):
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


def topdeg_dict(T, r):
    """[A(T)]_top: sector todo-positivo, alternante en las dos mitades."""
    D = {}
    n = 2 * r
    for a in itertools.permutations(range(r)):
        for b in itertools.permutations(range(r)):
            q = [0] * n
            e = [0] * r
            for i in range(r):
                q[i] = 2 * a[i]
                e[a[i]] += T[i]
            for i in range(r):
                q[r + i] = 2 * b[i] + 1
                e[b[i]] -= T[r + i]
            k = tuple(e)
            D[k] = D.get(k, 0) + perm_sign(q)
    return dict((k, v) for k, v in D.items() if v != 0)


def botdeg_dict(T, r):
    """[A(T)]_bot: emparejamiento de consecutivos, todos los signos +1 (permanente)."""
    D = {}
    n = 2 * r
    for s in itertools.permutations(range(r)):
        q = [0] * n
        e = [0] * r
        for i in range(r):
            q[2 * i] = 2 * s[i]
            q[2 * i + 1] = 2 * s[i] + 1
            e[s[i]] = T[2 * i] - T[2 * i + 1]
        D[tuple(e)] = D.get(tuple(e), 0) + perm_sign(q)
    return dict((k, v) for k, v in D.items() if v != 0)


def inv_of(T, r):
    """el atajo INV de second_stratum.py, aqui solo para contrastarlo (A1)."""
    H, L = T[:r], T[r:]
    alpha = [H[i] - (r - 1 - i) for i in range(r)]
    atil = tuple(a - alpha[-1] for a in alpha)
    Lstar = [L[0] - L[r - 1 - i] for i in range(r)]
    astar = tuple(Lstar[i] - (r - 1 - i) for i in range(r))
    return (H[-1] - L[0], tuple(sorted([atil, astar])))


def analyse(beta, t, r):
    N = len(beta)
    cl = defaultdict(list)
    for i, b in enumerate(beta):
        cl[b % t].append(i)
    if len(cl) < t:                       # hipotesis de ocupacion (O)
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    if not E:
        return None
    S = sorted((beta[i] for k in E for i in cl[k]), reverse=True)
    Sset = set(S)
    C = S[0] + S[-1]
    cond_i = set(C - v for v in S) == Sset

    keys = sorted(cl)
    tr = []
    for pick in itertools.product(*[cl[k] for k in keys]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)     # ya decreciente
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        sel = tuple(beta[i] for i in pick)                      # valor retirado por clase (orden keys)
        dif = tuple(T[2 * i] - T[2 * i + 1] for i in range(r))
        tr.append((sel, T, w, sum(T[:r]) - sum(T[r:]), sum(dif), tuple(sorted(dif, reverse=True))))

    Dmax = max(x[3] for x in tr)
    Dmin = min(x[4] for x in tr)
    G = [x for x in tr if x[3] == Dmax]
    Gb = [x for x in tr if x[4] == Dmin]

    top = {}
    for x in G:
        for k, v in topdeg_dict(list(x[1]), r).items():
            top[k] = top.get(k, 0) + x[2] * v
    top = dict((k, v) for k, v in top.items() if v != 0)
    bot = {}
    for x in Gb:
        for k, v in botdeg_dict(list(x[1]), r).items():
            bot[k] = bot.get(k, 0) + x[2] * v
    bot = dict((k, v) for k, v in bot.items() if v != 0)

    # el recuento con signo por multiconjunto de diferencias
    byD = defaultdict(int)
    membersD = defaultdict(list)
    for x in Gb:
        byD[x[5]] += x[2]
        membersD[x[5]].append(x)
    bot_zero_count = all(v == 0 for v in byD.values())

    # atajo INV, solo para A1
    top_zero_inv = (len(G) == 2 and inv_of(G[0][1], r) == inv_of(G[1][1], r)
                    and G[0][2] == -G[1][2])

    idxE = [i for i, k in enumerate(keys) if k in E]
    tied = set()
    g_com = None
    if len(G) == 2:
        tied = set(keys[i] for i in range(len(keys)) if G[0][0][i] != G[1][0][i])
        g_com = set(G[0][0][i] for i in idxE if G[0][0][i] == G[1][0][i])

    moved_bot = set(keys[i] for i in range(len(keys))
                    if len(set(x[0][i] for x in Gb)) > 1)

    return dict(t=t, r=r, beta=beta, S=S, C=C, e=len(E), E=E, cond_i=cond_i,
                Dmax=Dmax, Dmin=Dmin, nG=len(G), nGb=len(Gb), G=G, Gb=Gb,
                top=top, bot=bot, top_zero=(not top), bot_zero=(not bot),
                bot_zero_count=bot_zero_count, byD=dict(byD), membersD=dict(membersD),
                top_zero_inv=top_zero_inv, tied=tied, g_com=g_com, moved_bot=moved_bot,
                keys=keys)


def main():
    tot = 0
    a0 = a1 = a2 = a2n = a3 = a3n = 0
    c2_bad = c2_n = 0
    c5_closed_botzero = c5_closed_only = c5_botzero_only = 0
    per_cfg = {}
    target = []                      # top_zero y sigma_C(g_com) != g_com
    b0_bad = b0_n = 0                # (i) falsa  <=>  g_com no simetrico, dado top_zero
    hand = None

    for (t, r, M) in CONFIGS:
        N = t + 2 * r
        if M < N - 1:
            print("  cfg t=%d r=%d M=%d : SALTADA (M < N-1)" % (t, r, M))
            continue
        n_cfg = ntgt = 0
        for comb in itertools.combinations(range(M + 1), N):
            beta = tuple(sorted(comb, reverse=True))
            a = analyse(beta, t, r)
            if a is None:
                continue
            n_cfg += 1
            tot += 1

            # A0  el recuento con signo por D  ==  el desarrollo monomial
            if a['bot_zero_count'] != a['bot_zero']:
                a0 += 1
            # A1  el atajo INV+signo  ==  la expansion
            if a['top_zero_inv'] != a['top_zero']:
                a1 += 1
            # A2  sigma_C conserva degmin y D  (sobre los propios T, como conjuntos de enteros)
            # (sin tope: la primera version corria solo Gb[:2] y el tope no estaba escrito en
            #  ninguna parte -- el mismo defecto que depth.py de la noche del 12.)
            for x in a['Gb']:
                U = tuple(sorted((a['C'] - v for v in x[1]), reverse=True))
                dif = tuple(sorted((U[2 * i] - U[2 * i + 1] for i in range(r)), reverse=True))
                a2n += 1
                if dif != x[5] or sum(dif) != a['Dmin']:
                    a2 += 1
            # A3  la descomposicion de S
            if a['nG'] == 2:
                K = set(a['G'][0][1]) & set(a['G'][1][1])
                four = (set(a['G'][0][1]) ^ set(a['G'][1][1])) | (
                    set(a['G'][0][0]) ^ set(a['G'][1][0]))
                a3n += 1
                if set(a['S']) != K | a['g_com'] | four or len(a['g_com']) != a['e'] - 2:
                    a3 += 1
            # C2  senuelo: sumar todos los w de golpe
            c2_n += 1
            if (sum(x[2] for x in a['Gb']) == 0) != a['bot_zero']:
                c2_bad += 1
            # C5  clausura de G_bot bajo sigma_C  vs  bot = 0
            setsGb = set(x[1] for x in a['Gb'])
            closed = all(tuple(sorted((a['C'] - v for v in T), reverse=True)) in setsGb
                         for T in setsGb)
            if closed and a['bot_zero']:
                c5_closed_botzero += 1
            elif closed:
                c5_closed_only += 1
            elif a['bot_zero']:
                c5_botzero_only += 1

            if a['top_zero']:
                sym = (a['g_com'] is not None
                       and set(a['C'] - v for v in a['g_com']) == a['g_com'])
                b0_n += 1
                if sym == a['cond_i']:
                    pass
                else:
                    b0_bad += 1
                if not sym:
                    ntgt += 1
                    target.append(a)
                    if hand is None and t >= 4 and a['e'] > 2:
                        hand = a
        per_cfg[(t, r, M)] = (n_cfg, ntgt)
        print("  cfg t=%d r=%d M=%d : %d formas, %d en la poblacion objetivo" % (t, r, M, n_cfg, ntgt))
        sys.stdout.flush()

    print("")
    print("=" * 100)
    print("ACEPTACION -- si alguna de estas no es 0, lo demas no significa nada")
    print("=" * 100)
    print("  A0  [Phi]_bot por recuento con signo != por expansion : %d fallos de %d" % (a0, tot))
    print("  A1  atajo INV+signo != expansion para [Phi]_top       : %d discrepancias de %d" % (a1, tot))
    print("  A2  sigma_C conserva degmin y el multiconjunto D      : %d fallos de %d" % (a2, a2n))
    print("  A3  S = K |_| g_com |_| {x1,x2,y1,y2}, |g_com| = e-2  : %d fallos de %d" % (a3, a3n))
    fatal = a0 + a2 + a3
    print("")
    print("=" * 100)
    print("B0  dado [Phi]_top = 0:   (i) falsa   <=>   sigma_C(g_com) != g_com")
    print("=" * 100)
    print("  %d discrepancias de %d formas con [Phi]_top = 0" % (b0_bad, b0_n))
    print("")
    print("=" * 100)
    print("C4  POBLACION OBJETIVO  ([Phi]_top = 0  y  g_com no simetrico) : n = %d" % len(target))
    print("=" * 100)
    for k in sorted(per_cfg):
        print("      t=%d r=%d M=%d : %d de %d formas" % (k[0], k[1], k[2], per_cfg[k][1], per_cfg[k][0]))
    t2 = [a for a in target if a['t'] == 2]
    print("")
    print("  C3  NULO CONTROLADO: poblacion objetivo en t = 2 (el teorema publicado dice 0) : %d"
          % len(t2))

    if not target:
        print("")
        print("  poblacion objetivo VACIA -- H-A, H-B y H-C quedan SIN TESTAR.  No hay veredicto.")
        print("DONE")
        return 0 if fatal == 0 else 1

    # ---------------------------------------------------------------- H-A, H-B, H-C -------------
    ha = sum(1 for a in target if not a['bot_zero'])
    hb = hb_n = 0
    c1_bad = 0
    hc = hc_n = hc_single = 0
    depth = defaultdict(int)
    for a in target:
        asym = set(v for v in a['g_com'] if (a['C'] - v) not in a['g_com'])
        asym_cls = set(v % a['t'] for v in asym)
        hb_n += 1
        if a['moved_bot'] & asym_cls:
            hb += 1
        if a['moved_bot'] <= a['tied']:
            c1_bad += 1
        # H-C: la clase D que sobrevive, en que clases difieren sus miembros
        surv = [D for D, s in a['byD'].items() if s != 0]
        for D in surv:
            mem = a['membersD'][D]
            hc_n += 1
            diffc = set(a['keys'][i] for i in range(len(a['keys']))
                        if len(set(x[0][i] for x in mem)) > 1)
            # SIN el disyunto de conveniencia que tenia la primera version
            # ("o la clase D tiene un solo miembro y moved_bot no esta en tied"), que hacia H-C
            # mas facil de pasar sin decirlo.  Las clases D de un solo miembro se cuentan aparte:
            # ahi la pregunta "en que difieren sus miembros" no tiene sentido.
            if len(mem) == 1:
                hc_single += 1
            elif diffc & asym_cls:
                hc += 1
        depth[a['Dmax'] - a['Dmin']] += 1

    print("")
    print("=" * 100)
    print("H-A   top = 0  y  g_com no simetrico   =>   [Phi]_bot != 0")
    print("=" * 100)
    print("  %d de %d   (fallos: %d)" % (ha, len(target), len(target) - ha))
    print("")
    print("=" * 100)
    print("H-B   la variacion de G_bot toca las clases donde g_com pierde la simetria")
    print("=" * 100)
    print("  moved_bot & asym_cls != vacio : %d de %d" % (hb, hb_n))
    print("  C1 SENUELO  moved_bot contenido en tied (abajo mira lo mismo que arriba) : %d de %d"
          % (c1_bad, hb_n))
    print("")
    print("=" * 100)
    print("H-C   la clase D que SOBREVIVE esta separada por una clase de g_com")
    print("=" * 100)
    print("  %d de %d clases D supervivientes  (y %d de las %d tienen UN SOLO miembro, donde la"
          % (hc, hc_n, hc_single, hc_n))
    print("  pregunta no tiene sentido: no hay dos transversales que difieran en nada)")
    print("")
    print("  profundidad Dmax - Dmin en la poblacion objetivo:")
    for k in sorted(depth):
        print("      %s : %d" % (k, depth[k]))

    # ---------------------------------------------------------------- los fallos de H-A ---------
    fails = [a for a in target if a['bot_zero']]
    if fails:
        print("")
        print("=" * 100)
        print("LOS FALLOS DE H-A -- top = 0 Y bot = 0 Y sigma_C(S) != S.  Si son reales, el enunciado")
        print("'[Phi]_top = 0 y [Phi]_bot = 0  =>  sigma_C(S) = S'  es FALSO.")
        print("=" * 100)
        print("")
        print("     t   r  beta                              e  |G_bot|  Dmax  Dmin  moved_bot  tied"
              "  moved<=tied")
        print("  " + "-" * 96)
        for a in fails:
            print("  %4d %3d  %-32s %2d %7d %5d %5d  %-9s %-5s %s"
                  % (a['t'], a['r'], str(list(a['beta'])), a['e'], a['nGb'], a['Dmax'], a['Dmin'],
                     str(sorted(a['moved_bot'])), str(sorted(a['tied'])),
                     a['moved_bot'] <= a['tied']))
        nc1 = sum(1 for a in fails if a['moved_bot'] <= a['tied'])
        print("")
        print("  cruce con el senuelo C1: %d de los %d fallos tienen moved_bot contenido en tied"
              % (nc1, len(fails)))

    print("")
    print("=" * 100)
    print("SENUELOS Y CONTROLES SOBRE TODAS LAS FORMAS")
    print("=" * 100)
    print("  C2  sumar todos los w de G_bot sin repartir por D : %d discrepancias de %d"
          % (c2_bad, c2_n))
    print("  C5  clausura de G_bot bajo sigma_C   vs   [Phi]_bot = 0:")
    print("        cerrado y bot = 0 : %d      cerrado y bot != 0 : %d      bot = 0 sin clausura : %d"
          % (c5_closed_botzero, c5_closed_only, c5_botzero_only))

    # ---------------------------------------------------------------- el caso a mano ------------
    if hand is not None:
        a = hand
        print("")
        print("=" * 100)
        print("EL CASO A MANO")
        print("=" * 100)
        print("  t = %d, r = %d, beta = %s" % (a['t'], a['r'], list(a['beta'])))
        print("  S = %s     C = %d     e = %d     (i) concentrico: %s"
              % (a['S'], a['C'], a['e'], a['cond_i']))
        print("  clases de exceso E = %s     empatadas (tied) = %s" % (a['E'], sorted(a['tied'])))
        print("  g_com = %s     sigma_C(g_com) = %s"
              % (sorted(a['g_com'], reverse=True),
                 sorted((a['C'] - v for v in a['g_com']), reverse=True)))
        print("  |G| = %d  Dmax = %d      |G_bot| = %d  Dmin = %d"
              % (a['nG'], a['Dmax'], a['nGb'], a['Dmin']))
        print("  moved_bot = %s" % sorted(a['moved_bot']))
        print("")
        print("  los maximizadores:")
        for x in a['G']:
            print("      retira %s -> T = %s   w = %+d" % (list(x[0]), list(x[1]), x[2]))
        print("")
        print("  el estrato de abajo, por multiconjunto D:")
        for D in sorted(a['byD'], reverse=True):
            print("      D = %-22s suma w = %+d   miembros:" % (str(list(D)), a['byD'][D]))
            for x in a['membersD'][D]:
                print("          retira %-18s T = %-26s w = %+d"
                      % (str(list(x[0])), str(list(x[1])), x[2]))
        print("")
        print("  [Phi]_bot = 0 ? %s        [Phi]_top = 0 ? %s" % (a['bot_zero'], a['top_zero']))

    print("")
    if fatal:
        print("ACEPTACION FALLIDA -- veredicto suspendido.")
        return 1
    print("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
