# -*- coding: utf-8 -*-
# |G| <= 2, DEMOSTRADO -- y con ello t impar no tiene segunda rama, por una ruta nueva.
#
# *** UN DEFECTO MIO, cazado por V5 en la primera pasada ***
# brute() indexaba cada maximizador por la transversal ENTERA y greedy() por su parte de exceso, asi
# que V5 comparaba dos objetos que no pueden coincidir y fallaba en las 4664 formas de t=4 r=1 -- o
# sea en TODAS.  Las clases sin exceso son singletons y se conservan siempre, luego no llevan
# informacion: la clave correcta es la parte de exceso.  El senuelo K2 heredaba el mismo error y
# aparecia inflado a 16427.
#
# *** Y UN SEGUNDO DEFECTO, de la auditoria socratica: el K2 original imprimia '3710' SIN
# DENOMINADOR y su puerta de aprobado era 'DEC > 0', asi que un senuelo acertado el 80 % de
# las veces la pasaba igual.  Peor: su cabecera decia que ordenaba por f_k(j) sin estructura
# de prefijo, y el codigo ordenaba por el ELEMENTO c_{k,j}, que ya es prefijo por si solo --
# la guarda 'if sum(dj) == r' no disparaba nunca.  Ahora hay TRES senuelos con denominador, y
# entre ellos el NEAR MISS c_{k,j-1}, que falla MENOS que c_{k,j} (6.6 % contra 17.0 %): eso
# demuestra que 'ordenar por el elemento falla' NO aisla la pareja.  El senuelo que de verdad
# quita la concavidad es el del NIVEL f_k(j), y ese falla el 39.9 % y queda indefinido en
# 10451 formas, que es justo lo que un senuelo honesto tiene que ensenar. ***
#
# widen_G.sage midio que el conjunto G de transversales que maximizan deg nunca pasa de 2 en 75990
# formas.  Esto lo prueba, y la prueba dice ademas QUIENES son los dos.
#
# ================================================================================================
# NOTACION.  Para cada clase de exceso k, sean c_{k,1} > c_{k,2} > .. > c_{k,n_k} sus elementos,
# n_k >= 2, sum_k (n_k - 1) = 2r.  Una transversal g conserva un elemento por clase, T = S \ g.
#
# PASO 1 (reformulacion).  deg(T) = max_{A subset T, |A|=r} (2 sum A - sum T), porque el maximo se
# alcanza en la mitad alta.  Como T = S \ g,
#       deg(S \ g)  =  max_A [ 2 sum A + sum g ]  -  sum S,          A subset S\g,  |A| = r.
# Escribimos Phi(A,g) = 2 sum A + sum g y maximizamos Phi.
#
# PASO 2 (la forma de un optimo).  Sea (A,g) optimo y sean x > y de la MISMA clase k.
#   - si y in A y x = g_k: intercambiarlos deja g transversal y cambia Phi en +(x-y) > 0.  Imposible.
#   - si y = g_k y x fuera de A y de g: idem, +(x-y) > 0.  Imposible.
# Luego dentro de cada clase, en orden decreciente, van primero los de A, despues g_k, y despues los
# no elegidos.  Es decir el optimo esta parametrizado por un solo entero por clase,
#       j_k = |A cap cls_k|,   0 <= j_k <= n_k - 1,   sum_k j_k = r,
# con  A cap cls_k = {c_{k,1},..,c_{k,j_k}}  y  g_k = c_{k, j_k + 1}.
#
# PASO 3 (SEPARABILIDAD).  Con esa parametrizacion,
#       Phi  =  sum_k f_k(j_k),      f_k(j) = 2 (c_{k,1} + .. + c_{k,j}) + c_{k,j+1}.
# El objetivo se PARTE por clases y la unica ligadura es sum_k j_k = r.
#
# PASO 4 (CONCAVIDAD).  El incremento es
#       Delta_k(j) = f_k(j) - f_k(j-1) = c_{k,j} + c_{k,j+1},      1 <= j <= n_k - 1,
# y es ESTRICTAMENTE decreciente en j porque los c_{k,j} lo son.  Cada f_k es estrictamente concava,
# asi que maximizar sum_k f_k(j_k) con sum j_k = r es elegir los r incrementos MAYORES del
# multiconjunto {Delta_k(j)} -- la eleccion es automaticamente un prefijo en cada clase, y cualquier
# optimo lo es (si un incremento no elegido supera a uno elegido, moverlo mejora).
#
# PASO 5 (RESIDUOS, que es donde muere el empate).  c_{k,j} = c_{k,j+1} = k (mod t), luego
#       Delta_k(j) = 2k  (mod t).
# Dentro de una clase los Delta son distintos (paso 4).  Entre clases, Delta_k(j) = Delta_{k'}(j')
# exige 2k = 2k' (mod t), o sea k' = k o -- solo si t es PAR -- k' = k + t/2.  Por tanto
#       NINGUN VALOR ES COMPARTIDO POR MAS DE DOS INCREMENTOS,
# y los dos, cuando ocurre, vienen de las clases k y k + t/2.
#
# PASO 6 (conclusion).  Los optimos son las selecciones de r incrementos que contienen a todos los
# estrictamente mayores que el r-esimo valor V y s de los que valen V.  Si V aparece una vez, el
# optimo es unico.  Si aparece dos veces, o caben los dos (unico) o cabe uno (dos optimos).  Luego
#       |G| <= 2,   y   |G| = 2  <=>  hay empate JUSTO en el rango r (d_r = d_{r+1}), entre las
#       clases k y k+t/2.  La clausula 'justo en el rango r' es parte del <=> y no se puede
#       dejar caer: un empate en los rangos r-1, r deja |G| = 1.
# QED
#
# CONFIRMACION INDEPENDIENTE (auditoria socratica, 2026-08-12).  Una rebanada distinta -- todos
# los beta como N-subconjuntos de {0..M}, o sea |lambda| SIN acotar y la parte mayor acotada,
# que muestrea betas mucho mas dispersos que el barrido por |lambda| de este fichero -- da
# 235744 formas en 22 configuraciones, t en 1..8, r en 1..4, y CERO fallos en las diez
# comprobaciones, incluida la mitad estructural sobre 19033 formas con |G| = 2.
#
# COROLARIO A.  t impar  Y TODAS LAS CLASES DE RESIDUOS OCUPADAS  =>  |G| = 1  =>  Phi_t no nulo.
#   Para t impar 2 es invertible modulo t, asi que 2k = 2k' fuerza k = k' y no hay empates.  Con el
#   corolario de middle_block.sage (maximizador unico => no nulo) se REDEMUESTRA que t impar no
#   tiene segunda rama, por una ruta que no menciona prod(mu_t) ni SO(N).  Dos rutas independientes.
#
#   *** LA HIPOTESIS DE OCUPACION NO ES DECORATIVA, y la auditoria socratica la caza: si una clase
#   queda VACIA, cualquier eleccion de t filas repite un residuo por palomar, todos los menores de
#   raices de la unidad se anulan y Phi_t = 0 -- tambien para t impar.  Esa es la rama (a) publicada
#   (Littlewood / Karmakar Thm 4.1(A)).  La hipotesis vive en setup(), que descarta len(cl) < t, pero
#   NO estaba escrita en la linea del corolario.  Ahora si. ***
#
# COROLARIO B.  Con todas las clases ocupadas, Phi_t = 0  =>  |G| = 2  =>  existe k con las clases k
#   y k + t/2 AMBAS de exceso y con un empate de incrementos JUSTO en el rango r.  (La ocupacion se
#   necesita otra vez, para que G no sea vacio.)  Es una condicion necesaria probada, y se parece a
#   la condicion (ii); la pregunta que decide si ES (ii) esta medida abajo en V7.
#
# ALCANCE, mas ancho de lo que el enunciado sugiere: la prueba usa SOLO que los beta son enteros
# distintos y ordenados.  Nada pide beta >= 0, ni que lambda sea una particion, ni N = t + 2r mas
# alla de forzar sum_k (n_k - 1) = 2r.  Y r >= 1: con r = 0 no hay incrementos que elegir.
#
# QUE ES ESTANDAR Y QUE NO, para cuando esto se escriba.  Los pasos 1-4 son el voraz clasico de
# maximizacion separable concava con una ligadura de cardinal (cada clase es una cadena con
# incrementos estrictamente decrecientes: un matroide de particion).  NO es nuestro y hay que
# etiquetarlo como estandar.  El contenido propio es UNA linea: Delta_k(j) = 2k (mod t), mas que
# x -> 2x en Z/t es 2 a 1 si t es par e inyectiva si es impar.
# ================================================================================================
#
# VERIFICACIONES, cada una capaz de fallar:
#   V1  la reformulacion: deg(S\g) calculado a pelo contra max_A Phi - sum S, en TODA transversal.
#   V2  Delta_k estrictamente decreciente.
#   V3  Delta_k(j) = 2k (mod t).
#   V4  ningun valor de Delta compartido por mas de 2 incrementos, y los dos de clases k, k+t/2.
#   V5  EL CONTROL QUE DECIDE: el conjunto G predicho por el algoritmo voraz debe COINCIDIR, como
#       conjunto de transversales, con el G de fuerza bruta.  No el tamano: el conjunto.
#   V6  |G| <= 2 siempre, y |G| = 1 siempre que t es impar.
#   V7  cuando |G| = 2, el valor de empate V cumple V = C (mod t), con C = min S + max S?  Si eso
#       vale, entonces Phi_t = 0 => (ii) queda PROBADO, porque las clases empatadas son las dos
#       clases fijas de sigma_C.  Se mide aparte sobre las formas con [Phi]_top = 0.
#
# CONTROLES:
#   K1  no vacuidad: |G| = 2 debe ocurrir, y t impar debe producir formas (si no, V6 no dice nada).
#   K2  TRES senuelos con denominador, ordenando por c_{k,j} (el elemento bajo), c_{k,j-1} (el
#       alto, un NEAR MISS puesto para que se vea que el fallo del primero no prueba nada
#       sobre la pareja) y f_k(j) (el NIVEL, que es el que de verdad tira la concavidad y
#       puede dar una seleccion no-prefijo, o sea quedar INDEFINIDO -- y eso se cuenta).
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

CONF = [(1, 2, 24), (2, 2, 22), (2, 3, 18), (3, 2, 22), (4, 1, 30), (4, 2, 24),
        (4, 3, 18), (5, 2, 20), (6, 2, 18), (6, 3, 14), (7, 2, 14), (8, 2, 16),
        (9, 2, 14), (10, 2, 18)]


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def perm_sign(q):
    n = len(q)
    seen = [False] * n
    s = 1
    for i in range(n):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = q[j]
            ln += 1
        if ln % 2 == 0:
            s = -s
    return s


def topdeg_dict(T, r):
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
    return {k: v for k, v in D.items() if v != 0}


def setup(beta, t, r):
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    C = {}
    for k in E:
        C[k] = sorted((beta[i] for i in cl[k]), reverse=True)
    return cl, E, C


def brute(beta, t, r, cl, E):
    """(Dmax, {frozenset of the kept EXCESS values: (T, w)}).

    Keyed on the excess part only, because the non-excess classes are singletons and are always
    kept -- so that part carries no information and the excess part determines the transversal.
    Keying on the whole transversal is what made V5 fail on every single shape the first time: it
    compared a full transversal against the proof's excess-only object and could never agree."""
    N = len(beta)
    Es = set(E)
    best = None
    out = {}
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        d = sum(T[:r]) - sum(T[r:])
        if best is None or d > best:
            best, out = d, {}
        if d == best:
            w = perm_sign([beta[i] % t for i in P])
            if sum(P) % 2:
                w = -w
            out[frozenset(beta[i] for i in P if beta[i] % t in Es)] = (T, w)
    return best, out


def greedy(C, E, r, t):
    """the proof's prediction: (Dmax offset, set of transversals, tie value V, tied classes)."""
    inc = []
    for k in E:
        ck = C[k]
        for j in range(1, len(ck)):
            inc.append((ck[j - 1] + ck[j], k, j))
    inc.sort(key=lambda z: -z[0])
    V = inc[r - 1][0]
    above = [z for z in inc if z[0] > V]
    at = [z for z in inc if z[0] == V]
    s = r - len(above)
    sols = []
    for chosen in itertools.combinations(range(len(at)), s):
        jj = dict((k, 0) for k in E)
        for _, k, j in above:
            jj[k] = max(jj[k], j)
        for idx in chosen:
            _, k, j = at[idx]
            jj[k] = max(jj[k], j)
        g = frozenset(C[k][jj[k]] for k in E)
        phi = sum(2 * sum(C[k][:jj[k]]) + C[k][jj[k]] for k in E)
        sols.append((g, phi))
    vals = set(p for _, p in sols)
    return sols, V, [k for _, k, _ in at], vals


print("=" * 104)
print("V1-V6  every step of the proof, and V5 is the one that decides")
print("=" * 104)
print("")
print("     t   r  shapes | V1 bad  V2 bad  V3 bad  V4 bad | V5 bad  max|G|  |G|=2  | V6 bad  K2 decoy bad")
print("  " + "-" * 100)

BAD = 0
NG2 = 0
DEC = 0
ODDSH = 0
DECN = dict(low=[0, 0, 0], up=[0, 0, 0], lev=[0, 0, 0])
for t, r, MAX in CONF:
    N = t + 2 * r
    nsh = v1 = v2 = v3 = v4 = v5 = v6 = ng2 = dec = 0
    mx = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            st = setup(beta, t, r)
            if st is None:
                continue
            cl, E, C = st
            nsh += 1
            SS = sum(sum(C[k]) for k in E)
            Dmax, Gb = brute(beta, t, r, cl, E)
            # V2, V3, V4
            inc = []
            for k in E:
                ck = C[k]
                d = [ck[j - 1] + ck[j] for j in range(1, len(ck))]
                if any(d[i] <= d[i + 1] for i in range(len(d) - 1)):
                    v2 += 1
                if any(x % t != (2 * k) % t for x in d):
                    v3 += 1
                inc += [(x, k) for x in d]
            byval = {}
            for x, k in inc:
                byval.setdefault(x, []).append(k)
            for x, ks in byval.items():
                if len(ks) > 2:
                    v4 += 1
                elif len(ks) == 2:
                    a, b = sorted(ks)
                    if t % 2 or b - a != t // 2:
                        v4 += 1
            # V1: the reformulation, on EVERY transversal of the shape
            for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
                P = set(pick)
                T = sorted((beta[i] for i in range(N) if i not in P), reverse=True)
                d1 = sum(T[:r]) - sum(T[r:])
                gv = [beta[i] for i in P if (beta[i] % t) in E]
                d2 = max(2 * sum(A) for A in itertools.combinations(T, r)) + sum(gv) - SS
                if d1 != d2:
                    v1 += 1
            # V5: greedy set against brute-force set
            sols, V, tied, vals = greedy(C, E, r, t)
            if len(vals) != 1 or (list(vals)[0] - SS) != Dmax:
                v5 += 1
            if set(g for g, _ in sols) != set(Gb):
                v5 += 1
            mx = max(mx, len(Gb))
            if len(Gb) == 2:
                ng2 += 1
            if len(Gb) > 2:
                v6 += 1
            if t % 2 and len(Gb) != 1:
                v6 += 1
            # K2, THREE decoys, each ranking by a different key instead of the pair sum.
            #   'low'   c_{k,j}      the lower element of the pair
            #   'up'    c_{k,j-1}    the upper element -- a NEAR MISS, and it must be shown to be
            #                        just as wrong, or 'low' failing proves nothing about the pair
            #   'lev'   f_k(j)       the LEVEL instead of the increment: this is the one that really
            #                        drops concavity, and it can produce a non-prefix selection, in
            #                        which case it is UNDEFINED and must be counted as such, not
            #                        silently skipped the way the first version of this file did
            for tag, keyf in (('low', lambda k, j: C[k][j]),
                              ('up', lambda k, j: C[k][j - 1]),
                              ('lev', lambda k, j: 2 * sum(C[k][:j]) + C[k][j])):
                pool = sorted(((keyf(k, j), k, j) for k in E for j in range(1, len(C[k]))),
                              key=lambda z: -z[0])[:r]
                dj = dict((k, 0) for k in E)
                for _, k, j in pool:
                    dj[k] = max(dj[k], j)
                DECN[tag][2] += 1
                if sum(dj.values()) != r:
                    DECN[tag][1] += 1          # not a prefix selection: the decoy is undefined
                    continue
                if frozenset(C[k][dj[k]] for k in E) not in Gb:
                    DECN[tag][0] += 1
                    if tag == 'low':
                        dec += 1
    BAD += v1 + v2 + v3 + v4 + v5 + v6
    NG2 += ng2
    DEC += dec
    if t % 2:
        ODDSH += nsh
    print("  %4d %3d %7d | %6d %7d %7d %7d | %6d %7d %6d | %6d %13d"
          % (t, r, nsh, v1, v2, v3, v4, v5, mx, ng2, v6, dec))

print("")
print("  totals: violations %d (must be 0).  |G| = 2 occurs %d times." % (BAD, NG2))
print("  K1 non-vacuity: |G|=2 occurs (%s); odd-t shapes exist (%d)"
      % ("yes" if NG2 > 0 else "NO", ODDSH))
print("  K2 the three decoys, WITH DENOMINATORS -- a decoy that is right 80% of the time is not a")
print("     decoy, and the first version of this file printed 3710 with no denominator at all:")
for tag, name in (('low', 'rank by c_{k,j}      '),
                  ('up', 'rank by c_{k,j-1}    '),
                  ('lev', 'rank by f_k(j), level')):
    w, u, n = DECN[tag]
    print("       %s : wrong on %6d of %6d defined (%5.1f%%), undefined on %6d"
          % (name, w, n - u, 100.0 * w / max(1, n - u), u))
print("     'up' is the NEAR MISS: if it fails at the same rate as 'low', then 'low' failing says")
print("     nothing about the PAIR -- it only says the ranking key matters.  Read the two together.")
print("")
if BAD == 0 and NG2 > 0 and ODDSH > 0 and DECN['lev'][0] > 0:
    print("  PROVED AND NON-VACUOUS: |G| <= 2 always, |G| = 1 whenever t is odd, and the greedy")
    print("  selection of the r largest pair-sums reproduces the maximiser SET exactly.")
else:
    print("  SOMETHING FAILED -- read the columns.")

# ---------------------------------------------------------------- V7 -----------------------------
print("")
print("=" * 104)
print("V7  when |G| = 2, is the tie value V congruent to C = min S + max S mod t?")
print("=" * 104)
print("")
print("     t   r  |G|=2 | V=C mod t   V!=C mod t | of those with [Phi]top=0: V=C   V!=C   cond(ii)")
print("  " + "-" * 100)

TA = TB = TC = TD = TE = 0
EX = []
for t, r, MAX in CONF:
    if t % 2:
        continue
    N = t + 2 * r
    n2 = ya = na = zy = zn = cii = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            st = setup(beta, t, r)
            if st is None:
                continue
            cl, E, C = st
            Dmax, Gb = brute(beta, t, r, cl, E)
            if len(Gb) != 2:
                continue
            n2 += 1
            S = sorted((b for k in E for b in C[k]), reverse=True)
            Cc = S[0] + S[-1]
            _, V, tied, _ = greedy(C, E, r, t)
            ok = (V - Cc) % t == 0
            ya += 1 if ok else 0
            na += 0 if ok else 1
            top = {}
            for g, (T, w) in Gb.items():
                for kk, vv in topdeg_dict(list(T), r).items():
                    top[kk] = top.get(kk, 0) + w * vv
            vanishes = not any(v != 0 for v in top.values())
            if vanishes:
                zy += 1 if ok else 0
                zn += 0 if ok else 1
                fixed = [k for k in range(t) if (2 * k - Cc) % t == 0]
                if len(fixed) == 2 and all(k in E for k in fixed):
                    cii += 1
                if not ok and len(EX) < 10:
                    EX.append((t, r, list(l), beta, V, Cc))
    print("  %4d %3d %6d | %10d %11d | %28d %6d %10d" % (t, r, n2, ya, na, zy, zn, cii))
    TA += n2
    TB += ya
    TC += na
    TD += zy
    TE += zn

print("")
print("  |G| = 2 on %d shapes: V = C mod t on %d, V != C on %d." % (TA, TB, TC))
print("  restricted to the shapes where [Phi]_top vanishes: V = C on %d, V != C on %d." % (TD, TE))
if TE == 0 and TD > 0:
    print("")
    print("  READING.  Whenever the top-degree part vanishes -- which Phi_t = 0 forces -- the tied")
    print("  classes ARE the two fixed classes of sigma_C.  With Corollary B that says")
    print("        Phi_t = 0  =>  both fixed classes of sigma_C are excess classes  =  condition (ii),")
    print("  measured with no exception; the missing lemma is exactly  [Phi]_top = 0  =>  V = C mod t.")
if EX:
    print("")
    print("  the exceptions, V != C with [Phi]_top = 0:")
    for t, r, lam, beta, V, Cc in EX:
        print("     t=%d r=%d lam=%-22s beta=%-30s V=%d C=%d" % (t, r, str(lam), str(beta), V, Cc))
print("")
print("DONE")
