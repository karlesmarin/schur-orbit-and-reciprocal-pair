# -*- coding: utf-8 -*-
# LA CADENA CERRADA: Phi_t = 0 => (ii), salvo UN solo enunciado, y ya no es nuestro.
#
# Dos pasos que ayer eran "medidos" se prueban aqui, y con ellos el eslabon se reduce a A1 sola.
#
# ================================================================================================
# RECORDATORIO.  |G| <= 2 (probado).  Si |G| = 1, [Phi]_top = +-P(T) != 0.  Luego Phi_t = 0 fuerza
# |G| = 2, y los dos maximizadores difieren en un intercambio entre las clases k y k' = k+t/2:
#       T_A = K u {a1, a2},  T_B = K u {b1, b2},   a1 + b1 = a2 + b2 = V   (el empate),
# con a1, b1 en la clase k y a2, b2 en la k'.  Y  [Phi]_top = w_A P(T_A) + w_B P(T_B).
#
# A1 (medido exhaustivamente, link_P_separates): P(T) = +-P(T')  =>  T' = T + m  o  T' = c - T.
#
# ------------------------------------------------------------------------------------------------
# PASO 5, PROBADO.  La rama de traslacion IMPLICA la de reflexion.
#   Si T_B = T_A + m con m != 0, entonces |T_A cap (T_A+m)| = |K| = 2r-2, o sea T_A es union de
#   EXACTAMENTE DOS progresiones aritmeticas de paso m,
#         T_A = {p, p+m, .., p+(u-1)m}  u  {q, q+m, .., q+(v-1)m},   u+v = 2r,
#   y entonces  DA = T_A \ T_B = {p, q}  (los dos comienzos)  y  DB = {p+um, q+vm}  (los dos finales
#   desplazados).  El empate dice sigma_V(DA) = DB, y solo hay dos emparejamientos:
#     alpha)  V-p = p+um  y  V-q = q+vm.  Entonces V - T_A = {p+m..p+um} u {q+m..q+vm} = T_B.
#     beta )  V-p = q+vm  y  V-q = p+um.  Restando, um = vm, luego u = v, y otra vez V - T_A = T_B.
#   En ambos casos T_B = V - T_A.  QED
#
# PASO 6, PROBADO.  El centro de la reflexion es V, nunca el cruzado.
#   Si T_B = c - T_A entonces, como sigma_c es involucion que cruza los dos, sigma_c(K) = K y
#   sigma_c(DA) = DB.  Dos emparejamientos:
#     recto )  c - a1 = b1 y c - a2 = b2, o sea c = a1+b1 = V.
#     cruzado) c - a1 = b2 y c - a2 = b1.  Sumando los dos: a1+a2+b1+b2 = 2c, y el empate da
#              a1+a2+b1+b2 = 2V, luego c = V TAMBIEN.  Pero c = a1+b2 = k+k' = 2k + t/2 (mod t)
#              mientras V = 2k (mod t): exigiria t/2 = 0 (mod t), falso para t >= 2.
#   Luego el cruzado no existe y c = V.  QED
# ------------------------------------------------------------------------------------------------
#
# LA CADENA, ENTERA:
#   1. Phi_t = 0  =>  [Phi]_top = 0                                        trivial
#   2. |G| <= 2 y |G| = 1 => [Phi]_top != 0,  luego |G| = 2                PROBADO (proof_G_le_2)
#   3. [Phi]_top = 0  =>  P(T_A) = +- P(T_B)                               trivial (w = +-1)
#   4. P(T_A) = +-P(T_B)  =>  T_B = T_A + m  o  T_B = c - T_A              *** A1, solo MEDIDO ***
#   5. la rama traslacion implica la de reflexion                          PROBADO aqui
#   6. el centro es V                                                      PROBADO aqui
#   7. T_B = V - T_A  =>  V = C                                            PROBADO, ver abajo
#      *** LA LINEA QUE FALTABA, y la debo a la auditoria socratica: lemma_V_eq_C HIPOTETIZA
#      T_B = C - T_A con C el centro de la involucion y concluye V = C, o sea es el paso 6 otra
#      vez.  Lo que faltaba escribir es: sigma_V(T_A) = T_B => sigma_V(K) = K y sigma_V(DA) = DB
#      => sigma_V(S) = S => min S + max S = V.  Una linea, trivialmente cierta, nunca escrita --
#      y sin ella los pasos 6 y 7 parecen dos resultados independientes cuando son uno. ***
#   8. V = 2k (mod t) y V = C  =>  condicion (ii)                          PROBADO (paso 5 de |G|<=2)
# El unico hueco es A1, y A1 no habla de particiones ni de raices de la unidad: es la rigidez de un
# PRODUCTO DE DOS FUNCIONES DE SCHUR.  Si esta en la literatura, el recíproco-mitad-(ii) es teorema.
#
# VERIFICADO AQUI, cada cosa capaz de fallar:
#   S1  en cada forma con traslacion, T_A es union de EXACTAMENTE dos PA de paso m, DA son los dos
#       comienzos y DB los dos finales desplazados.  Es la sustancia del paso 5.
#   S2  y en cada una de ellas la reflexion T_B = V - T_A vale.  Paso 5, la conclusion.
#   S3  cual de los dos emparejamientos ocurre, alpha o beta.  Los dos son legitimos; si uno sale 0
#       hay que decirlo, no esconderlo.
#   S4  paso 6: no hay reflexion de centro distinto de V, y el cruzado no ocurre nunca.
#   S5  no vacuidad: debe haber traslaciones (si no, el paso 5 no dice nada) y reflexiones sin
#       traslacion (si no, la rama de traslacion seria la unica y el paso 5 seria todo el asunto).
#   S6  control forzado sobre el conjunto: [Phi]_top = 0  <=>  T_B = V - T_A, cero desacuerdos.
#
# F10, dicho y no escondido: en t = 2 con r >= 2 la rama de TRASLACION del paso 5 no se activa
# NUNCA -- 0 de 2588 formas con |G|=2 en r=2 y 0 de 780 en r=3.  En r=1 las 84 son traslaciones,
# pero para conjuntos de dos elementos traslacion y reflexion son la misma condicion, asi que no se
# testa nada.  O sea: en el caso que se cierra, el paso 5 tiene CERO apoyo empirico y descansa
# entero en el argumento en prosa.  El argumento es correcto y esta re-derivado por el auditor,
# pero la asimetria entre lo probado y lo medido hay que decirla.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

# F1, de la auditoria socratica: este fichero EXCLUIA t = 2.  Anadido.
CONF = [(2, 1, 34), (2, 2, 26), (2, 3, 20), (4, 1, 30), (4, 2, 24), (4, 3, 18),
        (6, 2, 18), (6, 3, 14), (8, 2, 16), (10, 2, 18)]


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


def P_of(T, r):
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


def chains(T, m):
    """decompose T into maximal arithmetic progressions of step m; return list of (start, length)."""
    Ts = set(T)
    out = []
    for x in sorted(Ts):
        if x - m not in Ts:
            ln = 1
            while x + ln * m in Ts:
                ln += 1
            out.append((x, ln))
    return out


def analyse(beta, t, r):
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    CC = dict((k, sorted((beta[i] for i in cl[k]), reverse=True)) for k in E)
    best = None
    G = {}
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        Pp = sorted(pick)
        Ps = set(Pp)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        d = sum(T[:r]) - sum(T[r:])
        if best is None or d > best:
            best, G = d, {}
        if d == best:
            w = perm_sign([beta[i] % t for i in Pp])
            if sum(Pp) % 2:
                w = -w
            G[frozenset(beta[i] for i in Pp)] = (T, w)
    inc = []
    for k in E:
        ck = CC[k]
        for j in range(1, len(ck)):
            inc.append((ck[j - 1] + ck[j], k))
    inc.sort(key=lambda z: -z[0])
    return G, inc[r - 1][0]


print("=" * 104)
print("S1-S6  step 5 (translation implies reflection) and step 6 (the centre is V)")
print("=" * 104)
print("")
print("     t   r  |G|=2 | transl | S1 bad  S2 bad  case a  case b | refl no transl | S4 bad | S6 bad")
print("  " + "-" * 100)

TOT = dict(g=0, tr=0, s1=0, s2=0, ca=0, cb=0, rn=0, s4=0, s6=0)
for t, r, MAX in CONF:
    N = t + 2 * r
    g2 = ntr = s1 = s2 = ca = cb = rn = s4 = s6 = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            an = analyse(beta, t, r)
            if an is None:
                continue
            G, V = an
            if len(G) != 2:
                continue
            g2 += 1
            (TA, wA), (TB, wB) = [G[g] for g in G]
            SA, SB = set(TA), set(TB)
            K = SA & SB
            DA = sorted(SA - K)
            DB = sorted(SB - K)
            refl = tuple(sorted((V - x for x in TA), reverse=True)) == tuple(sorted(TB, reverse=True))
            ms = [m for m in set(b - a for a in TA for b in TB)
                  if m != 0 and set(x + m for x in SA) == SB]
            if ms:
                ntr += 1
                m = ms[0]
                ch = chains(TA, m)
                starts = sorted(x for x, _ in ch)
                ends = sorted(x + ln * m for x, ln in ch)
                if len(ch) != 2 or starts != DA or ends != DB:
                    s1 += 1
                if not refl:
                    s2 += 1
                # which pairing: alpha keeps each chain in its own class, beta swaps the chains
                if len(ch) == 2:
                    p, u = ch[0]
                    if V - p == p + u * m:
                        ca += 1
                    else:
                        cb += 1
            elif refl:
                rn += 1
            # S4: a reflection about a centre other than V
            for c in set(x + y for x in TA for y in TB):
                if c != V and tuple(sorted((c - x for x in TA), reverse=True)) == \
                        tuple(sorted(TB, reverse=True)):
                    s4 += 1
                    break
            top = {}
            for TT_, ww in ((TA, wA), (TB, wB)):
                for kk, vv in P_of(list(TT_), r).items():
                    top[kk] = top.get(kk, 0) + ww * vv
            z = not any(v != 0 for v in top.values())
            if z != refl:
                s6 += 1
    print("  %4d %3d %6d | %6d | %6d %7d %7d %7d | %14d | %6d | %6d"
          % (t, r, g2, ntr, s1, s2, ca, cb, rn, s4, s6))
    for a, b in (('g', g2), ('tr', ntr), ('s1', s1), ('s2', s2), ('ca', ca), ('cb', cb),
                 ('rn', rn), ('s4', s4), ('s6', s6)):
        TOT[a] += b

print("")
print("  totals over %d shapes with |G| = 2:" % TOT['g'])
print("     translations: %d.  S1 (two APs, DA = starts, DB = shifted ends) bad: %d."
      % (TOT['tr'], TOT['s1']))
print("     S2 (every translation is also the reflection about V) bad: %d." % TOT['s2'])
print("     the pairing that occurs: alpha %d, beta %d." % (TOT['ca'], TOT['cb']))
print("     S4 (a reflection about a centre other than V) : %d" % TOT['s4'])
print("     S5 non-vacuity: translations %s, reflections WITHOUT translation %d."
      % ("yes" if TOT['tr'] else "NO", TOT['rn']))
print("     S6 ([Phi]_top = 0  <=>  T_B = V - T_A) disagreements: %d" % TOT['s6'])
print("")
if TOT['s1'] == TOT['s2'] == TOT['s4'] == TOT['s6'] == 0 and TOT['tr'] > 0 and TOT['rn'] > 0:
    print("  STEPS 5 AND 6 HOLD.  The chain now has exactly one unproved link, A1, and A1 is a")
    print("  statement about products of two Schur functions with no partitions and no roots of")
    print("  unity in it.  If A1 is in the literature, Phi_t = 0 => (ii) is a theorem.")
    for nm, n in (("alpha", TOT['ca']), ("beta", TOT['cb'])):
        if n == 0:
            print("")
            print("  NOTE, said rather than hidden: pairing %s never occurs in this range, so that" % nm)
            print("  branch of the step-5 proof is UNTESTED.  It is kept because the proof has to be")
            print("  exhaustive, not because a measurement supports it.  Only pairing %s is measured."
                  % ("beta" if nm == "alpha" else "alpha"))
else:
    print("  SOMETHING FAILED -- read the columns.")
print("")
print("DONE")
