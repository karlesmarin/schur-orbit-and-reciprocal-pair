# -*- coding: utf-8 -*-
# LA PROFUNDIDAD EN FORMA CERRADA, y la hipotesis exacta.  13 de agosto de 2026.
#
# QUE ES.  Consolida en un solo guion, con controles, lo que el 13 salio en scripts sueltos: la
# profundidad de anulacion NO hace falta medirla expandiendo estratos.  Es la diferencia de DOS
# maximizaciones voraces, y la hipotesis bajo la que eso vale es intrinseca.
#
# LA VARIABLE, y por que es esa.  Toda transversal elige UN elemento por clase de residuo, luego una
# clase con un solo elemento OBLIGA a elegirlo y ese elemento nunca queda en T.  Entonces
#
#     S = union de las clases de EXCESO  =  exactamente los elementos que una transversal puede dejar
#
# y  C = min S + max S  no se elige: una involucion que preserva S manda el minimo al maximo.  Esto
# NO es beta_0 + beta_{N-1}: coinciden solo cuando todas las clases tienen dos elementos, que es el
# accidente de t=4 r=2 (N=8, t=4, sum(n_k - 1) = 2r = 4).  Medir C sobre beta en vez de sobre S fue
# el error que me dio diez residuos falsos en t=8 r=3, y N2 lo enseña con su denominador.
#
# LAS DOS FORMULAS
#   F1   D     = sum_k c_{k,1} + (los r mayores Delta_k(j) = c_{k,j} + c_{k,j+1}) - sum(S)
#        es el voraz separable-concavo de la prueba de |G| <= 2 (12 de agosto), reusado tal cual.
#   F2   first = max{ deg(T) : T conserva EXACTAMENTE UNO de max S, min S }
#        sale del lema de la anatomia (N3): el estrato de arriba lo soportan bloques con LOS DOS
#        extremos, el primer estrato vivo con UNO.
#   y entonces   prof = F1 - F2   sin expandir ningun polinomio.
#
# LA HIPOTESIS
#   H*  los DOS maximizadores conservan los DOS extremos de S.
#   Es la condicion (ii) mirada desde los extremos: por |G| <= 2 los dos maximizadores difieren SOLO
#   en las dos clases empatadas, que son las clases fijas de sigma_C; un maximizador pierde un extremo
#   si y solo si ese extremo vive en una de ellas.
#
# COLUMNAS
#   C0  fatal: probe() contra scan(), y las tres formas testigo ya publicadas.
#   N1  Q1 (T subset S) y Q2 (C - min S = max S): las dos razones de que la variable sea S.
#   N2  H*  <->  F2 = first  <->  LA LEY del extremo fugitivo.  Tabla de contingencia con las cuatro
#       casillas, y el numero que decide es el de las DOS casillas de desacuerdo.
#       Y el SEÑUELO: la misma tabla con C leida sobre beta en vez de sobre S, que es como la lei yo.
#   N3  la anatomia: quien soporta el estrato de arriba y quien el primer estrato vivo.
#   N4  los dos regimenes del paso s, que solo dependen de s mod t.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python closed_form_prof.py

import itertools
import json
import os
import sys
from collections import Counter, defaultdict

from second_stratum import setup, all_transversals
from depth_histogram import stratify, stratum
from survivors_wide import scan

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
assert "def probe(" in _head and "def shapes_of_width(" in _head, "k_vs_m.py cambio de forma"
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]
shapes_of_width = _ns["shapes_of_width"]

CFG = [(4, 2, 26), (6, 3, 21), (8, 3, 22)]
SMIN = {4: 2, 6: 6, 8: 8}          # el paso del regimen robusto: s == 0 (mod t), y par
OUT_JSON = "closed_form_prof_RESULT.json"


def excess(beta, t):
    """(cl, E, S, C) con S la union de las clases de EXCESO y C = min S + max S -- como en probe()."""
    cl, E, Cd = setup(beta, t)
    S = sorted({v for k in E for v in Cd[k]})
    return cl, E, S, S[0] + S[-1]


def F1(beta, t, r):
    """D por el voraz.  No mira ningun polinomio."""
    cl, E, S, C = excess(beta, t)
    by = defaultdict(list)
    for b in beta:
        by[b % t].append(b)
    for k in by:
        by[k].sort(reverse=True)
    inc = []
    for k, v in by.items():
        for j in range(1, len(v)):
            inc.append(v[j - 1] + v[j])
    inc.sort(reverse=True)
    return sum(v[0] for v in by.values()) + sum(inc[:r]) - sum(beta)


def F2(beta, t, r, sobre_beta=False):
    """first por el voraz restringido: transversales que conservan EXACTAMENTE UN extremo.
    sobre_beta=True usa los extremos de beta en vez de los de S -- el SEÑUELO, que es mi error."""
    cl, E, S, C = excess(beta, t)
    hi, lo = (beta[0], beta[-1]) if sobre_beta else (S[-1], S[0])
    best = None
    for (_, T, _, d) in all_transversals(beta, cl, r, t):
        if ((hi in T) + (lo in T)) == 1:
            best = d if best is None else max(best, d)
    return best


def Hstar(beta, t, r, sobre_beta=False):
    cl, E, S, C = excess(beta, t)
    hi, lo = (beta[0], beta[-1]) if sobre_beta else (S[-1], S[0])
    tr = all_transversals(beta, cl, r, t)
    D = max(x[3] for x in tr)
    return all((hi in x[1]) and (lo in x[1]) for x in tr if x[3] == D)


def medido(beta, t, r):
    """(D, first) por expansion -- la verdad contra la que se contrastan las formulas."""
    cl, E, Cd = setup(beta, t)
    tr = all_transversals(beta, cl, r, t)
    D = max(x[3] for x in tr)
    B = stratify([(x[2], x[1]) for x in tr], r)
    for s_ in sorted(B, reverse=True):
        if stratum(B[s_], r):
            return D, s_
    return D, None


def ley(beta, t, r, s, J=6):
    """el generador sobre los extremos de S: prof tiene que crecer exactamente s por paso."""
    cl, E, S, C = excess(beta, t)
    hi, lo = S[-1], S[0]
    p0 = probe(beta, t, r)['prof']
    for j in range(1, J + 1):
        b = tuple(sorted([(x + s * j if x == hi else (x - s * j if x == lo else x)) for x in beta],
                         reverse=True))
        rec = probe(b, t, r)
        if rec is None or not rec['surv'] or rec['prof'] is None or rec['prof'] != p0 + s * j:
            return False
    return True


def survivors(t, r, Wmax):
    N = t + 2 * r
    out = []
    for W in range(N - 1, Wmax + 1):
        for beta in shapes_of_width(W, N):
            rec = probe(beta, t, r)
            if rec and rec['surv'] and rec['prof'] is not None:
                out.append(beta)
    return out


# ===================================================================== C0 ========================
print("=" * 116)
print("C0  ACEPTACION -- fatal")
print("=" * 116)
print("")
bad = 0
for (t, r, M) in [(4, 2, 15), (6, 3, 18)]:
    n_ref, cont_ref, sv_ref = scan(t, r, M)
    mine, mb = Counter(), []
    for comb in itertools.combinations(range(M + 1), t + 2 * r):
        beta = tuple(sorted(comb, reverse=True))
        rec = probe(beta, t, r, deep=False)
        if rec is None:
            continue
        mine[(rec['e'] == t, rec['surv'])] += 1
        if rec['surv']:
            mb.append(beta)
    ok = (sum(mine.values()) == n_ref and mine == cont_ref
          and sorted(mb) == sorted(x['beta'] for x in sv_ref))
    bad += not ok
    print("  C0a  probe() == scan()  t=%d r=%d M=%d : %s" % (t, r, M, "ok" if ok else "*** FALLA ***"))
for (beta, p) in [((18, 17, 11, 8, 7, 6, 1, 0), 6), ((38, 23, 21, 18, 17, 16, 15, 0), 16),
                  ((114, 61, 59, 56, 55, 54, 53, 0), 54)]:
    rec = probe(beta, 4, 2)
    ok = rec is not None and rec['surv'] and rec['prof'] == p
    bad += not ok
    print("  C0b  testigo %-40s prof %s/%d  %s"
          % (str(beta), rec['prof'] if rec else "-", p, "ok" if ok else "*** FALLA ***"))
print("")
if bad:
    print("  C0 FALLA -- el resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA")
RES = {}

SV = {}
for (t, r, Wmax) in CFG:
    SV[(t, r)] = survivors(t, r, Wmax)

# ===================================================================== N1 ========================
print("")
print("=" * 116)
print("N1  POR QUE LA VARIABLE ES S -- las dos razones, medidas")
print("=" * 116)
print("")
print("     Q1  T subset S : una clase de UN elemento obliga a elegirlo, luego nunca queda en T")
print("     Q2  C - min S = max S : una involucion que preserva S manda el minimo al maximo")
print("")
for (t, r, Wmax) in CFG:
    sv = SV[(t, r)]
    q1 = q2 = coin = 0
    for beta in sv:
        cl, E, S, C = excess(beta, t)
        Ss = set(S)
        q1 += all(set(x[1]) <= Ss for x in all_transversals(beta, cl, r, t))
        q2 += (C - S[0] == S[-1])
        coin += (S[-1] == beta[0] and S[0] == beta[-1])
    print("     t=%d r=%d (%3d formas) : Q1 %d/%d, Q2 %d/%d   |   S == beta en %d (%.0f%%) <- el accidente"
          % (t, r, len(sv), q1, len(sv), q2, len(sv), coin, 100.0 * coin / len(sv)))
    RES["N1_%d_%d" % (t, r)] = dict(n=len(sv), q1=q1, q2=q2, S_eq_beta=coin)
sys.stdout.flush()

# ===================================================================== N2 ========================
print("")
print("=" * 116)
print("N2  H*  <->  F2 = first  <->  LA LEY.  Y el SEÑUELO: la misma tabla con C leida sobre beta")
print("=" * 116)
for (t, r, Wmax) in CFG:
    sv = SV[(t, r)]
    tab = defaultdict(int)
    tabS = defaultdict(int)
    f1ok = 0
    for beta in sv:
        D, first = medido(beta, t, r)
        f1ok += (F1(beta, t, r) == D)
        tab[(Hstar(beta, t, r), F2(beta, t, r) == first, ley(beta, t, r, SMIN[t]))] += 1
        # SEÑUELO, arreglado.  La primera version comparaba H*(beta) contra F2(beta)==first, o sea
        # las dos mitades leidas sobre beta: cuando las dos se equivocan a la vez sale 0 desacuerdos
        # y el control NO PUEDE FALLAR.  Lo que hay que contar es cuantas veces la formula leida
        # sobre beta acierta el first MEDIDO, entre las formas donde H* (la buena) se cumple.
        if Hstar(beta, t, r):
            tabS[F2(beta, t, r, True) == first] += 1
    print("")
    print("  t=%d r=%d  (%d supervivientes, paso s=%d)" % (t, r, len(sv), SMIN[t]))
    print("     F1 = D (el voraz) : %d de %d" % (f1ok, len(sv)))
    print("       H*  | F2=first | LEY | n")
    for k in sorted(tab, reverse=True):
        print("       %-3s | %-8s | %-3s | %d" % tuple(["SI" if x else "no" for x in k] + [tab[k]]))
    mal = sum(v for k, v in tab.items() if not (k[0] == k[1] == k[2]))
    print("     desacuerdos entre las TRES: %d %s"
          % (mal, "<-- LAS TRES SON LA MISMA COSA" if mal == 0 else "*** NO son la misma ***"))
    malS = tabS[False]
    print("     SEÑUELO: F2 con los extremos de BETA acierta el first medido en %d de %d formas con H*"
          % (tabS[True], tabS[True] + malS))
    print("              %s"
          % ("coincide aqui porque S == beta en todas" if malS == 0
             else "<-- FALLA en %d: AQUI ES DONDE MI LECTURA FALLABA" % malS))
    RES["N2_%d_%d" % (t, r)] = dict(n=len(sv), F1=f1ok, desacuerdos=mal, senuelo=malS,
                                    tabla={str(k): v for k, v in tab.items()})
    sys.stdout.flush()

# ===================================================================== N3 ========================
print("")
print("=" * 116)
print("N3  LA ANATOMIA -- quien soporta cada extremo del intervalo de estratos")
print("=" * 116)
print("")
print("     un bloque de Laplace tiene grado 2*sum_S(T) - sum(T).  Bajo eps_s, el desplazamiento del")
print("     grado depende SOLO de cuantos extremos de S conserva T:  dos -> +-2s o 0 ;  uno -> +-s .")
print("")
for (t, r, Wmax) in CFG:
    cnt_top = Counter()
    cnt_first = Counter()
    for beta in SV[(t, r)][:40]:
        cl, E, S, C = excess(beta, t)
        hi, lo = S[-1], S[0]
        tr = all_transversals(beta, cl, r, t)
        D = max(x[3] for x in tr)
        B = stratify([(x[2], x[1]) for x in tr], r)
        first = next((s_ for s_ in sorted(B, reverse=True) if stratum(B[s_], r)), None)
        if first is None:
            continue
        cnt_top[tuple(sorted(Counter((hi in x[1]) + (lo in x[1])
                                     for x in tr if x[3] == D).items()))] += 1
        tipos = {(hi in set(T)) + (lo in set(T)) for (_, T, _, _) in B[first]}
        cnt_first[tuple(sorted(tipos))] += 1
    print("     t=%d r=%d : estrato de ARRIBA, extremos conservados por los maximizadores : %s"
          % (t, r, dict(cnt_top)))
    print("               primer estrato VIVO, tipos de bloque presentes                 : %s"
          % dict(cnt_first))
    RES["N3_%d_%d" % (t, r)] = dict(top=str(dict(cnt_top)), first=str(dict(cnt_first)))
sys.stdout.flush()

# ===================================================================== N4 ========================
print("")
print("=" * 116)
print("N4  LOS DOS REGIMENES DEL PASO -- solo dependen de s mod t")
print("=" * 116)
print("")
t, r, Wmax = 8, 3, 22
sv = SV[(8, 3)]
print("     t=8 r=3, %d supervivientes.  t/2 = 4." % len(sv))
print("       s   | s mod t | que le hace a las clases          | cumplen la ley")
print("     " + "-" * 88)
N4 = {}
for s in (4, 8, 12, 16):
    n = sum(ley(beta, t, r, s) for beta in sv)
    N4[s] = n
    print("     %3d  |   %2d    | %-32s | %d de %d"
          % (s, s % t, "los extremos a su clase PAREJA" if s % t == t // 2
             else ("los extremos a su PROPIA clase" if s % t == 0 else "otra"), n, len(sv)))
RES['N4'] = N4
print("")
print("     el regimen robusto es s == 0 (mod t): preserva la particion en clases ENTERA.")
print("     para t=6 el otro esta PROHIBIDO -- t/2 = 3 es impar y s tiene que ser par, porque para")
print("     t par Phi_t es homogeneo mod 2 en el grado (probado el 12).  Por eso alli s_min = 6 = t.")

# ===================================================================== N5 ========================
print("")
print("=" * 116)
print("N5  VEREDICTO")
print("=" * 116)
print("")
tot = sum(RES["N2_%d_%d" % (t, r)]['n'] for (t, r, _) in CFG)
des = sum(RES["N2_%d_%d" % (t, r)]['desacuerdos'] for (t, r, _) in CFG)
f1 = sum(RES["N2_%d_%d" % (t, r)]['F1'] for (t, r, _) in CFG)
print("     F1 = D  (voraz, sin polinomios)          : %d de %d" % (f1, tot))
print("     H* <-> F2=first <-> LEY                  : %d desacuerdos en %d formas" % (des, tot))
print("")
if des == 0 and f1 == tot:
    print("     prof = F1 - F2 bajo H*.  La profundidad deja de costar una expansion polinomica y")
    print("     pasa a costar dos ordenaciones.  Y H* es la condicion (ii) del 12 mirada desde los")
    print("     extremos de S, no una hipotesis nueva.")
else:
    print("     ALGO NO CUADRA -- ver las tablas de N2, que es para lo que estan.")
print("")
print("     ALCANCE: los supervivientes salen de barridos exhaustivos con tope de RELOJ.  't=8 r=2'")
print("     y 't=10 r=2' no aparecen porque no tienen NINGUN superviviente hasta su tope, y eso es un")
print("     null con alcance, no un 'no existen'.")
print("")
json.dump(RES, open(OUT_JSON, "w"), indent=1)
print("DONE")
