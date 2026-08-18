# -*- coding: utf-8 -*-
# SE MUEVE K CON M?  El barrido que hasta ahora no se podia hacer.  13 de agosto de 2026.
#
# POR QUE.  En dos dias han muerto CUATRO lecturas, las cuatro por el mismo mecanismo: medidas a M
# corto leidas como firma estructural.  La ultima es la profundidad: <= 4 a M <= 19, y a M = 19..21
# aparece 6.  El histograma de depth_histogram.py dice profundidad 4 en 108 formas y 6 en 16, nada
# mas abajo -- pero A ESE M.  La pregunta que decide si "K = 4" es un enunciado o el quinto artefacto
# de rango es si K se mueve al subir M.  Antes costaba 12 s por forma y no se podia preguntar; con el
# medidor por estratos cuesta 0,01 s.
#
# LA VARIABLE NO ES M, ES LA ANCHURA.  Todo lo que se mide es INVARIANTE POR TRASLACION de beta:
# deg_of es una diferencia de sumas de r elementos cada una, inv_of y los monomios T_a - T_b son
# diferencias, y la condicion (i) es concentricidad.  Lo unico que se mueve es w, y solo por un signo
# GLOBAL (los residuos rotan, y una rotacion de Z_t es una permutacion fija): un signo global no toca
# |G|, ni a[2] == -b[2], ni Delta = 0, ni la profundidad.  Luego barrer subconjuntos de {0..M} contaba
# cada forma M+1-W veces.  Aqui se barre por ANCHURA W = max beta - min beta, una vez cada forma, y
#
#     K(M) = max_{W <= M} K(W)        exacto, sin re-barrer
#
# asi que un solo barrido en W da la curva K(M) ENTERA.  Y el numero de formas de anchura W es
# C(W-1, N-2) en vez de C(M+1, N): el mismo alcance sale mucho mas barato.
#
# COLUMNAS
#   C0  ACEPTACION, fatal, dos controles independientes:
#       C0a  mi probe(beta) reproduce scan() de survivors_wide.py EXACTO -- objetivo, contingencia
#            y la lista de betas supervivientes -- en sus propias configuraciones.  Es el control de
#            que no he reescrito el criterio al reimplementar el bucle.
#       C0b  sum_W (M+1-W) * n_W  ==  el n de scan(t,r,M), objetivo Y supervivientes.  Es el control
#            de la invariancia por traslacion Y del recuento por anchura, a la vez.  Si la
#            invariancia fuera falsa, esta suma no cuadraria.
#   N1  por configuracion: W -> objetivo, supervivientes, histograma de profundidad, K(W), K_sop(W).
#       SE IMPRIME AL TERMINAR CADA W.  La corrida anterior se paro a mano y perdio la tabla; esta
#       deja datos validos aunque se pare.
#
#       *** LA COLUMNA K_sop SE AÑADIO EL 13 DE AGOSTO POR LA TARDE, Y ESTE ES EL ARREGLO DE UN
#       DEFECTO DE ESTE MISMO GUION. ***  probe() calcula desde el principio 'vac_cancelan' y
#       'vac_sin_soporte' -- los peldanos vacios que CANCELAN y los que no tienen ni un monomio -- y
#       esta tabla imprimia SOLO prof, que los suma.  La cabecera de measure() en depth_histogram.py
#       ya avisaba de que "no son lo mismo para quien quiera probar una cota".  Con la columna oculta,
#       la corrida de la madrugada leyo la recta K = W/2 - 3 como "no hay ninguna cota, se cae el
#       instrumento entero", cuando tres cuartas partes de ese crecimiento son HUECOS DE ESPECTRO.
#       El dato estaba calculado y archivado; no estaba a la vista.  Un instrumento que calcula una
#       distincion y no la enseña induce el error que su propio autor documento.
#
#       SALIDAS: k_vs_m_OUT.txt es la corrida parada a mano del 13 de madrugada, y k_vs_m_OUT2.txt la
#       completa de la manana -- LAS DOS de la version SIN esta columna.  k_vs_m_OUT3.txt es la de
#       esta version.  No se sobrescriben: son salidas de guiones distintos.
#   N2  K(M) acumulado, que es la respuesta.
#   N3  la anchura MINIMA en que aparece cada profundidad -- si 6 aparece en W=19 y 8 no aparece
#       hasta el tope, el tope es la respuesta honesta, no "no existe".
#
# DEFINICIONES IMPORTADAS, no reescritas: setup / all_transversals / inv_of de second_stratum.py,
# dim_gl / halves de dim_certificate.py, P_poly de second_vanishes.py, el medidor por estratos
# (stratify / stratum / measure) de depth_histogram.py, y scan de survivors_wide.py para C0.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python k_vs_m.py

import itertools
import json
import sys
import time
from collections import Counter, defaultdict

from second_stratum import setup, all_transversals, inv_of
from second_vanishes import P_poly
from dim_certificate import dim_gl, halves
from depth_histogram import measure
from survivors_wide import scan

# (t, r, W_max).  El tope de cada una esta puesto por reloj, no por teoria, y se dice en la salida.
SWEEP = [(4, 2, 36), (6, 2, 26), (8, 2, 24), (6, 3, 24), (8, 3, 24), (6, 4, 23), (10, 2, 23)]
OUT_JSON = "k_vs_m_SWEEP.json"


def probe(beta, t, r, deep=True):
    """None si no es poblacion objetivo; si no, el registro con la profundidad ya medida.
    Mismo criterio que survivors_wide.scan, y C0a lo comprueba contra el."""
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(beta, cl, r, t)
    D = max(x[3] for x in tr)
    G = [x for x in tr if x[3] == D]
    if len(G) != 2:
        return None
    a, b = G
    if not (inv_of(a[1], r) == inv_of(b[1], r) and a[2] == -b[2]):
        return None                                   # [Phi]_top != 0
    S = sorted({v for k in E for v in Cd[k]})
    C = S[0] + S[-1]
    if set(C - v for v in S) == set(S):
        return None                                   # (i) cierta
    rest = [x for x in tr if x[3] < D]
    if not rest:
        return None
    D2 = max(x[3] for x in rest)
    G2 = [x for x in rest if x[3] == D2]
    delta = 0
    for (_, T, w, _) in G2:
        at, ast = halves(T, r)
        delta += w * dim_gl(at) * dim_gl(ast)
    rec = dict(e=len(E), nG2=len(G2), D1=D, surv=(delta == 0))
    if not rec['surv'] or not deep:
        return rec
    acc = defaultdict(int)
    for (_, T, w, _) in G2:
        for k, v in P_poly(T, r).items():
            acc[k] += w * v
    rec['poly_zero'] = not any(acc.values())
    out = measure([(x[2], x[1]) for x in tr], r)
    if out is None:
        rec['prof'] = None                            # Phi_t == 0
    else:
        first, cancel, _, sup = out
        rec['prof'] = D - first
        empty = list(range(D, first, -2))
        rec['vac_cancelan'] = sum(1 for d in empty if d in sup)
        rec['vac_sin_soporte'] = len(empty) - rec['vac_cancelan']
    return rec


def shapes_of_width(W, N):
    """las formas de anchura exacta W: contienen 0 y W.  C(W-1, N-2) de ellas."""
    for mid in itertools.combinations(range(1, W), N - 2):
        yield (W,) + tuple(reversed(mid)) + (0,)


# ===================================================================== C0 ========================
print("=" * 108)
print("C0  ACEPTACION -- dos controles independientes, y el veredicto se suspende si falla uno")
print("=" * 108)
print("")
print("  C0a  probe() reproduce scan() de survivors_wide.py: objetivo, contingencia y las BETAS")
a_bad = 0
for (t, r, M) in [(4, 2, 15), (4, 3, 15), (6, 2, 17), (6, 3, 18)]:
    n_ref, cont_ref, sv_ref = scan(t, r, M)
    mine = Counter()
    mine_beta = []
    for comb in itertools.combinations(range(M + 1), t + 2 * r):
        beta = tuple(sorted(comb, reverse=True))
        rec = probe(beta, t, r, deep=False)
        if rec is None:
            continue
        mine[(rec['e'] == t, rec['surv'])] += 1
        if rec['surv']:
            mine_beta.append(beta)
    ok = (sum(mine.values()) == n_ref and mine == cont_ref
          and sorted(mine_beta) == sorted(x['beta'] for x in sv_ref))
    a_bad += not ok
    print("       t=%d r=%d M=%d : objetivo %d/%d, contingencia %s, betas %d/%d   %s"
          % (t, r, M, sum(mine.values()), n_ref, "igual" if mine == cont_ref else "DISTINTA",
             len(mine_beta), len(sv_ref), "ok" if ok else "*** FALLA ***"))
    sys.stdout.flush()
print("       C0a %s" % ("PASA" if not a_bad else "FALLA"))
print("")

print("  C0b  sum_W (M+1-W) * n_W  ==  n de scan(t,r,M)   [invariancia por traslacion + recuento]")
b_bad = 0
for (t, r, M) in [(4, 2, 19), (6, 3, 19), (8, 3, 21)]:
    N = t + 2 * r
    n_ref, _, sv_ref = scan(t, r, M)
    tot = surv = 0
    for W in range(N - 1, M + 1):
        nw = sw = 0
        for beta in shapes_of_width(W, N):
            rec = probe(beta, t, r, deep=False)
            if rec is None:
                continue
            nw += 1
            sw += rec['surv']
        tot += (M + 1 - W) * nw
        surv += (M + 1 - W) * sw
    ok = (tot == n_ref and surv == len(sv_ref))
    b_bad += not ok
    print("       t=%d r=%d M=%d : objetivo %d (scan %d), supervivientes %d (scan %d)   %s"
          % (t, r, M, tot, n_ref, surv, len(sv_ref), "ok" if ok else "*** FALLA ***"))
    sys.stdout.flush()
print("       C0b %s" % ("PASA" if not b_bad else "FALLA"))
print("")
if a_bad or b_bad:
    print("  C0 FALLA -- el barrido por anchura no mide lo mismo que scan.  El resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA -- barrer por anchura es barrer por M, deduplicado.")

# ===================================================================== N1 ========================
print("")
print("=" * 108)
print("N1  EL BARRIDO.  Cada linea se imprime al terminar su W: una parada a mano deja datos validos")
print("=" * 108)
DATA = defaultdict(dict)
for (t, r, Wmax) in SWEEP:
    N = t + 2 * r
    print("")
    print("  t=%d r=%d  (N=%d, tope W=%d puesto POR RELOJ, no por teoria)" % (t, r, N, Wmax))
    print("     W    formas  objetivo  surv |  histograma de profundidad        K(W)  K(<=W) | K_sop  los que CANCELAN")
    print("  " + "-" * 122)
    Kacc = 0
    Kacc_sop = 0
    t0 = time.time()
    for W in range(N - 1, Wmax + 1):
        n_sh = n_tg = 0
        prof = Counter()
        extra = Counter()
        can = Counter()          # histograma de vac_cancelan: los peldanos vacios que SI cancelan
        for beta in shapes_of_width(W, N):
            n_sh += 1
            rec = probe(beta, t, r)
            if rec is None:
                continue
            n_tg += 1
            if not rec['surv']:
                continue
            prof[rec['prof']] += 1
            if rec['prof'] is not None:
                extra['poly_bad'] += not rec['poly_zero']
                extra['sin_sop'] += rec['vac_sin_soporte'] > 0
                can[rec['vac_cancelan']] += 1
            else:
                extra['phi_zero'] += 1
        ns = sum(prof.values())
        K = max((k for k in prof if k is not None), default=0)
        K_sop = max(can, default=0)            # la distincion que este guion calculaba y no enseñaba
        nuevo = K > Kacc                      # una profundidad que no se habia visto a menor anchura
        nuevo_sop = K_sop > Kacc_sop
        Kacc = max(Kacc, K)
        Kacc_sop = max(Kacc_sop, K_sop)
        DATA[(t, r)][W] = dict(formas=n_sh, objetivo=n_tg, surv=ns,
                               prof={str(k): v for k, v in prof.items()},
                               K=K, K_sop=K_sop, can={str(k): v for k, v in can.items()},
                               extra=dict(extra))
        hist = " ".join("p%s:%d" % (k, v) for k, v in
                        sorted(prof.items(), key=lambda x: (x[0] is None, x[0]))) or "-"
        hcan = " ".join("c%d:%d" % (k, v) for k, v in sorted(can.items())) or "-"
        print("  %4d %9d %9d %5d |  %-32s %4d %6d | %5d  %-26s%s%s"
              % (W, n_sh, n_tg, ns, hist, K, Kacc, K_sop, hcan,
                 "   <-- K SUBE" if nuevo else "",
                 "  <-- K_sop SUBE" if nuevo_sop else ""))
        sys.stdout.flush()
    print("     %d anchuras en %.0f s" % (Wmax - N + 2, time.time() - t0))

json.dump({"%d_%d" % k: v for k, v in DATA.items()}, open(OUT_JSON, "w"), indent=1)
print("")
print("  barrido archivado en %s" % OUT_JSON)

# ===================================================================== N2, N3 ====================
print("")
print("=" * 108)
print("N2  K(M) -- LA RESPUESTA.  K(M) = max_{W <= M} K(W), exacto")
print("=" * 108)
print("")
allW = sorted({W for c in DATA for W in DATA[c]})
print("     M   | %s" % "  ".join("t=%d r=%d" % (t, r) for (t, r, _) in SWEEP))
print("  " + "-" * 96)
prev = None
for M in allW:
    row = []
    for (t, r, Wmax) in SWEEP:
        d = DATA[(t, r)]
        ks = [d[w]['K'] for w in d if w <= M]
        row.append("%7s" % (max(ks) if ks else "-") if M <= Wmax else "%7s" % ".")
    line = "  %4d  | %s" % (M, "  ".join(row))
    if line.split("|")[1] != (prev or ""):
        print(line + "   <-- cambia")
    else:
        print(line)
    prev = line.split("|")[1]

print("")
print("=" * 108)
print("N3  la anchura MINIMA en que aparece cada profundidad")
print("=" * 108)
print("")
first_at = {}
for (t, r, Wmax) in SWEEP:
    for W in sorted(DATA[(t, r)]):
        for k in DATA[(t, r)][W]['prof']:
            if k == "None":
                continue
            key = int(k)
            if key not in first_at or W < first_at[key][0]:
                first_at[key] = (W, t, r)
print("     profundidad | primera anchura | donde")
print("  " + "-" * 60)
for k in sorted(first_at):
    W, t, r = first_at[k]
    print("     %11d | %15d | t=%d r=%d" % (k, W, t, r))
print("")
Kglob = max(first_at) if first_at else 0
topes = ", ".join("t=%d r=%d: %d" % (t, r, W) for (t, r, W) in SWEEP)
print("  VEREDICTO")
print("     K maximo en TODO el barrido : %d" % Kglob)
Ksop_glob = max((DATA[(t, r)][W].get('K_sop', 0) for (t, r, _) in SWEEP for W in DATA[(t, r)]),
                default=0)
print("     K_sop maximo (los que CANCELAN) : %d" % Ksop_glob)
print("     *** LEER LAS DOS COLUMNAS, NO SOLO K. ***  K cuenta peldanos de la escalera COMPLETA")
print("     (D, D-2, D-4, ... se cuenten o no) y K_sop solo los que tienen soporte y aun asi se")
print("     anulan.  Un K grande con K_sop pequeño NO es profundidad de cancelacion: es que el")
print("     espectro se ahueca al ensanchar beta.  La corrida del 13 de madrugada leyo solo K y")
print("     concluyo 'no hay ninguna cota, se cae el instrumento entero'.  Era cierto sobre K y")
print("     falso sobre el objeto.  Ver LAW_RUNAWAY_EXTREMES.md.")
print("")
print("     profundidad 6 aparece por primera vez en W = %s"
      % (first_at.get(6, ("-",))[0] if 6 in first_at else "no aparece"))
print("     profundidad 8 : %s"
      % ("APARECE en W = %d (t=%d r=%d) -- K SE MUEVE CON M, y K=4 era el quinto artefacto de rango"
         % first_at[8] if 8 in first_at else
         "NO APARECE en ningun sitio del barrido"))
print("")
print("     Los topes son de RELOJ: %s." % topes)
print("     Un 'no aparece' hasta el tope NO es 'no existe': es un null con alcance, y el alcance")
print("     es ese.  Las cuatro lecturas muertas en dos dias murieron entre una y dos tallas por")
print("     encima de donde se habian medido.")
print("")
print("DONE")
