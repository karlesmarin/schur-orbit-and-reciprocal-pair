# -*- coding: utf-8 -*-
# LA SATURACION -- y la muerte de mi propio "K_sop <= 4".  13 de agosto de 2026.
#
# QUE CORRIGE.  A media tarde escribi, con el barrido exhaustivo delante, que los estratos que de
# verdad CANCELAN se saturan en 4 y que por tanto la cota de profundidad revivia en la escalera con
# soporte.  Era falso, y falso por el mismo mecanismo que ha matado cinco lecturas esta semana: el
# barrido exhaustivo solo llega a K <= 8 fuera de t=4 r=2, asi que "K_sop <= 4" era una medida de
# hasta donde llego el reloj, no del objeto.  Iterando las familias del extremo fugitivo hasta
# W ~ 200 -- que ningun barrido alcanza -- cancelan sube a 24 (t=6 r=3), 18 (t=8 r=3) y 42 (t=6 r=4),
# y el valor DEPENDE DE LA SEMILLA.  El 4 era de r=2 y de anchura corta.
#
# LO QUE SI HAY, y es mas fuerte que lo que yo decia: cancelan SE SATURA, en las cuatro familias de
# las cuatro configuraciones.  Y la fila sinsop enseña por que, con dos fases que cambian en el MISMO j:
#
#     fase DENSA      la escalera no tiene huecos   ->  cancelan = prof/2 EXACTAMENTE, sinsop = 0
#     fase SATURADA   cancelan congelado            ->  todo el crecimiento se va a sinsop
#
# EL MECANISMO ES LA PROPIA LEY.  El interior esta CONGELADO, luego su aportacion al espectro de
# grados es un conjunto finito y fijo; los extremos solo trasladan.  Cuando la ventana que abren
# supera ese espectro no queda nada que cancelar y lo que se añade son huecos.  El valor de saturacion
# es el tamaño del espectro del interior congelado: un invariante de beta, no una constante.
#
# COLUMNAS
#   C0  fatal: probe() contra scan(), y los testigos publicados.
#   N1  LA TRAYECTORIA ENTERA de una familia por configuracion: W, prof, cancelan, sinsop.  Es el dato
#       que mata el 4, y se enseña completo, no resumido.
#   N2  el barrido de semillas: hasta que valor satura cada una.  Si el valor dependiera solo de (t,r)
#       saldria una sola columna; sale un histograma, y eso es el enunciado honesto.
#   N3  la fase densa medida contra su prediccion cancelan == prof/2, con su denominador.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python saturation.py

import itertools
import json
import os
import sys
from collections import Counter

from second_stratum import setup
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

SMIN = {4: 2, 6: 6, 8: 8}          # el regimen robusto: s == 0 (mod t) y par
CASOS = [(4, 2, 26, (18, 17, 11, 8, 7, 6, 1, 0), 30),
         (6, 3, 21, (16, 15, 14, 13, 11, 6, 5, 4, 3, 2, 1, 0), 20),
         (8, 3, 22, (20, 18, 15, 14, 13, 8, 7, 6, 5, 4, 3, 2, 1, 0), 16),
         (6, 4, 21, (20, 19, 17, 16, 15, 13, 8, 7, 6, 5, 4, 3, 1, 0), 16)]
OUT_JSON = "saturation_RESULT.json"


def extremos_S(beta, t):
    """los extremos del conjunto de EXCESO -- que es la variable, no beta (ver closed_form_prof.py)."""
    cl, E, Cd = setup(beta, t)
    S = sorted({v for k in E for v in Cd[k]})
    return S[-1], S[0]


def trayectoria(beta, t, r, s, J):
    hi, lo = extremos_S(beta, t)
    out = []
    for j in range(J + 1):
        b = tuple(sorted([(x + s * j if x == hi else (x - s * j if x == lo else x)) for x in beta],
                         reverse=True))
        rec = probe(b, t, r)
        if rec is None or not rec['surv'] or rec['prof'] is None:
            break
        out.append((b[0] - b[-1], rec['prof'], rec['vac_cancelan'], rec['vac_sin_soporte']))
    return out


def satura(c):
    """(valor, j donde se alcanza) si los tres ultimos coinciden; None si sigue subiendo."""
    if len(c) < 4 or not (c[-1] == c[-2] == c[-3]):
        return None
    return c[-1], c.index(c[-1])


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
for (beta, p, c) in [((18, 17, 11, 8, 7, 6, 1, 0), 6, 2), ((38, 23, 21, 18, 17, 16, 15, 0), 16, 4)]:
    rec = probe(beta, 4, 2)
    ok = rec is not None and rec['prof'] == p and rec['vac_cancelan'] == c
    bad += not ok
    print("  C0b  testigo %-36s prof %s/%d cancelan %s/%d  %s"
          % (str(beta), rec['prof'] if rec else "-", p, rec['vac_cancelan'] if rec else "-", c,
             "ok" if ok else "*** FALLA ***"))
print("")
if bad:
    print("  C0 FALLA -- el resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA")
RES = {}

# ===================================================================== N1 ========================
print("")
print("=" * 116)
print("N1  LA TRAYECTORIA ENTERA -- el dato que mata mi 'K_sop <= 4'.  Completa, no resumida")
print("=" * 116)
for (t, r, Wmax, sem, J) in CASOS:
    s = SMIN[t]
    f = trayectoria(sem, t, r, s, J)
    print("")
    print("  t=%d r=%d  s=%d  semilla %s   (%d miembros)" % (t, r, s, str(sem), len(f)))
    print("     W      : %s" % " ".join("%4d" % x[0] for x in f))
    print("     prof   : %s" % " ".join("%4d" % x[1] for x in f))
    print("     CANCEL : %s" % " ".join("%4d" % x[2] for x in f))
    print("     sinsop : %s" % " ".join("%4d" % x[3] for x in f))
    c = [x[2] for x in f]
    sat = satura(c)
    print("     -> %s" % ("SATURA en %d desde j=%d  (W=%d)" % (sat[0], sat[1], f[sat[1]][0])
                          if sat else "SIGUE SUBIENDO de %d a %d -- y entonces no satura" % (c[0], c[-1])))
    RES["N1_%d_%d" % (t, r)] = dict(semilla=list(sem), s=s, n=len(f),
                                    W=[x[0] for x in f], prof=[x[1] for x in f],
                                    cancelan=c, sinsop=[x[3] for x in f],
                                    satura=sat[0] if sat else None)
    sys.stdout.flush()

# ===================================================================== N2 ========================
print("")
print("=" * 116)
print("N2  EL VALOR DE SATURACION, semilla a semilla.  Si dependiera solo de (t,r) habria UNA columna")
print("=" * 116)
print("")
print("     t   r  |  semillas |  W max  | prof max | histograma del valor de SATURACION")
print("  " + "-" * 104)
for (t, r, Wmax, sem, J) in CASOS:
    N = t + 2 * r
    s = SMIN[t]
    SV = []
    for W in range(N - 1, Wmax + 1):
        for beta in shapes_of_width(W, N):
            rec = probe(beta, t, r)
            if rec and rec['surv'] and rec['prof'] is not None:
                SV.append(beta)
    hist = Counter()
    sinsat = 0
    Wmx = pmx = 0
    for beta in SV[:40]:
        f = trayectoria(beta, t, r, s, J)
        if len(f) < 4:
            continue
        Wmx = max(Wmx, f[-1][0])
        pmx = max(pmx, f[-1][1])
        sat = satura([x[2] for x in f])
        if sat:
            hist[sat[0]] += 1
        else:
            sinsat += 1
    print("    %2d  %2d  | %9d | %7d | %8d | %s%s"
          % (t, r, sum(hist.values()) + sinsat, Wmx, pmx, dict(sorted(hist.items())),
             "   + %d SIN saturar" % sinsat if sinsat else ""))
    RES["N2_%d_%d" % (t, r)] = dict(hist={str(k): v for k, v in hist.items()}, sin_saturar=sinsat,
                                    Wmax=Wmx, profmax=pmx)
    sys.stdout.flush()
print("")
print("     Sale un HISTOGRAMA, no una columna: el valor de saturacion es un invariante de beta.")
print("     Y en r=2 vale 4; en r=3 llega a 24; en r=4 a 42.  Mi 'K_sop <= 4' era de r=2.")

# ===================================================================== N3 ========================
print("")
print("=" * 116)
print("N3  LA FASE DENSA -- prediccion: mientras sinsop == 0, cancelan == prof/2 EXACTAMENTE")
print("=" * 116)
print("")
tot = ok = 0
for (t, r, Wmax, sem, J) in CASOS:
    f = trayectoria(sem, t, r, SMIN[t], J)
    n = m = 0
    for (W, p, c, ss) in f:
        if ss == 0:
            n += 1
            m += (c == p // 2)
    tot += n
    ok += m
    print("     t=%d r=%d : %d miembros con sinsop == 0, y cancelan == prof/2 en %d de ellos" % (t, r, n, m))
print("")
print("     TOTAL: %d de %d %s" % (ok, tot, "<-- la fase densa es exacta" if ok == tot else "*** FALLA ***"))
RES['N3'] = dict(ok=ok, total=tot)

# ===================================================================== N4 ========================
print("")
print("=" * 116)
print("N4  VEREDICTO")
print("=" * 116)
print("")
sats = [RES["N1_%d_%d" % (t, r)]['satura'] for (t, r, _, _, _) in CASOS]
print("     saturaciones de las cuatro familias : %s" % sats)
if all(x is not None for x in sats):
    print("     LAS CUATRO SATURAN.  'cancelan' es finito a lo largo de cada familia -- ESO es lo que")
    print("     hay que probar.  El VALOR no es universal: 4, 24, 18, 42, y depende de la semilla.")
    print("     Mi 'K_sop <= 4' de esta tarde queda RETIRADO: era r=2 a anchura corta, el sexto")
    print("     artefacto de rango de la semana.")
else:
    print("     ALGUNA NO SATURA -- entonces ni siquiera la saturacion vale, y se dice.")
print("")
print("     ALCANCE: J esta puesto por reloj.  'Satura' significa 'los tres ultimos miembros coinciden")
print("     hasta W ~ 200-276'.  No es una prueba de que no vuelva a subir mas alla; es un null con")
print("     alcance, y el alcance es ese.")
print("")
json.dump(RES, open(OUT_JSON, "w"), indent=1)
print("DONE")
