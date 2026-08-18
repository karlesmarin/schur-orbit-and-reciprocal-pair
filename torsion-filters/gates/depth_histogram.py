# -*- coding: utf-8 -*-
# EL HISTOGRAMA DE PROFUNDIDADES SOBRE LOS 124 SUPERVIVIENTES.  13 de agosto de 2026.
#
# LA DEUDA.  survivors_wide.py se paro A MANO dentro de su tabla por superviviente (bloque t=6 r=4),
# y con ella se quedo sin imprimir el histograma de profundidades -- que es justo la pregunta abierta:
# la ruta de t >= 4 decia "Phi_t != 0 => alguno de los TRES PRIMEROS estratos es no nulo", ya sabemos
# que hay formas de profundidad 6 (depth6_check_OUT.txt), y falta saber HASTA DONDE BAJA.
#
# POR QUE SE PARO Y QUE SE HACE AQUI.  full_expansion recorre las (2r)! permutaciones POR TRANSVERSAL:
# con t=6 r=4 son 8! = 40 320 por cada una de ~130 transversales, o sea ~12 s por forma.  Pero para la
# profundidad no hace falta la expansion entera, solo el PRIMER ESTRATO NO NULO.  Y el grado total de
# un monomio depende UNICAMENTE de que r de los 2r indices van al lado positivo:
#
#     Phi_t = sum_T w(T) det M(T),   M(T)[a][2j] = x_j^{+T_a},  M(T)[a][2j+1] = x_j^{-T_a}
#
# y el desarrollo de Laplace por las columnas pares parte ese determinante en C(2r,r) bloques, uno por
# cada reparto S (los que van a "+"), cada uno de grado total FIJO  sum_{S} T - sum_{S^c} T.  Se agrupan
# los bloques de todas las transversales por grado, se baja estrato a estrato y SE PARA en el primero
# no nulo.  Con profundidad pequena se tocan 3-4 bloques en vez de los 70, y cada bloque cuesta r!+r!
# en vez de (2r)!.  Medido: 12 s -> 0,01 s por forma, factor ~1000.  No es una aproximacion: es la
# MISMA suma, reordenada -- y C0a lo comprueba diccionario a diccionario.
#
# COLUMNAS
#   C0  ACEPTACION, fatal.  Cuatro controles, todos pueden fallar:
#       C0a  la suma de TODOS los estratos == full_expansion, monomio a monomio, en formas al azar.
#       C0b  las filas ARCHIVADAS de survivors_wide_OUT.txt (D1, primer grado, profundidad) se
#            reproducen exactamente -- incluidas las de t=6 r=4, que son las que costaban minutos.
#       C0c  las 4 formas de depth6_check_OUT.txt: espectro de grados y estratos vacios por encima.
#       C0d  DECOY: formas con [Phi]_top != 0 tienen profundidad 0 por construccion.  Si el medidor
#            les encuentra profundidad, mide otra cosa.
#   N1  los 124 supervivientes, la tabla ENTERA que se quedo a medias.
#   N2  EL HISTOGRAMA de profundidades, global y por configuracion.  La deuda declarada.
#   N3  N4 del guion anterior, completado: Delta = 0 contra el POLINOMIO en los 124, no en 28.
#   N4  supervivientes con Phi_t == 0 (contraejemplo a la necesidad de (i), seria algo peor).
#   N5  VEREDICTO: el menor K con "alguno de los K primeros estratos es no nulo" sobre esta poblacion.
#
# DEFINICIONES IMPORTADAS, no reescritas: scan / full_expansion de survivors_wide.py (a quien se le
# anadio un guard __main__ -- cambio de pura indentacion, verificado por diff inverso), setup /
# all_transversals / perm_sign de second_stratum.py.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python depth_histogram.py

import json
import re
import random
import sys
import time
from collections import Counter, defaultdict
from itertools import combinations, permutations

from second_stratum import setup, all_transversals, perm_sign
from survivors_wide import scan, full_expansion, CONF_NEW

OUT_SURV = "depth_histogram_SURV.json"


# ===================================================================== EL MEDIDOR ================
def split_sign(S, Sc, r):
    """signo de la permutacion que manda S a las columnas pares y S^c a las impares."""
    q = [0] * (2 * r)
    for j in range(r):
        q[S[j]] = 2 * j
        q[Sc[j]] = 2 * j + 1
    return perm_sign(q)


def alt(vals, r):
    """el alternante: r! monomios con signo, clave = exponentes por variable."""
    return {tuple(vals[p[j]] for j in range(r)): perm_sign(p)
            for p in permutations(range(r))}


def stratify(tm, r):
    """{grado total: [(w, T, S, S^c)]} -- los C(2r,r) bloques de Laplace de cada transversal."""
    buckets = defaultdict(list)
    idx = tuple(range(2 * r))
    for w, T in tm:
        tot = sum(T)
        for S in combinations(idx, r):
            Sc = tuple(a for a in idx if a not in S)
            buckets[2 * sum(T[a] for a in S) - tot].append((w, T, S, Sc))
    return buckets


def stratum(bucket, r):
    """el estrato entero: {monomio: coeficiente}, ya sin los ceros."""
    acc = defaultdict(int)
    for (w, T, S, Sc) in bucket:
        A = alt([T[a] for a in S], r)
        B = alt([-T[a] for a in Sc], r)
        base = w * split_sign(S, Sc, r)
        for ka, ca in A.items():
            c = base * ca
            for kb, cb in B.items():
                acc[tuple(ka[j] + kb[j] for j in range(r))] += c * cb
    return {k: v for k, v in acc.items() if v}


def measure(tm, r, want=1):
    """None si Phi_t == 0; si no (primer grado no nulo, grados que SE CANCELAN por encima,
    grados no nulos hasta 'want', conjunto de grados con soporte).

    CONVENIO DE ESTRATO -- el de depth.py:205 y el de depth6_check, que es el de los enunciados:
    los estratos son D1, D1-2, D1-4, ... se cuenten o no.  Un grado por encima del primero no nulo
    puede estar vacio por DOS razones distintas, y aqui se separan porque no son lo mismo para
    quien quiera probar una cota: o tiene monomios que se cancelan, o no tiene monomio ninguno."""
    B = stratify(tm, r)
    first = None
    spectrum = []
    cancel = []
    for s in sorted(B, reverse=True):
        if stratum(B[s], r):
            if first is None:
                first = s
            spectrum.append(s)
            if len(spectrum) >= want:
                break
        elif first is None:
            cancel.append(s)
    if first is None:
        return None
    return first, cancel, spectrum, set(B)


def full_dict(tm, r):
    """la expansion entera POR ESTRATOS (para el control C0a)."""
    out = {}
    for s, bucket in stratify(tm, r).items():
        out.update(stratum(bucket, r))
    return out


def transversals_of(beta, t, r):
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(beta, cl, r, t)
    return tr, [(x[2], x[1]) for x in tr], max(x[3] for x in tr)


if __name__ == "__main__":
    # ===================================================================== C0 ========================
    print("=" * 108)
    print("C0  ACEPTACION -- cuatro controles, y el veredicto se suspende si falla cualquiera")
    print("=" * 108)
    print("")

    # ---- C0a  la suma de todos los estratos ES la expansion entera ---------------------------------
    print("  C0a  suma de TODOS los estratos  ==  full_expansion, monomio a monomio")
    random.seed(20260813)
    a_ok = a_bad = 0
    for (t, r, M, n) in [(4, 2, 14, 40), (4, 3, 15, 25), (6, 2, 13, 25), (6, 3, 15, 12)]:
        ok = bad = 0
        for _ in range(n):
            beta = tuple(sorted(random.sample(range(M + 1), t + 2 * r), reverse=True))
            got = transversals_of(beta, t, r)
            if got is None:
                continue
            _, tm, _ = got
            if full_expansion(tm, r) == full_dict(tm, r):
                ok += 1
            else:
                bad += 1
                print("       *** DISTINTOS *** t=%d r=%d beta=%s" % (t, r, list(beta)))
        print("       t=%d r=%d : %3d formas iguales, %d distintas" % (t, r, ok, bad))
        a_ok += ok
        a_bad += bad
        sys.stdout.flush()
    print("       C0a %s (%d formas)" % ("PASA" if a_bad == 0 and a_ok else "FALLA", a_ok))
    print("")

    # ---- C0b  las filas archivadas del guion interrumpido ------------------------------------------
    print("  C0b  las filas de survivors_wide_OUT.txt, reproducidas con el medidor rapido")
    ROW = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+"
                     r"\[([\d,\s]*)\]\s+(si|NO)\s+(\d+)\s+(\d+)\s+(\d+)\s*$")
    rows = []
    for line in open("survivors_wide_OUT.txt", encoding="utf-8"):
        m = ROW.match(line.rstrip("\n"))
        if m:
            rows.append(dict(t=int(m.group(1)), r=int(m.group(2)),
                             beta=tuple(int(v) for v in m.group(6).split(",")),
                             D1=int(m.group(8)), first=int(m.group(9)), prof=int(m.group(10))))
    print("       filas archivadas encontradas: %d" % len(rows))
    b_ok = b_bad = 0
    t0 = time.time()
    for row in rows:
        got = transversals_of(row['beta'], row['t'], row['r'])
        _, tm, D1 = got
        out = measure(tm, row['r'])
        first = None if out is None else out[0]
        if D1 == row['D1'] and first == row['first'] and D1 - first == row['prof']:
            b_ok += 1
        else:
            b_bad += 1
            print("       *** NO COINCIDE *** beta=%s  archivado D1=%d primer=%d prof=%d  medido D1=%d "
                  "primer=%s" % (list(row['beta']), row['D1'], row['first'], row['prof'], D1, first))
    by_cfg = Counter((x['t'], x['r']) for x in rows)
    print("       %d de %d reproducidas  (%s)  en %.1f s -- el guion viejo tardaba minutos por fila"
          % (b_ok, len(rows), ", ".join("t=%d r=%d: %d" % (k[0], k[1], v) for k, v in sorted(by_cfg.items())),
             time.time() - t0))
    print("       C0b %s" % ("PASA" if b_bad == 0 and b_ok else "FALLA"))
    print("")

    # ---- C0c  las 4 formas de profundidad 6, con su espectro ---------------------------------------
    print("  C0c  depth6_check_OUT.txt: espectro de grados y estratos vacios por encima")
    txt = open("depth6_check_OUT.txt", encoding="utf-8").read().split("\n")
    c_ok = c_bad = 0
    for i, line in enumerate(txt):
        m = re.match(r"^beta=\[([\d,\s]+)\].*Dmax=(\d+)\s+primer grado no nulo=(\d+)\s+PROFUNDIDAD=(\d+)",
                     line)
        if not m:
            continue
        beta = tuple(int(v) for v in m.group(1).split(","))
        Dmax, first_a, prof_a = int(m.group(2)), int(m.group(3)), int(m.group(4))
        spec_a = [int(v) for v in re.search(r"\[([\d,\s]+)\]", txt[i + 1]).group(1).split(",")]
        empty_a = [int(v) for v in re.search(r"\[([\d,\s]+)\]", txt[i + 2]).group(1).split(",")]
        _, tm, D1 = transversals_of(beta, 4, 2)
        first, cancel, spec, sup = measure(tm, 2, want=len(spec_a))
        empty = list(range(D1, first, -2))                 # el convenio de depth.py y depth6_check
        ok = (D1 == Dmax and first == first_a and D1 - first == prof_a
              and spec == spec_a and empty == empty_a)
        print("       beta=%-32s D1=%d primer=%d prof=%d  espectro %s  vacios %s   %s"
              % (str(list(beta)), D1, first, D1 - first, spec, empty, "ok" if ok else "*** FALLA ***"))
        print("            de esos vacios: %s se cancelan, %s no tienen ni un monomio"
              % ([d for d in empty if d in sup], [d for d in empty if d not in sup]))
        c_ok += ok
        c_bad += not ok
    print("       C0c %s (%d formas)" % ("PASA" if c_bad == 0 and c_ok else "FALLA", c_ok))
    print("")

    # ---- C0d  el senuelo: si [Phi]_top != 0, la profundidad es 0 -----------------------------------
    print("  C0d  DECOY: formas con [Phi]_top != 0 deben dar profundidad 0")
    from second_stratum import inv_of                                                    # noqa: E402
    d_n = d_bad = 0
    for comb in combinations(range(17), 8):
        beta = tuple(sorted(comb, reverse=True))
        got = transversals_of(beta, 4, 2)
        if got is None:
            continue
        tr, tm, D1 = got
        G = [x for x in tr if x[3] == D1]
        if len(G) == 2 and inv_of(G[0][1], 2) == inv_of(G[1][1], 2) and G[0][2] == -G[1][2]:
            continue                                     # [Phi]_top = 0, no es senuelo
        d_n += 1
        if d_n > 400:
            break
        out = measure(tm, 2)
        if out is None or out[0] != D1:
            d_bad += 1
    print("       %d formas con [Phi]_top != 0 : %d con profundidad != 0  (debe ser 0)"
          % (min(d_n, 400), d_bad))
    print("       C0d %s" % ("PASA" if d_bad == 0 and d_n else "FALLA"))
    print("")

    C0 = (a_bad == 0 and a_ok and b_bad == 0 and b_ok and c_bad == 0 and c_ok and d_bad == 0 and d_n)
    print("  C0 %s" % ("PASA -- el medidor rapido es la expansion entera, reordenada." if C0 else
                       "FALLA -- el medidor no es de fiar, el resto NO vale."))
    if not C0:
        print("DONE (veredicto suspendido)")
        raise SystemExit(1)

    # ===================================================================== N1 ========================
    print("")
    print("=" * 108)
    print("N1  LOS SUPERVIVIENTES, la tabla entera (survivors_wide.py se paro en la fila %d)" % len(rows))
    print("=" * 108)
    print("")
    surv = []
    t0 = time.time()
    for (t, r, M) in CONF_NEW:
        out = scan(t, r, M)
        if out is None:
            continue
        n, cont, sv = out
        surv += sv
        print("  barrido t=%d r=%d M=%d : %d objetivo, %d supervivientes   [%.0f s]"
              % (t, r, M, n, len(sv), time.time() - t0), flush=True)
    print("")
    print("  supervivientes totales: %d" % len(surv))
    print("")
    print("     t   r  e |G2|  beta                                                 poly=0   D1   1er "
          " prof  estrato   vacios: cancelan / sin soporte")
    print("  " + "-" * 104)
    n_zero = n_poly_bad = n_odd = 0
    recs = []
    t0 = time.time()
    for s in sorted(surv, key=lambda x: (x['t'], x['r'], x['beta'])):
        r = s['r']
        if not s['poly_zero']:
            n_poly_bad += 1
        out = measure(s['tm'], r)
        if out is None:
            n_zero += 1
            prof = k = n_can = n_abs = None
            off_abs = []
            line = "   Phi_t == 0 (!!)"
        else:
            first, cancel, _, sup = out
            prof = s['D1'] - first
            empty = list(range(s['D1'], first, -2))        # convenio de depth.py: paso 2, se cuenten o no
            n_can = sum(1 for d in empty if d in sup)
            n_abs = len(empty) - n_can
            off_abs = [s['D1'] - d for d in empty if d not in sup]      # a que altura esta el sin-soporte
            n_odd += (prof % 2)
            k = prof // 2 + 1
            line = "%5d %5d %8d %11d / %d" % (first, prof, k, n_can, n_abs)
        print("  %4d %3d %2d %4d  %-52s %6s %4d %s"
              % (s['t'], r, s['e'], s['nG2'], str(list(s['beta'])),
                 "si" if s['poly_zero'] else "NO", s['D1'], line))
        recs.append(dict(t=s['t'], r=r, beta=list(s['beta']), e=s['e'], nG2=s['nG2'],
                         poly_zero=s['poly_zero'], D1=s['D1'], prof=prof, estrato=k,
                         vac_cancelan=n_can, vac_sin_soporte=n_abs, vac_altura=off_abs))
        sys.stdout.flush()
    print("")
    print("  los %d supervivientes medidos en %.1f s" % (len(surv), time.time() - t0))

    json.dump(recs, open(OUT_SURV, "w", encoding="utf-8"), indent=1)
    print("  archivados en %s" % OUT_SURV)

    # ===================================================================== N2 ========================
    print("")
    print("=" * 108)
    print("N2  EL HISTOGRAMA DE PROFUNDIDADES  (la deuda declarada)")
    print("=" * 108)
    print("")
    prof = Counter(x['prof'] for x in recs if x['prof'] is not None)
    estr = Counter(x['estrato'] for x in recs if x['estrato'] is not None)
    tot = sum(prof.values())
    print("     profundidad D1 - (primer grado no nulo)   |   estrato no nulo (1 = el de arriba)")
    print("  " + "-" * 88)
    for k in sorted(prof):
        print("     %2d : %4d formas  %-30s |   #%d" % (k, prof[k], "#" * min(prof[k], 30), k // 2 + 1))
    print("")
    print("  parity check: supervivientes con profundidad IMPAR (los grados deberian ir de 2 en 2): %d"
          % n_odd)
    print("")
    print("  los estratos VACIOS por encima del primero no nulo, separados por causa:")
    print("     (un grado sin ni un monomio NO se puede probar no nulo -- para una cota por estratos")
    print("      no es lo mismo que uno donde hay cancelacion)")
    vc = Counter(x['vac_cancelan'] for x in recs if x['prof'] is not None)
    va = Counter(x['vac_sin_soporte'] for x in recs if x['prof'] is not None)
    print("     %-28s %s" % ("vacios POR CANCELACION:", dict(sorted(vc.items()))))
    print("     %-28s %s" % ("vacios SIN SOPORTE:", dict(sorted(va.items()))))
    print("")
    print("  por configuracion:")
    print("     t   r | %s" % "  ".join("prof %2d" % k for k in sorted(prof)))
    print("  " + "-" * 88)
    for cfg in sorted({(x['t'], x['r']) for x in recs}):
        c = Counter(x['prof'] for x in recs if (x['t'], x['r']) == cfg and x['prof'] is not None)
        print("  %4d %3d | %s" % (cfg[0], cfg[1],
                                 "  ".join("%7d" % c[k] for k in sorted(prof))))

    # ===================================================================== N3-N5 =====================
    print("")
    print("=" * 108)
    print("N3-N5  lo que el guion interrumpido no llego a decir")
    print("=" * 108)
    print("")
    print("  N3  Delta = 0 pero POLINOMIO no nulo (no serian supervivientes): %d de %d"
          % (n_poly_bad, len(surv)))
    print("  N4  supervivientes con Phi_t == 0 (contraejemplo a la necesidad de (i)): %d" % n_zero)
    print("")
    Kmax = max(estr) if estr else 0
    Pmax = max(prof) if prof else 0
    print("  N5  VEREDICTO -- hasta donde baja")
    print("      profundidad maxima observada : %d   (el enunciado de la noche del 12 decia <= 4)" % Pmax)
    print("      estrato no nulo mas profundo : %d   (el enunciado decia: alguno de los TRES primeros)"
          % Kmax)
    print("      formas por debajo del tercer estrato : %d de %d (%.0f%%)"
          % (sum(v for k, v in estr.items() if k > 3), tot,
             100.0 * sum(v for k, v in estr.items() if k > 3) / max(tot, 1)))
    print("")
    if Kmax > 3:
        print("      'Phi_t != 0 => alguno de los TRES PRIMEROS estratos es no nulo' ES FALSO en esta")
        print("      poblacion, y el menor K que la aguanta ENTERA es K = %d.  Pero K se ha movido con" % Kmax)
        print("      el rango en cada re-medida (4 artefactos de rango en dos dias), asi que K = %d es" % Kmax)
        print("      una COTA OBSERVADA A ESTE M, no un enunciado.")
    else:
        print("      el enunciado de los tres estratos SOBREVIVE a esta poblacion.")

    # ---- N6: el mismo veredicto en el OTRO convenio ------------------------------------------------
    print("")
    print("  N6  EL MISMO HISTOGRAMA CONTANDO SOLO ESTRATOS POBLADOS")
    print("      El convenio de arriba (el de depth.py y el de los enunciados) cuenta D1, D1-2, D1-4...")
    print("      se cuenten o no.  Pero un grado SIN NI UN MONOMIO no es un estrato de Phi_t que se")
    print("      cancele: no existe.  Contando solo los grados con soporte, el estrato no nulo es:")
    pob = Counter(x['vac_cancelan'] + 1 for x in recs if x['prof'] is not None)
    for k in sorted(pob):
        print("        estrato poblado #%d : %4d formas  %s" % (k, pob[k], "#" * min(pob[k], 30)))
    deep2 = [x for x in recs if x['vac_cancelan'] is not None and x['vac_cancelan'] + 1 > 3]
    print("")
    print("      EL ENUNCIADO CAE EN LOS DOS CONVENIOS, pero no con la misma poblacion:")
    print("        en grados            : %d de %d formas por debajo del tercer estrato"
          % (sum(v for k, v in estr.items() if k > 3), tot))
    print("        en estratos poblados : %d de %d" % (len(deep2), tot))
    print("      Las %d formas que lo rompen en el convenio ESTRICTO -- las mas profundas que tenemos:"
          % len(deep2))
    for x in deep2:
        print("        t=%d r=%d beta=%s  D1=%d prof=%d" % (x['t'], x['r'], x['beta'], x['D1'], x['prof']))
    print("")
    n6 = sum(1 for x in recs if x['prof'] == 6)
    sin = [x for x in recs if x['prof'] == 6 and x['vac_sin_soporte'] == 1]
    alturas = Counter(tuple(x['vac_altura']) for x in sin)
    print("      Y la anatomia de la profundidad 6: de las %d formas, %d tienen uno de los tres estratos"
          % (n6, len(sin)))
    print("      vacio PORQUE NO TIENE NI UN MONOMIO -- no porque se cancele -- y esta a altura D1-%s."
          % "/D1-".join(str(a[0]) for a in sorted(alturas)))
    print("      Solo las %d restantes cancelan de verdad los tres pisos." % (n6 - len(sin)))
    print("")
    print("DONE")
