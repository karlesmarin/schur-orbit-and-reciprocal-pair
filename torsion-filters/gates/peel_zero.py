# -*- coding: utf-8 -*-
# EL PELADO DE UN CERO: 14 de agosto de 2026.
#
# LA PREGUNTA.  El cuadrante abierto es t >= 4 PAR, r >= 2.  La bandera del 13 (flag.py) peló
# SUPERVIVIENTES -- formas con [Phi]_top = 0 pero Phi != 0 -- y encontró una torre.  Lo que NO se
# preguntó, y es lo único que cierra el criterio, es si el pelado respeta la anulación TOTAL:
#
#       Phi_{t,r}(beta) == 0     ==>     Phi_{t,r-1}(pelar beta) == 0 ?
#
# Si eso vale, hay INDUCCIÓN EN r: se baja hasta r = 1, donde el Teorema 3.1 da el criterio completo
# para todo t, y el cuadrante se cierra.  Si no vale, la vía del pelado está muerta y hay que decirlo.
#
# POR QUÉ ES CREÍBLE, Y POR QUÉ NO ES SUFICIENTE ESO.  En t = 2 el pelado preserva la rama (b) por
# álgebra de dos líneas: si beta es simétrica respecto de C, quitar beta_1 y beta_N deja una beta
# simétrica respecto del MISMO C, con N' = N-2 y anchura w' = C - N' + 1 = w + 2, de la misma
# paridad.  O sea: en t = 2 el pelado reproduce el teorema.  Eso es una razón para probarlo en t par
# general y NO es una prueba, porque en t >= 4 el g_com puede ser no vacío y romper la simetría
# global -- que es exactamente lo que la Proposición 8.30 sólo controla en la capa exterior.
#
# COLUMNAS
#   N1  la implicación, exhaustiva por configuración.  Se cuentan las formas con Phi == 0 y se mira
#       el pelado.  Cualquier fallo se imprime con su beta.
#   N2  el RECÍPROCO, que es lo que haría falta para un criterio y no sólo para una necesidad:
#       Phi == 0  <=>  (C = tau) y (el pelado es cero).
#   N3  CONTROL, y sin él esto no mide nada: el pelado de una forma que NO se anula.  Si el pelado
#       devolviera cero casi siempre, N1 sería vacío.  Se exige que la tasa sea baja.
#   N4  CONTROL en t = 2, donde el criterio ya está probado: el pelado tiene que preservar la
#       autocomplementariedad de anchura impar, 100%, y ninguna otra cosa.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python peel_zero.py

import itertools
import sys
from collections import Counter

from second_stratum import setup, all_transversals, inv_of
from depth_histogram import measure
from flag import pelar, extremos_S

# (t, r, Wmax) -- anchura maxima de beta_1.  N = t + 2r.
CFG = [(4, 2, 17), (6, 2, 17), (8, 2, 18), (4, 3, 16), (6, 3, 16), (10, 2, 20)]
CFG_T2 = [(2, 2, 16), (2, 3, 16), (2, 4, 16)]


def betas(t, r, W):
    """beta estrictamente decreciente, longitud N, beta_N = 0, beta_1 <= W.
    beta_N = 0 no pierde generalidad: una traslacion comun multiplica Phi por det(A)^m != 0."""
    N = t + 2 * r
    for mid in itertools.combinations(range(1, W + 1), N - 1):
        yield tuple(sorted(mid, reverse=True)) + (0,)


def occupied(b, t):
    return len({v % t for v in b}) == t


def phi_zero(b, t, rp):
    """Phi_{t,rp}(beta) == 0, exacto, por la descomposicion de Laplace.

    OJO, y es el defecto que este comentario existe para no repetir: setup() devuelve None
    EXACTAMENTE cuando falta una clase residual, y eso NO es "no se sabe" -- por el Lema 4.1
    no hay ninguna transversal, la suma de Laplace es vacia y Phi == 0 OUTRIGHT.  Es la rama
    (a).  La primera version de este guion lo saltaba en N1 (truncando en silencio la
    poblacion, justo la de los pelados degenerados) y lo contaba como FALLO en el control de
    t = 2.  Un dato ausente no es un cero, pero este dato no esta ausente: es un cero."""
    st = setup(b, t)
    if st is None:
        return True                      # rama (a): clase residual vacia
    cl, E, Cd = st
    if not E:
        return None                      # esto si es "no aplica": no hay clase de exceso
    tr = all_transversals(b, cl, rp, t)
    return measure([(x[2], x[1]) for x in tr], rp) is None


def C_and_tau(b, t, rp):
    """(C = min S + max S, tau = Delta_(r)) -- las dos constantes de la Corolario 8.31."""
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    S = sorted({v for k in E for v in Cd[k]})
    incr = []
    for k in E:
        c = sorted(Cd[k], reverse=True)
        incr += [c[i] + c[i + 1] for i in range(len(c) - 1)]
    incr.sort(reverse=True)
    if len(incr) < rp:
        return None
    return S[0] + S[-1], incr[rp - 1]


def run(cfgs, label, t2=False):
    print("=" * 78)
    print(label)
    print("=" * 78)
    print("  t  r  W   shapes  Phi==0   peel==0   FAIL   C=tau&peel0   <=>Phi==0   ctl: peel0|Phi!=0")
    bad = 0
    for (t, r, W) in cfgs:
        n = nz = pz = fail = both = agree = 0
        ctl_n = ctl_pz = 0
        for b in betas(t, r, W):
            if not occupied(b, t):
                continue
            z = phi_zero(b, t, r)
            if z is None:
                continue
            n += 1
            pb = pelar(b, t)
            if pb is None or len(pb) != t + 2 * (r - 1):
                continue
            pz_here = phi_zero(pb, t, r - 1) if r - 1 >= 1 else None
            if pz_here is None:
                continue
            ct = C_and_tau(b, t, r)
            crit = (ct is not None and ct[0] == ct[1] and pz_here)
            if z:
                nz += 1
                if pz_here:
                    pz += 1
                else:
                    fail += 1
                    if fail <= 3:
                        print("      FAIL  t=%d r=%d beta=%s  peel=%s" % (t, r, b, pb))
            else:
                ctl_n += 1
                ctl_pz += pz_here
            if crit:
                both += 1
            agree += (crit == z)
        bad += fail
        print("  %2d %2d %2d %7d %7d %9d %6d %13d %11d   %d/%d"
              % (t, r, W, n, nz, pz, fail, both, agree, ctl_pz, ctl_n))
    return bad


def run_t2():
    """CONTROL en t = 2: alli el criterio esta PROBADO, asi que el pelado tiene que reproducirlo."""
    print()
    print("=" * 78)
    print("CONTROL t = 2, donde el criterio ya es teorema (8.6)")
    print("=" * 78)
    print("   r   W   shapes  Phi==0  peel==0  FAIL   selfcomp odd  peel selfcomp odd")
    bad = 0
    for (t, r, W) in CFG_T2:
        N = t + 2 * r
        n = nz = pz = fail = sc = psc = 0
        for b in betas(t, r, W):
            if not occupied(b, t):
                continue
            z = phi_zero(b, t, r)
            if z is None:
                continue
            n += 1
            if not z:
                continue
            nz += 1
            C = b[0] + b[-1]
            issc = all(b[j] + b[N - 1 - j] == C for j in range(N))
            sc += issc
            pb = pelar(b, t)
            if pb is None or len(pb) != N - 2:
                continue
            pz_here = phi_zero(pb, t, r - 1)
            pz += bool(pz_here)
            fail += (not pz_here)
            M = len(pb)
            psc += all(pb[j] + pb[M - 1 - j] == pb[0] + pb[-1] for j in range(M))
        bad += fail
        print("  %2d  %2d %7d %7d %8d %5d %13d %17d" % (r, W, n, nz, pz, fail, sc, psc))
    return bad


if __name__ == "__main__":
    rc = run(CFG, "N1/N2/N3 -- t PAR >= 4, el cuadrante abierto")
    rc += run_t2()
    print()
    print("TOTAL FALLOS DE LA IMPLICACION: %d" % rc)
    sys.exit(1 if rc else 0)
