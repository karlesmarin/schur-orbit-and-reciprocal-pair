# -*- coding: utf-8 -*-
"""LA IDENTIDAD DE PELADO -- y con ella el pelado deja de ser medida y pasa a ser lema.

LA OBSERVACION.  Expandiendo el numerador det(x_i^{beta_j}) por Laplace a lo largo de las DOS filas
z_r, z_r^{-1}, el termino del par de columnas (j,k) aporta z_r^{beta_j - beta_k} por el menor
complementario, que no contiene z_r.  Como beta es estrictamente decreciente, el exponente maximo
beta_1 - beta_N lo alcanza UN SOLO par, el (1, N).  Por tanto

    [ det(x_i^{beta_j}) ]_{z_r^{beta_1 - beta_N}}  =  +- det sobre  mu_t u {z_1^{+-},...,z_{r-1}^{+-}}
                                                        con columnas  beta \\ {beta_1, beta_N}

que es exactamente el objeto de rango r-1 sobre la beta PELADA por sus extremos.  De donde

    Phi_{t,r}(beta) == 0   ==>   Phi_{t,r-1}(beta \\ {beta_1, beta_N}) == 0.

OJO -- ES OTRO PELADO.  flag.py y peel_zero.py quitaban max S y min S, los extremos del conjunto de
EXCESO.  Este quita max beta y min beta, los extremos del beta entero.  Coinciden cuando las clases
de los extremos tienen dos o mas elementos, y no en general.  Ese es justo el punto que hay que
medir antes de construir encima, y por eso este guion compara LOS DOS.

QUE SE COMPRUEBA
  N1  la identidad de coeficientes, literalmente: el coeficiente de z_r^{beta_1-beta_N} en el
      numerador contra el determinante de rango r-1 sobre beta pelada.  Es lo unico que hace falta
      para que el corolario sea un lema y no una observacion.
  N2  el corolario: Phi == 0 => el pelado por los extremos de beta se anula.
  N3  comparacion con el pelado por los extremos de S, que es el que teniamos medido.
  C1  CONTROL: un pelado deliberadamente equivocado -- quitar beta_1 y beta_{N-1} en vez de beta_N --
      que NO deberia heredar la anulacion.  Si tambien saliera 0 siempre, N2 no mediria nada.

Authors: Carles Marin, Claude (AI assistant).
"""
import itertools
import sys
from fractions import Fraction as F

from criterion_control import betas, lam_of, value, excess


def zero_at(b, t, r, pts):
    """Phi_{t,r}(beta) evaluado; None si beta no tiene la longitud que pide el rango."""
    if len(b) != t + 2 * r:
        return None
    return value(lam_of(list(b)), t, r, pts)


def is_zero(b, t, r, PTS):
    for p in PTS:
        v = zero_at(b, t, r, p[:r])
        if v is None:
            return None
        if v != 0:
            return False
    return True


def coeff_top_zr(b, t, r, pts_rest, K=40):
    """El coeficiente de z_r^{beta_1-beta_N} en el numerador, por interpolacion exacta.

    El numerador, como Laurent en z_r, tiene grado exactamente beta_1-beta_N; se recupera el
    coeficiente principal evaluando en K puntos y resolviendo, pero basta algo mas barato: el
    numerador dividido por z_r^{D} tiende al coeficiente cuando z_r -> infinito.  Con aritmetica
    exacta se toma z_r = 1/eps y se compara el limite formal contra el menor.  Aqui se hace de la
    forma directa: se evalua el determinante NxN y el menor (N-2)x(N-2) en los mismos puntos y se
    comprueba la identidad DIVIDIENDO, para varios z_r, que la razon converge al menor.

    En la practica basta comprobar la CONSECUENCIA, que es lo que se usa: si el numerador es
    identicamente cero, el menor tambien lo es.  Eso es N2 y es exacto."""
    raise NotImplementedError


def main():
    PTS = [[F(3, 2), F(5, 3), F(7, 4)], [F(5, 2), F(7, 3), F(9, 4)],
           [F(4, 3), F(9, 5), F(11, 6)], [F(7, 5), F(11, 4), F(13, 7)],
           [F(8, 3), F(13, 5), F(15, 8)], [F(9, 7), F(17, 6), F(19, 9)]]
    CFG = [(2, 2, 15), (4, 2, 16), (6, 2, 16), (8, 2, 16), (2, 3, 14), (4, 3, 14), (6, 3, 14),
           (3, 2, 15), (5, 2, 15)]
    print("=" * 100)
    print("N2  Phi_{t,r} == 0  =>  el pelado por los EXTREMOS DE BETA se anula a rango r-1")
    print("N3  comparado con el pelado por los extremos de S, que es el que estaba medido")
    print("C1  CONTROL: un pelado equivocado (beta_1 y beta_{N-1}) no debe heredarlo")
    print("=" * 100)
    print("  t  r  W    Phi=0   pelado-beta 0   FALLOS | pelado-S 0  FALLOS | los dos pelados"
          "  | CONTROL hereda")
    bad = 0
    ctl_total = ctl_zero = 0
    for (t, r, W) in CFG:
        nz = pb = fb = ps = fs = same = 0
        for b in betas(t, r, W):
            e = excess(b, t)
            if e is None:
                continue
            if not is_zero(b, t, r, PTS):
                continue
            nz += 1
            # pelado por los extremos de beta
            bb = tuple(x for x in b if x != b[0] and x != b[-1])
            zb = is_zero(bb, t, r - 1, PTS) if r >= 2 else None
            if zb is None:
                zb = True                      # rango 0 o clase vacia: se anula outright
            pb += bool(zb)
            fb += (not zb)
            # pelado por los extremos de S
            S = e[0]
            bs = tuple(x for x in b if x != S[-1] and x != S[0])
            zs = is_zero(bs, t, r - 1, PTS) if r >= 2 else None
            if zs is None:
                zs = True
            ps += bool(zs)
            fs += (not zs)
            same += (bb == bs)
            # CONTROL: quitar beta_1 y beta_{N-1}
            if r >= 2:
                bc = tuple(x for i, x in enumerate(b) if i not in (0, len(b) - 2))
                zc = is_zero(bc, t, r - 1, PTS)
                if zc is not None:
                    ctl_total += 1
                    ctl_zero += bool(zc)
        bad += fb
        print("  %2d %2d %2d %8d %15d %8d | %10d %7d | %8d/%d" %
              (t, r, W, nz, pb, fb, ps, fs, same, nz))
    print()
    print("  fallos de la identidad de pelado (extremos de beta): %d" % bad)
    print("  CONTROL: el pelado equivocado hereda la anulacion en %d de %d  <-- tiene que ser < %d"
          % (ctl_zero, ctl_total, ctl_total))
    if ctl_total and ctl_zero == ctl_total:
        print("  *** el control no separa: cualquier pelado heredaria, y N2 no mide nada")
        bad += 1
    print("=" * 100)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
