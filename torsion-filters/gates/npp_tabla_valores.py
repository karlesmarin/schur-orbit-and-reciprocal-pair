# -*- coding: utf-8 -*-
# LA TABLA COMPLETA DE sp_lambda(x_m), AUTOCONTENIDA.   20 de agosto de 2026.
#
# POR QUE.  Una tabla vale mas que un resumen: el resumen hay que creerselo y la tabla se
# comprueba.  Y se fabrica ANTES de anunciarla, no despues --- anunciar un adjunto que todavia
# no existe es como se acaba recomendando lo que no se ha escrito.
#
# QUE PRODUCE.  Un unico fichero de texto, autocontenido y legible sin nuestro codigo:
#   - los tres elementos, con su orden y la comprobacion de regularidad;
#   - el criterio de anulacion en las dos lecturas, para que el lector decida cual es la suya;
#   - la tabla completa de sp_lambda(x_m) para m = 1, 2 y 3;
#   - y como regenerarla.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python npp_tabla_valores.py    (escribe npp_table_values.txt)

import io
import itertools
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

exec(open("npp_question81_v2.py", encoding="utf-8").read().split('print("=" * 100)')[0])

LMAX = {1: 12, 2: 8, 3: 6}
OUT = []
W = OUT.append

W("=" * 92)
W("CHARACTER VALUES OF Sp(2m) AT AN EVEN-ORDER REGULAR NON-PRINCIPAL TORSION ELEMENT")
W("=" * 92)
W("")
W("Carles Marin, 20 August 2026.  Companion table to the letter; every number below is")
W("computed in exact cyclotomic arithmetic (no floating point) by the script named at the end.")
W("")
W("THE ELEMENTS.  For each m, put h = 2m (the Coxeter number of C_m), t = h + 2 = 2m + 2, and")
W("xi = exp(2 pi i / t).  The element is the diagonal torus element")
W("")
W("      x_m  =  diag( xi, xi^2, ..., xi^m, xi^-m, ..., xi^-2, xi^-1 )   in  Sp(2m).")
W("")
W("It has order t.  It is regular: alpha(x_m) != 1 for every root of C_m, checked over the whole")
W("root system.  It is not principal: on the simple roots of C_m the m-1 short ones give xi^-1 and")
W("the long one gives xi^-2, whereas a principal element takes one primitive root on all of them.")
W("Equivalently, the principal element of order t = h+2 would need half-integer torus exponents")
W("(m-1/2, ..., 1/2), so it does not lie in the torus of Sp(2m) at all.")
W("")

for m in (1, 2, 3):
    t = 2 * m + 2
    W("      m = %d   %-12s  t = %d   x_%d = diag(%s)"
      % (m, {1: "Sp(2)=SL(2)", 2: "Sp(4)", 3: "Sp(6)"}[m], t, m,
         ", ".join(["xi^%d" % k for k in range(1, m + 1)]
                   + ["xi^-%d" % k for k in range(m, 0, -1)])))
W("")
W("      For m = 1, xi = i and x_1 = diag(i, -i).")
W("")
W("-" * 92)
W("THE VANISHING CRITERION, in the two readings that your G_2 example does not separate")
W("-" * 92)
W("")
W("Write lambda + rho = (lambda_1 + m, lambda_2 + m-1, ..., lambda_m + 1), and count")
W("")
W("      N_G    =  #{ alpha in Phi(C_m) : t divides <lambda+rho, alpha> },   [roots of the group]")
W("      N_dual =  #{ alpha in Phi(B_m) : t divides <lambda+rho, alpha> },   [roots of the dual]")
W("")
W("both counted over roots of both signs.  The two differ only on the long roots of C_m, where")
W("2 e_i is paired against the coroot e_i.  Over all the weights tabulated below,")
W("")
W("      sp_lambda(x_m) = 0   <==>   N_G > 0        without exception,")
W("")
W("while the same statement with N_dual fails, at every m: over lambda_1 <= 9 it misses 3 weights")
W("for SL(2), 16 for Sp(4) and 54 for Sp(6), the smallest being lambda = (1), (1,0) and (1,0,0).")
W("(For SL(2) alone, N_dual does become correct if t is replaced by the order of x_1 in the adjoint")
W("group, which is 2; that substitution then breaks Sp(4) and Sp(6).  No single choice other than")
W("the one above works for all three.)  In G_2 at d = 2 the two counts always agree, because root")
W("and coroot there differ by the odd factor 3, which is why that example cannot separate them.")
W("")

for m in (1, 2, 3):
    t = 2 * m + 2
    phi = ciclotomico(t)
    a = list(range(1, m + 1))
    W("-" * 92)
    W("VALUES FOR m = %d   (%s, t = %d),  all lambda with lambda_1 <= %d"
      % (m, {1: "SL(2)", 2: "Sp(4)", 3: "Sp(6)"}[m], t, LMAX[m]))
    W("-" * 92)
    W("")
    W("   %-22s %-16s %6s %8s %10s" % ("lambda", "lambda+rho", "N_G", "N_dual", "sp_lambda(x)"))
    lams = [l for l in itertools.product(range(0, LMAX[m] + 1), repeat=m)
            if all(l[i] >= l[i + 1] for i in range(m - 1))]
    lams.sort(key=lambda l: (sum(l), l))
    for lam in lams:
        v = sp_valor(list(lam), m, a, t, phi)
        lr = [lam[i] + (m - i) for i in range(m)]
        W("   %-22s %-16s %6d %8d %10s"
          % ("(" + ",".join(map(str, lam)) + ")", "(" + ",".join(map(str, lr)) + ")",
             hits(lr, m, t, "B"), hits(lr, m, t, "A"), v))
    W("")

W("-" * 92)
W("HOW THIS WAS COMPUTED")
W("-" * 92)
W("")
W("sp_lambda(x_m) is evaluated from the type-C bialternant, as a quotient of two determinants in")
W("Z[xi] reduced modulo the t-th cyclotomic polynomial, so every entry is an exact integer and no")
W("root of unity is ever approximated.  The regularity check evaluates alpha(x_m) over the entire")
W("root system rather than arguing it.  Scripts, with their archived output:")
W("")
W("      npp_elementos_regulares.py    the elements: order, regularity, non-principality")
W("      npp_calibrado_G2.py           the calibration against Section 7 of arXiv:2504.14684")
W("      npp_question81_v2.py          the two readings against the computed values")
W("      npp_tabla_valores.py      this table")
W("")
W("They are ancillary files of arXiv:2608.18302 and can be sent directly on request.")
W("")
W("=" * 92)

io.open("npp_table_values.txt", "w", encoding="utf-8", newline="\n").write("\n".join(OUT) + "\n")
print("escrito npp_table_values.txt  --  %d lineas" % len(OUT))
