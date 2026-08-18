# -*- coding: utf-8 -*-
# ============================================================================================
#  VOLCADO DE DATOS de folding_t2.py -- una fila POR FORMA, no conteos.  14 de agosto de 2026.
#
#  POR QUE.  folding_t2_OUT.txt da el veredicto y los agregados; para presentar (y para que
#  cualquiera pueda re-derivar los conteos sin volver a correr nada) hacen falta los datos crudos.
#  Este guion NO recalcula nada distinto: reutiliza EXACTAMENTE las funciones de folding_t2.py,
#  cargando su cabecera, para que no pueda divergir de lo que se certifico.
#
#  SALIDAS
#    folding_t2_TABLE.tsv    una fila por lambda: r, lambda, beta, cero, signo, criterio, expansion
#    folding_t2_MIXED.txt    TODAS las expansiones con signos mezclados, completas (la refutacion)
#    folding_t2_ZEROS.txt    TODAS las formas con s_lambda = 0, con su beta y sus clases mod 2
#    folding_t2_DATA.json    lo mismo en JSON, para volver a contar sin re-correr
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python folding_t2_dump.py   (desde gates/, DESPUES de folding_t2.py)
# ============================================================================================

import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "folding_t2.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0 =")[0]
assert "def restriccion(" in _head and "def to_sp(" in _head, "folding_t2.py cambio de forma"
_ns = {"__name__": "folding_t2_preamble", "__file__": SRC}
exec(compile(_head, SRC, "exec"), _ns)
restriccion = _ns["restriccion"]
to_sp = _ns["to_sp"]
condiciones = _ns["condiciones"]
signo_de = _ns["signo_de"]
signo_predicho = _ns["signo_predicho"]
particiones = _ns["particiones"]

CFG = [(1, 14), (2, 8), (3, 5), (4, 2)]

tsv = open("folding_t2_TABLE.tsv", "w", encoding="utf-8")
tsv.write("r\tlambda\tbeta\tcero\tsigno\tsigno_predicho\tgenerico\tcriterio_i_ii\tn_mu\texpansion\n")
mix = open("folding_t2_MIXED.txt", "w", encoding="utf-8")
zer = open("folding_t2_ZEROS.txt", "w", encoding="utf-8")
mix.write("TODAS las expansiones con signos MEZCLADOS -- la refutacion de '+- combinacion no negativa'.\n")
mix.write("Formato:  r  lambda  beta  ->  suma de a_mu * sp_mu\n\n")
zer.write("TODAS las formas con s_lambda(1,-1,z^{+-1}) = 0, con sus dos clases mod 2 del beta-conjunto.\n")
zer.write("Formato:  r  lambda  beta  clase_par  clase_impar  (i)&(ii)\n\n")

DATA = {}
for r, LM in CFG:
    N = 2 * r + 2
    filas = []
    nmix = nzero = 0
    print("  r=%d ..." % r)
    sys.stdout.flush()
    for lam in particiones(N, LM):
        P, beta = restriccion(lam, r)
        cero = not P
        crit = condiciones(lam)
        coefs = {} if cero else to_sp(P, r)
        s = 0 if cero else signo_de(coefs)
        pred, gen = signo_predicho(beta, r)
        exp = sorted(((list(m), c) for m, c in coefs.items()),
                     key=lambda x: (-sum(x[0]), [-v for v in x[0]]))
        exps = " + ".join("%d*sp%s" % (c, tuple(m)) for m, c in exp) if exp else "0"
        tsv.write("%d\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%d\t%s\n"
                  % (r, "".join("%d," % x for x in lam)[:-1],
                     "".join("%d," % x for x in beta)[:-1],
                     1 if cero else 0, "" if cero else "%+d" % s,
                     "" if pred is None else "%+d" % pred, 1 if gen else 0,
                     "-" if crit is None else ("SI" if crit else "no"), len(exp), exps))
        if s == 0 and not cero:
            nmix += 1
            mix.write("r=%d  lambda=%s  beta=%s\n    %s\n" % (r, tuple(lam), tuple(beta), exps))
        if cero:
            nzero += 1
            par = sorted((b for b in beta if b % 2 == 0), reverse=True)
            imp = sorted((b for b in beta if b % 2 == 1), reverse=True)
            zer.write("r=%d  lambda=%-30s beta=%-32s par=%-20s impar=%-20s  %s\n"
                      % (r, str(tuple(lam)), str(tuple(beta)), str(par), str(imp),
                         "-" if crit is None else ("SI" if crit else "no")))
        filas.append({"lambda": lam, "beta": beta, "cero": cero, "signo": s,
                      "pred": pred, "generico": gen, "crit": crit,
                      "exp": [[list(m), c] for m, c in exp]})
    DATA["r%d" % r] = filas
    print("     %d formas, %d ceros, %d mezcladas" % (len(filas), nzero, nmix))
    sys.stdout.flush()

tsv.close()
mix.close()
zer.close()
json.dump(DATA, open("folding_t2_DATA.json", "w"), indent=0)

print("")
print("  escritos: folding_t2_TABLE.tsv, folding_t2_MIXED.txt, folding_t2_ZEROS.txt, folding_t2_DATA.json")
for f in ("folding_t2_TABLE.tsv", "folding_t2_MIXED.txt", "folding_t2_ZEROS.txt", "folding_t2_DATA.json"):
    print("     %-28s %d bytes" % (f, os.path.getsize(f)))
print("DONE")
