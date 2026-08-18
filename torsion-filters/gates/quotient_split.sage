# -*- coding: utf-8 -*-
# DENTRO DE UN SOLO CORE:  6 se anulan y 6 no.  ¿QUE LAS SEPARA?   15 de agosto de 2026.
#
# DE DONDE SALE.  core_conditional.sage:  a t=6, r=2, W=13 hay 36 t-cores y solo TRES tienen alguna
# anulacion.  El mejor es  core = (1,1,1):  12 formas ocupadas,  6 nulas y 6 no.
#
# POR QUE ESTE ES EL EXPERIMENTO.  core y quotient determinan lambda entre los dos.  Fijado el core,
# lo unico que queda variar es el QUOTIENT.  Y 6 contra 6, mismo core, misma paridad, es el tamaño en
# el que la diferencia esta OBLIGADA a verse -- o a no existir, que tambien es un resultado.
# Es [[stratify-and-stop-dont-expand]]: no ampliar la poblacion, estrechar el estrato.
#
# QUE SE IMPRIME.  Los 12 casos enteros, con su quotient componente a componente, para MIRARLOS.  No
# se busca una regla automatica sobre 12 objetos -- eso es ajustar ruido.  Se imprime el dato crudo y
# las diferencias agregadas, y cualquier regla que salga se prueba DESPUES sobre los otros dos cores
# y sobre las 479 no nulas.
#
# CONTROLES
#   C0  los tres cores, no solo el mejor.  Un patron que solo vale en (1,1,1) no es un patron.
#   C1  cualquier separador que aparezca se evalua sobre la poblacion ocupada ENTERA (491), no solo
#       sobre el estrato.  Se imprime su tasa de falsos positivos.
#   C2  no vacuidad: n impreso por estrato, y se dice explicitamente cuando un estrato es demasiado
#       pequeño para concluir.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage quotient_split.sage

import itertools, json, sys
from collections import defaultdict

load("pob_helper.py")

t, r, W = 6, 2, 13
N = t + 2 * r
d = list(range(N - 1, -1, -1))


# EL BETA-SET CANONICO ES beta, NO UNO RECALCULADO DESDE lambda.
# Error cazado por el propio dato: con el beta-set recomputado desde lambda con n = numero de partes
# no nulas, tres lambdas distintas -- (3,2,2,2), (1^9), (2,2,2,2,1) -- salian con el MISMO core y el
# MISMO quotient.  Imposible: core+quotient es biyectivo.  La causa es que el quotient solo es
# canonico con el beta-set de tamaño FIJO; con n variable las clases cambian de tamaño y se pierde la
# carga.  Y el beta-set de tamaño fijo ya lo teniamos: es beta, de tamaño N = t+2r, que es ademas
# exactamente el beta-set de la maquinaria de Ayyer-Kumari.  El CORE si es independiente de n
# (verificado a mano en los dos casos), el quotient no.
def core_quot(beta, tt):
    clases = defaultdict(list)
    for b in beta:
        clases[b % tt].append(b)
    quot = {}
    for q in range(tt):
        vs = sorted(clases[q], reverse=True)
        k = len(vs)
        quot[q] = tuple(x for x in ((vs[i] - q) // tt - (k - 1 - i) for i in range(k)) if x)
    nuevo = []
    for q in range(tt):
        nuevo += [q + tt * i for i in range(len(clases[q]))]
    nuevo.sort(reverse=True)
    core = tuple(x for x in (nuevo[i] - (len(nuevo) - 1 - i) for i in range(len(nuevo))) if x)
    return core, tuple(quot[q] for q in range(tt))


DATOS = []
for b in betas_py(t, r, W):
    if not occupied_py(b, t):
        continue
    z = phi_zero_py(b, t, r)
    if z is None:
        continue
    lam = tuple(x for x in (b[i] - d[i] for i in range(N)) if x != 0)
    c, q = core_quot(tuple(b), t)          # el beta-set canonico, tamaño N fijo
    tam = tuple(sum(1 for x in b if x % t == qq) for qq in range(t))
    DATOS.append({"beta": tuple(b), "lam": lam, "core": c, "quot": q, "tam": tam, "nula": bool(z)})

print("=" * 126)
print("DENTRO DE UN SOLO CORE  --  t=%d  r=%d  W=%d,  %d formas ocupadas, %d nulas"
      % (t, r, W, len(DATOS), sum(1 for x in DATOS if x["nula"])))
print("=" * 126)

CORES = [(1, 1, 1), (2, 1), (2, 2, 2, 1)]
for C in CORES:
    est = [x for x in DATOS if x["core"] == C]
    nl = [x for x in est if x["nula"]]
    no = [x for x in est if not x["nula"]]
    print("")
    print("-" * 126)
    print("CORE = %s     n = %d   (%d nulas / %d no nulas)%s"
          % (str(C), len(est), len(nl), len(no),
             "   *** estrato demasiado pequeño para concluir ***" if len(est) < 6 else ""))
    print("-" * 126)
    print("  nula | beta                             | lambda                    | tam clases | quotient (componente a componente)")
    for x in sorted(est, key=lambda y: (not y["nula"], y["beta"])):
        print("   %-3s | %-32s | %-25s | %-10s | %s"
              % ("SI" if x["nula"] else "no", str(x["beta"]), str(x["lam"]), str(x["tam"]),
                 "  ".join(str(c) if c else "-" for c in x["quot"])))
    if not nl or not no:
        continue
    # ---- diferencias agregadas, para MIRAR, no para ajustar --------------------------------
    def resumen(L):
        return {
            "|quot| total": sorted(sum(sum(c) for c in x["quot"]) for x in L),
            "componentes no vacias": sorted(sum(1 for c in x["quot"] if c) for x in L),
            "quot mas largo": sorted(max((len(c) for c in x["quot"]), default=0) for x in L),
            "posiciones ocupadas": sorted(tuple(i for i, c in enumerate(x["quot"]) if c) for x in L),
        }
    A, B = resumen(nl), resumen(no)
    print("")
    print("  %-24s | nulas%-30s | no nulas" % ("agregado", ""))
    for k in A:
        print("  %-24s | %-35s | %s" % (k, str(A[k])[:35], str(B[k])[:60]))
    sys.stdout.flush()

# ------------------------------------------------------------------ C1 --------------------------
print("")
print("=" * 126)
print("C1  CUALQUIER separador tiene que probarse sobre las %d ocupadas, no sobre el estrato." % len(DATOS))
print("=" * 126)
print("")
CAND = [
    ("suma del quotient impar",   lambda x: sum(sum(c) for c in x["quot"]) % 2 == 1),
    ("una sola componente no vacia", lambda x: sum(1 for c in x["quot"] if c) == 1),
    ("dos componentes no vacias",  lambda x: sum(1 for c in x["quot"] if c) == 2),
    ("componente 0 no vacia",      lambda x: bool(x["quot"][0])),
    ("quotient con parte >= 2",    lambda x: any(any(v >= 2 for v in c) for c in x["quot"])),
]
nulas_tot = [x for x in DATOS if x["nula"]]
print("  %-32s | cubre nulas | falsos + sobre las %d no nulas" % ("regla sobre el quotient", len(DATOS) - len(nulas_tot)))
print("  " + "-" * 100)
for nom, f in CAND:
    cn = sum(1 for x in nulas_tot if f(x))
    fp = sum(1 for x in DATOS if not x["nula"] and f(x))
    print("  %-32s | %5d/%-5d | %d %s"
          % (nom, cn, len(nulas_tot), fp,
             "*** CANDIDATA ***" if (cn == len(nulas_tot) and fp == 0) else ""))
print("")
json.dump([{"beta":[int(v) for v in x["beta"]],"lam":[int(v) for v in x["lam"]],"core":[int(v) for v in x["core"]],"quot":[[int(v) for v in c] for c in x["quot"]],"tam":[int(v) for v in x["tam"]],"nula":bool(x["nula"])} for x in DATOS],
          open("quotient_split_DUMP.json", "w"), indent=1)
print("  volcado en quotient_split_DUMP.json")
print("=" * 126)
print("DONE")
