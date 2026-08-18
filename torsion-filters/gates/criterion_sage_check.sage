# -*- coding: utf-8 -*-
# CONTROL EXTERNO del criterio del 14 de agosto, por VIA INDEPENDIENTE, con cascada y checkpoint.
#
# POR QUE.  Todo criterion_S.py sale de un solo instrumento -- setup(), all_transversals() y
# measure(): Laplace, estratos y voraz.  Si measure() tuviera un fallo, las 24 configuraciones con
# 0 discrepancias serian 24 acuerdos consigo mismo.  Aqui NO se usa nada de eso.
#
# LA VIA -- y OJO, no es la que decia este comentario.  La primera version expandia s_lambda en
# SUMAS DE POTENCIAS; se cambio a JACOBI-TRUDI por coste (ver el docstring de value()) y el rotulo
# se quedo atras.  Lo que corre es: se toma la generatriz de las h_k sobre el alfabeto entero,
#
#     H(u) = 1/(1-u^t) * prod_j 1/((1 - z_j u)(1 - u/z_j)),
#
# donde el factor 1/(1-u^t) ES la orbita colapsada, prod_{k<t}(1 - zeta^k u) = 1 - u^t, de modo que
# no aparece ninguna raiz de la unidad; y s_lambda = det(h_{lambda_i - i + j}), una NxN.
#
# Ni bialternante, ni expansion de Laplace, ni estratificacion por grados, ni greedy.  Es la ruta de
# funciones simetricas que ya usan excess_invariant.sage y sus hermanos.
#
# LA CASCADA -- y es lo que lo hace viable.  El primer intento de este control calculaba el
# determinante simbolico entero: 18 s por forma, 3-4 h de corrida.  El 99 % de la poblacion no se
# anula, y descartarlo es barato:
#
#     paso 1  evaluar en z_j = 1 sobre Q                  -> si != 0, NO se anula.  Definitivo.
#     paso 2  evaluar en 4 puntos racionales al azar      -> si alguno != 0, NO se anula.  Definitivo.
#     paso 3  solo los supervivientes pagan el anillo de Laurent exacto.
#
# Los pasos 1 y 2 sólo pueden dar negativos, y un negativo suyo es una PRUEBA de no anulacion; el
# paso 3 es el unico que puede afirmar el cero.  Esa asimetria es la que hace correcta la cascada.
# El patron es de paper/anc/selfcomp_law.sage, que ya lo hacia y no mire antes de escribir el mio.
#
# CHECKPOINT.  Cada forma se escribe a criterion_sage_CKPT.jsonl con flush, y al arrancar se lee y
# se salta lo ya hecho.  Una interrupcion cuesta la forma en curso, no la corrida.
#
# CONTROLES
#   K1  en r = 1 el instrumento tiene que reproducir el Teorema 3.1, que esta PROBADO.
#   K2  un señuelo que TIENE que fallar: pedir C - beta = beta en vez de C - S = S.
#   K3  la cascada tiene que ser correcta: se comprueba, sobre los que el paso 3 declara CERO, que
#       los pasos 1 y 2 tambien daban cero.  Un desacuerdo ahi seria un fallo del cedazo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         sagemath/sagemath:latest sage /work/criterion_sage_check.sage

import itertools, json, os, random, sys, time

CKPT = "/work/criterion_sage_CKPT.jsonl"
CFG = [(2, 1, 12), (3, 1, 12), (4, 1, 12), (5, 1, 11), (6, 1, 11),
       (2, 2, 13), (4, 2, 14), (6, 2, 14), (8, 2, 15), (3, 2, 13), (5, 2, 13),
       (2, 3, 13), (4, 3, 13), (6, 3, 14)]

Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()
random.seed(int(20260814))
_pcache = {}


def value(lam, t, r, zs, R):
    """s_lambda evaluado en mu_t union {z_j, 1/z_j}, por JACOBI-TRUDI.

    La primera version de este guion expandia s_lambda en sumas de potencias.  Es correcto y es
    inutilizable: para beta ancha, lambda llega a |lambda| = 50 y esa expansion es una suma sobre
    las 204226 particiones de 50.  El coste de Jacobi-Trudi depende de N y del grado, no de
    p(|lambda|), que es la razon por la que selfcomp_law.sage lo hace asi.

    Y la funcion generatriz se simplifica sola, porque la orbita entera colapsa:

        prod_{k<t} (1 - zeta^k u) = 1 - u^t,

    de modo que H(u) = 1/(1-u^t) * prod_j 1/((1 - z_j u)(1 - u/z_j)).  No aparece ninguna raiz de
    la unidad: todo el calculo vive sobre Q, o sobre el anillo de Laurent en el paso 3."""
    N = t + 2 * r
    L = list(lam) + [0] * (N - len(lam))
    D = int(L[0] + N + 1)
    PK = PowerSeriesRing(R, 'u', default_prec=D + 2)
    u = PK.gen()
    G = PK(1) / (1 - u ** t)
    for z in zs:
        G *= PK(1) / ((1 - R(z) * u) * (1 - u / R(z)))
    h = [G[k] for k in range(D + 1)]
    H = lambda k: R(0) if k < 0 else R(h[k])
    return matrix(R, N, N, lambda i, j: H(L[i] - i + j)).det()


def lam_of(b):
    N = len(b)
    return tuple(b[i] - (N - 1 - i) for i in range(N))


def excess(b, t):
    cl = {}
    for v in b:
        cl.setdefault(v % t, []).append(v)
    if len(cl) < t:
        return None
    S, incr = [], []
    for k, vs in cl.items():
        if len(vs) >= 2:
            c = sorted(vs, reverse=True)
            S += c
            incr += [c[i] + c[i + 1] for i in range(len(c) - 1)]
    if not S:
        return None
    S = sorted(S)
    return S, S[0] + S[-1], incr


def concentric(b, t):
    cl = {}
    for v in b:
        cl.setdefault(v % t, []).append(v)
    big = [sorted(v, reverse=True) for v in cl.values() if len(v) >= 2]
    if len(big) == 1:
        return False
    return big[0][0] + big[0][1] == big[1][0] + big[1][1]


def betas(t, r, W):
    N = t + 2 * r
    for mid in itertools.combinations(range(1, W + 1), N - 1):
        yield tuple(sorted(mid, reverse=True)) + (0,)


# ------------------------------------------------------------------ checkpoint
done = {}
if os.path.exists(CKPT):
    for line in open(CKPT):
        line = line.strip()
        if line:
            d = json.loads(line)
            done[d["k"]] = d
    print("checkpoint: %d formas ya calculadas" % len(done), flush=True)
FH = open(CKPT, "a")


def decide(lam, t, r):
    """la cascada.  Devuelve (cero?, paso que lo decidio, cero segun el cedazo barato)."""
    Q1 = [QQ(1)] * r
    if value(lam, t, r, Q1, QQ) != 0:
        return False, 1, False
    for _ in range(4):
        pts = [QQ(int(random.randint(int(3), int(200)))) / QQ(int(random.randint(int(3), int(200)))) for _ in range(r)]
        if value(lam, t, r, pts, QQ) != 0:
            return False, 2, False
    R = LaurentPolynomialRing(QQ, ['z%d' % j for j in range(r)])
    return value(lam, t, r, list(R.gens()), R) == 0, 3, True


print("=" * 96, flush=True)
print("CONTROL EXTERNO EN SAGE -- Jacobi-Trudi + cascada.  Sin Laplace, estratos ni voraz.", flush=True)
print("=" * 96, flush=True)
print("  t  r  W   formas  ceros  criterio  FP  FN | K1 r=1 | K2 señuelo FP FN | paso3  seg", flush=True)
bad = k1bad = decoy = k3bad = 0
TOTn = TOTz = 0
for (t, r, W) in CFG:
    t0 = time.time()
    n = nz = ncrit = fp = fn = dfp = dfn = k1 = deep = 0
    for b in betas(t, r, W):
        e = excess(b, t)
        if e is None:
            continue
        S, C, incr = e
        key = "%d,%d,%s" % (t, r, ",".join(map(str, b)))
        if key in done:
            z, step = done[key]["z"], done[key]["s"]
        else:
            z, step, cheap = decide(lam_of(b), t, r)
            FH.write(json.dumps({"k": key, "z": int(z), "s": int(step)}) + "\n")
            FH.flush()
        n += 1
        nz += bool(z)
        deep += (step == 3)
        crit = sorted(C - v for v in S) == S and (C in incr)
        ncrit += bool(crit)
        fp += bool(crit and not z)
        fn += bool(z and not crit)
        d = sorted(C - v for v in b) == sorted(b)
        dfp += bool(d and not z)
        dfn += bool(z and not d)
        if r == 1:
            k1 += bool(concentric(b, t) != z)
    bad += fp + fn
    k1bad += k1
    TOTn += n
    TOTz += nz
    if dfp or dfn:
        decoy += 1
    print("  %2d %2d %2d %8d %6d %9d %3d %3d | %6s | %14d %2d | %5d %5.0f"
          % (t, r, W, n, nz, ncrit, fp, fn,
             ("%d" % k1) if r == 1 else "-", dfp, dfn, deep, time.time() - t0), flush=True)

FH.close()

# ------------------------------------------------------------------ K3, el cedazo
# ESTE CONTROL FALTABA.  La cabecera lo prometia desde el primer dia y no estaba escrito; y tal
# como estaba enunciado -- "sobre los que el paso 3 declara CERO, comprobar que 1 y 2 daban cero" --
# NO PUEDE FALLAR: al paso 3 solo llegan los que dieron cero en 1 y en 2.  El control con contenido
# es el reciproco: coger formas que el cedazo DESCARTO en el paso 1 o 2 y pagarles el anillo de
# Laurent exacto.  Un cero ahi seria un falso descarte, o sea el cedazo comiendose un cero.
print(flush=True)
print("  K3  el cedazo: recalculo exacto sobre formas DESCARTADAS en los pasos 1-2", flush=True)
print("      t  r  muestra  ceros hallados (tiene que ser 0)", flush=True)
rnd3 = random.Random(int(20260815))
for (t, r, W) in CFG:
    pool = []
    for b in betas(t, r, W):
        if excess(b, t) is None:
            continue
        key = "%d,%d,%s" % (t, r, ",".join(map(str, b)))
        d = done.get(key)
        if d is not None and d["s"] in (1, 2) and d["z"] == 0:
            pool.append(b)
    if not pool:
        continue
    smp = rnd3.sample(pool, min(int(3), len(pool)))
    hit = 0
    for b in smp:
        R = LaurentPolynomialRing(QQ, ['z%d' % j for j in range(r)])
        if value(lam_of(b), t, r, list(R.gens()), R) == 0:
            hit += 1
            print("      *** FALSO DESCARTE t=%d r=%d beta=%s" % (t, r, b), flush=True)
    k3bad += hit
    print("      %2d %2d %8d %6d" % (t, r, len(smp), hit), flush=True)
if k3bad:
    bad += k3bad

print(flush=True)
print("=" * 96, flush=True)
print("  TOTAL: %d formas, %d ceros, %d configuraciones" % (TOTn, TOTz, len(CFG)), flush=True)
print("  fallos del criterio: %d" % bad, flush=True)
print("  K1  r=1 contra el Teorema 3.1 probado: %d desacuerdos" % k1bad, flush=True)
print("  K2  el señuelo (beta entero) falla en %d de %d configuraciones" % (decoy, len(CFG)), flush=True)
print("  K3  el cedazo: %d falsos descartes" % k3bad, flush=True)
if decoy == 0:
    print("  *** el señuelo no falla: el control no separa S de beta", flush=True)
    bad += 1
print("=" * 96, flush=True)
print("VEREDICTO: %s" % ("CONFIRMA" if bad + k1bad + k3bad == 0 else "DISCREPA"), flush=True)
