# -*- coding: utf-8 -*-
# ¿SON LOS DOS FILTROS LA MISMA FUNCION, TRASLADADA?   16 de agosto de 2026.
#
# DE DONDE SALE.  Cruzando las dos formulas de exponentes que el paper usa por separado:
#
#     tipo C, rango m :   x^C_i = xi^{ eta_i + m - i + 1 }          = xi^{c_i},  c = eta + rho_C
#     tipo B, rango m':   x^B_i = xi^{ 2 eta_i + 2(m'-i) + 1 }      = xi^{2u_i+1},  u_i = eta_i + m'-i
#
# El exponente de B es SIEMPRE impar porque rho_B es semientero.  Pero si 2 es invertible modulo t
# --que pasa exactamente cuando t es impar-- eso no es un obstaculo sino un CAMBIO DE COORDENADA.
# Con s = 2^{-1} = (t+1)/2:
#
#     2u+1 == 0        <=>  u == -s          <=>  c == 0            con c := u + s
#     2u_i+1 == 2u_j+1 <=>  u_i == u_j       <=>  c_i == c_j
#     2u_i+1 == -(2u_j+1) <=> u_i+u_j == -1  <=>  c_i + c_j == 0
#     2u+1 == -(2u+1)  <=>  2(2u+1) == 0     <=>  c == 0            (t impar: 2 invertible)
#
# Las CUATRO paredes de B van a las CUATRO paredes de C.  O sea, la conjetura:
#
#   (T1)  el lugar regular de tipo B a orden impar t es el de tipo C TRASLADADO por s en cada
#         coordenada;  en pesos, eta |--> eta + ((t-1)/2, ..., (t-1)/2).
#   (T2)  y si ademas los SIGNOS se corresponden con un signo global, los dos filtros son
#         literalmente la misma funcion sobre (Z/t)^m, y el filtro impar deja de necesitar
#         demostracion propia: es el par, movido.
#
# Esto NO es un ajuste: la traslacion es 2^{-1}, y existe exactamente cuando 2 es invertible.
# Si (T1) vale, la "segunda manifestacion de la paridad" del paper deja de ser una observacion
# sobre que paredes hay y pasa a ser un enunciado sobre la ACCION que las identifica.
#
# LO QUE SE MIDE, sobre TODOS los vectores de exponentes de (Z/t)^m, no sobre una muestra
#   T1  |{c regular en C}| contra |{u regular en B}|, y si el trasladado del segundo ES el primero.
#   T2  el signo: para cada u regular, signo_B(u) contra signo_C(u+s).  ¿constante el cociente?
#
# CONTROLES -- los tres pueden fallar, y dos TIENEN que fallar
#   C0  t PAR: no hay 2^{-1}, asi que (T1) tiene que FALLAR.  Si sale bien, el control no mide.
#   C1  DECOY: la misma traslacion con s' != 2^{-1} (todos los demas s' del rango) tiene que
#       fallar.  Si funcionara con cualquier s', la coincidencia es trivial.
#   C2  el conteo bruto |regular_B| == |regular_C| es condicion NECESARIA y se imprime aparte:
#       si los cardinales no cuadran ya no hay biyeccion posible, traslacion o no.
#   C3  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python filter_translate.py
import itertools
import json


def clases(v, t):
    """el multiconjunto {+-v_i} mod t, como lista de 2m elementos."""
    out = []
    for x in v:
        out.append(x % t)
        out.append((-x) % t)
    return out


def regular(v, t):
    """regular = ninguna clase nula y las 2m clases +-v_i distintas dos a dos."""
    cl = clases(v, t)
    if 0 in cl:
        return False
    return len(set(cl)) == len(cl)


def signo(v, t):
    """el signo del elemento del grupo de Weyl (permutaciones con signo) que lleva v al dominante.
       Cada coordenada se pliega al representante en (0, t/2) con su signo; luego se ordena."""
    val, sg = [], 1
    for x in v:
        e = x % t
        if e > t - e:
            e, s = t - e, -1
        else:
            s = 1
        val.append(e)
        sg *= s
    # signo de la permutacion que ordena val de forma decreciente (todos distintos si es regular)
    idx = sorted(range(len(val)), key=lambda i: -val[i])
    perm = list(idx)
    vis = [False] * len(perm)
    for i in range(len(perm)):
        if vis[i]:
            continue
        j, ciclo = i, 0
        while not vis[j]:
            vis[j] = True
            j = perm[j]
            ciclo += 1
        if ciclo % 2 == 0:
            sg = -sg
    return sg


print("=" * 118)
print("¿ES EL FILTRO IMPAR EL FILTRO PAR TRASLADADO POR 2^{-1}?")
print("=" * 118)

RES = []
for t in [3, 5, 7, 9, 11, 4, 6, 8]:
    par = (t % 2 == 0)
    for m in [1, 2, 3]:
        if t ** m > 200000:
            continue
        TODOS = list(itertools.product(range(t), repeat=m))
        # lugar regular en coordenadas de tipo C: el propio vector de exponentes
        regC = set(v for v in TODOS if regular(v, t))
        # lugar regular en coordenadas de tipo B: el vector u, con exponentes 2u+1
        regB = set(u for u in TODOS if regular([2 * x + 1 for x in u], t))
        s = None
        if not par:
            s = (t + 1) // 2                      # 2^{-1} mod t
        # T1 con el s bueno
        ok_s = None
        if s is not None:
            trans = set(tuple((x + s) % t for x in u) for u in regB)
            ok_s = (trans == regC)
        # C1 decoy: TODAS las demas traslaciones
        decoys_ok = []
        for sp in range(t):
            if s is not None and sp == s:
                continue
            trans = set(tuple((x + sp) % t for x in u) for u in regB)
            if trans == regC:
                decoys_ok.append(sp)
        # T2 signos, solo si T1 vale
        cocientes = set()
        if ok_s:
            for u in regB:
                c = tuple((x + s) % t for x in u)
                sb = signo([2 * x + 1 for x in u], t)
                sc = signo(c, t)
                cocientes.add(sb * sc)
        print("")
        print("  t=%2d %s  m=%d   |(Z/t)^m| = %d" % (t, "PAR " if par else "impar", m, len(TODOS)))
        print("     C2  |regular_B| = %-6d  |regular_C| = %-6d  %s"
              % (len(regB), len(regC), "cardinales iguales" if len(regB) == len(regC) else "!! DISTINTOS"))
        if s is not None:
            print("     T1  traslacion por s = 2^{-1} = %d : %s" % (s, "COINCIDE" if ok_s else "no coincide"))
        else:
            print("     C0  t par: no existe 2^{-1}, no hay traslacion candidata")
        print("     C1  decoy: otras traslaciones que tambien coinciden : %s"
              % (decoys_ok if decoys_ok else "ninguna"))
        if ok_s:
            print("     T2  signo_B(u) * signo_C(u+s) toma los valores %s   %s"
                  % (sorted(cocientes),
                     "-> los dos filtros son LA MISMA funcion, salvo un signo global"
                     if len(cocientes) == 1 else "-> los signos NO se corresponden"))
        RES.append({"t": int(t), "m": int(m), "par": bool(par),
                    "n_regB": len(regB), "n_regC": len(regC),
                    "s": (int(s) if s is not None else None),
                    "T1": (bool(ok_s) if ok_s is not None else None),
                    "decoys_que_coinciden": [int(x) for x in decoys_ok],
                    "cocientes_de_signo": sorted(int(x) for x in cocientes)})

print("")
print("=" * 118)
print("  LECTURA, escrita ANTES de correr:")
print("   * si T1 vale para todo t impar y NINGUN decoy coincide, y ademas falla para t par,")
print("     entonces el lugar del filtro impar es el del par trasladado por 2^{-1} -- y la razon")
print("     de que exista la traslacion es exactamente que 2 sea invertible.")
print("   * si ademas T2 da un solo valor, los dos filtros son la misma funcion sobre (Z/t)^m y")
print("     el filtro impar no necesita demostracion propia.")
print("   * si algun decoy coincide, la traslacion no distingue nada y esto no dice nada.")
json.dump(RES, open("filter_translate_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
