# -*- coding: utf-8 -*-
# LA TABLA DE LOS 16:  que pared mata al candidato, y quien sobrevive.   15 de agosto de 2026.
#
# DE DONDE SALE.  Vuelta 12, recomendacion final: no barrer mas beta, sino hacer UNA tabla sobre
# t=6,r=2 con las columnas
#
#     ( mu_candidate, eta_candidate, wall-type de eta, mu_max, eta_first surviving )
#
# y su prediccion bifurcada, que es lo que la hace util:
#
#   - si los 16 caen sobre LA MISMA especie de pared  =>  el (1,1) pertenece al filtro de torsion;
#   - si caen sobre paredes distintas y aun asi dan (1,1)  =>  el (1,1) NO es del filtro: es de la
#     regla de branching que transporta el defecto desde Sp_4 al factor libre Sp_4.
#
# Las tres especies para t=6 (m=2), con  a = (eta_1+2, eta_2+1):
#     W0   a_i = 0 mod 6           (pared a=0)
#     W3   a_i = 3 mod 6           (pared a=t/2)
#     WC   a_1 = +- a_2 mod 6      (colision entre columnas)
#
# LAS 16 SON LAS DE LA CARTA.  sixteen_betas.py las reconstruye y su C1 verifica que son las mismas
# que se mandaron en la vuelta 11 (16 fallos, Delta=(1,1) en los 16, poblacion critica n=299).
#
# CONTROLES
#   C0  FATAL, Y ES NUEVO.  mu_max calculado por el branching (esta maquinaria, Sage) tiene que
#       coincidir con v_1 -- el primer nivel no nulo de la expansion de Laplace (la otra maquinaria,
#       Python).  Son dos caminos completamente distintos al mismo peso.  Y Delta = v_max - mu_max
#       tiene que salir (1,1) otra vez, ahora por el camino simplectico.
#   C1  FATAL.  A_{mu_cand} = sum_eta B tau  tiene que dar EXACTAMENTE 0 en los 16: es la definicion
#       de que el candidato muere.  Si diera != 0, la poblacion no es la que creo.
#   C2  |A_{mu_max}| = 1 en los 16  (la (H) de su vuelta 12), medida aqui de nuevo.
#   C3  SEÑUELO.  Se cuenta tambien cuantos eta del bloque de mu_cand mueren por pared y cuantos
#       sobreviven-pero-se-cancelan.  Si TODOS murieran por pared, el (1,1) seria del filtro; si el
#       bloque tiene supervivientes que se cancelan entre si, NO lo es.  Las dos ramas se imprimen:
#       la tabla no puede "confirmar" por construccion.
#   C4  no vacuidad: n impreso siempre, y el reparto por especie de pared con sus 16 filas visibles.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage wall_table.sage

import itertools, json, sys, csv
from collections import defaultdict

t, r = 6, 2
m = (t - 2) // 2
R = r + m

# ------------------------------------------------------------------ Phi por el bialternante -----
def phi_bialternante(beta, tt, nvar):
    N = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    L = LaurentPolynomialRing(K, nvar, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(N - 1, -1, -1))
    def alt(expo):
        return matrix(L, N, N, lambda i, j: x[i] ** expo[j]).determinant()
    den = alt(delta)
    if den == 0:
        return None
    q = alt(list(beta)) / den
    try:
        q = L(q)
    except Exception:
        return "NO-POLINOMIO"
    out = {}
    for e, c in zip(q.exponents(), q.coefficients()):
        e = tuple(e) if hasattr(e, '__iter__') else (e,)
        if c != 0:
            out[e] = c
    return out

_SP = {}
def sp_char(mu, rr):
    key = (tuple(mu), rr)
    if key not in _SP:
        W = WeylCharacterRing("C%d" % rr)
        el = W(W.space().from_vector(vector(list(mu))))
        d = {}
        for wt, mult in el.weight_multiplicities().items():
            k = tuple(int(v) for v in wt.to_vector())
            d[k] = d.get(k, 0) + mult
        _SP[key] = d
    return _SP[key]

def sp_producto(eta, mu):
    a, b = sp_char(eta, m), sp_char(mu, r)
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            k = e1 + e2
            out[k] = out.get(k, 0) + c1 * c2
    return out

def pelar_branching(P, tope=8000):
    P = {e: c for e, c in P.items() if c != 0}
    out = {}
    for _ in range(tope):
        if not P:
            return out, {}
        dom = [e for e in P
               if list(e[:m]) == sorted(e[:m], reverse=True) and min(e[:m]) >= 0
               and list(e[m:]) == sorted(e[m:], reverse=True) and min(e[m:]) >= 0]
        if not dom:
            return out, P
        top = max(dom, key=lambda e: (sum(e), e))
        B = P[top]
        out[(tuple(top[:m]), tuple(top[m:]))] = out.get((tuple(top[:m]), tuple(top[m:])), 0) + B
        for k, v in sp_producto(tuple(top[:m]), tuple(top[m:])).items():
            nv = P.get(k, 0) - B * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return out, P

# ------------------------------------------------------------------ el filtro y sus paredes -----
def pared(eta):
    """(tau, especie).  Especies: '-' vive, 'W0', 'W3', 'WC'."""
    a = [eta[j] + (m - (j + 1) + 1) for j in range(m)]
    esp = []
    for v in a:
        c = v % t
        if c == 0:
            esp.append("W0")
        elif 2 * c == t:
            esp.append("W3")
    if esp:
        return 0, "+".join(sorted(set(esp)))
    cl = []
    sg = 1
    for v in a:
        c = v % t
        if c <= m:
            cl.append(c)
        else:
            cl.append(t - c); sg *= -1
    if len(set(cl)) != m:
        return 0, "WC"
    perm = [m - cl[j] for j in range(m)]
    inv = sum(1 for i in range(m) for j in range(i + 1, m) if perm[i] > perm[j])
    return sg * (-1) ** inv, "-"

# ------------------------------------------------------------------ EL DICCIONARIO ------------
# ERROR DE CONCEPTO CORREGIDO (15-ago, tras la primera corrida).  v(T) vive en coordenadas de
# NEWTON (exponentes del Laurent) y mu vive en coordenadas de PESO de Sp_{2r}.  Son dos reticulos
# distintos y compararlos directo da C0 fallando en las 16 y los bloques vacios.  El diccionario es
# nuestra propia ley del mu_max:
#
#     mu  =  exponente  -  sigma_r ,     sigma_r = (N-1, N-3, ..., N-2r+1)
#
# y sigma_r SE SIMPLIFICA, que es lo que la hace legible:
#
#     N - 2k + 1 = t + (2r - 2k + 1)   =>   sigma_r = (t-1)*(1,...,1) + 2*rho_{C_r}
#
# con 2*rho_{C_r} = (2r, 2r-2, ..., 2) la suma de las raices positivas de C_r.  O sea: el
# desplazamiento de la ley NO es una lista arbitraria -- es "suma de raices positivas del grupo
# libre" mas "un corrimiento uniforme t-1 que aporta el bloque de torsion".  Se verifica en C_dic.
def sigma_r():
    return tuple(N_total - 2 * k + 1 for k in range(1, r + 1))

def sigma_r_simplificada():
    dosrho = tuple(2 * (r - k + 1) for k in range(1, r + 1))
    return tuple((t - 1) + d for d in dosrho)

def a_peso(expo):
    s = sigma_r()
    return tuple(expo[i] - s[i] for i in range(r))


def dominantes_maximales(S):
    return [mu for mu in S
            if not any(nu != mu and all(sum(nu[:k + 1]) >= sum(mu[:k + 1]) for k in range(len(mu)))
                       for nu in S)]

# ================================================================== corrida =====================
D = json.load(open("sixteen_betas.json"))["t6_r2"]
N_total = len(D[0]["beta"])
print("=" * 132)
print("LA TABLA DE LOS 16   --   t=6, r=2,  branching  Sp_8 > Sp_4 x Sp_4,  filtro tau_6")
print("=" * 132)
print("")
print("  C_dic  EL DICCIONARIO Newton -> peso.   sigma_r = %s" % (sigma_r(),))
print("         forma simplificada  (t-1) + 2*rho_{C_r} = %s   ->   %s"
      % (sigma_r_simplificada(),
         "COINCIDE" if sigma_r() == sigma_r_simplificada() else "*** NO COINCIDE ***"))
print("")
print("  n = %d formas (las de la vuelta 11, reconstruidas y verificadas por sixteen_betas.py)" % len(D))
print("")
print("  #  beta                             | mu_cand | mu_max | Delta | A_max | eta+alto de mu_cand / pared | #eta en mu_cand: muertos/vivos | eta 1o superviviente de mu_max")
print("  " + "-" * 176)

esp_cuenta = defaultdict(int)
esp_alto = defaultdict(int)
malo0 = malo1 = malo2 = 0
vivos_en_cand = defaultdict(int)
FILAS = []
VOLCADO = []
for idx, rec in enumerate(D):
    b = tuple(rec["beta"])
    vmax = a_peso(tuple(rec["v_max"]))     # el candidato de Laplace, YA en coordenadas de peso
    v1 = a_peso(tuple(rec["v_1"]))         # el primer nivel no nulo, idem
    Psi = phi_bialternante(b, 2, R)
    if Psi in (None, "NO-POLINOMIO"):
        print("  %2d  %-32s | *** Psi no calculable ***" % (idx + 1, str(b)))
        continue
    Psi = {k: QQ(v) for k, v in Psi.items()}
    B, resto = pelar_branching(Psi)
    B = {k: v for k, v in B.items() if v != 0}
    A = defaultdict(lambda: 0)
    for (eta, mu), bb in B.items():
        A[mu] += bb * pared(eta)[0]
    A = {mu: a for mu, a in A.items() if a != 0}
    maxi = dominantes_maximales(list(A)) if A else []
    mu_max = maxi[0] if len(maxi) == 1 else (tuple(maxi) if maxi else None)

    # C0: mu_max (branching, Sage) contra v_1 (Laplace, Python)
    ok0 = (len(maxi) == 1 and maxi[0] == v1)
    malo0 += (not ok0)
    # C1: el candidato muere
    a_cand = sum(bb * pared(eta)[0] for (eta, mu), bb in B.items() if mu == vmax)
    malo1 += (a_cand != 0)
    # C2: |A_max| = 1
    amax = A.get(maxi[0]) if len(maxi) == 1 else None
    malo2 += (amax is None or abs(amax) != 1)

    bloque = sorted([eta for (eta, mu) in B if mu == vmax], key=lambda e: (sum(e), e), reverse=True)
    muertos = [e for e in bloque if pared(e)[0] == 0]
    vivos = [e for e in bloque if pared(e)[0] != 0]
    alto = bloque[0] if bloque else None
    esp = pared(alto)[1] if alto else "-"
    esp_alto[esp] += 1
    for e in muertos:
        esp_cuenta[pared(e)[1]] += 1
    vivos_en_cand[len(vivos)] += 1

    bl_max = sorted([eta for (eta, mu) in B if mu == maxi[0]], key=lambda e: (sum(e), e), reverse=True) if len(maxi) == 1 else []
    primero = next((e for e in bl_max if pared(e)[0] != 0), None)

    dl = tuple(vmax[i] - maxi[0][i] for i in range(r)) if len(maxi) == 1 else None
    print("  %2d  %-32s | %-7s | %-6s | %-5s | %-5s | %-10s %-15s | %6d / %-6d %s | %s"
          % (idx + 1, str(b), str(vmax), str(mu_max), str(dl), str(amax),
             str(alto), "[" + esp + "]", len(muertos), len(vivos),
             "" if ok0 else "*** C0 ***", str(primero)))
    sys.stdout.flush()
    FILAS.append((b, vmax, mu_max, dl, amax, alto, esp, len(muertos), len(vivos), primero))

    # ---- el volcado, que es lo que se presenta.  Bloque ENTERO, no el resumen. --------------
    def bloque_de(mu):
        # int() sobre cada campo: los Integer de Sage NO son serializables y la primera corrida
        # murio aqui DESPUES de imprimir la tabla entera.  [[save-the-outputs-not-just-the-scripts]]
        return [{"eta": [int(v) for v in e], "B": int(bb), "tau": int(pared(e)[0]),
                 "especie": pared(e)[1]}
                for (e, nu), bb in sorted(B.items(), key=lambda kv: (-sum(kv[0][0]), kv[0][0]))
                if nu == mu]
    VOLCADO.append({
        "beta": [int(v) for v in b],
        "v_max_laplace": [int(v) for v in vmax], "v_1_laplace": [int(v) for v in v1],
        "mu_max_branching": [int(v) for v in maxi[0]] if len(maxi) == 1 else None,
        "Delta": [int(v) for v in dl] if dl else None,
        "A_mu_max": int(amax) if amax is not None else None,
        "A_mu_cand": int(a_cand),
        "n_terminos_branching": len(B),
        "eta_alto_de_mu_cand": [int(v) for v in alto] if alto else None,
        "especie_del_alto": esp,
        "eta_primero_superviviente_de_mu_max": [int(v) for v in primero] if primero else None,
        "bloque_mu_cand": bloque_de(vmax),
        "bloque_mu_max": bloque_de(maxi[0]) if len(maxi) == 1 else [],
        "C0_mu_max_coincide_con_laplace": bool(ok0),
    })

print("  " + "-" * 176)
print("")
print("  C0  mu_max por branching == v_1 por Laplace, y Delta=(1,1) por el camino simplectico : %s"
      % ("PASA en %d/%d" % (len(D) - malo0, len(D)) if malo0 == 0 else "*** FALLA en %d ***" % malo0))
print("  C1  A_{mu_cand} = 0 en los 16 (el candidato muere)                                   : %s"
      % ("PASA" if malo1 == 0 else "*** FALLA en %d ***" % malo1))
print("  C2  |A_{mu_max}| = 1  --  la (H) de su vuelta 12                                     : %s"
      % ("PASA" if malo2 == 0 else "*** FALLA en %d ***" % malo2))
print("")
print("=" * 132)
print("LA PREGUNTA QUE LA TABLA CONTESTA")
print("=" * 132)
print("")
print("  especie de pared del eta MAS ALTO del bloque de mu_cand, sobre las %d:" % len(D))
for k, v in sorted(esp_alto.items(), key=lambda kv: -kv[1]):
    print("      %-10s  %2d de %d" % (k, v, len(D)))
print("")
print("  reparto de TODOS los eta muertos del bloque de mu_cand, por especie:")
for k, v in sorted(esp_cuenta.items(), key=lambda kv: -kv[1]):
    print("      %-10s  %4d" % (k, v))
print("")
print("  C3  eta SUPERVIVIENTES dentro del bloque de mu_cand (los que el filtro NO mata, y que por")
print("      tanto se cancelan entre si para dar A_cand = 0).  Histograma de cuantos hay por forma:")
for k in sorted(vivos_en_cand):
    print("      %2d supervivientes : %2d formas" % (k, vivos_en_cand[k]))
print("")
print("  LECTURA, y las dos ramas estaban escritas ANTES de correr:")
print("    * si el histograma de C3 esta concentrado en 0 supervivientes y una sola especie manda,")
print("      el (1,1) es del FILTRO DE TORSION;")
print("    * si hay formas con supervivientes que se cancelan, el (1,1) NO es del filtro: es de la")
print("      regla de branching, que es la segunda rama de su vuelta 12.")

# ================================================================== EL VOLCADO ==================
# [[save-the-outputs-not-just-the-scripts]], y su reciproco: la vuelta 11 mando una tabla cuyo guion
# no se guardo.  Aqui se guardan LAS DOS COSAS, y ademas los datos en formato presentable.
# int() TAMBIEN en config: la primera vez se corrigieron las filas y no esto, y el volcado murio
# tras escribir 22 bytes -- con la tabla ya impresa, o sea con la corrida pareciendo un exito.
json.dump({"config": {"t": int(t), "r": int(r), "m": int(m), "R": int(R), "n": int(len(D)),
                      "poblacion": "las 16 formas de t=6,r=2 con Delta!=0, ver sixteen_betas.py"},
           "controles": {"C0_fallos": int(malo0), "C1_fallos": int(malo1), "C2_fallos": int(malo2)},
           "especie_del_eta_alto": {str(k): int(v) for k, v in esp_alto.items()},
           "reparto_muertos_por_especie": {str(k): int(v) for k, v in esp_cuenta.items()},
           "histograma_supervivientes_en_mu_cand": {str(k): int(v) for k, v in vivos_en_cand.items()},
           "filas": VOLCADO},
          open("wall_table_DUMP.json", "w"), indent=1)

with open("wall_table.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["beta", "mu_cand", "mu_max", "Delta", "A_mu_max", "A_mu_cand",
                "eta_alto_de_mu_cand", "especie_pared", "n_eta_muertos", "n_eta_vivos",
                "eta_1o_superviviente_de_mu_max", "n_terminos_branching"])
    for f in VOLCADO:
        w.writerow([f["beta"], f["v_max_laplace"], f["mu_max_branching"], f["Delta"],
                    f["A_mu_max"], f["A_mu_cand"], f["eta_alto_de_mu_cand"], f["especie_del_alto"],
                    len([x for x in f["bloque_mu_cand"] if x["tau"] == 0]),
                    len([x for x in f["bloque_mu_cand"] if x["tau"] != 0]),
                    f["eta_primero_superviviente_de_mu_max"], f["n_terminos_branching"]])

print("")
print("  DATOS GUARDADOS, no solo la tabla impresa:")
print("      wall_table_DUMP.json   volcado completo -- cada bloque eta a eta, con B, tau y especie")
print("      wall_table.csv         las %d filas, para pegar en la carta" % len(VOLCADO))
print("      wall_table_OUT.txt     esta salida")
print("")
print("=" * 132)
print("DONE")
