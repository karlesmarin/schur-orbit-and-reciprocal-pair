# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# EL ORDEN Z[g_1..g_m] SIN LA EXPLOSION DE MONOMIOS DE SAGE.
#
# POR QUE EXISTE.  `L.order(gens)` llama a absolute_order_from_ring_generators, que construye
# monomials(gens, n) con n_i = grado del polinomio minimo de g_i (order.py:2762 en sage 10.9): son
# PROD n_i elementos del cuerpo como generadores de modulo.  Con m simetricas de grado d eso es d^m.
# (7,19): 9^7 = 4.8 millones -> OOM a 5 GB.  (6,23): 11^6 = 1.8 M termino; (6,29): 14^6 = 7.5 M murio.
# Los 137 de los barridos cel6/cel7/cel8 del 16-sep eran ESTO, no casillas aritmeticamente caras.
#
# LO QUE HACE.  Cierre de reticulos:  M_0 = Z.1 ,  M_{k+1} = HNF(M_k + sum_i g_i.M_k) , hasta punto fijo.
#   - el punto fijo contiene 1 y es estable por cada g_i  =>  contiene Z[g].
#   - todo elemento de cada M_k sale de 1 multiplicando por g's y sumando  =>  esta dentro de Z[g].
#   Luego el punto fijo ES Z[g].  Termina: cadena creciente de submodulos de Z[a] (g_i enteros).
#   Coste por paso: (1+m).d productos y una HNF de ((1+m).d) x d.
#
# COMO SE USA.  load("orden_generado.sage") y luego  A = og_orden(L, gens).  Los nombres llevan og_
# porque load() comparte el espacio de nombres con el guion que lo carga (el 16-sep un "prod = ..." pisado tumbo
# 15 de 16 casillas de un guion que lo cargaba).
#
# CONTROL.  orden_generado_control.sage compara og_orden contra L.order en casillas baratas.
#
# Authors: Carles Marin, Claude (AI assistant).
from sage.rings.number_field.order import absolute_order_from_module_generators as og_desde_modulo
def og_hnf(L, vs):
    M = matrix(QQ, [list(L(v)) for v in vs])
    N = lcm([QQ(t).denominator() for t in M.list()])
    H = (matrix(ZZ, N*M)).hermite_form(include_zero_rows=False)
    return H / N
def og_orden(L, gens, max_pasos=None):
    gens = [L(g) for g in gens]
    d = L.degree()
    if max_pasos is None:
        max_pasos = 10 * d + 200
    H = og_hnf(L, [L(1)])
    pasos = 0
    while True:
        base = [L(list(r)) for r in H.rows()]
        vs = base + [g * b for g in gens for b in base]
        H2 = og_hnf(L, vs)
        pasos += 1
        if H2 == H:
            break
        if pasos > max_pasos:
            raise RuntimeError("og_orden: sin punto fijo tras %d pasos" % pasos)
        H = H2
    if H.nrows() != d:
        raise ValueError("og_orden: los generadores no dan rango completo (%d de %d)" % (H.nrows(), d))
    return og_desde_modulo([L(list(r)) for r in H.rows()], check_integral=False, check_rank=False, check_is_ring=False)
