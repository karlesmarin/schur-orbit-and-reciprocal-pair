# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: datos del piso 1 por (n, p) para las tablas, reutilizando SIN modificar orden_por_caracteres.py (orden: e exacto
#      90/90) e indice_certificado.py (certificado (K) 90/90). Esos scripts ejecutan sus baterias al importarse; aqui se
#      carga solo el codigo de sus funciones (hasta la linea 'PAT = re.compile'), cada uno en su propio espacio de nombres.
#      piso1(n, p) -> dict: d, e (conductor; exacto si cotas iguales), e_inf, e_sup, v_p del indice (delta_A =
#      p d - 2 + Z(pn) + Z(n)), Z(pn), Z(n), certificado K (Delta = delta_seed - delta_A; 0 = graduado entero generado).
# FECHA: 1-oct-2026
import types


def _carga(ruta):
    src = open(ruta, encoding='utf-8').read()
    corte = src.index('PAT = re.compile')
    ns = {'__name__': 'carga_' + ruta.split('/')[-1][:-3]}
    exec(compile(src[:corte], ruta, 'exec'), ns)  # espacio de nombres propio (no compartimos globales)
    return ns


_OC = _carga('/work/orden_por_caracteres.py')
_OCmod = types.SimpleNamespace(B_prim=_OC['B_prim'], orden=_OC['orden'])
_src_ic = open('/work/indice_certificado.py', encoding='utf-8').read()
_src_ic = _src_ic.replace('import orden_por_caracteres as OC  # noqa: F401  (su bloque final tambien corre: 90 casos, rapido)',
                          '')
_ns_ic = {'__name__': 'carga_indice_certificado', 'OC': _OCmod}
exec(compile(_src_ic[:_src_ic.index('PAT = re.compile')], '/work/indice_certificado.py', 'exec'), _ns_ic)
assert 'certificado' in _ns_ic and 'soportes' in _ns_ic, "no se cargaron las funciones del certificado"


def piso1(n, p):
    lo, hi, d = _OC['orden'](n, p)
    d2, dA, ds, k1, k0 = _ns_ic['certificado'](n, p)
    S, llena, npar, nimp = _ns_ic['soportes'](n, p, 6 * p + 4)
    dim = (lambda k: npar if k % 2 == 0 else nimp)
    Zpn = sum(dim(p + j) - len(S[p + j]) for j in range(p - 1))
    Zn = dim(2 * p - 1) - len(S[2 * p - 1])
    return {'d': int(d), 'e': (int(lo) if lo == hi else None), 'e_inf': int(lo), 'e_sup': (int(hi) if hi else None),
            'vp_indice': int(dA), 'Z_pn': int(Zpn), 'Z_n': int(Zn),
            'cert_K': (None if ds is None else int(ds - dA))}
