# -*- coding: utf-8 -*-
# Envoltorio minimo para que un guion de Sage use las poblaciones de peel_zero sin que el
# preparser de Sage toque el modulo original.  No calcula nada: solo reexporta.
#
# Authors: Carles Marin, Claude (AI assistant).

from peel_zero import betas as _betas, occupied as _occupied, phi_zero as _phi_zero


def betas_py(t, r, W):
    return _betas(t, r, W)


def occupied_py(b, t):
    return _occupied(b, t)


def phi_zero_py(b, t, r):
    return _phi_zero(b, t, r)
