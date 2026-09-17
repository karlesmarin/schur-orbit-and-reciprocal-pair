#!/bin/bash
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# Casillas que deciden R2 (prediccion escrita antes de los datos), elegidas
# ANTES de calcularlas:
#   asesinas (dim F_p[<p>].S fuera de {1,f-1,f}): 19:156 14:116 15:93 ; dim 0 (S = 0 mod p): 10:63 16:99
#   controles dim = f (predicen c_p != P^2): 4:36 5:33 5:44 6:52 8:51 9:57 9:76 12:75 12:100 13:108
# Reanudable.  Lectura: python reparto_p2.py ; python perfil_W.py
#
#   MSYS_NO_PATHCONV=1 docker run -d --name r2 --memory=8g \
#     -v "$PWD:/work" -w /work sage-normaliz:local bash r2_barrido.sh
#
# Authors: Carles Marin, Claude (AI assistant).
cd /work
TOPE_SEG=2400 MODO=fichero bash correr_reanudable.sh h2_residuo.sage \
  "10:63,4:36,5:33,5:44,6:52,8:51,9:57,9:76,12:75,12:100,13:108,19:156,14:116,15:93,16:99" r2_OUT.txt
