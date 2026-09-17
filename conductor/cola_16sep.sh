#!/bin/bash
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# SEGUNDA SESION DE LOS BARRIDOS m = 6, 7, 8 (celdas_m6/7/8_OUT.txt), EN SERIE Y REANUDABLE.  Se
# relanza con la misma orden y sigue donde se quedo (correr_reanudable.sh + casillas_hechas.py).
#
# La primera sesion de estos logs la escribio correr_por_casilla.sh.  En ella algunas casillas
# murieron por memoria (137): eran los d^m monomios de L.order, que factorizar_conductor.sage ya no
# usa (ver orden_generado.sage).  Aqui las fallidas de la sesion VIEJA se reintentan, las hechas se
# saltan, y un fallo nuevo no se reintenta al relanzar.
#
# (La cola original de ese dia tenia antes dos pasos mas, sobre otros logs que este paquete no usa;
#  aqui solo queda el paso que escribio estos tres.)
#
# Uso:
#   MSYS_NO_PATHCONV=1 docker run -d --name cola --memory=8g \
#     -v "$PWD:/work" -w /work sage-normaliz:local \
#     bash cola_16sep.sh
# Progreso:  python casillas_hechas.py <log> --resumen
#
# Authors: Carles Marin, Claude (AI assistant).
cd /work
export TOPE_SEG="${TOPE_SEG:-1800}"
for M in 6 7 8; do
  case $M in
    6) LISTA="6:13,6:15,6:16,6:17,6:19,6:20,6:21,6:23,6:24,6:25,6:27,6:28,6:29,6:31,6:32,6:33,6:35,6:36,6:37,6:39,6:40,6:41,6:44,6:45,6:48,6:51,6:52,6:55,6:56,6:57,6:60,6:63,6:64,6:68,6:72,6:75,6:76,6:80,6:84,6:88,6:96,6:100,6:108,6:120,6:132" ;;
    7) LISTA="7:15,7:16,7:17,7:19,7:20,7:21,7:23,7:24,7:25,7:27,7:28,7:29,7:31,7:32,7:33,7:35,7:36,7:37,7:39,7:40,7:41,7:44,7:45,7:48,7:51,7:52,7:55,7:56,7:57,7:60,7:63,7:64,7:68,7:72,7:75,7:76,7:80,7:84,7:88,7:96,7:100,7:108,7:120,7:132" ;;
    8) LISTA="8:17,8:19,8:20,8:21,8:23,8:24,8:25,8:27,8:28,8:29,8:31,8:32,8:33,8:35,8:36,8:37,8:39,8:40,8:41,8:44,8:45,8:48,8:51,8:52,8:55,8:56,8:57,8:60,8:63,8:64,8:68,8:72,8:75,8:76,8:80,8:84,8:88,8:96,8:100,8:108,8:120,8:132" ;;
  esac
  REINTENTAR_FALLOS=viejos MODO=stdin bash correr_reanudable.sh factorizar_conductor.sage "$LISTA" "celdas_m${M}_OUT.txt"
done
echo "== cola terminada: $(date -u +%Y-%m-%dT%H:%M:%SZ) =="
