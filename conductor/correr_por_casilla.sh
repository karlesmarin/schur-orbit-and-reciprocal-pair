#!/bin/bash
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# UN PROCESO POR CASILLA.  El barrido en un solo proceso murio con OOM (137) tres veces, y las tres
# NO por una casilla cara: (5,44) revento dentro del barrido y sale suelta sin problema.  La memoria
# se ACUMULA entre casillas (pila de PARI, cacheo de cuerpos).  Un proceso por casilla la acota.
#
# Uso (dentro del contenedor):  bash correr_por_casilla.sh "<lista m:q,...>" <fichero de salida>
#
# Cada casilla se anade con >> y con su propia cabecera; el fichero es un LOG, no una tabla, y los
# lectores de filas (reparto.py) ya parsean por linea.
#
# Authors: Carles Marin, Claude (AI assistant).
LISTA="$1"
SALIDA="$2"
echo "== un proceso por casilla ==  lista: $LISTA" > "$SALIDA"
echo "== inicio: $(date -u +%Y-%m-%dT%H:%M:%SZ) ==" >> "$SALIDA"
IFS=',' read -ra CELDAS <<< "$LISTA"
for celda in "${CELDAS[@]}"; do
  # El salto de linea de delante NO es cosmetico: la salida de sage no termina en newline, asi que
  # sin el, el marcador queda PEGADO a mitad de la linea anterior y un `grep -c '^--- casilla'` lo
  # cuenta como 0.  El 16-sep eso me hizo leer "1 casilla hecha" en un barrido de 50 ya terminado.
  echo "" >> "$SALIDA"
  echo "--- casilla $celda ---" >> "$SALIDA"
  # timeout por casilla: una que se cuelgue no puede llevarse el resto del barrido por delante.
  CASILLAS="$celda" timeout "${TOPE_SEG:-900}" sage < factorizar_conductor.sage >> "$SALIDA" 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "*** casilla $celda TERMINO CON CODIGO $rc (137=OOM, 124=timeout) -- NO es un nulo ***" >> "$SALIDA"
  fi
done
echo "== fin: $(date -u +%Y-%m-%dT%H:%M:%SZ) ==" >> "$SALIDA"
