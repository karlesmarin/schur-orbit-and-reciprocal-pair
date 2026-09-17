#!/bin/bash
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# BARRIDO REANUDABLE, UN PROCESO POR CASILLA.  Se puede matar (Docker cae, OOM, Ctrl-C) y relanzar
# con LA MISMA orden: sigue donde se quedo.
#
#   - Cada casilla escribe en  $SALIDA.en_curso  y SOLO al terminar su bloque se anade al log con la
#     marca  "=== hecha m:q ==="  o  "=== fallo m:q rc=N ===".  Una casilla cortada no deja filas a
#     medias en el log, asi que al reanudar no hay filas duplicadas que los lectores contarian dos veces.
#   - Al arrancar, casillas_hechas.py lee el log (formato nuevo y el viejo de
#     correr_por_casilla.sh) y dice que saltar.  Las fallidas se saltan tambien: un OOM repetido
#     quema una hora.  Para reintentarlas:  REINTENTAR_FALLOS=1.
#   - El log NUNCA se trunca.  Lo unico que se reescribe es .en_curso.
#
# Uso (dentro del contenedor, en este directorio):
#   bash correr_reanudable.sh <guion.sage> "<m:q,m:q,...>" <log>
# Variables:  TOPE_SEG (por casilla, 900)   MODO=stdin|fichero (stdin: "sage < guion", como los
#             barridos viejos, para que el formato de filas que leen los guiones .py no cambie)
#             REINTENTAR_FALLOS=1 (todas) | viejos (solo las del formato viejo)
#
# Authors: Carles Marin, Claude (AI assistant).
GUION="$1"
LISTA="$2"
SALIDA="$3"
MODO="${MODO:-stdin}"
TOPE="${TOPE_SEG:-900}"
if [ -z "$GUION" ] || [ -z "$LISTA" ] || [ -z "$SALIDA" ]; then
  echo "uso: bash correr_reanudable.sh <guion.sage> \"<m:q,...>\" <log>" >&2
  exit 2
fi
PY="python3"
command -v python3 >/dev/null 2>&1 || PY="sage -python"
EXTRA=""
[ "${REINTENTAR_FALLOS:-0}" = "1" ] && EXTRA="--reintentar-fallos"
[ "${REINTENTAR_FALLOS:-0}" = "viejos" ] && EXTRA="--reintentar-fallos-viejos"
SALTAR="$($PY casillas_hechas.py "$SALIDA" --saltar $EXTRA | tr -d '\r')"
touch "$SALIDA"
echo "" >> "$SALIDA"
echo "== sesion: $(date -u +%Y-%m-%dT%H:%M:%SZ)  guion=$GUION  modo=$MODO  tope=${TOPE}s ==  lista: $LISTA" >> "$SALIDA"
IFS=',' read -ra CELDAS <<< "$LISTA"
total=0; saltadas=0; hechas=0; fallidas=0
for celda in "${CELDAS[@]}"; do
  # quitar espacios y \r: una lista pasada desde PowerShell traia "7:132\r" y la marca salia partida
  celda="$(echo "$celda" | tr -d ' \r\n\t')"
  [ -z "$celda" ] && continue
  total=$((total+1))
  if echo "$SALTAR" | grep -qxF "$celda"; then
    saltadas=$((saltadas+1))
    continue
  fi
  EN_CURSO="$SALIDA.en_curso"
  echo "--- casilla $celda ---  $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$EN_CURSO"
  if [ "$MODO" = "stdin" ]; then
    CASILLAS="$celda" timeout "$TOPE" sage < "$GUION" >> "$EN_CURSO" 2>&1
  else
    CASILLAS="$celda" timeout "$TOPE" sage "$GUION" >> "$EN_CURSO" 2>&1
  fi
  rc=$?
  # salto de linea delante de la marca: la salida de sage no acaba en newline
  echo "" >> "$EN_CURSO"
  if [ $rc -eq 0 ]; then
    echo "=== hecha $celda ===" >> "$EN_CURSO"
    hechas=$((hechas+1))
  else
    echo "=== fallo $celda rc=$rc ===  (137=OOM, 124=timeout) -- NO es un nulo" >> "$EN_CURSO"
    fallidas=$((fallidas+1))
  fi
  cat "$EN_CURSO" >> "$SALIDA"
  rm -f "$EN_CURSO"
done
echo "== fin de sesion: $(date -u +%Y-%m-%dT%H:%M:%SZ)  lista $total  saltadas $saltadas  hechas $hechas  fallidas $fallidas ==" >> "$SALIDA"
$PY casillas_hechas.py "$SALIDA" --resumen >> "$SALIDA"
