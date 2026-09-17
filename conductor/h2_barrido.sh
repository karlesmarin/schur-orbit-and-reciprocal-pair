#!/bin/bash
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# Barrido que prueba H2 / H2' (prediccion escrita antes de los datos).  Reanudable: relanzar con la misma orden.
#
#   MSYS_NO_PATHCONV=1 docker run -d --name h2 --memory=8g \
#     -v "$PWD:/work" -w /work sage-normaliz:local \
#     bash h2_barrido.sh
# Lectura:  python perfil_W.py ; python datos_indice.py
#
# Authors: Carles Marin, Claude (AI assistant).
cd /work
# h2_lista.txt es la lista con la que se corrio (m in {4,5,9,10} con 2m < q <= 132, y m in {6,7,8}
# con 133 <= q <= 160, grado 2 <= phi(q)/2 <= 30, sin las casillas ya calculadas en otros logs).
# Es la misma que la cabecera de la sesion en h2_OUT.txt.
LISTA="$(head -1 h2_lista.txt)"
if [ -z "$LISTA" ]; then
  echo "lista vacia" >&2
  exit 3
fi
TOPE_SEG=1200 MODO=fichero bash correr_reanudable.sh h2_residuo.sage "$LISTA" h2_OUT.txt
