"""Recuento sobre capa2_OUT.txt, bloque (B): en cuantos de los 38 pares se cumple W_{P+2} = W_P.

La version 1 de la nota decia "en los 38"; la salida archivada da otra cosa, y este guion imprime
la cifra que la version 2 cita, con las filas que fallan.  Solo lee capa2_OUT.txt.
"""
import re
from pathlib import Path

bloque = Path(__file__).with_name("capa2_OUT.txt").read_text(encoding="utf-8").split("=== (B)")[1]
fila = re.compile(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+\d+\s+\|\s+\[([^\]]*)\]")
total = cumple = cumple_p3 = total_p3 = 0
for linea in bloque.splitlines():
    m = fila.match(linea)
    if not m:
        continue
    mm, q, p, P = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))
    W = [int(x) for x in m.group(5).split(",")]
    total += 1
    total_p3 += p == 3
    if W[P + 2] == W[P]:
        cumple += 1
        cumple_p3 += p == 3
    else:
        print("no cumple: m=%d q=%d p=%d  W_P=%d  W_{P+2}=%d" % (mm, q, p, W[P], W[P + 2]))
print("W_{P+2} = W_P en %d de %d pares; con p=3: %d de %d" % (cumple, total, cumple_p3, total_p3))
