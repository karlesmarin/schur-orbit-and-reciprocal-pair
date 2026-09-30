"""Rangos completos de tres medidas del articulo del conductor (punto 5 del informe T75, 28-sep-2026).

Lee las salidas archivadas, no recalcula nada:
  capa2_OUT.txt  seccion (A): la segunda capa como cuadrado de la primera (238 pares, p impar, p no divide N);
  l2r_OUT.txt    END-TO-END: lema de duplicacion sobre las capas (304 casos, q' = 2r);
  medio_OUT.txt  seccion (B): medio periodo con 4 | q' (156 casos distintos en 158 lineas).
Para cada una imprime el recuento (y lo compara con el del texto) y el rango de cada parametro.
"""
import re
from pathlib import Path

AQUI = Path(__file__).resolve().parent
fallos = 0


def rango(nombre, valores):
    v = sorted(set(valores))
    if len(v) <= 16:
        return "%s en %s" % (nombre, v)
    return "%s en [%s, %s] (%d valores distintos)" % (nombre, v[0], v[-1], len(v))


def informe(titulo, filas, esperado, claves):
    global fallos
    distintas = sorted(set(filas))
    ok = len(distintas) == esperado
    fallos += not ok
    print("== %s: %d lineas, %d distintas (texto: %d) %s" % (titulo, len(filas), len(distintas), esperado,
                                                           "OK" if ok else "FALLA"))
    for i, k in enumerate(claves):
        if k:
            print("   ", rango(k, [f[i] for f in distintas]))
    return distintas


# (A) de capa2: m q p q' d W_1 W_2 dim(V*V) | ...
t = (AQUI / "capa2_OUT.txt").read_text(encoding="utf-8")
sec = t[:t.index("=== (B)")]
filas = [tuple(int(x) for x in m.groups()) for m in
         re.finditer(r"^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\d+\s+\|", sec, flags=re.M)]
d = informe("segunda capa como cuadrado", filas, 238, ["m", "q", "p", "q'", "d"])
print("    familias (q' divide 2m, 2m+1 o 2m+2):",
      {k: sum(1 for f in d if (2 * f[0] + k) % f[3] == 0) for k in (0, 1, 2)})
print("    v = v_p(q):", sorted({next(v for v in range(1, 40) if f[1] % f[2] ** (v + 1)) for f in d}))

# l2r END-TO-END
t = (AQUI / "l2r_OUT.txt").read_text(encoding="utf-8")
filas = [(int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(4), int(m.group(5)), int(m.group(6)),
          int(m.group(7)), int(m.group(8))) for m in
         re.finditer(r"^r=\s*(\d+) q'=(\d+) p=(\d+) m=(\S+)\s+N/B=(\d+) m=(\d+) v=(\d+) q=(\d+) \|", t, flags=re.M)]
d = informe("duplicacion, capas", filas, 304, ["r", "q'", "p", "familia", "N/B", "m", "v", "q"])

# medio (B)
t = (AQUI / "medio_OUT.txt").read_text(encoding="utf-8")
sec = t[t.index("=== (B)"):t.index("=== (C)")]
filas = [tuple(int(x) for x in m.groups()) for m in
         re.finditer(r"^q'=\s*(\d+) p=\s*(\d+) B=(\d+) m=\s*(\d+) v=(\d+) q=\s*(\d+) \|", sec, flags=re.M)]
d = informe("medio periodo, capas", filas, 156, ["q'", "p", "B", "m", "v", "q"])
print("    familias: m = Br", sum(1 for f in d if (2 * f[3]) % f[0] == 0 and (2 * f[3] // f[0]) % 2),
      "| m = Br-1", sum(1 for f in d if (2 * f[3] + 2) % f[0] == 0))
print("fallos:", fallos)
