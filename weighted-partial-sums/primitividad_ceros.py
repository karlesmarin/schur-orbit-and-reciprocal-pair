# -*- coding: utf-8 -*-
r"""primitividad_ceros.py -- el Corolario 26 sin la hipotesis de primitividad, contado a mano.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

LA ACUSACION (informe adversario del 19-sep-2026, F01).  El Teorema 25 da la formula cerrada
T = (f/2) B_{1,chi} (chi(2)^{-1} - 1) y anade PRIMITIVO para concluir B_{1,chi} != 0.  El
Corolario 26 quita esa hipotesis, se queda con "f >= 5 impar", y cuenta TODOS los caracteres
impares: 0 si -1 esta en <2>, y d/t si no.  Testigo propuesto: f=21, el caracter impar inducido
desde el modulo 3, con chi(2) = -1 y sin embargo T = 0.

ESTO NO SE JUZGA LEYENDO.  Se enumeran los caracteres impares modulo f de verdad, se evalua
T(m,chi) = suma_{c<=m} c chi(c) con m = (f-1)/2 exactamente en Z[zeta], y se compara el numero de
ceros con lo que el corolario predice.  Se separan ademas los primitivos, que es donde el teorema
si aplica: si la formula acierta en los primitivos y falla en los imprimitivos, la reparacion es
la hipotesis y no el enunciado.

Salida: gates/capas_tipo/primitividad_ceros_OUT.txt
"""
import sys
from itertools import product

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def factoriza(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def orden(a, n):
    k, x = 1, a % n
    while x != 1:
        x = x * a % n
        k += 1
    return k


def raiz_primitiva(p, e):
    """generador de (Z/p^e)^x para p impar."""
    m = p ** e
    phi = (p - 1) * p ** (e - 1)
    fac = list(factoriza(phi))
    for g in range(2, m):
        if g % p == 0:
            continue
        if all(pow(g, phi // q, m) != 1 for q in fac):
            return g
    raise RuntimeError("sin raiz primitiva")


class Caracteres(object):
    """Los caracteres de Dirichlet modulo f impar, como vectores de exponentes.

    (Z/f)^x = producto de ciclicos <g_i> de orden n_i (f impar, asi que cada p^e aporta uno).
    Un caracter es una tupla (a_1..a_k) con 0 <= a_i < n_i y chi(x) = zeta_{n_i}^{a_i * log_i(x)}.
    Se representa el VALOR como un racional del circulo: la suma de a_i*log_i / n_i modulo 1.
    Asi T(m,chi) se calcula como un vector de multiplicidades por raiz de unidad, exacto.
    """

    def __init__(self, f):
        self.f = f
        self.comp = []           # (modulo p^e, generador, orden)
        for p, e in factoriza(f).items():
            m = p ** e
            self.comp.append((m, raiz_primitiva(p, e), (p - 1) * p ** (e - 1)))
        # logaritmo discreto en cada componente
        self.log = []
        for m, g, n in self.comp:
            tabla, x = {}, 1
            for k in range(n):
                tabla[x] = k
                x = x * g % m
            self.log.append(tabla)
        self.ordenes = [n for _, _, n in self.comp]

    def todos(self):
        return product(*[range(n) for n in self.ordenes])

    def valor(self, a, x):
        """chi_a(x) como fraccion del circulo (num, den), o None si gcd(x,f) > 1."""
        from fractions import Fraction
        if any(x % m_p == 0 for m_p, _, _ in
               [(p, 0, 0) for p in factoriza(self.f)]):
            pass
        s = Fraction(0)
        for (m, g, n), tab, ai in zip(self.comp, self.log, a):
            xr = x % m
            if xr not in tab:
                return None
            s += Fraction(ai * tab[xr], n)
        return s % 1

    def impar(self, a):
        v = self.valor(a, self.f - 1)
        from fractions import Fraction
        return v == Fraction(1, 2)

    def conductor(self, a):
        """el conductor de chi_a: producto de p^e sobre las componentes no triviales, con el
        exponente minimo que la soporta.  Para f libre de cuadrados basta mirar si a_i = 0."""
        c = 1
        for (m, g, n), ai in zip(self.comp, a):
            if ai == 0:
                continue
            p = min(factoriza(m))
            e = 1
            # el orden de chi en esa componente
            from math import gcd
            o = n // gcd(ai, n)
            # chi es trivial en 1 + p^k Z sii o | p^{e-k}(p-1) ... para e=1 basta p
            mm, ee = p, 1
            while (p - 1) * p ** (ee - 1) % o != 0:
                ee += 1
                mm *= p
            c *= mm
        return c


def T(ch, a, m):
    """suma_{c=1..m} c chi_a(c), exacta: como diccionario raiz -> coeficiente entero."""
    from fractions import Fraction
    acc = {}
    for c in range(1, m + 1):
        v = ch.valor(a, c)
        if v is None:
            continue
        acc[v] = acc.get(v, 0) + c
    return acc


def es_cero(acc, ch):
    """La suma de coef_j * zeta^{e_j} es cero?  Se reduce con el ciclotomico del orden comun."""
    from fractions import Fraction
    if not acc:
        return True
    n = 1
    for v in acc:
        n = n * v.denominator // __import__("math").gcd(n, v.denominator)
    # coeficientes sobre Z[x]/(x^n - 1), reducidos modulo el n-esimo ciclotomico
    import sympy
    x = sympy.symbols("x")
    poly = sum(co * x ** int(v * n) for v, co in acc.items())
    return sympy.rem(sympy.Poly(poly, x), sympy.Poly(sympy.cyclotomic_poly(n, x), x)) == 0


def main():
    print("=" * 92)
    print("EL COROLARIO 26 SIN PRIMITIVIDAD -- conteo exacto de ceros")
    print("=" * 92)
    print("\n  f    d    t  predice  ceros  primit.  imprim.  veredicto")
    print("  " + "-" * 74)
    fallos, filas = 0, []
    for f in range(5, 122, 2):
        if f % 2 == 0:
            continue
        ch = Caracteres(f)
        phi = 1
        for n in ch.ordenes:
            phi *= n
        d = phi // 2
        t = orden(2, f) if __import__("math").gcd(2, f) == 1 else None
        if t is None:
            continue
        # -1 en <2>?
        en = any(pow(2, k, f) == f - 1 for k in range(1, t + 1))
        predice = 0 if en else d // t
        ceros_p, ceros_i = 0, 0
        m = (f - 1) // 2
        for a in ch.todos():
            if not ch.impar(a):
                continue
            if es_cero(T(ch, a, m), ch):
                if ch.conductor(a) == f:
                    ceros_p += 1
                else:
                    ceros_i += 1
        total = ceros_p + ceros_i
        ok = (total == predice)
        filas.append((f, d, t, predice, total, ceros_p, ceros_i, ok))
        if not ok:
            fallos += 1
        print("  %3d  %3d  %3d  %5d  %5d    %5d    %5d  %s" %
              (f, d, t, predice, total, ceros_p, ceros_i, "ok" if ok else "*** FALLA ***"))

    print("\n  " + "-" * 74)
    print("  %d de %d modulos impares con discrepancia" % (fallos, len(filas)))
    # EL DESGLOSE PRIMO/COMPUESTO SE IMPRIME, porque es la cifra que el papel cita.  Un numero que
    # el papel da y que ninguna salida archivada contiene es una afirmacion sin respaldo, aunque
    # yo lo haya calculado: lo que no esta en la salida no lo puede reencontrar nadie.
    pr = [r for r in filas if len(factoriza(r[0])) == 1 and max(factoriza(r[0]).values()) == 1]
    co = [r for r in filas if r not in pr]
    print("  desglose: %d primos y %d compuestos; la formula acierta en los %d primos"
          % (len(pr), len(co), len(pr)))
    print("            y falla en %d de los %d compuestos" % (len([r for r in co if not r[7]]),
                                                              len(co)))

    # la pregunta que decide la reparacion: falla alguna vez SOLO en los primitivos?
    solo_prim = [r for r in filas if r[5] != (0 if r[3] == 0 and r[6] == 0 else r[3] - r[6])]
    print("\n  LA REPARACION.  Si todos los ceros de mas son imprimitivos, la hipotesis que falta")
    print("  es la primitividad y el enunciado se salva restringiendola:")
    malos = [r for r in filas if not r[7]]
    print("    modulos que fallan: %s" % [r[0] for r in malos][:20])
    print("    de ellos, con algun cero PRIMITIVO de mas: %s" %
          [r[0] for r in malos if r[5] > r[3]] or "ninguno")
    print("    todos los modulos que fallan son compuestos: %s" %
          all(len(factoriza(r[0])) > 1 or max(factoriza(r[0]).values()) > 1 for r in malos))

    print("\n" + "=" * 92)
    print("TOTAL DE DISCREPANCIAS: %d" % fallos)
    print("=" * 92)
    return 1 if fallos else 0


if __name__ == "__main__":
    sys.exit(main())
