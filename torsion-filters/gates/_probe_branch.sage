# sonda: ¿existe la regla Sp_8 -> Sp_4 x Sp_4 en esta imagen?
for nombre in ("C4", "C3"):
    for sub in ("C2xC2", "C1xC2"):
        try:
            br = branching_rule(nombre, sub, "orthogonal_sum")
            W = WeylCharacterRing(nombre); X = WeylCharacterRing(sub)
            v = [0] * int(nombre[1])
            v[0] = 1
            el = W(W.space().from_vector(vector([Integer(u) for u in v])))
            print(nombre, "->", sub, ":", el.branch(X, rule=br))
        except Exception as e:
            print(nombre, "->", sub, "FALLA:", type(e).__name__, str(e)[:90])
print("DONE")
