import time
def one(t,r,b):
    N=t+2*r
    K = CyclotomicField(t) if t>2 else QQ
    zeta = K.zeta(t) if t>2 else K(-1)
    R = LaurentPolynomialRing(K, ['z%d'%j for j in range(r)])
    zs = R.gens()
    alpha=[K(zeta)**k for k in range(t)]+[x for j in range(r) for x in (zs[j],zs[j]**-1)]
    M=matrix(R,N,N,lambda i,j: alpha[i]**b[j])
    t0=time.time(); M.determinant()==0; return time.time()-t0
for (t,r,b) in [(2,1,(7,5,2,0)),(4,1,(9,7,4,2,1,0)),(2,2,(11,8,6,4,1,0)),
                (4,2,(12,10,8,6,3,2,1,0)),(6,2,(13,11,9,7,5,4,2,1,1,0)[:10]),
                (4,3,(13,11,9,7,6,4,3,2,1,0))]:
    N=t+2*r
    bb=tuple(sorted(set(b),reverse=True))[:N]
    if len(bb)<N: continue
    print("t=%d r=%d N=%d  %.3f s"%(t,r,N,one(t,r,bb)), flush=True)
