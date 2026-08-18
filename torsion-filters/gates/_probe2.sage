import time
def zero_uni(t,r,b):
    N=t+2*r; B=2*b[0]+1
    K = CyclotomicField(t) if t>2 else QQ
    zeta = K.zeta(t) if t>2 else K(-1)
    R.<w> = LaurentPolynomialRing(K)
    zs=[w**(B**j) for j in range(r)]
    alpha=[K(zeta)**k for k in range(t)]+[x for j in range(r) for x in (zs[j],zs[j]**-1)]
    M=matrix(R,N,N,lambda i,j: alpha[i]**b[j])
    return M.determinant()==0
for (t,r,b) in [(4,2,(12,10,8,6,3,2,1,0)),(6,2,(13,11,9,7,6,5,4,2,1,0)),
                (4,3,(13,11,9,7,6,4,3,2,1,0)),(8,3,(16,14,12,10,8,6,5,4,3,2,1,0,)+(0,)*2)]:
    N=t+2*r; bb=tuple(sorted(set(b),reverse=True))
    if len(bb)<N: print("skip t=%d r=%d"%(t,r)); continue
    bb=bb[:N]; t0=time.time(); z=zero_uni(t,r,bb)
    print("t=%d r=%d N=%d  %.3f s   zero=%s"%(t,r,N,time.time()-t0,z), flush=True)
