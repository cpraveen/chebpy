import numpy as np 
from spectral import chebpts,tocheb,fromcheb,diffcheb

for N in [8,16,32,64]: 
    x = chebpts(N) 
    ex = np.exp(x) 
    u = np.sin(ex) 
    dudx = ex*np.cos(ex)
    uc = tocheb(u,x)
    duc = diffcheb(uc) 
    du = fromcheb(duc,x) 
    err = np.max(np.abs(du-dudx))
    print("with N = %d error is %e" % (N, err))
