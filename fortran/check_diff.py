import numpy as np 
import spectral as cheb

for N in [8,16,32,64]: 
    x = cheb.chebpts(N) 
    ex = np.exp(x) 
    u = np.sin(ex) 
    dudx = ex*np.cos(ex)
    uc = cheb.tocheb(u,x)
    duc = cheb.diffcheb(uc) 
    du = cheb.fromcheb(duc,x) 
    err = np.max(np.abs(du-dudx))
    print("with N = %d error is %e" % (N, err))
