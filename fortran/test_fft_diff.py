import numpy as np
from scipy.fftpack import diff

def fd(u):
    """ Return 2*dx* finite-difference x-derivative of u. """
    ud = np.empty_like(u)
    ud[1:-1] = u[2: ] - u[ :-2]
    ud[0] = u[1] - u[-1]
    ud[-1] = u[0]-u[-2]
    return ud

for N in [4,8,16,32,64,128,256]:
    dx = 2.0*np.pi/N
    x = np.linspace(0,2.0*np.pi,N,endpoint=False)
    u = np.exp(np.sin(x))
    du_ex = np.cos(x)*u
    du_sp = diff(u)
    du_fd = fd(u)/(2.0*dx)
    err_sp = np.max(np.abs(du_sp-du_ex))
    err_fd = np.max(np.abs(du_fd-du_ex))
    print("N=%5d   err_sp=%12.4e   err_fd=%12.4e" % (N,err_sp,err_fd))
