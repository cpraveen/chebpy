'''
Solves viscous Burgers' equation
See: Stewart, Python for Scientists, Section 9.9
'''
import numpy as np
from scipy.integrate import odeint 
import matplotlib.pyplot as plt
import spectral as cheb

c, mu = 1.0, 0.1
N = 64
x = cheb.chebpts(N)
tau1, tau2 = N**2, N**2
t_initial, t_final = -2.0, 2.0

def u_exact(t,x):
    """ Exact kink solution of Burgers’ equation. """ 
    return c*(1.0 + np.tanh(c*(c*t - x)/(2*mu)))

def mu_ux_exact(t,x):
    """ Return mu*du/dx for exact solution. """ 
    arg = np.tanh(c*(c*t-x)/(2*mu))
    return 0.5*c*c*(arg*arg - 1.0)

def f(x):
    """ Return initial data. """ 
    return u_exact(t_initial,x)

def g1(t):
    """ Return function needed at left boundary. """ 
    return (u_exact(t,-1.0))**2 - mu_ux_exact(t,-1.0)

def g2(t):
    """ Return function needed at right boundary. """
    return mu_ux_exact(t,1.0)

def rhs(u, t):
    """ Return du/dt. """
    u_cheb    = cheb.tocheb(u,x)
    ux_cheb   = cheb.diffcheb(u_cheb)
    uxx_cheb  = cheb.diffcheb(ux_cheb)
    ux        = cheb.fromcheb(ux_cheb,x)
    uxx       = cheb.fromcheb(uxx_cheb,x)
    dudt      = -u*ux + mu*uxx
    dudt[0]  -= tau1*(u[0]**2 - mu*ux[0] - g1(t))
    dudt[-1] -= tau2*(mu*ux[-1] - g2(t))
    return dudt

t = np.linspace(t_initial,t_final,81)
u_initial = f(x)
sol = odeint(rhs,u_initial,t,rtol=10e-12,atol=1.0e-12,mxstep=5000)
xg,tg = np.meshgrid(x,t)
ueg = u_exact(tg,xg)
err = sol - ueg
print("With %d points error is %e" % (N,np.max(np.abs(err))))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xg,tg,sol,rstride=1,cstride=2,alpha=0.9)
ax.set_xlabel('x',style='italic')
ax.set_ylabel('t',style='italic')
ax.set_zlabel('u',style='italic')

plt.show()
