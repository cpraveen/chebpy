'''
Solve
   u_t + 2*pi*u_x = 0 for x in (0,2*pi)
with periodic bc.
'''
import numpy as np
from scipy.fftpack import diff
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def u_exact(t,x):
    """ Exact solution. """
    return np.exp(np.sin(x-2*np.pi*t))

def rhs(u, t):
    """ Return rhs. """
    return -2.0*np.pi*diff(u)

N = 32
x = np.linspace(0,2*np.pi,N,endpoint=False)
u0 = u_exact(0,x)
t_initial, t_final = 0.0, 1.0
t = np.linspace(t_initial,t_final,101)
sol = odeint(rhs,u0,t,mxstep=5000)

# Compute error at final time
u_ex = u_exact(t[-1],x)
err = np.abs(np.max(sol[-1,: ]-u_ex))
print("With %d Fourier nodes the final error = %g" % (N, err))

plt.figure()
plt.plot(x, sol[ 0,:], label="Initial")
plt.plot(x, sol[-1,:], label="Final")
plt.xlabel("x"); plt.ylabel("u")
plt.legend()

# Plot the surface
x_gr, t_gr = np.meshgrid(x, t)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(t_gr, x_gr, sol, vmin=sol.min())
ax.elev, ax.azim = 47, -137
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u')

plt.show()
