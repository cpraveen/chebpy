from numpy import sqrt,arange,diag,sort,argsort
from numpy.linalg import eig

# Return n-point gauss quadrature
#   x[n] : quadrature points in (-1,1)
#   w[n] : Quadrature weights, summing to 2
def gauss(n):
    beta = 0.5/sqrt(1 - 1/(2*arange(1,n))**2)
    T = diag(beta,-1) + diag(beta,+1)
    x, V = eig(T); i = argsort(x); x = sort(x)
    w = 2*V[0,i]**2
    return x, w
