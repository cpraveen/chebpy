from numpy import pi,arange,cos,zeros,ones,mod

# (N+1)-point Clenshaw-Curtiss quadrature
def clencurt(N):
    theta = pi*arange(N+1)/N
    x = -cos(theta)
    w = zeros(N+1)
    v = ones(N-1)
    if mod(N,2) == 0:
        w[0] = 1.0/(N**2-1)
        w[-1] = w[0]
        for k in range(1,N//2):
            v = v - 2*cos(2*k*theta[1:-1])/(4*k**2-1)
        v = v - cos(N*theta[1:-1])/(N**2-1)
    else:
        w[0] = 1.0/N**2
        w[-1]= w[0]
        for k in range(1,(N-1)//2+1):
            v = v - 2*cos(2*k*theta[1:-1])/(4*k**2-1)
    w[1:-1] = 2*v/N
    return x,w
