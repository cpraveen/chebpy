from numpy import pi,cos,arange,ones,tile,dot,eye,diag,linspace

def cheb(N):
    '''Chebushev polynomial differentiation matrix.
       Ref.: Trefethen's 'Spectral Methods in MATLAB' book.
    '''
    x      = cos(pi*arange(0,N+1)/N)
    if N%2 == 0:
        x[N/2] = 0.0 # only when N is even!
    c      = ones(N+1)
    c[0]   = 2.0
    c[N]   = 2.0
    c      = c * (-1.0)**arange(0,N+1)
    X      = tile(x.reshape(N+1,1), (1,N+1))
    dX     = X - X.T
    D      = dot(c.reshape(N+1,1),(1.0/c).reshape(1,N+1))
    D      = D / (dX+eye(N+1))
    D      = D - diag( D.T.sum(axis=0) )
    return D,x
