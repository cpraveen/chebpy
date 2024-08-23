C     Fourier and Chebyshev codes from
C        Fornberg: Practical Guide to Pseudospectral Methods

      SUBROUTINE FFT (A,B,IS,N,ID)
C-- +--------------------------------------------------------------
C-- | A CALL TO FFT REPLACES THE COMPLEX DATA VALUES A(J)+i B(J), 
C-- | J=0,1,... ,N-1 WITH THEIR TRANSFORM
C-- |
C-- |                            2 i ID PI K J / N
C-- |    SUM      (A(K) + iB(K))e ,                 J=0,1,.. .,N-1
C-- |  K=0..N-1
C-- |
C-- | INPUT AND OUTPUT PARAMETERS
C-- |    A    ARRAY A (0: *), REAL PART OF DATA/TRANSFORM
C-- |    B    ARRAY B (0: *), IMAGINARY PART OF DATA/TRANSFORM
C-- | INPUT PARAMETERS
C-- |    IS   SPACING BETWEEN CONSECUTIVE ELEMENTS IN A AND B
C-- |         (USE IS=+1 FOR ELEMENTS STORED CONSECUTIVELY)
C-- |    N    NUMBER OF DATA VALUES, MUST BE A POWER OF TWO
C-- |    ID   USE +1 OR -1 TO SPECIFY DIRECTION OF TRANSFORM
C-- +--------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(0:*),B(0:*)
Cf2py intent(in, out) A, B
Cf2py intent(in) IS, ID
Cf2py integer intent(hide), depend(A) :: N = len(A)
      J=0
C---  APPLY PERMUTATION MATRIX ----
      DO 20 I=0,(N-2)*IS,IS
         IF (I.LT.J) THEN 
            TR = A(J) 
            A(J) = A(I) 
            A(I) = TR
            TI = B(J) 
            B(J) = B(I) 
            B(I) = TI
         ENDIF
         K = IS*N/2
 10      IF (K.LE.J) THEN
            J = J - K
            K = K/2
            GOTO 10
         ENDIF
         J = J + K
 20   CONTINUE
C---  PERFORM THE LOG2 N MATRIX-VECTOR MULTIPLICATIONS ---
      S = 0.0D0
      C = -1.0D0
      L = IS
 30   LH = L
      L = L + L
      UR = 1.0D0
      UI = 0.0D0
      DO 50 J=0,LH-IS,IS
         DO 40 I=J,(N-1)*IS,L
            IP = I + LH
            TR = A(IP)*UR - B(IP)*UI
            TI = A(IP)*UI + B(IP)*UR
            A(IP) = A(I) - TR
            B(IP) = B(I) - TI
            A(I)  = A(I) + TR
            B(I)  = B(I) + TI
 40      CONTINUE
         TI = UR*S + UI*C
         UR = UR*C - UI*S
         UI = TI
 50   CONTINUE 
      S = SQRT(0.5D0*(1.0D0-C))*ID
      C = SQRT(0.5D0*(1.0D0+C))
      IF(L.LT.N*IS) GOTO 30
      RETURN
      END

      SUBROUTINE FCT (A,X,N,B)
C-- +--------------------------------------------------------------
C-- | A CALL TO FCT PLACES IN B(0:N) THE COSINE TRANSFORM OF THE 
C-- | VALUES IN A(0:N)
C-- |
C-- | B(J) = SUM C(K)*A(K)*COS(PI*K*J/N), J=0,1,...,N,  K=0,...,N
C-- |
C-- | WHERE
C-- |    C(K) = 1.0 FOR K=0,N
C-- |         = 2.0 FOR K=1,2,...,N-1
C-- |
C-- | INPUT PARAMETERS
C-- |    A   A(0:N)   ARRAY WITH INPUT DATA
C-- |    X   X(0:N)   ARRAY WITH CHEBYSHEV GRID POINT LOCATIONS
C-- |                 X(J) = -COS(PI*J/N), J=0,1,...,N
C-- |    N            SIZE OF TRANSFORM - MUST BE A POWER OF TWO
C-- |
C-- | OUTPUT PARAMETER
C-- |    B   B(0:N)   ARRAY FOR RECEIVING TRANSFORM COEFFICIENTS
C-- |                 (MAY BE IDENTICAL WITH THE ARRAY A)
C-- +--------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(0:*),X(0:*),B(0:*) 
Cf2py intent(in) A, X
Cf2py intent(out) B
Cf2py integer intent(hide), depend(A) :: N = len(A)
      N2 = N/2
      A0 = A(N2-1)+A(N2+1)
      A9 = A(1)
      DO 10 I=2,N2-2,2
         A0 = A0+A9+A(N+1-I)
         A1 = A( I+1)-A9
         A2 = A(N+1-I)-A(N-1-I) 
         A3 = A(I)+A(N-I)
         A4 = A(I)-A(N-I)
         A5 = A1-A2
         A6 = A1+A2
         A1 = X(N2-I)
         A2 = X(I)
         A7 = A1*A4+A2*A6
         A8 = A1*A6-A2*A4
         A9 = A(I+1)
         B(I ) = A3+A7
         B(N-I) = A3-A7
         B(I+1 ) = A8+A5 
         B(N+1-I) = A8-A5
 10   CONTINUE
      B(1) = A(0)-A(N)
      B(0) = A(0)+A(N)
      B(N2 ) = 2.D0*A(N2)
      B(N2+1) = 2.D0*(A9-A(N2+1)) 
      CALL FFT(B(0),B(1),2,N2,1) 
      A0 = 2.D0*A0
      B(N) = B(0)-A0 
      B(0) = B(0)+A0 
      DO 20 I=1,N2-1
         A1 = 0.5 D0 *(B(I)+B(N-I))
         A2 = 0.25D0/X(N2+I)*(B(I)-B(N-I)) 
         B(I ) = A1+A2
         B(N-I) = A1-A2
 20   CONTINUE
      RETURN 
      END

      SUBROUTINE TOFOUR(A,B,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(0:N-1),B(0:N-1)
Cf2py intent(in,out) A, B
Cf2py integer intent(hide), depend(A) :: N = len(A)
      CALL FFT(A,B,1,N,-1)
      A1 = 1.0D0/N
      DO 10 I=0,N-1
         A(I) = A(I)*A1
         B(I) = B(I)*A1
 10   CONTINUE
      RETURN
      END

      SUBROUTINE FROMFOUR(A,B,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(0:N-1),B(0:N-1)
Cf2py intent(in,out) A, B
Cf2py integer intent(hide), depend(A) :: N = len(A)
      CALL FFT(A,B,1,N,+1)
      RETURN
      END

      SUBROUTINE DIFFFOUR(A,B,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(0:N-1),B(0:N-1)
Cf2py intent(in,out) A, B
Cf2py integer intent(hide), depend(A) :: N = len(A)
      PI    = 4.0D0*ATAN(1.0D0)
      N2    = N/2
      A(0)  = 0.0D0
      B(0)  = 0.0D0
      A(N2) = 0.0D0
      B(N2) = 0.0D0
      DO 10 I=1,N2-1
         J    = N-I
         F    = I*PI
         A1   = A(I)*F
         A(I) = -B(I)*F
         B(I) = A1
         A1   = A(J)*F
         A(J) = B(J)*F
         B(J) = -A1
 10   CONTINUE
      RETURN
      END

      SUBROUTINE CHEBPTS(X, N)
C     CREATE N+1 CHEBYSHEV DATA POINTS IN X 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(0:N)
Cf2py intent(in) N
Cf2py intent(out) X
Cf2py depend(N) X
      PI = 4.D0*ATAN(1.D0) 
      DO 10 I=0,N
         X(I) = -COS(PI*I/N) 
 10   CONTINUE
      RETURN
      END

      SUBROUTINE FROMCHEB (A,X,N,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(0:N),X(0:N),B(0:N)
Cf2py intent(in) A, X
Cf2py intent(out) B
Cf2py integer intent(hide), depend(A) :: N = len(A) - 1
      B(0) = A(0)
      A1 = 0.5D0
      DO 10 I=1,N-1
         A1 = -A1
         B(I) = A1*A(I) 
 10   CONTINUE
      B(N) = A(N)
      CALL FCT(B,X,N,B)
      RETURN
      END
         
      SUBROUTINE TOCHEB (A,X,N,B) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION A(0:N),X(0:N),B(0:N)
Cf2py intent(in) A, X
Cf2py intent(out) B
Cf2py integer intent(hide), depend(A) :: N = len(A) - 1
      CALL FCT(A,X,N,B) 
      B1 = 0.5D0/N
      B(0) = B(0)*B1 
      B(N) = B(N)*B1
      B1 = 2.D0*B1 
      DO 10 I=1,N-1
         B1 = -B1
         B(I) = B(I)*B1
 10   CONTINUE
      RETURN 
      END

      SUBROUTINE DIFFCHEB (A,N,B) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION A(0:N),B(0:N)
Cf2py intent(in) A
Cf2py intent(out) B
Cf2py integer intent(hide), depend(A) :: N = len(A) - 1
      A1 = A(N)
      A2 = A(N-1)
      B(N) = 0.D0
      B(N-1) = 2.D0*N*A1
      A1 = A2
      A2 = A(N-2)
      B(N-2) = 2.D0*(N-1)*A1 
      DO 10 I=N-2,2,-1
         A1 = A2
         A2 = A(I-1)
         B(I-1) = B(I+1)+2.D0*I*A1
 10   CONTINUE
      B(0) = 0.5D0*B(2)+A2
      RETURN
      END
