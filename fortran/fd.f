C +----------------------------------------------------------------------------+
C | Taken from                                                                 |
C |   Fornberg: Practical Guide to Pseudospectral Methods                      |
C +----------------------------------------------------------------------------+

      SUBROUTINE WEIGHTS (XI,X,N,M,C)
C +----------------------------------------------------------------------------+
C | INPUT PARAMETERS:                                                          |
C |   XI  POINT AT WHICH THE APPROXIMATIONS ARE TO BE ACCURATE                 |
C |   X   X-COORDINATES FOR GRID POINTS, ARRAY DIMENSIONED X(0:N)              |
C |   N   THE GRID POINTS ARE AT X(0),X(1),...,X(N) (I.E. N+1 IN ALL)          |
C |   M   HIGHEST ORDER OF DERIVATIVE TO BE APPROXIMATED                       |
C |                                                                            |
C | OUTPUT PARAMETER:                                                          |
C |   C   WEIGHTS, ARRAY DIMENSIONED C(0:N,0:N,0:M)                            |
C |       ON RETURN, THE ELEMENT C(J,I,K) CONTAINS THE WEIGHT TO BE            |
C |       APPLIED AT X(J) WHEN THE K:TH DERIVATIVE IS APPROXIMATED             |
C |       BY A STENCIL EXTENDING OVER X(0),X(1),...,X(I)                       |
C +----------------------------------------------------------------------------+
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(0:N),C(0:N,0:N,0:M)
      C(0,0,0) = 1.0D0
      C1       = 1.0D0
      C4       = X(0)-XI
      DO 40 I=1,N
         MN = MIN(I,M)
         C2 = 1.0D0
         C5 = C4
         C4 = X(I)-XI
         DO 20 J=0,I-1
            C3 = X(I)-X(J)
            C2 = C2*C3
            IF(I.LE.M) C(J,I-1,I) = 0.0D0
            C(J,I,0) = C4*C(J,I-1,0)/C3
            DO 10 K=1,MN
               C(J,I,K) = (C4*C(J,I-1,K) - K*C(J,I-1,K-1))/C3
 10         CONTINUE
 20      CONTINUE
         C(I,I,0) = -C1*C5*C(I-1,I-1,0)/C2
         DO 30 K=1,MN
            C(I,I,K) = C1*(K*C(I-1,I-1,K-1) - C5*C(I-1,I-1,K))/C2
 30      CONTINUE
         C1 = C2
 40   CONTINUE
      RETURN
      END

      SUBROUTINE WEIGHTS1 (XI,X,N,M,C)
C +----------------------------------------------------------------------------+
C | INPUT PARAMETERS:                                                          |
C |   XI  POINT AT WHICH THE APPROXIMATIONS ARE TO BE ACCURATE                 |
C |   X   X-COORDINATES FOR GRID POINTS, ARRAY DIMENSIONED X(0:N)              |
C |   N   THE GRID POINTS ARE AT X(0),X(1),...,X(N) (I.E. N+1 IN ALL)          |
C |   M   HIGHEST ORDER OF DERIVATIVE TO BE APPROXIMATED                       |
C |                                                                            |
C | OUTPUT PARAMETER:                                                          |
C |   C   WEIGHTS, ARRAY DIMENSIONED C(0:N,0:M)                                |
C |       ON RETURN, THE ELEMENT C(J,K) CONTAINS THE WEIGHT TO BE              |
C |       APPLIED AT X(J) WHEN THE K:TH DERIVATIVE IS APPROXIMATED             |
C |       BY A STENCIL EXTENDING OVER X(0),X(1),...,X(N)                       |
C +----------------------------------------------------------------------------+
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(0:N),C(0:N,0:M)
      C1 = 1.0D0
      C4 = X(0)-XI
      DO 10 K=0,M
         DO 10 J=0,N
            C(J,K) = 0.0D0
 10   CONTINUE
      C(0,0) = 1.0D0
      DO 50 I=1,N
         MN = MIN(I,M)
         C2 = 1.0D0
         C5 = C4
         C4 = X(I)-XI
         DO 40 J=0,I-1
            C3 = X(I)-X(J)
            C2 = C2*C3
            DO 20 K=MN,1,-1
               C(I,K) = C1*(K*C(I-1,K-1) - C5*C(I-1,K))/C2
 20         CONTINUE
            C(I,0) = -C1*C5*C(I-1,0)/C2
            DO 30 K=MN,1,-1
               C(J,K) = (C4*C(J,K) - K*C(J,K-1))/C3
 30         CONTINUE
            C(J,0) = C4*C(J,0)/C3
 40      CONTINUE
         C1 = C2
 50   CONTINUE
      RETURN
      END
