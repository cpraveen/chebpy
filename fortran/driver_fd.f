      PROGRAM TEST
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(M=4,N=8)
      DIMENSION X(0:N),C(0:N,0:N,0:M),D(0:N,0:M)
      IO = 6
      DO 10 J=0,N
         X(J) = J
 10   CONTINUE
      CALL WEIGHTS (0.0D0,X,N,M,C)
      DO 30 K=1,M
         WRITE(IO,'(" Derivative order = ",I3)') K
         DO 20 I=1,N
            WRITE(IO,40) (C(J,I,K),J=0,I)
 20      CONTINUE
         WRITE(IO,*)
 30   CONTINUE
 40   FORMAT (1X,9F8.3)
      CALL WEIGHTS1 (0.0D0,X,N,M,D)
      DO 50 K=0,M
         WRITE(IO,40) (D(J,K),J=0,N)
 50   CONTINUE
      STOP
      END
