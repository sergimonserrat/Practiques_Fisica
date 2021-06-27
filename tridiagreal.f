c 
c Solves the problem T psi =R
c 
c where T is a tridiagonal matrix, A (lower), B (central), C (upper)
c  A   0 a1, a2, a3, ...., aIMAX
c  B   b1 b2, b3, b4, ...., bIMAX
c  C   c1, c2, c3, c4, ....,0 

        SUBROUTINE TRIDIAG(A,B,C,R,PSI,IMAX)
        IMPLICIT double precision (A-H,K,O-Z)
        IMPLICIT INTEGER (I-J , L-N)
        double precision  BET
        double precision  GAM(4001)
        double precision A(IMAX),B(IMAX),C(IMAX),R(IMAX),PSI(IMAX)

        IF(B(1).EQ.0.) PAUSE
        BET=B(1)
        PSI(1)=R(1)/BET
        DO 11 J=2,IMAX
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF(BET.EQ.0) PAUSE
        PSI(J)=(R(J)-A(J)*PSI(J-1))/BET
11      CONTINUE

        DO 12 J=IMAX-1,1,-1
        PSI(J)=PSI(J)-GAM(J+1)*PSI(J+1)
12      CONTINUE

       RETURN
       END

