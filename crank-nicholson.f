c     Aquest programa fa calor
      program cranknicholson
      implicit none
      integer nx,j,compta
      double precision h,Lx,T(0:140),k,delta(0:140)
      double precision T_nova(0:140),eps,alfa,delta_t
      double precision lower(0:140),upper(0:140),central(0:140)

      h=0.25d0
      Lx=35.d0
      nx=int(Lx/h)
      k=5.4d0

c     Condicions de contorn
      T(0)=15.d0
      T(nx)=55.d0

c     Distribucio inicial de temperatures
      do j=1,nx-1
       T(j)=30.d0
      enddo

c     Distribucio estacionaria
      call subroutine gauseidel(T,T_nova,h,nx,delta)
      eps=maxval(delta) 
      do while (eps.lt.(10.d-6))
       call subroutine gauseidel(T,T_nova,h,nx,delta)
       eps=maxval(delta)
      enddo
      open(1,file='guaseidel.dat')
      
      do j=0,nx
       write(1,*) dble(j)*h,T_nova(j)
      enddo

c     Evolucio temporal
      T=T_nova
      T(0)=0.d0
      T(nx)=0.d0
      delta_t=2.d-3
      alfa=(k*delta_t)/h
      do j=0,nx
      upper(j)=-alfa
      lower(j)=-alfa
      central(j)=(1.d0+2.d0*alfa)
      enddo
      upper(0)=0.d0
      lower(nx)=0.d0
      do 
      end program

c     **************  
c     *Gauss-Seidel*
c     **************
      subroutine gauseidel(u_old,u,h,nx,delta)
      integer nx
      double precision u_old(0:nx),u(0:nx),h
      double precision delta(0:nx)
      integer j
      u(0)=u_old(0)
      u(nx)=u_old(nx)
      delta(0)=0.d0
      delta(nx)=0.d0
      do j=1,nx-1
       u(j)=(1.d0/2.d0)*(u_old(j+1)+u_old(j-1)-2.d0*u_old(j)+
     + -(h**2)*ro(j*h))
       delta(j)=dabs(u(j)-u_old(j))
      enddo

      u_old=u
      return
      end subroutine

c     ******************
c     *Terme inhomogeni*
c     ******************
      double precision function ro(x)
      double precision x
      ro=4.9d0*dexp(-(dabs(x-6.d0))**2.d0)+
     +2.1d0*dexp(-((dabs(x-16.d0))**2.d0)/(0.3d0**2))
      return
      end function

c     ********************
c     *Matriu tridiagonal*
c     ********************

c     Solves the problem T psi =R

c     where T is a tridiagonal matrix, A (lower), B (central), C (upper)
c     A   0 a1, a2, a3, ...., aIMAX
c     B   b1 b2, b3, b4, ...., bIMAX
c     C   c1, c2, c3, c4, ....,0 

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
11    CONTINUE

      DO 12 J=IMAX-1,1,-1
       PSI(J)=PSI(J)-GAM(J+1)*PSI(J+1)
12    CONTINUE

       RETURN
       END