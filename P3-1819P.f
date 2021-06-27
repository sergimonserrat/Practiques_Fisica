c     Aquest programa calcula el periode del Big Ben
      program Bigben
      implicit none
      double precision C,w0,g,L,phi0,PI,phi,f,trapezis,simpson,h
      double precision T1,T2,T3,S1,S2,S3,eT1,eT2,eT3
      double precision eS1,eS2,eS3,Terror1,Serror1,Terror2,Serror2
      double precision Terror3,Serror3
      integer m,q
      common/PARAMETRE/phi0
      external f
      g=9.807d0
      L=4.d0
      w0=g/L
      C=4.d0/w0
      PI=dacos(-1.d0)
      h=(PI/2.d0)/2**m
      open(1,file='P3-1819P-res1.dat')
100   format(13(e24.18,1x),I2)
      do m=2,20
       h=(PI/2.d0)/2**m      
       do q=1,3
        if (q.eq.1) then
         phi0=0.d0
         T1=C*trapezis(0.d0,PI/2.d0,m,f)
         S1=C*simpson(0.d0,PI/2.d0,m,f)
         Terror1=C*trapezis(0.d0,PI/2.d0,20,f)
         Serror1=C*simpson(0.d0,PI/2.d0,20,f)
         eT1=dabs(Terror1-T1)
         eS1=dabs(Serror1-S1)
        endif
        if (q.eq.2) then 
         phi0=PI-1.d-1
         T2=C*trapezis(0.d0,PI/2.d0,m,f)
         S2=C*simpson(0.d0,PI/2.d0,m,f)
         Terror2=C*trapezis(0.d0,PI/2.d0,20,f)
         Serror2=C*simpson(0.d0,PI/2.d0,20,f)
         eT2=dabs(Terror2-T2)
         eS2=dabs(Serror2-S2)
        endif
        if (q.eq.3) then
         phi0=PI-1.d-3
         T3=C*trapezis(0.d0,PI/2.d0,m,f)
         S3=C*simpson(0.d0,PI/2.d0,m,f)
         Terror3=C*trapezis(0.d0,PI/2.d0,20,f)
         Serror3=C*simpson(0.d0,PI/2.d0,20,f)
         eT3=dabs(Terror3-T3)
         eS3=dabs(Serror3-S3)
        endif
       write(1,100) h,T1,T2,T3,S1,S2,S3,
     + eT1,eT2,eT3,eS1,eS2,eS3,m
       enddo
      enddo
      close(1)
      end program

c     **********************
c     *Metode dels trapezis*
c     **********************
      double precision function trapezis(a,b,m,fcn)
      double precision a,b,h,extrems,enmig,phi,fcn
      integer N,k
      N=2**m
      h=(b-a)/dble(N)
      
      enmig=0
      do k=1,N-1
       phi=a+(dble(k)*h)
       enmig=enmig+2*fcn(phi)
      enddo
      trapezis=(h/2.d0)*enmig
      return
      end function

c     *******************
c     *Metode de Simpson*
c     *******************
      double precision function simpson(a,b,m,fcn)
      double precision a,b,h,extrems,enmig,phi,fcn
      integer N,k,m
      N=2**m
      h=(b-a)/dble(N)
      extrems=fcn(a)+fcn(b)
      enmig=0
      do k=1,N-1
       phi=a+(dble(k)*h)
       if (mod(k,2).eq.1) then
        enmig=enmig+4*fcn(phi)
       endif
       if (mod(k,2).eq.0) then
        enmig=enmig+2*fcn(phi)
       endif
      enddo
      simpson=(h/3.d0)*(enmig+extrems)
      return
      end function

c     *******************
c     *Funcio a integrar*
c     *******************
      double precision function f(phi)
      double precision phi,mu,f0,phi0
      common/PARAMETRE/phi0
      mu=(dsin(phi0/2))**2
      f0=1/dsqrt(1-mu*(dsin(phi))**2)
      f=f0
      return
      end function
