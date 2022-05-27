c     This program uses montecarlo techniques
c     to generate random numbers according to
c     a probability distribution using rejection
c     methods.
      program montecarlo
      implicit none
      double precision sumaup,sumadown,ran2,cru(1:10000),nup,x
      double precision quadratup,quadratdown,varup,vardown
      double precision sigmaup,sigmadown,ndown,ud,L,PI,prob
      double precision I2,var,sigma,suma,quadrat,g,psi(1:3)
      double precision xnums(1:1000000)
      integer i,compta,j,idum,N
      external prob

      do i=1,3
       do j=1,i-1
        write(*,*) 'i',i,'j',j
       enddo
      enddo
c     Valors de variables
      idum=-16876996
      L=50.d0
c     Montecarlo cru
      open(1,file='P6-1819P-res.dat')
      compta=1
      do i=1,100
        sumaup=0.d0
        sumadown=0.d0
        quadratup=0.d0
        quadratdown=0.d0
        N=compta*100
        do j=1,N
         cru(j)=ran2(idum)
         x=cru(j)
         sumaup=sumaup+3.06d0*ud(x)
         sumadown=sumadown+5.11d0*ud(x)
         quadratup=quadratup+(3.06d0*ud(x))**2
         quadratdown=quadratup+(5.11d0*ud(x))**2
        enddo
        nup=sumaup/dble(N)
        ndown=sumadown/dble(N)
        varup=-(nup)**2+quadratup/dble(N)
        vardown=-(ndown)**2+quadratdown/dble(N)
        sigmaup=sqrt(varup/dble(N))
        sigmadown=sqrt(vardown/dble(N))
        write(1,*) N,nup,sigmaup,ndown,sigmadown
        compta=compta+1
      enddo
      write(1,*) 
      write(1,*)

c     Metode rebuig
      call acceptrebuig(1000000,xnums,0.d0,2*L,(1.d0/L),prob)
c     Fita 1/L ja que es el maxim de la prob(x)

c     Montecarlo amb prob(x)   
      compta=1
      do i=1,100
        suma=0.d0
        quadrat=0.d0
        N=compta*10000
        do j=1,N
         x=xnums(j)
         suma=suma+g(x)
         quadrat=quadrat+(g(x))**2
        enddo
        I2=suma/dble(N)
        var=-(I2)**2+quadrat/dble(N)
        sigma=sqrt(var/dble(N))
        write(1,*) N,I2,sigma
        compta=compta+1
      enddo
      write(1,*) 
      write(1,*) 

c     Montecarlo multidimensional
c      compta=1
c      do i=1,33
c        suma=0.d0
c        quadrat=0.d0
c        N=compta*10000
c        do j=1,N
c         do k=1,3
c          xvector(k)=xnums(j)
c         enddo
c         suma=suma+(psi(xvector)/prob(xvector(1))*prob(xvector(2)*etcetera 
c         quadrat=quadrat+(psi(xvector)/prob(xvector))**2  
c        enddo
c        I3=suma/dble(N)
c        var=-(I2)**2+quadrat/dble(N)
c        sigma=sqrt(var/dble(N))
c        write(1,*) N,I3,sigma
c        compta=compta+1
c      enddo
c      write(1,*) 
c      write(1,*) 

      end program

c     ****************
c     *Funcio up-down*
c     ****************
      double precision function ud(x)
      double precision x
      
      ud=x**(-0.2d0)*(1.d0-x)**4
      return
      end function

c     ***********************
c     *Densitat probabilitat*      
c     ***********************
      double precision function prob(x)
      double precision PI,x,L
      PI=dacos(-1.d0)
      L=50.d0
      prob=(1.d0/L)*dsin(PI*(x-2.d0*L)/(2.d0*L))**2
      return
      end function

c     *************
c     *Funcio g(x)*
c     *************
      double precision function g(x)
      double precision PI,x,L
      PI=dacos(-1.d0)
      L=50.d0
      g=dsin(PI*8.d0*(x-2.d0*L)/L)**2
      return
      end function

c     **************
c     *Funcio d'ona*
c     **************
      double precision function psi(x)
      double precision x(1:3),PI,L,t
      double precision productori,psi,sumatori
      integer m,k,j,i
      PI=dacos(-1.d0)
      L=50.d0
      t=PI/(2.d0*L)
      psi=1.d0
      productori=1.d0
      do k=1,3
       do j=1,k-1
        do i=1,3
         sumatori=0.d0
         do m=1,3
          sumatori=sumatori+dcos(t*x(m))
         enddo
         productori=productori*sumatori*dsin(t*x(i))
        enddo
        psi=psi*(dcos(t*x(j))-dcos(t*x(k))*productori
      enddo
      end function

c     ******************
c     *Metode de rebuig*
c     ******************
      subroutine acceptrebuig(ndat,xnums,a,b,M,fun)
      double precision a,b,M,xnums(1:ndat),r(1:2),x,p,suma,mitjana
      double precision suma2,ran2
      integer ndat,i,j,k,idum
      external fun
      idum=-16876996
      do j=1,ndat
1      do i=1,2
        r(i)=ran2(idum)
       enddo
       x=(b-a)*r(1)+a
       p=M*r(2)
       if (fun(x).gt.p) then
        xnums(j)=x
       else 
        goto 1
       endif
      enddo
      suma=0
      suma2=0
      do k=1,ndat
       suma=suma+xnums(k)
       suma2=suma2+(xnums(k))**2
      enddo
      mitjana=suma/ndat
      varianca=suma2/ndat
      desvest=varianca-(mitjana)**2
      return
      end subroutine

c     ********************
c     *Generador aleatori*
c     ********************
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     +   IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     +   IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        !Long period (> 2 x 10 18 ) random number generator of L'Ecuyer with Bays-Durham shuffle
        !and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
        !of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
        !alter idum between successive deviates in a sequence. RNMX should approximate the largest
        !floating value that is less than 1.
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then               !Initialize.
          idum=max(-idum,1)             !Be sure to prevent idum = 0.
          idum2=idum
          do j=NTAB+8,1,-1           !Load the shuffle table (after 8 warm-ups).
               k=idum/IQ1
               idum=IA1*(idum-k*IQ1)-k*IR1
               if (idum.lt.0) idum=idum+IM1
               if (j.le.NTAB) iv(j)=idum
          enddo 
          iy=iv(1)
      endif
      k=idum/IQ1                        !Start here when not initializing.
      idum=IA1*(idum-k*IQ1)-k*IR1       !Compute idum=mod(IA1*idum,IM1) without over-
      if (idum.lt.0) idum=idum+IM1      !flows by Schrage's method.
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2     !Compute idum2=mod(IA2*idum2,IM2) likewise.
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV                       !Will be in the range 1:NTAB.
      iy=iv(j)-idum2                    !Here idum is shuffled, idum and idum2 are com-
      iv(j)=idum                        !bined to generate output.
      if(iy.lt.1)iy=iy+IMM1
      ran2=dmin1(AM*dble(iy),RNMX)              !Because users don't expect endpoint values.
      return
      END FUNCTION
