c     Aquest programa no es el circuit de Monaco
      program montecarlo
      implicit none
      double precision ran2,f,xnums(1:1050000),PI,cru(1:150000),suma
      integer i,idum,compta,j
      external fun

c     Valor de les variables
      idum=-16876996
      PI=dacos(-1.d0)

c     Vector nul 
      do i=1,150000
       cru(j)=0
      enddo

c     Generar amb metode rebuig      
      open(1,file='rebuig.dat')
      call acceptrebuig(1050000,xnums,-PI,PI,0.4d0,fun)
      do i=1,1050000
       write(1,*) xnums(i)
      enddo
      close (1)
c      do i=1,1050000,2500
c       write(*,*) xnums(i)
c      enddo

c     Montecarlo cru   
      open(1,file='Montecarlocru.dat')
      compta=1
      do i=1,60
        suma=0
        do j=1,2500*compta
         cru(j)=ran2(idum)
         suma=suma+2*PI*sqrt(PI**2-((2*PI*cru(j))-PI)**2)
        enddo
        f=suma/(2500*compta)
        write(*,*) 2500*compta,f
        write(1,*) 2500*compta,f
        compta=compta+1
      enddo

      end program

c     ********
c     *Funcio*      
c     ********
      double precision function fun(x)
      double precision PI,x
      PI=dacos(-1.d0)
      fun=((5.d0/4.d0)*dexp(-dabs(x))*((dsin(x))**2))/(1.d0-dexp(-PI))
      return
      end function

c     ******************
c     *Metode de rebuig*
c     ******************
      subroutine acceptrebuig(ndat,xnums,a,b,M,fun)
      double precision a,b,M,xnums(1:ndat),r(1:2),x,p,suma,mitjana
      double precision suma2
      integer ndat,i,j,k,idum
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
      !com escriure en un fitxer desde subrutina
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