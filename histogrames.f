c     Aquest programa fa histogrames
      program histogrames
      implicit none
      double precision ran2,fun
      external fun
      fun=(6.d0/5.d0)(1-(x**2)/2)

      call acceptrebuig(20000,xnums,0.d0,1.d0,6.d0/5.d0,fun)
      call histograma(20000,xnums,0.d0,1.d0,30,xhis,vhis,errhis,boxsize,ierr)
      open(1,file='P5-1819P-res.dat')
      do
      enddo
      close(1)
      end program

c     ******************
c     *Metode de rebuig*
c     ******************
      subroutine acceptrebuig(ndat,xnums,a,b,M,fun)
      double precision a,b,M,xnums(1:ndat),r(1:2),x,p,suma,mitjana
      double precision suma2
      integer ndat,i,j,k

      do j=1,ndat
1      do i=1,2
        r(i)=ran2(-16876996)
       enddo
       x=(b-a)*r(1)+a
       p=M*r(2)
       if (fun(x).gt.p) then
        xnums(j)=x
       else goto 1
      enddo
      suma=0
      suma2
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
      
      *********************
      *Generador gaussiana*
      *********************
      subroutine subgauss(ndat,xmu,xsigma,xgaus)
      double precision radi,PI,xmu,xsigma,xgaus(1:ndat),r(1:2),phi
      integer ndat,i,j
      PI=dacos(-1.d0)
      do j=1,ndat
       do i=1,2
        r(i)=ran2(-16876996)
       enddo
       phi=r(1)*(2.d0*PI)
       radi=sqrt(-2.d0*dlog(r(2)))
       xgaus(j)=radi*dcos(phi) !modificar per tenir sigma i mu desitjades
      enddo
      end subroutine
      
      ***********************
      *Generador exponencial*
      ***********************
      subroutine subexpo(ndat,xlambda,xexpo)
      double precision xlamba,xexpo(1:ndat)
      integer ndat,i,j

      do i=1,ndat
       r=ran2(-16876996)
       xexpo(i)=xlambda*dlog(r)
      enddo
      **********************
      *Generador histograma*
      **********************
      subroutine histograma(ndat,xdata,xa,xb,nbox,xhis,vhis,errhis,
     +boxsize,ierr)
      implicit none
      integer ndat,nbox,ierr,i,icount
      double precision xdata(1:ndat),xa,xb,xhis(nbox),errhis(nbox)
      double precision boxsize
      
      if (xa.ge.xb) then
       ierr=1
       return
      endif
      
      boxsize=(xb-xa)/nbox

      icount=0

      do i=1,nbox
       vhis(i)=0
       errhis(i)=0
      enddo
      
      do i=1,ndat
       if (xdata(i).ge.xa.and.xdata(i).lt.xb) then
        ibox=int((xdata(i)-xa)/boxsize)+1
        if (ibox.eq.nbox+1) ibox=nbox
        vhis(ibox)=vhis(ibox)+1
        icount=icount+1
       endif
      enddo

      if (icount.eq.0) then
       ierr=2
       return
      endif

      ierr=0
      write(*,*) icount,'punts acceptats de',ndat,'punts'

      do i=1,nbox
       xhis(i)=xa+boxsize/2.d0+(i-1)*boxsize
       errhis(i)=sqrt(vhis(i)/icount*(1.d0-vhis(i)/icount))/boxsize
     + /sqrt(dble(icount))
       vhis(i)=vhis(i)/icount/boxsize
      enddo
      end subroutine

      ********************
      *Generador aleatori*
      ********************
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

