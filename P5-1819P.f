c     Aquest programa fa coses de nombres aleatoris
      program gaussianes
      implicit none
      double precision xgaus(1:110000),xhis(1:100),vhis(1:100)
      double precision errhis(1:100),boxsize, exacte(1:100),PI
      double precision xmol(1:220),ymol(1:220),suma,sumaquadrat
      double precision deltax(1:220),deltay(1:220),d,dt
      double precision varianca(1:250),t(1:250)
      integer i,ierr,comptador,j
      PI=dacos(-1.d0)
      d=2.62d-5
      dt=2.d-2
c     Exercici histograma
      call subgauss(110000,0.d0,1.d0,xgaus)
      open(1,file='P5-1819P-res.dat')
      open(2,file='P5-1819P-res2.dat')
c     Aquest codi comentat es per comprovar que els nombres es generen bé.
c      do i=1,110000 
c       write(1,*) i,xgaus(i) 
c      enddo
      call histograma(110000,xgaus,-4.5d0,4.5d0,100,xhis,vhis,errhis,
     +boxsize,ierr)
      do i=1,100
       exacte(i)=(1/sqrt(2*PI))*dexp(-(xhis(i)**2)/2)
       write(1,*) xhis(i),vhis(i),errhis(i),exacte(i)
      enddo
      write(1,*) 
      write(1,*) 
c     Exercici molecules
c     Coordenades inicials
      do i=1,220
       xmol(i)=0.d0
       ymol(i)=0.d0
      enddo
c     Trajectories
      comptador=1
      do j=1,250
       t(j)=j*dt
       suma=0.d0
       sumaquadrat=0.d0
       do i=1,220
        xmol(i)=xmol(i)+xgaus(comptador)*d*dt
        comptador=comptador+1
        ymol(i)=ymol(i)+xgaus(comptador)*d*dt
        comptador=comptador+1
        suma=suma+xmol(i)
        sumaquadrat=sumaquadrat+(xmol(i)**2)
       enddo
      varianca(j)=(sumaquadrat/220.d0)-(suma/220.d0)**2
      write(1,*) xmol(1),ymol(1),xmol(2),ymol(2),xmol(3),ymol(3),
     +xmol(4),ymol(4),xmol(5),ymol(5)
      write(2,*) t(j),varianca(j)
      enddo
        
      
      end program

c     *********************
c     *Generador gaussiana*
c     *********************
      subroutine subgauss(ndat,xmu,xsigma,xgaus)
      double precision ran2,radi,PI,xmu,xsigma,xgaus(1:ndat),r(1:2),phi
      integer ndat,i,j,idum
      idum=-16876996
      PI=dacos(-1.d0)
      do j=1,ndat
       do i=1,2
        r(i)=ran2(idum)
       enddo
       phi=r(1)*(2.d0*PI)
       radi=sqrt(-2.d0*dlog(r(2)))
       xgaus(j)=radi*dcos(phi)
      enddo
      end subroutine

c     **********************
c     *Generador histograma*
c     **********************
      subroutine histograma(ndat,xdata,xa,xb,nbox,xhis,vhis,errhis,
     +boxsize,ierr)
      implicit none
      integer ndat,nbox,ierr,i,icount,ibox
      double precision xdata(1:ndat),xa,xb,xhis(nbox),errhis(nbox)
      double precision boxsize,vhis(nbox)
      
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

