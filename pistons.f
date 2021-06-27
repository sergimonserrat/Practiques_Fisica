c     Aquest programa calcula coses de pistons
      program pistons
      implicit none
      double precision x(1:4),t(0:50),R,w,L,w0,radius
      integer k,j

      w(k)=w0*((dble(k)/2.d0)+13.d-1)
      w0=13.d-1
      L=250
      open(1,file='P2_1819P-res1.dat')
      do j=0,50
       call posicio(dble(j)*2.d-2,w0,x,L)
       write(1,*) j,x
      enddo
      end program 

c     *****************
c     *Calcul del radi*
c     *****************
      double precision function radius(k,L)
      double precision L,R
      integer k

      R=L-1.d0-3.d0*(dble(k)-1)
      radius=R
      return
      end function

c     **********************
c     *Calcul de la posicio*
c     **********************
      subroutine posicio(t,w0,x,L)
      double precision posi(1:4),t,w0,x(1:4),L,w,radius
      integer i
      w(i)=w0*((dble(i)/2.d0)+13.d-1)
      do i=1,4
       posi(i)=dsqrt(-((radius(i,L)**2.d0)*dsin(w(i)*t))+L**2.d0)+
     + radius(i,L)*dcos(w(i)*t)
       x(i)=posi(i)
      enddo
      return
      end subroutine
