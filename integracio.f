c     Ay caramba dijo el bar sinso
      program integracio
      implicit none
      double precision f0,f,x,R,trapezis,arg,b,a,enmig,extrems,h
      integer k,N
      external f_1
c      x0=0
c      y0=0
c      h=1000
       N=2**20
       R=42.325d0
       a=-R
       b=R
c      do i=0,200000000
c      x1=x0+1.d-2
c      if(mod(i,100000).eq.0) write(*,*) x1
c      x0=x1
c      enddo
c      do i=0,2000
c      y1=dble(h)*dble(i)
c      write(*,*) y1
c      y1=y0
c      enddo
      h=(b-a)/dble(N)
 
      enmig=0
      do k=1,N-1
      arg=a+(dble(k)*h)
      enmig=enmig+2*f_1(arg)
      if (k.eq.(N-1)) write(*,*) arg,f_1(arg),(h/2.d0)*enmig
      enddo
      write(*,*)
      end program
      
c     **********************
c     *Metode dels trapezis*
c     **********************
      double precision function trapezis(a,b,N,fcn)
      double precision a,b,h,extrems,enmig,arg,total
      integer N,k
      h=(b-a)/dble(N)
      
      enmig=0
      do k=1,N-1
      arg=a+(dble(k)*h)
      enmig=enmig+2*fcn(arg)
      enddo
      trapezis=enmig
      return
      end function

c     *******************
c     *Metode de Simpson*
c     *******************
      double precision function simpson(a,b,m,fcn)
      double precision a,b,h,extrems,enmig,arg,total
      integer N,k,m
      N=2**m
      h=(b-a)/dble(N)
      
      enmig=0
      do k=1,N-1
      arg=a+(dble(k)*h)
      if (mod(k,2).eq.0) then

      if (mod(k,2).eq.1) then

      arg=a+(dble(k)*h)
      enmig=enmig+2*fcn(arg)
      enddo
      trapezis=enmig
      return
      end function

c     **********
c     *Funcio 1*
c     **********
      double precision function f_1(x)
      double precision x,R,f0,f1
      R=42.325d0
      f0=-x/(R*(dsqrt(1.d0-((x/R)**2))))
      f1=dsqrt(1.d0+f0**2)
      f_1=f1
      return
      end function

