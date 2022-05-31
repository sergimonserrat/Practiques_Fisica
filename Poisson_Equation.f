c     This program uses the Gauss-Seidel method to solve
c     the Poisson equation for heat conduction.
      program poisson
      implicit none
      double precision h,T_nova(0:67,0:91),Lx,Ly,eps,T_new(0:67,0:91)
      double precision T_inicial(0:67,0:91),w,delta(0:6097)
      integer x_estudi,y_estudi
      integer j,n,nx,ny
      w=1.35d0
      h=0.5d0
      Lx=33.5d0
      Ly=45.5d0
      nx=int(Lx/h)
      ny=int(Ly/h)
c     Condicions de contorn
      do n=0,ny
       T_inicial(0,n)=17.d0
       T_inicial(nx,n)=25.3d0
      enddo
      do j=0,nx
       T_inicial(j,0)=0.5d0
       T_inicial(j,ny)=11.2d0
      enddo
c     Distribució inicial de temperatures
      do n=1,ny-1
       do j=1,nx-1
        T_inicial(j,n)=10.d0
c        write(*,*) T_inicial(j,n)
       enddo
      enddo
      open(1,file='P9-1819P-Gauss-Seidel.dat')
c     Gauss-Seidel
      eps=132465.d0
      x_estudi=int(7.5d0/h)
      y_estudi=int(23.5d0/h)
      do while (eps.gt.(1.d-6))
       call gauseidel(T_inicial,h,nx,ny,delta)!,T_new)
       write(1,*) T_inicial(x_estudi,y_estudi)
c       T_inicial=T_new
       eps=maxval(delta)
c       write(1,*) eps
      enddo
      write(1,*) 
      write(1,*)
      do n=0,ny
       do j=0,nx
        write(1,*) j*h,n*h,T_inicial(j,n)
       enddo
      enddo
      close(1)
      open(2,file='P9-1819P-sobrerelaxacio.dat')
c     Sobrerelaxacio
c      eps=132465.d0
c      do while (eps.gt.(10.d-6))
c       call sobrerelaxacio(T_nova,T_inicial,h,nx,ny,w,delta)
c        write(1,*) T_inicial(x_estudi,y_estudi)
c       eps=maxval(delta)
c      enddo
      
c      do n=0,ny
c       do j=0,nx
c        write(2,*) j*h,n*h,T_nova(j,n)
c       enddo
c      enddo
      close(2)
      end program
c     **************  
c     *Gauss-Seidel*
c     **************
      subroutine gauseidel(u_old,h,nx,ny,delta)!,u)
      integer nx,ny
      double precision ro(0:nx,0:ny),r2,r3,x,y,ro_2,ro_3,fogo2,fogo3
      double precision u_old(0:nx,0:ny),u(0:nx,0:ny),q,h
      double precision delta(0:nx*ny),m
      integer j,n,compta
      do n=0,ny
       do j=0,nx
       x=j*h
       y=n*h
       ro_2=10.d0
       ro_3=6.d0
       r2=dsqrt(((x-8.d0)**2)+((y-22.5d0)**2))
       r3=dsqrt(((x-22.d0)**2)+((y-10.d0)**2))
       fogo2=ro_2*dexp((-(r2-5.d0)**2)/((0.3d0)**2))
       fogo3=ro_3*dexp((-(r3-4.d0)**2)/((0.8d0)**2))
        if (((x.gt.(18.d0)).or.(x.lt.(22.d0))).and.
     +  ((y.gt.(29.d0)).or.(y.lt.(35.d0)))) then
         ro(j,n)=fogo2+fogo3+3.d0
        else
         ro(j,n)=fogo2+fogo3
        endif
       enddo
      enddo


      compta=0
      do n=1,ny-1
       do j=1,nx-1
        u(j,n)=(1.d0/4.d0)*(u_old(j+1,n)+u_old(j-1,n)+
     +  u_old(j,n+1)+
     +  u_old(j,n-1)+(h**2)*ro(j,n))
c        m=u(j,n)
        compta=compta+1
        delta(compta)=dabs(u(j,n)-u_old(j,n))
       enddo
      enddo
      u_old=u
      do n=1,ny-1
       do j=1,nx-1
c        write(*,*) u_old(j,n)
       enddo
      enddo
      return
      end subroutine

c     ****************
c     *Sobrerelaxacio*
c     ****************
c      subroutine sobrerelaxacio(nova,vella,h,nx,ny,w,delta)
c      integer nx,ny
c      double precision vella(nx,ny),h,nova(nx,ny),w,delta(0:nx*ny),ro
c      double precision delta1
c      integer j,n
c      do n=1,n-1
c       do j=1,j-1
c        delta1=(1.d0/4.d0)*(vella(j+1,n)+nova(j-1,n)+
c     +  vella(j,n+1)+nova(j,n-1)-4*vella(j,n)+
c     +  (h**2)*ro(dble(j*h),dble(n*h)))
c        nova(j,n)=vella(j,n)+w*delta1
c       enddo
c      enddo
c      return
c      end subroutine

c     ******************
c     *Terme inhomogeni*
c     ******************
c      double precision function ro(x,y)
c      double precision r2,r3,x,y,ro_2,ro_3,fogo2,fogo3
c      ro_2=10.d0
c      ro_3=6.d0
c      r2=dsqrt(((x-8.d0)**2)+((y-22.5d0)**2))
c      r3=dsqrt(((x-22.d0)**2)+((y-10.d0)**2))
c      fogo2=ro_2*dexp((-(r2-5.d0)**2)/((0.3d0)**2))
c      fogo3=ro_3*dexp((-(r3-4.d0)**2)/((0.8d0)**2))
c      if (((x.gt.(18.d0)).or.(x.lt.(22.d0))).and.
c     +((y.gt.(29.d0)).or.(y.lt.(35.d0)))) then
c        ro=fogo2+fogo3+3.d0
c      else
c        ro=fogo2+fogo3
c      endif 
      
c      return
c      end function
