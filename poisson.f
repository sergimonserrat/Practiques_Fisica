c     Aquest programa fa calor
      program poisson
      implicit none
      double precision h,T_antiga(j,n),T_nova
      integer j,n
      external ro

      h=0.5d0
      Lx=30.5d0
      Ly=18.5d0
      
      end program

c     ********
c     *Jacobi*
c     ********
      subroutine jacobi(vella,nova,h,j,n)
      double precision vella(j,n),h,nova(j,n)
      integer j,n
      do n=1,n-1
       do j=1,j-1
        nova(j,n)=(1.d0/4.d0)*(vella(j+1,n)+vella(j-1,n)+
     +  vella(j,n+1)+vella(j,n-1)+(h**2)*ro(j,n))
       enddo
      enddo
      return
      end subroutine


c     *************  
c     *Gauss-Seidel*
c     *************
      subroutine gauss-seidel(u,h,j,n)
      double precision u(j,n),h
      integer j,n
      do n=1,n-1
       do j=1,j-1
        u(j,n)=(1.d0/4.d0)*(u(j+1,n)+u(j-1,n)+u(j,n+1)+u(j,n-1)+(h**2)*ro(j,n))
       enddo
      enddo
      return
      end subroutine

c     ****************
c     *Sobrerelaxacio*
c     ****************
      subroutine sobrerelaxacio(nova,vella,h,j,n)
      double precision vella(j,n),h,nova(j,n)
      integer j,n
      do n=1,n-1
       do j=1,j-1
        nova(j,n)=(1.d0/4.d0)*(vella(j+1,n)+vella(j-1,n)+
     +  vella(j,n+1)+vella(j,n-1)+(h**2)*ro)
       enddo
      enddo
      return
      end subroutine
      
c     ******************
c     *Terme inhomogeni*
c     ******************
      double precision function ro(x,y)
      double precision r,x,y,ro_0
      ro_0=2.17d0
      r=dsqrt(((x-22.d0)**2)+((y-10.d0)**2))
      ro=ro_0*dexp((-(r-7.d0)**2)/((0.5d0)**2))
      return
      end function