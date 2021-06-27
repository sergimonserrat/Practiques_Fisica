c     Aquest programa calcula Van der Waals
      program vanderwaals
      implicit none
      double precision fun,x1,x2,e,xarrel
      integer n
      external fun
      
      call Bisection(-1,3,1.d-12,n,xarrel)
      write (*,*)

      end program

c     **********
c     *Biseccio*
c     **********
      subroutine Bisection(A,B,eps,fun,niter,xarrel)
      double precision A,B,C,eps,fun,xarrel,delta
      integer niter,n
      delta=2.d0*eps
      n=0
      C=(A+B)/2.d0
      do while(delta.gt.eps)
       n=n+1
       call fun(A,fu,dfu)
       call fun(B,fu,dfu)
       if ((fu(A)*fu(B)).lt.0) A=C
       if ((fu(A)*fu(B)).gt.0) B=C
       delta=dabs(fu(B)-fu(A))
      enddo
      niter=n
      xarrel=(A+B)/2.d0
      end subroutine

c     ****************
c     *Newton Raphson*
c     ****************
      subroutine NewtonRaphson(x0,eps,fun,niter,xarrel)
      double precision x0,eps,fun,xarrel,delta,x1
      integer niter,n
      delta=2.d0*eps
      n=0
      do while(delta.gt.eps)
       n=n+1
       call fun(x0,fu,dfu)
       x1=x0-(fu/dfu)
       delta=dabs(x1-x0)
       x0=x1
      enddo
      niter=n
      xarrel=x0
      end subroutine

c     ***********************************
c     *Avaluació i derivada de la funcio*
c     ***********************************
      subroutine fun(x,fu,dfu)
      double precision x,fu,dfu

      fu=(35.d0/16.d0)+(x/2.d0)-(61.d0/20.d0)*(x**2)+x**3
      dfu=(1.d0/2.d0)-(61.d0/10.d0)*x+3.d0*(x**2)
      return
      end subroutine

c     *******************
c     *Derivada centrada*
c     *******************
      subroutine derfun(ndat)
      end subroutine