      program main
      implicit none
      external f
      double precision f,int,legendre,x
      integer n,i


      do n=2,20
       call GaussLegendre(f,n,int)
       write(*,*) n,int
      enddo
      stop
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      double precision function f(x)
       double precision x,a0,a1,a2,a3,a4
       a0=0.5d0
       a1=0.24d0
       a2=-1.13d0
       a3=6.35d0
       a4=5.27d0
       f=a4*x**4+a3*x**3+a2*x**2+a1*x+a0
      return
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       subroutine GaussLegendre(f,n,integral)
	  implicit none
	  double precision x,x0,x1,e,delta,xi,y,legendre,zeros,weights,f
	  double precision integral
	  integer n,j,i,m
	  dimension zeros(1:1000)
	  dimension weights(1:1000)
	  parameter(delta=1.d-14)

c%%%%%%%%escombrat ràpid de l'interval [-1,1] per trobar els zeros de forma aproximada

        m=100*n
1	  y=1.d0
	  j=0
        if(mod(m,2).eq.0) m=m+1
	  do i=m,-m,-2
	   xi=dble(i)/dble(m)
c         write(*,*) xi,n,y,legendre(xi,n)
	   if((legendre(xi,n)*y).lt.0.d0)then
		j=j+1
		zeros(j)=xi
		y=legendre(xi,n)
	   else
	    y=legendre(xi,n)
	   endif
	  enddo

        if(j.ne.n)then
c         write(*,*) j,n,'increasing number of sampled points'
         m=2*m
         goto 1
        endif

c        do i=1,n
c	   write(*,*) 'initial zeros',zeros(i)
c	  enddo

c%%%%%%%Newton-Rapson partint dels zeros aproximats per trobar els zeros amb bona precisió

        do i=1,n
	   x0=zeros(i)
         e=2.d0*delta
	   do while(e.gt.delta)
          x1=x0-legendre(x0,n)*(1.d0-x0**2)
     +    /(dble(n)*(legendre(x0,n-1)-x0*legendre(x0,n)))
	    e=dabs(x1-x0)
          x0=x1
	   enddo
         zeros(i)=x1
	  enddo

c%%%%%%%càlcul dels pesos per fer la integral

        do i=1,n
         weights(i)=2.d0*(1.d0-zeros(i)**2)
     +   /(dble(n)*legendre(zeros(i),n-1))**2
	  enddo

c%%%%%avaluació de la integral

        integral=0.d0
	  do i=1,n
	   integral=integral+f(zeros(i))*weights(i)
	  enddo

	  return
	  end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function legendre(x,n)
	  double precision x,xi,legendre,legendre0,legendre1
	  integer n,i
	   
	   if(n.eq.0)then
	    legendre=1.d0
	   else if(n.eq.1) then
	    legendre=x
	   else
	    legendre0=1.d0
          legendre1=x
            do i=2,n
             xi=dble(i)
             legendre=((2.d0*xi-1.d0)*x*legendre1-
     +       (xi-1.d0)*legendre0)/xi
             legendre0=legendre1
             legendre1=legendre
             enddo
         endif
	  return
	  end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





