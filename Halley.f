c%%%%%This program implements the Halley method to
c%%%%%find zeros of continuous functions

      program halley
      implicit none
	  
	  double precision x,x0,x1,f,df,df2,e,delta,gap
	  integer n
	  
	  parameter(delta=1.d-14)
	  
c%%%%%%funcio de la que busquem els zeros%%%%%%%%	  

c	  f(x)=dcos(x)-x
c      df(x)=-dsin(x)-1.
c      f(x)=1.+dsin(x) 
c      df(x)=dcos(x)
c	  f(x)=datan(x)
c	  df(x)=1.d0/(1.d0+x**2)
c	  df2(x)=-2.d0*x/(1.d0+x**2)**2
      f(x)=5.d0*x*dexp(-x)-1.d0
      df(x)=5.d0*(1.d0-x)*dexp(-x)
      df2(x)=-5.d0*(2.d0-x)*dexp(-x)


c%%%%%%funcio de la que busquem els zeros%%%%%%%%	

c%%%%%%ens demana per pantalla el valor inicial%%%

1     write(*,*) 'entra el valor inicial'
	  read(*,*) x0

c%%%%%%ens demana per pantalla el valor inicial%%%

      open(1,file='Halley.dat',status='unknown')

c%%%%%comença la newton-raphson%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	  

      n=0
	  e=2.*delta
	  do while(e.gt.delta)
	   n=n+1
       x1=x0-2.d0*f(x0)*df(x0)/(2.d0*(df(x0))**2-f(x0)*df2(x0))

       write(*,*) 'la nova estimacio es',x1
	   e=dabs(x1-x0)
	   gap=dabs(f(x1)-f(x0))
	   write(*,*) 'amb un error relatiu de',e,n
c	   write(*,*) 'amb un gap de',gap
	   write(1,*) n,x0,x1,e,gap
	   x0=x1 		
	  enddo

      close(1)
      end

