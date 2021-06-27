c%%%%%Aquest programa implementa el metode de la biseccio per trobar
c%%%%%zeros de funcions continues. Garanteix que trobara al menys un
c%%%%%zero sempre que aquest no sigui simultaneament un extrem

      program biseccio
      implicit none
	  
	  double precision x,x0,x1,x2,f,e,delta,gap
	  integer n
	  
	  parameter(delta=1.d-15)
	  
c%%%%%%funcio de la que busquem els zeros%%%%%%%%	  

	  f(x)=dcos(x)-x
c      f(x)=exp(x)-x-1. 
c      f(x)=1.d0/(x-1.d0)

c%%%%%%funcio de la que busquem els zeros%%%%%%%%	

c%%%%%%ens demana per pantalla els limits inferior i superior%%%

1     write(*,*) 'entra el limit inferior'
	  read(*,*) x0
	  write(*,*) 'entra el limit superior'
	  read(*,*) x1
	  
	  if(((f(x0)*f(x1)).gt.0.))then
	   write(*,*) 'interval erroni'
	   goto 1
	  endif 
c%%%%%%ens demana per pantalla els limits inferior i superior%%%

       open(1,file='biseccio.dat',status='unknown')

c%%%%%comença la biseccio%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	  

      n=0
	  e=2.d0*delta
	  do while(e.gt.delta)
	   n=n+1
	   x2=(x0+x1)/2.
	    if((f(x2)*f(x1)).lt.0.) x0=x2
		if((f(x0)*f(x2)).lt.0.) x1=x2
       write(*,*) 'el nou interval es',x0,x1
	   e=dabs(x1-x0)/dabs(x0+x1)
	   gap=dabs(f(x1)-f(x0))
	   write(*,*) 'amb un error relatiu de',e,n
c	   write(*,*) 'amb un gap de',gap
	   write(1,*) n,x0,x1,e,gap
	  enddo

      close(1)
      stop
	  end
	  