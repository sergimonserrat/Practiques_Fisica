c     This program uses Euler methods to solve
c     the differential equation for a pendulum.
c     Then it computes the potential and kinetic
c     energy.
      program pendol
      implicit none
      double precision g,l,Tn,PI,phi,phi2,m,h,velphi,t,p
      double precision velphi2,phigrans,velgrans,phigrans2
      double precision velgrans2,epot,ekin,kin,pot
      double precision phipot,velkin,phipot2,velkin2
      double precision phiplus,phiplus2,velplus,velplus2
      double precision phiminus,phiminus2,velminus,velminus2
      integer N,i

c     Variables      
      PI=dacos(-1.d0)
      g=1.66d0
      l=1.05d0
      Tn=2.d0*PI/sqrt(g/l)
      N=1400
      m=0.95d0

c     Condicions inicials
      phi=0.075d0
      velphi=0.d0
      phigrans=PI-0.12d0 
      velgrans=0.d0
      phipot=PI-0.012d0
      velkin=0.1d0
      t=0.d0
      
      open(1,file='P7-1819P-resf.dat')
      write(1,*) "Petites oscilacions"
      write(1,*) "Euler simple"

c     Càlcul del pas
      h=6.d0*Tn/dble(N)
c     Mètode d'Euler
      do i=1,N
      t=h*dble(i)
      phi2=phi+h*velphi
      velphi2=velphi+h*(-(g/l)*dsin(phi))

      phigrans2=phigrans+h*velgrans
      velgrans2=velgrans+h*(-(g/l)*dsin(phigrans))

      phipot2=phipot+h*velkin
      velkin2=velkin+h*(-(g/l)*dsin(phipot))

      pot=epot(phipot,m,g,l)
      kin=ekin(velkin,m,l)
      write(1,*) t,phi2,velphi2,phigrans2,velgrans2,kin,pot+kin
      phi=phi2
      velphi=velphi2
      phigrans=phigrans2
      velgrans=velgrans2
      phipot=phipot2
      velkin=velkin2
      enddo
      write(1,*) 
      write(1,*)

c     Reset de condicions inicials
      phi=0.075d0
      velphi=0.d0
      phigrans=PI-0.12d0 
      velgrans=0.d0
      phipot=PI-0.012d0
      velkin=0.1d0
      phiplus=0.d0
      phiminus=0.d0
      velplus=2.d0*sqrt(g/l)+0.08d0
      velminus=2.d0*sqrt(g/l)-0.08d0
      t=0.d0
      
c     Predictor corrector
      write(1,*) "Predictor corrector"
      do i=1,N
      t=h*dble(i)
      phi2=phi+(h/2)*(velphi+(velphi+h*(-(g/l)*dsin(phi))))
      velphi2=velphi+(h/2)*(-(g/l)*dsin(phi)-(g/l)*dsin(phi+h*velphi))

      phigrans2=phigrans+
     +(h/2)*(velgrans+(velgrans+h*(-(g/l)*dsin(phigrans))))
      velgrans2=velgrans+
     +(h/2)*(-(g/l)*dsin(phigrans)-(g/l)*dsin(phigrans+h*velgrans))

      phipot2=phipot+(h/2)*(velkin+(velkin+h*(-(g/l)*dsin(phipot))))
      velkin2=velkin+
     +(h/2)*(-(g/l)*dsin(phipot)-(g/l)*dsin(phipot+h*velkin))

      phiplus2=phiplus+
     +(h/2)*(velplus+(velplus+h*(-(g/l)*dsin(phiplus))))
      velplus2=velplus+
     +(h/2)*(-(g/l)*dsin(phiplus)-(g/l)*dsin(phiplus+h*velplus))
      phiminus2=phiminus+
     +(h/2)*(velminus+(velminus+h*(-(g/l)*dsin(phiminus))))
      velminus2=velminus+
     +(h/2)*(-(g/l)*dsin(phiminus)-(g/l)*dsin(phiminus+h*velminus))

      pot=epot(phipot,m,g,l)
      kin=ekin(velkin,m,l)
      write(1,*) t,phi2,velphi2,phigrans2,velgrans2,kin,pot+
     +kin,phiplus2,velplus2,phiminus2,velminus2
      phi=phi2
      velphi=velphi2
      phigrans=phigrans2
      velgrans=velgrans2
      phipot=phipot2
      velkin=velkin2
      phiplus=phiplus2
      phiminus=phiminus2
      velplus=velplus2
      velminus=velminus2
      enddo 

      close(1)
      end program

c     *******************
c     *Energia potencial*
c     *******************
      double precision function epot(phi,m,g,l)
      double precision phi,m,g,l
      epot=-m*g*l*dcos(phi)
      return
      end function
c     ******************
c     *Energia cinetica*
c     ******************
      double precision function ekin(velphi,m,l)
      double precision velphi,m,l
      ekin=(m/2.d0)*(velphi**2)*(l**2)
      return
      end function
