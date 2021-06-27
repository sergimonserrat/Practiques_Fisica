c%%%%%%Aquest programa fa coses
       
       program fibonacci
       implicit none
       integer k,f,suma

1      write(*,*) 'introdueix el terme que vols calcular'
       read(*,*) k
       if(k.lt.3) go to 1
       if (k.gt.35) go to 1

       write(*,*) 'el terme',k,'val',f(k)
       write(*,*) 'la suma dels termes 4 a 32 es',suma(4,32)



       end program
c      *****************
c      *Serie fibonacci*
c      *****************
       integer function f(n)
       integer n,Pn,P0,P1

       P0=1
       P1=1

       do i=1,n-2
        Pn=P0+P1
        P0=P1
        P1=Pn
       enddo
       f=Pn
       return
       end function
c      ****************
c      *Suma de termes*
c      ****************
       integer function suma(N1,N2)
       integer f,N1,N2,S
       S=0
       do i=N1,N2
        S=S+f(i)
       enddo
       suma=S
       return
       end function