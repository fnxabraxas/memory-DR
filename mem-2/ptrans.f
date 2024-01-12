	program ptrans
c	R=1, P=0 is assumed
	implicit none
	integer, parameter :: nS1=16
	integer, parameter :: nS=2**nS1
	integer i,j,k,N,ii,jj
	double precision S,T,beta,pR,pS,pT,Wij,Wji,Wds,Wd(nS),maxdif,
     &	   A,Bij,Bji,pfixij,pfixji,prodij,prodji,difpbij,difpbji
	character*80 inpf, inpd, outf

	read(*,*) S,T
	read(*,*) N,beta
	read(*,*) inpd
	read(*,*) inpf
	read(*,*) outf
	maxdif=max(1.d0,S,T,0.d0)-min(1.d0,S,T,0.d0)
	beta=1.d0/maxdif ! no Fermi; just linear

	open(11,file=inpd,status='old')
	i=0
	do
	  read(11,*,end=50) pR,pS
	  i=i+1
	  Wd(i)=pR+pS*(S+T)
	enddo
  50	continue
	close(11)

	open(10,file=inpf,status='old')
	open(90,file=outf,status='unknown')
	i=1
	j=1
	do
	  read(10,*,end=51) pR,pS,pT
	  j=j+1
	  if (j.gt.nS) then
	    i=i+1
	    j=i+1
	    print*,i
	  endif
		print*, '     ',j
	  Wij=pR+pS*S+pT*T
	  Wji=pR+pS*T+pT*S
	  Wds=Wd(i)+Wd(j)
	  A=Wds-Wij-Wji
	  Bij=-Wds+N*(Wij-Wd(j))
	  Bji=-Wds+N*(Wji-Wd(i))

	  pfixij=0.d0
	  pfixji=0.d0
	  do ii=1,N-1
	    prodij=1.d0
	    prodji=1.d0
	    do jj=1,ii
	      difpbij=beta*(A*jj-Bij)/(N-1)
	      difpbji=beta*(A*jj-Bji)/(N-1)
	      !prodij=prodij*(1.d0+dexp(-difpbij))/(1.d0+dexp(difpbij))
	      !prodji=prodji*(1.d0+dexp(-difpbji))/(1.d0+dexp(difpbji))
	      prodij=prodij*(1.d0-difpbij)/(1.d0+difpbij)
	      prodji=prodji*(1.d0-difpbji)/(1.d0+difpbji)
	    enddo
	    pfixij=pfixij+prodij
	    pfixji=pfixji+prodji
	  enddo
	  pfixij=1.d0/(1.d0+pfixij)
	  pfixji=1.d0/(1.d0+pfixji)
	  write(90,'(2E16.6)') pfixij, pfixji
		if(i.eq.5) stop
	enddo
 51	continue
	close(10)
	close(90)
	

	stop
	end

