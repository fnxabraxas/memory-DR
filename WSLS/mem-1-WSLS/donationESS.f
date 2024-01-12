	program donationESS
c	Find b/c limit for the donation game
c	inpf_diag, outf, nF, inpf_1,inpf_2,...,inpf_nF
	implicit none
	integer, parameter :: nS1=4, nS=2**nS1
	double precision, parameter :: zero=1.d-10, supb=100.d0
	character*80 inpf, outf
	integer i,j,k, nF, kF, pI(nS1)
	double precision Wd(nS,2),pCC,pCD,pDC,Wds(nS),Tv(nS,2),T,
     &		pABcccd, pABccdc, pBAcccd, pBAccdc, num, den

	Tv(:,1)=-supb
	Tv(:,2)=supb
	read(*,*) inpf
	open(10,file=inpf,status='old')
	i=0
	do 
	  i=i+1
	  read(10,*,end=50) (Wd(i,j),j=1,2)
	  Wds(i)=Wd(i,1)+Wd(i,2)
	enddo
  50	continue
	close(10)
	if ((i-1).ne.nS) stop 'Error 1'
	read(*,*) outf
	open(90,file=outf,status='unknown')
	read(*,*) nF
	i=1
	do kF=1,nF
	  read(*,*) inpf
	  open(10,file=inpf,status='old')
	  do 
	    do j=i+1,nS
		
	      read(10,*,end=51) pCC,pCD,pDC
	      if(((Tv(i,1).lt.supb).and.(Tv(i,2).gt.-supb)).or.
     &		((Tv(j,1).lt.supb).and.(Tv(j,2).gt.-supb))) then

	       pABcccd=pCC+pCD
	       pABccdc=pCC+pDC
	       pBAcccd=pABccdc
	       pBAccdc=pABcccd

	       if((Tv(i,1).lt.supb).and.(Tv(i,2).gt.-supb)) then  ! A being stable against B
	      	num=Wds(i)-pBAcccd
	      	den=Wds(i)-pBAccdc
	      	if((dabs(num).lt.zero).and.(dabs(den).lt.zero)) then ! check second condition
c		  num=pABcccd-Wds(j)
c	      	  den=pABccdc-Wds(j)
c	      	  if((dabs(num).lt.zero).and.(dabs(den).lt.zero)) then
c		    Tv(i,1)=supb
c		    Tv(i,2)=-supb
c		  else
c		    T=num/den
c		    if (den.lt.0) then  ! change < >
c		      if((T.gt.-supb).and.(T.lt.Tv(i,2))) Tv(i,2)=T
c		    else
c		      if((T.lt.supb).and.(T.gt.Tv(i,1))) Tv(i,1)=T
c		    endif
c		  endif
	      	else
		  T=num/den
		  if (den.lt.0) then  ! change < >
		    if((T.lt.Tv(i,2))) Tv(i,2)=T
		  else
		    if((T.gt.Tv(i,1))) Tv(i,1)=T
		  endif
	      	endif
	       endif

	       if((Tv(j,1).lt.supb).and.(Tv(j,2).gt.-supb)) then  ! B being stable against A
	      	num=Wds(j)-pABcccd
	      	den=Wds(j)-pABccdc
	      	if((dabs(num).lt.zero).and.(dabs(den).lt.zero)) then ! check second condition
c		  num=pBAcccd-Wds(i)
c	      	  den=pBAccdc-Wds(i)
c	      	  if((dabs(num).lt.zero).and.(dabs(den).lt.zero)) then
c		    Tv(j,1)=supb
c		    Tv(j,2)=-supb
c		  else
c		    T=num/den
c		    if (den.lt.0) then  ! change < >
c		      if((T.gt.-supb).and.(T.lt.Tv(j,2))) Tv(j,2)=T
c		    else
c		      if((T.lt.supb).and.(T.gt.Tv(j,1))) Tv(j,1)=T
c		    endif
c		  endif
	      	else
		  T=num/den
		  if (den.lt.0) then  ! change < >
		    if((T.lt.Tv(j,2))) Tv(j,2)=T
		  else
		    if((T.gt.Tv(j,1))) Tv(j,1)=T
		  endif
	      	endif
	       endif

	      endif

	    enddo  ! j
	      !call s2binI(i-1,pI)
	      !write(*,'(I8,4I2)') i-1, (pI(k),k=1,nS1)
		print*,i
	    if ((Tv(i,1).lt.supb).and.(Tv(i,2).gt.-supb)) then
	      call s2binI(i-1,pI)
	      write(*,'(I8,4I2, 2F12.4)')
     &				i-1, (pI(k),k=1,nS1), Tv(i,1),Tv(i,2)
	      if (Tv(i,1).lt.(Tv(i,2)).and.(Tv(i,2).ge.0)) then
		write(90,'(I8,4I2,3F12.4)')
     &			i-1,(pI(k),k=1,nS1),Tv(i,1),Tv(i,2),Wds(i)
	      else
		!print*, i, 'No overlap', Tv(i,1),Tv(i,2)
	      endif
	    endif
	    i=i+1
	    if(i.gt.nS) goto 51
	  enddo
  51	  continue
	  i=i-1
	  print*,'End of file: ',i
	  close(10)
	enddo ! nF
	if (i.ne.nS) stop 'Error 2'

	close(90)

	stop
	end



	subroutine s2binI(Inum,pI)
	implicit none
	integer, parameter :: ns1=4
	integer Inum, pI(nS1), i
	double precision t1,t2
	pI=0
	t1=Inum
	do i=1,nS1
	  t2=int(t1/2.d0)
	  if(t2.ne.(t1/2.d0)) pI(i)=1
	  t1=t2
	enddo
	return
	end

