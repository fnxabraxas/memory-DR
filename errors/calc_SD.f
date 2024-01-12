	program calc_SD

	implicit none

	integer, parameter :: nS1=4, nS=nS1**2
	integer i,j,k, N
	double precision payM(nS,nS),b,beta,statD(nS),eps2,eps1,eps,
     &			coop(nS),CL
	character*80 outf

	write(*,*) 'Error action: '
	read(*,*) eps
	write(*,*) 'Errors remembering others and own '
	read(*,*) eps2, eps1
	write(*,*) 'b/c: '
	read(*,*) b
	print*,'Beta, N'
	read(*,*) beta, N
	print*,'Output file'
	read(*,*) outf

	payM=0.d0
	statD=0.d0

	call calcPAYM(payM,coop,nS1,eps,eps1,eps2,b)
			!do i=1,16
			!  print*,(payM(i,j),j=1,16)
			!enddo
	call SD(statD,payM,nS,N,beta)

	open(90,file=outf,status='unknown')
	CL=0.d0
	do i=1,nS
	  write(90,'(E12.4)') statD(i)
 	  CL=CL+coop(i)*statD(i)
	  !print*, coop(i),statD(i)
	enddo
	  write(90,'(E15.4)') CL
	close(90)


	stop
	end



	subroutine calcPAYM(payM,coop,nS1,eps,eps1,eps2,b)

	integer s(2), n, nS1
	double precision p(nS1,2),M(nS1,nS1),coop(nS1**2),
     &  v(nS1),Sp,Tp,Pp,payM(nS1**2,nS1**2),eps2,eps1,eps,b,c

	c=1.d0
	Rp=b-c
	Pp=0.d0
	Sp=-c
	Tp=b
	!Rp=1.d0
	!Pp=0.d0
	!Sp=-0.5d0
	!Tp=1.5d0	

	v=0.d0
	payM=0.d0
	p=0.d0
	do i=0,15
	  do j=0,15

		s(1)=i
		s(2)=j

		call s2bin(p,s,eps,eps1,eps2)

		call createM(p,M)

		call autovest(M,v,nS1)

		payM(i+1,j+1)=v(1)*Rp+v(2)*Sp+v(3)*Tp+v(4)*Pp
	
		if(i.eq.j) coop(i+1)=v(1)+v(2)

	  enddo
	enddo

      return
      end




      subroutine s2bin(pff,s,eps,eps1,eps2)
	integer s(2)
	double precision eps,eps1,eps2,pf(4),pff(4,2)
	integer it,jt
	pff=0.d0
	pf=0.d0
	do it=1,2
	  call dec2bin(pf,s(it),eps,eps1,eps2)
	  do jt=1,4
		pff(jt,it)=pf(jt)
	  enddo 
	enddo
      return
      end

	
      subroutine dec2bin(pf,dec,eps,eps1,eps2)
	integer dec
	double precision eps,eps1,eps2, p(4), pf(4)
	integer jt
	double precision tmp, tmp1
	p=0.d0+eps
	  tmp1=dec
	  do jt=1,3
	    tmp=int(tmp1/2.d0)
	    if(tmp.ne.(tmp1/2)) p(5-jt)=1.d0-eps
	    tmp1=tmp
	  enddo
	  if(tmp.eq.1) p(1)=1.d0-eps
	  pf(1)=(1.d0-eps2)*(1.d0-eps1)*p(1)
     &		+(1.d0-eps2)*eps1*p(3)
     &		+eps2*(1.d0-eps1)*p(2)
     &		+eps2*eps1*p(4)
	  pf(2)=(1.d0-eps2)*(1.d0-eps1)*p(2)
     &		+(1.d0-eps2)*eps1*p(4)
     &		+eps2*(1.d0-eps1)*p(1)
     &		+eps2*eps1*p(3)
	  pf(3)=(1.d0-eps2)*(1.d0-eps1)*p(3)
     &		+(1.d0-eps2)*eps1*p(1)
     &		+eps2*(1.d0-eps1)*p(4)
     &		+eps2*eps1*p(2)
	  pf(4)=(1.d0-eps2)*(1.d0-eps1)*p(4)
     &		+(1.d0-eps2)*eps1*p(2)
     &		+eps2*(1.d0-eps1)*p(3)
     &		+eps2*eps1*p(1)

      return
      end



      subroutine createM(p,M)

	double precision p(4,2), M(4,4)


	M(1,1)=p(1,1)*p(1,2)
	M(1,2)=p(1,1)*(1-p(1,2))
	M(1,3)=(1-p(1,1))*p(1,2)
	M(1,4)=(1-p(1,1))*(1-p(1,2))

	M(2,1)=p(2,1)*p(3,2)
	M(2,2)=p(2,1)*(1-p(3,2))
	M(2,3)=(1-p(2,1))*p(3,2)
	M(2,4)=(1-p(2,1))*(1-p(3,2))

	M(3,1)=p(3,1)*p(2,2)
	M(3,2)=p(3,1)*(1-p(2,2))
	M(3,3)=(1-p(3,1))*p(2,2)
	M(3,4)=(1-p(3,1))*(1-p(2,2))

	M(4,1)=p(4,1)*p(4,2)
	M(4,2)=p(4,1)*(1-p(4,2))
	M(4,3)=(1-p(4,1))*p(4,2)
	M(4,4)=(1-p(4,1))*(1-p(4,2))


      return
      end



      subroutine autovest(M,v,nS1)

	double precision M(nS1,nS1), v(nS1), Mt(nS1,nS1), 
     + A(nS1-1,nS1-1), b(nS1-1),invA(nS1-1,nS1-1), vt(nS1-1)
	integer errorflag, nS1

	errorflag=-99
	do it=1,4
	  do jt=1,4
		Mt(it,jt)=M(jt,it)
	  enddo
	enddo

	do it=1,4
	  Mt(it,it)=Mt(it,it)-1.d0
	enddo

	do it=1,3
	 do jt=1,3
	    A(it,jt)=Mt(it,jt)
	 enddo
	 b(it)=-Mt(it,4)
	enddo

	call FINDInv(A, invA, 3,3, errorflag)
	if(errorflag.ne.0) stop 'NO SE HA CALCULADO INVERSA'

	vt=matmul(invA,b)
	do it=1,3
	   v(it)=vt(it)/(vt(1)+vt(2)+vt(3)+1.d0)
	enddo
	v(4)=1.d0/(vt(1)+vt(2)+vt(3)+1.d0)


      return
      end


	subroutine SD(statD,payMAT,nS,N,beta)
	implicit none
	integer nS,N
	double precision beta,payMAT(nS,nS),statD(nS)
	integer i,j,is1,is2, errorflag
	double precision  rhoMAT(nS,nS), 
     +	   suma, MATT(nS,nS),MATinv(nS,nS),
     +		  pAA,pAB,pBB,pBA, pfix, MATTt(nS,nS)
	do is1=1,nS
	  do is2=1,nS
	    if(is1.ne.is2) then
		pAA=payMAT(is1,is1)
		pAB=payMAT(is1,is2)
		pBB=payMAT(is2,is2)
		pBA=payMAT(is2,is1)
		!print*,is1-1,is2-1,pAA,pAB,pBB,pBA
		rhoMAT(is2,is1)=pfix(pAA,pAB,pBB,pBA,N,beta)
c      			if((is2.eq.2).and.(is1.eq.9)) print*,
c     +				rhoMAT(is1,is2),pAA,pAB,pBB,pBA
	    endif
	  enddo
	enddo
	do i=1,nS
	  suma=0.d0
	  do j=1,nS
	   if(i.ne.j) suma=suma+rhoMAT(i,j)
	  enddo
	  rhoMAT(i,i)=-suma!/(nS-1)
	enddo
	MATT=0.d0
	do i=1,nS
	  do j=1,nS
	    MATT(i,j)=rhoMAT(j,i)
	  enddo
	enddo
c		do i=1,nS
c		write(*,'(20E10.2)') (payMAT(i,j),j=1,20)
c		enddo
c		print*,'----------------------------------'
c		do i=1,nS
c		write(*,'(16E10.2)') (rhoMAT(i,j),j=1,16)
c		enddo
c		stop	
	do j=1,nS
	  MATT(1,j)=1.d0
	enddo
		!do i=1,nS
	  	!do j=1,nS
	    	!  print*,i,j,MATT(i,j)
	  	!enddo
		!enddo
	call FINDInv(MATT, MATinv, nS,nS, errorflag)
		!do i=1,nS
	  	!do j=1,nS
	    	!  print*,i,j,MATinv(i,1)
	  	!enddo
		!enddo
	suma=0.d0
	do i=1,nS
	  statD(i)=MATinv(i,1)
	  suma=suma+statD(i)
	enddo
			!print*,suma
			!MATTt=matmul(MATT,MATinv)
		!do i=1,nS
	  	!do j=1,nS
	    	!  print*,i,j,MATTt(i,j)
	  	!enddo
		!enddo
			
	return
	end



	double precision function pfix(pAA,pAB,pBB,pBA,N,beta)
	implicit none
	integer N
	double precision pAA,pAB,pBB,pBA,beta
	integer i,j
	double precision payA(N),payB(N),prod,difpb
	do i=1,N
	  payA(i)=((i-1)*pAA+(N-i)*pAB)/(N-1)
	  payB(i)=(i*pBA+(N-i-1)*pBB)/(N-1)
	enddo
	pfix=0.d0
	do i=1,N-1
	  prod=1.d0
	  do j=1,i
	    difpb=beta*(payA(j)-payB(j))
	    prod=prod*(1.d0+dexp(-difpb))/(1.d0+dexp(difpb))
	  enddo
	  pfix=pfix+prod
	enddo
	pfix=1.d0/(1.d0+pfix)
	return
	end



     	SUBROUTINE FINDInv(matrix, inverse, n, nMAX, errorflag)
!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm           
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN) :: n, nMAX
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	double precision, INTENT(IN), DIMENSION(nMAX,nMAX) :: matrix  !Input matrix
	double precision, INTENT(OUT), DIMENSION(nMAX,nMAX) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	double precision :: m
	double precision, DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	!Test for invertibility
	!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1.d0
			Else
				augmatrix(i,j) = 0.d0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible_"
					inverse = 0.d0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible ."
			inverse = 0.d0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0

     	END SUBROUTINE FINDInv


