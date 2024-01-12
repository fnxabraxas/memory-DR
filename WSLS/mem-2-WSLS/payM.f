	program payM
	
	implicit none
	integer, parameter :: nSta=16, nS1=9
	integer, parameter :: nS=2**nS1
	character*80 outf
	integer i,j,k, di,st1,st2
	double precision eps,M(nSta,nSta), bpay(nSta), pR,pS,pT,pP,
     &		pI(nS1),pJ(nS1) !,sbinI(nS1),sbinJ(nS1)

	read(*,*) st1, st2, di
	read(*,*) eps,outf
	open(90,file=outf,status='unknown')

	if(di.eq.1) then
	do i=st1,st2	! Elements above the diagonal
	  call s2bin(i,eps,pI)
	  write(*,'(I12,9F8.4)'),i!,pI
	do j=i+1,nS-1
	  call s2bin(j,eps,pJ)
	  !write(*,'(A,I12,16F8.4)'),'     ',j!,pJ

	  call createM(pI,pJ,M)
	  call autovect(M,bpay)

	  pR=2.d0*bpay(1)+bpay(2)+bpay(3)+bpay(4)+bpay(5)
     &		+bpay(9)+bpay(13)
	  pS=bpay(2)+bpay(5)+2.d0*bpay(6)+bpay(7)+bpay(8)
     &		+bpay(10)+bpay(14)
	  pT=bpay(3)+bpay(7)+bpay(9)+bpay(10)+2.d0*bpay(11)
     &		+bpay(12)+bpay(15)
	  pP=bpay(4)+bpay(8)+bpay(12)+bpay(13)+bpay(14)
     &		+bpay(15)+2.d0*bpay(16)
	  write(90,'(3E16.6)') pR/2.d0,pS/2.d0,pT/2.d0
	enddo
	enddo

	else
	do i=0,nS-1 !46073,46073 !0,nS-1 ! Elements on the diagonal
	  call s2bin(i,eps,pI)

	  write(*,'(I12,16F8.4)'),i!,pI

	  call createM(pI,pI,M)
	  call autovect(M,bpay)

	  pR=2.d0*bpay(1)+bpay(2)+bpay(3)+bpay(4)+bpay(5)
     &		+bpay(9)+bpay(13)
	  pS=bpay(2)+bpay(5)+2.d0*bpay(6)+bpay(7)+bpay(8)
     &		+bpay(10)+bpay(14)
	  pT=bpay(3)+bpay(7)+bpay(9)+bpay(10)+2.d0*bpay(11)
     &		+bpay(12)+bpay(15)
	  pP=bpay(4)+bpay(8)+bpay(12)+bpay(13)+bpay(14)
     &		+bpay(15)+2.d0*bpay(16)
	  write(90,'(2E16.6)') pR/2.d0,pS/2.d0 !,pT/2.d0,pP/2.d0
	enddo
	endif

	close(90)
	stop
	end


	subroutine s2bin(Inum,eps,pI)
	implicit none
	integer, parameter :: nS1=9
	integer Inum, i
	double precision pI(nS1),eps, t1,t2
	pI=eps
	t1=Inum
	do i=1,nS1
	  t2=int(t1/2.d0)
	  if(t2.ne.(t1/2.d0)) pI(i)=1.d0-eps
	  t1=t2
	enddo
	return
	end


      	subroutine createM(pI,pJ,M)
	implicit none
	integer, parameter :: nSta=16, nS1=9
	double precision pI(nS1),pJ(nS1),M(nSta,nSta)
	integer i,j
	M=0.d0

	M(1,1)=pI(1)*pJ(1)
	M(1,5)=pI(1)*(1.d0-pJ(1))
	M(1,9)=(1.d0-pI(1))*pJ(1)
	M(1,13)=(1.d0-pI(1))*(1.d0-pJ(1))

	M(2,1)=pI(2)*pJ(4)
	M(2,5)=pI(2)*(1.d0-pJ(4))
	M(2,9)=(1.d0-pI(2))*pJ(4)
	M(2,13)=(1.d0-pI(2))*(1.d0-pJ(4))

	M(3,1)=pI(4)*pJ(2)
	M(3,5)=pI(4)*(1.d0-pJ(2))
	M(3,9)=(1.d0-pI(4))*pJ(2)
	M(3,13)=(1.d0-pI(4))*(1.d0-pJ(2))

	M(4,1)=pI(5)*pJ(5)
	M(4,5)=pI(5)*(1.d0-pJ(5))
	M(4,9)=(1.d0-pI(5))*pJ(5)
	M(4,13)=(1.d0-pI(5))*(1.d0-pJ(5))



	M(5,2)=pI(2)*pJ(4)
	M(5,6)=pI(2)*(1.d0-pJ(4))
	M(5,10)=(1.d0-pI(2))*pJ(4)
	M(5,14)=(1.d0-pI(2))*(1.d0-pJ(4))

	M(6,2)=pI(3)*pJ(7)
	M(6,6)=pI(3)*(1.d0-pJ(7))
	M(6,10)=(1.d0-pI(3))*pJ(7)
	M(6,14)=(1.d0-pI(3))*(1.d0-pJ(7))

	M(7,2)=pI(5)*pJ(5)
	M(7,6)=pI(5)*(1.d0-pJ(5))
	M(7,10)=(1.d0-pI(5))*pJ(5)
	M(7,14)=(1.d0-pI(5))*(1.d0-pJ(5))

	M(8,2)=pI(6)*pJ(8)
	M(8,6)=pI(6)*(1.d0-pJ(8))
	M(8,10)=(1.d0-pI(6))*pJ(8)
	M(8,14)=(1.d0-pI(6))*(1.d0-pJ(8))



	M(9,3)=pI(4)*pJ(2)
	M(9,7)=pI(4)*(1.d0-pJ(2))
	M(9,11)=(1.d0-pI(4))*pJ(2)
	M(9,15)=(1.d0-pI(4))*(1.d0-pJ(2))

	M(10,3)=pI(5)*pJ(5)
	M(10,7)=pI(5)*(1.d0-pJ(5))
	M(10,11)=(1.d0-pI(5))*pJ(5)
	M(10,15)=(1.d0-pI(5))*(1.d0-pJ(5))

	M(11,3)=pI(7)*pJ(3)
	M(11,7)=pI(7)*(1.d0-pJ(3))
	M(11,11)=(1.d0-pI(7))*pJ(3)
	M(11,15)=(1.d0-pI(7))*(1.d0-pJ(3))

	M(12,3)=pI(8)*pJ(6)
	M(12,7)=pI(8)*(1.d0-pJ(6))
	M(12,11)=(1.d0-pI(8))*pJ(6)
	M(12,15)=(1.d0-pI(8))*(1.d0-pJ(6))



	M(13,4)=pI(5)*pJ(5)
	M(13,8)=pI(5)*(1.d0-pJ(5))
	M(13,12)=(1.d0-pI(5))*pJ(5)
	M(13,16)=(1.d0-pI(5))*(1.d0-pJ(5))

	M(14,4)=pI(6)*pJ(8)
	M(14,8)=pI(6)*(1.d0-pJ(8))
	M(14,12)=(1.d0-pI(6))*pJ(8)
	M(14,16)=(1.d0-pI(6))*(1.d0-pJ(8))

	M(15,4)=pI(8)*pJ(6)
	M(15,8)=pI(8)*(1.d0-pJ(6))
	M(15,12)=(1.d0-pI(8))*pJ(6)
	M(15,16)=(1.d0-pI(8))*(1.d0-pJ(6))

	M(16,4)=pI(9)*pJ(9)
	M(16,8)=pI(9)*(1.d0-pJ(9))
	M(16,12)=(1.d0-pI(9))*pJ(9)
	M(16,16)=(1.d0-pI(9))*(1.d0-pJ(9))


      	return
      	end


	
	subroutine autovect(M,bpay)
	implicit none
	integer, parameter :: nSta=16
	double precision M(nSta,nSta),bpay(nSta),
     &			M2(nSta,nSta),Minv(nSta,nSta)
	double precision suma
	integer i,j, errorflag
	M2=transpose(M)
	do i=1,nSta
	  M2(i,i)=M(i,i)-1.d0
	  M2(nSta,i)=1.d0
	enddo
	call FINDInvt(M2,Minv,nSta,nSta,errorflag)
	if(errorflag.ne.0) stop 'NO SE HA CALCULADO INVERSA'
	suma=0.d0
	do i=1,nSta
	  bpay(i)=Minv(i,nSta)
	  suma=suma+bpay(i)
	enddo
	do i=1,nSta
	  bpay(i)=bpay(i)/suma
	enddo
		!do i=1,16
		!write(*,'(16E11.4)') (M(i,j),j=1,16)
		!enddo
		!print*, 'bpay: ',bpay
		!print*,M(1,1)
	return
	end



     	SUBROUTINE FINDInvt(matrix, inverse, n, nMAX, errorflag)
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

     	END SUBROUTINE FINDinvt



