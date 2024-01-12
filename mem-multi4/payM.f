	program payM
	
	implicit none
	integer, parameter :: nS1=2**4, nStt=2*4
	integer, parameter :: nS=2**(2*4)
	character*80 outf
	integer i,j,k, di,st1,st2, ii
	double precision eps,M(nS1,nS1),bpay(nS1), pwr,pw0,pD(4),pC(4),
     &		pI(nStt),pJ(nStt), coopR !,sbinI(nS1),sbinJ(nS1)

	read(*,*) st1, st2, di
	read(*,*) eps,outf
	open(90,file=outf,status='unknown')

	if(di.eq.1) then
	do i=st1,st2	! Elements outside of the diagonal
	  call s2bin(i,eps,pI)
	  write(*,'(I12,16F8.4)'),i!,pI
	do j=0,nS-1
	 if(i.ne.j) then
	  call s2bin(j,eps,pJ)
	  !write(*,'(A,I12,16F8.4)'),'     ',j!,pJ

	  call createM(pI,pJ,M)
	  call autovect(M,bpay)

	  pC(0+1)=bpay(8)
	  pC(1+1)=bpay(5)+bpay(6)+bpay(7)
	  pC(2+1)=bpay(2)+bpay(3)+bpay(4)
	  pC(3+1)=bpay(1)
	  pD(1+1)=bpay(5+8)+bpay(6+8)+bpay(7+8)
	  pD(2+1)=bpay(2+8)+bpay(3+8)+bpay(4+8)
	  pD(3+1)=bpay(1+8)	
	  pwr=pC(0+1)+2.d0*pC(1+1)+3.d0*pC(2+1)+4.d0*pC(3+1)
     &			+1.d0*pD(1+1)+2.d0*pD(2+1)+3.d0*pD(3+1)
	  pw0=pC(0+1)+pC(1+1)+pC(2+1)+pC(3+1)

	  write(90,'(2E22.12)') pwr,pw0
	 endif
	enddo
	enddo

	else
	do i=0,nS-1 ! Elements on the diagonal
	  call s2bin(i,eps,pI)
	  !write(*,'(I12,16F8.4)'),i,pI

	  call createM(pI,pI,M)
	  call autovect(M,bpay)

	  pC(0+1)=bpay(8)
	  pC(1+1)=bpay(5)+bpay(6)+bpay(7)
	  pC(2+1)=bpay(2)+bpay(3)+bpay(4)
	  pC(3+1)=bpay(1)
	  pD(0+1)=bpay(8+8)
	  pD(1+1)=bpay(5+8)+bpay(6+8)+bpay(7+8)
	  pD(2+1)=bpay(2+8)+bpay(3+8)+bpay(4+8)
	  pD(3+1)=bpay(1+8)	
	  pwr=pC(0+1)+2.d0*pC(1+1)+3.d0*pC(2+1)+4.d0*pC(3+1)
     &			+1.d0*pD(1+1)+2.d0*pD(2+1)+3.d0*pD(3+1)
	  pw0=pC(0+1)+pC(1+1)+pC(2+1)+pC(3+1)

c	  coopR=pI(4)*pC(0+1)+pI(3)*pC(1+1)+pI(2)*pC(2+1)+pI(1)*pC(3+1)
c     &		+pI(8)*pD(0+1)+pI(7)*pD(1+1)+pI(6)*pD(2+1)+pI(5)*pD(3+1)

c	     if(i.eq.240) then
c		print*,pI
c		print*
c		print*,bpay
c		print*
c		print*,pC
c		print*,pD
c		print*, pwr,pw0,coopR
c		stop
c	      endif
	  write(90,'(2E22.12)') pwr,pw0
	enddo
	endif

	close(90)
	stop
	end


	subroutine s2bin(Inum,eps,pI)
	implicit none
	integer, parameter :: nStt=2*4
	integer Inum, i
	double precision pI(nStt),eps, t1,t2
	pI=eps
	t1=Inum
	do i=1,nStt
	  t2=int(t1/2.d0)
	  if(t2.ne.(t1/2.d0)) pI(i)=1.d0-eps
	  t1=t2
	enddo
	return
	end


      	subroutine createM(pI,pJ,M)
	implicit none
	integer, parameter :: nS1=2**4, nStt=2*4
	double precision pI(nStt),pJ(nStt),M(nS1,nS1),pf
	integer i,j, sj,tj, si,sii,ti
	M=0.d0

	do tj=0,1
	  sj=tj*8
	do ti=0,1
	  si=ti*8
	  sii=ti*4

	pf=1.d0*tj+(1.d0-2.d0*tj)*pI(1+sii)
	M(1+si,1+sj)=pf*(pJ(1+ti)**3.d0)
	M(1+si,2+sj)=pf*(pJ(1+ti)**2.d0)*(1.d0-pJ(1+ti))
	M(1+si,3+sj)=M(1+si,2+sj)
	M(1+si,4+sj)=M(1+si,2+sj)
	M(1+si,5+sj)=pf*pJ(1+ti)*((1.d0-pJ(1+ti))**2.d0)
	M(1+si,6+sj)=M(1+si,5+sj)
	M(1+si,7+sj)=M(1+si,5+sj)
	M(1+si,8+sj)=pf*((1.d0-pJ(1+ti))**3.d0)

	pf=1.d0*tj+(1.d0-2.d0*tj)*pI(2+sii)
	M(2+si,1+sj)=pf*(pJ(2+ti)**2.d0)*pJ(5+ti)  ! XCCD -> YCCC
	M(2+si,2+sj)=pf*(pJ(2+ti)**2)*(1.d0-pJ(5+ti)) ! XCCD -> YCCD
	M(2+si,3+sj)=pf*pJ(2+ti)*(1.d0-pJ(2+ti))*pJ(5+ti) ! XCCD -> YCDC
	M(2+si,4+sj)=M(2+si,3+sj)
	M(2+si,5+sj)=pf*pJ(2+ti)*(1.d0-pJ(2+ti))*(1.d0-pJ(5+ti)) ! XCCD -> YCDD
	M(2+si,6+sj)=M(2+si,5+sj)
	M(2+si,7+sj)=pf*((1.d0-pJ(2+ti))**2.d0)*pJ(5+ti) ! XCCD -> YDDC
	M(2+si,8+sj)=pf*((1.d0-pJ(2+ti))**2.d0)*(1.d0-pJ(5+ti)) ! XCCD -> YDDD

	M(3+si,1+sj)=M(2+si,1+sj)
	M(3+si,2+sj)=M(2+si,3+sj)
	M(3+si,3+sj)=M(2+si,2+sj)
	M(3+si,4+sj)=M(2+si,3+sj)
	M(3+si,5+sj)=M(2+si,5+sj)
	M(3+si,6+sj)=M(2+si,7+sj)
	M(3+si,7+sj)=M(2+si,5+sj)
	M(3+si,8+sj)=M(2+si,8+sj)

	M(4+si,1+sj)=M(2+si,1+sj)
	M(4+si,2+sj)=M(2+si,3+sj)
	M(4+si,3+sj)=M(2+si,3+sj)
	M(4+si,4+sj)=M(2+si,2+sj)
	M(4+si,5+sj)=M(2+si,7+sj)
	M(4+si,6+sj)=M(2+si,5+sj)
	M(4+si,7+sj)=M(2+si,5+sj)
	M(4+si,8+sj)=M(2+si,8+sj)

	pf=1.d0*tj+(1.d0-2.d0*tj)*pI(3+sii)
	M(5+si,1+sj)=pf*pJ(3+ti)*pJ(6+ti)**2.d0 ! XCDD -> YCCC
	M(5+si,2+sj)=pf*pJ(3+ti)*pJ(6+ti)*(1.d0-pJ(6+ti)) ! XCDD -> YCCD
	M(5+si,3+sj)=M(5+si,2+sj)
	M(5+si,4+sj)=pf*(1.d0-pJ(3+ti))*pJ(6+ti)**2.d0 ! XCDD -> YDCC
	M(5+si,5+sj)=pf*pJ(3+ti)*(1.d0-pJ(6+ti))**2.d0 ! XCDD -> YCDD
	M(5+si,6+sj)=pf*(1.d0-pJ(3+ti))*pJ(6+ti)*(1.d0-pJ(6+ti)) ! XCDD -> YDCD
	M(5+si,7+sj)=M(5+si,6+sj)
	M(5+si,8+sj)=pf*(1.d0-pJ(3+ti))*((1.d0-pJ(6+ti))**2.d0) ! XCDD -> YDDD	

	M(6+si,1+sj)=M(5+si,1+sj)
	M(6+si,2+sj)=M(5+si,2+sj)
	M(6+si,3+sj)=M(5+si,4+sj)
	M(6+si,4+sj)=M(5+si,2+sj)
	M(6+si,5+sj)=M(5+si,6+sj)
	M(6+si,6+sj)=M(5+si,5+sj)
	M(6+si,7+sj)=M(5+si,6+sj)
	M(6+si,8+sj)=M(5+si,8+sj)

	M(7+si,1+sj)=M(5+si,1+sj)
	M(7+si,2+sj)=M(5+si,4+sj)
	M(7+si,3+sj)=M(5+si,2+sj)
	M(7+si,4+sj)=M(5+si,2+sj)
	M(7+si,5+sj)=M(5+si,6+sj)
	M(7+si,6+sj)=M(5+si,6+sj)
	M(7+si,7+sj)=M(5+si,5+sj)
	M(7+si,8+sj)=M(5+si,8+sj)

	pf=1.d0*tj+(1.d0-2.d0*tj)*pI(4+sii)
	M(8+si,1+sj)=pf*pJ(7+ti)**3
	M(8+si,2+sj)=pf*(pJ(7+ti)**2.d0)*(1.d0-pJ(7+ti))
	M(8+si,3+sj)=M(8+si,2+sj)
	M(8+si,4+sj)=M(8+si,2+sj)
	M(8+si,5+sj)=pf*pJ(7+ti)*((1.d0-pJ(7+ti))**2.d0)
	M(8+si,6+sj)=M(8+si,5+sj)
	M(8+si,7+sj)=M(8+si,5+sj)
	M(8+si,8+sj)=pf*((1.d0-pJ(7+ti))**3.d0)

	enddo
	enddo

      	return
      	end


	
	subroutine autovect(M,bpay)
	implicit none
	integer, parameter :: nS1=2**4
	double precision M(nS1,nS1),bpay(nS1),M2(nS1,nS1),Minv(nS1,nS1)
	double precision suma
	integer i,j, errorflag
	M2=transpose(M)
	do i=1,nS1
	  M2(i,i)=M(i,i)-1.d0
	  M2(nS1,i)=1.d0
	enddo
	call FINDInvt(M2,Minv,nS1,nS1,errorflag)
	if(errorflag.ne.0) stop 'NO SE HA CALCULADO INVERSA'
	suma=0.d0
	do i=1,nS1
	  bpay(i)=Minv(i,nS1)
	  suma=suma+bpay(i)
	enddo
	do i=1,nS1
	  bpay(i)=bpay(i)/suma
	enddo
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



