SUBROUTINE dpstrf( N, NNZ, NZ, A, INFO )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DPSTRF computes the Cholesky factorization of a real symmetric
  !  positive definite matrix A stored in skyline format.
  !
  !  The factorization has the form
  !     A = C * C'
  !  where C is a matrix stored in skyline format.
  !
  !  Arguments
  !  =========
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  NNZ     (input) INTEGER
  !          The number of non-zeros of the matrix A.  NNZ >= 0
  !
  !  NZ      (input) INTEGER array, dimension (N)
  !          The non-zeros of each row/column.
  !
  !  A       (input/output) DOUBLE PRECISION array, dimension (NNZ)
  !          On entry, the symmetric matrix A stored in skyline format.  
  !
  !          On exit, if INFO = 0, the factor C from the Cholesky
  !          factorization A = C*C' stored in skyline format.
  !
  !  INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order k is not
  !               positive definite, and the factorization could not be
  !               !ompleted.
  !
  !  =====================================================================
  !
  !  -- Written on 2-March-2008.
  !     Fernando Mut, George Mason University.
  !
  use def_kintyp, only : ip,rp
  implicit none
  INTEGER(ip) :: INFO,N,NNZ
  INTEGER(ip) :: NZ(N)
  real(rp)    :: A(NNZ)
  INTEGER(ip) :: I,J,JJ,K,IJ,NIJ
  INTEGER(ip) :: I0,J0,K0
  INTEGER(ip) :: ICNT,JCNT
  real(rp)    :: AJJ,AJJINV,TEMP

  INFO = 0
  IF( N<0 ) THEN
     INFO = -1
  ELSE IF( NNZ<0 ) THEN
     INFO = -2
  END IF
  IF( INFO/=0 ) THEN
     RETURN
  END IF
  !
  ! Quick return if possible.
  !
  IF( N==0.OR.NNZ==0 ) RETURN
  !
  ! Compute the Cholesky factorization A = C*C'.
  !
  JCNT = 0
  !     
  DO J = 1, N
     !     
     ! Compute C(J,J) and test for non-positive-definiteness.
     !     
     JJ = JCNT + NZ(J)      ! diagonal
     K0 = JCNT + 1
     !     
     TEMP = 0.0_rp
     DO K = K0, JJ-1
        TEMP = TEMP + A(K)*A(K)
     end do
     !     
     AJJ = A(JJ) - TEMP
     !     
     IF ( AJJ<=0.0_rp ) THEN
        A(JJ) = AJJ
        GO TO 90
     ENDIF
     !     
     AJJ    = SQRT(AJJ)
     A(JJ)  = AJJ
     AJJINV = 1.0_rp / AJJ
     !     
     ! Compute elements J+1:N of column J.
     !     
     ICNT = JCNT + NZ(J)
     !     
     DO I = J+1, N
        !     
        IJ = I - J
        IF ( IJ>=NZ(I) ) GO TO 49
        !     
        NIJ = NZ(I) - NZ(J)
        IF ( IJ<NIJ ) THEN
           I0 = ICNT + NIJ - IJ
           J0 = JCNT
        ELSE
           I0 = ICNT
           J0 = JCNT - NIJ + IJ
        ENDIF
        !     
        TEMP = 0.0_rp
        DO K = 1, JJ-J0-1
           J0 = J0 + 1
           I0 = I0 + 1
           TEMP = TEMP + A(J0)*A(I0)
        end do
        !     
        I0 = I0 + 1
        A(I0) = (A(I0) - TEMP) * AJJINV
        !     
49      CONTINUE
        !    
        ICNT = ICNT + NZ(I)
        !     
     end do
     !     
     JCNT = JCNT + NZ(J)
     !     
  end do
  !     
  GO TO 100
  !     
90 CONTINUE
  INFO = J
  !     
100 CONTINUE

END SUBROUTINE dpstrf
