SUBROUTINE DPSTRS( N, NNZ, NZ, NRHS, A, B, LDB, INFO )
  !
  !  Purpose
  !  =======
  !
  !  DPOTRS solves a system of linear equations A*X = B with a symmetric
  !  positive definite matrix A stored in skyline format using the Cholesky
  !  factorization A = C**T*C computed by DPSTRF.
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
  !  A       (input) DOUBLE PRECISION array, dimension (NNZ)
  !          The triangular factor C stored in skyline format from the 
  !          !holesky factorization A = L*L**T, as computed by DPSTRF.
  !
  !  NRHS    (input) INTEGER
  !          The number of right hand sides, i.e., the number of columns
  !          of the matrix B.  NRHS >= 0.
  !
  !  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
  !          On entry, the right hand side matrix B.
  !          On exit, the solution matrix X.
  !
  !  LDB     (input) INTEGER
  !          The leading dimension of the array B.  LDB >= max(1,N).
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !
  !  =====================================================================
  !
  !  -- Written on 2-March-2008.
  !     Fernando Mut, George Mason University.
  !
  use def_kintyp, only : ip,rp
  implicit none
  INTEGER(ip) :: INFO, LDB, N, NNZ, NRHS
  INTEGER(ip) :: NZ( N )
  real(rp)    :: A( NNZ ),B( LDB, * )
  INTEGER(ip) :: I,J, K
  INTEGER(ip) :: K0, NZJ
  INTEGER(ip) :: KCNT, JCNT
  real(rp)    :: TEMP
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( N<0 ) THEN
     INFO = -1
  ELSE IF( NNZ<0 ) THEN
     INFO = -2
  ELSE IF( NRHS<0 ) THEN
     INFO = -4
  ELSE IF( LDB<MAX( 1_ip, N ) ) THEN
     INFO = -7
  END IF
  IF( INFO/=0 ) THEN
     RETURN
  END IF
  !
  ! Quick return if possible.
  !
  IF( N==0 .OR. NNZ==0 .OR. NRHS==0 ) RETURN
  !
  !     Loop over NRHS
  !
  DO I = 1, NRHS
     !     
     ! Forward substitution.
     !     
     JCNT = 0
     !     
     DO J = 1, N
        !     
        K0 = J - NZ(J) + 1
        !     
        TEMP = B(J,I)
        DO K = K0, J-1
           !     
           JCNT = JCNT + 1
           TEMP = TEMP - A(JCNT)*B(K,I)
           !     
        end do
        !     
        JCNT = JCNT + 1
        B(J,I) = TEMP/A(JCNT)
        !     
     end do
     !     
     !     Backward substitution.
     !     
     NZJ = 0
     !     
     DO J = N, 1, -1
        !     
        JCNT = JCNT - NZJ
        TEMP = B(J,I) / A(JCNT)
        B(J,I) = TEMP
        !     
        NZJ = NZ(J)
        K0 = J - NZJ + 1
        !
        KCNT = JCNT - NZJ
        DO K = K0, J-1
           !     
           KCNT = KCNT + 1
           B(K,I) = B(K,I) - A(KCNT)*TEMP
           !     
        end do
        !     
     end do
     !     
  end do

END SUBROUTINE DPSTRS
