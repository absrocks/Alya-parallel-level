subroutine chofac( N, NNZ, NZ, A, INFO )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  chol_fact computes the Cholesky factorization of a real symmetric
  !  positive definite matrix A stored in skyline format.
  !
  !  The factorization has the form
  !     A = LL^T
  !  where L is a matrix stored in skyline format.
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
  !  NZ      (input) INTEGER array, dimension (N+1)
  !          The non-zeros of each row/column.
  !
  !  A       (input/output) DOUBLE PRECISION array, dimension (NNZ)
  !          On entry, the symmetri! matrix A stored in skyline format.  
  !
  !          On exit, if INFO = 0, the factor C from the Cholesky
  !          factorization A = C*C' stored in skyline format.
  !
  !  INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order k is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !  =====================================================================
  !
  !
  !     .. Local Scalars ..
  use def_kintyp, only :  ip,rp
  implicit none
  integer(ip)          :: INFO, N, NNZ
  integer(ip)          :: NZ( N+1 )
  real(rp)             :: A( NNZ )
  INTEGER(ip)          :: I, J,  K
  INTEGER(ip)          :: I0, J0
  INTEGER(ip)          :: IPOS,IDIF,IROW,JCOL
  real(rp)             :: TEMP, RDIAG
  !     ..
  !     .. Intrinsi! Functions ..
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  if( N<0 ) then
     write(*,*)'N <0' 
     INFO = -1
  else if( NNZ<0 ) then
     write(*,*)'NNZ <0' 
     INFO = -2
  end if
  if( INFO.NE.0 ) then
     return
  end if
  !
  !     Quick return if possible.
  !
  if( N==0.or.NNZ==0 )  return
  !
  !     Compute the Cholesky factorization A = LL^T.
  !
  !
  !     Initialize pointer
  !
  IPOS = 2
  !
  !     Compute first pivot
  !
  if(A(1)<0.0_rp)then
     write(*,*) 'first pivot negative'
     INFO=1
     return
  endif
  ! 
  A(1) = SQRT ( A(1))
  !    
  do  I = 2, N
     !
     !     -----SET TO 0 DIAGONAL TERM
     !
     RDIAG = 0.0_rp
     !
     !     -----COLUMN NUMBER OF THE FIRST NON ZERO IN ROW I
     !
     JCOL=I-(NZ(I+1)-NZ(I))+1
     !     
     !     -----Compute elements JCOL:(I-1) of line I.
     !     
     do  J = JCOL , I-1
        !
        !     -----ROW NUMBER OF FIRST NON ZERO IN COLUMN J
        !
        IROW=J-(NZ(J+1)-NZ(J))+1
        !
        !     -----Check for non zero bounds to setup the pointers for the scalar product
        !
        if(JCOL>IROW)then
           IDIF = J-JCOL
           I0 = NZ(I)
           J0 = NZ(J+1)-IDIF-1
        else
           IDIF = J-IROW
           I0 = IPOS-IDIF
           J0 = NZ(J)
        endif
        !
        !     -----compute scalar product
        !
        TEMP = 0.0_rp
        do K= 1,IDIF
           TEMP = TEMP + A( I0 )*A( J0 )
           I0 = I0 + 1
           J0 = J0 + 1
        enddo
        !
        A(IPOS)=(A(IPOS)-TEMP)/A(NZ(J+1)-1)
        !
        !     -----ACCUMULATE IN RDIAG
        !
        RDIAG=RDIAG+A(IPOS)*A(IPOS)
        !
        !     -----MOVE POINTER
        !
        IPOS=IPOS+1
        !
     enddo
     !
     !     -----COMPUTE THE DIAGONAL
     !

     A(IPOS) = A(IPOS)-RDIAG
     if(A(IPOS)<0.0_rp)then
        write(*,*)'Pivot negative in line ',I
        INFO=I
        return
     endif
     !
     A(IPOS)=SQRT(A(IPOS)) 
     IPOS=IPOS+1
     !
  enddo
  !
  return
end subroutine chofac
!
!    ------------------------------------------------------------------------------------------------------------
!

subroutine chosol( N, NNZ, NZ, NRHS, A, B, LDB, INFO )
  !     ..
  !
  !  Purpose
  !  =======
  !!  cho_solve solves a system of linear equations A*X = B with a symmetric
  !  positive definite matrix A stored in skyline format using the Cholesky
  !  factorization A = LL^T computed by chol_fact.
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
  !          Cholesky factorization A = L*L**T, as computed by chol_fact.
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
  !
  use def_kintyp, only :  ip,rp
  implicit none
  INTEGER(ip)          :: INFO, LDB, N, NNZ, NRHS
  INTEGER(ip)          :: NZ( N+1 )
  real(rp)             :: A( NNZ ), B( LDB, * )
  INTEGER(ip)          :: J, K
  INTEGER(ip)          :: K0
  INTEGER(ip)          :: JCNT, I
  real(rp)             :: TEMP
  !     ..
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  if( N<0 ) then
     write(*,*)'N <0' 
     INFO = -1
  else if( NNZ<0 ) then
     write(*,*)'NNZ <0' 
     INFO = -2
  else if( NRHS<0 ) then
     write(*,*)'NRHS <0' 
     INFO = -4
  else if( LDB<MAX( 1_ip, N ) ) then
     write(*,*)'LDB < MAX(1,N)' 
     INFO = -7
  end if
  if( INFO.NE.0 ) then
     return
  end if
  !
  !     Quick return if possible.
  !
  if( N==0 .OR. NNZ==0 .OR. NRHS==0 )  return
  !
  !     Loop over NRHS
  !
  do  I = 1, NRHS
     !     
     !     Forward substitution.
     !    
     B(1,I)=B(1,I)/A(1)
     JCNT = 2 
     !
     do  J = 2, N
        !     
        K0 = J - (NZ(J+1)-JCNT) + 1
        !     
        TEMP = B(J,I)
        do  K = K0, J-1
           !     
           TEMP = TEMP - A(JCNT)*B(K,I)
           JCNT = JCNT + 1
           !     
        enddo
        !     
        B(J,I) = TEMP /A(JCNT)
        JCNT=JCNT+1 
        !     
     enddo
     !     
     !     Backward substitution.
     !     
     do  J = N, 1, -1
        !     
        JCNT=JCNT-1
        TEMP = B(J,I) / A(JCNT)
        B(J,I) = TEMP
        !     
        K0 = J - (NZ(J+1)-NZ(J)) + 1
        !
        do  K =  J-1, K0, -1
           !     
           JCNT = JCNT - 1
           B(K,I) = B(K,I) - A(JCNT)*TEMP
           !     
        enddo
        !     
     enddo
     !     
  enddo
  !
  return
end subroutine chosol
!
!     ---------------------------------------------------------------------------------------------------------
!
subroutine lufact( N, NNZ, NZ, A, IDIAG, INFO )
  !
  !
  !     .. Scalar Arguments ..
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  LU_fact computes the LU Doolittle factorization of a real
  !  matrix A stored in nonsymmetri! skyline format.
  !
  !  The factorization has the form
  !     A = LU
  !  where L has 1 in its diagonal.
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
  !  NZ      (input) INTEGER array, dimension (N+1)
  !          The pointer to A.
  ! 
  !  IDIAG    The position of the diagonal in each row/column 
  !
  !  A       (input/output) DOUBLE PRECISION array, dimension (NNZ)
  !          On entry, the real matrix A stored in nonsymmetri! skyline format.  
  !
  !          On exit, if INFO = 0, the factor LU from the LU
  !          factorization A = LU stored in non symmetri! skyline format.
  !
  !  INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order k is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !  =====================================================================
  !
  !
  use def_kintyp, only :  ip,rp
  implicit none
  INTEGER(ip)          :: INFO, N, NNZ
  INTEGER(ip)          :: NZ( N+1 ), IDIAG( N )
  real(rp)             :: A( NNZ )
  INTEGER(ip)          :: I, J, K
  INTEGER(ip)          :: I0, J0
  INTEGER(ip)          :: IDIF, IROW, IIDIAG, JCOL, IPOS
  real(rp)             :: TEMP,TOL
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  TOL  = 1.0e-14_rp
  INFO = 0
  if( N<0 ) then
     INFO = -1
  else if( NNZ<0 ) then
     INFO = -2
  end if
  if( INFO.NE.0 ) then
     return
  end if
  !
  !     Quick return if possible.
  !
  if( N==0.OR.NNZ==0 )  return
  !
  !     Compute the LU factorization A = LU.
  !     Initially L(1,1)=1 and U(1,1)=A(1,1) 
  !     
  !
  !     Compute first pivot
  !
  if(abs(A(1))<TOL)then
     write(*,*) 'first pivot null'
     INFO=1
     return
  endif
  !
  !     Set IPOS
  !
  IPOS=2
  !     
  !     Outer loop on the order of the matrix
  !
  do I = 2, N
     !     
     !     Compute the lower part L(I,J) 
     !
     !     -----COLUMN NUMBER OF THE FIRST NON ZERO IN ROW I
     !
     JCOL=I-(IDIAG(I)-NZ(I))
     !
     !     -----Fill row I
     !   
     do J=JCOL,I-1
        !
        !     -----ROW NUMBER OF FIRST NON ZERO IN COLUMN J
        !
        IROW=J-(NZ(J+1)-1-IDIAG(J))    
        !
        !     -----Check for non zero bounds to setup the pointers for the scalar product
        !
        if(JCOL>IROW)then
           IDIF = J-JCOL
           I0 = NZ(I)
           J0 = NZ(J+1)-IDIF
        else
           IDIF = J-IROW
           I0 = IPOS-IDIF
           J0 = IDIAG(J)+1 
        endif
        !
        !     -----compute scalar product
        !
        TEMP = 0.0_rp  
        do K= 1,IDIF
           TEMP = TEMP + A( I0 )*A( J0 )
           I0 = I0 + 1
           J0 = J0 + 1
        enddo
        !
        A(IPOS)=(A(IPOS)-TEMP)/A(IDIAG(J))
        IPOS=IPOS+1 
        !
     enddo
     !
     !     -----JUMP THE DIAGONAL
     !
     IPOS=IPOS+1
     !
     !     -----COMPUTE THE UPPER PART OF THE FACTORIZATION
     !
     !
     !     -----ROW NUMBER OF THE FIRST NON ZERO IN COLUMN I
     !
     IROW=I-(NZ(I+1)-1-IDIAG(I))
     !
     !     -----FILL COLUMN I
     !
     do  J = IROW,I-1
        !
        !     -----COLUMN NUMBER OF FIRST NON ZERO IN ROW J
        !
        JCOL=J-(IDIAG(J)-NZ(J))    
        !
        !     -----check for non zero bounds to setup the pointers
        !
        if(JCOL>IROW)then
           IDIF = J - JCOL
           I0 = NZ(J)
           J0 = IPOS-IDIF
        else
           IDIF = J-IROW
           I0 = IDIAG(J)-IDIF
           J0 = IDIAG(I)+1
        endif
        !
        !     -----compute scalar product
        !
        TEMP = 0.0_rp  
        do  K= 1,IDIF
           TEMP = TEMP + A( I0 )*A( J0 )
           I0 = I0 + 1
           J0 = J0 + 1
        enddo
        !
        A(IPOS) = A(IPOS)-TEMP
        IPOS = IPOS + 1
        !
     enddo
     !
     !
     !     -----DIAGONAL TERM
     !
     !
     !     -----COLUMN NUMBER OF FIRST NON ZERO IN ROW I
     !
     JCOL=I-(IDIAG(I)-NZ(I))    
     !
     !     -----check for non zero bounds to setup the pointers
     !
     if(JCOL>IROW)then
        IDIF = I - JCOL
        I0 = NZ(I)
        J0 = IPOS-IDIF
     else
        IDIF = I-IROW
        I0 = IDIAG(I)-IDIF
        J0 = IDIAG(I)+1
     endif
     !
     !     -----compute scalar product
     !
     TEMP = 0.0_rp 
     do  K= 1,IDIF
        TEMP = TEMP + A( I0 )*A( J0 )
        I0 = I0 + 1
        J0 = J0 + 1
     enddo
     !
     IIDIAG=IDIAG(I)
     A(IIDIAG) = A(IIDIAG)-TEMP
     !
     !     Compute first pivot
     !
     if(abs(A(IIDIAG))<TOL)then
        write(*,*) 'pivot null line',I
        INFO=I
        return
     endif
     !
  enddo
  !     
  return
end subroutine lufact
!
!     -------------------------------------------------------------------------------------------------------------
!

subroutine lusolv( N, NNZ, NZ, NRHS, A, B, LDB, IDIAG, INFO )
  !
  !     .. Scalar Arguments ..
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  LU_solve solves a system of linear equations A*X = B with a real non  symmetric
  !  matrix A stored in skyline format using the LU
  !  factorization A = LU computed by LU_fact.
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
  !          LU factorization A = LU as computed by LU_fact. L contains 1
  !          on its diagonal. The skyline format stores the rows until the diagonal
  !          and then the columns 
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
  use def_kintyp, only :  ip,rp
  implicit none
  INTEGER(ip)          :: INFO, LDB, N, NNZ, NRHS
  INTEGER(ip)          :: NZ( N+1 ), IDIAG ( N )
  real(rp)             :: A( NNZ ), B( LDB, NRHS )
  INTEGER(ip)          :: J, K
  INTEGER(ip)          :: K0
  INTEGER(ip)          :: JCNT, IIDIAG, I
  real(rp)             :: TEMP
  !     ..
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  if( N < 0 ) then
     INFO = -1
  else if( NNZ < 0 ) then
     INFO = -2
  else if( NRHS < 0 ) then
     INFO = -4
  else if( LDB < MAX( 1_ip, N ) ) then
     INFO = -7
  end if
  if( INFO /= 0 ) then
     return
  end if
  !
  !     Quick return if possible.
  !
  if( N==0 .OR. NNZ==0 .OR. NRHS==0 )   return
  !
  !     Loop over NRHS
  !
  do  I = 1, NRHS
     !     
     !     Forward substitution. L has 1 on its diagonal
     !     
     do  J = 2, N
        !
        JCNT = NZ(J)
        K0   = J - (IDIAG(J)-JCNT)  
        !     
        TEMP = B(J,I)
        do  K = K0, J-1
           !     
           TEMP = TEMP - A(JCNT)*B(K,I) 
           JCNT = JCNT + 1
           !     
        enddo
        !
        B(J,I) = TEMP
        !
     enddo
     !     
     !     Backward substitution.
     !     
     do  J = N, 1, -1
        !    
        IIDIAG = IDIAG (J) 
        TEMP   = B(J,I) / A(IIDIAG)
        B(J,I) = TEMP
        !
        K0   = J - (NZ(J+1)-1-IIDIAG) 
        JCNT = NZ(J+1)
        !
        do  K =  J-1, K0, -1
           !     
           JCNT   = JCNT - 1
           B(K,I) = B(K,I) - A(JCNT)*TEMP
           !     
        enddo
        !     
     enddo
     !     
  enddo
  !
  return
end subroutine lusolv

