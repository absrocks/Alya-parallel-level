SUBROUTINE hlm_findinv(matrix,inverse,n,errorflag)

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_findinv.f90
  ! NAME
  !    hlm_findinv
  ! DESCRIPTION
  !    This routine finds the inverse of a square real matrix
  ! INPUT ARGUMENTS
  !    MATRIX ...... Input square matrix
  !    N ........... Matrix dimension
  !    ERRORFLAG ... Return error status: -1 for error, 0 for normal
  ! OUTPUT ARGUMENTS
  !    INVERSE ..... Inverted matrix
  ! USES
  ! USED BY
  !    hlm_wls
  !-----------------------------------------------------------------------

  USE def_master

  IMPLICIT NONE

  !Declarations
  INTEGER(ip), INTENT(IN)  :: n
  INTEGER(ip), INTENT(OUT) :: errorflag  
  REAL(rp), INTENT(IN)     :: matrix(n,n)   
  REAL(rp), INTENT(OUT)    :: inverse(n,n) 

  LOGICAL     :: FLAG = .TRUE.
  INTEGER(ip) :: i, j, k, l
  REAL(rp)    :: m
  REAL(rp)    :: augmatrix(n,2*n)       !Augmented matrix

  !Augment input matrix with an identity matrix
  DO i = 1, n
      DO j = 1, 2*n
          IF (j <= n ) THEN
              augmatrix(i,j) = matrix(i,j)
          ELSE IF ((i+n) == j) THEN
              augmatrix(i,j) = 1
          ELSE
              augmatrix(i,j) = 0
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
                  WRITE(*,*) 'MATRIX IS NON-INVERTIBLE 1!'
                  inverse = 0
                  errorflag = -1
                  RETURN
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
          WRITE(*,*) 'MATRIX IS NON-INVERTIBLE 2!'
          inverse = 0
          errorflag = -1
          RETURN
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
              augmatrix(i,j) = augmatrix(i,j) - augmatrix(k+1,j) * m
          END DO
      END DO
  END DO

  !Store answer
  DO i =1, n
      DO j = 1, n
          inverse(i,j) = augmatrix(i,j+n)
      END DO
  END DO
  errorflag = 0

END SUBROUTINE hlm_findinv


SUBROUTINE hlm_invcmplx(matrix,inverse,n,errorflag,needsymm)

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_findinv.f90
  ! NAME
  !    hlm_invcmplx
  ! DESCRIPTION
  !    This routine finds the inverse of a square complex matrix
  ! INPUT ARGUMENTS
  !    MATRIX ...... Input square matrix
  !    N ........... Matrix dimension
  !    ERRORFLAG ... Return error status: -1 for error, 0 for normal
  !    NEEDSYMM  ... Make the resulting matrix symmetric (often this routine returns non-symm inverses for symm matrices)
  ! OUTPUT ARGUMENTS
  !    INVERSE ..... Inverted matrix
  ! USES
  ! USED BY
  !    AINV preconditioner
  !-----------------------------------------------------------------------

  USE def_master

  IMPLICIT NONE

  !Declarations
  INTEGER(ip), INTENT(IN)  :: n,needsymm
  INTEGER(ip), INTENT(OUT) :: errorflag  
  COMPLEX(rp), INTENT(IN)  :: matrix(n,n)   
  COMPLEX(rp), INTENT(OUT) :: inverse(n,n) 
        
  LOGICAL     :: FLAG = .TRUE.
  INTEGER(ip) :: i, j, k, l
  COMPLEX(rp) :: m
  COMPLEX(rp) :: augmatrix(n,2*n)       !Augmented matrix

  !Augment input matrix with an identity matrix
  DO i = 1, n
      DO j = 1, 2_ip*n
          IF (j <= n ) THEN
              augmatrix(i,j) = matrix(i,j)
          ELSE IF ((i+n) == j) THEN
              augmatrix(i,j) = cmplx(1.0_rp,0.0_rp,kind=rp)
          ELSE
              augmatrix(i,j) = cmplx(0.0_rp,0.0_rp,kind=rp)
          ENDIF
      END DO
  END DO

 !Reduce augmented matrix to upper traingular form
  DO k =1, n-1
      IF (augmatrix(k,k) == cmplx(0.0_rp,0.0_rp,kind=rp)) THEN
          FLAG = .FALSE.
          DO i = k+1, n
              IF (augmatrix(i,k) /= cmplx(0.0_rp,0.0_rp,kind=rp)) THEN
                  DO j = 1,2_ip*n
                      augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                  END DO
                  FLAG = .TRUE.
                  EXIT
              ENDIF
              IF (FLAG .EQV. .FALSE.) THEN
                  WRITE(*,*) 'COMPLEX MATRIX IS NON-INVERTIBLE 1!'
                  inverse = 0
                  errorflag = -1
                  RETURN
              ENDIF
          END DO
      ENDIF
      DO j = k+1, n
          m = augmatrix(j,k)/augmatrix(k,k)
          DO i = k, 2_ip*n
              augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
          END DO
      END DO
  END DO

  !Test for invertibility
  DO i = 1, n
      IF (augmatrix(i,i) == cmplx(0.0_rp,0.0_rp,kind=rp)) THEN
          WRITE(*,*) 'COMPLEX MATRIX IS NON-INVERTIBLE 2!'
          inverse = 0
          errorflag = -1
          RETURN
      ENDIF
  END DO

  !Make diagonal elements as 1
  DO i = 1, n
      m = augmatrix(i,i)
      DO j = i, 2_ip*n
          augmatrix(i,j) = (augmatrix(i,j) / m)
      END DO
  END DO

  !Reduced right side half of augmented matrix to identity matrix
  DO k = n-1, 1, -1
      DO i = 1, k
          m = augmatrix(i,k+1)
          DO j = k, 2_ip*n
              augmatrix(i,j) = augmatrix(i,j) - augmatrix(k+1,j) * m
          END DO
      END DO
  END DO

  if (needsymm == 1) then
      DO i = 1, n
          inverse(i,i) = augmatrix(i,i+n)
          DO j = i+1, n
              inverse(i,j) = augmatrix(i,j+n)
              inverse(j,i) = augmatrix(i,j+n)
          END DO
      END DO
  else
      DO i = 1, n
          DO j = 1, n
              inverse(i,j) = augmatrix(i,j+n)
          END DO
      END DO
  end if
  errorflag = 0

END SUBROUTINE hlm_invcmplx

