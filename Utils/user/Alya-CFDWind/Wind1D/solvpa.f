      subroutine solvpa(iamat, jamat, amatr, rhsid, unkno, nequa, nrhsi,
     .                  nthre, outso, lusol)
c****************************************************************************
c
c**** This routine solves the system of equations A*x=b using Pardiso Solver
c**** n = number of equations
c**** iamat, jamat, amatr  = CSR coefficients of matrix A
c**** x = unknown
c**** b = right hand side
c**** called by
c     pardis
c****************************************************************************
!       use def_solver, only:cpu_solve
       IMPLICIT NONE

        INTEGER, INTENT(IN)  :: nequa, nthre, nrhsi, lusol
        INTEGER, INTENT(IN)  :: iamat(nequa+1) 
        INTEGER, INTENT(IN)  :: jamat(iamat(nequa+1)-1)
        REAL*8,  INTENT(IN)  :: amatr(iamat(nequa+1)-1) 
        REAL*8,  INTENT(IN)  :: rhsid(nequa)
        REAL*8, INTENT(out)  :: unkno(nequa)   
        real*8             :: cpuin,cpuen    
        logical*1            :: outso
c        integer nthreads, OMP_GET_NUM_THREADS, OMP_GET_MAX_THREADS, maxthr
        
        !PARDISO CONTROL PARAMETERS
C..     Internal solver memory pointer for 64-bit architectures
       INTEGER*8 pt(64)
C..     Internal solver memory pointer for 32-bit architectures
C        INTEGER*4 pt(64)
C..     This is OK in both cases.
C..     INTEGER*8 pt(64)
        INTEGER maxfc, mnumm, mtype, phase, error, msglvl, idum 
        INTEGER iparm(64), memor
        REAL*8 ddum

        DATA maxfc /1/, mnumm /1/
        
        call cpu_time(cpuin)
        
C
C  .. Setup Pardiso control parameters and initialize the solvers     
C     internal adress pointers. This is only necessary for the FIRST   
C     call of the PARDISO solver.                                     
C
      mtype     = 1      ! matrix type = structurally symmetric     
C      mtype     = 11     ! for unsymmetric matrix 
      
      call pardisoinit(pt, mtype, iparm)

C  .. Numbers of Processors ( value of OMP_NUM_THREADS )
      iparm(3) = 8 ! number of threads
c      NTHREADS= OMP_GET_NUM_THREADS()
c      maxthr=OMP_GET_MAX_THREADS()
      
    
C..   Reordering and Symbolic Factorization, This step also allocates
C     all memory that is necessary for the factorization 
C
      phase     = 11      ! only reordering and symbolic factorization
      msglvl    = 0       ! with statistical information
               
      CALL pardiso (pt, maxfc, mnumm, mtype, phase, nequa, amatr, iamat
     1 , jamat,idum, nrhsi, iparm, msglvl, ddum, ddum, error)
     
      
      IF (error .NE. 0) THEN
        WRITE(LUSOL,*) 'The following pardiso ERROR was detected: ',
     1 error     
        write(*,*) 'SOLVPA: ERROR DURING PARDISO SYMBOLIC 
     1 FACTORIZATION'
      END IF
      memor = iparm(16) !permanent memory symbolic factorization


C.. Factorization.
      phase     = 22  ! only factorization
                 
                 
      CALL pardiso(pt, maxfc, mnumm, mtype, phase, nequa, amatr, iamat,
     1 jamat, idum, nrhsi, iparm, msglvl, ddum, ddum, error) 

       IF (error .NE. 0) THEN
        
        WRITE(LUSOL,*) 'The following ERROR was detected in pardiso 
     1 solver: ', error
        write(*,*)  'SOLVPA: ERROR DETECTED IN PARDISO SOLVER'
        
      END IF

C.. Back substitution and iterative refinement
      phase     = 33  ! only factorization
      iparm(8)  = 1   ! max numbers of iterative refinement steps
     
      CALL pardiso (pt, maxfc, mnumm, mtype, phase, nequa, amatr, iamat,
     1 jamat, idum, nrhsi, iparm, msglvl, rhsid, unkno, error) 
     
      IF (error .NE. 0) THEN
        write(*,*) 'SOLVPA: ERROR DETECTED IN PARDISO BACK
     1 SUBSTITUTION'
        WRITE(LUSOL,*) 'The following ERROR was detected in pardiso
     1 solver: ', error
        STOP
      END IF

!      iparm(17): total memory compsumption of the solver for the factorization and solution phases.
      if(outso) write(lusol,110)  float(iparm(17)+memor )/1000 
      
C.. Termination and release of memory
      phase     = -1           ! release internal memory
      CALL pardiso (pt, maxfc, mnumm, mtype, phase, nequa, ddum, idum,
     1 idum, idum, nrhsi, iparm, msglvl, ddum, ddum, error)

  110 format(11x,'*** VOLATILE MEMORY (MBYTES): ',f12.3)
!      call cputim(cpuen)
  !    cpu_solve = cpu_solve + cpuen - cpuin
  
      END
