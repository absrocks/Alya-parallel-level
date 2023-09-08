FUNCTION Handfp (signum, excnum)
  !DEC$ ATTRIBUTES C :: handfp
  USE DFLIB
  INTEGER(2) signum, excnum
  INTEGER(4) eptrs

  WRITE(*,*) 'In signal handler for SIG$FPE'
  WRITE(*,*) 'signum = ', signum
  WRITE(*,*) 'exception = ', excnum

  EPTRS = GETEXCEPTIONPTRSQQ()
  CALL TRACEBACKQQ("Application SIGFPE error!",USER_EXIT_CODE=-1,EPTR=EPTRS)

  SELECT CASE(excnum)
  CASE( FPE$INVALID )
     STOP ' Floating point exception: Invalid number'
  CASE( FPE$DENORMAL )
     STOP ' Floating point exception: Denormalized number'
  CASE( FPE$ZERODIVIDE )
     STOP ' Floating point exception: Zero divide'
  CASE( FPE$OVERFLOW )
     STOP ' Floating point exception: Overflow'
  CASE( FPE$UNDERFLOW )
     STOP ' Floating point exception: Underflow'
  CASE( FPE$INEXACT )
     STOP ' Floating point exception: Inexact precision'
  CASE DEFAULT
     STOP ' Floating point exception: Non-IEEE type'
  END SELECT
  handfp = 1

END FUNCTION Handfp
