!------------------------------------------------------------------------
!> @addtogroup Function
!> @{
!> @name    ToolBox for space/time functions
!> @file    mod_ker_space_time_function.f90
!> @author  Guillaume Houzeaux
!> @date    22/02/2013
!> @brief   ToolBox for space/time functions
!> @details Allocate memory, parse formulas, etc.
!>          This module is based in fparser module from Roland Schmehl.
!>
!>          \verbatim
!>          !------- -------- --------- --------- --------- --------- --------- --------- -------
!>          ! Fortran 90 function parser v1.0
!>          !------- -------- --------- --------- --------- --------- --------- --------- -------
!>          !
!>          ! This public domain function parser module is intended for applications
!>          ! where a set of mathematical expressions is specified at runtime and is
!>          ! then evaluated for a large number of variable values. This is done by
!>          ! compiling the set of function strings into byte code, which is interpreted
!>          ! very efficiently for the various variable values.
!>          !
!>          ! The source code is available from:
!>          ! http://www.its.uni-karlsruhe.de/~schmehl/opensource/fparser-v1.0.tar.gz
!>          !
!>          ! Please send comments, corrections or questions to the author:
!>          ! Roland Schmehl <Roland.Schmehl@mach.uni-karlsruhe.de>
!>          !
!>          !------- -------- --------- --------- --------- --------- --------- --------- -------
!>          ! The function parser concept is based on a C++ class library written by Warp
!>          ! <warp@iki.fi> available from:
!>          ! http://www.students.tut.fi/~warp/FunctionParser/fparser.zip
!>          !------- -------- --------- --------- --------- --------- --------- --------- -------
!>          \endverbatim
!>
!>          To add a function xxx, modify the following:
!>          1. Define variable number Cxxx
!>          2. Add function name in Funcs 'xxx  '
!>          3. Code your function case (Cxxx)
!> 
!> @{
!------------------------------------------------------------------------

module mod_ker_space_time_function

  use def_kintyp, only : ip,rp,lg
  use def_kermod, only : number_space_time_function,space_time_function
  implicit none
  private
  save
  integer,     parameter  :: is = selected_int_kind(1) ! Data type of bytecode
  integer(ip), private    :: EvalErrType               ! =0: no error occured, >0: evaluation error
  integer(is), parameter  :: cImmed   =  1,        &
       &                     fneg     =  2,        &
       &                     fadd     =  3,        &
       &                     fsub     =  4,        &
       &                     fmul     =  5,        &
       &                     fdiv     =  6,        &
       &                     fpow     =  7,        &
       &                     fAbs     =  8,        &
       &                     fExp     =  9,        &
       &                     fLog10   = 10,        &
       &                     fLog     = 11,        &
       &                     fSqrt    = 12,        &
       &                     fsinh    = 13,        &
       &                     fcosh    = 14,        &
       &                     ftanh    = 15,        &
       &                     fsin     = 16,        &
       &                     fcos     = 17,        &
       &                     ftan     = 18,        &
       &                     fasin    = 19,        &
       &                     facos    = 20,        &
       &                     fatan    = 21,        &
       &                     fpi      = 22,        &
       &                     fpos     = 23,        &
       &                     VarBegin = 24
  character(len=1), dimension(fadd:fpow),  parameter :: Ops = (/ &
       &                     '+',    &
       &                     '-',    &
       &                     '*',    &
       &                     '/',    &
       &                     '^' /)
  character(len=5), dimension(fAbs:fpos),  parameter :: Funcs = (/ &
       &                     'abs  ', &
       &                     'exp  ', & 
       &                     'log10', &
       &                     'log  ', &
       &                     'sqrt ', &
       &                     'sinh ', &
       &                     'cosh ', &
       &                     'tanh ', &
       &                     'sin  ', &
       &                     'cos  ', &
       &                     'tan  ', &
       &                     'asin ', &
       &                     'acos ', & 
       &                     'atan ', & 
       &                     'pi   ', & 
       &                     'pos  ' /)
    character (len=52),                   parameter :: calpha = 'abcdefghijklmnopqrstuvwxyz'// &
         'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  type tComp
     integer(is), dimension(:),  pointer     :: ByteCode
     integer(ip)                             :: ByteCodeSize
     real(rp),    dimension(:),  pointer     :: Immed
     integer(ip)                             :: ImmedSize
     real(rp),    dimension(:),  pointer     :: Stack
     integer(ip)                             :: StackSize
     integer(ip)                             :: StackPtr
  end type tComp
  type (tComp),   dimension(:),  pointer     :: Comp              ! Bytecode
  integer(ip),    dimension(:),  allocatable :: ipos              ! Associates function strings 
  !
  ! Variable names and values
  !
  integer(ip),    parameter                  :: numbervariables=4
  character(10)                              :: variablenames(numbervariables)      
  real(rp)                                   :: variablesvalues(numbervariables)
  integer(ip),    pointer                    :: function_position(:)
  integer(ip)                                :: ker_n_functions 
  !
  ! Interface for space/time functions evaluation
  !
  interface ker_space_time_function
     module procedure ker_space_time_function_scalar,   &
          &           ker_space_time_function_scalar_1, &
          &           ker_space_time_function_vector,   &
          &           ker_space_time_function_vector_1
  end interface

  public :: ker_n_functions 
  public :: ker_init_space_time_function
  public :: ker_space_time_function
  public :: space_time_function_number

contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Initialize space/time functions
  !> @details Allocate memory for FUNCTION_NUMBER functions
  !
  !----------------------------------------------------------------------

  subroutine ker_init_space_time_function()
    integer(ip) :: ifunc
    integer(ip) :: idime
    integer(ip) :: function_number
    !
    ! Variables
    !
    variablenames(1) = 'x'
    variablenames(2) = 'y'
    variablenames(3) = 'z'
    variablenames(4) = 't'
    !
    ! Allocate 
    !
    if( number_space_time_function > 0 ) &
         allocate(function_position(number_space_time_function))
    !
    ! FUNCTION_NUMBER= number of functions
    !
    function_number = 0
    do ifunc = 1,number_space_time_function
       function_position(ifunc) = function_number
       function_number = function_number + size(space_time_function(ifunc) % expression)
    end do
    !
    ! Initialize function parser for function_numberfunctions
    !
    if( function_number > 0 ) call initf(function_number) 
    !
    ! Parse and bytecompile ifunc-th function string
    !
    function_number = 0 
    do ifunc = 1,number_space_time_function
       do idime = 1,size(space_time_function(ifunc) % expression)
          function_number = function_number + 1
          call parsef(function_number,space_time_function(ifunc) % expression(idime), variablenames)
       end do
    end do

   ker_n_functions = function_number  

  end subroutine ker_init_space_time_function

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Get a function number
  !> @details Get a space and time function number given the function
  !>          name
  !
  !----------------------------------------------------------------------

  function space_time_function_number(wfname)
    integer(ip)                :: space_time_function_number
    character(*),  intent(in)  :: wfname
    integer(ip)                :: ifunc

    space_time_function_number = 0
    do ifunc = 1,number_space_time_function
       if( trim(wfname) == trim(space_time_function(ifunc) % name) ) then
          space_time_function_number = ifunc
       end if
    end do
    if( space_time_function_number == 0 ) &
         call runend('SPACE TIME FUNCTION '//trim(wfname)//' DOES NOT EXIST')

  end function space_time_function_number

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Evaluate a space/time function
  !> @details Evaluate multidimensional arrays using space/time functions
  !
  !----------------------------------------------------------------------

  subroutine ker_space_time_function_vector(ifunc,x,y,z,t,xvalu)
    integer(ip),  intent(in)  :: ifunc
    real(rp),     intent(in)  :: x
    real(rp),     intent(in)  :: y
    real(rp),     intent(in)  :: z
    real(rp),     intent(in)  :: t
    real(rp),     intent(out) :: xvalu(:)
    integer(ip)               :: idime,nsize,nexpr
    integer(ip)               :: current_function

    variablesvalues(1) = x
    variablesvalues(2) = y
    variablesvalues(3) = z
    variablesvalues(4) = t
    !
    ! Look for function position 
    !
    current_function = function_position(ifunc)
    nsize = size(xvalu)
    nexpr = size(space_time_function(ifunc) % expression)
   
    if( nexpr == 1 ) then
       !
       ! Only one function for all degrees of freedom
       !
       current_function = current_function + 1
       do idime = 1,nsize
          xvalu(idime) = evalf(current_function,variablesvalues)
          if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()
       end do
    else if( nexpr == nsize ) then
       !
       ! One function per degree of freedom
       !
       do idime = 1,nsize
          current_function = current_function + 1
          xvalu(idime) = evalf(current_function,variablesvalues)
          if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()
       end do
    else
       !
       ! Wrong combination
       !
       write(*,*) 'Error space/time function dimension'
    end if

  end subroutine ker_space_time_function_vector

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Evaluate a space/time function
  !> @details Evaluate multidimensional arrays using space/time functions
  !
  !----------------------------------------------------------------------

  subroutine ker_space_time_function_vector_1(ifunc,x,t,xvalu)
    integer(ip),  intent(in)  :: ifunc
    real(rp),     intent(in)  :: x(:)
    real(rp),     intent(in)  :: t
    real(rp),     intent(out) :: xvalu(:)
    integer(ip)               :: idime,nsize,nexpr
    integer(ip)               :: current_function,ndime,ndim2,ndim3

    ndime = size(x)
    ndim2 = min(2_ip,ndime)
    ndim3 = min(3_ip,ndime)
    variablesvalues(1) = x(1)
    variablesvalues(2) = x(ndim2)
    variablesvalues(3) = x(ndim3)
    variablesvalues(4) = t
    !
    ! Look for function position 
    !
    current_function = function_position(ifunc)
    nsize = size(xvalu)
    nexpr = size(space_time_function(ifunc) % expression)
   
    if( nexpr == 1 ) then
       !
       ! Only one function for all degrees of freedom
       !
       current_function = current_function + 1
       do idime = 1,nsize
          xvalu(idime) = evalf(current_function,variablesvalues)
          if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()
       end do
    else if( nexpr == nsize ) then
       !
       ! One function per degree of freedom
       !
       do idime = 1,nsize
          current_function = current_function + 1
          xvalu(idime) = evalf(current_function,variablesvalues)
          if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()
       end do
    else
       !
       ! Wrong combination
       !
       write(*,*) 'Error space/time function dimension'
    end if

  end subroutine ker_space_time_function_vector_1

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Evaluate a space/time function
  !> @details Evaluate a scalar using space/time functions
  !
  !----------------------------------------------------------------------

  subroutine ker_space_time_function_scalar(ifunc,x,y,z,t,xvalu)
    integer(ip),  intent(in)  :: ifunc
    real(rp),     intent(in)  :: x
    real(rp),     intent(in)  :: y
    real(rp),     intent(in)  :: z
    real(rp),     intent(in)  :: t
    real(rp),     intent(out) :: xvalu
    integer(ip)               :: current_function

    variablesvalues(1) = x
    variablesvalues(2) = y
    variablesvalues(3) = z
    variablesvalues(4) = t

    current_function = function_position(ifunc) + 1
    xvalu = evalf(current_function,variablesvalues)
    if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()

  end subroutine ker_space_time_function_scalar

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Evaluate a space/time function
  !> @details Evaluate a scalar using space/time functions
  !
  !----------------------------------------------------------------------

  subroutine ker_space_time_function_scalar_1(ifunc,x,t,xvalu)
    integer(ip),  intent(in)  :: ifunc
    real(rp),     intent(in)  :: x(:)
    real(rp),     intent(in)  :: t
    real(rp),     intent(out) :: xvalu
    integer(ip)               :: current_function,ndime,ndim2,ndim3

    ndime = size(x)
    ndim2 = min(2_ip,ndime)
    ndim3 = min(3_ip,ndime)
    variablesvalues(1) = x(1)
    variablesvalues(2) = x(ndim2)
    variablesvalues(3) = x(ndim3)
    variablesvalues(4) = t

    current_function = function_position(ifunc) + 1
    xvalu = evalf(current_function,variablesvalues)
    if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()

  end subroutine ker_space_time_function_scalar_1

  subroutine initf (n)

    !--------------------------------------------------------------------
    !
    ! Initialize function parser for n functions
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip), intent(in) :: n                                 ! Number of functions
    integer(ip)             :: i
    
    allocate (Comp(n))
    do i=1,n
       nullify (Comp(i)%ByteCode,Comp(i)%Immed,Comp(i)%Stack)
    end do
  end subroutine initf
 
  subroutine parsef (i, FuncStr, Var)

    !--------------------------------------------------------------------
    !
    ! Parse ith function string FuncStr and compile it into bytecode
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),                     intent(in) :: i         ! Function identifier
    character (len=*),               intent(in) :: FuncStr   ! Function string
    character (len=*), dimension(:), intent(in) :: Var       ! Array with variable names
    character (len=len(FuncStr))                :: Func      ! Function string, local use

    if (i < 1 .or. i > size(Comp)) then
       write(*,*) '*** Parser error: Function number ',i,' out of range'
       stop
    end if
    allocate (ipos(len_trim(FuncStr)))                       ! Char. positions in orig. string
    Func = FuncStr                                           ! Local copy of function string
    call Replace ('**','^ ',Func)                            ! Exponent into 1-Char. format
    call RemoveSpaces (Func)                                 ! Condense function string
    call CheckSyntax (Func,FuncStr,Var)
    deallocate (ipos)
    call Compile (i,Func,Var)                                ! Compile into bytecode

  end subroutine parsef
  
  function evalf (i, Val) result (res)

    !--------------------------------------------------------------------
    !
    ! Evaluate bytecode of ith function for the values passed in array Val(:)
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),            intent(in) :: i                  ! Function identifier
    real(rp), dimension(:), intent(in) :: Val                ! Variable values
    real(rp)                           :: res                ! Result
    integer(ip)                        :: IK                 ! Instruction pointer
    integer(ip)                        :: DP                 ! Data pointer
    integer(ip)                        :: SP                 ! Stack pointer
    real(rp),               parameter  :: zero = 0.0_rp
    
    DP = 1
    SP = 0

    do IK = 1,Comp(i) % ByteCodeSize
       select case (Comp(i) % ByteCode(ik))

       case (cImmed)
          SP=SP+1; Comp(i)%Stack(SP) = Comp(i)%Immed(DP); DP=DP+1

       case   (fneg)
          !
          ! -x
          !
          Comp(i) % Stack(SP)   = -Comp(i) % Stack(SP)

       case   (fadd)
          !
          ! x + y
          !
          Comp(i) % Stack(SP-1) =  Comp(i) % Stack(SP-1) + Comp(i) % Stack(SP); SP=SP-1

       case   (fsub)
          !
          ! x - y
          !
          Comp(i) % Stack(SP-1) =  Comp(i) % Stack(SP-1) - Comp(i) % Stack(SP); SP=SP-1

       case   (fmul)
          !
          ! x * y
          !
          Comp(i) % Stack(SP-1) =  Comp(i) % Stack(SP-1) * Comp(i) % Stack(SP); SP=SP-1

       case   (fdiv)
          !
          ! x * y
          !
          if (Comp(i) % Stack(SP)==0.0_rp) then
             EvalErrType=1
             res=zero
             return
          end if
          Comp(i) % Stack(SP-1)  = Comp(i) % Stack(SP-1)/Comp(i) % Stack(SP); SP=SP-1

       case   (fpow)
          !
          ! x^y
          !
          Comp(i) % Stack(SP-1) = Comp(i) % Stack(SP-1)**Comp(i) % Stack(SP); SP=SP-1

       case   (fAbs)
          !
          ! |x|
          !
          Comp(i) % Stack(SP)   = abs(Comp(i) % Stack(SP))

       case   (fExp)
          !
          ! exp(x)
          !
          Comp(i) % Stack(SP)   = exp(Comp(i) % Stack(SP))

       case (fLog10)
          !
          ! log10(x)
          !
          if (Comp(i) % Stack(SP) <= 0.0_rp) then; EvalErrType=3; res=zero; return; endif
          Comp(i) % Stack(SP)   = log10(Comp(i) % Stack(SP))

       case   (fLog)
          !
          ! log(x)
          !
          if (Comp(i) % Stack(SP) <= 0.0_rp) then; EvalErrType=3; res=zero; return; endif
          Comp(i) % Stack(SP)   = log(Comp(i) % Stack(SP))

       case  (fSqrt)
          !
          ! sqrt(x)
          !
          if (Comp(i) % Stack(SP) <  0.0_rp) then; EvalErrType=3; res=zero; return; endif
          Comp(i) % Stack(SP)   = sqrt(Comp(i) % Stack(SP)) 

       case  (fsinh)
          !
          ! sinh(x)
          !
          Comp(i) % Stack(SP)   = sinh( Comp(i) % Stack(SP) )

       case  (fcosh)
          !
          ! cosh(x)
          !
          Comp(i) % Stack(SP)   = cosh( Comp(i) % Stack(SP) )

       case  (ftanh)
          !
          ! tanh(x)
          !
          Comp(i) % Stack(SP)   = tanh( Comp(i) % Stack(SP) )

       case   (fsin)
          !
          ! sin(x)
          !
          Comp(i) % Stack(SP)   = sin(  Comp(i) % Stack(SP) )

       case   (fcos)
          !
          ! cos(x)
          !
          Comp(i) % Stack(SP)   = cos(  Comp(i) % Stack(SP) )

       case   (ftan)
          !
          ! tan(x)
          !
         Comp(i) % Stack(SP)   = tan(  Comp(i) % Stack(SP) )

       case  (fasin)
          !
          ! asin(x)
          !
          if ((Comp(i)%Stack(SP)<-1.0_rp).or.(Comp(i)%Stack(SP)>1.0_rp)) then
          EvalErrType=4; res=zero; return; endif
          Comp(i)%Stack(SP)   = asin(Comp(i)%Stack(SP))

       case  (facos)
          !
          ! acos(x)
          !
          if ((Comp(i)%Stack(SP)<-10._rp).or.(Comp(i)%Stack(SP)>1.0_rp)) then
          EvalErrType=4; res=zero; return; endif
          Comp(i)%Stack(SP)   = acos(Comp(i)%Stack(SP))

       case  (fatan)
          !
          ! atan(x)
          !
          Comp(i)%Stack(SP)   = atan(Comp(i)%Stack(SP))

       case  (fpi)
          !
          ! pi(x) = pi * x
          !
          Comp(i)%Stack(SP)   = Comp(i)%Stack(SP) * 3.141592653589793238462643383279502884197_rp

       case  (fpos)
          !
          ! pos(x) = max(x,0), positive
          !
          Comp(i)%Stack(SP)   = max(Comp(i)%Stack(SP),epsilon(1.0_rp))

       case  DEFAULT
          SP=SP+1; Comp(i)%Stack(SP)=Val(int(Comp(i)%ByteCode(ik),ip)-int(VarBegin,ip)+1)

       end select 
    end do
    EvalErrType = 0
    res = Comp(i)%Stack(1)
  end function evalf
 
  subroutine CheckSyntax (Func,FuncStr,Var)

    !--------------------------------------------------------------------
    !
    ! Check syntax of function string,  returns 0 if syntax is ok
    !
    !--------------------------------------------------------------------

    implicit none
    character (len=*),               intent(in) :: Func      ! Function string without spaces
    character (len=*),               intent(in) :: FuncStr   ! Original function string
    character (len=*), dimension(:), intent(in) :: Var       ! Array with variable names
    integer(is)                                 :: n
    character (len=1)                           :: c
    real(rp)                                    :: r
    logical(lg)                                 :: err
    integer(ip)                                 :: ParCnt    ! Parenthesis counter
    integer(ip)                                 :: j,ib,in,lFunc

    j = 1
    ParCnt = 0
    lFunc = len_trim(Func)
    step: do
       if (j > lFunc) call ParseErrMsg (j, FuncStr)
       c = Func(j:j)
       !
       ! Check for valid operand (must appear)
       !
       if (c == '-' .or. c == '+') then                      ! Check for leading - or +
          j = j+1
          if (j > lFunc) call ParseErrMsg (j, FuncStr, 'Missing operand')
          c = Func(j:j)
          if (any(c == Ops)) call ParseErrMsg (j, FuncStr, 'Multiple operators')
       end if
       n = MathFunctionIndex (Func(j:))
       if (n > 0_1) then                                       ! Check for math function
          j = j+len_trim(Funcs(n))
          if (j > lFunc) call ParseErrMsg (j, FuncStr, 'Missing function argument')
          c = Func(j:j)
          if (c /= '(') call ParseErrMsg (j, FuncStr, 'Missing opening parenthesis')
       end if
       if (c == '(') then                                    ! Check for opening parenthesis
          ParCnt = ParCnt+1
          j = j+1
          cycle step
       end if

       if (scan(c,'0123456789.') > 0) then                   ! Check for number
          r = RealNum (Func(j:),ib,in,err)
          if (err) call ParseErrMsg (j, FuncStr, 'Invalid number format: '//Func(j+ib-1:j+in-2))
          j = j+in-1
          if (j > lFunc) exit
          c = Func(j:j)
       !else if (scan(c,'pi') > 0) then                       ! Special numbers
       !   r = 3.1_rp
       !   j = j+2    
       !   if (j > lFunc) exit
       !   c = Func(j:j)
       else                                                  ! Check for variable
          n = VariableIndex (Func(j:),Var,ib,in)
          if (n == 0_1) call ParseErrMsg (j, FuncStr, 'Invalid element: '//Func(j+ib-1:j+in-2))
          j = j+in-1
          if (j > lFunc) exit
          c = Func(j:j)
       end if
       do while (c == ')')                                   ! Check for closing parenthesis
          ParCnt = ParCnt-1
          if (ParCnt < 0) call ParseErrMsg (j, FuncStr, 'Mismatched parenthesis')
          if (Func(j-1:j-1) == '(') call ParseErrMsg (j-1, FuncStr, 'Empty parentheses')
          j = j+1
          if (j > lFunc) exit
          c = Func(j:j)
       end do
       !
       ! Now, we have a legal operand: A legal operator or end of string must follow
       !
       if (j > lFunc) exit
       if (any(c == Ops)) then                               ! Check for multiple operators
          if (j+1 > lFunc) call ParseErrMsg (j, FuncStr)
          if (any(Func(j+1:j+1) == Ops)) call ParseErrMsg (j+1, FuncStr, 'Multiple operators')
       else                                                  ! Check for next operand
          call ParseErrMsg (j, FuncStr, 'Missing operator')
       end if
       !
       ! Now, we have an operand and an operator: the next loop will check for another
       ! operand (must appear)
       !
       j = j+1
    end do step
    if (ParCnt > 0) call ParseErrMsg (j, FuncStr, 'Missing )')
  end subroutine CheckSyntax
 
  function EvalErrMsg () result (msg)

    !--------------------------------------------------------------------
    !
    ! Return error message
    !
    !--------------------------------------------------------------------

    implicit none
    character (len=32), dimension(4) :: m(4)
    character (len=len(m))           :: msg

    m(1) = 'Division by zero                '
    m(2) = 'Argument of SQRT negative       '
    m(3) = 'Argument of LOG negative        '
    m(4) = 'Argument of ASIN or ACOS illegal'
    if (EvalErrType < 1 .or. EvalErrType > size(m)) then
       msg = ''
    else
       msg = m(EvalErrType)
    endif
  end function EvalErrMsg
  
  subroutine ParseErrMsg (j, FuncStr, Msg)

    !--------------------------------------------------------------------
    !
    ! Print error message and terminate program
    !
    !--------------------------------------------------------------------
   
    implicit none
    integer(ip),                 intent(in) :: j
    character(len=*),            intent(in) :: FuncStr       ! Original function string
    character(len=*), optional,  intent(in) :: Msg
    integer(ip)                             :: k
    
    if (present(Msg)) then
       write(*,*) '*** Error in syntax of function string: '//Msg
    else
       write(*,*) '*** Error in syntax of function string:'
    endif
    write(*,*)
    write(*,'(A)') ' '//FuncStr
    do k=1,ipos(j)
       write(*,'(A)',ADVANCE='NO') ' '                       ! Advance to the jth position
    end do
    write(*,'(A)') '?'
    stop
  end subroutine ParseErrMsg
  
  function OperatorIndex (c) result (n)

    !--------------------------------------------------------------------
    !
    ! Return operator index
    !
    !--------------------------------------------------------------------

    implicit none
    character (len=1), intent(in) :: c
    integer(is)                   :: n,j

    n = 0
    do j=fadd,fpow
       if (c == Ops(j)) then
          n = j
          exit
       end if
    end do
  end function OperatorIndex
  
  function MathFunctionIndex (str) result (n)

    !--------------------------------------------------------------------
    !
    ! Return index of math function beginnig at 1st position of string str
    !
    !--------------------------------------------------------------------

    implicit none
    character(len=*), intent(in) :: str
    integer(is)                  :: n,j
    integer(ip)                  :: k
    character(len=len(Funcs))    :: fun

    n = 0
    do j=fAbs,VarBegin-1                                     ! Check all math functions
       k = min(len_trim(Funcs(j)), len(str))
       call LowCase (str(1:k), fun)
       if (fun == Funcs(j)) then                             ! Compare lower case letters
          n = j                                              ! Found a matching function
          exit
       end if
    end do
  end function MathFunctionIndex
  
  function VariableIndex (str, Var, ibegin, inext) result (n)

    !--------------------------------------------------------------------
    !
    ! Return index of variable at begin of string str (returns 0 if no variable found)
    !
    !--------------------------------------------------------------------

    implicit none
    character(len=*),               intent(in)  :: str       ! String
    character(len=*), dimension(:), intent(in)  :: Var       ! Array with variable names
    integer(is)                                 :: n         ! Index of variable
    integer(ip), optional,          intent(out) :: ibegin    ! Start position of variable name
    integer(ip), optional,          intent(out) :: inext     ! Position of character after name
    integer(ip)                                 :: j,ib
    integer(ip)                                 :: in,lstr

    n = 0
    !print*, "IN VariableIndex:", str, Var!, ibegin, inext
    lstr = len_trim(str)
    if( lstr > 0 ) then
       do ib = 1,lstr                                        ! Search for first character in str
          if( str(ib:ib) /= ' ' ) exit                       ! When lstr>0 at least 1 char in str
       end do
       do in = ib,lstr                                       ! Search for name terminators
          if( scan(str(in:in),'+-*/^) ') > 0 ) exit
       end do
       do j = 1,size(Var)
          if( str(ib:in-1) == Var(j) ) then
             n = int(j,is)                                   ! Variable name found
             exit
          end if
       end do
    end if
    if( present(ibegin) ) ibegin = ib
    if( present(inext)  ) inext  = in
  end function VariableIndex
  
  subroutine RemoveSpaces (str)

    !--------------------------------------------------------------------
    !
    ! Remove Spaces from string, remember positions of characters in old string
    !
    !--------------------------------------------------------------------

    implicit none
    character(len=*), intent(inout) :: str
    integer(ip)                     :: k,lstr

    lstr = len_trim(str)
    ipos = (/ (k,k=1,lstr) /)
    k = 1
    do while( str(k:lstr) /= ' ' )
       if( str(k:k) == ' ' ) then
          str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
          ipos(k:lstr) = (/ ipos(k+1:lstr), 0_ip /)             ! Move 1 element to left
          k = k-1
       end if
       k = k+1
    end do
  end subroutine RemoveSpaces
  
  subroutine Replace (ca,cb,str)

    !--------------------------------------------------------------------
    !
    ! Replace ALL appearances of character set ca in string str by character set cb
    !
    !--------------------------------------------------------------------

    implicit none
    character(len=*),       intent(in)    :: ca
    character(len=len(ca)), intent(in)    :: cb                ! LEN(ca) must be LEN(cb)
    character(len=*),       intent(inout) :: str
    integer(ip)                           :: j,lca
    
    lca = len(ca)
    do j = 1,len_trim(str)-lca+1
       if( str(j:j+lca-1) == ca ) str(j:j+lca-1) = cb
    end do
  end subroutine Replace
  
  subroutine Compile (i, F, Var)

    !--------------------------------------------------------------------
    !
    ! Compile i-th function string F into bytecode
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),                    intent(in) :: i         ! Function identifier
    character(len=*),               intent(in) :: F         ! Function string
    character(len=*), dimension(:), intent(in) :: Var       ! Array with variable names
    integer(ip)                                :: istat,lentri

    !print*, "IN Compile", i, F, Var
    if (associated(Comp(i)%ByteCode)) &
         deallocate ( Comp(i)%ByteCode, Comp(i)%Immed, Comp(i)%Stack )
    Comp(i) % ByteCodeSize = 0
    Comp(i) % ImmedSize    = 0
    Comp(i) % StackSize    = 0
    Comp(i) % StackPtr     = 0
    lentri= len_trim(F)
    call CompileSubstr (i,F,1_ip,lentri,Var)               ! Compile string to determine size
    allocate ( &
         Comp(i) % ByteCode(Comp(i) % ByteCodeSize), &
         Comp(i) % Immed(Comp(i) % ImmedSize),       &
         Comp(i) % Stack(Comp(i) % StackSize),       &
         STAT = istat                            )
    if (istat /= 0) then
       write(*,*) '*** Parser error: Memmory allocation for byte code failed'
       stop
    else
       Comp(i) % ByteCodeSize = 0
       Comp(i) % ImmedSize    = 0
       Comp(i) % StackSize    = 0
       Comp(i) % StackPtr     = 0
       call CompileSubstr (i,F,1_ip,lentri,Var)            ! Compile string into bytecode
    end if
    
  end subroutine Compile
  
  subroutine AddCompiledByte (i, b)

    !--------------------------------------------------------------------
    !
    ! Add compiled byte to bytecode
    !
    !--------------------------------------------------------------------
    implicit none
    integer(ip),     intent(in) :: i                             ! Function identifier
    integer(is),     intent(in) :: b                             ! Value of byte to be added

    Comp(i) % ByteCodeSize = Comp(i) % ByteCodeSize + 1
    if (associated(Comp(i) % ByteCode)) Comp(i) % ByteCode(Comp(i) % ByteCodeSize) = b

  end subroutine AddCompiledByte
 
  function MathItemIndex (i, F, Var) result (n)

    !--------------------------------------------------------------------
    !
    ! Return math item index, if item is real number, enter it into Comp-structure
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),                     intent(in) :: i         ! Function identifier
    character(len=*),                intent(in) :: F         ! Function substring
    character(len=*), dimension(:),  intent(in) :: Var       ! Array with variable names
    integer(is)                                 :: n         ! Byte value of math item
    
    n = 0
    !print*, "IN MathItemIndex", i, F, Var
    if (scan(F(1:1),'0123456789.') > 0) then                 ! Check for begin of a number
       Comp(i) % ImmedSize = Comp(i) % ImmedSize + 1
       if (associated(Comp(i) % Immed)) Comp(i) % Immed(Comp(i) % ImmedSize) = RealNum (F)
       n = cImmed
    else                                                     ! Check for a variable
       !print*, "B4 VariableIndex ", F, Var
       n = VariableIndex (F, Var)
       if (n > 0_1) n = int(VarBegin,is)+n-1_is
    end if
  end function MathItemIndex
  
  function CompletelyEnclosed (F, b, e) result (res)

    !--------------------------------------------------------------------
    !
    ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
    !
    !--------------------------------------------------------------------

    implicit none
    character(len=*), intent(in) :: F                       ! Function substring
    integer(ip),      intent(in) :: b,e                     ! First and last pos. of substring
    logical(lg)                  :: res
    integer(ip)                  :: j,k

    res=.false.
    if (F(b:b) == '(' .and. F(e:e) == ')') then
       k = 0
       do j=b+1,e-1
          if     (F(j:j) == '(') then
             k = k+1
          elseif (F(j:j) == ')') then
             k = k-1
          end if
          if (k < 0) exit
       end do
       if (k == 0) res=.true.                                ! All opened parenthesis closed
    end if

  end function CompletelyEnclosed
  
  recursive subroutine CompileSubstr (i, F, b, e, Var)

    !--------------------------------------------------------------------
    !
    ! Compile i-th function string F into bytecode
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),                     intent(in) :: i         ! Function identifier
    character (len=*),               intent(in) :: F         ! Function substring
    integer(ip),                     intent(in) :: b,e       ! Begin and end position substring
    character (len=*), dimension(:), intent(in) :: Var       ! Array with variable names
    integer(is)                                 :: n
    integer(ip)                                 :: b2,j,k,io
    !
    ! Check for special cases of substring
    !
    !print*, "IN CompileSubstr", i, F, b, e, Var
    if     (F(b:b) == '+') then                              ! Case 1: F(b:e) = '+...'
       !      WRITE(*,*)'1. F(b:e) = "+..."'
       call CompileSubstr (i, F, b+1, e, Var)
       return
    elseif (CompletelyEnclosed (F, b, e)) then               ! Case 2: F(b:e) = '(...)'
       !      WRITE(*,*)'2. F(b:e) = "(...)"'
       call CompileSubstr (i, F, b+1, e-1, Var)
       return
    elseif (scan(F(b:b),calpha) > 0) then
       n = MathFunctionIndex (F(b:e))
       if (n > 0_1) then
          b2 = b+index(F(b:e),'(')-1
          if (CompletelyEnclosed(F, b2, e)) then             ! Case 3: F(b:e) = 'fcn(...)'
             !            WRITE(*,*)'3. F(b:e) = "fcn(...)"'
             call CompileSubstr(i, F, b2+1, e-1, Var)
             call AddCompiledByte (i, n)
             return
          end if
       end if
    elseif (F(b:b) == '-') then
       if (CompletelyEnclosed (F, b+1, e)) then              ! Case 4: F(b:e) = '-(...)'
          !         WRITE(*,*)'4. F(b:e) = "-(...)"'
          call CompileSubstr (i, F, b+2, e-1, Var)
          call AddCompiledByte (i, fneg)
          return
       elseif (scan(F(b+1:b+1),calpha) > 0) then
          n = MathFunctionIndex (F(b+1:e))
          if (n > 0_1) then
             b2 = b+index(F(b+1:e),'(')
             if (CompletelyEnclosed(F, b2, e)) then          ! Case 5: F(b:e) = '-fcn(...)'
                !               WRITE(*,*)'5. F(b:e) = "-fcn(...)"'
                call CompileSubstr(i, F, b2+1, e-1, Var)
                call AddCompiledByte (i, n)
                call AddCompiledByte (i, fneg)
                return
             end if
          end if
       endif
    end if
    !
    ! Check for operator in substring: check only base level (k=0), exclude expr. in ()
    !
    do io=fadd,fpow                                          ! Increasing priority +-*/^
       k = 0
       do j=e,b,-1
          if     (F(j:j) == ')') then
             k = k+1
          elseif (F(j:j) == '(') then
             k = k-1
          end if
          if (k == 0 .and. F(j:j) == Ops(io) .and. IsBinaryOp (j, F)) then
             if (any(F(j:j) == Ops(fmul:fpow)) .and. F(b:b) == '-') then ! Case 6: F(b:e) = '-...Op...' with Op > -
                !               WRITE(*,*)'6. F(b:e) = "-...Op..." with Op > -'
                call CompileSubstr (i, F, b+1, e, Var)
                call AddCompiledByte (i, fneg)
                return
             else                                                        ! Case 7: F(b:e) = '...BinOp...'
                !               WRITE(*,*)'7. Binary operator',F(j:j)
                call CompileSubstr (i, F, b, j-1, Var)
                call CompileSubstr (i, F, j+1, e, Var)
                call AddCompiledByte (i, OperatorIndex(Ops(io)))
                Comp(i)%StackPtr = Comp(i)%StackPtr - 1
                return
             end if
          end if
       end do
    end do
    !
    ! Check for remaining items, i.e. variables or explicit numbers
    !
    b2 = b
    if (F(b:b) == '-') b2 = b2+1
    !print*, "B4 MathItemIndex", i, F(b2:e), Var
    n = MathItemIndex(i, F(b2:e), Var)
    !   WRITE(*,*)'8. AddCompiledByte ',n
    call AddCompiledByte (i, n)
    Comp(i)%StackPtr = Comp(i)%StackPtr + 1
    if (Comp(i)%StackPtr > Comp(i)%StackSize) Comp(i)%StackSize = Comp(i)%StackSize + 1
    if (b2 > b) call AddCompiledByte (i, fneg)

  end subroutine CompileSubstr
  
  function IsBinaryOp (j, F) result (res)

    !--------------------------------------------------------------------
    !
    ! Check if operator F(j:j) in string F is binary operator
    ! Special cases already covered elsewhere:              (that is corrected in v1.1)
    ! - operator character F(j:j) is first character of string (j=1)
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),       intent(in) :: j                       ! Position of Operator
    character (len=*), intent(in) :: F                       ! String
    logical(lg)                   :: res                     ! Result
    integer(ip)                   :: k
    logical(lg)                   :: Dflag,Pflag

    res=.true.
    if (F(j:j) == '+' .or. F(j:j) == '-') then               ! Plus or minus sign:
       if (j == 1) then                                      ! - leading unary operator ?
          res = .false.
       elseif (scan(F(j-1:j-1),'+-*/^(') > 0) then           ! - other unary operator ?
          res = .false.
       elseif (scan(F(j+1:j+1),'0123456789') > 0 .and. &     ! - in exponent of real number ?
            scan(F(j-1:j-1),'eEdD')       > 0) then
          Dflag=.false.; Pflag=.false.
          k = j-1
          do while (k > 1)                                   !   step to the left in mantissa
             k = k-1
             if     (scan(F(k:k),'0123456789') > 0) then
                Dflag=.true.
             elseif (F(k:k) == '.') then
                if (Pflag) then
                   exit                                      !   * EXIT: 2nd appearance of '.'
                else
                   Pflag=.true.                              !   * mark 1st appearance of '.'
                endif
             else
                exit                                         !   * all other characters
             end if
          end do
          if (Dflag .and. (k == 1 .or. scan(F(k:k),'+-*/^(') > 0)) res = .false.
       end if
    end if
  end function IsBinaryOp
  
  function RealNum (str, ibegin, inext, error) result (res)

    !--------------------------------------------------------------------
    !
    ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
    !
    !--------------------------------------------------------------------

    implicit none
    character (len=*),     intent(in)  :: str                    ! String
    real(rp)                           :: res                    ! Real number
    integer(ip), optional, intent(out) :: ibegin                 ! Start position of real number
    integer(ip), optional, intent(out) :: inext                  ! 1st character after real number
    logical(lg), optional, intent(out) :: error                  ! Error flag
    integer(ip)                        :: ib,in,istat
    logical(lg)                        :: Bflag,               & ! .T. at begin of number in str
         &                                InMan,               & ! .T. in mantissa of number
         &                                Pflag,               & ! .T. after 1st '.' encountered
         &                                Eflag,               & ! .T. at exponent identifier 'eEdD'
         &                                InExp,               & ! .T. in exponent of number
         &                                DInMan,              & ! .T. if at least 1 digit in mant.
         &                                DInExp,              & ! .T. if at least 1 digit in exp.
         &                                err                    ! Local error flag

    Bflag=.true.; InMan=.false.; Pflag=.false.; Eflag=.false.; InExp=.false.
    DInMan=.false.; DInExp=.false.
    ib   = 1
    in   = 1
    do while (in <= len_trim(str))
       select case (str(in:in))
       case (' ')                                            ! Only leading blanks permitted
          ib = ib+1
          if (InMan .or. Eflag .or. InExp) exit
       case ('+','-')                                        ! Permitted only
          if     (Bflag) then
             InMan=.true.; Bflag=.false.                     ! - at beginning of mantissa
          elseif (Eflag) then
             InExp=.true.; Eflag=.false.                     ! - at beginning of exponent
          else
             exit                                            ! - otherwise STOP
          endif
       case ('0':'9')                                        ! Mark
          if     (Bflag) then
             InMan=.true.; Bflag=.false.                     ! - beginning of mantissa
          elseif (Eflag) then
             InExp=.true.; Eflag=.false.                     ! - beginning of exponent
          endif
          if (InMan) DInMan=.true.                           ! Mantissa contains digit
          if (InExp) DInExp=.true.                           ! Exponent contains digit
       case ('.')
          if     (Bflag) then
             Pflag=.true.                                    ! - mark 1st appearance of '.'
             InMan=.true.; Bflag=.false.                     !   mark beginning of mantissa
          elseif (InMan .and..not.Pflag) then
             Pflag=.true.                                    ! - mark 1st appearance of '.'
          else
             exit                                            ! - otherwise STOP
          end if
       case ('e','E','d','D')                                ! Permitted only
          if (InMan) then
             Eflag=.true.; InMan=.false.                     ! - following mantissa
          else
             exit                                            ! - otherwise STOP
          endif
       case DEFAULT
          exit                                               ! STOP at all other characters
       end select
       in = in+1
    end do
    err = (ib > in-1) .or. (.not.DInMan) .or.  ((Eflag.or.InExp).and..not.DInExp)
    if (err) then
       res = 0.0_rp
    else
       read(str(ib:in-1),*,IOSTAT=istat) res
       err = istat /= 0
    end if
    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in
    if (present(error))  error  = err

  end function RealNum
 
  subroutine LowCase (str1, str2)

    !--------------------------------------------------------------------
    !
    ! Transform upper case letters in str1 into lower case letters, result is str2
    !
    !--------------------------------------------------------------------
    implicit none
    character(len=*),   intent(in)  :: str1
    character(len=*),   intent(out) :: str2
    integer(ip)                     :: j,k
    character(len=*),   parameter   :: lc = 'abcdefghijklmnopqrstuvwxyz'
    character(len=*),   parameter   :: uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    str2 = str1
    do j = 1,len_trim(str1)
       k = index(uc,str1(j:j))
       if( k > 0 ) str2(j:j) = lc(k:k)
    end do

  end subroutine LowCase

end module mod_ker_space_time_function
