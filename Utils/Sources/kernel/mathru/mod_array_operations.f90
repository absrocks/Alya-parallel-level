!-----------------------------------------------------------------------
!
!> @addtogroup Array_Operations_Toolbox
!> Toolbox for array operations, like axpy, dot product, imitating BLAS.
!> @{
!> @name    ToolBox for array operations
!> @file    mod_array_operations.f90
!> @date    22/05/2015
!> @author  Guillaume Houzeaux
!> @brief   Array operations ("a la" BLAS)
!> @details The following subroutines are available:
!>          \verbatim
!>          AXPY ........ y = y + alpha * x
!>          AXPBY ....... y = beta * y + alpha * x
!>          COPY ........ y = x
!>          AXOZ ........ y = alpha * x (*,/,+,-) z
!>          CONST ....... y = alpha
!>          \endverbatim
!>          The following functions are available:
!>          \verbatim
!>          NORM2 ....... alpha = sqrt(x.x)
!>          DOT ......... alpha = x.y
!>          \endverbatim
!
!-----------------------------------------------------------------------

module mod_array_operations

  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : INOTMASTER
  use mod_communications, only : PAR_SUM
  implicit none
  private 

  interface array_operations_axpy
     module procedure array_operations_axpy_1     , &
          &           array_operations_axpy_2     , &
          &           array_operations_axpy_1n
  end interface array_operations_axpy
  interface array_operations_axpby
     module procedure array_operations_axpby_1    , &
          &           array_operations_axpby_2    , &
          &           array_operations_axpby_3    , &
          &           array_operations_axpby_1n   , &
          &           array_operations_axpby_2n   , &
          &           array_operations_axpby_12n  , &
          &           array_operations_axpby_21n 
  end interface array_operations_axpby

  interface array_operations_copy
     module procedure array_operations_copy_11  , &
          &           array_operations_copy_22  , &
          &           array_operations_copy_33  , &
          &           array_operations_copy_12  , &
          &           array_operations_copy_21  , &
          &           array_operations_copy_23  , &
          &           array_operations_copy_p1  , &
          &           array_operations_copy_p2    
  end interface array_operations_copy

  interface array_operations_axoz
     module procedure array_operations_axoz_11  , &
          &           array_operations_axoz_p1
  end interface array_operations_axoz

  interface array_operations_const
     module procedure array_operations_const_p1  , &
          &           array_operations_const_p2
  end interface array_operations_const

  interface array_operations_norm2
     module procedure array_operations_norm2_p1  , &
          &           array_operations_norm2_p2
  end interface array_operations_norm2

  interface array_operations_dot
     module procedure array_operations_dot_p1    , &
          &           array_operations_dot_p2
  end interface array_operations_dot

  interface array_operations_initialization
     module procedure array_operations_initialization_1, &
          &           array_operations_initialization_2
  end interface array_operations_initialization

  public :: array_operations_axpy
  public :: array_operations_axpby
  public :: array_operations_copy
  public :: array_operations_axoz
  public :: array_operations_const
  public :: array_operations_norm2
  public :: array_operations_dot
  public :: array_operations_initialization

contains

  !-----------------------------------------------------------------------
  !
  !> @brief   x = 0
  !> @details Array initialization
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_initialization_1(xx)
    real(rp),   pointer :: xx(:)
    integer(ip)         :: ii

    do ii = 1,size(xx,KIND=ip)
       xx(ii) = 0.0_rp
    end do

  end subroutine array_operations_initialization_1

  subroutine array_operations_initialization_2(xx)
    real(rp),   pointer :: xx(:,:)
    integer(ip)         :: ii,jj

    do ii = 1,size(xx,2,KIND=ip)
       do jj = 1,size(xx,1,KIND=ip)
          xx(jj,ii) = 0.0_rp
       end do
    end do

  end subroutine array_operations_initialization_2

  !-----------------------------------------------------------------------
  !
  !> @brief   y = y + alpha*x
  !> @details Array operation using or not using OMP: y = y + alpha*x
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_axpy_generic(nn,alpha,xx,yy,use_omp)
    integer(ip), intent(in)    :: nn
    real(rp),    intent(in)    :: alpha
    real(rp),    intent(in)    :: xx(*)
    real(rp),    intent(inout) :: yy(*)
    logical(lg), intent(in)    :: use_omp
    integer(ip)                :: ii

    if( INOTMASTER ) then
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC)  &
          !$OMP DEFAULT  ( NONE )              &
          !$OMP SHARED   ( xx, yy, alpha, nn ) &
          !$OMP PRIVATE  ( ii )            
          do ii = 1,nn
             yy(ii) = yy(ii) + alpha * xx(ii)
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,nn
             yy(ii) = yy(ii) + alpha * xx(ii)
          end do
       end if
    end if

  end subroutine array_operations_axpy_generic

  subroutine array_operations_axpy_1(alpha,xx,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),     pointer,  intent(in)    :: xx(:)
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       nn = size(xx,KIND=ip)       
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpy_generic(nn,alpha,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpy_1

  subroutine array_operations_axpy_1n(nn,alpha,xx,yy,worder)
    integer(ip),            intent(in)    :: nn
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpy_generic(nn,alpha,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpy_1n

  subroutine array_operations_axpy_2(alpha,xx,yy,worder)
    real(rp),     pointer,  intent(in)    :: xx(:,:)
    real(rp),     pointer,  intent(inout) :: yy(:,:)
    real(rp),               intent(in)    :: alpha
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       nn = size(xx,1,KIND=ip)*size(xx,2,KIND=ip)       
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) .or. size(xx,2,KIND=ip) /= size(yy,2,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpy_generic(nn,alpha,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpy_2

  !-----------------------------------------------------------------------
  !
  !> @brief   y = alpha*x + beta*y
  !> @details Array operation using or not using OMP: y = alpha*x + beta*y
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_axpby_generic(nn,alpha,beta,xx,yy,use_omp)
    integer(ip), intent(in)    :: nn
    real(rp),    intent(in)    :: alpha
    real(rp),    intent(in)    :: beta
    real(rp),    intent(in)    :: xx(*)
    real(rp),    intent(inout) :: yy(*)
    logical(lg), intent(in)    :: use_omp
    integer(ip)                :: ii

    if( use_omp ) then
       !$OMP PARALLEL DO SCHEDULE (STATIC)        &
       !$OMP DEFAULT  ( NONE )                    &
       !$OMP SHARED   ( xx, yy, alpha, beta, nn ) &
       !$OMP PRIVATE  ( ii )            
       do ii = 1,nn
          yy(ii) = alpha * xx(ii) + beta * yy(ii)
       end do
       !$OMP END PARALLEL DO
    else
       do ii = 1,nn
          yy(ii) = alpha * xx(ii) + beta * yy(ii)
       end do
    end if

  end subroutine array_operations_axpby_generic

  subroutine array_operations_axpby_1(alpha,beta,xx,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),     pointer,  intent(in)    :: xx(:)
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       nn = size(xx,KIND=ip)       
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(nn,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_1

  subroutine array_operations_axpby_2(alpha,beta,xx,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),     pointer,  intent(in)    :: xx(:,:)
    real(rp),     pointer,  intent(inout) :: yy(:,:)
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       nn = size(xx,1,KIND=ip)*size(xx,2,KIND=ip)       
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) .or. size(xx,2,KIND=ip) /= size(yy,2,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(nn,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_2

  subroutine array_operations_axpby_3(alpha,beta,xx,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),     pointer,  intent(in)    :: xx(:,:,:)
    real(rp),     pointer,  intent(inout) :: yy(:,:,:)
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       nn = size(xx,1,KIND=ip)*size(xx,2,KIND=ip)       
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) .or. size(xx,2,KIND=ip) /= size(yy,2,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(nn,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_3

  subroutine array_operations_axpby_1n(n1,alpha,beta,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(n1,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_1n

  subroutine array_operations_axpby_2n(n1,n2,alpha,beta,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),               intent(in)    :: xx(n1,*)
    real(rp),               intent(inout) :: yy(n1,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(n1*n2,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_2n

  subroutine array_operations_axpby_12n(n1,n2,alpha,beta,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(n2,*)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(n1*n2,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_12n

  subroutine array_operations_axpby_21n(n1,n2,alpha,beta,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(n1,*)
    real(rp),               intent(inout) :: yy(*)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(n1*n2,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_21n

  !-----------------------------------------------------------------------
  !
  !> @brief   y = x
  !> @details Array operation using or not using OMP: y = x
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_copy_generic(nn,xx,yy,use_omp)
    integer(ip), intent(in)    :: nn
    real(rp),    intent(in)    :: xx(*)
    real(rp),    intent(inout) :: yy(*)
    logical(lg), intent(in)    :: use_omp
    integer(ip)                :: ii

    if( use_omp ) then
       !$OMP PARALLEL DO SCHEDULE (STATIC)  &
       !$OMP DEFAULT  ( NONE )              &
       !$OMP SHARED   ( xx, yy, nn )        &
       !$OMP PRIVATE  ( ii )            
       do ii = 1,nn
          yy(ii) = xx(ii)
       end do
       !$OMP END PARALLEL DO
    else
       do ii = 1,nn
          yy(ii) = xx(ii)
       end do
    end if

  end subroutine array_operations_copy_generic

  subroutine array_operations_copy_11(n1,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_11

  subroutine array_operations_copy_22(n1,n2,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(n1,*)
    real(rp),               intent(out)   :: yy(n1,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_22

  subroutine array_operations_copy_33(n1,n2,n3,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    integer(ip),            intent(in)    :: n3
    real(rp),               intent(in)    :: xx(n1,n2,*)
    real(rp),               intent(inout) :: yy(n1,n2,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2*n3,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_33

  subroutine array_operations_copy_12(n1,n2,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(n1,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_12

  subroutine array_operations_copy_21(n1,n2,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(n1,*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_21

  subroutine array_operations_copy_23(n1,n2,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(n1,n2,*)
    real(rp),               intent(inout) :: yy(n1,n2,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_23

  subroutine array_operations_copy_p1(xx,yy,worder)
    real(rp),     pointer,  intent(in)    :: xx(:)
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp
    integer(ip)                           :: nn

    if( INOTMASTER ) then
       nn = size(xx,1,KIND=ip)
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(nn,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_p1

  subroutine array_operations_copy_p2(xx,yy,worder)
    real(rp),     pointer,   intent(in)    :: xx(:,:)
    real(rp),     pointer,   intent(inout) :: yy(:,:)
    character(*), optional,  intent(in)    :: worder
    logical(lg)                            :: use_omp
    integer(ip)                            :: n1,n2

    if( INOTMASTER ) then
       n1 = size(xx,1,KIND=ip)
       n2 = size(xx,2,KIND=ip)
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_p2

  !-----------------------------------------------------------------------
  !
  !> @brief   y = alpha*x (*,/,+,-) z
  !> @details Array operation using or not using OMP: y = alpha*x (*,/,+,-) z
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_axoz_generic(nn,alpha,woperation,xx,zz,yy,use_omp)
    integer(ip),  intent(in)    :: nn
    real(rp),     intent(in)    :: alpha
    character(1), intent(in)    :: woperation
    real(rp),     intent(in)    :: xx(*)
    real(rp),     intent(in)    :: zz(*)
    real(rp),     intent(inout) :: yy(*)
    logical(lg),  intent(in)    :: use_omp
    integer(ip)                 :: ii

    if( INOTMASTER ) then

       if( woperation == '*' ) then
          !
          ! y = alpha * x * z
          !
          if( use_omp ) then
             !$OMP PARALLEL DO SCHEDULE (STATIC)      &
             !$OMP DEFAULT  ( NONE )                  &
             !$OMP SHARED   ( xx, zz, yy, alpha, nn ) &
             !$OMP PRIVATE  ( ii )            
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) * zz(ii)
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) * zz(ii)
             end do
          end if

       else if( woperation == '/' ) then
          !
          ! y = alpha * x / z
          !
          if( use_omp ) then
             !$OMP PARALLEL DO SCHEDULE (STATIC)      &
             !$OMP DEFAULT  ( NONE )                  &
             !$OMP SHARED   ( xx, zz, yy, alpha, nn ) &
             !$OMP PRIVATE  ( ii )            
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) / zz(ii)
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) / zz(ii)
             end do
          end if

       else if( woperation == '+' ) then
          !
          ! y = alpha * x + z
          !
          if( use_omp ) then
             !$OMP PARALLEL DO SCHEDULE (STATIC)      &
             !$OMP DEFAULT  ( NONE )                  &
             !$OMP SHARED   ( xx, zz, yy, alpha, nn ) &
             !$OMP PRIVATE  ( ii )            
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) + zz(ii)
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) + zz(ii)
             end do
          end if

       else if( woperation == '-' ) then
          !
          ! y = alpha * x - z
          !
          if( use_omp ) then
             !$OMP PARALLEL DO SCHEDULE (STATIC)      &
             !$OMP DEFAULT  ( NONE )                  &
             !$OMP SHARED   ( xx, zz, yy, alpha, nn ) &
             !$OMP PRIVATE  ( ii )            
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) - zz(ii)
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) - zz(ii)
             end do
          end if

       end if

    end if

  end subroutine array_operations_axoz_generic

  subroutine array_operations_axoz_11(nn,alpha,woperation,xx,zz,yy,worder)
    integer(ip),            intent(in)    :: nn
    real(rp),               intent(in)    :: alpha
    character(1),           intent(in)    :: woperation
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(in)    :: zz(*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axoz_generic(nn,alpha,woperation,xx,zz,yy,use_omp)
    end if

  end subroutine array_operations_axoz_11

  subroutine array_operations_axoz_p1(alpha,woperation,xx,zz,yy,worder)
    real(rp),               intent(in)    :: alpha
    character(1),           intent(in)    :: woperation
    real(rp),     pointer,  intent(in)    :: xx(:)
    real(rp),     pointer,  intent(in)    :: zz(:)
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp
    integer(ip)                           :: nn

    if( INOTMASTER ) then
       nn = size(xx,KIND=ip)
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axoz_generic(nn,alpha,woperation,xx,zz,yy,use_omp)
    end if

  end subroutine array_operations_axoz_p1

  subroutine array_operations_const_p1(alpha,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp
    integer(ip)                           :: ii

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC) &
          !$OMP DEFAULT  ( NONE )             &
          !$OMP SHARED   ( yy, alpha )        &
          !$OMP PRIVATE  ( ii )            
          do ii = 1,size(yy,KIND=ip)
             yy(ii) = alpha 
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,size(yy,KIND=ip)
             yy(ii) = alpha 
          end do
       end if
    end if

  end subroutine array_operations_const_p1

  subroutine array_operations_const_p2(alpha,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),     pointer,  intent(inout) :: yy(:,:)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp
    integer(ip)                           :: ii,jj,n1,n2

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       n1 = size(yy,1,KIND=ip)
       n2 = size(yy,2,KIND=ip)
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC)  &
          !$OMP DEFAULT  ( NONE )              &
          !$OMP SHARED   ( yy, alpha, n1, n2 ) &
          !$OMP PRIVATE  ( ii, jj )            
          do jj = 1,n2
             do ii = 1,n1
                yy(ii,jj) = alpha
             end do
          end do
          !$OMP END PARALLEL DO
       else
          do jj = 1,n2
             do ii = 1,n1
                yy(ii,jj) = alpha
             end do
          end do
       end if
    end if

  end subroutine array_operations_const_p2

  subroutine array_operations_norm2_generic(nn,yy,yynorm2,use_omp)
    integer(ip), intent(in)  :: nn
    real(rp),    intent(in)  :: yy(*)
    real(rp),    intent(out) :: yynorm2
    logical(lg), intent(in)  :: use_omp
    integer(ip)              :: ii

    yynorm2 = 0.0_rp

    if( INOTMASTER ) then    
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC) &
          !$OMP DEFAULT  ( NONE )             &
          !$OMP PRIVATE  ( ii )               &
          !$OMP SHARED   ( yy, nn )           &
          !$OMP REDUCTION (+:yynorm2)  
          do ii = 1,nn
             yynorm2 = yynorm2 + yy(ii) * yy(ii)
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,nn
             yynorm2 = yynorm2 + yy(ii) * yy(ii)
          end do
       end if
    end if
    call PAR_SUM(yynorm2,'IN MY CODE')
    yynorm2 = sqrt(yynorm2)

  end subroutine array_operations_norm2_generic

  function array_operations_norm2_p1(yy,worder)
    real(rp),     pointer,  intent(in) :: yy(:) 
    character(*), optional, intent(in) :: worder
    real(rp)                           :: array_operations_norm2_p1
    logical(lg)                        :: use_omp
    integer(ip)                        :: nn

    if( present(worder) ) then
       if( worder == 'USE OPENMP' ) then
          use_omp = .true.
       else if( worder == 'DO NOT USE OPENMP' ) then
          use_omp = .false.
       end if
    else
       use_omp = .true.
    end if
    nn = size(yy,KIND=ip)
    call array_operations_norm2_generic(nn,yy,array_operations_norm2_p1,use_omp)

  end function array_operations_norm2_p1

  function array_operations_norm2_p2(yy,worder)
    real(rp),     pointer,  intent(in) :: yy(:,:) 
    character(*), optional, intent(in) :: worder
    real(rp)                           :: array_operations_norm2_p2
    logical(lg)                        :: use_omp
    integer(ip)                        :: nn

    if( present(worder) ) then
       if( worder == 'USE OPENMP' ) then
          use_omp = .true.
       else if( worder == 'DO NOT USE OPENMP' ) then
          use_omp = .false.
       end if
    else
       use_omp = .true.
    end if
    nn = size(yy,KIND=ip)
    call array_operations_norm2_generic(nn,yy,array_operations_norm2_p2,use_omp)

  end function array_operations_norm2_p2

  subroutine array_operations_dot_generic(nn,xx,yy,xdoty,use_omp)
    integer(ip), intent(in)  :: nn
    real(rp),    intent(in)  :: xx(*)
    real(rp),    intent(in)  :: yy(*)
    real(rp),    intent(out) :: xdoty
    logical(lg), intent(in)  :: use_omp
    integer(ip)              :: ii

    xdoty = 0.0_rp

    if( INOTMASTER ) then    
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC) &
          !$OMP DEFAULT  ( NONE )             &
          !$OMP PRIVATE  ( ii )               &
          !$OMP SHARED   ( yy, nn )           &
          !$OMP REDUCTION (+:xdoty)  
          do ii = 1,nn
             xdoty = xdoty + xx(ii) * yy(ii)
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,nn
             xdoty = xdoty + xx(ii) * yy(ii)
          end do
       end if
    end if
    call PAR_SUM(xdoty,'IN MY CODE')

  end subroutine array_operations_dot_generic

  function array_operations_dot_p1(xx,yy,worder)
    real(rp),     pointer,  intent(in) :: xx(:) 
    real(rp),     pointer,  intent(in) :: yy(:) 
    character(*), optional, intent(in) :: worder
    real(rp)                           :: array_operations_dot_p1
    logical(lg)                        :: use_omp
    integer(ip)                        :: nn

    if( present(worder) ) then
       if( worder == 'USE OPENMP' ) then
          use_omp = .true.
       else if( worder == 'DO NOT USE OPENMP' ) then
          use_omp = .false.
       end if
    else
       use_omp = .true.
    end if
    nn = size(yy,KIND=ip)
    call array_operations_dot_generic(nn,xx,yy,array_operations_dot_p1,use_omp)

  end function array_operations_dot_p1

  function array_operations_dot_p2(xx,yy,worder)
    real(rp),     pointer,  intent(in) :: xx(:,:) 
    real(rp),     pointer,  intent(in) :: yy(:,:) 
    character(*), optional, intent(in) :: worder
    real(rp)                           :: array_operations_dot_p2
    logical(lg)                        :: use_omp
    integer(ip)                        :: nn

    if( present(worder) ) then
       if( worder == 'USE OPENMP' ) then
          use_omp = .true.
       else if( worder == 'DO NOT USE OPENMP' ) then
          use_omp = .false.
       end if
    else
       use_omp = .true.
    end if
    nn = size(yy,KIND=ip)
    call array_operations_dot_generic(nn,xx,yy,array_operations_dot_p2,use_omp)

  end function array_operations_dot_p2

end module mod_array_operations
!> @}
