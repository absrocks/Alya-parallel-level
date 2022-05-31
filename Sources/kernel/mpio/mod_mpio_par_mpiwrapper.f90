!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_mpiwrapper.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO Wrapper
!> @details This module is a wrapper of MPI-IO functions.
!>          It adds error management and makes type management transparent to the user.
!> @}
!-----------------------------------------------------------------------


module mod_mpio_par_mpiwrapper

  use def_kintyp,                     only : ip,rp,lg,r1p
  use def_master,                     only : comm_data_par
  use def_master,                     only : intost, retost, kfl_outfo
  use def_master,                     only : IPARSLAVE, IMASTER, INOTMASTER
  use mod_memory,                     only : memory_alloca, memory_deallo
  use def_domain,                     only : ndime
  use def_mpio,                       only : mpio_memor

#ifdef MPI_OFF
#define MPI_OFFSET_KIND 8
#endif

#ifndef MPI_OFF
  use mod_communications,             only : PAR_DEFINE_COMMUNICATOR
#endif

  implicit none

  private
#ifndef MPI_OFF
  include 'mpif.h'
#endif

#ifndef MPI_OFF
  integer(ip)                             ::  ierr                                ! error flag
  integer                                 ::  info=MPI_INFO_NULL                  ! info flag
  integer, dimension (MPI_STATUS_SIZE)    ::  status                              ! MPI status
#endif

  character(150)                          ::  wherein_world="IN THE WORLD"                    !
  character(150)                          ::  wherein_code="IN MY CODE"                       !
  character(150)                          ::  wherein_p="IN MPIO WITH MASTER"     !
  character(150)                          ::  wherein_pwm="IN MPIO"               !
  character(23)                           ::  vacal = "mod_mpio_par_mpiwrapper"
  interface PAR_FILE_SET_VIEW
     module procedure & 
          PAR_FILE_SET_VIEW_INT4_0  , &
          PAR_FILE_SET_VIEW_INT8_0  , &
          PAR_FILE_SET_VIEW_REAL4_0 , &
          PAR_FILE_SET_VIEW_REAL8_0 , &
          PAR_FILE_SET_VIEW_INT4_V  , &
          PAR_FILE_SET_VIEW_INT4_M  , &
          PAR_FILE_SET_VIEW_INT8_V  , &
          PAR_FILE_SET_VIEW_INT8_M  , &
          PAR_FILE_SET_VIEW_REAL4_V , &
          PAR_FILE_SET_VIEW_REAL4_M , &
#ifndef __PGI
          PAR_FILE_SET_VIEW_REAL16_V, &
          PAR_FILE_SET_VIEW_REAL16_M, &
#endif
          PAR_FILE_SET_VIEW_REAL8_V , &
          PAR_FILE_SET_VIEW_REAL8_M
  end interface PAR_FILE_SET_VIEW


  interface PAR_FILE_READ
     module procedure                        PAR_FILE_READ_INT4,&
          PAR_FILE_READ_INT4_V,&
          PAR_FILE_READ_INT4_M,&
          PAR_FILE_READ_INT8,&
          PAR_FILE_READ_INT8_V,&
          PAR_FILE_READ_INT8_M,&
          PAR_FILE_READ_REAL8,&
          PAR_FILE_READ_REAL8_V,&
          PAR_FILE_READ_REAL8_M,&
          PAR_FILE_READ_REAL4,&
          PAR_FILE_READ_REAL4_V,&
          PAR_FILE_READ_REAL4_M,&
#ifndef __PGI
          PAR_FILE_READ_REAL16,&
          PAR_FILE_READ_REAL16_V,&
          PAR_FILE_READ_REAL16_M,&
#endif
          PAR_FILE_READ_CHAR8,&
          PAR_FILE_READ_CHAR8_P
  end interface PAR_FILE_READ

  interface PAR_FILE_WRITE
     module procedure                        PAR_FILE_WRITE_INT4,&
          PAR_FILE_WRITE_INT4_V,&
          PAR_FILE_WRITE_INT4_M,&
          PAR_FILE_WRITE_INT8,&
          PAR_FILE_WRITE_INT8_V,&
          PAR_FILE_WRITE_INT8_M,&
          PAR_FILE_WRITE_REAL8,&
          PAR_FILE_WRITE_REAL8_V,&
          PAR_FILE_WRITE_REAL8_M,&
          PAR_FILE_WRITE_REAL4,&
          PAR_FILE_WRITE_REAL4_V,&
          PAR_FILE_WRITE_REAL4_M,&
#ifndef __PGI
          PAR_FILE_WRITE_REAL16,&
          PAR_FILE_WRITE_REAL16_V,&
          PAR_FILE_WRITE_REAL16_M,&
#endif
          PAR_FILE_WRITE_CHAR8,&
          PAR_FILE_WRITE_CHAR8_P
  end interface PAR_FILE_WRITE

  interface PAR_FILE_READ_ALL
     module procedure                        PAR_FILE_READ_ALL_INT4_V,&
          PAR_FILE_READ_ALL_INT4_M,&
          PAR_FILE_READ_ALL_INT8_V,&
          PAR_FILE_READ_ALL_INT8_M,&
          PAR_FILE_READ_ALL_REAL4_V,&
          PAR_FILE_READ_ALL_REAL4_M,&
#ifndef __PGI
          PAR_FILE_READ_ALL_REAL16_V,&
          PAR_FILE_READ_ALL_REAL16_M,&
#endif
          PAR_FILE_READ_ALL_REAL8_V,&
          PAR_FILE_READ_ALL_REAL8_M
  end interface PAR_FILE_READ_ALL

  interface PAR_FILE_WRITE_ALL
     module procedure                        PAR_FILE_WRITE_ALL_INT4_V,&
          PAR_FILE_WRITE_ALL_INT4_M,&
          PAR_FILE_WRITE_ALL_INT8_V,&
          PAR_FILE_WRITE_ALL_INT8_M,&
          PAR_FILE_WRITE_ALL_REAL4_V,&
          PAR_FILE_WRITE_ALL_REAL4_M,&
#ifndef __PGI
          PAR_FILE_WRITE_ALL_REAL16_V,&
          PAR_FILE_WRITE_ALL_REAL16_M,&
#endif
          PAR_FILE_WRITE_ALL_REAL8_V,&
          PAR_FILE_WRITE_ALL_REAL8_M
  end interface PAR_FILE_WRITE_ALL

  interface PAR_FILE_WRITE_ALL_BEGIN
     module procedure                        PAR_FILE_WRITE_ALL_BEGIN_INT4_V,&
          PAR_FILE_WRITE_ALL_BEGIN_INT4_M,&
          PAR_FILE_WRITE_ALL_BEGIN_INT8_V,&
          PAR_FILE_WRITE_ALL_BEGIN_INT8_M,&
          PAR_FILE_WRITE_ALL_BEGIN_REAL4_V,&
          PAR_FILE_WRITE_ALL_BEGIN_REAL4_M,&
#ifndef __PGI
          PAR_FILE_WRITE_ALL_BEGIN_REAL16_V,&
          PAR_FILE_WRITE_ALL_BEGIN_REAL16_M,&
#endif
          PAR_FILE_WRITE_ALL_BEGIN_REAL8_V,&
          PAR_FILE_WRITE_ALL_BEGIN_REAL8_M
  end interface PAR_FILE_WRITE_ALL_BEGIN

  interface PAR_FILE_WRITE_ALL_END
     module procedure                        PAR_FILE_WRITE_ALL_END_INT4_V,&
          PAR_FILE_WRITE_ALL_END_INT4_M,&
          PAR_FILE_WRITE_ALL_END_INT8_V,&
          PAR_FILE_WRITE_ALL_END_INT8_M,&
          PAR_FILE_WRITE_ALL_END_REAL4_V,&
          PAR_FILE_WRITE_ALL_END_REAL4_M,&
#ifndef __PGI
          PAR_FILE_WRITE_ALL_END_REAL16_V,&
          PAR_FILE_WRITE_ALL_END_REAL16_M,&
#endif
          PAR_FILE_WRITE_ALL_END_REAL8_V,&
          PAR_FILE_WRITE_ALL_END_REAL8_M
  end interface PAR_FILE_WRITE_ALL_END



  integer(ip)                             ::  i



  public                                  ::  PAR_FILE_OPEN_READ,&
       PAR_FILE_OPEN_WRITE,&
       PAR_INFO_CREATE,&
       PAR_INFO_FREE,&
       PAR_INFO_SET,&
       PAR_FILE_CLOSE,&
       PAR_FILE_SET_VIEW,&
       PAR_FILE_SEEK_SET,&
       PAR_FILE_SET_SIZE,&
       PAR_FILE_GET_SIZE,&
       PAR_FILE_READ,&
       PAR_FILE_WRITE,&
       PAR_FILE_READ_ALL,&
       PAR_FILE_WRITE_ALL,&
       PAR_FILE_WRITE_ALL_BEGIN,&
       PAR_FILE_WRITE_ALL_END


contains

  subroutine PAR_INFO_CREATE()
#ifndef MPI_OFF
    call MPI_Info_create(info, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("MPI info could not be created. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_INFO_CREATE

  subroutine PAR_INFO_SET(field, value)
    character(*),          intent(in)              :: field, value
#ifndef MPI_OFF
    call MPI_Info_set(info, field, value , ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("MPI info could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_INFO_SET

  subroutine PAR_INFO_FREE()
#ifndef MPI_OFF
    call MPI_Info_free(info, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("MPI info could not be set free. MPI error code: "//trim(intost(ierr)))
    end if
    info=MPI_INFO_NULL
#endif
  end subroutine PAR_INFO_FREE

  !-----------------------------------------------------------------------
  !> 
  !> @author  dosimont and houzeaux
  !> @date    2020-05-14
  !> @brief   Open file for reading
  !> @details Open file for reading
  !> 
  !-----------------------------------------------------------------------

  subroutine PAR_FILE_OPEN_READ(fh, filename, wherein, PAR_COMM)

    integer(ip),           intent(inout)           :: fh
    character(*),          intent(in)              :: filename
    character(*),          optional, intent(in)    :: wherein
    integer(ip),           optional, intent(in)    :: PAR_COMM
    type(comm_data_par),   pointer                 :: commu
    integer(4)                                     :: PAR_COMM_TO_USE
#ifndef MPI_OFF
    fh=0
    ierr=0

    if( present(PAR_COMM) ) then
       PAR_COMM_TO_USE = int(PAR_COMM,KIND=4)
    else
       call par_my_com(PAR_COMM_TO_USE, commu, wherein)
    end if
    call MPI_FILE_OPEN(PAR_COMM_TO_USE, filename, MPI_MODE_RDONLY, info, fh, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File "//trim(filename)//" could not be open. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_OPEN_READ

  !-----------------------------------------------------------------------
  !> 
  !> @author  dosimont and houzeaux
  !> @date    2020-05-14
  !> @brief   Open file for writting
  !> @details Open file for writting
  !> 
  !-----------------------------------------------------------------------

  subroutine PAR_FILE_OPEN_WRITE(fh, filename, wherein, PAR_COMM)
    use def_master 
    integer(ip),           intent(inout)           :: fh
    character(*),          intent(in)              :: filename
    character(*),          optional, intent(in)    :: wherein
    integer(ip),           optional, intent(in)    :: PAR_COMM
    type(comm_data_par),   pointer                 :: commu
    integer(4)                                     :: PAR_COMM_TO_USE
#ifndef MPI_OFF
    fh=0
    ierr=0
    if( present(PAR_COMM) ) then
       PAR_COMM_TO_USE = int(PAR_COMM,KIND=4)
    else
       call par_my_com(PAR_COMM_TO_USE, commu, wherein)
    end if
    call MPI_FILE_OPEN(PAR_COMM_TO_USE, filename, MPI_MODE_CREATE+MPI_MODE_RDWR, info, fh, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File "//trim(filename)//" could not be open nor created. MPI error code: "//trim(intost(ierr))//&
            '. Execute ulimit -n and increase value or run with SYNCHRONOUS: ON')
    end if
#endif
  end subroutine PAR_FILE_OPEN_WRITE

  subroutine PAR_FILE_CLOSE(fh)
    integer(ip),           intent(inout)           :: fh
    type(comm_data_par),   pointer                 :: commu
    integer(4)                                     :: PAR_COMM_TO_USE
#ifndef MPI_OFF
    ierr=0
    call MPI_FILE_CLOSE(fh, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be closed. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_CLOSE

  subroutine PAR_FILE_SET_VIEW_INT4_0(fh, buf, offset, count)
    integer(ip),              intent(in) :: fh
    integer(4),               intent(in) :: buf(*)
    integer(MPI_OFFSET_KIND), intent(in) :: offset
    integer(4),               intent(in) :: count
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_INT4_0

  subroutine PAR_FILE_SET_VIEW_INT8_0(fh, buf, offset, count)
    integer(ip),              intent(in) :: fh
    integer(8),               intent(in) :: buf(*)
    integer(MPI_OFFSET_KIND), intent(in) :: offset
    integer(8),               intent(in) :: count
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_INTEGER8, MPI_INTEGER8, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_INT8_0

  subroutine PAR_FILE_SET_VIEW_REAL4_0(fh, buf, offset, count)
    integer(ip),              intent(in) :: fh
    real(4),                  intent(in) :: buf(*)
    integer(MPI_OFFSET_KIND), intent(in) :: offset
    integer(ip),              intent(in) :: count
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_REAL4_0

  subroutine PAR_FILE_SET_VIEW_REAL8_0(fh, buf, offset, count)
    integer(ip),              intent(in) :: fh
    real(8),                  intent(in) :: buf(*)
    integer(MPI_OFFSET_KIND), intent(in) :: offset
    integer(ip),              intent(in) :: count
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_REAL8_0


  subroutine PAR_FILE_SET_VIEW_INT4_V(fh, buf, offset)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(in)              :: buf(:)
    integer(MPI_OFFSET_KIND),           intent(in)              :: offset
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_INT4_V

  subroutine PAR_FILE_SET_VIEW_INT4_M(fh, buf, offset)
    integer(ip),            intent(in)              :: fh
    integer(4), pointer,   intent(in)              :: buf(:,:)
    integer(MPI_OFFSET_KIND),           intent(in)              :: offset
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_INT4_M

  subroutine PAR_FILE_SET_VIEW_INT8_V(fh, buf, offset)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(in)              :: buf(:)
    integer(MPI_OFFSET_KIND),           intent(in)              :: offset
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_INTEGER8, MPI_INTEGER8, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_INT8_V

  subroutine PAR_FILE_SET_VIEW_INT8_M(fh, buf, offset)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(in)              :: buf(:,:)
    integer(MPI_OFFSET_KIND),           intent(in)              :: offset
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_INTEGER8, MPI_INTEGER8, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_INT8_M

  subroutine PAR_FILE_SET_VIEW_REAL8_V(fh, buf, offset)
    integer(ip),               intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:)
    integer(MPI_OFFSET_KIND),               intent(in)              :: offset
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_REAL8_V

  subroutine PAR_FILE_SET_VIEW_REAL8_M(fh, buf, offset)
    integer(ip),               intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:,:)
    integer(MPI_OFFSET_KIND),               intent(in)              :: offset
#ifndef MPI_OFF
    call MPI_FILE_SET_VIEW(fh, int(offset,MPI_OFFSET_KIND), MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" set view could not be set. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_VIEW_REAL8_M

  subroutine PAR_FILE_SET_VIEW_REAL4_V(fh, buf, offset)
    integer(ip),               intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:)
    integer(MPI_OFFSET_KIND),               intent(in)              :: offset
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_SET_VIEW_REAL4_V

  subroutine PAR_FILE_SET_VIEW_REAL4_M(fh, buf, offset)
    integer(ip),               intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:,:)
    integer(MPI_OFFSET_KIND),               intent(in)              :: offset
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_SET_VIEW_REAL4_M

#ifndef __PGI
  subroutine PAR_FILE_SET_VIEW_REAL16_V(fh, buf, offset)
    integer(ip),               intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:)
    integer(MPI_OFFSET_KIND),               intent(in)              :: offset
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_SET_VIEW_REAL16_V

  subroutine PAR_FILE_SET_VIEW_REAL16_M(fh, buf, offset)
    integer(ip),               intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:,:)
    integer(MPI_OFFSET_KIND),               intent(in)              :: offset
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_SET_VIEW_REAL16_M
#endif

  subroutine PAR_FILE_READ_INT4(fh, buf)
    integer(ip),           intent(in)              :: fh
    integer(4),            intent(inout)           :: buf
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, 1_4, MPI_INTEGER4, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_INT4

  subroutine PAR_FILE_READ_INT4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_INT4_V

  subroutine PAR_FILE_READ_INT4_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_INT4_M

  subroutine PAR_FILE_READ_INT8(fh, buf)
    integer(ip),           intent(in)              :: fh
    integer(8),            intent(inout)           :: buf
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, 1, MPI_INTEGER8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_INT8

  subroutine PAR_FILE_READ_INT8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_INT8_V

  subroutine PAR_FILE_READ_INT8_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_INT8_M

  subroutine PAR_FILE_READ_REAL8(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(8),               intent(inout)           :: buf
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, 1, MPI_REAL8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_REAL8

  subroutine PAR_FILE_READ_REAL8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, nsize, MPI_REAL8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_REAL8_V

  subroutine PAR_FILE_READ_REAL8_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, nsize, MPI_REAL8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_REAL8_M

  subroutine PAR_FILE_READ_REAL4(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(4),               intent(inout)           :: buf
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_REAL4

  subroutine PAR_FILE_READ_REAL4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_REAL4_V

  subroutine PAR_FILE_READ_REAL4_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_REAL4_M

#ifndef __PGI
  subroutine PAR_FILE_READ_REAL16(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(16),               intent(inout)           :: buf
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_REAL16

  subroutine PAR_FILE_READ_REAL16_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_REAL16_V

  subroutine PAR_FILE_READ_REAL16_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_REAL16_M
#endif

  subroutine PAR_FILE_READ_CHAR8(fh, buf)
    integer(ip),           intent(in)              :: fh
    character(8),          intent(inout)           :: buf
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, 8, MPI_CHARACTER, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_CHAR8

  subroutine PAR_FILE_READ_CHAR8_P(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    character(8), dimension(*), intent(inout)      :: buf
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_READ(fh, buf, 8*nsize, MPI_CHARACTER, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be read. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_CHAR8_P

  subroutine PAR_FILE_READ_ALL_INT4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
    integer(4), pointer                            :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_4)
       call MPI_FILE_READ_ALL(fh, dummy, 0_ip, MPI_INTEGER4, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_READ_ALL(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be READ_ALL. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_ALL_INT4_V

  subroutine PAR_FILE_READ_ALL_INT4_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    integer(4), pointer                            :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_4, 1_4)
       call MPI_FILE_READ_ALL(fh, dummy, 0, MPI_INTEGER4, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_READ_ALL(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be READ_ALL. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_ALL_INT4_M

  subroutine PAR_FILE_READ_ALL_INT8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
    integer(8), pointer                            :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_8)
       call MPI_FILE_READ_ALL(fh, dummy, 0, MPI_INTEGER8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_READ_ALL(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be READ_ALL. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_ALL_INT8_V

  subroutine PAR_FILE_READ_ALL_INT8_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    integer(8), pointer                            :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_8, 1_8)
       call MPI_FILE_READ_ALL(fh, dummy, 0, MPI_INTEGER8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_READ_ALL(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be READ_ALL. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_ALL_INT8_M

  subroutine PAR_FILE_READ_ALL_REAL8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
    real(8), pointer                               :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_ip)
       call MPI_FILE_READ_ALL(fh, dummy, 0, MPI_REAL8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_READ_ALL(fh, buf, nsize, MPI_REAL8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be READ_ALL. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_ALL_REAL8_V

  subroutine PAR_FILE_READ_ALL_REAL8_M(fh, buf, nsize)
    use def_master
    use mod_memory
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    real(8), pointer                               :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF

    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_ip, 1_ip)
       call MPI_FILE_READ_ALL(fh, dummy, 0, MPI_REAL8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_READ_ALL(fh, buf, nsize, MPI_REAL8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be READ_ALL. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_READ_ALL_REAL8_M

  subroutine PAR_FILE_READ_ALL_REAL4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
    real(4), pointer                               :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_ALL_REAL4_V

  subroutine PAR_FILE_READ_ALL_REAL4_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    real(4), pointer                               :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_ALL_REAL4_M

#ifndef __PGI
  subroutine PAR_FILE_READ_ALL_REAL16_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(inout)           :: buf(:)
    integer(ip),           intent(in)              :: nsize
    real(16), pointer                               :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_ALL_REAL16_V

  subroutine PAR_FILE_READ_ALL_REAL16_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(inout)           :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    real(16), pointer                               :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_READ_ALL_REAL16_M
#endif

  subroutine PAR_FILE_WRITE_INT4(fh, buf)
    integer(ip),           intent(in)              :: fh
    integer(4),            intent(in)              :: buf
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, 1, MPI_INTEGER4, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_INT4

  subroutine PAR_FILE_WRITE_INT4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_INT4_V

  subroutine PAR_FILE_WRITE_INT4_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_INT4_M

  subroutine PAR_FILE_WRITE_INT8(fh, buf)
    integer(ip),           intent(in)              :: fh
    integer(8),            intent(in)              :: buf
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, 1, MPI_INTEGER8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_INT8

  subroutine PAR_FILE_WRITE_INT8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_INT8_V

  subroutine PAR_FILE_WRITE_INT8_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_INT8_M

  subroutine PAR_FILE_WRITE_REAL8(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(8),               intent(in)              :: buf
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, 1, MPI_REAL8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_REAL8

  subroutine PAR_FILE_WRITE_REAL8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, nsize, MPI_REAL8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_REAL8_V

  subroutine PAR_FILE_WRITE_REAL8_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, nsize, MPI_REAL8, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_REAL8_M

  subroutine PAR_FILE_WRITE_REAL4(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(4),               intent(in)              :: buf
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_REAL4

  subroutine PAR_FILE_WRITE_REAL4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_REAL4_V

  subroutine PAR_FILE_WRITE_REAL4_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_REAL4_M

#ifndef __PGI
  subroutine PAR_FILE_WRITE_REAL16(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(16),               intent(in)              :: buf
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_REAL16

  subroutine PAR_FILE_WRITE_REAL16_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_REAL16_V

  subroutine PAR_FILE_WRITE_REAL16_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_REAL16_M
#endif

  subroutine PAR_FILE_WRITE_CHAR8(fh, buf)
    integer(ip),           intent(in)              :: fh
    character(8),          intent(in)              :: buf
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, 8, MPI_CHARACTER, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_CHAR8

  subroutine PAR_FILE_WRITE_CHAR8_P(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    character(8), dimension(*), intent(in)         :: buf
    integer(ip),           intent(in)              :: nsize
#ifndef MPI_OFF
    call MPI_FILE_WRITE(fh, buf, 8*nsize, MPI_CHARACTER, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_CHAR8_P

  subroutine PAR_FILE_WRITE_ALL_INT4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    integer(4), pointer                            :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_4)
       call MPI_FILE_WRITE_ALL(fh, dummy, 0, MPI_INTEGER4, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_INT4_V

  subroutine PAR_FILE_WRITE_ALL_INT4_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    integer(4), pointer                            :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_4, 1_4)
       call MPI_FILE_WRITE_ALL(fh, dummy, 0, MPI_INTEGER4, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_INT4_M

  subroutine PAR_FILE_WRITE_ALL_INT8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    integer(8), pointer                            :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_8)
       call MPI_FILE_WRITE_ALL(fh, dummy, 0, MPI_INTEGER8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_INT8_V

  subroutine PAR_FILE_WRITE_ALL_INT8_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    integer(8), pointer                            :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_8, 1_8)
       call MPI_FILE_WRITE_ALL(fh, dummy, 0, MPI_INTEGER8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_INT8_M

  subroutine PAR_FILE_WRITE_ALL_REAL8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    real(8), pointer                               :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_ip)
       call MPI_FILE_WRITE_ALL(fh, dummy, 0, MPI_REAL8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL(fh, buf, nsize, MPI_REAL8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_REAL8_V

  subroutine PAR_FILE_WRITE_ALL_REAL8_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    real(8), pointer                               :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    if (nsize==0 .or. .not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_ip, 1_ip)
       call MPI_FILE_WRITE_ALL(fh, dummy, 0, MPI_REAL8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL(fh, buf, nsize, MPI_REAL8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_REAL8_M

  subroutine PAR_FILE_WRITE_ALL_REAL4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    real(4), pointer                               :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_REAL4_V

  subroutine PAR_FILE_WRITE_ALL_REAL4_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    real(4), pointer                               :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_REAL4_M

#ifndef __PGI
  subroutine PAR_FILE_WRITE_ALL_REAL16_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    real(16), pointer                               :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_REAL16_V

  subroutine PAR_FILE_WRITE_ALL_REAL16_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    real(16), pointer                               :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_REAL16_M
#endif

  subroutine PAR_FILE_WRITE_ALL_BEGIN_INT4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    integer(4), pointer                            :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    if (.not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_4)
       call MPI_FILE_WRITE_ALL_BEGIN(fh, dummy, 0, MPI_INTEGER4, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL_BEGIN(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_INT4_V

  subroutine PAR_FILE_WRITE_ALL_BEGIN_INT4_M(fh, buf, nsize)
    integer(ip),               intent(in)          :: fh
    integer(4), pointer,       intent(in)          :: buf(:,:)
    integer(ip),               intent(in)          :: nsize
    integer(4), pointer                            :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    if (.not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_4, 1_4)
       call MPI_FILE_WRITE_ALL_BEGIN(fh, dummy, 0, MPI_INTEGER4, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL_BEGIN(fh, buf, nsize, MPI_INTEGER4, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_INT4_M

  subroutine PAR_FILE_WRITE_ALL_BEGIN_INT8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    integer(8), pointer                            :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    if (.not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_8)
       call MPI_FILE_WRITE_ALL_BEGIN(fh, dummy, 0, MPI_INTEGER8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL_BEGIN(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_INT8_V

  subroutine PAR_FILE_WRITE_ALL_BEGIN_INT8_M(fh, buf, nsize)
    integer(ip),               intent(in)          :: fh
    integer(8), pointer,       intent(in)          :: buf(:,:)
    integer(ip),               intent(in)          :: nsize
    integer(8), pointer                            :: dummy(:,:)
    nullify(dummy)
#ifndef MPI_OFF
    if (.not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_8, 1_8)
       call MPI_FILE_WRITE_ALL_BEGIN(fh, dummy, 0, MPI_INTEGER8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL_BEGIN(fh, buf, nsize, MPI_INTEGER8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_INT8_M

  subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL8_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    real(8), pointer                               :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    if (.not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_ip)
       call MPI_FILE_WRITE_ALL_BEGIN(fh, dummy, 0, MPI_REAL8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL_BEGIN(fh, buf, nsize, MPI_REAL8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL8_V

  subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL8_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    real(8), pointer                               :: dummy(:,:)
#ifndef MPI_OFF
    if (.not.associated(buf)) then
       call memory_alloca(mpio_memor,'dummy',vacal,  dummy, 1_ip, 1_ip)
       call MPI_FILE_WRITE_ALL_BEGIN(fh, dummy, 0, MPI_REAL8, status, ierr)
       call memory_deallo(mpio_memor,'dummy',vacal,  dummy)
    else
       call MPI_FILE_WRITE_ALL_BEGIN(fh, buf, nsize, MPI_REAL8, status, ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL8_M

  subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL4_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    real(4), pointer                               :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL4_V

  subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL4_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    real(4), pointer                               :: dummy(:,:)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL4_M

#ifndef __PGI
  subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL16_V(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:)
    integer(ip),           intent(in)              :: nsize
    real(16), pointer                               :: dummy(:)
    nullify(dummy)
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL16_V

  subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL16_M(fh, buf, nsize)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:,:)
    integer(ip),           intent(in)              :: nsize
    real(16), pointer                               :: dummy(:,:)
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_BEGIN_REAL16_M
#endif

  subroutine PAR_FILE_WRITE_ALL_END_INT4_V(fh, buf)
    integer(ip),           intent(in)              :: fh
    integer(4), pointer,   intent(in)              :: buf(:)
#ifndef MPI_OFF
    call MPI_FILE_WRITE_ALL_END(fh, buf, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_INT4_V

  subroutine PAR_FILE_WRITE_ALL_END_INT4_M(fh, buf)
    integer(ip),            intent(in)             :: fh
    integer(4), pointer,   intent(in)              :: buf(:,:)
#ifndef MPI_OFF
    call MPI_FILE_WRITE_ALL_END(fh, buf, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_INT4_M


  subroutine PAR_FILE_WRITE_ALL_END_INT8_V(fh, buf)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(in)              :: buf(:)
#ifndef MPI_OFF
    call MPI_FILE_WRITE_ALL_END(fh, buf, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_INT8_V

  subroutine PAR_FILE_WRITE_ALL_END_INT8_M(fh, buf)
    integer(ip),           intent(in)              :: fh
    integer(8), pointer,   intent(in)              :: buf(:,:)
#ifndef MPI_OFF
    call MPI_FILE_WRITE_ALL_END(fh, buf, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_INT8_M

  subroutine PAR_FILE_WRITE_ALL_END_REAL8_V(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:)
#ifndef MPI_OFF
    call MPI_FILE_WRITE_ALL_END(fh, buf, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_REAL8_V

  subroutine PAR_FILE_WRITE_ALL_END_REAL8_M(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(8), pointer,      intent(in)              :: buf(:,:)
#ifndef MPI_OFF
    call MPI_FILE_WRITE_ALL_END(fh, buf, status, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be written. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_REAL8_M

  subroutine PAR_FILE_WRITE_ALL_END_REAL4_V(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_REAL4_V

  subroutine PAR_FILE_WRITE_ALL_END_REAL4_M(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(4), pointer,      intent(in)              :: buf(:,:)
#ifndef MPI_OFF
    call runend("MPIO: RP 4 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_REAL4_M

#ifndef __PGI
  subroutine PAR_FILE_WRITE_ALL_END_REAL16_V(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:)
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_REAL16_V

  subroutine PAR_FILE_WRITE_ALL_END_REAL16_M(fh, buf)
    integer(ip),           intent(in)              :: fh
    real(16), pointer,      intent(in)              :: buf(:,:)
#ifndef MPI_OFF
    call runend("MPIO: RP 16 IS NOT IMPLEMENTED")
#endif
  end subroutine PAR_FILE_WRITE_ALL_END_REAL16_M
#endif

  subroutine PAR_FILE_SEEK_SET(fh, offset)
    integer(ip),           intent(in)              :: fh
    integer(MPI_OFFSET_KIND),           intent(in) :: offset
#ifndef MPI_OFF
    ierr=0
    call MPI_FILE_SEEK(fh, int(offset,MPI_OFFSET_KIND), MPI_SEEK_SET, ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" offset could not be modified (File seek). MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SEEK_SET

  subroutine PAR_FILE_SET_SIZE(fh, nsize)
    integer(ip),           intent(in)              :: fh
    integer(MPI_OFFSET_KIND),           intent(in) :: nsize
#ifndef MPI_OFF
    ierr=0
    call MPI_FILE_SET_SIZE(fh, int(nsize,MPI_OFFSET_KIND), ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" could not be resized. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_SET_SIZE

  subroutine PAR_FILE_GET_SIZE(fh, nsize)
    integer(ip),           intent(in)              :: fh
    integer(MPI_OFFSET_KIND),           intent(inout):: nsize
#ifndef MPI_OFF
    ierr=0
    call MPI_FILE_GET_SIZE(fh, int(nsize,MPI_OFFSET_KIND), ierr)
    if (ierr/=MPI_SUCCESS) then
       call runend("File descriptor "//trim(intost(fh))//" size unavailable. MPI error code: "//trim(intost(ierr)))
    end if
#endif
  end subroutine PAR_FILE_GET_SIZE

  subroutine par_my_com(PAR_COMM_TO_USE, commu, wherein)
    character(*),          optional, intent(in)     :: wherein
    type(comm_data_par),   pointer, intent(inout)   :: commu
    integer(4),            intent(inout)            :: PAR_COMM_TO_USE
#ifndef MPI_OFF
    if(present(wherein)) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR(wherein_p,PAR_COMM_TO_USE,commu)
    end if
#endif
  end subroutine par_my_com



end module mod_mpio_par_mpiwrapper
