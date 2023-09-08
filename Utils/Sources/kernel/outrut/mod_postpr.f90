!-----------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @name    Toolbox for pstprocess
!> @file    mod_postpr.f90
!> @author  houzeaux
!> @date    2018-02-19
!> @brief   postprocess
!> @details postprocess
!> @{
!-----------------------------------------------------------------------

module mod_postpr

  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  use def_postpr
  use mod_mpio_seq_log
  use mod_mpio_postpr
  use mod_permut
  use def_inpout
  use def_mpio
  use mod_outfor,   only : outfor
  use mod_messages, only : livinf
  use mod_iofile,   only : iofile_flush_unit

  use mod_communications
  implicit none

  private
  !
  ! Variables fo headers
  !
  integer(ip), parameter :: nwwww=10_ip
  character(5)           :: wwwww(nwwww)
  character(8)           :: wwww8(nwwww)
  integer(4)             :: iiiii(10)=0_4
  real(8)                :: rrrrr(10)=0.0_8
  !
  ! vtk variable
  !
  integer(ip)            :: success
  integer(ip)            :: mastervtk = 666_ip,opart

  character(5)           :: bench_rw
  integer(ip)            :: kfl_perme
  integer(ip)            :: kfl_permb

  integer(ip)            :: bridge_ip(1)
  real(rp)               :: bridge_rp(1)

  interface postpr
     module procedure &
          postpr_real_scalar,postpr_real_vector,postpr_int_scalar,postpr_int_vector,&
          posr3p,posvex
  end interface postpr

  interface postpr_right_now
     module procedure &
          postpr_right_now_i1,&
          postpr_right_now_r1,&
          postpr_right_now_r2,&
          postpr_right_now_r0
  end interface postpr_right_now

  public :: postpr                        ! Postprocess a variable
  public :: postpr_right_now              ! Postprocess a variable on the fly
  public :: postpr_right_now_r0           ! Postprocess a variable on the fly
  public :: postpr_at_current_time_step   ! Postprocess at current step
  public :: postpr_initialization         ! Initialize postprocess variables
  public :: postpr_real_vector
  public :: postpr_real_scalar

contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Initialize module variables
  !> @details Initialize variables and nullify arrays of this module
  !>
  !----------------------------------------------------------------------

  subroutine postpr_initialization()

    iiiii = 0_4
    rrrrr = 0.0_8
    wwwww = ''
    wwww8 = ''
    nullify(gescp)
    nullify(gevep)
    nullify(giscp)
    nullify(givep)

  end subroutine postpr_initialization

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess "sobre la marcha"
  !> @details Postprocess a variable on the fly
  !>
  !----------------------------------------------------------------------

  subroutine postpr_right_now_i1(wopo1,wopo2,wopo3,iscalar)
    character(5),         intent(in) :: wopo1,wopo2,wopo3
    integer(ip), pointer, intent(in) :: iscalar(:)
    character(5)                     :: wopos_loc(3)
    integer(ip)                      :: ipoin,nbound

    if( INOTMASTER ) then
       if( wopo3 == 'NELEM') then
          nbound = nelem
       else
          nbound = npoin
       end if
       call postpr_memory(0_ip,nbound,0_ip)
       do ipoin = 1,nbound
          gescp(ipoin) = real(iscalar(ipoin),rp)
       end do
    end if
    wopos_loc(1) = wopo1
    wopos_loc(2) = wopo2
    wopos_loc(3) = wopo3
    call postpr(gescp,wopos_loc,0_ip,0.0_rp)
    if( INOTMASTER ) call postpr_memory(2_ip,nbound,0_ip)

  end subroutine postpr_right_now_i1

  subroutine postpr_right_now_r1(wopo1,wopo2,wopo3,rvector,pleng_opt)
    character(5),         intent(in)             :: wopo1,wopo2,wopo3
    real(rp),    pointer, intent(inout)          :: rvector(:)
    integer(ip),          intent(in),   optional :: pleng_opt
    character(5)                                 :: wopos_loc(3)
    integer(ip)                                  :: ipoin,nbound

    wopos_loc(1) = wopo1
    wopos_loc(2) = wopo2
    wopos_loc(3) = wopo3
    call postpr(rvector,wopos_loc,0_ip,0.0_rp,pleng_opt=pleng_opt)

  end subroutine postpr_right_now_r1

  subroutine postpr_right_now_r0(wopo1,wopo2,wopo3,kpoin,rvector,TIME_STEP_OUTPUT,TIME_OUTPUT)
    integer(ip),  intent(in)             :: kpoin
    character(5), intent(in)             :: wopo1,wopo2,wopo3
    real(rp),     intent(inout)          :: rvector(*)
    real(rp),     intent(in),   optional :: TIME_OUTPUT
    integer(ip),  intent(in),   optional :: TIME_STEP_OUTPUT
    character(5)                         :: wopos_loc(3)
    integer(ip)                          :: ipoin,nbound
    real(rp)                             :: my_time
    integer(ip)                          :: my_time_step
    real(rp),     pointer                :: generic_scalar(:)

    nullify(generic_scalar)

    wopos_loc(1) = wopo1
    wopos_loc(2) = wopo2
    wopos_loc(3) = wopo3
    if( INOTMASTER ) then
       allocate(generic_scalar(npoin))
       do ipoin = 1,npoin
          generic_scalar(ipoin) = rvector(ipoin)
       end do
    end if
    if( present(TIME_OUTPUT) ) then
       my_time = TIME_OUTPUT
    else
       my_time = 0.0_rp
    end if
    if( present(TIME_STEP_OUTPUT) ) then
       my_time_step = TIME_STEP_OUTPUT
    else
       my_time_step = 0
    end if

    call postpr(generic_scalar,wopos_loc,my_time_step,my_time,TAG1=TIME_STEP_OUTPUT)

    if( INOTMASTER ) deallocate(generic_scalar)

  end subroutine postpr_right_now_r0

  subroutine postpr_right_now_r2(wopo1,wopo2,wopo3,rvector)
    character(5),         intent(in) :: wopo1,wopo2,wopo3
    real(rp),    pointer, intent(in) :: rvector(:,:)
    character(5)                     :: wopos_loc(3)
    integer(ip)                      :: ipoin,nbound
    real(rp),    pointer             :: my_vect(:,:)

    nullify(my_vect)
    if( INOTMASTER ) then
       if( wopo3 == 'NELEM') then
          nbound = nelem
       else
          nbound = npoin
       end if
       allocate(my_vect(ndime,nbound))
       do ipoin = 1,nbound
          my_vect(1:ndime,ipoin) = rvector(1:ndime,ipoin)
       end do
    end if
    wopos_loc(1) = wopo1
    wopos_loc(2) = wopo2
    wopos_loc(3) = wopo3
    call postpr(my_vect,wopos_loc,0_ip,0.0_rp)
    if( INOTMASTER ) deallocate(my_vect)

  end subroutine postpr_right_now_r2

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess current time step
  !> @details Postprocess current time step
  !>
  !----------------------------------------------------------------------

  function postpr_at_current_time_step()
    integer(ip) :: postpr_at_current_time_step
    integer(ip) :: imodu,ivari,itime,iok
    logical(lg) :: if_already

    postpr_at_current_time_step = 0

    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 ) then
          do ivari = 1,nvarp

             if_already = .false.

             if( ittyp == ITASK_ENDRUN ) then
                iok = 1
             else
                iok = 0
             end if

             if( &
                  ittim >= momod(imodu) % postp(1) % npp_inits .and. &
                  iok == 0 .and. &
                  momod(imodu) % postp(1) % npp_stepi(ivari) > 0 ) then
                if( mod(ittim, momod(imodu) % postp(1) % npp_stepi(ivari)) == 0 ) then
                   postpr_at_current_time_step = postpr_at_current_time_step + 1
                   if_already = .true.
                end if
             end if
             !
             ! At a given time
             !
             if( iok == 0 .and. ( .not. if_already ) ) then
                do itime = 1,nvart
                   if(   abs( momod(imodu) % postp(1) % pos_times(itime,ivari)-cutim) < (0.5_rp*dtime) .and. &
                        &     momod(imodu) % postp(1) % pos_times(itime,ivari)        > 0.0_rp) then
                      postpr_at_current_time_step = postpr_at_current_time_step + 1
                      if_already = .true.
                   end if
                end do
             end if
             !
             ! At a given time period
             !
             if( cutim >= momod(imodu) % postp(1) % pos_tinit .and. iok == 0 &
                  .and. ( .not. if_already ) ) then
                if(    abs(momod(imodu) % postp(1) % pos_times(1,ivari)-cutim) < (0.6_rp*dtime).and.&
                     &     momod(imodu) % postp(1) % pos_perio(ivari)          > 0.0_rp) then
                   postpr_at_current_time_step = postpr_at_current_time_step + 1
                   if_already = .true.
                end if
             end if

          end do
       end if
    end do
  end function postpr_at_current_time_step

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of real scalar variables
  !> @details Bridge to the postprocess of a real variables
  !>
  !----------------------------------------------------------------------

  subroutine postpr_real_scalar(bridge,wopos,itste,ttime,worig,pleng_opt,TAG1,TAG2)

    implicit none
    character(*), intent(in)             :: wopos(*)
    real(rp),     intent(inout), pointer :: bridge(:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip)                          :: pdime
    real(rp)                             :: dummr(2)

    pdime = 1
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if
    !
    !
    !
    if (postpr_is_mpio(wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)) then
       gescar_mpio    => bridge

       call posmpio_real_v()


    else if( kfl_outfo == 50 ) then
       !
       ! HDFPOS
       !
       wopos_hdf(1) =  wopos(1)
       wopos_hdf(2) =  wopos(2)
       gesca_hdf    => bridge
       call Hdfpos(1_ip)
       !
       ! VTK ( only write postpr not restart, restart is alya)
       !
    else if(( kfl_outfo == 40 .OR. kfl_outfo == 41).AND. kfl_reawr == 0) then
       if( IMASTER ) then
          call geovtk(dummr,wopos,itste,ttime,pdime)
       else if ( ISLAVE ) then
          call geovtk(bridge,wopos,itste,ttime,pdime)
       else if( ISEQUEN ) then
          call geovtk(bridge,wopos,itste,ttime,pdime)
       endif
    else
       !
       ! Alya
       !
       if ((kfl_reawr/=1)) then
          bench_rw='WRITE'
       else
          bench_rw='READ '
       end if
       call start_timer(wopos(1), bench_rw, 'ALYA', 1_ip, 1_ip)
       if( IMASTER ) then
          call postpr_postprocess_rp(bridge_rp,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)
       else
          if( associated(bridge) ) then
             call postpr_postprocess_rp(bridge,    wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)
          else
             call postpr_postprocess_rp(bridge_rp,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)
          end if
       end if
       call end_timer()

    end if

  end subroutine postpr_real_scalar

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of real vector variables
  !> @details Bridge to the postprocess of a real variables
  !>          When kdime is present, kdime rules over ndime.
  !>
  !----------------------------------------------------------------------

  subroutine postpr_real_vector(bridge,wopos,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2)

    implicit none
    character(*), intent(in)             :: wopos(*)
    real(rp),     intent(inout), pointer :: bridge(:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip)                          :: pdime
    real(rp)                             :: dummr(2,2)

    if( present(kdime) ) then
       pdime = kdime
    else
       pdime = ndime
    end if
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if

    if (postpr_is_mpio(wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)) then
       gevecr_mpio    => bridge

       call posmpio_real_m()

    else if( kfl_outfo == 50 ) then
       !
       ! HDFPOS
       !
       wopos_hdf(1) =  wopos(1)
       wopos_hdf(2) =  wopos(2)
       gevec_hdf    => bridge
       !
       ! VORTEX EXTRACTION HDF5
       !
       if( wopos(3) == 'VORTX' ) then
          call Hdfpos(11_ip)
       else
          call Hdfpos(2_ip)
       end if
       !
       ! VTK ( only write postpr not restart, restart is alya format)
       !
    else if(( kfl_outfo == 40 .OR. kfl_outfo == 41).AND. kfl_reawr == 0) then
       if( IMASTER ) then
          call geovtk(dummr,wopos,itste,ttime,pdime)
       else if ( ISLAVE ) then
          call geovtk(bridge,wopos,itste,ttime,pdime)
       else if( ISEQUEN ) then
          call geovtk(bridge,wopos,itste,ttime,pdime)
       endif
    else
       !
       ! Alya
       !
       if ((kfl_reawr/=1)) then
          bench_rw='WRITE'
       else
          bench_rw='READ '
       end if
       call start_timer(wopos(1), bench_rw, 'ALYA', 1_ip, 1_ip)
       if( IMASTER ) then
          call postpr_postprocess_rp(bridge_rp,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)
       else
          if( associated(bridge) ) then
             call postpr_postprocess_rp(bridge,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)
          else
             call postpr_postprocess_rp(bridge_rp,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)
          end if
       end if
       call end_timer()

    end if

  end subroutine postpr_real_vector

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of integer scalar variables
  !> @details Bridge to the postprocess of a integer variables
  !>
  !----------------------------------------------------------------------

  subroutine postpr_int_scalar(bridge,wopos,itste,ttime,worig,TAG1,TAG2)

    implicit none
    character(*), intent(in)             :: wopos(*)
    integer(ip),  intent(inout), pointer :: bridge(:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip)                          :: pdime
    integer(ip)                          :: dummi(2)

    pdime = 1
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if

    if (postpr_is_mpio(wopos,itste,ttime,pdime,0_ip,TAG1,TAG2)) then
       gescai_mpio    => bridge

       call posmpio_int_v()

    else if( kfl_outfo == 50 ) then
       !
       ! HDFPOS
       !
       call runend('NOT CODED')
    else if( kfl_outfo == 40 .AND. kfl_outfo == 41 ) then
       !
       ! VTK
       !
       call runend('NOT CODED')
    else
       !
       ! Alya
       !
       if ((kfl_reawr/=1)) then
          bench_rw='WRITE'
       else
          bench_rw='READ '
       end if
       call start_timer(wopos(1), bench_rw, 'ALYA', 1_ip, 1_ip)
       if( IMASTER ) then
          call postpr_postprocess_ip(bridge_ip,wopos,itste,ttime,pdime,TAG1,TAG2)
       else
          if( associated(bridge) ) then
             call postpr_postprocess_ip(bridge,wopos,itste,ttime,pdime,TAG1,TAG2)
          else
             call postpr_postprocess_ip(bridge_ip,wopos,itste,ttime,pdime,TAG1,TAG2)
          end if
       end if
       call end_timer()
    end if

#ifdef EVENT
    call mpitrace_user_function(0)
#endif

  end subroutine postpr_int_scalar

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of integer vector variables
  !> @details Bridge to the postprocess of a integer variables
  !>          When kdime is present, kdime rules over ndime.
  !>
  !----------------------------------------------------------------------

  subroutine postpr_int_vector(bridge,wopos,itste,ttime,kdime,worig,TAG1,TAG2)

    implicit none
    character(*), intent(in)             :: wopos(*)
    integer(ip),  intent(inout), pointer :: bridge(:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip)                          :: pdime
    integer(ip)                          :: dummi(2,2)

    if( present(kdime) ) then
       pdime = kdime
    else
       pdime = ndime
    end if
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if
    if (postpr_is_mpio(wopos,itste,ttime,pdime,0_ip,TAG1,TAG2)) then
       geveci_mpio    => bridge

       call posmpio_int_m()

    else if( kfl_outfo == 50 ) then
       !
       ! HDFPOS
       !
       call runend('NOT CODED')
    else if( kfl_outfo == 40 .AND. kfl_outfo == 41 ) then
       !
       ! VTK
       !
       call runend('NOT CODED')
    else
       !
       ! Alya
       !
       if ((kfl_reawr/=1)) then
          bench_rw='WRITE'
       else
          bench_rw='READ '
       end if
       call start_timer(wopos(1), bench_rw, 'ALYA', 1_ip, 1_ip)
       if( IMASTER ) then
          call postpr_postprocess_ip(bridge_ip,wopos,itste,ttime,pdime,TAG1,TAG2)
       else
          if( associated(bridge) ) then
             call postpr_postprocess_ip(bridge,wopos,itste,ttime,pdime,TAG1,TAG2)
          else
             call postpr_postprocess_ip(bridge_ip,wopos,itste,ttime,pdime,TAG1,TAG2)
          end if
       end if
       call end_timer()
    end if

  end subroutine postpr_int_vector

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of r3p type
  !> @details Bridge to the postprocess of a r3p types
  !>          When kdime is present, kdime rules over ndime.
  !>
  !----------------------------------------------------------------------

  subroutine posr3p(bridge,wopos,itste,ttime,kdime)

    implicit none
    character(*), intent(in)              :: wopos(*)
    type(r3p),    intent(inout), pointer  :: bridge(:)
    integer(ip),  intent(in)              :: itste
    real(rp),     intent(in)              :: ttime
    integer(ip),  intent(in),    optional :: kdime
    integer(ip)                           :: ileng,pleng
    integer(ip)                           :: idim1,idim2,ptota,ielty,ii
    integer(ip)                           :: pdim1,pdim2,itota,pdime

    if( present(kdime) ) then
       pdime = kdime
    else
       pdime = 1
    end if

    if ( kfl_outfo == 50 ) then
       !
       ! HDFPOS
       !
       wopos_hdf(1) =  wopos(1)
       wopos_hdf(2) =  wopos(2)
       ger3p_hdf   => bridge

       call Hdfpos(3_ip)
    else
       !
       ! Type of postprocess
       !
       wwwww(1) = 'ALYA '
       wwwww(2) = 'V0003'
       wwwww(3) = wopos(1)
       wwwww(4) = wopos(2)
       wwwww(5) = 'R3P  '
       wwwww(6) = 'REAL '
       wwwww(7) = '8BYTE'
       if( wopos(2) == 'R3PVE' ) then
          iiiii(1) = int(ndime,4)
       else
          iiiii(1) = 1_4
       end if
       iiiii(2) = 0_4
       iiiii(4) = int(itste,4)
       rrrrr(1) = ttime
       wopos_pos(1) = wopos(1)
       parii    = NELEM_TYPE
       pleng    = nelem
       if( IMASTER ) pleng = nelem_total

       iiiii(2) = int(pleng,4)
       !
       ! Output solution
       !
       call postpr_open_file()
       !call opfpos(1_ip)
       if( kfl_reawr < 0 ) then
          kfl_reawr = abs(kfl_reawr)
          return
       end if
 
       call postpr_header()

       if( ISEQUEN ) then

          if( kfl_reawr /= 1 ) then
             ptota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                ptota = ptota + pdim1 * pdim2
             end do
             write(lun_postp) ( ngaus(ielty),ielty=iesta_dom,iesto_dom)
             write(lun_postp) pleng
             write(lun_postp) ptota
             write(lun_postp) (((bridge(ileng)%a(idim1,idim2,pdime),idim1=1,size(bridge(ileng)%a,1,kind=ip)),&
                  idim2=1,size(bridge(ileng)%a,2,kind=ip)),ileng=1,pleng)
          else
             read(lun_postp)  ( ii,ielty=iesta_dom,iesto_dom)
             read(lun_postp)  pleng
             read(lun_postp)  ptota
             read(lun_postp)  (((bridge(ileng)%a(idim1,idim2,pdime),idim1=1,size(bridge(ileng)%a,1,kind=ip)),&
                  idim2=1,size(bridge(ileng)%a,2,kind=ip)),ileng=1,pleng)
          end if

       else if( IMASTER ) then

          if( kfl_reawr /= 1 ) then
             write(lun_postp) ( ngaus(ielty),ielty=iesta_dom,iesto_dom)
          else if( kfl_reawr == 1 ) then
             read(lun_postp) ( ii,ielty=iesta_dom,iesto_dom)
          end if

          do kfl_desti_par = 1,npart

             if( kfl_reawr /= 1 ) then
                !
                ! - Master writes postpro without filter
                ! - Master writes restart
                !
                pleng = nelem_par(kfl_desti_par)
                call pararr('RCV',0_ip,1_ip,ptota)
                write(lun_postp) pleng
                write(lun_postp) ptota
                call postpr_memory(0_ip,ptota,0_ip)
                call pararr('RCV',0_ip,ptota,gescp)
                write(lun_postp) ( gescp(ileng), ileng = 1,ptota )
                call postpr_memory(2_ip,ptota,0_ip)

             else if( kfl_reawr == 1 ) then
                !
                ! - Master reads restart and sends to slave
                !
                pleng = nelem_par(kfl_desti_par)
                read(lun_postp) pleng
                read(lun_postp) ptota
                call postpr_memory(0_ip,ptota,0_ip)
                read(lun_postp) ( gescp(ileng), ileng = 1,ptota )
                call pararr('SND',0_ip,ptota,gescp)
                call postpr_memory(2_ip,ptota,0_ip)

             end if

          end do

       else if( ISLAVE ) then

          kfl_desti_par = 0

          if( kfl_reawr /= 1 ) then
             !
             ! Slaves without filter
             !
             ptota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                ptota = ptota + pdim1 * pdim2
             end do
             call parari('SND',0_ip,1_ip,ptota)
             call postpr_memory(0_ip,ptota,0_ip)
             itota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                do idim2=1,pdim2
                   do idim1=1,pdim1
                      itota = itota + 1
                      gescp(itota) = bridge(ileng)%a(idim1,idim2,pdime)
                   end do
                end do
             end do
             call pararr('SND',0_ip,ptota,gescp)
             call postpr_memory(2_ip,ptota,0_ip)

          else if( kfl_reawr == 1 ) then
             !
             ! Slaves receive
             !
             ptota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                ptota = ptota + pdim1 * pdim2
             end do
             call postpr_memory(0_ip,ptota,0_ip)
             call pararr('RCV',0_ip,ptota,gescp)
             itota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                do idim2=1,pdim2
                   do idim1=1,pdim1
                      itota = itota + 1
                      bridge(ileng)%a(idim1,idim2,pdime) = gescp(itota)
                   end do
                end do
             end do
             call postpr_memory(2_ip,ptota,0_ip)

          end if

       end if

       call postpr_close_file()
       !call opfpos(2_ip)

    end if

#ifdef EVENT
    call mpitrace_user_function(0)
#endif

  end subroutine posr3p

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of complex vector variables
  !> @details Bridge to the postprocess of a complex variables
  !>          When kdime is present, kdime rules over ndime.
  !>
  !----------------------------------------------------------------------

  subroutine posvex(bridge,wopos,itste,ttime,kdime,worig)

    implicit none
    character(*), intent(in)             :: wopos(*)
    complex(rp),  intent(inout), pointer :: bridge(:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip)                          :: pdime

    if( present(kdime) ) then
       pdime = kdime
    else
       pdime = ndime
    end if
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if

    if( kfl_outfo == 50 ) then
       !
       ! HDFPOS
       !
       call runend('POSVEX: NOT CODED')
       !wopos_hdf(1) =  wopos(1)
       !wopos_hdf(2) =  wopos(2)
       !gevec_hdf    => bridge
       !if( wopos(3) == 'VORTX' ) then
       !   call Hdfpos(11_ip)
       !else
       !   call Hdfpos(2_ip)
       !end if

    else
       !
       ! Alya
       !
       !if( IMASTER ) then
       !   call posrex( dummr,wopos,itste,ttime,pdime)
       !else
       !   call posrex(bridge,wopos,itste,ttime,pdime)
       !end if

    end if

  end subroutine posvex

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess a real array
  !> @details Postprocess a real array
  !>
  !----------------------------------------------------------------------

  subroutine postpr_postprocess_rp(bridge,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)

    use mod_communications , only : PAR_SEND_RECEIVE
    use mod_communications , only : PAR_SUM

    character(*), intent(in)           :: wopos(*)
    real(rp),     intent(inout)        :: bridge(*)
    integer(ip),  intent(in)           :: itste
    real(rp),     intent(in)           :: ttime
    integer(ip),  intent(in)           :: pdime
    integer(ip),  intent(in), optional :: pleng_opt
    integer(ip),  intent(in), optional :: TAG1
    integer(ip),  intent(in), optional :: TAG2
    integer(ip)                        :: ileng,idime,pleng,kpoin,ipoin,dummi(2)
    integer(ip)                        :: imesh,idofn,qpoin,pdofn,pleng_total
    integer(ip)                        :: pleng_com(2),tag1_opt,tag2_opt

    !
    ! Options
    !
    tag1_opt = 0
    tag2_opt = 0
    if( present(TAG1) ) tag1_opt = TAG1
    if( present(TAG2) ) tag2_opt = TAG2
    !
    ! Type of postprocess
    !
    wwwww(1) = 'ALYA '
    wwwww(2) = 'V0003'
    wwwww(3) = wopos(1)
    wwwww(4) = wopos(2)
    wwwww(5) = wopos(3)
    iiiii(5) = int(tag1_opt,4)
    iiiii(6) = int(tag2_opt,4)
    wwwww(6) = 'REAL '
    wwwww(7) = '8BYTE'
    iiiii(1) = int(pdime,4)
    iiiii(2) = 0_4
    iiiii(4) = int(itste,4)
    rrrrr(1) = ttime

    wopos_pos(1) = wopos(1)

    if( kfl_reawr == 0 ) then
       imesh = kfl_posdi
    else
       imesh = ndivi
    end if
    !
    ! Map the result onto original mesh:
    ! - Only for NPOIN type results
    !
    kfl_permu = 0
    kfl_perme = 0
    kfl_permb = 0

    if( kfl_origi == 1 .and. ( kfl_posdi /= ndivi ) .and. wopos(3) == 'NPOIN' ) kfl_permu = 1
    if( kfl_origi == 1 .and. ( kfl_posdi /= ndivi ) .and. wopos(3) == 'NELEM' ) kfl_perme = 1
    if( kfl_origi == 1 .and. ( kfl_posdi /= ndivi ) .and. wopos(3) == 'NBOUN' ) kfl_permb = 1

    if( wopos(3) == 'NPOIN' ) then

       parii = NPOIN_TYPE
       if( associated(meshe) .and. size(meshe,kind=ip) >= imesh ) then
          pleng = meshe(imesh) % npoin
          if( IMASTER ) pleng = meshe(imesh) % npoin_total
       else
          pleng = npoin
          if( IMASTER ) pleng = npoin_total
       end if

    else if( wopos(3) == 'NELEM' ) then

       parii = NELEM_TYPE
       pleng = meshe(imesh) % nelem
       if( IMASTER ) pleng = meshe(imesh) % nelem_total

    else if( wopos(3) == 'NBOUN' ) then

       parii = NBOUN_TYPE
       pleng = meshe(imesh) % nboun
       if( IMASTER ) pleng = meshe(imesh) % nboun_total

    else if( wopos(3) == 'VORTX' ) then
       call vortwr()
       ! write that file in vu format : call vortvu()
       return

    else if( wopos(3) == 'WHATE' ) then

       parii = NBOUN_TYPE
       if( present(pleng_opt) ) then
          pleng = pleng_opt
       else
          call runend('MOD_POSTPR - POSTPR_POSTPROCESS_RP: A SIZE SHOULD BE GIVEN')
       end if
       if( IMASTER ) then
          pleng_total = 0
       else
          pleng_total = pleng
       end if
       call PAR_SUM(pleng_total)
       if( IMASTER ) pleng = pleng_total

    else
       call runend('MOD_POSTPR - POSTPR_POSTPROCESS_RP: UNDEFINED POSTPROCESS')
    end if
    iiiii(2) = int(pleng,4)
    !
    ! Get filter
    !
    call fildef(3_ip)
    !
    ! Output solution
    !
    call postpr_open_file(TAG1,TAG2)

    if( kfl_reawr < 0 ) then
       kfl_reawr = abs(kfl_reawr)
       return
    end if
    !
    ! Write header
    !
    call postpr_header()

    if( ISEQUEN ) then

       !-----------------------------------------------------------------
       !
       ! SEQUENTIAL
       !
       !-----------------------------------------------------------------

       if( kfl_reawr == 0 .and. kfl_filte == 0 ) then
          !
          ! Postprocess without filter
          !
          write(lun_postp) pleng
          write(lun_postp) ( bridge(ileng), ileng = 1,pdime*pleng )

       else if( kfl_reawr == 0 .and. kfl_filte /= 0 ) then
          !
          ! Postprocess with filter
          !
          call postpr_memory(0_ip,pdime*meshe(imesh) % npoin,0_ip)
          call postpr_memory(1_ip,meshe(imesh) % npoin,0_ip)

          pleng = 0
          do qpoin = 1,meshe(imesh) % npoin

             if( kfl_permu == 1 ) then
                ipoin = lpmsh(qpoin)
             else
                ipoin = qpoin
             end if
             if( kfl_perme == 1 ) then
                ipoin = lemsh(qpoin)
             else
                ipoin = qpoin
             end if
             if( kfl_permb == 1 ) then
                ipoin = lbmsh(qpoin)
             else
                ipoin = qpoin
             end if

             if( gefil(ipoin) > 0 ) then
                pleng = pleng + 1
                giscp(pleng) = ipoin
                pdofn = (pleng-1) * pdime
                idofn = (ipoin-1) * pdime
                do idime = 1,pdime
                   pdofn = pdofn + 1
                   idofn = idofn + 1
                   gescp(pdofn) = bridge(idofn)
                end do
             end if

          end do

          write(lun_postp) pleng
          write(lun_postp) ( giscp(ileng), ileng = 1,pleng )
          write(lun_postp) ( gescp(ileng), ileng = 1,pdime*pleng )

          call postpr_memory(3_ip,meshe(imesh) % npoin,0_ip)
          call postpr_memory(2_ip,pdime*meshe(imesh) % npoin,0_ip)

       else if( kfl_reawr == 2 ) then
          !
          ! Preliminary
          !
          write(lun_postp) pleng
          write(lun_postp) ( bridge(ileng), ileng = 1,pdime*pleng )

       else if( kfl_reawr == 1 ) then
          !
          ! Restart
          !
          read(lun_postp) pleng
          read(lun_postp)  ( bridge(ileng), ileng = 1,pdime*pleng )

       end if

    else if( IMASTER ) then

       !-----------------------------------------------------------------
       !
       ! MASTER
       !
       !-----------------------------------------------------------------

       do kfl_desti_par = 1,npart

          if( kfl_filte == 0 .and. kfl_reawr /= 1 ) then
             !
             ! - Master writes postpro without filter
             ! - Master writes restart
             !
             if(      wopos(3) == 'NPOIN' ) then
                if( associated(meshe) .and. size(meshe,kind=ip) >= imesh ) then
                   pleng = meshe(imesh) % npoin_par(kfl_desti_par)
                else
                   pleng = npoin_par(kfl_desti_par)
                end if

             else if( wopos(3) == 'NELEM' ) then
                pleng = meshe(imesh) % nelem_par(kfl_desti_par)
             else if( wopos(3) == 'NBOUN' ) then
                pleng = meshe(imesh) % nboun_par(kfl_desti_par)
             else if( wopos(3) == 'WHATE' ) then
                call PAR_SEND_RECEIVE(0_ip,1_ip,dummi,pleng_com,'IN MY CODE',kfl_desti_par)
                pleng = pleng_com(1)
             end if

             call postpr_memory(0_ip,max(1_ip,pdime*pleng),0_ip)

             call pararr('RCV',0_ip,pdime*pleng,gescp)

             write(lun_postp) pleng
             write(lun_postp) ( gescp(ileng), ileng = 1,pdime*pleng )

             call postpr_memory(2_ip,max(1_ip,pdime*pleng),0_ip)

          else if( kfl_reawr == 1 ) then
             !
             ! - Master reads restart and sends to slave
             !
             if(      wopos(3) == 'NPOIN' ) then
                pleng = npoin_par(kfl_desti_par)
             else if( wopos(3) == 'NELEM' ) then
                pleng = nelem_par(kfl_desti_par)
             else if( wopos(3) == 'NBOUN' ) then
                pleng = nboun_par(kfl_desti_par)
             else if( wopos(3) == 'WHATE' ) then
                call PAR_SEND_RECEIVE(0_ip,1_ip,dummi,pleng_com,'IN MY CODE',kfl_desti_par)
                pleng = pleng_com(1)
             end if
             call postpr_memory(0_ip,max(1_ip,pdime*pleng),0_ip)
             read(lun_postp) pleng
             read(lun_postp) ( gescp(ileng), ileng = 1,pdime*pleng )
             call pararr('SND',0_ip,pdime*pleng,gescp)
             call postpr_memory(2_ip,max(1_ip,pdime*pleng),0_ip)

          else
             !
             ! - Master writes postpro with filter
             !
             call parari('RCV',0_ip,1_ip,pleng)
             if( pleng > 0 ) then
                call postpr_memory(1_ip,pleng,0_ip)
                call parari('RCV',0_ip,pleng,giscp)
                call postpr_memory(0_ip,pdime*pleng,0_ip)
                call pararr('RCV',0_ip,pdime*pleng,gescp)
                write(lun_postp) pleng
                write(lun_postp) ( giscp(ileng), ileng = 1,pleng )
                write(lun_postp) ( gescp(ileng), ileng = 1,pdime*pleng )
                call postpr_memory(2_ip,pdime*pleng,0_ip)
                call postpr_memory(3_ip,      pleng,0_ip)
             else
                write(lun_postp) pleng
             end if
          end if
       end do

    else if( ISLAVE ) then

       !-----------------------------------------------------------------
       !
       ! SLAVES
       !
       !-----------------------------------------------------------------

       kfl_desti_par = 0

       if( wopos(3) == 'WHATE' ) then
          pleng_com(1) = pleng
          call PAR_SEND_RECEIVE(1_ip,0_ip,pleng_com,dummi,'IN MY CODE',kfl_desti_par)
       end if

       if( kfl_reawr == 0 .and. kfl_filte == 0 ) then
          !
          ! Slaves without filter for postprocess
          !
          if(      kfl_permu == 1 ) then

             call postpr_memory(0_ip,pdime,pleng)
             call permut_rp2(pdime,pleng,lpmsh,bridge,gevep)
             call pararr('SND',0_ip,pdime*pleng,gevep)
             call postpr_memory(2_ip,pdime,pleng)

          else if( kfl_perme == 1 ) then

             call postpr_memory(0_ip,pdime,pleng)
             call permut_rp2(pdime,pleng,lemsh,bridge,gevep)
             call pararr('SND',0_ip,pdime*pleng,gevep)
             call postpr_memory(2_ip,pdime,pleng)

          else if( kfl_permb == 1 ) then

             call postpr_memory(0_ip,pdime,pleng)
             call permut_rp2(pdime,pleng,lbmsh,bridge,gevep)
             call pararr('SND',0_ip,pdime*pleng,gevep)
             call postpr_memory(2_ip,pdime,pleng)

          else

             call pararr('SND',0_ip,pdime*pleng,bridge)

          end if

       else if( kfl_reawr == 0 .and. kfl_filte /= 0 ) then
          !
          ! Slaves with filter for postprocess
          !
          call postpr_memory(0_ip,pdime,meshe(imesh) % npoin)
          call postpr_memory(1_ip,meshe(imesh) % npoin,0_ip)

          pleng = 0
          do qpoin = 1,meshe(imesh) % npoin
             if( kfl_permu == 1 ) then
                ipoin = lpmsh(qpoin)
             else if( kfl_perme == 1 ) then
                ipoin = lemsh(qpoin)
             else
                ipoin = qpoin
             end if
             if( gefil(ipoin) /= 0 ) then
                kpoin = (ipoin-1) * pdime
                pleng = pleng + 1
                giscp(pleng) = qpoin
                do idime = 1,pdime
                   kpoin = kpoin + 1
                   gevep(idime,pleng) = bridge(kpoin)
                end do
             end if
          end do

          call parari('SND',0_ip, 1_ip,pleng)
          if( pleng > 0 ) then
             call parari('SND',0_ip,pleng,giscp)
             call pararr('SND',0_ip,pdime*pleng,gevep)
          end if

          call postpr_memory(3_ip,meshe(imesh) % npoin,0_ip)
          call postpr_memory(2_ip,pdime,meshe(imesh) % npoin)

       else if( kfl_reawr == 1 ) then
          !
          ! Slaves receive for restart
          !
          call pararr('RCV',0_ip,pdime*pleng,bridge)

       else if( kfl_reawr == 2 ) then
          !
          ! Slaves send for preliminary
          !
          call pararr('SND',0_ip,pdime*pleng,bridge)

       end if

    end if

    call postpr_close_file()

  end subroutine postpr_postprocess_rp

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess an integer array
  !> @details Postprocess an integer array
  !>
  !----------------------------------------------------------------------

  subroutine postpr_postprocess_ip(bridge,wopos,itste,ttime,pdime,TAG1,TAG2)

    character(*), intent(in)           :: wopos(*)
    integer(ip),  intent(inout)        :: bridge(*)
    integer(ip),  intent(in)           :: itste
    real(rp),     intent(in)           :: ttime
    integer(ip),  intent(in)           :: pdime
    integer(ip),  intent(in), optional :: TAG1
    integer(ip),  intent(in), optional :: TAG2
    integer(ip)                        :: ileng,idime,pleng,kpoin,ipoin,ppoin
    integer(ip)                        :: imesh,kinteger,tag1_opt,tag2_opt
    integer(ip),  pointer              :: gisc2(:)
    !
    ! Options
    !
    tag1_opt = 0
    tag2_opt = 0
    if( present(TAG1) ) tag1_opt = TAG1
    if( present(TAG2) ) tag2_opt = TAG2
    !
    ! Header
    !
    kinteger = kind(pdime)
    wwwww(1) = 'ALYA '
    wwwww(2) = 'V0003'
    wwwww(3) = wopos(1)
    wwwww(4) = wopos(2)
    wwwww(5) = wopos(3)
    iiiii(5) = int(tag1_opt,4)
    iiiii(6) = int(tag2_opt,4)
    wwwww(6) = 'INTEG'
    wwwww(7) = '4BYTE'
    iiiii(1) = int(pdime,4)
    iiiii(2) = 0_4
    iiiii(4) = int(itste,4)
    rrrrr(1) = ttime

    if( kinteger == 8_ip ) wwwww(7) = '8BYTE'
    wopos_pos(1) = wopos(1)
    !
    ! Mesh multiplication?
    !
    if( kfl_reawr == 0 ) then
       imesh = kfl_posdi
    else
       imesh = ndivi
    end if
    !
    ! Type
    !
    if( wopos(3) == 'NPOIN' ) then

       parii = NPOIN_TYPE
       pleng = meshe(imesh) % npoin
       if( IMASTER ) pleng = meshe(imesh) % npoin_total

    else if( wopos(3) == 'NELEM' ) then

       parii = NELEM_TYPE
       pleng = meshe(imesh) % nelem
       if( IMASTER ) pleng = meshe(imesh) % nelem_total

    else if( wopos(3) == 'NBOUN' ) then

       parii = NBOUN_TYPE
       pleng =  meshe(imesh) % nboun
       if( IMASTER ) pleng = meshe(imesh) % nboun_total

    else

       call runend('MOD_POSTPR - POSTPR_POSTPROCESS_IP: UNDEFINED POSTPROCESS')

    end if
    iiiii(2) = int(pleng,4)
    !
    ! Get filter
    !
    call fildef(3_ip)
    !
    ! Output solution
    !
    call postpr_open_file()
    if( kfl_reawr < 0 ) then
       kfl_reawr = abs(kfl_reawr)
       return
    end if
    !
    ! Read or write header
    !
    call postpr_header()

    if( ISEQUEN ) then

       if( kfl_filte == 0 .and. kfl_reawr /= 1 ) then
          write(lun_postp) pleng
          write(lun_postp) ( bridge(ileng), ileng = 1,pdime*pleng )
       else if( kfl_reawr == 1 ) then
          read(lun_postp) pleng
          read(lun_postp)  ( bridge(ileng), ileng = 1,pdime*pleng )
       else if( kfl_filte /= 0 ) then
          call runend('POSTPR_POSTPROCESS_IP: FILTER NOT CODED')
       end if

    else if( IMASTER ) then

       do kfl_desti_par = 1,npart

          if( kfl_filte == 0 ) then
             !
             ! Master without filter
             !
             if(      wopos(3) == 'NPOIN' ) then
                pleng = meshe(imesh) % npoin_par(kfl_desti_par)
             else if( wopos(3) == 'NELEM' ) then
                pleng = meshe(imesh) % nelem_par(kfl_desti_par)
             else if( wopos(3) == 'NBOUN' ) then
                pleng = meshe(imesh) % nboun_par(kfl_desti_par)
             end if
             call postpr_memory(1_ip,pdime*pleng,0_ip)
             if( pleng > 0 ) call parari('RCV',0_ip,pdime*pleng,giscp)
             write(lun_postp) pleng
             write(lun_postp) ( giscp(ileng), ileng = 1,pdime*pleng )
             call postpr_memory(3_ip,pdime*pleng,0_ip)

          else
             !
             ! Master with filter
             !
             call parari('RCV',0_ip,1_ip,pleng)
             if( pleng > 0 ) then
                allocate(gisc2(pleng))
                call parari('RCV',0_ip,pleng,gisc2)
                call postpr_memory(1_ip,pdime*pleng,0_ip)
                call parari('RCV',0_ip,pdime*pleng,giscp)
                write(lun_postp) pleng
                write(lun_postp) ( gisc2(ileng), ileng = 1,pleng )
                write(lun_postp) ( giscp(ileng), ileng = 1,pdime*pleng )
                call postpr_memory(3_ip,pdime*pleng,0_ip)
                deallocate(gisc2)
             else
                write(lun_postp) pleng
             end if
          end if
       end do

    else if( ISLAVE ) then

       kfl_desti_par = 0
       if( kfl_filte == 0 ) then
          !
          ! Slaves without filter
          !
          if( pleng > 0 ) call parari('SND',0_ip,pdime*pleng,bridge)
       else
          !
          ! Salves with filter
          !
          call postpr_memory(1_ip,pdime*npoin,0_ip)
          pleng = 0
          loop_ppoin: do ppoin = 1,npoin
             ipoin = gefil(ppoin)
             if( ipoin /= 0 ) then
                kpoin = (ipoin-1) * pdime
                do idime = 1,pdime
                   pleng = pleng + 1
                   kpoin = kpoin + 1
                   giscp(pleng) = bridge(kpoin)
                end do
             else
                exit loop_ppoin
             end if
          end do loop_ppoin
          pleng = pleng / pdime
          call parari('SND',0_ip, 1_ip,pleng)
          if( pleng > 0 ) then
             call parari('SND',0_ip,pleng,gefil)
             call parari('SND',0_ip,pdime*pleng,giscp)
          end if
          call postpr_memory(3_ip,pdime*npoin,0_ip)
       end if

    end if

    call postpr_close_file()

  end subroutine postpr_postprocess_ip

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-02-23
  !> @brief   Header
  !> @details Write/read header of a postprocess file
  !>
  !-----------------------------------------------------------------------

  subroutine postpr_header()

    integer(ip) :: ii
    integer(4)  :: ihead

    npari = 0
    nparr = 0
    nparc = 0
    nparx = 0
    pardi = int(iiiii(1),ip)

    if( INOTSLAVE ) then
       if( IMASTER ) then
          wwwww(8) = 'PARAL'
          iiiii(3) = int(npart,4)
       else
          wwwww(8) = 'SEQUE'
          iiiii(3) = int(1,4)
       end if
       if( kfl_filte == 0 ) then
          wwwww(9) = 'NOFIL'
       else
          wwwww(9) = 'FILTE'
       end if

       if( kfl_reawr == 1 ) then
          read(lun_postp)  ihead
          if( ihead /= 1234 ) call runend('MOD_POSTPR: ALYA FILE HAS A WRONG BINARY FORMAT')
          read(lun_postp)  wwww8(1)
          read(lun_postp)  wwww8(2)
          read(lun_postp)  wwww8(3)
          read(lun_postp)  wwww8(4)
          read(lun_postp)  wwww8(5)
          read(lun_postp)  wwww8(6)
          read(lun_postp)  wwww8(7)
          read(lun_postp)  wwww8(8)
          read(lun_postp)  wwww8(9)
          read(lun_postp)  iiiii(1)
          read(lun_postp)  iiiii(2)
          read(lun_postp)  iiiii(3)
          read(lun_postp)  iiiii(4)
          if( wwww8(2) == 'V0003' ) then
             read(lun_postp)  iiiii(5)
             read(lun_postp)  iiiii(6)
          end if
          read(lun_postp)  rrrrr(1)
          do ii = 1,nwwww
             wwwww(ii)=wwww8(ii)(1:5)
          end do
       else
          wwww8(1) = 'ALYAPOST'
          do ii = 2,nwwww
             wwww8(ii)=wwwww(ii)//'   '
          end do
          ihead = 1234_4
          write(lun_postp) ihead
          write(lun_postp) wwww8(1)
          write(lun_postp) wwww8(2)
          write(lun_postp) wwww8(3)
          write(lun_postp) wwww8(4)
          write(lun_postp) wwww8(5)
          write(lun_postp) wwww8(6)
          write(lun_postp) wwww8(7)
          write(lun_postp) wwww8(8)
          write(lun_postp) wwww8(9)
          write(lun_postp) iiiii(1)
          write(lun_postp) iiiii(2)
          write(lun_postp) iiiii(3)
          write(lun_postp) iiiii(4)
          write(lun_postp) iiiii(5)
          write(lun_postp) iiiii(6)
          write(lun_postp) rrrrr(1)
          if( kfl_reawr == 0 ) then
             if( ncoun_pos == 0 ) then
                ncoun_pos = ncoun_pos + 1
                write( lun_pos02,'(a,1x,i9,1x,e12.5)') 'START',iiiii(4),rrrrr(1)
             end if
             write( lun_pos02,'(a)') wwwww(3)(1:5)
             call iofile_flush_unit(lun_pos02)
          end if

       end if
    end if

  end subroutine postpr_header

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-02-19
  !> @brief   Allocate and deallocate memory for generic scalar and vector arrays
  !> @details Allocate and deallocate memory for generic scalar and vector arrays:
  !>          ITASK=0 ... Allocate memory
  !>          ITASK=2 ... Deallocate memory
  !>
  !-----------------------------------------------------------------------

  subroutine postpr_memory(itask,ndim1,ndim2)

    use def_parame
    use def_master
    use mod_memory, only : memory_alloca
    use mod_memory, only : memory_deallo
    implicit none
    integer(ip), intent(in) :: itask,ndim1,ndim2
    integer(4)              :: istat

    select case(itask)

    case(0_ip)
       !
       ! Allocate memory for real
       !
       if(ndim1>0.and.ndim2==0) then
          call memory_alloca(memke,'GESCP','postpr_memory',gescp,ndim1)
       else if(ndim1>0.and.ndim2>0) then
          call memory_alloca(memke,'GEVEP','postpr_memory',gevep,ndim1,ndim2)
       else if(ndim1<0.and.ndim2>0) then
          call runend('POSTPR_MEMORY: NOT CODED')
          !allocate(getep(-ndim1,-ndim1,ndim2),stat=istat)
          !call memchk(zero,istat,memke,'GETEP','postpr_memory',getep)
       end if

    case(1_ip)
       !
       ! Allocate memory for integer
       !
       if(ndim1/=0.and.ndim2==0) then
          call memory_alloca(memke,'GISCP','postpr_memory',giscp,ndim1)
       else if(ndim1/=0.and.ndim2/=0) then
          call memory_alloca(memke,'GIVEP','postpr_memory',givep,ndim1,ndim2)
       end if

    case(2_ip)
       !
       ! Deallocate memory for real
       !
       if(ndim1>0.and.ndim2==0) then
          call memory_deallo(memke,'GESCP','postpr_memory',gescp)
       else if(ndim1>0.and.ndim2>0) then
          call memory_deallo(memke,'GEVEP','postpr_memory',gevep)
       else if(ndim1<0.and.ndim2>0) then
          call memory_deallo(memke,'GETEP','postpr_memory',getep)
       end if

    case(3_ip)
       !
       ! Deallocate memory for integer
       !
       if(ndim1/=0.and.ndim2==0) then
          call memory_deallo(memke,'GISCP','postpr_memory',giscp)
       else if(ndim1/=0.and.ndim2/=0) then
          call memory_deallo(memke,'GIVEP','postpr_memory',givep)
       end if

    case(4_ip)
       !
       ! Allocate memory for r3p
       !
       call runend('POSTPR_MEMORY: NOT PRGRAMMED')

    case(5_ip)
       !
       ! Deallocate memory for r3p
       !
       call runend('POSTPR_MEMORY: NOT PRGRAMMED')

    end select

  end subroutine postpr_memory

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-02-19
  !> @brief   Open file
  !> @details Open postprocess/restart file
  !>
  !-----------------------------------------------------------------------

  subroutine postpr_open_file(TAG1,TAG2)

    use mod_iofile
    use mod_opfpos

    integer(ip),  intent(in), optional  :: TAG1
    integer(ip),  intent(in), optional  :: TAG2
    character(200)                      :: filsa
    logical                             :: mesh
    character(200)                         :: messa
    call opfpos_name(wopos_pos(1), ".alyabin", mesh, TAG1, TAG2)
    if (mesh) then
       if(INOTSLAVE) then
          call iofile(0_ip,lun_postp,fil_postp,'ALYA POSTPROCESS FILE','unknown','unformatted')
       end if
    else if( kfl_reawr == 0 ) then
       !
       ! Open postprocess file name
       !
       if( INOTSLAVE ) then
          call iofile(0_ip,lun_postp,fil_postp,'ALYA POSTPROCESS FILE','unknown','unformatted')
          write(lun_pos01,'(a)') trim(fil_postp)
          call iofile_flush_unit(lun_pos01)
       end if

    else if( kfl_reawr == 1 ) then
       !
       ! Open restart file name for reading
       !
       if( INOTSLAVE )  then
          call iofile(4_ip,lun_postp,fil_postp,'ALYA RESTART FILE','unknown','unformatted')
       end if
       call parari('BCT',0_ip,1_ip,kfl_reawr)
       if( kfl_reawr < 0 ) then
          messa = 'CANNOT OPEN RESTART FILE: '//trim(fil_postp)
          call outfor(2_ip,0_ip,trim(messa))
          call livinf(0_ip,trim(messa),0_ip)
       end if
       if( kfl_reawr == 1 .and. INOTSLAVE ) then
          call iofile(0_ip,lun_postp,fil_postp,'ALYA RESTART FILE','unknown','unformatted')
       end if

    else if( kfl_reawr == 2 .and. INOTSLAVE ) then
       !
       ! Open restart file name for writing
       !
       call iofile(0_ip,lun_postp,fil_postp,'ALYA RESTART FILE','unknown','unformatted')
    end if

  end subroutine postpr_open_file

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-02-19
  !> @brief   Close file
  !> @details Close postprocess/restart file
  !>
  !-----------------------------------------------------------------------

  subroutine postpr_close_file()

    use mod_iofile

    if( INOTSLAVE ) then
       if( kfl_reawr == 0 ) then
          call iofile(2_ip,lun_postp,' ','ALYA POSTPROCESS FILE')
       else
          call iofile(2_ip,lun_postp,' ','ALYA RESTART FILE')
       end if
    end if

  end subroutine postpr_close_file

  function postpr_is_mpio(wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2) result(ret)
    character(*), intent(in)            :: wopos(*)
    integer(ip) , intent(in)            :: itste
    real(rp)    , intent(in)            :: ttime
    integer(ip),  intent(in)            :: pdime
    integer(ip),  intent(in), optional  :: pleng_opt
    integer(ip),  intent(in), optional  :: TAG1
    integer(ip),  intent(in), optional  :: TAG2
    logical(lg)                         :: ret
    ret=.false.
    if ((kfl_mpio_post > IO_CLASSIC .AND. kfl_reawr == 0) .OR. (kfl_mpio_rst > IO_CLASSIC .AND. kfl_reawr /= 0)) then
       wopos_mpio(1) =  wopos(1)
       wopos_mpio(2) =  wopos(2)
       wopos_mpio(3) =  wopos(3)
       pdime_mpio = pdime
       ttime_mpio = ttime
       itste_mpio = itste
       if ( present(TAG1) ) then
          tag1_mpio  = TAG1
       else
          tag1_mpio = -1
       end if
       if ( present(TAG2) ) then
          tag2_mpio  = TAG2
       else
          tag2_mpio = -1
       end if
       if( present(pleng_opt) ) then
          pleng_mpio = pleng_opt
       else
          pleng_mpio = 0
       end if
       ret=.true.
    end if
  end function postpr_is_mpio

end module mod_postpr



!> @}

