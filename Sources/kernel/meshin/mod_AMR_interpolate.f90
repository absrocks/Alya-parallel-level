!-----------------------------------------------------------------------
!> @addtogroup AMR
!> @{
!> @file    mod_AMR_interpolate.f90
!> @author  houzeaux
!> @date    2020-06-16
!> @brief   AMR interpolation
!> @details Interpolation for AMR
!-----------------------------------------------------------------------

module mod_AMR_interpolate

  use def_kintyp_basic,   only : ip,rp,lg
  use def_domain,         only : memor_dom
  use mod_memory_basic,   only : memory_alloca
  use mod_memory_basic,   only : memory_deallo
  use mod_memory_basic,   only : memory_size
  use mod_memory_basic,   only : memory_copy
  use def_coupli,         only : typ_color_coupling
  use mod_strings,        only : integer_to_string
  use mod_interpolation,  only : COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_messages,       only : messages_live
  use mod_communications, only : PAR_MAX
  use mod_type,           only : my_type
  use def_AMR
  use def_master
  implicit none
  private

  interface AMR_interpolate
     module procedure &
          AMR_interpolate_IP_1,&
          AMR_interpolate_IP_2,&
          AMR_interpolate_RP_1,&
          AMR_interpolate_RP_2,&
          AMR_interpolate_RP_3
  end interface AMR_interpolate
  
  type(typ_color_coupling), pointer :: current_coupling
  integer(ip)                       :: nenti_new
  integer(8)                        :: memor_loc(2)
  character(50)                     :: my_variable_name
  logical(lg)                       :: if_verbose

  public :: AMR_interpolate

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-07
  !> @brief   Interopolation
  !> @details Interpolation functions
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_interpolate_IP_1(xx,wtype,posit,memor,variable_name,coupling,VERBOSE)

    integer(ip),              pointer,           intent(inout) :: xx(:)
    character(*),                                intent(in)    :: wtype
    integer(ip),                       optional, intent(in)    :: posit
    integer(8),                        optional, intent(inout) :: memor(2)
    character(*),                      optional, intent(in)    :: variable_name
    type(typ_color_coupling), target,  optional, intent(in)    :: coupling
    logical(lg),                       optional, intent(in)    :: VERBOSE
    integer(ip),              pointer                          :: xx_target(:)
    integer(ip)                                                :: ii,ndofn,jj

    nullify(xx_target)

    ii = memory_size(xx)
    call PAR_MAX(ii)
    if( ii <= 0 ) return

    call AMR_options(wtype,memor,variable_name,coupling,VERBOSE,my_type(xx))
    
    if( if_verbose ) call messages_live('INTERPOLATE '//trim(variable_name))
    
    call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target,nenti_new)    
    call COU_GET_INTERPOLATE_POINTS_VALUES(xx,xx_target,current_coupling)

    call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx)
    call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx,nenti_new)

    if( associated(current_coupling % commd % lrecv_perm) ) then
       do ii = 1,nenti_new
          jj = current_coupling % commd % lrecv_perm(ii)
          xx(jj) = xx_target(ii)
       end do
    else
       do ii = 1,nenti_new
          xx(ii) = xx_target(ii)
       end do       
    end if
    
    call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target)

    if( present(memor) ) then
       memor     = memor_loc
    else
       memor_dom = memor_loc 
    end if

  end subroutine AMR_interpolate_IP_1

  subroutine AMR_interpolate_IP_2(xx,wtype,posit,memor,variable_name,coupling,VERBOSE)

    integer(ip),              pointer,           intent(inout) :: xx(:,:)
    character(*),                                intent(in)    :: wtype
    integer(ip),                       optional, intent(in)    :: posit
    integer(8),                        optional, intent(inout) :: memor(2)
    character(*),                      optional, intent(in)    :: variable_name
    type(typ_color_coupling), target,  optional, intent(in)    :: coupling
    logical(lg),                       optional, intent(in)    :: VERBOSE
    integer(ip),              pointer                          :: xx_target(:,:)
    integer(ip)                                                :: ii,ndofn,jj

    nullify(xx_target)

    ii = memory_size(xx)
    call PAR_MAX(ii)
    if( ii <= 0 ) return

    ndofn = memory_size(xx,1_ip)
    call AMR_options(wtype,memor,variable_name,coupling,VERBOSE,my_type(xx))
    
    if( if_verbose ) call messages_live('INTERPOLATE '//trim(variable_name))

    call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target,ndofn,nenti_new)    
    call COU_GET_INTERPOLATE_POINTS_VALUES(xx,xx_target,current_coupling)
 
    call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx)
    call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx,ndofn,nenti_new)

    if( posit == 1) call runend('IP_2 COUPLING NOT CODED FOR AMR')
    
    if( associated(current_coupling % commd % lrecv_perm) ) then
       do ii = 1,nenti_new
          jj = current_coupling % commd % lrecv_perm(ii)
          xx(:,jj) = xx_target(:,ii)
       end do
    else
       do ii = 1,nenti_new
          xx(:,ii) = xx_target(:,ii)
       end do
    end if
    
    call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target)

    if( present(memor) ) then
       memor     = memor_loc
    else
       memor_dom = memor_loc 
    end if

  end subroutine AMR_interpolate_IP_2

  subroutine AMR_interpolate_RP_3(xx,wtype,posit,memor,variable_name,coupling,VERBOSE)

    real(rp),                 pointer,           intent(inout) :: xx(:,:,:)
    character(*),                                intent(in)    :: wtype
    integer(ip),                       optional, intent(in)    :: posit
    integer(8),                        optional, intent(inout) :: memor(2)
    character(*),                      optional, intent(in)    :: variable_name
    type(typ_color_coupling), target,  optional, intent(in)    :: coupling
    logical(lg),                       optional, intent(in)    :: VERBOSE
    real(rp),                 pointer                          :: xx_target(:,:)
    real(rp),                 pointer                          :: xx_cpy(:,:,:),xx2(:,:)
    integer(ip)                                                :: ii,jj,ndofn,ndim3,idim3
    
    nullify(xx_target)
    nullify(xx_cpy)

    ii = memory_size(xx)
    call PAR_MAX(ii)
    if( ii <= 0 ) return

    ndofn = memory_size(xx,1_ip)
    ndim3 = memory_size(xx,3_ip)
    
    call AMR_options(wtype,memor,variable_name,coupling,VERBOSE,my_type(xx))
    
    if( if_verbose ) call messages_live('INTERPOLATE '//trim(variable_name))
    
    if( posit == 2 ) then
       
       call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target,ndofn,nenti_new)    
       call memory_copy  (memor_loc,'XX_CPY'              ,'AMR_interpolate',xx,xx_cpy) 
       call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx)
       call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx,ndofn,nenti_new,ndim3)
       if( .not. associated(xx_cpy) ) then
          call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_cpy,1_ip,1_ip,ndim3)
       end if
       !
       ! Coupling
       !
       do idim3 = 1,ndim3
          xx2 => xx_cpy(:,:,idim3)
          call COU_GET_INTERPOLATE_POINTS_VALUES(xx2,xx_target,current_coupling)
          if( associated(current_coupling % commd % lrecv_perm) ) then
             do ii = 1,nenti_new
                jj = current_coupling % commd % lrecv_perm(ii)
                xx(:,jj,idim3) = xx_target(:,ii)
             end do
          else
             do ii = 1,nenti_new
                xx(:,ii,idim3) = xx_target(:,ii)
             end do
          end if          
       end do

       call memory_deallo(memor_loc,'XX_CPY'              ,'AMR_interpolate',xx_cpy)
       call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target)
 
    else
       call runend('AMR_interpolate_RP_3: NOT CODED')
    end if

    if( present(memor) ) then
       memor     = memor_loc
    else
       memor_dom = memor_loc 
    end if

  end subroutine AMR_interpolate_RP_3

  subroutine AMR_interpolate_RP_2(xx,wtype,posit,memor,variable_name,coupling,VERBOSE)

    real(rp),                 pointer,           intent(inout) :: xx(:,:)
    character(*),                                intent(in)    :: wtype
    integer(ip),                       optional, intent(in)    :: posit
    integer(8),                        optional, intent(inout) :: memor(2)
    character(*),                      optional, intent(in)    :: variable_name
    type(typ_color_coupling), target,  optional, intent(in)    :: coupling
    logical(lg),                       optional, intent(in)    :: VERBOSE
    real(rp),                 pointer                          :: xx_target(:,:)
    integer(ip)                                                :: ii,ndofn,idofn,jj,nenti
    real(rp),                 pointer                          :: yy(:,:)

    nullify(xx_target,yy)

    ii = memory_size(xx)
    call PAR_MAX(ii)
    if( ii <= 0 ) return

    ndofn = memory_size(xx,1_ip)
    
    call AMR_options(wtype,memor,variable_name,coupling,VERBOSE,my_type(xx))        
    if( if_verbose ) call messages_live('INTERPOLATE '//trim(variable_name))    
    !
    ! Coupling
    !
    if( posit == 2_ip ) then

       call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target,ndofn,nenti_new)    
       call COU_GET_INTERPOLATE_POINTS_VALUES(xx,xx_target,current_coupling)
       
       call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx)
       call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx,ndofn,nenti_new)
       
       if( associated(current_coupling % commd % lrecv_perm) ) then
          do ii = 1,nenti_new
             jj = current_coupling % commd % lrecv_perm(ii)
             xx(:,jj) = xx_target(:,ii)
          end do
       else
          do ii = 1,nenti_new
             xx(:,ii) = xx_target(:,ii)
          end do
       end if       

    else if( posit == 1_ip ) then
       !
       ! TEMPE(:,NDOFN)
       !
       nenti = memory_size(xx,1_ip)
       ndofn = memory_size(xx,2_ip)
       call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target,ndofn,nenti_new)       
       call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',yy,ndofn,nenti)

       do ii = 1,nenti
          do idofn = 1,ndofn
             yy(idofn,ii) = xx(ii,idofn)
          end do
       end do

       call COU_GET_INTERPOLATE_POINTS_VALUES(yy,xx_target,current_coupling)

       call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',yy)
       call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx)
       call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx,nenti_new,ndofn)

       if( associated(current_coupling % commd % lrecv_perm) ) then
          do ii = 1,nenti_new
             jj = current_coupling % commd % lrecv_perm(ii)
             do idofn = 1,ndofn
                xx(jj,idofn) = xx_target(idofn,ii)
             end do
          end do
       else          
          do ii = 1,nenti_new          
             do idofn = 1,ndofn
                xx(ii,idofn) = xx_target(idofn,ii)
             end do
          end do
       end if

    else
       call runend('MOD_AMR: INTERPOLATION NOT CODED FOR POSIT='//integer_to_string(posit)//' FOR VARIABLE '//trim(variable_name))
    end if

    call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target)

    if( present(memor) ) then
       memor     = memor_loc
    else
       memor_dom = memor_loc 
    end if

  end subroutine AMR_interpolate_RP_2

  subroutine AMR_interpolate_RP_1(xx,wtype,posit,memor,variable_name,coupling,VERBOSE)

    real(rp),                 pointer,           intent(inout) :: xx(:)
    character(*),                                intent(in)    :: wtype
    integer(ip),                       optional, intent(in)    :: posit
    integer(8),                        optional, intent(inout) :: memor(2)
    character(*),                      optional, intent(in)    :: variable_name
    type(typ_color_coupling), target,  optional, intent(in)    :: coupling
    logical(lg),                       optional, intent(in)    :: VERBOSE
    real(rp),                 pointer                          :: xx_target(:)
    integer(ip)                                                :: ii,ndofn,jj

    nullify(xx_target)

    ii = memory_size(xx)
    call PAR_MAX(ii)
    if( ii <= 0 ) return

    call AMR_options(wtype,memor,variable_name,coupling,VERBOSE,my_type(xx))            
    if( if_verbose ) call messages_live('INTERPOLATE '//trim(variable_name))

    call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target,nenti_new)    
    call COU_GET_INTERPOLATE_POINTS_VALUES(xx,xx_target,current_coupling)

    call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx)
    call memory_alloca(memor_loc,trim(my_variable_name),'AMR_interpolate',xx,nenti_new)

    if( associated(current_coupling % commd % lrecv_perm) ) then
       do ii = 1,nenti_new
          jj = current_coupling % commd % lrecv_perm(ii)
          xx(jj) = xx_target(ii)
       end do
    else
       do ii = 1,nenti_new
          xx(ii) = xx_target(ii)
       end do
    end if
    
    call memory_deallo(memor_loc,trim(my_variable_name),'AMR_interpolate',xx_target)

    if( present(memor) ) then
       memor     = memor_loc
    else
       memor_dom = memor_loc 
    end if

  end subroutine AMR_interpolate_RP_1

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-16
  !> @brief   Define local variables
  !> @details Define local variables to avoid repeating the same all the
  !>          time
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_options(wtype,memor,variable_name,coupling,VERBOSE,current_type)

    character(*),                                intent(in)    :: wtype
    integer(8),                        optional, intent(inout) :: memor(2)
    character(*),                      optional, intent(in)    :: variable_name
    type(typ_color_coupling), target,  optional, intent(in)    :: coupling
    logical(lg),                       optional, intent(in)    :: VERBOSE
    character(len=*),                  optional, intent(in)    :: current_type
    
    if(      wtype == 'NELEM' ) then
       nenti_new = nelem_new
    else if( wtype == 'NPOIN' ) then
       nenti_new = npoin_new
    else if( wtype == 'NBOUN' ) then
       nenti_new = nboun_new
    end if
    
    if( present(coupling) ) then
       current_coupling => coupling
    else if( wtype == 'NELEM' ) then
       current_coupling => coupling_AMR_nelem
    else if( wtype == 'NPOIN' ) then
       current_coupling => coupling_AMR_npoin      
    else if( wtype == 'NBOUN' ) then
       current_coupling => coupling_AMR_nboun             
    end if

    if( present(current_type) ) then
       if( current_type == 'integer' .and. wtype == 'NPOIN' .and. .not. present(coupling) ) &
            current_coupling => coupling_AMR_npoin_nearest
    end if
    
    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if

    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = memor_dom
    end if

    if( present(VERBOSE) ) then
       if_verbose = VERBOSE
    else
       if_verbose = .true.
    end if
    
  end subroutine AMR_options
  
end module mod_AMR_interpolate
!> @}
