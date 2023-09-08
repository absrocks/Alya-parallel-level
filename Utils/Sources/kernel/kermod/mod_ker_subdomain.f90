  !-----------------------------------------------------------------------
  !> @addtogroup Kermod
  !> @{
  !> @file    mod_ker_subdomain.f90
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Subdomain
  !> @details Subdomain treatment: motion, etc.
  !>
  !-----------------------------------------------------------------------

module mod_ker_subdomain

  use def_kintyp,   only : ip,rp
  use def_parame,   only : pi
  use def_master,   only : INOTMASTER
  use def_master,   only : ID_KERMOD
  use def_master,   only : mem_modul
  use def_domain,   only : nsubd
  use mod_memory,   only : memory_copy
  use mod_memory,   only : memory_alloca
  use mod_memory,   only : memory_deallo
  use mod_maths,    only : maths_normalize_vector
  use mod_messages, only : messages_live
  use def_parame,   only : pi
!!  use Begste,       only : lnsub
  use def_kermod
  
  implicit none
  integer(ip), parameter :: ALE              = 1
  integer(ip), parameter :: TOTAL_LAGRANGIAN = 2
  real(rp),    pointer   :: coord_initial(:,:)
  integer(ip), pointer   :: mask(:)
  
  private

  public :: ker_subdomain_initialization
  public :: ker_subdomain_read_data
  public :: ker_subdomain_function_name_to_number
  public :: ker_subdomain_update_coordinates
  public :: ker_subdomain_parall
  public :: ker_subdomain_motion_exists
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   If motion exists
  !> @details If motion exists
  !> 
  !-----------------------------------------------------------------------

  function ker_subdomain_motion_exists()

    logical(lg) :: ker_subdomain_motion_exists
    
    ker_subdomain_motion_exists = .false.
    if( associated(subdomain) ) then
       if( maxval(subdomain(:) % kfl_formulation) > 0 ) then
          ker_subdomain_motion_exists = .true.
       end if
    end if

  end function ker_subdomain_motion_exists
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Send data to slaves
  !> @details Send data to slaves
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_subdomain_parall()

    integer(ip) :: ii,isubd

    if( associated(subdomain) ) then
       do isubd = 1,nsubd
          call iexcha(subdomain(isubd) % kfl_formulation)
          call iexcha(subdomain(isubd) % rotation_number)
          call cexcha(5_ip,subdomain(isubd) % rotation_name)
          call iexcha(subdomain(isubd) % omega_number)
          call cexcha(5_ip,subdomain(isubd) % omega_name)
          call rexcha(subdomain(isubd) % omega)
          do ii = 1,3
             call rexcha(subdomain(isubd) % rotation_axis(ii))
          end do
          do ii = 1,3
             call rexcha(subdomain(isubd) % rotation_center(ii))
          end do
          call rexcha(subdomain(isubd) % rotation_angle)
          call iexcha(subdomain(isubd) % translation_number)
          call cexcha(5_ip,subdomain(isubd) % translation_name)
          do ii = 1,3
             call rexcha(subdomain(isubd) % translation(ii))
          end do
       end do
    end if

  end subroutine ker_subdomain_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Initialization of structure
  !> @details Initialization of structure
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_subdomain_initialization()
    
    integer(ip) :: isubd

    
    nullify(coord_initial)
    nullify(mask)
    allocate(subdomain(nsubd))
    if( associated(subdomain) ) then
       do isubd = 1,nsubd
          subdomain(isubd) % kfl_formulation    = 0
          subdomain(isubd) % rotation_number    = 0
          subdomain(isubd) % rotation_name      = 'NONE'
          subdomain(isubd) % omega_number       = 0
          subdomain(isubd) % omega_name         = 'NONE'
          subdomain(isubd) % omega              = 0.0_rp
          subdomain(isubd) % rotation_axis      = 0.0_rp
          subdomain(isubd) % rotation_center    = 0.0_rp
          subdomain(isubd) % rotation_angle     = 0.0_rp
          subdomain(isubd) % translation_number = 0
          subdomain(isubd) % translation_name   = 'NONE'
          subdomain(isubd) % translation        = 0.0_rp
       end do
    end if 

  end subroutine ker_subdomain_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Read data
  !> @details Read subdomain data
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_subdomain_read_data()

    use def_inpout
    use mod_ecoute, only : ecoute
    integer(ip) :: isubd
  
    isubd = getint('SUBDO',1_ip,'#SUBDOMAIN')
    if( isubd < 0 .or. isubd > nsubd ) call runend('KER_SUBDOMAIN_READ_DATA: WRONG SUBDOMAIN NUMBER')
    call ecoute('ker_readat')
    do while( words(1) /= 'ENDSU' ) 
       if( words(1) == 'FORMU' ) then
          if( words(2) == 'ALE  ' ) then
             subdomain(isubd) % kfl_formulation = ALE
          else if( words(2) == 'TOTAL' ) then
             subdomain(isubd) % kfl_formulation = TOTAL_LAGRANGIAN
          end if
       else if( words(1) == 'DISPL' ) then
          if( words(2) == 'SPACE' ) then
              subdomain(isubd) % translation_name = getcha('SPACE','NULL ','#Space/time Function name')
          end if
       else if( words(1) == 'ROTAT' ) then
          if( words(2) == 'SPACE' ) then
              subdomain(isubd) % rotation_name = getcha('SPACE','NULL ','#Space/time Function name')
          else if( words(2) == 'OMEGA' ) then
              subdomain(isubd) % omega_name = getcha('OMEGA','NULL ','#Space/time Function name')
          else if( words(2) == 'VELOC' ) then
             subdomain(isubd) % omega = param(2)
          else if( words(2) == 'AXIS ' ) then
             subdomain(isubd) % rotation_axis(1:3) = param(2:4)
             call maths_normalize_vector(3_ip,subdomain(isubd) % rotation_axis)
          else if( words(2) == 'CENTE' ) then
             subdomain(isubd) % rotation_center(1:3) = param(2:4)
          end if
       end if
       call ecoute('ker_readat')
    end do
    
  end subroutine ker_subdomain_read_data

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Function name to function number
  !> @details Function name to function number
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_subdomain_function_name_to_number()

    use mod_ker_space_time_function, only : space_time_function_number
    integer(ip) :: isubd

    if( associated(subdomain) ) then
       do isubd = 1,size(subdomain)
          if( trim(subdomain(isubd) % rotation_name) /= 'NONE' ) &
               subdomain(isubd) % rotation_number    = space_time_function_number(subdomain(isubd) % rotation_name)
          if( trim(subdomain(isubd) % omega_name) /= 'NONE' ) &
               subdomain(isubd) % omega_number    = space_time_function_number(subdomain(isubd) % omega_name)
          if( trim(subdomain(isubd) % translation_name) /= 'NONE' ) &
               subdomain(isubd) % translation_number = space_time_function_number(subdomain(isubd) % translation_name)
       end do
    end if
    
  end subroutine ker_subdomain_function_name_to_number
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Update subdomain coordinates
  !> @details Update subdomain coordinates
  !>          See https://en.wikipedia.org/wiki/Rotation_matrix
  !>
  !-----------------------------------------------------------------------

  subroutine ker_subdomain_update_coordinates() 

    use def_master,                  only : cutim
    use def_master,                  only : dtinv
    use def_master,                  only : dispm
    use def_master,                  only : velom
    use def_coupli,                  only : mcoup
    use def_coupli,                  only : coupling_type
    use mod_coupling_memory
    use def_domain,                  only : coord,npoin_2,npoin
    use def_domain,                  only : lnods,lnnod
    use def_domain,                  only : ndime,nelem
    use def_domain,                  only : npoin,lesub
    use mod_ker_space_time_function, only : ker_space_time_function
    use mod_elsest,                  only : elsest_deallocate
    use mod_elsest,                  only : elsest_initialization
    use mod_par_bin_structure,       only : par_bin_structure
    use mod_par_bin_structure,       only : par_bin_structure_deallocate
    use mod_communications,          only : PAR_MAX
    
    integer(ip)          :: ipoin,ielem,isubd,icoup
    integer(ip)          :: ifunc_translation,kk
    integer(ip)          :: ifunc_rotation
    integer(ip)          :: ifunc_omega
    real(rp)             :: rotma(3,3)
    real(rp)             :: coord0(3),x0(3),coord1(3)
    real(rp)             :: w_vec(3),r_vec(3)
    integer(ip)          :: istat


    !if( ittim < 100 ) then
    !   return
    !else if( ittim == 100 ) then
    !   sttim = cutim-dtime
    !end if
    
    if ( ker_subdomain_motion_exists() ) then

       if( INOTMASTER ) then
          
          if( .not. associated(coord_initial) ) then
             call memory_copy(mem_modul(1:2,ID_KERMOD),'COORD_INITIAL','ker_subdomain_update_coordinates',coord,coord_initial,'DO_NOT_DEALLOCATE')
          end if
          
          coord0   = 0.0_rp
          x0       = 0.0_rp
          velom    = 0.0_rp
          dispm    = 0.0_rp
          
          if( .not. associated(mask) ) then
             call memory_alloca(mem_modul(1:2,ID_KERMOD),'mask','ker_subdomain_update_coordinates',mask,npoin_2)
             mask = 0
             do isubd = 1,nsubd
                if( subdomain(isubd) % kfl_formulation > 0 ) then
                   do ielem = 1,nelem
                      if( lesub(ielem) == isubd ) & 
                         mask(lnods(1:lnnod(ielem),ielem)) = isubd
                   end do
                end if
             end do
          end if

          do isubd = 1,nsubd

             if( subdomain(isubd) % kfl_formulation > 0 ) then

                ifunc_rotation = subdomain(isubd) % rotation_number
                
                if( ifunc_rotation > 0 ) then
                   call ker_space_time_function(&
                        ifunc_rotation,x0(1),x0(2),x0(3),cutim,subdomain(isubd) % rotation_angle)
                else 
                   subdomain(isubd) % rotation_angle = subdomain(isubd) % omega * cutim
                end if
                call rotmat(ndime,subdomain(isubd) % rotation_angle,subdomain(isubd) % rotation_axis,rotma(:,:))

                ifunc_omega = subdomain(isubd) % omega_number
                if( ifunc_omega > 0 ) then
                   call ker_space_time_function(&
                        ifunc_omega,x0(1),x0(2),x0(3),cutim,subdomain(isubd) % omega)
                end if
                w_vec(1:ndime) = subdomain(isubd) % rotation_axis * subdomain(isubd) % omega

                ifunc_translation = subdomain(isubd) % translation_number
                if( ifunc_translation > 0 ) then
                   call ker_space_time_function(&
                        ifunc_translation,x0(1),x0(2),x0(3),cutim,subdomain(isubd) % translation(1:ndime))
                end if

                do ipoin = 1,npoin
                   if( mask(ipoin) == isubd ) then
!!                   if( lnsub(ipoin) == isubd ) then
                      coord0(1:ndime)      = coord_initial(1:ndime,ipoin)
                      r_vec(1:ndime)       = coord0 - subdomain(isubd) % rotation_center
                      coord1               = matmul(rotma,r_vec) + subdomain(isubd) % rotation_center
                      coord1(1:ndime)      = coord1(1:ndime) + subdomain(isubd) % translation(1:ndime)
                      call vecpro (w_vec(1:ndime), r_vec(1:ndime), velom(1:ndime, ipoin), ndime)
                      dispm(1:ndime,ipoin,1) = ( coord1(1:ndime) - coord_initial(1:ndime,ipoin) )
                      coord(1:ndime,ipoin) = coord1(1:ndime)
                      velom(1:ndime,ipoin) = velom(1:ndime,ipoin) + subdomain(isubd) % translation(1:ndime) * dtinv
                   end if 
                end do
             end if
          end do
          !
          ! Update wet point coordinates
          !
          do icoup = 1,mcoup
             do kk = 1,coupling_type(icoup) % wet % npoin_wet
                ipoin = coupling_type(icoup) % wet % lpoin_wet(kk)
                coupling_type(icoup) % wet % coord_wet(1:ndime,kk) = coord(1:ndime,ipoin)
             end do
          end do
         
       end if
    !
    ! Recompute couplings
    !
       call messages_live('MOVING SUBDOMAIN UPDATE','START SECTION')
       if( INOTMASTER ) then
          call elsest_deallocate(ielse)
          call elsest_initialization()
       end if
       call elsini()
       call par_bin_structure_deallocate()
       call par_bin_structure()
       istat = 1_ip
       do while(istat /= 0_ip)
          do icoup = 1,mcoup
             call COU_DEALLOCATE_SINGLE_COUPLING(coupling_type(icoup),WET=.false.)
             call COU_INITIALIZATION_SINGLE_COUPLING(coupling_type(icoup),WET=.false.,INPUT= .false.)
          end do
          call cou_initialize_coupling(istat)
          call PAR_MAX(istat)
       end do
       call domarr(2_ip)
       call messages_live('MOVING SUBDOMAIN UPDATE','END SECTION')
    end if
    
  end subroutine ker_subdomain_update_coordinates
  
end module mod_ker_subdomain
!> @}
