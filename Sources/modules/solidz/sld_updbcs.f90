!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_updbcs.f90
!> @author  Solidz
!> @date
!> @brief   Boundary conditions update
!> @details This routine updates the boundary conditions:
!>          0. At the beginning of the run
!>          1. Before a time step begins
!>          2. Before a global iteration begins
!>          3. Before an inner iteration begins
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_updbcs(itask)

  use def_kintyp,                  only : ip, rp
  use def_master,                  only : INOTEMPTY, INOTMASTER
  use def_master,                  only : FUNCTION_TIME
  use def_master,                  only : FUNCTION_MODULE
  use def_master,                  only : FUNCTION_SPACE_TIME
  use def_master,                  only : cutim, timei, timef
  use def_master,                  only : ITER_K, ITER_AUX, TIME_N
  use def_master,                  only : displ
  use def_domain,                  only : ndime, npoin, kfl_codno
  use def_domain,                  only : lpoty, coord
  use def_kermod,                  only : number_space_time_function
  use mod_ker_space_time_function, only : ker_space_time_function
  use def_solidz,                  only : bvess_sld
  use def_solidz,                  only : kfl_conbc_sld, kfl_funno_sld
  use def_solidz,                  only : nfunc_sld
  use def_solidz,                  only : kfl_local_sld, kfl_conta_stent
  use def_solidz,                  only : kfl_windk_sld, ptota_sld
  use def_solidz,                  only : bodyf_sld, kfl_bodyf_sld
  use def_solidz,                  only : conbo_sld, contactbou_sld
  use mod_sld_csys,                only : sld_csys_build_jacrot
  use mod_sld_stent,               only : sld_calculation_contact_stent
  use def_master,                  only : ITASK_TURNON, ITASK_BEGSTE
  use def_master,                  only : ITASK_BEGITE, ITASK_INNITE
  use def_domain,                  only : ndime, npoin
  use def_domain,                  only : lpoty
  use mod_ker_space_time_function, only : ker_space_time_function
  use def_solidz,                  only : bvess_sld
  use def_solidz,                  only : kfl_conbc_sld, kfl_funno_sld
  use def_solidz,                  only : nfunc_sld, kfl_funtn_sld
  use def_solidz,                  only : kfl_windk_sld, ptota_sld
  use def_solidz,                  only : bodyf_sld, kfl_bodyf_sld
  use def_solidz,                  only : kfl_fixno_sld
  use def_solidz,                  only : kfl_contn_stent
  use mod_ker_functions,           only : ker_functions
  use mod_sld_stent,               only : sld_set_boundaries
  use mod_sld_csys,                only : sld_csys_build_jacrot

  implicit none

  integer(ip), intent(in) :: itask      !< What to do
  integer(ip)             :: ipoin,ibopo,ifunc,itype,idime
  real(rp)                :: dinew(ndime),fonew(ndime)

  select case (itask)

  case( ITASK_TURNON )

  case( ITASK_BEGSTE )

     !-------------------------------------------------------------------
     !
     ! Before a time step
     !
     !-------------------------------------------------------------------

 
     if ( INOTEMPTY ) then

        !
        ! Get rotation matrix
        !
        !call sld_csys_build_jacrot 

        !
        ! Update boundary conditions (Transient)
        !
        if ( kfl_conbc_sld == 0 ) then ! Non-constant (Transient)

           do ipoin = 1, npoin
              ibopo = lpoty(ipoin)

              dinew(:) = 0.0_rp
              fonew(:) = 0.0_rp
              ifunc    = kfl_funno_sld(ipoin)
              itype    = kfl_funtn_sld(ipoin)

              if ( itype == FUNCTION_MODULE ) then
                 !
                 ! Solidz Functions
                 !
                 ! Displacements
                 call sld_funbou( 2_ip, ipoin, dinew)
                 !
                 ! Forces
                 if ( kfl_bodyf_sld > 0 ) then
                    call sld_funbou( 3_ip, ipoin,fonew)
                    bodyf_sld(1:ndime) = fonew(1:ndime)
                 end if

              else if ( itype == FUNCTION_SPACE_TIME ) then

                 call ker_functions(ipoin,ifunc,itype,bvess_sld(:,ipoin,2),dinew)

              else
                 !
                 ! Codes
                 !
                 dinew(1:ndime) = bvess_sld(1:ndime,ipoin,ITER_AUX)*cutim/(timef-timei)

              end if
              !
              ! Update Dircihlet displacement
              !
              do idime = 1,ndime
                 if( kfl_fixno_sld(idime,ipoin) > 0 ) then
                    bvess_sld(idime,ipoin,ITER_K) = dinew(idime)
                 end if
              end do
              !
              ! STent computations
              !
 !             if (kfl_conta_stent /= 0_ip) then
                   if (kfl_funno_sld(ipoin) < 0 .and. kfl_contn_stent(ipoin)==1_ip) then
                      ifunc = -kfl_funno_sld(ipoin)
                          
                          call ker_space_time_function(&
                              ifunc,coord(1,ipoin),coord(2,ipoin)+displ(2,ipoin,TIME_N),coord(ndime,ipoin),cutim,dinew(1:ndime))
                          call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))

                       if (kfl_conta_stent == 3_ip) then 
                         if ( (coord(2,ipoin)+displ(2,ipoin,TIME_N) <= 0.00_rp) .and. (coord(2,ipoin)+displ(2,ipoin,TIME_N))>-4.0_rp ) then
                            call ker_space_time_function(&
                                ifunc,coord(1,ipoin),coord(2,ipoin)+displ(2,ipoin,TIME_N),coord(ndime,ipoin),cutim,dinew(1:ndime))
                            call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))
                         else if ( (coord(2,ipoin)+displ(2,ipoin,TIME_N))<=-4.0_rp ) then
                            call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))
                         end if
!                    else if (kfl_conta_stent == 4_ip) then
!                        if ( ((coord(2,ipoin) - cutim*12.5_rp) <= 0.0_rp) .and.  ((coord(2,ipoin) - cutim*12.5_rp) > -4.0_rp) ) then
!                           call ker_space_time_function(&
!                                ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,dinew(1:ndime))
!                           call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))
!                        else if ( ((coord(2,ipoin) - cutim*12.5_rp) <= -4.0_rp) ) then 
!                           call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(1:ndime)) 
!                        end if
                    end if
                   !  if ( ibopo > 0 ) then
                   !    call sld_calculation_contact_stent(ipoin,ibopo,ifunc,dinew(:))
                   !  else
                   !    call runend('WE ARE TRYING CONTACT WITHOUT A BOUNDARY IBOPO')
                   !  end if
 !               end if
             end if    
           end do
      
        end if

     end if

  case( ITASK_BEGITE )

     !-------------------------------------------------------------------
     !
     ! Before a global iteration
     !
     !-------------------------------------------------------------------

  case( ITASK_INNITE )

     !-------------------------------------------------------------------
     !
     ! Before an inner iteration
     !
     !-------------------------------------------------------------------

  end select

end subroutine sld_updbcs
