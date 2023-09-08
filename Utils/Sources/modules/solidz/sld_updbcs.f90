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
  use def_master,                  only : INOTMASTER
  use def_master,                  only : cutim, timei, timef
  use def_master,                  only : ITER_K, ITER_AUX, TIME_N
  use def_domain,                  only : ndime, npoin
  use def_domain,                  only : lpoty, coord
  use def_kermod,                  only : number_space_time_function
  use mod_ker_space_time_function, only : ker_space_time_function
  use def_solidz,                  only : bvess_sld
  use def_solidz,                  only : kfl_conbc_sld, kfl_funno_sld
  use def_solidz,                  only : nfunc_sld
  use def_solidz,                  only : kfl_local_sld
  use def_solidz,                  only : kfl_windk_sld, ptota_sld
  use def_solidz,                  only : bodyf_sld, kfl_bodyf_sld
  use mod_sld_csys,                only : sld_csys_build_jacrot

  implicit none

  integer(ip), intent(in) :: itask      !< What to do
  integer(ip)             :: ipoin,ibopo,ifunc
  real(rp)                :: dinew(ndime),fonew(ndime)

  select case (itask)

  case( 0_ip )

     !-------------------------------------------------------------------
     !
     ! At the beginning of the run
     !
     !-------------------------------------------------------------------

     if ( INOTMASTER ) then
        !
        ! Build rotation matrix for local axes prescription
        !
        if ( kfl_local_sld /= 0_ip ) call sld_csys_build_jacrot()

     end if

  case( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Before a time step
     !
     !-------------------------------------------------------------------

     if ( INOTMASTER ) then
        !
        ! Compute endo pressure
        !
        if ( nfunc_sld > 0 ) then
           do ifunc = 1, kfl_windk_sld ! how many cavities to compute the windkessel pressure
              call sld_funbou(1_ip,-ifunc,ptota_sld(ifunc))
           end do
        end if
        !
        ! Update boundary conditions (Transient)
        !
        if ( kfl_conbc_sld == 0 ) then ! Non-constant (Transient)

           do ipoin = 1, npoin
              ibopo = lpoty(ipoin)

              dinew(:) = 0.0_rp
              fonew(:) = 0.0_rp
              if ( kfl_funno_sld(ipoin) > 0 ) then
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

              else if ( kfl_funno_sld(ipoin) < 0 .and. number_space_time_function > 0 ) then
                 !
                 ! Space/Time Functions
                 !
                 ifunc = -kfl_funno_sld(ipoin)
                 call ker_space_time_function(&
                      ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,dinew(1:ndime))

              else
                 !
                 ! Codes
                 !
                 dinew(1:ndime) = bvess_sld(1:ndime,ipoin,ITER_AUX)*cutim/(timef-timei)

              end if
              !
              ! Update dirichlet displacement
              !
              bvess_sld(1:ndime,ipoin,ITER_K) = dinew(1:ndime)

           end do

        end if

     end if

  case( 2_ip )

     !-------------------------------------------------------------------
     !
     ! Before a global iteration
     !
     !-------------------------------------------------------------------

  case( 3_ip )

     !-------------------------------------------------------------------
     !
     ! Before an inner iteration
     !
     !-------------------------------------------------------------------

  end select

end subroutine sld_updbcs
