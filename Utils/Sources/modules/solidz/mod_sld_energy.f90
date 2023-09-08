!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_energy.f90
!> @author  Gerard Guillamet
!> @date    February, 2018
!>           - Subroutine written
!> @brief   Calculation of energies and energetic balance checks
!> @details Calculation of energies and energetic balance checks
!>
!>          \verbatim
!>          The following energies are implemented:
!>           - ALLIE: All internal energy
!>           - ALLWK: All external work
!>           - ALLKE: All kinetic energy
!>           - ETOTA: Total energy
!>
!>          Energy balance:
!>           - Energetic balance according to Belytschko
!>          \endverbatim
!>
!>          References:\n
!>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
!>          Nonlinear Finite elements for Continua and Structures.
!>
!> @todo    To do list::\n
!>           - Energetic balance for subdomains (locally, see Belytschko)
!>           - Include energy from natural forces (fsi problems)
!>           - Include strain energy
!>           - Include artificial energy from rayleigh damping
!>
!> @}
!------------------------------------------------------------------------

module mod_sld_energy

  use def_kintyp,         only : ip, rp
  use def_domain,         only : ndime, npoin
  use def_master,         only : INOTMASTER, ITER_K, TIME_N
  use def_solidz,         only : fintt_sld, fextt_sld
  use def_solidz,         only : allie_sld, allwk_sld, allke_sld, etota_sld

  implicit none

  private

  public :: sld_energy
  public :: sld_updene

contains

  !-----------------------------------------------------------------------
  !>
  !> @brief   This routine calculates energies and perform a an energy
  !>          balance check.
  !> @details
  !>          \verbatim
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine sld_energy()

    use def_master,         only : ittim,rhsid
    use mod_operations,     only : operations_parallel_dot_product
    use def_solidz,         only : SLD_DYNAMIC_PROBLEM
    use def_solidz,         only : kfl_timei_sld
    use def_solidz,         only : veloc_sld, ddisp_sld
    use def_solidz,         only : vmass_sld
    use def_solidz,         only : eener_sld

    implicit none

    integer(ip)                 :: ipoin,idime
    real(rp)                    :: allwk, allie, allke          ! Energies
    real(rp)                    :: eener
    real(rp)                    :: fintsum(ndime,npoin), fextsum(ndime,npoin)
    real(rp)                    :: rhstota(ndime,npoin)
    real(rp)                    :: velosqu(ndime,npoin)
    real(rp)                    :: vmassto(ndime,npoin)

    !
    ! Preliminary calculations
    !
    if ( INOTMASTER ) then

       ! Forces summation between current and previous iteration
       fintsum(1:ndime,1:npoin) = fintt_sld(1:ndime,1:npoin,ITER_K) + fintt_sld(1:ndime,1:npoin,TIME_N)
       fextsum(1:ndime,1:npoin) = fextt_sld(1:ndime,1:npoin,ITER_K) + fextt_sld(1:ndime,1:npoin,TIME_N)
       ! Square of the velocity and lumped mass matrix
       if (kfl_timei_sld == SLD_DYNAMIC_PROBLEM) then
          ! Square of velocity
          velosqu(1:ndime,1:npoin) = veloc_sld(1:ndime,1:npoin,ITER_K)*veloc_sld(1:ndime,1:npoin,ITER_K)
          ! Lumped mass matrix
          do ipoin = 1,npoin
             do idime=1,ndime
                vmassto(idime,ipoin) = vmass_sld(ipoin)
             end do
          end do
       end if
       ! rhs in array form
       rhstota = reshape(rhsid, (/ndime,npoin/))

    end if

    !
    ! All energies (whole model)
    !
    ! All internal energy
    call operations_parallel_dot_product(ndime,fintsum,ddisp_sld(:,:,ITER_K),allie,'IN MY CODE')
    if (ittim == 1_ip) then
       allie_sld(ITER_K) = allie_sld(TIME_N) + abs(allie)
    else
       allie_sld(ITER_K) = allie_sld(TIME_N) + 0.5_rp*abs(allie)
    end if
    ! All work (external energy)
    call operations_parallel_dot_product(ndime,fextsum,ddisp_sld(:,:,ITER_K),allwk,'IN MY CODE')
    if (ittim == 1_ip) then
       allwk_sld(ITER_K) = allwk_sld(TIME_N) + abs(allwk)
    else
       allwk_sld(ITER_K) = allwk_sld(TIME_N) + 0.5_rp*abs(allwk)
    end if
    ! All kinetic energy
    if (kfl_timei_sld == SLD_DYNAMIC_PROBLEM) then
       call operations_parallel_dot_product(ndime,velosqu,vmassto,allke,'IN MY CODE')
       allke_sld(ITER_K) = 0.5_rp*allke
    end if

    !
    ! Total energy
    !
    etota_sld(ITER_K) = allie_sld(ITER_K) + allke_sld(ITER_K) - allwk_sld(ITER_K)

    !
    ! Error in energy flow
    !
    call operations_parallel_dot_product(ndime,rhstota,ddisp_sld(:,:,ITER_K),eener,'IN MY CODE')
    eener_sld = abs(eener)

    !
    ! Energy balance
    !
    if (etota_sld(ITER_K) <= 1.0e-2_rp*max(allie_sld(ITER_K),allwk_sld(ITER_K),allke_sld(ITER_K)) ) then
       !print*,'Energy balance OK'
    end if

  end subroutine sld_energy

  !-----------------------------------------------------------------------
  !>
  !> @brief   This routine performs several updates for energies and
  !>          variables related to them.
  !> @details
  !>          \verbatim
  !>          ITASK = 1 ... Initial guess outer iterations
  !>                  2 ... Initial guess inner iterations
  !>                  3 ... Updates (end interation)
  !>                  4 ... Updates (end iteration converged)
  !>                  5 ... End time step
  !>                  6 ... Initial guess for restart
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine sld_updene(itask)

    use def_master,           only : ITER_AUX
    use def_solidz,           only : nprev_sld

    implicit none

    integer(ip), intent(in)       :: itask !< Update variables at selected case

    select case (itask)

    case(1_ip)

       !------------------------------------------------------------------
       !
       !  Outer iterations (begin time step)
       !
       !------------------------------------------------------------------
       ! Force vectors
       if (INOTMASTER) then
          fintt_sld(1:ndime,1:npoin,ITER_AUX) = fintt_sld(1:ndime,1:npoin,nprev_sld)
          fextt_sld(1:ndime,1:npoin,ITER_AUX) = fextt_sld(1:ndime,1:npoin,nprev_sld)
       end if
       ! Energies
       allie_sld(ITER_AUX) = allie_sld(nprev_sld)
       allwk_sld(ITER_AUX) = allwk_sld(nprev_sld)
       allke_sld(ITER_AUX) = allke_sld(nprev_sld)
       etota_sld(ITER_AUX) = etota_sld(nprev_sld)

    case(2_ip)

       !------------------------------------------------------------------
       !
       !  Inner iterations (begin iterations)
       !
       !------------------------------------------------------------------
       ! Force vectors
       if (INOTMASTER) then
          fintt_sld(1:ndime,1:npoin,ITER_K) = fintt_sld(1:ndime,1:npoin,ITER_AUX)
          fextt_sld(1:ndime,1:npoin,ITER_K) = fextt_sld(1:ndime,1:npoin,ITER_AUX)
       end if
       ! Energies
       allie_sld(ITER_K) = allie_sld(ITER_AUX)
       allwk_sld(ITER_K) = allwk_sld(ITER_AUX)
       allke_sld(ITER_K) = allke_sld(ITER_AUX)
       etota_sld(ITER_K) = etota_sld(ITER_AUX)

    case(3_ip)

    case(4_ip)

       !------------------------------------------------------------------
       !
       !  Updates (end iterations converged)
       !
       !------------------------------------------------------------------
       ! Force vectors
       if (INOTMASTER) then
          fintt_sld(1:ndime,1:npoin,ITER_AUX) = fintt_sld(1:ndime,1:npoin,ITER_K)
          fextt_sld(1:ndime,1:npoin,ITER_AUX) = fextt_sld(1:ndime,1:npoin,ITER_K)
       end if
       ! Energies
       allie_sld(ITER_AUX) = allie_sld(ITER_K)
       allwk_sld(ITER_AUX) = allwk_sld(ITER_K)
       allke_sld(ITER_AUX) = allke_sld(ITER_K)
       etota_sld(ITER_AUX) = etota_sld(ITER_K)

    case(5_ip)

       !------------------------------------------------------------------
       !
       !  End Time Step
       !
       !------------------------------------------------------------------
       ! Force vectors
       if (INOTMASTER) then
          fintt_sld(1:ndime,1:npoin,TIME_N) = fintt_sld(1:ndime,1:npoin,ITER_K)
          fextt_sld(1:ndime,1:npoin,TIME_N) = fextt_sld(1:ndime,1:npoin,ITER_K)
       end if
       ! Energies
       allie_sld(TIME_N) = allie_sld(ITER_K)
       allwk_sld(TIME_N) = allwk_sld(ITER_K)
       allke_sld(TIME_N) = allke_sld(ITER_K)
       etota_sld(TIME_N) = etota_sld(ITER_K)

     case(6_ip)

       !------------------------------------------------------------------
       !
       !  Initial guess after reading restart
       !
       !------------------------------------------------------------------
       ! Force vectors
       if (INOTMASTER) then
          fintt_sld(1:ndime,1:npoin,ITER_K) = fintt_sld(1:ndime,1:npoin,nprev_sld)
          fextt_sld(1:ndime,1:npoin,ITER_K) = fextt_sld(1:ndime,1:npoin,nprev_sld)
       end if
       ! Energies
       allie_sld(ITER_K) = allie_sld(nprev_sld)
       allwk_sld(ITER_K) = allwk_sld(nprev_sld)
       allke_sld(ITER_K) = allke_sld(nprev_sld)
       etota_sld(ITER_K) = etota_sld(nprev_sld)

    end select

  end subroutine sld_updene

end module mod_sld_energy

