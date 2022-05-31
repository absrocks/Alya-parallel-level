!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_doiter.f90
!> @author  Solidz Team
!> @date
!> @brief   This routine controls the internal loop of the module.
!> @details
!>          \verbatim
!>          Iterations depend on the integration scheme:
!>          Implicit scheme
!>            - Iterates until convergence is achieved
!>          Explicit scheme
!>            - Perform only one iteration (no convergence)
!>          Runge-Kutta scheme
!>            - Perform only one iteration (no convergence)
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_doiter

  use def_kintyp,  only : ip
  use def_master,  only : ITASK_ENDINN, ITASK_ENDITE
  use def_domain,  only : kfl_elcoh
  use def_solidz,  only : kfl_rigid_sld
  use def_solidz,  only : kfl_stead_sld, kfl_goite_sld
  use def_solidz,  only : kfl_conta_sld, SLD_PDN_UNILATERAL
  use def_solidz,  only : release_nodes
  use mod_sld_rbo, only : sld_rbo_doiter
#ifdef COMMDOM
#if COMMDOM==2
  use mod_sld_pdn_contact, only : commdom_sld_nodes_release
#endif
#endif

  implicit none

  if( kfl_stead_sld == 0_ip ) then

     if ( kfl_rigid_sld == 0_ip ) then

        !-------------------------------------------------------------------
        !
        ! Deformable body
        !
        !-------------------------------------------------------------------
        !
        ! Begin iteration
        !
        call sld_begite()
        !
        ! Iterate
        !
        do while( kfl_goite_sld == 1_ip ) ! if not converged
           call sld_solite()
           call sld_endite(ITASK_ENDINN)
#ifdef COMMDOM
#if COMMDOM==2
           call commdom_sld_nodes_release()
#endif
#endif
        end do
        !
        ! End iteration
        !
        call sld_endite(ITASK_ENDITE)
        !
        ! Update traction cohesive elements (XFM/CZM)
        !
        if( kfl_elcoh /= 0 ) call sld_uptcoh()
        !
        ! PDN contact
        !
        if( kfl_conta_sld == SLD_PDN_UNILATERAL .and. associated(release_nodes) ) then
           deallocate(release_nodes)
        end if

     else if ( kfl_rigid_sld == 1_ip ) then

        !-------------------------------------------------------------------
        !
        ! Rigid body
        !
        !-------------------------------------------------------------------

        call sld_rbo_doiter()

     end if

  end if

end subroutine sld_doiter
