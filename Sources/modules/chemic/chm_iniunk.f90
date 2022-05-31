subroutine chm_iniunk()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_iniunk
  ! NAME 
  !    chm_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the concentrations.
  ! USED BY
  !    chm_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ker_proper
  use mod_chm_finiteRate,            only : chm_getProp_finiteRate
  use mod_chm_droplets
  use mod_chm_sectional_soot_model,  only : chm_validation_test_ssm 
  implicit none
  integer(ip) :: iclas,ipoin,icomp

  if( kfl_rstar == 0 ) then

     !
     ! Load initial conditions
     !
     if (kfl_model_chm /= 4) then

        if(kfl_field_chm(1) /= 0_ip) call chm_fields()  ! Initialization by fields activated
        
        do ipoin = 1,npoin
           do iclas = 1,nclas_chm
              conce(ipoin,iclas,1) = bvess_chm(iclas,ipoin)
           end do
        end do
        !
        ! Apply Dirichlet boundary conditions on (:,1)
        !
        call chm_updbcs(ITASK_INIUNK)
        
        !
        ! Droplet identification and output
        !
        if ( kfl_droplet_id_chm /= 0 ) then
           call chm_droplet_id 
           call chm_droplet_output
        end if

        call chm_updunk(ITASK_INIUNK)

     end if

  end if

  !
  ! Update chemistry model
  !
  call chm_update_model()

  !
  ! Soot model validation test
  !
  if (kfl_soot_chm < 0) &
     call chm_validation_test_ssm()

end subroutine chm_iniunk
