subroutine nsa_timste
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_timste
  ! NAME 
  !    nsa_timste
  ! DESCRIPTION
  !    This routine computes the time step
  ! USES
  !    nsa_iniunk
  !    nsa_updtss
  !    nsa_updbcs
  !    nsa_updunk
  ! USED BY
  !    Nastal
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_nastal
  implicit none

  !
  ! Time step size. 
  !
  if(kfl_stead_nsa/=1) then
     if (kfl_nacdr_nsa == 0) then
        call nsa_updtss(DT_PHYSICAL)
     else
        ! explicit convection-diffusion-reaction equation
        call nsa_updtss_cdr
     end if
  end if

end subroutine nsa_timste
