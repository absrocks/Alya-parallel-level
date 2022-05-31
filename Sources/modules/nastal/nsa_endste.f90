subroutine nsa_endste
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_endste
  ! NAME 
  !    nsa_endste
  ! DESCRIPTION
  !    This routine ends the time step
  ! USES
  !    nsa_cvgunk
  !    nsa_updunk
  !    nsa_output
  ! USED BY
  !    Nastal
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_nastal
  implicit none
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if(kfl_stead_nsa==0.and.kfl_timei_nsa==1) then
     call nsa_cvgunk(three)
     call nsa_updunk(five)  
  end if
  !
  ! Compute averaged variables
  !
  call nsa_averag()
  !
  ! Write restart file
  !
  call nsa_restar(two)
  !
  ! If not steady, go on
  !
  if(kfl_stead_nsa==0) then
     if(kfl_timei_nsa==1) kfl_gotim = 1
  end if

  !
  ! If the initial velocity was zero, now is ok, so go back to preconditioning when needed
  !
  if (kfl_zevel_nsa(1) == 1) then
     kfl_zevel_nsa(1) = 0  
     kfl_lopre_nsa    = kfl_zevel_nsa(2)
  end if
  

end subroutine nsa_endste
