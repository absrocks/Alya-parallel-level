subroutine chm_matrix()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_matrix
  ! NAME 
  !    chm_marix
  ! DESCRIPTION
  !    This routine assembles matrix and RHS
  ! USES
  !    chm_matrix
  !    Soldir
  !    Solite
  ! USED BY
  !    chm_solite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_chemic
  use def_domain
  use def_solver
  use mod_postpr
  use mod_messages, only : livinf
  implicit none
  real(rp) :: time1,time2,time3,dtinv_tmp,dtmin,xnorm
  !
  ! Solver initializations and preconditioners
  !
  call inisol()
  solve_sol(1) % xdiag = 1.0_rp/dtinv_chm
  call chm_massma(pmatr)
  !
  ! Save values
  !
  if( kfl_assem_chm >= 2 ) then
     dtinv_tmp = dtinv_chm
     dtinv_chm = 0.0_rp
  end if
  call cputim(time1)
  !
  ! Element assembly
  !
  if( INOTMASTER .and. kfl_assem_chm /= 4 ) then

     if( kfl_model_chm == 4 ) then
        !
        ! Combustion - Standard models
        !      
        call chm_elmcom(1_ip) 

     else if( kfl_model_chm == 5 ) then
        !
        ! CFI Combustion model (Univ. Twente)
        !      
        call chm_elmcfi(1_ip)

     end if
  else
     dtinv_chm = 1.0e6_rp
  end if

  if (kfl_reset /= -1) then
     !
     ! Look for minimum over subdomains
     !
     call pararr('MIN',0_ip,1_ip,dtcri_chm)
     !
     !
     if (dtcri_chm /= 0.0) dtinv_chm = 1.0_rp/(dtcri_chm*safet_chm)
     !
     ! If changing dtinv is necessary, activate reset mode
     if (dtinv_chm > reset_factor * dtinv) then
        dtinv     = dtinv_chm
        kfl_reset = 1
        if( INOTSLAVE ) call livinf(-12_ip,'REQUESTED RESET OF TIME STEP', 1_ip)
     endif
  endif
  !
  ! Boundary assembly
  !
  call cputim(time2)
  if( INOTMASTER .and. kfl_assem_chm /= 4 ) call chm_bouope() 
  call cputim(time3)
  !
  ! Recover old values
  !
  if( kfl_assem_chm >= 2 ) then 
     dtinv_chm = dtinv_tmp
  end if
  !
  ! Split assembly
  !
  if( solve_sol(1)%kfl_algso == 9 .or. &
      solve_sol(1)%kfl_algso == 10 ) then
     dtinv_tmp = dtinv_chm
     dtinv_chm = 0.0_rp     
  end if
  if( INOTMASTER ) call chm_splass()
  if( solve_sol(1)%kfl_algso == 9 .or. &
      solve_sol(1)%kfl_algso == 10 ) then
     dtinv_chm = dtinv_tmp
  end if

  cputi_chm(1) = cputi_chm(1) + (time2-time1)
  cputi_chm(2) = cputi_chm(2) + (time3-time2)

  cpu_modul(CPU_ASSEMBLY,modul) = cpu_modul(CPU_ASSEMBLY,modul) + time3 - time1

end subroutine chm_matrix

