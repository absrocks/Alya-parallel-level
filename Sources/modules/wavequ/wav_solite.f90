subroutine wav_solite
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_solite
  ! NAME 
  !    wav_solite
  ! DESCRIPTION
  !    This routine solves an iteration of the temperature equations.
  ! USES
  !    wav_matrix
  !    Soldir
  !    Solite
  ! USED BY
  !    wav_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  implicit none
  real(rp) :: cpu_refe1,cpu_refe2

  call cputim(cpu_refe1)
  !
  ! Update inner iteration counter and write headings in the solver
  ! file.
  !
  itinn(modul) = itinn(modul) + 1

  if(kfl_paral/=0) then
     !
     ! Construct the system matrix and right-hand-side
     !
     call wav_matrix()
     !
     ! Solve the algebraic system
     !  
     if(kfl_timet_wav==1) then
        call wav_explic()
     else
        call solver(rhsid,unkno,amatr,pmatr)
     end if
  end if

  call cputim(cpu_refe2)
  cpu_modul(3,modul) =  cpu_modul(3,modul) + cpu_refe2 - cpu_refe1 

end subroutine wav_solite
