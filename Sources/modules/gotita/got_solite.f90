subroutine got_solite()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_solite
  ! NAME 
  !    got_solite
  ! DESCRIPTION
  !    This routine solves an iteration for the incomcdropible NS equations
  !    using a monolitic scheme.
  ! USES
  !    got_ifconf
  !    got_matrix
  !    Soldir
  !    Solite
  !    got_rotunk
  ! USED BY
  !    got_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  implicit none
  integer(ip) :: kfl_linea_old
  real(rp)    :: cpu_refe1,cpu_refe2

  call cputim(cpu_refe1)
  !
  ! Update inner iteration counter and write headings in the solver file.
  !
  itinn(modul) = itinn(modul) + 1
  ittot_got    = ittot_got + 1
  !
  ! Linearization
  !
  kfl_linea_old=kfl_linea_got
  if(ittot_got<=npica_got) kfl_linea_got=1
  !
  ! Update bc for alpha
  !
  call got_updfix()
  !
  ! Solve system
  !
  if(kfl_algor_got==4) then
     !
     ! Block Gauss-Seidel scheme
     !
     call got_solbgs()
  else
     !
     ! Monolithic scheme
     !
     call got_solmon()
  end if
  !
  ! Linearization
  !
  kfl_linea_got=kfl_linea_old

  call cputim(cpu_refe2)
  cpu_modul(3,modul) =  cpu_modul(3,modul) + cpu_refe2 - cpu_refe1 

end subroutine got_solite
