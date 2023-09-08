subroutine nsa_solite
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_solite
  ! NAME 
  !    nsa_solite
  ! DESCRIPTION
  !    This routine is the bridge to the different iteration solution schemes
  ! USES
  !    nsa_...
  ! USED BY
  !    nsa_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastal
  use mod_nsa_euler, only: EULER, nsa_euler_calc
  implicit none
  real(rp)    :: cpu_refe1,cpu_refe2

  call cputim(cpu_refe1)  
  !
  ! Update inner iteration counter and write headings in the solver file.
  !
  itinn(modul) = itinn(modul) + 1

  !  
  ! Clean rshax the first itinn
  !
  call nsa_uptimi(one)

  !
  ! Update boundary conditions
  !
  if(euler_nsa) then  
     call nsa_euler_calc( EULER )
  else
     call nsa_updbcs(three)
     !
     ! Time advance 
     !
     call nsa_gocons 
  endif

  call cputim(cpu_refe2)
  cpu_modul(3,modul) =  cpu_modul(3,modul) + cpu_refe2 - cpu_refe1 

end subroutine nsa_solite
