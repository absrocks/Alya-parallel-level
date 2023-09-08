subroutine par_senbcs()
  !------------------------------------------------------------------------
  !****f* Parall/par_senbcs
  ! NAME
  !    par_senbcs
  ! DESCRIPTION
  !    Send boundary conditions to slaves  
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_parall
  use def_domain
  use def_master
  use mod_memchk
  implicit none  

  strin = 'par_senbcs'
  strre = 'par_senbcs'
  strch = 'par_senbcs'

  !----------------------------------------------------------------------
  !
  ! Codes
  !
  !----------------------------------------------------------------------
  !
  ! Slaves allocate memory
  !
  if( ISLAVE ) then
     nboun_2 = nboun
     call membcs(1_ip)
     call membcs(2_ip)
  end if
  !
  ! Gather code arrays
  !
  if( kfl_icodn > 0 ) call par_parari('GAT',NPOIN_TYPE,mcono*npoin,kfl_codno)
  if( kfl_icodb > 0 ) call par_parari('GAT',NBOUN_TYPE,nboun,kfl_codbo)
  !
  ! Master deallocates memory
  !
  if( IMASTER ) then
     call membcs(-1_ip)
     call membcs(-2_ip)
  end if

end subroutine par_senbcs
