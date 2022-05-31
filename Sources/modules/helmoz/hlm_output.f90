subroutine hlm_output()

  !-----------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_output.f90
  ! NAME 
  !    hlm_output
  ! DESCRIPTION
  !    This routine outputs and postprocesses a solution.
  !    Also, it calculates values of field vectors from FE-computed potentials.
  ! USES
  !    hlm_memmas
  !    hlm_outvar
  !    hlm_fivecs
  ! USED BY
  !    hlm_iniunk
  !    hlm_endite
  !    hlm_endste
  !-----------------------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain
  use def_helmoz
  use mod_postpr
  use mod_iofile

  implicit none

  integer(ip) :: ivari,ivarp


  !if(ittyp == ITASK_ENDRUN) then

  !This routine is called three times from different routines and field values should be calculated only once,
  !after FE computation of potentials (cntrl_hlm == 1_ip)
  cntrl_hlm = 0_ip

  !Master allocates memory for problem unknowns
  !if (IMASTER) call hlm_memmas()       

  !Output and postprocess solution of FE computation
  do ivarp = 1,nvarp
     ivari = ivarp 
     call posdef(11_ip,ivari) 
     call hlm_outvar(ivari)      
  enddo

  !Calculate values of field vectors from FE-computed potentials
  if (cntrl_hlm == 1_ip .and. INOTSLAVE) call hlm_fivecs()       

  !end if

end subroutine hlm_output
