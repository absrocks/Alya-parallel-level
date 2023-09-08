subroutine chm_begite()
  !-----------------------------------------------------------------------
  !****f* partis/chm_begite
  ! NAME 
  !    chm_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the transport
  !    equation
  ! USES
  !    chm_tittim
  !    chm_updbcs
  !    chm_inisol
  !    chm_updunk
  ! USED BY
  !    chm_doiter
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_chemic
  use def_solver
  use def_kermod
  use mod_ker_proper 
  use mod_messages, only : livinf
  implicit none
  !
  ! Initializations
  !
  kfl_goite_chm = 1 
  itinn(modul)  = 0
  kfl_under_chm = 0
  kfl_overs_chm = 0

  if(itcou==1) call chm_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Update boundary conditions
  !
  call chm_updbcs(2_ip)
  !
  ! In combustion, update initial values of specific heat and viscosity with temperature coming from outside
  !
  if (kfl_model_chm == 4 ) then
     call chm_upcpmu(1_ip) 
     if (kfl_stagg_chm /= 0) then 
        call livinf(164_ip,' ',1_ip)
        if (kfl_stagg_chm == 1) then 
           call livinf(59_ip,'EVOLVING CHEMICAL SOURCE TERMS',7_ip)
        else
           call livinf(59_ip,'EVOLVING ALL SOURCE TERMS',7_ip)
        endif
        call chm_updrea(1.0_rp/dtinv) ! Evolve reactions, if staggered scheme
        call livinf(15_ip,' ',modul)
     endif

  endif
  !
  ! Obtain the initial guess for inner iterations
  !
  call chm_updunk(2_ip)
  !
  !
  if (kfl_model_chm == 4 ) then
     call chm_upwmea(3_ip)            ! Mean molecular weight
     call ker_updpro()                ! We need to update the global viscosity,density, and Cp
     call chm_omegak(1_ip,nspec_chm)  ! Update all the mass source terms in combustion
  elseif (kfl_model_chm == 5 ) then
     call chm_upwmea(3_ip)                          ! wmean(:,2) => wmean(:,1) 
  endif
  
end subroutine chm_begite
    
