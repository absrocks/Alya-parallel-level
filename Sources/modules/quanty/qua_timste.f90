subroutine qua_timste()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_timste
  ! NAME 
  !    qua_timste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the shrodinger
  !    equation      
  ! USES
  !    qua_updtss
  !    qua_outlat
  !    qua_updtss
  ! USED BY
  !    Quanty
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  implicit none

  if(ittim==0) then
     
	 call qua_output(-1_ip)
	 
	 ! Aca generola sol inicial 
	 ! desde las phi que yo le indique en la entrada
     ! pero solo cuando ittim==0 La guarda en rhoon(*,ncomp_qua) ???
     
	 ! a partir de N y L
     ! calcular degeneracion
	 !         multiplicidad
	 ! phion teoricas con Z y tipo de atomo
	 !
	 ! y finalmente rhoon

  endif

  ! Time step size 
  !
!  if(kfl_stead_qua/=1) call qua_updtss()     ! fijo
 
end subroutine qua_timste

