subroutine qua_endite(itask)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_endite
  ! NAME 
  !    qua_endite
  ! DESCRIPTION
  !    This routine checks convergence and performs updates of the
  !    solution  at:
  !    - itask=1 The end of an internal iteration
  !    - itask=2 The end of the internal loop iteration
  ! USES
  !    qua_cvgunk
  !    qua_updunk
  !    qua_output
  ! USED BY
  !    qua_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_quanty
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: itask

  select case(itask)

  case(1_ip)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || T(n,i,j) - T(n,i,j-1)|| / ||T(n,i,j)||) and update unknowns:
     !  T(n,i,j-1) <-- T(n,i,j) 
     !
     if(kfl_dftgs_qua==1 .or. kfl_alele_qua ==1) then
        call qua_updunk(4_ip)  ! calcula nueva rho(2) en base a los autoestados solucion
     endif

     call qua_cvgunk(1_ip) ! Residuo 

     if(kfl_dftgs_qua==1 .or. kfl_alele_qua ==1) then
        call qua_updunk(3_ip)  ! actualizo rho(1)
        call qua_integr(1_ip)
     end if

  case(2_ip)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || T(n,i,*) - T(n,i-1,*)|| / ||T(n,i,*)||) and update unknowns:
     !  T(n,i-1,*) <-- T(n,i,*) 
     !
     call livinf(16_ip,' ',itinn(modul))
     call qua_cvgunk(2_ip) ! Residual: ||TEMPE(:,2)-TEMPE(:,1)||
     !     call qua_updunk(NO_HAY) ! Update:   rho
     call qua_output(2_ip)   ! sale el calculo de energia!

  end select

end subroutine qua_endite
