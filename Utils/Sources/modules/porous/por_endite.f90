!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_endite.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Checks convergence and performs updates
!> @details Checks convergence and performs updates
!!    - itask=1 The end of an internal iteration
!!    - itask=2 The end of the internal loop iteration
!> @} 
!------------------------------------------------------------------------
subroutine por_endite(itask)
  use      def_parame
  use      def_master
  use      def_porous
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: itask

  select case(itask)

  case(1)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || P(n,i,j) - P(n,i,j-1)|| / ||P(n,i,j)||) and update unknowns:
     !  P(n,i,j-1) <-- P(n,i,j) 
     !
     call por_cvgunk(1_ip) ! Residual:   ||UNKNO(:)-PRESS(:,1)||
     call por_updunk(3_ip) ! Update:     PRESS(:,1)=UNKNO
     !
     ! Solve Subgrid scale equation
     !
!     call por_solsgs()   ! for the moment no sgs

  case(2)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || P&S(n,i,*) - P&S(n,i-1,*)|| / ||P&S(n,i,*)||) and update unknowns:
     !  P&S(n,i-1,*) <-- P&S(n,i,*) 
     !
     call livinf(16_ip,' ',itinn(modul))
     call por_cvgunk(2_ip) ! Residual: ||PRESS(:,2)-PRESS(:,1)||     !also SATUR
     call por_updunk(4_ip) ! Update:   PRESS(:,2) = PRESS(:,1)

  case(3)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || S(n,i,j) - S(n,i,j-1)|| / ||S(n,i,j)||) and update unknowns:
     !  S(n,i,j-1) <-- S(n,i,j) 
     !
     call por_cvgunk(4_ip) ! Residual:   ||UNKNO(:)-SATUR(:,1)|| 
     call por_updunk(7_ip) ! Update:     SATUR(:,1)=UNKNO
     !
     ! Solve Subgrid scale equation
     !
!     call por_solsgs()   ! for the moment no sgs

  end select

end subroutine por_endite
