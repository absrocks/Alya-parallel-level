!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_begite.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Starts an internal iteration for the porous equation
!> @details Starts an internal iteration for the porous equation
!> @} 
!------------------------------------------------------------------------
subroutine por_begite 
  use      def_parame
  use      def_master
  use      def_domain
  use      def_porous
  use mod_messages, only : livinf
  implicit none
  !
  ! Initializations
  !
  kfl_goite_por = 1 
  itinn(modul)  = 0
  if(itcou==1) call por_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Update boundary conditions
  !
  call por_updbcs(two)
  !
  ! Set up the solver parameters for the porous equation
  !
  ! call por_inisol()   ! por_inisol must not be called he but in some more inner subroutine where kprsa_por is already defined 
  !
  ! Obtain the initial guess for inner iterations
  !
  call por_updunk(two)
  !
  ! Coupling
  !
!  call por_coupli(ITASK_BEGITE)   !for the moment we do not have coupling in porous

end subroutine por_begite
    
