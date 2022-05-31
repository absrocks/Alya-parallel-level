!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_begite.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Starts each newton (inner) iteration 
!> @details Starts each newton (inner) iteration
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_begite
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  use mod_messages, only : livinf

  implicit none
  !
  ! Initializations
  !
  kfl_goite_nsa = 1
  itinn(modul)  = 0
  if(miinn_nsa==0) kfl_goite_nsa=0 
  if(itcou==1)  call nsa_tistep  
  !  call livinf(15_ip,' ',modul)
  !
  ! 
  !
!!$  PRUEBA PARA EL PSEUDO: NO CALCULAR EL DT Y USAR EN LAS ITERACIONES DE NEWTON EL DEL PASO DE TIEMPO
!!$  if (kfl_pseud_nsa == 1) then
!!$     call nsa_updtss(two)
!!$     ! dtinv_nsa = dtinv_nsa * safet_nsa / safet_pseud_nsa
!!$  end if
  ! 
  ! Update boundary conditions
  !
  call nsa_updbcs(two)
  !
  ! Obtain the initial guess for inner iterations: y(2) <- y(1)
  !
  call nsa_updunk(two)  ! u(,ITER_K) <-- u(,ITER_AUX)


end subroutine nsa_begite
