!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_endite.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Check convergence and updates unknowns
!> @details Check convergence and updates unknowns
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_endite(itask)
  use      def_master
  use      def_parame
  use      def_nastal
  use      def_solver
  use mod_messages, only : livinf
  implicit none
  integer(ip)      :: itask !> subroutine's calling point
  character(300)   :: messa
  integer(ip)      :: maxiter

  select case(itask)

  case(one)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||) and update unknowns:
     !  u(n,i,j-1) <-- u(n,i,j)
     !
     call nsa_cvgunk(one)

     call nsa_updunk(three)   !  u(,ITER_K) <-- unkno

!!!!! PROBANDO ELMOPERATIONS, QUE NO USA SETVAR:
     if (kfl_mod_elmop_nsa == 1)  call nsa_computephysical(one) 

     ! Compute mean velocities only if total pressure conditions are present
     if (nuinlet_nsa > 0) call nsa_evaluate_uinlet(1_ip,ITER_K)


!!$  PRUEBA PARA EL PSEUDO: NO CALCULAR EL DT Y USAR EN LAS ITERACIONES DE NEWTON EL DEL PASO DE TIEMPO
!!$     if (kfl_goite_nsa == 1 .and. kfl_pseud_nsa == 1) then 
!!$        call nsa_updtss(two)
!!$       ! if (dtinv > dtinv_nsa) call runend('NSA_ENDITE: THE PHYSICAL TIME IS SMALLER THAN THE PSEUDO-TIME')
!!$     end if

     call livinf(56_ip,' ',modul)
!     call livinf(16_ip,' ',itinn(modul))
     maxiter= solve_sol(1) % miter
!!     messa = &
!!          ' (SUBIT: '//trim(intost(itinn(modul)))//'/'//trim(intost(miinn_nsa))//' IT: '//trim(intost(iters))//'/'//trim(intost(maxiter))//')'
     messa = &
          ' (SUBIT: '//trim(intost(itinn(modul)))//'/'//trim(intost(miinn_nsa))//' IT: '//trim(intost(last_iters_nsa))//'/'//trim(intost(maxiter))//')'
     call livinf(-3_ip,messa,one)

  case(two)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||) and update unknowns:
     !  u(n,i-1,*) <-- u(n,i,*)
     !
!!     call livinf(16_ip,' ',itinn(modul))
     call nsa_cvgunk(two)
     call nsa_updunk(four)    !  u(,ITER_AUX) <-- u(,ITER_K)
     ! ****************************
     ! Tests for FSI force exchange
     ! call nsa_fsiexch(1_ip)
     ! call nsa_fsiexch(2_ip)
     ! ****************************     
  end select

end subroutine nsa_endite
