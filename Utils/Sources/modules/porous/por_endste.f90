!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_endite.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Ends a time step of the porous equation.
!> @details Ends a time step of the porous equation.
!> @} 
!------------------------------------------------------------------------
subroutine por_endste()
  use def_parame
  use def_master
  use def_porous
  implicit none
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if(kfl_stead_por==0.and.kfl_timei_por(1)==1.and.kfl_timei_por(2)==1) then    ! missing check if pressure or satur here
!     call por_cvgunk(three)
     call por_updunk( five)          ! for Crank-Nicholson or BDF2
  end if
  !
  ! Write restart file
  !
!  call por_restar(2_ip)  ! restart not ready
  !
  ! If not steady, go on
  !
!  if(kfl_stead_por==0.and.kfl_timei_por==1.and.kfl_conve(modul)==1) 
   kfl_gotim = 1  ! for the moment the time step loop never stops!!!!!!!!!!!!!!!!

end subroutine por_endste
