!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_output.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Output and postprocess of solution
!> @details Output and postprocess of solution
!> @} 
!------------------------------------------------------------------------
subroutine por_output()
  use def_parame
  use def_master
  use def_domain
  use def_porous
  use mod_iofile
  implicit none
  integer(ip) :: ivari,ivarp
  !
  ! Initial solution, end of a time step and and of run
  !
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     call por_outvar(ivari)
  end do

  !!call por_outwel()

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Calculations on sets
     !
     call por_outset()

  else if( ittyp == ITASK_ENDRUN ) then
     !
     ! End of the run
     !

  end if

end subroutine por_output
