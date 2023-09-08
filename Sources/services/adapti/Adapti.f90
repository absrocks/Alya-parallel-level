subroutine Adapti(order)
!-----------------------------------------------------------------------
!****f* adapti/Adapti
! NAME 
!    Adapti
! DESCRIPTION
!    This routine is the bridge for Adapti service
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  implicit none
  integer(ip), intent(in) :: order
  real(rp)                :: time1,time2

  servi=9
  call cputim(time1)
  
  select case (order)

  case(ITASK_REAPRO) 
     call ada_reapro 
  case(ITASK_TURNON)
     if(kfl_servi(servi)/=0) then
        call ada_turnon     
        call ada_modmsh(ITASK_TURNON)
     end if
  case(ITASK_ENDSTE)
     if(kfl_servi(servi)/=0) then
        call ada_modmsh(ITASK_ENDSTE)
        call ada_outmsh  
        call ada_outerr
        call ada_prject
     end if
  end select

  call cputim(time2)
  cpu_servi(1,servi) =  cpu_servi(1,servi) + time2-time1


end subroutine adapti
