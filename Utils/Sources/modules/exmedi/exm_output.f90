subroutine exm_output()
!-----------------------------------------------------------------------
!****f* Exmedi/exm_output
! NAME 
!    exm_output
! DESCRIPTION
! USES
!    output
!    postpr
!    exm_outrep
! USED BY
!    exm_endste (itask=1)
!    exm_turnof (itask=2)
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain

  use      def_exmedi

  use      mod_postpr
  use      mod_iofile
 
  implicit none
  integer(ip) :: ivari,ivarp
  !
  ! Initial solution, end of a time step and and of run
  !
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     call exm_outvar(ivari)
  end do
  !
  ! At the end of a time step
  !
  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Postprocess on sets
     !
     call exm_outset()
     !
     ! Postprocess on witness points
     !
     call exm_outwit()

  end if

end subroutine exm_output
