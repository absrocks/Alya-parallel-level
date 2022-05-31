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
  use      mod_output_postprocess
 
  implicit none
  external    :: exm_outvar
  integer(ip) :: ivari,ivarp
  !
  ! Initial solution, end of a time step and and of run
  !
  call output_postprocess_variables(exm_outvar)
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
