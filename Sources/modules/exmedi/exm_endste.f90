subroutine exm_endste
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_endste
  ! NAME 
  !    exm_endste
  ! DESCRIPTION
  !    This routine ends a time step of the incompressible NS equations.
  ! USES
  !    exm_output
  ! USED BY
  !    exm_endste (itask=1)
  !    exm_turnof (itask=2)
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain

  use      def_exmedi

  implicit none


  if(kfl_timei_exm==1) then
     call exm_updunk(ITASK_ENDSTE)         
     call exm_upcell(ITASK_ENDSTE)         
     call exm_isochr( 1_ip)
  end if

  if(kfl_timei_exm==1) kfl_gotim = 1
  
end subroutine exm_endste
