!> opt_parall.f90
!> @file opt_parall.f90 
!> @fn opt_parall 
!> Text automatically added for recognition of doxygen parser
!>
subroutine opt_parall(itask)
  !-----------------------------------------------------------------------
  !****f* Optsol/opt_parall
  ! NAME
  !    opt_parall
  ! DESCRIPTION
  !    This routine is a bridge to Parall service  
  ! USED BY
  !    opt_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  implicit none
  integer(ip), intent(in) :: itask

  if(kfl_paral>=0) then

     select case(itask)

     case(-1)
        !
        ! Exchange kfl_ndvars_opt (number of design variables) always using MPI
        !
        call opt_sendat(-1_ip)

     case(1)
        !
        ! Exchange data read in opt_reanut always using MPI
        !
        call opt_sendat(1_ip)

     end select

  end if

end subroutine opt_parall
