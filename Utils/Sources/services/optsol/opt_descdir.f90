subroutine opt_descdir
  !------------------------------------------------------------------------
  !****f* Optsol/opt_descdir
  ! NAME
  !    opt_descdir
  ! DESCRIPTION
  ! OUTPUT
  ! USED BY
  !    Optsol
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_optsol
  use      def_inpout
  use      def_solver
  implicit none

  integer(ip) :: indvars

  ! Currently, we have implemented only the Steepest descent direction
  if(kfl_curlin_opt==1) then
     do indvars=1,kfl_ndvars_opt
        descdir(indvars)= 0.0_rp  
        if (diffj_illum(indvars)/=0.0_rp)then
           descdir(indvars)= -1.0_rp * diffj(indvars) 
        end if
     end do
  end if

end subroutine opt_descdir
