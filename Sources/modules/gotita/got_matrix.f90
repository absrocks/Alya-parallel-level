subroutine got_matrix(itask)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_matrix
  ! NAME 
  !    got_matrix
  ! DESCRIPTION
  !    This routine computes elemental matrix and RHS
  ! USES
  !    got_elmope
  ! USED BY
  !    got_solite
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_solver
  use def_gotita
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp

  if(kfl_paral/=0) then
     ! 
     ! Initializations
     !
     solve_sol => solve(ivari_got:)
     call inisol()
     tamin_got= 1.0e16
     tamax_got=-1.0e16
     dimin_got= 1.0e16
     dimax_got=-1.0e16
     !
     ! Smooth diffusion
     !
     call got_diffus()
     !
     ! Element assembly
     !
     if(kfl_probl_got==1) then
        call got_elmope(itask)
     else if(kfl_probl_got==2) then
        call got_elmop2()
     else if(kfl_probl_got==3) then
        call got_elmop3()
     end if
  end if

end subroutine got_matrix
