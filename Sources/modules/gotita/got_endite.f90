subroutine got_endite(itask)
!-----------------------------------------------------------------------
!****f* Gotita/got_endite
! NAME 
!    got_endite
! DESCRIPTION
!    This routine checks convergence and updates unknowns at:
!    - itask=1 The end of an internal iteration
!    - itask=2 The end of the internal loop iteration
! USES
!    got_cvgunk
!    got_updunk
! USED BY
!    got_doiter
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_gotita
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: itask

  select case(itask)

  case(one)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||) and update unknowns:
     !  u(n,i,j-1) <-- u(n,i,j)
     !
     call got_cvgunk(  one)
     call got_updunk(three)
     !
     ! Output solution
     !
     call got_output(  two)

  case(two)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||) and update unknowns:
     !  u(n,i-1,*) <-- u(n,i,*)
     !
     call livinf(16_ip,' ',itinn(modul))
     call got_cvgunk( two)
     call got_updunk(four)
     
  end select

end subroutine got_endite
