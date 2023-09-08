subroutine rad_endite(itask)
!-----------------------------------------------------------------------
!****f* Radiat/rad_endite
! NAME 
!    rad_endite
! DESCRIPTION
!    This routine checks convergence and performs updates of the
!    radiation  at:
!    - itask=1 The end of an internal iteration
!    - itask=2 The end of the internal loop iteration
! USES
!    rad_cvgunk
!    rad_updunk
! USED BY
!    rad_doiter
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_radiat
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: itask

  select case(itask)

  case(1)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || G(n,i,j) - G(n,i,j-1)|| / ||G(n,i,j)||) and update unknowns:
     !  G(n,i,j-1) <-- G(n,i,j) 
     !
     call rad_cvgunk(1_ip) ! Residual:   ||UNKNO(:)-RADIA(:,1)||
     call rad_updunk(3_ip) ! Update:     RADIA(:,1)=UNKNO
     call rad_updunk(7_ip) ! Update RADSO
     !
     ! Solve Subgrid scale equation
     !
     call rad_solsgs()

  case(2)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || G(n,i,*) - G(n,i-1,*)|| / ||G(n,i,*)||) and update unknowns:
     !  G(n,i-1,*) <-- G(n,i,*) 
     !
     call livinf(16_ip,' ',itinn(modul))
     call rad_cvgunk(2_ip) ! Residual: ||RADAV(:,2)-RADAV(:,1)||
     call rad_updunk(4_ip) ! Update:   RADAV(:,2) = RADAV(:,1)


  end select

end subroutine rad_endite
