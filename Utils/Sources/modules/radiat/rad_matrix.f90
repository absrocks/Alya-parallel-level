subroutine rad_matrix()
  !------------------------------------------------------------------------
  !
  ! Compute elemental matrix and RHS of the following problem:
  !   _                             _
  !  |                             |
  !  |  rho*cp/(theta*dt) T*v dw + |  rho*cp*[u.grad(T)] dw
  ! _|W                           _|W
  !
  !    _                         _               _
  !   |                         |               |
  ! + |  k*grad(T).grad(v) dw + |  ar*Tr*v ds = |  Q*v dw
  !  _|W                       _|S             _|
  !
  !   _                               _
  !  |                               |
  !  |  rho*cp/(theta*dt) T^n*v dw - |  (qr-ar*Tr) v ds
  ! _|W                             _|S
  ! 
  ! All terms are calculated at current time step, i.e. n+theta
  !
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_radiat
  use def_domain

  implicit none

  !
  ! Initializations
  !
  call inisol()
  resgs_rad(1) = 0.0_rp
  resgs_rad(2) = 0.0_rp

  if( INOTMASTER ) then 
     !
     ! Element assembly
     !
     call rad_elmope(1_ip) 
     !
     ! Boundary assembly
     !
     call rad_bouope()
     !
     ! Immersed boundary
     !
     !call rad_bounib()

  end if

end subroutine rad_matrix

