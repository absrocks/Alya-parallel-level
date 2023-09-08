subroutine qua_matrix()
  !------------------------------------------------------------------------
  !
  ! Compute elemental matrix and RHS of the following problem:
  !
  !    _                                  _                    _
  !   |                                  |                    |
  ! + |  (-0.5)*grad(Phi).grad(Phi) dw + |  V(x,y,z)*Phi dw = |  Phi dw
  !  _|W                                _|W                  _|
  !
  ! with results into a matrix system in the way:
  !
  !      K*Phi = landa*M*Phi 
  !
  ! All terms are calculated at current time step, i.e. n+theta
  !
  !------------------------------------------------------------------------
  use def_master
  implicit none
  !
  ! Initializations
  !
  if( INOTMASTER ) then

     call inieig()
     !
     ! Element assembly
     !
     call qua_elmope(0_ip)

  endif

end subroutine qua_matrix
