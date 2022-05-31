subroutine elsest_octpoi(imesh)
  !
  ! Points to current mesh structure
  !
  use def_elsest
  implicit none
  integer(ip), intent(in) :: imesh

  iallo     => oct_struc(imesh)%iallo
  iwhat     => oct_struc(imesh)%iwhat
  tree_root => oct_struc(imesh)%tree_root
  memor     => oct_struc(imesh)%memor
  limit     => oct_struc(imesh)%limit
  divmax    => oct_struc(imesh)%divmax
  memax     => oct_struc(imesh)%memax
  kstat     => oct_struc(imesh)%kstat
  ksear     => oct_struc(imesh)%ksear
  kfirs     => oct_struc(imesh)%kfirs
  kseco     => oct_struc(imesh)%kseco
  comin     => oct_struc(imesh)%comin
  comax     => oct_struc(imesh)%comax
  elcod     => oct_struc(imesh)%elcod
  lboel     => oct_struc(imesh)%lboel
  cputi     => oct_struc(imesh)%cputi

end subroutine elsest_octpoi
