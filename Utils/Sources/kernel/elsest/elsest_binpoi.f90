subroutine elsest_binpoi(imesh)
  !
  ! Points to current mesh structure
  !
  use def_elsest
  implicit none
  integer(ip), intent(in) :: imesh

  nboxe => bin_struc(imesh)%nboxe
  nboxx => bin_struc(imesh)%nboxx
  dataf => bin_struc(imesh)%dataf
  iallo => bin_struc(imesh)%iallo
  lboel => bin_struc(imesh)%lboel
  pboel => bin_struc(imesh)%pboel
  tboel => bin_struc(imesh)%tboel
  memor => bin_struc(imesh)%memor
  memax => bin_struc(imesh)%memax
  kstat => bin_struc(imesh)%kstat
  ksear => bin_struc(imesh)%ksear
  kfirs => bin_struc(imesh)%kfirs
  kseco => bin_struc(imesh)%kseco
  comin => bin_struc(imesh)%comin
  comax => bin_struc(imesh)%comax
  delta => bin_struc(imesh)%delta
  elcod => bin_struc(imesh)%elcod
  cputi => bin_struc(imesh)%cputi

end subroutine elsest_binpoi
