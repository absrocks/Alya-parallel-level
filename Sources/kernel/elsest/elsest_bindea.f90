subroutine elsest_bindea(ithre,imesh,lmesh)
  !
  ! Bin deallocation
  !
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), intent(in) :: ithre,imesh,lmesh
  integer(ip)             :: istat

  if(imesh/=lmesh) call elsest_binpoi(imesh)
  !
  ! It has not been allocated
  !
  if(  bin_struc(imesh)%iallo == 0 ) return

  bin_struc(imesh)%iallo=0
  if(bin_struc(imesh)%dataf==0) then
     call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'TBOEL','elsest_bindea',bin_struc(imesh)%tboel)
     deallocate(bin_struc(imesh)%tboel,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'TBOEL','elsest_bindea',0_ip)
  else if(bin_struc(imesh)%dataf==1) then
     call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'PBOEL','elsest_bindea',bin_struc(imesh)%pboel)
     deallocate(bin_struc(imesh)%pboel,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'PBOEL','elsest_bindea',0_ip)
  end if
  call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'ELCOD','elsest_bindea',bin_struc(imesh)%elcod)
  deallocate(bin_struc(imesh)%elcod,stat=istat)
  if(istat/=0) call elsest_memerr(2_ip,'ELCOD','elsest_bindea',0_ip)
  !
  ! Deallocate arrays for memory and cpu statistics
  !
  if(associated(bin_struc(imesh)%cputi)) then
     deallocate(bin_struc(imesh)%cputi,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'CPUTI','elsest_deallo',0_ip)
  end if

  if(associated(bin_struc(imesh)%memor)) then
     deallocate(bin_struc(imesh)%memor,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'MEMOR','elsest_deallo',0_ip)
  end if

  if(associated(bin_struc(imesh)%kstat)) then
     deallocate(bin_struc(imesh)%kstat,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KSTAT','elsest_alloc',0_ip)
  end if

  if(associated(bin_struc(imesh)%ksear)) then
     deallocate(bin_struc(imesh)%ksear,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KSEAR','elsest_alloc',0_ip)
  end if

  if(associated(bin_struc(imesh)%kfirs)) then
     deallocate(bin_struc(imesh)%kfirs,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KFIRS','elsest_alloc',0_ip)
  end if

  if(associated(bin_struc(imesh)%kseco)) then
     deallocate(bin_struc(imesh)%kseco,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KSECO','elsest_alloc',0_ip)
  end if

end subroutine elsest_bindea
