subroutine elsest_octdea(ithre,imesh,lmesh)
  !
  ! Oct deallocation
  !
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), intent(in) :: ithre,imesh,lmesh
  integer(ip)             :: istat

  if(oct_struc(imesh)%iallo/=0) then

     !call elsest_octpoi(imesh)

     oct_struc(imesh)%iallo =  0
     current(ithre)%o       => oct_struc(imesh)%tree_root
     memor                  => oct_struc(imesh)%memor
     call elsest_octdes(ithre)
     !
     ! Deallocate ELCOD memory
     !
     call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'ELCOD','elsest_octdea',oct_struc(imesh)%elcod)
     deallocate(oct_struc(imesh)%elcod,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'ELCOD','elsest_octdea',0_ip)
     !
     ! Deallocate arrays for memory and cpu statistics
     !
     deallocate(oct_struc(imesh)%cputi,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'CPUTI','elsest_deallo',0_ip)
     deallocate(oct_struc(imesh)%memor,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'MEMOR','elsest_deallo',0_ip)
     deallocate(oct_struc(imesh)%kstat,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KSTAT','elsest_alloc',0_ip)
     deallocate(oct_struc(imesh)%ksear,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KSEAR','elsest_alloc',0_ip)
     deallocate(oct_struc(imesh)%kfirs,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KFIRS','elsest_alloc',0_ip)
     deallocate(oct_struc(imesh)%kseco,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KSECO','elsest_alloc',0_ip)

  end if

end subroutine elsest_octdea
