function elsest_resvec(arr_to_resize,currsize1,desirsize1,ithre,memnum)
  !
  !  Resize a given 1D vector
  !
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), pointer     :: arr_to_resize(:), elsest_resvec(:)
  integer(ip), intent(in)  :: currsize1,desirsize1,ithre,memnum
  integer(ip)              :: i,istat
  allocate(elsest_resvec(desirsize1),stat=istat)
!*OMP CRITICAL(memorycheck)
  call elsest_memchk(0_ip,ithre,istat,memor(memnum,ithre),'RESIZE_V','elsest_resvec',elsest_resvec)
!*OMP END CRITICAL(memorycheck)
     do i=1,currsize1
     elsest_resvec(i) = arr_to_resize(i)
     end do
!*OMP CRITICAL(memorycheck)
  call elsest_memchk(2_ip,ithre,istat,memor(memnum,ithre),'RESIZE_V','elsest_resvec',arr_to_resize)
!*OMP END CRITICAL(memorycheck)
  deallocate(arr_to_resize,stat=istat)
  if(istat/=0) call elsest_memerr(2_ip,'RESIZE_V','elsest_resvec',0_ip)
  arr_to_resize => elsest_resvec
  return
end function
!
