function elsest_resmat(arr_to_resize,currsize1,currsize2,desirsize1,desirsize2,ithre,memnum)
  !
  !  Resize a given 2D matrix
  !
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), pointer     :: arr_to_resize(:,:), elsest_resmat(:,:)
  integer(ip), intent(in)  :: currsize1,currsize2,desirsize1,desirsize2,ithre,memnum
  integer(ip)              :: i,j,istat
  allocate(elsest_resmat(desirsize1,desirsize2),stat=istat)
  !*OMP CRITICAL(memorycheck)
  call elsest_memchk(0_ip,ithre,istat,memor(memnum,ithre),'RESIZE_M','elsest_resmat',elsest_resmat)
  !*OMP END CRITICAL(memorycheck)
  do j=1,currsize2
     do i=1,currsize1
        elsest_resmat(i,j) = arr_to_resize(i,j)
     end do
  end do
  !*OMP CRITICAL(memorycheck)
  call elsest_memchk(2_ip,ithre,istat,memor(memnum,ithre),'RESIZE_M','elsest_resmat',arr_to_resize)
  !*OMP END CRITICAL(memorycheck)
  deallocate(arr_to_resize,stat=istat)
  if(istat/=0) call elsest_memerr(2_ip,'RESIZE_M','elsest_resmat',0_ip)
  arr_to_resize => elsest_resmat
  return

end function elsest_resmat
