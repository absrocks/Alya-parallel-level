subroutine alyafaces
  
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain 
  use mod_memchk
  use mod_htable
  implicit none

  integer(4)                :: istat
  integer(ip)               :: &
       ielty,pface,ielem

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! llamadas a cosas de alya
  
  call elmtyp

  if (.not.associated(ltype)) allocate(ltype(        nelem))
  if (.not.associated(lelch)) allocate(lelch(        nelem))
  if (.not.associated(lpoty)) allocate(lpoty(        npoin))

  if (associated(nepoi)) deallocate(nepoi)
  if (associated(pelpo)) deallocate(pelpo)
  if (associated(lelpo)) deallocate(lelpo)
  if (associated(lelfa)) deallocate(lelfa)
  if (associated(facel)) deallocate(facel)
  if (associated(lfacg)) deallocate(lfacg)
  if (associated(lelbf)) deallocate(lelbf)

  do ielem=1,nelem
     call fintyp(ndime,lnnod(ielem),ielty)
     lexis(ielty)= 1
     ltype(ielem)= ielty
  end do

  INOTMASTER = .true.
  INOTSLAVE  = .true.

  call cderda
  call connpo

  nelem_2 = nelem
  pelpo_2 => pelpo
  lelpo_2 => lelpo

  kfl_lface= 1
  kfl_lelbf= 0

  lpoty = 1   ! esto es necesario

  call lgface


end subroutine alyafaces
