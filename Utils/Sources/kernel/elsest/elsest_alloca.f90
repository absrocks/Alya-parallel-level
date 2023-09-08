subroutine elsest_alloca(itask,ipara)
  !
  ! Allocate memory for bin and/or quad/oct structures
  !
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), intent(in)    :: itask
  integer(ip), intent(inout) :: ipara(*)
  integer(ip)                :: nmesh,istat
#ifdef _OPENMPPPPPPPPP
  integer :: OMP_GET_MAX_THREADS
#endif

#ifdef _OPENMPPPPPPPPP
  nthre=OMP_GET_MAX_THREADS()
#else
  nthre=1
#endif
  !
  ! Allocate bin structure
  !
  if( (itask==0 .or. itask==1 ) .and. kfl_memor(1) == 0 ) then
     nmesh        = ipara(5)
     kfl_memor(1) = 1
     allocate(bin_struc(nmesh),stat=istat)
     if( istat /= 0 ) call elsest_memerr(0_ip,'BIN_STRUC','elsest_alloc',0_ip)
     bin_struc(1:nmesh) = bin_struc_init
  end if
  !
  ! Allocate oct structure
  !
  if( (itask==0 .or. itask==2 ) .and. kfl_memor(2) == 0 ) then
     nmesh        = ipara(5)
     kfl_memor(2) = 1
     allocate(oct_struc(nmesh),stat=istat)
     if( istat /= 0 ) call elsest_memerr(0_ip,'OCT_STRUC','elsest_alloc',0_ip)
     oct_struc(1:nmesh) = oct_struc_init
  end if
  !
  ! Current
  !
  if( kfl_memor(3) == 0 ) then
     kfl_memor(3) = 1
     allocate(current(nthre),stat=istat)
  end if

end subroutine elsest_alloca
