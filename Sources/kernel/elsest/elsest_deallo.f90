subroutine elsest_deallo()
  !
  ! Deallocate memory for bin and quad/oct structures
  !
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip) :: istat
  !
  ! Deallocate Bin structure
  !
  if( kfl_memor(1)==1 ) then
     kfl_memor(1)=0
     deallocate(bin_struc,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'BIN_STRUC','elsest_deallo',0_ip)
  end if
  !
  ! Deallocate Quad/Oct structure
  !
  if( kfl_memor(2)==1 ) then
     kfl_memor(2)=0
     deallocate(oct_struc,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'OCT_STRUC','elsest_deallo',0_ip)
  end if
  !
  ! Cuurent
  !
  if( kfl_memor(3)==1 ) then
     kfl_memor(3)=0
     deallocate(current,stat=istat)
  end if


end subroutine elsest_deallo
