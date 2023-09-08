subroutine elsest_octdes(ithre)
  !------------------------------------------------------------------------
  !****f* elsest/elsest_octdes
  ! NAME 
  !    elsest_octdes
  ! DESCRIPTION
  !    Deallocate the quad/oct-tree structure
  ! USES
  ! USED BY
  !    elsest_octdea
  !***
  !------------------------------------------------------------------------
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip),intent(in) :: ithre
  integer(ip)            :: istat
  logical(lg)            :: conti

  conti=.true.
  do while(conti)
     !
     ! First go to deepest level in first branch
     !
     do while(current(ithre)%o%whoiam==0)
        current(ithre)%o=>current(ithre)%o%children(1)
     end do
     !
     ! Deallocate list of elements
     !
     if(current(ithre)%o%whoiam==1) then
        call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'CURRENT%ELEMS','elsest_octdes',current(ithre)%o%elems)
        deallocate(current(ithre)%o%elems,stat=istat) 
        if(istat/=0) call elsest_memerr(2_ip,'CURRENT%ELEMS','elsest_octdes',0_ip)
     end if

     if(current(ithre)%o%childid<divmax.and.current(ithre)%o%childid/=0) then
        !
        ! I'm not the last child neither the Padrino
        !
        current(ithre)%o => current(ithre)%o%parent%children(current(ithre)%o%childid+1)

     else if(current(ithre)%o%childid==divmax) then
        !
        ! I'm the last child
        !
        current(ithre)%o => current(ithre)%o%parent
        deallocate(current(ithre)%o%children,stat=istat)
        if(istat/=0) call elsest_memerr(2_ip,'CURRENT_CHILDREN','elsest_octdes',0_ip)
        current(ithre)%o%whoiam=3

     else if(current(ithre)%o%id==0) then
        !
        ! I'm the Padrino: end of deallocation
        !
        deallocate(current(ithre)%o,stat=istat)
        if(istat/=0) call elsest_memerr(2_ip,'CURRENT','elsest_octdes',0_ip)
        conti=.false.

     end if

  end do

end subroutine elsest_octdes
