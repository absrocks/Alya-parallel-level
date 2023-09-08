recursive subroutine elsest_octdep(itask,ithre,ndime,ipoin,ielem)
  !
  ! Descend down the hierarchy for post-process
  !
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), intent(in)    :: itask,ithre,ndime
  integer(ip), intent(inout) :: ipoin,ielem
  integer(ip)                :: inode,lnods(8)
  !
  ! First go to deepest level in first branch
  !
  do while(current(ithre)%o%whoiam==0)
     current(ithre)%o=>current(ithre)%o%children(1)
  end do
  !
  ! Current bin has elements
  !
  if(current(ithre)%o%whoiam/=0) then
     if(itask==1) then
        if(ndime==2) then
           ipoin = ipoin+1
           write(iunit(2),*) ipoin,current(ithre)%o%minc(1),current(ithre)%o%minc(2)
           ipoin = ipoin+1
           write(iunit(2),*) ipoin,current(ithre)%o%maxc(1),current(ithre)%o%minc(2)
           ipoin = ipoin+1
           write(iunit(2),*) ipoin,current(ithre)%o%maxc(1),current(ithre)%o%maxc(2)
           ipoin = ipoin+1
           write(iunit(2),*) ipoin,current(ithre)%o%minc(1),current(ithre)%o%maxc(2)
        end if
     else if(itask==2) then
        if(ndime==2) then
           inode = 0
           ielem = ielem+1
           ipoin = ipoin+1
           inode = inode+1
           lnods(inode)=ipoin
           ipoin = ipoin+1
           inode = inode+1
           lnods(inode)=ipoin
           ipoin = ipoin+1
           inode = inode+1
           lnods(inode)=ipoin
           ipoin = ipoin+1
           inode = inode+1
           lnods(inode)=ipoin
           write(iunit(2),'(6(1x,i9))')  ielem,lnods(1),lnods(2),lnods(3),lnods(4),current(ithre)%o%level
        end if
     else if(itask==3) then
        ielem = ielem+1
        write(iunit(3),*) ielem,current(ithre)%o%npoinbox
     else if(itask==4) then
        ielem = ielem+1
        write(iunit(3),*) ielem,current(ithre)%o%nelembox   
     else if(itask==5) then
        ielem = ielem+1
        write(iunit(3),*) ielem,current(ithre)%o%level
     end if
  end if

  if(current(ithre)%o%childid < divmax .and. current(ithre)%o%childid /=0) then
     !
     ! I'm not the last child neither the Padrino
     !
     current(ithre)%o => current(ithre)%o%parent%children(current(ithre)%o%childid+1)

  else if(current(ithre)%o%childid==divmax) then
     !
     ! I'm the last child of this generation: postprocess 
     !
     do while(current(ithre)%o%id>0)
        if(current(ithre)%o%parent%id == 0) then
           return
        else
           if(current(ithre)%o%parent%childid /=divmax) then
              current(ithre)%o => current(ithre)%o%parent%parent%children(current(ithre)%o%parent%childid+1)
              exit
           else 
              current(ithre)%o => current(ithre)%o%parent
           end if
        end if
     end do

  else if(current(ithre)%o%id==0) then
     !
     ! I'm the Padrino
     !
     return
  end if

  call elsest_octdep(itask,ithre,ndime,ipoin,ielem)

end subroutine elsest_octdep
