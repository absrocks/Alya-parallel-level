subroutine renum1(lpont,nodad,nnode,npoin,nelem,nposi)

!-----------------------------------------------------------------------
!
! Computes the connections between degrees of freedom (stored as 
! a linked list)
!
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip) :: nnode,npoin,nelem,nposi
  integer(ip) :: lpont(nnode,nelem), nodad(1:*)
  integer(ip) :: ipoin,jpoin,kpoin,inode,jnode,nfree
  integer(ip) :: nadro,nadre,ielem

  do ipoin = 1,npoin
     nodad(ipoin) = 0
  end do
  nfree = npoin+1
  
  do ielem = 1,nelem
     do inode = 1,nnode
        ipoin = lpont(inode,ielem)
        do jnode = 1,nnode
           if(jnode/=inode) then
              jpoin = lpont(jnode,ielem)
              nadre = nodad(ipoin)
              nadro = ipoin
              do while (nadre>0)
                 kpoin = nodad(nadre)
                 if(kpoin==jpoin) go to 10
                 nadro = nadre + 1
                 nadre = nodad(nadro)
              end do
              nodad(nadro) = nfree
              nodad(nfree) = jpoin
              nodad(nfree+1) = 0
              nfree = nfree + 2
10            continue
           end if
        end do
     end do
  end do
  nposi = nfree
  
end subroutine renum1
