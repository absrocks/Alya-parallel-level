subroutine ada_outgid
!-----------------------------------------------------------------------
!****f* adapti/ada_outgid
! NAME 
!    ada_outgid
! DESCRIPTION
!    This routine writes the post.msh file for gid
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_inpout
  use      def_adapti
  implicit none
  integer(ip)              :: &
       lutmp,ielem,ipoin,inode,itrce,iolne,lnofa(mnode),jelem,jolne

  !
  ! Write again the original gid file, but in the adapti gid one
  !
  call ada_openfi(3)
  lutmp = lun_outpu_dom
  lun_outpu_dom = lun_ougid_ada
  !call geogid

  !
  ! Reopen it to append new mesh elements
  !
  call ada_openfi(4)
  
  if (nimmo_ada == 0) then
     !
     ! No immersed objects around here
     !
     write(lun_ougid_ada,'(a)') 'MESH refined_mesh dimension 2 Elemtype Triangle Nnode 3'
     write(lun_ougid_ada,'(a)') 'coordinates'  
     do ipoin= 1,npnew_ada
        write(lun_ougid_ada,*) ipoin+npoin, conew_ada(1:ndime,ipoin)
     end do
     write(lun_ougid_ada,'(a)') 'end coordinates'     
     write(lun_ougid_ada,'(a)') 'elements'
     itrce = 3
     do ielem = 1,nenew_ada
        write(lun_ougid_ada,*) ielem, (lnnew_ada(inode,ielem),inode=1,mnode), itrce     
     end do
     write(lun_ougid_ada,'(a)') 'end elements'

  else if (nimmo_ada > 0) then
     !
     ! Immersed objects are present
     !
     !
     ! Reopen it to append immersed body and sausage
     !
     
     write(lun_ougid_ada,'(a)') 'MESH immersed_body dimension 2 Elemtype Linear Nnode  2'
     write(lun_ougid_ada,'(a)') 'coordinates'  
     do ipoin= 1,nposu_ada
        write(lun_ougid_ada,*) ipoin+npoin, conew_ada(1:ndime,ipoin)
     end do
     write(lun_ougid_ada,'(a)') 'end coordinates'
     write(lun_ougid_ada,'(a)') 'elements'
     itrce = 3
     do ielem = 1,nelsu_ada
        write(lun_ougid_ada,*) ielem+nelem, (npoin+lnosu_ada(inode,ielem),inode=1,nnodb_ada), itrce     
     end do
     write(lun_ougid_ada,'(a)') 'end elements'
     write(lun_ougid_ada,'(a)') 'MESH sausage dimension 2 Elemtype Triangle Nnode  3'
     write(lun_ougid_ada,'(a)') 'coordinates'  
     write(lun_ougid_ada,'(a)') 'end coordinates'
     write(lun_ougid_ada,'(a)') 'elements'
     itrce = 6
     do ielem = 1,nelem
        if (lenew_ada(ielem)%neseg < 0) then
           write(lun_ougid_ada,*) ielem, (lnods(inode,ielem),inode=1,mnode), itrce     
           !        write(6,*) ielem, (lnods(inode,ielem),inode=1,mnode), itrce     
        end if
     end do
     write(lun_ougid_ada,'(a)') 'end elements'
     
     
     !
     ! write adapti domain output (without the sausage)
     !
     
     write(lun_oudom_ada,'(a)') 'MESH sausage dimension 2 Elemtype Triangle Nnode  3'
     write(lun_oudom_ada,    *) '# npoco_ada :',npoco_ada
     write(lun_oudom_ada,    *) '# nelco_ada :',nelco_ada
     write(lun_oudom_ada,    *) '# npopa_ada :',npopa_ada
     write(lun_oudom_ada,'(a)') 'coordinates'
     do iolne=1,npoin
        ipoin= lolne_ada(4,iolne)
        if (iolne <= npopa_ada) then
           jolne= lolne_ada(2,iolne)
           ipoin= lolne_ada(4,jolne)
        end if
        if (ipoin > 0) write (lun_oudom_ada,*) iolne, coord(1:ndime,ipoin)
     end do
     write(lun_oudom_ada,'(a)') 'end coordinates'
     write(lun_oudom_ada,'(a)') 'elements'
     itrce = 3
     jelem = 0
     do ielem=1,nelem
        if (lenew_ada(ielem)%neseg == 0) then
           do inode=1,mnode
              ipoin= lnods(inode,ielem)
              lnofa(inode)= lolne_ada(3,ipoin)
              if (lnofa(inode) <= npopa_ada) then
                 lnofa(inode) =  lolne_ada(1,lnofa(inode))
              end if
           end do
           jelem = jelem + 1
           write(lun_oudom_ada,*) jelem, (lnofa(inode),inode=1,mnode), itrce             
        end if
     end do
     write(lun_oudom_ada,'(a)') 'end elements'

  end if
  
     

end subroutine ada_outgid
