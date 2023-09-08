subroutine ada_outmgi
!-----------------------------------------------------------------------
!****f* adapti/ada_outmgi
! NAME 
!    ada_outmgi
! DESCRIPTION
!    This routine dumps the file with the HSMB and the immersed boundary
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_adapti
  implicit none
  integer(ip)   :: &
       ielem,ipoin,iolne,jolne,ielsu,ielad,inode
  integer(ip)   :: &
       lnofa(mnode)
  character(10)       :: &
       mycha

  
  if (ndime == 2) then
     !
     ! ONLY 2D IS PROGRAMMED!!!!!!!!!
     !
    
     ! Triangle mesh generator format:
     !

     fil_oumgi_ada = &
          adjustl(trim(fil_oumgi_ada))//'.poly'

     call ada_openfi(2)

     !
     ! write mesh generator input file
     !

     write(lun_oumgi_ada,*) '# File produced by adapti service'
     write(lun_oumgi_ada,*) '# Use these options: '
     write(lun_oumgi_ada,*) '#     triangle -pPBYq32.5aXXX ',adjustl(trim(fil_oumgi_ada))
     write(lun_oumgi_ada,*) '#     where XXX is the element size:',sizet_ada
     write(lun_oumgi_ada,*) '# File produced by adapti service -pPBYq32.5'
     write(lun_oumgi_ada,*) nposu_ada, ' 2  0  0 '
     do iolne=1,npopa_ada
        ipoin= lolne_ada(2,iolne)
        write (lun_oumgi_ada,*) iolne, conew_ada(1:ndime,ipoin) , '  1.0 '
     end do
     do ipoin=npopa_ada+1, nposu_ada
        write (lun_oumgi_ada,*) ipoin, conew_ada(1:ndime,ipoin) , '  2.0 '
     end do
     write(lun_oumgi_ada,*) nelsu_ada, ' 0 '
     do ielsu=1,nelpa_ada
        do inode=1,nnodb_ada
           lnofa(inode) = lolne_ada(1,lnoad_ada(inode,ielsu))
        end do
        write(lun_oumgi_ada,*) ielsu,lnofa(1:nnodb_ada)
     end do
     do ielsu=nelpa_ada+1,nelsu_ada
        write(lun_oumgi_ada,*) ielsu,lnosu_ada(1:nnodb_ada,ielsu)
     end do
     write(lun_oumgi_ada,*) '    0 '

  
  else
     !
     ! ONLY 2D IS PROGRAMMED!!!!!!!!!
     !

  end if

  
end subroutine ada_outmgi
