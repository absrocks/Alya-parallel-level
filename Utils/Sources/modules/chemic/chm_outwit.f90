subroutine chm_outwit()
  !------------------------------------------------------------------------
  !****f* chemic/chm_outwit
  ! NAME 
  !    chm_outwit
  ! DESCRIPTION
  !    Output chemic species on witness points
  ! USES
  ! USED BY
  !    chm_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_chemic
  implicit none
  integer(ip)       :: iwitn,ielem,inode,pnode,pelty,ipoin,iclas,ivawi,dummi

  if( nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0 ) then
     !
     ! Results on witness points
     !
     witne => postp(1) % witne
  
     if( INOTMASTER ) then

       do iwitn = 1,nwitn
          ielem = lewit(iwitn)
          if( ielem > 0 ) then
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             do ivawi = 1,postp(1) % nvawi
                witne(ivawi,iwitn) = 0.0_rp
             end do 
      
             if( postp(1) % npp_witne(1) == 1 ) then
                do iclas = 1,nclas_chm
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      witne(iclas,iwitn) = witne(iclas,iwitn) + shwit(inode,iwitn) * conce(ipoin,iclas,1)
                   end do
                end do
             end if             

          end if

       end do
     end if
     !
     ! Parall
     !
     call posdef(24_ip,dummi)

  end if

end subroutine chm_outwit

