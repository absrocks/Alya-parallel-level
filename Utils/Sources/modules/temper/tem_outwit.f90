subroutine tem_outwit()
  !------------------------------------------------------------------------
  !****f* Temper/tem_outwit
  ! NAME 
  !    tem_outwit
  ! DESCRIPTION
  !    Output results on witness points
  ! USES
  ! USED BY
  !    tem_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_temper
  implicit none
  integer(ip) :: iwitn,ielem,inode,pnode,pelty,ipoin,ivawi,dummi

  if( nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0 ) then
     !
     ! Results on witness points
     !
     witne => postp(1) % witne

     if( INOTMASTER ) then 

        do iwitn = 1, nwitn
           ielem = lewit(iwitn)
           if( ielem > 0 ) then
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              do ivawi = 1,postp(1) % nvawi
                 witne(ivawi,iwitn) = 0.0_rp
              end do

              if( postp(1) % npp_witne(1) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(1,iwitn) = witne(1,iwitn) + shwit(inode,iwitn) * tempe(ipoin,1)
                 end do
              end if

              if( postp(1) % npp_witne(2) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(2,iwitn) = witne(2,iwitn) + dewit(1,inode,iwitn) * tempe(ipoin,1)
                 end do
              end if

              if( postp(1) % npp_witne(3) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(3,iwitn) = witne(3,iwitn) + dewit(2,inode,iwitn) * tempe(ipoin,1)
                 end do
              end if
 
              if( postp(1) % npp_witne(4) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(4,iwitn) = witne(4,iwitn) + dewit(3,inode,iwitn) * tempe(ipoin,1)
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

end subroutine tem_outwit

