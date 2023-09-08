subroutine rad_outwit()
  !------------------------------------------------------------------------
  !****f* Radiat/rad_outwit
  ! NAME 
  !    rad_outwit
  ! DESCRIPTION
  !    Output results on witness points
  ! USES
  ! USED BY
  !    rad_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_radiat
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
                    witne(1,iwitn) = witne(1,iwitn) + shwit(inode,iwitn) * radav_rad(ipoin,1)
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

end subroutine rad_outwit

