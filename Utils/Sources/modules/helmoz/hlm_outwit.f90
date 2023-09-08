 subroutine hlm_outwit()
  !------------------------------------------------------------------------
  !****f* Helmoz/hlm_outwit
  ! NAME 
  !    hlm_outwit
  ! DESCRIPTION
  !    Output results on witness points
  ! USES
  ! USED BY
  !    hlm_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_helmoz
  implicit none
  integer(ip) :: iwitn,ielem,inode,pnode,pelty,ipoin,ivawi,dummi

  if( nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0 ) then

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
                    witne(1,iwitn) = witne(1,iwitn) + shwit(inode,iwitn) * real(elefi_hlm(1,ipoin))
                 end do
              end if

              if( postp(1) % npp_witne(2) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(2,iwitn) = witne(2,iwitn) + shwit(inode,iwitn) * real(elefi_hlm(2,ipoin))
                 end do
              end if

              if( postp(1) % npp_witne(3) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(3,iwitn) = witne(3,iwitn) + shwit(inode,iwitn) * real(elefi_hlm(ndime,ipoin))
                 end do
              end if

              if( postp(1) % npp_witne(4) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(4,iwitn) = witne(4,iwitn) + shwit(inode,iwitn) * aimag(elefi_hlm(1,ipoin))
                 end do
              end if

              if( postp(1) % npp_witne(5) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(5,iwitn) = witne(5,iwitn) + shwit(inode,iwitn) * aimag(elefi_hlm(2,ipoin))
                 end do
              end if

              if( postp(1) % npp_witne(6) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(6,iwitn) = witne(6,iwitn) + shwit(inode,iwitn) * aimag(elefi_hlm(ndime,ipoin))
                 end do
              end if

             if( postp(1) % npp_witne(7) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(7,iwitn) = witne(7,iwitn) + shwit(inode,iwitn) * real(magfi_hlm(1,ipoin))
                 end do
              end if

              if( postp(1) % npp_witne(8) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(8,iwitn) = witne(8,iwitn) + shwit(inode,iwitn) * real(magfi_hlm(2,ipoin))
                 end do
              end if

              if( postp(1) % npp_witne(9) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(9,iwitn) = witne(9,iwitn) + shwit(inode,iwitn) * real(magfi_hlm(ndime,ipoin))
                 end do
              end if 

             if( postp(1) % npp_witne(10) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(10,iwitn) = witne(10,iwitn) + shwit(inode,iwitn) * aimag(magfi_hlm(1,ipoin))
                 end do
              end if

              if( postp(1) % npp_witne(11) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(11,iwitn) = witne(11,iwitn) + shwit(inode,iwitn) * aimag(magfi_hlm(2,ipoin))
                 end do
              end if

              if( postp(1) % npp_witne(12) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(12,iwitn) = witne(12,iwitn) + shwit(inode,iwitn) * aimag(magfi_hlm(ndime,ipoin))
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

end subroutine hlm_outwit

