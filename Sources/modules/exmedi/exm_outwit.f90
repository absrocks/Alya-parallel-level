!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_outwit.f90
!> @author  Mariano Vazquez
!> @brief   Output witness points
!> @date    16/11/1966
!> @details Output witness points
!> @} 
!-----------------------------------------------------------------------
subroutine exm_outwit()

  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_exmedi
  implicit none
  integer(ip) :: iwitn,ielem,inode,pnode,pelty,ipoin,ifwit, ii

  ! Check: write witness or not?
  ifwit= 0
  call posdef(5_ip,ifwit)

  if (ifwit == 1) then
  
     if( INOTMASTER ) then 

        do iwitn = 1, nwitn
           ielem = lewit(iwitn)
           if( ielem > 0 ) then

               pelty = ltype(ielem)
               pnode = nnode(pelty)
               !
               ! Intra
               !
               if( postp(1) % npp_witne(1) == 1 ) then
                  do inode = 1,pnode
                     ipoin = lnods(inode,ielem)
                     witne(1,iwitn) = witne(1,iwitn) &
                           + shwit(inode,iwitn) *  elmag(ipoin,1) * postp(1) % witne_dt(1)
                  end do
               end if
               !
               ! Coord x
               !
               if( postp(1) % npp_witne(2) == 1 ) then
                  do inode = 1,pnode
                     ipoin = lnods(inode,ielem)
                     witne(2,iwitn) = witne(2,iwitn) &
                           + shwit(inode,iwitn) *  coord(1,ipoin) * postp(1) % witne_dt(2)
                  end do
               end if
               !
               ! VCONC
               !
               do ii = 1,size(vconc, 1_ip, KIND=ip) 
                  if( postp(1) % npp_witne(2 + ii) == 1 ) then
                     do inode = 1,pnode
                        ipoin = lnods(inode,ielem)
                        witne( 2+ii, iwitn ) = witne( 2+ii, iwitn ) &
                              + shwit(inode,iwitn) *  vconc(ii,ipoin,1) * postp(1) % witne_dt( 2+ii )
                     end do
                  end if
               end do

               !
               ! VICEL
               !
               do ii = 1,size(vicel_exm, 1_ip, KIND=ip) 
                  if( postp(1) % npp_witne(13 + ii) == 1 ) then
                     do inode = 1,pnode
                        ipoin = lnods(inode,ielem)
                        witne( 13+ii, iwitn ) = witne( 13+ii, iwitn ) &
                              + shwit(inode,iwitn) *  vicel_exm(ii,ipoin,1) * postp(1) % witne_dt( 13+ii )
                     end do
                  end if
               end do



           end if

        end do

        
     end if

  end if


end subroutine exm_outwit
