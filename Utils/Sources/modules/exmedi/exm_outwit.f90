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
  integer(ip) :: iwitn,ielem,inode,pnode,pelty,ipoin,ivawi,dummi,ifwit,idime
  real(rp)    :: fb(3),st(6),stfib,fmod

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
              do ivawi = 1,postp(1) % nvawi
                 witne(ivawi,iwitn) = 0.0_rp
              end do
              !
              ! Intra
              !
              if( postp(1) % npp_witne(1) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(1,iwitn) = witne(1,iwitn) &
                         + shwit(inode,iwitn) *  elmag(ipoin,1) 
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


end subroutine exm_outwit
