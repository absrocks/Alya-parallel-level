!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_elmgat.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Gather operations
!> @details Gather operations
!> @} 
!------------------------------------------------------------------------
subroutine por_elmgat(&
     pnode,lnods,elpre,elswa,elcod,elvel,kfl_gavel)
 
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode,coord
  use def_master, only     :  press,wasat,veloc
  use def_porous, only     :  kfl_timei_por,kfl_tisch_por,kprsa_por
  implicit none
  integer(ip), intent(in)  :: pnode,kfl_gavel
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: elpre(pnode,*)
  real(rp),    intent(out) :: elswa(pnode,*)
  real(rp),    intent(out) :: elcod(ndime,mnode)
  real(rp),    intent(out) :: elvel(ndime,mnode)
  integer(ip)              :: inode,ipoin,idime
  !
  ! Current pressure, saturation and coordinates
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     elpre(inode,1) = press(ipoin,1)
     elswa(inode,1) = wasat(ipoin,1)
     do idime=1,ndime
        elcod(idime,inode) = coord(idime,ipoin)
     end do  
  end do
  if ( kfl_gavel == 1_ip ) then
     do inode=1,pnode
        ipoin=lnods(inode)
        do idime=1,ndime
           elvel(idime,inode) = veloc(idime,ipoin,1)
        end do
     end do
  end if
  !
  ! Time integration
  !
  if( kfl_timei_por(kprsa_por) == 1 ) then
     do inode=1,pnode
        ipoin=lnods(inode)
        elpre(inode,2) = press(ipoin,3)
        elswa(inode,2) = wasat(ipoin,3)
     end do
     if(kfl_tisch_por==2) then
        call runend('POR_ELMGAT:kfl_tisch_por==2 not ready')
     end if
  end if 

end subroutine por_elmgat
