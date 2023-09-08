!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_difwel.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Add diffusion at the wells 
!> @details Add diffusion at the wells - for the moment a constant value given at por.dat file
!> @} 
!------------------------------------------------------------------------
subroutine por_difwel(&
     pnode,pgaus,lnods,gpdif_aux)
!
  use def_kintyp, only       :  ip,rp
  use def_porous, only       :  difwe_por,iwell_por
  implicit none 
  integer(ip), intent(in)    :: pnode                        !> number of nodes
  integer(ip), intent(in)    :: pgaus                        !> number of gauss points
  integer(ip), intent(in)    :: lnods(pnode)                 !> conectivity
  real(rp),    intent(out)   :: gpdif_aux(pgaus)             !> Added diffusion
 

  integer(ip)                :: inode,igaus,iauxi

  iauxi = 0_ip
  do inode = 1,pnode
     if ( iwell_por( lnods(inode) ) /= 0 ) iauxi = 1_ip
  end do  

  if ( iauxi == 1_ip) then
     do igaus = 1,pgaus
        gpdif_aux(igaus) = gpdif_aux(igaus) + difwe_por
     end do
  end if

end subroutine por_difwel
