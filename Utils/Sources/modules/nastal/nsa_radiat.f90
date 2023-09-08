!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_radiat.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Nastal coupling to a radiation heat source  
!> @details Nastal coupling to a radiation heat source  
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_radiat()
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  rhsid, radso, kfl_coupl,ID_NASTAL,ID_RADIAT
  use def_domain, only       :  npoin, vmass,mgaus,mnode
  use def_nastal

  implicit none 

  integer(ip)                :: pnode,inode,pgaus,igaus,ipoin,kpoin,ielem
  integer(ip)                :: pelty
  real(rp)                   :: gpsha(mnode,mgaus)
  real(rp)                   :: gprad(mgaus)
  !
  ! Coupling with RADIATion
  !
  if (kfl_coupl(ID_NASTAL,ID_RADIAT) >= 1 ) then
     do ipoin = 1,npoin
        rhsid(ipoin) = rhsid(ipoin)-radso(ipoin,1)*vmass(ipoin)
     end do
  end if

end subroutine nsa_radiat
