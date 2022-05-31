!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    rhsmod.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Correct nodal vectors for periodicity and parallelization  
!> @details Correct nodal vectors for periodicity and parallelization according to:  
!>          NVBAR > 0  ... Array of size XARRA( NBVAR,NPOIN)
!>          NVBAR < 0  ... Array of size XARRA(-NBVAR,NBOPO)
!> @} 
!-----------------------------------------------------------------------
subroutine rhsmod(nbvar,xarra)
  use def_kintyp, only               :  ip,rp
  use def_master, only               :  INOTMASTER,NPOIN_REAL_2DIM,parr1,&
       &                                parr2,NBOPO_REAL_2DIM,NPOIN_TYPE,&
       &                                ISLAVE
  use def_domain, only               :  lpoty,lhang,nhang,npoin,nbopo
  use mod_couplings
  use def_master
  implicit none
  integer(ip), intent(in)            :: nbvar                  !< run through either npoin or nbopo
  real(rp),    intent(inout), target :: xarra(abs(nbvar),*)    !< in-out corrected vector
  integer(ip)                        :: ii,kk,nbva2,ihang,dummi
  real(rp)                           :: dummr

  if( INOTEMPTY ) then

     if( nbvar > 0 ) then
        !
        ! Parall service: exchange result between slaves
        !
        call pararr('SLX',NPOIN_TYPE,npoin*nbvar,xarra)

     else if( nbvar < 0 ) then
        !
        ! Parall service: exchange result between slaves
        !
        call runend('RHSMOD: DEPRECATED OPTION')
        
     end if

  end if

end subroutine rhsmod
