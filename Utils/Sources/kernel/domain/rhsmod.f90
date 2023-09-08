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
  implicit none
  integer(ip), intent(in)            :: nbvar                  !< run through either npoin or nbopo
  real(rp),    intent(inout), target :: xarra(abs(nbvar),*)    !< in-out corrected vector
  integer(ip)                        :: ii,kk,nbva2,ihang,dummi
  real(rp)                           :: dummr

  if( INOTMASTER ) then

     if( nbvar > 0 ) then
        !
        ! Periodicity
        ! 
        call presol(0_ip,nbvar,nbvar,xarra,dummr,dummr)
        !
        ! hanging nodes
        !
        if( nhang > 0 ) then
           !do ihang = 1,nhang
           !   ii = lhang(1,ihang)
           !   do kk = 1,nbvar
           !      xarra(kk,ii) = 0.0_rp
           !   end do
           !end do
        end if
        !
        ! Parall service: exchange result between slaves
        !
        if( ISLAVE .and. npoin > 0 ) then
           call pararr('SLX',NPOIN_TYPE,npoin*nbvar,xarra)
        end if
        !
        ! Periodicity: recover solution on slave nodes
        ! 
        call presol(2_ip,nbvar,nbvar,dummr,dummr,xarra)

     else if( nbvar < 0 ) then
        nbva2 = -nbvar
        !
        ! Parall service: exchange result between slaves
        !
        if( ISLAVE ) then
           call vocabu(NBOPO_REAL_2DIM,nbva2,0_ip)
           parr2 => xarra(1:nbva2,1:nbopo)
           call Parall(400_ip)
           nullify(parr2)
        end if
     end if

  end if

end subroutine rhsmod
