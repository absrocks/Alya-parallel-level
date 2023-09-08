!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_elmope.f90
!> @author  Guillaume Houzeaux
!> @brief   Coupling with Partis module
!> @details Remove momentum from the Navier-Stokes equations
!>          \verbatim
!>          x_p ......... particle position
!>          Ff .......... force fluid on paticle
!>          Navier-Stokes = - \int_W (Ff/Vp).v.delta(x-x_p) dw,
!>                        = - (Ff(xp)/Vp).v(xp)*Vp
!>          So that the nodeal force Fi is:
!>                     Fi = - Fp * Ni
!>          \endverbatim
!> @}
!------------------------------------------------------------------------

subroutine nsi_partis()
  use def_master
  use def_domain
  use def_elmtyp
  use def_nastin
  implicit none
  integer(ip) :: idime,idofn,ipoin

  do ipoin = 1,npoin
     do idime = 1,ndime
        idofn = (ipoin-1)*ndime+idime
        if( kfl_fixno_nsi(idime,ipoin) <= 0 ) then
           rhsid(idofn) = rhsid(idofn) + momen(idime,ipoin)
        end if
     end do
     momen(1:ndime,ipoin) = 0.0_rp
  end do

end subroutine nsi_partis
