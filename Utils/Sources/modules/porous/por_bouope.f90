!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_bouope.f90
!> @date    29/07/2013
!> @author  Herbert Owen
!> @brief   Compute boundary contribution to matrix and RHS for the pressure
!> @details Compute boundary contribution to matrix and RHS for the pressure
!>          for the moment we do not need anything in the boundary - see lyx
!> @} 
!------------------------------------------------------------------------
subroutine por_bouope()
  use def_parame
  use def_master
  use def_kermod
  use mod_ker_proper 
  use def_domain
  use def_porous
  implicit none
  real(rp)    :: elrhs(mnode)
  real(rp)    :: baloc(ndime,ndime)
  real(rp)    :: gbpre(mgaub)
  real(rp)    :: gbswa(mgaub)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: bocod(ndime,mnodb)
  real(rp)    :: bopre(mnodb)
  real(rp)    :: boswa(mnodb)
  integer(ip) :: ielem,kboun,igaub,iboun,pblty
  integer(ip) :: pnodb,pmate,pnode,pelty,pgaus
  real(rp)    :: eucta,gbsur,auxii,auxi2,lam_o
  real(rp)    :: lam_w,dummr,deltp,xx_o,xx_w
  real(rp)    :: Bo,Bw,kro,krw
  !
  ! Loop over boundaries  
  !


end subroutine por_bouope
