subroutine exm_ionicunew(ipoin,xlian,xpion,xpfhn,xpfen,xpten)

!-----------------------------------------------------------------------
!
! This routine computes the ionic current contribution
!
!
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_solver

  use      def_exmedi

  implicit none
  integer(ip) :: ipoin,kmodel
  real(rp)    :: xmeps, xchic ,xmalp,xmgam,taufi,xpfhn,xpfen,xpten,xpion,xlian

  xpfhn = 1.0_rp
  xpfen = 0.0_rp
  xpten = 0.0_rp



  kmodel= kfl_cellmod(nodemat(ipoin))


  if (kmodel == 0) then                       !no model
     xmeps = 0.0_rp
     xmalp = 0.0_rp
     xmgam = 0.0_rp
     taufi = 0.0_rp
     xpion = 0.0_rp

  else if (kmodel == 11) then          !fitzhugh nagumo
     
     xmeps = xmopa_exm( 1,nodemat(ipoin))           ! epsilon for the rcp
     xchic = xmopa_exm( 2,nodemat(ipoin))
     xmalp = xmopa_exm( 3,nodemat(ipoin))           ! alpha for the Iion
     xmgam = xmopa_exm( 4,nodemat(ipoin))           ! gamma for the rcp
     taufi = xmopa_exm( 5,nodemat(ipoin))           ! C for the Iion

     xpion =        - taufi * xlian * (xlian - 1.0_rp) * (xlian - xmalp)

     xpfhn = 1.0_rp + taufi * (xlian - 1.0_rp) * (xlian - xmalp) * dtime / xmccm_exm 

  else
     call runend('EXM_IONICUnew: IONIC CURRENT FOR SELECTED MODEL NOT PROGRAMMED')

  end if
  
end subroutine exm_ionicunew
