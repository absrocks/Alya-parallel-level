!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_setvar.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute derived fields
!> @details Compute derived fields
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_evaluate_uinlet(itask,icomp)
  use      def_master
  use      def_domain
  use      def_parame
  use      def_kermod
  use      mod_communications, only : PAR_SUM
  use      def_nastal

  implicit none
  integer(ip)       :: itask,icomp,ipoin,kpoin,icoun,kfixi
  real(rp)          :: vmodu,totvolu


  if (itask == 0) then

     nuinlet_nsa = 0

     if (INOTMASTER) then
        ! only count on "my nodes", done in two loops
        do ipoin=1,npoi1
           if (kfl_fixno_nsa(ndime+1,ipoin) == 3) then
              nuinlet_nsa = nuinlet_nsa + 1
           end if
        end do
        do ipoin=npoi2,npoi3
           if (kfl_fixno_nsa(ndime+1,ipoin) == 3) then
              nuinlet_nsa = nuinlet_nsa + 1
           end if
        end do
     end if

     call PAR_SUM(nuinlet_nsa)

  else if (itask == 1) then

     uinlet_nsa = 0.0_rp
     totvolu = 0.0_rp

     if (INOTMASTER) then
        
        ! only count on "my nodes", done in two loops
        do ipoin=1,npoi1        
           if (kfl_fixno_nsa(ndime+1,ipoin)  == 3) then
              vmodu = veloc(1,ipoin,icomp) * veloc(1,ipoin,icomp) + veloc(2,ipoin,icomp) * veloc(2,ipoin,icomp) 
              if (ndime == 3 ) vmodu = vmodu + veloc(3,ipoin,icomp) * veloc(3,ipoin,icomp)
              vmodu = sqrt(vmodu)
              uinlet_nsa = uinlet_nsa + vmass(ipoin)*vmodu
              totvolu = totvolu + vmass(ipoin)
           end if
        end do
        do ipoin=npoi2,npoi3        
           if (kfl_fixno_nsa(ndime+1,ipoin)  == 3) then
              vmodu = veloc(1,ipoin,icomp) * veloc(1,ipoin,icomp) + veloc(2,ipoin,icomp) * veloc(2,ipoin,icomp) 
              if (ndime == 3 ) vmodu = vmodu + veloc(3,ipoin,icomp) * veloc(3,ipoin,icomp)
              vmodu = sqrt(vmodu)
              uinlet_nsa = uinlet_nsa + vmass(ipoin)*vmodu
              totvolu = totvolu + vmass(ipoin)
           end if
        end do
     end if

     call PAR_SUM(uinlet_nsa)
     call PAR_SUM(totvolu)

     uinlet_nsa = uinlet_nsa / totvolu
     uinlet_nsa = sqrt(uinlet_nsa)

  end if

end subroutine nsa_evaluate_uinlet
