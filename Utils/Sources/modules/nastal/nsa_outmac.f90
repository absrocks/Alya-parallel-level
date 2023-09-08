subroutine nsa_outmac(ktask,jpoin,icomp,xmach)
!-----------------------------------------------------------------------
!****f* Nastal/nsa_outmac
! NAME 
!    nsa_setmac
! DESCRIPTION
!    This routine computes mach number, storing it in unkno(1:npoin)
! USED BY
!    
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame
  use      def_kermod
  use      mod_ker_proper

  use      def_nastal

  implicit none
  integer(ip) :: ktask,icomp,ipoin,idime,jpoin,dummi
  real(rp)    :: vmodu,sound,xmach,dummy(ndime,ndime),adgam

  xmach = 0.0_rp

  if(kfl_paral/=0) then

     if (ktask == zero) then                                            

        if (kfl_prope /= 0 ) then
!           wmean = mowei_nsa
           !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa,dummi,dummi,dummy(:,1),dummy(:,1))
           call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa)
        else
           wmean = mowei_nsa
           shecp_nsa = cpcoe_nsa
        endif

        do ipoin=1,npoin
   
           ! Thermodynamic properties
           !
           adgam = shecp_nsa(ipoin) / (shecp_nsa(ipoin) - runiv_nsa/ wmean(ipoin,1))

           ! Compute the field mach values for node ktask and component icomp

           sound = sqrt(adgam * press(ipoin,icomp) / densi(ipoin,icomp))
!           sound= sofac_nsa * sqrt(tempe(ipoin,icomp))        
!           sound= sqrt((adgam_nsa - 1.0_rp)*cpcoe_nsa*tempe(ipoin,icomp))        
           if (kfl_relat_nsa == 1) then
              sound= sqrt((1.0_rp - 1.0_rp / entha(ipoin,icomp))*(adgam-1.0_rp))              
           end if
           vmodu = veloc(1,ipoin,icomp)*veloc(1,ipoin,icomp) + veloc(2,ipoin,icomp)*veloc(2,ipoin,icomp) 
           if (ndime == 3) vmodu = vmodu + veloc(3,ipoin,icomp)*veloc(3,ipoin,icomp) 
           vmodu = sqrt(vmodu)

           xmach = vmodu/sound

           if (ktask == 0) unkno(ipoin) = xmach

        end do

     end if

  end if

end subroutine nsa_outmac
