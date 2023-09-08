subroutine nsa_lawvis(ktask,icomp,xvisc,xtemp,dvite)
!-----------------------------------------------------------------------
!****f* Nastal/nsa_lawvis
! NAME 
!    nsa_lawvis
! DESCRIPTION
!    This routine computes the viscosity from a given viscosity law 
!
!    Calculation of visco and dvis/dtem in a given point.
!
!    Non dimensionalized law (nsa_lawvis=1):
!            
!        mu = mu_ref T^1.5 (1 + c)/(T + c)
!
!        c = 110K/Tphys_inf = 198.6R/Tphys_inf
!
!    Here, sutha = mu_ref (1 + c)  
!           suthb = c
!           T = T/Tref
!
!    Non dimensionalized restricted law (nsa_lawvis=2):
!            
!        mu = mu_ref T^0.768
!
!    Here, sutha = mu_ref   
!           suthb = c
!           T = T/Tref
!      
! USED BY
!    nsa_updtss
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame
  use      def_nastal

  implicit none
  integer(ip),   intent(in) :: ktask,icomp
  real(rp),      intent(out):: xvisc,dvite
  real(rp),    intent(inout):: xtemp
  integer(ip)               :: ipoin,kpoin

  xvisc= 0.0_rp
  dvite= 0.0_rp

  if (ktask == zero) then

     ! Compute the whole viscosity field

     if (lawvi_nsa==0) then
        do ipoin = 1,npoin
           visco(ipoin,icomp)=  visco_nsa
        end do
     else if (lawvi_nsa==1) then
        do ipoin = 1,npoin
           xtemp= tempe(ipoin,icomp)
           visco(ipoin,icomp)=  (vispa_nsa(1) * xtemp ** vispa_nsa(2)) / (xtemp+vispa_nsa(3))
        end do
     else if (lawvi_nsa==2) then
        do ipoin = 1,npoin
           xtemp= tempe(ipoin,icomp) / tempe_nsa
           visco(ipoin,icomp)=  visco_nsa * xtemp ** 0.76_rp
        end do        
     end if
         

  else if (ktask > zero) then

     ! Compute the viscosity field value for ktask node
     
     ipoin = ktask

     if (lawvi_nsa==0) then
        visco(ipoin,icomp)=  visco_nsa
     else if (lawvi_nsa==1) then               ! power law
        xtemp= tempe(ipoin,icomp)
        visco(ipoin,icomp)=  vispa_nsa(1) * xtemp ** vispa_nsa(2)
        dvite=  vispa_nsa(1)*vispa_nsa(2)*xtemp**(vispa_nsa(2)-1.0_rp)        
     else if (lawvi_nsa==2) then               ! sutherland law
        xtemp= tempe(ipoin,icomp)
        visco(ipoin,icomp)=  (vispa_nsa(1) * xtemp ** 1.5_rp) / (xtemp+vispa_nsa(2))
        dvite=  visco(ipoin,icomp)*( 1.5_rp/xtemp - 1.0_rp/(xtemp+vispa_nsa(2)))
     end if

  else if (ktask == -1) then

     ! Compute the viscosity xvisc using the arguments
     
     if (lawvi_nsa==0) then
        xvisc=  visco_nsa
     else if (lawvi_nsa==1) then
        xvisc=  vispa_nsa(1) * xtemp ** vispa_nsa(2)
        dvite=  vispa_nsa(1)*vispa_nsa(2)*xtemp**(vispa_nsa(2)-1.0_rp)        
     else if (lawvi_nsa==2) then
        xvisc= (vispa_nsa(1) * xtemp ** 1.5_rp) / (xtemp + vispa_nsa(2))
        dvite=  xvisc*( 1.5_rp/xtemp - 1.0_rp/(xtemp+vispa_nsa(2)))
     end if

  end if

end subroutine nsa_lawvis
