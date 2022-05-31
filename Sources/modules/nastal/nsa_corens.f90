!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_corens.f90
!> @author  Mariano Vazquez
!> @date    19/03/2014
!> @brief   Compute navier-stokes nodal residuals 
!> @details Compute navier-stokes nodal residuals 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_corens(inew,iold) 
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  
  implicit none
  integer(ip)    :: inew   !> current step label
  integer(ip)    :: iold   !> old step label
  integer(ip)    :: kpoin,ipoin,idime
  real(rp)       :: va,vd,va_plus,vd_plus,debug_resi(3)
  
  if (INOTMASTER) then
     do ipoin = 1,npoin
        va_plus = 0.0_rp
        vd_plus = 0.0_rp
        do idime=1,ndime
           vd = umome(idime,ipoin,inew) - umome(idime,ipoin,iold)
           va = umome(idime,ipoin,inew)
           va_plus = va_plus + va * va
           vd_plus = vd_plus + vd * vd
        end do
        crens_nsa(1,ipoin) = sqrt(vd_plus)
        if (sqrt(va_plus) > zeror) then
           crens_nsa(1,ipoin) = crens_nsa(1,ipoin) / sqrt(va_plus)
        end if

        vd = densi(ipoin,inew) - densi(ipoin,iold)
        va = densi(ipoin,inew)
        va_plus = va * va
        vd_plus = vd * vd
        crens_nsa(2,ipoin) = sqrt(vd_plus)
        if (sqrt(va_plus) > zeror) then
           crens_nsa(2,ipoin) = crens_nsa(2,ipoin) / sqrt(va_plus)
        end if

        vd = energ(ipoin,inew) - energ(ipoin,iold)
        va = energ(ipoin,1)
        va_plus = va * va
        vd_plus = vd * vd
        crens_nsa(3,ipoin) = sqrt(vd_plus)
        if (sqrt(va_plus) > zeror) then
           crens_nsa(3,ipoin) = crens_nsa(3,ipoin) / sqrt(va_plus)
        end if

     end do
  end if

end subroutine nsa_corens
