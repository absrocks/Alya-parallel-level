!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_residu.f90
!> @author  Mariano Vazquez
!> @date    04/12/2013
!> @brief   Compute convergence residuals and convergence fields
!> @details Compute convergence residuals and convergence fields 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_residu(rinsa) 
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  
  implicit none
  integer(ip)    :: kpoin,ipoin,idime
  real(rp)       :: rinsa(4)           !> convergence residuals
  real(rp)       :: va,vd,va_plus,vd_plus,resid(3,2),debug_resi(3)
  

  ! OJO QUE FALTA PARAL!!!

  if (INOTMASTER) then
     rinsa= 0.0_rp
     resid= 0.0_rp
     debug_resi= 0.0_rp
     do ipoin = 1,npoin
        va_plus = 0.0_rp
        vd_plus = 0.0_rp
        do idime=1,ndime
           vd = umome(idime,ipoin,1) - umome(idime,ipoin,3)
           va = umome(idime,ipoin,1)
           resid(1,1)= resid(1,1)+vd*vd
           resid(1,2)= resid(1,2)+va*va
           va_plus = va_plus + va * va
           vd_plus = vd_plus + vd * vd
        end do
        crens_nsa(1,ipoin) = sqrt(vd_plus)
!        if (sqrt(va_plus) > zeror) then
!           crens_nsa(1,ipoin) = crens_nsa(1,ipoin) / sqrt(va_plus)
!        end if
!        crens_nsa(1,ipoin) = crens_nsa(1,ipoin) / resou_nsa(1)

        vd = densi(ipoin,1) - densi(ipoin,3)
        va = densi(ipoin,1)
        resid(2,1)= resid(2,1)+vd*vd
        resid(2,2)= resid(2,2)+va*va
        va_plus = va * va
        vd_plus = vd * vd
        crens_nsa(2,ipoin) = sqrt(vd_plus)
!        if (sqrt(va_plus) > zeror) then
!           crens_nsa(2,ipoin) = crens_nsa(2,ipoin) / sqrt(va_plus)
!        end if
!        crens_nsa(2,ipoin) = crens_nsa(2,ipoin) / resou_nsa(2)

        vd = energ(ipoin,1) - energ(ipoin,3)
        va = energ(ipoin,1)
        resid(3,1)= resid(3,1)+vd*vd
        resid(3,2)= resid(3,2)+va*va
        va_plus = va * va
        vd_plus = vd * vd
        crens_nsa(3,ipoin) = sqrt(vd_plus)
!        if (sqrt(va_plus) > zeror) then
!           crens_nsa(3,ipoin) = crens_nsa(3,ipoin) / sqrt(va_plus)
!        end if
!        crens_nsa(3,ipoin) = crens_nsa(3,ipoin) / resou_nsa(3)

     end do

     if(resid(1,2)>zeror) then
        rinsa(1) = sqrt(resid(1,1)/resid(1,2))
     else
        rinsa(1) = sqrt(resid(1,1))
     end if
     
     if(resid(2,2)>zeror) then
        rinsa(2) = sqrt(resid(2,1)/resid(2,2))
     else
        rinsa(2) = sqrt(resid(2,1))
     end if
     
     if(resid(3,2)>zeror) then
        rinsa(3) = sqrt(resid(3,1)/resid(3,2))
     else
        rinsa(3) = sqrt(resid(3,1))
     end if

     rinsa(4) = rinsa(1) + rinsa(2) + rinsa(3)

     rinsa(1)= rinsa(1)/resou_nsa(1)
     rinsa(2)= rinsa(2)/resou_nsa(2)
     rinsa(3)= rinsa(3)/resou_nsa(3)
     rinsa(4)= rinsa(4)/resou_nsa(4)


  end if

end subroutine nsa_residu
