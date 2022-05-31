!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_ohaland_calcium.f90
!> @author  Jazmin Aguado-Sierra and Francesc Levrero-Florencio
!> @brief   Single cell run for Initial condition setup for Ohara-Rudy 2011 heterogeneous model
!> @date   16/NOV/1966
!> @details Runs a single cell simulation at th given frequency and pathologic conditions \n
!!    It performs single cell runs under normal, heart failure or drugs \n
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_ohaland_calcium(kmcmdn,cmdnmax,trpnmax,vnsr,vmyo,vss,acap,farad,dtimeEP,ipoin,cai,sac,volt)
!subroutine exm_ohaland_calcium(kmcmdn,kmtrpn,cmdnmax,trpnmax,vnsr,vmyo,vss,acap,farad,dtimeEP,ipoin,cai,sac,volt)

  use      def_parame
  use      def_master
  use      def_elmtyp
  use      def_domain
  use      def_exmedi
  use      mod_exm_oharaequations
  use      mod_exm_sld_eccoupling, only: troponin_ecc

  implicit none

  real(rp), intent(in) :: kmcmdn,cmdnmax,trpnmax,vnsr,vmyo,vss,acap,farad,dtimeEP,volt
  !real(rp), intent(in) :: kmtrpn
  real(rp), intent(out) :: cai, sac
  real(rp) :: vaux1,vaux3,rhsx1,rhsx2,rhsx,val0,k_TRPN,n_TRPN,Ca50,Ca50_ref,beta_1,lambda,  &
   C_tens(3,3),dCaTRPN,lambda0,gsac,esac

  integer(ip) , intent(in) :: ipoin
  integer(ip) :: idime,jdime,kdime

  ! Declare some parameters of the Land model
  k_TRPN = 0.1_rp
  n_TRPN = 2.0_rp
  Ca50_ref = 0.805_rp
  beta_1 = -2.4_rp

  ! Recover the right Cauchy-Green strain tensor
  C_tens = 0.0_rp
  do idime = 1,ndime
    do jdime = 1,ndime
      do kdime = 1,ndime
        C_tens(idime,jdime) = C_tens(idime,jdime)+gdepo(kdime,idime,ipoin)*                   &
         gdepo(kdime,jdime,ipoin)
      end do
    end do
  end do

  ! Recover the stretch in the fibre direction
  lambda = 0.0_rp
  do idime = 1,ndime
    do jdime = 1,ndime
      lambda = lambda+fiber(idime,ipoin)*C_tens(idime,jdime)*fiber(jdime,ipoin)
    end do
  end do
  lambda = SQRT(lambda)

  ! Calculate calcium sensitivity
  lambda0 = lambda
  lambda = MIN(1.2_rp,lambda)
  Ca50 = Ca50_ref+beta_1*(lambda-1.0_rp)

  ! Calcium calcium in the coupled model
  vaux1 = (kmcmdn+vconc(1,ipoin,2)) ** 2.0_rp
  !vaux2 = (kmtrpn+vconc(1,ipoin,2)) ** 2.0_rp
  !vaux3 = 1.0_rp/(1.0_rp + (cmdnmax*kmcmdn/vaux1) + (trpnmax*kmtrpn/vaux2))
  vaux3 = 1.0_rp/(1.0_rp + (cmdnmax*kmcmdn/vaux1))
  rhsx1 = vicel_exm(14,ipoin,1) + vicel_exm(13,ipoin,1) - (2.0_rp*vicel_exm(8,ipoin,1))
  !rhsx2 = -(vicel_exm(18,ipoin,1)*vnsr/vmyo) + (vicel_exm(15,ipoin,1)*vss/vmyo)
  dCaTRPN = troponin_ecc(ipoin,1) - troponin_ecc(ipoin,2)
  rhsx2 = -(vicel_exm(18,ipoin,1)*vnsr/vmyo) + (vicel_exm(15,ipoin,1)*vss/vmyo) - trpnmax*dCaTRPN
  rhsx = vaux3 *(-(rhsx1 *acap/(2.0_rp*farad*vmyo)) + rhsx2)
  val0 = vconc(1,ipoin,2)
  cai = val0 + dtimeEP*rhsx

  ! Update troponin
  troponin_ecc(ipoin,2) = troponin_ecc(ipoin,1)

  ! Calculate SAC current
  if (lambda0 >= 1.0_rp) then
    gsac = 0.1_rp
    esac = 0.0_rp
    sac = gsac*((lambda0-1.0_rp)/(1.1_rp-1.0_rp))*(volt-esac)
  else
    sac = 0.0_rp  
  end if

end subroutine exm_ohaland_calcium
