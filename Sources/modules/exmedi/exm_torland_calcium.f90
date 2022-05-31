!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_torland_calcium.f90
!> @author  Zhinuo Jenny Wang
!> @brief   Single cell run for Initial condition setup for ToR-ORd heterogeneous model
!> @date   16/NOV/1966
!> @details Runs a single cell simulation at th given frequency and pathologic conditions \n
!!    It performs single cell runs under normal, heart failure or drugs \n
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_torland_calcium(ipoin,sac,dCaTRPN)

  use      def_parame
  use      def_master
  use      def_elmtyp
  use      def_domain
  use      def_exmedi
  use      mod_exm_oharaequations
  use      mod_exm_sld_eccoupling, only: state_ecc, troponin_ecc

  implicit none

  real(rp), intent(out) :: sac,dCaTRPN
  real(rp) :: vaux1,vaux2,vaux3,rhsx1,rhsx2,rhsx,val0,k_TRPN,n_TRPN,Ca50,Ca50_ref,beta_1,lambda,  &
   C_tens(3,3),lambda0,gsac,esac

  integer(ip) , intent(in) :: ipoin
  integer(ip) :: idime,jdime,kdime

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

  ! Calcium calcium in the coupled model
  dCaTRPN = troponin_ecc(ipoin,1) - troponin_ecc(ipoin,2)

  ! Update troponin
  troponin_ecc(ipoin,2) = troponin_ecc(ipoin,1)



  ! Calculate SAC current
  !if (lambda0 >= 1.0_rp) then
  !  gsac = 0.1_rp
  !  esac = 0.0_rp
  !  sac = gsac*((lambda0-1.0_rp)/(1.1_rp-1.0_rp))*(volt-esac)
  !else
    sac = 0.0_rp  
  !end if

end subroutine exm_torland_calcium
