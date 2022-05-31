!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_inihet.f90
!> @author  Jazmin Aguado-Sierra
!> @brief   Initial condition setup for TenTuscher-Panfilov 2006 heterogeneous model
!> @details Obtains initial conditions of Normal or Heart Failure cell at  70 bpm or 857 ms \n
!!    Otherwise, it obtains the initial conditions from single cell runs in ONECEL.F90 AND OCEIHE.F90 \n
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_inihet(ipoin,mat)
!subroutine exm_inihet(kmodel_ipoin,ipoin,mat)

  use      def_master
  use      def_elmtyp
  use      def_domain

  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: ipoin, mat
  integer(ip)             :: ituss_exm
  real(rp)    ::  vaux1, vaux2, vaux3, rhsx
  !integer(ip) ::  mat

  if(kfl_timei_exm==1) then 
     !
     ! Load initial conditions for the potentials
     !
     !mat=nodemat(ipoin)
     ituss_exm = int(celty_exm(1,ipoin),KIND=ip)
     !write(*,*) ituss_exm
     if(ituss_exm == EXM_CELLTYPE_EPI) then   !epicardial
        elmag(ipoin,1:3) = vminimate_exm(3,mat)          ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        vauxi_exm(1,ipoin,1:3) = vauin_exm(1,3,mat)          !Variable m
        vauxi_exm(2,ipoin,1:3) = vauin_exm(2,3,mat)          !Variable h 
        vauxi_exm(3,ipoin,1:3) = vauin_exm(3,3,mat)          !Variable j
        vauxi_exm(4,ipoin,1:3) = vauin_exm(4,3,mat)          !Variable d
        vauxi_exm(5,ipoin,1:3) = vauin_exm(5,3,mat)          !Variable f
        vauxi_exm(12,ipoin,1:3) = vauin_exm(12,3,mat)        !Variable f2
        vauxi_exm(6,ipoin,1:3) = vauin_exm(6,3,mat)          ! Variable fCa
        vauxi_exm(7,ipoin,1:3) = vauin_exm(7,3,mat)          !Variable r
        vauxi_exm(8,ipoin,1:3) = vauin_exm(8,3,mat)          !Variable s
        vauxi_exm(9,ipoin,1:3) = vauin_exm(9,3,mat)          !Variable xs 
        vauxi_exm(10,ipoin,1:3) = vauin_exm(10,3,mat)        !Variable xr1
        vauxi_exm(11,ipoin,1:3) = vauin_exm(11,3,mat)        !Variable xr2 
        vconc(1,ipoin,1:3) = vcoin_exm(1,3,mat)              !Variable Cai
        vconc(2,ipoin,1:3) = vcoin_exm(2,3,mat)              !Variable CaSR
        vconc(3,ipoin,1:3) = vcoin_exm(3,3,mat)              !Variable Nai
        vconc(4,ipoin,1:3) = vcoin_exm(4,3,mat)              !Variable Ki
        vconc(5,ipoin,1:3) = vcoin_exm(5,3,mat)              !Variable CaSS
        vconc(9,ipoin,1:3) = vcoin_exm(9,3,mat)              !Variable Rprime
     else if(ituss_exm == EXM_CELLTYPE_ENDO) then  !endocardial
        elmag(ipoin,1:3) = vminimate_exm(1,mat)          ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        vauxi_exm(1,ipoin,1:3) = vauin_exm(1,1,mat)          !Variable m
        vauxi_exm(2,ipoin,1:3) = vauin_exm(2,1,mat)          !Variables h 
        vauxi_exm(3,ipoin,1:3) = vauin_exm(3,1,mat)          !Variable j
        vauxi_exm(4,ipoin,1:3) = vauin_exm(4,1,mat)          !Variable d
        vauxi_exm(5,ipoin,1:3) = vauin_exm(5,1,mat)          !Variable f
        vauxi_exm(12,ipoin,1:3) = vauin_exm(12,1,mat)        !Variable f2
        vauxi_exm(6,ipoin,1:3) = vauin_exm(6,1,mat)          !Variable fCa
        vauxi_exm(7,ipoin,1:3) = vauin_exm(7,1,mat)          !Variable r
        vauxi_exm(8,ipoin,1:3) = vauin_exm(8,1,mat)          !Variable s
        vauxi_exm(9,ipoin,1:3) = vauin_exm(9,1,mat)          !Variable xs 
        vauxi_exm(10,ipoin,1:3) = vauin_exm(10,1,mat)        !Variable xr1
        vauxi_exm(11,ipoin,1:3) = vauin_exm(11,1,mat)        !Variable xr2 
        vconc(1,ipoin,1:3) = vcoin_exm(1,1,mat)              !Variable Cai
        vconc(2,ipoin,1:3) = vcoin_exm(2,1,mat)              !Variable CaSR
        vconc(3,ipoin,1:3) = vcoin_exm(3,1,mat)              !Variable Nai
        vconc(4,ipoin,1:3) = vcoin_exm(4,1,mat)              !Variable Ki
        vconc(5,ipoin,1:3) = vcoin_exm(5,1,mat)              !Variable CaSS
        vconc(9,ipoin,1:3) = vcoin_exm(9,1,mat)              !Variable Rprime
     else if(ituss_exm == EXM_CELLTYPE_MID) then  !mid
        elmag(ipoin,1:3) = vminimate_exm(2,mat)          ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        vauxi_exm(1,ipoin,1:3) = vauin_exm(1,2,mat)          !Variable m
        vauxi_exm(2,ipoin,1:3) = vauin_exm(2,2,mat)          !Variables h 
        vauxi_exm(3,ipoin,1:3) = vauin_exm(3,2,mat)          !Variable j
        vauxi_exm(4,ipoin,1:3) = vauin_exm(4,2,mat)          !Variable d
        vauxi_exm(5,ipoin,1:3) = vauin_exm(5,2,mat)          !Variable f
        vauxi_exm(12,ipoin,1:3) = vauin_exm(12,2,mat)        !Variable f2
        vauxi_exm(6,ipoin,1:3) = vauin_exm(6,2,mat)          ! Variable fCa
        vauxi_exm(7,ipoin,1:3) = vauin_exm(7,2,mat)          !Variable r
        vauxi_exm(8,ipoin,1:3) = vauin_exm(8,2,mat)          !Variable s
        vauxi_exm(9,ipoin,1:3) = vauin_exm(9,2,mat)          !Variable xs 
        vauxi_exm(10,ipoin,1:3) = vauin_exm(10,2,mat)        !Variable xr1
        vauxi_exm(11,ipoin,1:3) = vauin_exm(11,2,mat)        !Variable xr2 
        vconc(1,ipoin,1:3) = vcoin_exm(1,2,mat)              !Variable Cai
        vconc(2,ipoin,1:3) = vcoin_exm(2,2,mat)              !Variable CaSR
        vconc(3,ipoin,1:3) = vcoin_exm(3,2,mat)              !Variable Nai
        vconc(4,ipoin,1:3) = vcoin_exm(4,2,mat)              !Variable Ki
        vconc(5,ipoin,1:3) = vcoin_exm(5,2,mat)              !Variable CaSS
        vconc(9,ipoin,1:3) = vcoin_exm(9,2,mat)              !Variable Rprime
     end if
     ! Variable of activation O (related with I_rel current)
     ! nconc = 10
     ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
     !  kcasr = max_sr - (max_sr - min_sr)/(1+(pow((EC/Ca_SR), 2)))
     vaux1 = 1.0_rp + ((1.5_rp / vconc(2,ipoin,2))*(1.5_rp / vconc(2,ipoin,2)))  !conc2 = CaSRfree = CaSR
     vaux1 = 1.0_rp / vaux1
     vaux1 = 2.5_rp - (1.5_rp * vaux1)  !kCaSR
     !     O = ( k1*pow(Ca_ss, 2)*R_prime)/(k3+ k1*pow(Ca_ss, 2))      
     vaux3 = 0.15_rp / vaux1   !K1
     vaux2 = 0.045_rp * vaux1  !K2
     rhsx = 1.0_rp / (0.06_rp + (vaux3*vconc(5,ipoin,2)*vconc(5,ipoin,2)))      !O for calcium dynamics
     rhsx = vaux3 * vconc(5,ipoin,2) * vconc(5,ipoin,2) * vconc(9,ipoin,2) * rhsx !O
     vconc(10,ipoin,1:3) = rhsx 
  end if

  !write(991,*) vconc(1,:,1)  !JAS
end subroutine exm_inihet
