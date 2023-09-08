module mod_exm_oharaequations

  use def_kintyp, only : ip,rp

  implicit none

  private :: eq_inf_coefs,eq_tau_coefs
  
  public

  integer(ip), parameter :: &
     OHR_M_INF             =1,  OHR_H_INF       =2,  OHR_H_CAMK_INF=3,  OHR_M_L_INF    =4,  &
     OHR_H_L_INF           =5,  OHR_H_L_CAMK_INF=6,  OHR_A_INF     =7,  OHR_INV_TAU_A_1=8,  &
     OHR_INV_TAU_A_2       =9,  OHR_I_INF       =10, OHR_A_I_FAST  =11, OHR_A_CAMK_INF =12, &
     OHR_DELTA_CAMK_RECOVER=13, OHR_D_INF       =14, OHR_F_INF     =15, OHR_A_F_CA_FAST=16, &
     OHR_X_R_INF           =17, OHR_A_XR_FAST   =18, OHR_R_KR_1    =19, OHR_R_KR_2     =20, &
     OHR_X_S1_INF          =21, OHR_X_KB        =22, OHR_DELTA_EPI =23
  integer(ip), parameter :: EXM_OHR_INF_EQUATIONS = 23

  integer(ip), parameter :: &
     OHR_TAU_M      =1,  OHR_TAU_H_FAST =2,  OHR_TAU_H_SLOW        =3,  OHR_TAU_J        =4, &
     OHR_TAU_I_FAST =5,  OHR_TAU_I_SLOW =6,  OHR_DELTA_CAMK_DEVELOP=7,  OHR_TAU_D        =8, &
     OHR_TAU_F_FAST =9,  OHR_TAU_F_SLOW =10, OHR_TAU_F_CA_FAST     =11, OHR_TAU_F_CA_SLOW=12, &
     OHR_TAU_XR_FAST=13, OHR_TAU_XR_SLOW=14, OHR_TAU_X_S1          =15, OHR_TAU_X_S2     =16, &
     OHR_TAU_X_K1   =17
  integer(ip), parameter :: EXM_OHR_TAU_EQUATIONS = 17
  

  real(rp)  :: eq_inf_coefs(EXM_OHR_INF_EQUATIONS,4) = reshape( (/ &
     0.0_rp,  1.0_rp,             -1.0_rp / 9.871_rp,    39.57_rp,   & ! M_INF
     0.0_rp,  1.0_rp,              1.0_rp / 6.086_rp,     82.9_rp,   & ! H_INF !Passini et al. 6.22_rp,78.5_rp, (original: 6.086_rp,    82.9_rp,)
     0.0_rp,  1.0_rp,              1.0_rp / 6.086_rp,     89.1_rp,   & ! H_CAMK_INF !Passini et al. 6.22_rp,84.7_rp, (original: 6.086_rp,    89.1_rp,) 
     0.0_rp,  1.0_rp,             -1.0_rp / 5.264_rp,    42.85_rp,   & ! M_L_INF
     0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    87.61_rp,   & ! H_L_INF
     0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    93.81_rp,   & ! H_L_CAMK_INF
     0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -14.34_rp,   & ! A_INF
     0.0_rp,  1.0_rp / 1.2089_rp, -1.0_rp / 29.3814_rp, -18.4099_rp, & ! INV_TAU_A_1
     0.0_rp,  3.5_rp,              1.0_rp / 29.3814_rp,  100.0_rp,   & ! INV_TAU_A_2
     0.0_rp,  1.0_rp,              1.0_rp / 5.711_rp,    43.94_rp,   & ! I_INF
     0.0_rp,  1.0_rp,              1.0_rp / 151.2_rp,   -213.6_rp,   & ! A_I_FAST
     0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -24.34_rp,   & ! A_CAMK_INF
     1.0_rp, -0.5_rp,              1.0_rp / 20.0_rp,     70.0_rp,    & ! DELTA_CAMK_RECOVER
     0.0_rp,  1.0_rp,             -1.0_rp / 4.230_rp,    3.940_rp,   & ! D_INF
     0.0_rp,  1.0_rp,              1.0_rp / 3.696_rp,    19.58_rp,   & ! F_INF
     0.3_rp,  0.6_rp,              1.0_rp / 10.0_rp,    -10.0_rp,    & ! A_F_CA_FAST
     0.0_rp,  1.0_rp,             -1.0_rp / 6.789_rp,    8.337_rp,   & ! X_R_INF
     0.0_rp,  1.0_rp,              1.0_rp / 38.21_rp,    54.81_rp,   & ! A_XR_FAST
     0.0_rp,  1.0_rp,              1.0_rp / 75.0_rp,     55.0_rp,    & ! R_KR_1
     0.0_rp,  1.0_rp,              1.0_rp / 30.0_rp,    -10.0_rp,    & ! R_KR_2
     0.0_rp,  1.0_rp,             -1.0_rp / 8.932_rp,    11.60_rp,   & ! X_S1_INF
     0.0_rp,  1.0_rp,             -1.0_rp / 18.34_rp,   -14.48_rp,   & ! X_KB
     1.0_rp, -0.95_rp,             1.0_rp / 5.0_rp,      70.0_rp     & ! DELTA_EPI
     /), (/ EXM_OHR_INF_EQUATIONS,4_ip /), ORDER= (/ 2_ip,1_ip /) )
 

  real(rp)  :: eq_tau_coefs(EXM_OHR_TAU_EQUATIONS,8) = reshape( (/ &
     0.0_rp,    1.0_rp,    6.765_rp,      1.0_rp / 34.77_rp,  11.64_rp, 8.552_rp,     -1.0_rp / 5.955_rp,   77.42_rp,  & ! TAU_M
     0.0_rp,    1.0_rp,    1.432E-05_rp, -1.0_rp / 6.285_rp,  1.196_rp, 6.149_rp,      1.0_rp / 20.27_rp,   0.5096_rp,      & ! TAU_H_FAST !Ohara Original
     0.0_rp,    1.0_rp,    0.009794_rp,  -1.0_rp / 28.05_rp,  17.95_rp, 0.3343_rp,     1.0_rp / 56.66_rp,   5.730_rp,  & ! TAU_H_SLOW
     2.038_rp,  1.0_rp,    0.02136_rp,   -1.0_rp / 8.281_rp,  100.6_rp, 0.3052_rp,     1.0_rp / 38.45_rp,   0.9941_rp,     & ! TAU_J  !Ohara original
     4.562_rp,  1.0_rp,    0.3933_rp,    -1.0_rp / 100.0_rp,  100.0_rp, 0.08004_rp,    1.0_rp / 16.59_rp,   50.0_rp,   & ! TAU_I_FAST
     23.62_rp,  1.0_rp,    0.001416_rp,  -1.0_rp / 59.05_rp,  96.52_rp, 1.7808E-8_rp,  1.0_rp / 8.079_rp,   114.1_rp,  & ! TAU_I_SLOW
     1.354_rp,  1.0E-3_rp, 1.0_rp,        1.0_rp / 15.89_rp, -167.4_rp, 1.0_rp,       -1.0_rp / 0.2154_rp, -12.23_rp,  & ! DELTA_CAMK_DEVELOP
     0.6_rp,    1.0_rp,    1.0_rp,       -0.05_rp,            6.0_rp,   1.0_rp,        0.09_rp,             14.0_rp,   & ! TAU_D
     7.0_rp,    1.0_rp,    0.0045_rp,    -1.0_rp / 10.0_rp,   20.0_rp,  0.0045_rp,     1.0_rp / 10.0_rp,    20.0_rp,   & ! TAU_F_FAST
     1000.0_rp, 1.0_rp,    0.000035_rp,  -1.0_rp / 4.0_rp,    5.0_rp,   0.000035_rp,   1.0_rp / 6.0_rp,     5.0_rp,    & ! TAU_F_SLOW
     7.0_rp,    1.0_rp,    0.04_rp,      -1.0_rp / 7.0_rp,   -4.0_rp,   0.04_rp,       1.0_rp / 7.0_rp,    -4.0_rp,    & ! TAU_F_CA_FAST
     100.0_rp,  1.0_rp,    0.00012_rp,   -1.0_rp / 3.0_rp,    0.0_rp,   0.00012_rp,    1.0_rp / 7.0_rp,     0.0_rp,    & ! TAU_F_CA_SLOW
     12.98_rp,  1.0_rp,    0.3652_rp,     1.0_rp / 3.869_rp, -31.66_rp, 4.123E-5_rp,  -1.0_rp / 20.38_rp,  -47.78_rp,  & ! TAU_XR_FAST
     1.865_rp,  1.0_rp,    0.06629_rp,    1.0_rp / 7.355_rp, -34.70_rp, 1.128E-5_rp,  -1.0_rp / 25.94_rp,  -29.74_rp,  & ! TAU_XR_SLOW
     817.3_rp,  1.0_rp,    2.326E-4_rp,   1.0_rp / 17.80_rp,  48.28_rp, 0.001292_rp,  -1.0_rp / 230.0_rp,   210.0_rp,  & ! TAU_X_S1
     0.0_rp,    1.0_rp,    0.01_rp,       1.0_rp / 20.0_rp,  -50.0_rp,  0.0193_rp,    -1.0_rp / 31.0_rp,    66.54_rp,  & ! TAU_X_S2
     0.0_rp,    122.2_rp,  1.0_rp,       -1.0_rp / 20.36_rp,  127.2_rp, 1.0_rp,        1.0_rp / 69.33_rp,   236.8_rp   & ! TAU_X_K1
     /), (/ EXM_OHR_TAU_EQUATIONS, 8_ip /), ORDER = (/ 2_ip, 1_ip /) )

contains

  subroutine exm_ohara_infs(v, res,ipoin,kfl_paral)

     ! Computes for each i in [1..neqs]:
     ! res(i) = a(i) + b(i) / ( 1 + exp(c(i) ( v + d(i) ) ) )
     real(rp), intent(in)    :: v
     real(rp), intent(out)   :: res(EXM_OHR_INF_EQUATIONS)
     integer(ip) :: kfl_paral,ipoin

     integer(ip) :: ieq

     do ieq = 1,EXM_OHR_INF_EQUATIONS
        res(ieq) = eq_inf_coefs(ieq,1) &
                 + eq_inf_coefs(ieq,2) / ( 1.0_rp + exp( eq_inf_coefs(ieq,3) * ( v + eq_inf_coefs(ieq,4) ) ) )
     end do

  end subroutine exm_ohara_infs


  subroutine exm_ohara_taus(v, res)

     ! Computes for each i in [1..neqs]:
     ! res(i) = a(i) + b(i) / ( c(i) exp( d(i) (v + e(i)) ) + f(i) exp( g(i) (v + h(i)) ))
     real(rp), intent(in)    :: v
     real(rp), intent(out)   :: res(EXM_OHR_TAU_EQUATIONS)

     integer(ip) :: ieq

     do ieq = 1,EXM_OHR_TAU_EQUATIONS
        res(ieq) = eq_tau_coefs(ieq,1) &
                 + eq_tau_coefs(ieq,2) / ( eq_tau_coefs(ieq,3) * exp( eq_tau_coefs(ieq,4) * (v + eq_tau_coefs(ieq,5)) ) +&
                                           eq_tau_coefs(ieq,6) * exp( eq_tau_coefs(ieq,7) * (v + eq_tau_coefs(ieq,8)) ) )
     end do
  end subroutine exm_ohara_taus
  
  subroutine exm_passini_ohara_table()
  
  implicit none
  real(rp)    :: inf1, inf2, inf3, inf4, tau1, tau2, tau3, tau4, tau5, tau6
  real(rp)    :: tau7, tau8, tau9, tau10, tau11, tau12, tau13  

  inf1=6.22_rp   !Hinf gate coef 3
  inf2=78.5_rp   !Hinf gate coef 4
  inf3=6.22_rp   !Hsp gate coef 3
  inf4=84.7_rp   !Hsp gate coef 4
  tau1=3.686E-06_rp !Hf tau coef 3
  tau2=7.8579_rp    !Hf tau coef 4
  tau3=3.8875_rp    !Hf tau coef 5
  tau4=16.0_rp      !Hf tau coef 6
  tau5=9.1843_rp    !Hf tau coef 7
  tau6=0.4963_rp    !Hf tau coef 8           
  tau7=4.859_rp     !J tau coef 1
  tau8=0.8628_rp    !J tau coef 3
  tau9=7.6005_rp    !J tau coef 4
  tau10=116.7258_rp !J tau coef 5
  tau11=1.1096_rp   !J tau coef 6
  tau12=9.0358_rp   !J tau coef 7
  tau13=6.2719_rp   !J tau coef 8    
  
  eq_inf_coefs = reshape( (/ &
     0.0_rp,  1.0_rp,             -1.0_rp / 9.871_rp,    39.57_rp,   & ! M_INF
     0.0_rp,  1.0_rp,              1.0_rp / inf1,        inf2,       & ! H_INF !Passini et al. 6.22_rp,78.5_rp, (original: 6.086_rp,    82.9_rp,)
     0.0_rp,  1.0_rp,              1.0_rp / inf3,        inf4,       & ! H_CAMK_INF !Passini et al. 6.22_rp,84.7_rp, (original: 6.086_rp,    89.1_rp,) 
     0.0_rp,  1.0_rp,             -1.0_rp / 5.264_rp,    42.85_rp,   & ! M_L_INF
     0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    87.61_rp,   & ! H_L_INF
     0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    93.81_rp,   & ! H_L_CAMK_INF
     0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -14.34_rp,   & ! A_INF
     0.0_rp,  1.0_rp / 1.2089_rp, -1.0_rp / 29.3814_rp, -18.4099_rp, & ! INV_TAU_A_1
     0.0_rp,  3.5_rp,              1.0_rp / 29.3814_rp,  100.0_rp,   & ! INV_TAU_A_2
     0.0_rp,  1.0_rp,              1.0_rp / 5.711_rp,    43.94_rp,   & ! I_INF
     0.0_rp,  1.0_rp,              1.0_rp / 151.2_rp,   -213.6_rp,   & ! A_I_FAST
     0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -24.34_rp,   & ! A_CAMK_INF
     1.0_rp, -0.5_rp,              1.0_rp / 20.0_rp,     70.0_rp,    & ! DELTA_CAMK_RECOVER
     0.0_rp,  1.0_rp,             -1.0_rp / 4.230_rp,    3.940_rp,   & ! D_INF
     0.0_rp,  1.0_rp,              1.0_rp / 3.696_rp,    19.58_rp,   & ! F_INF
     0.3_rp,  0.6_rp,              1.0_rp / 10.0_rp,    -10.0_rp,    & ! A_F_CA_FAST
     0.0_rp,  1.0_rp,             -1.0_rp / 6.789_rp,    8.337_rp,   & ! X_R_INF
     0.0_rp,  1.0_rp,              1.0_rp / 38.21_rp,    54.81_rp,   & ! A_XR_FAST
     0.0_rp,  1.0_rp,              1.0_rp / 75.0_rp,     55.0_rp,    & ! R_KR_1
     0.0_rp,  1.0_rp,              1.0_rp / 30.0_rp,    -10.0_rp,    & ! R_KR_2
     0.0_rp,  1.0_rp,             -1.0_rp / 8.932_rp,    11.60_rp,   & ! X_S1_INF
     0.0_rp,  1.0_rp,             -1.0_rp / 18.34_rp,   -14.48_rp,   & ! X_KB
     1.0_rp, -0.95_rp,             1.0_rp / 5.0_rp,      70.0_rp     & ! DELTA_EPI
     /), (/ EXM_OHR_INF_EQUATIONS,4_ip /), ORDER= (/ 2_ip,1_ip /) )
 

  eq_tau_coefs = reshape( (/ &
     0.0_rp,    1.0_rp,    6.765_rp,      1.0_rp / 34.77_rp,  11.64_rp, 8.552_rp,     -1.0_rp / 5.955_rp,   77.42_rp,  & ! TAU_M
     0.0_rp,    1.0_rp,    tau1,         -1.0_rp / tau2,      tau3,     tau4,          1.0_rp / tau5,       tau6,      & ! TAU_H_FAST !Ohara Original
     0.0_rp,    1.0_rp,    0.009794_rp,  -1.0_rp / 28.05_rp,  17.95_rp, 0.3343_rp,     1.0_rp / 56.66_rp,   5.730_rp,  & ! TAU_H_SLOW
     tau7,      1.0_rp,    tau8,         -1.0_rp / tau9,      tau10,    tau11,         1.0_rp / tau12,      tau13,     & ! TAU_J  !Ohara original
     4.562_rp,  1.0_rp,    0.3933_rp,    -1.0_rp / 100.0_rp,  100.0_rp, 0.08004_rp,    1.0_rp / 16.59_rp,   50.0_rp,   & ! TAU_I_FAST
     23.62_rp,  1.0_rp,    0.001416_rp,  -1.0_rp / 59.05_rp,  96.52_rp, 1.7808E-8_rp,  1.0_rp / 8.079_rp,   114.1_rp,  & ! TAU_I_SLOW
     1.354_rp,  1.0E-3_rp, 1.0_rp,        1.0_rp / 15.89_rp, -167.4_rp, 1.0_rp,       -1.0_rp / 0.2154_rp, -12.23_rp,  & ! DELTA_CAMK_DEVELOP
     0.6_rp,    1.0_rp,    1.0_rp,       -0.05_rp,            6.0_rp,   1.0_rp,        0.09_rp,             14.0_rp,   & ! TAU_D
     7.0_rp,    1.0_rp,    0.0045_rp,    -1.0_rp / 10.0_rp,   20.0_rp,  0.0045_rp,     1.0_rp / 10.0_rp,    20.0_rp,   & ! TAU_F_FAST
     1000.0_rp, 1.0_rp,    0.000035_rp,  -1.0_rp / 4.0_rp,    5.0_rp,   0.000035_rp,   1.0_rp / 6.0_rp,     5.0_rp,    & ! TAU_F_SLOW
     7.0_rp,    1.0_rp,    0.04_rp,      -1.0_rp / 7.0_rp,   -4.0_rp,   0.04_rp,       1.0_rp / 7.0_rp,    -4.0_rp,    & ! TAU_F_CA_FAST
     100.0_rp,  1.0_rp,    0.00012_rp,   -1.0_rp / 3.0_rp,    0.0_rp,   0.00012_rp,    1.0_rp / 7.0_rp,     0.0_rp,    & ! TAU_F_CA_SLOW
     12.98_rp,  1.0_rp,    0.3652_rp,     1.0_rp / 3.869_rp, -31.66_rp, 4.123E-5_rp,  -1.0_rp / 20.38_rp,  -47.78_rp,  & ! TAU_XR_FAST
     1.865_rp,  1.0_rp,    0.06629_rp,    1.0_rp / 7.355_rp, -34.70_rp, 1.128E-5_rp,  -1.0_rp / 25.94_rp,  -29.74_rp,  & ! TAU_XR_SLOW
     817.3_rp,  1.0_rp,    2.326E-4_rp,   1.0_rp / 17.80_rp,  48.28_rp, 0.001292_rp,  -1.0_rp / 230.0_rp,   210.0_rp,  & ! TAU_X_S1
     0.0_rp,    1.0_rp,    0.01_rp,       1.0_rp / 20.0_rp,  -50.0_rp,  0.0193_rp,    -1.0_rp / 31.0_rp,    66.54_rp,  & ! TAU_X_S2
     0.0_rp,    122.2_rp,  1.0_rp,       -1.0_rp / 20.36_rp,  127.2_rp, 1.0_rp,        1.0_rp / 69.33_rp,   236.8_rp   & ! TAU_X_K1
     /), (/ EXM_OHR_TAU_EQUATIONS, 8_ip /), ORDER = (/ 2_ip, 1_ip /) )
     
  end subroutine exm_passini_ohara_table

end module mod_exm_oharaequations

