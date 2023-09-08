!########################### CONSTANTS #####################
module mod_exm_hipscventriconst
!! Constants
  use def_kintyp, only : ip,rp

  implicit none

!ionic concentrations
real(rp), parameter :: &
Na_o = 151.0,  K_o = 5.4,  Ca_o = 1.8,  K_i = 150.0

!!Ventricular
!Cell size
real(rp), parameter :: &
C_m = 98.7109e-12,    V_c = 8800.0,    V_sr = 583.73, &
!Maximum conductances and currents
g_Na = 3.6712302e3,    g_to = 29.9038,    g_K1 = 28.1492,    P_NaK = 1.841424,    K_NaCa = 4900.0,    Vmax_up = 0.56064


!Maximum conductances and currents
real(rp), parameter :: &
g_CaL = 8.635702e-5,  g_Kr = 29.8667,  g_Ks = 2.041,  g_f = 30.10312,  a_rel = 16.464, &
b_rel = 0.25,  c_rel = 8.232,  V_leak = 4.4444e-4,  g_pCa = 0.4125,  g_bNa = 0.9,  g_bCa = 0.69264

!Other constants
real(rp), parameter :: &
Buf_c = 0.25,  Buf_sr = 10.0,  K_Buf_c = 0.001,  K_Buf_sr = 0.3,  K_up = 0.00025, &
K_pCa = 0.0005,  F = 96485.3415,  Rcons = 8.314472,  Temp = 310.0,  L_0 = 0.025,  P_kna = 0.03, &
K_sat = 0.1,  Km_Ca = 1.38,  Km_Na_i = 87.5,  Alpha = 2.8571432,  Gamma = 0.35, &
K_mNa = 40.0,  K_mk = 1.0
end module mod_exm_hipscventriconst

!#####################################################################################
