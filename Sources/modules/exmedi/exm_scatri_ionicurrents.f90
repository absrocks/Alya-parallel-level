
!!#############################
    !voltage: elmag
    
    !Ca_i:vconc(1,:)
    !Ca_SR:vconc(2,:)
    !Na_i:vconc(3,:)
    
    !h:vauxi_exm(1,:)  ;  j:vauxi_exm(2,:)  ;  !m:vauxi_exm(3,:)
    !d:vauxi_exm(4,:)  ;  f_Ca:vauxi_exm(5,:)  ;  !f1:vauxi_exm(6,:)
    !f2:vauxi_exm(7,:)  ;  r:vauxi_exm(8,:)  ;  !q:vauxi_exm(9,:)
    !xr1:vauxi_exm(10,:)  ;  xr2:vauxi_exm(11,:)  ;  !Xs:vauxi_exm(12,:)
    !Xf:vauxi_exm(13,:)  ;  g:vauxi_exm(14,:)
    
    !I_Na:vicel_exm(1,:)  ;  I_Cal:vicel_exm(2,:)  ;  !I_to:vicel_exm(3,:)
    !I_kr:vicel_exm(4,:)  ;  I_ks:vicel_exm(5,:)  ;  !I_k1:vicel_exm(6,:)
    !I_f:vicel_exm(7,:)  ;  I_NaK:vicel_exm(8,:)  ;  !I_NaCa:vicel_exm(9,:)
    !I_rel:vicel_exm(10,:)  ;  I_up:vicel_exm(11,:)  ;  !I_leak:vicel_exm(12,:)
    !I_pCa:vicel_exm(13,:)  ;  I_bNa:vicel_exm(14,:)  ;  I_bCa:vicel_exm(15,:)
!!#############################
    
!! main code
subroutine exm_scatri_ionicurrents(ipoin,xioni,dioni)

  use      def_parame
  use      def_master
  use      def_elmtyp
  use      def_domain
  use      def_exmedi
  !use      mod_exm_hipscatriaconst
  implicit none

  integer(ip), intent(in) :: ipoin !< node
  real(rp), intent(out) :: xioni   !< current
  real(rp), intent(out)   :: dioni !< current derivative
  
  real(rp) :: tt, h, j, m, d, f_Ca, f1, f2, r, q, xr1, xr2, xs, xf, g, Ca_i, Ca_SR, Na_i, V
  
  real(rp) ::  func_CaiA, func_CaSRA, func_dA, func_f1A, &
  func_f2A, func_fCaA, func_XsA, func_gA, func_hA, func_jA, func_mA, &
  func_NaiA, func_qA, func_rA, func_VA, func_XfA, func_xr1A, func_xr2A
  
  real(rp) :: E_Na, E_K, E_Ks, E_f, E_Ca, Alpha_K1, beta_K1, x_K1_inf

  !! !! defining functions and Solving the system of ODE using Runge-kutta of 4th order

  real(rp) :: df1 = -1.0

  !! RK constants declaration
  real(rp) :: K1h, K1j, K1m, K1d, K1f_Ca, K1f1, K1f2, &
  K1r, K1q, K1xr1, K1xr2, K1Xs, K1Xf, K1g, K1Ca_i, K1Ca_SR, K1Na_i, K1V

  real(rp) :: K2h, K2j, K2m, K2d, K2f_Ca, K2f1, K2f2, &
  K2r, K2q, K2xr1, K2xr2, K2Xs, K2Xf, K2g, K2Ca_i, K2Ca_SR, K2Na_i, K2V

  real(rp) :: K3h, K3j, K3m, K3d, K3f_Ca, K3f1, K3f2, &
  K3r, K3q, K3xr1, K3xr2, K3Xs, K3Xf, K3g, K3Ca_i, K3Ca_SR, K3Na_i, K3V

  real(rp) :: K4h, K4j, K4m, K4d, K4f_Ca, K4f1, K4f2, &
  K4r, K4q, K4xr1, K4xr2, K4Xs, K4Xf, K4g, K4Ca_i, K4Ca_SR, K4Na_i, K4V

  real(rp) :: Na_o, K_o, Ca_o, K_i
  real(rp) :: C_m, V_c, V_sr,  g_Na, g_to, g_K1, P_NaK, K_NaCa, Vmax_up
  real(rp) :: g_CaL, g_Kr, g_Ks, g_f, a_rel, b_rel, c_rel, V_leak, g_pCa, g_bNa, g_bCa
  real(rp) :: Buf_c, Buf_sr, K_Buf_c, K_Buf_sr, K_up, K_pCa, F, Rcons, Temp, L_0, P_kna
  real(rp) :: K_sat, Km_Ca, Km_Na_i, Alpha, Gamma, K_mNa, K_mk
!Maximum conductances and current



!ionic concentrations
Na_o = 151.0_rp  
K_o = 5.4_rp  
Ca_o = 1.8_rp  
K_i = 150.0_rp

!!Atrial
!Cell size
C_m = 78.6671e-12_rp    
V_c = 7012.0_rp    
V_sr = 465.20_rp
!Maximum conductances and currents
g_Na = 6.646185_rp*(10.0_rp**3.0_rp)   
g_to = 59.8077_rp    
g_K1 = 19.1925_rp    
P_NaK = 1.4731392_rp    
K_NaCa = 2450.0_rp    
Vmax_up = 0.22_rp
!Maximum conductances and currents
g_CaL = 8.635702e-5_rp  
g_Kr = 29.8667_rp  
g_Ks = 2.041_rp  
g_f = 30.10312_rp  
a_rel = 16.464_rp
b_rel = 0.25_rp  
c_rel = 8.232_rp  
V_leak = 4.4444e-4_rp  
g_pCa = 0.4125_rp  
g_bNa = 0.9_rp  
g_bCa = 0.69264_rp
!Other constants
Buf_c = 0.25_rp  
Buf_sr = 10.0_rp  
K_Buf_c = 0.001_rp  
K_Buf_sr = 0.3_rp  
K_up = 0.00025_rp
K_pCa = 0.0005_rp
F = 96485.3415_rp  
Rcons = 8.314472_rp  
Temp = 310.0_rp  
L_0 = 0.025_rp  
P_kna = 0.03_rp
K_sat = 0.1_rp  
Km_Ca = 1.38_rp  
Km_Na_i = 87.5_rp  
Alpha = 2.8571432_rp  
Gamma = 0.35_rp
K_mNa = 40.0_rp  
K_mk = 1.0_rp

  !! State variables declaration
  V = elmag(ipoin,ITER_K)/1000.0_rp
  tt = cutim
  h = vauxi_exm(1,ipoin,2)
  j = vauxi_exm(2,ipoin,2)
  m = vauxi_exm(3,ipoin,2)
  d = vauxi_exm(4,ipoin,2)
  f_Ca = vauxi_exm(5,ipoin,2)
  f1 = vauxi_exm(6,ipoin,2)
  f2 = vauxi_exm(7,ipoin,2)
  r = vauxi_exm(8,ipoin,2)
  Q = vauxi_exm(9,ipoin,2)
  xr1 = vauxi_exm(10,ipoin,2)
  xr2 = vauxi_exm(11,ipoin,2)
  xs = vauxi_exm(12,ipoin,2)
  xf = vauxi_exm(13,ipoin,2)
  g = vauxi_exm(14,ipoin,2)
  Ca_i = vconc(1,ipoin,2)
  Ca_SR = vconc(2,ipoin,2)
  Na_i = vconc(3,ipoin,2)
  
      !!Currents
    E_Na = (Rcons*Temp/F)*log(Na_o/Na_i)
    E_K = (Rcons*Temp/F)*log(K_o/K_i)
    E_Ks = (Rcons*Temp/F)*log((K_o+P_kna*Na_o)/(K_i+P_kna*Na_i))
    E_f = -0.017
    E_Ca = (0.5*Rcons*Temp/F)*log(Ca_o/Ca_i)
    Alpha_K1 = 3.91/(1+exp(0.5942*(V*1000.0-E_K*1000.0-200.0)))
    beta_K1 = (-1.509*exp(0.0002*(V*1000.0-E_K*1000.0+100.0))+exp(0.5886*(V*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(V*1000.0-E_K*1000.0)))
    x_K1_inf = Alpha_K1/(Alpha_K1+beta_K1)
    
    vicel_exm(1,ipoin,1) = g_Na*(m**3.0)*h*j*(V-E_Na)
    vicel_exm(2,ipoin,1) = (g_CaL*4.0*V*F**2.0)*(Ca_i*exp(2.0*V*F/(Rcons*Temp))-0.341*Ca_o)*d*f1*f2*f_Ca/(Rcons*Temp*(exp(2.0*V*F/(Rcons*Temp))-1.0))  
    vicel_exm(3,ipoin,1) = g_to*r*q*(V-E_K)
    vicel_exm(4,ipoin,1) = g_Kr*sqrt(K_o/5.4)*xr1*xr2*(V-E_K)
    vicel_exm(5,ipoin,1) = g_Ks*Xs**2.0*(1.0+0.6/(1.0+(3.8*1e-5/Ca_i)**1.4))*(V-E_Ks)
    vicel_exm(6,ipoin,1) = g_K1*x_K1_inf*sqrt(K_o/5.4)*(V-E_K)
    vicel_exm(7,ipoin,1) = g_f*Xf*(V-E_f)
    vicel_exm(8,ipoin,1) = ((P_NaK*K_o*Na_i/(K_o+K_mk))/(Na_i+K_mNa))/(1.0+0.1245*exp(-0.1*V*F/(Rcons*Temp))+0.353*exp(-V*F/(Rcons*Temp)))
    vicel_exm(9,ipoin,1) = K_NaCa*(exp(Gamma*V*F/(Rcons*Temp))*Na_i**3.0*Ca_o-exp((Gamma-1.0)*V*F/(Rcons*Temp))*Na_o**3.0*Ca_i*Alpha)/((Km_Na_i**3.0+Na_o**3.0)*(Km_Ca+Ca_o)*(1.0+K_sat*exp((Gamma-1.0)*V*F/(Rcons*Temp))))
    vicel_exm(10,ipoin,1) = ((a_rel*(Ca_SR**2.0)/(b_rel**2.0+Ca_SR**2.0))+c_rel)*d*g*(0.0556)
    vicel_exm(11,ipoin,1) = Vmax_up/(1.0+K_up**2.0/Ca_i**2.0)
    vicel_exm(12,ipoin,1) = V_leak*(Ca_SR-Ca_i)
    vicel_exm(13,ipoin,1) = g_pCa*Ca_i/(Ca_i+K_pCa)
    vicel_exm(14,ipoin,1) = g_bNa*(V-E_Na)
    vicel_exm(15,ipoin,1) = g_bCa*(V-E_Ca)
  
  !! Applying Runge-Kutta
  !open(unit=1, file='voltage.dat')
    
    !! K1
    K1h = dtime*func_hA(tt,V,h)
    K1j = dtime*func_jA(tt,V,j)
    K1m = dtime*func_mA(tt,V,m)
    K1d = dtime*func_dA(tt,V,d)
    K1f_Ca = dtime*func_fCaA(tt,f_Ca,Ca_i,V)
    K1f1 = dtime*func_f1A(tt,V,f1,Ca_i,df1)
    K1f2 = dtime*func_f2A(tt,V,f2)
    K1r = dtime*func_rA(tt,V,r)
    K1q = dtime*func_qA(tt,V,q)
    K1xr1 = dtime*func_xr1A(tt,V,xr1)
    K1xr2 = dtime*func_xr2A(tt,V,xr2)
    K1Xs = dtime*func_XsA(tt,V,xs)
    K1Xf = dtime*func_XfA(tt,V,xf)
    K1g = dtime*func_gA(tt,Ca_i,g,V)
    K1Ca_i = dtime*func_CaiA(tt,Ca_i,Ca_SR,d,g,V,Na_i,f1,f2,f_Ca)
    K1Ca_SR = dtime*func_CaSRA(tt,Ca_i,Ca_SR,d,g)
    K1Na_i = dtime*func_NaiA(tt,m,h,j,V,Na_i,Ca_i)
    K1V = dtime*func_VA(tt,V,r,q,xr1,xr2,xs,Ca_i,d,f1,f2,f_Ca,Na_i,m,h,j,xf)

    !! K2
    K2h =  dtime*func_hA(tt+dtime/2.0,  V+K1V/2.0,  h+K1h/2.0)
    K2j =  dtime*func_jA(tt+dtime/2.0,  V+K1V/2.0,  j+K1j/2.0)
    K2m =  dtime*func_mA(tt+dtime/2.0,  V+K1V/2.0,  m+K1m/2.0)
    K2d =  dtime*func_dA(tt+dtime/2.0,  V+K1V/2.0,  d+K1d/2.0)
    K2f_Ca =  dtime*func_fCaA(tt+dtime/2.0,  f_Ca+K1f_Ca/2.0,  Ca_i+K1Ca_i/2.0 ,V+K1V/2.0)
    K2f1 = dtime*func_f1A(tt+dtime/2.0,  V+K1V/2.0,  f1+K1f1/2.0  ,Ca_i+K1f1/2.0 ,df1)
    K2f2 = dtime*func_f2A(tt+dtime/2.0,  V+K1V/2.0,  f2+K1f2/2.0)
    K2r = dtime*func_rA(tt+dtime/2.0,  V+K1V/2.0,  r+K1r/2.0)
    K2q = dtime*func_qA(tt+dtime/2.0,  V+K1V/2.0,  q+K1q/2.0)
    K2xr1 = dtime*func_xr1A(tt+dtime/2.0,V+K1V/2.0,  xr1+K1xr1/2.0)
    K2xr2 = dtime*func_xr2A(tt+dtime/2.0,  V+K1V/2.0,  xr2+K1xr2/2.0)
    K2Xs = dtime*func_XsA(tt+dtime/2.0,  V+K1V/2.0,  xs+K1Xs/2.0)
    K2Xf = dtime*func_XfA(tt+dtime/2.0,  V+K1V/2.0,  xf+K1Xf/2.0)
    K2g = dtime*func_gA(tt+dtime/2.0,  Ca_i+K1Ca_i/2.0,  g+K1g/2.0 ,V+K1V/2.0)
    K2Ca_i = dtime*func_CaiA(tt+dtime/2.0,  Ca_i+K1Ca_i/2.0,  Ca_SR+K1Ca_SR/2.0,  d+K1d/2.0,  g+K1g/2.0,  V+K1V/2.0,  Na_i+K1Na_i/2.0,  f1+K1f1/2.0  ,f2+K1f2/2.0  ,f_Ca+K1f_Ca/2.0)
    K2Ca_SR = dtime*func_CaSRA(tt+dtime/2.0,  Ca_i+K1Ca_i/2.0,  Ca_SR+K1Ca_SR/2.0,  d+K1d/2.0,  g+K1g/2.0)
    K2Na_i = dtime*func_NaiA(tt+dtime/2.0,  m+K1m/2.0,  h+K1h/2.0,  j+K1j/2.0,  V+K1V/2.0,  Na_i+K1Na_i/2.0,  Ca_i+K1Ca_i/2.0)
    K2V = dtime*func_VA(tt+dtime/2.0,  V+K1V/2.0,  r+K1r/2.0,  q+K1q/2.0,  xr1+K1xr1/2.0,  xr2+K1xr2/2.0,  xs+K1Xs/2.0,  Ca_i+K1Ca_i/2.0,  d+K1d/2.0,  f1+K1f1/2.0,  f2+K1f2/2.0,  f_Ca+K1f_Ca/2.0,  Na_i+K1Na_i/2.0,  m+K1m/2.0,  h+K1h/2.0,  j+K1j/2.0  ,xf+K1Xf/2.0)

    !! K3
    K3h = dtime*func_hA(tt+dtime/2.0,  V+K2V/2.0,  h+K2h/2.0)
    K3j = dtime*func_jA(tt+dtime/2.0,  V+K2V/2.0,  j+K2j/2.0)
    K3m =     dtime*func_mA(tt+dtime/2.0,  V+K2V/2.0,  m+K2m/2.0)
    K3d =     dtime*func_dA(tt+dtime/2.0,  V+K2V/2.0,  d+K2d/2.0)
    K3f_Ca =  dtime*func_fCaA(tt+dtime/2.0,  f_Ca+K2f_Ca/2.0,  Ca_i+K2Ca_i/2.0 ,V+K2V/2.0)
    K3f1 = dtime*func_f1A(tt+dtime/2.0,  V+K2V/2.0,  f1+K2f1/2.0, Ca_i+K2f1/2.0,df1)
    K3f2 =    dtime*func_f2A(tt+dtime/2.0,  V+K2V/2.0,  f2+K2f2/2.0)
    K3r =     dtime*func_rA(tt+dtime/2.0,  V+K2V/2.0,  r+K2r/2.0)
    K3q =     dtime*func_qA(tt+dtime/2.0,  V+K2V/2.0,  q+K2q/2.0)
    K3xr1 =   dtime*func_xr1A(tt+dtime/2.0,V+K2V/2.0,  xr1+K2xr1/2.0)
    K3xr2 =  dtime* func_xr2A(tt+dtime/2.0,  V+K2V/2.0,  xr2+K2xr2/2.0)
    K3Xs =    dtime*func_XsA(tt+dtime/2.0,  V+K2V/2.0,  xs+K2Xs/2.0)
    K3Xf =    dtime*func_XfA(tt+dtime/2.0,  V+K2V/2.0,  xf+K2Xf/2.0)
    K3g =     dtime*func_gA(tt+dtime/2.0,  Ca_i+K2Ca_i/2.0,  g+K2g/2.0 ,V+K2V/2.0)
    K3Ca_i =  dtime*func_CaiA(tt+dtime/2.0,  Ca_i+K2Ca_i/2.0,  Ca_SR+K2Ca_SR/2.0,  d+K2d/2.0,  g+K2g/2.0,  V+K2V/2.0,  Na_i+K2Na_i/2.0,  f1+K2f1/2.0  ,f2+K2f2/2.0  ,f_Ca+K2f_Ca/2.0)
    K3Ca_SR = dtime*func_CaSRA(tt+dtime/2.0,  Ca_i+K2Ca_i/2.0,  Ca_SR+K2Ca_SR/2.0,  d+K2d/2.0,  g+K2g/2.0)
    K3Na_i =  dtime*func_NaiA(tt+dtime/2.0,  m+K2m/2.0,  h+K2h/2.0,  j+K2j/2.0,  V+K2V/2.0,  Na_i+K2Na_i/2.0,  Ca_i+K2Ca_i/2.0)
    K3V =     dtime*func_VA(tt+dtime/2.0,  V+K2V/2.0,  r+K2r/2.0,  q+K2q/2.0,  xr1+K2xr1/2.0,  xr2+K2xr2/2.0,  xs+K2Xs/2.0,  Ca_i+K2Ca_i/2.0,  d+K2d/2.0,  f1+K2f1/2.0,  f2+K2f2/2.0,  f_Ca+K2f_Ca/2.0,  Na_i+K2Na_i/2.0,  m+K2m/2.0,  h+K2h/2.0,  j+K2j/2.0  ,xf+K2Xf/2.0)

    !! K4
    K4h = dtime*func_hA(tt+dtime,  V+K3V,  h+K3h)
    K4j = dtime*func_jA(tt+dtime,  V+K3V,  j+K3j)
    K4m = dtime*func_mA(tt+dtime,  V+K3V,  m+K3m)
    K4d = dtime*func_dA(tt+dtime,  V+K3V,  d+K3d)
    K4f_Ca =  dtime*func_fCaA(tt+dtime,  f_Ca+K3f_Ca,  Ca_i+K3Ca_i ,V+K3V)
    K4f1 = dtime*func_f1A(tt+dtime,  V+K3V,  f1+K3f1, Ca_i+K3f1 ,df1)
    K4f2 =    dtime*func_f2A(tt+dtime,  V+K3V,  f2+K3f2)
    K4r =     dtime*func_rA(tt+dtime,  V+K3V,  r+K3r)
    K4q =     dtime*func_qA(tt+dtime,  V+K3V,  q+K3q)
    K4xr1 =   dtime*func_xr1A(tt+dtime, V+K3V,  xr1+K3xr1)
    K4xr2 =   dtime*func_xr2A(tt+dtime,  V+K3V,  xr2+K3xr2)
    K4Xs =    dtime*func_XsA(tt+dtime,  V+K3V,  xs+K3Xs)
    K4Xf =    dtime*func_XfA(tt+dtime,  V+K3V,  xf+K3Xf)
    K4g =     dtime*func_gA(tt+dtime,  Ca_i+K3Ca_i,  g+K3g ,V+K3V)
    K4Ca_i =  dtime*func_CaiA(tt+dtime,  Ca_i+K3Ca_i,  Ca_SR+K3Ca_SR,  d+K3d,  g+K3g,  V+K3V,  Na_i+K3Na_i,  f1+K3f1  ,f2+K3f2  ,f_Ca+K3f_Ca)
    K4Ca_SR = dtime*func_CaSRA(tt+dtime,  Ca_i+K3Ca_i,  Ca_SR+K3Ca_SR,  d+K3d,  g+K3g)
    K4Na_i =  dtime*func_NaiA(tt+dtime,  m+K3m,  h+K3h,  j+K3j,  V+K3V,  Na_i+K3Na_i,  Ca_i+K3Ca_i)
    K4V =     dtime*func_VA(tt+dtime,  V+K3V,  r+K3r,  q+K3q,  xr1+K3xr1,  xr2+K3xr2,  xs+K3Xs,  Ca_i+K3Ca_i,  d+K3d,  f1+K3f1,  f2+K3f2,  f_Ca+K3f_Ca,  Na_i+K3Na_i,  m+K3m,  h+K3h,  j+K3j  ,xf+K3Xf)

    !!Gates, Concentrations and Voltage
    vauxi_exm(1,ipoin,1) = h + (K1h + 2.0*K2h + 2.0*K3h + K4h)*(1.0/6.0)
    vauxi_exm(2,ipoin,1) = j + (K1j + 2.0*K2j + 2.0*K3j + K4j)*(1.0/6.0)
    vauxi_exm(3,ipoin,1) = m + (K1m + 2.0*K2m + 2.0*K3m + K4m)*(1.0/6.0)
    vauxi_exm(4,ipoin,1) = d + (K1d + 2.0*K2d + 2.0*K3d + K4d)*(1.0/6.0)
    vauxi_exm(5,ipoin,1) = f_Ca + (K1f_Ca + 2.0*K2f_Ca + 2.0*K3f_Ca + K4f_Ca)*(1.0/6.0)
    vauxi_exm(6,ipoin,1) = f1 + (K1f1 + 2.0*K2f1 + 2.0*K3f1 + K4f1)*(1.0/6.0)
    vauxi_exm(7,ipoin,1) = f2 + (K1f2 + 2.0*K2f2 + 2.0*K3f2 + K4f2)*(1.0/6.0)
    vauxi_exm(8,ipoin,1) = r + (K1r + 2.0*K2r + 2.0*K3r + K4r)*(1.0/6.0)
    vauxi_exm(9,ipoin,1) = q + (K1q + 2.0*K2q + 2.0*K3q + K4q)*(1.0/6.0)
    vauxi_exm(10,ipoin,1) = xr1 + (K1xr1 + 2.0*K2xr1 + 2.0*K3xr1 + K4xr1)*(1.0/6.0)
    vauxi_exm(11,ipoin,1) = xr2 + (K1xr2 + 2.0*K2xr2 + 2.0*K3xr2 + K4xr2)*(1.0/6.0)
    vauxi_exm(12,ipoin,1) = xs + (K1Xs + 2.0*K2Xs + 2.0*K3Xs + K4Xs)*(1.0/6.0)
    vauxi_exm(13,ipoin,1) = xf + (K1Xf + 2.0*K2Xf + 2.0*K3Xf + K4Xf)*(1.0/6.0)
    vauxi_exm(14,ipoin,1) = g + (K1g + 2.0*K2g + 2.0*K3g + K4g)*(1.0/6.0)
    vconc(1,ipoin,1) = Ca_i + (K1Ca_i + 2.0*K2Ca_i + 2.0*K3Ca_i + K4Ca_i)*(1.0/6.0)
    vconc(2,ipoin,1) = Ca_SR + (K1Ca_SR + 2.0*K2Ca_SR + 2.0*K3Ca_SR + K4Ca_SR)*(1.0/6.0)
    vconc(3,ipoin,1) = Na_i + (K1Na_i + 2.0*K2Na_i + 2.0*K3Na_i + K4Na_i)*(1.0/6.0)
    xioni=-func_VA(tt,V,r,q,xr1,xr2,xs,Ca_i,d,f1,f2,f_Ca,Na_i,m,h,j,xf)
    !xioni = V + (K1V + 2.0*K2V + 2.0*K3V + K4V)*(1.0/6.0)
    df1 = (vauxi_exm(6,ipoin,1)-vauxi_exm(6,ipoin,2))/dtime
    dioni = 0.0
    xioni = xioni * 1000.0_rp
    !write(*,*) df1
    
    vconc(:,ipoin,3)=vconc(:,ipoin,2)         
    vconc(:,ipoin,2)=vconc(:,ipoin,1) 
    vauxi_exm(:,ipoin,2)=vauxi_exm(:,ipoin,1)   

  
  !open(unit=1, file='voltage_Atrial.dat')
  !do z=1,counter
  !    write(1,*) tspan(z), elmag(z)
  !end do
  !close(unit=1)

end subroutine exm_scatri_ionicurrents


!################################# functions are defined here#########################
!#####################################################################################
real(rp) function func_CaiA(tt,Ca_i,Ca_SR,d,g,V,Na_i,f1,f2,f_Ca)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, Ca_i, Ca_SR, d, g, V, Na_i, f1, f2, f_Ca
  real(rp) :: Ca_ibufc, E_Ca, I_CaL, I_NaCa, I_bCa, I_pCa, I_rel, I_up, I_leak

  !! Reversal potentials
  E_Ca = (0.5*Rcons*Temp/F)*log(Ca_o/Ca_i)

  !! L-type Ca2+ current, I_CaL
  I_CaL = (g_CaL*4.0*V*F**2.0)*(Ca_i*exp(2.0*V*F/(Rcons*Temp))-0.341*Ca_o)*d*f1*f2*f_Ca/(Rcons*Temp*(exp(2.0*V*F/(Rcons*Temp))-1.0))  
  !! Na+/Ca2+ exchanger current, I_NaCa
  I_NaCa = K_NaCa*(exp(Gamma*V*F/(Rcons*Temp))*Na_i**3.0*Ca_o-exp((Gamma-1.0)*V*F/(Rcons*Temp))*Na_o**3.0*Ca_i*Alpha)/((Km_Na_i**3.0+Na_o**3.0)*(Km_Ca+Ca_o)*(1.0+K_sat*exp((Gamma-1.0)*V*F/(Rcons*Temp))))
    
  !! background currents
  I_bCa = g_bCa*(V-E_Ca)

  !! Ca2+ pump current, I_pCa
  I_pCa = g_pCa*Ca_i/(Ca_i+K_pCa)
  !! Ca2+ dynamics
  !Atrial
  I_rel = ((a_rel*(Ca_SR**2.0)/(b_rel**2.0+Ca_SR**2.0))+c_rel)*d*g*(0.0556)


  I_up = Vmax_up/(1.0+K_up**2.0/Ca_i**2.0)
  I_leak = V_leak*(Ca_SR-Ca_i)

  Ca_ibufc = 1.0/(1.0+Buf_c*K_Buf_c/(Ca_i+K_Buf_c)**2.0)
  func_CaiA = Ca_ibufc*(I_leak-I_up+I_rel-(I_CaL+I_bCa+I_pCa-2.0*I_NaCa)*C_m/(2.0*V_c*F*1.0e-18))
end function

!#####################################################################################
real(rp) function func_CaSRA(tt,Ca_i,Ca_SR,d,g)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, Ca_i, Ca_SR, d, g
  real(rp) :: Ca_srbufsr, I_rel, I_up, I_leak

  !! Ca2+ dynamics
  !Atrial
  I_rel = ((a_rel*(Ca_SR**2.0)/(b_rel**2.0+Ca_SR**2.0))+c_rel)*d*g*(0.0556)

  I_up = Vmax_up/(1.0+K_up**2.0/Ca_i**2.0)
  I_leak = V_leak*(Ca_SR-Ca_i)

  !!
  Ca_srbufsr = 1.0/(1.0+Buf_sr*K_Buf_sr/(Ca_SR+K_Buf_sr)**2.0)
  func_CaSRA = Ca_srbufsr*V_c*(I_up-(I_rel+I_leak))/V_sr

end function

!#####################################################################################
real(rp) function func_dA(tt,V,d)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, d

  real(rp) :: d_inf, Alpha_d, beta_d, Gamma_d, taw_d
  !I_CaL, d gate
  !Atrial
  d_inf = 1.0/(1.0+exp(-(V*1000.0+5.986)/7.0))

  Alpha_d = (1.4/(1.0+exp((-35.0-V*1000.0)/13.0)))+0.25
  beta_d = 1.4/(1.0+exp((V*1000.0+5.0)/5.0))
  Gamma_d = 1.0/(1.0+exp((50.0-V*1000.0)/20.0))
  taw_d = (Alpha_d*beta_d+Gamma_d)/1000.0

  func_dA = (d_inf-d)/taw_d;
end function 

!#####################################################################################
real(rp) function  func_f1A(tt,V,f1,Ca_i,df1)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, f1, Ca_i, df1
  
  real(rp) :: f1_inf, taw_f1

  !! I_CaL, f1 gate
  !Atrial
  f1_inf = 1.0/(1.0+exp((V*1000.0+25.226)/3.0))

  if (df1>0.0) then
    taw_f1 = ((1102.5*(exp(-(((V*1000.0+27.0)**2.0)/15.0)**2.0))+200.0/(1.0+exp((13.0-V*1000.0)/10.0))+180.0/(1.0+exp((30.0+V*1000.0)/10.0))+20.0)*(1.0+1433.0*(Ca_i-50.0*1e-6)))/1000.0
  else
    taw_f1 = ((1102.58*(exp(-(((V*1000.0+27.0)**2.0)/15.0)**2.0))+200.0/(1.0+exp((13.0-V*1000.0)/10.0))+180.0/(1.0+exp((30.0+V*1000.0)/10.0))+20.0))/1000.0
  end if

  func_f1A = (f1_inf-f1)/taw_f1;

end function


!#####################################################################################
real(rp) function func_f2A(tt,V,f2)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, f2
  
  real(rp) :: f2_inf, taw_f2

  !! I_CaL, f2 gate
  !Atrial
  f2_inf = (0.67/(1.0+exp((V*1000.0+31.226)/4.0)))+0.33
  taw_f2 = ((600.0*(exp((-(V*1000.0+25.0)**2.0)/170.0))+31.0/(1.0+exp((25.0-V*1000.0)/10.0))+16.0/(1.0+exp((30.0+V*1000.0)/10.0)))*2.0)/1000.0
  func_f2A = (f2_inf-f2)/taw_f2;
end function

!#####################################################################################
! From here on, functions V = A (HISPCM)
!#####################################################################################

real(rp) function func_fCaA(tt,f_Ca,Ca_i,V)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, f_Ca, Ca_i, V
  
  real(rp) :: Alpha_fCa, beta_fCa, Gamma_fCa, f_Ca_inf, constf_Ca
  real(rp), parameter :: taw_fCa = 0.002
  
  !! I_CaL, f_Ca gate
  Alpha_fCa = 1.0/(1.0+(Ca_i/0.0006)**8.0)
  beta_fCa = 0.1/(1.0+exp((Ca_i-0.0009)/0.0001))
  Gamma_fCa = 0.3/(1.0+exp((Ca_i-0.00075)/0.0008))
  f_Ca_inf = (Alpha_fCa+beta_fCa+Gamma_fCa)/1.3156
  
  if ((f_Ca_inf>f_Ca) .and. (V>-0.06)) then
    constf_Ca = 0.0
  else
    constf_Ca = 1.0
  end if

  func_fCaA = constf_Ca*(f_Ca_inf-f_Ca)/taw_fCa
end function

!#####################################################################################
real(rp) function  func_XsA(tt,V,Xs)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, Xs
  
  real(rp) :: xs_inf, Alpha_xs, beta_xs, taw_xs

  !I_Ks, Xs gate
  xs_inf = 1.0/(1.0+exp((-V*1000.0-20.0)/16.0))
  Alpha_xs =1100.0/sqrt(1.0+exp((-V*1000.0-10.0)/6.0))
  beta_xs = 1.0/(1.0+exp((V*1000.0-60.0)/20.0))
  taw_xs = (Alpha_xs*beta_xs)/1000.0

  func_XsA = (xs_inf-Xs)/taw_xs
end function

!#####################################################################################
real(rp) function func_gA(tt,Ca_i,g,V)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, Ca_i, g, V
  real(rp) :: g_inf
  real(rp), parameter :: taw_g = 0.002
  real(rp) :: constii

  if (Ca_i<=0.00035) then
    g_inf = 1.0/(1.0+(Ca_i/0.00035)**6.0)
  else
    g_inf = 1.0/(1.0+(Ca_i/0.00035)**16.0)
  end if

  if ((g_inf>g) .and. (V>-0.06)) then
    constii = 0.0
  else
    constii = 1.0
  end if

  func_gA = constii*(g_inf-g)/taw_g
end function

!#####################################################################################
real(rp) function func_hA(tt,V,h)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, h
  
  real(rp) :: h_inf, Alpha_h, beta_h, taw_h

  !! I_Na, h gate
  h_inf = 1.0/(sqrt(1.0+exp((V*1000.0+72.1)/5.7)))

  if (V<-0.04) then
    Alpha_h = 0.057*exp(-(V*1000.0+80.0)/6.8)
    beta_h = 2.7*exp(0.079*V*1000.0)+3.1*(10.0**5.0)*exp(0.3485*V*1000.0)
  else
    Alpha_h = 0.0
    beta_h = 0.77/(0.13*(1.0+exp((1000.0*V+10.66)/(-11.1))))
  end if
  if (V<-0.04) then
    taw_h = 1.5/((Alpha_h+beta_h)*1000.0)
  else
    taw_h = 2.542/1000.0
  end if
  func_hA = (h_inf-h)/taw_h
end function

!#####################################################################################
real(rp) function func_jA(tt,V,j)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, j
  
  real(rp) :: j_inf, Alpha_j, beta_j, taw_j

  !! I_Na, j gate
  j_inf = 1.0/(sqrt(1.0+exp((1000.0*V+72.1)/5.7)))

  if (V<-0.04) then
    Alpha_j = (-25428.0*exp(0.2444*V*1000.0)-6.948*(1e-6)*exp(-0.04391*V*1000.0))*(V*1000.0+37.78)/(1.0+exp(0.311*(V*1000.0+79.23)))
    beta_j =  0.02424*exp(-0.01052*V*1000.0)/(1.0+exp(-0.1378*(V*1000.0+40.14)))
  else
    Alpha_j = 0.0
    beta_j =  0.6*exp(0.057*V*1000.0)/(1.0+exp(-0.1*(V*1000.0+32.0)))
  end if

  taw_j = 7.0/((Alpha_j+beta_j)*1000.0)
  func_jA = (j_inf-j)/taw_j
end function

!#####################################################################################
real(rp) function func_mA(tt,V,m)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, m
  
  real(rp) :: m_inf, Alpha_m, beta_m, taw_m

  !! I_Na, m gate
  m_inf =  (1.0+exp((-34.1-V*1000.0)/5.9))**(-1.0/3.0)
  Alpha_m =  1.0/(1.0+exp((-60.0-V*1000.0)/5.0))
  beta_m =  0.1/(1.0+exp((V*1000.0+35.0)/5.0))+0.1/(1.0+exp((V*1000.0-50.0)/200.0))
  taw_m = (Alpha_m*beta_m)/1000.0

  func_mA = (m_inf-m)/taw_m
end function


!#####################################################################################
real(rp) function func_NaiA(tt,m,h,j,V,Na_i,Ca_i)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, m, h, j, V, Na_i, Ca_i
  real(rp) :: E_Na, I_bNa, I_Na, I_NaK, I_NaCa
 
  !! Reversal potentials
  E_Na = (Rcons*Temp/F)*log(Na_o/Na_i)

  !! background currents
  I_bNa = g_bNa*(V-E_Na)
  I_Na = g_Na*(m**3.0)*h*j*(V-E_Na)

  !! Na+/K+ pump current, I_NaK
  I_NaK = ((P_NaK*K_o*Na_i/(K_o+K_mk))/(Na_i+K_mNa))/(1.0+0.1245*exp(-0.1*V*F/(Rcons*Temp))+0.353*exp(-V*F/(Rcons*Temp)))

  !! Na+/Ca2+ exchanger current, I_NaCa
  I_NaCa = K_NaCa*(exp(Gamma*V*F/(Rcons*Temp))*(Na_i**3.0)*Ca_o-exp((Gamma-1.0)*V*F/(Rcons*Temp))*(Na_o**3.0)*Ca_i*Alpha)/((Km_Na_i**3.0+Na_o**3.0)*(Km_Ca+Ca_o)*(1.0+K_sat*exp((Gamma-1.0)*V*F/(Rcons*Temp))))

  !! sodium dynamics
  func_NaiA = -C_m*(I_Na+I_bNa+3.0*I_NaK+3.0*I_NaCa)/(F*V_c*1e-18)
end function

!#####################################################################################
real(rp) function func_qA(tt,V,q)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, q
  
  real(rp) :: q_inf, taw_q

  !! I_to, q gate
  q_inf = 1.0/(1.0+exp((V*1000.0+53.0)/13.0))
  taw_q = (39.102/(0.57*exp(-0.08*(V*1000.0+44.0))+0.065*exp(0.1*(V*1000.0+45.93)))+6.06)/1000.0

  func_qA = (q_inf-q)/taw_q
end function

!#####################################################################################
real(rp) function func_rA(tt,V,r)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, r
  
  real(rp) :: r_inf, taw_r

  !I_to, r gate
  r_inf = 1.0/(1.0+exp((22.3-V*1000.0)/18.75))
  taw_r = (14.405516/(1.037*exp(0.09*(V*1000.0+30.61))+0.369*exp(-0.12*(V*1000.0+23.84)))+2.75352)/1000.0

  func_rA = (r_inf-r)/taw_r
end function

!#####################################################################################
real(rp) function func_VA(tt,V,r,q,xr1,xr2,Xs,Ca_i,d,f1,f2,f_Ca,Na_i,m,h,j,Xf)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, r, q, xr1, xr2, &
  Xs, Ca_i, d, f1, f2, f_Ca, Na_i, m, h, j, Xf
  
  real(rp) :: Alpha_K1, beta_K1, x_K1_inf, E_Na, E_K, E_Ks, E_Ca, E_f, &
      I_CaL, I_to, I_Kr, I_Ks, I_K1, I_f, I_NaK, I_NaCa, I_bNa, I_bCa, I_Na, &
      I_pCa

  !! Reversal potentials
  E_Na = (Rcons*Temp/F)*log(Na_o/Na_i)
  E_K = (Rcons*Temp/F)*log(K_o/K_i)
  E_Ks = (Rcons*Temp/F)*log((K_o+P_kna*Na_o)/(K_i+P_kna*Na_i))
  E_Ca = (0.5*Rcons*Temp/F)*log(Ca_o/Ca_i)
  E_f = -0.017


  !! L-type Ca2+ current, I_CaL
I_CaL = (g_CaL*4.0*V*F**2.0)*(Ca_i*exp(2.0*V*F/(Rcons*Temp))-0.341*Ca_o)*d*f1*f2*f_Ca/(Rcons*Temp*(exp(2.0*V*F/(Rcons*Temp))-1.0))
  !!
  !! Transient outward current, I_to
  I_to = g_to*r*q*(V-E_K)

  !! Rapid delayed rectifier K+ current, I_Kr
  I_Kr = g_Kr*sqrt(K_o/5.4)*xr1*xr2*(V-E_K)
  !! Slow delayed rectifier K+ current, I_Ks
  I_Ks = g_Ks*Xs**2.0*(1.0+0.6/(1.0+(3.8*1e-5/Ca_i)**1.4))*(V-E_Ks)

  !! Inward rectifier K+ current, I_K1
  Alpha_K1 = 3.91/(1+exp(0.5942*(V*1000.0-E_K*1000.0-200.0)))
  beta_K1 = (-1.509*exp(0.0002*(V*1000.0-E_K*1000.0+100.0))+exp(0.5886*(V*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(V*1000.0-E_K*1000.0)))
  x_K1_inf = Alpha_K1/(Alpha_K1+beta_K1)
  I_K1 = g_K1*x_K1_inf*sqrt(K_o/5.4)*(V-E_K)

  !! Hyperpolarization activated funny current, I_f
  I_f = g_f*Xf*(V-E_f)

  !! Na+/K+ pump current, I_NaK
  I_NaK = ((P_NaK*K_o*Na_i/(K_o+K_mk))/(Na_i+K_mNa))/(1.0+0.1245*exp(-0.1*V*F/(Rcons*Temp))+0.353*exp(-V*F/(Rcons*Temp)));

  !! Na+/Ca2+ exchanger current, I_NaCa
I_NaCa = K_NaCa*(exp(Gamma*V*F/(Rcons*Temp))*Na_i**3.0*Ca_o-exp((Gamma-1.0)*V*F/(Rcons*Temp))*Na_o**3.0*Ca_i*Alpha)/((Km_Na_i**3.0+Na_o**3.0)*(Km_Ca+Ca_o)*(1.0+K_sat*exp((Gamma-1.0)*V*F/(Rcons*Temp))))

  !! background currents
  I_bNa = g_bNa*(V-E_Na)
  I_bCa = g_bCa*(V-E_Ca)
  I_Na = g_Na*m**3.0*h*j*(V-E_Na)

  !! Ca2+ pump current, I_pCa
  I_pCa = g_pCa*Ca_i/(Ca_i+K_pCa)

  !! Membrane potential
  func_VA = -(I_K1+I_to+I_Kr+I_Ks+I_CaL+I_NaK+I_Na+I_NaCa+I_pCa+I_f+I_bNa+I_bCa)
end function

!#####################################################################################
real(rp) function func_XfA(tt,V,Xf)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, Xf
  
  real(rp) :: xf_inf, taw_f

  !I_f, X_f gate
  xf_inf = 1.0/(1.0+exp((V*1000.0+77.85)/5.0))
  taw_f = 1900.0/((1.0+exp((V*1000.0+15.0)/10.0))*1000.0)

  func_XfA = (xf_inf-Xf)/taw_f
end function

!#####################################################################################
real(rp) function func_xr1A(tt,V,xr1)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, xr1

  real(rp), parameter :: Qcons = 2.3
  real(rp) :: V_half, xr1_inf, Alpha_xr1, beta_xr1, taw_xr1

  V_half = 1000.0*(-(Rcons*Temp/(F*Qcons))*log(((1.0+Ca_o/2.6)**4.0)/(L_0*(1.0+Ca_o/0.58)**4.0))-0.019)
  xr1_inf = 1.0/(1.0+exp((V_half - V*1000.0)/4.9))
  Alpha_xr1 = 450.0/(1.0+exp((-45.0-V*1000.0)/10.0))
  beta_xr1 = 6.0/(1.0+exp((V*1000.0+30.0)/11.5))
  taw_xr1 = (Alpha_xr1*beta_xr1)/1000.0

  func_xr1A = (xr1_inf-xr1)/taw_xr1
end function

!#####################################################################################
real(rp) function func_xr2A(tt,V,xr2)
  use mod_exm_hipscatriaconst

  real(rp), intent(in) :: tt, V, xr2
  
real(rp) :: xr2_inf, Alpha_xr2, beta_xr2, taw_xr2

  !! I_Kr, xr2 gate
  xr2_inf = 1.0/(1.0+exp((V*1000.0+88.0)/50.0))
  Alpha_xr2 = 3.0/(1.0+exp((-V*1000.0-60.0)/20.0))
  beta_xr2 = 1.12/(1.0+exp((V*1000.0-60.0)/20.0))
  taw_xr2 = (Alpha_xr2*beta_xr2)/1000.0

  func_xr2A = (xr2_inf-xr2)/taw_xr2
end function
