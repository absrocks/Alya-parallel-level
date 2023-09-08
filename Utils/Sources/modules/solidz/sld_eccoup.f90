!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_eccoup.f90
!> @author  Mariano Vazquez and Francesc Levrero-Florencio
!> @date    16/11/1966
!> @brief   Calculate the active (2nd P-K and cauchy) active stress 
!> @details Calculate the active (2nd P-K and cauchy) active stress 
!> @} 
!-----------------------------------------------------------------------
subroutine sld_eccoup(&
     pgaus,pmate,igaus,gpdet,gptlo,lamda,nfibt,nshtt,normt,statvar,&
     gpigd,castr_active,gpstr_active,gpcac,tactf,ielem,jac,&
     param)

  !-----------------------------------------------------------------------
  !
  !    GPSTR_active ... 2nd P-K Active Stress tensor ....................S_a
  !    CASTR_active ... Cauchy Stress tensor .................... sigma_a
  !
  ! USES
  !
  ! USED BY
  !    sld_stress_model_134
  !    sld_stress_model_xxx  
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip,rp
  use def_domain, only : ndime
  use def_master, only : cutim,kfl_cellmod,CELL_TT2006_EXMEDI,CELL_OHARA_EXMEDI,                  &
   CELL_FITZHUGH_EXMEDI,CELL_NOMOD_EXMEDI,ittim,dtime,kfl_eccty, cell_ca0_ecc
  use def_solidz           

  implicit none
  integer(ip), intent(in) :: pgaus,pmate,igaus,ielem
  integer(ip) :: idime,jdime,kdime,ldime,inewton,icomp,max_nr

  real(rp), intent(in) :: gptlo(pgaus),lamda(3),gpigd(ndime,ndime,pgaus),gpdet(pgaus),            &
   nfibt(ndime,pgaus),nshtt(ndime,pgaus),normt(ndime,pgaus),gpcac

  real(rp), intent(inout) :: statvar(6,2)

  real(rp), intent(out) :: castr_active(ndime,ndime),gpstr_active(ndime,ndime),tactf,Jac(6,6),    &
   param(25)
  
  real(rp) :: ca0,cam,thau,cca2p,n,beta,Cal50,tmaxi,tacts,tactn,tacis,k,fporf_t(ndime,ndime),     &
   spors_t(ndime,ndime),nporn_t(ndime,ndime),nfibe_length,lambda_aux,ydot(6),I_6(6,6),res(6),     &
   sol(6),lambda_rate,ort_act_k1,lambda_prev,lambda_prev_prev,dtime_aux,S,W,CaTRPN,B,zeta_s,      &
   k_TRPN,n_TRPN,Ca50_ref,n_Tm,TRPN_50,k_uw,k_ws,r_s,r_w,y_s,y_w,phi,A_eff,beta_0,beta_1,T_ref,   &
   k_wu,k_su,A_s,c_s,c_w,y_su,y_wu,Ca50,k_b,k_u,k_u_2,h,norm_res,res_tol,gpcac_aux,g,Q1,Q2,TRPN,  &
   XB,Fn,n_xb,A1,A2,alpha_1,alpha_2,k_xb,perm_tot,hprime,Q,a_vel,d_mxb,Jac_inv(6,6),term,         &
   ort_act_k2,zeta_w,prestress

  !constants for active stress (units: CGS)

  !!!
  ! Hunter-McCulloch ECC model
  !!!
  if (kfl_eccty(pmate) == 1) then

     ! Ca0 is not necesarily zero for the EP model.
     ! An average of the three Ca0 concentrations is set for the whole tissue.
     ! in this way we avoid the gather operation.
     ! The 1000 factor is for unit conversion.
     ca0=sum(cell_ca0_ecc(:,pmate))*1000.0_rp/3.0_rp

     cam = 1.0_rp !mmole

     !Default value: 0.06_rp sec (human) ... !0.01268 sec Rat
     thau = timec_sld(pmate) 

     !tmaxi= 1.0e6_rp !0.2e6_rp !barye !! 1.0e-3!1.0e-3 !stress in msec, cm
     tmaxi = 1.0e6_rp*cocof_sld(pmate)
     !  JAZZY DID THIS TO SPEED THINGS UP AND DEBUG!!!  
     if (kfl_coupt_sld(pmate) == 11) then 
        tmaxi = tmaxi*covar_sld(ielem)
     end if
     !!  END 
     beta = 1.45_rp  !!1.45_rp  changed to value from paper: JAZ 2.0
     n = hillc_sld(pmate)
     k = 0.3_rp

     castr_active = 0.0_rp
     gpstr_active = 0.0_rp

     ! * * * * * *    
     ! Active stress T(lamda, [Ca++])m x m  (Hunter p 688)
     ! * * * * * * 

     !calculate [Ca2+] 
     if ((kfl_cellmod(pmate) == CELL_TT2006_EXMEDI) .or.                                          &
      (kfl_cellmod(pmate) == CELL_OHARA_EXMEDI)) then   !! Ten Tuscher calcula el Ca2+ dins del propi model
        !
        ! Ten Tuscher model
        !
        cca2p = gpcac    !*1000  Changed by JAZ, was wrong!!
        !C50
        Cal50 = cal50_sld(pmate)!0.6 !value by Hung 
        ! Cal50= 4.35_rp/sqrt(exp(4.75_rp*((1.9*lamda(1))-1.58_rp))-1.0_rp) 
     else if (kfl_cellmod(pmate) == CELL_FITZHUGH_EXMEDI) then  !! FHN
        !
        ! Otherwise
        !
        cca2p= ca0 + (cam-ca0)*(gptlo(igaus)/thau)*exp(1.0_rp-(gptlo(igaus)/thau))
        !C50
        Cal50 = cal50_sld(pmate)!0.6 !value by Hung 
        !     Cal50= 0.5_rp!0.6 !value by Hung 
        ! Cal50= 4.35_rp/sqrt(exp(4.75_rp*((1.9*lamda(1))-1.58_rp))-1.0_rp) 
     else if (kfl_cellmod(pmate) == CELL_NOMOD_EXMEDI)  then  !
        !
     else 
        call runend ("SLD_ECCOUP: CELL MODEL NOT INCLUDED IN THE COUPLING WITH SOLIDZ MODULE")
     end if

     !write(990,*) cca2p
     !  if (tactf == 1) then  !Hunter model for active force (PAULA)
     !T, active tension in ff

     cca2p=cca2p-ca0
     tactf = (cca2p**n/(cca2p**n+Cal50**n))*tmaxi*(1.0_rp+beta*(lamda(1)-1.0_rp))
     !  tacts = 0.0_rp!tactf!0.0!k*tactf*(lamda(2)/lamda(1))
     !  tactn = 0.0_rp!tactf!0.0!k*tactf*(lamda(3)/lamda(1))
     !  else                  !Rice myofilament model (PAULA) ------> if (tactf == 1) aún por mejorar!
     !Cal50 here should be 0.5
     !     call sld_myofil(cca2p, Factive)
     !     tactf = Factive            !Gamma for tensile active force is necessary in this case? (PAULA)
     !  end if

  else if (kfl_eccty(pmate) == 2) then

     ! RICE: TO BE PROGRAMMED

  else if (kfl_eccty(pmate) == 3 .or. kfl_eccty(pmate)== 4) then

    ! ----------------------------------------------------------------------
    ! This is the Electromechanical coupling model of Land et al (2017), 
    ! solved with Forward Euler.
    !
    ! The system consists of six ODEs, which are
    !  dS/dt      = k_ws*W-k_su*S-y_su*S
    !  dW/dt      = k_uw*U-k_wu*W-k_ws*W-y_wu*W
    !  dCaTRPN/dt = k_TRPN*(((Ca_i/Ca_T50)^n_TRPN)*(1-CaTRPN)-CaTRPN)
    !  dB/dt      = k_b*CaTRPN^(-n_tm/2)*U-k_u*CaTRPN^(n_tm/2)*B
    !  dzeta_s/dt = A_s*(dlambda/dt)-c_s*zeta_s
    !  dzeta_w/dt = A_w*(dlambda/dt)-c_w*zeta_w
    !
    ! With an active tension (and its length-dependence) of
    !  T_a = h*(T_ref/r_s)*(S*(zeta_s+1)+W*zeta_w)
    !  h   = 0                             if lambda<=-1/(2*beta_0-1.87)
    !      = 1+2*beta_0*lambda-1.87*beta_0 if (lambda>(-1/(2*beta_0-1.87)))and(lambda<=0.87)
    !      = 1+beta_0*lambda-beta_0        if (lambda>0.87)and(lambda<=1.2)
    !      = 1+0.2*beta_0                  if lambda>1.2
    !
    ! and a calcium half-activation point of
    !  Ca_T50 = Ca_ref+beta_1*(lambda-1) if lambda<=1.2
    !         = Ca_ref+0.2*beta-1        if lambda>1.2
    !
    ! where the variables are defined as
    !  lambda = ||F*f_0|| (f_0 is the fibre direction in reference configuration)
    !  U      = 1-B-S-W
    !  y_wu   = y_w*|zeta_w|
    !  y_su   = y_s*(-zeta_s-1) if zeta_s+1<0
    !         = y_s*zeta_s      if zeta_s+1>1
    !         = 0               otherwise
    !  A_s    = A_w = A_eff*r_s/((1-r_s)*r_w+r_s)
    !  k_wu   = k_uw*r_w*(1/r_w-1)-k_ws
    !  k_su   = k_ws*r_w*(1/r_s-1)
    !  k_b    = k_u*(CaTRPN^n_tm)/(1-r_s-(1-r_s)r_w)
    !  c_w    = phi*k_uw*U/W = phi*k_uw*((1-r_s)*(1-r_w))/((1-r_s)*r_w)
    !  c_s    = phi*k_ws*W/S = phi*k_ws*((1-r_s)*r_w)/r_s
    !
    ! and the parameters for these variables are given by
    !  k_TRPN          = 0.1   (1/ms)
    !  n_TRPN          = 2
    !  [Ca^2+]^ref_T50 = 2.5   (muM)
    !  k_u             = 1     (1/ms)
    !  n_tm            = 5
    !  TRPN_50         = 0.35
    !  k_uw            = 0.182 (1/ms)
    !  k_ws            = 0.012 (1/ms)
    !  r_w             = 0.5
    !  r_s             = 0.25
    !  y_s             = 0.0085
    !  y_w             = 0.615
    !  phi             = 2.23
    !  A_eff           = 25
    !  beta_0          = 2.3
    !  beta_1          = -2.4
    !  T_ref           = 120   (kPa)
    ! ----------------------------------------------------------------------

    ! Time increment
    dtime_aux = dtime

    ! Declare parameters of the model
    ! ----------------------------------------------------------------------
    k_TRPN = 100.0_rp
    n_TRPN = 2.0_rp
    Ca50_ref = 0.805_rp
    k_u = 1000.0_rp
    n_tm = 5.0_rp
    TRPN_50 = 0.35_rp
    k_uw = 182.0_rp
    k_ws = 12.0_rp
    r_w = 0.5_rp
    r_s = 0.25_rp
    y_s = 0.0085_rp
    y_w = 0.615_rp
    phi = 2.23_rp
    A_eff = 25.0_rp
    beta_0 = 2.3_rp
    beta_1 = -2.4_rp
    T_ref = 120.0_rp
    k_wu = k_uw*((1.0_rp/r_w)-1.0_rp)-k_ws
    k_su = k_ws*((1.0_rp/r_s)-1.0_rp)*r_w
    A_s = (A_eff*r_s)/((1.0_rp-r_s)*r_w+r_s)
    c_s = phi*k_ws*(1.0_rp-r_s)*r_w/r_s
    c_w = phi*k_uw*((1.0_rp-r_s)*(1.0_rp-r_w))/((1.0_rp-r_s)*r_w)
    k_u_2 = (k_u*(TRPN_50**n_tm))/(1.0_rp-r_s-(1.0_rp-r_s)*r_w)

    param(1) = k_TRPN
    param(2) = n_TRPN
    param(3) = Ca50_ref
    param(4) = k_u
    param(5) = n_tm
    param(6) = TRPN_50
    param(7) = k_uw
    param(8) = k_ws
    param(9) = r_w
    param(10) = r_s
    param(11) = y_s
    param(12) = y_w
    param(13) = phi
    param(14) = A_eff
    param(15) = beta_0
    param(16) = beta_1
    param(17) = T_ref*10000.0_rp ! Convert kPa to baryes
    param(18) = k_wu
    param(19) = k_su
    param(20) = A_s
    param(21) = c_s
    param(22) = c_w
    param(23) = k_u_2
    ! ----------------------------------------------------------------------

    ! Retrieve previous values of lambda for lambda and lambda rate calculations
    lambda_prev = exm_lambda_sld(ielem,igaus,1)

    ! Retrieve lambda
    lambda_aux = MIN(1.2_rp,lamda(1))
    !lambda_aux = lamda(1)

    ! Find the lambda rate through finite differences
    lambda_rate = (lamda(1)-lambda_prev)/dtime_aux

    ! Calculate Ca50
    Ca50 = Ca50_ref+beta_1*(lambda_aux-1.0_rp)
    !Ca50 = 5.4725_rp*EXP(-1.897_rp*lambda_aux)

    ! Calculate h
    h = MAX(0.0_rp,1.0_rp+beta_0*(lambda_aux+MIN(0.87_rp,lambda_aux)-1.87_rp))
    !if (lambda_aux < 0.6_rp)  then
    !  h = 0.0_rp
    !else
    !  h = 2.0_rp*(lambda_aux**2.0_rp)-2.4_rp*lambda_aux+0.72_rp
    !end if

    ! Loop for Newton-Raphson (Backward Euler)
    max_nr = 100
    Jac = 0.0_rp
    I_6 = 0.0_rp
    do icomp = 1,6
      I_6(icomp,icomp) = 1.0_rp
    end do
    res_tol = 0.0000000001_rp
    statvar(1:6,2) = statvar(1:6,1)

    do inewton = 1,max_nr

      ! State variables equal to the converged ones of the previous time step
      !S = MAX(0.0_rp,statvar(1,2))
      !W = MAX(0.0_rp,statvar(2,2))
      !CaTRPN = MAX(0.0_rp,statvar(3,2))
      !B = statvar(4,2)
      !zeta_s = statvar(5,2)
      !zeta_w = statvar(6,2)
      S = statvar(1,2)
      W = statvar(2,2)
      CaTRPN = statvar(3,2)
      B = statvar(4,2)
      zeta_s = statvar(5,2)
      zeta_w = statvar(6,2)

      sol(1) = S
      sol(2) = W
      sol(3) = CaTRPN
      sol(4) = B
      sol(5) = zeta_s
      sol(6) = zeta_w

      ! State variable dependent parameters (regularised parameters)
      if ((zeta_s+1.0_rp) < 0.0_rp) then
         y_su = y_s*(-zeta_s-1.0_rp)
      elseif ((zeta_s+1.0_rp) > 1.0_rp) then
         y_su = y_s*zeta_s
      else
         y_su = 0.0_rp
      end if
      !y_su = 7.887E-4*(zeta_s**2.0_rp)+6.3531E-4*zeta_s

      y_wu = y_w*ABS(zeta_w)
      !y_wu = 0.0577*(zeta_w**2.0_rp)

      ! Define the RHS of the ODE system
      ydot(1) = k_ws*W-k_su*S-y_su*S
      ydot(2) = k_uw*(1.0_rp-B-S-W)-k_wu*W-k_ws*W-y_wu*W
      ydot(3) = k_TRPN*(((gpcac/Ca50)**n_TRPN)*(1.0_rp-CaTRPN)-CaTRPN)
      ydot(4) = k_u_2*MIN((CaTRPN**(-n_tm/2.0_rp)),100.0_rp)*(1.0_rp-B-S-W)-k_u*                  &
       (CaTRPN**(n_tm/2.0_rp))*B
      ydot(5) = A_s*lambda_rate-c_s*zeta_s
      ydot(6) = A_s*lambda_rate-c_w*zeta_w

      ! Define the residuals
      res(1:6) = sol(1:6)-statvar(1:6,1)-dtime_aux*ydot(1:6)

      ! Define the Jacobian
      ! ======================================================================
      Jac = 0.0_rp

      Jac(1,1) = -k_su-y_su
      Jac(1,2) = k_ws

      if ((zeta_s+1.0_rp) < 0.0_rp) then
        Jac(1,5) = -y_s
      elseif ((zeta_s+1.0_rp) > 1.0_rp) then
        Jac(1,5) = y_s
      else
        Jac(1,5) = 0.0_rp
      end if
      !Jac(1,5) = 2.0_rp*7.887E-4*zeta_s+6.3531E-4
       
      Jac(2,1) = -k_uw
      Jac(2,2) = -k_uw-k_wu-k_ws-y_wu
      Jac(2,4) = -k_uw

      if (zeta_w < 0.0_rp) then
        Jac(2,6) = -y_w
      else
        Jac(2,6) = y_w
      end if
      !Jac(2,6) = 2.0_rp*0.0577*zeta_w

      Jac(3,3) = -k_TRPN*(((gpcac/Ca50)**n_TRPN)+1.0_rp)

      if ((CaTRPN**(-n_tm/2.0_rp)) <= 100.0_rp) then
        term = CaTRPN**(-n_tm/2.0_rp)
      else
        term = 100.0_rp
      end if
      Jac(4,1) = -k_u_2*term

      Jac(4,2) = Jac(4,1)

      if ((CaTRPN**(-n_tm/2.0_rp)) <= 100.0_rp) then
        Jac(4,3) = (-n_tm/2.0_rp)*(k_u_2*(1.0_rp-B-S-W)*(CaTRPN**(-(n_tm/2.0_rp)-1.0_rp))+        &
         k_u*B*(CaTRPN**((n_tm/2.0_rp)-1.0_rp)))
      else
        Jac(4,3) = (-n_tm/2.0_rp)*k_u*B*(CaTRPN**((n_tm/2.0_rp)-1.0_rp))
      end if

      Jac(4,4) = -k_u_2*term-k_u*(CaTRPN**(n_tm/2.0_rp))
      Jac(5,5) = -c_s
      Jac(6,6) = -c_w

      Jac = I_6-dtime_aux*Jac
      ! ======================================================================

      ! Compute the inverse of the Jacobian
      Jac_inv = Jac
      call invert(Jac_inv,6,6)

      ! Solution of the state variables
      statvar(1:6,2) = statvar(1:6,2)-MATMUL(Jac_inv(1:6,1:6),res(1:6))

      ! Check for convergence
      norm_res = 0.0_rp
      do icomp = 1,6
        norm_res = norm_res+(res(icomp)**2.0_rp)
      end do
      norm_res = SQRT(norm_res)

      if (norm_res < res_tol) then
        exit
      end if
    end do
    ! ----------------------------------------------------------------------  

    ! Calculate the active tension (tactf)
    ! ----------------------------------------------------------------------  
    tactf = h*(T_ref/r_s)*(statvar(1,2)*(statvar(5,2)+1.0_rp)+statvar(2,2)*statvar(6,2))

    ! Convert the active tension to CGS units (from kPa)
    tactf = tactf*10000.0_rp
    ! ----------------------------------------------------------------------  

  end if

  tacts = 0.0_rp
  tactn = 0.0_rp
  tacis = 0.0_rp  ! for transversally isotropic coupling

  if (kfl_coupt_sld(pmate) == 2) then
    tacis = tactf*trans_sld(1,pmate)
    tactf = tactf*(1.0_rp-trans_sld(1,pmate))
  end if

  tactv_sld= tactv_sld + tactf/real(pgaus,rp)

  ! Orthotropic activation (according to Rossi et al., 2014)
  if( kfl_fiber_sld == 2 .or. kfl_fiber_sld == 3) then
    ort_act_k1 = ortk1_sld(pmate)
    ort_act_k2 = ortk2_sld(pmate)

    param(24) = ort_act_k1
    param(25) = ort_act_k2

    tacts = ort_act_k1*tactf
    tactn = ort_act_k2*tactf
  end if

  ! Pre-stress
  !prestress = 0.0_rp
  !if (cutim <= preti_sld(pmate)) then
  !  prestress = prest_sld(pmate) * (cutim / preti_sld(pmate))
  !else
  !  prestress = prest_sld(pmate)
  !end if
  ! tactf = tactf + prest_sld(1)

  !f_true x f_true (outer products [3x3]) 
  fporf_t = 0.0_rp
  spors_t = 0.0_rp
  nporn_t = 0.0_rp    
  nfibe_length = 0.0_rp
  do idime = 1,ndime
     nfibe_length = nfibe_length+nfibt(idime,igaus)*nfibt(idime,igaus) 
     do jdime = 1,ndime
        fporf_t(idime,jdime) = nfibt(idime,igaus)*nfibt(jdime,igaus) 
        spors_t(idime,jdime) = nshtt(idime,igaus)*nshtt(jdime,igaus) 
        nporn_t(idime,jdime) = normt(idime,igaus)*normt(jdime,igaus)
    end do
  end do
  nfibe_length = sqrt(nfibe_length) 

  ! Cauchy active stress tensor (castr_active)
  do idime = 1,ndime
    do jdime = 1,ndime          
      castr_active(idime,jdime) = (1.0_rp/gpdet(igaus))*((fporf_t(idime,jdime)*tactf)+            &
       (spors_t(idime,jdime)*tacts)+(nporn_t(idime,jdime)*tactn))
    end do
  end do

  if (kfl_coupt_sld(pmate) == 2) then
    ! add isotropy in the transversally isotropic coupling
    do jdime = 1,ndime          
      castr_active(idime,idime) = castr_active(idime,idime)+tacis
    end do
  end if

101 format (9(e13.6,' ')) 

end subroutine sld_eccoup
