!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_myofil.f90
!> @author  
!> @date    19/12/2014
!> @brief   Computing the active Force. Rice myofilament model. 
!> @details Active force computed in the myofilament model of Rice. No internal cell shortening Isosarcometric case.
!> @} 
!-----------------------------------------------------------------------
subroutine sld_myofil(&
     pgaus,pmate,igaus,gpdet,gptlo,lamda,nfibt,nshtt,normt,nfibe_length,&
     gpigd,castr_active,gpstr_active,gpcac,tactf,ielem)

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_solidz
  use def_master, only       :  dtime,cutim,kfl_cellmod,&
       CELL_TT2006_EXMEDI,CELL_OHARA_EXMEDI,CELL_FITZHUGH_EXMEDI,CELL_NOMOD_EXMEDI,ittim,dtinv

  implicit none

  integer(ip), intent(in)    :: pgaus,pmate,igaus,ielem  
  real(rp),    intent(in)    :: gptlo(pgaus),lamda(3),gpigd(ndime,ndime,pgaus),gpdet(pgaus),&
       nfibt(ndime,pgaus),nshtt(ndime,pgaus),normt(ndime,pgaus),gpcac
  real(rp),    intent(out)   :: castr_active(ndime,ndime),gpstr_active(ndime,ndime),tactf
  real(rp)                   :: ca0,cam,thau,cca2p,n,beta,Cal50,tmaxi,tacts,tactn,tacis
  real(rp)                   :: fporf_t(ndime,ndime),spors_t(ndime,ndime),nporn_t(ndime,ndime)
  real(rp)                   :: Factive, koffmod, nfibe_length
  integer(ip)                :: idime,jdime,kdime,ldime
  real(rp)                   :: Force, Fpreload, Fafterload, Fpassive, Fnordv, dintf
  real(rp)                   :: SLmax, SLmin, lenthin, lenthick, lenhbare, SL0, lencoeff1, lencoeff2
  real(rp)                   :: Qkon, Qkoff, Qknp, Qkpn, Qfapp, Qgapp
  real(rp)                   :: Qhf, Qhb, Qgxb, QT, TempK
  real(rp)                   :: kon, koffL, koffH, knp, kpn
  real(rp)                   :: konT, koffLT, koffHT, knpT, kpnT, kxb
  real(rp)                   :: permtot, invpermtot, Tropreg
  real(rp)                   :: fapp, gapp, gslmod, hf, hfmdc, hb, gxb, xbmodsp, perm50
  real(rp)                   :: gappslmod, hfmod, hbmod, gxbmod, hbmdc
  real(rp)                   :: fappT, gappT, hfT, hbT, gxbT
  real(rp)                   :: sigmap, sigman
  real(rp)                   :: x0, phi
  real(rp)                   :: SLrest, PConti, PExpti, PConcol, PExpCol, Ftitin, Fcol
  real(rp)                   :: SLcol, mass, visc, KSE
  real(rp)                   :: sovrze, sovrcle, lensovr, SOVFthin, SOVFthick
  real(rp)                   :: XBdutyprer, XBdutypostr, XBMPrer, XBMPostr
  real(rp)                   :: t, h, N0XB 
  real(rp)                   :: nperm
  integer(ip)                :: j, k, i!, l 


  castr_active=0.0_rp
  gpstr_active=0.0_rp
  ca0  = 0.0_rp!0.01 !mmole
  cam  = 1.0_rp !mmole
  !Default value: 0.06_rp sec (human) ... !0.01268 sec Rat
  thau = timec_sld(pmate) 


  if ((kfl_cellmod(pmate) == CELL_TT2006_EXMEDI).or. &
       (kfl_cellmod(pmate) == CELL_OHARA_EXMEDI)) then   !! Ten Tuscher calcula el Ca2+ dins del propi model
     !
     ! Cell model
     !
     cca2p=gpcac  !
     !C50
     Cal50= cal50_sld(pmate)!0.6 !value by Hung 
     ! Cal50= 4.35_rp/sqrt(exp(4.75_rp*((1.9*lamda(1))-1.58_rp))-1.0_rp) 
  else if (kfl_cellmod(pmate) == CELL_FITZHUGH_EXMEDI)  then  !! FHN
     cca2p= ca0 + (cam-ca0)*(gptlo(igaus)/thau)*exp(1.0_rp-(gptlo(igaus)/thau))
  else if (kfl_cellmod(pmate) == CELL_NOMOD_EXMEDI)  then  !
     !
  else 
     call runend ("SLD_MYOFIL: CELL MODEL NOT INCLUDED IN THE COUPLING WITH SOLIDZ MODULE")
  end if

  !--------------------------------------------------------------
  !             Initial Conditions are HARD CODED IN SLD_INIVAR   
  !--------------------------------------------------------------

  !iy(1) = 0.0_rp  !xXBprer
  !iy(2) = 0.007_rp  * 0.0001_rp !xXBPostr (in cm)
  !iy(3) = 1.0_rp  !N
  !iy(4) = 2.0579e-7_rp  !XBprer
  !iy(5) = 8.7571e-7_rp  !XBpostr
  !iy(6) = (1.0_rp-1.8e-7_rp) !N_NoXB - Outside of single overlap
  !iy(7) = 1.8e-7_rp  !P_NoXB - Outside of single overlap
  !iy(8) = 0.0196_rp  !TRPNCaL
  !iy(9) = 0.1667_rp  !TRPNCaH
  !iy(10) =  1.85_rp* 0.0001_rp      ! SL
  !iy(11) = 0.0_rp !intF  integral of force
  !iy(12)= 0.0_rp ! P  probability 


  !--------------------------------------------------------------
  !            CONSTANTS                                         
  !--------------------------------------------------------------

  kon          = 50.0_rp      !(uM^-1 s^-1)
  koffL        = 250.0_rp     !(s^-1)
  koffH        = 25.0_rp      !(s^-1)
  perm50       = 0.6738_rp
  !     perm50       = 0.55
  knp          = 500.0_rp      !(1/s) In the paper: 50 Unitless
  kpn          = 50.0_rp     !(1/s)  In the paper: 500 Unitless


  fapp         = 500.0_rp    !(1/s)
  gapp         = 70.0_rp     !(1/s)
  gslmod       = 6.0_rp      
  hf           = 2000.0_rp   !(1/s)
  hfmdc        = 5.0_rp
  hb           = 400.0_rp    !(1/s)
  hbmdc        = 0.0_rp
  gxb          = 70.0_rp     !(1/s)
  kxb          = 1200000.0_rp    !(Dynes/cm^2)     
  sigmap       = 8.0_rp
  sigman       = 1.0_rp
  xbmodsp      = 0.2_rp
  koffmod      = 1.0_rp        ! mod to change species 

  !! Temperature Dependence (degrees Kelvin)
  TempK        = 310.0_rp  !273 + 37
  QT           = (TempK-310_rp)/10_rp
  Qkon         = 1.5_rp**QT
  Qkoff        = 1.3_rp**QT
  Qknp         = 1.6_rp**QT
  Qkpn         = 1.6_rp**QT
  Qfapp        = 6.25_rp**QT
  Qgapp        = 2.5_rp**QT
  Qhf          = 6.25_rp**QT
  Qhb          = 6.25_rp**QT
  Qgxb         = 6.25_rp**QT


  x0           = 0.007_rp * 0.0001_rp    !(cm)
  phi          = 2.0_rp  !xPsi


  SLmax         = 2.4_rp * 0.0001_rp     !(cm)
  SLmin         = 1.4_rp * 0.0001_rp     !(cm)
  lenthick      = 1.65_rp * 0.0001_rp    !(cm)
  lenthin       = 1.2_rp * 0.0001_rp     !(cm)
  lenhbare      = 0.1_rp * 0.0001_rp     !(cm)
  SL0           = 1.85_rp * 0.0001_rp     !(cm)


  SLrest        = 1.85_rp * 0.0001_rp     !(cm)
  SLcol         = 2.25_rp * 0.0001_rp     !(cm)
  PConti        = 0.002_rp     !(Unit normalized force)
  PExpti        = 10.0_rp      
  PConcol       = 0.02_rp      !(Unit normalized force)
  PExpcol       = 70.0_rp
  nperm         = 8.1997_rp
  !    nperm         = 7.5

  mass          = 0.2_rp !([norm Force s^2]/cm)
  visc          = 30.0_rp   !(norm Force s/cm)
  KSE           = 10000.0_rp     !(Unit normalised force)/cm 

  !------------------------------------------------------------
  !         PARAMETERS THAT DEPEND ON OTHER CONSTANTS  
  !------------------------------------------------------------


  !--------------------------------------------------------------
  !        SOLVING THE SYSTEM OF ODE's. EULER.
  !--------------------------------------------------------------

  !y = iy  myofi_sld(9,3)
  t = cutim
  h=  dtime


  sovrze       = min(0.5_rp*lenthick,0.5_rp*myofi_sld(10,1))
  sovrcle      = max(myofi_sld(10,1)/2_rp-(myofi_sld(10,1)-lenthin),lenhbare*0.5_rp)
  lensovr      = sovrze-sovrcle
  SOVFthick    = (2.0_rp*lensovr)/(lenthick-lenhbare)
  SOVFthin     = lensovr/lenthin

  ! Compute combined Ca binding to high-(w/XB) and low-(no XB) sites
  Tropreg      = (1.0_rp-SOVFthin)*myofi_sld(8,1) + SOVFthin*myofi_sld(9,1)
  permtot      = sqrt(1.0_rp/(1.0_rp+(perm50/Tropreg)**nperm))
  invpermtot   = min(1.0_rp/permtot,100.0_rp)

  !Adjustments for Ca activation, temperature, SL, stress and strain
  konT         = kon*Qkon
  koffLT       = koffL*Qkoff*koffmod
  koffHT       = koffH*Qkoff*koffmod
  knpT         = knp*permtot*Qknp
  kpnT         = kpn*invpermtot*Qkpn
  fappT        = fapp*xbmodsp*Qfapp     
  gappslmod    = 1_rp + (1_rp-SOVFthick)*gslmod
  gappT        = gapp*gappslmod*xbmodsp*Qgapp
  hfmod        = exp(SIGN(myofi_sld(1,1),-1.0_rp)*hfmdc*(myofi_sld(1,1)/x0)*(myofi_sld(1,1)/x0))
  !hfmod        = exp((SIGN(0.5_rp,myofi_sld(1,1))-SIGN(0.5_rp,-myofi_sld(1,1)))*hfmdc*) 
  hbmod        = exp(SIGN((myofi_sld(1,2)-x0),1.0_rp)*hbmdc*(myofi_sld(2,1)/x0)*(myofi_sld(2,1)/x0)) 
  hfT          = hf*hfmod*xbmodsp*Qhf
  hbT          = hb*hbmod*xbmodsp*Qhb
  if ( (x0-myofi_sld(2,1)) >= 0.0_rp )   then
     gxbmod        = exp(sigman*((x0-myofi_sld(2,1))/x0)*((x0-myofi_sld(2,1))/x0))
  else 
     gxbmod        = exp(sigmap*((x0-myofi_sld(2,1))/x0)*((x0-myofi_sld(2,1))/x0))
  end if
  gxbT         = gxb*gxbmod*xbmodsp*Qgxb

  !Regulation and crossbridge cycling state derivatives
  myofi_sld(12,1) = 1_rp - myofi_sld(3,1) - myofi_sld(4,1) - myofi_sld(5,1) !P

  myofi_sld(8,3) =cca2p*konT*(1-myofi_sld(8,1))- koffLT*myofi_sld(8,1)  !dTRPNCaL
  myofi_sld(9,3) =cca2p*konT*(1-myofi_sld(9,1))- koffHT*myofi_sld(9,1)  !dTRPNCaH

  myofi_sld(6,3) = -knpT*myofi_sld(6,1)+kpnT*myofi_sld(7,1)          !dN_NoXB - Outside of single overlap
  myofi_sld(7,3) = knpT*myofi_sld(6,1)-kpnT*myofi_sld(7,1)           !dP_NoXB - Outside of single overlap
  myofi_sld(3,3) = -knpT*myofi_sld(3,1)+ kpnT * myofi_sld(12,1)       !dN

  myofi_sld(4,3) = fappT*myofi_sld(12,1) - gappT*myofi_sld(4,1)- hfT*myofi_sld(4,1) + hbT*myofi_sld(5,1)   !dXBprer
  myofi_sld(5,3) = hfT*myofi_sld(4,1)-(hbT+gxbT)*myofi_sld(5,1)     !dXBpostr
  myofi_sld(12,3) = -(myofi_sld(3,3)+myofi_sld(4,3)+myofi_sld(5,3)) !dP

  ! steady-state fractions in XBprer and XBpostr using King-Altman rule 
  XBMPrer      = (fapp*hb+fapp*gxb)/(gxb*hf+fapp*hf+gapp*hb+gapp*gxb+fapp*hb+fapp*gxb)  !SSXBprer
  XBMPostr     = fapp*hf/(gxb*hf+gapp*hb+fapp*hf+gapp*gxb+fapp*hb+fapp*gxb)  !SSXBpostr
  ! % normalization for scaling active and passive force (maximal force)    
  Fnordv = kxb*x0*XBMPostr!

  ! Calculate Forces (active, passive, preload, afterload)
  Force   = (kxb*SOVFthick*(myofi_sld(1,1)*myofi_sld(4,1)+myofi_sld(2,1)*myofi_sld(5,1)))     
  Factive = Force/Fnordv     
  if (myofi_sld(10,1) >= SLrest) then 
     Ftitin       = PConti*(exp(PExpti*(myofi_sld(10,1)-SLrest))-1.0_rp)
  else if (myofi_sld(10,1) < SLrest) then
     Ftitin       = -PConti*(exp(PExpti*(SLrest-myofi_sld(10,1)))-1.0_rp)
  end if

  if (myofi_sld(10,1) >= SLcol) then
     Fcol = PConcol*(exp(PExpcol*(myofi_sld(10,1)-SLcol))-1)
  else if (myofi_sld(10,1) < SLcol) then
     Fcol         = 0.0_rp
  end if

  Fpassive     = Ftitin + Fcol
  Fpreload     = 0.0_rp   ! Applied force that could induce a longer sarcomere length than resting length 
  !preload = sign(SLset-SLrest)*PCon_t*(exp(PExp_t*abs(SLset-SLrest))-1);    
  Fafterload   = 0.0_rp !KSE*(myofi_sld(10,1)-SLrest)  
  !print*, 'Fafterload', Fafterload, 'Fpassive', Fpassive
  myofi_sld(11,3) = -Fpassive+Fpreload-Factive+Fafterload     !dintf
  !change in SL    
  lencoeff1 = 1.0_rp
  if (myofi_sld(10,1) < SLmin) lencoeff1=0
  lencoeff2 = 1.0_rp
  if (SLmax < myofi_sld(10,1)) lencoeff2=0     
  myofi_sld(10,3) = lencoeff1 * lencoeff2  * (myofi_sld(11,1)+(SLrest-(myofi_sld(10,1)* nfibe_length))*visc)/mass  !dSL   

  ! Mean strain of strongly-bound states due to SL motion and XB cycling     
  XBdutyprer   = (fappT*hbT+fappT*gxbT)/(gxbT*hfT+fappT*hfT+gappT*hbT+gappT*gxbT+fappT*hbT+fappT*gxbT)
  XBdutypostr  =(fappT*hfT)/(gxbT*hfT+fappT*hfT+gappT*hbT+gappT*gxbT+fappT*hbT+fappT*gxbT)

  myofi_sld(1,3) = myofi_sld(10,3)*0.5_rp + phi/XBdutyprer*(-myofi_sld(1,1)*fappT+(myofi_sld(2,1)-x0-myofi_sld(1,1))*hbT) !dxXBprer
  myofi_sld(2,3) = myofi_sld(10,3)*0.5_rp + phi/XBdutypostr*(x0+myofi_sld(1,1)-myofi_sld(2,1))*hfT  !dxXBPostr

  !!N0XB    = 1_rp - myofi_sld(7,1) dunno what this is

  ! Ca buffering by low-affinity troponin C (LTRPNCa)
  !FrSBXB    = (XBpostr+XBprer)/(SSXBpostr + SSXBprer);
  !dFrSBXB   = (dXBpostr+dXBprer)/(SSXBpostr + SSXBprer);

  !dsovr_ze  = -dSL/2*heav(len_thick-SL);
  !dsovr_cle = -dSL/2*heav((2*len_thin-SL)-len_hbare);
  !dlen_sovr = dsovr_ze-dsovr_cle;
  !dSOVFThin = dlen_sovr/len_thin;
  !dSOVFThick= 2*dlen_sovr/(len_thick-len_hbare);

  !TropTot = Trop_conc*((1-SOVFThin)*TRPNCaL + ...
  !  SOVFThin*(FrSBXB*TRPNCaH+(1-FrSBXB)*TRPNCaL));
  !dTropTot= Trop_conc*(-dSOVFThin*TRPNCaL+(1-SOVFThin)*dTRPNCaL + ...
  !  dSOVFThin*(FrSBXB*TRPNCaH+(1-FrSBXB)*TRPNCaL) + ...
  !  SOVFThin*(dFrSBXB*TRPNCaH+FrSBXB*dTRPNCaH-dFrSBXB*TRPNCaL+...
  !  (1-FrSBXB)*dTRPNCaL));

  !dforce = kxb*dSOVFThick*(xXBpostr*XBpostr+xXBprer*XBprer) + ...
  !  kxb*SOVFThick*(dxXBpostr*XBpostr+xXBpostr*dXBpostr + ...
  !  dxXBprer*XBprer+xXBprer*dXBprer);

  do j = 1, 12         
     myofi_sld(j,2) = myofi_sld(j,1) + h*myofi_sld(j,3)         
  end do

  myofi_sld(:,1) = myofi_sld(:,2)
  t = t + h

  tactf = -myofi_sld(11,1) * cocof_sld(pmate) * 1.0e9_rp !example to scale

  !       t = t + h(l)
  tacts = 0.0_rp
  tactn = 0.0_rp
  tacis = 0.0_rp  ! for transversally isotropic coupling

  !if (kfl_coupt_sld(pmate) == 2) then
  !   tacis = tactf * trans_sld(1,pmate)
  !   tactf = tactf * (1.0_rp - trans_sld(1,pmate))
  !end if


  tactv_sld= tactv_sld + tactf/real(pgaus,rp)


  !  if (ielem == 1253) write(6,*) 'poto',ielem,nfibe_length,tactf
  !  if (ielem ==    1) write(6,*) 'poto',ielem,nfibe_length,tactf

  ! Cauchy active stress tensor (castr_active)
  do idime=1,ndime
     do jdime=1,ndime          
        castr_active(idime,jdime) = (fporf_t(idime,jdime)*tactf) +&
             (spors_t(idime,jdime)*tacts) +&
             (nporn_t(idime,jdime)*tactn)
     end do
  end do

  !     write(100, *) t, y(8)!, h 
  !     write(200, *) t, y(9)!, h 
  !     write(300, *) t, y(1)!, h 
  !     write(400, *) t, y(2)!, h
  !     write(500, *) t, y(3)!, h 
  !     write(600, *) t, y(4)!, h
  !     write(700, *) t, y(5)!, h 
  !     write(800, *) t, y(6)!, h 
  !     write(900, *) t, y(7)!, h
  !     write(1000,*) t, P
  !     write(2000,*) t, Factive
  !    write(4000,*) t, myofi_sld(10,1)
  !    write(6000,*) t, tactf

  !  end do





end subroutine sld_myofil
