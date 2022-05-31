!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_oceohr.f90
!> @author  Jazmin Aguado-Sierra
!> @brief   Single cell run for Initial condition setup for Ohara-Rudy 2011 heterogeneous model
!> @date   16/NOV/1966
!> @details Runs a single cell simulation at th given frequency and pathologic conditions \n
!!    It performs single cell runs under normal, heart failure or drugs \n
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_oharaf_ionicurrents(ipoin,xioni,dioni,cai)

  use      def_parame
  use      def_master
  use      def_elmtyp
  use      def_domain
  use      def_exmedi
  use      mod_exm_oharaequations

  ! definition of variables
  implicit none
  integer(ip), intent(in) :: ipoin !< node
  real(rp), intent(out) :: xioni   !< current
  real(rp), intent(out)   :: dioni !< current derivative
  real(rp), intent(out)   :: cai !< current calcium
  integer(ip) :: i,n,imate,ituss_exm
  real(rp)    :: vinf, xitaux, vaux0, vaux1, vaux2, vaux3, a2bas
  real(rp)    ::  vffrt, vfrt, ena, ek, eks, bt, a_rel, btp, a_relp, ccm, Inet
  real(rp)    :: pkna, farad,  v1, v2
  real(rp)    ::  a1, a2, a3, a4, b1, b2, b3, b4, k1p, k2p, k3p, k4p, k1m, k2m, k3m, k4m
  real(rp)    :: h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, x1, x2, x3, x4
  real(rp)    :: k1, k2, k3, k4, k5, k6, k7, k8, e1, e2, e3, e4, k3pp, k4pp
  real(rp)    :: knai,knao
  real(rp)    :: gna, gnal, gto, gks, gk1, gncx, gkb, gkr, gpca, phical, phicana, phicak, nao
  real(rp)    :: ahf, ahs, hp, hh, finalp, ii, aff, afs, afcaf, afcas, fcap
  real(rp)    :: zca, pca, pcap, pcak, pcanap, pcakp, ksca, kasymm
  real(rp)    :: wna, wca, wnaca, kcaon, kcaoff, qna, qca, hca, hna, kmcaact, allo, zna, jncxna, jncxca
  real(rp)    :: delta, kki, kko, mgadp, mgatp, kmgatp, ep, khp, knap, kxkur, pp, zk
  real(rp)    :: jnakna, jnakk, pnak, xkb, pnab, ipp, aif, ais
  real(rp)    :: fp, ko, axrf, axrs, xr, rkr, rk1, fjupp, pcab, pcana
  real(rp)    :: cao, ff, fca, ficalp, finap,  fjrelp
  real(rp)    :: rgas, temp, jupnp, jupp, fitop, hbig, camko, kmcam
  real(rp)    :: cass, cansr, cajsr, nass, nai, ki, kss, camkt
  real(rp)    :: rhsx1, rhsx2, rhsx, val0, bslmax, bsrmax,cmdnmax, csqnmax
  real(rp)    :: leng, rad, vcell, ageo, acap, vmyo, vnsr, vjsr, vss, kmcamk, kmcmdn
  real(rp)    :: acamk, bcamk, trpnmax, kmcsqn, kmtrpn, kmbsl, kmbsr, INaCa
  real(rp)    :: dosis1, dosis2, dosis3, i50cal1, i50na1, i50kr1, gcaldrug1, gkrdrug1, gnadrug1, dtimeEP
  real(rp)    :: dosis4, i50k11, gk1drug1,dosis5, dosis6, dosis7, dosis8, i50cal2, i50na2, i50kr2
  real(rp)    :: i50k12, gk1drug2, gcaldrug2, gkrdrug2, gnadrug2, dosis9, dosis10, i50nal1, i50nal2
  real(rp)    :: gnaldrug1, gnaldrug2, dosis11, dosis12, gksdrug1, gksdrug2, i50ks1, i50ks2   
  real(rp)    :: sac
  real(rp)    :: inf_results(EXM_OHR_INF_EQUATIONS), tau_results(EXM_OHR_TAU_EQUATIONS)
  real(rp), parameter :: KNA1 = 15.0_rp
  real(rp), parameter :: KNA2 = 5.0_rp
  real(rp), parameter :: KNA3 = 88.12_rp
  real(rp), parameter :: INV_KNA1 = 1.0_rp / KNA1
  real(rp), parameter :: INV_KNA2 = 1.0_rp / KNA2
  real(rp), parameter :: INV_KNA3 = 1.0_rp / KNA3
  real(rp), parameter :: KNAI0 = 9.073_rp
  real(rp), parameter :: KNAO0 = 27.78_rp
  logical :: flag_land


  gna  = 0.0_rp
  gnal = 0.0_rp
  gk1  = 0.0_rp
  gna  = 0.0_rp
  gkr  = 0.0_rp
  pca  = 0.0_rp


  if(INOTMASTER) then

     sac=0.0_rp !initialize sac current
     flag_land = .false.

     dtimeEP=dtime * 1000.0_rp
     if (kfl_atbhe_exm == 0_ip) then
        a2bas = 1.0_rp
     else 
        a2bas = atbhe_exm(1,ipoin)
     end if

     
     if (kfl_cellmod(nodemat(ipoin)) == 5_ip) then

        ituss_exm = int(celty_exm(1,ipoin))  

        n = nodemat(ipoin)
  
        if (n==0) then
            call runend('Elements with material 0 are present')
        end if

        do imate= 1,nmate_exm
           if (kfl_eccty(imate) == 4_ip) flag_land = .TRUE.
        end do
        
!!!!  DRUG DEFINITION TO IKR, INA OR ICAL
  !!!!  DRUG DEFINITION TO IKR, INA OR ICAL
        if(kfl_drugsmate_exm(n)==1_ip)then
           dosis1 = drugdmate_exm(1,n)
           i50cal1 = drugdmate_exm(2,n)
           dosis2 = drugdmate_exm(3,n)
           i50kr1 = drugdmate_exm(4,n)
           dosis3 = drugdmate_exm(5,n)
           i50na1 = drugdmate_exm(6,n)
           dosis4 = drugdmate_exm(7,n)
           i50k11 = drugdmate_exm(8,n)
           dosis5 = drugdmate_exm(9,n)
           i50nal1 = drugdmate_exm(10,n)
           dosis6 = drugdmate_exm(11,n)
           i50ks1 = drugdmate_exm(12,n)
           !a second drug dose that may affect the same channels    
           dosis7 = drugdmate_exm(13,n)
           i50cal2 = drugdmate_exm(14,n)
           dosis8 = drugdmate_exm(15,n)
           i50kr2 = drugdmate_exm(16,n)   
           dosis9 = drugdmate_exm(17,n)
           i50na2 = drugdmate_exm(18,n)   
           dosis10 = drugdmate_exm(19,n)
           i50k12 = drugdmate_exm(20,n)
           dosis11 = drugdmate_exm(21,n)
           i50nal2 = drugdmate_exm(22,n)   
           dosis12 = drugdmate_exm(23,n)
           i50ks2 = drugdmate_exm(24,n)
        
           gcaldrug1 = 1.0_rp / (1.0_rp + (dosis1/i50cal1))
           gkrdrug1 = 1.0_rp / (1.0_rp + (dosis2/i50kr1))
           gnadrug1 = 1.0_rp / (1.0_rp + (dosis3/i50na1))
           gk1drug1 = 1.0_rp / (1.0_rp + (dosis4/i50k11))
           gnaldrug1 = 1.0_rp / (1.0_rp + (dosis5/i50nal1))
           gksdrug1 = 1.0_rp / (1.0_rp + (dosis6/i50ks1))
        
           gcaldrug2 = 1.0_rp / (1.0_rp + (dosis7/i50cal2))
           gkrdrug2 = 1.0_rp / (1.0_rp + (dosis8/i50kr2))
           gnadrug2 = 1.0_rp / (1.0_rp + (dosis9/i50na2))
           gk1drug2 = 1.0_rp / (1.0_rp + (dosis10/i50k12))
           gnaldrug2 = 1.0_rp / (1.0_rp + (dosis11/i50nal2))  
           gksdrug2 = 1.0_rp / (1.0_rp + (dosis12/i50ks2))  

           if ((dosis1  > 0.00001_rp) .and. (dosis7 > 0.00001_rp))  then 
           pca = 0.0001_rp * (1.0_rp - ((1.0_rp - gcaldrug1) + (1.0_rp - gcaldrug2))) 
           else if ((dosis1 > 0.00001_rp) .and. (dosis7 <= 0.00001_rp)) then
             pca = 0.0001_rp * gcaldrug1
           else if ((dosis7 > 0.00001_rp) .and. (dosis1 <= 0.00001_rp)) then
             pca = 0.0001_rp * gcaldrug2
           else if ((dosis1 <= 0.00001_rp) .and. (dosis7 <= 0.00001_rp)) then
             pca = 0.0001_rp
           end if   
           if ((dosis2 > 0.00001_rp) .and. (dosis8 > 0.00001_rp))  then 
           gkr = 0.046_rp * (1.0_rp - ((1.0_rp - gkrdrug1) + (1.0_rp - gkrdrug2))) 
           else if ((dosis2  > 0.00001_rp) .and. (dosis8 <= 0.00001_rp))  then
             gkr = 0.046_rp * gkrdrug1
           else if ((dosis8  > 0.00001_rp) .and. (dosis2 <= 0.00001_rp)) then
             gkr = 0.046_rp * gkrdrug2
           else if ((dosis2 <= 0.00001_rp) .and. (dosis8 <= 0.00001_rp)) then
             gkr = 0.046_rp
           end if       
           if ((dosis3  > 0.00001_rp) .and. (dosis9 > 0.00001_rp))  then 
            gna = 75.0_rp * (1.0_rp - ((1.0_rp - gnadrug1) + (1.0_rp - gnadrug2))) 
            else if ((dosis3 > 0.00001_rp) .and. (dosis9 <= 0.00001_rp)) then
              gna = 75.0_rp * gnadrug1 
            else if ((dosis9 > 0.00001_rp)  .and. (dosis3 <= 0.00001_rp)) then
              gna = 75.0_rp * gnadrug2 
            else if ((dosis3 <= 0.00001_rp) .and. (dosis9 <= 0.00001_rp)) then
              gna = 75.0_rp  
           end if
           if ((dosis4 > 0.00001_rp) .and. (dosis10 > 0.00001_rp))  then 
           gk1 = 0.1908_rp * (1.0_rp - ((1.0_rp - gk1drug1) + (1.0_rp - gk1drug2)))  
           else if ((dosis4 > 0.00001_rp) .and. (dosis10 <= 0.00001_rp)) then
             gk1 = 0.1908_rp * gk1drug1
           else if ((dosis10 > 0.00001_rp) .and. (dosis4 <= 0.00001_rp)) then
             gk1 = 0.1908_rp * gk1drug2
           else if ((dosis4 <= 0.00001_rp) .and. (dosis10 <= 0.00001_rp)) then
             gk1 = 0.1908_rp
           end if 
           if ((dosis5 > 0.00001_rp) .and. (dosis11 > 0.00001_rp))  then 
           gnal = 0.0075_rp * (1.0_rp - ((1.0_rp - gnaldrug1) + (1.0_rp - gnaldrug2)))  
           else if ((dosis5 > 0.00001_rp) .and. (dosis11 <= 0.00001_rp)) then
             gnal = 0.0075_rp * gnaldrug1
           else if ((dosis11 > 0.00001_rp) .and. (dosis5 <= 0.00001_rp)) then
             gnal = 0.0075_rp * gnaldrug2
           else if ((dosis5 <= 0.00001_rp) .and. (dosis11 <= 0.00001_rp)) then    
             gnal = 0.0075_rp
           end if 
           if ((dosis6 > 0.00001_rp) .and. (dosis12 > 0.00001_rp))  then 
           gks = 0.0034_rp * (1.0_rp - ((1.0_rp - gksdrug1) + (1.0_rp - gksdrug2))) 
           else if ((dosis6  > 0.00001_rp) .and. (dosis12 <= 0.00001_rp))  then
             gks = 0.0034_rp * gksdrug1
           else if ((dosis6  > 0.00001_rp) .and. (dosis12 <= 0.00001_rp)) then
             gks = 0.0034_rp * gksdrug2
           else if ((dosis6 <= 0.00001_rp) .and. (dosis12 <= 0.00001_rp)) then   
             gks = 0.0034_rp
           end if      
        else
           gna = 75.0_rp
           pca = 0.0001_rp 
           gkr = 0.046_rp
           gk1 = 0.1908_rp
           gnal = 0.0075_rp
           gks = 0.0034_rp
        endif 
        !extracellular ionic concentrations
        nao = 140.0_rp * ttparmate_exm(1,13,n)
        cao = 1.8_rp * ttparmate_exm(2,13,n)
        ko = 5.4_rp * ttparmate_exm(3,13,n)

        !physical constants
        rgas = 8314.0_rp
        temp = 310.0_rp
        farad = 96485.0_rp 
        !cell geometry
        leng = 0.01_rp
        rad = 0.0011_rp
        vcell = 1000.0_rp*3.14_rp*rad*rad*leng
        ageo = 2.0_rp*3.14_rp*rad*rad + 2.0_rp*3.14_rp*rad*leng
        acap = 2.0_rp*ageo
        vmyo = 0.68_rp*vcell
        vnsr = 0.0552_rp*vcell
        vjsr = 0.0048_rp*vcell
        vss = 0.02_rp*vcell
        !%camk constants
        kmcamk = 0.15_rp
        acamk = 0.05_rp
        bcamk = 0.00068_rp
        camko = 0.05_rp
        kmcam = 0.0015_rp
        !Current constants:
        pkna = 0.01833_rp
!!!!!!   calculate ina
        ahf = 0.99_rp
        ahs = 1.0_rp - ahf
        !gnal = 0.0075_rp
        gto = 0.02_rp !* gitodrug !this is no longer available for now
!!! calculate ical
!!!!!! calculate ff
        aff = 0.6_rp
        afs = 1.0_rp - aff  
!!!!! calculate icana and icak
        zca = 2.0_rp
        !!  inal  
!!!  calculate IKr
!!!  calculate IKs
        gks = gks * a2bas
!!!! calculate IK1  

        gncx = 0.0008_rp
        pnak = 30.0_rp
        !!calculate IKb
        gkb = 0.003_rp
        !%calcium buffer constants
        cmdnmax = 0.05_rp

        if(ituss_exm == EXM_CELLTYPE_EPI) then !!epi
           gnal = gnal * 0.6_rp*ttparmate_exm(3,6,n)
           gto = gto*4.0_rp*ttparmate_exm(3,1,n)
           pca = pca*1.2_rp*ttparmate_exm(3,9,n)  !epi
           gkr = gkr * 1.3_rp*ttparmate_exm(3,4,n)
           gk1 = gk1 * 1.2_rp*ttparmate_exm(3,3,n)
           gncx = gncx*1.1_rp*ttparmate_exm(3,7,n)
           pnak = pnak*0.9_rp*ttparmate_exm(3,10,n)
           gkb = gkb*0.6_rp*ttparmate_exm(3,8,n)
           cmdnmax = cmdnmax*1.3_rp*ttparmate_exm(3,11,n) 
           gks = gks*1.4_rp*ttparmate_exm(3,2,n)
           gna = gna*ttparmate_exm(3,5,n)
        else if(ituss_exm == EXM_CELLTYPE_MID) then !!MID
           gto = gto*4.0_rp*ttparmate_exm(2,1,n)
           pca = pca*2.5_rp*ttparmate_exm(2,9,n)  !mid 
           gkr = gkr * 0.8_rp*ttparmate_exm(2,4,n)
           gks = gks*ttparmate_exm(2,2,n)     
           gk1 = gk1 * 1.3_rp*ttparmate_exm(2,3,n)
           gncx = gncx*1.4_rp*ttparmate_exm(2,7,n)
           pnak = pnak*0.7_rp*ttparmate_exm(2,10,n)
           gna = gna*ttparmate_exm(2,5,n)
           gnal = gnal*ttparmate_exm(2,6,n)
           gkb = gkb*ttparmate_exm(2,8,n)
           cmdnmax = cmdnmax*ttparmate_exm(2,11,n)    
           if(kfl_inaga_exm(n) == 1_ip) then
              pca=(pca*1.8_rp)/2.5_rp   !mid !Minchole modification
           end if   
        else if(ituss_exm == EXM_CELLTYPE_ENDO) then  !ENDO
           gto = gto*ttparmate_exm(1,1,n)
           pca = pca*ttparmate_exm(1,9,n)  !endo
           gkr = gkr * ttparmate_exm(1,4,n)
           gks = gks* ttparmate_exm(1,2,n)     
           gk1 = gk1 * ttparmate_exm(1,3,n)
           gncx = gncx*ttparmate_exm(1,7,n)
           pnak = pnak*ttparmate_exm(1,10,n)
           gna = gna*ttparmate_exm(1,5,n)
           gnal = gnal*ttparmate_exm(1,6,n)
           gkb = gkb*ttparmate_exm(1,8,n)
           cmdnmax = cmdnmax*ttparmate_exm(1,11,n)  
        end if
        
        
        pcap = 1.1_rp * pca
        pcana = 0.00125_rp * pca
        pcak = 0.0003574_rp * pca
        pcanap = 0.00125_rp * pcap
        pcakp = 0.0003574_rp * pcap
        !%calculate inaca_i
        kasymm = 12.5_rp
        wna = 60000.0_rp
        wca = 60000.0_rp
        wnaca = 5000.0_rp
        kcaon = 1500000.0_rp
        kcaoff = 5000.0_rp
        qna = 0.5224_rp
        qca = 0.1670_rp
        kmcaact = 0.000150_rp !% and to calculate inaca_ss
        !%calculate inak
        k1p = 949.5_rp
        k1m = 182.4_rp
        k2p = 687.2_rp
        k2m = 39.4_rp
        k3m = 79300.0_rp
        k4m = 40.0_rp
        delta = -0.1550_rp
        kki = 0.5_rp
        kko = 0.3582_rp
        mgadp = 0.05_rp
        mgatp = 9.8_rp
        kmgatp = 0.0000001698_rp 
        hbig = 0.00000010_rp
        ep = 4.2_rp
        khp = 0.0000001698_rp 
        knap = 224.0_rp
        kxkur = 292.0_rp
        !%calculate inab
        pnab = 0.000000000375_rp
        !%calculate icab
        pcab = 0.000000025_rp
        !%calculate ipca
        gpca = 0.0005_rp
        kmcmdn = 0.00238_rp
        trpnmax = 0.07_rp
        kmtrpn = 0.0005_rp
        bsrmax = 0.047_rp
        kmbsr = 0.00087_rp
        bslmax = 1.124_rp
        kmbsl = 0.0087_rp
        csqnmax = 10.0_rp
        kmcsqn = 0.8_rp

        ccm = 1.0_rp


        !xistim = appfi_exm(ipoin)
        vicel_exm(23,ipoin,1) = 0.0_rp!xistim

        !% START OF CEIOHR CURRENT CALCULATIONS
        vaux1 = rgas*temp/farad
        ena = vaux1 * log(nao/vconc(5,ipoin,2))
        ek = vaux1 * log(ko/vconc(3,ipoin,2))
        eks = vaux1 *log((ko + pkna*nao)/(vconc(3,ipoin,2) + pkna*vconc(5,ipoin,2)))

        !%convenient shorthand calculations
        vffrt = elmag(ipoin,ITER_K)*farad*farad / (rgas*temp)
        vfrt = elmag(ipoin,ITER_K)*farad / (rgas*temp)

!!! First calculate the auxiliaries, Currents, then concentrations and then calculate new voltage,
        !!CHECK CAMKA
!!! Update CaMKa
        vaux2 = 1.0_rp / (1.0_rp + (kmcam/vconc(6,ipoin,2)))
        vicel_exm(26,ipoin,1) = vaux2 * camko* (1.0_rp - vconc(11,ipoin,2))
        !vicel_exm(26,2) = vicel_exm(26,1)
        vicel_exm(22,ipoin,1) = vicel_exm(26,ipoin,1) + vconc(11,ipoin,2)  ! camka
        !vicel_exm(22,2) = vicel_exm(22,1)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
!!! calculate exponential equations
        if (kfl_inaga_exm(n) == 1_ip) then !Passini
                call exm_ohara_infs(elmag(ipoin,ITER_K),inf_results, OHR_CELLMODEL_OHARA_PASSINI)
                call exm_ohara_taus(elmag(ipoin,ITER_K),tau_results, OHR_CELLMODEL_OHARA_PASSINI)
        else !normal o'hara
                call exm_ohara_infs(elmag(ipoin,ITER_K),inf_results, OHR_CELLMODEL_OHARA_RUDY)
                call exm_ohara_taus(elmag(ipoin,ITER_K),tau_results, OHR_CELLMODEL_OHARA_RUDY)
        end if

        !call exm_ohara_updateaux(inf_results,tau_results,elmag(ipoin,ITER_K),ko,vconc(6,ipoin,2),dtimeEP,vauxi_exm(:,ipoin,2),vauxi_exm(:,ipoin,1))
        call exm_ohara_updateaux(ipoin,inf_results,tau_results,elmag(ipoin,ITER_K),ko,vconc(6,ipoin,2),dtimeEP, ituss_exm)
        !! CALCULATE CONCENTRATIONS


        !%update camk
        rhsx1 = acamk * vicel_exm(26,ipoin,1) * (vicel_exm(26,ipoin,1) + vconc(11,ipoin,2))
        rhsx =  rhsx1 - (bcamk*vconc(11,ipoin,2))
        val0  =  vconc(11,ipoin,2)
        camkt = val0 + dtimeEP*rhsx
        !k1 = rhsx
        !k2 = rhsx + 0.5_rp * dtimeEP* k1
        !k3 = rhsx + 0.5_rp * dtimeEP * k2
        !k4 = rhsx + dtimeEP * k3
        !camkt = val0 + ((dtimeEP / (6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
        vconc(11,ipoin,1) = camkt      ! value of camkt concentration


        !%calculate ryanodione receptor calcium induced calcium release from the jsr
        bt = 4.75_rp
        a_rel = 0.5_rp * bt
        vaux1 = (1.5_rp/vconc(8,ipoin,2)) ** 8.0_rp
        vinf = a_rel * (-vicel_exm(4,ipoin,1)) / (1.0_rp + vaux1)
        if(ituss_exm == EXM_CELLTYPE_MID) then
           vinf = vinf * 1.7_rp
        end if
        xitaux = bt / (1.0_rp + 0.0123_rp/vconc(8,ipoin,2))

        if(xitaux < 0.001_rp) then  !fixed bug was wrongly 0.005
           xitaux = 0.001_rp
        end if
        vaux0 = vconc(9,ipoin,2)
        !vaux1 = (vinf - vaux0) * xitaux;
        !vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux;
        !vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux;
        !vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux;
        !vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4);
        !vaux3 = vaux5
        vaux3 = vinf - (vinf - vaux0) * safe_exp(-dtimeEP/xitaux)
        vconc(9,ipoin,1) = vaux3      ! value of variable jrelnp

        !!  JRELP                
        btp = 1.25_rp * bt
        a_relp = 0.5_rp * btp
        vaux1 = (1.5_rp/vconc(8,ipoin,2)) ** 8.0_rp 
        vinf = a_relp * (-vicel_exm(4,ipoin,1)) / (1.0_rp + vaux1)
        if(ituss_exm == EXM_CELLTYPE_MID) then
           vinf = vinf * 1.7_rp
        end if
        xitaux = btp / (1.0_rp+ (0.0123_rp/vconc(8,ipoin,2)))

        if(xitaux < 0.001_rp) then  !fixed bug was wrongly 0.005
           xitaux = 0.001_rp 
        end if
        vaux0 = vconc(10,ipoin,2)
        !vaux1 = (vinf - vaux0) * xitaux;
        !vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux;
        !vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux;
        !vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux;
        !vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
        !vaux3 = vaux5
        vaux3 = vinf - (vinf - vaux0) * safe_exp(-dtimeEP/xitaux)    
        vconc(10,ipoin,1) = vaux3      ! value of jrelp


!!!!!  CONCENTRATIONS!!!!   

        !%update intracellular concentrations, using buffers for cai, cass, cajsr
        !! calculate na current
        vaux1 =  vicel_exm(1,ipoin,1) + vicel_exm(2,ipoin,1) + vicel_exm(12,ipoin,1)
        vaux2 = 3.0_rp*vicel_exm(8,ipoin,1) + 3.0_rp*vicel_exm(10,ipoin,1)
        vaux3 = -(vaux1+vaux2) * acap/(farad*vmyo)
        rhsx = vaux3 + (vicel_exm(16,ipoin,1)*vss/vmyo)
        val0 = vconc(5,ipoin,2)
        nai = val0 + dtimeEP*rhsx
        !k1 = rhsx 
        !k2 = rhsx + 0.5_rp * dtimeEP* k1
        !k3 = rhsx + 0.5_rp * dtimeEP * k2 
        !k4 = rhsx + dtimeEP * k3 
        !nai = val0 + (dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
        vconc(5,ipoin,1) = nai

!!! calculate na current in subspace ss
        vaux1 = (vicel_exm(24,ipoin,1) + 3.0_rp * vicel_exm(9,ipoin,1)) * acap / (farad*vss)
        rhsx =  -vaux1 - vicel_exm(16,ipoin,1)
        val0 = vconc(2,ipoin,2)
        nass = val0 + dtimeEP*rhsx
        !k1 = rhsx 
        !k2 = rhsx + 0.5_rp * dtimeEP* k1
        !k3 = rhsx + 0.5_rp * dtimeEP * k2 
        !k4 = rhsx + dtimeEP * k3 
        !nass = val0 + (dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
        vconc(2,ipoin,1) = nass             !!! nass

!!! calculate k current
        vaux1 = vicel_exm(3,ipoin,1) + vicel_exm(5,ipoin,1) + vicel_exm(6,ipoin,1) + vicel_exm(23,ipoin,1)
        vaux2 = vicel_exm(7,ipoin,1) + vicel_exm(11,ipoin,1)  - (2.0_rp*vicel_exm(10,ipoin,1)) 
        vaux3 = (vicel_exm(17,ipoin,1)*vss/vmyo)
        rhsx = -((vaux1+ vaux2) * acap/(farad*vmyo)) + vaux3
        val0 = vconc(3,ipoin,2)
        ki = val0 + dtimeEP*rhsx
        !k1 = rhsx 
        !k2 = rhsx + 0.5_rp * dtimeEP* k1
        !k3 = rhsx + 0.5_rp * dtimeEP * k2 
        !k4 = rhsx + dtimeEP * k3 
        !ki = val0 + (dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4) 
        vconc(3,ipoin,1) = ki             !!! ki  

!!!!  calculate k current in the subspace ss               
        rhsx = -(vicel_exm(25,ipoin,1)*acap/(farad*vss)) - vicel_exm(17,ipoin,1)
        val0 = vconc(4,ipoin,2)
        kss = val0 + dtimeEP*rhsx
        !k1 = rhsx 
        !k2 = rhsx + 0.5_rp * dtimeEP* k1
        !k3 = rhsx + 0.5_rp * dtimeEP * k2 
        !k4 = rhsx + dtimeEP * k3 
        !kss = val0 + (dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
        vconc(4,ipoin,1) = kss             !!! kss

!!!  calculate ca current (Cai)

        if (flag_land) then

          ! Solve the coupled ORd and Land model
          !call exm_ohaland_calcium(kmcmdn,kmtrpn,cmdnmax,trpnmax,vnsr,vmyo,vss,acap,farad,dtimeEP,ipoin,cai,sac,elmag(ipoin,ITER_K))
          ! kmtrpn is not used inside
          call exm_ohaland_calcium(kmcmdn,cmdnmax,trpnmax,vnsr,vmyo,vss,acap,farad,dtimeEP,ipoin,cai,sac,elmag(ipoin,ITER_K))
          vconc(1,ipoin,1) = cai

        else
          vaux1 = (kmcmdn+vconc(1,ipoin,2)) * (kmcmdn+vconc(1,ipoin,2))
          vaux2 = (kmtrpn+vconc(1,ipoin,2)) * (kmtrpn+vconc(1,ipoin,2))
          vaux3 = 1.0_rp/(1.0_rp + (cmdnmax*kmcmdn/vaux1) + (trpnmax*kmtrpn/vaux2))
          rhsx1 = vicel_exm(14,ipoin,1)+vicel_exm(13,ipoin,1) - (2.0_rp*vicel_exm(8,ipoin,1))
          rhsx2 = -(vicel_exm(18,ipoin,1)*vnsr/vmyo) + (vicel_exm(15,ipoin,1)*vss/vmyo)
          rhsx = vaux3 *(-(rhsx1 *acap/(2.0_rp*farad*vmyo)) + rhsx2 )
          val0 = vconc(1,ipoin,2)
          cai = val0 + dtimeEP*rhsx
          !k1 = rhsx 
          !k2 = rhsx + 0.5_rp * dtimeEP* k1
          !k3 = rhsx + 0.5_rp * dtimeEP * k2 
          !k4 = rhsx + dtimeEP * k3 
          !cai = val0 + (dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
          vconc(1,ipoin,1) = cai             !!! cai
        end if

!!!! calculate ca current in the subspace ss
        vaux1 = (kmbsr+vconc(6,ipoin,2)) * (kmbsr+vconc(6,ipoin,2))
        vaux2 = (kmbsl+vconc(6,ipoin,2)) * (kmbsl+vconc(6,ipoin,2))
        vaux3 = 1.0_rp / (1.0_rp + (bsrmax*kmbsr/vaux1) + (bslmax*kmbsl/vaux2)) !Bcass
        rhsx1 = vicel_exm(4,ipoin,1) - (2.0_rp*vicel_exm(9,ipoin,1))
        rhsx2 = vicel_exm(21,ipoin,1)*vjsr/vss
        rhsx = vaux3 *(-rhsx1 *acap/(2.0_rp*farad*vss) + rhsx2 - vicel_exm(15,ipoin,1))
        val0 = vconc(6,ipoin,2)
        cass = val0 + dtimeEP*rhsx
        !k1 = rhsx 
        !k2 = rhsx + 0.5_rp * dtimeEP* k1
        !k3 = rhsx + 0.5_rp * dtimeEP * k2 
        !k4 = rhsx + dtimeEP * k3 
        !cass = val0 + (dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
        vconc(6,ipoin,1) = cass             !!! cass

!!! Update CaMKa
        vaux2 = 1.0_rp / (1.0_rp + (kmcam/vconc(6,ipoin,1)))
        vicel_exm(26,ipoin,1) = vaux2 * camko* (1.0_rp - vconc(11,ipoin,1))
        vicel_exm(22,ipoin,1) = vicel_exm(26,ipoin,1) + vconc(11,ipoin,1)  ! camka

        fjrelp= 1.0_rp / (1.0_rp + (kmcamk/vicel_exm(22,ipoin,1)))
        vicel_exm(21,ipoin,1)= (1.0_rp - fjrelp) * vconc(9,ipoin,1) + (fjrelp * vconc(10,ipoin,1)) !!! jrel


        !%calculate diffusion fluxes
        vicel_exm(16,ipoin,1)=(vconc(2,ipoin,1)-vconc(5,ipoin,1))*0.5_rp    !!! jdiffna
        vicel_exm(17,ipoin,1)=(vconc(4,ipoin,1)-vconc(3,ipoin,1))*0.5_rp    !!! jdiffk
        vicel_exm(15,ipoin,1)=(vconc(6,ipoin,1)-vconc(1,ipoin,1))*5.0_rp    !!! jdiff
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

!!!! calculate ca current in the sarcoplasmic reticulum nsr
        rhsx1 = vicel_exm(20,ipoin,1)*vjsr / vnsr
        rhsx = vicel_exm(18,ipoin,1) - rhsx1
        val0 = vconc(7,ipoin,2)
        cansr = val0 + dtimeEP*rhsx
        !k1 = rhsx 
        !k2 = rhsx + 0.5_rp * dtimeEP* k1
        !k3 = rhsx + 0.5_rp * dtimeEP * k2 
        !k4 = rhsx + dtimeEP * k3 
        !cansr = val0 + (dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
        vconc(7,ipoin,1) = cansr             !!! cansr

        !%calculate serca pump, ca uptake flux
        jupnp = 0.004375_rp * vconc(1,ipoin,1) / (vconc(1,ipoin,1) + 0.00092_rp)
        vaux1 = 1.0_rp / (vconc(1,ipoin,1) + 0.00092_rp - 0.00017_rp)
        jupp = 2.75_rp * 0.004375_rp * vconc(1,ipoin,1) * vaux1
        if (ituss_exm == EXM_CELLTYPE_EPI) then
            jupnp = jupnp*1.3_rp*ttparmate_exm(3, 12, n)
            jupp = jupp*1.3_rp*ttparmate_exm(3, 12, n)
        else if (ituss_exm == EXM_CELLTYPE_ENDO) then  !fixed this input 
            jupnp = jupnp*ttparmate_exm(1, 12, n)
            jupp = jupp*ttparmate_exm(1, 12, n)
        else if (ituss_exm == EXM_CELLTYPE_MID) then
            jupnp = jupnp*ttparmate_exm(2, 12, n)
            jupp = jupp*ttparmate_exm(2, 12, n)
        end if
        fjupp = (1.0_rp / (1.0_rp + kmcamk/vicel_exm(22,ipoin,1)))
        vicel_exm(19,ipoin,1) = 0.0039375_rp*vconc(7,ipoin,1)/15.0_rp   !!!! jleak
        !vicel_exm(19,3) = vicel_exm(19,2)
        vicel_exm(18,ipoin,1) = (1.0_rp-fjupp) * jupnp + fjupp * jupp - vicel_exm(19,ipoin,1)  !!! jup

!!!! calculate ca current in the junctional sarcoplasmic reticulum jsr
        vaux1 = (kmcsqn+vconc(8,ipoin,2)) * (kmcsqn+vconc(8,ipoin,2))
        vaux3 = csqnmax * kmcsqn / vaux1
        vaux2 = 1.0_rp / (1.0_rp + vaux3) 
        rhsx = vaux2 *(vicel_exm(20,ipoin,1)-vicel_exm(21,ipoin,1))
        val0 = vconc(8,ipoin,2)
        cajsr = val0 + dtimeEP*rhsx
        !k1 = rhsx 
        !k2 = rhsx + 0.5_rp * dtimeEP* k1
        !k3 = rhsx + 0.5_rp * dtimeEP * k2 
        !k4 = rhsx + dtimeEP * k3 
        !cajsr = val0 + ((dtimeEP / (6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
        vconc(8,ipoin,1) = cajsr             !!! cajsr

        !%calculate tranlocation flux
        vicel_exm(20,ipoin,1) = (vconc(7,ipoin,1)-vconc(8,ipoin,1)) * 0.01_rp  !!!! jtr          

        !Calculate INA
!!!!  hp
        hp = (ahf * vauxi_exm(2,ipoin,1)) + (ahs * vauxi_exm(5,ipoin,1))
        vaux1 = ahf * vauxi_exm(2,ipoin,1)
        hh = vaux1 + (ahs * vauxi_exm(3,ipoin,1))  !h
        finap = 1.0_rp / (1.0_rp + (kmcamk/vicel_exm(22,ipoin,1)))
        vaux1 = gna*(elmag(ipoin,ITER_K)-ena)
        vaux1 = vaux1 * vauxi_exm(1,ipoin,1) * vauxi_exm(1,ipoin,1) * vauxi_exm(1,ipoin,1)
        vaux2 = (1.0_rp-finap) * hh * vauxi_exm(4,ipoin,1)
        vaux2 = vaux2 + (finap * hp * vauxi_exm(6,ipoin,1))
        vicel_exm(1,ipoin,1) = vaux1 * vaux2  !!!!  ina

!!!  Calculate inal current
        finalp = 1.0_rp + (kmcamk/vicel_exm(22,ipoin,1))
        finalp = 1.0_rp / finalp
        vaux1 = gnal*(elmag(ipoin,ITER_K)-ena)*vauxi_exm(7,ipoin,1)
        vaux2 = (1.0_rp-finalp)*vauxi_exm(8,ipoin,1) + (finalp*vauxi_exm(9,ipoin,1))
        vicel_exm(2,ipoin,1)= vaux1 * vaux2

!!!!!!  Ito current
!!!! calculate aif, ais and ii
        aif = inf_results(OHR_A_I_FAST)
        ais = 1.0_rp - aif
        ii = aif*vauxi_exm(11,ipoin,1) + (ais*vauxi_exm(12,ipoin,1))
        ipp = aif*vauxi_exm(14,ipoin,1) + (ais*vauxi_exm(15,ipoin,1))
        vaux1 = kmcamk/vicel_exm(22,ipoin,1)
        fitop = 1.0_rp/(1.0_rp + vaux1)
        vaux2 = (fitop*vauxi_exm(13,ipoin,1) * ipp)
        vaux1 = ((1.0_rp-fitop) * vauxi_exm(10,ipoin,1) * ii ) + vaux2
        vaux1 = gto * (elmag(ipoin,ITER_K)-ek) * vaux1
        vicel_exm(3,ipoin,1) = vaux1        !!! ito

!!! calculate fca
        afcaf = inf_results(OHR_A_F_CA_FAST)
        afcas = 1.0_rp - afcaf
        fca = (afcaf*vauxi_exm(19,ipoin,1)) + (afcas*vauxi_exm(20,ipoin,1))
!!!!  calculate fcap
        fcap = afcaf*vauxi_exm(24,ipoin,1) + (afcas*vauxi_exm(20,ipoin,1))
!!!!! calculate fp
        fp = aff*vauxi_exm(23,ipoin,1) + (afs*vauxi_exm(18,ipoin,1))
        vaux2 = safe_exp(2.0_rp*vfrt)
        vaux1 = safe_exp(1.0_rp*vfrt)
        vaux3 = 1.0_rp / (vaux2-1.0_rp)
        phical = 4.0_rp * vffrt*((vconc(6,ipoin,1)*vaux2) - (0.341_rp*cao)) * vaux3
        vaux3 = 1.0_rp / (vaux1 -1.0_rp)
        phicana = 1.0_rp * vffrt*((0.75_rp*vconc(2,ipoin,1)*vaux1) - (0.75_rp*nao)) * vaux3
        phicak = 1.0_rp * vffrt*((0.75_rp*vconc(4,ipoin,1)*vaux1) - (0.75_rp*ko)) * vaux3
        ficalp = 1.0_rp / (1.0_rp + (kmcamk/vicel_exm(22,ipoin,1)))

!!!!!! calculate ff  = F
        ff = (aff * vauxi_exm(17,ipoin,1)) + (afs * vauxi_exm(18,ipoin,1))   !value of variable ff

!!!! calculate ical current
        vaux1 = (fp * (1.0_rp-vauxi_exm(22,ipoin,1)) + (vauxi_exm(21,ipoin,1)*fcap*vauxi_exm(22,ipoin,1)))
        vaux1 = ficalp * pcap * phical * vauxi_exm(16,ipoin,1) * vaux1
        vaux2 = ff * (1.0_rp-vauxi_exm(22,ipoin,1)) + (vauxi_exm(21,ipoin,1) * fca * vauxi_exm(22,ipoin,1))
        vaux3 = (1.0_rp-ficalp) * pca * phical * vauxi_exm(16,ipoin,1)
        vicel_exm(4,ipoin,1)= (vaux3 * vaux2) + vaux1  !!! ical

!!!! calculate icana current
        vaux1 = (fp*(1.0_rp-vauxi_exm(22,ipoin,1)) + (vauxi_exm(21,ipoin,1)*fcap*vauxi_exm(22,ipoin,1)))
        vaux1 = vaux1 * ficalp*pcanap*phicana*vauxi_exm(16,ipoin,1)
        vaux2 = (ff*(1.0_rp-vauxi_exm(22,ipoin,1))) + (vauxi_exm(21,ipoin,1)*fca*vauxi_exm(22,ipoin,1))
        vaux2 = vaux2 * (1.0_rp-ficalp)*pcana*phicana*vauxi_exm(16,ipoin,1)
        vicel_exm(24,ipoin,1)= vaux2 + vaux1  !! icana

!!!!! calculate icak
        vaux1 = fp * (1.0_rp-vauxi_exm(22,ipoin,1)) + (vauxi_exm(21,ipoin,1) * fcap * vauxi_exm(22,ipoin,1))
        vaux1 = vaux1 * ficalp * pcakp * phicak * vauxi_exm(16,ipoin,1)
        vaux2 = ff * (1.0_rp-vauxi_exm(22,ipoin,1)) + (vauxi_exm(21,ipoin,1) * fca * vauxi_exm(22,ipoin,1))
        vaux2 = vaux2 * (1.0_rp-ficalp) * pcak * phicak * vauxi_exm(16,ipoin,1)
        vicel_exm(25,ipoin,1) = vaux2 + vaux1   !!! icak

        !%calculate ikr CURRENT
        axrf = inf_results(OHR_A_XR_FAST)
        axrs = 1.0_rp - axrf
        xr = axrf*vauxi_exm(25,ipoin,1) + (axrs*vauxi_exm(26,ipoin,1))
        rkr = inf_results(OHR_R_KR_1) * inf_results(OHR_R_KR_2)
        vaux1 = sqrt(ko/5.4_rp) * (elmag(ipoin,ITER_K)-ek)
        vicel_exm(5,ipoin,1) = gkr * vaux1 * xr * rkr   !!! ikr

!!!! calculate iks CURRENTS
        ksca = (0.000038_rp/vconc(1,ipoin,1)) ** (1.4_rp) !!(7.0_rp/5.0_rp)
        ksca = 1.0_rp + (0.6_rp / (1.0_rp + ksca))
        vaux1 = gks * ksca * vauxi_exm(27,ipoin,1) * vauxi_exm(28,ipoin,1)
        vicel_exm(6,ipoin,1) =  vaux1 * (elmag(ipoin,ITER_K)-eks)

!!!! calculate ik1 CURRENT
        rk1 = 1.0_rp + safe_exp((elmag(ipoin,ITER_K) + 105.8_rp - 2.6_rp * ko) / 9.493_rp)
        rk1 = 1.0_rp / rk1
        vaux1 = gk1 * sqrt(ko) * rk1 * vauxi_exm(29,ipoin,1)
        vicel_exm(7,ipoin,1) = vaux1 * (elmag(ipoin,ITER_K)-ek) !!! ik1

        !%calculate INaCa CURRENT
        hca = safe_exp(qca * vfrt)
        hna = safe_exp(qna * vfrt)
        h1 = 1.0_rp + (vconc(5,ipoin,1)*INV_KNA3) * (1.0_rp +hna)
        h2 = (vconc(5,ipoin,1)*hna) / (KNA3*h1)
        h3 = 1.0_rp/h1
        h4 = 1.0_rp + vconc(5,ipoin,1) * INV_KNA1 * (1.0_rp + (vconc(5,ipoin,1) * INV_KNA2))
        h5 = vconc(5,ipoin,1)*vconc(5,ipoin,1) / (h4*KNA1*KNA2)
        h6 = 1.0_rp/h4
        h7 = 1.0_rp + nao*INV_KNA3*(1.0_rp + (1.0_rp/hna))
        h8 = nao / (KNA3*hna*h7)
        h9 = 1.0_rp / h7
        h10 = kasymm + 1.0_rp + nao*INV_KNA1 * (1.0_rp + (nao*INV_KNA2))
        h11 = nao*nao / (h10*KNA1*KNA2)
        h12 = 1.0_rp / h10
        k1 = h12*cao*kcaon
        k2 = kcaoff
        k3p = h9*wca
        k3pp = h8*wnaca
        k3 = k3p + k3pp
        k4p = h3*wca / hca
        k4pp = h2*wnaca
        k4 =k4p+k4pp
        k5 = kcaoff
        k6 = h6*vconc(1,ipoin,1)*kcaon
        k7 = h5*h2*wna
        k8 = h8*h11*wna
        x1 = (k2*k4*(k7+k6)) + (k5*k7*(k2+k3))
        x2 = (k1*k7*(k4+k5)) + (k4*k6*(k1+k8))
        x3 = (k1*k3*(k7+k6)) + (k8*k6*(k2+k3))
        x4 = (k2*k8*(k4+k5)) + (k3*k5*(k1+k8))
        e1 = x1 / (x1+x2+x3+x4)
        e2 = x2 / (x1+x2+x3+x4)
        e3 = x3 / (x1+x2+x3+x4)
        e4 = x4 / (x1+x2+x3+x4)
        allo = (kmcaact/vconc(1,ipoin,1)) * (kmcaact/vconc(1,ipoin,1))
        allo = 1.0_rp / (1.0_rp + allo)
        zna = 1.0_rp
        jncxna = (3.0_rp*(e4*k7-e1*k8)) + (e3*k4pp) - (e2*k3pp)
        jncxca = (e2*k2) - (e1*k1)
        vaux1 = (zna*jncxna) + (zca*jncxca)
        vicel_exm(8,ipoin,1) = 0.8_rp * gncx * allo * vaux1  !!! inaca_i Current

        !%calculate inaca_ss
        h1 = 1.0_rp + vconc(2,ipoin,1)*INV_KNA3 * (1.0_rp + hna)
        h2 = (vconc(2,ipoin,1)*hna) / (KNA3*h1)
        h3 = 1.0_rp/h1
        h4 = 1.0_rp + vconc(2,ipoin,1)*INV_KNA1 * (1.0_rp + (vconc(2,ipoin,1)*INV_KNA2))
        h5 = vconc(2,ipoin,1)*vconc(2,ipoin,1) / (h4*KNA1*KNA2)
        h6 = 1.0_rp / h4
        h7 = 1.0_rp + (nao*INV_KNA3) *(1.0_rp + (1.0_rp/hna))
        h8 = nao / (KNA3*hna*h7)
        h9 = 1.0_rp / h7
        h10 = kasymm + 1.0_rp + nao*INV_KNA1 * (1.0_rp + (nao*INV_KNA2))
        h11 = nao*nao / (h10*KNA1*KNA2)
        h12 = 1.0_rp / h10
        k1 = h12*cao*kcaon
        k2 = kcaoff
        k3p = h9*wca
        k3pp = h8*wnaca
        k3 = k3p + k3pp
        k4p = h3*wca / hca
        k4pp = h2*wnaca
        k4 = k4p + k4pp
        k5 = kcaoff
        k6 = h6*vconc(6,ipoin,1)*kcaon
        k7 = h5*h2*wna
        k8 = h8*h11*wna
        x1 = (k2*k4*(k7+k6)) + (k5*k7*(k2+k3))
        x2 = (k1*k7*(k4+k5)) + (k4*k6*(k1+k8))
        x3 = (k1*k3*(k7+k6)) + (k8*k6*(k2+k3))
        x4 = (k2*k8*(k4+k5)) + (k3*k5*(k1+k8))
        e1 = x1 / (x1+x2+x3+x4)
        e2 = x2 / (x1+x2+x3+x4)
        e3 = x3 / (x1+x2+x3+x4)
        e4 = x4 / (x1+x2+x3+x4)
        allo = (kmcaact/vconc(6,ipoin,1)) * (kmcaact/vconc(6,ipoin,1))
        allo = 1.0_rp / (1.0_rp + allo)
        jncxna = (3.0_rp*(e4*k7-e1*k8)) + (e3*k4pp) - (e2*k3pp)
        jncxca = (e2*k2) - (e1*k1)
        vicel_exm(9,ipoin,1) = 0.2_rp * gncx * allo * ((zna*jncxna) + (zca*jncxca))   !!! inaca_ss

        INaCa = vicel_exm(8,ipoin,1)+vicel_exm(9,ipoin,1)

        !%calculate inak
        k3p = 1899.0_rp
        k4p = 639.0_rp
        knai = KNAI0 * safe_exp(delta*vfrt/3.0_rp)
        knao = KNAO0 * safe_exp((1.0_rp-delta)*vfrt/3.0_rp)
        pp = ep / (1.0_rp + hbig/khp + (vconc(5,ipoin,1)/knap) + (vconc(3,ipoin,1)/kxkur))
        vaux1 = (vconc(5,ipoin,1)/knai) * (vconc(5,ipoin,1)/knai) * (vconc(5,ipoin,1)/knai)
        vaux2 = (1.0_rp + (vconc(5,ipoin,1)/knai)) * (1.0_rp + (vconc(5,ipoin,1)/knai)) * (1.0_rp + (vconc(5,ipoin,1)/knai))
        vaux3 = (1.0_rp + (vconc(3,ipoin,1)/kki)) * (1.0_rp + (vconc(3,ipoin,1)/kki))
        a1 = (k1p*vaux1) / (vaux2 + vaux3 - 1.0_rp)
        b1 = k1m * mgadp
        a2 = k2p
        vaux1 = (nao/knao) * (nao/knao) * (nao/knao)
        vaux2 = (1.0_rp + (nao/knao)) * (1.0_rp + (nao/knao)) * (1.0_rp + (nao/knao))
        vaux3 = (1.0_rp + (ko/kko)) * (1.0_rp + (ko/kko))
        b2 = (k2m*vaux1)/(vaux2 + vaux3 -1.0_rp)
        vaux1 = (ko/kko) * (ko/kko)
        a3 = (k3p * vaux1) / (vaux2 + vaux3 - 1.0_rp)
        b3 = (k3m*pp*hbig) / (1.0_rp + (mgatp/kmgatp))
        a4 = (k4p*mgatp/kmgatp) / (1.0_rp + (mgatp/kmgatp))
        vaux1 = (vconc(3,ipoin,1)/kki) * (vconc(3,ipoin,1)/kki)
        vaux2 = (1.0_rp + (vconc(5,ipoin,1)/knai)) * (1.0_rp + (vconc(5,ipoin,1)/knai)) * (1.0_rp + (vconc(5,ipoin,1)/knai))
        vaux3 = (1.0_rp + (vconc(3,ipoin,1)/kki)) * (1.0_rp + (vconc(3,ipoin,1)/kki))
        b4 = (k4m*vaux1)/(vaux2 + vaux3 -1.0_rp)
        x1 = (a4*a1*a2) + (b2*b4*b3) + (a2*b4*b3) + (b3*a1*a2)
        x2 = (b2*b1*b4) + (a1*a2*a3) + (a3*b1*b4) + (a2*a3*b4)
        x3 = (a2*a3*a4) + (b3*b2*b1) + (b2*b1*a4) + (a3*a4*b1)
        x4 = (b4*b3*b2) + (a3*a4*a1) + (b2*a4*a1) + (b3*b2*a1)
        e1 = x1/(x1+x2+x3+x4)
        e2 = x2/(x1+x2+x3+x4)
        e3 = x3/(x1+x2+x3+x4)
        e4 = x4/(x1+x2+x3+x4)
        zk = 1.0_rp
        jnakna = 3.0_rp*((e1*a3)-(e2*b3))
        jnakk = 2.0_rp*((e4*b1)-(e3*a1))
        vicel_exm(10,ipoin,1) = pnak*(zna*jnakna + zk*jnakk)  !!!! inak

        !%calculate ikb CURRENT
        xkb = inf_results(OHR_X_KB)
        vicel_exm(11,ipoin,1) = gkb*xkb*(elmag(ipoin,ITER_K)-ek)  !!!! ikb

        !%calculate inab CURRENT
        vaux1 = pnab*vffrt*(vconc(5,ipoin,1)*safe_exp(vfrt) - nao)
        vicel_exm(12,ipoin,1) =  vaux1 / (safe_exp(vfrt)-1.0_rp)

        !%calculate icab CURRENT
        vaux1 = (vconc(1,ipoin,1)*safe_exp(2.0_rp*vfrt) - (0.341_rp*cao))/(safe_exp(2.0_rp*vfrt)-1.0_rp)
        vicel_exm(13,ipoin,1) = pcab * 4.0_rp*  vffrt * vaux1   !!!!icab

        !%calculate ipca CURRENT
        vicel_exm(14,ipoin,1) = gpca * vconc(1,ipoin,1) / (0.0005_rp+vconc(1,ipoin,1))  !!! ipca


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !!  VOLTAGE INTEGRATION
        v1 = 0.0_rp
        v2 = 0.0_rp
        Inet = 0.0_rp 
        do i=1,14
           v1= v1 + vicel_exm(i,ipoin,1)
        end do
        do i=23,25
           v2= v2 + vicel_exm(i,ipoin,1)
        end do
        do i=2,7
           Inet= Inet + vicel_exm(i,ipoin,1)
        end do
        xioni=v1+v2+sac
        qneto_exm(ipoin)= qneto_exm(ipoin)+ (Inet*dtimeEP)

        !!  Runge Kutta, now done outside
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        vconc(:,ipoin,3)=vconc(:,ipoin,2)         
        vconc(:,ipoin,2)=vconc(:,ipoin,1) 
        vauxi_exm(:,ipoin,2)=vauxi_exm(:,ipoin,1)         
        !vauxi_exm(:,ipoin,2)=vauxi_exm(:,ipoin,1)                 

        !end do
        dioni = 0.0_rp
     end if
  end if

end subroutine exm_oharaf_ionicurrents


!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_ohara_updateaux.f90
!> @author  Mariano Vazquez
!> @brief   Updating cell model variables
!> @date   16/NOV/1966
!> @details Updating cell model variables
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_ohara_updateaux(ipoin,inf_results,tau_results,elmag_local,ko,conc6,dtimeEP,ituss_exm)
   use def_kintyp, only : ip,rp
   use def_exmedi
   use mod_exm_oharaequations
   implicit none

   integer(ip), intent(in) :: ipoin, ituss_exm !< node
   real(rp), intent(in)  :: inf_results(EXM_OHR_INF_EQUATIONS) !< inf
   real(rp), intent(in)  :: tau_results(EXM_OHR_TAU_EQUATIONS) !< tau
   real(rp), intent(in)  :: elmag_local !< elmag
   real(rp), intent(in)  :: ko !< ko
   real(rp), intent(in)  :: conc6 !< conc6
   real(rp), intent(in)  :: dtimeEP !< dtimeEP

   integer(ip) :: iauxi
   real(rp) :: infs(nauxi_exm),taus(nauxi_exm)
   real(rp) :: anca,dnca,km2n,vaux1,vaux4

   real(rp), parameter   :: KMN = 0.002_rp
   real(rp), parameter   :: K2N = 1000.0_rp

   infs(1) = inf_results(OHR_M_INF)
   taus(1) = tau_results(OHR_TAU_M)

   infs(2) = inf_results(OHR_H_INF)
   taus(2) = tau_results(OHR_TAU_H_FAST)

   infs(3) = inf_results(OHR_H_INF)
   taus(3) = tau_results(OHR_TAU_H_SLOW)

   infs(4) = inf_results(OHR_H_INF)  ! j_inf = h_inf
   taus(4) = tau_results(OHR_TAU_J)

   infs(5) = inf_results(OHR_H_CAMK_INF)
   taus(5) = 3.0_rp * tau_results(OHR_TAU_H_SLOW)

   infs(6) = inf_results(OHR_H_INF)  ! j_camk_inf = j_inf = h_inf
   taus(6) = 1.46_rp * tau_results(OHR_TAU_J)

   infs(7) = inf_results(OHR_M_L_INF)
   taus(7) = tau_results(OHR_TAU_M)  ! tau_m_l = tau_m

   infs(8) = inf_results(OHR_H_L_INF)
   taus(8) = 200.0_rp

   infs(9) = inf_results(OHR_H_L_CAMK_INF)
   taus(9) = 3.0_rp * 200.0_rp

   ! NOTE: OHR_INV_TAU_A_1/2 are infs, not taus!
   infs(10) = inf_results(OHR_A_INF)
   taus(10) = 1.0515_rp / (inf_results(OHR_INV_TAU_A_1) + inf_results(OHR_INV_TAU_A_2))

   infs(11) = inf_results(OHR_I_INF)
   taus(11) = tau_results(OHR_TAU_I_FAST)
   if(ituss_exm == EXM_CELLTYPE_EPI) then
      taus(11) = taus(11) * inf_results(OHR_DELTA_EPI)
   end if

   infs(12) = inf_results(OHR_I_INF)
   taus(12) = tau_results(OHR_TAU_I_SLOW)
   if(ituss_exm == EXM_CELLTYPE_EPI) then
      taus(12) = taus(12) * inf_results(OHR_DELTA_EPI)
   end if

   infs(13) = inf_results(OHR_A_CAMK_INF)
   taus(13) = taus(10) ! tau_a_camk_inf = tau_a

   infs(14) = inf_results(OHR_I_INF) ! i_camk_inf = i_inf
   taus(14) = tau_results(OHR_DELTA_CAMK_DEVELOP) * inf_results(OHR_DELTA_CAMK_RECOVER) * taus(11)

   infs(15) = inf_results(OHR_I_INF) ! i_camk_inf = i_inf
   taus(15) = tau_results(OHR_DELTA_CAMK_DEVELOP) * inf_results(OHR_DELTA_CAMK_RECOVER) * taus(12)

   infs(16) = inf_results(OHR_D_INF)
   taus(16) = tau_results(OHR_TAU_D)

   infs(17) = inf_results(OHR_F_INF)
   taus(17) = tau_results(OHR_TAU_F_FAST)

   infs(18) = inf_results(OHR_F_INF)
   taus(18) = tau_results(OHR_TAU_F_SLOW)

   infs(19) = inf_results(OHR_F_INF) ! f_ca_inf = f_inf
   taus(19) = tau_results(OHR_TAU_F_CA_FAST)

   infs(20) = inf_results(OHR_F_INF) ! f_ca_inf = f_inf
   taus(20) = tau_results(OHR_TAU_F_CA_SLOW)

   infs(21) = inf_results(OHR_F_INF) ! j_ca_inf = f_ca_inf = f_inf
   taus(21) = 75.0_rp

   infs(23) = inf_results(OHR_F_INF) ! f_camk_inf = f_inf
   taus(23) = 2.5_rp * tau_results(OHR_TAU_F_FAST)

   infs(24) = inf_results(OHR_F_INF) ! f_ca_camk_inf = f_inf
   taus(24) = 2.5_rp * tau_results(OHR_TAU_F_CA_FAST)

   infs(25) = inf_results(OHR_X_R_INF)
   taus(25) = tau_results(OHR_TAU_XR_FAST)

   infs(26) = inf_results(OHR_X_R_INF)
   taus(26) = tau_results(OHR_TAU_XR_SLOW)

   infs(27) = inf_results(OHR_X_S1_INF)
   taus(27) = tau_results(OHR_TAU_X_S1)

   infs(28) = inf_results(OHR_X_S1_INF)
   taus(28) = tau_results(OHR_TAU_X_S2)

   infs(29) = 1.0_rp / (1.0_rp + safe_exp(-(elmag_local + 2.5538_rp * ko + 144.59_rp) / (1.5692_rp * ko + 3.8115_rp)))
   taus(29) = tau_results(OHR_TAU_X_K1)

   do iauxi=1,21
      vauxi_exm(iauxi,ipoin,1) = infs(iauxi) - (infs(iauxi) - vauxi_exm(iauxi,ipoin,2)) * safe_exp(-dtimeEP/taus(iauxi))
   end do
   do iauxi=23,29
      vauxi_exm(iauxi,ipoin,1) = infs(iauxi) - (infs(iauxi) - vauxi_exm(iauxi,ipoin,2)) * safe_exp(-dtimeEP/taus(iauxi))
   end do
        
   !!!!! calculate icana and icak
   km2n = vauxi_exm(21,ipoin,1)
   vaux1 = KMN / conc6
   vaux4 = (1.0_rp+vaux1) *(1.0_rp+vaux1) * (1.0_rp+vaux1) * (1.0_rp+vaux1)
   anca = 1.0_rp / ((K2N / km2n) + vaux4)
   dnca = (anca*k2n/km2n)
   vauxi_exm(22,ipoin,1) = dnca - (dnca - vauxi_exm(22,ipoin,2)) * safe_exp(-dtimeEP*km2n)

end subroutine exm_ohara_updateaux



