!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_oceohr.f90
!> @author  Jazmin Aguado-Sierra
!> @brief   Single cell run for Initial condition setup for Ohara-Rudy 2011 heterogeneous model
!> @details Runs a single cell simulation at th given frequency and pathologic conditions \n
!!    It performs single cell runs under normal, heart failure or drugs \n
!> @} 
!!-----------------------------------------------------------------------
logical function exm_oceohr(mat) 

  use      def_parame
  use      def_master
  use      def_elmtyp
  use      def_domain
  use      def_exmedi
  
    ! definition of variables
  implicit none
  integer(ip), intent(in) :: mat

  integer(ip) :: nbeat, j, i, itim, imate, nsamples_per_beat
  real(rp)    :: vinf, xitaux, vaux0, vaux1, vaux2, vaux3, vaux4, vaux5
  real(rp)    ::  vffrt, vfrt, ena, ek, eks, bt, a_rel, btp, a_relp, ccm
  real(rp)    :: delepi, ta, devel, recov, fss, batim, dur, aux1, atim
  real(rp)    :: pkna, farad, iss,  tif, tis, tj, v1, v2
  real(rp)    ::  a1, a2, a3, a4, b1, b2, b3, b4, k1p, k2p, k3p, k4p, k1m, k2m, k3m, k4m
  real(rp)    :: h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, x1, x2, x3, x4
  real(rp)    ::  k1, k2, k3, k4, k5, k6, k7, k8, e1, e2, e3, e4, k3pp, k4pp, knai, knao0
  real(rp)    :: gna, gnal, gto, gks, gk1, gncx, gkb, gkr, gpca, phical, phicana, phicak, nao
  real(rp)    :: ahf, ahs, hp, hh, finalp, ii, aff, afs, afcaf, afcas, fcap, kmn, k2n, km2n
  real(rp)    :: jss, ths, tm, zca, pca, pcap, pcak, pcanap, pcakp, ksca, kna1, kna2, kna3, kasymm
  real(rp)    :: wna, wca, wnaca, kcaon, kcaoff, qna, qca, hca, hna, kmcaact, allo, zna, jncxna, jncxca
  real(rp)    :: knai0, knao, delta, kki, kko, mgadp, mgatp, kmgatp, ep, khp, knap, kxkur, pp, zk
  real(rp)    :: jnakna, jnakk, pnak, xkb, pnab, ipp, aif, ais
  real(rp)    :: tff, tfcaf, fp, ko, axrf, axrs, xr, rkr, rk1, fjupp, pcab, pcana
  real(rp)    :: anca, cao, dnca, ff, fca, ficalp, finap,  fjrelp
  real(rp)    :: rgas, temp, jupnp, jupp, fitop, hbig, camko, kmcam
  real(rp)    :: cai, cass, cansr, cajsr, nass, nai, ki, kss, camkt
  real(rp)    :: rhsx1, rhsx2, rhsx, val0, bslmax, bsrmax,cmdnmax, csqnmax
  real(rp)    :: leng, rad, vcell, ageo, acap, vmyo, vnsr, vjsr, vss, kmcamk, kmcmdn
  real(rp)    :: acamk, bcamk, trpnmax, kmcsqn, kmtrpn, kmbsl, kmbsr, INaCa
  real(rp)    :: dosis1, dosis2, dosis3, i50cal, i50na, i50kr, gcaldrug, gkrdrug, gnadrug, dtimon
  real(rp)    :: dosis4, i50ito, gitodrug
      
  !arrays to store voltages of the last 2 beats
  !voltages_last will always point ot the last elemnt of the array, 
  !when reached the end of the array, the voltages_last restarts from the beginning
  real(rp), dimension (:), allocatable :: ss_voltages_beat_current, ss_voltages_beat_previous
  
  !number of voltages actually saved in the array, needed to verify that there was no error
  integer(ip) :: ss_voltage_current_id
  real(rp) :: ss_voltages_epsilon ! threshold to declare steady statewhen comparing two beats
  real(rp) :: ss_rmse, ss_min_voltage_threshold, ss_min_voltage
  integer(ip) :: ss_rmse_counter, ss_rmse_counter_max !the rmse has to maintain low value for ss_rmse_counter_max number of consecutive beats


  
  !!!!  DRUG DEFINITION TO IKR, INA OR ICAL
 if(kfl_drugsmate_exm(mat)==1_ip)then
    dosis1 = drugdmate_exm(1,mat)
    i50cal = drugdmate_exm(2,mat)
    dosis2 = drugdmate_exm(3,mat)
    i50kr = drugdmate_exm(4,mat)
    dosis3 = drugdmate_exm(5,mat)
    i50na = drugdmate_exm(6,mat)
    dosis4 = drugdmate_exm(7,mat)
    i50ito = drugdmate_exm(8,mat)
    gcaldrug = 1.0_rp / (1.0_rp + (dosis1/i50cal))
    gkrdrug = 1.0_rp / (1.0_rp + (dosis2/i50kr))
    gnadrug = 1.0_rp / (1.0_rp + (dosis3/i50na))
    gitodrug = 1.0_rp / (1.0_rp + (dosis4/i50ito))
 else
    gcaldrug = 1.0_rp
    gkrdrug = 1.0_rp
    gnadrug = 1.0_rp
    gitodrug = 1.0_rp
 endif 
  !extracellular ionic concentrations
  nao = 140.0_rp
  cao = 1.8_rp
  ko = 5.4_rp

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
  gna = 75.0_rp * gnadrug
  gnal = 0.0075_rp
  gto = 0.02_rp * gitodrug !
!!! calculate ical
!!!!!! calculate ff
  aff = 0.6_rp
  afs = 1.0_rp - aff  
!!!!! calculate icana and icak
  kmn = 0.002_rp
  k2n = 1000.0_rp
  zca = 2.0_rp
  pca = 0.0001_rp * gcaldrug !endo 
!!  inal  
!!!  calculate IKr
  gkr = 0.046_rp * gkrdrug
!!!  calculate IKs
  gks = 0.0034_rp
 !!!! calculate IK1  
  gk1 = 0.1908_rp
  gncx = 0.0008_rp
  pnak = 30.0_rp
  !!calculate IKb
  gkb = 0.003_rp
  !%calcium buffer constants
  cmdnmax = 0.05_rp
  if(ituss_exm == 3) then !!epi
      gnal = gnal * 0.6_rp*ttparmate_exm(3,6,mat)
      gto = gto*4.0_rp*ttparmate_exm(3,1,mat)
      pca = pca*1.2_rp*ttparmate_exm(3,9,mat)  !epi
      gkr = gkr * 1.3_rp*ttparmate_exm(3,4,mat)
      gk1 = gk1 * 1.2_rp*ttparmate_exm(3,3,mat)
      gncx = gncx*1.1_rp*ttparmate_exm(3,7,mat)
      pnak = pnak*0.9_rp*ttparmate_exm(3,10,mat)
      gkb = gkb*0.6_rp*ttparmate_exm(3,8,mat)
      cmdnmax = cmdnmax*1.3_rp*ttparmate_exm(3,11,mat) 
      gks = gks*ttparmate_exm(3,2,mat)
      gna = gna*ttparmate_exm(3,5,mat) 
      gkb = gkb*ttparmate_exm(3,8,mat)
   else if(ituss_exm == 2) then !!MID
      gto = gto*4.0_rp*ttparmate_exm(2,1,mat)
      pca = pca*2.5_rp*ttparmate_exm(2,9,mat)  !mid 
      gkr = gkr * 0.8_rp*ttparmate_exm(2,4,mat)
      gks = gks*1.4_rp*ttparmate_exm(2,2,mat)     
      gk1 = gk1 * 1.3_rp*ttparmate_exm(2,3,mat)
      gncx = gncx*1.4_rp*ttparmate_exm(2,7,mat)
      pnak = pnak*0.7_rp*ttparmate_exm(2,10,mat)
      gna = gna*ttparmate_exm(2,5,mat) 
      gnal = gnal*ttparmate_exm(2,6,mat)
      gkb = gkb*ttparmate_exm(2,8,mat) 
      cmdnmax = cmdnmax*ttparmate_exm(3,11,mat)      
      if(kfl_inaga_exm(mat) == 1_ip) then
         pca=(pca*1.8_rp)/2.5_rp   !mid !Minchole modification
      end if
  else if(ituss_exm == 1) then  !ENDO
      gto = gto*ttparmate_exm(1,1,mat)
      pca = pca*1.0_rp*ttparmate_exm(1,9,mat)  !endo
      gkr = gkr * ttparmate_exm(1,4,mat)
      gks = gks* ttparmate_exm(1,2,mat)     
      gk1 = gk1 * ttparmate_exm(1,3,mat)
      gncx = gncx*ttparmate_exm(1,7,mat)
      pnak = pnak*ttparmate_exm(1,10,mat)
      gna = gna*ttparmate_exm(1,5,mat) 
      gnal = gnal*ttparmate_exm(1,6,mat)
      gkb = gkb*ttparmate_exm(1,8,mat)
      cmdnmax = cmdnmax*ttparmate_exm(3,11,mat)  
   end if
  pcap = 1.1_rp * pca
  pcana = 0.00125_rp * pca
  pcak = 0.0003574_rp * pca
  pcanap = 0.00125_rp * pcap
  pcakp = 0.0003574_rp * pcap
!%calculate inaca_i
  kna1 = 15.0_rp
  kna2 = 5.0_rp
  kna3 = 88.12_rp
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
  knai0 = 9.073_rp
  knao0 = 27.78_rp
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
  
  dtimon = 0.005_rp   !0.0018_rp  0.0015_rp
  nbeat = 0
  itim = 0
  batim = 0
  j = 0
  ccm = 1.0_rp
  atim = 0.0_rp

  !----------------------------------------------
  !
  ! begin: Initialize variables for steady state detection
  !
  !----------------------------------------------
  !allocate the arrays to store two full beats
  nsamples_per_beat = nint(moneclmate_exm(2,mat) / dtimon)
  allocate ( ss_voltages_beat_current(nsamples_per_beat) )   
  allocate ( ss_voltages_beat_previous(nsamples_per_beat) ) 
  ss_voltages_beat_current = 0.0_rp
  ss_voltages_beat_previous = 0.0_rp
  ss_voltages_epsilon = 1.0_rp !this is maximum allowed  mean voltage difference
  ss_rmse_counter = 0_ip
  ss_rmse_counter_max = 3_ip
  ss_min_voltage_threshold = -75.0_rp !the min volytage within the beat has to be lower than this, to consider that we are in steady state
  !----------------------------------------------
  !
  ! end: Initialize variables for steady state detection
  !
  !----------------------------------------------

  exm_oceohr = .FALSE.


    do while (batim < (moneclmate_exm(1,mat) * moneclmate_exm(2,mat)))
        aux1 = 0.0_rp
        atim = dtimon * aux1
        

        !for steady state verification 
        ss_voltage_current_id = 1_ip


        do while (atim < moneclmate_exm(2,mat))
        !if (atim < moneclmate_exm(2) )then            
            aux1 = aux1 + 1.0_rp
            atim = dtimon * aux1
            dur = 0.5_rp ! duration time (ms)
            if (atim < dur) then 
            ! External stimulus current I_stim [pA/pF]
                 viclo_exm(23,1) = -80.0_rp
                 !viclo_exm(23,3) = -80.0_rp     ! Value of I_stim current [pA/pF]
            else
                viclo_exm(23,1) = 0.0_rp
                !viclo_exm(23,3) = 0.0_rp
            end if
                !% START OF CEIOHR CURRENT CALCULATIONS
                vaux1 = rgas*temp/farad
                ena = vaux1 * log(nao/vcolo_exm(5,2))
                ek = vaux1 * log(ko/vcolo_exm(3,2))
                eks = vaux1 *log((ko + pkna*nao)/(vcolo_exm(3,2) + pkna*vcolo_exm(5,2)))
           
                !%convenient shorthand calculations
                vffrt = elmlo_exm(2)*farad*farad / (rgas*temp)
                vfrt = elmlo_exm(2)*farad / (rgas*temp)
                
           !!! First calculate the auxiliaries, Currents, then concentrations and then calculate new voltage,
            
           !!! Update CaMKa
                vaux2 = 1.0_rp / (1.0_rp + (kmcam/vcolo_exm(6,2)))
                viclo_exm(26,1) = vaux2 * camko* (1.0_rp - vcolo_exm(11,2))
                !viclo_exm(26,2) = viclo_exm(26,1)
                viclo_exm(22,1) = viclo_exm(26,1) + vcolo_exm(11,2)  ! camka
                !viclo_exm(22,2) = viclo_exm(22,1)
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
           !!!! CALCULATE GATES
           !!!!!!   calculate INA
                !calculate m gate
                vinf = exp(-(elmlo_exm(2) + 39.57_rp) / 9.871_rp)
                vinf = 1.0_rp / (1.0_rp + vinf)
                vaux1 = 8.552_rp * exp(-(elmlo_exm(2) + 77.42_rp) / 5.955_rp)
                xitaux = vaux1 + (6.765_rp * exp((elmlo_exm(2) + 11.64_rp) / 34.77_rp))
                tm = 1.0_rp/ xitaux
                vaux0 = vaulo_exm(1,2)
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/tm)    
                vaulo_exm(1,1) = vaux3 ! value of variable m
                
           !!!!!!! hf
                if(kfl_inaga_exm(mat) == 1_ip) then
                   vinf = exp((elmlo_exm(2) + 78.50_rp) / 6.22_rp) !modified from Passini et al. with LMS no reg(original: exp((elmlo_exm(2) + 82.90_rp) / 6.086_rp))
                   vinf = 1.0_rp / (1.0_rp + vinf)
                   vaux1 = exp(-(elmlo_exm(2) + 3.8875_rp) / 7.8579_rp) !modified from Passini et al. with LMS no reg(original: exp(-(elmlo_exm(2) + 1.196_rp) / 6.285_rp))
                   vaux2 = exp((elmlo_exm(2) + 0.4963_rp) / 9.1843_rp) !modified from Passini et al. with LMS no reg(original: exp((elmlo_exm(2) + 0.5096_rp) / 20.27_rp))
                   xitaux = (0.000003686_rp * vaux1) + (16.0_rp * vaux2) !modified from Passini et al. with LMS no reg(original: (0.00001432_rp * vaux1) + (6.149_rp * vaux2) )
                else
                   vinf = exp((elmlo_exm(2) + 82.90_rp) / 6.086_rp) !original O'Hara-Rudy formulation
                   vinf = 1.0_rp / (1.0_rp + vinf)
                   vaux1 = exp(-(elmlo_exm(2) + 1.196_rp) / 6.285_rp) !original O'Hara-Rudy formulation
                   vaux2 = exp((elmlo_exm(2) + 0.5096_rp) / 20.27_rp) !origina O'Hara-Rudy formulation
                   xitaux = (0.00001432_rp * vaux1) + (6.149_rp * vaux2) !original O'Hara-Rudy formulation
                end if
                xitaux = 1.0_rp / xitaux
                vaux0 = vaulo_exm(2,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(2,1) = vaux3      ! value of variable hf
           
                jss = vinf    
           !!!! hs
                vaux1 = exp(-(elmlo_exm(2) + 17.95_rp) / 28.05_rp)
                vaux2 = exp((elmlo_exm(2) + 5.730_rp) / 56.66_rp)
                xitaux = (0.009794_rp * vaux1) + (0.3343_rp * vaux2) 
                ths = 1.0_rp / xitaux
                vaux0 = vaulo_exm(3,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/ths)    
                vaulo_exm(3,1) = vaux3      ! value of variable hs    
           
           !!!!!  j
                if(kfl_inaga_exm(mat) == 1_ip) then
                   vaux1 = exp(-(elmlo_exm(2) + 116.7258_rp) / 7.6005_rp) !modified from Passini et al. with LMS no reg(original: exp(-(v+100.6)/8.281))
                   vaux2 = exp((elmlo_exm(2) + 6.2719_rp) / 9.0358_rp) !modified from Passini et al. with LMS no reg(original:  exp((elmlo_exm(2) + 0.9941_rp) / 38.45_rp))
                   xitaux = (0.8628_rp * vaux1) + (1.1096_rp * vaux2) !modified from Passini et al. with LMS no reg(original: (0.02136_rp * vaux1) + (0.3052_rp * vaux2))
                   tj = 4.8590_rp + (1.0_rp / xitaux) !modified from Passini et al. with LMS no reg(original: 2.038_rp + (1.0_rp / xitaux))
                else
                   vaux1 = exp(-(elmlo_exm(2) + 100.6_rp) / 8.281_rp) !original O'Hara-Rudy formulation
                   vaux2 = exp((elmlo_exm(2) + 0.9941_rp) / 38.45_rp) !original O'Hara-Rudy formulation
                   xitaux = (0.02136_rp * vaux1) + (0.3052_rp * vaux2) !original O'Hara-Rudy formulation
                   tj = 2.038_rp + (1.0_rp / xitaux) !original O'Hara-Rudy formulation     
                end if
                vaux0 = vaulo_exm(4,2)  
                vaux3 = jss - (jss - vaux0) * exp(-dtimon/tj)    
                vaulo_exm(4,1) = vaux3      ! value of variable j       
           
           !!!!  hsp
                if(kfl_inaga_exm(mat) == 1_ip) then
                   vinf = 1.0_rp / (1.0_rp + exp((elmlo_exm(2) + 84.7_rp) / 6.22_rp)) !modified from Passini et al. (original:   1.0_rp / (1.0_rp + exp((elmlo_exm(2) + 89.1_rp) / 6.086_rp)))
                else
                   vinf = 1.0_rp / (1.0_rp + exp((elmlo_exm(2) + 89.1_rp) / 6.086_rp)) !original:   1.0_rp / (1.0_rp + exp((elmlo_exm(2) + 89.1_rp) / 6.086_rp)))
                end if
                xitaux = 3.0_rp * ths
                vaux0 = vaulo_exm(5,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(5,1) = vaux3      ! value of variable hsp           
                
           !!! jp  
                xitaux = 1.46_rp * tj           
                vaux0 = vaulo_exm(6,2)  
                vaux3 = jss - (jss - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(6,1) = vaux3      ! value of variable jp             
           
                !Calculate INA
                !!!!  hp
                hp = (ahf * vaulo_exm(2,1)) + (ahs * vaulo_exm(5,1)) 
                vaux1 = ahf * vaulo_exm(2,1)
                hh = vaux1 + (ahs * vaulo_exm(3,1))  !h
                finap = 1.0_rp / (1.0_rp + (kmcamk/viclo_exm(22,1)))
                vaux1 = gna*(elmlo_exm(2)-ena) 
                vaux1 = vaux1 * vaulo_exm(1,1) * vaulo_exm(1,1) * vaulo_exm(1,1)
                vaux2 = (1.0_rp-finap) * hh * vaulo_exm(4,1)
                vaux2 = vaux2 + (finap * hp * vaulo_exm(6,1)) 
                viclo_exm(1,1) = vaux1 * vaux2  !!!!  ina               
                      
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
                !%calculate inal
           !!!!!!  ml
                vinf = 1.0_rp + exp(-(elmlo_exm(2) + 42.85_rp) / 5.264_rp)
                vinf = 1.0_rp / vinf
                vaux0 = vaulo_exm(7,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/tm)    
                vaulo_exm(7,1) = vaux3      ! value of variable ml   
           
           !!!!! hl
                vinf = 1.0_rp + exp((elmlo_exm(2) + 87.61_rp) / 7.488_rp)
                vinf = 1.0_rp / vinf
                xitaux = 200.0_rp
                vaux0 = vaulo_exm(8,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(8,1) = vaux3      ! value of variable hl      
           
           !!!!  hlp
                vinf = 1.0_rp + exp((elmlo_exm(2) + 93.81_rp) / 7.488_rp)
                vinf = 1.0_rp / vinf
                xitaux = 3.0_rp * 200.0_rp          
                vaux0 = vaulo_exm(9,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(9,1) = vaux3      ! value of variable hlp                     

             !!!  Calculate inal current            
                finalp = 1.0_rp + (kmcamk/viclo_exm(22,1))
                finalp = 1.0_rp / finalp
                vaux1 = gnal*(elmlo_exm(2)-ena)*vaulo_exm(7,1)
                vaux2 = (1.0_rp-finalp)*vaulo_exm(8,1) + (finalp*vaulo_exm(9,1))
                viclo_exm(2,1)= vaux1 * vaux2
          
           !!!!! calculate ito
           !!!!!!  calculate variable a
                vinf = 1.0_rp + exp(-(elmlo_exm(2) - 14.34_rp) / 14.82_rp)
                vinf = 1.0_rp / vinf
                vaux1 =  1.0_rp / (1.2089_rp * (1.0_rp + exp(-(elmlo_exm(2) - 18.4099_rp) / 29.3814_rp)))
                vaux2 = 3.5_rp / (1.0_rp + exp((elmlo_exm(2) + 100.0_rp) / 29.3814_rp))
                ta =  1.0515_rp / (vaux1 + vaux2)
                vaux0 = vaulo_exm(10,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/ta)    
                vaulo_exm(10,1) = vaux3      ! value of variable a
           
           !!!!! calculate ifast  iss
                vinf = 1.0_rp + exp((elmlo_exm(2) + 43.94_rp) / 5.711_rp)
                iss = 1.0_rp / vinf
                if(ituss_exm == 3) then
                   delepi = (1.0_rp + exp((elmlo_exm(2) + 70.0_rp) / 5.0_rp))
                   delepi = (0.95_rp / delepi)
                   delepi = 1.0_rp - delepi
                else
                   delepi = 1.0_rp
                end if           
                vaux1 = 0.08004_rp * exp((elmlo_exm(2) + 50.0_rp) / 16.59_rp) !tif
                vaux2 = 0.3933_rp * exp(-(elmlo_exm(2) + 100.0_rp) / 100.0_rp)
                vaux3 = vaux1 + vaux2
                xitaux = 4.562_rp + (1.0_rp/ vaux3 )
                xitaux = delepi * xitaux
                tif = xitaux
                vaux0 = vaulo_exm(11,2)  
                vaux3 = iss - (iss - vaux0) * exp(-dtimon/tif)    
                vaulo_exm(11,1) = vaux3      ! value of variable ifast
           
           !!!!!!  calculate islow  tis
                vaux1 = 0.001416_rp * exp(-(elmlo_exm(2) + 96.52_rp) / 59.05_rp)
                vaux2 = 0.00000001780_rp * exp((elmlo_exm(2) + 114.1_rp) / 8.079_rp)
                vaux3 = vaux1 + vaux2
                xitaux = 23.62_rp + (1.0_rp / vaux3)
                xitaux = delepi * xitaux
                tis = xitaux         
                vaux0 = vaulo_exm(12,2)  
                vaux3 = iss - (iss - vaux0) * exp(-dtimon/tis)    
                vaulo_exm(12,1) = vaux3      ! value of variable islow    

           !!!!  calculate ap  (assp)
                vinf = 1.0_rp + exp(-(elmlo_exm(2)-24.34_rp)/14.82_rp)
                vinf = 1.0_rp / vinf
                vaux0 = vaulo_exm(13,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/ta)    
                vaulo_exm(13,1) = vaux3      ! value of variable ap   
           
           !!!!!!!!!! calculate ifp (dti_develop)
                vaux1 = exp((elmlo_exm(2)-167.4_rp)/15.89_rp)
                vaux2 = exp(-(elmlo_exm(2)-12.23_rp)/0.2154_rp)
                devel =  vaux1 + vaux2 
                devel = 1.354_rp + (0.0001_rp / devel)
           
                vaux1 = exp((elmlo_exm(2)+70.0_rp)/20.0_rp)
                vaux1 = 1.0_rp + vaux1
                recov = 1.0_rp - (0.5_rp / vaux1)
                xitaux = devel * recov * tif
                vaux0 = vaulo_exm(14,2)  
                vaux3 = iss - (iss - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(14,1) = vaux3      ! value of variable ifp
           
           !!!!! calculate isp
                xitaux = devel * recov * tis    
                vaux0 = vaulo_exm(15,2)  
                vaux3 = iss - (iss - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(15,1) = vaux3      ! value of variable isp      
           
                !!!!!!  Ito current
                !!!! calculate aif, ais and ii
                vaux1 = exp((elmlo_exm(2)-213.6_rp)/151.2_rp)
                aif = 1.0_rp / (1.0_rp + vaux1)
                ais = 1.0_rp - aif
                ii = aif*vaulo_exm(11,1) + (ais*vaulo_exm(12,1))
                ipp = aif*vaulo_exm(14,1) + (ais*vaulo_exm(15,1))
                vaux1 = kmcamk/viclo_exm(22,1)
                fitop = 1.0_rp/(1.0_rp + vaux1)
                vaux2 = (fitop*vaulo_exm(13,1) * ipp)
                vaux1 = ((1.0_rp-fitop) * vaulo_exm(10,1) * ii ) + vaux2
                vaux1 = gto * (elmlo_exm(2)-ek) * vaux1
                viclo_exm(3,1) = vaux1        !!! ito           
                    
           !!! calculate ical, icana, icak
          !%calculate ical
           !!!!  calculate gate d
                vaux1 = exp(-(elmlo_exm(2) + 3.940_rp) / 4.230_rp)
                vinf = 1.0_rp / (1.0_rp + vaux1)
                vaux1 =  exp(0.09_rp*(elmlo_exm(2) + 14.0_rp))
                vaux2 = exp(-0.05_rp*(elmlo_exm(2) + 6.0_rp))
                xitaux = vaux1 + vaux2
                xitaux = 0.6_rp + (1.0_rp/xitaux)
                vaux0 = vaulo_exm(16,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(16,1) = vaux3      ! value of gate d
           
           !!!! calculate gate ff
                vaux1 = exp((elmlo_exm(2) + 19.58_rp) / 3.696_rp)
                vinf = 1.0_rp/(1.0_rp + vaux1)
                fss = vinf
                vaux1 = 0.0045_rp * exp((elmlo_exm(2) + 20.0_rp) / 10.0_rp)
                vaux2 = 0.0045_rp * exp(-(elmlo_exm(2) + 20.0_rp) / 10.0_rp) 
                tff = 7.0_rp + (1.0_rp / (vaux1 + vaux2))
                vaux0 = vaulo_exm(17,2)  
                vaux3 = fss - (fss - vaux0) * exp(-dtimon/tff)    
                vaulo_exm(17,1) = vaux3      ! value of gate ff
           
           !!!!!!  calculate gate fs
                vaux1 =  0.000035_rp * exp((elmlo_exm(2)+5.0_rp)/6.0_rp)
                vaux2 = 0.000035_rp * exp(-(elmlo_exm(2)+5.0_rp)/4.0_rp) 
                xitaux = 1000.0_rp + (1.0_rp / (vaux1 + vaux2))          
                vaux0 = vaulo_exm(18,2)  
                vaux3 = fss - (fss - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(18,1) = vaux3      ! value of variable fs
           
           !!!!! calculate fcass
                !!fcass=fss
                vaux1 = (0.04_rp*exp((elmlo_exm(2)-4.0_rp)/7.0_rp))
                vaux2 = (0.04_rp*exp(-(elmlo_exm(2)-4.0_rp)/7.0_rp))  
                tfcaf = 7.0_rp + (1.0_rp / (vaux1 + vaux2))           
                vaux0 = vaulo_exm(19,2)  
                vaux3 = fss - (fss - vaux0) * exp(-dtimon/tfcaf)    
                vaulo_exm(19,1) = vaux3      ! value of variable fcaf   
           
           !!!! calculate gate fcas
                vaux1 = (0.00012_rp * exp(elmlo_exm(2)/7.0_rp))
                vaux2 = (0.00012_rp * exp(-elmlo_exm(2)/3.0_rp))  
                xitaux = 100.0_rp + (1.0_rp / (vaux1 + vaux2))
                vaux0 = vaulo_exm(20,2)  
                vaux3 = fss - (fss - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(20,1) = vaux3      ! value of variable fcas                       

           !!! calculate gate jca
                xitaux = 75.0_rp           
                vaux0 = vaulo_exm(21,2)  
                vaux3 = fss - (fss - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(21,1) = vaux3      ! value of variable jca  
           
           !!! calculate gate ffp 
                xitaux = 2.5_rp * tff     
                vaux0 = vaulo_exm(23,2)  
                vaux3 = fss - (fss - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(23,1) = vaux3      ! value of variable ffp   
                
           !!!! calculate gate fcafp
                xitaux = 2.5_rp * tfcaf
                vaux0 = vaulo_exm(24,2)  
                vaux3 = fss - (fss - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(24,1) = vaux3      ! value of variable ffp   
              
           !!!!! calculate icana and icak
                km2n = vaulo_exm(21,1) * 1.0_rp
                vaux1 = kmn / vcolo_exm(6,2)
                vaux4 = (1.0_rp+vaux1) *(1.0_rp+vaux1) * (1.0_rp+vaux1) * (1.0_rp+vaux1)
                vaux2 = (k2n / km2n) +vaux4
                anca = 1.0_rp / vaux2
                dnca = (anca*k2n/km2n)               
                vaux0 = vaulo_exm(22,2)                  
                vaux3 = dnca - (dnca - vaux0) * exp(-dtimon*km2n)    
                vaulo_exm(22,1) = vaux3

           !!!!!! calculate ff  = F
                ff = (aff * vaulo_exm(17,1)) + (afs * vaulo_exm(18,1))   !value of variable ff
               !!! calculate fca
                vaux1 = 1.0_rp + exp((elmlo_exm(2)-10.0_rp)/10.0_rp)
                afcaf = 0.3_rp + (0.6_rp/vaux1)
                afcas = 1.0_rp - afcaf    
                fca = (afcaf*vaulo_exm(19,1)) + (afcas*vaulo_exm(20,1))
           !!!!  calculate fcap
                fcap = afcaf*vaulo_exm(24,1) + (afcas*vaulo_exm(20,1))
           !!!!! calculate fp
                fp = aff*vaulo_exm(23,1) + (afs*vaulo_exm(18,1))                
                vaux2 = exp(2.0_rp*vfrt)
                vaux1 = exp(1.0_rp*vfrt)
                vaux3 = 1.0_rp / (vaux2-1.0_rp)
                phical = 4.0_rp * vffrt*((vcolo_exm(6,2)*vaux2) - (0.341_rp*cao)) * vaux3
                vaux3 = 1.0_rp / (vaux1 -1.0_rp)
                phicana = 1.0_rp * vffrt*((0.75_rp*vcolo_exm(2,2)*vaux1) - (0.75_rp*nao)) * vaux3
                phicak = 1.0_rp * vffrt*((0.75_rp*vcolo_exm(4,2)*vaux1) - (0.75_rp*ko)) * vaux3
                ficalp = 1.0_rp / (1.0_rp + (kmcamk/viclo_exm(22,1)))
           
           !!!! calculate ical current
                vaux1 = (fp * (1.0_rp-vaulo_exm(22,1)) + (vaulo_exm(21,1)*fcap*vaulo_exm(22,1)))
                vaux1 = ficalp * pcap * phical * vaulo_exm(16,1) * vaux1
                vaux2 = ff * (1.0_rp-vaulo_exm(22,1)) + (vaulo_exm(21,1) * fca * vaulo_exm(22,1))
                vaux3 = (1.0_rp-ficalp) * pca * phical * vaulo_exm(16,1)
                viclo_exm(4,1)= (vaux3 * vaux2) + vaux1  !!! ical
           
           !!!! calculate icana current
                vaux1 = (fp*(1.0_rp-vaulo_exm(22,1)) + (vaulo_exm(21,1)*fcap*vaulo_exm(22,1)))
                vaux1 = vaux1 * ficalp*pcanap*phicana*vaulo_exm(16,1)
                vaux2 = (ff*(1.0_rp-vaulo_exm(22,1))) + (vaulo_exm(21,1)*fca*vaulo_exm(22,1))
                vaux2 = vaux2 * (1.0_rp-ficalp)*pcana*phicana*vaulo_exm(16,1)
                viclo_exm(24,1)= vaux2 + vaux1  !! icana
           
           !!!!! calculate icak
                vaux1 = fp * (1.0_rp-vaulo_exm(22,1)) + (vaulo_exm(21,1) * fcap * vaulo_exm(22,1))
                vaux1 = vaux1 * ficalp * pcakp * phicak * vaulo_exm(16,1)
                vaux2 = ff * (1.0_rp-vaulo_exm(22,1)) + (vaulo_exm(21,1) * fca * vaulo_exm(22,1))
                vaux2 = vaux2 * (1.0_rp-ficalp) * pcak * phicak * vaulo_exm(16,1)
                viclo_exm(25,1) = vaux2 + vaux1   !!! icak
                 
           !!!!  calculate ikr GATES  XRF
                vinf = exp(-(elmlo_exm(2) + 8.337_rp) / 6.789_rp)
                vinf = 1.0_rp / (1.0_rp + vinf)
                vaux1 = 0.00004123_rp * exp(-(elmlo_exm(2) - 47.78_rp) / 20.38_rp)
                vaux2 = 0.3652_rp * exp((elmlo_exm(2) - 31.66_rp) / 3.869_rp) 
                xitaux = 12.98_rp + (1.0_rp / (vaux1 + vaux2))          
                vaux0 = vaulo_exm(25,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(25,1) = vaux3      ! value of variable xrf 
           !! GATE XRS
                vaux1 = 0.06629_rp * exp((elmlo_exm(2) - 34.70_rp) / 7.355_rp)
                vaux2 = 0.00001128_rp * exp(-(elmlo_exm(2) - 29.74_rp) / 25.94_rp)
                xitaux = 1.865_rp + (1.0_rp / (vaux1 + vaux2))          
                vaux0 = vaulo_exm(26,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(26,1) = vaux3      ! value of variable xrs
 
                !%calculate ikr CURRENT               
                vaux1 = exp((elmlo_exm(2)+54.81_rp)/38.21_rp)
                axrf = 1.0_rp / (1.0_rp + vaux1)
                axrs = 1.0_rp - axrf    
                xr = axrf*vaulo_exm(25,1) + (axrs*vaulo_exm(26,1))
                vaux1 = 1.0_rp + exp((elmlo_exm(2)+55.0_rp)/75.0_rp)
                vaux2 = 1.0_rp + exp((elmlo_exm(2)-10.0_rp)/30.0_rp)
                vaux3 = vaux1*vaux2                
                rkr = 1.0_rp / vaux3
                vaux1 = sqrt(ko/5.4_rp) * (elmlo_exm(2)-ek)
                viclo_exm(5,1) = gkr * vaux1 * xr * rkr   !!! ikr 
          
           !!! %calculate iks GATES
                vinf =  exp(-(elmlo_exm(2) + 11.60_rp)/8.932_rp)
                vinf = 1.0_rp / (1.0_rp + vinf)
                vaux1 = 0.0002326_rp * exp((elmlo_exm(2) + 48.28_rp)/17.80_rp)
                vaux2 = 0.001292_rp * exp(-(elmlo_exm(2) + 210.0_rp)/230.0_rp)
                xitaux = 817.3_rp + (1.0_rp/(vaux1 + vaux2))          
                vaux0 = vaulo_exm(27,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(27,1) = vaux3      ! value of variable xs1
           
                !vinf = xs2ss = xs1ss
                vaux1 = 0.01_rp * exp((elmlo_exm(2) - 50.0_rp) / 20.0_rp)
                vaux2 =  0.0193_rp * exp(-(elmlo_exm(2) + 66.54_rp) / 31.0_rp)
                xitaux = 1.0_rp / (vaux1 + vaux2) 
                vaux0 = vaulo_exm(28,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(28,1) = vaux3      ! value of variable xs2
               
                !!!! calculate iks CURRENTS
                ksca = (0.000038_rp/vcolo_exm(1,2)) ** (1.4_rp) !!(7.0_rp/5.0_rp)
                ksca = 1.0_rp + (0.6_rp / (1.0_rp + ksca))
                vaux1 = gks * ksca * vaulo_exm(27,1) * vaulo_exm(28,1)
                viclo_exm(6,1) =  vaux1 * (elmlo_exm(2)-eks)
            
           !!!!  calculate gate xk1 GATES
                vaux1 = 1.0_rp / (1.5692_rp * ko + 3.8115_rp)
                vaux2 = elmlo_exm(2) + 2.5538_rp * ko + 144.59_rp
                vinf = exp(-vaux2 * vaux1 )
                vinf = 1.0_rp / (1.0_rp + vinf)
                vaux1 = exp(-(elmlo_exm(2) + 127.2_rp) / 20.36_rp)
                vaux2 = exp((elmlo_exm(2) + 236.8_rp) / 69.33_rp)
                xitaux = 122.2_rp / (vaux1 + vaux2)
                vaux0 = vaulo_exm(29,2)  
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vaulo_exm(29,1) = vaux3      ! value of variable xk1    
 
                !!!! calculate ik1 CURRENT
                rk1 = 1.0_rp + exp((elmlo_exm(2) + 105.8_rp - 2.6_rp * ko) / 9.493_rp)
                rk1 = 1.0_rp / rk1
                vaux1 = gk1 * sqrt(ko) * rk1 * vaulo_exm(29,1) 
                viclo_exm(7,1) = vaux1 * (elmlo_exm(2)-ek) !!! ik1 
  
                !%calculate INaCa CURRENT	
                hca = exp(qca * vfrt)
                hna = exp(qna * vfrt)
                h1 = 1.0_rp + (vcolo_exm(5,2)/kna3) * (1.0_rp +hna)
                h2 = (vcolo_exm(5,2)*hna) / (kna3*h1)
                h3 = 1.0_rp/h1
                h4 = 1.0_rp + vcolo_exm(5,2) / kna1 * (1.0_rp + (vcolo_exm(5,2)/kna2))
                h5 = vcolo_exm(5,2)*vcolo_exm(5,2) / (h4*kna1*kna2)
                h6 = 1.0_rp/h4
                h7 = 1.0_rp + nao/kna3*(1.0_rp + (1.0_rp/hna))
                h8 = nao / (kna3*hna*h7)
                h9 = 1.0_rp / h7
                h10 = kasymm + 1.0_rp + nao/kna1 * (1.0_rp + (nao/kna2))
                h11 = nao*nao / (h10*kna1*kna2)
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
                k6 = h6*vcolo_exm(1,2)*kcaon
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
                allo = (kmcaact/vcolo_exm(1,2)) * (kmcaact/vcolo_exm(1,2))
                allo = 1.0_rp / (1.0_rp + allo)
                zna = 1.0_rp
                jncxna = (3.0_rp*(e4*k7-e1*k8)) + (e3*k4pp) - (e2*k3pp)
                jncxca = (e2*k2) - (e1*k1)
                vaux1 = (zna*jncxna) + (zca*jncxca)
                viclo_exm(8,1) = 0.8_rp * gncx * allo * vaux1  !!! inaca_i Current
  
                !%calculate inaca_ss
                h1 = 1.0_rp + vcolo_exm(2,2)/kna3 * (1.0_rp + hna)
                h2 = (vcolo_exm(2,2)*hna) / (kna3*h1)
                h3 = 1.0_rp/h1
                h4 = 1.0_rp + vcolo_exm(2,2)/ kna1 * (1.0_rp + (vcolo_exm(2,2)/kna2))
                h5 = vcolo_exm(2,2)*vcolo_exm(2,2) / (h4*kna1*kna2)
                h6 = 1.0_rp / h4
                h7 = 1.0_rp + (nao/kna3) *(1.0_rp + (1.0_rp/hna))
                h8 = nao / (kna3*hna*h7)
                h9 = 1.0_rp / h7
                h10 = kasymm + 1.0_rp + nao/kna1 * (1.0_rp + (nao/kna2))
                h11 = nao*nao / (h10*kna1*kna2)
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
                k6 = h6*vcolo_exm(6,2)*kcaon
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
                allo = (kmcaact/vcolo_exm(6,2)) * (kmcaact/vcolo_exm(6,2))
                allo = 1.0_rp / (1.0_rp + allo)
                jncxna = (3.0_rp*(e4*k7-e1*k8)) + (e3*k4pp) - (e2*k3pp)
                jncxca = (e2*k2) - (e1*k1)
                viclo_exm(9,1) = 0.2_rp * gncx * allo * ((zna*jncxna) + (zca*jncxca))   !!! inaca_ss
             
                INaCa = viclo_exm(8,1)+viclo_exm(9,1)
  
                !%calculate inak
                k3p = 1899.0_rp
                k4p = 639.0_rp
                knai = knai0 * exp(delta*vfrt/3.0_rp)
                knao = knao0 * exp((1.0_rp-delta)*vfrt/3.0_rp)
                pp = ep / (1.0_rp + hbig/khp + (vcolo_exm(5,2)/knap) + (vcolo_exm(3,2)/kxkur))
                vaux1 = (vcolo_exm(5,2)/knai) * (vcolo_exm(5,2)/knai) * (vcolo_exm(5,2)/knai)
                vaux2 = (1.0_rp + (vcolo_exm(5,2)/knai)) * (1.0_rp + (vcolo_exm(5,2)/knai)) * (1.0_rp + (vcolo_exm(5,2)/knai))
                vaux3 = (1.0_rp + (vcolo_exm(3,2)/kki)) * (1.0_rp + (vcolo_exm(3,2)/kki))
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
                vaux1 = (vcolo_exm(3,2)/kki) * (vcolo_exm(3,2)/kki)
                vaux2 = (1.0_rp + (vcolo_exm(5,2)/knai)) * (1.0_rp + (vcolo_exm(5,2)/knai)) * (1.0_rp + (vcolo_exm(5,2)/knai))
                vaux3 = (1.0_rp + (vcolo_exm(3,2)/kki)) * (1.0_rp + (vcolo_exm(3,2)/kki))
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
                viclo_exm(10,1) = pnak*(zna*jnakna + zk*jnakk)  !!!! inak
           
                !%calculate ikb CURRENT
                xkb = 1.0_rp / (1.0_rp + exp(-(elmlo_exm(2)-14.48_rp)/18.34_rp))
                viclo_exm(11,1) = gkb*xkb*(elmlo_exm(2)-ek)  !!!! ikb
                
                !%calculate inab CURRENT
                vaux1 = pnab*vffrt*(vcolo_exm(5,2)*exp(vfrt) - nao)
                viclo_exm(12,1) =  vaux1 / (exp(vfrt)-1.0_rp)
           
                !%calculate icab CURRENT
                vaux1 = (vcolo_exm(1,2)*exp(2.0_rp*vfrt) - (0.341_rp*cao))/(exp(2.0_rp*vfrt)-1.0_rp)
                viclo_exm(13,1) = pcab * 4.0_rp*  vffrt * vaux1   !!!!icab
           
                !%calculate ipca CURRENT
                viclo_exm(14,1) = gpca * vcolo_exm(1,2) / (0.0005_rp+vcolo_exm(1,2))  !!! ipca
              
              
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              !!  VOLTAGE INTEGRATION    
                v1 = 0.0_rp
                v2 = 0.0_rp
                do i=1,14
                    v1= v1 + viclo_exm(i,1)
                end do
                do i=23,25
                    v2= v2 + viclo_exm(i,1)
                end do
                elmlo_exm(1) =  -(v1+v2) ! *farad/acap  
                elmlo_exm(1) = elmlo_exm(2) + dtimon*elmlo_exm(1)
            !!  Runge Kutta
                !k1 = elmlo_exm(1) 
                !k2 = elmlo_exm(1) + 0.5_rp * dtimon * k1
                !k3 = elmlo_exm(1) + 0.5_rp * dtimon * k2 
                !k4 = elmlo_exm(1) + dtimon * k3 
                !elmlo_exm(1) = elmlo_exm(2) + (dtimon/6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
                elmlo_exm(3) = elmlo_exm(2)
                elmlo_exm(2) = elmlo_exm(1)  
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

           !!! Update CaMKa
                vaux2 = 1.0_rp / (1.0_rp + (kmcam/vcolo_exm(6,2)))
                viclo_exm(26,1) = vaux2 * camko* (1.0_rp - vcolo_exm(11,2))
                viclo_exm(22,1) = viclo_exm(26,1) + vcolo_exm(11,2)  ! camka
                 
                !%update camk
                rhsx1 = acamk * viclo_exm(26,1) * (viclo_exm(26,1) + vcolo_exm(11,2))
                rhsx =  rhsx1 - (bcamk*vcolo_exm(11,2))
                val0  =  vcolo_exm(11,3)
                camkt = val0 + dtimon*rhsx
                !k1 = rhsx 
                !k2 = rhsx + 0.5_rp * dtimon* k1
                !k3 = rhsx + 0.5_rp * dtimon * k2 
                !k4 = rhsx + dtimon * k3 
                !camkt = val0 + ((dtimon / (6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
                vcolo_exm(11,1) = camkt      ! value of camkt concentration

                !%calculate diffusion fluxes
                viclo_exm(16,1)=(vcolo_exm(2,2)-vcolo_exm(5,2))/2.0_rp    !!! jdiffna
                viclo_exm(17,1)=(vcolo_exm(4,2)-vcolo_exm(3,2))/2.0_rp    !!! jdiffk
                viclo_exm(15,1)=(vcolo_exm(6,2)-vcolo_exm(1,2))/0.2_rp    !!! jdiff
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
                
                !%calculate ryanodione receptor calcium induced calcium release from the jsr
                bt = 4.75_rp
                a_rel = 0.5_rp * bt
                vaux1 = (1.5_rp/vcolo_exm(8,2)) ** 8.0_rp 
                vinf = a_rel * (-viclo_exm(4,1)) / (1.0_rp + vaux1)
                if(ituss_exm == 2) then
                    vinf = vinf * 1.7_rp
                end if
                xitaux = bt / (1.0_rp + 0.0123_rp/vcolo_exm(8,2))
                
                if(xitaux < 0.005_rp) then
                   xitaux = 0.005_rp
                end if
                vaux0 = vcolo_exm(9,3)
                !vaux1 = (vinf - vaux0) * xitaux;
                !vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
                !vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
                !vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
                !vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
                !vaux3 = vaux5
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vcolo_exm(9,1) = vaux3      ! value of variable jrelnp   

                !!  JRELP                
                btp = 1.25_rp * bt
                a_relp = 0.5_rp * btp
                vaux1 = (1.5_rp/vcolo_exm(8,2)) ** 8.0_rp 
                vinf = a_relp * (-viclo_exm(4,1)) / (1.0_rp + vaux1)
                if(ituss_exm == 2) then
                    vinf = vinf * 1.7_rp
                end if
                xitaux = btp / (1.0_rp+ (0.0123_rp/vcolo_exm(8,2)))
                
                if(xitaux < 0.005_rp) then
                   xitaux = 0.005_rp 
                end if
                vaux0 = vcolo_exm(10,3)
                !vaux1 = (vinf - vaux0) * xitaux;
                !vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
                !vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
                !vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
                !vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
                !vaux3 = vaux5
                vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/xitaux)    
                vcolo_exm(10,1) = vaux3      ! value of jrelp
           
                fjrelp= 1.0_rp / (1.0_rp + (kmcamk/viclo_exm(22,1)))
                viclo_exm(21,1)= (1.0_rp - fjrelp) * vcolo_exm(9,1) + (fjrelp * vcolo_exm(10,1)) !!! jrel
             
                !%calculate serca pump, ca uptake flux
                jupnp = (0.004375_rp * vcolo_exm(1,2) / (vcolo_exm(1,2) + 0.00092_rp)) * ttparmate_exm(3,12,mat)
                vaux1 = 1.0_rp / (vcolo_exm(1,2) + 0.00092_rp - 0.00017_rp)
                jupp = 2.75_rp * 0.004375_rp * vcolo_exm(1,2) * vaux1 * ttparmate_exm(3,12,mat)
                if(ituss_exm == 3) then
                   jupnp = jupnp * 1.3_rp * ttparmate_exm(3,12,mat)
                   jupp = jupp * 1.3_rp * ttparmate_exm(3,12,mat)
                end if
                fjupp = (1.0_rp / (1.0_rp + kmcamk/viclo_exm(22,1)))
                viclo_exm(19,1) = 0.0039375_rp*vcolo_exm(7,2)/15.0_rp   !!!! jleak
                !viclo_exm(19,3) = viclo_exm(19,2)
                viclo_exm(18,1) = (1.0_rp-fjupp) * jupnp + fjupp * jupp - viclo_exm(19,1)  !!! jup
           
                !%calculate tranlocation flux
                viclo_exm(20,1) = (vcolo_exm(7,2)-vcolo_exm(8,2)) / 100.0_rp  !!!! jtr
               
                        !!!!!  CONCENTRATIONS!!!!   
                
                !%update intracellular concentrations, using buffers for cai, cass, cajsr
                !! calculate na current
                vaux1 =  viclo_exm(1,1) + viclo_exm(2,1) + viclo_exm(12,1)
                vaux2 = 3.0_rp*viclo_exm(8,1) + 3.0_rp*viclo_exm(10,1)
                vaux3 = -(vaux1+vaux2) * acap/(farad*vmyo)
                rhsx = vaux3 + (viclo_exm(16,1)*vss/vmyo)
                val0 = vcolo_exm(5,3)
                nai = val0 + dtimon*rhsx
                !k1 = rhsx 
                !k2 = rhsx + 0.5_rp * dtimon* k1
                !k3 = rhsx + 0.5_rp * dtimon * k2 
                !k4 = rhsx + dtimon * k3 
                !nai = val0 + (dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
                vcolo_exm(5,1) = nai
                
                !!! calculate na current in subspace ss
                vaux1 = (viclo_exm(24,1) + 3.0_rp * viclo_exm(9,1)) * acap / (farad*vss)
                rhsx =  -vaux1 - viclo_exm(16,1)
                val0 = vcolo_exm(2,3)
                nass = val0 + dtimon*rhsx
                !k1 = rhsx 
                !k2 = rhsx + 0.5_rp * dtimon* k1
                !k3 = rhsx + 0.5_rp * dtimon * k2 
                !k4 = rhsx + dtimon * k3 
                !nass = val0 + (dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
                vcolo_exm(2,1) = nass             !!! nass
             
                 !!! calculate k current
                vaux1 = viclo_exm(3,1) + viclo_exm(5,1) + viclo_exm(6,1) + viclo_exm(23,1)
                vaux2 = viclo_exm(7,1) + viclo_exm(11,1)  - (2.0_rp*viclo_exm(10,1)) 
                vaux3 = (viclo_exm(17,1)*vss/vmyo)
                rhsx = -((vaux1+ vaux2) * acap/(farad*vmyo)) + vaux3
                val0 = vcolo_exm(3,3)
                ki = val0 + dtimon*rhsx
                !k1 = rhsx 
                !k2 = rhsx + 0.5_rp * dtimon* k1
                !k3 = rhsx + 0.5_rp * dtimon * k2 
                !k4 = rhsx + dtimon * k3 
                !ki = val0 + (dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4) 
                vcolo_exm(3,1) = ki             !!! ki  
                
                !!!!  calculate k current in the subspace ss               
                rhsx = -(viclo_exm(25,1)*acap/(farad*vss)) - viclo_exm(17,1)
                val0 = vcolo_exm(4,3)
                kss = val0 + dtimon*rhsx
                !k1 = rhsx 
                !k2 = rhsx + 0.5_rp * dtimon* k1
                !k3 = rhsx + 0.5_rp * dtimon * k2 
                !k4 = rhsx + dtimon * k3 
                !kss = val0 + (dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
                vcolo_exm(4,1) = kss             !!! kss
                
                !!!  calculate ca current (Cai)
                vaux1 = (kmcmdn+vcolo_exm(1,2)) * (kmcmdn+vcolo_exm(1,2))
                vaux2 = (kmtrpn+vcolo_exm(1,2)) * (kmtrpn+vcolo_exm(1,2))
                vaux3 = 1.0_rp/(1.0_rp + (cmdnmax*kmcmdn/vaux1) + (trpnmax*kmtrpn/vaux2))
                rhsx1 = viclo_exm(14,1)+viclo_exm(13,1) - (2.0_rp*viclo_exm(8,1))
                rhsx2 = -(viclo_exm(18,1)*vnsr/vmyo) + (viclo_exm(15,1)*vss/vmyo)
                rhsx = vaux3 *(-(rhsx1 *acap/(2.0_rp*farad*vmyo)) + rhsx2 )
                val0 = vcolo_exm(1,3)
                cai = val0 + dtimon*rhsx
                !k1 = rhsx 
                !k2 = rhsx + 0.5_rp * dtimon* k1
                !k3 = rhsx + 0.5_rp * dtimon * k2 
                !k4 = rhsx + dtimon * k3 
                !cai = val0 + (dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
                vcolo_exm(1,1) = cai             !!! cai
            
                !!!! calculate ca current in the subspace ss
                vaux1 = (kmbsr+vcolo_exm(6,2)) * (kmbsr+vcolo_exm(6,2))
                vaux2 = (kmbsl+vcolo_exm(6,2)) * (kmbsl+vcolo_exm(6,2))
                vaux3 = 1.0_rp / (1.0_rp + (bsrmax*kmbsr/vaux1) + (bslmax*kmbsl/vaux2)) !Bcass
                rhsx1 = viclo_exm(4,1) - (2.0_rp*viclo_exm(9,1))
                rhsx2 = viclo_exm(21,1)*vjsr/vss
                rhsx = vaux3 *(-rhsx1 *acap/(2.0_rp*farad*vss) + rhsx2 - viclo_exm(15,1))
                val0 = vcolo_exm(6,3)
                cass = val0 + dtimon*rhsx
                !k1 = rhsx 
                !k2 = rhsx + 0.5_rp * dtimon* k1
                !k3 = rhsx + 0.5_rp * dtimon * k2 
                !k4 = rhsx + dtimon * k3 
                !cass = val0 + (dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
                vcolo_exm(6,1) = cass             !!! cass
                
                !!!! calculate ca current in the sarcoplasmic reticulum nsr
                rhsx1 = viclo_exm(20,1)*vjsr / vnsr
                rhsx = viclo_exm(18,1) - rhsx1
                val0 = vcolo_exm(7,3)
                cansr = val0 + dtimon*rhsx
                !k1 = rhsx 
                !k2 = rhsx + 0.5_rp * dtimon* k1
                !k3 = rhsx + 0.5_rp * dtimon * k2 
                !k4 = rhsx + dtimon * k3 
                !cansr = val0 + (dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
                vcolo_exm(7,1) = cansr             !!! cansr
                
                !!!! calculate ca current in the junctional sarcoplasmic reticulum jsr
                vaux1 = (kmcsqn+vcolo_exm(8,2)) * (kmcsqn+vcolo_exm(8,2))
                vaux3 = csqnmax * kmcsqn / vaux1
                vaux2 = 1.0_rp / (1.0_rp + vaux3) 
                rhsx = vaux2 *(viclo_exm(20,1)-viclo_exm(21,1))
                val0 = vcolo_exm(8,3)
                cajsr = val0 + dtimon*rhsx
                !k1 = rhsx 
                !k2 = rhsx + 0.5_rp * dtimon* k1
                !k3 = rhsx + 0.5_rp * dtimon * k2 
                !k4 = rhsx + dtimon * k3 
                !cajsr = val0 + ((dtimon / (6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
                vcolo_exm(8,1) = cajsr             !!! cajsr
                
                vcolo_exm(:,2)=vcolo_exm(:,1)
                vcolo_exm(:,3)=vcolo_exm(:,2)                 
                vaulo_exm(:,2)=vaulo_exm(:,1)                 
                !vaulo_exm(:,2)=vaulo_exm(:,1)
            
                itim = itim + 1
                if (j>=200) then  !it was 10 for saving
                 !   write(995,400) vaulo_exm(1:29,1)
                    write(996,*) vcolo_exm(1,2)  !atim
                 !   write(994,300) viclo_exm(1:26,1)
                    write(997,*) elmlo_exm(2)
                    j = 0
                end if
                    j = j + 1
                !viclo_exm(:,2)=viclo_exm(:,1)

                
                ss_voltages_beat_current(ss_voltage_current_id) = elmlo_exm(2)
                ss_voltage_current_id = ss_voltage_current_id+1_ip

        end do


        batim = itim * dtimon

        

        !--------------------------------------------
        !
        !  begin: check if we reached steady state. If so, break      
        !
        !----------------------------------------
        if (ss_voltage_current_id-1 .NE. nsamples_per_beat) then
            print *, 'exm_oceohr.f90: The beat length was estimated incorrectly'
        end if


        ss_min_voltage = minval(ss_voltages_beat_current)
        !print *, 'Beat ',itim/nsamples_per_beat,': min voltage: ',ss_min_voltage
        if ( ss_min_voltage .LE. ss_min_voltage_threshold ) then

            ss_rmse = sum( (ss_voltages_beat_current - ss_voltages_beat_previous)**2/nsamples_per_beat )
            !print *, 'Beat ',itim/nsamples_per_beat,': rmse between beats: ',ss_rmse
            if ( ss_rmse < ss_voltages_epsilon ) then
                ss_rmse_counter = ss_rmse_counter+1
            else
                ss_rmse_counter = 0_ip
            end if

            if (ss_rmse_counter .GE. ss_rmse_counter_max ) then
                exm_oceohr = .TRUE.
                exit
            end if

        end if
        
        !otherwise copy the beat and go to the beginning
        ss_voltages_beat_previous = ss_voltages_beat_current
        !--------------------------------------------
        !
        !  end: check if we reached steady state        
        !
        !----------------------------------------



    end do
    !write(993,*) vcolo_exm(:,2)  !JAS
    !write(994,*) vaulo_exm(:,2)
     !write(995,*) elmlo_exm(2)
  
    write(997,*) "999999"  !separator between cells

    deallocate ( ss_voltages_beat_current )   
    deallocate ( ss_voltages_beat_previous ) 

        if(ituss_exm==3) then  !epi
            vauin_exm(:,3,mat)=vaulo_exm(:,1)
            vcoin_exm(:,3,mat) = vcolo_exm(:,1)
            vminimate_exm(3,mat) = elmlo_exm(2)
            !vauin_exm(:,1,3)=vaulo_exm(:,1)
            !vcoin_exm(:,1,3) = vcolo_exm(:,1)
            !vminimate_exm(3) = elmlo_exm(1)
        else if (ituss_exm==1) then  !endo
            vauin_exm(:,1,mat)=vaulo_exm(:,1)
            vcoin_exm(:,1,mat) = vcolo_exm(:,1)
            vminimate_exm(1,mat) = elmlo_exm(2)
            !vauin_exm(:,2,3)=vaulo_exm(:,1)
            !vcoin_exm(:,2,3) = vcolo_exm(:,1)
            !vminimate_exm(3) = elmlo_exm(1)
        else if (ituss_exm==2) then  !mid
            vauin_exm(:,2,mat)=vaulo_exm(:,1)
            vcoin_exm(:,2,mat) = vcolo_exm(:,1)
            vminimate_exm(2,mat) = elmlo_exm(2)
            !vauin_exm(:,3,3)=vaulo_exm(:,1)
            !vcoin_exm(:,3,3) = vcolo_exm(:,1)
            !vminimate_exm(1,3) = elmlo_exm(1)
        end if
  
  200 format(11(f14.8))
  300 format(26(f14.8))
  400 format(29(f14.8))
  !end if

end function exm_oceohr
