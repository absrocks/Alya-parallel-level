!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_iniohr.f90
!> @author  Jazmin Aguado-Sierra
!> @brief   Initial condition setup for Ohara-Rudy heterogeneous model
!> @details Obtains initial conditions of Normal or Heart Failure cell at  70 bpm or 857 ms \n
!!   C++ Implementation of the O'Hara-Rudy dynamic (ORd) model for the \n
!!   undiseased human ventricular action potential and calcium transient \n
!!  \n
!!   The ORd model is described in the article "Simulation of the Undiseased \n
!!   Human Cardiac Ventricular Action Potential: Model Formulation and \n
!!   Experimental Validation" \n
!!   by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy \n
!!   \n
!!   The article and supplemental materails are freely available in the \n
!!   Open Access jounal PLoS Computational Biology \n
!!   Link to Article: \n
!!   http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061 \n
!!    \n
!!   Email: tom.ohara@gmail.com / rudy@wustl.edu \n
!!   Web: http://rudylab.wustl.edu \n
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_iniohr(ipoin,mat)
!subroutine exm_iniohr(kmodel_ipoin,ipoin,mat)

  use      def_master
  use      def_domain
  use      def_elmtyp
  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: ipoin, mat
  integer(ip)   :: ituss_exm
  !ituss = 0_ip  !!%endo = 1, epi = 0, M = 2  


  if(kfl_timei_exm==1_ip) then 
     !mat=kmodel_ipoin-3    
     ituss_exm = int(celty_exm(1,ipoin))

     if(ituss_exm == EXM_CELLTYPE_EPI) then   !epicardial
        elmag(ipoin,1:3) = vminimate_exm(3,mat)     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        elmlo_exm(1:2) = vminimate_exm(3,mat)
        vicel_exm(1:26,ipoin,1) = 0.0_rp
        vconc(1:11,ipoin,1) = vcoin_exm(1:11,3,mat) !8.17202753912090_rp      ! nai=2
        vconc(1:11,ipoin,2) = vcoin_exm(1:11,3,mat) !8.17202753912090_rp      ! nai=2
        vconc(1:11,ipoin,3) = vcoin_exm(1:11,3,mat) !8.17202753912090_rp      ! nai=2
        vauxi_exm(1:29,ipoin,1) = vauin_exm(1:29,3,mat) 
        vauxi_exm(1:29,ipoin,2) = vauin_exm(1:29,3,mat) 
        vauxi_exm(1:29,ipoin,3) = vauin_exm(1:29,3,mat) 
        !           elmag(1,ipoin,1:4) = vminimate_exm(3,mat) * 0.000003335640952_rp    ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        !            vconc(2,ipoin,1:3) = 8.17212071948942_rp       !nass=3);
        !            vconc(3,ipoin,1:3) = 143.675184333376_rp       !ki=4);
        !            vconc(4,ipoin,1:3) = 143.675147146842_rp       !kss=5);
        !            vconc(5,ipoin,1:3) = 8.03149767106260e-05_rp       !cai=6);
        !            vconc(6,ipoin,1:3) = 7.96912949697674e-05_rp       !cass=7);
        !            vconc(7,ipoin,1:3) = 2.14811455007091_rp       !cansr=8);
        !            vconc(8,ipoin,1:3) = 2.03335899236798_rp      !cajsr=9);
        !            vauxi_exm(1,ipoin,1:3) = 0.00739719746920272_rp       !m=10);
        !            vauxi_exm(2,ipoin,1:3) = 0.695621622011335_rp       !hf=11);
        !            vauxi_exm(3,ipoin,1:3) = 0.695601842634086_rp       !hs=12);
        !            vauxi_exm(4,ipoin,1:3) = 0.695486248719023_rp       !j=13);
        !            vauxi_exm(5,ipoin,1:3) = 0.452023628358454_rp       !hsp=14);
        !            vauxi_exm(6,ipoin,1:3) = 0.695403157533235_rp       !jp=15);
        !            vauxi_exm(7,ipoin,1:3) = 0.000190839777466418_rp       !mL=16);
        !            vauxi_exm(8,ipoin,1:3) = 0.493606704642336_rp       !hL=17);
        !            vauxi_exm(9,ipoin,1:3) = 0.264304293390731_rp       !hLp=18);
        !            vauxi_exm(10,ipoin,1:3) = 0.00100594231451985_rp      !a=19);
        !            vauxi_exm(11,ipoin,1:3) = 0.999548606668578_rp      !iF=20);
        !            vauxi_exm(12,ipoin,1:3) = 0.999488774162635_rp      !iS=21);
        !            vauxi_exm(13,ipoin,1:3) = 0.000512555980943569_rp      !ap=22);
        !            vauxi_exm(14,ipoin,1:3) = 0.999548607287668_rp      !iFp=23);
        !            vauxi_exm(15,ipoin,1:3) = 0.999488774162635_rp      !iSp=24);
        !            vauxi_exm(16,ipoin,1:3) = 2.38076098345898e-09_rp      !d=25);
        !            vauxi_exm(17,ipoin,1:3) = 0.999999990696210_rp      !ff=26);
        !            vauxi_exm(18,ipoin,1:3) = 0.904906458666787_rp      !fs=27);
        !            vauxi_exm(19,ipoin,1:3) = 0.999999990696060_rp      !fcaf=28);
        !            vauxi_exm(20,ipoin,1:3) = 0.999581201974281_rp      !fcas=29);
        !            vauxi_exm(21,ipoin,1:3) = 0.999903346883777_rp      !jca=30);
        !            vauxi_exm(22,ipoin,1:3) = 0.00215555277945401_rp      !nca=31);
        !            vauxi_exm(23,ipoin,1:3) = 0.999999990680285_rp      !ffp=32);
        !            vauxi_exm(24,ipoin,1:3) = 0.999999990692529_rp      !fcafp=33);
        !            vauxi_exm(25,ipoin,1:3) = 8.64222375034682e-06_rp      !xrf=34);
        !            vauxi_exm(26,ipoin,1:3) = 0.487585264457487_rp      !xrs=35);
        !            vauxi_exm(27,ipoin,1:3) = 0.276203479404767_rp      !xs1=36);
        !            vauxi_exm(28,ipoin,1:3) = 0.000194412216700766_rp      !xs2=37);
        !            vauxi_exm(29,ipoin,1:3) = 0.996778581263402_rp      !xk1=38);
        !            vconc(9,ipoin,1:3) = 4.74946280300893e-07_rp       !Jrelnp=39);
        !            vconc(10,ipoin,1:3) = 5.93539009244893e-07_rp      !Jrelp=40);
        !            vconc(11,ipoin,1:3) = 0.0228529042639590_rp      !CaMKt=41);
     else if(ituss_exm == EXM_CELLTYPE_ENDO) then   !endocardial
        elmag(ipoin,1:3) = vminimate_exm(1,mat)     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        elmlo_exm(1:2) = vminimate_exm(1,mat)
        vicel_exm(1:26,ipoin,1) = 0.0_rp
        vconc(1:11,ipoin,1) = vcoin_exm(1:11,1,mat)      ! nai=2
        vconc(1:11,ipoin,2) = vcoin_exm(1:11,1,mat)       ! nai=2
        vconc(1:11,ipoin,3) = vcoin_exm(1:11,1,mat)      ! nai=2
        vauxi_exm(1:29,ipoin,1) = vauin_exm(1:29,1,mat) 
        vauxi_exm(1:29,ipoin,2) = vauin_exm(1:29,1,mat) 
        vauxi_exm(1:29,ipoin,3) = vauin_exm(1:29,1,mat) 
        !           elmag(1,ipoin,1:4) = vminimate_exm(1,mat) * 0.000003335640952_rp   
        !            vconc(1,ipoin,1:3) = 7.50131413867337_rp     ! nai=2
        !            vconc(2,ipoin,1:3) = 7.50141023194658_rp       !nass=3);
        !            vconc(3,ipoin,1:3) = 144.397336592521_rp       !ki=4);
        !            vconc(4,ipoin,1:3) = 144.397306654274_rp       !kss=5);
        !            vconc(5,ipoin,1:3) = 9.06297090590454e-05_rp       !cai=6);
        !            vconc(6,ipoin,1:3) = 8.99057612476716e-05_rp       !cass=7);
        !            vconc(7,ipoin,1:3) = 1.70096618787449_rp       !cansr=8);
        !            vconc(8,ipoin,1:3) = 1.61469442244241_rp       !cajsr=9);
        !            vauxi_exm(1,ipoin,1:3) = 0.00735179196722269_rp       !m=10);
        !            vauxi_exm(2,ipoin,1:3) = 0.697747459518797_rp       !hf=11);
        !            vauxi_exm(3,ipoin,1:3) = 0.697717242697596_rp       !hs=12);
        !            vauxi_exm(4,ipoin,1:3) = 0.697546140963364_rp       !j=13);
        !            vauxi_exm(5,ipoin,1:3) = 0.454478178913152_rp       !hsp=14);
        !            vauxi_exm(6,ipoin,1:3) = 0.697428694436028_rp       !jp=15);
        !            vauxi_exm(7,ipoin,1:3) = 0.000188633292218496_rp       !mL=16);
        !            vauxi_exm(8,ipoin,1:3) = 0.490629978954784_rp       !hL=17);
        !            vauxi_exm(9,ipoin,1:3) = 0.256189072459875_rp       !hLp=18);
        !            vauxi_exm(10,ipoin,1:3) = 0.00100180193761673_rp      !a=19);
        !            vauxi_exm(11,ipoin,1:3) = 0.999553321574707_rp      !iF=20);
        !            vauxi_exm(12,ipoin,1:3) = 0.519939585076356_rp      !iS=21);
        !            vauxi_exm(13,ipoin,1:3) = 0.000510445304511332_rp      !ap=22);
        !            vauxi_exm(14,ipoin,1:3) = 0.999553334882239_rp      !iFp=23);
        !            vauxi_exm(15,ipoin,1:3) = 0.519939585076356_rp      !iSp=24);
        !            vauxi_exm(16,ipoin,1:3) = 2.34656800943690e-09_rp      !d=25);
        !            vauxi_exm(17,ipoin,1:3) = 0.999999990847537_rp      !ff=26);
        !            vauxi_exm(18,ipoin,1:3) = 0.890968078383048_rp      !fs=27);
        !            vauxi_exm(19,ipoin,1:3) = 0.999999990847864_rp      !fcaf=28);
        !            vauxi_exm(20,ipoin,1:3) = 0.999331407376516_rp      !fcas=29);
        !            vauxi_exm(21,ipoin,1:3) = 0.999861393178778_rp      !jca=30);
        !            vauxi_exm(22,ipoin,1:3) = 0.00342031822852804_rp      !nca=31);
        !            vauxi_exm(23,ipoin,1:3) = 0.999999990811414_rp      !ffp=32);
        !            vauxi_exm(24,ipoin,1:3) = 0.999999990842929_rp      !fcafp=33);
        !            vauxi_exm(25,ipoin,1:3) = 9.06339108941520e-06_rp      !xrf=34);
        !            vauxi_exm(26,ipoin,1:3) = 0.519363745312448_rp      !xrs=35);
        !            vauxi_exm(27,ipoin,1:3) = 0.305861220188224_rp      !xs1=36);
        !            vauxi_exm(28,ipoin,1:3) = 0.000193120615511300_rp      !xs2=37);
        !            vauxi_exm(29,ipoin,1:3) = 0.996762767851278_rp      !xk1=38);
        !            vconc(9,ipoin,1:3) =  2.70257677876248e-07_rp      !Jrelnp=39);
        !            vconc(10,ipoin,1:3) = 3.37551159684025e-07_rp      !Jrelp=40);
        !            vconc(11,ipoin,1:3) = 0.0172511826393971_rp      !CaMKt=41);
     else if(ituss_exm == EXM_CELLTYPE_MID) then   !MIDmyocardial
        elmag(ipoin,1:3) = vminimate_exm(2,mat)     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        elmlo_exm(1:2) = vminimate_exm(2,mat)
        vicel_exm(1:26,ipoin,1) = 0.0_rp
        vconc(1:11,ipoin,1) = vcoin_exm(1:11,2,mat)       ! nai=2
        vconc(1:11,ipoin,2) = vcoin_exm(1:11,2,mat)       ! nai=2
        vconc(1:11,ipoin,3) = vcoin_exm(1:11,2,mat)       ! nai=2
        vauxi_exm(1:29,ipoin,1) = vauin_exm(1:29,2,mat) 
        vauxi_exm(1:29,ipoin,2) = vauin_exm(1:29,2,mat) 
        vauxi_exm(1:29,ipoin,3) = vauin_exm(1:29,2,mat) 
        !          elmag(1,ipoin,1:4) = vminimate_exm(2,mat) * 0.000003335640952_rp   
        !            vconc(1,ipoin,1:3) =  9.31318571775255_rp     ! nai=2
        !            vconc(2,ipoin,1:3) =  9.31335418275460_rp      !nass=3);
        !            vconc(3,ipoin,1:3) =  142.448454117416_rp      !ki=4);
        !            vconc(4,ipoin,1:3) =  142.448413204153_rp      !kss=5);
        !            vconc(5,ipoin,1:3) = 0.000109610482790345_rp       !cai=6);
        !            vconc(6,ipoin,1:3) = 0.000106921402857184_rp       !cass=7);
        !            vconc(7,ipoin,1:3) = 2.55909692264644_rp       !cansr=8);
        !            vconc(8,ipoin,1:3) = 2.46618898848524_rp       !cajsr=9);          
        !            vauxi_exm(1,ipoin,1:3) = 0.00762063517901599_rp       !m=10);
        !            vauxi_exm(2,ipoin,1:3) = 0.685227812017847_rp       !hf=11);
        !            vauxi_exm(3,ipoin,1:3) = 0.685186564524828_rp       !hs=12);
        !            vauxi_exm(4,ipoin,1:3) = 0.684946001497933_rp       !j=13);
        !            vauxi_exm(5,ipoin,1:3) = 0.439930036176900_rp       !hsp=14);
        !            vauxi_exm(6,ipoin,1:3) = 0.684770613692175_rp       !jp=15);
        !            vauxi_exm(7,ipoin,1:3) = 0.000201874904059089_rp       !mL=16);
        !            vauxi_exm(8,ipoin,1:3) = 0.469593913079625_rp       !hL=17);
        !            vauxi_exm(9,ipoin,1:3) = 0.232197586816947_rp       !hLp=18);
        !            vauxi_exm(10,ipoin,1:3) = 0.00102621879122940_rp      !a=19);
        !            vauxi_exm(11,ipoin,1:3) = 0.999524469433189_rp      !iF=20);
        !            vauxi_exm(12,ipoin,1:3) = 0.484655616738848_rp      !iS=21);
        !            vauxi_exm(13,ipoin,1:3) = 0.000522892623132084_rp      !ap=22);
        !            vauxi_exm(14,ipoin,1:3) = 0.999524488104748_rp      !iFp=23);
        !            vauxi_exm(15,ipoin,1:3) = 0.484655616738848_rp      !iSp=24);
        !            vauxi_exm(16,ipoin,1:3) = 2.55334758464344e-09_rp      !d=25);
        !            vauxi_exm(17,ipoin,1:3) = 0.999999989917362_rp      !ff=26);
        !            vauxi_exm(18,ipoin,1:3) = 0.850786415635368_rp      !fs=27);
        !            vauxi_exm(19,ipoin,1:3) = 0.999999989918100_rp      !fcaf=28);
        !            vauxi_exm(20,ipoin,1:3) = 0.998660902959982_rp      !fcas=29);
        !            vauxi_exm(21,ipoin,1:3) = 0.999737086062351_rp      !jca=30);
        !            vauxi_exm(22,ipoin,1:3) = 0.00660100554600241_rp      !nca=31);
        !            vauxi_exm(23,ipoin,1:3) = 0.999999989887965_rp      !ffp=32);
        !            vauxi_exm(24,ipoin,1:3) = 0.999999989911158_rp      !fcafp=33);
        !            vauxi_exm(25,ipoin,1:3) = 1.29078830800004e-05_rp      !xrf=34);
        !            vauxi_exm(26,ipoin,1:3) = 0.556418397118542_rp      !xrs=35);
        !            vauxi_exm(27,ipoin,1:3) = 0.365366368645388_rp      !xs1=36);
        !            vauxi_exm(28,ipoin,1:3) = 0.000201042604770464_rp      !xs2=37);
        !            vauxi_exm(29,ipoin,1:3) = 0.996855464208211_rp      !xk1=38);
        !            vconc(9,ipoin,1:3) =  1.89377036465905e-06_rp      !Jrelnp=39);
        !            vconc(10,ipoin,1:3) = 2.36711601348107e-06_rp      !Jrelp=40);
        !            vconc(11,ipoin,1:3) = 0.0406722957165916_rp       !CaMKt=41);
     end if
  end if
  !write(996,*) vauxi_exm(:,1,:)
end subroutine exm_iniohr
