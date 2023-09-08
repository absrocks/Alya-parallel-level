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
subroutine exm_iniohr(kmodel_ipoin,ipoin,mat)

  use      def_master
  use      def_domain
  use      def_elmtyp
  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: kmodel_ipoin, ipoin, mat
  !integer(ip)   :: mat
  !ituss = 0_ip  !!%endo = 1, epi = 0, M = 2  


  if(kfl_timei_exm==1_ip) then 
     !mat=kmodel_ipoin-3    
     ituss_exm = int(celty_exm(1,ipoin))

     if(ituss_exm == 3_ip) then   !epicardial
        elmag(ipoin,1:3) = vminimate_exm(3,mat)     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        elmlo_exm(1:2) = vminimate_exm(3,mat)
        vicel_exm(1:26,ipoin,1:3) = 0.0_rp
        vconc(1:11,ipoin,1) = vcoin_exm(1:11,3,mat) !8.17202753912090      ! nai=2
        vconc(1:11,ipoin,2) = vcoin_exm(1:11,3,mat) !8.17202753912090      ! nai=2
        vconc(1:11,ipoin,3) = vcoin_exm(1:11,3,mat) !8.17202753912090      ! nai=2
        vauxi_exm(1:29,ipoin,1) = vauin_exm(1:29,3,mat) 
        vauxi_exm(1:29,ipoin,2) = vauin_exm(1:29,3,mat) 
        vauxi_exm(1:29,ipoin,3) = vauin_exm(1:29,3,mat) 
        !           elmag(1,ipoin,1:4) = vminimate_exm(3,mat) * 0.000003335640952_rp    ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        !            vconc(2,ipoin,1:3) = 8.17212071948942       !nass=3);
        !            vconc(3,ipoin,1:3) = 143.675184333376       !ki=4);
        !            vconc(4,ipoin,1:3) = 143.675147146842       !kss=5);
        !            vconc(5,ipoin,1:3) = 8.03149767106260e-05       !cai=6);
        !            vconc(6,ipoin,1:3) = 7.96912949697674e-05       !cass=7);
        !            vconc(7,ipoin,1:3) = 2.14811455007091       !cansr=8);
        !            vconc(8,ipoin,1:3) = 2.03335899236798      !cajsr=9);
        !            vauxi_exm(1,ipoin,1:3) = 0.00739719746920272       !m=10);
        !            vauxi_exm(2,ipoin,1:3) = 0.695621622011335       !hf=11);
        !            vauxi_exm(3,ipoin,1:3) = 0.695601842634086       !hs=12);
        !            vauxi_exm(4,ipoin,1:3) = 0.695486248719023       !j=13);
        !            vauxi_exm(5,ipoin,1:3) = 0.452023628358454       !hsp=14);
        !            vauxi_exm(6,ipoin,1:3) = 0.695403157533235       !jp=15);
        !            vauxi_exm(7,ipoin,1:3) = 0.000190839777466418       !mL=16);
        !            vauxi_exm(8,ipoin,1:3) = 0.493606704642336       !hL=17);
        !            vauxi_exm(9,ipoin,1:3) = 0.264304293390731       !hLp=18);
        !            vauxi_exm(10,ipoin,1:3) = 0.00100594231451985      !a=19);
        !            vauxi_exm(11,ipoin,1:3) = 0.999548606668578      !iF=20);
        !            vauxi_exm(12,ipoin,1:3) = 0.999488774162635      !iS=21);
        !            vauxi_exm(13,ipoin,1:3) = 0.000512555980943569      !ap=22);
        !            vauxi_exm(14,ipoin,1:3) = 0.999548607287668      !iFp=23);
        !            vauxi_exm(15,ipoin,1:3) = 0.999488774162635      !iSp=24);
        !            vauxi_exm(16,ipoin,1:3) = 2.38076098345898e-09      !d=25);
        !            vauxi_exm(17,ipoin,1:3) = 0.999999990696210      !ff=26);
        !            vauxi_exm(18,ipoin,1:3) = 0.904906458666787      !fs=27);
        !            vauxi_exm(19,ipoin,1:3) = 0.999999990696060      !fcaf=28);
        !            vauxi_exm(20,ipoin,1:3) = 0.999581201974281      !fcas=29);
        !            vauxi_exm(21,ipoin,1:3) = 0.999903346883777      !jca=30);
        !            vauxi_exm(22,ipoin,1:3) = 0.00215555277945401      !nca=31);
        !            vauxi_exm(23,ipoin,1:3) = 0.999999990680285      !ffp=32);
        !            vauxi_exm(24,ipoin,1:3) = 0.999999990692529      !fcafp=33);
        !            vauxi_exm(25,ipoin,1:3) = 8.64222375034682e-06      !xrf=34);
        !            vauxi_exm(26,ipoin,1:3) = 0.487585264457487      !xrs=35);
        !            vauxi_exm(27,ipoin,1:3) = 0.276203479404767      !xs1=36);
        !            vauxi_exm(28,ipoin,1:3) = 0.000194412216700766      !xs2=37);
        !            vauxi_exm(29,ipoin,1:3) = 0.996778581263402      !xk1=38);
        !            vconc(9,ipoin,1:3) = 4.74946280300893e-07       !Jrelnp=39);
        !            vconc(10,ipoin,1:3) = 5.93539009244893e-07      !Jrelp=40);
        !            vconc(11,ipoin,1:3) = 0.0228529042639590      !CaMKt=41);
     else if(ituss_exm == 1_ip) then   !endocardial
        elmag(ipoin,1:3) = vminimate_exm(1,mat)     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        elmlo_exm(1:2) = vminimate_exm(1,mat)
        vicel_exm(1:26,ipoin,1:3) = 0.0_rp
        vconc(1:11,ipoin,1) = vcoin_exm(1:11,1,mat)      ! nai=2
        vconc(1:11,ipoin,2) = vcoin_exm(1:11,1,mat)       ! nai=2
        vconc(1:11,ipoin,3) = vcoin_exm(1:11,1,mat)      ! nai=2
        vauxi_exm(1:29,ipoin,1) = vauin_exm(1:29,1,mat) 
        vauxi_exm(1:29,ipoin,2) = vauin_exm(1:29,1,mat) 
        vauxi_exm(1:29,ipoin,3) = vauin_exm(1:29,1,mat) 
        !           elmag(1,ipoin,1:4) = vminimate_exm(1,mat) * 0.000003335640952_rp   
        !            vconc(1,ipoin,1:3) = 7.50131413867337     ! nai=2
        !            vconc(2,ipoin,1:3) = 7.50141023194658       !nass=3);
        !            vconc(3,ipoin,1:3) = 144.397336592521       !ki=4);
        !            vconc(4,ipoin,1:3) = 144.397306654274       !kss=5);
        !            vconc(5,ipoin,1:3) = 9.06297090590454e-05       !cai=6);
        !            vconc(6,ipoin,1:3) = 8.99057612476716e-05       !cass=7);
        !            vconc(7,ipoin,1:3) = 1.70096618787449       !cansr=8);
        !            vconc(8,ipoin,1:3) = 1.61469442244241       !cajsr=9);
        !            vauxi_exm(1,ipoin,1:3) = 0.00735179196722269       !m=10);
        !            vauxi_exm(2,ipoin,1:3) = 0.697747459518797       !hf=11);
        !            vauxi_exm(3,ipoin,1:3) = 0.697717242697596       !hs=12);
        !            vauxi_exm(4,ipoin,1:3) = 0.697546140963364       !j=13);
        !            vauxi_exm(5,ipoin,1:3) = 0.454478178913152       !hsp=14);
        !            vauxi_exm(6,ipoin,1:3) = 0.697428694436028       !jp=15);
        !            vauxi_exm(7,ipoin,1:3) = 0.000188633292218496       !mL=16);
        !            vauxi_exm(8,ipoin,1:3) = 0.490629978954784       !hL=17);
        !            vauxi_exm(9,ipoin,1:3) = 0.256189072459875       !hLp=18);
        !            vauxi_exm(10,ipoin,1:3) = 0.00100180193761673      !a=19);
        !            vauxi_exm(11,ipoin,1:3) = 0.999553321574707      !iF=20);
        !            vauxi_exm(12,ipoin,1:3) = 0.519939585076356      !iS=21);
        !            vauxi_exm(13,ipoin,1:3) = 0.000510445304511332      !ap=22);
        !            vauxi_exm(14,ipoin,1:3) = 0.999553334882239      !iFp=23);
        !            vauxi_exm(15,ipoin,1:3) = 0.519939585076356      !iSp=24);
        !            vauxi_exm(16,ipoin,1:3) = 2.34656800943690e-09      !d=25);
        !            vauxi_exm(17,ipoin,1:3) = 0.999999990847537      !ff=26);
        !            vauxi_exm(18,ipoin,1:3) = 0.890968078383048      !fs=27);
        !            vauxi_exm(19,ipoin,1:3) = 0.999999990847864      !fcaf=28);
        !            vauxi_exm(20,ipoin,1:3) = 0.999331407376516      !fcas=29);
        !            vauxi_exm(21,ipoin,1:3) = 0.999861393178778      !jca=30);
        !            vauxi_exm(22,ipoin,1:3) = 0.00342031822852804      !nca=31);
        !            vauxi_exm(23,ipoin,1:3) = 0.999999990811414      !ffp=32);
        !            vauxi_exm(24,ipoin,1:3) = 0.999999990842929      !fcafp=33);
        !            vauxi_exm(25,ipoin,1:3) = 9.06339108941520e-06      !xrf=34);
        !            vauxi_exm(26,ipoin,1:3) = 0.519363745312448      !xrs=35);
        !            vauxi_exm(27,ipoin,1:3) = 0.305861220188224      !xs1=36);
        !            vauxi_exm(28,ipoin,1:3) = 0.000193120615511300      !xs2=37);
        !            vauxi_exm(29,ipoin,1:3) = 0.996762767851278      !xk1=38);
        !            vconc(9,ipoin,1:3) =  2.70257677876248e-07      !Jrelnp=39);
        !            vconc(10,ipoin,1:3) = 3.37551159684025e-07      !Jrelp=40);
        !            vconc(11,ipoin,1:3) = 0.0172511826393971      !CaMKt=41);
     else if(ituss_exm == 2_ip) then   !MIDmyocardial
        elmag(ipoin,1:3) = vminimate_exm(2,mat)     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
        elmlo_exm(1:2) = vminimate_exm(2,mat)
        vicel_exm(1:26,ipoin,1:3) = 0.0_rp
        vconc(1:11,ipoin,1) = vcoin_exm(1:11,2,mat)       ! nai=2
        vconc(1:11,ipoin,2) = vcoin_exm(1:11,2,mat)       ! nai=2
        vconc(1:11,ipoin,3) = vcoin_exm(1:11,2,mat)       ! nai=2
        vauxi_exm(1:29,ipoin,1) = vauin_exm(1:29,2,mat) 
        vauxi_exm(1:29,ipoin,2) = vauin_exm(1:29,2,mat) 
        vauxi_exm(1:29,ipoin,3) = vauin_exm(1:29,2,mat) 
        !          elmag(1,ipoin,1:4) = vminimate_exm(2,mat) * 0.000003335640952_rp   
        !            vconc(1,ipoin,1:3) =  9.31318571775255     ! nai=2
        !            vconc(2,ipoin,1:3) =  9.31335418275460      !nass=3);
        !            vconc(3,ipoin,1:3) =  142.448454117416      !ki=4);
        !            vconc(4,ipoin,1:3) =  142.448413204153      !kss=5);
        !            vconc(5,ipoin,1:3) = 0.000109610482790345       !cai=6);
        !            vconc(6,ipoin,1:3) = 0.000106921402857184       !cass=7);
        !            vconc(7,ipoin,1:3) = 2.55909692264644       !cansr=8);
        !            vconc(8,ipoin,1:3) = 2.46618898848524       !cajsr=9);          
        !            vauxi_exm(1,ipoin,1:3) = 0.00762063517901599       !m=10);
        !            vauxi_exm(2,ipoin,1:3) = 0.685227812017847       !hf=11);
        !            vauxi_exm(3,ipoin,1:3) = 0.685186564524828       !hs=12);
        !            vauxi_exm(4,ipoin,1:3) = 0.684946001497933       !j=13);
        !            vauxi_exm(5,ipoin,1:3) = 0.439930036176900       !hsp=14);
        !            vauxi_exm(6,ipoin,1:3) = 0.684770613692175       !jp=15);
        !            vauxi_exm(7,ipoin,1:3) = 0.000201874904059089       !mL=16);
        !            vauxi_exm(8,ipoin,1:3) = 0.469593913079625       !hL=17);
        !            vauxi_exm(9,ipoin,1:3) = 0.232197586816947       !hLp=18);
        !            vauxi_exm(10,ipoin,1:3) = 0.00102621879122940      !a=19);
        !            vauxi_exm(11,ipoin,1:3) = 0.999524469433189      !iF=20);
        !            vauxi_exm(12,ipoin,1:3) = 0.484655616738848      !iS=21);
        !            vauxi_exm(13,ipoin,1:3) = 0.000522892623132084      !ap=22);
        !            vauxi_exm(14,ipoin,1:3) = 0.999524488104748      !iFp=23);
        !            vauxi_exm(15,ipoin,1:3) = 0.484655616738848      !iSp=24);
        !            vauxi_exm(16,ipoin,1:3) = 2.55334758464344e-09      !d=25);
        !            vauxi_exm(17,ipoin,1:3) = 0.999999989917362      !ff=26);
        !            vauxi_exm(18,ipoin,1:3) = 0.850786415635368      !fs=27);
        !            vauxi_exm(19,ipoin,1:3) = 0.999999989918100      !fcaf=28);
        !            vauxi_exm(20,ipoin,1:3) = 0.998660902959982      !fcas=29);
        !            vauxi_exm(21,ipoin,1:3) = 0.999737086062351      !jca=30);
        !            vauxi_exm(22,ipoin,1:3) = 0.00660100554600241      !nca=31);
        !            vauxi_exm(23,ipoin,1:3) = 0.999999989887965      !ffp=32);
        !            vauxi_exm(24,ipoin,1:3) = 0.999999989911158      !fcafp=33);
        !            vauxi_exm(25,ipoin,1:3) = 1.29078830800004e-05      !xrf=34);
        !            vauxi_exm(26,ipoin,1:3) = 0.556418397118542      !xrs=35);
        !            vauxi_exm(27,ipoin,1:3) = 0.365366368645388      !xs1=36);
        !            vauxi_exm(28,ipoin,1:3) = 0.000201042604770464      !xs2=37);
        !            vauxi_exm(29,ipoin,1:3) = 0.996855464208211      !xk1=38);
        !            vconc(9,ipoin,1:3) =  1.89377036465905e-06      !Jrelnp=39);
        !            vconc(10,ipoin,1:3) = 2.36711601348107e-06      !Jrelp=40);
        !            vconc(11,ipoin,1:3) = 0.0406722957165916       !CaMKt=41);
     end if
  end if
  !write(996,*) vauxi_exm(:,1,:)
end subroutine exm_iniohr
