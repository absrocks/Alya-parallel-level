!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_oneohr.f90
!> @date    12/04/2013
!> @author  Jazmin Aguado-Sierra
!> @brief   Sets initial conditions for Ohara-Rudy 2011 model
!> @details Runs for changes of cycle length and drug administration \n
!> @}
!------------------------------------------------------------------------
logical function exm_oneohr(mat) 
  
  use def_master
  use def_exmedi

  implicit none
  integer(ip), intent(in) :: mat
  logical     :: exm_oceohr, exm_oceola, flag_land
         
  exm_oneohr = .TRUE.  

  flag_land = .FALSE.
   
     if (kfl_eccty(mat) == 3_ip .or. kfl_eccty(mat)==4_ip ) flag_land = .TRUE.    
     
     
  if (((coupling('SOLIDZ','EXMEDI') >= 1_ip) .or. (coupling('EXMEDI','SOLIDZ') >= 1_ip)) .and. flag_land) then

     elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
     !elmlo_exm(2) = -87.99_rp
     ituss_exm = 3    !epicardial
     viclo_exm(1:26,1) = 0.0_rp      
     vcolo_exm(1,1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
     vcolo_exm(2,1:3) = 7.0_rp !8.17212071948942       !nass=3);
     vcolo_exm(3,1:3) = 145.0_rp !143.675184333376       !ki=4);
     vcolo_exm(4,1:3) = 145.0_rp  !143.675147146842       !kss=5);
     vcolo_exm(5,1:3) = 7.0_rp  !8.17202753912090      ! nai=2
     vcolo_exm(6,1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
     vcolo_exm(7,1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
     vcolo_exm(8,1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
     vaulo_exm(1,1:2) = 0.0_rp !0.00739719746920272       !m=10);
     vaulo_exm(2,1:2) = 1.0_rp  !0.695621622011335       !hf=11);
     vaulo_exm(3,1:2) = 1.0_rp  !0.695601842634086       !hs=12);
     vaulo_exm(4,1:2) = 1.0_rp  !0.695486248719023       !j=13);
     vaulo_exm(5,1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
     vaulo_exm(6,1:2) = 1.0_rp  !0.695403157533235       !jp=15);
     vaulo_exm(7,1:2) = 0.0_rp !0.000190839777466418       !mL=16);
     vaulo_exm(8,1:2) = 1.0_rp  !0.493606704642336       !hL=17);
     vaulo_exm(9,1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
     vaulo_exm(10,1:2) = 0.0_rp !0.00100594231451985      !a=19);
     vaulo_exm(11,1:2) = 1.0_rp  !0.999548606668578      !iF=20);
     vaulo_exm(12,1:2) = 1.0_rp  !0.999488774162635      !iS=21);
     vaulo_exm(13,1:2) = 0.0_rp !0.000512555980943569      !ap=22);
     vaulo_exm(14,1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
     vaulo_exm(15,1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
     vaulo_exm(16,1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
     vaulo_exm(17,1:2) = 1.0_rp  !0.999999990696210      !ff=26);
     vaulo_exm(18,1:2) = 1.0_rp  !0.904906458666787      !fs=27);
     vaulo_exm(19,1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
     vaulo_exm(20,1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
     vaulo_exm(21,1:2) = 1.0_rp  !0.999903346883777      !jca=30);
     vaulo_exm(22,1:2) = 0.0_rp !0.00215555277945401      !nca=31);
     vaulo_exm(23,1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
     vaulo_exm(24,1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
     vaulo_exm(25,1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
     vaulo_exm(26,1:2) = 0.0_rp !0.487585264457487      !xrs=35);
     vaulo_exm(27,1:2) = 0.0_rp !0.276203479404767      !xs1=36);
     vaulo_exm(28,1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
     vaulo_exm(29,1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
     vcolo_exm(9,1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
     vcolo_exm(10,1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
     vcolo_exm(11,1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);  

     !endocardium
     ituss_exm = 1  !endocardial
     elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
     !elmlo_exm(2) = -87.99_rp
     viclo_exm(1:26,1) = 0.0_rp      
     vcolo_exm(1,1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
     vcolo_exm(2,1:3) = 7.0_rp !8.17212071948942       !nass=3);
     vcolo_exm(3,1:3) = 145.0_rp !143.675184333376       !ki=4);
     vcolo_exm(4,1:3) = 145.0_rp  !143.675147146842       !kss=5);
     vcolo_exm(5,1:3) = 7.0_rp  !8.17202753912090      ! nai=2
     vcolo_exm(6,1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
     vcolo_exm(7,1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
     vcolo_exm(8,1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
     vaulo_exm(1,1:2) = 0.0_rp !0.00739719746920272       !m=10);
     vaulo_exm(2,1:2) = 1.0_rp  !0.695621622011335       !hf=11);
     vaulo_exm(3,1:2) = 1.0_rp  !0.695601842634086       !hs=12);
     vaulo_exm(4,1:2) = 1.0_rp  !0.695486248719023       !j=13);
     vaulo_exm(5,1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
     vaulo_exm(6,1:2) = 1.0_rp  !0.695403157533235       !jp=15);
     vaulo_exm(7,1:2) = 0.0_rp !0.000190839777466418       !mL=16);
     vaulo_exm(8,1:2) = 1.0_rp  !0.493606704642336       !hL=17);
     vaulo_exm(9,1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
     vaulo_exm(10,1:2) = 0.0_rp !0.00100594231451985      !a=19);
     vaulo_exm(11,1:2) = 1.0_rp  !0.999548606668578      !iF=20);
     vaulo_exm(12,1:2) = 1.0_rp  !0.999488774162635      !iS=21);
     vaulo_exm(13,1:2) = 0.0_rp !0.000512555980943569      !ap=22);
     vaulo_exm(14,1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
     vaulo_exm(15,1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
     vaulo_exm(16,1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
     vaulo_exm(17,1:2) = 1.0_rp  !0.999999990696210      !ff=26);
     vaulo_exm(18,1:2) = 1.0_rp  !0.904906458666787      !fs=27);
     vaulo_exm(19,1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
     vaulo_exm(20,1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
     vaulo_exm(21,1:2) = 1.0_rp  !0.999903346883777      !jca=30);
     vaulo_exm(22,1:2) = 0.0_rp !0.00215555277945401      !nca=31);
     vaulo_exm(23,1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
     vaulo_exm(24,1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
     vaulo_exm(25,1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
     vaulo_exm(26,1:2) = 0.0_rp !0.487585264457487      !xrs=35);
     vaulo_exm(27,1:2) = 0.0_rp !0.276203479404767      !xs1=36);
     vaulo_exm(28,1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
     vaulo_exm(29,1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
     vcolo_exm(9,1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
     vcolo_exm(10,1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
     vcolo_exm(11,1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);    
    
     ituss_exm = 2  !MIDmyocardial
     elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
     !elmlo_exm(2) = -87.99_rp
     viclo_exm(1:26,1) = 0.0_rp      
     vcolo_exm(1,1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
     vcolo_exm(2,1:3) = 7.0_rp !8.17212071948942       !nass=3);
     vcolo_exm(3,1:3) = 145.0_rp !143.675184333376       !ki=4);
     vcolo_exm(4,1:3) = 145.0_rp  !143.675147146842       !kss=5);
     vcolo_exm(5,1:3) = 7.0_rp  !8.17202753912090      ! nai=2
     vcolo_exm(6,1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
     vcolo_exm(7,1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
     vcolo_exm(8,1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
     vaulo_exm(1,1:2) = 0.0_rp !0.00739719746920272       !m=10);
     vaulo_exm(2,1:2) = 1.0_rp  !0.695621622011335       !hf=11);
     vaulo_exm(3,1:2) = 1.0_rp  !0.695601842634086       !hs=12);
     vaulo_exm(4,1:2) = 1.0_rp  !0.695486248719023       !j=13);
     vaulo_exm(5,1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
     vaulo_exm(6,1:2) = 1.0_rp  !0.695403157533235       !jp=15);
     vaulo_exm(7,1:2) = 0.0_rp !0.000190839777466418       !mL=16);
     vaulo_exm(8,1:2) = 1.0_rp  !0.493606704642336       !hL=17);
     vaulo_exm(9,1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
     vaulo_exm(10,1:2) = 0.0_rp !0.00100594231451985      !a=19);
     vaulo_exm(11,1:2) = 1.0_rp  !0.999548606668578      !iF=20);
     vaulo_exm(12,1:2) = 1.0_rp  !0.999488774162635      !iS=21);
     vaulo_exm(13,1:2) = 0.0_rp !0.000512555980943569      !ap=22);
     vaulo_exm(14,1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
     vaulo_exm(15,1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
     vaulo_exm(16,1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
     vaulo_exm(17,1:2) = 1.0_rp  !0.999999990696210      !ff=26);
     vaulo_exm(18,1:2) = 1.0_rp  !0.904906458666787      !fs=27);
     vaulo_exm(19,1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
     vaulo_exm(20,1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
     vaulo_exm(21,1:2) = 1.0_rp  !0.999903346883777      !jca=30);
     vaulo_exm(22,1:2) = 0.0_rp !0.00215555277945401      !nca=31);
     vaulo_exm(23,1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
     vaulo_exm(24,1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
     vaulo_exm(25,1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
     vaulo_exm(26,1:2) = 0.0_rp !0.487585264457487      !xrs=35);
     vaulo_exm(27,1:2) = 0.0_rp !0.276203479404767      !xs1=36);
     vaulo_exm(28,1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
     vaulo_exm(29,1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
     vcolo_exm(9,1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
     vcolo_exm(10,1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
     vcolo_exm(11,1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);  

     if (exm_oneohr) then
        exm_oneohr = exm_oceola(mat)
     end if
  end if

  if (kfl_hfmodmate_exm(mat) == 0_ip) then
    
    if ((moneclmate_exm(2,mat) == 1000) .and. (kfl_drugsmate_exm(mat) == 0_ip)) then  !if it's Normal and at 1000 bpm
            !
            vminimate_exm(3,mat) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
    
            ituss_exm = 3_ip    !epicardial
            viclo_exm(1:26,1) = 0.0_rp      
            vcoin_exm(1,3,mat) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2,3,mat) = 7.0_rp !8.17212071948942       !nass=3);
            vcoin_exm(3,3,mat) = 145.0_rp !143.675184333376       !ki=4);
            vcoin_exm(4,3,mat) = 145.0_rp  !143.675147146842       !kss=5);
            vcoin_exm(5,3,mat) = 7.0_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6,3,mat) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7,3,mat) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8,3,mat) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1,3,mat) = 0.0_rp !0.00739719746920272       !m=10);
            vauin_exm(2,3,mat) = 1.0_rp  !0.695621622011335       !hf=11);
            vauin_exm(3,3,mat) = 1.0_rp  !0.695601842634086       !hs=12);
            vauin_exm(4,3,mat) = 1.0_rp  !0.695486248719023       !j=13);
            vauin_exm(5,3,mat) = 1.0_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6,3,mat) = 1.0_rp  !0.695403157533235       !jp=15);
            vauin_exm(7,3,mat) = 0.0_rp !0.000190839777466418       !mL=16);
            vauin_exm(8,3,mat) = 1.0_rp  !0.493606704642336       !hL=17);
            vauin_exm(9,3,mat) = 1.0_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10,3,mat) = 0.0_rp !0.00100594231451985      !a=19);
            vauin_exm(11,3,mat) = 1.0_rp  !0.999548606668578      !iF=20);
            vauin_exm(12,3,mat) = 1.0_rp  !0.999488774162635      !iS=21);
            vauin_exm(13,3,mat) = 0.0_rp !0.000512555980943569      !ap=22);
            vauin_exm(14,3,mat) = 1.0_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15,3,mat) = 1.0_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16,3,mat) = 0.0_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17,3,mat) = 1.0_rp  !0.999999990696210      !ff=26);
            vauin_exm(18,3,mat) = 1.0_rp  !0.904906458666787      !fs=27);
            vauin_exm(19,3,mat) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20,3,mat) = 1.0_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21,3,mat) = 1.0_rp  !0.999903346883777      !jca=30);
            vauin_exm(22,3,mat) = 0.0_rp !0.00215555277945401      !nca=31);
            vauin_exm(23,3,mat) = 1.0_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24,3,mat) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25,3,mat) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26,3,mat) = 0.0_rp !0.487585264457487      !xrs=35);
            vauin_exm(27,3,mat) = 0.0_rp !0.276203479404767      !xs1=36);
            vauin_exm(28,3,mat) = 0.0_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29,3,mat) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9,3,mat) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10,3,mat) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11,3,mat) = 0.0_rp !0.0228529042639590      !CaMKt=41);     
                  
            !endocardium
            ituss_exm = 1_ip  !endocardial
            vminimate_exm(3,mat) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26,1) = 0.0_rp      
            vcoin_exm(1,1,mat) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2,1,mat) = 7.0_rp !8.17212071948942       !nass=3);
            vcoin_exm(3,1,mat) = 145.0_rp !143.675184333376       !ki=4);
            vcoin_exm(4,1,mat) = 145.0_rp  !143.675147146842       !kss=5);
            vcoin_exm(5,1,mat) = 7.0_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6,1,mat) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7,1,mat) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8,1,mat) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1,1,mat) = 0.0_rp !0.00739719746920272       !m=10);
            vauin_exm(2,1,mat) = 1.0_rp  !0.695621622011335       !hf=11);
            vauin_exm(3,1,mat) = 1.0_rp  !0.695601842634086       !hs=12);
            vauin_exm(4,1,mat) = 1.0_rp  !0.695486248719023       !j=13);
            vauin_exm(5,1,mat) = 1.0_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6,1,mat) = 1.0_rp  !0.695403157533235       !jp=15);
            vauin_exm(7,1,mat) = 0.0_rp !0.000190839777466418       !mL=16);
            vauin_exm(8,1,mat) = 1.0_rp  !0.493606704642336       !hL=17);
            vauin_exm(9,1,mat) = 1.0_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10,1,mat) = 0.0_rp !0.00100594231451985      !a=19);
            vauin_exm(11,1,mat) = 1.0_rp  !0.999548606668578      !iF=20);
            vauin_exm(12,1,mat) = 1.0_rp  !0.999488774162635      !iS=21);
            vauin_exm(13,1,mat) = 0.0_rp !0.000512555980943569      !ap=22);
            vauin_exm(14,1,mat) = 1.0_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15,1,mat) = 1.0_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16,1,mat) = 0.0_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17,1,mat) = 1.0_rp  !0.999999990696210      !ff=26);
            vauin_exm(18,1,mat) = 1.0_rp  !0.904906458666787      !fs=27);
            vauin_exm(19,1,mat) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20,1,mat) = 1.0_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21,1,mat) = 1.0_rp  !0.999903346883777      !jca=30);
            vauin_exm(22,1,mat) = 0.0_rp !0.00215555277945401      !nca=31);
            vauin_exm(23,1,mat) = 1.0_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24,1,mat) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25,1,mat) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26,1,mat) = 0.0_rp !0.487585264457487      !xrs=35);
            vauin_exm(27,1,mat) = 0.0_rp !0.276203479404767      !xs1=36);
            vauin_exm(28,1,mat) = 0.0_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29,1,mat) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9,1,mat) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10,1,mat) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11,1,mat) = 0.0_rp !0.0228529042639590      !CaMKt=41);                   
    
            ituss_exm = 2_ip  !MIDmyocardial
            vminimate_exm(2,mat) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26,1) = 0.0_rp      
            vcoin_exm(1,2,mat) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2,2,mat) = 7.0_rp !8.17212071948942       !nass=3);
            vcoin_exm(3,2,mat) = 145.0_rp !143.675184333376       !ki=4);
            vcoin_exm(4,2,mat) = 145.0_rp  !143.675147146842       !kss=5);
            vcoin_exm(5,2,mat) = 7.0_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6,2,mat) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7,2,mat) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8,2,mat) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1,2,mat) = 0.0_rp !0.00739719746920272       !m=10);
            vauin_exm(2,2,mat) = 1.0_rp  !0.695621622011335       !hf=11);
            vauin_exm(3,2,mat) = 1.0_rp  !0.695601842634086       !hs=12);
            vauin_exm(4,2,mat) = 1.0_rp  !0.695486248719023       !j=13);
            vauin_exm(5,2,mat) = 1.0_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6,2,mat) = 1.0_rp  !0.695403157533235       !jp=15);
            vauin_exm(7,2,mat) = 0.0_rp !0.000190839777466418       !mL=16);
            vauin_exm(8,2,mat) = 1.0_rp  !0.493606704642336       !hL=17);
            vauin_exm(9,2,mat) = 1.0_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10,2,mat) = 0.0_rp !0.00100594231451985      !a=19);
            vauin_exm(11,2,mat) = 1.0_rp  !0.999548606668578      !iF=20);
            vauin_exm(12,2,mat) = 1.0_rp  !0.999488774162635      !iS=21);
            vauin_exm(13,2,mat) = 0.0_rp !0.000512555980943569      !ap=22);
            vauin_exm(14,2,mat) = 1.0_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15,2,mat) = 1.0_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16,2,mat) = 0.0_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17,2,mat) = 1.0_rp  !0.999999990696210      !ff=26);
            vauin_exm(18,2,mat) = 1.0_rp  !0.904906458666787      !fs=27);
            vauin_exm(19,2,mat) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20,2,mat) = 1.0_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21,2,mat) = 1.0_rp  !0.999903346883777      !jca=30);
            vauin_exm(22,2,mat) = 0.0_rp !0.00215555277945401      !nca=31);
            vauin_exm(23,2,mat) = 1.0_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24,2,mat) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25,2,mat) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26,2,mat) = 0.0_rp !0.487585264457487      !xrs=35);
            vauin_exm(27,2,mat) = 0.0_rp !0.276203479404767      !xs1=36);
            vauin_exm(28,2,mat) = 0.0_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29,2,mat) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9,2,mat) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10,2,mat) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11,2,mat) = 0.0_rp !0.0228529042639590      !CaMKt=41);  
    
    
       else if((moneclmate_exm(2,mat) == 857) .and. (kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_inaga_exm(mat) == 1_ip)) then   !Normal at 70 bpm
    
            vminimate_exm(3,mat) = -87.99_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !ituss_exm = 3  !epicardial
             !epicardium   
            vcoin_exm(1,3,mat) = 0.0000803149767106260_rp       !cai=6);
            vcoin_exm(2,3,mat) = 8.17212071948942_rp       !nass=3);
            vcoin_exm(3,3,mat) = 143.675184333376_rp       !ki=4);
            vcoin_exm(4,3,mat) = 143.675147146842_rp       !kss=5);
            vcoin_exm(5,3,mat) = 8.17202753912090_rp      ! nai=2
            vcoin_exm(6,3,mat) = 0.0000796912949697674_rp       !cass=7);
            vcoin_exm(7,3,mat) = 2.14811455007091_rp       !cansr=8);
            vcoin_exm(8,3,mat) = 2.03335899236798_rp      !cajsr=9);
            vauin_exm(1,3,mat) = 0.00739719746920272_rp       !m=10);
            vauin_exm(2,3,mat) = 0.695621622011335_rp       !hf=11);
            vauin_exm(3,3,mat) = 0.695601842634086_rp       !hs=12);
            vauin_exm(4,3,mat) = 0.695486248719023_rp       !j=13);
            vauin_exm(5,3,mat) = 0.452023628358454_rp       !hsp=14);
            vauin_exm(6,3,mat) = 0.695403157533235_rp       !jp=15);
            vauin_exm(7,3,mat) = 0.000190839777466418_rp       !mL=16);
            vauin_exm(8,3,mat) = 0.493606704642336_rp       !hL=17);
            vauin_exm(9,3,mat) = 0.264304293390731_rp       !hLp=18);
            vauin_exm(10,3,mat) = 0.00100594231451985_rp      !a=19);
            vauin_exm(11,3,mat) = 0.999548606668578_rp      !iF=20);
            vauin_exm(12,3,mat) = 0.999488774162635_rp      !iS=21);
            vauin_exm(13,3,mat) = 0.000512555980943569_rp      !ap=22);
            vauin_exm(14,3,mat) = 0.999548607287668_rp      !iFp=23);
            vauin_exm(15,3,mat) = 0.999488774162635_rp      !iSp=24);
            vauin_exm(16,3,mat) = 0.00000000238076098345898_rp      !d=25);
            vauin_exm(17,3,mat) = 0.999999990696210_rp      !ff=26);
            vauin_exm(18,3,mat) = 0.904906458666787_rp      !fs=27);
            vauin_exm(19,3,mat) = 0.999999990696060_rp      !fcaf=28);
            vauin_exm(20,3,mat) = 0.999581201974281_rp      !fcas=29);
            vauin_exm(21,3,mat) = 0.999903346883777_rp      !jca=30);
            vauin_exm(22,3,mat) = 0.00215555277945401_rp      !nca=31);
            vauin_exm(23,3,mat) = 0.999999990680285_rp      !ffp=32);
            vauin_exm(24,3,mat) = 0.999999990692529_rp      !fcafp=33);
            vauin_exm(25,3,mat) = 0.00000864222375034682_rp     !xrf=34);
            vauin_exm(26,3,mat) = 0.487585264457487_rp      !xrs=35);
            vauin_exm(27,3,mat) = 0.276203479404767_rp      !xs1=36);
            vauin_exm(28,3,mat) = 0.000194412216700766_rp     !xs2=37);
            vauin_exm(29,3,mat) = 0.996778581263402_rp      !xk1=38);
            vcoin_exm(9,3,mat) = 0.000000474946280300893_rp       !Jrelnp=39);
            vcoin_exm(10,3,mat) = 0.000000593539009244893_rp      !Jrelp=40);
            vcoin_exm(11,3,mat) = 0.0228529042639590_rp      !CaMKt=41);
    
             !endocardium
            vminimate_exm(1,mat) = -87.99_rp
            !ituss_exm = 1  !endocardial
            vcoin_exm(1,1,mat) = 0.0000906297090590454_rp       !cai=6); 
            vcoin_exm(2,1,mat) = 7.50141023194658_rp       !nass=3);
            vcoin_exm(3,1,mat) = 144.397336592521_rp       !ki=4);
            vcoin_exm(4,1,mat) = 144.397306654274_rp       !kss=5);
            vcoin_exm(5,1,mat) = 7.50131413867337_rp     ! nai=2
            vcoin_exm(6,1,mat) = 0.0000899057612476716_rp       !cass=7);
            vcoin_exm(7,1,mat) = 1.70096618787449_rp       !cansr=8);
            vcoin_exm(8,1,mat) = 1.61469442244241_rp      !cajsr=9);
            vauin_exm(1,1,mat) = 0.00735179196722269_rp       !m=10);
            vauin_exm(2,1,mat) = 0.697747459518797_rp       !hf=11);
            vauin_exm(3,1,mat) = 0.697717242697596_rp       !hs=12);
            vauin_exm(4,1,mat) = 0.697546140963364_rp       !j=13);
            vauin_exm(5,1,mat) = 0.454478178913152_rp       !hsp=14);
            vauin_exm(6,1,mat) = 0.697428694436028_rp       !jp=15);
            vauin_exm(7,1,mat) = 0.000188633292218496_rp       !mL=16);
            vauin_exm(8,1,mat) = 0.490629978954784_rp       !hL=17);
            vauin_exm(9,1,mat) = 0.256189072459875_rp       !hLp=18);
            vauin_exm(10,1,mat) = 0.00100180193761673_rp      !a=19);
            vauin_exm(11,1,mat) = 0.999553321574707_rp      !iF=20);
            vauin_exm(12,1,mat) = 0.519939585076356_rp      !iS=21);
            vauin_exm(13,1,mat) = 0.000510445304511332_rp      !ap=22);
            vauin_exm(14,1,mat) = 0.999553334882239_rp      !iFp=23);
            vauin_exm(15,1,mat) = 0.519939585076356_rp      !iSp=24);
            vauin_exm(16,1,mat) = 0.00000000234656800943690_rp      !d=25);
            vauin_exm(17,1,mat) = 0.999999990847537_rp      !ff=26);
            vauin_exm(18,1,mat) = 0.890968078383048_rp      !fs=27);
            vauin_exm(19,1,mat) = 0.999999990847864_rp      !fcaf=28);
            vauin_exm(20,1,mat) = 0.999331407376516_rp      !fcas=29);
            vauin_exm(21,1,mat) = 0.999861393178778_rp      !jca=30);
            vauin_exm(22,1,mat) = 0.00342031822852804_rp      !nca=31);
            vauin_exm(23,1,mat) = 0.999999990811414_rp      !ffp=32);
            vauin_exm(24,1,mat) = 0.999999990842929_rp      !fcafp=33);
            vauin_exm(25,1,mat) = 0.00000906339108941520_rp      !xrf=34);
            vauin_exm(26,1,mat) = 0.519363745312448_rp      !xrs=35);
            vauin_exm(27,1,mat) = 0.305861220188224_rp      !xs1=36);
            vauin_exm(28,1,mat) = 0.000193120615511300_rp      !xs2=37);
            vauin_exm(29,1,mat) = 0.996762767851278_rp      !xk1=38);
            vcoin_exm(9,1,mat) =  0.000000270257677876248_rp      !Jrelnp=39);
            vcoin_exm(10,1,mat) = 0.000000337551159684025_rp      !Jrelp=40);
            vcoin_exm(11,1,mat) = 0.0172511826393971_rp      !CaMKt=41);
             
             !midmyocardium
            vminimate_exm = -87.99_rp
            !ituss_exm = 2  !midmyocardial
            vcoin_exm(1,2,mat) =  0.000109610482790345_rp       !cai=6);
            vcoin_exm(2,2,mat) =  9.31335418275460_rp      !nass=3);
            vcoin_exm(3,2,mat) =  142.448454117416_rp      !ki=4);
            vcoin_exm(4,2,mat) =  142.448413204153_rp      !kss=5);
            vcoin_exm(5,2,mat) =  9.31318571775255_rp     ! nai=2
            vcoin_exm(6,2,mat) = 0.000106921402857184_rp       !cass=7);
            vcoin_exm(7,2,mat) = 2.55909692264644_rp       !cansr=8);
            vcoin_exm(8,2,mat) = 2.46618898848524_rp       !cajsr=9);          
            vauin_exm(1,2,mat) = 0.00762063517901599_rp       !m=10);
            vauin_exm(2,2,mat) = 0.685227812017847_rp       !hf=11);
            vauin_exm(3,2,mat) = 0.685186564524828_rp       !hs=12);
            vauin_exm(4,2,mat) = 0.684946001497933_rp       !j=13);
            vauin_exm(5,2,mat) = 0.439930036176900_rp       !hsp=14);
            vauin_exm(6,2,mat) = 0.684770613692175_rp       !jp=15);
            vauin_exm(7,2,mat) = 0.000201874904059089_rp       !mL=16);
            vauin_exm(8,2,mat) = 0.469593913079625_rp       !hL=17);
            vauin_exm(9,2,mat) = 0.232197586816947_rp       !hLp=18);
            vauin_exm(10,2,mat) = 0.00102621879122940_rp      !a=19);
            vauin_exm(11,2,mat) = 0.999524469433189_rp      !iF=20);
            vauin_exm(12,2,mat) = 0.484655616738848_rp      !iS=21);
            vauin_exm(13,2,mat) = 0.000522892623132084_rp      !ap=22);
            vauin_exm(14,2,mat) = 0.999524488104748_rp      !iFp=23);
            vauin_exm(15,2,mat) = 0.484655616738848_rp      !iSp=24);
            vauin_exm(16,2,mat) = 0.00000000255334758464344_rp      !d=25);
            vauin_exm(17,2,mat) = 0.999999989917362_rp      !ff=26);
            vauin_exm(18,2,mat) = 0.850786415635368_rp      !fs=27);
            vauin_exm(19,2,mat) = 0.999999989918100_rp      !fcaf=28);
            vauin_exm(20,2,mat) = 0.998660902959982_rp      !fcas=29);
            vauin_exm(21,2,mat) = 0.999737086062351_rp      !jca=30);
            vauin_exm(22,2,mat) = 0.00660100554600241_rp      !nca=31);
            vauin_exm(23,2,mat) = 0.999999989887965_rp      !ffp=32);
            vauin_exm(24,2,mat) = 0.999999989911158_rp      !fcafp=33);
            vauin_exm(25,2,mat) = 0.0000129078830800004_rp      !xrf=34);
            vauin_exm(26,2,mat) = 0.556418397118542_rp      !xrs=35);
            vauin_exm(27,2,mat) = 0.365366368645388_rp      !xs1=36);
            vauin_exm(28,2,mat) = 0.000201042604770464_rp      !xs2=37);
            vauin_exm(29,2,mat) = 0.996855464208211_rp      !xk1=38);
            vcoin_exm(9,2,mat) =  0.00000189377036465905_rp      !Jrelnp=39);
            vcoin_exm(10,2,mat) = 0.00000236711601348107_rp      !Jrelp=40);
            vcoin_exm(11,2,mat) = 0.0406722957165916_rp       !CaMKt=41);
             
               
       !else if (kfl_hfmod_exm(mat) == 1 .and. kfl_hfmodmate_exm(mat) == 0 ) then  !if it's a Normal PIG            
       else if((moneclmate_exm(2,mat) == 400) .and. (kfl_hfmod_exm(mat) == 1_ip) .and. (kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_inaga_exm(mat) == 1_ip)) then   !Normal PIG AT 400 BCL
        
            vminimate_exm(1,mat) = -87.8061985708231_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL

            ituss_exm = 3_ip    !epicardial
            viclo_exm(1:26,1) = 0.0_rp      
            
            vcoin_exm(1,3,mat) = 0.0001040361496276368_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2,3,mat) =  7.43617875723091_rp !8.17212071948942       !nass=3);
            vcoin_exm(3,3,mat) = 144.311817237457_rp !143.675184333376       !ki=4);
            vcoin_exm(4,3,mat) = 144.311793923281_rp  !143.675147146842       !kss=5);
            vcoin_exm(5,3,mat) = 7.43604263653445_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6,3,mat) = 0.0001102064090489676_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7,3,mat) = 2.63411603951980_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8,3,mat) = 1.76343231449702_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1,3,mat) = 0.007490096291770637_rp !0.00739719746920272       !m=10);
            vauin_exm(2,3,mat) = 0.816838169651186_rp  !0.695621622011335       !hf=11);
            vauin_exm(3,3,mat) = 0.816869319742196_rp  !0.695601842634086       !hs=12);
            vauin_exm(4,3,mat) = 0.779416631991491_rp  !0.695486248719023       !j=13);
            vauin_exm(5,3,mat) = 0.621614033150460_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6,3,mat) = 0.719030208978245_rp  !0.695403157533235       !jp=15);
            vauin_exm(7,3,mat) = 0.0001953921132645023_rp !0.000190839777466418       !mL=16);
            vauin_exm(8,3,mat) = 0.350601487895652_rp  !0.493606704642336       !hL=17);
            vauin_exm(9,3,mat) = 0.165483509753837_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10,3,mat) = 0.001014440468937102_rp !0.00100594231451985      !a=19);
            vauin_exm(11,3,mat) = 0.999538652553908_rp  !0.999548606668578      !iF=20);
            vauin_exm(12,3,mat) = 0.940646485484524_rp  !0.999488774162635      !iS=21);
            vauin_exm(13,3,mat) = 0.0005168881869164542_rp !0.000512555980943569      !ap=22);
            vauin_exm(14,3,mat) = 0.999538658574815_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15,3,mat) = 0.961396539930704_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16,3,mat) = 0.000000002451977924935544_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17,3,mat) = 0.999999990358952_rp  !0.999999990696210      !ff=26);
            vauin_exm(18,3,mat) = 0.721766724570291_rp  !0.904906458666787      !fs=27);
            vauin_exm(19,3,mat) = 0.999999990360122_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20,3,mat) = 0.952750811437167_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21,3,mat) = 0.962655802194759_rp  !0.999903346883777      !jca=30);
            vauin_exm(22,3,mat) = 0.007755537033307213_rp !0.00215555277945401      !nca=31);
            vauin_exm(23,3,mat) = 0.999989171225269_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24,3,mat) = 0.999998595329089_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25,3,mat) = 0.01177184211893142_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26,3,mat) = 0.803851729756987_rp !0.487585264457487      !xrs=35);
            vauin_exm(27,3,mat) = 0.525297568917089_rp !0.276203479404767      !xs1=36);
            vauin_exm(28,3,mat) = 0.0004266215773384047_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29,3,mat) = 0.996814746325558_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9,3,mat) = 0.0000003867443373607881_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10,3,mat) = 0.0000004808818699484394_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11,3,mat) = 0.08284541646773119_rp !0.0228529042639590      !CaMKt=41);     
                        
            vminimate_exm(1,mat) = -87.8077947685269_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            ituss_exm = 1    !endocardial
            
            viclo_exm(1:26,1) = 0.0_rp      
            vcoin_exm(1,1,mat) = 0.0001244664562667345_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2,1,mat) = 7.38274584190492_rp !8.17212071948942       !nass=3);
            vcoin_exm(3,1,mat) = 144.434867863316_rp !143.675184333376       !ki=4);
            vcoin_exm(4,1,mat) = 144.434843290996_rp  !143.675147146842       !kss=5);
            vcoin_exm(5,1,mat) = 7.38257910315418_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6,1,mat) = 0.0001303227834520188_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7,1,mat) = 2.22033410012224_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8,1,mat) = 1.57032538537506_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1,1,mat) = 0.007488900815703561_rp !0.00739719746920272       !m=10);
            vauin_exm(2,1,mat) = 0.816818226826766_rp  !0.695621622011335       !hf=11);
            vauin_exm(3,1,mat) = 0.816860080147452_rp  !0.695601842634086       !hs=12);
            vauin_exm(4,1,mat) = 0.760420946701455_rp  !0.695486248719023       !j=13);
            vauin_exm(5,1,mat) = 0.621427482142547_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6,1,mat) = 0.686823279932452_rp  !0.695403157533235       !jp=15);
            vauin_exm(7,1,mat) = 0.0001953332078784692_rp !0.000190839777466418       !mL=16);
            vauin_exm(8,1,mat) = 0.321919485372494_rp  !0.493606704642336       !hL=17);
            vauin_exm(9,1,mat) = 0.146893625330722_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10,1,mat) = 0.001014350908419775_rp !0.00100594231451985      !a=19);
            vauin_exm(11,1,mat) = 0.999537929745094_rp  !0.999548606668578      !iF=20);
            vauin_exm(12,1,mat) = 0.199370793235339_rp  !0.999488774162635      !iS=21);
            vauin_exm(13,1,mat) = 0.0005168425303908772_rp !0.000512555980943569      !ap=22);
            vauin_exm(14,1,mat) = 0.999538048856435_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15,1,mat) = 0.226043556528992_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16,1,mat) = 0.000000002451190519566430_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17,1,mat) = 0.999999990351502_rp  !0.999999990696210      !ff=26);
            vauin_exm(18,1,mat) = 0.713425153962173_rp  !0.904906458666787      !fs=27);
            vauin_exm(19,1,mat) = 0.999999990357108_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20,1,mat) = 0.948912597173054_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21,1,mat) = 0.957191372519062_rp  !0.999903346883777      !jca=30);
            vauin_exm(22,1,mat) = 0.01457185320720413_rp !0.00215555277945401      !nca=31);
            vauin_exm(23,1,mat) = 0.999971016040121_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24,1,mat) = 0.999997537179088_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25,1,mat) = 0.01854399770068590_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26,1,mat) = 0.819260557541050_rp !0.487585264457487      !xrs=35);
            vauin_exm(27,1,mat) = 0.535044673733069_rp !0.276203479404767      !xs1=36);
            vauin_exm(28,1,mat) = 0.0006762955384707180_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29,1,mat) = 0.996815803757701_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9,1,mat) = 0.0000002397313133386565_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10,1,mat) = 0.0000002975619496753238_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11,1,mat) = 0.08074947458467882_rp !0.0228529042639590      !CaMKt=41);     
                        
            vminimate_exm(2,mat) = -87.5165638807769_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            ituss_exm = 2    !midmyocardial
            
            viclo_exm(1:26,1) = 0.0_rp      
            vcoin_exm(1,2,mat) = 0.0001302993381963230_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2,2,mat) = 7.78650116481851_rp !8.17212071948942       !nass=3);
            vcoin_exm(3,2,mat) = 143.967943434319_rp !143.675184333376       !ki=4);
            vcoin_exm(4,2,mat) = 143.967928583526_rp  !143.675147146842       !kss=5);
            vcoin_exm(5,2,mat) = 7.78627290424874_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6,2,mat) = 0.0001367372627130561_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7,2,mat) = 2.59618964681759_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8,2,mat) = 1.50264239254053_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1,2,mat) = 0.007711425764228608_rp !0.00739719746920272       !m=10);
            vauin_exm(2,2,mat) = 0.809593790171627_rp  !0.695621622011335       !hf=11);
            vauin_exm(3,2,mat) = 0.809668554567417_rp  !0.695601842634086       !hs=12);
            vauin_exm(4,2,mat) = 0.676310283476841_rp  !0.695486248719023       !j=13);
            vauin_exm(5,2,mat) = 0.609473347484487_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6,2,mat) = 0.574825632378327_rp  !0.695403157533235       !jp=15);
            vauin_exm(7,2,mat) = 0.0002064428815447646_rp !0.000190839777466418       !mL=16);
            vauin_exm(8,2,mat) = 0.246316732862273_rp  !0.493606704642336       !hL=17);
            vauin_exm(9,2,mat) = 0.103887440257557_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10,2,mat) = 0.001034493006853912_rp !0.00100594231451985      !a=19);
            vauin_exm(11,2,mat) = 0.999513339328340_rp  !0.999548606668578      !iF=20);
            vauin_exm(12,2,mat) = 0.151520865587831_rp  !0.999488774162635      !iS=21);
            vauin_exm(13,2,mat) = 0.0005271107526648554_rp !0.000512555980943569      !ap=22);
            vauin_exm(14,2,mat) = 0.999513525196830_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15,2,mat) = 0.167105642226589_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16,2,mat) = 0.000000002626133449123401_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17,2,mat) = 0.999999986012299_rp  !0.999999990696210      !ff=26);
            vauin_exm(18,2,mat) = 0.611718123469599_rp  !0.904906458666787      !fs=27);
            vauin_exm(19,2,mat) = 0.999999989545206_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20,2,mat) = 0.904565662128277_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21,2,mat) = 0.922416557341182_rp  !0.999903346883777      !jca=30);
            vauin_exm(22,2,mat) = 0.01810224279643378_rp !0.00215555277945401      !nca=31);
            vauin_exm(23,2,mat) = 0.999584611307330_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24,2,mat) = 0.999972928171199_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25,2,mat) = 0.05599708391131620_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26,2,mat) = 0.860854175657345_rp !0.487585264457487      !xrs=35);
            vauin_exm(27,2,mat) = 0.619347616204376_rp !0.276203479404767      !xs1=36);
            vauin_exm(28,2,mat) = 0.003773361894746919_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29,2,mat) = 0.996893681111765_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9,2,mat) = 0.0000006261445113237194_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10,2,mat) = 0.0000007713396255442226_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11,2,mat) = 0.133755388160624_rp !0.0228529042639590      !CaMKt=41);     

       else if((moneclmate_exm(2,mat) == 600) .and. (kfl_hfmod_exm(mat) == 2_ip) .and.(kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_inaga_exm(mat)==1_ip)) then  !if it's a Normal MALE HUMAN  AT 600BCL
                    vminimate_exm(3,mat) = -87.7717_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
                    !ituss_exm = 3  !epicardial
                    !EPICARDIUM
                    vcoin_exm(1,3,mat) = 0.0001096_rp       !cai=6);
                    vcoin_exm(2,3,mat) = 8.8338_rp       !nass=3);
                    vcoin_exm(3,3,mat) = 142.6926_rp       !ki=4);
                    vcoin_exm(4,3,mat) = 142.6925_rp       !kss=5);
                    vcoin_exm(5,3,mat) = 8.8337_rp      ! nai=2
                    vcoin_exm(6,3,mat) = 0.00010776_rp       !cass=7);
                    vcoin_exm(7,3,mat) = 3.6249_rp       !cansr=8);
                    vcoin_exm(8,3,mat) = 3.3755_rp      !cajsr=9);
                    vauin_exm(1,3,mat) = 0.0075_rp       !m=10);
                    vauin_exm(2,3,mat) = 0.6901_rp       !hf=11);
                    vauin_exm(3,3,mat) = 0.69_rp       !hs=12);
                    vauin_exm(4,3,mat) = 0.6898_rp       !j=13);
                    vauin_exm(5,3,mat) = 0.4455_rp       !hsp=14);
                    vauin_exm(6,3,mat) = 0.6896_rp       !jp=15);
                    vauin_exm(7,3,mat) = 0.00019668_rp       !mL=16);
                    vauin_exm(8,3,mat) = 0.4566_rp       !hL=17);
                    vauin_exm(9,3,mat) = 0.2334_rp       !hLp=18);
                    vauin_exm(10,3,mat) = 0.001_rp      !a=19);
                    vauin_exm(11,3,mat) = 0.9995_rp      !iF=20);
                    vauin_exm(12,3,mat) = 0.9973_rp      !iS=21);
                    vauin_exm(13,3,mat) = 0.00051807_rp      !ap=22);
                    vauin_exm(14,3,mat) = 0.9995_rp      !iFp=23);
                    vauin_exm(15,3,mat) = 0.9986_rp      !iSp=24);
                    vauin_exm(16,3,mat) = 0.0000000024717_rp      !d=25);
                    vauin_exm(17,3,mat) = 1_rp      !ff=26);
                    vauin_exm(18,3,mat) = 0.8579_rp      !fs=27);
                    vauin_exm(19,3,mat) = 1_rp      !fcaf=28);
                    vauin_exm(20,3,mat) = 0.9957_rp      !fcas=29);
                    vauin_exm(21,3,mat) = 0.9978_rp      !jca=30);
                    vauin_exm(22,3,mat) = 0.0068_rp      !nca=31);
                    vauin_exm(23,3,mat) = 1_rp      !ffp=32);
                    vauin_exm(24,3,mat) = 1_rp      !fcafp=33);
                    vauin_exm(25,3,mat) = 0.00010734_rp     !xrf=34);
                    vauin_exm(26,3,mat) = 0.6305_rp      !xrs=35);
                    vauin_exm(27,3,mat) = 0.3536_rp      !xs1=36);
                    vauin_exm(28,3,mat) = 0.00019802_rp     !xs2=37);
                    vauin_exm(29,3,mat) = 0.9968_rp      !xk1=38);
                    vcoin_exm(9,3,mat) = 0.00000053112_rp       !Jrelnp=39);
                    vcoin_exm(10,3,mat) = 0.00000066389_rp      !Jrelp=40);
                    vcoin_exm(11,3,mat) = 0.0481_rp      !CaMKt=41);

                    !ENDOCARDIUM
                    vminimate_exm(1,mat) = -87.9538_rp
                    !ituss_exm = 1  !endocardial
                    vcoin_exm(1,1,mat) = 0.000107_rp       !cai=6);
                    vcoin_exm(2,1,mat) = 8.0504_rp       !nass=3);
                    vcoin_exm(3,1,mat) = 143.7827_rp       !ki=4);
                    vcoin_exm(4,1,mat) = 143.7826_rp       !kss=5);
                    vcoin_exm(5,1,mat) = 8.0503_rp     ! nai=2
                    vcoin_exm(6,1,mat) = 0.00010827_rp       !cass=7);
                    vcoin_exm(7,1,mat) = 2.0238_rp       !cansr=8);
                    vcoin_exm(8,1,mat) = 1.7311_rp      !cajsr=9);
                    vauin_exm(1,1,mat) = 0.0074_rp       !m=10);
                    vauin_exm(2,1,mat) = 0.6964_rp       !hf=11);
                    vauin_exm(3,1,mat) = 0.6963_rp       !hs=12);
                    vauin_exm(4,1,mat) = 0.6958_rp       !j=13);
                    vauin_exm(5,1,mat) = 0.4527_rp       !hsp=14);
                    vauin_exm(6,1,mat) = 0.6952_rp       !jp=15);
                    vauin_exm(7,1,mat) = 0.00018999_rp       !mL=16);
                    vauin_exm(8,1,mat) = 0.4425_rp       !hL=17);
                    vauin_exm(9,1,mat) = 0.2167_rp       !hLp=18);
                    vauin_exm(10,1,mat) = 0.001_rp      !a=19);
                    vauin_exm(11,1,mat) = 0.9996_rp      !iF=20);
                    vauin_exm(12,1,mat) = 0.3632_rp      !iS=21);
                    vauin_exm(13,1,mat) = 0.00051176_rp      !ap=22);
                    vauin_exm(14,1,mat) = 0.9996_rp      !iFp=23);
                    vauin_exm(15,1,mat) = 0.4086_rp      !iSp=24);
                    vauin_exm(16,1,mat) = 0.0000000023677_rp      !d=25);
                    vauin_exm(17,1,mat) = 1_rp      !ff=26);
                    vauin_exm(18,1,mat) = 0.833_rp      !fs=27);
                    vauin_exm(19,1,mat) = 1_rp      !fcaf=28);
                    vauin_exm(20,1,mat) = 0.833_rp      !fcas=29);
                    vauin_exm(21,1,mat) = 0.9967_rp      !jca=30);
                    vauin_exm(22,1,mat) = 0.007_rp      !nca=31);
                    vauin_exm(23,1,mat) = 1_rp      !ffp=32);
                    vauin_exm(24,1,mat) = 1_rp      !fcafp=33);
                    vauin_exm(25,1,mat) = 0.00024763_rp      !xrf=34);
                    vauin_exm(26,1,mat) = 0.6669_rp      !xrs=35);
                    vauin_exm(27,1,mat) = 0.3946_rp      !xs1=36);
                    vauin_exm(28,1,mat) = 0.00019447_rp      !xs2=37);
                    vauin_exm(29,1,mat) = 0.9968_rp      !xk1=38);
                    vcoin_exm(9,1,mat) =  0.00000031607_rp      !Jrelnp=39);
                    vcoin_exm(10,1,mat) = 0.00000039438_rp      !Jrelp=40);
                    vcoin_exm(11,1,mat) = 0.0407_rp      !CaMKt=41);

                    !MIDMYOCARDIUM
                    vminimate_exm(2,mat) = -87.5215_rp
                    !ituss_exm = 2  !midmyocardial
                    vcoin_exm(1,2,mat) =  0.00013215_rp       !cai=6);
                    vcoin_exm(2,2,mat) =  10.0687_rp      !nass=3);
                    vcoin_exm(3,2,mat) =  141.4514_rp      !ki=4);
                    vcoin_exm(4,2,mat) =  141.4514_rp      !kss=5);
                    vcoin_exm(5,2,mat) =  10.0684_rp     ! nai=2
                    vcoin_exm(6,2,mat) = 0.00013224_rp       !cass=7);
                    vcoin_exm(7,2,mat) = 4.0357_rp       !cansr=8);
                    vcoin_exm(8,2,mat) = 3.6872_rp       !cajsr=9);
                    vauin_exm(1,2,mat) = 0.0077_rp       !m=10);
                    vauin_exm(2,2,mat) = 0.6812_rp       !hf=11);
                    vauin_exm(3,2,mat) = 0.6811_rp       !hs=12);
                    vauin_exm(4,2,mat) = 0.68_rp       !j=13);
                    vauin_exm(5,2,mat) = 0.4349_rp       !hsp=14);
                    vauin_exm(6,2,mat) = 0.6786_rp       !jp=15);
                    vauin_exm(7,2,mat) = 0.00020625_rp       !mL=16);
                    vauin_exm(8,2,mat) = 0.3975_rp       !hL=17);
                    vauin_exm(9,2,mat) = 0.1798_rp       !hLp=18);
                    vauin_exm(10,2,mat) = 0.001_rp      !a=19);
                    vauin_exm(11,2,mat) = 0.9995_rp      !iF=20);
                    vauin_exm(12,2,mat) = 0.3166_rp      !iS=21);
                    vauin_exm(13,2,mat) = 0.0005269_rp      !ap=22);
                    vauin_exm(14,2,mat) = 0.9995_rp      !iFp=23);
                    vauin_exm(15,2,mat) = 0.3501_rp      !iSp=24);
                    vauin_exm(16,2,mat) = 0.0000000026226_rp      !d=25);
                    vauin_exm(17,2,mat) = 1_rp      !ff=26);
                    vauin_exm(18,2,mat) = 0.769_rp      !fs=27);
                    vauin_exm(19,2,mat) = 1_rp      !fcaf=28);
                    vauin_exm(20,2,mat) = 0.9864_rp      !fcas=29);
                    vauin_exm(21,2,mat) = 0.994_rp      !jca=30);
                    vauin_exm(22,2,mat) = 0.0148_rp      !nca=31);
                    vauin_exm(23,2,mat) = 1_rp      !ffp=32);
                    vauin_exm(24,2,mat) = 1_rp      !fcafp=33);
                    vauin_exm(25,2,mat) = 0.0009926_rp      !xrf=34);
                    vauin_exm(26,2,mat) = 0.7143_rp      !xrs=35);
                    vauin_exm(27,2,mat) = 0.4686_rp      !xs1=36);
                    vauin_exm(28,2,mat) = 0.0002076_rp      !xs2=37);
                    vauin_exm(29,2,mat) = 0.9969_rp      !xk1=38);
                    vcoin_exm(9,2,mat) =  0.0000019734_rp      !Jrelnp=39);
                    vcoin_exm(10,2,mat) = 0.0000024672_rp      !Jrelp=40);
                    vcoin_exm(11,2,mat) = 0.1839_rp       !CaMKt=41);

       else if((moneclmate_exm(2,mat) == 600) .and. (kfl_hfmod_exm(mat) == 3_ip) .and.(kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_inaga_exm(mat)==1_ip)) then  !if it's a Normal FEMALE HUMAN at 600BCL
                    vminimate_exm(3,mat) = -87.5541_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
                    !ituss_exm = 3  !epicardial
                    !EPICARDIUM
                    vcoin_exm(1,3,mat) = 0.00012597_rp       !cai=6);
                    vcoin_exm(2,3,mat) = 9.3243_rp       !nass=3);
                    vcoin_exm(3,3,mat) = 141.9779_rp       !ki=4);
                    vcoin_exm(4,3,mat) = 141.9778_rp       !kss=5);
                    vcoin_exm(5,3,mat) = 9.3242_rp      ! nai=2
                    vcoin_exm(6,3,mat) = 0.00012154_rp       !cass=7);
                    vcoin_exm(7,3,mat) = 5.2532_rp       !cansr=8);
                    vcoin_exm(8,3,mat) = 5.0894_rp      !cajsr=9);
                    vauin_exm(1,3,mat) = 0.0077_rp       !m=10);
                    vauin_exm(2,3,mat) = 0.6824_rp       !hf=11);
                    vauin_exm(3,3,mat) = 0.6824_rp       !hs=12);
                    vauin_exm(4,3,mat) = 0.6823_rp       !j=13);
                    vauin_exm(5,3,mat) = 0.4368_rp       !hsp=14);
                    vauin_exm(6,3,mat) = 0.6821_rp       !jp=15);
                    vauin_exm(7,3,mat) = 0.00020497_rp       !mL=16);
                    vauin_exm(8,3,mat) = 0.4403_rp       !hL=17);
                    vauin_exm(9,3,mat) = 0.2178_rp       !hLp=18);
                    vauin_exm(10,3,mat) = 0.001_rp      !a=19);
                    vauin_exm(11,3,mat) = 0.9995_rp      !iF=20);
                    vauin_exm(12,3,mat) = 0.9957_rp      !iS=21);
                    vauin_exm(13,3,mat) = 0.00052572_rp      !ap=22);
                    vauin_exm(14,3,mat) = 0.9995_rp      !iFp=23);
                    vauin_exm(15,3,mat) = 0.9978_rp      !iSp=24);
                    vauin_exm(16,3,mat) = 0.0000000026021_rp      !d=25);
                    vauin_exm(17,3,mat) = 1_rp      !ff=26);
                    vauin_exm(18,3,mat) = 0.8391_rp      !fs=27);
                    vauin_exm(19,3,mat) = 1_rp      !fcaf=28);
                    vauin_exm(20,3,mat) = 0.9941_rp      !fcas=29);
                    vauin_exm(21,3,mat) = 0.9973_rp      !jca=30);
                    vauin_exm(22,3,mat) = 0.0107_rp      !nca=31);
                    vauin_exm(23,3,mat) = 1_rp      !ffp=32);
                    vauin_exm(24,3,mat) = 1_rp      !fcafp=33);
                    vauin_exm(25,3,mat) = 0.00019016_rp     !xrf=34);
                    vauin_exm(26,3,mat) = 0.6589_rp      !xrs=35);
                    vauin_exm(27,3,mat) = 0.3798_rp      !xs1=36);
                    vauin_exm(28,3,mat) = 0.00020289_rp     !xs2=37);
                    vauin_exm(29,3,mat) = 0.9969_rp      !xk1=38);
                    vcoin_exm(9,3,mat) = 0.000000563_rp       !Jrelnp=39);
                    vcoin_exm(10,3,mat) = 0.00000070371_rp      !Jrelp=40);
                    vcoin_exm(11,3,mat) = 0.0771_rp      !CaMKt=41);

                    !ENDOCARDIUM
                    vminimate_exm(1,mat) = -87.8796_rp
                    !ituss_exm = 1  !endocardial
                    vcoin_exm(1,1,mat) = 0.00010362_rp       !cai=6);
                    vcoin_exm(2,1,mat) = 8.4849_rp       !nass=3);
                    vcoin_exm(3,1,mat) = 143.335_rp       !ki=4);
                    vcoin_exm(4,1,mat) = 143.335_rp       !kss=5);
                    vcoin_exm(5,1,mat) = 8.4848_rp     ! nai=2
                    vcoin_exm(6,1,mat) = 0.00010455_rp       !cass=7);
                    vcoin_exm(7,1,mat) = 2.4087_rp       !cansr=8);
                    vcoin_exm(8,1,mat) = 2.0558_rp      !cajsr=9);
                    vauin_exm(1,1,mat) = 0.0074_rp       !m=10);
                    vauin_exm(2,1,mat) = 0.6939_rp       !hf=11);
                    vauin_exm(3,1,mat) = 0.6938_rp       !hs=12);
                    vauin_exm(4,1,mat) = 0.6932_rp       !j=13);
                    vauin_exm(5,1,mat) = 0.4497_rp       !hsp=14);
                    vauin_exm(6,1,mat) = 0.6923_rp       !jp=15);
                    vauin_exm(7,1,mat) = 0.00019267_rp       !mL=16);
                    vauin_exm(8,1,mat) = 0.422_rp       !hL=17);
                    vauin_exm(9,1,mat) = 0.1981_rp       !hLp=18);
                    vauin_exm(10,1,mat) = 0.001_rp      !a=19);
                    vauin_exm(11,1,mat) = 0.9995_rp      !iF=20);
                    vauin_exm(12,1,mat) = 0.3378_rp      !iS=21);
                    vauin_exm(13,1,mat) = 0.00051432_rp      !ap=22);
                    vauin_exm(14,1,mat) = 0.9995_rp      !iFp=23);
                    vauin_exm(15,1,mat) = 0.3769_rp      !iSp=24);
                    vauin_exm(16,1,mat) = 0.0000000024096_rp      !d=25);
                    vauin_exm(17,1,mat) = 1_rp      !ff=26);
                    vauin_exm(18,1,mat) = 0.8081_rp      !fs=27);
                    vauin_exm(19,1,mat) = 1_rp      !fcaf=28);
                    vauin_exm(20,1,mat) = 0.9906_rp      !fcas=29);
                    vauin_exm(21,1,mat) = 0.9955_rp      !jca=30);
                    vauin_exm(22,1,mat) = 0.0061_rp      !nca=31);
                    vauin_exm(23,1,mat) = 1_rp      !ffp=32);
                    vauin_exm(24,1,mat) = 1_rp      !fcafp=33);
                    vauin_exm(25,1,mat) = 0.00053579_rp      !xrf=34);
                    vauin_exm(26,1,mat) = 0.6919_rp      !xrs=35);
                    vauin_exm(27,1,mat) = 0.4271_rp      !xs1=36);
                    vauin_exm(28,1,mat) = 0.00019696_rp      !xs2=37);
                    vauin_exm(29,1,mat) = 0.9968_rp      !xk1=38);
                    vcoin_exm(9,1,mat) =  0.00000039242_rp      !Jrelnp=39);
                    vcoin_exm(10,1,mat) = 0.00000049021_rp      !Jrelp=40);
                    vcoin_exm(11,1,mat) = 0.0502_rp      !CaMKt=41);

                    !MIDMYOCARDIUM
                    vminimate_exm(2,mat) = -87.5541_rp
                    !ituss_exm = 2  !midmyocardial
                    vcoin_exm(1,2,mat) =  0.00012597_rp       !cai=6);
                    vcoin_exm(2,2,mat) =  9.3243_rp      !nass=3);
                    vcoin_exm(3,2,mat) =  141.9779_rp      !ki=4);
                    vcoin_exm(4,2,mat) =  141.9778_rp      !kss=5);
                    vcoin_exm(5,2,mat) =  9.3242_rp     ! nai=2
                    vcoin_exm(6,2,mat) = 0.00012154_rp       !cass=7);
                    vcoin_exm(7,2,mat) = 5.2532_rp       !cansr=8);
                    vcoin_exm(8,2,mat) = 5.0894_rp       !cajsr=9);
                    vauin_exm(1,2,mat) = 0.0077_rp       !m=10);
                    vauin_exm(2,2,mat) = 0.6824_rp       !hf=11);
                    vauin_exm(3,2,mat) = 0.6824_rp       !hs=12);
                    vauin_exm(4,2,mat) = 0.6823_rp       !j=13);
                    vauin_exm(5,2,mat) = 0.4368_rp       !hsp=14);
                    vauin_exm(6,2,mat) = 0.6821_rp       !jp=15);
                    vauin_exm(7,2,mat) = 0.00020497_rp       !mL=16);
                    vauin_exm(8,2,mat) = 0.4403_rp       !hL=17);
                    vauin_exm(9,2,mat) = 0.2178_rp       !hLp=18);
                    vauin_exm(10,2,mat) = 0.001_rp      !a=19);
                    vauin_exm(11,2,mat) = 0.9995_rp      !iF=20);
                    vauin_exm(12,2,mat) = 0.9957_rp      !iS=21);
                    vauin_exm(13,2,mat) = 0.00052572_rp      !ap=22);
                    vauin_exm(14,2,mat) = 0.9995_rp      !iFp=23);
                    vauin_exm(15,2,mat) = 0.9978_rp      !iSp=24);
                    vauin_exm(16,2,mat) = 0.0000000026021_rp      !d=25);
                    vauin_exm(17,2,mat) = 1_rp      !ff=26);
                    vauin_exm(18,2,mat) = 0.8391_rp      !fs=27);
                    vauin_exm(19,2,mat) = 1_rp      !fcaf=28);
                    vauin_exm(20,2,mat) = 0.9941_rp      !fcas=29);
                    vauin_exm(21,2,mat) = 0.9973_rp      !jca=30);
                    vauin_exm(22,2,mat) = 0.0107_rp      !nca=31);
                    vauin_exm(23,2,mat) = 1_rp      !ffp=32);
                    vauin_exm(24,2,mat) = 1_rp      !fcafp=33);
                    vauin_exm(25,2,mat) = 0.00019016_rp      !xrf=34);
                    vauin_exm(26,2,mat) = 0.6589_rp      !xrs=35);
                    vauin_exm(27,2,mat) = 0.3798_rp      !xs1=36);
                    vauin_exm(28,2,mat) = 0.00020289_rp      !xs2=37);
                    vauin_exm(29,2,mat) = 0.9969_rp      !xk1=38);
                    vcoin_exm(9,2,mat) =  0.000000563_rp      !Jrelnp=39);
                    vcoin_exm(10,2,mat) = 0.00000070371_rp      !Jrelp=40);
                    vcoin_exm(11,2,mat) = 0.0771_rp       !CaMKt=41);       
    
       elseif ((moneclmate_exm(2,mat) /= 1000) .and. (kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_hfmod_exm(mat) /= 4_ip) ) then
           !.or. (kfl_drugsmate_exm(mat) == 1) Normal with different heart rate .and. (moneclmate_exm(2,mat) /= 857).and. (moneclmate_exm(2,mat) /= 600
            
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
    
            ituss_exm = 3_ip    !epicardial
            viclo_exm(1:26,1) = 0.0_rp      
            vcolo_exm(1,1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2,1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3,1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4,1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5,1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6,1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7,1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8,1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1,1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2,1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3,1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4,1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5,1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6,1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7,1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8,1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9,1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10,1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11,1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12,1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13,1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14,1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15,1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16,1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17,1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18,1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19,1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20,1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21,1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22,1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23,1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24,1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25,1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26,1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27,1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28,1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29,1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9,1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10,1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11,1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);     
            if (exm_oneohr) then
                exm_oneohr = exm_oceohr(mat)
            end if
                  
            !endocardium
            ituss_exm = 1_ip  !endocardial
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26,1) = 0.0_rp      
            vcolo_exm(1,1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2,1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3,1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4,1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5,1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6,1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7,1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8,1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1,1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2,1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3,1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4,1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5,1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6,1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7,1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8,1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9,1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10,1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11,1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12,1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13,1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14,1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15,1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16,1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17,1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18,1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19,1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20,1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21,1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22,1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23,1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24,1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25,1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26,1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27,1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28,1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29,1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9,1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10,1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11,1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);    
            if (exm_oneohr) then
                exm_oneohr = exm_oceohr(mat)
            end if
    
            ituss_exm = 2_ip  !MIDmyocardial
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26,1) = 0.0_rp      
            vcolo_exm(1,1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2,1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3,1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4,1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5,1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6,1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7,1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8,1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1,1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2,1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3,1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4,1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5,1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6,1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7,1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8,1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9,1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10,1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11,1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12,1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13,1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14,1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15,1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16,1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17,1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18,1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19,1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20,1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21,1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22,1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23,1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24,1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25,1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26,1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27,1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28,1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29,1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9,1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10,1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11,1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);    
            if (exm_oneohr) then
                exm_oneohr = exm_oceohr(mat)
            end if

       end if      

  else if (kfl_hfmodmate_exm(mat) == 1_ip) then  !if it's Modified and Initial conditions need to be calculated                    
     
    if ((kfl_hfmod_exm(mat) /= 4_ip) )  then      !.and. (moneclmate_exm(2,mat) /= 857) .and. (moneclmate_exm(2,mat) /= 1000).or. (kfl_drugsmate_exm(mat) == 1)Normal with different heart rate
            
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
    
            ituss_exm = 3_ip    !epicardial
            viclo_exm(1:26,1) = 0.0_rp      
            vcolo_exm(1,1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2,1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3,1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4,1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5,1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6,1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7,1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8,1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1,1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2,1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3,1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4,1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5,1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6,1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7,1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8,1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9,1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10,1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11,1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12,1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13,1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14,1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15,1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16,1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17,1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18,1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19,1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20,1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21,1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22,1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23,1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24,1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25,1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26,1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27,1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28,1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29,1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9,1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10,1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11,1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);     
            if (exm_oneohr) then
                exm_oneohr = exm_oceohr(mat)
            end if
                  
            !endocardium
            ituss_exm = 1_ip  !endocardial
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26,1) = 0.0_rp      
            vcolo_exm(1,1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2,1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3,1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4,1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5,1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6,1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7,1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8,1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1,1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2,1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3,1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4,1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5,1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6,1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7,1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8,1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9,1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10,1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11,1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12,1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13,1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14,1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15,1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16,1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17,1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18,1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19,1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20,1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21,1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22,1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23,1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24,1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25,1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26,1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27,1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28,1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29,1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9,1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10,1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11,1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);    
            if (exm_oneohr) then
                exm_oneohr = exm_oceohr(mat)
            end if
    
            ituss_exm = 2_ip  !MIDmyocardial
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26,1) = 0.0_rp      
            vcolo_exm(1,1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2,1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3,1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4,1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5,1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6,1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7,1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8,1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1,1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2,1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3,1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4,1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5,1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6,1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7,1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8,1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9,1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10,1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11,1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12,1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13,1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14,1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15,1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16,1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17,1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18,1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19,1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20,1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21,1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22,1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23,1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24,1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25,1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26,1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27,1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28,1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29,1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9,1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10,1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11,1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);    
            if (exm_oneohr) then
                exm_oneohr = exm_oceohr(mat)
            end if

    else if ((kfl_hfmod_exm(mat) == 4_ip) .and. (kfl_inaga_exm(mat) == 1_ip)) then  !if it's a Modified PIG @ 400 BCL with INa modified channel from Passini
       !if( (kfl_drugsmate_exm(mat) == 0) then   !no drugs
        
            elmlo_exm(1:2) = -87.8061985708231_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL

            ituss_exm = 3_ip    !epicardial
            viclo_exm(1:26,1) = 0.0_rp      
            vcolo_exm(1,1:3) = 0.0001040361496276368_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2,1:3) =  7.43617875723091_rp !8.17212071948942       !nass=3);
            vcolo_exm(3,1:3) = 144.311817237457_rp !143.675184333376       !ki=4);
            vcolo_exm(4,1:3) = 144.311793923281_rp  !143.675147146842       !kss=5);
            vcolo_exm(5,1:3) = 7.43604263653445_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6,1:3) = 0.0001102064090489676_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7,1:3) = 2.63411603951980_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8,1:3) = 1.76343231449702_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1,1:2) = 0.007490096291770637_rp !0.00739719746920272       !m=10);
            vaulo_exm(2,1:2) = 0.816838169651186_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3,1:2) = 0.816869319742196_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4,1:2) = 0.779416631991491_rp  !0.695486248719023       !j=13);
            vaulo_exm(5,1:2) = 0.621614033150460_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6,1:2) = 0.719030208978245_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7,1:2) = 0.0001953921132645023_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8,1:2) = 0.350601487895652_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9,1:2) = 0.165483509753837_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10,1:2) = 0.001014440468937102_rp !0.00100594231451985      !a=19);
            vaulo_exm(11,1:2) = 0.999538652553908_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12,1:2) = 0.940646485484524_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13,1:2) = 0.0005168881869164542_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14,1:2) = 0.999538658574815_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15,1:2) = 0.961396539930704_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16,1:2) = 0.000000002451977924935544_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17,1:2) = 0.999999990358952_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18,1:2) = 0.721766724570291_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19,1:2) = 0.999999990360122_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20,1:2) = 0.952750811437167_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21,1:2) = 0.962655802194759_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22,1:2) = 0.007755537033307213_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23,1:2) = 0.999989171225269_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24,1:2) = 0.999998595329089_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25,1:2) = 0.01177184211893142_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26,1:2) = 0.803851729756987_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27,1:2) = 0.525297568917089_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28,1:2) = 0.0004266215773384047_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29,1:2) = 0.996814746325558_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9,1:3) = 0.0000003867443373607881_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10,1:3) = 0.0000004808818699484394_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11,1:3) = 0.08284541646773119_rp !0.0228529042639590      !CaMKt=41);     
            if (exm_oneohr) then
                exm_oneohr = exm_oceohr(mat)
            end if
                        
            elmlo_exm(1:2) = -87.8077947685269_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            ituss_exm = 1_ip    !endocardium
            
            viclo_exm(1:26,1) = 0.0_rp      
            vcolo_exm(1,1:3) = 0.0001244664562667345_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2,1:3) = 7.38274584190492_rp !8.17212071948942       !nass=3);
            vcolo_exm(3,1:3) = 144.434867863316_rp !143.675184333376       !ki=4);
            vcolo_exm(4,1:3) = 144.434843290996_rp  !143.675147146842       !kss=5);
            vcolo_exm(5,1:3) = 7.38257910315418_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6,1:3) = 0.0001303227834520188_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7,1:3) = 2.22033410012224_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8,1:3) = 1.57032538537506_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1,1:2) = 0.007488900815703561_rp !0.00739719746920272       !m=10);
            vaulo_exm(2,1:2) = 0.816818226826766_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3,1:2) = 0.816860080147452_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4,1:2) = 0.760420946701455_rp  !0.695486248719023       !j=13);
            vaulo_exm(5,1:2) = 0.621427482142547_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6,1:2) = 0.686823279932452_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7,1:2) = 0.0001953332078784692_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8,1:2) = 0.321919485372494_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9,1:2) = 0.146893625330722_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10,1:2) = 0.001014350908419775_rp !0.00100594231451985      !a=19);
            vaulo_exm(11,1:2) = 0.999537929745094_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12,1:2) = 0.199370793235339_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13,1:2) = 0.0005168425303908772_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14,1:2) = 0.999538048856435_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15,1:2) = 0.226043556528992_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16,1:2) = 0.000000002451190519566430_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17,1:2) = 0.999999990351502_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18,1:2) = 0.713425153962173_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19,1:2) = 0.999999990357108_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20,1:2) = 0.948912597173054_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21,1:2) = 0.957191372519062_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22,1:2) = 0.01457185320720413_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23,1:2) = 0.999971016040121_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24,1:2) = 0.999997537179088_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25,1:2) = 0.01854399770068590_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26,1:2) = 0.819260557541050_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27,1:2) = 0.535044673733069_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28,1:2) = 0.0006762955384707180_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29,1:2) = 0.996815803757701_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9,1:3) = 2.397313133386565E-007_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10,1:3) = 2.975619496753238E-007_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11,1:3) = 0.08074947458467882_rp !0.0228529042639590      !CaMKt=41);     
            if (exm_oneohr) then
                exm_oneohr = exm_oceohr(mat)
            end if
                        
            elmlo_exm(1:2) = -87.5165638807769_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            ituss_exm = 2_ip    !midmyocardial
            
            viclo_exm(1:26,1) = 0.0_rp      
            vcolo_exm(1,1:3) = 0.0001302993381963230_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2,1:3) = 7.78650116481851_rp !8.17212071948942       !nass=3);
            vcolo_exm(3,1:3) = 143.967943434319_rp !143.675184333376       !ki=4);
            vcolo_exm(4,1:3) = 143.967928583526_rp  !143.675147146842       !kss=5);
            vcolo_exm(5,1:3) = 7.78627290424874_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6,1:3) = 0.0001367372627130561_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7,1:3) = 2.59618964681759_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8,1:3) = 1.50264239254053_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1,1:2) = 0.007711425764228608_rp !0.00739719746920272       !m=10);
            vaulo_exm(2,1:2) = 0.809593790171627_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3,1:2) = 0.809668554567417_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4,1:2) = 0.676310283476841_rp  !0.695486248719023       !j=13);
            vaulo_exm(5,1:2) = 0.609473347484487_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6,1:2) = 0.574825632378327_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7,1:2) = 0.0002064428815447646_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8,1:2) = 0.246316732862273_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9,1:2) = 0.103887440257557_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10,1:2) = 0.001034493006853912_rp !0.00100594231451985      !a=19);
            vaulo_exm(11,1:2) = 0.999513339328340_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12,1:2) = 0.151520865587831_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13,1:2) = 0.0005271107526648554_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14,1:2) = 0.999513525196830_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15,1:2) = 0.167105642226589_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16,1:2) = 0.000000002626133449123401_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17,1:2) = 0.999999986012299_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18,1:2) = 0.611718123469599_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19,1:2) = 0.999999989545206_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20,1:2) = 0.904565662128277_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21,1:2) = 0.922416557341182_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22,1:2) = 0.01810224279643378_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23,1:2) = 0.999584611307330_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24,1:2) = 0.999972928171199_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25,1:2) = 0.05599708391131620_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26,1:2) = 0.860854175657345_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27,1:2) = 0.619347616204376_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28,1:2) = 0.003773361894746919_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29,1:2) = 0.996893681111765_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9,1:3) = 0.0000006261445113237194_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10,1:3) = 0.0000007713396255442226_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11,1:3) = 0.133755388160624_rp !0.0228529042639590      !CaMKt=41);     
            if (exm_oneohr) then
                exm_oneohr = exm_oceohr(mat)
            end if
    end if
  end if

end function exm_oneohr
