!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exm_oneohr.f90
!> @author  mixed
!> @date    2019-11-16
!> @brief   mod_exm_oneohr
!> @details mod_exm_oneohr
!> @}
!-----------------------------------------------------------------------
module mod_exm_oneohr

   use def_master
   use def_exmedi
   use mod_messages, only: messages_live
   use mod_exm_oceohr,      only : exm_oceohr, exm_oceola
   use mod_exm_sld_eccoupling, only :  has_exmsld_coupling


   implicit none

   integer(ip), parameter  :: &
      EXM_ONEOHR_COUPLED_LAND = 0_ip, &
      !precalculated hardcoded initial conditions
      EXM_ONEOHR_NORMAL_1000BCL = 1_ip, &             !ohara-rudy,    MYOCYTE=normal,        CYCLELENGTH=1000 (genderless)
      EXM_ONEOHR_NORMAL_70BPM = 2_ip, &               !ohara-passini, MYOCYTE=normal,        CYCLELENGTH=857
      EXM_ONEOHR_NORMAL_PIG_400BCL = 3_ip, &          !ohara-passini, MYOCYTE=pig,           CYCLELENGTH=400
      EXM_ONEOHR_NORMAL_HUMAN_MALE_600BCL = 4_ip, &   !ohara-passini, MYOCYTE=normal male,   CYCLELENGTH=600
      EXM_ONEOHR_NORMAL_HUMAN_FEMALE_600BCL = 5_ip, & !ohara-passini, MYOCYTE=normal female, CYCLELENGTH=600
      !forces recalculation of the initial conditions
      EXM_ONEOHR_PIG_NODRUGS = 6_ip, & !pig
      EXM_ONEOHR_NORMAL_MODIFIED = 7_ip, & !modified
      EXM_ONEOHR_MODPIG_400BCL_MODINA = 8_ip !modified pig, same as EXM_ONEOHR_NORMAL_PIG_400BCL with modified Ina

   integer(ip), parameter  :: &
      EXM_ONEOHR_ERROR_CONVERGED = 0_ip, &
      EXM_ONEOHR_ERROR_NOTCONVERGED = 1_ip, &
      EXM_ONEOHR_ERROR_NOTINITIALIZED = 2_ip

   public :: exm_oneohr
   private :: exm_init_voltages, localvoltages2globalvolatges

contains

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
! icelltype -- epi, endo, mid (1-3)
! Returns:
! success_status == 0_ip -- sucess, >0 -- error, see EXM_ONEOHR_ERROR* for the meaning
! land_variables:
!     land_variables(1) = S
!     land_variables(2) = W
!     land_variables(3) = CaTRPN
!     land_variables(4) = B
!     land_variables(5) = zeta_s
!     land_variables(6) = zeta_w
   subroutine exm_oneohr(mat, icelltype, land_variables, success_status, ohara_stats)

      implicit none
      integer(ip), intent(in) :: mat, icelltype
      integer(ip), intent(out) :: success_status
      integer(ip)              :: success_status_tmp
      logical     :: flag_land, initialized
      real(rp), dimension(1:6), intent(out) ::land_variables

      real(rp), intent(inout)    :: ohara_stats(3_ip)  !saves ohara stats: number of beats, tolerance

      success_status = 0_ip

      flag_land = .FALSE.
      initialized = .FALSE. !flag to check if any of the IF was executed

      if (kfl_eccty(mat) == 3_ip .or. kfl_eccty(mat) == 4_ip) flag_land = .TRUE.

      if ( has_exmsld_coupling() .and. flag_land) then

         call exm_init_voltages(EXM_ONEOHR_COUPLED_LAND, mat, icelltype)

         if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
            call exm_oceola(mat, icelltype, land_variables, success_status_tmp, ohara_stats)
            success_status = success_status + success_status_tmp
         end if

         call localvoltages2globalvolatges(mat, icelltype)
         initialized = .TRUE.
      end if

      if (kfl_hfmodmate_exm(mat) == 0_ip) then

         if ((moneclmate_exm(2, mat) == 1000) .and. (kfl_drugsmate_exm(mat) == 0_ip)) then  !if it's Normal and at 1000 bpm
            call exm_init_voltages(EXM_ONEOHR_NORMAL_1000BCL, mat, icelltype)
            initialized = .TRUE.

         else if ((moneclmate_exm(2, mat) == 857) .and. (kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_inaga_exm(mat) == 1_ip)) then   !Normal at 70 bpm
            call exm_init_voltages(EXM_ONEOHR_NORMAL_70BPM, mat, icelltype)
            initialized = .TRUE.

            !else if (kfl_hfmod_exm(mat) == 1 .and. kfl_hfmodmate_exm(mat) == 0 ) then  !if it's a Normal PIG
        else if((moneclmate_exm(2,mat) == 400) .and. (kfl_hfmod_exm(mat) == 1_ip) .and. (kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_inaga_exm(mat) == 1_ip)) then   !Normal PIG AT 400 BCL
            call exm_init_voltages(EXM_ONEOHR_NORMAL_PIG_400BCL, mat, icelltype)
            initialized = .TRUE.

        else if((moneclmate_exm(2,mat) == 600) .and. (kfl_hfmod_exm(mat) == 2_ip) .and.(kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_inaga_exm(mat)==1_ip)) then  !if it's a Normal MALE HUMAN  AT 600BCL
            call exm_init_voltages(EXM_ONEOHR_NORMAL_HUMAN_MALE_600BCL, mat, icelltype)
            call localvoltages2globalvolatges(mat, icelltype)
            initialized = .TRUE.

        else if((moneclmate_exm(2,mat) == 600) .and. (kfl_hfmod_exm(mat) == 3_ip) .and.(kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_inaga_exm(mat)==1_ip)) then  !if it's a Normal FEMALE HUMAN at 600BCL
            call exm_init_voltages(EXM_ONEOHR_NORMAL_HUMAN_FEMALE_600BCL, mat, icelltype)
            call localvoltages2globalvolatges(mat, icelltype)
            initialized = .TRUE.

         elseif ((moneclmate_exm(2, mat) /= 1000) .and. (kfl_drugsmate_exm(mat) == 0_ip) .and. (kfl_hfmod_exm(mat) /= 4_ip)) then
            !.or. (kfl_drugsmate_exm(mat) == 1) Normal with different heart rate .and. (moneclmate_exm(2,mat) /= 857).and. (moneclmate_exm(2,mat) /= 600
            call exm_init_voltages(EXM_ONEOHR_PIG_NODRUGS, mat, icelltype)

            if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
               call exm_oceohr(mat, icelltype, success_status_tmp, ohara_stats)
               success_status = success_status + success_status_tmp
            end if

            call localvoltages2globalvolatges(mat, icelltype)
            initialized = .TRUE.

         end if

      else if (kfl_hfmodmate_exm(mat) == 1_ip) then  !if it's Modified and Initial conditions need to be calculated

         if ((kfl_hfmod_exm(mat) /= 4_ip)) then      !.and. (moneclmate_exm(2,mat) /= 857) .and. (moneclmate_exm(2,mat) /= 1000).or. (kfl_drugsmate_exm(mat) == 1)Normal with different heart rate

            call exm_init_voltages(EXM_ONEOHR_NORMAL_MODIFIED, mat, icelltype)

            if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
               call exm_oceohr(mat, icelltype, success_status_tmp, ohara_stats)
               success_status = success_status + success_status_tmp
            end if

            call localvoltages2globalvolatges(mat, icelltype)
            initialized = .TRUE.

         else if ((kfl_hfmod_exm(mat) == 4_ip) .and. (kfl_inaga_exm(mat) == 1_ip)) then  !if it's a Modified PIG @ 400 BCL with INa modified channel from Passini

            call exm_init_voltages(EXM_ONEOHR_MODPIG_400BCL_MODINA, mat, icelltype)

            if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
               call exm_oceohr(mat, icelltype, success_status_tmp, ohara_stats)
               success_status = success_status + success_status_tmp
            end if

            call localvoltages2globalvolatges(mat, icelltype)
            initialized = .TRUE.

         end if
      end if

      if (success_status > 0_ip) then
         success_status = EXM_ONEOHR_ERROR_NOTCONVERGED
      end if

      if (.NOT. initialized) then
         success_status = EXM_ONEOHR_ERROR_NOTINITIALIZED
      end if

   end subroutine exm_oneohr

   subroutine localvoltages2globalvolatges(mat, ituss_exm)
      implicit none

      integer(ip), intent(in) :: mat, ituss_exm

      if (ituss_exm == EXM_CELLTYPE_EPI) then  !epi
         vauin_exm(:, 3, mat) = vaulo_exm(:, 1)
         vcoin_exm(:, 3, mat) = vcolo_exm(:, 1)
         vminimate_exm(3, mat) = elmlo_exm(2)
         !vauin_exm(:,1,3)=vaulo_exm(:,1)
         !vcoin_exm(:,1,3) = vcolo_exm(:,1)
         !vminimate_exm(3) = elmlo_exm(1)
      else if (ituss_exm == EXM_CELLTYPE_ENDO) then  !endo
         vauin_exm(:, 1, mat) = vaulo_exm(:, 1)
         vcoin_exm(:, 1, mat) = vcolo_exm(:, 1)
         vminimate_exm(1, mat) = elmlo_exm(2)
         !vauin_exm(:,2,3)=vaulo_exm(:,1)
         !vcoin_exm(:,2,3) = vcolo_exm(:,1)
         !vminimate_exm(3) = elmlo_exm(1)
      else if (ituss_exm == EXM_CELLTYPE_MID) then  !mid
         vauin_exm(:, 2, mat) = vaulo_exm(:, 1)
         vcoin_exm(:, 2, mat) = vcolo_exm(:, 1)
         vminimate_exm(2, mat) = elmlo_exm(2)
         !vauin_exm(:,3,3)=vaulo_exm(:,1)
         !vcoin_exm(:,3,3) = vcolo_exm(:,1)
         !vminimate_exm(1,3) = elmlo_exm(1)
      end if
   end subroutine

   subroutine exm_init_voltages(condition_type, mat, celltype)
      !------------------------------------------------------------------------
      !> @addtogroup Exmedi
      !> @{
      !> @file    exm_onecel.f90
      !> @date    03/10/2012
      !> @author  Mariano Vazquez
      !> @brief   One cell simulations to obtain initial conditions
      !> @details Runs for changes of cycle length and drug administration
      !> @}
      !------------------------------------------------------------------------

      implicit none

      integer(ip), intent(in) :: condition_type, mat, celltype

      select case (condition_type)

      case (EXM_ONEOHR_PIG_NODRUGS)
         select case (celltype)
         case (EXM_CELLTYPE_EPI)
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp

            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_ENDO)
            !endocardium
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_MID)
            !MIDmyocardial
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case default
            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
         end select

      case (EXM_ONEOHR_NORMAL_HUMAN_FEMALE_600BCL)
         select case (celltype)
         case (EXM_CELLTYPE_ENDO) !1     
            elmlo_exm(1:2) =   -87.664902277317367_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) =    1.1455758338192559E-004_rp
            vcolo_exm(2, 1:3) =    8.8679915666602636_rp
            vcolo_exm(3, 1:3) =    142.60572798713852_rp
            vcolo_exm(4, 1:3) =    142.60568727724311_rp
            vcolo_exm(5, 1:3) =    8.8678427399199471_rp
            vcolo_exm(6, 1:3) =    1.1587317953261430E-004_rp
            vcolo_exm(7, 1:3) =    3.1369878324189133_rp
            vcolo_exm(8, 1:3) =    2.7944805944253157_rp
            vcolo_exm(9, 1:3) =    7.0702906508163200E-007_rp
            vcolo_exm(10, 1:3) =    8.8375189003228634E-007_rp
            vcolo_exm(11, 1:3) =    7.4513326488172979E-002_rp
            vcolo_exm(12, 1:3) =    0.0000000000000000_rp
            vaulo_exm(1, 1:2) =    7.5972512633682120E-003_rp
            vaulo_exm(2, 1:2) =   0.81351764109286584_rp
            vaulo_exm(3, 1:2) =   0.81353058798150524_rp
            vaulo_exm(4, 1:2) =   0.80819140775794773_rp
            vaulo_exm(5, 1:2) =   0.61670107627034498_rp
            vaulo_exm(6, 1:2) =   0.78906362841872835_rp
            vaulo_exm(7, 1:2) =    2.0070623504288491E-004_rp
            vaulo_exm(8, 1:2) =   0.41228785270410750_rp
            vaulo_exm(9, 1:2) =   0.19128492327140989_rp
            vaulo_exm(10, 1:2) =    1.0241152945904416E-003_rp
            vaulo_exm(11, 1:2) =   0.99952689235110237_rp
            vaulo_exm(12, 1:2) =   0.32466518417485002_rp
            vaulo_exm(13, 1:2) =    5.2182028268431757E-004_rp
            vaulo_exm(14, 1:2) =   0.99952692611944138_rp
            vaulo_exm(15, 1:2) =   0.36076385857210003_rp
            vaulo_exm(16, 1:2) =    2.5350254970419230E-009_rp
            vaulo_exm(17, 1:2) =   0.99999998999729311_rp
            vaulo_exm(18, 1:2) =   0.77728679107429577_rp
            vaulo_exm(19, 1:2) =   0.99999998999759832_rp
            vaulo_exm(20, 1:2) =   0.98697365138801163_rp
            vaulo_exm(21, 1:2) =   0.99421320450114170_rp
            vaulo_exm(22, 1:2) =    9.0113765701275959E-003_rp
            vaulo_exm(23, 1:2) =   0.99999997975008115_rp
            vaulo_exm(24, 1:2) =   0.99999998956811254_rp
            vaulo_exm(25, 1:2) =    7.3225421282804371E-004_rp
            vaulo_exm(26, 1:2) =   0.70374586064651934_rp
            vaulo_exm(27, 1:2) =   0.46041257143396341_rp
            vaulo_exm(28, 1:2) =    2.0270889409482307E-004_rp
            vaulo_exm(29, 1:2) =   0.99684836820524714_rp
            vaulo_exm(30, 1:2) =    0.0000000000000000_rp         
         case (EXM_CELLTYPE_MID)  !2 
            elmlo_exm(1:2) =   -87.435195284762415_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) =    1.1890657170525055E-004_rp
            vcolo_exm(2, 1:3) =    9.7706696760397378_rp
            vcolo_exm(3, 1:3) =    141.66483952218840_rp
            vcolo_exm(4, 1:3) =    141.66479112719344_rp
            vcolo_exm(5, 1:3) =    9.7704648272293344_rp
            vcolo_exm(6, 1:3) =    1.1898783498170927E-004_rp
            vcolo_exm(7, 1:3) =    2.8021201980086876_rp
            vcolo_exm(8, 1:3) =    2.4236072150570278_rp
            vcolo_exm(9, 1:3) =    1.4071451498516887E-006_rp
            vcolo_exm(10, 1:3) =    1.7586168867274698E-006_rp
            vcolo_exm(11, 1:3) =    7.8962722885925360E-002_rp
            vcolo_exm(12, 1:3) =    0.0000000000000000_rp
            vaulo_exm(1, 1:2) =    7.7747291376501189E-003_rp
            vaulo_exm(2, 1:2) =   0.80784540409976191_rp
            vaulo_exm(3, 1:2) =   0.80786105547820475_rp
            vaulo_exm(4, 1:2) =   0.80265157891836159_rp
            vaulo_exm(5, 1:2) =   0.60792453114967615_rp
            vaulo_exm(6, 1:2) =   0.78403552237781493_rp
            vaulo_exm(7, 1:2) =    2.0965657866126548E-004_rp
            vaulo_exm(8, 1:2) =   0.41110392493267284_rp
            vaulo_exm(9, 1:2) =   0.19178946886523357_rp
            vaulo_exm(10, 1:2) =    1.0400968637168772E-003_rp
            vaulo_exm(11, 1:2) =   0.99950747422536579_rp
            vaulo_exm(12, 1:2) =   0.32788959685532149_rp
            vaulo_exm(13, 1:2) =    5.2996757393124953E-004_rp
            vaulo_exm(14, 1:2) =   0.99950751019835571_rp
            vaulo_exm(15, 1:2) =   0.36509196261881782_rp
            vaulo_exm(16, 1:2) =    2.6765011517336192E-009_rp
            vaulo_exm(17, 1:2) =   0.99999998935554446_rp
            vaulo_exm(18, 1:2) =   0.78131589613574182_rp
            vaulo_exm(19, 1:2) =   0.99999998935589263_rp
            vaulo_exm(20, 1:2) =   0.98788614956431220_rp
            vaulo_exm(21, 1:2) =   0.99431744834580349_rp
            vaulo_exm(22, 1:2) =    9.9542465724058414E-003_rp
            vaulo_exm(23, 1:2) =   0.99999998133288193_rp
            vaulo_exm(24, 1:2) =   0.99999998895728170_rp
            vaulo_exm(25, 1:2) =    6.6213956799707088E-004_rp
            vaulo_exm(26, 1:2) =   0.69998596202985708_rp
            vaulo_exm(27, 1:2) =   0.45745549722925427_rp
            vaulo_exm(28, 1:2) =    2.0754361583136973E-004_rp
            vaulo_exm(29, 1:2) =   0.99690661867951635_rp
            vaulo_exm(30, 1:2) =    0.0000000000000000_rp               
         case (EXM_CELLTYPE_EPI)  !3 
            elmlo_exm(1:2) =   -87.397959426295145_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) =    8.5433179446188953E-005_rp
            vcolo_exm(2, 1:3) =    9.7033871421276938_rp
            vcolo_exm(3, 1:3) =    140.63418994101886_rp
            vcolo_exm(4, 1:3) =    140.63414316826149_rp
            vcolo_exm(5, 1:3) =    9.7032770674922091_rp
            vcolo_exm(6, 1:3) =    8.3969893961209863E-005_rp
            vcolo_exm(7, 1:3) =    8.7338011361931187_rp
            vcolo_exm(8, 1:3) =    8.5090253044155926_rp
            vcolo_exm(9, 1:3) =    9.4348216776708040E-007_rp
            vcolo_exm(10, 1:3) =    1.1792687199169907E-006_rp
            vcolo_exm(11, 1:3) =   0.23304109107375579_rp
            vcolo_exm(12, 1:3) =    0.0000000000000000_rp
            vaulo_exm(1, 1:2) =    7.8038762211742538E-003_rp
            vaulo_exm(2, 1:2) =   0.80697892118129511_rp
            vaulo_exm(3, 1:2) =   0.80698048877554496_rp
            vaulo_exm(4, 1:2) =   0.80427282308871484_rp
            vaulo_exm(5, 1:2) =   0.60674144785235762_rp
            vaulo_exm(6, 1:2) =   0.79118681882009789_rp
            vaulo_exm(7, 1:2) =    2.1114419773127173E-004_rp
            vaulo_exm(8, 1:2) =   0.42822745461488432_rp
            vaulo_exm(9, 1:2) =   0.20738103535881042_rp
            vaulo_exm(10, 1:2) =    1.0426904028929566E-003_rp
            vaulo_exm(11, 1:2) =   0.99950453328143118_rp
            vaulo_exm(12, 1:2) =   0.99424423356861691_rp
            vaulo_exm(13, 1:2) =    5.3128975397284521E-004_rp
            vaulo_exm(14, 1:2) =   0.99950453350151214_rp
            vaulo_exm(15, 1:2) =   0.99701906127760564_rp
            vaulo_exm(16, 1:2) =    2.7000139951840619E-009_rp
            vaulo_exm(17, 1:2) =   0.99999998925616540_rp
            vaulo_exm(18, 1:2) =   0.81022448363846356_rp
            vaulo_exm(19, 1:2) =   0.99999998925619915_rp
            vaulo_exm(20, 1:2) =   0.99151690213156662_rp
            vaulo_exm(21, 1:2) =   0.99630964501327157_rp
            vaulo_exm(22, 1:2) =    2.6414450407093136E-003_rp
            vaulo_exm(23, 1:2) =   0.99999998823676040_rp
            vaulo_exm(24, 1:2) =   0.99999998919167776_rp
            vaulo_exm(25, 1:2) =    3.0381511578187663E-004_rp
            vaulo_exm(26, 1:2) =   0.67685433163727815_rp
            vaulo_exm(27, 1:2) =   0.41836968833467597_rp
            vaulo_exm(28, 1:2) =    2.0676096444089193E-004_rp
            vaulo_exm(29, 1:2) =   0.99691451138133580_rp
            vaulo_exm(30, 1:2) =    0.0000000000000000_rp                      
         case default
            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
         end select

      case (EXM_ONEOHR_NORMAL_HUMAN_MALE_600BCL)
         select case (celltype)
         case (EXM_CELLTYPE_ENDO) !1 
            elmlo_exm(1:2) =   -87.840677417569623_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) =    1.0642364619202648E-004_rp
            vcolo_exm(2, 1:3) =    7.9772077683146829_rp
            vcolo_exm(3, 1:3) =    143.86444519547052_rp
            vcolo_exm(4, 1:3) =    143.86440869454694_rp
            vcolo_exm(5, 1:3) =    7.9770757936138228_rp
            vcolo_exm(6, 1:3) =    1.0772077419226267E-004_rp
            vcolo_exm(7, 1:3) =    1.9994980946461740_rp
            vcolo_exm(8, 1:3) =    1.7061582377871256_rp
            vcolo_exm(9, 1:3) =    3.1469432182109148E-007_rp
            vcolo_exm(10, 1:3) =    3.9256396771398118E-007_rp
            vcolo_exm(11, 1:3) =    4.0081050009406739E-002_rp
            vcolo_exm(12, 1:3) =    0.0000000000000000_rp
            vaulo_exm(1, 1:2) =    7.4641620415121365E-003_rp
            vaulo_exm(2, 1:2) =   0.81778093747987457_rp
            vaulo_exm(3, 1:2) =   0.81778996374016499_rp
            vaulo_exm(4, 1:2) =   0.81525534740532923_rp
            vaulo_exm(5, 1:2) =   0.62341439722652625_rp
            vaulo_exm(6, 1:2) =   0.80380141909773406_rp
            vaulo_exm(7, 1:2) =    1.9411613137063226E-004_rp
            vaulo_exm(8, 1:2) =   0.43868600477080355_rp
            vaulo_exm(9, 1:2) =   0.21407703205078979_rp
            vaulo_exm(10, 1:2) =    1.0120484465101275E-003_rp
            vaulo_exm(11, 1:2) =   0.99954127668039483_rp
            vaulo_exm(12, 1:2) =   0.35456205922719886_rp
            vaulo_exm(13, 1:2) =    5.1566877333477810E-004_rp
            vaulo_exm(14, 1:2) =   0.99954130276542819_rp
            vaulo_exm(15, 1:2) =   0.39798819179101697_rp
            vaulo_exm(16, 1:2) =    2.4318143636470063E-009_rp
            vaulo_exm(17, 1:2) =   0.99999999046347776_rp
            vaulo_exm(18, 1:2) =   0.81829227729265408_rp
            vaulo_exm(19, 1:2) =   0.99999999046370502_rp
            vaulo_exm(20, 1:2) =   0.99159904776531749_rp
            vaulo_exm(21, 1:2) =   0.99607918538745932_rp
            vaulo_exm(22, 1:2) =    6.8365035674538763E-003_rp
            vaulo_exm(23, 1:2) =   0.99999998908994181_rp
            vaulo_exm(24, 1:2) =   0.99999999037131027_rp
            vaulo_exm(25, 1:2) =    2.9908608383753262E-004_rp
            vaulo_exm(26, 1:2) =   0.67433044628262184_rp
            vaulo_exm(27, 1:2) =   0.41460644563979088_rp
            vaulo_exm(28, 1:2) =    1.9695576450942084E-004_rp
            vaulo_exm(29, 1:2) =   0.99680275895078962_rp
            vaulo_exm(30, 1:2) =    0.0000000000000000_rp            
         case (EXM_CELLTYPE_MID)  !2    
            elmlo_exm(1:2) =   -87.435195284762415_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) =    1.1890657170525055E-004_rp
            vcolo_exm(2, 1:3) =    9.7706696760397378_rp
            vcolo_exm(3, 1:3) =    141.66483952218840_rp
            vcolo_exm(4, 1:3) =    141.66479112719344_rp
            vcolo_exm(5, 1:3) =    9.7704648272293344_rp
            vcolo_exm(6, 1:3) =    1.1898783498170927E-004_rp
            vcolo_exm(7, 1:3) =    2.8021201980086876_rp
            vcolo_exm(8, 1:3) =    2.4236072150570278_rp
            vcolo_exm(9, 1:3) =    1.4071451498516887E-006_rp
            vcolo_exm(10, 1:3) =    1.7586168867274698E-006_rp
            vcolo_exm(11, 1:3) =    7.8962722885925360E-002_rp
            vcolo_exm(12, 1:3) =    0.0000000000000000_rp
            vaulo_exm(1, 1:2) =    7.7747291376501189E-003_rp
            vaulo_exm(2, 1:2) =   0.80784540409976191_rp
            vaulo_exm(3, 1:2) =   0.80786105547820475_rp
            vaulo_exm(4, 1:2) =   0.80265157891836159_rp
            vaulo_exm(5, 1:2) =   0.60792453114967615_rp
            vaulo_exm(6, 1:2) =   0.78403552237781493_rp
            vaulo_exm(7, 1:2) =    2.0965657866126548E-004_rp
            vaulo_exm(8, 1:2) =   0.41110392493267284_rp
            vaulo_exm(9, 1:2) =   0.19178946886523357_rp
            vaulo_exm(10, 1:2) =    1.0400968637168772E-003_rp
            vaulo_exm(11, 1:2) =   0.99950747422536579_rp
            vaulo_exm(12, 1:2) =   0.32788959685532149_rp
            vaulo_exm(13, 1:2) =    5.2996757393124953E-004_rp
            vaulo_exm(14, 1:2) =   0.99950751019835571_rp
            vaulo_exm(15, 1:2) =   0.36509196261881782_rp
            vaulo_exm(16, 1:2) =    2.6765011517336192E-009_rp
            vaulo_exm(17, 1:2) =   0.99999998935554446_rp
            vaulo_exm(18, 1:2) =   0.78131589613574182_rp
            vaulo_exm(19, 1:2) =   0.99999998935589263_rp
            vaulo_exm(20, 1:2) =   0.98788614956431220_rp
            vaulo_exm(21, 1:2) =   0.99431744834580349_rp
            vaulo_exm(22, 1:2) =    9.9542465724058414E-003_rp
            vaulo_exm(23, 1:2) =   0.99999998133288193_rp
            vaulo_exm(24, 1:2) =   0.99999998895728170_rp
            vaulo_exm(25, 1:2) =    6.6213956799707088E-004_rp
            vaulo_exm(26, 1:2) =   0.69998596202985708_rp
            vaulo_exm(27, 1:2) =   0.45745549722925427_rp
            vaulo_exm(28, 1:2) =    2.0754361583136973E-004_rp
            vaulo_exm(29, 1:2) =   0.99690661867951635_rp
            vaulo_exm(30, 1:2) =    0.0000000000000000_rp               
         case (EXM_CELLTYPE_EPI)  !3 
            elmlo_exm(1:2) =   -87.708936866906583_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) =    8.0028671995560340E-005_rp
            vcolo_exm(2, 1:3) =    8.6136328499931931_rp
            vcolo_exm(3, 1:3) =    142.87867795811988_rp
            vcolo_exm(4, 1:3) =    142.87863578222584_rp
            vcolo_exm(5, 1:3) =    8.6135333361932869_rp
            vcolo_exm(6, 1:3) =    8.0347466740107238E-005_rp
            vcolo_exm(7, 1:3) =    3.3957368575898395_rp
            vcolo_exm(8, 1:3) =    3.1113557277334025_rp
            vcolo_exm(9, 1:3) =    4.7491800922715200E-007_rp
            vcolo_exm(10, 1:3) =    5.9362015321547805E-007_rp
            vcolo_exm(11, 1:3) =    3.9107271067814031E-002_rp
            vcolo_exm(12, 1:3) =    0.0000000000000000_rp
            vaulo_exm(1, 1:2) =    7.5636859241911452E-003_rp
            vaulo_exm(2, 1:2) =   0.81463650323285075_rp
            vaulo_exm(3, 1:2) =   0.81463987634369694_rp
            vaulo_exm(4, 1:2) =   0.81361739286268686_rp
            vaulo_exm(5, 1:2) =   0.61856297179956743_rp
            vaulo_exm(6, 1:2) =   0.80730247582796777_rp
            vaulo_exm(7, 1:2) =    1.9903435708260969E-004_rp
            vaulo_exm(8, 1:2) =   0.45892752188085889_rp
            vaulo_exm(9, 1:2) =   0.23701050019239650_rp
            vaulo_exm(10, 1:2) =    1.0210647092482490E-003_rp
            vaulo_exm(11, 1:2) =   0.99953077560081349_rp
            vaulo_exm(12, 1:2) =   0.99757663568230093_rp
            vaulo_exm(13, 1:2) =    5.2026513031970978E-004_rp
            vaulo_exm(14, 1:2) =   0.99953077622510722_rp
            vaulo_exm(15, 1:2) =   0.99877342602002361_rp
            vaulo_exm(16, 1:2) =    2.5086639233833685E-009_rp
            vaulo_exm(17, 1:2) =   0.99999999012175589_rp
            vaulo_exm(18, 1:2) =   0.86178725892506836_rp
            vaulo_exm(19, 1:2) =   0.99999999012183993_rp
            vaulo_exm(20, 1:2) =   0.99587532884666397_rp
            vaulo_exm(21, 1:2) =   0.99788907700093699_rp
            vaulo_exm(22, 1:2) =    2.2320899641454662E-003_rp
            vaulo_exm(23, 1:2) =   0.99999999004661755_rp
            vaulo_exm(24, 1:2) =   0.99999999011105045_rp
            vaulo_exm(25, 1:2) =    9.4095695945113766E-005_rp
            vaulo_exm(26, 1:2) =   0.62481556900203306_rp
            vaulo_exm(27, 1:2) =   0.34825026454480218_rp
            vaulo_exm(28, 1:2) =    1.9934081572285325E-004_rp
            vaulo_exm(29, 1:2) =   0.99683592429184664_rp
            vaulo_exm(30, 1:2) =    0.0000000000000000_rp                   
         case default
            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
         end select

      case (EXM_ONEOHR_NORMAL_PIG_400BCL)
         select case (celltype)
         case (EXM_CELLTYPE_EPI)
            vminimate_exm(3, mat) = -87.8061985708231_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            viclo_exm(1:26, 1) = 0.0_rp
            vcoin_exm(1, 3, mat) = 0.0001040361496276368_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2, 3, mat) = 7.43617875723091_rp !8.17212071948942       !nass=3);
            vcoin_exm(3, 3, mat) = 144.311817237457_rp !143.675184333376       !ki=4);
            vcoin_exm(4, 3, mat) = 144.311793923281_rp  !143.675147146842       !kss=5);
            vcoin_exm(5, 3, mat) = 7.43604263653445_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6, 3, mat) = 0.0001102064090489676_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7, 3, mat) = 2.63411603951980_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8, 3, mat) = 1.76343231449702_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1, 3, mat) = 0.007490096291770637_rp !0.00739719746920272       !m=10);
            vauin_exm(2, 3, mat) = 0.816838169651186_rp  !0.695621622011335       !hf=11);
            vauin_exm(3, 3, mat) = 0.816869319742196_rp  !0.695601842634086       !hs=12);
            vauin_exm(4, 3, mat) = 0.779416631991491_rp  !0.695486248719023       !j=13);
            vauin_exm(5, 3, mat) = 0.621614033150460_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6, 3, mat) = 0.719030208978245_rp  !0.695403157533235       !jp=15);
            vauin_exm(7, 3, mat) = 0.0001953921132645023_rp !0.000190839777466418       !mL=16);
            vauin_exm(8, 3, mat) = 0.350601487895652_rp  !0.493606704642336       !hL=17);
            vauin_exm(9, 3, mat) = 0.165483509753837_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10, 3, mat) = 0.001014440468937102_rp !0.00100594231451985      !a=19);
            vauin_exm(11, 3, mat) = 0.999538652553908_rp  !0.999548606668578      !iF=20);
            vauin_exm(12, 3, mat) = 0.940646485484524_rp  !0.999488774162635      !iS=21);
            vauin_exm(13, 3, mat) = 0.0005168881869164542_rp !0.000512555980943569      !ap=22);
            vauin_exm(14, 3, mat) = 0.999538658574815_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15, 3, mat) = 0.961396539930704_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16, 3, mat) = 0.000000002451977924935544_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17, 3, mat) = 0.999999990358952_rp  !0.999999990696210      !ff=26);
            vauin_exm(18, 3, mat) = 0.721766724570291_rp  !0.904906458666787      !fs=27);
            vauin_exm(19, 3, mat) = 0.999999990360122_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20, 3, mat) = 0.952750811437167_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21, 3, mat) = 0.962655802194759_rp  !0.999903346883777      !jca=30);
            vauin_exm(22, 3, mat) = 0.007755537033307213_rp !0.00215555277945401      !nca=31);
            vauin_exm(23, 3, mat) = 0.999989171225269_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24, 3, mat) = 0.999998595329089_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25, 3, mat) = 0.01177184211893142_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26, 3, mat) = 0.803851729756987_rp !0.487585264457487      !xrs=35);
            vauin_exm(27, 3, mat) = 0.525297568917089_rp !0.276203479404767      !xs1=36);
            vauin_exm(28, 3, mat) = 0.0004266215773384047_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29, 3, mat) = 0.996814746325558_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9, 3, mat) = 0.0000003867443373607881_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10, 3, mat) = 0.0000004808818699484394_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11, 3, mat) = 0.08284541646773119_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_ENDO)
            !endocardial
            vminimate_exm(1, mat) = -87.8077947685269_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            viclo_exm(1:26, 1) = 0.0_rp
            vcoin_exm(1, 1, mat) = 0.0001244664562667345_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2, 1, mat) = 7.38274584190492_rp !8.17212071948942       !nass=3);
            vcoin_exm(3, 1, mat) = 144.434867863316_rp !143.675184333376       !ki=4);
            vcoin_exm(4, 1, mat) = 144.434843290996_rp  !143.675147146842       !kss=5);
            vcoin_exm(5, 1, mat) = 7.38257910315418_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6, 1, mat) = 0.0001303227834520188_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7, 1, mat) = 2.22033410012224_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8, 1, mat) = 1.57032538537506_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1, 1, mat) = 0.007488900815703561_rp !0.00739719746920272       !m=10);
            vauin_exm(2, 1, mat) = 0.816818226826766_rp  !0.695621622011335       !hf=11);
            vauin_exm(3, 1, mat) = 0.816860080147452_rp  !0.695601842634086       !hs=12);
            vauin_exm(4, 1, mat) = 0.760420946701455_rp  !0.695486248719023       !j=13);
            vauin_exm(5, 1, mat) = 0.621427482142547_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6, 1, mat) = 0.686823279932452_rp  !0.695403157533235       !jp=15);
            vauin_exm(7, 1, mat) = 0.0001953332078784692_rp !0.000190839777466418       !mL=16);
            vauin_exm(8, 1, mat) = 0.321919485372494_rp  !0.493606704642336       !hL=17);
            vauin_exm(9, 1, mat) = 0.146893625330722_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10, 1, mat) = 0.001014350908419775_rp !0.00100594231451985      !a=19);
            vauin_exm(11, 1, mat) = 0.999537929745094_rp  !0.999548606668578      !iF=20);
            vauin_exm(12, 1, mat) = 0.199370793235339_rp  !0.999488774162635      !iS=21);
            vauin_exm(13, 1, mat) = 0.0005168425303908772_rp !0.000512555980943569      !ap=22);
            vauin_exm(14, 1, mat) = 0.999538048856435_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15, 1, mat) = 0.226043556528992_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16, 1, mat) = 0.000000002451190519566430_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17, 1, mat) = 0.999999990351502_rp  !0.999999990696210      !ff=26);
            vauin_exm(18, 1, mat) = 0.713425153962173_rp  !0.904906458666787      !fs=27);
            vauin_exm(19, 1, mat) = 0.999999990357108_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20, 1, mat) = 0.948912597173054_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21, 1, mat) = 0.957191372519062_rp  !0.999903346883777      !jca=30);
            vauin_exm(22, 1, mat) = 0.01457185320720413_rp !0.00215555277945401      !nca=31);
            vauin_exm(23, 1, mat) = 0.999971016040121_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24, 1, mat) = 0.999997537179088_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25, 1, mat) = 0.01854399770068590_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26, 1, mat) = 0.819260557541050_rp !0.487585264457487      !xrs=35);
            vauin_exm(27, 1, mat) = 0.535044673733069_rp !0.276203479404767      !xs1=36);
            vauin_exm(28, 1, mat) = 0.0006762955384707180_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29, 1, mat) = 0.996815803757701_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9, 1, mat) = 0.0000002397313133386565_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10, 1, mat) = 0.0000002975619496753238_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11, 1, mat) = 0.08074947458467882_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_MID)
            vminimate_exm(2, mat) = -87.5165638807769_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !midmyocardial
            viclo_exm(1:26, 1) = 0.0_rp
            vcoin_exm(1, 2, mat) = 0.0001302993381963230_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2, 2, mat) = 7.78650116481851_rp !8.17212071948942       !nass=3);
            vcoin_exm(3, 2, mat) = 143.967943434319_rp !143.675184333376       !ki=4);
            vcoin_exm(4, 2, mat) = 143.967928583526_rp  !143.675147146842       !kss=5);
            vcoin_exm(5, 2, mat) = 7.78627290424874_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6, 2, mat) = 0.0001367372627130561_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7, 2, mat) = 2.59618964681759_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8, 2, mat) = 1.50264239254053_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1, 2, mat) = 0.007711425764228608_rp !0.00739719746920272       !m=10);
            vauin_exm(2, 2, mat) = 0.809593790171627_rp  !0.695621622011335       !hf=11);
            vauin_exm(3, 2, mat) = 0.809668554567417_rp  !0.695601842634086       !hs=12);
            vauin_exm(4, 2, mat) = 0.676310283476841_rp  !0.695486248719023       !j=13);
            vauin_exm(5, 2, mat) = 0.609473347484487_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6, 2, mat) = 0.574825632378327_rp  !0.695403157533235       !jp=15);
            vauin_exm(7, 2, mat) = 0.0002064428815447646_rp !0.000190839777466418       !mL=16);
            vauin_exm(8, 2, mat) = 0.246316732862273_rp  !0.493606704642336       !hL=17);
            vauin_exm(9, 2, mat) = 0.103887440257557_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10, 2, mat) = 0.001034493006853912_rp !0.00100594231451985      !a=19);
            vauin_exm(11, 2, mat) = 0.999513339328340_rp  !0.999548606668578      !iF=20);
            vauin_exm(12, 2, mat) = 0.151520865587831_rp  !0.999488774162635      !iS=21);
            vauin_exm(13, 2, mat) = 0.0005271107526648554_rp !0.000512555980943569      !ap=22);
            vauin_exm(14, 2, mat) = 0.999513525196830_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15, 2, mat) = 0.167105642226589_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16, 2, mat) = 0.000000002626133449123401_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17, 2, mat) = 0.999999986012299_rp  !0.999999990696210      !ff=26);
            vauin_exm(18, 2, mat) = 0.611718123469599_rp  !0.904906458666787      !fs=27);
            vauin_exm(19, 2, mat) = 0.999999989545206_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20, 2, mat) = 0.904565662128277_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21, 2, mat) = 0.922416557341182_rp  !0.999903346883777      !jca=30);
            vauin_exm(22, 2, mat) = 0.01810224279643378_rp !0.00215555277945401      !nca=31);
            vauin_exm(23, 2, mat) = 0.999584611307330_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24, 2, mat) = 0.999972928171199_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25, 2, mat) = 0.05599708391131620_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26, 2, mat) = 0.860854175657345_rp !0.487585264457487      !xrs=35);
            vauin_exm(27, 2, mat) = 0.619347616204376_rp !0.276203479404767      !xs1=36);
            vauin_exm(28, 2, mat) = 0.003773361894746919_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29, 2, mat) = 0.996893681111765_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9, 2, mat) = 0.0000006261445113237194_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10, 2, mat) = 0.0000007713396255442226_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11, 2, mat) = 0.133755388160624_rp !0.0228529042639590      !CaMKt=41);
         case default
            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
         end select

      case (EXM_ONEOHR_NORMAL_70BPM)
         select case (celltype)
         case (EXM_CELLTYPE_ENDO)
            !endocardium
            vminimate_exm(1, mat) = -87.99_rp
            vcoin_exm(1, 1, mat) = 0.0000906297090590454_rp       !cai=6);
            vcoin_exm(2, 1, mat) = 7.50141023194658_rp       !nass=3);
            vcoin_exm(3, 1, mat) = 144.397336592521_rp       !ki=4);
            vcoin_exm(4, 1, mat) = 144.397306654274_rp       !kss=5);
            vcoin_exm(5, 1, mat) = 7.50131413867337_rp     ! nai=2
            vcoin_exm(6, 1, mat) = 0.0000899057612476716_rp       !cass=7);
            vcoin_exm(7, 1, mat) = 1.70096618787449_rp       !cansr=8);
            vcoin_exm(8, 1, mat) = 1.61469442244241_rp      !cajsr=9);
            vauin_exm(1, 1, mat) = 0.00735179196722269_rp       !m=10);
            vauin_exm(2, 1, mat) = 0.697747459518797_rp       !hf=11);
            vauin_exm(3, 1, mat) = 0.697717242697596_rp       !hs=12);
            vauin_exm(4, 1, mat) = 0.697546140963364_rp       !j=13);
            vauin_exm(5, 1, mat) = 0.454478178913152_rp       !hsp=14);
            vauin_exm(6, 1, mat) = 0.697428694436028_rp       !jp=15);
            vauin_exm(7, 1, mat) = 0.000188633292218496_rp       !mL=16);
            vauin_exm(8, 1, mat) = 0.490629978954784_rp       !hL=17);
            vauin_exm(9, 1, mat) = 0.256189072459875_rp       !hLp=18);
            vauin_exm(10, 1, mat) = 0.00100180193761673_rp      !a=19);
            vauin_exm(11, 1, mat) = 0.999553321574707_rp      !iF=20);
            vauin_exm(12, 1, mat) = 0.519939585076356_rp      !iS=21);
            vauin_exm(13, 1, mat) = 0.000510445304511332_rp      !ap=22);
            vauin_exm(14, 1, mat) = 0.999553334882239_rp      !iFp=23);
            vauin_exm(15, 1, mat) = 0.519939585076356_rp      !iSp=24);
            vauin_exm(16, 1, mat) = 0.00000000234656800943690_rp      !d=25);
            vauin_exm(17, 1, mat) = 0.999999990847537_rp      !ff=26);
            vauin_exm(18, 1, mat) = 0.890968078383048_rp      !fs=27);
            vauin_exm(19, 1, mat) = 0.999999990847864_rp      !fcaf=28);
            vauin_exm(20, 1, mat) = 0.999331407376516_rp      !fcas=29);
            vauin_exm(21, 1, mat) = 0.999861393178778_rp      !jca=30);
            vauin_exm(22, 1, mat) = 0.00342031822852804_rp      !nca=31);
            vauin_exm(23, 1, mat) = 0.999999990811414_rp      !ffp=32);
            vauin_exm(24, 1, mat) = 0.999999990842929_rp      !fcafp=33);
            vauin_exm(25, 1, mat) = 0.00000906339108941520_rp      !xrf=34);
            vauin_exm(26, 1, mat) = 0.519363745312448_rp      !xrs=35);
            vauin_exm(27, 1, mat) = 0.305861220188224_rp      !xs1=36);
            vauin_exm(28, 1, mat) = 0.000193120615511300_rp      !xs2=37);
            vauin_exm(29, 1, mat) = 0.996762767851278_rp      !xk1=38);
            vcoin_exm(9, 1, mat) = 0.000000270257677876248_rp      !Jrelnp=39);
            vcoin_exm(10, 1, mat) = 0.000000337551159684025_rp      !Jrelp=40);
            vcoin_exm(11, 1, mat) = 0.0172511826393971_rp      !CaMKt=41);
         case (EXM_CELLTYPE_EPI)
            vminimate_exm(3, mat) = -87.99_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !epicardium
            vcoin_exm(1, 3, mat) = 0.0000803149767106260_rp       !cai=6);
            vcoin_exm(2, 3, mat) = 8.17212071948942_rp       !nass=3);
            vcoin_exm(3, 3, mat) = 143.675184333376_rp       !ki=4);
            vcoin_exm(4, 3, mat) = 143.675147146842_rp       !kss=5);
            vcoin_exm(5, 3, mat) = 8.17202753912090_rp      ! nai=2
            vcoin_exm(6, 3, mat) = 0.0000796912949697674_rp       !cass=7);
            vcoin_exm(7, 3, mat) = 2.14811455007091_rp       !cansr=8);
            vcoin_exm(8, 3, mat) = 2.03335899236798_rp      !cajsr=9);
            vauin_exm(1, 3, mat) = 0.00739719746920272_rp       !m=10);
            vauin_exm(2, 3, mat) = 0.695621622011335_rp       !hf=11);
            vauin_exm(3, 3, mat) = 0.695601842634086_rp       !hs=12);
            vauin_exm(4, 3, mat) = 0.695486248719023_rp       !j=13);
            vauin_exm(5, 3, mat) = 0.452023628358454_rp       !hsp=14);
            vauin_exm(6, 3, mat) = 0.695403157533235_rp       !jp=15);
            vauin_exm(7, 3, mat) = 0.000190839777466418_rp       !mL=16);
            vauin_exm(8, 3, mat) = 0.493606704642336_rp       !hL=17);
            vauin_exm(9, 3, mat) = 0.264304293390731_rp       !hLp=18);
            vauin_exm(10, 3, mat) = 0.00100594231451985_rp      !a=19);
            vauin_exm(11, 3, mat) = 0.999548606668578_rp      !iF=20);
            vauin_exm(12, 3, mat) = 0.999488774162635_rp      !iS=21);
            vauin_exm(13, 3, mat) = 0.000512555980943569_rp      !ap=22);
            vauin_exm(14, 3, mat) = 0.999548607287668_rp      !iFp=23);
            vauin_exm(15, 3, mat) = 0.999488774162635_rp      !iSp=24);
            vauin_exm(16, 3, mat) = 0.00000000238076098345898_rp      !d=25);
            vauin_exm(17, 3, mat) = 0.999999990696210_rp      !ff=26);
            vauin_exm(18, 3, mat) = 0.904906458666787_rp      !fs=27);
            vauin_exm(19, 3, mat) = 0.999999990696060_rp      !fcaf=28);
            vauin_exm(20, 3, mat) = 0.999581201974281_rp      !fcas=29);
            vauin_exm(21, 3, mat) = 0.999903346883777_rp      !jca=30);
            vauin_exm(22, 3, mat) = 0.00215555277945401_rp      !nca=31);
            vauin_exm(23, 3, mat) = 0.999999990680285_rp      !ffp=32);
            vauin_exm(24, 3, mat) = 0.999999990692529_rp      !fcafp=33);
            vauin_exm(25, 3, mat) = 0.00000864222375034682_rp     !xrf=34);
            vauin_exm(26, 3, mat) = 0.487585264457487_rp      !xrs=35);
            vauin_exm(27, 3, mat) = 0.276203479404767_rp      !xs1=36);
            vauin_exm(28, 3, mat) = 0.000194412216700766_rp     !xs2=37);
            vauin_exm(29, 3, mat) = 0.996778581263402_rp      !xk1=38);
            vcoin_exm(9, 3, mat) = 0.000000474946280300893_rp       !Jrelnp=39);
            vcoin_exm(10, 3, mat) = 0.000000593539009244893_rp      !Jrelp=40);
            vcoin_exm(11, 3, mat) = 0.0228529042639590_rp      !CaMKt=41);
         case (EXM_CELLTYPE_MID)
            !midmyocardium
            vminimate_exm(2, mat) = -87.99_rp
            vcoin_exm(1, 2, mat) = 0.000109610482790345_rp       !cai=6);
            vcoin_exm(2, 2, mat) = 9.31335418275460_rp      !nass=3);
            vcoin_exm(3, 2, mat) = 142.448454117416_rp      !ki=4);
            vcoin_exm(4, 2, mat) = 142.448413204153_rp      !kss=5);
            vcoin_exm(5, 2, mat) = 9.31318571775255_rp     ! nai=2
            vcoin_exm(6, 2, mat) = 0.000106921402857184_rp       !cass=7);
            vcoin_exm(7, 2, mat) = 2.55909692264644_rp       !cansr=8);
            vcoin_exm(8, 2, mat) = 2.46618898848524_rp       !cajsr=9);
            vauin_exm(1, 2, mat) = 0.00762063517901599_rp       !m=10);
            vauin_exm(2, 2, mat) = 0.685227812017847_rp       !hf=11);
            vauin_exm(3, 2, mat) = 0.685186564524828_rp       !hs=12);
            vauin_exm(4, 2, mat) = 0.684946001497933_rp       !j=13);
            vauin_exm(5, 2, mat) = 0.439930036176900_rp       !hspexm_oneohr=15);
            vauin_exm(7, 2, mat) = 0.000201874904059089_rp       !mL=16);
            vauin_exm(8, 2, mat) = 0.469593913079625_rp       !hL=17);
            vauin_exm(9, 2, mat) = 0.232197586816947_rp       !hLp=18);
            vauin_exm(10, 2, mat) = 0.00102621879122940_rp      !a=19);
            vauin_exm(11, 2, mat) = 0.999524469433189_rp      !iF=20);
            vauin_exm(12, 2, mat) = 0.484655616738848_rp      !iS=21);
            vauin_exm(13, 2, mat) = 0.000522892623132084_rp      !ap=22);
            vauin_exm(14, 2, mat) = 0.999524488104748_rp      !iFp=23);
            vauin_exm(15, 2, mat) = 0.484655616738848_rp      !iSp=24);
            vauin_exm(16, 2, mat) = 0.00000000255334758464344_rp      !d=25);
            vauin_exm(17, 2, mat) = 0.999999989917362_rp      !ff=26);
            vauin_exm(18, 2, mat) = 0.850786415635368_rp      !fs=27);
            vauin_exm(19, 2, mat) = 0.999999989918100_rp      !fcaf=28);
            vauin_exm(20, 2, mat) = 0.998660902959982_rp      !fcas=29);
            vauin_exm(21, 2, mat) = 0.999737086062351_rp      !jca=30);
            vauin_exm(22, 2, mat) = 0.00660100554600241_rp      !nca=31);
            vauin_exm(23, 2, mat) = 0.999999989887965_rp      !ffp=32);
            vauin_exm(24, 2, mat) = 0.999999989911158_rp      !fcafp=33);
            vauin_exm(25, 2, mat) = 0.0000129078830800004_rp      !xrf=34);
            vauin_exm(26, 2, mat) = 0.556418397118542_rp      !xrs=35);
            vauin_exm(27, 2, mat) = 0.365366368645388_rp      !xs1=36);
            vauin_exm(28, 2, mat) = 0.000201042604770464_rp      !xs2=37);
            vauin_exm(29, 2, mat) = 0.996855464208211_rp      !xk1=38);
            vcoin_exm(9, 2, mat) = 0.00000189377036465905_rp      !Jrelnp=39);
            vcoin_exm(10, 2, mat) = 0.00000236711601348107_rp      !Jrelp=40);
            vcoin_exm(11, 2, mat) = 0.0406722957165916_rp       !CaMKt=41);
         case default
            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
         end select

      case (EXM_ONEOHR_NORMAL_1000BCL)
         select case (celltype)
         case (EXM_CELLTYPE_ENDO)
            !endocardium
            vminimate_exm(1, mat) = -87.0_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcoin_exm(1, 1, mat) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2, 1, mat) = 7.0_rp !8.17212071948942       !nass=3);
            vcoin_exm(3, 1, mat) = 145.0_rp !143.675184333376       !ki=4);
            vcoin_exm(4, 1, mat) = 145.0_rp  !143.675147146842       !kss=5);
            vcoin_exm(5, 1, mat) = 7.0_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6, 1, mat) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7, 1, mat) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8, 1, mat) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1, 1, mat) = 0.0_rp !0.00739719746920272       !m=10);
            vauin_exm(2, 1, mat) = 1.0_rp  !0.695621622011335       !hf=11);
            vauin_exm(3, 1, mat) = 1.0_rp  !0.695601842634086       !hs=12);
            vauin_exm(4, 1, mat) = 1.0_rp  !0.695486248719023       !j=13);
            vauin_exm(5, 1, mat) = 1.0_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6, 1, mat) = 1.0_rp  !0.695403157533235       !jp=15);
            vauin_exm(7, 1, mat) = 0.0_rp !0.000190839777466418       !mL=16);
            vauin_exm(8, 1, mat) = 1.0_rp  !0.493606704642336       !hL=17);
            vauin_exm(9, 1, mat) = 1.0_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10, 1, mat) = 0.0_rp !0.00100594231451985      !a=19);
            vauin_exm(11, 1, mat) = 1.0_rp  !0.999548606668578      !iF=20);
            vauin_exm(12, 1, mat) = 1.0_rp  !0.999488774162635      !iS=21);
            vauin_exm(13, 1, mat) = 0.0_rp !0.000512555980943569      !ap=22);
            vauin_exm(14, 1, mat) = 1.0_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15, 1, mat) = 1.0_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16, 1, mat) = 0.0_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17, 1, mat) = 1.0_rp  !0.999999990696210      !ff=26);
            vauin_exm(18, 1, mat) = 1.0_rp  !0.904906458666787      !fs=27);
            vauin_exm(19, 1, mat) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20, 1, mat) = 1.0_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21, 1, mat) = 1.0_rp  !0.999903346883777      !jca=30);
            vauin_exm(22, 1, mat) = 0.0_rp !0.00215555277945401      !nca=31);
            vauin_exm(23, 1, mat) = 1.0_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24, 1, mat) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25, 1, mat) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26, 1, mat) = 0.0_rp !0.487585264457487      !xrs=35);
            vauin_exm(27, 1, mat) = 0.0_rp !0.276203479404767      !xs1=36);
            vauin_exm(28, 1, mat) = 0.0_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29, 1, mat) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9, 1, mat) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10, 1, mat) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11, 1, mat) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_EPI)
            vminimate_exm(3, mat) = -87.0_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcoin_exm(1, 3, mat) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2, 3, mat) = 7.0_rp !8.17212071948942       !nass=3);
            vcoin_exm(3, 3, mat) = 145.0_rp !143.675184333376       !ki=4);
            vcoin_exm(4, 3, mat) = 145.0_rp  !143.675147146842       !kss=5);
            vcoin_exm(5, 3, mat) = 7.0_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6, 3, mat) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7, 3, mat) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8, 3, mat) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1, 3, mat) = 0.0_rp !0.00739719746920272       !m=10);
            vauin_exm(2, 3, mat) = 1.0_rp  !0.695621622011335       !hf=11);
            vauin_exm(3, 3, mat) = 1.0_rp  !0.695601842634086       !hs=12);
            vauin_exm(4, 3, mat) = 1.0_rp  !0.695486248719023       !j=13);
            vauin_exm(5, 3, mat) = 1.0_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6, 3, mat) = 1.0_rp  !0.695403157533235       !jp=15);
            vauin_exm(7, 3, mat) = 0.0_rp !0.000190839777466418       !mL=16);
            vauin_exm(8, 3, mat) = 1.0_rp  !0.493606704642336       !hL=17);
            vauin_exm(9, 3, mat) = 1.0_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10, 3, mat) = 0.0_rp !0.00100594231451985      !a=19);
            vauin_exm(11, 3, mat) = 1.0_rp  !0.999548606668578      !iF=20);
            vauin_exm(12, 3, mat) = 1.0_rp  !0.999488774162635      !iS=21);
            vauin_exm(13, 3, mat) = 0.0_rp !0.000512555980943569      !ap=22);
            vauin_exm(14, 3, mat) = 1.0_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15, 3, mat) = 1.0_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16, 3, mat) = 0.0_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17, 3, mat) = 1.0_rp  !0.999999990696210      !ff=26);
            vauin_exm(18, 3, mat) = 1.0_rp  !0.904906458666787      !fs=27);
            vauin_exm(19, 3, mat) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20, 3, mat) = 1.0_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21, 3, mat) = 1.0_rp  !0.999903346883777      !jca=30);
            vauin_exm(22, 3, mat) = 0.0_rp !0.00215555277945401      !nca=31);
            vauin_exm(23, 3, mat) = 1.0_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24, 3, mat) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25, 3, mat) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26, 3, mat) = 0.0_rp !0.487585264457487      !xrs=35);
            vauin_exm(27, 3, mat) = 0.0_rp !0.276203479404767      !xs1=36);
            vauin_exm(28, 3, mat) = 0.0_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29, 3, mat) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9, 3, mat) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10, 3, mat) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11, 3, mat) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_MID)
            vminimate_exm(2, mat) = -87.0_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcoin_exm(1, 2, mat) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcoin_exm(2, 2, mat) = 7.0_rp !8.17212071948942       !nass=3);
            vcoin_exm(3, 2, mat) = 145.0_rp !143.675184333376       !ki=4);
            vcoin_exm(4, 2, mat) = 145.0_rp  !143.675147146842       !kss=5);
            vcoin_exm(5, 2, mat) = 7.0_rp  !8.17202753912090      ! nai=2
            vcoin_exm(6, 2, mat) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcoin_exm(7, 2, mat) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcoin_exm(8, 2, mat) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vauin_exm(1, 2, mat) = 0.0_rp !0.00739719746920272       !m=10);
            vauin_exm(2, 2, mat) = 1.0_rp  !0.695621622011335       !hf=11);
            vauin_exm(3, 2, mat) = 1.0_rp  !0.695601842634086       !hs=12);
            vauin_exm(4, 2, mat) = 1.0_rp  !0.695486248719023       !j=13);
            vauin_exm(5, 2, mat) = 1.0_rp  !0.452023628358454       !hsp=14);
            vauin_exm(6, 2, mat) = 1.0_rp  !0.695403157533235       !jp=15);
            vauin_exm(7, 2, mat) = 0.0_rp !0.000190839777466418       !mL=16);
            vauin_exm(8, 2, mat) = 1.0_rp  !0.493606704642336       !hL=17);
            vauin_exm(9, 2, mat) = 1.0_rp  !0.264304293390731       !hLp=18);
            vauin_exm(10, 2, mat) = 0.0_rp !0.00100594231451985      !a=19);
            vauin_exm(11, 2, mat) = 1.0_rp  !0.999548606668578      !iF=20);
            vauin_exm(12, 2, mat) = 1.0_rp  !0.999488774162635      !iS=21);
            vauin_exm(13, 2, mat) = 0.0_rp !0.000512555980943569      !ap=22);
            vauin_exm(14, 2, mat) = 1.0_rp  !0.999548607287668      !iFp=23);
            vauin_exm(15, 2, mat) = 1.0_rp  !0.999488774162635      !iSp=24);
            vauin_exm(16, 2, mat) = 0.0_rp !2.38076098345898e-09      !d=25);
            vauin_exm(17, 2, mat) = 1.0_rp  !0.999999990696210      !ff=26);
            vauin_exm(18, 2, mat) = 1.0_rp  !0.904906458666787      !fs=27);
            vauin_exm(19, 2, mat) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vauin_exm(20, 2, mat) = 1.0_rp  !0.999581201974281      !fcas=29);
            vauin_exm(21, 2, mat) = 1.0_rp  !0.999903346883777      !jca=30);
            vauin_exm(22, 2, mat) = 0.0_rp !0.00215555277945401      !nca=31);
            vauin_exm(23, 2, mat) = 1.0_rp  !0.999999990680285      !ffp=32);
            vauin_exm(24, 2, mat) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vauin_exm(25, 2, mat) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vauin_exm(26, 2, mat) = 0.0_rp !0.487585264457487      !xrs=35);
            vauin_exm(27, 2, mat) = 0.0_rp !0.276203479404767      !xs1=36);
            vauin_exm(28, 2, mat) = 0.0_rp !0.000194412216700766      !xs2=37);
            vauin_exm(29, 2, mat) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcoin_exm(9, 2, mat) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcoin_exm(10, 2, mat) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcoin_exm(11, 2, mat) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case default
            call runend("EXM_INIT_VOLTAGES: Undefined cell type.")
         end select

      case (EXM_ONEOHR_COUPLED_LAND)
         select case (celltype)
         case (EXM_CELLTYPE_ENDO)
            !endocardium
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_EPI)
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_MID)
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case default
            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
         end select

      case (EXM_ONEOHR_NORMAL_MODIFIED)
         select case (celltype)
         case (EXM_CELLTYPE_ENDO)
            !endocardium
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_EPI)
            !epicardial
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_MID)
            !MIDmyocardial
            elmlo_exm(1:2) = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !elmlo_exm(2) = -87.99_rp
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.0_rp !0.0228529042639590      !CaMKt=41);
         case default
            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
         end select
      case (EXM_ONEOHR_MODPIG_400BCL_MODINA)
         select case (celltype)
         case (EXM_CELLTYPE_ENDO)
            !endocardium
            elmlo_exm(1:2) = -87.8077947685269_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001244664562667345_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.38274584190492_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 144.434867863316_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 144.434843290996_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.38257910315418_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001303227834520188_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 2.22033410012224_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.57032538537506_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.007488900815703561_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 0.816818226826766_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 0.816860080147452_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 0.760420946701455_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 0.621427482142547_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 0.686823279932452_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0001953332078784692_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 0.321919485372494_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 0.146893625330722_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.001014350908419775_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 0.999537929745094_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 0.199370793235339_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0005168425303908772_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 0.999538048856435_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 0.226043556528992_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.000000002451190519566430_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 0.999999990351502_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 0.713425153962173_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 0.999999990357108_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 0.948912597173054_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 0.957191372519062_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.01457185320720413_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 0.999971016040121_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 0.999997537179088_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.01854399770068590_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.819260557541050_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.535044673733069_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0006762955384707180_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 0.996815803757701_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 2.397313133386565E-007_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 2.975619496753238E-007_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.08074947458467882_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_EPI)
            !epicardial
            elmlo_exm(1:2) = -87.8061985708231_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001040361496276368_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.43617875723091_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 144.311817237457_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 144.311793923281_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.43604263653445_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001102064090489676_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 2.63411603951980_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.76343231449702_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.007490096291770637_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 0.816838169651186_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 0.816869319742196_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 0.779416631991491_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 0.621614033150460_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 0.719030208978245_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0001953921132645023_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 0.350601487895652_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 0.165483509753837_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.001014440468937102_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 0.999538652553908_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 0.940646485484524_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0005168881869164542_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 0.999538658574815_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 0.961396539930704_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.000000002451977924935544_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 0.999999990358952_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 0.721766724570291_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 0.999999990360122_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 0.952750811437167_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 0.962655802194759_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.007755537033307213_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 0.999989171225269_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 0.999998595329089_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.01177184211893142_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.803851729756987_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.525297568917089_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0004266215773384047_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 0.996814746325558_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0000003867443373607881_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0000004808818699484394_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.08284541646773119_rp !0.0228529042639590      !CaMKt=41);
         case (EXM_CELLTYPE_MID)
            !midmyocardial
            elmlo_exm(1:2) = -87.5165638807769_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            viclo_exm(1:26, 1) = 0.0_rp
            vcolo_exm(1, 1:3) = 0.0001302993381963230_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:3) = 7.78650116481851_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:3) = 143.967943434319_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:3) = 143.967928583526_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:3) = 7.78627290424874_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:3) = 0.0001367372627130561_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:3) = 2.59618964681759_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:3) = 1.50264239254053_rp  !2.03335899236798      !cajsr=9);
            vaulo_exm(1, 1:2) = 0.007711425764228608_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 0.809593790171627_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 0.809668554567417_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 0.676310283476841_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 0.609473347484487_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 0.574825632378327_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0002064428815447646_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 0.246316732862273_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 0.103887440257557_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.001034493006853912_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 0.999513339328340_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 0.151520865587831_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0005271107526648554_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 0.999513525196830_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 0.167105642226589_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.000000002626133449123401_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 0.999999986012299_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 0.611718123469599_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 0.999999989545206_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 0.904565662128277_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 0.922416557341182_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.01810224279643378_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 0.999584611307330_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 0.999972928171199_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.05599708391131620_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.860854175657345_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.619347616204376_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.003773361894746919_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 0.996893681111765_rp  !0.996778581263402      !xk1=38);
            vcolo_exm(9, 1:3) = 0.0000006261445113237194_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:3) = 0.0000007713396255442226_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:3) = 0.133755388160624_rp !0.0228529042639590      !CaMKt=41);
         case default
            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
         end select
      case default
         call runend("EXM_INIT_VOLATGES: Undefined condition type.")
      end select
   end subroutine exm_init_voltages

end module mod_exm_oneohr
