!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exm_onetor.f90
!> @author  mixed
!> @date    2019-11-16
!> @brief   mod_exm_onetor
!> @details mod_exm_onetor
!> @}
!-----------------------------------------------------------------------
module mod_exm_onetor

   use def_master
   use def_exmedi
   use mod_messages, only: messages_live
   use mod_exm_ocetor, only: exm_ocetor, exm_ocetla

   implicit none

   integer(ip), parameter :: &
      EXM_ONETOR_ERROR_CONVERGED = 0_ip, &
      EXM_ONETOR_ERROR_NOTCONVERGED = 1_ip, &
      EXM_ONETOR_ERROR_NOTINITIALIZED = 2_ip

   public :: exm_onetor
   private :: exm_init_voltages, localvoltages2globalvoltages

contains

!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_onetor.f90
!> @date    12/04/2013
!> @author  Zhinuo Jenny Wang
!> @brief   Sets initial conditions for ToR-ORd 2019 model
!> @}
!------------------------------------------------------------------------

   subroutine exm_onetor(mat, icelltype, land_variables, success_status, torord_stats)

      implicit none
      integer(ip), intent(in)  :: mat, icelltype
      integer(ip), intent(out) :: success_status
      integer(ip)              :: success_status_tmp
      logical                  :: flag_land, initialized
      real(rp), dimension(1:7), intent(out) :: land_variables

      real(rp), intent(inout)  :: torord_stats(3_ip) ! saves ToR-ORd states: number of beats, tolerance

      success_status = 0_ip

      flag_land = .FALSE.
      initialized = .FALSE. !flag to check if any of the IF was executed

      if (kfl_eccty(mat) == 3_ip .or. kfl_eccty(mat) == 4_ip) flag_land = .TRUE.

      ! No stored initial states at the moment, steady-state needs to be calculated for every simulation.
      call exm_init_voltages(mat, icelltype)
      if (((coupling('SOLIDZ', 'EXMEDI') >= 1_ip) .or. (coupling('EXMEDI', 'SOLIDZ') >= 1_ip)) .and. flag_land) then

         if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
            call exm_ocetla(mat, icelltype, land_variables, success_status_tmp, torord_stats)
            success_status = success_status + success_status_tmp
         end if
      else
         if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
            call exm_ocetor(mat, icelltype, success_status_tmp, torord_stats)
            success_status = success_status + success_status_tmp
         end if
      end if

      call localvoltages2globalvoltages(mat, icelltype)
      initialized = .TRUE.

      if (success_status > 0_ip) then
         success_status = EXM_ONETOR_ERROR_NOTCONVERGED
      end if

      if (.NOT. initialized) then
         success_status = EXM_ONETOR_ERROR_NOTINITIALIZED
      end if
  end subroutine exm_onetor

  subroutine localvoltages2globalvoltages(mat, ituss_exm)
     implicit none

     integer(ip), intent(in) :: mat, ituss_exm

     if (ituss_exm == EXM_CELLTYPE_EPI) then ! epi
        vauin_exm(:, 3, mat) = vaulo_exm(:, 1)
        vcoin_exm(:, 3, mat) = vcolo_exm(:, 1)
        vminimate_exm(3, mat) = elmlo_exm(2)
     else if (ituss_exm == EXM_CELLTYPE_ENDO) then ! endo
        vauin_exm(:, 1, mat) = vaulo_exm(:, 1)
        vcoin_exm(:, 1, mat) = vcolo_exm(:, 1)
        vminimate_exm(1, mat) = elmlo_exm(2)
     else if (ituss_exm == EXM_CELLTYPE_MID) then ! mid
        vauin_exm(:, 2, mat) = vaulo_exm(:, 1)
        vcoin_exm(:, 2, mat) = vcolo_exm(:, 1)
        vminimate_exm(2, mat) = elmlo_exm(2)
     end if
  end subroutine

  subroutine exm_init_voltages(mat, celltype)
     implicit none

     integer(ip), intent(in) :: mat, celltype

     vminimate_exm(celltype,mat) = -88.76_rp
 
     if (celltype == EXM_CELLTYPE_EPI) then
        ! Intracellular concentrations
        vcoin_exm(1,3,mat) = 12.8363_rp   ! [Nai] = mM
        vcoin_exm(2,3,mat) = 12.8366_rp   ! [Nass] = mM
        vcoin_exm(3,3,mat) = 142.6951_rp  ! [Ki] = mM
        vcoin_exm(4,3,mat) = 142.6951_rp  ! [Kss] = mM
        vcoin_exm(5,3,mat) = 6.6309e-5_rp ! [Cai] = mM
        vcoin_exm(6,3,mat) = 5.7672e-5_rp ! [Cass] = mM
        vcoin_exm(7,3,mat) = 1.8119_rp    ! [Cansr] = mM
        vcoin_exm(8,3,mat) = 1.8102_rp    ! [Cajsr] = mM
        vcoin_exm(9,3,mat) = 24.0_rp      ! [Cli] = mM
        vcoin_exm(10,3,mat) = 0.0129_rp   ! [CaMKt] = mM
        vcoin_exm(11,3,mat) = 0.15_rp     ! [KmCaMK] = mM
        vcoin_exm(12,3,mat) = 0.0015_rp   ! [KmCaM] = mM
        vcoin_exm(13,3,mat) = 0.0_rp      ! [Jrel_np] = mM/ms
        vcoin_exm(14,3,mat) = 0.0_rp      ! [Jrel_p] = mM/ms
     else if (celltype == EXM_CELLTYPE_MID) then
        vcoin_exm(1,2,mat) = 15.0038_rp   ! [Nai] = mM
        vcoin_exm(2,2,mat) = 15.0043_rp   ! [Nass] = mM
        vcoin_exm(3,2,mat) = 143.0403_rp  ! [Ki] = mM
        vcoin_exm(4,2,mat) = 143.0402_rp  ! [Kss] = mM
        vcoin_exm(5,2,mat) = 8.166e-05_rp ! [Cai] = mM
        vcoin_exm(6,2,mat) = 6.5781e-05_rp ! [Cass] = mM
        vcoin_exm(7,2,mat) = 1.9557_rp    ! [Cansr] = mM
        vcoin_exm(8,2,mat) = 1.9593_rp    ! [Cajsr] = mM
        vcoin_exm(9,2,mat) = 24.0_rp      ! [Cli] = mM
        vcoin_exm(10,2,mat) = 0.0192_rp   ! [CaMKt] = mM
        vcoin_exm(11,2,mat) = 0.15_rp     ! [KmCaMK] = mM
        vcoin_exm(12,2,mat) = 0.0015_rp   ! [KmCaM] = mM
        vcoin_exm(13,2,mat) = 0.0_rp      ! [Jrel_np] = mM/ms
        vcoin_exm(14,2,mat) = 0.0_rp      ! [Jrel_p] = mM/ms
     else if (celltype == EXM_CELLTYPE_ENDO) then 
        vcoin_exm(1,1,mat) = 12.1025_rp   ! [Nai] = mM
        vcoin_exm(2,1,mat) = 12.1029_rp   ! [Nass] = mM
        vcoin_exm(3,1,mat) = 142.3002_rp  ! [Ki] = mM
        vcoin_exm(4,1,mat) = 142.3002_rp  ! [Kss] = mM
        vcoin_exm(5,1,mat) = 7.453481e-05_rp ! [Cai] = mM
        vcoin_exm(6,1,mat) = 7.0305e-5_rp ! [Cass] = mM
        vcoin_exm(7,1,mat) = 1.5211_rp    ! [Cansr] = mM
        vcoin_exm(8,1,mat) = 1.5214_rp    ! [Cajsr] = mM
        vcoin_exm(9,1,mat) = 24.0_rp      ! [Cli] = mM
        vcoin_exm(10,1,mat) = 0.0111_rp   ! [CaMKt] = mM
        vcoin_exm(11,1,mat) = 0.15_rp     ! [KmCaMK] = mM
        vcoin_exm(12,1,mat) = 0.0015_rp   ! [KmCaM] = mM
        vcoin_exm(13,1,mat) = 0.0_rp      ! [Jrel_np] = mM/ms
        vcoin_exm(14,1,mat) = 0.0_rp      ! [Jrel_p] = mM/ms
     end if

     ! Gating state variables
     if (celltype == EXM_CELLTYPE_EPI) then 
        vauin_exm(1,3,mat) = 0.00074303_rp   ! m
        vauin_exm(2,3,mat) = 0.8360_rp      ! h
        vauin_exm(3,3,mat) = 0.8359_rp      ! j
        vauin_exm(4,3,mat) = 0.6828_rp      ! hp
        vauin_exm(5,3,mat) = 0.8357_rp      ! jp
        vauin_exm(6,3,mat) = 0.00015166_rp      ! mL
        vauin_exm(7,3,mat) = 0.5401_rp      ! hL
        vauin_exm(8,3,mat) = 0.3034_rp      ! hLp
        vauin_exm(9,3,mat) = 0.00092716_rp      ! a
        vauin_exm(10,3,mat) = 0.9996_rp        ! i_F
        vauin_exm(11,3,mat) = 0.9996_rp        ! iS
        vauin_exm(12,3,mat) = 0.0004724_rp     ! ap
        vauin_exm(13,3,mat) = 0.9996_rp        ! iFp
        vauin_exm(14,3,mat) = 0.9996_rp        ! iSp
        vauin_exm(15,3,mat) = 0.0_rp                ! d
        vauin_exm(16,3,mat) = 1.0_rp                   ! ff
        vauin_exm(17,3,mat) = 0.9485_rp        ! fs
        vauin_exm(18,3,mat) = 1.0_rp                   ! fcaf
        vauin_exm(19,3,mat) = 0.9999_rp        ! fcas
        vauin_exm(20,3,mat) = 1.0_rp                   ! jca
        vauin_exm(21,3,mat) = 1.0_rp                   ! ffp
        vauin_exm(22,3,mat) = 1.0_rp                   ! fcafp
        vauin_exm(23,3,mat) = 0.00030853_rp    ! nca_ss
        vauin_exm(24,3,mat) = 0.00053006_rp    ! nca_i
        vauin_exm(25,3,mat) = 0.00067941_rp   ! C1
        vauin_exm(26,3,mat) = 0.00082869_rp    ! C2
        vauin_exm(27,3,mat) = 0.9982_rp                   ! C3
        vauin_exm(28,3,mat) = 0.00027561_rp     ! O_IKr
        vauin_exm(29,3,mat) = 9.5416e-06_rp    ! I_IKr
        vauin_exm(30,3,mat) = 0.2309_rp        ! xs1
        vauin_exm(31,3,mat) = 0.00016975_rp     ! xs2
     else if (celltype == EXM_CELLTYPE_MID) then
        vauin_exm(1,2,mat) = 0.00073818_rp   ! m
        vauin_exm(2,2,mat) = 0.8365_rp      ! h
        vauin_exm(3,2,mat) = 0.8363_rp      ! j
        vauin_exm(4,2,mat) = 0.6838_rp      ! hp
        vauin_exm(5,2,mat) = 0.8358_rp      ! jp
        vauin_exm(6,2,mat) = 0.00015079_rp      ! mL
        vauin_exm(7,2,mat) = 0.5327_rp      ! hL
        vauin_exm(8,2,mat) = 0.2834_rp      ! hLp
        vauin_exm(9,2,mat) = 0.00092527_rp      ! a
        vauin_exm(10,2,mat) = 0.9996_rp        ! i_F
        vauin_exm(11,2,mat) = 0.5671_rp        ! iS
        vauin_exm(12,2,mat) = 0.00047143_rp     ! ap
        vauin_exm(13,2,mat) = 0.9996_rp        ! iFp
        vauin_exm(14,2,mat) = 0.6261_rp        ! iSp
        vauin_exm(15,2,mat) = 0.0_rp                ! d
        vauin_exm(16,2,mat) = 1.0_rp                   ! ff
        vauin_exm(17,2,mat) = 0.92_rp        ! fs
        vauin_exm(18,2,mat) = 1.0_rp                   ! fcaf
        vauin_exm(19,2,mat) = 0.9998_rp        ! fcas
        vauin_exm(20,2,mat) = 1.0_rp                   ! jca
        vauin_exm(21,2,mat) = 1.0_rp                   ! ffp
        vauin_exm(22,2,mat) = 1.0_rp                   ! fcafp
        vauin_exm(23,2,mat) = 0.00051399_rp    ! nca_ss
        vauin_exm(24,2,mat) = 0.0012_rp    ! nca_i
        vauin_exm(25,2,mat) = 0.00069560_rp   ! C1
        vauin_exm(26,2,mat) = 0.00082672_rp    ! C2
        vauin_exm(27,2,mat) = 0.9979_rp                   ! C3
        vauin_exm(28,2,mat) = 0.00054206_rp     ! O_IKr
        vauin_exm(29,2,mat) = 1.8784e-05_rp    ! I_IKr
        vauin_exm(30,2,mat) = 0.2653_rp        ! xs1
        vauin_exm(31,2,mat) = 0.00016921_rp     ! xs2
     else if(celltype == EXM_CELLTYPE_ENDO) then 
        vauin_exm(1,1,mat) = 8.0572e-4_rp   ! m
        vauin_exm(2,1,mat) = 0.8286_rp      ! h
        vauin_exm(3,1,mat) = 0.8284_rp      ! j
        vauin_exm(4,1,mat) = 0.6707_rp      ! hp
        vauin_exm(5,1,mat) = 0.8281_rp      ! jp
        vauin_exm(6,1,mat) = 1.629e-4_rp      ! mL
        vauin_exm(7,1,mat) = 0.5255_rp      ! hL
        vauin_exm(8,1,mat) = 0.2872_rp      ! hLp
        vauin_exm(9,1,mat) = 9.5098e-4_rp      ! a
        vauin_exm(10,1,mat) = 0.9996_rp        ! i_F
        vauin_exm(11,1,mat) = 0.5936_rp        ! iS
        vauin_exm(12,1,mat) = 4.8454e-4_rp     ! ap
        vauin_exm(13,1,mat) = 0.9996_rp        ! iFp
        vauin_exm(14,1,mat) = 0.6538_rp        ! iSp
        vauin_exm(15,1,mat) = 8.1084e-9_rp                ! d
        vauin_exm(16,1,mat) = 1.0_rp                   ! ff
        vauin_exm(17,1,mat) = 0.939_rp        ! fs
        vauin_exm(18,1,mat) = 1.0_rp                   ! fcaf
        vauin_exm(19,1,mat) = 0.9999_rp        ! fcas
        vauin_exm(20,1,mat) = 1.0_rp                   ! jca
        vauin_exm(21,1,mat) = 1.0_rp                   ! ffp
        vauin_exm(22,1,mat) = 1.0_rp                   ! fcafp
        vauin_exm(23,1,mat) = 6.6462e-4_rp    ! nca_ss
        vauin_exm(24,1,mat) = 0.0012_rp    ! nca_i
        vauin_exm(25,1,mat) = 7.0344e-4_rp   ! C1
        vauin_exm(26,1,mat) = 8.5109e-4_rp   ! C2
        vauin_exm(27,1,mat) = 0.9981_rp                   ! C3
        vauin_exm(28,1,mat) = 3.7585e-4_rp     ! O_IKr
        vauin_exm(29,1,mat) = 1.3289e-5_rp    ! I_IKr
        vauin_exm(30,1,mat) = 0.248_rp        ! xs1
        vauin_exm(31,1,mat) = 1.7707e-4_rp     ! xs2
     end if

     ! Initialise single cell variables
     viclo_exm(:,1:2) = 0.0_rp
     vcolo_exm(:,1) = vcoin_exm(:, celltype, mat)
     vcolo_exm(:,2) = vcoin_exm(:, celltype, mat)
     elmlo_exm(1:2) = vminimate_exm(celltype,mat)
     vaulo_exm(:,1) = vauin_exm(:, celltype, mat)
     vaulo_exm(:,2) = vauin_exm(:, celltype, mat)
  end subroutine exm_init_voltages

end module mod_exm_onetor
