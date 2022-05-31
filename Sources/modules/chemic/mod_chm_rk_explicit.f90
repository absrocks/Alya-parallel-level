module mod_chm_rk_explicit

  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_gradie
  use mod_ker_proper
  use mod_solver,          only : solver_lumped_mass_system
  use mod_solver,          only : solver_explicit
  use mod_memory,          only : memory_alloca
  use mod_memory,          only : memory_deallo
  use mod_chm_finiteRate,  only : chm_IntegrateSource_finiteRate
  use mod_chm_entropy

  implicit none

  real(rp),   pointer :: bt(:)         ! Energy RHS
  real(rp),   pointer :: tt(:)         ! Unknown
  real(rp),   pointer :: tt0(:,:)      ! Unknown
  real(rp),   pointer :: rt(:,:,:)     ! Unknown residuals
  real(rp),   pointer :: aux_rt(:,:)   ! Unknown residuals
  real(rp),   pointer :: aux_tt(:,:)   ! Unknown residuals
  real(rp),   pointer :: Mass(:)       ! lumped mass
  real(rp),   pointer :: Tim(:,:,:)    ! projection dt / (rho)
  real(rp)            :: dt(3)

  real (rp) :: a(5)
  real (rp) :: b(5,5)
  real (rp) :: sigma(5)
  
  type(soltyp), pointer  :: solve1(:)
  type(soltyp), pointer  :: solve2(:)
  real(rp),     pointer  :: x(:,:)
  
  private

  public :: chm_rk_explicit_solution
  public :: chm_rk_explicit_initialization
  public :: chm_rk_explicit_memory

contains
  
  subroutine chm_rk_explicit_initialization()

    nullify(rt)
    nullify(Tim)
    nullify(tt0)

    call chm_entropy_initialization()

  end subroutine chm_rk_explicit_initialization

  subroutine chm_rk_explicit_memory()

    bt       => rhsid
    tt       => unkno
    Mass     => vmass

    if( INOTEMPTY ) then
       if (kfl_model_chm == 4) then
          call memory_alloca(mem_modul(1:2,modul),'RT' ,'chm_rk_explicit_allocate',rt,nvar_CMC_chm,npoin,5_ip)
          call memory_alloca(mem_modul(1:2,modul),'AUX_RT','chm_rk_explicit_allocate',aux_rt,nvar_CMC_chm,npoin)
          call memory_alloca(mem_modul(1:2,modul),'AUX_TT','chm_rk_explicit_allocate',aux_tt,nvar_CMC_chm,npoin)
          call memory_alloca(mem_modul(1:2,modul),'TIM','chm_rk_explicit_allocate',Tim,nvar_CMC_chm,npoin,4_ip)
          call memory_alloca(mem_modul(1:2,modul),'TT0','chm_rk_explicit_allocate',tt0,nvar_CMC_chm,npoin)
       else
          call memory_alloca(mem_modul(1:2,modul),'RT' ,'chm_rk_explicit_allocate',rt,nclas_chm,npoin,5_ip)
          call memory_alloca(mem_modul(1:2,modul),'AUX_RT','chm_rk_explicit_allocate',aux_rt,nclas_chm,npoin)
          call memory_alloca(mem_modul(1:2,modul),'AUX_TT','chm_rk_explicit_allocate',aux_tt,nclas_chm,npoin)
          call memory_alloca(mem_modul(1:2,modul),'TIM','chm_rk_explicit_allocate',Tim,nclas_chm,npoin,4_ip)
          call memory_alloca(mem_modul(1:2,modul),'TT0','chm_rk_explicit_allocate',tt0,nclas_chm,npoin)
       end if
    else
       call memory_alloca(mem_modul(1:2,modul),'RT' ,'chm_rk_explicit_allocate',rt,1_ip,1_ip,5_ip)
       call memory_alloca(mem_modul(1:2,modul),'AUX_RT','chm_rk_explicit_allocate',aux_rt,1_ip,1_ip)
       call memory_alloca(mem_modul(1:2,modul),'TIM','chm_rk_explicit_allocate',Tim,1_ip,1_ip,4_ip)
       call memory_alloca(mem_modul(1:2,modul),'TT0','chm_rk_explicit_allocate',tt0,1_ip,1_ip)
    end if

    solve1 => solve(1:)
    solve2 => solve(2:)
    x      => conce(:,:,1)

    dt    = 1.0_rp
    a     = 0.0_rp
    b     = 0.0_rp
    sigma = 0.0_rp

    if(kfl_tiacc_chm == 1) then !RK1
       a(4) = 1.0_rp
       sigma(4) = 1.0_rp
    else if(kfl_tiacc_chm == 2) then !Heun’s method or RK2 just for validation 
       a(4) = 1.0_rp

       b(4,3) = 1.0_rp

       sigma(3) = 0.5_rp
       sigma(4) = 0.5_rp
    else if(kfl_tiacc_chm == 3) then ! standard 3 order RK
       a(3) = 1.0_rp
       a(4) = 0.5_rp  

       b(3,2) = 1.0_rp
       b(4,2) = 1.0_rp/4.0_rp
       b(4,3) = 1.0_rp/4.0_rp

       sigma(2) = 1.0_rp/6.0_rp 
       sigma(3) = 1.0_rp/6.0_rp
       sigma(4) = 2.0_rp/3.0_rp
    else if(kfl_tiacc_chm == 4) then ! standard 4 order RK
       a(2) = 0.5_rp
       a(3) = 0.5_rp
       a(4) = 1.0_rp

       b(2,1) = 0.5_rp 
       b(3,2) = 0.5_rp
       b(4,3) = 1.0_rp

       sigma(1) = 1.0_rp/6.0_rp 
       sigma(2) = 1.0_rp/3.0_rp 
       sigma(3) = 1.0_rp/3.0_rp 
       sigma(4) = 1.0_rp/6.0_rp 
    end if

    call chm_entropy_memory()

  end subroutine chm_rk_explicit_memory


  subroutine chm_multi_step_fs_eval(iclai,iclaf,istep)
    integer(ip) , intent(in) :: iclai
    integer(ip) , intent(in) :: iclaf
    integer(ip) , intent(in) :: istep
    integer(ip)              :: ipoin,iclass,kpoin,iclaf_gas,iclai_liq

    dt_rho_chm = 0.0_rp

    if ( kfl_spray_chm /= 0_ip ) & 
       dt_chm = 0.0_rp

    call chm_matrix()

    call solver_lumped_mass_system(1_ip,dt_rho_chm)

    if ( kfl_spray_chm /= 0_ip ) & 
         call solver_lumped_mass_system(1_ip,dt_chm)   
    
    !
    ! Coupling with partis
    !
    call chm_partis()

    !
    ! Solve explicit system
    !
    if ( kfl_spray_chm /= 0_ip ) then

       iclaf_gas = iclaf     - 2
       iclai_liq = iclaf_gas + 1

       do ipoin = 1,npoin
          do iclass = iclai,iclaf_gas 
             Tim(iclass,ipoin,1) = dt_rho_chm(ipoin)/dt(1)
          end do
          do iclass = iclai_liq,iclaf
             Tim(iclass,ipoin,1) = dt_chm(ipoin)/dt(1)
          end do
       end do

    else
       do ipoin = 1,npoin
          do iclass = iclai,iclaf
             Tim(iclass,ipoin,1) = dt_rho_chm(ipoin)/dt(1)
          end do
       end do
    end if

    kpoin = 0
    do ipoin = 1,npoin
       do iclass = iclai,iclaf
          kpoin             = kpoin + 1
          aux_rt(iclass,ipoin) = bt(kpoin)
       end do
    end do

    if (kfl_model_chm == 4) then
       call rhsmod(nvar_CMC_chm,aux_rt)
    else
       call rhsmod(nclas_chm,aux_rt)
    end if

    if( INOTMASTER ) then 
       do ipoin = 1,npoin
          do iclass = iclai,iclaf
             rt(iclass,ipoin,istep) = aux_rt(iclass,ipoin)
          end do
       end do
    endif
    
    if(istep == 4_ip) then
       call chm_multi_step_fs_solution_sij(iclai,iclaf,4_ip,sigma,1.0_rp)
    else
       call chm_multi_step_fs_solution_sij(iclai,iclaf,istep,b(istep+1,:),a(istep+1))
    endif

  end subroutine chm_multi_step_fs_eval


 subroutine chm_multi_step_fs_solution_sij(iclai,iclaf,istep,weight1,weight2)
   integer(ip) , intent(in) :: iclai
   integer(ip) , intent(in) :: iclaf
   integer(ip) , intent(in) :: istep
   real(rp) ,    intent(in) :: weight1(5)
   real(rp) ,    intent(in) :: weight2
   integer(ip)       :: ipoin,jstep,iclass,kpoin
   real(rp) :: aux2  

   if( INOTMASTER ) then

      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf             
            aux2 = 0.0_rp
            kpoin     = kpoin + 1
            do jstep=1,istep
               aux2   = aux2 + rt(iclass,ipoin,jstep)*weight1(jstep)
            end do
            !tt(kpoin) = tt0(iclass,ipoin) + ( Tim(iclass,ipoin,1)*aux2*dt(1) ) / Mass(ipoin)
            aux_tt(iclass,ipoin) = Tim(iclass,ipoin,1)*aux2*dt(1)
         end do
      end do
      !
      ! Solve system
      !
      call solver_explicit(solve1,aux_tt,EXCHANGE=.false.,solve_consistent=solve2)
      
      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
            kpoin = kpoin+1
            tt(kpoin) = tt0(iclass,ipoin) + aux_tt(iclass,ipoin)
         end do
      end do
      !
      ! Dirichet bbcc 
      !
      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
            kpoin = kpoin + 1
            if( kfl_fixno_chm(iclass,ipoin) > 0 ) &
                 tt(kpoin) = bvess_chm(iclass,ipoin)
         end do
      end do

   end if

  end subroutine chm_multi_step_fs_solution_sij


 subroutine chm_rk_explicit_solution()

   use mod_chm_operations_CMC,     only : chm_calc_diff_condVar_mixfraction_CMC, &
                                          chm_integrate_chem_source_CMC, &
                                          chm_calc_temp_CMC, &
                                          chm_local2global_CMC, chm_global2local_CMC, &
                                          chm_cp_k_gauss_CMC


   integer(ip) :: ipoin,kpoin
   integer(ip), save :: iter = 0
   real(rp)          :: dt_split

   iter = iter + 1

   !
   ! Update inner iteration counter
   !
   itinn(modul) = itinn(modul) + 1
   ittot_chm    = ittot_chm + 1
   rtpts_chm     =  0.0_rp
   comin_chm     =  1.0e9_rp
   comax_chm     = -1.0e9_rp
   kfl_goit2_chm =  0_ip 

   !
   ! Allocate memory 
   !
   bt       => rhsid
   tt       => unkno
   Mass     => vmass

   !
   ! Assemble equations: 
   !
   dt(3) = dt(2)
   dt(2) = dt(1)
   dt(1) = 1.0_rp/dtinv

   if (kfl_model_chm == 4) then
      if(kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0) then
         do ipoin=1,npoin
            veloc_CFD_chm(:,ipoin) = advec(:,ipoin,1)
         end do
      end if
      call chm_calc_diff_condVar_mixfraction_CMC   ! For CMC compute second deriv. in mixt. frac.
   end if

   if (kfl_model_chm == 4) then
      iclai_chm = 1
      iclaf_chm = nvar_CMC_chm
   else
      iclai_chm = 1
      iclaf_chm = nclas_chm
   end if

   !
   ! Update boundary conditions
   !
   if (kfl_start_CMC_chm == 0)  call chm_updbcs(1)

   if( INOTMASTER) then
      !
      ! Integrate source terms
      !
      if(kfl_split_chm > 1_ip) then
         dt_split = dt(1) / dble(kfl_split_chm)
         if (kfl_model_chm == 3) then
             call chm_IntegrateSource_finiteRate(dt_split)
         else if (kfl_model_chm == 4) then
             call chm_integrate_chem_source_CMC(dt_split)
         end if
      end if
   end if


   mixt_fr: do imixf_rk = 2,nZ_CMC_chm-1  ! If we are using finite rate nZ_CMC_chm is fixed at 3 -> we enter in the loop just once

      if( INOTMASTER) then
         !
         ! Define unkno(ipoin) from conce(:,iclas,1) for all iclas or 
         ! from Yk_CMC_chm and enthalp_CMC_chm for CMC model
         !
         if (kfl_model_chm == 4) then
            call chm_global2local_CMC(imixf_rk)
            call chm_updbcs(2)
         end if
         call chm_updunk(ITASK_BEGINN)  

         if ( kfl_model_chm == 4 ) then
            ! Compute properties
            call chm_cp_k_gauss_CMC
         end if

         kpoin = 0
         do ipoin = 1,npoin
            do iclas_chm = iclai_chm, iclaf_chm
               kpoin                = kpoin + 1
               tt0(iclas_chm,ipoin) = tt(kpoin)
            end do
         end do

      end if

      !
      ! Entropy prediction:
      !
      if(kfl_entropy_chm == 1_ip) call chm_entropy_solution()
      !!!!!!!!!!!!!!! VER SI AQUÍ HAY QUE CAMBIAR ALGO PARA EL CMC


      if(     kfl_tiacc_chm == 1) then

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,4_ip)

      else if(     kfl_tiacc_chm == 2) then

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,3_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,4_ip)

      else if(kfl_tiacc_chm == 3) then  

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,2_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)
         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,3_ip)

         !
         ! Clipping + update conce(:,:,1)
         !

         call chm_endite(ITASK_INNITE)

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,4_ip)

      else ! order 4

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,1_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,2_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,3_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,4_ip)

      end if

      call chm_endite(ITASK_ENDINN)
      ! Note: if enthalpy is not transported it does not change with time -> enthalpy
      ! matrix already filled during the initialization

      if (kfl_model_chm == 4) then
         call chm_local2global_CMC(imixf_rk)
         call chm_calc_temp_CMC(imixf_rk)
      end if

   end do mixt_fr

 end subroutine chm_rk_explicit_solution

end module mod_chm_rk_explicit
 
