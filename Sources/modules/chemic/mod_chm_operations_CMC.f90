module mod_chm_operations_CMC

#ifdef CANTERA
  use cantera
#endif

  use def_kintyp,              only   : ip,rp
  use def_domain,              only   : ndime,npoin,coord,mgaus,mnode
  use def_chemic,              only   : nclas_chm
  use def_master,              only   : inotmaster
  use mod_memory,              only   : memory_alloca, memory_deallo
  use mod_memory,              only   : memory_deallo
  use def_master,              only   : mem_modul,modul
  use mod_chm_finiteRate,      only   : chm_elmprc_finiteRate

  implicit none

  private

  public :: chm_inert_eq_CMC
  public :: chm_compute_initial_fields_CMC
  public :: chm_initialization_domain_CMC
  public :: chm_bc_type_species_CMC
  public :: chm_save_unconditional_fields_bc_CMC
  public :: chm_get_cond_fields_bc_Dir_CMC
  public :: chm_data_CFD_to_CMC_domain_CMC
  public :: chm_data_CMC_to_CFD_domain_CMC
  public :: chm_local2global_CMC
  public :: chm_global2local_CMC
  public :: chm_element_operations_CMC
  public :: chm_AMC_generate_Z_S_vectors_CMC
  public :: chm_AMC_integrals_CMC
  public :: compute_Xnormalized_profile_CMC
  public :: chm_calc_diff_condVar_mixfraction_CMC
  public :: chm_calc_temp_CMC
  public :: chm_calc_densi_visco_gauss_CMC
  public :: chm_integrate_flow_var_points_CMC
  public :: chm_integrate_flow_var_gauss_CMC
  public :: chm_rho_visco_nodal_project_CMC
  public :: chm_integrate_chem_source_CMC
  public :: chm_updtcc_CMC
  public :: chm_cp_k_gauss_CMC
  public :: chm_limit_Yk_CMC
  public :: chm_heatRelease_integral_CMC
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! I N I T I A L   S U B R O U T I N E S !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_inert_eq_CMC
  
    !-----------------------------------------------------------------------
    !****f* Chemic/mod_chm_operations_CMC/chm_inert_eq_CMC
    ! NAME
    !    chm_inert_eq_CMC
    ! DESCRIPTION
    !    Compute the inert and equilbrium to initialize CMC fields for the
    !    first time.
    ! USES
    ! USED BY
    !    chm_reaphy
    !    chm_initialization_domain_CMC
    !*** 
    !-----------------------------------------------------------------------
  
    use def_parame
    use def_master,     only :  prthe, speci
    use def_kermod,     only :  gasco
    use def_chemic,     only :  nZ_CMC_chm, react_scal_bc_CMC_chm, rscal_inert_CMC_chm, &
                                rscal_equil_CMC_chm, Z_CMC_chm, Zs_CMC_chm, T_bc_CMC_chm, &
                                W_k, mechanism_path, nsize_mech_name, temp_inert_CMC_chm, &
                                temp_equil_CMC_chm
    use mod_physics,    only :  physics_H_2_TCp, physics_T_2_HCp

    implicit none
    integer(ip)              :: imixf, iclas, iter, ivalu
    real(rp)                 :: ratio_mf, h_mixt, Tini, Taux
    real(rp)                 :: coeff_cp_k_aux(nclas_chm,15) ! Memory still not assigned to coeff_cp_k
    real(rp)                 :: cploc(6,2), aux_cp_lt(8), aux_cp_ht(8)
    real(rp)                 :: aux_cp
    character(len=20)        :: names_rscal(nclas_chm+3), number_str
  
#ifdef CANTERA
    call cantera_initialization(1_ip,mechanism_path,nsize_mech_name) ! Create gas object to be used in next functions
    call cantera_coeff(mechanism_path,nsize_mech_name,coeff_cp_k_aux, W_k)

    ! Get enthalpies at boundaries from temperature and composition
    do iter = 1,2
       ! chm_calc_h_from_TY_CMC cannot be called since coeff_cp_k has not memory assigned yet
       call cantera_alya_cp(nclas_chm, coeff_cp_k_aux, react_scal_bc_CMC_chm(2_ip:nclas_chm+1_ip,iter), &
                 W_k, aux_cp_lt(1:8), aux_cp_ht(1:8))
       do ivalu = 1, 6
          cploc(ivalu,1) = aux_cp_lt(ivalu)
          cploc(ivalu,2) = aux_cp_ht(ivalu)
       end do
       call physics_T_2_HCp(T_bc_CMC_chm(iter), cploc, h_mixt, aux_cp)

       react_scal_bc_CMC_chm(1,iter) = h_mixt
    end do
#endif


    ! Get inert profiles
    do iclas = 1,nclas_chm+1_ip
       rscal_inert_CMC_chm(1,iclas)          = react_scal_bc_CMC_chm(iclas,1)
       rscal_inert_CMC_chm(nZ_CMC_chm,iclas) = react_scal_bc_CMC_chm(iclas,2)
    end do
  
    do imixf = 2,nZ_CMC_chm-1_ip
       ratio_mf = Z_CMC_chm(imixf) / Zs_CMC_chm
       do iclas = 1,nclas_chm+1_ip
          rscal_inert_CMC_chm(imixf,iclas) = react_scal_bc_CMC_chm(iclas,1) + &
               (react_scal_bc_CMC_chm(iclas,2) - react_scal_bc_CMC_chm(iclas,1)) * ratio_mf
       end do
    end do

    ! Get inert temperature from enthalpy and mass fractions
#ifdef CANTERA
    do imixf  = 1,nZ_CMC_chm
       ! Cantera is called to get an initial temperature to start iterating from
       ! Cantera and Alya values are very close but for consistency we take as the good
       ! one the one from Alya
       call cantera_get_temperature_from_hy(prthe(1), rscal_inert_CMC_chm(imixf,1), &
               rscal_inert_CMC_chm(imixf,2_ip:nclas_chm+1_ip), Tini)

       call cantera_alya_cp(nclas_chm, coeff_cp_k_aux, rscal_inert_CMC_chm(imixf,2_ip:nclas_chm+1_ip), &
              W_k, aux_cp_lt(1:8), aux_cp_ht(1:8))
       do ivalu = 1, 6
          cploc(ivalu,1) = aux_cp_lt(ivalu)
          cploc(ivalu,2) = aux_cp_ht(ivalu)
       end do
       Taux = Tini
       call physics_H_2_TCp(rscal_inert_CMC_chm(imixf,1), cploc, Taux, aux_cp)

       temp_inert_CMC_chm(imixf) = Taux 
    end do
#endif
 
    !
    ! Get equilibrium profiles
    !
    temp_equil_CMC_chm  = temp_inert_CMC_chm
    rscal_equil_CMC_chm = rscal_inert_CMC_chm  ! Enthalpy is the same in inert and equilibrium conditions. Species mass fractions
                                               ! are overwritten in the following
#ifdef CANTERA
    do imixf = 2,nZ_CMC_chm-1_ip
       ! Get chemical equilibrium
       call cantera_equilibrium_from_hy(prthe(1), rscal_equil_CMC_chm(imixf,1), &
               rscal_equil_CMC_chm(imixf,2_ip:nclas_chm+1_ip), Tini)

       ! Cantera and Alya values are very close but for consistency we take as the good
       ! one the one from Alya
       call cantera_alya_cp(nclas_chm, coeff_cp_k_aux, rscal_equil_CMC_chm(imixf,2_ip:nclas_chm+1_ip), &
              W_k, aux_cp_lt(1:8), aux_cp_ht(1:8))
       do ivalu = 1, 6
          cploc(ivalu,1) = aux_cp_lt(ivalu)
          cploc(ivalu,2) = aux_cp_ht(ivalu)
       end do
       Taux = Tini
       call physics_H_2_TCp(rscal_equil_CMC_chm(imixf,1), cploc, Taux, aux_cp)

       temp_equil_CMC_chm(imixf) = Taux
    end do
#endif

    !
    ! Write files for inert and equilibrium solutions
    !

    ! Get the names of the reactive scalars in a vector for the header
    names_rscal(1) = 'MIXT_FRAC [-]:1'
    names_rscal(2) = 'TEMPERATURE [K]:2'
    names_rscal(3) = 'ENTHALPY [J/kg]:3'
    do iclas = 1,nclas_chm
       write(number_str,fmt="(i4)") (iclas+3)
       names_rscal(iclas+3) = 'Y' // trim(speci(iclas)%name) // ' [-]:' // number_str
    end do

    ! Write inert mixture file
    open(unit=1, file='inert_solution.log', status="replace", action="write")
    write(unit=1,fmt="(300(a20))") names_rscal
    do imixf = 1,nZ_CMC_chm
       write(unit=1,fmt="(300(e12.6,8x))") Z_CMC_chm(imixf), temp_inert_CMC_chm(imixf), &
            rscal_inert_CMC_chm(imixf,1_ip:nclas_chm+1_ip)
    end do
    close(unit=1)
    print*, '--| ALYA     CHEMIC: INERT SOLUTION WRITTEN FOR CMC MODEL'  

    ! Write equilibrium mixture file
    open(unit=1, file='equilibrium_solution.log', status="replace", action="write")
    write(unit=1,fmt="(300(a20))") names_rscal
    do imixf = 1,nZ_CMC_chm
       write(unit=1,fmt="(300(e12.6,8x))") Z_CMC_chm(imixf), temp_equil_CMC_chm(imixf), &
            rscal_equil_CMC_chm(imixf,1_ip:nclas_chm+1_ip)
    end do
    close(unit=1)
    print*, '--| ALYA     CHEMIC: EQUILIBRIUM SOLUTION WRITTEN FOR CMC MODEL'
  
  end subroutine chm_inert_eq_CMC


  subroutine chm_compute_initial_fields_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_compute_initial_fields_CMC
     ! NAME 
     !    chm_compute_initial_fields_CMC
     ! DESCRIPTION
     !    CMC fields initialization for bvess_CMC_chm
     ! USED BY
     !    chm_iniunk
     !***
     !-----------------------------------------------------------------------
     use def_parame
     use def_master
     use def_chemic,     only :  nZ_CMC_chm, rscal_inert_CMC_chm, rscal_equil_CMC_chm, &
                                 kfl_weigh_in_eq_CMC_chm, bvess_chm, bvess_CMC_chm, &
                                 kfl_solve_enth_CMC_chm, nvar_CMC_chm, Zavg_CFD_chm, &
                                 Zvar_CFD_chm

     implicit none
     integer(ip)              :: ipoin,imixf,iclas
     real(rp)                 :: aux_dif
     real(rp)                 :: Yk_poin(nZ_CMC_chm,nclas_chm)

     ! Initialization of CMC domain based on a coefficient field and inert and equilibrium solutions
       
     if ( kfl_rstar == 0_ip .and. kfl_weigh_in_eq_CMC_chm == 1_ip) then ! Initial fields defined by a coefficient alpha

        print*, '--| ALYA     CHEMIC: COMPUTING INITIAL FIELDS FOR CMC MODEL'

        ! At this point the only field that is filled in bvess_chm is the first
        ! one with the alpha values

        ! Fill the inner nodes of bvess_CMC_chm matrix (physical + mixt. frac. spaces)
        do ipoin = 1,npoin
           do iclas = 1,nclas_chm
                 ! Mixture fraction boundaries
                 bvess_CMC_chm(1,ipoin,iclas)          = rscal_inert_CMC_chm(1,iclas+1_ip)
                 bvess_CMC_chm(nZ_CMC_chm,ipoin,iclas) = rscal_inert_CMC_chm(nZ_CMC_chm,iclas+1_ip)

                 ! Intermediate mixture fractions
                 do imixf = 2,nZ_CMC_chm-1_ip
                    aux_dif = rscal_equil_CMC_chm(imixf,iclas+1_ip) - rscal_inert_CMC_chm(imixf,iclas+1_ip)
                    bvess_CMC_chm(imixf,ipoin,iclas) = rscal_inert_CMC_chm(imixf,iclas+1_ip) + &
                        bvess_chm(nvar_CMC_chm+1_ip,ipoin) * aux_dif
                 end do
           end do
        end do

        if (kfl_solve_enth_CMC_chm /= 0_ip) then
           do ipoin = 1,npoin
              do imixf = 1, nZ_CMC_chm
                 do iclas = 1, nclas_chm
                    Yk_poin(imixf,iclas) = bvess_CMC_chm(imixf,ipoin,iclas)
                 end do
              end do

              call chm_get_cond_enthalpy_average_CMC(bvess_chm(nvar_CMC_chm,ipoin), &
                       Yk_poin(1:nZ_CMC_chm,1:nclas_chm), Zavg_CFD_chm(ipoin), &
                       Zvar_CFD_chm(ipoin), bvess_chm(nvar_CMC_chm+1_ip,ipoin), &
                       bvess_CMC_chm(1:nZ_CMC_chm,ipoin,nvar_CMC_chm))
           end do
        end if

     else ! Restart or initialize at t=0 all the files provided when previous option is not used
     ! Make possible that if one species does not appear assign a null field to that species
     !!!!!!! FILL

        call runend('Restart option not implemented yet')
        print*, '--| ALYA     CHEMIC: READING INITIAL FIELDS FOR CMC MODEL'
        
     end if

  end subroutine chm_compute_initial_fields_CMC


  subroutine chm_initialization_domain_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_initialization_domain_CMC
     ! NAME 
     !    chm_initialization_domain_CMC
     ! DESCRIPTION
     !    CMC fields initialization
     ! USED BY
     !    chm_iniunk
     !***
     !-----------------------------------------------------------------------
     use def_parame
     use def_master
     use def_chemic,     only :  nZ_CMC_chm, rscal_inert_CMC_chm, rscal_equil_CMC_chm, &
                                 enthalp_CMC_chm, Yk_CMC_chm, temp_CMC_chm, &
                                 bvess_CMC_chm, kfl_solve_enth_CMC_chm, nvar_CMC_chm

     implicit none
     integer(ip)              :: ipoin,imixf,iclas

     ! Fill matrices Yk_CMC_chm and enthalp_CMC_chm
     do ipoin = 1,npoin
        do iclas = 1,nclas_chm
           do imixf = 1,nZ_CMC_chm
              Yk_CMC_chm(imixf,ipoin,iclas) = bvess_CMC_chm(imixf,ipoin,iclas)
           end do
        end do
     end do


     if (kfl_solve_enth_CMC_chm /= 0_ip) then
        do ipoin = 1,npoin
           do imixf = 1,nZ_CMC_chm
              enthalp_CMC_chm(imixf,ipoin) = bvess_CMC_chm(imixf,ipoin,nvar_CMC_chm)
           end do
        end do
     else
        do imixf = 1,nZ_CMC_chm
           enthalp_CMC_chm(imixf,1) = rscal_inert_CMC_chm(imixf,1)
        end do
     end if

#ifdef CANTERA
     ! Compute conditional temperature, viscosity and specific volume at nodes
     do imixf = 1,nZ_CMC_chm

        ! Assign to conce and therm
        do ipoin = 1,npoin
           if (kfl_solve_enth_CMC_chm /= 0_ip) then
              therm(ipoin,1) = enthalp_CMC_chm(imixf,ipoin)
           else
              therm(ipoin,1) = enthalp_CMC_chm(imixf,1)
           end if

           do iclas = 1,nclas_chm
              conce(ipoin,iclas,1) = Yk_CMC_chm(imixf,ipoin,iclas)
           end do

           call cantera_get_temperature_from_hy(prthe(1), therm(ipoin,1), conce(ipoin,1:nclas_chm,1), &
                     temp_CMC_chm(imixf,ipoin)) ! The value for temperature obtained here is used
           ! as initial value in chm_calc_temp_CMC
        end do
        
        call chm_calc_temp_CMC(imixf)
        
        ! Assign to tempe
        do ipoin = 1, npoin
           tempe(ipoin,1) = temp_CMC_chm(imixf,ipoin)
        end do
        
        call chm_calc_densi_visco_gauss_CMC(imixf)
     end do
#endif

     ! Integrations
     call chm_integrate_flow_var_points_CMC
     call chm_integrate_flow_var_gauss_CMC
     call chm_rho_visco_nodal_project_CMC

  end subroutine chm_initialization_domain_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! B O U N D A R Y   C O N D I T I O N S !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_bc_type_species_CMC
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_bc_type_species_CMC
     ! NAME                                                                   
     !    chm_bc_type_species_CMC
     ! DESCRIPTION                                                            
     !    Idenitify if all the species have the same type of boundary conditions
     !    for the physical points of the mesh
     ! USED BY                                                                
     !    chm_iniunk
     !***                                                                     
     !-----------------------------------------------------------------------

     use def_chemic,             only :  kfl_bc_type_spec_CMC_chm, kfl_fixno_chm

     implicit none
     integer(ip)                      :: ipoin, iclas, type_clas1

     do ipoin = 1, npoin
        type_clas1 = kfl_fixno_chm(1,ipoin)
        do iclas = 2, nclas_chm
           if (kfl_fixno_chm(iclas,ipoin) /= type_clas1) then
              kfl_bc_type_spec_CMC_chm(ipoin) = 0
           end if
        end do
     end do

  end subroutine chm_bc_type_species_CMC


  subroutine chm_save_unconditional_fields_bc_CMC
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_save_unconditional_fields_bc_CMC
     ! NAME                                                                   
     !    chm_save_unconditional_fields_bc_CMC
     ! DESCRIPTION                                                            
     !    Assign bvess_chm values to bvess_ufield_CMC_chm when starting
     !    the simulation.
     ! USED BY                                                                
     !    chm_iniunk
     !***                                                                     
     !-----------------------------------------------------------------------

     use def_chemic,             only :  bvess_chm, bvess_ufield_CMC_chm, &
                                         nvar_CMC_chm

     implicit none
     integer(ip)                      :: ipoin, iclas
     
     do ipoin = 1, npoin
        do iclas = 1, nvar_CMC_chm
           bvess_ufield_CMC_chm(ipoin,iclas) = bvess_chm(iclas,ipoin)
        end do
     end do

  end subroutine chm_save_unconditional_fields_bc_CMC


  subroutine chm_get_cond_species_all_average_CMC(ipoin,Yk_avg,Zavg,Zvar,Yk_prof)
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_get_cond_species_all_average_CMC
     ! NAME                                                                   
     !    chm_get_cond_species_all_average_CMC
     ! DESCRIPTION                                                            
     !    Find conditional species profiles along mixture fraction from average
     !    values. Do this operation for all the species at the same time.
     ! USED BY                                                                
     !    chm_updbcs
     !***                                                                     
     !-----------------------------------------------------------------------

     use def_chemic,             only :  nZ_CMC_chm, rscal_inert_CMC_chm, &
                                         rscal_equil_CMC_chm, Yk_CMC_chm

     implicit none
     integer(ip), intent(in)          :: ipoin
     real(rp), intent(in)             :: Yk_avg(nclas_chm)
     real(rp), intent(in)             :: Zavg
     real(rp), intent(in)             :: Zvar
     real(rp), intent(out)            :: Yk_prof(nZ_CMC_chm,nclas_chm)
     integer(ip)                      :: iclas, imixf
     real(rp)                         :: inert_eq_int(2_ip*nclas_chm)
     real(rp)                         :: alpha
     real(rp)                         :: Yk_sum
     real(rp)                         :: diff_aux

     Yk_prof(1:nZ_CMC_chm,1:nclas_chm) = 0.0_rp

     ! Integrate inert and equilibrium profiles
     call chm_mxt_fr_integr_previous_steps_CMC(2_ip*nclas_chm, &
            (/rscal_inert_CMC_chm(1:nZ_CMC_chm,2:nclas_chm+1_ip), &
            rscal_equil_CMC_chm(1:nZ_CMC_chm,2:nclas_chm+1_ip)/), &
            Zavg, Zvar, inert_eq_int(1:2_ip*nclas_chm))

     do iclas = 1, nclas_chm
        diff_aux = inert_eq_int(iclas+nclas_chm) - inert_eq_int(iclas)
        if (diff_aux == 0.0_rp) then
           ! Take the profile from previous time step
           Yk_prof(1:nZ_CMC_chm,iclas) = Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
        else
           alpha = (Yk_avg(iclas) - inert_eq_int(iclas)) / diff_aux
           Yk_prof(1:nZ_CMC_chm,iclas) = rscal_inert_CMC_chm(1:nZ_CMC_chm,iclas+1_ip) + &
              alpha * (rscal_equil_CMC_chm(1:nZ_CMC_chm,iclas+1_ip) - & 
              rscal_inert_CMC_chm(1:nZ_CMC_chm,iclas+1_ip))
        end if
     end do

     ! Normalize in case mass fractions do not sum 1
     do imixf = 2, nZ_CMC_chm-1_ip
        Yk_sum = 0.0_rp
        do iclas = 1, nclas_chm
           Yk_sum = Yk_sum + Yk_prof(imixf,iclas)
        end do
        if (Yk_sum /= 0.0_rp) then
           Yk_prof(imixf,1:nclas_chm) = Yk_prof(imixf,1:nclas_chm) / Yk_sum
        end if
     end do

  end subroutine chm_get_cond_species_all_average_CMC


  
  subroutine chm_get_cond_species_average_CMC(iclas,ipoin,Yk_avg,Zavg,Zvar,Yk_prof)
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_get_cond_field_average_CMC
     ! NAME                                                                   
     !    chm_get_cond_field_average_CMC
     ! DESCRIPTION                                                            
     !    Find conditional species profiles along mixture fraction from average
     !    values. Do this operation for just one single species.
     ! USED BY                                                                
     !    chm_updbcs
     !***                                                                     
     !-----------------------------------------------------------------------

     use def_chemic,             only :  nZ_CMC_chm, rscal_inert_CMC_chm, &
                                         rscal_equil_CMC_chm, Yk_CMC_chm

     implicit none
     integer(ip), intent(in)          :: iclas
     integer(ip), intent(in)          :: ipoin
     real(rp), intent(in)             :: Yk_avg
     real(rp), intent(in)             :: Zavg
     real(rp), intent(in)             :: Zvar
     real(rp), intent(out)            :: Yk_prof(nZ_CMC_chm)
     integer(ip)                      :: imixf
     real(rp)                         :: inert_eq_int(2)
     real(rp)                         :: alpha
     real(rp)                         :: diff_aux

     Yk_prof(1:nZ_CMC_chm) = 0.0_rp

     ! Integrate inert and equilibrium profiles
     call chm_mxt_fr_integr_previous_steps_CMC(2_ip*nclas_chm, &
            (/rscal_inert_CMC_chm(:,iclas+1_ip), rscal_equil_CMC_chm(:,iclas+1_ip)/), &  ! Enthalpy is the first entry
            Zavg, Zvar, inert_eq_int(1:2))

     diff_aux = inert_eq_int(2) - inert_eq_int(1)
     if (diff_aux == 0) then
        ! Take the profile from previous time step
        Yk_prof(1:nZ_CMC_chm) = Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
     else
        alpha = (Yk_avg - inert_eq_int(1)) / diff_aux
        Yk_prof(1:nZ_CMC_chm) = rscal_inert_CMC_chm(1:nZ_CMC_chm,iclas+1_ip) + &
           alpha * (rscal_equil_CMC_chm(1:nZ_CMC_chm,iclas+1_ip) - rscal_inert_CMC_chm(1:nZ_CMC_chm,iclas+1_ip))
     end if

  end subroutine chm_get_cond_species_average_CMC


  subroutine chm_get_cond_enthalpy_average_CMC(temp_uncond,Yk_avg,Zavg,Zvar,alpha,enthalp_prof)
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_get_cond_enthalpy_average_CMC
     ! NAME                                                                   
     !    chm_get_cond_enthalpy_average_CMC
     ! DESCRIPTION                                                            
     !    Find conditional enthalpy profile along mixture fraction from average
     !    temperature and conditional mass fractions.
     ! USED BY                                                                
     !    chm_compute_initial_fields_CMC
     !***                                                                     
     !-----------------------------------------------------------------------

     use def_master,                 only :  prthe
     use def_kermod,                 only :  gasco
     use def_chemic,                 only :  nZ_CMC_chm, temp_inert_CMC_chm, &
                                             temp_equil_CMC_chm

     implicit none
     real(rp), parameter                  :: Tlim = 200.0_rp
     real(rp), intent(in)                 :: temp_uncond
     real(rp), intent(in)                 :: Yk_avg(nZ_CMC_chm,nclas_chm)
     real(rp), intent(in)                 :: Zavg
     real(rp), intent(in)                 :: Zvar
     real(rp), intent(in)                 :: alpha
     real(rp), intent(out)                :: enthalp_prof(nZ_CMC_chm)
     integer(ip)                          :: imixf, iclas
     real(rp)                             :: temp_prof(nZ_CMC_chm)
     real(rp)                             :: temp_adbiab_cond(nZ_CMC_chm,1)
     real(rp)                             :: temp_int(1)
     real(rp)                             :: factor_T
     real(rp)                             :: h_aux

     temp_adbiab_cond(1:nZ_CMC_chm,1) = temp_inert_CMC_chm(1:nZ_CMC_chm) + &
         alpha * (temp_equil_CMC_chm(1:nZ_CMC_chm) - temp_inert_CMC_chm(1:nZ_CMC_chm))

     ! Integrate temperature profile with adiabatic condition
     call chm_mxt_fr_integr_previous_steps_CMC(1_ip, &
            temp_adbiab_cond, Zavg, Zvar, temp_int)

     factor_T = temp_uncond / temp_int(1)

     temp_prof(1:nZ_CMC_chm) = factor_T * temp_adbiab_cond(1:nZ_CMC_chm,1)

     do imixf = 1, nZ_CMC_chm
        if (temp_prof(imixf) < Tlim) &
           call runend('CHEMIC OPERATIONS MOD. CMC: Conditional temperature lower than the minimum possible temperature')
     end do
     
     ! Get enthalpy
     do imixf = 1, nZ_CMC_chm
        call chm_calc_h_from_TY_CMC(Yk_avg(imixf,1:nclas_chm), temp_prof(imixf), &
                h_aux)
        enthalp_prof(imixf) = h_aux
     end do

  end subroutine chm_get_cond_enthalpy_average_CMC


  subroutine chm_get_cond_fields_bc_Dir_CMC
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_get_cond_fields_bc_Dir_CMC
     ! NAME                                                                   
     !    chm_get_cond_fields_bc_Dir_CMC
     ! DESCRIPTION                                                            
     !    Compute conditional fields for points that belong to the boundary
     !    with Dirichlet boundary condition.
     ! USED BY                                                                
     !    chm_updbcs
     !***                                                                     
     !-----------------------------------------------------------------------

     use def_chemic,                only :  kfl_bc_type_spec_CMC_chm, kfl_fixno_chm, &
                                            bvess_ufield_CMC_chm, Zavg_CFD_chm, Zvar_CFD_chm, &
                                            bvess_CMC_chm, Yk_CMC_chm, enthalp_CMC_chm, &
                                            kfl_solve_enth_CMC_chm, rscal_inert_CMC_chm, &
                                            rscal_equil_CMC_chm, nZ_CMC_chm, nvar_CMC_chm, &
                                            index_N2

     implicit none
     real(rp), parameter                 :: small = 1.0e-5_rp
     integer(ip)                         :: ipoin, iclas, imixf, kfl_dir
     real(rp)                            :: Yk_sum, alpha_avg, Yk_sum_aux, diff_aux
     real(rp)                            :: Yk_cond_bc_Dir(nZ_CMC_chm,nclas_chm)
     real(rp)                            :: Yk_poin(nZ_CMC_chm,nclas_chm)
     real(rp)                            :: inert_interm_eq_int(3_ip*nclas_chm)

     Yk_cond_bc_Dir(1:nZ_CMC_chm,1:nclas_chm) = 0.0_rp

     do ipoin = 1,npoin
        if (kfl_bc_type_spec_CMC_chm(ipoin) == 1_ip) then
           if( kfl_fixno_chm(1,ipoin) > 0_ip ) then
              call chm_get_cond_species_all_average_CMC(ipoin,bvess_ufield_CMC_chm(ipoin,1:nclas_chm), &
                     Zavg_CFD_chm(ipoin), Zvar_CFD_chm(ipoin), Yk_cond_bc_Dir(1:nZ_CMC_chm,1:nclas_chm))
              do imixf = 1,nZ_CMC_chm
                 do iclas = 1, nclas_chm
                    bvess_CMC_chm(imixf,ipoin,iclas) = Yk_cond_bc_Dir(imixf,iclas)
                    Yk_CMC_chm(imixf,ipoin,iclas)    = bvess_CMC_chm(imixf,ipoin,iclas)
                 end do   
              end do
           end if
     
        else
           ! Probably never used
           kfl_dir = 0_ip
           do iclas = 1, nclas_chm
              if( kfl_fixno_chm(iclas,ipoin) > 0_ip ) then
                 call chm_get_cond_species_average_CMC(iclas,ipoin,bvess_ufield_CMC_chm(ipoin,iclas), &
                        Zavg_CFD_chm(ipoin), Zvar_CFD_chm(ipoin), Yk_cond_bc_Dir(1:nZ_CMC_chm,iclas))
                 kfl_dir = 1_ip
              end if
           end do

           if (kfl_dir == 1_ip) then
              do iclas = 1, nclas_chm
                 if( kfl_fixno_chm(iclas,ipoin) == 0_ip ) then
                    ! In case for some species we have Dirichlet and for other species Neumann
                    ! (kfl_dir). For Neumann assign values from previous time step
                    Yk_cond_bc_Dir(1:nZ_CMC_chm,iclas) = Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
                 end if
              end do

              ! Normalize in case mass fractions do not sum 1
              do imixf = 2, nZ_CMC_chm-1_ip
                 Yk_sum = 0.0_rp
                 do iclas = 1, nclas_chm
                    Yk_sum = Yk_sum + Yk_cond_bc_Dir(imixf,iclas)
                 end do
                 if (Yk_sum /= 0.0_rp) then
                    Yk_cond_bc_Dir(imixf,1:nclas_chm) = Yk_cond_bc_Dir(imixf,1:nclas_chm) / Yk_sum
                 end if
              end do

              ! Assign to bvess_CMC_chm
              do imixf = 1,nZ_CMC_chm
                 do iclas = 1, nclas_chm
                    bvess_CMC_chm(imixf,ipoin,iclas) = Yk_cond_bc_Dir(imixf,iclas)
                    Yk_CMC_chm(imixf,ipoin,iclas)    = bvess_CMC_chm(imixf,ipoin,iclas)
                 end do
              end do
           end if

        end if
     end do

     if (kfl_solve_enth_CMC_chm /= 0_ip) then
        do ipoin = 1, npoin
           if( kfl_fixno_chm(nvar_CMC_chm,ipoin) > 0_ip ) then
              alpha_avg  = 0.0_rp
              Yk_sum_aux = 0.0_rp

              do imixf = 1, nZ_CMC_chm
                 do iclas = 1, nclas_chm
                    Yk_poin(imixf,iclas) = Yk_CMC_chm(imixf,ipoin,iclas)
                 end do
              end do

              ! Get a parameter that weighes inert and equilibrium solutions
              call chm_mxt_fr_integr_previous_steps_CMC(3_ip*nclas_chm, &
                   (/rscal_inert_CMC_chm(:,2:nclas_chm+1_ip), Yk_poin(1:nZ_CMC_chm,1:nclas_chm), &
                   rscal_equil_CMC_chm(:,2:nclas_chm+1_ip)/), &
                   Zavg_CFD_chm(ipoin), Zvar_CFD_chm(ipoin), inert_interm_eq_int(1:3_ip*nclas_chm))

              do iclas = 1, nclas_chm
                 if (iclas /= index_N2) then  ! N2 is the same in inert and equilibrium conditions except
                    ! for nitrogen oxides which are very low in mass. Better not to account for N2 in the
                    ! following average since YN2 is very large and can lead to spurious results
                    diff_aux = inert_interm_eq_int(iclas+2_ip*nclas_chm) - inert_interm_eq_int(iclas)

                    if (abs(diff_aux) >= small) then
                       alpha_avg = alpha_avg + inert_interm_eq_int(iclas+nclas_chm) * &
                          (inert_interm_eq_int(iclas+nclas_chm) - inert_interm_eq_int(iclas)) / diff_aux
                       Yk_sum_aux = Yk_sum_aux + inert_interm_eq_int(iclas+nclas_chm)
                    end if

                 end if
              end do
 
              alpha_avg = alpha_avg / Yk_sum_aux

              call chm_get_cond_enthalpy_average_CMC(bvess_ufield_CMC_chm(ipoin,nvar_CMC_chm), &
                       Yk_poin(1:nZ_CMC_chm,1:nclas_chm), Zavg_CFD_chm(ipoin), &
                       Zvar_CFD_chm(ipoin), alpha_avg, &
                       bvess_CMC_chm(1:nZ_CMC_chm,ipoin,nvar_CMC_chm))

              do imixf = 1, nZ_CMC_chm
                 do iclas = 1, nclas_chm
                    enthalp_CMC_chm(imixf,ipoin) = bvess_CMC_chm(imixf,ipoin,iclas)
                 end do
              end do
           end if
        end do
     end if

  end subroutine chm_get_cond_fields_bc_Dir_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! S U B R O U T I N E S   F O R   C O M M U N I C A T I O N !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_data_CFD_to_CMC_domain_CMC                                
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_data_CFD_to_CMC_domain_CMC
     ! NAME                                                                   
     !    chm_data_CFD_to_CMC_domain_CMC
     ! DESCRIPTION                                                            
     !
     ! USED BY                                                                
     !    
     !***                                                                     
     !----------------------------------------------------------------------- 
     use def_parame                                                           
     use def_master                                                           
     use def_domain                                                           
                                                                              
     !!DMM 1/ Send CFD variables: Z, Zv, grad_Z, nu_t, u, N                   
     !! Add a small quantity to Zvar to avoid numerical problems. Do this considering
     !! a segregation factor
                                                                              
                                                                              
  end subroutine chm_data_CFD_to_CMC_domain_CMC


  subroutine chm_data_CMC_to_CFD_domain_CMC
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_data_CMC_to_CFD_domain_CMC
     ! NAME                                                                   
     !    chm_data_CMC_to_CFD_domain_CMC
     ! DESCRIPTION                                                            
     !    This routine computes the integrated density and laminar viscosity,
     !    interpolates these CMC fields in CMC mesh into CFD mesh and sends them
     !    to Alya CFD execution.
     ! USED BY                                                                
     !
     !***                                                                     
     !----------------------------------------------------------------------- 
     use def_parame
     use def_master
     use def_domain
     use def_kermod
     use def_chemic

     implicit none

     !
     ! INTERPOLATIONS. Not available yet
     !

     !
     ! Send quantities to CFD
     !

  end subroutine chm_data_CMC_to_CFD_domain_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! T R A N S F E R   L O C A L   <--->   G L O B A L !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_global2local_CMC(imixf)
     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/chm_global2local_CMC
     ! NAME 
     !    chm_global2local_CMC
     ! DESCRIPTION
     !    This routine transfers information from Yk_CMC_chm, enthalp_CMC_chm,
     !    etc. to conce, therm, etc.
     ! USED BY
     !***
     !-----------------------------------------------------------------------

     use def_master,         only : therm, tempe, conce
     use def_chemic,         only : enthalp_CMC_chm, temp_CMC_chm, Yk_CMC_chm, &
                                    kfl_solve_enth_CMC_chm

     implicit none
     integer(ip)                 :: ipoin, iclas, imixf

     if( INOTMASTER ) then
        ! Assign conce, therm and tempe
        do ipoin = 1,npoin
           tempe(ipoin,1) = temp_CMC_chm(imixf,ipoin)  ! Used to obtain properties (chm_cp_k_gauss_CMC)
           do iclas = 1,nclas_chm
              conce(ipoin,iclas,1) = Yk_CMC_chm(imixf,ipoin,iclas)
           end do
        end do

        if (kfl_solve_enth_CMC_chm == 0_ip) then
           therm(1:npoin,1) = enthalp_CMC_chm(imixf,1)
        else
           do ipoin = 1,npoin
              therm(ipoin,1) = enthalp_CMC_chm(imixf,ipoin)
           end do
        end if
     end if

  end subroutine chm_global2local_CMC


  subroutine chm_local2global_CMC(imixf)
     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/chm_local2global_CMC
     ! NAME 
     !    chm_local2global_CMC
     ! DESCRIPTION
     !    This routine transfers information from conce, therm, etc. to Yk_CMC_chm,
     !    enthalp_CMC_chm, etc.
     ! USED BY
     !***
     !-----------------------------------------------------------------------

     use def_master,         only : therm, tempe, conce
     use def_chemic,         only : enthalp_CMC_chm, temp_CMC_chm, Yk_CMC_chm, &
                                    kfl_solve_enth_CMC_chm

     implicit none
     integer(ip)                 :: ipoin, iclas, imixf

     if( INOTMASTER ) then
        ! Assign conce, therm and tempe
        do ipoin = 1,npoin
           do iclas = 1,nclas_chm
              Yk_CMC_chm(imixf,ipoin,iclas) = conce(ipoin,iclas,1)
           end do
        end do

        if (kfl_solve_enth_CMC_chm /= 0_ip) then
           do ipoin = 1,npoin
              enthalp_CMC_chm(imixf,ipoin) = therm(ipoin,1)
           end do
        end if
     end if

  end subroutine chm_local2global_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! T E R M S   M O D E L L I N G   A N D   A S S E M B L Y !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_element_operations_CMC(order,pnode,pgaus,list_elements,imixf)
     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/chm_element_operations_CMC
     ! NAME 
     !    chm_element_operations_CMC
     ! DESCRIPTION
     !    This routine computes the conditioned variables not solved in CMC,
     !    that is, the terms that appear in the CMC transport equations and
     !    have to be modelled for a given mixture fraction level. Then, it
     !    assembles the terms in order to solve CMC transport
     !    equations for all the species and one mixture fraction level.
     ! USED BY
     !    chm_elmope_all
     !***
     !-----------------------------------------------------------------------
     use def_parame
     use def_master
     use def_domain
     use mod_ker_proper
     use def_kermod
     use def_chemic,             only : Z_CMC_chm, Zs_CMC_chm, &
                                        Xnormalized_prof_CMC_chm, dt_rho_chm, &
                                        condu_gp_CMC_chm, sphec_gp_CMC_chm, &
                                        spvol_gp_CMC_chm, &
                                        ADR_chm, kfl_taust_chm, kfl_advec_chm, &
                                        kfl_ellen_chm, kfl_entropy_chm, nvar_CMC_chm, &
                                        kfl_solve_enth_CMC_chm, Le_k

     use mod_matrix
     use def_solver

     use mod_ADR,                only : ADR_element_assembly
     use mod_ADR,                only : ADR_projections_and_sgs_assembly
     use mod_ADR,                only : ADR_add_sgs_or_bubble
     use mod_ADR,                only : ELEMENT_ASSEMBLY                 ! 1
     use mod_ADR,                only : PROJECTIONS_AND_SGS_ASSEMBLY     ! 4
     use mod_ADR,                only : BUBBLE_ASSEMBLY                  ! 5
     use mod_ADR,                only : mreac_adr
     use mod_chm_entropy,        only : chm_entropy_viscosity

     implicit none
     integer(ip), intent(in)          :: order             ! =1 defaul or =2 compute SGS only
     integer(ip), intent(in)          :: pnode             ! Number of nodes
     integer(ip), intent(in)          :: pgaus             ! Number of Gauss points
     integer(ip), intent(in), pointer :: list_elements(:)  ! List of elements
     integer(ip), intent(in)          :: imixf             ! Level of mixture fraction

     integer(ip)              :: iclas,kelem,ielem,igaus,dummi
     integer(ip)              :: pelty,porde,ptopo,plapl,izmat,izrhs

     real(rp)                 :: PDF_val
     real(rp)                 :: chale(3),chave(3),hleng(3),tragl(9)
     real(rp)                 :: dummr(mgaus*ndime), kcp_la

     ! Matrices for elements (values at nodes)
     real(rp)                 :: elmat(mnode,mnode)                     ! Elemental matrix contribution
     real(rp)                 :: elrhs(mnode)                           ! Right hand side
     real(rp)                 :: elcod(ndime,mnode)                     ! Coordinates
     real(rp)                 :: elvel_CFD(ndime,mnode)                 ! Velocity from CFD
     real(rp)                 :: elZavg_CFD(mnode)                      ! Average mixture fraction from CFD
     real(rp)                 :: elZvar_CFD(mnode)                      ! Mixture fraction variance from CFD
     real(rp)                 :: elXtot_CFD(mnode)                      ! Total scalar dissipation rate from CFD
     real(rp)                 :: elZgrad_CFD(ndime,mnode)               ! Average mixture fraction gradient from CFD
     real(rp)                 :: elturb_dif_CFD(mnode)                  ! Mass turbulent diffusion coefficient from CFD
     real(rp)                 :: elderiv2_unkno_CMC(mnode,nvar_CMC_chm) ! Second derivative in mixture fraction direction for Yk and species
     real(rp)                 :: elunkno_CMC(mnode,nvar_CMC_chm)        ! Yk and enthalpy at iZ and nodes

     ! Matrices for elements (values at Gaussian points)
     real(rp)                 :: gpvol(mgaus)                           ! |J|*w
     real(rp)                 :: gpcar(ndime,mnode,mgaus)               ! dNk/dxj
     real(rp)                 :: gphes(ntens,mnode,mgaus)               ! d2Nk/dxidxj
     real(rp)                 :: gpdif(mgaus,nvar_CMC_chm)              ! D_k
     real(rp)                 :: gprea(mgaus,mreac_adr)                 ! r
     real(rp)                 :: gpgrd(ndime,mgaus)                     ! grad(k) = grad(D_k)
     real(rp)                 :: gpvel_CFD(ndime,mgaus)                 ! Velocity
     real(rp)                 :: gpZavg_CFD(mgaus)                      ! Average mixture fraction
     real(rp)                 :: gpZvar_CFD(mgaus)                      ! Mixture fraction variance
     real(rp)                 :: gpXtot_CFD(mgaus)                      ! Total scalar dissipation rate
     real(rp)                 :: gpZgrad_CFD(ndime,mgaus)               ! Average mixture fraction gradient
     real(rp)                 :: gpturb_dif_CFD(mgaus)                  ! Mass turbulent diffusion coefficient
     real(rp)                 :: gplam_dif_CMC(mgaus,nvar_CMC_chm)      ! Mass laminar diffusion coefficient
     real(rp)                 :: del_gpdif(mgaus,nvar_CMC_chm)          ! change in diffusivity from entropy stable stabilization
     real(rp)                 :: gp_densi_PDF_CMC(mgaus)                ! Uncond. rho * probability density function
     real(rp)                 :: gp_veloc_CMC_chm(1:ndime,mgaus)        ! Conditional velocity
     real(rp)                 :: gp_diff_phys_spc(mgaus,nvar_CMC_chm)   ! Diffusion in physical space
     real(rp)                 :: gp_diff_mf(mgaus,nvar_CMC_chm)         ! Diffusion in mixture fraction direction
     real(rp)                 :: gpderiv2_unkno_CMC(mgaus,nvar_CMC_chm) ! Second derivative in mixture fraction direction for Yk
     real(rp)                 :: gpXcond_CMC                            ! Conditional scalar dissipation rate
     real(rp)                 :: gpPDF_param(7)                         ! PDF parameters at Gaussian points
     real(rp)                 :: factor0_vel(ndime)                     ! Constant in the conditional velocity model
     real(rp)                 :: factor1_vel(ndime)                     ! Slope in the conditional velocity model
     real(rp)                 :: gpX0_CMC                               ! Modulator scalar dissipation rate for AMC model
     real(rp)                 :: gpunkno_CMC(mgaus,nvar_CMC_chm)        ! Yk and enthalpy at iZ and Gaussian points



     ! Loop over all the elements
     elements: do kelem = 1, size(list_elements)
        ielem = list_elements(kelem)
        if( ielem > 0 ) then
           !
           ! Step 1: find values at Gaussian points for: veloc_CFD_chm, Zavg_CFD_chm,
           ! Zvar_CFD_chm, Xtot_CFD_chm, grad_Zavg_CFD_chm, turbulent diffusion coefficient, density
           !

           !
           ! Element dimensions
           !
           pelty = ltype(ielem)
           porde = lorde(pelty)
           ptopo = ltopo(pelty)

           !
           ! Initialization of variables in Gaussian points
           !

           elmat(1:mnode,1:mnode)                     = 0.0_rp
           elrhs(1:mnode)                             = 0.0_rp
           elcod(1:ndime,1:mnode)                     = 0.0_rp
           elvel_CFD(1:ndime,1:mnode)                 = 0.0_rp
           elZavg_CFD(1:mnode)                        = 0.0_rp
           elZvar_CFD(1:mnode)                        = 0.0_rp
           elXtot_CFD(1:mnode)                        = 0.0_rp
           elZgrad_CFD(1:ndime,1:mnode)               = 0.0_rp
           elturb_dif_CFD(1:mnode)                    = 0.0_rp
           elderiv2_unkno_CMC(1:mnode,1:nvar_CMC_chm) = 0.0_rp
           elunkno_CMC(1:mnode,1:nvar_CMC_chm)        = 0.0_rp
           gpvol(1:mgaus)                             = 0.0_rp
           gpcar(1:ndime,1:mnode,1:mgaus)             = 0.0_rp
           gphes(1:ntens,1:mnode,1:mgaus)             = 0.0_rp
           gpvel_CFD(1:ndime,1:mgaus)                 = 0.0_rp
           gpZavg_CFD(1:mgaus)                        = 0.0_rp
           gpZvar_CFD(1:mgaus)                        = 0.0_rp
           gpXtot_CFD(1:mgaus)                        = 0.0_rp
           gpZgrad_CFD(1:ndime,1:mgaus)               = 0.0_rp
           gplam_dif_CMC(1:mgaus,1:nvar_CMC_chm)      = 0.0_rp
           gpturb_dif_CFD(1:mgaus)                    = 0.0_rp
           gp_densi_PDF_CMC                           = 0.0_rp
           gp_veloc_CMC_chm(1:ndime,1:mgaus)          = 0.0_rp
           gpdif(1:mgaus,1:nvar_CMC_chm)              = 0.0_rp
           gp_diff_mf(1:mgaus,1:nvar_CMC_chm)         = 0.0_rp
           gp_diff_phys_spc(1:mgaus,1:nvar_CMC_chm)   = 0.0_rp
           gpderiv2_unkno_CMC(1:mgaus,1:nvar_CMC_chm) = 0.0_rp
           gpXcond_CMC                                = 0.0_rp
           gpX0_CMC                                   = 0.0_rp
           gpPDF_param(1:7)                           = 0.0_rp
           gpPDF_param(4)                             = Zs_CMC_chm
           gprea(1:mgaus,1:mreac_adr)                 = 0.0_rp
           gpgrd(1:ndime,1:mgaus)                     = 0.0_rp
           gpunkno_CMC(1:mgaus,1:nvar_CMC_chm)        = 0.0_rp
           del_gpdif(1:mgaus,1:nvar_CMC_chm)          = 0.0_rp
           factor0_vel(1:ndime)                       = 0.0_rp
           factor1_vel(1:ndime)                       = 0.0_rp
           plapl                                      = 0_ip

           !
           ! Gather values at the element. The turbulent diffusion coef. got in elturb_dif_CFD is for mass
           !
            call chm_elmgac_CMC(&
                   pnode, lnods(1:pnode,ielem), elcod, elvel_CFD, elZavg_CFD, elZvar_CFD,&
                   elXtot_CFD, elZgrad_CFD, elturb_dif_CFD, elderiv2_unkno_CMC, elunkno_CMC, imixf)

           !
           ! CHALE, HLENG and TRAGL 
           !
           if( kfl_taust_chm /= 0 ) then
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
              call elmchl(&
                   tragl,hleng,elcod,dummr,chave,chale,pnode,porde,hnatu(pelty),&
                   kfl_advec_chm,kfl_ellen_chm)
           else
              plapl = 0
           end if

           !
           ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
           !
           call elmcar(&
               pnode, pgaus, plapl, elmar(pelty)%weigp, elmar(pelty)%shape,&
               elmar(pelty)%deriv, elmar(pelty)%heslo, elcod, gpvol, gpcar,&
               gphes, ielem)

           !
           ! Send quantities to Gaussian points
           !
           call chm_elmpre_CMC(&
                   pnode, pgaus, elmar(pelty)%shape, elvel_CFD, elZavg_CFD, &
                   elZvar_CFD, elXtot_CFD, elZgrad_CFD, elturb_dif_CFD, elderiv2_unkno_CMC, &
                   elunkno_CMC, gpvel_CFD, gpZavg_CFD, gpZvar_CFD, gpXtot_CFD, gpZgrad_CFD, &
                   gpturb_dif_CFD, gpderiv2_unkno_CMC, gpunkno_CMC)

           !
           ! Compute the laminar diffusion coefficient D at Gaussian points
           !

           do igaus = 1, pgaus
              kcp_la = condu_gp_CMC_chm(ielem)%a(igaus,1,1) * spvol_gp_CMC_chm(ielem)%a(imixf,igaus,1) / &
                                 sphec_gp_CMC_chm(ielem) % a(igaus,1,1)
              do iclas = 1,nclas_chm
                 gpdif(igaus,iclas) = kcp_la / Le_k(iclas)  ! It contains rho*D and not only D
              end do
              if (kfl_solve_enth_CMC_chm /= 0_ip) then  ! In case enthalpy is transported
                 gpdif(igaus,nvar_CMC_chm) = kcp_la
              end if
           end do

           !
           ! Entropy stable viscosity. !!!!!!!!!!!!!!!!NEED TO BE INCLUDED?
           !
           !!!!!!!!!!!!!!!!!! CHECK THE FOLLOWING LINES
           if(kfl_entropy_chm == 1_ip) then
              del_gpdif = 0.0_rp
              do iclas = 1,nvar_CMC_chm
              !   call chm_entropy_viscosity(ielem,pnode,pgaus,1_ip,pgaus,iclas,elmar(pelty)%shape,gpcar,elvel,gpden,hleng,gpvol,del_gpdif)
              enddo

              !   do iclas = iclai_chm,iclaf_chm
              !      gpdif(1:pgaus,iclas) = gpdif(1:pgaus,iclas) + del_gpdif(1:pgaus,iclas) 
              !   enddo ! At this point gpdif contains the laminar and entropy viscosity contributions
           end if


           gauss: do igaus = 1, pgaus
              !
              ! Find values at Gaussian points
              !
              call chm_find_PDF_parameters_CMC(gpZavg_CFD(igaus), gpZvar_CFD(igaus), gpPDF_param(1:7))

              call chm_find_veloc_parameters_CMC(gpvel_CFD(1:ndime,igaus), gpZavg_CFD(igaus), &
                        gpZvar_CFD(igaus), gpZgrad_CFD(1:ndime,igaus), gpturb_dif_CFD(igaus), &
                        factor0_vel(1:ndime), factor1_vel(1:ndime))

              ! Scalar dissipation rate
              call chm_find_scalar_dissip_rate_X0_CMC(gpZavg_CFD(igaus), gpZvar_CFD(igaus), gpXtot_CFD(igaus), gpX0_CMC)

              ! Find pdf value
              PDF_val = gpPDF_param(3) * Z_CMC_chm(imixf)**(gpPDF_param(1)-1.0_rp) * &
                            (gpPDF_param(4) - Z_CMC_chm(imixf))**(gpPDF_param(2)-1.0_rp)

              ! Find density times PDF
              gp_densi_PDF_CMC(igaus) = densi_gp(ielem) % a(igaus,1,1) * PDF_val


              ! Diffusion in physical space
              !!!!! MODIFY? THE FOLLOWING CODE DEPENDS ON IF WE HAVE rho*D_turb or D_turb
              do iclas = 1,nvar_CMC_chm
                 gp_diff_phys_spc(igaus,iclas) = gp_densi_PDF_CMC(igaus) * &
                                                  (gpdif(igaus,iclas) + gpturb_dif_CFD(igaus))
              end do

              ! Conditional velocity
              gp_veloc_CMC_chm(1:ndime,igaus) = factor0_vel(1:ndime) + factor1_vel(1:ndime) * Z_CMC_chm(imixf)

              !
              ! Diffusive term in mixture fraction direction
              !
              ! Density * pdf * scalar dissip. term * d2Y_k / dZ2
              gp_diff_mf(igaus,1:nvar_CMC_chm) = gp_densi_PDF_CMC(igaus) * gpX0_CMC * &
                                                   Xnormalized_prof_CMC_chm(imixf) * gpderiv2_unkno_CMC(igaus,1:nvar_CMC_chm)

           end do gauss


           ! Once values at all Gaussian points are computed, do assembly, if required
           !
           ! Projections of rho/dt and 1/dt
           !
           call chm_rhodt(  &
                pnode,pgaus,porde,lnods(:,ielem),elmar(pelty)%shape,gpvol, &
                gp_densi_PDF_CMC,dt_rho_chm)

           !
           ! Assemble matrix
           !  
           izmat = 1
           izrhs = 1


           ASSEMBLY_ICLAS: do iclas = 1,nvar_CMC_chm
              ! 
              ! Assemble equation for iclas
              !
              if( order == ELEMENT_ASSEMBLY ) then

                 call ADR_element_assembly(&
                      ielem,pnode,pgaus,elcod,elmar(pelty)%shape,gpcar,elmar(pelty)%deriv,gphes,gpvol,chale,&
                      elmar(pelty)%shape_bub,elmar(pelty)%deriv_bub,ADR_chm(iclas),cutim, &
                      gp_densi_PDF_CMC(1:pgaus),gp_veloc_CMC_chm(1:ndime,1:pgaus), gp_diff_phys_spc(1:pgaus,iclas), &
                      gpgrd,gprea, gp_diff_mf(1:pgaus,iclas), gpunkno_CMC(1:pgaus,iclas), &
                      elunkno_CMC(1:pnode,iclas),elmat,elrhs)

                 
                 !
                 ! Prescribe Dirichlet boundary conditions
                 !
                 if( solve(1) % kfl_iffix == 0 ) &
                      call chm_elmdir(&
                      iclas,pnode,lnods(1,ielem),elmat,elrhs)

                 !
                 ! Call solver
                 !
                 call matrix_assexp(solve(1)%ndofn,1_ip,pnode,npoin,lnods(1:pnode,ielem),elrhs,elmat, &
                                    elunkno_CMC(1:pnode,iclas),rhsid,iclas)

                 izrhs = izrhs + npoin                                       !solve(1)%nzrhs
                 izmat = izmat + solve(1)%nzmat/nvar_CMC_chm**2_ip
              end if

           end do ASSEMBLY_ICLAS

        end if
     end do elements

  end subroutine chm_element_operations_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! P O I N T S   T R A N S F E R !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_elmgac_CMC(&
                 pnode,lnods,elcod,elvel_CFD,elZavg_CFD,elZvar_CFD,&
                 elXtot_CFD,elZgrad_CFD,elturb_dif_CFD,&
                 elderiv2_unkno_CMC, elunkno_CMC, imixf)

     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_elmgac_CMC
     ! NAME 
     !    chm_elmgac_CMC
     ! DESCRIPTION
     !    Gather operations for the set of not solved variables that do not 
     !    depend on the mixture fraction for CMC model.
     ! USES
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !------------------------------------------------------------------------

     use def_chemic,      only: veloc_CFD_chm, Zavg_CFD_chm, Zvar_CFD_chm, &
                                Xtot_CFD_chm, grad_Zavg_CFD_chm, &
                                visco_turb_CFD_chm, Yk_CMC_chm, deriv2_Yk_CMC_chm, &
                                deriv2_enthalp_CMC_chm, enthalp_CMC_chm, nvar_CMC_chm, &
                                kfl_solve_enth_CMC_chm, diffu_chm

     implicit none
     integer(ip), intent(in)  :: pnode
     integer(ip), intent(in)  :: lnods(pnode)
     integer(ip), intent(in)  :: imixf
     real(rp),    intent(out) :: elcod(ndime,pnode)
     real(rp),    intent(out) :: elvel_CFD(ndime,pnode)
     real(rp),    intent(out) :: elZavg_CFD(pnode)
     real(rp),    intent(out) :: elZvar_CFD(pnode)
     real(rp),    intent(out) :: elXtot_CFD(pnode)
     real(rp),    intent(out) :: elZgrad_CFD(ndime,pnode)
     real(rp),    intent(out) :: elturb_dif_CFD(pnode)
     real(rp),    intent(out) :: elderiv2_unkno_CMC(pnode,nvar_CMC_chm)
     real(rp),    intent(out) :: elunkno_CMC(pnode,nvar_CMC_chm)

     integer(ip)              :: inode, ipoin

     !
     ! Initialization
     !
     elcod(1:ndime,1:pnode)                     = 0.0_rp
     elvel_CFD(1:ndime,1:pnode)                 = 0.0_rp
     elZavg_CFD(1:pnode)                        = 0.0_rp
     elZvar_CFD(1:pnode)                        = 0.0_rp
     elXtot_CFD(1:pnode)                        = 0.0_rp
     elZgrad_CFD(1:ndime,1:pnode)               = 0.0_rp
     elturb_dif_CFD(1:pnode)                    = 0.0_rp
     elderiv2_unkno_CMC(1:pnode,1:nvar_CMC_chm) = 0.0_rp
     elunkno_CMC(1:pnode,1:nvar_CMC_chm)        = 0.0_rp
 
     !
     ! Values transference to elemental matrices
     !
     do inode = 1,pnode
        ipoin                                  = lnods(inode)
        elcod(1:ndime,inode)                   = coord(1:ndime,ipoin)
        elvel_CFD(1:ndime,inode)               = veloc_CFD_chm(1:ndime,ipoin)
        elZavg_CFD(inode)                      = Zavg_CFD_chm(ipoin)
        elZvar_CFD(inode)                      = Zvar_CFD_chm(ipoin)
        elXtot_CFD(inode)                      = Xtot_CFD_chm(ipoin)
        elZgrad_CFD(1:ndime,inode)             = grad_Zavg_CFD_chm(1:ndime,ipoin)
        elturb_dif_CFD(inode)                  = visco_turb_CFD_chm(ipoin) / diffu_chm(1,1)
        elderiv2_unkno_CMC(inode,1:nclas_chm)  = deriv2_Yk_CMC_chm(imixf,ipoin,1:nclas_chm)
        elunkno_CMC(inode,1:nclas_chm)         = Yk_CMC_chm(imixf,ipoin,1:nclas_chm)
     end do

     if (kfl_solve_enth_CMC_chm /= 0_ip) then   ! In case enthalpy is transported
        do inode = 1,pnode
           ipoin                                  = lnods(inode)
           elderiv2_unkno_CMC(inode,nvar_CMC_chm) = deriv2_enthalp_CMC_chm(imixf,ipoin)
           elunkno_CMC(inode,nvar_CMC_chm)        = enthalp_CMC_chm(imixf,ipoin)
        end do
     end if

  end subroutine chm_elmgac_CMC



  subroutine chm_elmpre_CMC(&
                   pnode, pgaus, gpsha, elvel_CFD, elZavg_CFD, elZvar_CFD, elXtot_CFD, &
                   elZgrad_CFD, elturb_dif_CFD, elderiv2_unkno_CMC, elunkno_CMC, &
                   gpvel_CFD, gpZavg_CFD, gpZvar_CFD, gpXtot_CFD, gpZgrad_CFD, &
                   gpturb_dif_CFD, gpderiv2_unkno_CMC, gpunkno_CMC)

     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_elmpre_CMC
     ! NAME 
     !    chm_elmpre_CMC
     ! DESCRIPTION
     !    Transfer elemental values to Gaussian points for the set of not 
     !    solved variables that do not depend on the mixture fraction for
     !    CMC model.
     ! USES
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !------------------------------------------------------------------------ 

     use def_chemic,     only :  nvar_CMC_chm

     implicit none
     integer(ip), intent(in)  :: pnode, pgaus
     real(rp),    intent(in)  :: gpsha(pnode,pgaus)
     real(rp),    intent(in)  :: elvel_CFD(ndime,pnode)
     real(rp),    intent(in)  :: elZavg_CFD(pnode)
     real(rp),    intent(in)  :: elZvar_CFD(pnode)
     real(rp),    intent(in)  :: elXtot_CFD(pnode)
     real(rp),    intent(in)  :: elZgrad_CFD(ndime,pnode)
     real(rp),    intent(in)  :: elturb_dif_CFD(pnode)
     real(rp),    intent(in)  :: elderiv2_unkno_CMC(pnode,nvar_CMC_chm)
     real(rp),    intent(in)  :: elunkno_CMC(pnode,nvar_CMC_chm)
     real(rp),    intent(out) :: gpvel_CFD(ndime,pgaus)
     real(rp),    intent(out) :: gpZavg_CFD(pgaus)
     real(rp),    intent(out) :: gpZvar_CFD(pgaus)
     real(rp),    intent(out) :: gpXtot_CFD(pgaus)
     real(rp),    intent(out) :: gpZgrad_CFD(ndime,pgaus)
     real(rp),    intent(out) :: gpturb_dif_CFD(pgaus)
     real(rp),    intent(out) :: gpderiv2_unkno_CMC(pgaus,nvar_CMC_chm)
     real(rp),    intent(out) :: gpunkno_CMC(pgaus,nvar_CMC_chm)

     integer(ip)              :: inode, igaus

     do igaus = 1,pgaus
        do inode = 1,pnode
           gpvel_CFD(1:ndime,igaus)                 = gpvel_CFD(1:ndime,igaus) &
                                                       + gpsha(inode,igaus) * elvel_CFD(1:ndime,inode)
           gpZavg_CFD(igaus)                        = gpZavg_CFD(igaus) &
                                                       + gpsha(inode,igaus) * elZavg_CFD(inode)
           gpZvar_CFD(igaus)                        = gpZvar_CFD(igaus) &
                                                       + gpsha(inode,igaus) * elZvar_CFD(inode)
           gpXtot_CFD(igaus)                        = gpXtot_CFD(igaus) &
                                                       + gpsha(inode,igaus) * elXtot_CFD(inode)
           gpZgrad_CFD(1:ndime,igaus)               = gpZgrad_CFD(1:ndime,igaus) &
                                                       + gpsha(inode,igaus) * elZgrad_CFD(1:ndime,inode)
           gpturb_dif_CFD(igaus)                    = gpturb_dif_CFD(igaus) &
                                                       + gpsha(inode,igaus) * elturb_dif_CFD(inode)
           gpderiv2_unkno_CMC(igaus,1:nvar_CMC_chm) = gpderiv2_unkno_CMC(igaus,1:nvar_CMC_chm) &
                                                       + gpsha(inode,igaus) * elderiv2_unkno_CMC(inode,1:nvar_CMC_chm)
           gpunkno_CMC(igaus,1:nvar_CMC_chm)        = gpunkno_CMC(igaus,1:nvar_CMC_chm) &
                                                       + gpsha(inode,igaus) * elunkno_CMC(inode,1:nvar_CMC_chm)
        end do
     end do

  end subroutine chm_elmpre_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! T E R M S   M O D E L L I N G !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_find_veloc_parameters_CMC(gpvel_CFD, gpZavg_CFD, gpZvar_CFD, &
                        gpZgrad_CFD, gpturb_dif_CFD, factor0, factor1)

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/chm_calc_conditional_veloc_CMC
     ! NAME 
     !    chm_calc_conditional_veloc_CMC
     ! DESCRIPTION
     !    It computes the conditional velocity assuming a normal joint distribution
     !    between mixture fraction and velocity.
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,      only: Zs_CMC_chm

     implicit none
     real(rp),    parameter   :: S_threshold = 1.0e-2_rp
     real(rp),    intent(in)  :: gpvel_CFD(ndime)
     real(rp),    intent(in)  :: gpZavg_CFD
     real(rp),    intent(in)  :: gpZvar_CFD
     real(rp),    intent(in)  :: gpZgrad_CFD(ndime)
     real(rp),    intent(in)  :: gpturb_dif_CFD
     real(rp),    intent(out) :: factor0(ndime), factor1(ndime)

     real(rp)                 :: aux(ndime), S

     S = gpZvar_CFD / (gpZavg_CFD*(Zs_CMC_chm-gpZavg_CFD))

     if (S <= S_threshold) then
        ! In this case we take <v|mixf> = <v> because if Zvar is very small grad(Z) will be
        ! very small but the ratio can be misleading
        factor0(1:ndime) = gpvel_CFD(1:ndime)
        factor1(1:ndime) = 0.0_rp
     else
        aux(1:ndime)     = gpturb_dif_CFD * gpZgrad_CFD(1:ndime) / gpZvar_CFD
        factor0(1:ndime) = gpvel_CFD(1:ndime) + aux(1:ndime) *  gpZavg_CFD
        factor1(1:ndime) = - aux(1:ndime)
     end if

  end subroutine chm_find_veloc_parameters_CMC


  subroutine chm_find_scalar_dissip_rate_X0_CMC(Zavg, Zvar, Xtot, X0)

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/chm_calc_scalar_dissip_rate_CMC
     ! NAME 
     !    chm_calc_scalar_dissip_rate_CMC
     ! DESCRIPTION
     !    It computes the scalar dissipation rate at Gaussian points.
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,      only: Zs_CMC_chm, Z_AMC_CMC_chm, S_AMC_CMC_chm, &
                                nZ_AMC_CMC_chm, nS_AMC_CMC_chm, Smax_AMC_CMC_chm, &
                                Xintegrated_table_AMC_CMC_chm

     implicit none
     real(rp), parameter      :: small = 1.0e-4_rp

     real(rp),    intent(in)  :: Zavg
     real(rp),    intent(in)  :: Zvar
     real(rp),    intent(in)  :: Xtot
     real(rp),    intent(out) :: X0

     integer(ip)              :: index_mf1, index_mf2, index_S1, index_S2
     real(rp)                 :: S, denominator

     if (Zavg <= small .or. Zavg >= (Zs_CMC_chm-small)) then
        ! Both Xtot and denominator will be very small so to avoid numerical noise 0 is assigned
        X0 = 0.0_rp
     else
        S = Zvar / (Zavg * (Zs_CMC_chm - Zavg))

        call get_index_vector(0.0_rp, Zs_CMC_chm, nZ_AMC_CMC_chm, 1.0_rp, Zavg, index_mf1)
        index_mf2 = index_mf1 + 1_ip

        if (S >= Smax_AMC_CMC_chm) then
           index_S2 = nS_AMC_CMC_chm
           index_S1 = index_S2 - 1_ip
        else
           call get_index_vector(0.0_rp, Smax_AMC_CMC_chm, nS_AMC_CMC_chm, 1.0_rp, S, index_S1)
           index_S2 = index_S1 + 1_ip
        end if

       ! Bilinear interpolation for AMC denominator
        call bilinear_intepolation((/Z_AMC_CMC_chm(index_mf1), Z_AMC_CMC_chm(index_mf2), S_AMC_CMC_chm(index_S1), S_AMC_CMC_chm(index_S2)/), &
               (/Xintegrated_table_AMC_CMC_chm(index_mf1,index_S1), Xintegrated_table_AMC_CMC_chm(index_mf1, index_S2), &
               Xintegrated_table_AMC_CMC_chm(index_mf2,index_S1), Xintegrated_table_AMC_CMC_chm(index_mf2,index_S2)/), &
               Zavg, S, denominator)

        if (denominator > small) then  ! The maximum value for the normalized X profile is 1, 
                                       ! then small is a good estimator to discern non-meaningful values
           !!!!!! CHECK NUMBER 2
           X0 = Xtot / (2.0_rp * denominator)
        else
           X0 = 0.0_rp
        end if
     end if

  end subroutine chm_find_scalar_dissip_rate_X0_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!
  !!! A M C   M O D E L !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_AMC_generate_Z_S_vectors_CMC

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/chm_AMC_generate_Z_S_vectors_CMC
     ! NAME 
     !    chm_AMC_generate_Z_S_vectors_CMC
     ! DESCRIPTION
     !    It generates the vectors for mixture fraction and segregation factor
     !    that define the table where the integrals for AMC model are saved.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,     only:  nZ_AMC_CMC_chm, nS_AMC_CMC_chm, Smax_AMC_CMC_chm, &
                                Zs_CMC_chm, Z_AMC_CMC_chm, S_AMC_CMC_chm

     implicit none
     real(rp)                :: delta_Z, delta_S
     integer(ip)             :: iZ, iS

     delta_Z = Zs_CMC_chm / (nZ_AMC_CMC_chm-1)
     Z_AMC_CMC_chm(nZ_AMC_CMC_chm) = Zs_CMC_chm
     do iZ = 2,nZ_AMC_CMC_chm-1_ip
        Z_AMC_CMC_chm(iZ) = (iZ-1) * delta_Z
     end do

     delta_S = Smax_AMC_CMC_chm / (nS_AMC_CMC_chm-1)
     S_AMC_CMC_chm(nS_AMC_CMC_chm) = Smax_AMC_CMC_chm
     do iS = 2,nS_AMC_CMC_chm-1_ip
        S_AMC_CMC_chm(iS) = (iS-1) * delta_S
     end do

  end subroutine chm_AMC_generate_Z_S_vectors_CMC


  subroutine chm_AMC_integrals_CMC

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/chm_AMC_integrals_CMC
     ! NAME 
     !    chm_AMC_integrals_CMC
     ! DESCRIPTION
     !    It generates a matrix that contains the integrals for AMC model.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,     only:  nZ_AMC_CMC_chm, nS_AMC_CMC_chm, Zs_CMC_chm, &
                                Xintegrated_table_AMC_CMC_chm, Xnormalized_prof_CMC_chm, &
                                Z_AMC_CMC_chm, S_AMC_CMC_chm, S_threshold

     implicit none
     integer(ip)             :: iZ, iS
     real(rp)                :: Zvar, PDF_param(7), aux(1)
     character(len=20)       :: names_rscal(3)

     PDF_param(4) = Zs_CMC_chm
 
     do iZ = 2, nZ_AMC_CMC_chm-1_ip
        do iS = 1, nS_AMC_CMC_chm
           if (S_AMC_CMC_chm(iS) <= S_threshold) then
              ! Direct interpolation
              call chm_mixture_fraction_interpolation_CMC(Xnormalized_prof_CMC_chm, 1_ip, &
                               Z_AMC_CMC_chm(iZ), aux(1))
           else
              ! Perform integration
              Zvar = S_AMC_CMC_chm(iS) * Z_AMC_CMC_chm(iZ) * (Zs_CMC_chm - Z_AMC_CMC_chm(iZ))
              call chm_find_PDF_parameters_CMC(Z_AMC_CMC_chm(iZ), Zvar, PDF_param)
              call chm_mixtFraction_integration_CMC(Xnormalized_prof_CMC_chm, 1_ip, &
                               PDF_param, aux(1))
           end if
           Xintegrated_table_AMC_CMC_chm(iZ,iS) = aux(1)
        end do
     end do

     ! Write file with integrated normalized profiles

     ! Get the names in a vector for the header
     names_rscal(1) = 'MF_AVG [-]:1'
     names_rscal(2) = 'S [-]:2'
     names_rscal(3) = 'X_INTEGRAL [-]:3'

     ! Write inert mixture file
     open(unit=1, file='Xintegrals.log', status="replace", action="write")
     write(unit=1,fmt="(300(a20))") names_rscal
     do iZ = 1, nZ_AMC_CMC_chm
        do iS = 1, nS_AMC_CMC_chm
           write(unit=1,fmt="(3(e12.6,8x))") Z_AMC_CMC_chm(iZ), S_AMC_CMC_chm(iS), Xintegrated_table_AMC_CMC_chm(iZ,iS)
        end do
     end do
     close(unit=1)

  end subroutine chm_AMC_integrals_CMC


  subroutine get_index_vector(bound_min, bound_max, n_vec, expon, val, index_v)

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/get_index_vector
     ! NAME 
     !    get_index_vector
     ! DESCRIPTION
     !    Given a value it finds the rounded down index for which the corresponding
     !    value in a vector defined by a potential distribution is closest.
     ! USED BY
     !    chm_calc_scalar_dissip_rate_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)          :: n_vec
     integer(ip), intent(out)         :: index_v
     real(rp), intent(in)             :: bound_min, bound_max, expon, val

     if (expon == 1.0_rp) then  ! Separate from others just to avoid numerical noise when doing 1.0_rp/1.0_rp
        index_v = floor( 1.0_rp + (n_vec-1) * ((val - bound_min) / (bound_max - bound_min)) )
     else
        index_v = floor( 1.0_rp + (n_vec-1) * ((val - bound_min) / (bound_max - bound_min))**(1.0_rp/expon) )
     end if

  end subroutine get_index_vector


  subroutine bilinear_intepolation(coordinates, heights, x, y, z)

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/bilinear_intepolation
     ! NAME 
     !    bilinear_intepolation
     ! DESCRIPTION
     !    Do a bilinear interpolation.
     ! USED BY
     !    chm_calc_scalar_dissip_rate_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     real(rp), intent(in)       :: coordinates(4), heights(4), x, y
     real(rp), intent(out)      :: z
     real(rp)                   :: z1, z2

     z1 = heights(1) + (heights(3) - heights(1)) * (x-coordinates(1)) / (coordinates(2) - coordinates(1))
     z2 = heights(2) + (heights(4) - heights(2)) * (x-coordinates(1)) / (coordinates(2) - coordinates(1))
     z  = z1 + (z2 -z1) * (y-coordinates(3)) / (coordinates(4) - coordinates(3))

  end subroutine bilinear_intepolation


  subroutine compute_Xnormalized_profile_CMC

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/compute_Xnormalized_profile_CMC
     ! NAME 
     !    compute_Xnormalized_profile_CMC
     ! DESCRIPTION
     !    Compute normalized scalar dissipation rate profile for AMC model from
     !    CMC model.
     ! USED BY
     !    chm_reaphy
     !***
     !-----------------------------------------------------------------------

     use def_chemic,         only:  Xnormalized_prof_CMC_chm, Z_CMC_chm, nZ_CMC_chm, &
                                    Zs_CMC_chm

     implicit none
     real(rp), parameter         :: small = 1.0e-8_rp

     integer(ip)                 :: iZ
     real(rp)                    :: aux(nZ_CMC_chm), aux2(nZ_CMC_chm)

     aux = 2.0_rp * Z_CMC_chm/Zs_CMC_chm - 1.0_rp
     aux2(1:nZ_CMC_chm) = 0.0_rp

     do iZ = 2,nZ_CMC_chm-1
        call compute_erfinv(aux(iZ),aux2(iZ))
        Xnormalized_prof_CMC_chm(iZ) = exp(-2.0_rp * aux2(iZ)**2.0_rp) * Zs_CMC_chm**2.0_rp
        if (Xnormalized_prof_CMC_chm(iZ) < small) then
           Xnormalized_prof_CMC_chm(iZ) = 0.0_rp
        end if
     end do

  end subroutine compute_Xnormalized_profile_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! M I X T.   F R A C T.  D I F F U S I O N !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_calc_diff_condVar_mixfraction_CMC

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/chm_calc_diff_condVar_mixfraction_CMC
     ! NAME 
     !    chm_calc_diff_condVar_mixfraction_CMC
     ! DESCRIPTION
     !    It computes the second derivative of the conditioned variables to be
     !    solved.
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !-----------------------------------------------------------------------
     use def_chemic,         only:  nZ_CMC_chm, diff_Z_CMC_chm, Yk_CMC_chm, &
                                    enthalp_CMC_chm, deriv2_Yk_CMC_chm, deriv2_enthalp_CMC_chm, &
                                    kfl_solve_enth_CMC_chm

     implicit none
     integer(ip)                 :: imixf, ipoin, iclas
     real(rp)                    :: numerator, denominator_inv(nZ_CMC_chm-2_ip), aux_var

     !
     ! This second derivative is of second order if the nodes are equally spaced and
     ! first order otherwise.
     !

     if( INOTMASTER ) then
        do imixf = 1, nZ_CMC_chm-2_ip
           aux_var = diff_Z_CMC_chm(imixf) * diff_Z_CMC_chm(imixf+1_ip) &
                   * (diff_Z_CMC_chm(imixf) + diff_Z_CMC_chm(imixf+1_ip))
           denominator_inv(imixf) = 1.0_rp / aux_var
        end do
  
        do ipoin = 1, npoin
           do imixf = 2, nZ_CMC_chm-1_ip
              do iclas = 1, nclas_chm
                 numerator = Yk_CMC_chm(imixf+1_ip,ipoin,iclas) * diff_Z_CMC_chm(imixf-1_ip) &
                            + Yk_CMC_chm(imixf-1_ip,ipoin,iclas) * diff_Z_CMC_chm(imixf) &
                            - Yk_CMC_chm(imixf,ipoin,iclas) * (diff_Z_CMC_chm(imixf-1)+diff_Z_CMC_chm(imixf))
                 deriv2_Yk_CMC_chm(imixf,ipoin,iclas) = 2.0_rp * numerator * denominator_inv(imixf-1_ip)
              end do
           end do
        end do

        if (kfl_solve_enth_CMC_chm /= 0_ip) then
           do ipoin = 1, npoin
              do imixf = 2, nZ_CMC_chm-1_ip
                 numerator = enthalp_CMC_chm(imixf+1_ip,ipoin) * diff_Z_CMC_chm(imixf-1_ip) &
                               + enthalp_CMC_chm(imixf-1_ip,ipoin) * diff_Z_CMC_chm(imixf) &
                               - enthalp_CMC_chm(imixf,ipoin) * (diff_Z_CMC_chm(imixf-1_ip)+diff_Z_CMC_chm(imixf))
                 deriv2_enthalp_CMC_chm(imixf,ipoin) = 2.0_rp * numerator * denominator_inv(imixf-1_ip)
              end do
           end do
        end if
     end if
  end subroutine chm_calc_diff_condVar_mixfraction_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! C O M P U T E   T I M E   S T E P !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_updtcc_CMC(dtmin)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_updtcc_CMC
     ! NAME 
     !    chm_updtcc_CMC
     ! DESCRIPTION
     !    This routine computes the critical time step for CMC model.
     ! USED BY
     !    chm_updtcc_CMC
     !***
     !-----------------------------------------------------------------------
     use def_parame
     use def_master
     use def_domain
     use def_chemic,         only: nZ_CMC_chm, kfl_advec_chm, kfl_ellen_chm, &
                                   condu_gp_CMC_chm, sphec_gp_CMC_chm, &
                                   spvol_gp_CMC_chm, &
                                   ADR_chm, Zs_CMC_chm, Z_CMC_chm, nvar_CMC_chm, &
                                   kfl_solve_enth_CMC_chm, Le_k
     use mod_ker_proper 
     use def_kermod
     use mod_ADR,            only : ADR_critical_time_step
     use mod_ADR,            only : mreac_adr
     use mod_ADR,            only : FROM_CRITICAL
     use mod_communications, only : PAR_MIN

     implicit none 
     real(rp),   intent(inout)    :: dtmin
     real(rp)                     :: dtcri(2)
     integer(ip)                  :: ielem,iclas,imixf,igaus                   ! Indices and dimensions
     integer(ip)                  :: pelty,pnode,dummi
     integer(ip)                  :: pgaus,plapl,porde,ptopo

     real(rp)                     :: PDF_val
     real(rp)                     :: dummr(mgaus*ndime), kcp_la
     real(rp)                     :: chale(3),chave(3),hleng(3),tragl(9)


     ! Matrices for elements (values at nodes)
     real(rp)                     :: elcod(ndime,mnode)                     ! Coordinates
     real(rp)                     :: elvel_CFD(ndime,mnode)                 ! Velocity from CFD
     real(rp)                     :: elZavg_CFD(mnode)                      ! Average mixture fraction from CFD
     real(rp)                     :: elZvar_CFD(mnode)                      ! Mixture fraction variance from CFD
     real(rp)                     :: elXtot_CFD(mnode)                      ! Total scalar dissipation rate from CFD
     real(rp)                     :: elZgrad_CFD(ndime,mnode)               ! Average mixture fraction gradient from CFD
     real(rp)                     :: elturb_dif_CFD(mnode)                  ! Mass turbulent diffusion coefficient from CFD
     real(rp)                     :: eldummy_CMC(mnode,nvar_CMC_chm)
     real(rp)                     :: eldumm2_CMC(mnode,nvar_CMC_chm)

     ! Matrices for elements (values at Gaussian points)
     real(rp)                     :: gpvol(mgaus)                           ! |J|*w
     real(rp)                     :: gprea(mgaus,mreac_adr)                 ! r
     real(rp)                     :: gpcar(ndime,mnode,mgaus)               ! dNk/dxj
     real(rp)                     :: gphes(ntens,mnode,mgaus)               ! dNk/dxidxj
     real(rp)                     :: gpdif(mgaus,nvar_CMC_chm)              ! D_k
     real(rp)                     :: gpvel_CFD(ndime,mgaus)                 ! Velocity
     real(rp)                     :: gpZavg_CFD(mgaus)                      ! Average mixture fraction
     real(rp)                     :: gpZvar_CFD(mgaus)                      ! Mixture fraction variance
     real(rp)                     :: gpXtot_CFD(mgaus)                      ! Total scalar dissipation rate
     real(rp)                     :: gpZgrad_CFD(ndime,mgaus)               ! Average mixture fraction gradient
     real(rp)                     :: gpturb_dif_CFD(mgaus)                  ! Mass turbulent diffusion coefficient
     real(rp)                     :: gplam_dif_CMC(mgaus,nvar_CMC_chm)      ! Mass laminar diffusion coefficient
     real(rp)                     :: gp_densi_PDF_CMC(mgaus)                ! Uncond. rho * probability density function
     real(rp)                     :: gp_veloc_CMC_chm(1:ndime,mgaus)        ! Conditional velocity
     real(rp)                     :: gp_diff_phys_spc(mgaus,nvar_CMC_chm)   ! Diffusion in physical space
     real(rp)                     :: gpdummy_CMC(mgaus,nvar_CMC_chm)
     real(rp)                     :: gpdumm2_CMC(mgaus,nvar_CMC_chm)
     real(rp)                     :: gpPDF_param(7)                         ! PDF parameters at Gaussian points
     real(rp)                     :: factor0_vel(ndime)                     ! Constant in the conditional velocity model
     real(rp)                     :: factor1_vel(ndime)                     ! Slope in the conditional velocity model


     if( INOTMASTER ) then
        mixt_fr: do imixf = 2, nZ_CMC_chm-1_ip
   
           call chm_global2local_CMC(imixf)

           call chm_cp_k_gauss_CMC

           elements: do ielem = 1,nelem
    
              !
              ! Element dimensions
              !
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              porde = lorde(pelty)
              ptopo = ltopo(pelty)

              !
              ! Initialization
              !

              elcod(1:ndime,1:mnode)                     = 0.0_rp
              elvel_CFD(1:ndime,1:mnode)                 = 0.0_rp
              elZavg_CFD(1:mnode)                        = 0.0_rp
              elZvar_CFD(1:mnode)                        = 0.0_rp
              elXtot_CFD(1:mnode)                        = 0.0_rp
              elZgrad_CFD(1:ndime,1:mnode)               = 0.0_rp
              elturb_dif_CFD(1:mnode)                    = 0.0_rp
              eldummy_CMC(1:mnode,1:nvar_CMC_chm)        = 0.0_rp
              eldumm2_CMC(1:mnode,1:nvar_CMC_chm)        = 0.0_rp
              gpvol(1:mgaus)                             = 0.0_rp
              gprea(1:mgaus,1:mreac_adr)                 = 0.0_rp
              gpcar(1:ndime,1:mnode,1:mgaus)             = 0.0_rp
              gphes(1:ntens,1:mnode,1:mgaus)             = 0.0_rp
              gpdif(1:mgaus,1:nvar_CMC_chm)              = 0.0_rp
              gp_diff_phys_spc(1:mgaus,1:nvar_CMC_chm)   = 0.0_rp
              gpdummy_CMC(1:mgaus,1:nvar_CMC_chm)        = 0.0_rp
              gpdumm2_CMC(1:mgaus,1:nvar_CMC_chm)        = 0.0_rp
              gpvel_CFD(1:ndime,1:mgaus)                 = 0.0_rp
              gpZavg_CFD(1:mgaus)                        = 0.0_rp
              gpZvar_CFD(1:mgaus)                        = 0.0_rp
              gpXtot_CFD(1:mgaus)                        = 0.0_rp
              gpZgrad_CFD(1:ndime,1:mgaus)               = 0.0_rp
              gplam_dif_CMC(1:mgaus,1:nvar_CMC_chm)      = 0.0_rp
              gpturb_dif_CFD(1:mgaus)                    = 0.0_rp
              gp_densi_PDF_CMC                           = 0.0_rp
              gp_veloc_CMC_chm(1:ndime,1:mgaus)          = 0.0_rp
              gpPDF_param(1:7)                           = 0.0_rp
              gpPDF_param(4)                             = Zs_CMC_chm
              factor0_vel(1:ndime)                       = 0.0_rp
              factor1_vel(1:ndime)                       = 0.0_rp

              !
              ! Gather values at the element
              !
               call chm_elmgac_CMC(&
                      pnode, lnods(1:pnode,ielem), elcod, elvel_CFD, elZavg_CFD, elZvar_CFD,&
                      elXtot_CFD, elZgrad_CFD, elturb_dif_CFD, eldummy_CMC, eldumm2_CMC, imixf)

              !
              ! CHALE, HLENG and TRAGL 
              !
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
                   
              call elmchl(&
                   tragl,hleng,elcod,dummr,chave,chale,pnode,porde,hnatu(pelty),&
                   kfl_advec_chm,kfl_ellen_chm)
      
              !
              ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
              !
              call elmcar(&
                   pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                   elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
                   gphes,ielem)

              !
              ! Send quantities to Gaussian points
              !
              call chm_elmpre_CMC(&
                      pnode, pgaus, elmar(pelty)%shape, elvel_CFD, elZavg_CFD, &
                      elZvar_CFD, elXtot_CFD, elZgrad_CFD, elturb_dif_CFD, eldummy_CMC, &
                      eldumm2_CMC, gpvel_CFD, gpZavg_CFD, gpZvar_CFD, gpXtot_CFD, gpZgrad_CFD, &
                      gpturb_dif_CFD, gpdummy_CMC, gpdumm2_CMC)

              ! Compute the laminar diffusion coefficient D at Gaussian points

              do igaus = 1, pgaus
                 kcp_la = condu_gp_CMC_chm(ielem)%a(igaus,1,1) * spvol_gp_CMC_chm(ielem)%a(imixf,igaus,1) / &
                                 sphec_gp_CMC_chm(ielem) % a(igaus,1,1)
                 do iclas = 1,nclas_chm
                    gpdif(igaus,iclas) = kcp_la / Le_k(iclas)
                 end do
                 if (kfl_solve_enth_CMC_chm /= 0_ip) then  ! In case enthalpy is transported
                    gpdif(igaus,nvar_CMC_chm) = kcp_la
                 end if
              end do

              gauss: do igaus = 1, pgaus
                 !
                 ! Find values at Gaussian points
                 !
                 call chm_find_PDF_parameters_CMC(gpZavg_CFD(igaus), gpZvar_CFD(igaus), gpPDF_param(1:7))

                 call chm_find_veloc_parameters_CMC(gpvel_CFD(1:ndime,igaus), gpZavg_CFD(igaus), &
                           gpZvar_CFD(igaus), gpZgrad_CFD(1:ndime,igaus), gpturb_dif_CFD(igaus), &
                           factor0_vel(1:ndime), factor1_vel(1:ndime))

                 !
                 ! Compute all the conditional values for the not solved variables
                 !

                 ! Find pdf value
                 PDF_val = gpPDF_param(3) * Z_CMC_chm(imixf)**(gpPDF_param(1)-1.0_rp) * &
                               (gpPDF_param(4) - Z_CMC_chm(imixf))**(gpPDF_param(2)-1.0_rp)

                 ! Find density times PDF
                 gp_densi_PDF_CMC(igaus) = densi_gp(ielem) % a(igaus,1,1) * PDF_val

                 ! Diffusion in physical space
                 !!!!! MODIFY? THE FOLLOWING CODE DEPENDS ON IF WE HAVE rho*D_turb or D_turb
                 do iclas = 1,nvar_CMC_chm
                    gp_diff_phys_spc(igaus,iclas) = gp_densi_PDF_CMC(igaus) * &
                                                     (gpdif(igaus,iclas) + gpturb_dif_CFD(igaus))
                 end do

                 ! Conditional velocity
                 gp_veloc_CMC_chm(1:ndime,igaus) = factor0_vel(1:ndime) + factor1_vel(1:ndime) * Z_CMC_chm(imixf)

              end do gauss

              do iclas = 1, nvar_CMC_chm
                 ! Compute time-step

                 call ADR_critical_time_step(ADR_chm(iclas),gp_densi_PDF_CMC(1:pgaus), &
                        gp_veloc_CMC_chm(1:ndime,1:pgaus),gp_diff_phys_spc(1:pgaus,iclas), &
                        gprea,dtcri,chale(1),chale(2))
                 ! Take minimum time step
                 dtmin = min(dtmin,dtcri(1))
              end do

           end do elements

        end do mixt_fr
   
     end if
     !
     ! Look for minimum over subdomains
     !
     call PAR_MIN(dtmin,'IN MY CODE')

  end subroutine chm_updtcc_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! C O M P U T E   V A R I A B L E S   O F   I N T E R E S T !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_calc_temp_CMC(imixf)
    !-----------------------------------------------------------------------
    !****f* chemic/chm_calc_temp_CMC
    ! NAME 
    !    chm_calc_temp_CMC
    ! DESCRIPTION
    !    It computes conditional values for temperature at nodes.
    ! USES
    ! USED BY
    !-----------------------------------------------------------------------
    use def_master,      only : therm, conce
    use def_chemic,      only : temp_CMC_chm

    implicit none

    integer(ip), intent(in)  :: imixf
    integer(ip)              :: ipoin
    real(rp)                 :: T, Tini

    if (INOTMASTER) then
       !
       ! Loop over points
       !

       points: do ipoin = 1, npoin
          ! Find conditional temperature
          Tini = temp_CMC_chm(imixf,ipoin)  ! Initial temperature from which start iteration
          T    = Tini                       ! Assignement to not leave an empty value
          call chm_calc_T_from_hY_CMC(conce(ipoin,1:nclas_chm,1), therm(ipoin,1), &
                 Tini, T)
          temp_CMC_chm(imixf,ipoin) = T
       end do points
    end if
  end subroutine chm_calc_temp_CMC



  subroutine chm_calc_T_from_hY_CMC(Yk,h,Tini,T)
    !-----------------------------------------------------------------------
    !****f* chemic/chm_calc_T_from_hY_CMC
    ! NAME 
    !    chm_calc_T_from_hY_CMC
    ! DESCRIPTION
    !    It computes temperature from enthalpy and species mass fractions.
    ! USES
    ! USED BY
    !-----------------------------------------------------------------------
    use def_chemic,      only : coeff_cp_k, W_k
#ifdef CANTERA
    use def_chemic,      only : gas_chm
#endif
    use mod_physics,     only : physics_H_2_TCp

    implicit none

    real(rp), intent(in)     :: Yk(nclas_chm), h, Tini
    real(rp), intent(out)    :: T
    integer(ip)              :: ivalu
    real(rp)                 :: cploc(6,2), aux_cp_lt(8), aux_cp_ht(8)
    real(rp)                 :: aux_cp


#ifdef CANTERA
    ! Find temperature from enthalpy and species mass fractions
    call cantera_alya_cp(nclas_chm,coeff_cp_k, Yk(1:nclas_chm), &
              W_k, aux_cp_lt(1:8), aux_cp_ht(1:8))
    do ivalu = 1, 6
       cploc(ivalu,1) = aux_cp_lt(ivalu)
       cploc(ivalu,2) = aux_cp_ht(ivalu)
    end do
    T = Tini
    call physics_H_2_TCp(h, cploc, T, aux_cp)
#endif
  end subroutine chm_calc_T_from_hY_CMC



  subroutine chm_calc_h_from_TY_CMC(Yk,T,h)
    !-----------------------------------------------------------------------
    !****f* chemic/chm_calc_h_from_TY_CMC
    ! NAME 
    !    chm_calc_h_from_TY_CMC
    ! DESCRIPTION
    !    It computes enthalpy from temperature and species mass fractions.
    ! USES
    ! USED BY
    !-----------------------------------------------------------------------
    use def_chemic,      only : coeff_cp_k, W_k
#ifdef CANTERA
    use def_chemic,      only : gas_chm
#endif
    use mod_physics,     only : physics_T_2_HCp

    implicit none

    real(rp), intent(in)     :: Yk(nclas_chm), T
    real(rp), intent(out)    :: h
    integer(ip)              :: ivalu
    real(rp)                 :: cploc(6,2), aux_cp_lt(8), aux_cp_ht(8)
    real(rp)                 :: aux_cp


#ifdef CANTERA
    ! Find enthalpy from temperature and species mass fractions
    call cantera_alya_cp(nclas_chm,coeff_cp_k, Yk(1:nclas_chm), &
              W_k, aux_cp_lt(1:8), aux_cp_ht(1:8))
    do ivalu = 1, 6
       cploc(ivalu,1) = aux_cp_lt(ivalu)
       cploc(ivalu,2) = aux_cp_ht(ivalu)
    end do
    call physics_T_2_HCp(T, cploc, h, aux_cp)
#endif
  end subroutine chm_calc_h_from_TY_CMC


  subroutine chm_calc_cp_k_W_nodes_CMC(imixf,ipoin,Yk,cp,condu,wmean)
    !-----------------------------------------------------------------------
    !****f* chemic/chm_calc_cp_k_W_nodes_CMC
    ! NAME 
    !    chm_calc_cp_k_W_nodes_CMC
    ! DESCRIPTION
    !    It computes conditional values for specific heat, conductivity and
    !    molecular weight at a given mixture fraction and physical point.
    ! USES
    ! USED BY
    !-----------------------------------------------------------------------
    use def_master,      only : prthe
    use def_kermod,      only : gasco
    use def_chemic,      only : coeff_cp_k, W_k, temp_CMC_chm
#ifdef CANTERA
    use def_chemic,      only : gas_chm
#endif
    use mod_physics,     only : physics_T_2_HCp

    implicit none

    integer(ip), intent(in)  :: imixf, ipoin
    real(rp),    intent(in)  :: Yk(nclas_chm)
    real(rp),    intent(out) :: cp,condu,wmean
    integer(ip)              :: ivalu
    real(rp)                 :: cploc(6,2), aux_cp_lt(8), aux_cp_ht(8)
    real(rp)                 :: aux_h, aux_T

#ifdef CANTERA
    ! Compute cp
    call cantera_alya_cp(nclas_chm,coeff_cp_k, Yk(1:nclas_chm), &
              W_k, aux_cp_lt(1:8), aux_cp_ht(1:8))
    do ivalu = 1, 6
       cploc(ivalu,1) = aux_cp_lt(ivalu)
       cploc(ivalu,2) = aux_cp_ht(ivalu)
    end do
    aux_T = temp_CMC_chm(imixf,ipoin)
    call physics_T_2_HCp(aux_T,cploc,aux_h,cp)

    ! Compute molecular weight
    call compute_molWeight_mixt_CMC_chm(Yk(1:nclas_chm),wmean)

    ! Compute conductivity
    call setState_TPX(gas_chm,temp_CMC_chm(imixf,ipoin),prthe(1),Yk(1:nclas_chm))
    call setMassFractions(gas_chm,Yk(1:nclas_chm))
    condu = thermalConductivity(gas_chm)
#endif

  end subroutine chm_calc_cp_k_W_nodes_CMC



  subroutine chm_calc_densi_visco_gauss_CMC(imixf)
    !-----------------------------------------------------------------------
    !****f* chemic/chm_calc_densi_visco_gauss_CMC
    ! NAME 
    !    chm_calc_densi_visco_gauss_CMC
    ! DESCRIPTION
    !    It computes conditional values for density and viscosity at Gauss 
    !    points at a given mixture fraction level.
    !    Variables transfered to the CFD (which are deemed critical) are com-
    !    puted at Gauss points for the sake of accuracy instead of nodes.
    ! USES
    ! USED BY
    !-----------------------------------------------------------------------
    use def_master,      only : prthe
    use def_kermod,      only : gasco
    use def_domain,      only : nelem,ltype,nnode,ltypb,nboun,lnodb,&
                                ngaus,llapl,lorde,ltopo,elmar,lnods                               
    use def_chemic,      only : visco_gp_CMC_chm, spvol_gp_CMC_chm
#ifdef CANTERA
    use def_chemic,      only : gas_chm
#endif

    implicit none

    integer(ip), intent(in) :: imixf
    integer(ip)             :: ielem,igaus,inode,iclas
    integer(ip)             :: pblty,pnodb
    integer(ip)             :: pelty,pnode,lnods_loc(mnode)
    integer(ip)             :: pgaus
    real(rp)                :: elcon(mnode,nclas_chm)
    real(rp)                :: elcod(ndime,mnode)
    real(rp)                :: elh(mnode)
    real(rp)                :: eltem(mnode)
    real(rp)                :: gpcon(mgaus,nclas_chm)
    real(rp)                :: gph(mgaus)
    real(rp)                :: gptem(mgaus)
    real(rp)                :: aux_wmean_mixf

    if (INOTMASTER) then
#ifdef CANTERA
       !
       ! Loop over elements to find visco_gp_CMC_chm y spvol_gp_CMC_chm
       !
       elements: do ielem = 1, nelem
          !
          ! Element dimensions
          !
          pelty = ltype(ielem)

          if( pelty > 0 ) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)

              !
              ! Gather all
              !
              lnods_loc(1:pnode) = lnods(1:pnode,ielem)
              call chm_gatherProp_CMC( &
                       pnode,lnods_loc,elcod,elcon(1:pnode,1:nclas_chm),elh,eltem)

              !
              ! Initialization variables
              !
              gpcon(1:pgaus,1:nclas_chm) = 0.0_rp
              gph(1:mgaus)               = 0.0_rp
              gptem(1:mgaus)             = 0.0_rp

              !
              ! Species mass fraction Y_k at Gauss points
              !
              do iclas = 1,nclas_chm
                 do igaus = 1,pgaus
                    do inode = 1,pnode
                       gpcon(igaus,iclas) = gpcon(igaus,iclas)&
                                            + elmar(pelty)%shape(inode,igaus) * elcon(inode,iclas)
                    end do
                 end do
              end do

              !
              ! Enthalpy and temperature at Gauss points
              !
              do igaus = 1,pgaus
                 do inode = 1,pnode
                    gph(igaus)   = gph(igaus) &
                                         + elmar(pelty)%shape(inode,igaus) * elh(inode)
                    gptem(igaus) = gptem(igaus) &
                                         + elmar(pelty)%shape(inode,igaus) * eltem(inode)
                 end do
              end do

              !
              ! Compute transport properties
              ! Cantera properties
              !
              gauss: do igaus = 1,pgaus
                 call setState_TPX(gas_chm,gptem(igaus),prthe(1),gpcon(igaus,1:nclas_chm))
                 call setMassFractions(gas_chm,gpcon(igaus,1:nclas_chm))

                 visco_gp_CMC_chm(ielem) % a(imixf,igaus,1) = viscosity(gas_chm)
                 
                 call compute_molWeight_mixt_CMC_chm(gpcon(igaus,1:nclas_chm),aux_wmean_mixf)
                 spvol_gp_CMC_chm(ielem) % a(imixf,igaus,1) = gptem(igaus) * gasco / &
                                                              aux_wmean_mixf / prthe(1)
              end do gauss
          end if
       end do elements
#endif

    end if

  end subroutine chm_calc_densi_visco_gauss_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! I N T E G R A T I O N S !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_integrate_flow_var_points_CMC

     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_integrate_flow_var_points_CMC
     ! NAME 
     !    chm_integrate_flow_var_points_CMC
     ! DESCRIPTION
     !    This routine integrates density, viscosity, heat release and if required
     !    mass fractions, enthalpy and temperature at nodes.
     ! USED BY
     !    
     !***
     !-----------------------------------------------------------------------

     use def_domain
     use def_master,              only :  postp, wmean, sphek, condk, ittim
     use def_chemic,              only :  hrr_CMC_chm, nZ_CMC_chm, temp_int_CMC_chm, &
                                          Yk_CMC_chm, enthalp_CMC_chm, temp_CMC_chm, &
                                          hrr_chm, Yk_int_CMC_chm, enthalp_int_CMC_chm, &
                                          src_Yk_CMC_chm, src_Yk_int_CMC_chm, &
                                          Zavg_CFD_chm, Zvar_CFD_chm, nvar_CMC_chm, &
                                          kfl_solve_enth_CMC_chm

     implicit none
     integer(ip)                       :: ipoin, ielem, imixf, iclas
     real(rp)                          :: cp_Z(nZ_CMC_chm)
     real(rp)                          :: condu_Z(nZ_CMC_chm)
     real(rp)                          :: wmean_Z(nZ_CMC_chm)
     real(rp), allocatable             :: aux_integrals_CMC_chm(:)
     real(rp), allocatable             :: aux_var(:,:)


     if( INOTMASTER ) then
        ! Compute heat release at nodes
        call chm_heatRelease_field_CMC

        print*, '--| ALYA     CHEMIC: VARIABLES INTEGRATION AT NODES'

        !!!!!!! ÉSTE ES EL CORRECTO
        !if( mod(ittim, postp(1) % npp_stepi(56) ) == 0 ) then   ! 56 is the position for conditional mass fractions

        !!!!!!! PROVISIONAL -> ELIMINAR
        if( mod(ittim, postp(1) % npp_stepi(57,0) ) == 0_ip ) then   ! 57 is for unconditional mass fractions
           ! When writing integrate all variables

           allocate(aux_integrals_CMC_chm(2_ip*nclas_chm+6_ip))
           allocate(aux_var(nZ_CMC_chm,2_ip*nclas_chm+1_ip))

           do ipoin = 1,npoin

              ! Compute cp, conductivity and molecular weight at nodes
              do imixf = 1, nZ_CMC_chm
                 call chm_calc_cp_k_W_nodes_CMC(imixf,ipoin, &
                      Yk_CMC_chm(imixf,ipoin,1:nclas_chm),cp_Z(imixf),&
                      condu_Z(imixf),wmean_Z(imixf))
              end do

              ! Integrate all the variables
              do iclas = 1, nclas_chm
                 aux_var(1:nZ_CMC_chm,iclas)           = Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
                 aux_var(1:nZ_CMC_chm,iclas+nclas_chm) = src_Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
              end do

              if (kfl_solve_enth_CMC_chm == 0) then
                 aux_var(1:nZ_CMC_chm,2_ip*nclas_chm+1_ip) = enthalp_CMC_chm(1:nZ_CMC_chm,1)
              else
                 aux_var(1:nZ_CMC_chm,2_ip*nclas_chm+1_ip) = enthalp_CMC_chm(1:nZ_CMC_chm,ipoin)
              end if

              call chm_mxt_fr_integr_previous_steps_CMC(2_ip*nclas_chm+6_ip, &
                     (/temp_CMC_chm(:,ipoin), cp_Z(:), condu_Z(:), wmean_Z(:), &
                     hrr_CMC_chm(:,ipoin), aux_var(:,:) /), &
                     Zavg_CFD_chm(ipoin), Zvar_CFD_chm(ipoin), &
                     aux_integrals_CMC_chm(:))

              temp_int_CMC_chm(ipoin)               = aux_integrals_CMC_chm(1)
              sphek(ipoin,1)                        = aux_integrals_CMC_chm(2)
              condk(ipoin,1)                        = aux_integrals_CMC_chm(3)
              wmean(ipoin,1)                        = aux_integrals_CMC_chm(4)
              hrr_chm(ipoin)                        = aux_integrals_CMC_chm(5)
              Yk_int_CMC_chm(ipoin,1:nclas_chm)     = aux_integrals_CMC_chm(6:nclas_chm+5)
              src_Yk_int_CMC_chm(ipoin,1:nclas_chm) = aux_integrals_CMC_chm(nclas_chm+6:2*nclas_chm+5)
              enthalp_int_CMC_chm(ipoin)            = aux_integrals_CMC_chm(2*nclas_chm+6)

           end do

           deallocate(aux_integrals_CMC_chm)
           deallocate(aux_var)
        else

           allocate(aux_integrals_CMC_chm(5))

           do ipoin = 1,npoin
              ! Compute cp, conductivity and molecular weight at nodes
              do imixf = 1, nZ_CMC_chm
                 call chm_calc_cp_k_W_nodes_CMC(imixf,ipoin, &
                      Yk_CMC_chm(imixf,ipoin,1:nclas_chm),cp_Z(imixf),&
                      condu_Z(imixf),wmean_Z(imixf))
              end do        

              ! Integrate temperature, specific heat, conductivity, molecular weight and heat release
              call chm_mxt_fr_integr_previous_steps_CMC(5_ip, &
                   (/temp_CMC_chm(:,ipoin), cp_Z(:), condu_Z(:), wmean_Z(:), &
                   hrr_CMC_chm(:,ipoin)/), &
                   Zavg_CFD_chm(ipoin), Zvar_CFD_chm(ipoin), &
                   aux_integrals_CMC_chm(:))

              temp_int_CMC_chm(ipoin) = aux_integrals_CMC_chm(1)
              sphek(ipoin,1)          = aux_integrals_CMC_chm(2)
              condk(ipoin,1)          = aux_integrals_CMC_chm(3)
              wmean(ipoin,1)          = aux_integrals_CMC_chm(4)
              hrr_chm(ipoin)          = aux_integrals_CMC_chm(5)

           end do

           deallocate(aux_integrals_CMC_chm)
        end if
     end if

  end subroutine chm_integrate_flow_var_points_CMC



  subroutine chm_integrate_flow_var_gauss_CMC

     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_integrate_flow_var_gauss_CMC
     ! NAME 
     !    chm_integrate_flow_var_gauss_CMC
     ! DESCRIPTION
     !    This routine integrates density, viscosity, heat release and if required
     !    mass fractions, enthalpy and temperature at Gauss points.
     ! USED BY
     !    
     !***
     !-----------------------------------------------------------------------

     use def_master,              only :  densi_gp, visco_gp
     use def_domain,              only :  nelem,ltype,nnode,ltypb,nboun,lnodb,&
                                          ngaus,llapl,lorde,ltopo,elmar,lnods
     use def_chemic,              only :  nZ_CMC_chm, W_k, spvol_gp_CMC_chm, &
                                          visco_gp_CMC_chm

     implicit none
     integer(ip)                       :: ielem, pelty, pnode, pgaus
     integer(ip)                       :: imixf, igaus, lnods_loc(mnode)
     real(rp)                          :: elvisco(mnode)
     real(rp)                          :: elspvol(mnode)
     real(rp)                          :: elZavg(mnode)
     real(rp)                          :: elZvar(mnode)
 
     real(rp)                          :: gpZavg(mgaus)
     real(rp)                          :: gpZvar(mgaus)
     real(rp)                          :: gpvisco(nZ_CMC_chm,mgaus)
     real(rp)                          :: gpspvol(nZ_CMC_chm,mgaus)
     real(rp)                          :: gp_aux(2)

     print*, '--| ALYA     CHEMIC: DENSITY AND VISCOSITY INTEGRATION AT GAUSS POINTS'

     if (INOTMASTER) then
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          !
          ! Element dimensions
          !
          pelty = ltype(ielem)

          if( pelty > 0 ) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)

              ! Initialization
              gpZavg(1:mgaus)               = 0.0_rp
              gpZvar(1:mgaus)               = 0.0_rp
              gpvisco(1:nZ_CMC_chm,1:mgaus) = 0.0_rp
              gpspvol(1:nZ_CMC_chm,1:mgaus) = 0.0_rp

              !
              ! Gather all
              !
              lnods_loc(1:pnode) = lnods(1:pnode,ielem)
              call chm_gatherIntegral_CMC_chm(pnode,lnods,elZavg,elZvar)

              !
              ! Send quantities to Gaussian points
              !
              call chm_elmpreIntegral_CMC(pnode,pgaus,elmar(pelty)%shape, &
                      elZavg,elZvar,gpZavg,gpZvar)

              do imixf = 1,nZ_CMC_chm
                 do igaus = 1,pgaus
                    gpspvol(imixf,igaus) = spvol_gp_CMC_chm(ielem) % a(imixf,igaus,1)
                    gpvisco(imixf,igaus) = visco_gp_CMC_chm(ielem) % a(imixf,igaus,1)
                 end do
              end do

              do igaus = 1,pgaus
                 ! Integrations in mixture fraction
                 call chm_mxt_fr_integr_previous_steps_CMC(2_ip,(/gpspvol(:,igaus), gpvisco(:,igaus)/), &
                     gpZavg(igaus),gpZvar(igaus), gp_aux)
                 densi_gp(ielem) % a(igaus,1,1) = 1.0_rp / gp_aux(1)
                 visco_gp(ielem) % a(igaus,1,1) = gp_aux(2)
              end do

          end if
       end do elements

     end if

  end subroutine chm_integrate_flow_var_gauss_CMC


  subroutine chm_gatherIntegral_CMC_chm(pnode,lnods,elZavg,elZvar)

     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_gatherIntegral_CMC_chm
     ! NAME 
     !    chm_gatherIntegral_CMC_chm
     ! DESCRIPTION
     !    Gather values for chm_integrate_flow_var_gauss_CMC.
     ! USED BY
     !    chm_integrate_flow_var_gauss_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,              only :  Zavg_CFD_chm, Zvar_CFD_chm

     implicit none
     integer(ip), intent(in)  :: pnode, lnods(pnode)
     real(rp),    intent(out) :: elZavg(pnode)
     real(rp),    intent(out) :: elZvar(pnode)

     integer(ip)              :: inode,ipoin

     !
     ! Initialization
     !
     elZavg(1:pnode) = 0.0_rp
     elZvar(1:pnode) = 0.0_rp

     do inode = 1,pnode
        ipoin = lnods(inode)
        elZavg(inode) = Zavg_CFD_chm(ipoin)
        elZvar(inode) = Zvar_CFD_chm(ipoin)
     end do
  end subroutine chm_gatherIntegral_CMC_chm


  subroutine chm_elmpreIntegral_CMC(pnode,pgaus,gpsha,elZavg,elZvar,gpZavg,gpZvar)

     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_elmpreIntegral_CMC
     ! NAME 
     !    chm_elmpreIntegral_CMC
     ! DESCRIPTION
     !    Transfer elemental values to Gaussian points for the set of not 
     !    solved variables that do not depend on the mixture fraction.
     ! USES
     ! USED BY
     !    chm_integrate_flow_var_gauss_CMC
     !***
     !------------------------------------------------------------------------ 

     use def_chemic,             only :  nZ_CMC_chm

     implicit none
     integer(ip), intent(in)  :: pnode, pgaus
     real(rp),    intent(in)  :: gpsha(pnode,pgaus)
     real(rp),    intent(in)  :: elZavg(pnode)
     real(rp),    intent(in)  :: elZvar(pnode)
     real(rp),    intent(out) :: gpZavg(pgaus)
     real(rp),    intent(out) :: gpZvar(pgaus)

     integer(ip)              :: igaus, inode


     do igaus = 1,pgaus
        do inode = 1,pnode
           gpZavg(igaus)  = gpZavg(igaus) + gpsha(inode,igaus) * elZavg(inode)
           gpZvar(igaus)  = gpZvar(igaus) + gpsha(inode,igaus) * elZvar(inode)
        end do
     end do

  end subroutine chm_elmpreIntegral_CMC


  subroutine chm_rho_visco_nodal_project_CMC
  !------------------------------------------------------------------------
  ! NAME 
  !****f* Chemic/mod_chm_operations_CMC/chm_rho_visco_nodal_project_CMC
  ! DESCRIPTION
  !    Projection of unconditional density and viscosity
  ! USES
  ! USED BY
  !    
  !***
  !------------------------------------------------------------------------

  use def_parame
  use def_elmtyp
  use def_master
  use def_domain,             only :  nelem,ltype,nnode,ngaus,lorde,elmar,lnods, &
                                      ntens
  use def_chemic,             only :  densi_int_CMC_chm, visco_lam_int_CMC_chm
  
  implicit none
  integer(ip)                      :: ielem,inode,jnode,ipoin,igaus
  integer(ip)                      :: pnode, pgaus, porde, aux_plapl,pelty
  integer(ip)                      :: lnods_loc(mnode)
  real(rp)                         :: elcod(ndime,mnode)
  real(rp)                         :: aux_gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
  real(rp)                         :: aux_gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj  
  real(rp)                         :: gpvol(mgaus)
  real(rp)                         :: fact_densi, fact_visco
  real(rp)                         :: eldensi(mnode)
  real(rp)                         :: elvisco(mnode)


  if (INOTMASTER) then
       !
       ! Loop over elements
       !
       elements: do ielem = 1, nelem
          !
          ! Element dimensions
          !
          pelty = ltype(ielem)

          if( pelty > 0 ) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              porde = lorde(pelty)

              !
              ! Initialization
              !
              gpvol     = 0.0_rp
              aux_gpcar = 0.0_rp
              aux_gphes = 0.0_rp
              eldensi(pnode) = 0.0_rp
              elvisco(pnode) = 0.0_rp

              !
              ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
              !
              call elmcar(&
                   pnode,pgaus,aux_plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                   elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,aux_gpcar,&
                   aux_gphes,ielem)

              if( porde == 1 ) then
                 !
                 ! Element assembly
                 !
                 do igaus = 1,pgaus
                    fact_densi = gpvol(igaus) * densi_gp(ielem) % a(igaus,1,1)
                    fact_visco = gpvol(igaus) * visco_gp(ielem) % a(igaus,1,1)
                    do inode = 1,pnode
                       eldensi(inode) = eldensi(inode) + &
                            elmar(pelty)%shape(inode,igaus) * fact_densi
                       elvisco(inode) = elvisco(inode) + &
                            elmar(pelty)%shape(inode,igaus) * fact_visco
                    end do
                 end do
                 !
                 ! Nodal projection
                 !
                 lnods_loc(1:pnode) = lnods(1:pnode,ielem)
                 do inode = 1,pnode
                    ipoin = lnods_loc(inode)
                    densi_int_CMC_chm(ipoin)     = densi_int_CMC_chm(ipoin) + eldensi(inode)
                    visco_lam_int_CMC_chm(ipoin) = visco_lam_int_CMC_chm(ipoin) + elvisco(inode)
                 end do
              else
                 call runend('CHEMIC OPERATIONS MOD. CMC: projections for density and viscosity not implemented for orders higher than 1')
              end if

          end if
       end do elements
  end if

  !if( porde == 1 ) then
  !else
  !   
  !   do inode=1,pnode
  !      do jnode=1,pnode
  !         elmat(inode,jnode)=0.0_rp
  !      end do
  !   end do

  !   do igaus=1,pgaus
  !      do inode=1,pnode
  !         fact=gpvol(igaus)*gpsha(inode,igaus)/(gpden(igaus)*dtinv)
  !         do jnode=1,pnode
  !            elmat(inode,jnode)=elmat(inode,jnode) +fact*gpsha(jnode,igaus)
  !         end do     
  !      end do
  !   end do

  !   trace  = 0.0_rp
  !   elmass = 0.0_rp
  !   do inode = 1,pnode                       
  !      trace = trace + elmat(inode,inode)
  !      do jnode = 1,pnode                       
  !         elmass = elmass + elmat(inode,jnode)
  !      end do
  !   end do

  !   !
  !   ! Nodal projection
  !   !
  !   do inode = 1,pnode
  !      ipoin = lnods_loc(inode)
  !      dt_rho_chm_loc(ipoin) = dt_rho_chm_loc(ipoin) + elmat(inode,inode)*(elmass/trace)
  !   end do
  !end if

  end subroutine chm_rho_visco_nodal_project_CMC


  subroutine chm_mxt_fr_integr_previous_steps_CMC(nfield, fields_Z, &
                            Zavg, Zvar, int_fields)
 
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_mxt_fr_integr_previous_steps_CMC
     ! NAME 
     !    chm_mxt_fr_integr_previous_steps_CMC
     ! DESCRIPTION
     !    This routine takes the fields in mixture fraction
     !    for a given domain, does some preliminary actions and performs the
     !    corresponding integrals for each field.
     !
     !    - nfield: number of fields to be integrated.
     !    - fields_Z(nZ_CMC_chm,nfield): matrix containing all the fields for each
     !    mixture fraction.
     !    - Zavg: field of averaged/filtered mixture fraction.
     !    - Zvar: field of averaged/filtered mixture fraction variance.
     !    - int_fields(nfield): integrated fields.
     ! USED BY
     !    chm_integrate_flow_var_points_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,          only :  nZ_CMC_chm, S_threshold, Zs_CMC_chm

     implicit none
     integer(ip), intent(in)       :: nfield
     real(rp),    intent(in)       :: fields_Z(nZ_CMC_chm, nfield), Zavg, Zvar
     real(rp),    intent(out)      :: int_fields(nfield)

     real(rp)                      :: PDF_param(7)
     real(rp)                      :: S

     PDF_param(4) = Zs_CMC_chm

     S = Zvar / (Zavg*(Zs_CMC_chm - Zavg))
     if (S <= S_threshold) then
        ! Direct interpolation
        call chm_mixture_fraction_interpolation_CMC(fields_Z(1:nZ_CMC_chm, 1:nfield), &
                nfield, Zavg, int_fields(1:nfield))
     else
        ! Perform integration
        call chm_find_PDF_parameters_CMC(Zavg, Zvar, PDF_param)
        call chm_mixtFraction_integration_CMC(fields_Z(1:nZ_CMC_chm, 1:nfield), &
                nfield, PDF_param, int_fields(1:nfield))
     end if

  end subroutine chm_mxt_fr_integr_previous_steps_CMC



  subroutine chm_mixtFraction_integration_CMC(fields_Z, nfield, PDF_param, int_fields)

     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_mixtFraction_integration_CMC
     ! NAME 
     !    chm_mixtFraction_integration_CMC
     ! DESCRIPTION
     !    This routine computes the integrals for the PDF/FPDF times each 
     !    variable profile along mixture fraction. It only computes such 
     !    integral for one pair of averaged/filtered mixture fraction and variance
     !    but for all the fields.
     !
     !    - fields_Z(nmixt,nfield): fields along mixture fraction (point fixed 
     !    in physical space).
     !    - nfield: number of fields to be integrated.
     !    - Z_CMC_chm: vector of mixture fractions.
     !    - nZ_CMC_chm: number of mixture fractions.
     !    - PDF_param: defining parameters for the PDF/FPDF.
     !    - p_cut: vector with the positions of Z_CMC_chm where divide the integral.
     !    - int_fields(nfield): integrated fields.
     !
     !    Actions:
     !    1. Evaluation of the regularized incomplete beta function for coefficients
     !    pairs (alfa, beta) and (alfa+1, beta).
     !    2. Compute the integral from the field values and the incomplete beta
     !    functions considering the constant and linear contributions from the fields
     !    (piecewise linear interpolation).
     !    The algorithm is not iterative.
     ! USED BY
     !    chm_mxt_fr_integr_previous_steps_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,          only :  nZ_CMC_chm, Z_CMC_chm

     implicit none
     integer(ip), intent(in)       :: nfield
     real(rp),    intent(in)       :: fields_Z(nZ_CMC_chm,nfield), PDF_param(7)
     real(rp),    intent(out)      :: int_fields(nfield)
     
     integer(ip)                   :: imixf, ifield
     real(rp)                      :: incompl_beta(nZ_CMC_chm,2), ratio_f, diffZ_inv
     real(rp)                      :: factor_0, factor_1
     
     int_fields(1:nfield) = 0.0_rp

     incompl_beta(1,1:2) = 0.0_rp       ! Since Z_CMC_chm(1) = 0
     incompl_beta(nZ_CMC_chm,1:2) = 1.0_rp  ! This is B(x,y)/B(x,y)=1
     do imixf=2, nZ_CMC_chm-1
        call incob(PDF_param(1), PDF_param(2), Z_CMC_chm(imixf)/PDF_param(4), PDF_param(6),incompl_beta(imixf,1))
        call incob(PDF_param(1)+1.0_rp, PDF_param(2), Z_CMC_chm(imixf)/PDF_param(4), PDF_param(7), incompl_beta(imixf,2))
     end do

     do imixf = 1, nZ_CMC_chm-1
        diffZ_inv = 1.0_rp / (Z_CMC_chm(imixf+1) - Z_CMC_chm(imixf))
        factor_0 = incompl_beta(imixf+1,1) - incompl_beta(imixf,1)
        factor_1 = PDF_param(5) * (incompl_beta(imixf+1,2) - incompl_beta(imixf,2))

        do ifield = 1, nfield
           ratio_f = (fields_Z(imixf+1,ifield) - fields_Z(imixf,ifield)) * diffZ_inv

           ! Contribution from the constant part
           int_fields(ifield) = int_fields(ifield) + (fields_Z(imixf,ifield) - ratio_f*Z_CMC_chm(imixf)) * factor_0

           ! Contribution from the linear part
           int_fields(ifield) = int_fields(ifield) + ratio_f * factor_1
        end do
     end do


  end subroutine chm_mixtFraction_integration_CMC



  subroutine chm_mixture_fraction_interpolation_CMC(fields_Z, nfield, Zavg, int_fields)

     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_mixture_fraction_interpolation_CMC
     ! NAME 
     !    chm_mixture_fraction_interpolation_CMC
     ! DESCRIPTION
     !    This routine interpolates the fields at a given mixture fraction value.
     !
     !    - fields_Z(nmixt,nfield): fields along mixture fraction (point fixed
     !    in physical space).
     !    - nfield: number of fields to be integrated.
     !    - Z_CMC_chm: vector of mixture fractions.
     !    - nZ_CMC_chm: number of mixture fractions.
     !    - Zavg: averaged/filtered mixture fraction.
     !    - int_fields(nfield): interpolated fields.
     ! USED BY
     !    chm_mxt_fr_integr_previous_steps_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,          only :  nZ_CMC_chm, Z_CMC_chm

     implicit none
     integer(ip), intent(in)       :: nfield
     real(rp),    intent(in)       :: fields_Z(nZ_CMC_chm,nfield), Zavg
     real(rp),    intent(out)      :: int_fields(nfield)

     integer(ip)                   :: pos_aux1, pos_aux2, ifield
     real(rp)                      :: fact, diff_Z_int, diff_Z


     call find_pos_min_diff(Z_CMC_chm, nZ_CMC_chm, Zavg, pos_aux1)
     if (Zavg == Z_CMC_chm(pos_aux1)) then
        int_fields(1:nfield) = fields_Z(pos_aux1, 1:nfield)
     else
        if(Zavg > Z_CMC_chm(pos_aux1)) then
           pos_aux2 = pos_aux1 + 1_ip
        else
           pos_aux2 = pos_aux1
           pos_aux1 = pos_aux1 - 1_ip
        end if
        diff_Z_int = Z_CMC_chm(pos_aux2) - Z_CMC_chm(pos_aux1)
        diff_Z = Zavg - Z_CMC_chm(pos_aux1)
        do ifield = 1, nfield
           fact = (fields_Z(pos_aux2, ifield) - fields_Z(pos_aux1, ifield)) / diff_Z_int
           int_fields(ifield) = fields_Z(pos_aux1, ifield) + diff_Z * fact
        end do
     end if

  end subroutine chm_mixture_fraction_interpolation_CMC



  subroutine chm_find_PDF_parameters_CMC(Zavg, Zvar_avg, PDF_param)

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/chm_find_PDF_parameters_CMC
     ! NAME 
     !    chm_find_PDF_parameters_CMC
     ! DESCRIPTION
     !    It computes the PDF parameters for a beta function: alpha, beta,
     !    Zs^(1-alfa-beta)/B(alfa,beta), Zs, Zs*B(alfa+1,beta)/B(alfa,beta),
     !    B(alfa,beta), B(alfa+1,beta)
     ! USED BY
     !    chm_mxt_fr_integr_previous_steps_CMC, chm_elmpre_pdf_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     real(rp), intent(in)          :: Zavg, Zvar_avg
     real(rp), intent(inout)       :: PDF_param(7)

     PDF_param(1) = ( Zavg*(PDF_param(4)-Zavg)/Zvar_avg - 1.0_rp ) * Zavg / PDF_param(4)
     PDF_param(2) = PDF_param(1) * (PDF_param(4)/Zavg - 1.0_rp)
     call beta_func(PDF_param(1), PDF_param(2), PDF_param(6))
     PDF_param(3) = PDF_param(4)**(1.0_rp - PDF_param(1) - PDF_param(2)) / PDF_param(6)
     call beta_func(PDF_param(1)+1.0_rp, PDF_param(2), PDF_param(7))
     PDF_param(5) = PDF_param(4) * PDF_param(7) / PDF_param(6)

  end subroutine chm_find_PDF_parameters_CMC



  subroutine find_pos_min_diff(vector, N, val, pos)

     !----------------------------------------------------------------------- 
     !****f* Chemic/find_pos_min_diff
     ! NAME 
     !    find_pos_min_diff
     ! DESCRIPTION
     !    Find the integer pos such that minimizes abs(vector(i)-val) i=1,...,N.
     !    It there are several indexes for which the difference is equal and
     !    minimum it takes the lowest index.
     ! USED BY
     !    chm_mxt_fr_integr_previous_steps_CMC, chm_mixture_fraction_interpolation_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)       :: N
     real(rp),    intent(in)       :: vector(N), val
     integer(ip), intent(out)      :: pos

     integer(ip)                   :: i
     real(rp)                      :: v_aux(N)

     v_aux(1:N) = abs(vector(1:N) - val)
     pos = 1_ip
     do i = 2, N
        if ( v_aux(i) < v_aux(pos) )  pos = i
     end do

  end subroutine find_pos_min_diff


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! V A R I A B L E S   O F   I N T E R E S T !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_cp_k_gauss_CMC
    !-----------------------------------------------------------------------
    !****f* chemic/getProp_CMC
    ! NAME 
    !    getProp_CMC
    ! DESCRIPTION
    !    Compute properties for CMC combustion model
    ! USES
    ! USED BY
    !
    !-----------------------------------------------------------------------
    use def_master,       only : prthe,tempe, conce,sphec, sphec_gp
    use def_kermod,       only : gasco
                                   
    use def_domain,       only : nelem,ltype,nnode,ltypb,nboun,lnodb,&
                                 ngaus,llapl,lorde,ltopo,elmar,lnods
    use def_chemic,       only : coeff_cp_k, W_k, condu_gp_CMC_chm, &
                                 sphec_gp_CMC_chm

#ifdef CANTERA
    use def_chemic,       only : gas_chm
#endif
      use mod_physics,    only : physics_T_2_HCp

    implicit none
    real(rp), parameter      :: Tmin = 200.0_rp, Tmax = 3000.0_rp
    integer(ip)              :: iboun,pblty,inodb,pnodb,ipoin
    integer(ip)              :: ielem,igaus,inode,iclas,ivalu
    integer(ip)              :: pelty,pnode,lnods_loc(mnode)
    integer(ip)              :: pgaus,plapl,porde,ptopo
    real(rp)                 :: elcon(mnode,nclas_chm)
    real(rp)                 :: elcod(ndime,mnode)
    real(rp)                 :: elh(mnode)
    real(rp)                 :: eltem(mnode)
    real(rp)                 :: gpcon(mgaus,nclas_chm)
    real(rp)                 :: gph(mgaus)
    real(rp)                 :: gptem(mgaus)
    real(rp)                 :: dummr
    real(rp)                 :: cploc(6,2)
    real(rp)                 :: sphec_gp_lt_aux(6)
    real(rp)                 :: sphec_gp_ht_aux(6)
    real(rp)                 :: aux_h

    if (INOTMASTER) then
       !
       ! Loop over elements
       !

       elements: do ielem = 1,nelem
          !
          ! Element dimensions
          !
          pelty = ltype(ielem)

          if( pelty > 0 ) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              plapl = llapl(pelty) 
              porde = lorde(pelty)
              ptopo = ltopo(pelty)
              
              !
              ! Gather all
              !
              lnods_loc(1:pnode) = lnods(1:pnode,ielem)
              call chm_gatherProp_CMC( &
                       pnode,lnods_loc,elcod,elcon(1:pnode,1:nclas_chm),elh,eltem)

              !
              ! Initialization variables
              !
              gpcon                          = 0.0_rp 
              gph                            = 0.0_rp 
              gptem                          = 0.0_rp 
              condu_gp_CMC_chm(ielem) % a    = 0.0_rp
              sphec_gp_CMC_chm(ielem) % a    = 0.0_rp
 
              !
              ! Species mass fraction Y_k at Gauss points
              !
              do iclas = 1,nclas_chm
                 do igaus = 1,pgaus
                    do inode = 1,pnode
                       gpcon(igaus,iclas) = gpcon(igaus,iclas)&
                                            + elmar(pelty)%shape(inode,igaus) * elcon(inode,iclas)
                    end do
                 end do
              end do 

              !
              ! Enthalpy and temperature at Gauss points
              !
              do igaus = 1,pgaus
                 do inode = 1,pnode
                    gph(igaus)   = gph(igaus) &
                                         + elmar(pelty)%shape(inode,igaus) * elh(inode)
                    gptem(igaus) = gptem(igaus) &
                                         + elmar(pelty)%shape(inode,igaus) * eltem(inode)
                 end do
              end do
              
              !
              ! Compute transport properties
              ! Cantera properties
              !
#ifdef CANTERA
              do igaus = 1,pgaus
                 call setState_TPX(gas_chm,gptem(igaus),prthe(1),gpcon(igaus,1:nclas_chm))
                 call setMassFractions(gas_chm,gpcon(igaus,1:nclas_chm))
                 
                 call cantera_alya_cp(nclas_chm,coeff_cp_k,gpcon(igaus,1:nclas_chm), &
                                   W_k,sphec_gp_lt_aux, sphec_gp_ht_aux)

                 cploc(1:6,1) = sphec_gp_lt_aux(1:6)
                 cploc(1:6,2) = sphec_gp_ht_aux(1:6)
                 call physics_T_2_HCp(gptem(igaus), cploc, aux_h, &
                       sphec_gp_CMC_chm(ielem) % a(igaus,1,1))

                 condu_gp_CMC_chm(ielem) % a(igaus,1,1)   = thermalConductivity(gas_chm) 

              end do
#endif
          end if
       end do elements 

       !
       ! Cp coefficients on boundary for BC's
       !

#ifdef CANTERA
       boundaries: do iboun = 1,nboun
          pblty = ltypb(iboun)
          pnodb = nnode(pblty)
          do inodb = 1,pnodb
             ipoin = lnodb(inodb,iboun)


             call setState_TPX(gas_chm, max(Tmin,min(Tmax,tempe(ipoin,1))), &
                               prthe(1), conce(ipoin,1:nclas_chm,1))
             call setMassFractions(gas_chm,conce(ipoin,1:nclas_chm,1))
             
             call cantera_alya_cp(nclas_chm,coeff_cp_k,conce(ipoin,1:nclas_chm,1), &
                                      W_k,sphec(ipoin,1:6,1),sphec(ipoin,1:6,2))

          end do
       end do boundaries
#endif

    end if

  end subroutine chm_cp_k_gauss_CMC


  subroutine chm_gatherProp_CMC( &
                 pnode,lnods,elcod,elcon,elh,eltem)
     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_element_operations/chm_gatherProp_CMC
     ! NAME 
     !    chm_gatherProp_CMC
     ! DESCRIPTION
     !    Gather operations for the combustion models
     ! USES
     ! USED BY
     !    chm_cp_k_gauss_CMC
     !***
     !------------------------------------------------------------------------ 
     use def_master, only     :  conce,therm,tempe

     implicit none 
     real(rp), parameter      :: Tmin = 200.0_rp, Tmax = 3000.0_rp
     integer(ip), intent(in)  :: pnode
     integer(ip), intent(in)  :: lnods(pnode)
     real(rp),    intent(out) :: elcod(ndime,pnode)
     real(rp),    intent(out) :: elcon(pnode,nclas_chm)
     real(rp),    intent(out) :: elh(pnode)
     real(rp),    intent(out) :: eltem(pnode)

     integer(ip)              :: inode,ipoin,iclas,idime

     !
     ! Initialization
     !
     elh(1:pnode)               = 0.0_rp
     eltem(1:pnode)             = 0.0_rp
     elcod(1:ndime,1:pnode)     = 0.0_rp
     elcon(1:pnode,1:nclas_chm) = 0.0_rp

     !
     ! Concentration and coordinates
     !
     do inode=1,pnode
        ipoin=lnods(inode)
        do iclas=1,nclas_chm
           elcon(inode,iclas) = conce(ipoin,iclas,1)
        end do

        do idime=1,ndime
           elcod(idime,inode)   = coord(idime,ipoin)
        end do

        elh(inode)   = therm(ipoin,1)
        eltem(inode) = max(Tmin, min(Tmax, tempe(ipoin,1)))
     end do

  end subroutine chm_gatherProp_CMC


  subroutine compute_molWeight_mixt_CMC_chm(Yk,wmean_mixf)
     !----------------------------------------------------------------------- 
     !****f* Chemic/compute_molWeight_mixt_CMC_chm
     ! NAME                                                                   
     !    compute_molWeight_mixt_CMC_chm
     ! DESCRIPTION                                                            
     !    This routine computes the molecular weight for all the mixture
     !    fraction levels for a given physical point.
     ! USED BY                                                                
     !    
     !***                                                                     
     !-----------------------------------------------------------------------

     use def_chemic,             only :  nclas_chm, W_k

     implicit none
     real(rp), intent(in)             :: Yk(nclas_chm)
     real(rp), intent(out)            :: wmean_mixf
     integer(ip)                      :: iclas

     wmean_mixf = 0.0_rp
     do iclas = 1, nclas_chm
        wmean_mixf =  wmean_mixf + Yk(iclas) / W_k(iclas)
     end do
     wmean_mixf = 1.0_rp / wmean_mixf

  end subroutine compute_molWeight_mixt_CMC_chm



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! C H E M I C A L   C A L C U L A T I O N S !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_integrate_chem_source_CMC(dt)

     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_integrate_chem_source_CMC
     ! NAME 
     !    chm_integrate_chem_source_CMC
     ! DESCRIPTION
     !    Compute chemical source terms for CMC transport equations.
     ! USED BY
     !    
     !***
     !-----------------------------------------------------------------------

     use def_master,          only :  prthe
     use def_kermod,          only :  gasco
     use def_chemic,          only :  Yk_CMC_chm, src_Yk_CMC_chm, temp_CMC_chm, &
                                      nZ_CMC_chm
     use mod_physics,         only :  physics_T_2_HCp  

     implicit none
     real(rp), parameter           :: temp_min = 500.0_rp
     real(rp), intent(in)          :: dt
     integer(ip)                   :: imixf, ipoin
     real(rp)                      :: T_next, Yk_next(nclas_chm), dt_inv
     real(rp)                      :: Yk_aux(nclas_chm)
     real(rp)                      :: wmean_mixf
     real(rp)                      :: densi_t0, densi_next


     if( INOTMASTER ) then
        print*, '--| ALYA     CHEMIC: COMPUTING CHEMICAL SOURCE TERMS...'

        dt_inv = 1.0_rp / dt
#ifdef CANTERA
        do imixf = 2, nZ_CMC_chm-1
           !!!!print*, '--| ALYA     CHEMIC: FOR MIXTURE FRACTION LEVEL ', imixf
           do ipoin = 1, npoin

              if (temp_CMC_chm(imixf,ipoin) > temp_min)  then
                 T_next               = temp_CMC_chm(imixf,ipoin)
                 Yk_next(1:nclas_chm) = Yk_CMC_chm(imixf,ipoin,1:nclas_chm)
                 call cantera_integrate(T_next, prthe(1), Yk_next(1:nclas_chm), dt)

                 ! Compute density
                 Yk_aux(1:nclas_chm) = Yk_CMC_chm(imixf,ipoin,1:nclas_chm)
                 call compute_molWeight_mixt_CMC_chm(Yk_aux(:),wmean_mixf)
                 densi_t0 = prthe(1) * wmean_mixf / (temp_CMC_chm(imixf,ipoin) * gasco)
                 call compute_molWeight_mixt_CMC_chm(Yk_next(:),wmean_mixf)
                 densi_next = prthe(1) * wmean_mixf / (T_next * gasco)

                 ! Compute chemical source terms and update mass fractions
                 src_Yk_CMC_chm(imixf,ipoin,1:nclas_chm) = &
                    (Yk_next(1:nclas_chm) - Yk_CMC_chm(imixf,ipoin,1:nclas_chm)) * &
                    dt_inv * 0.5_rp * (densi_t0 + densi_next)
                 temp_CMC_chm(imixf,ipoin)           = T_next
                 Yk_CMC_chm(imixf,ipoin,1:nclas_chm) = Yk_next(1:nclas_chm)
              else
                 src_Yk_CMC_chm(imixf,ipoin,1:nclas_chm) = 0.0_rp
              end if

           end do

        end do
#endif
     end if

  end subroutine chm_integrate_chem_source_CMC



  subroutine chm_heatRelease_field_CMC

     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_integrate_chem_source_CMC
     ! NAME 
     !    chm_integrate_chem_source_CMC
     ! DESCRIPTION
     !    Compute chemical source terms for CMC transport equations.
     !    heat_release = sum_j=1^N(par_mh_j * omega^molar_j)
     !    omega^molar_j = omega^mass_j / W_j = dY_j/dt * rho * volume. Substituting:
     !    heat_release = density * volume * sum_j=1^N(par_mh_j * dY_j/dt / W_j)
     !    However, the heat release per volume unit is obtained at each cell so volume is omitted.
     ! USED BY
     !    
     !***
     !-----------------------------------------------------------------------

     use def_master,                 only : prthe
     use def_chemic,                 only : src_Yk_CMC_chm, Yk_CMC_chm, temp_CMC_chm, &
                                            hrr_CMC_chm, nZ_CMC_chm, W_k

     implicit none
     integer(ip)                         :: iclas, ipoin, iZ
     real(rp)                            :: par_mh(nclas_chm) ! Partial molar enthalpy


     if (INOTMASTER) then

        hrr_CMC_chm(1:nZ_CMC_chm,1:npoin) = 0.0_rp
#ifdef CANTERA
        do iZ = 1, nZ_CMC_chm
           do ipoin= 1, npoin
              ! Get molar enthalpies
              call cantera_partial_molar_enthalpies(nclas_chm,temp_CMC_chm(iZ,ipoin),prthe(1), &
                        Yk_CMC_chm(iZ,ipoin,1:nclas_chm),par_mh(1:nclas_chm))
              do iclas= 1, nclas_chm
                 hrr_CMC_chm(iZ,ipoin) = hrr_CMC_chm(iZ,ipoin) - &
                      par_mh(iclas) * src_Yk_CMC_chm(iZ,ipoin,iclas) / W_k(iclas)
              end do
           end do
        end do
#endif
     end if

  end subroutine chm_heatRelease_field_CMC



  subroutine chm_heatRelease_integral_CMC

     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_integrate_chem_source_CMC
     ! NAME 
     !    chm_integrate_chem_source_CMC
     ! DESCRIPTION
     !    Compute chemical source terms for CMC transport equations.
     !    heat_release = sum_j=1^N(par_mh_j * omega^molar_j)
     !    omega^molar_j = omega^mass_j / W_j = dY_j/dt * rho * volume. Substituting:
     !    heat_release = density * volume * sum_j=1^N(par_mh_j * dY_j/dt / W_j)
     !    However, the heat release per volume unit is obtained at each cell so volume is omitted.
     ! USED BY
     !    
     !***
     !-----------------------------------------------------------------------

     use def_domain
     use def_master
     use def_chemic,                only :  hrr_chm, hrr_int_chm

     implicit none
     integer(ip)                         :: ipoin, ielem, igaus, inode
     integer(ip)                         :: pnode, pgaus, pelty
     real(rp)                            :: gpvol, gpdet, gphre
     real(rp)                            :: elhre(mnode), elcod(ndime,mnode)
     real(rp)                            :: gpcar(ndime,mnode,mgaus)
     real(rp)                            :: xjaci(ndime,ndime),xjacm(ndime,ndime)


     if (INOTMASTER) then
        !
        ! Compute heat release integrated over the whole domain
        ! Transfer heat release field from nodes to Gaussian points and do the volumetric integral
        hrr_int_chm = 0.0_rp
        elements: do ielem = 1, nelem
           pelty = ltype(ielem)
           if (pelty > 0) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)

              !
              ! Gather operations
              ! 
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 elhre(inode) = hrr_chm(ipoin)
                 elcod(1:ndime,inode) = coord(1:ndime,ipoin)
              end do

              !
              ! 1st and 2nd order Cartesian derivatives, and dV:=GPVOL=|J|*wg
              !
              do igaus = 1, pgaus
                 call elmder(&
                       pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                       elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
                  gpvol = elmar(pelty)%weigp(igaus) * gpdet               ! |J|*wg
                  gphre = 0.0_rp
                  do inode = 1, pnode
                     gphre = gphre + elmar(pelty) % shape(inode,igaus) * elhre(inode)
                  end do
                  hrr_int_chm = hrr_int_chm + gphre * gpvol
              end do
           end if
        end do elements

     else
        hrr_int_chm = 0.0_rp
     end if

    call pararr('SUM',0_ip,1_ip,hrr_int_chm)

  end subroutine chm_heatRelease_integral_CMC


  subroutine chm_limit_Yk_CMC
     !----------------------------------------------------------------------- 
     !****f* Chemic/chm_limit_Yk_CMC
     ! NAME 
     !    chm_limit_Yk_CMC
     ! DESCRIPTION
     !    Limit the values of mass fractions.
     ! USED BY
     !    
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only :  nZ_CMC_chm, nclas_chm, Yk_CMC_chm

     implicit none
     integer(ip)                      :: imixf, ipoin, iclas
     real(rp)                         :: alpha, inv_alpha

     if (INOTMASTER) then
        do imixf = 1, nZ_CMC_chm
           do ipoin = 1, npoin
              alpha = 0.0_rp
              do iclas = 1, nclas_chm
                 Yk_CMC_chm(imixf,ipoin,iclas) = max(0.0_rp, min(Yk_CMC_chm(imixf,ipoin,iclas),1.0_rp))
                 alpha = alpha + Yk_CMC_chm(imixf,ipoin,iclas)
              end do
              inv_alpha = 1.0_rp / alpha
              Yk_CMC_chm(imixf,ipoin,1:nclas_chm) = inv_alpha * Yk_CMC_chm(imixf,ipoin,1:nclas_chm)
           end do
        end do
     end if

  end subroutine chm_limit_Yk_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! M A T H E M A T I C A L   F U N C T I O N S !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!! SHOULDN'T THESE FUNCTIONS BE IN THE KERNEL?


  subroutine incob ( a, b, x, bt, bix )

     !*****************************************************************************
     !
     !  INCOB computes the regularized incomplete beta function Ix(a,b). Hence, if we
     !  denote the incomplete beta function as
     !
     !  B(x; a, b) = int_0^x (t^(a-1) * (1-t)^(b-1) * dt)
     !
     !  and the beta function as
     !
     !  B(a,b) = int_0^1 (t^(a-1) * (1-t)^(b-1) * dt)
     !
     !  then the regularized incomplete beta function Ix(a,b) is
     !
     !  Ix(a,b) = B(x; a, b) / B(a,b)
     !
     !
     !  Licensing:
     !
     !    This routine is copyrighted by Shanjie Zhang and Jianming Jin. However, 
     !    they give permission to incorporate this routine into a user program 
     !    provided that the copyright is acknowledged.
     !
     !  Modified:
     !
     !    22 July 2012
     !
     !  Author:
     !
     !    Shanjie Zhang, Jianming Jin.
     !
     !  Reference:
     !
     !    Shanjie Zhang, Jianming Jin,
     !    Computation of Special Functions,
     !    Wiley, 1996,
     !    ISBN: 0-471-11963-6,
     !    LC: QA351.C45.
     !
     !  Parameters:
     !    a, b: parameters of the beta function
     !    bt: value of B(a,b)
     !    x: value for which evaluate the incomplete beta function
     !    bix: value of the incomplete beta function

     implicit none
     real(rp), intent(in)     :: a, b, x, bt
     real(rp), intent(out)    :: bix

     integer(ip)              :: k
     real(rp)                 :: dk(51), fk(51), s0
     real(rp)                 :: t1, t2, ta, tb


     ! ¿CAMBIAR LA PRECISIÓN DE LOS DOUBLE?

     s0 = ( a + 1.0D+00 ) / ( a + b + 2.0D+00 )

     if ( x <= s0 ) then

       do k = 1, 20
         dk(2*k) = k * ( b - k ) * x / &
           ( a + 2.0D+00 * k - 1.0D+00 ) / ( a + 2.0D+00 * k )
       end do

       do k = 0, 20
         dk(2*k+1) = - ( a + k ) * ( a + b + k ) * x &
           / ( a + 2.0D+00 * k ) / ( a + 2.0D+00 * k + 1.0D+00 )
       end do

       t1 = 0.0D+00
       do k = 20, 1, -1
         t1 = dk(k) / ( 1.0D+00 + t1 )
       end do
       ta = 1.0D+00 / ( 1.0D+00 + t1 )
       bix = x ** a * ( 1.0D+00 - x ) ** b / ( a * bt ) * ta

     else

       do k = 1, 20
         fk(2*k) = k * ( a - k ) * ( 1.0D+00 - x ) &
           / ( b + 2.0D+00 * k - 1.0D+00 ) / ( b + 2.0D+00 * k )
       end do

       do k = 0,20
         fk(2*k+1) = - ( b + k ) * ( a + b + k ) * ( 1.0D+00 - x ) &
           / ( b + 2.0D+00 * k ) / ( b + 2.0D+00 * k + 1.0D+00 )
       end do

       t2 = 0.0D+00
       do k = 20, 1, -1
         t2 = fk(k) / ( 1.0D+00 + t2 )
       end do
       tb = 1.0D+00 / ( 1.0D+00 + t2 )
       bix = 1.0D+00 - x ** a * ( 1.0D+00 - x ) ** b / ( b * bt ) * tb

     end if

     return

  end subroutine incob



  subroutine beta_func(x, y, beta_f)

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/beta_func
     ! NAME 
     !    beta_func
     ! DESCRIPTION
     !    It computes beta function. All the values are assumed to be strictly
     !    positive.
     ! USED BY
     !    chm_find_PDF_parameters_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     real(rp), intent(in)            :: x, y
     real(rp), intent(out)           :: beta_f

     beta_f = gamma(x) * gamma(y) / gamma(x+y)

  end subroutine beta_func


  
  subroutine compute_erfinv(y0,x)

     !----------------------------------------------------------------------- 
     !****f* Chemic/mod_chm_operations_CMC/compute_erfinv
     ! NAME 
     !    compute_erfinv
     ! DESCRIPTION
     !    Compute inverse error function by a bisection method.
     ! USED BY
     !    compute_Xnormalized_profile_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     real(rp), parameter             :: small=1.0e-8_rp, big=1.0e20_rp

     real(rp), intent(in)            :: y0
     real(rp), intent(out)           :: x

     real(rp)                        :: y, xmiddle, x_aux, residual1, residual2, xmin, xmax
     logical                         :: negative, convergence
  
     xmin        = 0.0_rp
     xmax        = 4.0_rp
     negative    = .false.
     convergence = .false.  

     if (y0 < 0.0_rp) then
        negative = .true.
        y = -y0
     else
        y = y0
     end if
  
     if (y >= erf(xmax)) then
        x = big
     else if (y==0.0_rp) then
        x = 0.0_rp
     else
        do while (.not. convergence)
           xmiddle = (xmin + xmax) / 2.0_rp
           residual1 =  y - erf(xmiddle)
  
           if (residual1 == 0.0_rp) then
              x = xmiddle
              convergence = .true.
           else
              if (residual1 < 0.0_rp) then
                 xmax = xmiddle
                 residual2 = y-erf(xmin)
                 x_aux = xmin
              else
                 xmin = xmiddle
                 residual2 = y-erf(xmax)
                 x_aux = xmax
              end if
  
              if (abs(residual1) < small) then
                 x = xmiddle
                 convergence = .true.
              else if (abs(residual2) < small) then
                 x = x_aux
                 convergence = .true.
              end if
           end if
  
        end do
     end if
  
     if (negative) then
        x = -x
     end if

  end subroutine compute_erfinv


end module mod_chm_operations_CMC
