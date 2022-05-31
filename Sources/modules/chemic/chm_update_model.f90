subroutine chm_update_model()
  
  use def_master
  use def_chemic
  use mod_ker_proper
  use mod_chm_finiteRate,         only : chm_getProp_finiteRate
  use mod_chm_operations_CMC,     only : chm_save_unconditional_fields_bc_CMC, &
                                         chm_bc_type_species_CMC, &
                                         chm_compute_initial_fields_CMC, &
                                         chm_initialization_domain_CMC

  implicit none
  integer(ip) :: ii
  
  if (kfl_model_chm == 1_ip) then
     
     !================!
     ! FLAMELET MODEL !
     !================!

     !
     ! Scaling reaction progress Yc from c for initialization
     !
     if ( kfl_premix_chm == 0_ip .and. kfl_fields_scale_chm /= 0_ip ) call chm_scale_Yc()

     !
     ! Read flamelet table for gas phase
     !
     if (kfl_spray_chm == 0 .or. ( kfl_spray_chm /= 0 .and. kfl_premix_chm == 0)) then
        if (kfl_lookg_chm > 0) then
           call chm_gp_reatab()
        else
           call chm_reatab()
        endif

        !
        ! Update a few times to get the scalar dissipation right
        !
        if (kfl_ufpv_chm > 0) then
           do ii=1,4
              call ker_updpro() 
              call chm_post_scalar_dissipation_rate(47_ip)
              call chm_post_scalar_dissipation_rate(49_ip)
              if ( kfl_fields_scale_chm == 1_ip ) call chm_scale_Yc()
              if (kfl_lookg_chm > 0) then
                 call chm_gp_reatab()
              else
                 call chm_reatab()
              endif
           enddo
        endif
     end if

  elseif (kfl_model_chm == 3_ip) then
     
     !=============================!
     ! FINITE RATE CHEMISTRY MODEL !
     !=============================!

#ifdef CANTERA
     !
     ! NASA polinomial coefficients and molecular weights for individual species
     !
     call cantera_initialization(1_ip,mechanism_path,nsize_mech_name)
     call cantera_coeff(mechanism_path,nsize_mech_name,coeff_cp_k, W_k)
     if (kfl_pfa_chm == 1_ip) then 
        call cantera_reduced(Red_spec)
     end if
#endif
     call chm_getProp_finiteRate()


  elseif (kfl_model_chm == 4_ip) then
     if( INOTMASTER ) then
        !==================================!
        ! CONDITIONAL MOMENT CLOSURE MODEL !
        !==================================!

        !!!!!!!!! FALTA LA PARTE DEL RESTART

#ifdef CANTERA
        !
        ! NASA polinomial coefficients and molecular weights for individual species
        !
        call cantera_initialization(1_ip,mechanism_path,nsize_mech_name)
        call cantera_coeff(mechanism_path,nsize_mech_name,coeff_cp_k, W_k)
#endif
        call chm_bc_type_species_CMC
        call chm_save_unconditional_fields_bc_CMC
        call chm_fields  ! Initialization by fields activated
        call chm_compute_initial_fields_CMC
        call chm_initialization_domain_CMC

     end if
  endif

end subroutine chm_update_model
