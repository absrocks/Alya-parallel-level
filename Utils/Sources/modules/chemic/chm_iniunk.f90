subroutine chm_iniunk()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_iniunk
  ! NAME 
  !    chm_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the concentrations.
  ! USED BY
  !    chm_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ker_proper
  implicit none
  integer(ip) :: iclas,ipoin,icomp


  avtim_chm = 0.0_rp  ! Accumulated time for time-averaging variables

  if( kfl_rstar == 0 ) then

     if( INOTMASTER ) then
        !
        ! Load initial conditions
        !
        if(kfl_field_chm(1) /= 0_ip) call chm_fields  ! Initialization by fields activated

        !
        ! Initial conditions come from a function
        !
        if(kfl_meteo_chm == -1) call chm_initial_conditions

        icomp = min(3_ip,ncomp_chm)
        do ipoin = 1,npoin
           do iclas = 1,nclas_chm
              if( kfl_fixno_chm(iclas,ipoin) >= 0 ) then
                 conce(ipoin,iclas,icomp) = bvess_chm(iclas,ipoin)
              end if
              conce(ipoin,iclas,1) = conce(ipoin,iclas,icomp)
           end do
        end do

        !
        ! Load external fields
        ! 
        call chm_updfie()

        !
        ! Reading temperature from fields for initilization
        !
        if (kfl_field_chm(2) /= 0_ip) then
           do ipoin=1,npoin
              tempe_chm(ipoin) = xfiel(kfl_field_chm(2))%a(1,ipoin,1)
           end do
        end if

        if (kfl_model_chm==4) then

           call chm_upwmea(1_ip) ! Init mean molecular weight
           call chm_upcpmu(1_ip) ! Init enthalpy and Cp
           call ker_updpro() 

        elseif (kfl_model_chm==5) then

           !
           ! Initialization of hot gases in the reacting layer
           !
           if (kfl_uncfi_chm == 1_ip) call chm_scale_cfi()
           !
           !
           !
           call chm_reatab()     ! Properties update from the table
           call chm_nsacfi()     ! send chemical heat to gauss points as needed by nastal
           call chm_upwmea(4_ip) ! wmean(ipoin,1) ==> wmean(ipoin,2) 
           call chm_upwmea(5_ip) ! wmean(ipoin,1) ==> wmean(ipoin,3)

        endif

     end if

  else
     !
     ! Read restart file
     !
     call chm_restar(1_ip)
     call chm_updunk(6_ip)
   
     avtim_chm = cutim ! Accumulated time for time-averaging variables

     if( INOTMASTER ) then
        !
        ! Load external fields
        ! 
        call chm_updfie()
        !
        ! Init secondary variables for export
        !
        if (kfl_model_chm==4) then
           call chm_upwmea(1_ip) ! Init mean molecular weight
           call chm_upcpmu(1_ip) ! Init enthalpy and Cp
           !! For some reason Kernel didn´t notice properties were updated, this is nasty patch for a bug
           call ker_updpro() 
        elseif (kfl_model_chm==5) then
           call chm_reatab()
           call chm_upwmea(4_ip) ! wmean(ipoin,1) ==> wmean(ipoin,2) 
           call chm_upwmea(7_ip) ! wmean(ipoin,1) ==> wmean(ipoin,3) 
           call chm_upwmea(5_ip) ! wmean(ipoin,3) ==> wmean(ipoin,4)
        endif
     endif

     if (kfl_model_chm==4) then
        call chm_upwmea(1_ip) ! Init mean molecular weight
        call chm_upcpmu(1_ip) ! Init enthalpy and Cp
        !! For some reason Kernel didn´t notice properties were updated, this is nasty patch for a bug
        call ker_updpro() 
     elseif (kfl_model_chm==5) then
        call chm_reatab()
        call chm_upwmea(4_ip) ! wmean(ipoin,1) ==> wmean(ipoin,2) 
        call chm_upwmea(7_ip) ! wmean(ipoin,1) ==> wmean(ipoin,3) 
        call chm_upwmea(5_ip) ! wmean(ipoin,3) ==> wmean(ipoin,4)
     endif

  end if

  !
  ! Mechano-biological model
  !
  if( kfl_model_chm == 3 .and. wprob_chm == 'OSTE1' ) then
     call chm_proads()
  end if

  !-------------------------------------------------------------------
  !
  ! Interpolation from coarse to fine mesh
  !
  !-------------------------------------------------------------------

  if(kfl_meshi_chm /= 0_ip) call chm_coarfine(1_ip)


end subroutine chm_iniunk
