subroutine pts_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Partis/pts_inivar
  ! NAME 
  !    pts_inivar
  ! DESCRIPTION
  !    Initial variables
  ! USES
  ! USED BY
  !    pts_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod,              only: gasco
  use def_domain
  use def_partis
  use def_solver 
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_communications,      only : PAR_SUM
  use mod_pts_injection,       only : pts_injection_initialization
  use mod_pts_parallelization, only : pts_parallelization_initialization
  use mod_physics,             only : physics_initialize_liquid
  use mod_physics,             only : universal_gas_constant 
  use mod_interp_tab,          only : fw_lookup 
  use mod_arrays,              only : arrays_register
      
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: itype,ii,ivari,ivar2
  real(rp)                :: y_c_eq, y_c_0
  real(rp), pointer       :: retva(:)
  real(rp), pointer       :: control(:)
  real(rp), pointer       :: tab_conce(:)

  select case( itask )

  case( 0_ip )
     !
     ! Postprocess
     !
     call arrays_register( 12_ip,(/'MOMSK','VECTO','NPOIN','PRIMA'/),momentum_sink,ENTITY_POSITION=2_ip)
     call arrays_register( 11_ip,(/'HEASK','SCALA','NPOIN','PRIMA'/),heat_sink    ,ENTITY_POSITION=1_ip)
     call arrays_register( 10_ip,(/'MASSK','SCALA','NPOIN','PRIMA'/),mass_sink    ,ENTITY_POSITION=1_ip)
     
     call arrays_register( 13_ip,(/'DEPOE','SCALA','NELEM','PRIMA'/),depoe_pts,    ENTITY_POSITION=2_ip)
     call arrays_register(  6_ip,(/'DEPOB','SCALA','NBOUN','PRIMA'/),depob_pts,    ENTITY_POSITION=2_ip)     
     call arrays_register(  5_ip,(/'RESID','SCALA','NELEM','PRIMA'/),resid_pts,    ENTITY_POSITION=2_ip)

     call arrays_register( 14_ip,(/'NRESI','SCALA','NPOIN','SECON'/),              ENTITY_POSITION=2_ip,COMPONENT_POSITION=1_ip)
     call arrays_register(  1_ip,(/'PARTI','SCALA','NPOIN','SECON'/),              ENTITY_POSITION=2_ip,COMPONENT_POSITION=1_ip)

     call arrays_register( 15_ip,(/'AVMIE','SCALA','NPOIN','PRIMA'/),avg_mie_pts  ,ENTITY_POSITION=1_ip)

     postp(1) % wopos (1, 2) = 'DEPOS'
     postp(1) % wopos (1, 3) = 'SLIPW'
     postp(1) % wopos (1, 4) = 'FRICT'
     postp(1) % wopos (1, 7) = 'DEPOT'
     postp(1) % wopos (1, 8) = 'BOUNC'
     postp(1) % wopos (1, 9) = 'TEMPE'

     postp(1) % wopos (2, 2) = 'MULTI'
     postp(1) % wopos (2, 3) = 'SCALA'
     postp(1) % wopos (2, 4) = 'SCALA'
     postp(1) % wopos (2, 7) = 'SCALA'
     postp(1) % wopos (2, 8) = 'SCALA'
     postp(1) % wopos (2, 9) = 'SCALA'

     postp(1) % wopos (3, 2) = 'NBOUN'
     !
     ! Witness geometry variables
     !
     postp(1) % wowig  ( 1)  = 'NUMBE' ! Number of particles     
     postp(1) % wowig  ( 2)  = 'DIA10' ! Diameter 
     postp(1) % wowig  ( 3)  = 'DIA32' ! SMD
     postp(1) % wowig  ( 4)  = 'TEMPE' ! Temperature
     postp(1) % wowig  (5:7) = 'VELOC' ! Velocity
     postp(1) % wowig  ( 8)  = 'DIA20' ! Diameter square
     postp(1) % wowig  ( 9)  = 'TEMP2' ! Temperature square
     postp(1) % wowig (10:12)= 'VELO2' ! Velocity square
     postp(1) % wowig  (13)  = 'MASS ' ! Mass
     postp(1) % wowig (14:16)= 'MFLUX' ! Mass flux
     !
     ! Format
     !
     kfl_ptsres_binary     = 0       ! pts.res saved as text by default     
     kfl_ptsres_split      = 0       ! pts.res not split by default 
     !
     ! Others
     !
     nlagr                 = 0       ! Absolute number of lagrangian particles
     nlacc_pts             = 0       ! Number of accumulated particles
     kfl_depos_pts         = 0       ! Deposition map on mesh is not postprocess
     kfl_timei             = 1       ! Problem is always transient
     particles_sent        = 0       ! # of particles received
     particles_recv        = 0       ! # of particles sent
     comm_loops_pts        = 0       ! # communication loops 
     kfl_resid_pts         = 0       ! Residence time not required
     nlapr                 = 0       ! # of properties
     kfl_injec             = 0       ! No injection
     nlagr_local_pts       = 0       ! Position taken in type
     nlagr_free_pts        = 0       ! Free positions taken in type
     kfl_slip_wall_pts     = 0       ! Slip wall
     kfl_bouncing_wall_pts = 0       ! Bouncing wall
     nvarp_pts             = 0       ! Number of deposition variables
     nvard_pts             = 0       ! Number of deposition variables
     param_veloc_pts       = 0.0_rp  ! Parameters for velocity injection
     param_tempe_pts       = 0.0_rp  ! Parameters for temperature injection
     param_mass_pts        = 0.0_rp  ! Parameters for mass injection
 
     nullify(permu_nlagr_pts)        ! Permutation
     nullify(depoe_pts)              ! Deposition map over elements
     nullify(depob_pts)              ! Deposition map over boundaries
     nullify(defor_pts)              ! Deformation tensor
     nullify(hleng_pts)              ! Element characteristics length
     nullify(lboue_pts)
     nullify(kfl_fixbo_pts)     
     nullify(bvnat_pts)  
     nullify(tbcod_pts)
     nullify(leleboun_pts) 
     nullify(bouno_pts)  
     nullify(walld_slip_pts)
     nullify(kfl_fixno_walld_slip_pts)
     nullify(friction_pts)
     nullify(walld_bouncing_pts)
     nullify(kfl_fixno_walld_bouncing_pts)
     nullify(resid_pts)
     nullify(postprocess_list_pts)
     nullify(deposition_list_pts)
     nullify(avg_mie_pts)
     !
     ! Postprocess variables
     !     
     postprocess_name_pts     = 'NULL'
     
     postprocess_name_pts( 1) = 'T'       ! t              
     postprocess_name_pts( 2) = 'ILAGR'   ! ilagr          
     postprocess_name_pts( 3) = 'ITYPE'   ! itype          
     postprocess_name_pts( 4) = 'IELEM'   ! ielem          
     postprocess_name_pts( 5) = 'EXIST'   ! kfl_exist      
     postprocess_name_pts( 6) = 'ITTIM'   ! ittim          
     postprocess_name_pts( 7) = 'SET  '   ! boundary_set     
     postprocess_name_pts( 8) = 'COORX'   ! coord(1)       
     postprocess_name_pts( 9) = 'COORY'   ! coord(2)       
     postprocess_name_pts(10) = 'COORZ'   ! coord(3)       
     postprocess_name_pts(11) = 'VELOX'   ! veloc(1)       
     postprocess_name_pts(12) = 'VELOY'   ! veloc(2)       
     postprocess_name_pts(13) = 'VELOZ'   ! veloc(3)       
     postprocess_name_pts(14) = 'ACCEX'   ! accel(1)       
     postprocess_name_pts(15) = 'ACCEY'   ! accel(2)       
     postprocess_name_pts(16) = 'ACCEZ'   ! accel(3)       
     postprocess_name_pts(17) = 'DTK  '   ! dt_k           
     postprocess_name_pts(18) = 'CD   '   ! Cd             
     postprocess_name_pts(19) = 'STK1 '   ! Stk(1)         
     postprocess_name_pts(20) = 'STK2 '   ! Stk(2)         
     postprocess_name_pts(21) = 'VF1  '   ! v_fluid_k(1)   
     postprocess_name_pts(22) = 'VF1  '   ! v_fluid_k(2)   
     postprocess_name_pts(23) = 'VF1  '   ! v_fluid_k(3)   
     postprocess_name_pts(24) = 'ACCDX'   ! acced(1)       
     postprocess_name_pts(25) = 'ACCDY'   ! acced(2)       
     postprocess_name_pts(26) = 'ACCDZ'   ! acced(3)       
     postprocess_name_pts(27) = 'ACCFX'   ! accee(1)       
     postprocess_name_pts(28) = 'ACCFY'   ! accee(2)       
     postprocess_name_pts(29) = 'ACCFZ'   ! accee(3)       
     postprocess_name_pts(30) = 'ACCGX'   ! acceg(1)       
     postprocess_name_pts(31) = 'ACCGY'   ! acceg(2)       
     postprocess_name_pts(32) = 'ACCGZ'   ! acceg(3)       
     postprocess_name_pts(33) = 'STRET'   ! stret          
     postprocess_name_pts(34) = 'DTKM1'   ! dt_km1         
     postprocess_name_pts(35) = 'DTKM2'   ! dt_km2         
     postprocess_name_pts(36) = 'DTG  '   ! dtg            
     postprocess_name_pts(37) = 'VFM1X'   ! v_fluid_km1(1) 
     postprocess_name_pts(38) = 'VFM1Y'   ! v_fluid_km1(2) 
     postprocess_name_pts(39) = 'VFM1Z'   ! v_fluid_km1(3) 
     postprocess_name_pts(40) = 'VFM2X'   ! v_fluid_km2(1) 
     postprocess_name_pts(41) = 'VFM2Y'   ! v_fluid_km2(2) 
     postprocess_name_pts(42) = 'VFM2Z'   ! v_fluid_km2(3) 
     postprocess_name_pts(43) = 'TINJE'   ! t_inject       
     postprocess_name_pts(44) = 'COKX '   ! coord_k(1)     
     postprocess_name_pts(45) = 'COKY '   ! coord_k(2)     
     postprocess_name_pts(46) = 'COKZ '   ! coord_k(3)     
     postprocess_name_pts(47) = 'COK1X'   ! coord_km1(1)   
     postprocess_name_pts(48) = 'COK1Y'   ! coord_km1(2)   
     postprocess_name_pts(49) = 'COK1Z'   ! coord_km1(3)     
     postprocess_name_pts(50) = 'DISTA'   ! dista          
     postprocess_name_pts(51) = 'COO1D'   ! coord1d        
     postprocess_name_pts(52) = 'SIGN '   ! sign           
     postprocess_name_pts(53) = 'TEMPK'   ! tempe_k           
     postprocess_name_pts(54) = 'TEKM1'   ! tempe_km1     
     postprocess_name_pts(55) = 'MASSK'   ! mass_k      
     postprocess_name_pts(56) = 'MAKM1'   ! mass_km1
     postprocess_name_pts(57) = 'MPIRA'   ! mpi_rank
     postprocess_name_pts(58) = 'DIAMK'   ! diam_k      
     postprocess_name_pts(59) = 'TEMPF'   ! Temp_fluid_k
     postprocess_name_pts(60) = 'YVAPF'   ! Yvap_fluid_k
     !                    
     ! Solver
     !
     call soldef(-2_ip)                     ! Allocate memory
     solve(1) % kfl_solve = 1               ! Slip wall solver
     solve(1) % wprob     = 'SLIP_WALL'     ! Name
     solve(1) % kfl_iffix = 2               ! Fix Dirichlet with value given by initial solution

     solve(2) % kfl_solve = 1               ! Bouncing wall solver
     solve(2) % wprob     = 'BOUNCING_WALL' ! Name
     solve(2) % kfl_iffix = 2               ! Fix Dirichlet with value given by initial solution
     !
     ! Others
     ! 
     nlagr_existing_pts      = 0
     nlagr_non_migrating_pts = 0
     nlagr_going_out_pts     = 0
     nlagr_zero_time_pts     = 0
     nlagr_deposited_pts     = 0
     nlagr_hits_wall_pts     = 0
     
     !
     ! Injection
     !  
     call pts_injection_initialization()
   
  case( 1_ip )
     !
     ! After reading the data
     !
     !
     ! Number of used types
     !
     ntyla_pts = 0
     number_types_pts = 0
     do itype = 1,mtyla
        if( parttyp(itype) % kfl_exist == 1 ) then
           ntyla_pts = max(ntyla_pts,itype)
           number_types_pts = number_types_pts + 1           
        end if
     end do


     !
     ! Initialize liquid structure
     !
     do itype = 1,mtyla
         if (( parttyp(itype) % kfl_exist == 1 ) .and. &
             ( parttyp(itype) % kfl_therm /= 0 )) then
             
             !
             ! Deafult pressure 
             !  
             if (prthe(1) < 1000.0_rp) prthe = 101325.0_rp

             if (parttyp(itype) % liq % name == 'USER ') then
                 call physics_initialize_liquid(parttyp(itype) % liq, 298.15_rp, prthe(1), & 
                                     parttyp(itype) % L_vapor, parttyp(itype) % T_boiling, &
                                     parttyp(itype) % Cp, parttyp(itype) % denpa,          & 
                                     parttyp(itype) % w)
             else
                 !
                 ! Using lookup tables
                 !
                 call physics_initialize_liquid(parttyp(itype) % liq, 298.15_rp, prthe(1))
                 
             endif 
             


             if ( parttyp(itype) % kfl_therm == 2 ) then
                 !
                 ! Initialize cp coefficents of fuel
                 !
                 !
                 ! Allocate control variables 
                 !
                 nullify(control) 
                 nullify(retva) 
                 nullify(tab_conce) 

                 call memory_alloca(mem_modul(1:2,modul),'CONTROL'  ,'pts_inivar',control  ,parttyp(itype) % table_fw % main_table % ndim)
                 call memory_alloca(mem_modul(1:2,modul),'RETVA'    ,'pts_inivar',retva    ,parttyp(itype) % table_fw % main_table % nvar)
                 call memory_alloca(mem_modul(1:2,modul),'TAB_CONCE','pts_inivar',tab_conce,parttyp(itype) % table_fw % main_table % ndim)

                 !
                 ! Lookup seen state ASSUMES ADAIBATIC FOR NOW
                 !        
                 control = 0.0_rp
                 do ii = 1, parttyp(itype) % table_fw % main_table % ndim
                    select case (parttyp(itype) % table_fw % main_table % coords(ii) % name)
                    case ('CMEAN','C    ')
                        control(ii) = 0.0_rp
                    case ('CVAR ')
                        control(ii) = 0.0_rp
                    case ('CHIST')
                        control(ii) = 0.0_rp 
                    case ('ZMEAN','Z    ')
                        control(ii) = 1.0_rp
                    case ('ZVAR ')
                        control(ii) = 0.0_rp
                    case ('IMEAN','I    ')
                        !
                        ! At this point the enthalpy shouldn't matter, 
                        ! because the Cp coefficeints of the fule should not change
                        ! with enthalpy
                        !
                        control(ii) = 0.0_rp
                    end select
                 enddo
                 call fw_lookup( control, tab_conce, parttyp(itype) % table_fw, retva )
                 
                 do ii = 1,6 
                   !
                   ! First is low temperature, second is high temperature in table
                   !
                   parttyp(itype) % cpcoef_v_chm(ii,1) = retva(4+ii)
                   parttyp(itype) % cpcoef_v_chm(ii,2) = retva(10+ii)
                 end do
                 !call physics_T_2_HCp(T_mean, cpcoef, H_mean, sphfl)
                 
                 call memory_deallo(mem_modul(1:2,modul),'CONTROL'  ,'pts_inivar',control  )
                 call memory_deallo(mem_modul(1:2,modul),'RETVA'    ,'pts_inivar',retva    )
                 call memory_deallo(mem_modul(1:2,modul),'TAB_CONCE','pts_inivar',tab_conce)
              endif
         endif
     enddo
     !
     ! Postprocess
     !
     nvarp_pts = count(postprocess_var_pts)
     call memory_alloca(mem_modul(1:2,modul),'POSTPROCESS_LIST_PTS','pts_output',postprocess_list_pts,nvarp_pts)
     ivar2 = 0
     do ivari = 1,mvarp_pts
        if( postprocess_var_pts(ivari) ) then
           ivar2 = ivar2 + 1
           postprocess_list_pts(ivar2) = ivari
        end if
     end do
     nvard_pts = count(deposition_var_pts)
     call memory_alloca(mem_modul(1:2,modul),'DEPOSITION_LIST_PTS','pts_output',deposition_list_pts,nvard_pts)
     ivar2 = 0
     do ivari = 1,mvarp_pts
        if( deposition_var_pts(ivari) ) then
           ivar2 = ivar2 + 1
           deposition_list_pts(ivar2) = ivari
        end if
     end do

  case( 2_ip )
      !
      ! Injection and parallelization
      !
!!!      call pts_injection_initialization()
      call pts_parallelization_initialization()

   case ( 3_ip )
      !
     ! Default postprocess
     !
     postprocess_var_pts(pts_name_to_variable_number(    'T')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('ILAGR')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('ITYPE')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('COORX')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('COORY')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('COORZ')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('VELOX')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('VELOY')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('VELOZ')) = .true.
     postprocess_var_pts(pts_name_to_variable_number(  'DTK')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('EXIST')) = .true.
     postprocess_var_pts(pts_name_to_variable_number(   'CD')) = .true.

     if( kfl_thermo_pts /= 0 ) then
        postprocess_var_pts(pts_name_to_variable_number('TEMPK')) = .true.
        postprocess_var_pts(pts_name_to_variable_number('MASSK')) = .true.
     end if

     deposition_var_pts( pts_name_to_variable_number(    'T')) = .true.
     deposition_var_pts( pts_name_to_variable_number('ILAGR')) = .true.
     deposition_var_pts( pts_name_to_variable_number('ITYPE')) = .true.
     deposition_var_pts( pts_name_to_variable_number('COORX')) = .true.
     deposition_var_pts( pts_name_to_variable_number('COORY')) = .true.
     deposition_var_pts( pts_name_to_variable_number('COORZ')) = .true.
     deposition_var_pts( pts_name_to_variable_number('EXIST')) = .true.
     deposition_var_pts( pts_name_to_variable_number(  'SET')) = .true.
    
     !postprocess_var_pts( 1: 9) = .true.    ! t,ilagr,itype,coord,veloc,dt_k
     !postprocess_var_pts(14:16) = .true.    ! dt_k,kfl_exist,Cd     
     !deposition_var_pts(  1: 6) = .true.    ! t,ilagr,itype,coord
     !deposition_var_pts(    15) = .true.    ! kfl_exist
     !deposition_var_pts(    31) = .true.    ! boundary_set

  end select

end subroutine pts_inivar
