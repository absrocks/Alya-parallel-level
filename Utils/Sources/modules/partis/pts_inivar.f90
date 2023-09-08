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
  use def_domain
  use def_partis
  use def_solver 
  use mod_memory,              only : memory_alloca
  use mod_communications,      only : PAR_SUM
  use mod_pts_injection,       only : pts_injection_initialization
  use mod_pts_parallelization, only : pts_parallelization_initialization
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: itype

  select case( itask )

  case( 0_ip )
     !
     ! Postprocess
     !
     postp(1) % wopos (1, 1) = 'PARTI'
     postp(1) % wopos (1, 2) = 'DEPOS'
     postp(1) % wopos (1, 3) = 'SLIPW'
     postp(1) % wopos (1, 4) = 'FRICT'
     postp(1) % wopos (1, 5) = 'RESID'
     postp(1) % wopos (1, 6) = 'DEPOB'
     postp(1) % wopos (1, 7) = 'DEPOT'
     postp(1) % wopos (1, 8) = 'BOUNC'

     postp(1) % wopos (2, 1) = 'MULTI'
     postp(1) % wopos (2, 2) = 'MULTI'
     postp(1) % wopos (2, 3) = 'SCALA'
     postp(1) % wopos (2, 4) = 'SCALA'
     postp(1) % wopos (2, 5) = 'MULTI'
     postp(1) % wopos (2, 6) = 'SCALA'
     postp(1) % wopos (2, 7) = 'SCALA'
     postp(1) % wopos (2, 8) = 'SCALA'

     postp(1) % wopos (3, 2) = 'NBOUN'
     postp(1) % wopos (3, 6) = 'NBOUN'

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
     !
     ! Initialization particle
     !
     lagrtyp_init % ilagr          = 0
     lagrtyp_init % itype          = 1
     lagrtyp_init % kfl_exist      = 0
     lagrtyp_init % ielem          = 0
     lagrtyp_init % ittim          = 0
     lagrtyp_init % boundary_set   = 0
     lagrtyp_init % iboun          = 0
     lagrtyp_init % coord          = 0.0_rp
     lagrtyp_init % veloc          = 0.0_rp
     lagrtyp_init % accel          = 0.0_rp
     lagrtyp_init % coord_k        = 0.0_rp
     lagrtyp_init % coord_km1      = 0.0_rp
     lagrtyp_init % v_fluid_k      = 0.0_rp
     lagrtyp_init % v_fluid_km1    = 0.0_rp
     lagrtyp_init % v_fluid_km2    = 0.0_rp
     lagrtyp_init % acced          = 0.0_rp
     lagrtyp_init % accee          = 0.0_rp
     lagrtyp_init % acceg          = 0.0_rp 
     lagrtyp_init % stret          = 1.0_rp
     lagrtyp_init % t_inject       = 0.0_rp
     lagrtyp_init % t              = 0.0_rp
     lagrtyp_init % dt_k           = dtime
     lagrtyp_init % dt_km1         = dtime
     lagrtyp_init % dt_km2         = dtime
     lagrtyp_init % dtg            = dtime
     lagrtyp_init % Cd             = 0.0_rp
     lagrtyp_init % Re             = 0.0_rp
     lagrtyp_init % Stk            = 0.0_rp
     lagrtyp_init % dista          = 0.0_rp
     lagrtyp_init % coord1d        = 0.0_rp
     lagrtyp_init % sign           = 1.0_rp
     lagrtyp_init % prope(1:mlapr) = parttyp(1) % prope(1:mlapr)
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

     deposition_var_pts( pts_name_to_variable_number(    'T')) = .true.
     deposition_var_pts( pts_name_to_variable_number('ILAGR')) = .true.
     deposition_var_pts( pts_name_to_variable_number('ITYPE')) = .true.
     deposition_var_pts( pts_name_to_variable_number('COORX')) = .true.
     deposition_var_pts( pts_name_to_variable_number('COORY')) = .true.
     deposition_var_pts( pts_name_to_variable_number('COORZ')) = .true.
     deposition_var_pts( pts_name_to_variable_number('EXIST')) = .true.
     deposition_var_pts( pts_name_to_variable_number(  'SET')) = .true.
     !
     ! Injection and parallelization
     !
     call pts_injection_initialization()
     call pts_parallelization_initialization()
     
     !postprocess_var_pts( 1: 9) = .true.    ! t,ilagr,itype,coord,veloc,dt_k
     !postprocess_var_pts(14:16) = .true.    ! dt_k,kfl_exist,Cd     
     !deposition_var_pts(  1: 6) = .true.    ! t,ilagr,itype,coord
     !deposition_var_pts(    15) = .true.    ! kfl_exist
     !deposition_var_pts(    31) = .true.    ! boundary_set
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
  case( 2_ip )
     !after allocating variables
     if(INOTMASTER) then
         if ( veloc_field_id>0 ) then
            !velocity, if preloaded velocity is used
            !setting veloc field since advec is not yet pointing to veloc
            if(kfl_field(4,veloc_field_id) == 1) then !only one timestep, these will stay constant in time
                veloc(:,:,1) = xfiel(veloc_field_id) % a(:,:,1) !first timestep
                veloc(:,:,3) = xfiel(veloc_field_id) % a(:,:,1) !previous timestep = first timestep
            else
                !set advec(:,:,1) to the last time step (N-1th timestep, since Nth timestep = 1st), since it will be copied 
                !to advec(:,:,3) in pts_begste and then advec(:,:,1) will be set to the first timestep
                !veloc(:,:,1) = xfiel(veloc_field_id) % a( :, :, kfl_field(4,veloc_field_id)-1 )
                veloc(:,:,1) = xfiel(veloc_field_id) % a( :, :, 1 )  !should be equivalent to the above
            end if
         end if
    end if
  end select

end subroutine pts_inivar
