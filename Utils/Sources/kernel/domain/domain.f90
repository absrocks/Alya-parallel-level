!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    domain.f90
!> @author  Guillaume Houzeaux
!> @brief   Domain construction
!> @details Perform the operations
!>          needed to build up the domain data for the run.
!>          All the arrays computed here only depend on the kernel
!>          requierements. Other arrays required e.g. by the modules are
!>          computed later on in Turnon.
!>          -  Subdivide mesh
!>          - Define domain variables
!>          - Compute shape functions & derivatives
!>
!> @}
!-----------------------------------------------------------------------
subroutine domain()

  use def_kintyp,                only : ip,rp
  use def_master,                only : cpu_start
  use def_master,                only : CPU_MESH_MULTIPLICATION
  use def_master,                only : CPU_CONSTRUCT_DOMAIN
  use def_master,                only : INOTMASTER,INOTSLAVE
  use def_master,                only : IPARALL
  use def_master,                only : kfl_algor_msh
  use def_kermod,                only : kfl_elm_graph
  use def_kermod,                only : ndivi,kfl_graph
  use def_kermod,                only : kfl_coo,kfl_ell
  use def_domain,                only : nzdom,mepoi
  use def_domain,                only : r_dom,c_dom
  use def_domain,                only : elmar,meshe
  use def_domain,                only : lelpo,pelpo
  use def_domain,                only : ompss_domains 
  use def_domain,                only : ompss_boundaries 
  use mod_parall,                only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,                only : par_hybrid
  use mod_parall,                only : PAR_OMPSS
  use mod_parall,                only : par_omp_nelem_chunk
  use mod_parall,                only : par_omp_nboun_chunk
  use mod_parall_openmp,         only : parall_openmp_chunk_sizes
  use mod_parall_openmp,         only : parall_openmp_coloring
  use mod_graphs,                only : graphs_csr_to_coo
  use mod_graphs,                only : graphs_csr_to_ell
  use mod_graphs,                only : graphs_element_element_graph
  use mod_parall_openmp,         only : parall_openmp_partition_and_adjacency_ompss
  use mod_ghost_geometry,        only : par_ghost_geometry
  use mod_ADR,                   only : ADR_load_mesh
  use mod_par_additional_arrays, only : par_multiplicity_ownership
  use mod_par_additional_arrays, only : par_matrix_exchange_on_interface_nodes
  use mod_par_additional_arrays, only : par_matrix_w_halos_exchange_on_interface_nodes
  use mod_par_additional_arrays, only : par_node_number_in_owner
  use mod_par_additional_arrays, only : par_full_row_communication_arrays
  use mod_par_additional_arrays, only : par_matrix_computational_halo_exchange
  use mod_par_output_partition,  only : par_output_partition
  use mod_par_output_partition,  only : par_output_global_matrix
  use mod_unity_tests,           only : unity_tests_integration_rules
  use mod_unity_tests,           only : unity_tests_check_halos
  use mod_commdom_driver,        only : CNT_CPLNG,commdom_driver_init,commdom_driver_set_mesh !< 2016MAR30
  use mod_measurements,          only : measurements_set_function                             !< 2017JAN07
  use mod_mesh_type,             only : mesh_type_allocate_initialize
  use mod_mesh_type,             only : mesh_type_save_original_mesh
  use mod_mesh_type,             only : mesh_type_update_last_mesh
  use mod_mesh_multiplication,   only : mesh_multiplication
  use mod_mesh_multiplication,   only : mesh_multiplication_node_codes
  use mod_renumbering,           only : renumbering_temporal_locality
  use mod_messages,              only : messages_live
  use mod_element_data_base,     only : element_data_base_save
  use mod_par_bin_structure,     only : par_bin_structure
  use mod_communications,        only : PAR_BARRIER
  use mod_materials,             only : materials_on_nodes
  use mod_periodicity,           only : periodicity_enhance_node_graph
  use mod_outfor,                only : outfor
  use mod_output,                only : output_domain

  use def_domain
  use def_master
  implicit none
  real(rp) :: time1,time2,time3
 
  integer(ip) :: ipoin,ielem,inode
  
  !----------------------------------------------------------------------
  !
  ! From now on, Master does not have any mesh info :o(
  !
  !----------------------------------------------------------------------
  !  
  ! First things to do before continuing...
  !
  call domvar(1_ip)                                       ! LFACE
  call domvar(2_ip)                                       ! LNUTY, LTYPF, NNODF, LPERI
  !
  ! Mesh multiplication
  !
  call cputim(time1) 

  call mesh_type_allocate_initialize()                    ! MESHE
  call mesh_type_save_original_mesh()                     ! MESHE(0)

  call mesh_multiplication()

  call mesh_type_update_last_mesh()                       ! MESHE(NDIVI) =>

  call cputim(time2) 
  cpu_start(CPU_MESH_MULTIPLICATION) = time2 - time1
  !
  ! Renumber nodes to have own nodes
  ! (interior+own boundary) first
  !
  call par_renumber_nodes()                               ! All nodal arrays read in reageo
  !
  ! Global numbering
  !
  call par_global_numbering()
  call PAR_BARRIER()

  !----------------------------------------------------------------------
  !
  ! Coupling initialization
  !
  !----------------------------------------------------------------------
  !
  ! Zones and subdomains
  !
  call in_zone_in_subd()
  !
  ! Compute parallel communicators
  !
  call commdom_driver_init( CNT_CPLNG )                   !< 2016ABRIL04
  call par_color_communicators()
  !
  ! Bin structure for partitions
  !
  call PAR_BARRIER()
  call messages_live('PARALL: COMPUTE BIN STRUCTURE FOR PARTITIONS')
  call par_bin_structure()

  !----------------------------------------------------------------------
  !
  ! Other mesh arrays
  !
  !----------------------------------------------------------------------
  !
  ! Compute shape functions & derivatives
  !
  call cshder(3_ip)                                       ! SHAPE, DERIV, HESLO, WEIGP...116
  !
  ! Turn on the remeshing algorithm
  !
  call PAR_BARRIER()
  if( kfl_algor_msh > 0 ) call Newmsh(0_ip)
  !
  ! Create mesh graph
  !
  call PAR_BARRIER()
  call domgra(2_ip)                                       ! R_DOM, C_DOM, NZDOM, LELPO, PELPO
  !
  ! Extended graph
  !
  if( kfl_graph == 1 ) call par_extgra()                  ! R_DOM, C_DOM, NZDOM
  !
  ! Point to mesh type
  !
  meshe(ndivi) % r_dom => r_dom
  meshe(ndivi) % c_dom => c_dom
  meshe(ndivi) % nzdom =  nzdom 
  !
  ! Groups of deflated: 1 group / subdomain
  !
  call grodom(2_ip)                                       ! LGROU_DOM
  !
  ! Renumber elements and check spatial/temporal locality
  !
  call PAR_BARRIER()
  call renelm()                                           ! LTYPE, LELCH, LNNOD, LESUB, LMATE, LNODS LEINV_LOC, LBOEL, LESET, XFIEL, LELPO, PELPO
  !call renumbering_temporal_locality(meshe(ndivi))
  !
  ! Domain variables depending on mesh
  !
  call PAR_BARRIER()
  call domvar(3_ip)                                       ! LNNOB, LGAUS, LNOCH, LBONO, LMAST
  !
  ! Save element data base
  !
  call element_data_base_save()                           ! GPCAR, GPVOL, HLENG, TRAGL, GPHES
  !
  ! Edges connectivity, graphs, parallel data structure, etc
  !
  call PAR_BARRIER()
  call edge_data_structures()                             ! LEDGS, LNNED, EDGE_TO_NODE, LOCAL_TO_GLOBAL_EDGE...

  !----------------------------------------------------------------------
  !
  ! HYBRID PARALLELIZATION
  !
  !----------------------------------------------------------------------
  !
  ! Output partition info 
  !
  call par_output_info_partition()
  !
  ! Color elements for OMP
  !
  call parall_openmp_coloring(meshe(ndivi),mepoi,pelpo,lelpo)
  call parall_openmp_coloring(meshe(ndivi),ON_BOUNDARIES=.true.)
  !
  ! Size of blocks and chunk sizes for OpenMP
  !
  call parall_openmp_chunk_sizes(meshe(ndivi))            ! PAR_OMP_NELEM_CHUNK, PAR_OMP_NPOIN_CHUNK, PAR_OMP_NBOUN_CHUNK
  !
  ! OMPSS multidependencies: Local mesh partition
  !
  if( INOTMASTER .and. par_hybrid == PAR_OMPSS ) then
     call parall_openmp_partition_and_adjacency_ompss(&        
          par_omp_nelem_chunk,meshe(ndivi),ompss_domains)
     call parall_openmp_partition_and_adjacency_ompss(&         
          par_omp_nboun_chunk,meshe(ndivi),ompss_boundaries,ON_BOUNDARIES=.true.)
  end if

  !----------------------------------------------------------------------
  !
  ! Other mesh stuffs
  !
  !----------------------------------------------------------------------
  !
  ! Modify node graph to account for periodicity
  !
  call periodicity_enhance_node_graph('LOCAL GRAPH')      ! R_DOM, C_DOM, NZDOM
  !
  ! Element to csr
  !
  call elecsr()
  !
  ! Compute boundary connectivity
  !
  call setlbe()                                           ! LBOEL
  !
  ! Check mesh
  !
  call mescek(1_ip)
  !
  ! Materials
  !
  call materials_on_nodes()                               ! LMATN
  !
  ! Prepare voxels postprocess
  !
  call prevox()
  !
  ! Parall: Fringe geometry
  !
  call par_ghost_geometry()
  !
  ! Full row graph for own nodes including halo nodes     ! R_DOM_OWN, C_DOM_OWN
  !
  call messages_live('COMPUTE GRAPH FOR OWN NODES (FOR FULL ROW MATRIX FORMAT)')
  call full_row_graph(meshe(ndivi))
  !call periodicity_enhance_node_graph('FULL ROWS GRAPH')

  !----------------------------------------------------------------------
  !
  ! Parallel arrays
  !
  !----------------------------------------------------------------------
  !
  ! Parall communication arrays depending on zones
  !
  call Parall(608_ip)
  call par_zone_communication_arrays()                    ! PAR_COMM_COLOR_ARRAY
  !
  ! Some parallel addition arrays
  !
  if( IPARALL ) then
     call par_multiplicity_ownership(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1))             ! Multiplicity and ownership
     call par_matrix_exchange_on_interface_nodes(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1)) ! Must be done after par_extgra
     call par_node_number_in_owner(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1))               ! Obtains COMM % node_number_in_owner
     call par_full_row_communication_arrays(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1))      ! Communicator for full row exchange
#if defined(PARAL_AGMG) || defined(INC_PSBLAS)
     call par_matrix_w_halos_exchange_on_interface_nodes(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1))
#else
     call par_matrix_computational_halo_exchange(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1))
#endif
  end if
  !
  ! Initialize Elsest
  !
  call PAR_BARRIER()
  call elsini()
  !
  ! Witness point information
  !
  call PAR_BARRIER()
  call witnes()                                           ! LEWIT, SHWIT

  !----------------------------------------------------------------------
  !
  ! Coupling computations
  !
  !----------------------------------------------------------------------

  call PAR_BARRIER()
#ifndef COMMDOM  
  !
  call cou_turnon()
  call cou_define_wet_geometry()                          ! COUPLING_TYPE
  call cou_initialize_coupling()                          ! COUPLING_TYPE
  call cou_inivar(3_ip)
  call cou_communication_arrays()
  call cou_output()
#else
  !
  !!call commdom_driver_set_mesh( CNT_CPLNG )               !< 2016ABR12
  call measurements_set_function( commdom_driver_set_mesh, CNT_CPLNG, -1, -ID_KERMOD )    !< 2017JAN07
  !
#endif
  !
  ! Output partitioning and coupling
  !  
  call par_output_partition(1_ip)
  call par_output_global_matrix()
  !
  ! Compute mesh dependent variables
  !
  call domarr(1_ip)                                       ! VMASS, VMASC, EXNOR, YWALP, YWALB, WALLD, SKCOS, COORD
  !
  ! Operations on boundary conditions
  !
  call opebcs(1_ip)
  !
  ! Modify code if mesh subdivision has been used
  !
  call mesh_multiplication_node_codes()
  !
  ! Output domain mesh
  !
  call output_domain()
  !
  ! Geometrical normals and check codes
  !
  call opebcs(2_ip)
  call chkcod()
  !
  ! Edge boundary conditions
  !
  call edge_boundary_conditions()
  !
  ! Load mesh array and parameters for ADR toolbox
  !
  call ADR_load_mesh(meshe(ndivi:ndivi),elmar)
  !
  ! Graph in COO format
  !
  if( INOTMASTER .and. kfl_coo /= 0 ) then
     call graphs_csr_to_coo(&
          meshe(ndivi) % npoin,1_ip,meshe(ndivi) % r_dom   ,meshe(ndivi) % c_dom,&
          meshe(ndivi) % nzdom,meshe(ndivi) % coo_rows,meshe(ndivi) % coo_cols)
  end if
  !
  ! Graph in ELL format
  !
  if( INOTMASTER .and. kfl_ell /= 0 ) then
     call graphs_csr_to_ell(&
          meshe(ndivi) % npoin,1_ip,meshe(ndivi) % r_dom   ,meshe(ndivi) % c_dom,&
          meshe(ndivi) % nzdom_ell,meshe(ndivi) % ell_cols)
  end if
  !
  ! Element graph with halos
  !
  if( INOTMASTER .and. kfl_elm_graph == 1 ) then
     call messages_live('COMPUTE ELEMENT-ELEMENT GRAPH INCLUDING HALO')
     call graphs_element_element_graph(meshe(ndivi),'SHARING FACES','INCLUDING HALO','INCLUDING DIAGONAL')
  end if

  !----------------------------------------------------------------------
  !
  ! Unity tests
  !
  !----------------------------------------------------------------------

  if( INOTSLAVE ) call unity_tests_integration_rules()
  if( IPARALL )   call unity_tests_check_halos()
  !
  ! Close domain unit
  !
  call opfdom(3_ip)

  call cputim(time3)
  cpu_start(CPU_CONSTRUCT_DOMAIN) = time3 - time2

end subroutine domain 

subroutine cacapopo()

  use def_master
  use def_domain
  use mod_maths
  use mod_matrix
  use mod_parall
  implicit none
  integer(ip)          :: lsize,ii,nn,nz,iboun
  integer(ip), pointer :: kpoin(:)
  real(rp),    pointer :: xcoor(:,:)
  integer(ip), pointer :: ia(:)
  integer(ip), pointer :: ja(:)
  real(rp),    pointer :: aa(:)
  real(rp),    pointer :: diag(:)
  
  do iboun = 1,nboun
     if(lbinv_loc(iboun)==6) then
        print*,'a=',lninv_loc(lnodb(:,iboun))
        print*,'b=',xfiel(3) % a(:,iboun,1)
     end if
  end do
  call runend('O.K.!')
  if(inotmaster) then
     print*,kfl_paral,leinv_loc(1)
  end if
  call runend('O.K.!')

  goto 10
  !
  ! Unity test for matrix scaling
  !
  !                   +-          -+
  !                   |  1   2   3 |
  ! Original matrix = |  4   5   6 |
  !                   |  7   8   9 |
  !                   +-          -+
  !
  !                   +-          -+
  !                   |  2   0   0 |
  ! Scaling =         |  0   4   0 |
  !                   |  0   0   8 |
  !                   +-          -+
  !
  !                   +-          -+
  !                   |  2   4   6 |
  ! Left scaling =    | 16  20  24 |
  !                   | 56  64  72 |
  !                   +-          -+
  !
  !                   +-          -+
  !                   |  2   8  24 |
  ! Right scaling =   |  8  20  48 |
  !                   | 14  32  72 |
  !                   +-          -+
  !
  nn = 3
  nz = 9
  allocate(ia(nn+1),ja(nz))
  allocate(aa(nz))
  allocate(diag(nn))
  ia   = (/ 1,4,7,10 /)
  ja   = (/ 1,2,3,1,2,3,1,2,3 /)
  aa   = (/ 1.0_rp,2.0_rp,3.0_rp,4.0_rp,5.0_rp,6.0_rp,7.0_rp,8.0_rp,9.0_rp /)
  diag = (/ 2.0_rp,4.0_rp,8.0_rp /)

  call matrix_scaling_CSR(1_ip,nn,1_ip,ia,ja,aa,diag,LEFT_SCALING=.true.,RIGHT_SCALING=.false.)

  do ii = 1,nn     
     write(*,*) aa(ia(ii):ia(ii+1)-1)
  end do
  
  stop
  !
  ! Unity test
  !  
  if( INOTMASTER ) then
     lsize = 8
     allocate(kpoin(lsize))
     allocate(xcoor(ndime,lsize))
     do ii = 1,lsize
        kpoin(ii) = ii
     end do
     xcoor(1:3,1) = (/ -1.0_rp,  2.0_rp,  3.0_rp /)
     xcoor(1:3,2) = (/ -1.0_rp, -2.0_rp,  4.0_rp /)
     xcoor(1:3,3) = (/  0.0_rp, -2.0_rp, -7.0_rp /)
     xcoor(1:3,4) = (/ -1.0_rp, -2.0_rp,  2.0_rp /)
     xcoor(1:3,5) = (/ -1.0_rp, -2.0_rp, -1.0_rp /)
     xcoor(1:3,6) = (/  0.0_rp, -2.0_rp, -8.0_rp /)
     xcoor(1:3,7) = (/ -3.0_rp,  2.0_rp, 12.0_rp /)
     xcoor(1:3,8) = (/  0.0_rp, -2.0_rp,  8.0_rp /)
     call maths_geometrical_sort_using_coordinates(2_ip,ndime,lsize,xcoor,kpoin)
     do ii = 1,lsize
        print*,kpoin(ii),xcoor(1:3,ii)
     end do
  end if

10 continue
  if( INOTMASTER ) THEN
     !ii = PAR_GLOBAL_TO_LOCAL_NODE(103)
     !if( ii /= 0 ) print*,'popo=',kfl_paral,lninv_loc(ii)
     !ii = PAR_GLOBAL_TO_LOCAL_NODE(121)
     !if( ii /= 0 ) print*,'popo=',kfl_paral,lninv_loc(ii)
     !ii = PAR_GLOBAL_TO_LOCAL_NODE(4)
     !if( ii /= 0 ) print*,'popo=',kfl_paral,lninv_loc(ii)
  end if
  
  call runend('O.K.!')
end subroutine cacapopo
