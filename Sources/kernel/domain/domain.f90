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
subroutine Domain()

  use def_kintyp,                  only : ip,rp
  use def_master,                  only : cpu_start
  use def_master,                  only : CPU_MESH_MULTIPLICATION
  use def_master,                  only : CPU_CONSTRUCT_DOMAIN
  use def_master,                  only : INOTMASTER,INOTSLAVE
  use def_master,                  only : IPARALL
  use def_master,                  only : ID_KERMOD
  use def_master,                  only : kfl_algor_msh
  use def_kermod,                  only : kfl_elm_graph
  use def_kermod,                  only : ndivi
  use def_kermod,                  only : kfl_coo,kfl_ell
  use def_domain,                  only : nzdom,mepoi
  use def_domain,                  only : r_dom,c_dom
  use def_domain,                  only : elmar,meshe
  use def_domain,                  only : lelpo,pelpo
  use def_domain,                  only : ompss_domains
  use def_domain,                  only : ompss_boundaries
  use mod_parall,                  only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,                  only : par_hybrid
  use mod_parall,                  only : PAR_OMPSS
  use mod_parall,                  only : par_omp_nelem_chunk
  use mod_parall,                  only : par_omp_nboun_chunk
  use mod_par_tools,               only : par_tools_global_numbering
  use mod_parall_openmp,           only : parall_openmp_chunk_sizes
  use mod_parall_openmp,           only : parall_openmp_coloring
  use mod_graphs,                  only : graphs_csr_to_coo
  use mod_graphs,                  only : graphs_csr_to_ell
  use mod_graphs,                  only : graphs_element_element_graph
  use mod_parall_openmp,           only : parall_openmp_partition_and_adjacency_ompss
  use mod_ghost_geometry,          only : par_ghost_geometry
  use mod_ADR,                     only : ADR_load_mesh
  use mod_par_additional_arrays,   only : par_multiplicity_ownership
  use mod_par_additional_arrays,   only : par_matrix_exchange_on_interface_nodes
  use mod_par_additional_arrays,   only : par_matrix_w_halos_exchange_on_interface_nodes
  use mod_par_additional_arrays,   only : par_node_number_in_owner
  use mod_par_additional_arrays,   only : par_full_row_communication_arrays
  use mod_par_additional_arrays,   only : par_matrix_computational_halo_exchange
  use mod_par_output_partition,    only : par_output_partition
  use mod_par_output_partition,    only : par_output_global_matrix
  use mod_unity_tests,             only : unity_tests_integration_rules
  use mod_unity_tests,             only : unity_tests_check_halos
  use mod_commdom_driver,          only : CNT_CPLNG,commdom_driver_init,commdom_driver_set_mesh !< 2016MAR30
  use mod_measurements,            only : measurements_set_function                             !< 2017JAN07
  use mod_mesh_type,               only : mesh_type_allocate_initialize
  use mod_mesh_type,               only : mesh_type_save_original_mesh
  use mod_mesh_type,               only : mesh_type_update_last_mesh
  use mod_mesh_multiplication,     only : mesh_multiplication_node_codes
  use mod_renumbering,             only : renumbering_temporal_locality
  use mod_messages,                only : messages_live
  use mod_element_data_base,       only : element_data_base_save
  use mod_par_bin_structure,       only : par_bin_structure
  use mod_communications,          only : PAR_BARRIER
  use mod_materials,               only : materials_on_nodes
  use mod_outfor,                  only : outfor
  use mod_output,                  only : output_domain
  use mod_par_color_communicators, only : par_color_communicators_allocate
  use mod_par_color_communicators, only : par_color_communicators
  use mod_witness,                 only : witness_point
  use mod_witness,                 only : witness_mesh_create
  use mod_finite_volume,           only : finite_volume_arrays
  use mod_renumbering_nodes,       only : renumbering_nodes

  implicit none
  real(rp) :: time1,time2
  integer(ip) :: istat

  !----------------------------------------------------------------------
  !
  ! From now on, Master does not have any mesh info :o(
  !
  !----------------------------------------------------------------------
 
  call cputim(time1)
  call messages_live('CONSTRUCT DOMAIN','START SECTION')
  !
  ! Some domain variables
  !
  call domvar(4_ip)                                       ! LNNOD, LNNOB, LGAUS
  !
  ! Renumber nodes to have own nodes
  ! (interior+own boundary) first
  !
  call renumbering_nodes()                                ! All nodal arrays read in reageo
  !
  ! Point mesh type
  !
  call mesh_type_update_last_mesh()                       ! MESHE(NDIVI) =>
  !
  ! Global numbering
  !
  call par_tools_global_numbering()
  call PAR_BARRIER()
  
  !----------------------------------------------------------------------
  !
  ! Color communicators initialization
  !
  !----------------------------------------------------------------------
  !
  ! Allocate
  !
  call par_color_communicators_allocate()
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
  call par_extgra()                                       ! R_DOM, C_DOM, NZDOM
  !
  ! Point to mesh type
  !
  meshe(ndivi) % r_dom => r_dom
  meshe(ndivi) % c_dom => c_dom
  meshe(ndivi) % nzdom =  nzdom
  !
  ! Groups
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
  !
  ! Finite volume arrays
  !
  call finite_volume_arrays(meshe(ndivi))

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
  ! Element to csr
  !
  call elecsr()                                           ! LEZDO
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
  ! Ghost geometry (halos)
  !
  call PAR_BARRIER()
  call par_ghost_geometry() 
  !
  ! Full row graph for own nodes including halo nodes     ! R_DOM_OWN, C_DOM_OWN
  !
  call messages_live('COMPUTE GRAPH FOR OWN NODES (FOR FULL ROW MATRIX FORMAT)')
  call PAR_BARRIER()
  call full_row_graph(meshe(ndivi))

  !----------------------------------------------------------------------
  !
  ! Parallel arrays
  !
  !----------------------------------------------------------------------
  !
  ! Parall communication arrays depending on zones
  !
  call par_zone_communication_arrays()                    ! PAR_COMM_COLOR_ARRAY
  !
  ! Some parallel addition arrays
  !
  if( IPARALL ) then
     call par_multiplicity_ownership(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1),COMM_NAME='COMMD')             ! Multiplicity and ownership
     call par_matrix_exchange_on_interface_nodes(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1),COMM_NAME='COMMD') ! Must be done after par_extgra
     call par_node_number_in_owner(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1),COMM_NAME='COMMD')               ! Obtains COMM % node_number_in_owner
     call par_full_row_communication_arrays(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1),COMM_NAME='COMMD')      ! Communicator for full row exchange
#if defined(PARAL_AGMG) || defined(INC_PSBLAS)
     call par_matrix_w_halos_exchange_on_interface_nodes(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1))
#else
     call par_matrix_computational_halo_exchange(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1),COMM_NAME='COMMD')
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
  call witness_point()                                    ! LEWIT, SHWIT
  call witness_mesh_create()                              ! Witness meshes

  !----------------------------------------------------------------------
  !
  ! Coupling computations
  !
  !----------------------------------------------------------------------

  call PAR_BARRIER()

#ifndef COMMDOM

  call cou_inivar(1_ip)
  call cou_define_wet_geometry()                          ! COUPLING_TYPE
  istat = 1_ip
  call cou_initialize_coupling(istat)                     ! COUPLING_TYPE
  call cou_inivar(2_ip)
  call cou_communication_arrays()
  call cou_output()

#else

  call measurements_set_function( commdom_driver_set_mesh, CNT_CPLNG, -1, -ID_KERMOD )

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

  call cputim(time2)
  cpu_start(CPU_CONSTRUCT_DOMAIN) = cpu_start(CPU_CONSTRUCT_DOMAIN) + time2 - time1

  call messages_live('CONSTRUCT DOMAIN','END SECTION')


  block
    use mod_witness
    use def_kintyp_mesh_basic
    use def_kintyp_mesh
    use def_domain
    type(mesh_type_basic) :: mesh_out
    real(rp)              :: param(4)
    param(1) = 0.5_rp
    param(2) = 0.5_rp
    param(3) = 0.5_rp
    if( ndime == 2 ) then
       param(3) = -param(1)*0.5_rp-param(2)*0.5_rp
    else
       param(4) = -param(1)*0.5_rp-param(2)*0.5_rp-param(3)*0.5_rp
    end if
    !call witness_mesh_plane(mesh_out,meshe(ndivi),param)
  end block
  
end subroutine Domain
