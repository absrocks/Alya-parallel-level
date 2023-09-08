module def_adapti

  !-----------------------------------------------------------------------
  !****f* adapti/def_adapti
  ! NAME
  !    def_adapti
  ! DESCRIPTION
  !    Heading for the ADAPTI service
  ! USED BY
  !    Almost all service subroutines
  !***
  !-----------------------------------------------------------------------
  use def_kintyp

  !------------------------------------------------------------------------
  ! Types
  !------------------------------------------------------------------------

  type eleix_ada
     !
     ! This type defines an indexed list for the new elements
     !
     integer(ip)                             :: neseg   ! # of new segments
     integer(ip)                             :: knext   ! next element
     integer(ip)                             :: kprev   ! previous element
     integer(ip)                             :: kfaie   ! segmented face
     integer(ip)                             :: knoie   ! opposed inode
     integer(ip),    dimension(:,:), pointer :: lnose   ! their lnods
     integer(ip),    dimension(:)  , pointer :: ltyse   ! their types
  end type eleix_ada

  
  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip), parameter   :: &
       lun_pdata_ada = 5901,  &    ! 
       lun_outpu_ada = 5902,  &    ! 
       lun_oumgi_ada = 5903,  &    ! 
       lun_oudom_ada = 5904,  &    ! 
       lun_ouerr_ada = 5905,  &    ! 
       lun_ougid_ada = 5906,  &    ! 
       lun_aldom_ada = 5907,  &    ! 
       lun_algeo_ada = 5908,  &    ! 
       lun_alfix_ada = 5909,  &    ! 
       lun_chpoi_ada = 5910

  !------------------------------------------------------------------------
  ! File names
  !------------------------------------------------------------------------

  character(150)           :: &
       fil_rstar_ada,         &     ! Restart file
       fil_oumgi_ada,         &     ! Mesh generator input file
       fil_oudom_ada,         &     ! 
       fil_ouerr_ada,         &     ! 
       fil_aldom_ada,         &     ! 
       fil_algeo_ada,         &     ! 
       fil_alfix_ada,         &     ! 
       fil_chpoi_ada,         &     ! 
       fil_ougid_ada                ! 

  !------------------------------------------------------------------------
  ! Reastr
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_remst_ada,         & 
       kfl_redgr_ada,         & 
       kfl_inter_ada,         & 
       kfl_chpoi_ada,         & 
       kfl_dista_ada,         & 
       kfl_prpro_ada,         & 
       kfl_ougid_ada,         & 
       kfl_oudom_ada,         & 
       kfl_livei_ada,         & 
       nimmo_ada,             &    ! Number of immersed objects
       npoii_ada(10),         &     
       nelei_ada(10),         &     
       nelto_ada,             &     
       npoto_ada,             &     
       ndimi_ada,             &     
       npnew_ada,             &     
       nenew_ada,             &     
       nelsu_ada,             &     
       nposu_ada,             &     
       nelpa_ada,             &     
       npopa_ada,             &     
       nechi_ada,             &     
       nbnew_ada,             &     
       nnodb_ada,             &     
       nnoco_ada,             &     
       nnode_ada,             &     
       nnodi_ada,             &
       npoco_ada,             &
       nelco_ada,             &
       peltb_ada,             &
       pelte_ada,             &
       npori_ada,             &
       neori_ada
       
  real(rp)                 :: &
       sizet_ada,             &             
       sizem_ada(3,10),       &
       vloca_ada(3,10),       &
       vorig_ada(3,10),       &
       vscal_ada(3,10),       &
       vroti_ada(3,10)

  integer(ip) , pointer    :: &
       lnori_ada(:,:),        &
       lsori_ada(:)  ,        &
       lnnew_ada(:,:),        &
       lnoed_ada(:,:,:),      &
       lnosu_ada(:,:),        &
       lnoad_ada(:,:),        &
       lnodi_ada(:,:,:),        &
       lfaso_ada(:,:)  ,        &
       lsofa_ada(:)    ,        &
       lnseg_ada(:)    ,        &
       lnobi_ada(:,:,:),        &
       lelim_ada(:,:),          &
       lolne_ada(:,:),          &
       gpoic_ada(:,:,:)
  real(rp)    , pointer    :: &
       coori_ada(:,:,:),      &
       conew_ada(:,:),        &
       gshac_ada(:,:,:)
  real(rp)    , pointer    :: &
       venew_ada(:,:),        &
       scnew_ada(:)

  type(eleix_ada),      pointer :: &
       lenew_ada(:)


  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------


end module def_adapti
