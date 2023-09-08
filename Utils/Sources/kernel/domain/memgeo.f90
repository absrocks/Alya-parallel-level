subroutine memgeo(itask)
  !-----------------------------------------------------------------------
  !****f* Domain/memgeo
  ! NAME
  !    memgeo
  ! DESCRIPTION
  !    Allocate/Deallocate the geometry arrays 
  !    ITASK=1 ... Allocate memory
  !    ITASK=2 ... Deallocate memory
  ! OUTPUT
  ! USED BY
  !    reageo
  !    sengeo
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_kermod, only     :  kfl_noslw_ker
  use def_inpout
  use def_mpio, only       :  mpio_flag_geometry_export,PAR_MPIO_ON
  use mod_memchk
  use mod_memory
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: izone,nbou1,ipoin,kpoin,ifiel
  integer(ip)             :: inode,ielem,inodb,iboun
  integer(ip)             :: nsteps

  nbou1 = max(1_ip,nboun)

  select case(itask)

  case(1_ip)
     !
     ! Allocate memory of arrays read in reageo (deallocate in case = -1 )
     !
     call memory_alloca(memor_dom,'LTYPE'    ,'memgeo' , ltype     , nelem   )
     call memory_alloca(memor_dom,'LELCH'    ,'memgeo' , lelch     , nelem   )
     call memory_alloca(memor_dom,'LNODS'    ,'memgeo' , lnods     , mnode   , nelem )
     call memory_alloca(memor_dom,'LESUB'    ,'memgeo' , lesub     , nelem   )
     call memory_alloca(memor_dom,'LMATE'    ,'memgeo' , lmate     , nelem )
     !call memory_alloca(memor_dom,'LEINV_LOC','memgeo' , leinv_loc , nelem, 'IDENTITY') 
     
     call memory_alloca(memor_dom,'COORD'    ,'memgeo' , coord     , ndime   , npoin )
     call memory_alloca(memor_dom,'LNOCH'    ,'memgeo' , lnoch     , npoin   )
     call memory_alloca(memor_dom,'LMAST'    ,'memgeo' , lmast     , npoin   )
     !call memory_alloca(memor_dom,'LNINV_LOC','memgeo' , lninv_loc , npoin , 'IDENTITY')

     call memory_alloca(memor_dom,'LNODB'    ,'memgeo' , lnodb     , mnodb   , nbou1 )
     call memory_alloca(memor_dom,'LTYPB'    ,'memgeo' , ltypb     , nbou1   )
     call memory_alloca(memor_dom,'LBOCH'    ,'memgeo' , lboch     , nbou1   )
     call memory_alloca(memor_dom,'LELBO'    ,'memgeo' , lelbo     , nbou1   )
     !call memory_alloca(memor_dom,'LBINV_LOC','memgeo' , lbinv_loc , nbou1 , 'IDENTITY')
     
     if( ISEQUEN ) then
        call memory_alloca(memor_dom,'LNINV_LOC','memgeo',lninv_loc,npoin,'IDENTITY')
        call memory_alloca(memor_dom,'LEINV_LOC','memgeo',leinv_loc,nelem,'IDENTITY')
        call memory_alloca(memor_dom,'LBINV_LOC','memgeo',lbinv_loc,nboun,'IDENTITY')
     else
        !
        ! In parallel, these arrays are computed in the associated modules
        !
     end if

     
     if( nelem > 0 ) lesub(1:nelem) = 1
     if( nelem > 0 ) lmate(1:nelem) = 1

  case(-1_ip)
     !
     ! Deallocate memory of arrays read in reageo allocate in case = 1
     !
     call memory_deallo(memor_dom,'LTYPE'    ,'memgeo' , ltype )
     call memory_deallo(memor_dom,'LELCH'    ,'memgeo' , lelch )
     call memory_deallo(memor_dom,'LNODS'    ,'memgeo' , lnods )
     call memory_deallo(memor_dom,'LESUB'    ,'memgeo' , lesub )
     call memory_deallo(memor_dom,'LMATE'    ,'memgeo' , lmate )
     !call memory_deallo(memor_dom,'LEINV_LOC','memgeo' , leinv_loc )

     call memory_deallo(memor_dom,'COORD'    ,'memgeo' , coord )
     call memory_deallo(memor_dom,'LNOCH'    ,'memgeo' , lnoch )
     call memory_deallo(memor_dom,'LMAST'    ,'memgeo' , lmast )
     !call memory_deallo(memor_dom,'LNINV_LOC','memgeo' , lninv_loc )

     call memory_deallo(memor_dom,'LNODB'    ,'memgeo' , lnodb )
     call memory_deallo(memor_dom,'LBOEL'    ,'memgeo' , lboel )
     call memory_deallo(memor_dom,'LTYPB'    ,'memgeo' , ltypb )
     call memory_deallo(memor_dom,'LBOCH'    ,'memgeo' , lboch )
     call memory_deallo(memor_dom,'LELBO'    ,'memgeo' , lelbo )
     !call memory_deallo(memor_dom,'LBINV_LOC','memgeo' , lbinv_loc )
     
     if( ISEQUEN ) then
        call memory_deallo(memor_dom,'LNINV_LOC','memgeo',lninv_loc)
        call memory_deallo(memor_dom,'LEINV_LOC','memgeo',leinv_loc)
        call memory_deallo(memor_dom,'LBINV_LOC','memgeo',lbinv_loc)
     else
        !
        !
        !
     end if
     
  case(-2_ip)
     !
     ! Deallocate memory: Mesh arrays
     !
     call memory_deallo(memor_dom,'LTYPE'    ,'memgeo' , ltype )
     call memory_deallo(memor_dom,'LNNOD'    ,'memgeo' , lnnod )
     call memory_deallo(memor_dom,'LELCH'    ,'memgeo' , lelch )
     call memory_deallo(memor_dom,'LNODS'    ,'memgeo' , lnods )
     call memory_deallo(memor_dom,'LESUB'    ,'memgeo' , lesub )
     call memory_deallo(memor_dom,'LMATE'    ,'memgeo' , lmate )
     call memory_deallo(memor_dom,'LEINV_LOC','memgeo' , leinv_loc )

     call memory_deallo(memor_dom,'COORD'    ,'memgeo' , coord )
     call memory_deallo(memor_dom,'LNOCH'    ,'memgeo' , lnoch )
     call memory_deallo(memor_dom,'LMAST'    ,'memgeo' , lmast )
     call memory_deallo(memor_dom,'LNINV_LOC','memgeo' , lninv_loc )

     call memory_deallo(memor_dom,'LNODB'    ,'memgeo' , lnodb )
     call memory_deallo(memor_dom,'LTYPB'    ,'memgeo' , ltypb )
     call memory_deallo(memor_dom,'LBOCH'    ,'memgeo' , lboch )
     call memory_deallo(memor_dom,'LELBO'    ,'memgeo' , lelbo )
     call memory_deallo(memor_dom,'LBINV_LOC','memgeo' , lbinv_loc )
     !
     ! Dealloctae R_DOM and C_DOM to minimum size (used in deflated CG)
     !
     call memory_deallo(memor_dom,'R_DOM','memgeo' , r_dom )
     call memory_deallo(memor_dom,'C_DOM','memgeo' , c_dom )
     call memory_alloca(memor_dom,'R_DOM','memgeo' , r_dom , 1_ip )
     call memory_alloca(memor_dom,'C_DOM','memgeo' , c_dom , 1_ip )
     !
     ! Deallocate LELPO and PELPO (used in parall)
     !
     call memory_deallo(memor_dom,'PELPO','memgeo' , pelpo )
     call memory_deallo(memor_dom,'LELPO','memgeo' , lelpo )
     !
     ! Deallocate memory
     !
     call memose(10_ip)

  case(3_ip)
     !
     ! Deallocate memory: Mesh arrays
     !
     call runend('QUE PASA AQUI')
     
  case(4_ip)

     call memory_deallo(memor_dom,'COORD','memgeo' , coord )

  case(5_ip)
     !
     ! LBOEL
     !
     call memory_alloca(memor_dom,'LBOEL'    ,'memgeo' , lboel     , mnodb   , nbou1 )

  case(-5_ip)
     !
     ! LBOEL
     !
     call memory_deallo(memor_dom,'LBOEL'    ,'memgeo' , lboel   )

  case(14_ip)
     !
     ! Allocate SKCOS
     !
     call memory_alloca(memor_dom,'SKCOS','memgeo' , skcos , ndime , ndime , nbopo )

  case(15_ip)
     !
     ! Allocate EXNOR
     !
     call memory_alloca(memor_dom,'EXNOR','memgeo' , exnor , ndime , ndime , nbopo )

  case(-15_ip)
     !
     ! Allocate EXNOR
     !
     call memory_deallo(memor_dom,'EXNOR','memgeo' , exnor )

  case( 19_ip)
     !
     ! Allocate  boundary graph: R_BOU and C_BOU
     !
     if( kfl_crbou == 0 ) then
        call memory_alloca(memor_dom,'R_BOU','memgeo' , r_bou , npoin+1 )
        call memory_alloca(memor_dom,'C_BOU','memgeo' , c_bou , nzbou   )
        kfl_crbou = 1
     end if

  case(-19_ip)
     !
     ! Deallocate boundary graph: R_BOU and C_BOU
     !
     if( kfl_crbou == 1 ) then
        call memory_deallo(memor_dom,'R_BOU','memgeo' , r_bou )
        call memory_deallo(memor_dom,'C_BOU','memgeo' , c_bou )
        kfl_crbou = 0
     end if

  case(24_ip)
     !
     ! Arrays for variable wall distance on boundaries and boundary nodes
     !
     call memory_alloca(memor_dom,'YWALB','memgeo' , ywalb , nboun )
     call memory_alloca(memor_dom,'YWALP','memgeo' , ywalp , nbopo )
     if ( kfl_noslw_ker /= 0_ip) call memory_alloca(memor_dom,'YWALE','memgeo' , ywale , nelem )

  case(-24_ip)
     !
     ! Arrays for variable wall distance
     !
     call memory_deallo(memor_dom,'YWALB','memgeo' , ywalb )
     call memory_deallo(memor_dom,'YWALP','memgeo' , ywalp )
     if ( kfl_noslw_ker /= 0_ip) call memory_deallo(memor_dom,'YWALE','memgeo' , ywale )

  case(-25_ip)
     !
     ! Deallocate boundary
     !
     call memory_deallo(memor_dom,'LNODB','memgeo' , lnodb )
     call memory_deallo(memor_dom,'LBOEL','memgeo' , lboel )
     call memory_deallo(memor_dom,'LTYPB','memgeo' , ltypb )
     call memory_deallo(memor_dom,'LBOCH','memgeo' , lboch )
     call memory_deallo(memor_dom,'LELBO','memgeo' , lelbo )

  case( 25_ip)
     !
     ! Allocate boundary
     !
     call memory_alloca(memor_dom,'LNODB','memgeo' , lnodb , mnodb , nbou1 )
     call memory_alloca(memor_dom,'LBOEL','memgeo' , lboel , mnodb , nbou1 )
     call memory_alloca(memor_dom,'LTYPB','memgeo' , ltypb , nbou1 )
     call memory_alloca(memor_dom,'LBOCH','memgeo' , lboch , nbou1 )
     call memory_alloca(memor_dom,'LELBO','memgeo' , lelbo , nbou1 )

  case ( 27_ip )
     !
     ! Groups for deflated CG
     !
     call memory_alloca(memor_dom,'LGROU_DOM','memgeo' , lgrou_dom , npoin )

  case ( -27_ip )
     !
     ! Groups for deflated CG
     !
     call memory_deallo(memor_dom,'LGROU_DOM','memgeo' , lgrou_dom )

  case ( 28_ip )
     !
     ! Fields: allocate
     !
     call memory_alloca(memor_dom,'XFIEL'      ,'memgeo' , xfiel       , nfiel )

  case ( -28_ip )
     !
     ! Fields: deallocate
     !
     call memory_deallo(memor_dom,'XFIEL'      ,'memgeo' , xfiel     )

  case ( 29_ip )
     !
     ! Fields values: allocate
     !../../Sources/kernel/mpio/def_mpio.f90
     ! Master should not allocate these.
     !if (INOTMASTER) then
         if ( (kfl_field(6,igene) == 1) .AND. (mpio_flag_geometry_export /= PAR_MPIO_ON) ) then !ondemand 
            nsteps = nsteps_fiel_ondemand
         else
            nsteps = kfl_field(4,igene)
         end if

         if(      kfl_field(2,igene) == NELEM_TYPE ) then
            call memory_alloca(memor_dom,'XFIEL%A','memgeo' , xfiel(igene) % a , kfl_field(1,igene) , nelem , nsteps )
         else if( kfl_field(2,igene) == NPOIN_TYPE ) then                                                                        
            call memory_alloca(memor_dom,'XFIEL%A','memgeo' , xfiel(igene) % a , kfl_field(1,igene) , npoin , nsteps )
         else if( kfl_field(2,igene) == NBOUN_TYPE ) then                                                                        
            call memory_alloca(memor_dom,'XFIEL%A','memgeo' , xfiel(igene) % a , kfl_field(1,igene) , nbou1 , nsteps )
         end if
     !end if
  case ( -29_ip )
     !
     ! Fields values: deallocate
     !
     call memory_deallo(memor_dom,'XFIEL%A','memgeo'      , xfiel(igene) % a )

  case ( 30_ip )
     !
     ! Fields values: allocate
     !
     if(      kfl_field(2,igene) == NELEM_TYPE ) then
        call memory_alloca(memor_dom,'XFIEL%A','memgeo' , xfiel(igene) % a , kfl_field(1,igene) , nelem , kfl_field(4,igene) )
     else if( kfl_field(2,igene) == NPOIN_TYPE ) then                                                                        
        call memory_alloca(memor_dom,'XFIEL%A','memgeo' , xfiel(igene) % a , kfl_field(1,igene) , npoin , kfl_field(4,igene) )
     else if( kfl_field(2,igene) == NBOUN_TYPE ) then                                                                        
        call memory_alloca(memor_dom,'XFIEL%A','memgeo' , xfiel(igene) % a , kfl_field(1,igene) , nbou1 , kfl_field(4,igene) )
     end if

  case ( -30_ip )
     !
     ! Fields values: deallocate
     !
     call memory_deallo(memor_dom,'XFIEL%A','memgeo' , xfiel(igene) % a )

  case (  32_ip )
     !
     ! VMASS
     !
     call memory_alloca(memor_dom,'VMASS','memgeo' , vmass , npoin )

  case ( -32_ip )
     !
     ! VMASS
     !
     call memory_deallo(memor_dom,'VMASS','memgeo' , vmass )

  case (  33_ip )
     !
     ! VMASC
     !
     call memory_alloca(memor_dom,'VMASC','memgeo' , vmasc , npoin )

  case ( -33_ip )
     !
     ! VMASC
     !
     call memory_deallo(memor_dom,'VMASC','memgeo' , vmasc )

  case( 34_ip )
     !
     ! Allocate materials 
     !
     call memory_alloca(memor_dom,'LMATE','memgeo' , lmate , nelem )

  case ( -34_ip )
     !
     ! Deallocate materials LMATE
     !
     call memory_deallo(memor_dom,'LMATE','memgeo' , lmate )

  case( 35_ip )
     !
     ! Allocate materials: NODEMAT, LMATN and NMATN
     !
     call memory_alloca(memor_dom,'LMATN'  ,'memgeo' , lmatn , nmate )
     call memory_alloca(memor_dom,'NMATN'  ,'memgeo' , nmatn , nmate )
     call memory_alloca(memor_dom,'NODEMAT','MEMGEO' , nodemat,npoin )
     if( npoin > 0 ) nodemat = -1_ip


  case ( -35_ip )
     !
     ! Deallocate materials: LMATN and NMATN
     !
     call memory_deallo(memor_dom,'LMATN','memgeo' , lmatn )
     call memory_deallo(memor_dom,'NMATN','memgeo' , nmatn )

  case( 36_ip)
     !
     ! LMATN(IMATE) % L
     !
     call memory_alloca(memor_dom,'LMATN(IMATE)%L','memgeo' , lmatn(igene)%l , igen2 )

  case(-36_ip)
     !
     ! LMATN(IMATE) % L
     !
     call memory_deallo(memor_dom,'LMATN(IMATE)%L','memgeo' , lmatn(igene)%l )

  case( 40_ip)
     !
     ! LFCNT and LNCNT: Allocate lfcnt
     !
     call memory_alloca(memor_dom,'LFCNT','memgeo' , lfcnt , 5_ip , necnt )
     call memory_alloca(memor_dom,'LNCNT','memgeo' , lncnt , 2_ip , nncnt )

  case(-40_ip)
     !
     ! LFCNT and LNCNT: Dellocate lfcnt
     !
     call memory_deallo(memor_dom,'LFCNT','memgeo' , lfcnt )
     call memory_deallo(memor_dom,'LNCNT','memgeo' , lncnt )

  case( 41_ip)
     !
     ! PELEL
     !
     kfl_pelel = 1
     call memory_alloca(memor_dom,'PELEL','memgeo' , pelel , nelem+1_ip )

  case(-41_ip)
     !
     ! PELEL
     !
     kfl_pelel = 0
     call memory_deallo(memor_dom,'PELEL','memgeo' , pelel )

  case( 43_ip)
     !
     ! LELEL
     !
     call memory_alloca(memor_dom,'LELEL','memgeo' , lelel , nedge )

  case(-43_ip)
     !
     ! LELEL
     !
     call memory_deallo(memor_dom,'LELEL','memgeo' , lelel )

  case( 53_ip)
     !
     ! LPERI: periodicity
     !
     call memory_alloca(memor_dom,'LPERI','memgeo' , lperi , 2_ip , nperi )

  case(-53_ip)
     !
     ! LPERI: periodicity
     !
     call memory_deallo(memor_dom,'LPERI','memgeo' , lperi )

  case( 54_ip)
     !
     ! LEZDO and LBZDO
     !
     call memory_alloca(memor_dom,'LEZDO','memgeo' , lezdo , mnode , mnode , nelem )
     call memory_alloca(memor_dom,'LBZDO','memgeo' , lbzdo , mnodb , mnodb , nbou1 )

  case(-54_ip)
     !
     ! LEZDO and LBZDO
     !
     call memory_deallo(memor_dom,'LEZDO','memgeo' , lezdo )
     call memory_deallo(memor_dom,'LBZDO','memgeo' , lbzdo )

  case ( 59_ip )
     !
     ! Which zones and subdomain I'm in   
     !     
     call memory_alloca(memor_dom,'I_AM_IN_ZONE','memgeo',I_AM_IN_ZONE,nzone+1_ip,'INITIALIZE',lboun=0_ip)
     call memory_alloca(memor_dom,'I_AM_IN_SUBD','memgeo',I_AM_IN_SUBD,nsubd+1_ip,'INITIALIZE',lboun=0_ip)

  case (-59_ip )
     !
     ! Which zones and subdomain I'm in   
     !     
     call memory_deallo(memor_dom,'I_AM_IN_ZONE','memgeo',I_AM_IN_ZONE)
     call memory_deallo(memor_dom,'I_AM_IN_SUBD','memgeo',I_AM_IN_SUBD)

  case ( 60_ip )
     !
     ! List of boundary nodes
     !     
     call memory_alloca(memor_dom,'LBONO','memgeo',lbono,nbono)

  case (-60_ip )
     !
     ! List of boundary nodes
     !     
     call memory_deallo(memor_dom,'LBONO','memgeo',lbono)

  case ( 61_ip )
     !
     ! Number of nodes per boundary
     !
     call memory_alloca(memor_dom,'LNNOB','memgeo' , lnnob ,  nboun )

  case (-61_ip )
     !
     ! Number of nodes per boundary
     !
     call memory_deallo(memor_dom,'LNNOB','memgeo' , lnnob )

  case( 62_ip)
     !
     ! Allocate list of masters
     !
     call memory_alloca(memor_dom,'LMAST','memgeo' , lmast , npoin )

  case(-62_ip)
     !
     ! Deallocate list of masters
     !
     call memory_deallo(memor_dom,'LMAST','memgeo' , lmast )

  case (  63_ip )
     !
     ! lumma - for use in dual time step preconditioner
     !
     call memory_alloca(memor_dom,'LUMMA','memgeo' , lumma , npoin )      !nzdom * ndime * ndime    ! If I want full matrix

  case ( -63_ip )
     !
     ! lumma
     !
     call memory_deallo(memor_dom,'LUMMA','memgeo' , lumma )

  case ( 64_ip )
     !
     ! Number of Gauss points per elements
     !
     call memory_alloca(memor_dom,'LGAUS','memgeo' , lgaus ,  nelem )

  case (-64_ip )
     !
     ! Number of Gauss points per elements
     !
     call memory_deallo(memor_dom,'LGAUS','memgeo' , lgaus )

  case ( 65_ip )
     !
     ! lnnod
     !
     call memory_alloca(memor_dom,'LNNOD','memgeo' , lnnod ,  nelem   )
     
  case (-65_ip )
     !
     ! lnnod
     !
     call memory_deallo(memor_dom,'LNNOD','memgeo' , lnnod )

  case ( 66_ip )
     !
     ! TIME FIELD
     !
     call memory_alloca(memor_dom,'TIME_FIELD%A','memgeo' , time_field(igene) % a , kfl_field(4,igene) )

  case (-66_ip )
     !
     ! TIME FIELD
     !
     call memory_deallo(memor_dom,'TIME_FIELD%A','memgeo' , time_field(igene) % a )
    
  CASE DEFAULT
     call runend("Undefined task in memgeo()",itask)
  end select

end subroutine memgeo
