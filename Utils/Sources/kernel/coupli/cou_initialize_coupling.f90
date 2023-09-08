!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    cou_initialize_coupling.f90
!> @author  Guillaume Houzeaux
!> @date    03/03/2014
!> @brief   Read coupling data
!> @details Initialize coupling data
!> @} 
!-----------------------------------------------------------------------

subroutine cou_initialize_coupling()
 
  use def_kintyp,         only :  ip,rp,i1p,r2p,lg,comm_data_par
  use def_master,         only :  kfl_timin,kfl_paral
  use def_master,         only :  INOTMASTER
  use def_master,         only :  intost
  use def_domain,         only :  ndime,npoin
  use def_domain,         only :  coord
  use def_domain,         only :  nelem
  use def_domain,         only :  lesub
  use def_coupli,         only :  typ_color_coupling
  use def_coupli,         only :  mcoup
  use def_coupli,         only :  kfl_graph_cou
  use def_coupli,         only :  PROJECTION
  use def_coupli,         only :  STRESS_PROJECTION
  use def_coupli,         only :  TRANSPOSE_MIRROR
  use def_coupli,         only :  coupling_type
  use def_coupli,         only :  memor_cou
  use def_coupli,         only :  BETWEEN_SUBDOMAINS
  use def_coupli,         only :  BOUNDARY_INTERPOLATION
  use def_coupli,         only :  nboun_cou
  use def_coupli,         only :  lnodb_cou
  use def_coupli,         only :  ltypb_cou
  use def_coupli,         only :  lnnob_cou
  use def_coupli,         only :  lelbo_cou
  use def_coupli,         only :  IMPOSE_ZERO             
  use def_coupli,         only :  DONT_DO_ANYTHING        
  use def_coupli,         only :  STOP_ALYA               
  use def_coupli,         only :  STOP_ALYA_WITH_WARNINGS 
  use def_coupli,         only :  kfl_lost_wet_point_cou
  use mod_parall,         only :  color_target
  use mod_parall,         only :  color_source
  use mod_communications, only :  PAR_BARRIER
  use mod_communications, only :  PAR_SUM
  use mod_communications, only :  PAR_MAX
  use mod_communications, only :  PAR_GATHER
  use mod_communications, only :  PAR_GATHERV
  use mod_communications, only :  PAR_ALLGATHER
  use mod_communications, only :  PAR_SEND_RECEIVE_TO_ALL
  use mod_interpolation,  only :  COU_SOURCE_INTERFACE_MASS_MATRIX
  use mod_interpolation,  only :  COU_TARGET_INTERFACE_MASS_MATRIX
  use mod_interpolation,  only :  COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_couplings,      only :  COU_INIT_INTERPOLATE_POINTS_VALUES
  use mod_couplings,      only :  I_AM_IN_COUPLING
  use mod_kdtree,         only :  kdtree_construct
  use mod_memory,         only :  memory_alloca
  use mod_memory,         only :  memory_deallo
  use mod_maths,          only :  maths_heap_sort
  use mod_messages,       only :  livinf
  use mod_couplings_communications, only :  COU_GENERATE_LOCAL_TRANSMISSION_MATRICES
  use mod_couplings_communications, only :  COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES
  use mod_couplings_communications, only :  COU_PARALLELIZE_TRANSMISSION_MATRICES
  use mod_coupling_memory,          only :  COU_DEALLOCATE_GEOMETRY   
  use mod_parall,                   only :  PAR_DEALLOCATE_COMMUNICATION_ARRAY
  use mod_coupling_memory,          only :  COU_DEALLOCATE_TRANSMISSION          
  use mod_coupling_memory,          only :  COU_DEALLOCATE_SINGLE_COUPLING
  use mod_coupling_memory,          only :  COU_INITIALIZATION_SINGLE_COUPLING
    
  use mod_coupling_timer !< 2017ABR07

  implicit none
  integer(ip)                 :: icoup
  integer(ip)                 :: ipoin,kpoin,jcoup,zone_source
  integer(ip)                 :: subdomain_source,ierro,code_target
  integer(ip)                 :: pnodb,kboun,iboun,ielem
  logical(lg),    pointer     :: lesou(:)
  integer(ip),    pointer     :: lbsou(:)
  real(rp)                    :: time1,time2
  character(100), PARAMETER   :: vacal = "cou_initialize_coupling"

  nullify(lesou)
  nullify(lbsou)

  if( mcoup > 0 ) then

     call livinf(0_ip,'COUPLI: INITIALIZE COUPLING',0_ip)
 
     !-------------------------------------------------------------------
     !     
     ! Loop over couplings
     !
     !-------------------------------------------------------------------
     do icoup = 1,mcoup
        
        if( I_AM_IN_COUPLING(icoup) ) then

           code_target      = coupling_type(icoup) % code_target
           zone_source      = coupling_type(icoup) % zone_source
           subdomain_source = coupling_type(icoup) % subdomain_source
           color_target     = coupling_type(icoup) % color_target
           color_source     = coupling_type(icoup) % color_source
           jcoup            = coupling_type(icoup) % mirror_coupling 
           !
           ! Compute KDTree of the sources if necessary
           !
           if( kfl_timin == 1 ) call PAR_BARRIER('IN CURRENT COUPLING')
           call cputim(time1) 

           if( INOTMASTER ) then

              if(    coupling_type(icoup) % itype == STRESS_PROJECTION .or. &
                   & coupling_type(icoup) % itype == PROJECTION        ) then
                 !
                 ! Get mirror kdtree
                 !
                 if( jcoup /= 0 ) then
                    call kdtree_construct(&
                         nboun_cou,npoin,lnodb_cou,ltypb_cou,coord,coupling_type(icoup) % geome % kdtree,&
                         coupling_type(jcoup) % wet % lboun_wet)
                 end if

              else if( coupling_type(icoup) % itype == BOUNDARY_INTERPOLATION ) then
                 
                 !if( jcoup /= 0 ) then
                 !   call kdtree_construct(&
                 !      nboun_cou,npoin,lnodb_cou,ltypb_cou,coord,coupling_type(icoup) % geome % kdtree,&
                 !      coupling_type(jcoup) % wet % lboun_wet)
                 !else
                    !
                    ! Consider only boundaries involved as a source
                    !
                    call memory_alloca(memor_cou,'LESOU',vacal,lesou,nelem)
                    if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
                       do ielem = 1,nelem
                          if( lesub(ielem) == coupling_type(icoup) % subdomain_source ) lesou(ielem) = .true.
                       end do
                    else         
                       do ielem = 1,nelem
                          lesou(ielem) = .true.
                       end do
                    end if
                    kboun = 0
                    do iboun = 1,nboun_cou
                       pnodb = lnnob_cou(iboun)
                       ielem = lelbo_cou(iboun)
                       if( lesou(ielem) ) kboun = kboun + 1
                    end do
                    call memory_alloca(memor_cou,'LBSOU',vacal,lbsou,kboun)
                    kboun = 0
                    do iboun = 1,nboun_cou
                       pnodb = lnnob_cou(iboun)
                       ielem = lelbo_cou(iboun)
                       if( lesou(ielem) ) then
                          kboun = kboun + 1
                          lbsou(kboun) = iboun
                       end if
                    end do
                    call kdtree_construct(nboun_cou,npoin,lnodb_cou,ltypb_cou,coord,coupling_type(icoup) % geome % kdtree,lbsou)
                    call memory_deallo(memor_cou,'LBSOU',vacal,lbsou)
                    call memory_deallo(memor_cou,'LESOU',vacal,lesou)
                 !endif
              end if
           end if

           call cputim(time2) 
           coupling_type(icoup) % cputim(2) = coupling_type(icoup) % cputim(2) + time2 - time1
           !
           ! Get the CPUs in charge of my wet points
           !
           call livinf(0_ip,'COUPLI: DEFINE INTERPOLATION OF WET NODES FOR COUPLING '//trim(intost(icoup)),0_ip)
           if(.not. coupling_type(icoup) % itype == TRANSPOSE_MIRROR ) then 
              !
              call coupling_timer_set_clik()                !< 2017ABR07
              call COU_INIT_INTERPOLATE_POINTS_VALUES(coupling_type(icoup) % wet % coord_wet,color_target,&
                 color_source,coupling_type(icoup))              
              call coupling_timer_set_clak( 2_ip, 'SETMES') !< 2017ABR07
              !
              ! Check if all points have been found
              !
              ierro = 0
              if( coupling_type(icoup) % kfl_lost_wet_points == 0 ) then
                 do kpoin = 1,coupling_type(icoup) % wet % number_wet_points
                    ipoin = coupling_type(icoup) % geome % status(kpoin)
                    if( ipoin == 0 ) then
                       if( INOTMASTER .and. kfl_lost_wet_point_cou == STOP_ALYA_WITH_WARNINGS ) &
                          write(*,*) icoup,coupling_type(icoup) % wet % coord_wet(1:ndime,kpoin)
                       ierro = ierro + 1
                    end if
                 end do
              end if
           end if
           if( kfl_lost_wet_point_cou == STOP_ALYA_WITH_WARNINGS .or. kfl_lost_wet_point_cou == STOP_ALYA ) then
              call PAR_MAX(ierro,'IN CURRENT COUPLING')
              if( ierro > 0 ) call runend('COU_INITIALIZE_COUPLING: SOME WET NODES ARE LOST')
           end if
        end if

     ! call PAR_SCHEDULING(coupling_type(icoup))
     end do
     !
     ! Compute mass matrix if required
     !
     if( INOTMASTER ) then
        do icoup = 1,mcoup
           if( coupling_type(icoup) % itype == STRESS_PROJECTION ) then
              jcoup = coupling_type(icoup) % mirror_coupling 
              if( jcoup /= 0 ) then
                 call COU_SOURCE_INTERFACE_MASS_MATRIX(coupling_type(icoup),coupling_type(jcoup))
              end if
           else if( coupling_type(icoup) % itype == PROJECTION ) then
              call COU_TARGET_INTERFACE_MASS_MATRIX(coupling_type(icoup))
           end if
        end do
     end if
     
     !-------------------------------------------------------------------
     !     
     ! Generater Transmission Matrices
     !
     !-------------------------------------------------------------------
     ! 1. Generate local transmission matrices
     do icoup = 1,mcoup
        if( I_AM_IN_COUPLING(icoup) .and. (.not. coupling_type(icoup) % itype == TRANSPOSE_MIRROR ) ) then
           call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling_type(icoup))
        end if
     end do
     ! 2. Generate trans. matrices for transpose couplings
     do icoup = 1,mcoup
        if( I_AM_IN_COUPLING(icoup) .and. coupling_type(icoup) % itype == TRANSPOSE_MIRROR)  then
           call COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES(coupling_type(icoup))
        end if
     end do
     ! 3. Expand matrices for parallel executions
     do icoup = 1,mcoup
        if( I_AM_IN_COUPLING(icoup) ) then
           call COU_PARALLELIZE_TRANSMISSION_MATRICES(coupling_type(icoup))
        end if
     end do
  end if

  if(kfl_graph_cou == 1_ip .and. mcoup>0) then

     call cou_elements_graph()
     call memgeo(-41_ip)
     call memgeo(-43_ip)

  endif

!call runend('O.K.!')
end subroutine cou_initialize_coupling

!----------------------------------------------------------------------
!>
!> @author  Ricard Borrell
!> @date    13/06/2017
!> @brief   Element graph of coupled problelms 
!> @details Generates the graph (base on nodal connectivity) of the coupled 
!>          problem. Two elements of different codes are coupled if any of 
!>          their nodes are connected through the transmission matrices.  
!>          Actually this tool does not work properly if combined with mesh 
!>          multiplication, because the "internal" graph (lelel, pelel)  is 
!>          evaluated before carring out the multiplicaiton, and the "external"
!>          couplings after
!>
!>   OUTPUT (saved in case.cou.dat_graph) 
!>
!>   NELEM_COU:  #elmements of the graph representing overall coupled problem
!>   NEDGE_COU:  #edges of the graph
!>   PELEL_COU:  Pointers to the LELEL_COU array (of size NELEM_COU+1)
!>   PELEL_COU:  Adjacencies (edges) of each element
!>   LELEW_COU:  Edges weights
!>
!----------------------------------------------------------------------
subroutine cou_elements_graph()

  use def_kintyp,         only :  ip,rp,i1p,r2p,lg,comm_data_par
  use def_master,         only :  namda
  use def_master,         only :  IMASTER
  use def_master,         only :  current_code
  use def_master,         only :  leinv_loc
  use def_domain,         only :  npoin
  use def_domain,         only :  nelem
  use def_domain,         only :  pelpo
  use def_domain,         only :  lelpo
  use def_domain,         only :  pelel
  use def_domain,         only :  lelel
  use def_domain,         only :  nedge
  use def_coupli,         only :  typ_color_coupling
  use def_coupli,         only :  mcoup
  use def_coupli,         only :  coupling_type
  use def_coupli,         only :  memor_cou
  use mod_memory,         only :  memory_alloca
  use mod_memory,         only :  memory_deallo
  use mod_parall,         only :  color_target
  use mod_parall,         only :  color_source
  use mod_parall,         only :  PAR_MY_WORLD_RANK 
  use mod_parall,         only :  PAR_WORLD_SIZE
  use mod_communications, only :  PAR_SUM
  use mod_communications, only :  PAR_MAX
  use mod_communications, only :  PAR_GATHER
  use mod_communications, only :  PAR_GATHERV
  use mod_communications, only :  PAR_ALLGATHER
  use mod_communications, only :  PAR_SEND_RECEIVE_TO_ALL
  use mod_couplings,      only :  I_AM_IN_COUPLING
  use mod_maths,          only :  maths_heap_sort

  implicit none

  logical(lg)                 :: ICOUMASTER
  logical(lg)                 :: auxl
  character(150)              :: fname       !output file name
  integer(ip)                 :: u=189962    !output file unit     
  integer(ip)                 :: icoup,ipoin,irank,icode
  integer(ip)                 :: ielms,inewe,ineig,itmat
  integer(ip)                 :: icont,jcont,ia,ja
  integer(ip)                 :: ncode,nnewe,nneig
  integer(ip)                 :: nelem_world, nelem_code
  integer(ip)                 :: nedge_world, nedge_code
  integer(ip)                 :: npoin_recv,npelel,nlelel
  integer(ip)                 :: melpo
  integer(ip)                 :: cuelm
  integer(ip)                 :: sumnew
  integer(ip)                 :: auxi,auxj
  integer(ip)                 :: iaux,kaux
  integer(ip)                 :: auxs(2)
  integer(ip),pointer         :: nelco(:)          ! #elements per core
  integer(ip),pointer         :: nselc(:)          ! sum of #elements in "previous" codes
  integer(ip),pointer         :: lcora(:)          ! code associated to each rank
  real(rp),pointer            :: lnels(:,:)        ! list number of elements per point in source
  real(rp),pointer            :: lpels(:,:)        ! list elements indices per point source
  type(r2p),pointer           :: lnelt(:)          ! list number of elements per point in target
  type(r2p),pointer           :: lpelt(:)          ! list elements indices per point target
  integer(ip),pointer         :: lnewe(:,:)
  integer(4),pointer          :: lnnewe(:)         ! number of new edges generated on each rank
  integer(ip),pointer         :: lnewem(:,:)       ! list of new edged stored in COUMASTER
  type(i1p),pointer           :: extco(:)          ! external couplings
  integer(4),pointer          :: lnpel(:)          ! list of pelel sizes
  integer(4),pointer          :: lnlel(:)          ! list of lelel sizes
  integer(ip),pointer         :: pelel_buf(:)      ! buffer used by COUMASTER to receive pelels
  integer(ip),pointer         :: lelel_buf(:)      ! buffer used by COUMASTER to receive lelels  
  integer(ip),pointer         :: lelel_aux(:)
  type(typ_color_coupling)    :: coupling
  character(100), PARAMETER :: vacal = " cou_elements_graph "
  !
  !   Outputs (written to disc)
  !
  integer(ip)                 :: nelem_cou 
  integer(ip)                 :: nedge_cou
  integer(ip),pointer         :: pelel_cou(:)
  integer(ip),pointer         :: lelel_cou(:)
  integer(ip),pointer         :: lelew_cou(:)
  

  !
  !   Initializations
  !
  ICOUMASTER = .false.
  if(PAR_MY_WORLD_RANK == 0_ip) ICOUMASTER = .true.
  nullify(nelco,nselc,lcora,lnnewe,lnpel,lnlel,lnewe,lnewem,extco,pelel_cou)
  nullify(lelel_buf,pelel_buf,lelel_aux,lpels,lelel_cou,lnels,lnelt,lelew_cou)
  nullify(lpelt)
  sumnew = 0_ip
  !
  ! Evaluate: nedge_world and nelem_world
  !
  auxs(1)=nelem
  auxs(2)=nedge
  call PAR_SUM(2_ip, auxs,"IN THE WORLD","INCLUDE MASTER")
  nelem_world=auxs(1)
  nedge_world=auxs(2)
  !
  ! Evaluate: nedge_code and nelem_code
  !
  auxs(1)=nelem
  auxs(2)=nedge
  call PAR_SUM(2_ip, auxs,"IN MY CODE","INCLUDE MASTER")
  nelem_code=auxs(1)
  nedge_code=auxs(2)
  !
  ! Evaluate number of codes (ncode) elements per code (nelco(:))
  ! sum of elements in previous codes (nselc(:)) and the 
  ! code to which each rank is associated (lcora(:))
  !
  ncode = current_code
  call PAR_MAX(ncode,"IN THE WORLD")

  call memory_alloca(memor_cou,'nelco',vacal,nelco,ncode)
  call memory_alloca(memor_cou,'nselc',vacal,nselc,ncode)
  call memory_alloca(memor_cou,'lcora',vacal,lcora,PAR_WORLD_SIZE)

  if(IMASTER) then
     nelco(current_code) = nelem_code
  endif
  call PAR_SUM(nelco,"IN THE WORLD","INCLUDE MASTER")
  auxi = 0
  do icode = 1,ncode
     nselc(icode) = auxi
     auxi = auxi + nelco(icode)
  enddo
  if(auxi /= nelem_world) then
     call runend("cou_initialize_coupling: something wrong at summing up code elements")
  endif
  call PAR_ALLGATHER(current_code,lcora,1_4,"IN THE WORLD")
  !
  ! Send through the coupling commuincaiton, #elements (lnels(:)) 
  ! and elements containing each source point (lpels(:)). Those are 
  ! received in lnelt(:) and lpelt(:), respectively. Used melpo, maximum
  ! elements containing any point of the coupling
  !
  if(.not.IMASTER) then
     melpo = 0_ip 
     call memory_alloca(memor_cou,'lnels',vacal,lnels,1_ip,npoin)
     call memory_alloca(memor_cou,'lnelt',vacal,lnelt,mcoup)
     do ipoin = 1,npoin
        lnels(1_ip,ipoin) = real(pelpo(ipoin+1)-pelpo(ipoin),rp)
     enddo
     do icoup = 1,mcoup

        if( I_AM_IN_COUPLING(icoup) ) then

           coupling         = coupling_type(icoup)
           color_source     = coupling % color_source
           color_target     = coupling  % color_target

           npoin_recv       = coupling % commd % lrecv_dim
           nullify(lnelt(icoup) % a)
           call memory_alloca(memor_cou,'lnelt % a',vacal,lnelt(icoup) % a,1_ip,max(1_ip,npoin_recv))

           call PAR_SEND_RECEIVE_TO_ALL(1_ip,lnels,lnelt(icoup) % a,coupling  % commd&
              ,'ASYNCHRONOUS',coupling % commd % lsend_perm)
           melpo = max(melpo,int(maxval(lnelt(icoup) % a(1_ip,:)),ip))

        endif
     enddo
  endif

  call PAR_MAX(melpo,"IN THE WORLD")

  if(.not.IMASTER) then
     call memory_alloca(memor_cou,'lpels',vacal,lpels,melpo,npoin)
     call memory_alloca(memor_cou,'lpelt',vacal,lpelt,mcoup)
     lpels = -1_rp
     do ipoin = 1,npoin
        icont = 1_ip
        do ielms = pelpo(ipoin),pelpo(ipoin+1)-1_ip
           if(icont <= melpo) then
              lpels(icont,ipoin) = real(leinv_loc(lelpo(ielms)),rp) + real(nselc(lcora(PAR_MY_WORLD_RANK+1)),rp)
              icont = icont + 1_ip
           endif
        enddo
     enddo
     do icoup = 1,mcoup

        if( I_AM_IN_COUPLING(icoup) ) then

           coupling         = coupling_type(icoup)
           color_source     = coupling % color_source
           color_target     = coupling % color_target

           npoin_recv       = coupling % commd % lrecv_dim
           nullify(lpelt(icoup) % a)
           call memory_alloca(memor_cou,'lpelt % a',vacal,lpelt(icoup) % a,melpo,max(1_ip,npoin_recv))

           call PAR_SEND_RECEIVE_TO_ALL(mcoup,lpels,lpelt(icoup) % a,coupling % commd&
              ,'ASYNCHRONOUS',coupling % commd % lsend_perm)

        endif
     enddo
  endif
  !
  ! Count new edges (nnewe) and store them in lnewe(:,:)
  ! Note: for each new edge (i,j), the oposite edge (j,i) is also included
  !
  nnewe = 0_ip
  do icoup = 1,mcoup
     coupling         = coupling_type(icoup)
     if(associated(coupling % ltransmat_target)) then
        do ineig = 1,coupling % commd % nneig
           if(associated(coupling % ltransmat_target(ineig) % iA)) then

              do itmat = 1, size(coupling % ltransmat_target(ineig) % jA)

                 ja = coupling % ltransmat_target(ineig) % jA(itmat)
                 ia = coupling % wet % lpoin_wet( coupling % ltransmat_target(ineig) % iA(itmat))
                 nnewe = nnewe + int(lnelt(icoup) % a(1, ja) * lnels(1,ia),ip)

              enddo
           end if
        end do
     end if
  enddo

  call memory_alloca(memor_cou,'lnewe',vacal,lnewe,2_ip,2_ip*nnewe)
  inewe = 1_ip
  do icoup = 1,mcoup
     coupling         = coupling_type(icoup)
     if(associated(coupling % ltransmat_target)) then
        do ineig = 1,coupling % commd % nneig
           if(associated(coupling % ltransmat_target(ineig) % iA)) then

              do itmat = 1, size(coupling % ltransmat_target(ineig) % jA)

                 ja = coupling % ltransmat_target(ineig) % jA(itmat)
                 ia = coupling % wet % lpoin_wet( coupling % ltransmat_target(ineig) % iA(itmat))

                 do icont= 1,int(lnels(1,ia),ip)
                    do jcont = 1,int(lnelt(icoup) % a(1_ip,ja),ip) 
                       lnewe(1_ip,inewe) = int(lpels(icont,ia),ip)
                       lnewe(2_ip,inewe) = int(lpelt(icoup) % a(jcont, ja),ip)
                       inewe = inewe + 1_ip
                       lnewe(1_ip,inewe) = int(lpelt(icoup) % a(jcont, ja),ip)
                       lnewe(2_ip,inewe) = int(lpels(icont,ia),ip)
                       inewe = inewe + 1_ip
                    enddo
                 enddo

              enddo
           end if
        end do
     end if
  enddo
  if(inewe-1_ip /= nnewe*2_ip) call runend("Error in contruction of new edges of coupled graph")

  !
  ! Send all "external" couplings edges to COUMASTER
  !

  if(ICOUMASTER) call memory_alloca(memor_cou,'lnnewe',vacal,lnnewe,PAR_WORLD_SIZE)
  call PAR_GATHER(int(2*nnewe,4),lnnewe,"IN THE WORLD")

  if(ICOUMASTER)then
     nnewe = 0_ip
     nnewe = nnewe + sum(lnnewe)
     lnnewe(:) = 2_ip*lnnewe(:)
     call memory_alloca(memor_cou,'lnewem',vacal,lnewem,2_ip,nnewe)
  endif
  call PAR_GATHERV(lnewe,lnewem,lnnewe,"IN THE WORLD")

  !
  ! Store "external" couplings in extco double pointer:
  !
  !   extco[i] -> allocated in case there are couplings for element i
  !   extco[i][1] -> #edges
  !   extco[i][2] ... extco[i][exco[i][1]+1] -> edges with element i as initial point
  !

  if(ICOUMASTER) then
     !
     ! Order new edges
     !
     call maths_heap_sort(2_ip,nnewe,ivin=lnewem(1,:),ivo1=lnewem(2,:))
     !
     ! Allocate structure to estore graph of external couplings (count and allocate)  
     !
     call memory_alloca(memor_cou,'extco',vacal,extco,nelem_world)
     nneig = 0_ip
     cuelm = lnewem(1,1)
     do inewe=1,nnewe
        if(lnewem(1_ip,inewe)==cuelm) then
           nneig = nneig + 1_ip
        else
           call memory_alloca(memor_cou,'extco(cuelm) % l',vacal,extco(cuelm) % l,nneig+1)
           extco(cuelm) % l(1_ip) = 0_ip 
           cuelm = lnewem(1_ip,inewe)
           nneig = 1_ip
        endif
     enddo
     call memory_alloca(memor_cou,'extco(cuelm) % l',vacal,extco(cuelm) % l,nneig+1)
     extco(cuelm) % l(1_ip) = 0_ip 
     !
     ! Store couplings in allocated structure   
     !
     sumnew = 0_ip
     do inewe=1,nnewe
        cuelm = lnewem(1_ip,inewe)
        auxl = .false.
        do icont = 1,extco(cuelm) % l(1)
           if(extco(cuelm) % l(icont+1) == lnewem(2_ip,inewe)) then
              auxl = .true.
           endif 
        enddo
        if(.not.auxl) then
           extco(cuelm) % l(1) = extco(cuelm) % l(1) + 1_ip
           extco(cuelm) % l(extco(cuelm) % l(1)+1) = lnewem(2_ip,inewe)
           sumnew = sumnew + 1_ip
        endif
     enddo

  endif
  !
  ! Send all "pelels" to COUMASTER (this could be optimized, 
  !                                 with less comms)
  !
  if(ICOUMASTER) call memory_alloca(memor_cou,'lnpel',vacal,lnpel,PAR_WORLD_SIZE)
  call PAR_GATHER(int(size(pelel),4),lnpel,"IN THE WORLD")

  if(ICOUMASTER) then
     jcont = 1_ip
     do icont = 1,PAR_WORLD_SIZE
        if(lnpel(icont) == 1_ip) then
           lnpel(icont) = 0_ip
        else
           if(jcont /= lcora(icont)) then
              call runend("Codes disorderd across ranks. Code not ready for this")
           endif
           jcont = jcont + 1_ip
        endif
     enddo
     npelel = 0_ip
     npelel = npelel + sum(lnpel)
     call memory_alloca(memor_cou,'pelel_buf',vacal,pelel_buf,npelel)
  endif
  call PAR_GATHERV(pelel,pelel_buf,lnpel,"IN THE WORLD")
  !
  ! Send all "lelels" to COUMASTER 
  !
  if(ICOUMASTER) then
     call memory_alloca(memor_cou,'lnlel',vacal,lnlel,PAR_WORLD_SIZE)
     lnlel = 0_ip
     auxi = 0_ip
     do icont = 1,PAR_WORLD_SIZE
        if(lnpel(icont) > 0_ip) then
           auxi = auxi + lnpel(icont)
           lnlel(icont) = pelel_buf(auxi)-1_ip
        endif
     enddo
     nlelel = 0_ip
     nlelel = nlelel + sum(lnlel)
     call memory_alloca(memor_cou,'lelel_buf',vacal,lelel_buf,nlelel)
  endif 
  if(size(lelel) > 1_ip) then
     call memory_alloca(memor_cou,'lele_aux',vacal,lelel_aux,int(size(lelel),ip))
     lelel_aux(1:size(lelel)) = lelel(1:size(lelel)) + nselc(lcora(PAR_MY_WORLD_RANK+1))
  endif
  call PAR_GATHERV(lelel_aux,lelel_buf,lnlel,"IN THE WORLD")
  !
  ! Generate pelel_cou and lelel_cou
  !
  if(ICOUMASTER) then

     call memory_alloca(memor_cou,'pelel_cou',vacal,pelel_cou,npelel-ncode+1)
     call memory_alloca(memor_cou,'lelel_cou',vacal,lelel_cou,size(lelel_buf)+sumnew)
     call memory_alloca(memor_cou,'lelew_cou',vacal,lelew_cou,size(lelel_buf)+sumnew)

     pelel_cou(1) = 1_ip
     iaux = 1_ip
     kaux = 2_ip
     auxi = 0_ip
     auxj = 1_ip

     do irank = 1_ip,PAR_WORLD_SIZE
        do icont = auxi+2_ip,auxi+lnpel(irank)
           pelel_cou(kaux) = pelel_cou(kaux-1)+pelel_buf(icont)-pelel_buf(icont-1)
           do jcont = pelel_cou(kaux-1),pelel_cou(kaux)-1_ip
              lelel_cou(auxj) = lelel_buf(iaux)
              lelew_cou(auxj) = 1_ip
              iaux = iaux + 1_ip
              auxj = auxj + 1_ip
           enddo
           if(associated(extco(kaux-1) % l)) then
              pelel_cou(kaux) = pelel_cou(kaux) + extco(kaux-1_ip) % l(1_ip)
              do jcont = 2_ip,extco(kaux-1) % l(1) + 1_ip
                 lelel_cou(auxj) = extco(kaux-1) % l(jcont)
                 lelew_cou(auxj) = 2_ip
                 auxj = auxj + 1_ip
              enddo
           endif
           kaux = kaux + 1_ip              
        enddo
        auxi = auxi + lnpel(irank)
     enddo

     nelem_cou = nelem_world
     nedge_cou = size(lelel_cou)
     !
     ! Writte to disc
     !
     fname = adjustl(trim(namda))//".cou.dat_graph"
     open(u,file=fname,action="write")
     write(u,"(I6,I6,A6)") nelem_cou,nedge_cou,"001"
     do icont = 1,nelem_cou
        do jcont = pelel_cou(icont),pelel_cou(icont+1)-1_ip
           write(u,"(I6,I6)",advance='no') lelel_cou(jcont),lelew_cou(jcont)
        enddo
        write(u,*)
     enddo
     close(u)

       
  endif
  !
  ! Deallo memory
  !   
  if(associated(nelco)) call memory_deallo(memor_cou,'nelco',vacal,nelco)
  if(associated(nselc)) call memory_deallo(memor_cou,'nselc',vacal,nselc)
  if(associated(lcora)) call memory_deallo(memor_cou,'lcora',vacal,lcora)
  if(associated(lnels)) call memory_deallo(memor_cou,'lnels',vacal,lnels)
  if(associated(lpels)) call memory_deallo(memor_cou,'lpels',vacal,lpels)
  if(associated(lnelt)) call memory_deallo(memor_cou,'lnelt',vacal,lnelt)
  if(associated(lpelt)) call memory_deallo(memor_cou,'lpelt',vacal,lpelt)
  if(associated(lnewe)) call memory_deallo(memor_cou,'lnewe',vacal,lnewe)
  if(associated(lnnewe)) call memory_deallo(memor_cou,'lnnewe',vacal,lnnewe)
  if(associated(lnewem)) call memory_deallo(memor_cou,'lnewem',vacal,lnewem)
  if(associated(extco)) call memory_deallo(memor_cou,'extco',vacal,extco)
  if(associated(lnpel)) call memory_deallo(memor_cou,'lnpel',vacal,lnpel)
  if(associated(lnlel)) call memory_deallo(memor_cou,'lnlel',vacal,lnlel)
  if(associated(pelel_buf)) call memory_deallo(memor_cou,'pelel_buf',vacal,pelel_buf)
  if(associated(lelel_buf)) call memory_deallo(memor_cou,'lelel_buf',vacal,lelel_buf)
  if(associated(lelel_aux)) call memory_deallo(memor_cou,'lelel_aux',vacal,lelel_aux)
  if(associated(pelel_cou)) call memory_deallo(memor_cou,'pelel_cou',vacal,pelel_cou)
  if(associated(lelel_cou)) call memory_deallo(memor_cou,'lelel_cou',vacal,lelel_cou)
  if(associated(lelew_cou)) call memory_deallo(memor_cou,'lelew_cou',vacal,lelew_cou)

end subroutine cou_elements_graph

