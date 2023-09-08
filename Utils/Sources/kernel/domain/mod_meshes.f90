!-----------------------------------------------------------------------
!
!> @defgroup Meshes_Toolbox
!> Toolbox for meshes manipulations, like creating submeshes,
!> surface mesh, etc.
!> @{
!> @name    ToolBox for meshes operations
!> @file    mod_meshes.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for meshes operations
!> @details ToolBox for meshes operations
!
!-----------------------------------------------------------------------

module mod_meshes

  use def_kintyp, only : ip,rp,lg,i1p,i1pp
  use def_elmtyp, only : ELFEM
  use def_elmtyp, only : BAR02,TRI03
  use def_master, only : INOTMASTER,kfl_paral,leinv_loc
  use def_domain, only : memor_dom
  use def_domain, only : mesh_type
  use def_domain, only : mnode
  use mod_graphs, only : graphs_elepoi
  use mod_graphs, only : graphs_dealep
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_std
  implicit none
  private
  real(rp)   :: epsil = epsilon(1.0_rp)
  
  !interface meshes_boundary_mesh
  !   module procedure meshes_boundary_mesh_type,&
  !        &           meshes_boundary_mesh_arrays
  !end interface meshes_boundary_mesh

  public :: meshes_submesh
  public :: meshes_surface_from_nodal_array
  public :: meshes_surface_from_nodal_array_deallocate
  public :: meshes_list_boundary_nodes
  public :: meshes_list_boundary_elements

contains 

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute a submesh
  !> @details Compute a submesh given a permutation array
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------
 
  subroutine meshes_submesh(                                                    &
       & ndime,    npoin    ,nelem    ,coord    ,lnods    ,ltype    ,lnnod,     &
       &           npoin_sub,nelem_sub,coord_sub,lnods_sub,ltype_sub,lnnod_sub, &
       &           lperm   )

    integer(ip),          intent(in)  :: ndime
    integer(ip),          intent(in)  :: nelem
    integer(ip),          intent(in)  :: npoin
    real(rp),    pointer, intent(in)  :: coord(:,:)
    integer(ip), pointer, intent(in)  :: lnods(:,:)
    integer(ip), pointer, intent(in)  :: ltype(:)
    integer(ip), pointer, intent(in)  :: lnnod(:)
    integer(ip),          intent(out) :: nelem_sub
    integer(ip),          intent(out) :: npoin_sub
    real(rp),    pointer, intent(inout) :: coord_sub(:,:)
    integer(ip), pointer, intent(inout) :: lnods_sub(:,:)
    integer(ip), pointer, intent(inout) :: ltype_sub(:)
    integer(ip), pointer, intent(inout) :: lnnod_sub(:)
    integer(ip), pointer, intent(inout) :: lperm(:)
    integer(ip)                       :: ielem,inode,ipoin
    integer(ip)                       :: pelty,pnode,mnode_loc
    integer(ip)                       :: kelem,ipoin_sub
    integer(ip), pointer              :: lpoin(:)

    mnode_loc = size(lnods,1,KIND=ip)
    npoin_sub = 0
    if( .not. associated(lperm) ) then
       nelem_sub = 0
    else
       nelem_sub = size(lperm,KIND=ip)
    end if

    if( nelem_sub /= 0 ) then

       call memory_alloca(memor_dom,'LNODS_SUB','meshes_submesh',lnods_sub,mnode_loc,nelem_sub)
       call memory_alloca(memor_dom,'LTYPE_SUB','meshes_submesh',ltype_sub,          nelem_sub)
       call memory_alloca(memor_dom,'LNNOD_SUB','meshes_submesh',lnnod_sub,          nelem_sub)

       allocate(lpoin(npoin))
       do ipoin = 1,npoin
          lpoin(ipoin) = 0
       end do
       !
       ! Renumber nodes
       !
       do kelem = 1,nelem_sub
          ielem = lperm(kelem)
          pelty = abs(ltype(ielem))
          pnode = lnnod(ielem)
          ltype_sub(kelem) = pelty
          lnnod_sub(kelem) = pnode
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             if( lpoin(ipoin) == 0 ) then
                npoin_sub = npoin_sub + 1
                lpoin(ipoin) = npoin_sub
             end if
          end do
       end do
       !
       ! Copy connectivity and coordinates
       !
       call memory_alloca(memor_dom,'COORD_SUB','meshes_submesh',coord_sub,ndime,npoin_sub)

       do kelem = 1,nelem_sub
          ielem = lperm(kelem)
          pelty = abs(ltype(ielem))
          pnode = lnnod(ielem)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             lnods_sub(inode,kelem) = lpoin(ipoin)
          end do
       end do
       do ipoin = 1,npoin
          ipoin_sub = lpoin(ipoin)
          if( ipoin_sub /= 0 ) &
               coord_sub(1:ndime,ipoin_sub) = coord(1:ndime,ipoin) 
       end do

       deallocate(lpoin)

    else

       npoin_sub = 0
       nelem_sub = 0

    end if

  end subroutine meshes_submesh

  !-----------------------------------------------------------------------
  !
  !> @brief   Construct a surface mesh
  !> @details Construct a surface mesh from a nodal array xarra
  !> @date    29/09/2015
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine meshes_surface_from_nodal_array_deallocate(&
       lnodb_sur,coord_sur,ltypb_sur,lelem_sur)
    implicit none
    integer(ip),     pointer,           intent(inout) :: lnodb_sur(:,:)
    real(rp),        pointer,           intent(inout) :: coord_sur(:,:)
    integer(ip),     pointer, optional, intent(inout) :: ltypb_sur(:)
    integer(ip),     pointer, optional, intent(inout) :: lelem_sur(:)

    if( INOTMASTER ) then
       call memory_deallo(memor_dom,'LNODB_SUR','meshes_surface_from_nodal_array',lnodb_sur)
       call memory_deallo(memor_dom,'COORD_SUR','meshes_surface_from_nodal_array',coord_sur)
       if( present(ltypb_sur) ) call memory_deallo(memor_dom,'LTYPB_SUR','meshes_surface_from_nodal_array',ltypb_sur)
       if( present(lelem_sur) ) call memory_deallo(memor_dom,'LELEM_SUR','meshes_surface_from_nodal_array',lelem_sur)
    end if

  end subroutine meshes_surface_from_nodal_array_deallocate

  !-----------------------------------------------------------------------
  !
  !> @brief   Construct a surface mesh
  !> @details Construct a surface mesh from a nodal array xarra
  !> @date    29/09/2015
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine meshes_surface_from_nodal_array(&
       xarra,meshe,npoin_sur,nboun_sur,lnodb_sur,coord_sur,ltypb_sur,lelem_sur)
    implicit none
    real(rp),        pointer,           intent(in)  :: xarra(:,:)                 !< 
    type(mesh_type),                    intent(in)  :: meshe                      !< Mesh type
    integer(ip),                        intent(out) :: npoin_sur
    integer(ip),                        intent(out) :: nboun_sur
    integer(ip),     pointer,           intent(inout) :: lnodb_sur(:,:)
    real(rp),        pointer,           intent(inout) :: coord_sur(:,:)
    integer(ip),     pointer, optional, intent(inout) :: ltypb_sur(:)
    integer(ip),     pointer, optional, intent(inout) :: lelem_sur(:)
    integer(ip)                                     :: ielem,idime,icomp          ! Indices and dimensions
    integer(ip)                                     :: pelty,pnode,p1,p2,p3
    integer(ip)                                     :: inode,jnode,ipoin
    integer(ip)                                     :: compt,compli,compg,compi
    integer(ip)                                     :: ii,ij,inod1,inod2          ! element crossed by interface = 1 
    integer(ip)                                     :: signn,signp,sigtn,sigtp    ! sign inside the element 
    integer(ip)                                     :: tetra(4_ip,6_ip)
    integer(ip)                                     :: tetrp(4_ip,3_ip)
    integer(ip)                                     :: tetpy(4_ip,2_ip)
    integer(ip)                                     :: ipasq,ipasp,ipapy,itype
    integer(ip)                                     :: pnodb
    real(rp)                                        :: elcod(3,mnode)
    real(rp)                                        :: elarr(mnode),inter(3,4)
    real(rp)                                        :: l1,lp,x1,y1,x2,y2,z1,z2 
    real(rp)                                        :: gpgrl(3),gpcar(3,mnode)
    real(rp)                                        :: xjaci(9),xjacm(9),gpdet
    real(rp)                                        :: vec(3,3),coord_xxx(3)
    real(rp)                                        :: xx,nx,ny,nz
    real(rp)                                        :: norma(3),vecto(3)
    !
    ! Mesh arrays required
    !
    integer(ip)                                     :: nelem
    integer(ip)                                     :: npoin
    integer(ip)                                     :: ndime
    integer(ip),      pointer                       :: lnods(:,:)
    integer(ip),      pointer                       :: lnnod(:)
    integer(ip),      pointer                       :: lelch(:)
    real(rp),         pointer                       :: coord(:,:)

    if( INOTMASTER ) then
       !
       ! Mesh arrays
       !
       nelem =  meshe % nelem 
       npoin =  meshe % npoin 
       ndime =  meshe % ndime 
       lnods => meshe % lnods
       lnnod => meshe % lnnod
       lelch => meshe % lelch
       coord => meshe % coord
       compt =  0_ip

       if( ndime == 2 ) then
          !
          ! Loop over the elements to construct discrete interface 
          ! collection of segments (2d case)
          !
          pnodb = 2
          do ielem = 1,nelem
             if( lelch(ielem) == ELFEM ) then
                pnode = lnnod(ielem)
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   elarr(inode) = xarra(ipoin,1)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do

                compli = 0_ip
                do inode = 1,pnode
                   if( inode == pnode ) then
                      jnode = 1 
                   else
                      jnode = inode + 1
                   endif
                   !
                   ! Count interface points 
                   !
                   if( elarr(inode)*elarr(jnode) < 0.0_rp ) then
                      compli = compli + 1
                      if( compli == 3 ) compli = 1         ! Count interface points 
                      if( compli == 1 ) compt = compt + 1 ! Count interface segments
                   endif
                end do
             end if
          end do
          !
          ! Registering of discrete interface vectors size
          !
          nboun_sur = compt
          npoin_sur = ndime * nboun_sur
          compt     = 0_ip
          compg     = 0_ip
          !
          ! Memory allocation for discrete interface vectors
          ! connectivity (lnodb_sur) and points coordinates (coord_sur)
          !
          call memory_deallo(memor_dom,'LNODB_SUR','meshes_surface_from_nodal_array',lnodb_sur)
          call memory_deallo(memor_dom,'COORD_SUR','meshes_surface_from_nodal_array',coord_sur)
          if( present(ltypb_sur) ) call memory_deallo(memor_dom,'LTYPB_SUR','meshes_surface_from_nodal_array',ltypb_sur)
          if( present(lelem_sur) ) call memory_deallo(memor_dom,'LELEM_SUR','meshes_surface_from_nodal_array',lelem_sur)

          call memory_alloca(memor_dom,'LNODB_SUR','meshes_surface_from_nodal_array',lnodb_sur,pnodb,nboun_sur)
          call memory_alloca(memor_dom,'COORD_SUR','meshes_surface_from_nodal_array',coord_sur,ndime,npoin_sur)
          if( present(ltypb_sur) ) call memory_alloca(memor_dom,'LTYPB_SUR','meshes_surface_from_nodal_array',ltypb_sur,nboun_sur)
          if( present(lelem_sur) ) call memory_alloca(memor_dom,'LELEM_SUR','meshes_surface_from_nodal_array',lelem_sur,nboun_sur)
          !
          ! filling of discrete interface vectors 
          ! connectivity (lnodb_sur) and points coordinates (coord_sur)
          !
          do ielem=1,nelem
             if( lelch(ielem) == ELFEM ) then
                pnode = lnnod(ielem)
                do inode = 1,pnode
                   ipoin                = lnods(inode,ielem)
                   elarr(inode)         = xarra(ipoin,1)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do

                compli = 0_ip

                do inode = 1,pnode
                   if( inode == pnode ) then
                      jnode = 1 
                   else
                      jnode = inode + 1
                   endif
                   if( elarr(inode)*elarr(jnode) < 0.0_rp ) then
                      compg = compg + 1
                      compli = compli + 1
                      if( compli == 3 ) compli = 1
                      if( compli == 1 ) compt = compt + 1
                      !
                      ! Compute the intersection of the elements with the surface 
                      !
                      l1 = abs(elarr(inode))
                      lp = abs(elarr(inode)-elarr(jnode))
                      x1 = elcod(1,inode)
                      x2 = elcod(1,jnode)
                      y1 = elcod(2,inode)
                      y2 = elcod(2,jnode)
                      lnodb_sur(compli,compt) = compg
                      coord_sur(1,compg)     =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                      coord_sur(2,compg)     =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp
                      if( present(ltypb_sur) ) ltypb_sur(compt) = BAR02
                      if( present(lelem_sur) ) lelem_sur(compt) = ielem
                   end if
                end do
                !
                ! Orientate surface
                !
                if( compli == 2 ) then
                   x1       = coord_sur(1,compg-1)
                   y1       = coord_sur(2,compg-1)
                   x2       = coord_sur(1,compg)
                   y2       = coord_sur(2,compg)
                   norma(1) = y1-y2
                   norma(2) = x2-x1            
                   inode    = 1
                   vecto(1) = (elcod(1,inode)-0.5_rp*(x1+x2))
                   vecto(2) = (elcod(2,inode)-0.5_rp*(y1+y2))
                   l1       = norma(1)*vecto(1)+norma(2)*vecto(2)
                   if( l1 > 0.0_rp ) then
                      if( elarr(inode) >= 0.0_rp ) then
                         lnodb_sur(1,compt) = compg
                         lnodb_sur(2,compt) = compg-1
                      end if
                   else                      
                      if( elarr(inode) < 0.0_rp ) then
                         lnodb_sur(1,compt) = compg
                         lnodb_sur(2,compt) = compg-1
                      end if
                   end if
                end if

             end if
          end do

       else

          call runend('NOT CODED')

       end if

    end if

  end subroutine meshes_surface_from_nodal_array
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-31
  !> @brief   List of boundary nodes
  !> @details List of boundary noed of a mesh. No parallel exchange!
  !> 
  !-----------------------------------------------------------------------

  subroutine meshes_list_boundary_nodes(&
       nelem,npoin,lnods,ltype,&
       number_boundary_nodes,&
       list_boundary_nodes)

    use mod_maths,  only : maths_heap_sort
    use mod_elmgeo, only : element_type

    integer(ip),          intent(in)    :: nelem
    integer(ip),          intent(in)    :: npoin
    integer(ip), pointer, intent(in)    :: lnods(:,:)
    integer(ip), pointer, intent(in)    :: ltype(:)
    integer(ip),          intent(out)   :: number_boundary_nodes
    integer(ip), pointer, intent(inout) :: list_boundary_nodes(:)

    integer(ip)                         :: ielty,ielem,iface,inodf,ilist,ifacg
    integer(ip)                         :: inode,jelem,jface,jelty,ipoin,pnodf
    integer(ip)                         :: pepoi,ielpo,knode,mface,mnodf,ii,mnode_loc
    logical(lg)                         :: equal_faces  

    integer(ip)                         :: mepoi_loc
    integer(ip), pointer                :: lelpo_loc(:)
    integer(ip), pointer                :: pelpo_loc(:)
    integer(ip), pointer                :: lnnod_loc(:)

    integer(ip)                         :: nfacg
    integer(ip), pointer                :: facel(:,:,:)
    type(i1p),   pointer                :: lelfa(:)

    logical(lg), pointer                :: boundary_nodes(:)

    nullify(lelpo_loc)
    nullify(pelpo_loc)
    nullify(lnnod_loc)

    nullify(facel)
    nullify(lelfa)

    nullify(boundary_nodes)

    !----------------------------------------------------------------------
    !
    ! LELPO_LOC,PELPO_LOC: node-to-element graph
    !
    !----------------------------------------------------------------------

    if( INOTMASTER .and. nelem > 0 ) then

       call memory_alloca(memor_dom,'LNNOD_LOC','meshes_list_boundary_nodes',lnnod_loc,nelem)
       mnode_loc = size(lnods,1,KIND=ip)
       mface = 0
       mnodf = 0
       do ielem = 1,nelem
          ielty            = abs(ltype(ielem))
          mface            = max(mface,element_type(ielty) % number_faces)
          mnodf            = max(mnodf,element_type(ielty) % max_face_nodes)
          lnnod_loc(ielem) = element_type(ielty) % number_nodes       
       end do
       call graphs_elepoi(&
            npoin,nelem,mnode_loc,lnods,lnnod_loc,mepoi_loc,pelpo_loc,lelpo_loc)
       call memory_deallo(memor_dom,'LNNOD_LOC','meshes_list_boundary_nodes',lnnod_loc)

       !----------------------------------------------------------------------
       !
       ! LELFA: List of global faces
       !
       !----------------------------------------------------------------------
       !
       ! Allocate memory for lelfa, FACES, CFAEL AND NNODG
       !
       call memory_alloca(memor_dom,'LELFA','meshes_list_boundary_nodes',lelfa,nelem)
       call memory_alloca(memor_dom,'FACEL','meshes_list_boundary_nodes',facel,mnodf+1_ip,mface,nelem)
       !
       ! Construct and sort FACES
       !
       do ielem = 1,nelem                                          
          ielty = abs(ltype(ielem))
          call memory_alloca(memor_dom,'LELFA(IELEM)%L','meshes_list_boundary_nodes',lelfa(ielem)%l,element_type(ielty) % number_faces)
          do iface = 1,element_type(ielty) % number_faces
             pnodf = element_type(ielty) % node_faces(iface) 
             do inodf = 1,pnodf 
                inode = element_type(ielty) % list_faces(inodf,iface) 
                facel(inodf,iface,ielem) = lnods(inode,ielem)
             end do
             facel(mnodf+1,iface,ielem) = 1
             call maths_heap_sort(2_ip,pnodf,facel(:,iface,ielem))
          end do
       end do
       !
       ! Compute FACES
       !
       nfacg=0_ip
       do ielem = 1,nelem                                            ! Compare the faces and 
          ielty = abs(ltype(ielem))                                  ! eliminate the repited faces
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 ) then
                nfacg = nfacg + 1
                ipoin = facel(1,iface,ielem)
                ielpo = pelpo_loc(ipoin)-1
                do while( ielpo < pelpo_loc(ipoin+1)-1 )
                   ielpo = ielpo + 1
                   jelem = lelpo_loc(ielpo)
                   if( jelem /= ielem ) then
                      jelty = abs(ltype(jelem))                      ! eliminate the repited faces
                      jface = 0
                      do while( jface < element_type(jelty) % number_faces )
                         jface = jface + 1
                         if( facel(mnodf+1,jface,jelem) > 0 ) then
                            equal_faces = .true.
                            inodf = 0
                            do while( equal_faces .and. inodf < element_type(jelty) % node_faces(jface) )
                               inodf = inodf + 1 
                               if( facel(inodf,iface,ielem) /= facel(inodf,jface,jelem) ) equal_faces = .false.
                            end do
                            if( equal_faces ) then
                               facel(mnodf+1,iface,ielem) =  jelem                              ! Keep IELEM face
                               facel(mnodf+1,jface,jelem) = -ielem                              ! Elminate JELEM face
                               facel(      1,iface,ielem) = -jface                              ! Remember IFACE face
                               jface                      =  element_type(jelty) % number_faces ! Exit JFACE do
                               ielpo                      =  pelpo_loc(ipoin+1)                 ! Exit JELEM do  
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end do

       call memory_alloca(memor_dom,'BOUNDARY_NODES','meshes_list_boundary_nodes',boundary_nodes,npoin)

       do ielem = 1,nelem                      
          ielty = abs(ltype(ielem))                    
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
                pnodf = element_type(ielty) % node_faces(iface)
                do inodf = 1,pnodf
                   inode = element_type(ielty) % list_faces(inodf,iface) 
                   ipoin = lnods(inode,ielem)
                   boundary_nodes(ipoin) = .true.
                end do
             end if
          end do
       end do
       number_boundary_nodes = count(boundary_nodes,KIND=ip)
       if( associated(list_boundary_nodes) ) call runend('BOUNDARY NODE LIST ALREADY ASSOCIATED')
       call memory_alloca(memor_dom,'LELFA','meshes_list_boundary_nodes',list_boundary_nodes,number_boundary_nodes)
       if( number_boundary_nodes > 0 ) then
          ii = 0
          do ipoin = 1,npoin
             if( boundary_nodes(ipoin) ) then
                ii = ii + 1
                list_boundary_nodes(ii) = ipoin
             end if
          end do
       end if
       !
       ! Deallocate memory
       !
       call memory_deallo(memor_dom,'FACEL','meshes_list_boundary_nodes',facel)
       call memory_deallo(memor_dom,'LELFA','meshes_list_boundary_nodes',lelfa)
       call memory_deallo(memor_dom,'LELFA','meshes_list_boundary_nodes',boundary_nodes)
       call graphs_dealep(pelpo_loc,lelpo_loc)

    end if

  end subroutine meshes_list_boundary_nodes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-31
  !> @brief   List of boundary nodes
  !> @details List of boundary noed of a mesh. No parallel exchange!
  !> 
  !-----------------------------------------------------------------------

  subroutine meshes_list_boundary_elements(&
       nelem,npoin,lnods,ltype,&
       number_boundary_elements,&
       list_boundary_elements,&
       MEMORY_COUNTER)

    use mod_maths,  only : maths_heap_sort
    use mod_elmgeo, only : element_type

    integer(ip),          intent(in)    :: nelem
    integer(ip),          intent(in)    :: npoin
    integer(ip), pointer, intent(in)    :: lnods(:,:)
    integer(ip), pointer, intent(in)    :: ltype(:)
    integer(ip),          intent(out)   :: number_boundary_elements
    integer(ip), pointer, intent(inout) :: list_boundary_elements(:)
    integer(8),  optional               :: MEMORY_COUNTER(2)
    
    integer(ip)                         :: ielty,ielem,iface,inodf,ilist,ifacg
    integer(ip)                         :: inode,jelem,jface,jelty,ipoin,pnodf
    integer(ip)                         :: pepoi,ielpo,knode,mface,mnodf,ii,mnode_loc
    integer(8)                          :: memor_loc(2)
    logical(lg)                         :: equal_faces  

    integer(ip)                         :: mepoi_loc
    integer(ip), pointer                :: lelpo_loc(:)
    integer(ip), pointer                :: pelpo_loc(:)
    integer(ip), pointer                :: lnnod_loc(:)

    integer(ip)                         :: nfacg
    integer(ip), pointer                :: facel(:,:,:)
    type(i1p),   pointer                :: lelfa(:)

    logical(lg), pointer                :: boundary_elements(:)

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = memor_dom
    end if
    
    nullify(lelpo_loc)
    nullify(pelpo_loc)
    nullify(lnnod_loc)

    nullify(facel)
    nullify(lelfa)

    nullify(boundary_elements)

    !----------------------------------------------------------------------
    !
    ! LELPO_LOC,PELPO_LOC: node-to-element graph
    !
    !----------------------------------------------------------------------

    if( INOTMASTER ) then

       call memory_alloca(memor_loc,'LNNOD_LOC','meshes_list_boundary_elements',lnnod_loc,nelem)
       mnode_loc = size(lnods,1)
       mface = 0
       mnodf = 0
       do ielem = 1,nelem
          ielty            = abs(ltype(ielem))
          mface            = max(mface,element_type(ielty) % number_faces)
          mnodf            = max(mnodf,element_type(ielty) % max_face_nodes)
          lnnod_loc(ielem) = element_type(ielty) % number_nodes       
       end do
       call graphs_elepoi(&
            npoin,nelem,mnode_loc,lnods,lnnod_loc,mepoi_loc,pelpo_loc,lelpo_loc)
       call memory_deallo(memor_loc,'LNNOD_LOC','meshes_list_boundary_elements',lnnod_loc)

       !----------------------------------------------------------------------
       !
       ! LELFA: List of global faces
       !
       !----------------------------------------------------------------------
       !
       ! Allocate memory for lelfa, FACES, CFAEL AND NNODG
       !
       call memory_alloca(memor_loc,'LELFA','meshes_list_boundary_elements',lelfa,nelem)
       call memory_alloca(memor_loc,'FACEL','meshes_list_boundary_elements',facel,mnodf+1_ip,mface,nelem)
       !
       ! Construct and sort FACES
       !
       do ielem = 1,nelem                                          
          ielty = abs(ltype(ielem))
          call memory_alloca(memor_loc,'LELFA(IELEM)%L','meshes_list_boundary_elements',lelfa(ielem)%l,element_type(ielty) % number_faces)
          do iface = 1,element_type(ielty) % number_faces
             pnodf = element_type(ielty) % node_faces(iface) 
             do inodf = 1,pnodf 
                inode = element_type(ielty) % list_faces(inodf,iface) 
                facel(inodf,iface,ielem) = lnods(inode,ielem)
             end do
             facel(mnodf+1,iface,ielem) = 1
             call maths_heap_sort(2_ip,pnodf,facel(:,iface,ielem))
          end do
       end do
       !
       ! Compute FACES
       !
       nfacg=0_ip
       do ielem = 1,nelem                                            ! Compare the faces and 
          ielty = abs(ltype(ielem))                                  ! eliminate the repited faces
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 ) then
                nfacg = nfacg + 1
                ipoin = facel(1,iface,ielem)
                ielpo = pelpo_loc(ipoin)-1
                do while( ielpo < pelpo_loc(ipoin+1)-1 )
                   ielpo = ielpo + 1
                   jelem = lelpo_loc(ielpo)
                   if( jelem /= ielem ) then
                      jelty = abs(ltype(jelem))                      ! eliminate the repited faces
                      jface = 0
                      do while( jface < element_type(jelty) % number_faces )
                         jface = jface + 1
                         if( facel(mnodf+1,jface,jelem) > 0 ) then
                            equal_faces = .true.
                            inodf = 0
                            do while( equal_faces .and. inodf < element_type(jelty) % node_faces(jface) )
                               inodf = inodf + 1 
                               if( facel(inodf,iface,ielem) /= facel(inodf,jface,jelem) ) equal_faces = .false.
                            end do
                            if( equal_faces ) then
                               facel(mnodf+1,iface,ielem) =  jelem                              ! Keep IELEM face
                               facel(mnodf+1,jface,jelem) = -ielem                              ! Elminate JELEM face
                               facel(      1,iface,ielem) = -jface                              ! Remember IFACE face
                               jface                      =  element_type(jelty) % number_faces ! Exit JFACE do
                               ielpo                      =  pelpo_loc(ipoin+1)                 ! Exit JELEM do  
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end do       
       
       call memory_alloca(memor_loc,'BOUNDARY_ELEMENTS','meshes_list_boundary_elements',boundary_elements,nelem)

       do ielem = 1,nelem                      
          ielty = abs(ltype(ielem))                    
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
                boundary_elements(ielem) = .true.
             end if
          end do
       end do
       number_boundary_elements = count(boundary_elements)
       if( associated(list_boundary_elements) ) call runend('BOUNDARY NODE LIST ALREADY ASSOCIATED')
       call memory_alloca(memor_loc,'LELFA','meshes_list_boundary_elements',list_boundary_elements,number_boundary_elements)
       if( number_boundary_elements > 0 ) then          
          ii = 0
          do ielem = 1,nelem
             if( boundary_elements(ielem) ) then
                ii = ii + 1
                list_boundary_elements(ii) = ielem
             end if
          end do
       end if
       !
       ! Deallocate memory
       !
       call memory_deallo(memor_loc,'FACEL','meshes_list_boundary_elements',facel)
       call memory_deallo(memor_loc,'LELFA','meshes_list_boundary_elements',lelfa)
       call memory_deallo(memor_loc,'LELFA','meshes_list_boundary_elements',boundary_elements)
       call graphs_dealep(pelpo_loc,lelpo_loc)

    end if
   
    if( present(MEMORY_COUNTER) ) then
       MEMORY_COUNTER = memor_loc
    else
       memor_dom = memor_loc 
    end if

  end subroutine meshes_list_boundary_elements

end module mod_meshes
!> @}
