!-----------------------------------------------------------------------
!> @defgroup Output_Toolbox
!> Toolbox for output of meshes, element matrices
!> @{
!> @name    ToolBox for output
!> @file    mod_output.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for output
!> @details ToolBox for output, mainly for debugging
!
!-----------------------------------------------------------------------

module mod_output

  use def_kintyp, only : ip,rp,lg
  use def_domain, only : cenam,nnode,cetop
  use def_domain, only : memor_dom,nnode
  use def_domain, only : mesh_type
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_std
  implicit none

  private

  public :: output_mesh_gid_format
  public :: output_element_system
  public :: output_domain

contains

  !-----------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    18/09/2012
  !> @brief   Output boudnary mesh
  !> @details Output boudnary mesh in STL format
  !-----------------------------------------------------------------------

  subroutine output_stl()

    use def_kintyp,         only :  ip,rp,i1p
    use def_kintyp,         only :  comm_data_par
    use def_master,         only :  INOTMASTER,lninv_loc
    use def_master,         only :  INOTSLAVE
    use def_master,         only :  lun_tempo
    use def_master,         only :  namda
    use def_master,         only :  zeror
    use def_master,         only :  kfl_paral
    use def_master,         only :  intost
    use def_kermod,         only :  kfl_oustl
    use def_domain,         only :  cenam,ndime,nboun
    use def_domain,         only :  ltypb,lnnob,coord,npoin
    use def_domain,         only :  lnodb

    use def_domain,         only :  nelem,mnode,mnodb,nelty
    use def_domain,         only :  lnods,lnnod,ltype,nnodf
    use def_domain,         only :  nface,pelpo,lelpo,lface
    use def_domain,         only :  mface,ltypf

    use def_elmtyp,         only :  TRI03,QUA04
    use mod_communications, only :  PAR_DEFINE_COMMUNICATOR
    use mod_communications, only :  PAR_GATHER
    use mod_communications, only :  PAR_GATHERV
    use mod_parall,         only :  PAR_CODE_SIZE
    use mod_iofile,         only :  iofile
    use mod_graphs,         only :  graphs_list_faces
    use mod_graphs,         only :  graphs_deallocate_list_faces
    use mod_boundary_coarsening, only: boundary_coarsening,boundary_coarsening_graph
    use mod_messages,       only : livinf
    use mod_iofile,         only : iofile_flush_unit
    implicit none
    type(comm_data_par), pointer :: commu
    integer(ip)                  :: idime,istl
    integer(ip)                  :: pblty,iboun,inodb,ipoin
    integer(ip)                  :: ipart,ipoi1,ipoi2,ipoi3
    integer(ip)                  :: ifacg,ielem,iface,ielty
    integer(4)                   :: nstl4,nstl4_tot
    real(rp)                     :: vec(3,3),vnor,rboun
    integer(4),          pointer :: nstl4_gat(:)
    real(rp),            pointer :: xstl_gat(:)
    real(rp),            pointer :: xstl(:,:)
    character(150)               :: fil_tempo

    integer(ip)                  :: nfacg
    integer(ip),         pointer :: lfacg(:,:)
    type(i1p),           pointer :: lelfa(:)

    integer(ip)                  :: nboun_coarse          !< Number of boundaries
    integer(ip)                  :: npoin_coarse          !< Number of nodes
    integer(ip),         pointer :: lnodb_coarse(:,:)     !< Boundary connectivity
    real(rp),            pointer :: coord_coarse(:,:)     !< Boundary connectivity

    integer(ip)                  :: nbounThreshold
    integer(ip)                  :: coarseningPoints
    integer(ip)                  :: maxBoundElems

    select case ( kfl_oustl )

    case ( 1_ip )

       !-------------------------------------------------------------------
       !
       ! Master outputs the geometrical STL using lnodb
       !
       !-------------------------------------------------------------------

       if( ndime == 3 ) then

          call livinf(0_ip,'OUTPUT BOUNDARY MESH IN STL FORMAT',0_ip)
          nullify(nstl4_gat)
          nullify(xstl_gat)
          nullify(xstl)
          nstl4 = 0
          !call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
          if( INOTSLAVE ) then
             !
             ! Open file
             !
             fil_tempo = adjustl(trim(namda))//'.stl'
             call iofile(0_ip,lun_tempo,fil_tempo,'SYSTEM INFO')
          end if
          !
          ! Gather all STL
          !
          allocate( nstl4_gat(0:PAR_CODE_SIZE-1) )
          !
          ! Workers convert Alya format to STL
          !
          if( INOTMASTER ) then
             do iboun = 1,nboun
                pblty = ltypb(iboun)
                if( pblty == TRI03 ) then
                   nstl4  = nstl4 + 1
                else if( pblty == QUA04 ) then
                   nstl4  = nstl4 + 2
                end if
             end do
             allocate( xstl(9,nstl4) )
             nstl4 = 0
             do iboun = 1,nboun
                pblty = ltypb(iboun)
                if( pblty == TRI03 ) then
                   nstl4  = nstl4 + 1
                   idime = 1
                   do inodb = 1,lnnob(iboun)
                      ipoin = lnodb(inodb,iboun)
                      xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoin)
                      idime = idime + 3
                   end do
                else if( pblty == QUA04 ) then
                   nstl4  = nstl4 + 1
                   idime = 1
                   ipoi1 = lnodb(1,iboun)
                   ipoi2 = lnodb(2,iboun)
                   ipoi3 = lnodb(3,iboun)
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi1) ; idime = idime + 3
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi2) ; idime = idime + 3
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi3)
                   nstl4  = nstl4 + 1
                   idime = 1
                   ipoi1 = lnodb(1,iboun)
                   ipoi2 = lnodb(3,iboun)
                   ipoi3 = lnodb(4,iboun)
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi1) ; idime = idime + 3
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi2) ; idime = idime + 3
                   xstl(idime:idime+2,nstl4) = coord(1:ndime,ipoi3)
                else
                   call runend('OUTSTL: CANNOT GENERATE STL FOR THIS BOUNDARY TYPE '//cenam(pblty))
                end if
             end do
          end if
          !
          ! Gather all STL
          !
          nstl4 = 9_4*nstl4
          call PAR_GATHER(nstl4,nstl4_gat,'IN MY CODE')

          if( INOTSLAVE ) then
             nstl4_tot = 0
             do ipart = 0,PAR_CODE_SIZE-1
                nstl4_tot = nstl4_tot + nstl4_gat(ipart)
             end do
             allocate( xstl_gat(nstl4_tot) )
          end if

          call PAR_GATHERV(xstl,xstl_gat,nstl4_gat,'IN MY CODE')
          !
          ! Master outputs STL
          !
          if( INOTSLAVE ) then
             write(lun_tempo,'(a)') 'solid mesh_boundary'
             nstl4 = 0
             do ipart = 0,PAR_CODE_SIZE-1
                do istl = 1,nstl4_gat(ipart)/9
                   vec(1,1) = xstl_gat(nstl4+4) - xstl_gat(nstl4+1)
                   vec(2,1) = xstl_gat(nstl4+5) - xstl_gat(nstl4+2)
                   vec(3,1) = xstl_gat(nstl4+6) - xstl_gat(nstl4+3)
                   vec(1,2) = xstl_gat(nstl4+7) - xstl_gat(nstl4+1)
                   vec(2,2) = xstl_gat(nstl4+8) - xstl_gat(nstl4+2)
                   vec(3,2) = xstl_gat(nstl4+9) - xstl_gat(nstl4+3)
                   call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                   call vecnor(vec(1,3),ndime,vnor,2_ip)
                   vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                   write(lun_tempo,1)     'facet normal',vec(1:3,3)
                   write(lun_tempo,'(a)') 'outer loop'
                   write(lun_tempo,1)     'vertex',xstl_gat(nstl4+1:nstl4+3)
                   write(lun_tempo,1)     'vertex',xstl_gat(nstl4+4:nstl4+6)
                   write(lun_tempo,1)     'vertex',xstl_gat(nstl4+7:nstl4+9)
                   write(lun_tempo,'(a)') 'endloop'
                   write(lun_tempo,'(a)') 'endfacet'
                   nstl4 = nstl4 + 9
                end do
             end do
             write(lun_tempo,'(a)') 'endsolid mesh_boundary'
             deallocate( nstl4_gat )
             deallocate( xstl_gat  )
             call iofile_flush_unit(lun_tempo)
             close(lun_tempo)
          else
             deallocate( xstl )
          end if

       end if

    case ( 2_ip )

       !-------------------------------------------------------------------
       !
       ! Each slaves output its STL including subdomain interfaces
       !
       !-------------------------------------------------------------------

       if( INOTMASTER ) then
          nullify(lfacg)
          nullify(lelfa)
          call graphs_list_faces(&
               nelem,mnode,mnodb,nelty,mface,lnods,lnnod,ltype,&
               nnodf,lface,nface,pelpo,lelpo,nfacg,lfacg,lelfa)

          nbounThreshold = 1000
          rboun = real(nboun,rp)
          !maxBoundElems = 100000
          !coarseningPoints = maxBoundElems*(tanh( nboun/dfloat(maxBoundElems) ) +1)/2.0
          if(nboun>nbounThreshold) then!!!!
             !
             if(nboun<10000) then
                coarseningPoints = 2 + 1*(int(tanh(rboun/1000.0_rp-1.0_rp),ip)  +1)/2
             else if(nboun<100000) then
                coarseningPoints = 4 + 2*(int(tanh(rboun/10000.0_rp-1.0_rp),ip) +1)/2
             else
                coarseningPoints = 6 + 4*(int(tanh(rboun/100000.0_rp-1.0_rp),ip)+1)/2
             end if

             !-***default call
             !            call boundary_coarsening_graph(&
             !                nelem,mnode,mnodb,nelty,mface,lnods,lnnod,ltype,&
             !                nnodf,lface,nface,pelpo,lelpo,nfacg,lfacg,lelfa,&
             !                ndime,npoin,coord,&
             !                nboun_coarse,npoin_coarse,lnodb_coarse,coord_coarse)!Default:'divide',10
             !-*** call with optional arguments
             !-** 'divide' requires the factor to compute npoin_coarse=npoin/factor
             call boundary_coarsening_graph(&
                  nelem,mnode,mnodb,nelty,mface,lnods,lnnod,ltype,&
                  nnodf,lface,nface,pelpo,lelpo,nfacg,lfacg,lelfa,&
                  ndime,npoin,coord,&
                  nboun_coarse,npoin_coarse,lnodb_coarse,coord_coarse,&
                  'divide',coarseningPoints)
             !-** 'npoin_coarse' provides npoin_coarse (approximately)
             !            call boundary_coarsening_graph(&
             !                nelem,mnode,mnodb,nelty,mface,lnods,lnnod,ltype,&
             !                nnodf,lface,nface,pelpo,lelpo,nfacg,lfacg,lelfa,&
             !                ndime,npoin,coord,&
             !                nboun_coarse,npoin_coarse,lnodb_coarse,coord_coarse,&
             !                'npoin_coarse',coarseningPoints)

             fil_tempo = adjustl(trim(namda))//'-'//trim(intost(kfl_paral))//'.stl'
             call iofile(0_ip,lun_tempo,fil_tempo,'SYSTEM INFO')
             write(lun_tempo,'(a)') 'solid subdomain_boundary'
             do iboun = 1,nboun_coarse
                ipoi1    = lnodb_coarse(1,iboun)
                ipoi2    = lnodb_coarse(2,iboun)
                ipoi3    = lnodb_coarse(3,iboun)
                vec(1,1) = coord_coarse(1,ipoi2) - coord_coarse(1,ipoi1)
                vec(2,1) = coord_coarse(2,ipoi2) - coord_coarse(2,ipoi1)
                vec(3,1) = coord_coarse(3,ipoi2) - coord_coarse(3,ipoi1)
                vec(1,2) = coord_coarse(1,ipoi3) - coord_coarse(1,ipoi1)
                vec(2,2) = coord_coarse(2,ipoi3) - coord_coarse(2,ipoi1)
                vec(3,2) = coord_coarse(3,ipoi3) - coord_coarse(3,ipoi1)
                call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                call vecnor(vec(1,3),ndime,vnor,2_ip)
                vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                write(lun_tempo,1)     'facet normal',vec(1:3,3)
                write(lun_tempo,'(a)') 'outer loop'
                write(lun_tempo,1)     'vertex',coord_coarse(1:3,ipoi1)
                write(lun_tempo,1)     'vertex',coord_coarse(1:3,ipoi2)
                write(lun_tempo,1)     'vertex',coord_coarse(1:3,ipoi3)
                write(lun_tempo,'(a)') 'endloop'
                write(lun_tempo,'(a)') 'endfacet'
             end do
             deallocate( coord_coarse )
             deallocate( lnodb_coarse )

          else!!!!

             fil_tempo = adjustl(trim(namda))//'-'//trim(intost(kfl_paral))//'.stl'
             call iofile(0_ip,lun_tempo,fil_tempo,'SYSTEM INFO')
             write(lun_tempo,'(a)') 'solid subdomain_boundary'
             do ifacg = 1,nfacg
                if( lfacg(2,ifacg) == 0 ) then
                   iface = lfacg(3,ifacg)
                   ielem = lfacg(1,ifacg)
                   ielty = ltype(ielem)
                   pblty = ltypf(ielty) % l(iface)
                   if( pblty == TRI03 ) then
                      ipoi1    = lnods(lface(ielty) % l(1,iface),ielem)
                      ipoi2    = lnods(lface(ielty) % l(2,iface),ielem)
                      ipoi3    = lnods(lface(ielty) % l(3,iface),ielem)
                      vec(1,1) = coord(1,ipoi2) - coord(1,ipoi1)
                      vec(2,1) = coord(2,ipoi2) - coord(2,ipoi1)
                      vec(3,1) = coord(3,ipoi2) - coord(3,ipoi1)
                      vec(1,2) = coord(1,ipoi3) - coord(1,ipoi1)
                      vec(2,2) = coord(2,ipoi3) - coord(2,ipoi1)
                      vec(3,2) = coord(3,ipoi3) - coord(3,ipoi1)
                      call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                      call vecnor(vec(1,3),ndime,vnor,2_ip)
                      vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                      write(lun_tempo,1)     'facet normal',vec(1:3,3)
                      write(lun_tempo,'(a)') 'outer loop'
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi1)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi2)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi3)
                      write(lun_tempo,'(a)') 'endloop'
                      write(lun_tempo,'(a)') 'endfacet'
                   else if( pblty == QUA04 ) then
                      ipoi1    = lnods(lface(ielty) % l(1,iface),ielem)
                      ipoi2    = lnods(lface(ielty) % l(2,iface),ielem)
                      ipoi3    = lnods(lface(ielty) % l(3,iface),ielem)
                      vec(1,1) = coord(1,ipoi2) - coord(1,ipoi1)
                      vec(2,1) = coord(2,ipoi2) - coord(2,ipoi1)
                      vec(3,1) = coord(3,ipoi2) - coord(3,ipoi1)
                      vec(1,2) = coord(1,ipoi3) - coord(1,ipoi1)
                      vec(2,2) = coord(2,ipoi3) - coord(2,ipoi1)
                      vec(3,2) = coord(3,ipoi3) - coord(3,ipoi1)
                      call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                      call vecnor(vec(1,3),ndime,vnor,2_ip)
                      vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                      write(lun_tempo,1)     'facet normal',vec(1:3,3)
                      write(lun_tempo,'(a)') 'outer loop'
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi1)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi2)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi3)
                      write(lun_tempo,'(a)') 'endloop'
                      write(lun_tempo,'(a)') 'endfacet'
                      ipoi1    = lnods(lface(ielty) % l(1,iface),ielem)
                      ipoi2    = lnods(lface(ielty) % l(3,iface),ielem)
                      ipoi3    = lnods(lface(ielty) % l(4,iface),ielem)
                      vec(1,1) = coord(1,ipoi2) - coord(1,ipoi1)
                      vec(2,1) = coord(2,ipoi2) - coord(2,ipoi1)
                      vec(3,1) = coord(3,ipoi2) - coord(3,ipoi1)
                      vec(1,2) = coord(1,ipoi3) - coord(1,ipoi1)
                      vec(2,2) = coord(2,ipoi3) - coord(2,ipoi1)
                      vec(3,2) = coord(3,ipoi3) - coord(3,ipoi1)
                      call vecpro(vec(1,1),vec(1,2),vec(1,3),ndime)
                      call vecnor(vec(1,3),ndime,vnor,2_ip)
                      vec(1:3,3) = vec(1:3,3) / ( vnor + zeror )
                      write(lun_tempo,1)     'facet normal',vec(1:3,3)
                      write(lun_tempo,'(a)') 'outer loop'
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi1)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi2)
                      write(lun_tempo,1)     'vertex',coord(1:3,ipoi3)
                      write(lun_tempo,'(a)') 'endloop'
                      write(lun_tempo,'(a)') 'endfacet'
                   else
                      call runend('OUTSTL: CANNOT GENERATE STL FOR THIS BOUNDARY TYPE '//cenam(pblty))
                   end if
                end if
             end do

          end if!!!!!

          write(lun_tempo,'(a)') 'endsolid subdomain_boundary'
          flush(lun_tempo)
          close(lun_tempo)
          call graphs_deallocate_list_faces(lfacg,lelfa)
       end if

    end select

1   format(a,3(1x,e13.6),1x,i7)

  end subroutine output_stl


  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    05/01/2018
  !> @brief   Mesh in Gid format
  !> @details Output a mesh or a submesh in GiD format
  !>          1. Offsets can be used to call this subroutine recursively,
  !>             in order to append meshes: NODE_OFFSET, ELEMENT_OFFSET
  !>          2. To output a submesh, you can mark nodes with
  !>             MARK_NPOIN_OPT or elements with mark_nelem_opt
  !>          3. Boundary output can be disabled using OUTPUT_BOUNDARY='OFF'
  !
  !-----------------------------------------------------------------------

  subroutine output_mesh_gid_format(&
       meshe,title,lunit,mark_npoin_opt,mark_nelem_opt,&
       RESULT_INT,NODE_OFFSET,ELEMENT_OFFSET,OUTPUT_BOUNDARY,&
       MARK_ELEMENTS)

    type(mesh_type), intent(inout)                 :: meshe               !< Mesh type
    character(*),    intent(in)                    :: title               !< Title of the mesh
    integer(ip),     intent(in)                    :: lunit               !< Output unit
    logical(lg),     intent(in), pointer, optional :: mark_npoin_opt(:)   !< If a submesh is required
    logical(lg),     intent(in), pointer, optional :: mark_nelem_opt(:)   !< If a submesh is required
    integer(ip),     intent(in), pointer, optional :: RESULT_INT(:)       !< Integer results
    integer(ip),     intent(in),          optional :: NODE_OFFSET         !< Offset for node numbering
    integer(ip),     intent(in),          optional :: ELEMENT_OFFSET      !< Offset for element numbering
    character(*),    intent(in),          optional :: OUTPUT_BOUNDARY     !< Output boundary mesh
    integer(ip),     intent(in), pointer, optional :: MARK_ELEMENTS(:)    !< Mark elements

    integer(ip)                                    :: iesto,iesta,ibsta
    integer(ip)                                    :: ifirs,inode
    integer(ip)                                    :: ielty,pnodb
    integer(ip)                                    :: mnode,npoin,nelem
    integer(ip)                                    :: nboun,iboun
    integer(ip)                                    :: mnodb
    integer(ip)                                    :: ndime,ipoin,ielem
    integer(ip)                                    :: ipoin_offset
    integer(ip)                                    :: ielem_offset
    character(150)                                 :: dumml
    integer(ip),   pointer                         :: lnods(:,:)
    integer(ip),   pointer                         :: ltype(:)
    integer(ip),   pointer                         :: lnodb(:,:)
    integer(ip),   pointer                         :: ltypb(:)
    real(rp),      pointer                         :: coord(:,:)
    logical(lg),   pointer                         :: mark_nelem(:)
    logical(lg),   pointer                         :: mark_nboun(:)
    logical(lg)                                    :: output_type
    logical(lg)                                    :: boundary
    integer(ip),   pointer                         :: permu_npoin(:)
    !
    ! Nullify
    !
    nullify(lnods)
    nullify(ltype)
    nullify(lnodb)
    nullify(ltypb)
    nullify(coord)
    nullify(mark_nelem)
    nullify(mark_nboun)
    nullify(permu_npoin)
    !
    ! Output boundary
    !
    boundary = .true.
    if( present(OUTPUT_BOUNDARY) ) then
       if( trim(OUTPUT_BOUNDARY) == 'OFF' .or. trim(OUTPUT_BOUNDARY) == 'NO' ) then
          boundary = .false.
       end if
    end if
    !
    ! Submesh dimensions
    !
    mnode = meshe % mnode
    mnodb = meshe % mnodb
    ndime = meshe % ndime
    !
    ! Offsets
    !
    if( present(NODE_OFFSET) ) then
       ipoin_offset = NODE_OFFSET
    else
       ipoin_offset = 0
    end if
    if( present(ELEMENT_OFFSET) ) then
       ielem_offset = ELEMENT_OFFSET
    else
       ielem_offset = 0
    end if

    !--------------------------------------------------------------------
    !
    ! Compute submeshes
    !
    !--------------------------------------------------------------------

    if( present(mark_npoin_opt) ) then
       !
       ! Submesh where nodes are marked: mark associated elements and boundaries
       !
       call memory_alloca(memor_dom,'MARK_NELEM','output_mesh_gid_format',mark_nelem,meshe % nelem)
       call memory_alloca(memor_dom,'MARK_NBOUN','output_mesh_gid_format',mark_nboun,meshe % nboun)
       do ielem = 1,meshe % nelem
          loop_inode: do inode = 1,meshe % lnnod(ielem)
             ipoin = meshe % lnods(inode,ielem)
             if( mark_npoin_opt(ipoin) ) then
                mark_nelem(ielem) = .true.
                exit loop_inode
             end if
          end do loop_inode
       end do
       do iboun = 1,meshe % nboun
          pnodb = meshe % lnnob(iboun)
          ielem = meshe % lelbo(iboun)
          if( mark_nelem(ielem) ) mark_nboun(iboun) = .true.
       end do
    else if( present(mark_nelem_opt) ) then
       mark_nelem => mark_nelem_opt
    end if

    if( present(mark_npoin_opt) .or. present(mark_nelem_opt) ) then
       !
       ! Submesh where elements are marked
       !
       call memory_alloca(memor_dom,'PERMU_NPOIN','output_mesh_gid_format',permu_npoin,meshe % npoin)

       nelem = count(mark_nelem(1:meshe % nelem),KIND=ip)
       nboun = count(mark_nboun(1:meshe % nboun),KIND=ip)
       call memory_alloca(memor_dom,'LNODS','output_mesh_gid_format',lnods,mnode,nelem)
       call memory_alloca(memor_dom,'LTYPE','output_mesh_gid_format',ltype,nelem)
       call memory_alloca(memor_dom,'LNODB','output_mesh_gid_format',lnodb,mnodb,nboun)
       call memory_alloca(memor_dom,'LTYPB','output_mesh_gid_format',ltypb,nboun)
       nelem = 0
       do ielem = 1,meshe % nelem
          if( mark_nelem(ielem) ) then
             nelem = nelem + 1
             lnods(:,nelem) = meshe % lnods(:,ielem)
             ltype(nelem)   = meshe % ltype(ielem)
          end if
       end do
       nboun = 0
       do iboun = 1,meshe % nboun
          if( mark_nboun(iboun) ) then
             nboun = nboun + 1
             lnodb(:,nboun) = meshe % lnodb(:,iboun)
             ltypb(nboun)   = meshe % ltypb(iboun)
          end if
       end do

       npoin = 0
       do ielem = 1,meshe % nelem
          if( mark_nelem(ielem) ) then
             do inode = 1,meshe % lnnod(ielem)
                ipoin = meshe % lnods(inode,ielem)
                permu_npoin(ipoin) = 1
             end do
          end if
       end do
       npoin = count(permu_npoin==1,KIND=ip)
       call memory_alloca(memor_dom,'COORD','output_mesh_gid_format',coord,ndime,npoin)

       npoin = 0
       do ipoin = 1,meshe % npoin
          if( permu_npoin(ipoin) == 1 ) then
             npoin = npoin + 1
             permu_npoin(ipoin) = npoin
             coord(1:ndime,npoin) = meshe % coord(1:ndime,ipoin)
          end if
       end do
       do ielem = 1,nelem
          do inode = 1,mnode
             ipoin = lnods(inode,ielem)
             if( ipoin > 0 ) lnods(inode,ielem) = permu_npoin(ipoin)
          end do
       end do
       do iboun = 1,nboun
          do inode = 1,mnodb
             ipoin = lnodb(inode,iboun)
             if( ipoin > 0 ) lnodb(inode,iboun) = permu_npoin(ipoin)
          end do
       end do

    else
       !
       ! Whole mesh
       !
       npoin =  meshe % npoin
       nelem =  meshe % nelem
       nboun =  meshe % nboun
       lnods => meshe % lnods
       ltypb => meshe % ltypb
       lnodb => meshe % lnodb
       ltype => meshe % ltype
       coord => meshe % coord

    end if
    !
    ! Element range
    !
    if(      ndime == 1 ) then
       ibsta =  1
       iesta =  2
       iesto =  9
    else if( ndime == 2 ) then
       ibsta =  2
       iesta = 10
       iesto = 29
    else if( ndime == 3 ) then
       ibsta = 10
       iesta = 30
       iesto = 50
    end if

    !--------------------------------------------------------------------
    !
    ! Output mesh
    !
    !--------------------------------------------------------------------

    ifirs = 0
    do ielty = ibsta,iesto

       output_type = .false.
       if( any(abs(ltype)==ielty) ) output_type = .true.
       if( nboun > 0 ) then
          if( any(abs(ltypb)==ielty) ) output_type = .true.
       end if

       if( output_type ) then

          dumml = adjustl(trim(title)//'_'//cenam(ielty))
          !
          ! Header
          !
          write(lunit,1)&
               adjustl(trim(dumml)),ndime,&
               adjustl(trim(cetop(ielty))),nnode(ielty)
          !
          ! Coordinates
          !
          if( ifirs == 0 .and. npoin > 0 ) then
             ifirs = 1
             write(lunit,2) 'coordinates'
             if( ndime == 1 ) then
                do ipoin = 1,npoin
                   write(lunit,3) ipoin+ipoin_offset,coord(1,ipoin),0.0_rp
                end do
             else
                do ipoin = 1,npoin
                   write(lunit,3) ipoin+ipoin_offset,coord(1:ndime,ipoin)
                end do
             end if
             write(lunit,2) 'end coordinates'
          end if

          write(lunit,2) 'elements'

          if( ielty >= iesta ) then
             !
             ! Element connectivity
             !
             if( present(MARK_ELEMENTS) ) then
                do ielem = 1,nelem
                   if( abs(ltype(ielem)) == ielty ) then
                      write(lunit,4) ielem+ielem_offset,(lnods(inode,ielem),inode=1,nnode(ielty),MARK_ELEMENTS(ielem))
                   end if
                end do
             else
                do ielem = 1,nelem
                   if( abs(ltype(ielem)) == ielty ) then
                      write(lunit,4) ielem+ielem_offset,(lnods(inode,ielem),inode=1,nnode(ielty))
                   end if
                end do
             end if

          else if( boundary ) then
             !
             ! Boundary connectivity
             !
             do iboun = 1,nboun
                if( abs(ltypb(iboun)) == ielty ) then
                   write(lunit,4) iboun+ielem_offset,(lnodb(inode,iboun),inode=1,nnode(ielty))
                end if
             end do
          end if
          write(lunit,2) 'end elements'
          write(lunit,2) ''

       end if

    end do
    !
    ! Results
    !
    if( present(RESULT_INT) ) then
       if( present(mark_npoin_opt) .or. present(mark_nelem_opt) ) then
          call output_result_gid_format(lunit,meshe % npoin,RESULT_INT,'GENERIC',permu_npoin)
       else
          call output_result_gid_format(lunit,meshe % npoin,RESULT_INT,'GENERIC')
       end if
    end if
    !
    ! Deallocate
    !
    call memory_deallo(memor_dom,'PERMU_NPOIN','output_mesh_gid_format',permu_npoin)

    if( present(mark_npoin_opt) ) then

       call memory_deallo(memor_dom,'MARK_NELEM','output_mesh_gid_format',mark_nelem)
       call memory_deallo(memor_dom,'MARK_NBOUN','output_mesh_gid_format',mark_nboun)

    end if

    if( present(mark_npoin_opt) .or. present(mark_nelem_opt) ) then

       call memory_deallo(memor_dom,'LNODS','output_mesh_gid_format',lnods)
       call memory_deallo(memor_dom,'LTYPE','output_mesh_gid_format',ltype)
       call memory_deallo(memor_dom,'LNODB','output_mesh_gid_format',lnodb)
       call memory_deallo(memor_dom,'LTYPB','output_mesh_gid_format',ltypb)
       call memory_deallo(memor_dom,'COORD','output_mesh_gid_format',coord)

    end if

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i9, 3(1x,e16.8e3))
4   format(i9,50(1x,i9))

  end subroutine output_mesh_gid_format

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    05/01/2018
  !> @brief   Result in Gid format
  !> @details Output a result with possible permutation
  !
  !-----------------------------------------------------------------------

  subroutine output_result_gid_format(lunit,npoin,RESULT_INT,NAME,permu_npoin)

    integer(ip),     intent(in)                    :: lunit
    integer(ip),     intent(in)                    :: npoin
    integer(ip),     intent(in), pointer           :: RESULT_INT(:)       !< Integer results
    character(*),    intent(in),          optional :: NAME
    integer(ip),     intent(in), pointer, optional :: permu_npoin(:)
    integer(ip)                                    :: ipoin,kpoin,lmax
    character(10)                                  :: my_name

    write(lunit,1) 'GiD Post Results File 1.0'
    write(lunit,1) ' '

    if( present(NAME) ) then
       lmax = min(10,len(NAME))
       my_name = trim(NAME(1:lmax))
    else
       my_name = 'GENERIC'
    end if
    write(lunit,2) trim(my_name),'ALYA',1.0_rp,'Scalar'
    write(lunit,3) trim(my_name)

    write(lunit,1) 'Values'

    if( present(permu_npoin) ) then
       do ipoin = 1,npoin
          if( permu_npoin(ipoin) > 0 ) then
             kpoin = permu_npoin(ipoin)
             write(lunit,6) kpoin,RESULT_INT(kpoin)
          end if
       end do
    else
       do ipoin = 1,npoin
          write(lunit,6) ipoin,RESULT_INT(ipoin)
       end do
    end if

    write(lunit,1) 'Values'

1   format(a)
2   format('Result ',a,' ',a,' ',e15.8,' ',a,' OnNodes')
3   format('ComponentNames ',a)
6   format(i9, 3(1x,i8))

  end subroutine output_result_gid_format

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux version of nastin (adapted by ECR)
  !> @brief   Output system
  !> @details Output elemental system
  !>
  !-----------------------------------------------------------------------

  subroutine output_element_system(kdime,pnode,elmat,elrhs)

    use def_master, only :  intost

    implicit none

    integer(ip), intent(in)           :: kdime
    integer(ip), intent(in)           :: pnode
    ! Element matrices
    real(rp),    intent(in)           :: elmat(kdime*pnode,kdime*pnode)
    real(rp),    intent(in)           :: elrhs(kdime,pnode)

    integer(ip),   parameter          :: lreal=9
    character(13)                     :: FMT1
    integer(ip)                       :: inode,idime,idofn,jnode,jdime,jdofn
    character(3)                      :: vanam

    !
    ! Format
    !
    FMT1 = '(1x,e'//trim(intost(lreal))//'.'//trim(intost(lreal-7))
    !
    ! First line
    !
    write(99,'(a)',advance='no') '+-'
    do idofn = 1,kdime*pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do
    write(99,'(a)',advance='no') ' '
    do idofn = 1,pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do

    write(99,'(a)',advance='no') '-+'
    write(99,'(a)',advance='no') '  +-   -+'

    write(99,'(a)',advance='no') '     +-'
    do idime = 1,lreal
       write(99,'(a)',advance='no') ' '
    end do
    write(99,'(a)',advance='no') '-+'

    write(99,*)

    do inode = 1,pnode
       do idime = 1,kdime
          !
          ! Auu
          !
          write(99,'(a)',advance='no') '|'
          idofn = (inode-1)*kdime+idime
          do jnode = 1,pnode
             do jdime = 1,kdime
                jdofn = (jnode-1)*kdime+jdime
                write(99,FMT1,advance='no') elmat(idofn,jdofn)
             end do
          end do
          !
          ! ui
          !
          write(99,'(a)',advance='no') ' |'
          if( idime == 1 ) then
             vanam = 'u'
          else if( idime == 2 ) then
             vanam = 'v'
          else if( idime == 3 ) then
             vanam = 'w'
          end if
          vanam = trim(vanam) // trim(intost(inode))
          write(99,'(a)',advance='no') '  | '//vanam//' |'
          !
          ! bu
          !
          write(99,'(a)',advance='no') '     |'
          write(99,FMT1,advance='no') elrhs(idime,inode)
          write(99,'(a)',advance='no') ' |'

          write(99,*)
       end do
    end do
    !
    ! Last line
    !
    write(99,'(a)',advance='no') '+-'
    do idofn = 1,pnode*kdime
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do

    write(99,'(a)',advance='no') ' '
    do idofn = 1,pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do


    write(99,'(a)',advance='no') '-+'
    write(99,'(a)',advance='no') '  +-   -+'

    write(99,'(a)',advance='no') '     +-'
    do idime = 1,lreal
       write(99,'(a)',advance='no') ' '
    end do
    write(99,'(a)',advance='no') '-+'

    write(99,*)

  end subroutine output_element_system

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-05-02
  !> @brief   Output domain
  !> @details Output the mesh
  !>
  !-----------------------------------------------------------------------

  subroutine output_domain(CURRENT_MESH,ONLY_MESH)

    use def_parame
    use def_elmtyp
    use def_master
    use def_kermod
    use def_domain
    use def_mpio
    use mod_postpr
    use mod_memory
    use mod_messages, only : livinf
    use def_mpio,     only : mpio_flag_geometry_export, PAR_MPIO_ON
    use mod_iofile,   only : iofile_flush_unit
    implicit none

    integer(ip), optional, intent(in) :: CURRENT_MESH
    logical(lg), optional, intent(in) :: ONLY_MESH
    integer(ip)                       :: ipart,ifiel,kdime,istep,idivi
    logical(lg)                       :: output_boundaries
    real(rp),    pointer              :: dumm2(:,:),dumm1(:)
    logical(lg)                       :: if_only_mesh
    logical(lg)                       :: mpio_light
    logical(lg)                       :: mpio_not_light
    !
    ! Options
    !
    if( present(ONLY_MESH) ) then
       if_only_mesh = ONLY_MESH
    else
       if_only_mesh = .false.
    end if
    if( present(CURRENT_MESH) ) then
       idivi = CURRENT_MESH
    else
       idivi = kfl_posdi
    end if
    mpio_light = kfl_mpio_post > IO_CLASSIC .and. mpio_flag_post_light == PAR_MPIO_ON
    mpio_not_light = kfl_mpio_post > IO_CLASSIC .and. mpio_flag_post_light == PAR_MPIO_OFF
    nullify(dumm2,dumm1)

    if( maxval(kfl_oumes) == 1 .or. kfl_mpio_post > IO_CLASSIC ) then

       !-----------------------------------------------------------------
       !
       ! Output mesh information
       !
       !-----------------------------------------------------------------

       if( .not. if_only_mesh ) then
          call livinf(0_ip,'OUTPUT MESH',0_ip)
          if( IMASTER ) then
             write(lun_pos00,1) npart
             do ipart = 1,npart
                write(lun_pos00,2) ipart,&
                     &             meshe(idivi) % nelem_par(ipart),&
                     &             meshe(idivi) % npoin_par(ipart),&
                     &             meshe(idivi) % nboun_par(ipart)
             end do
          else if( ISEQUEN ) then
             write(lun_pos00,1) 1_ip
             write(lun_pos00,2) 1_ip,&
                  &          meshe(idivi) % nelem,&
                  &          meshe(idivi) % npoin,&
                  &          meshe(idivi) % nboun
          end if
          if( INOTSLAVE ) call iofile_flush_unit(lun_pos00)
       end if

       !-----------------------------------------------------------------
       !
       ! Output basic element, boundary and node arrays
       !
       !-----------------------------------------------------------------

       if( kfl_oumes(1) == 1  .or. kfl_mpio_post > IO_CLASSIC ) then

          call postpr(meshe(idivi) % coord,    postp(1) % wopos(:,16),ittim,cutim,ndime,'ORIGI')    ! COORD
          call postpr(meshe(idivi) % lnods,    postp(1) % wopos(:,15),ittim,cutim,mnode,'ORIGI')    ! LNODS
          call postpr(meshe(idivi) % ltype,    postp(1) % wopos(:,17),ittim,cutim,      'ORIGI')    ! LTYPE
          call postpr(meshe(idivi) % lninv_loc,postp(1) % wopos(:,18),ittim,cutim,      'ORIGI')    ! LNINV_LOC
          call postpr(meshe(idivi) % leinv_loc,postp(1) % wopos(:,23),ittim,cutim,      'ORIGI')    ! LEINV_LOC
          if ( .not. mpio_light ) then
            call postpr(meshe(idivi) % lelch,    postp(1) % wopos(:,19),ittim,cutim,      'ORIGI')    ! LELCH
            call postpr(meshe(idivi) % lnoch,    postp(1) % wopos(:,27),ittim,cutim,      'ORIGI')    ! LNOCH
            call postpr(meshe(idivi) % lesub,    postp(1) % wopos(:,22),ittim,cutim,      'ORIGI')    ! LESUB
            call postpr(meshe(idivi) % lmate,    postp(1) % wopos(:,24),ittim,cutim,      'ORIGI')    ! LMATE
            call postpr(meshe(idivi) % lmast,    postp(1) % wopos(:,75),ittim,cutim,      'ORIGI')    ! LMAST
          end if

          output_boundaries = .false.
          if( ( nboun_total > 0  .or. mpio_not_light ) ) output_boundaries = .true.

          if( output_boundaries ) then
             call postpr(meshe(idivi) % lnodb,    postp(1) % wopos(:,20),ittim,cutim,mnodb,'ORIGI') ! LNODB
             call postpr(meshe(idivi) % ltypb,    postp(1) % wopos(:,21),ittim,cutim,      'ORIGI') ! LTYPB
             call postpr(meshe(idivi) % lboch,    postp(1) % wopos(:,26),ittim,cutim,      'ORIGI') ! LBOCH
             call postpr(meshe(idivi) % lelbo,    postp(1) % wopos(:,33),ittim,cutim,      'ORIGI') ! LELBO
             call postpr(meshe(idivi) % lbinv_loc,postp(1) % wopos(:,25),ittim,cutim,      'ORIGI') ! LBINV_LOC
          end if
       end if

       !-----------------------------------------------------------------
       !
       ! Output SETS field
       !
       !-----------------------------------------------------------------

       if( kfl_oumes(2) == 1  .or. mpio_not_light) then
          if( neset > 0 ) call postpr(leset,postp(1) % wopos(:,30),ittim,cutim,'ORIGI')                          ! LESET
          if( nbset > 0 .and. output_boundaries ) call postpr(lbset,postp(1) % wopos(:,31),ittim,cutim,'ORIGI')  ! LBSET
       end if

       !-----------------------------------------------------------------
       !
       ! Output BOUNDARY_CONDITIONS
       !
       !-----------------------------------------------------------------

       if( kfl_oumes(3) == 1  .or. mpio_not_light) then
          if( kfl_icodb > 0 .and. output_boundaries ) call postpr(kfl_codbo,postp(1) % wopos(:,28),ittim,cutim,      'ORIGI') ! CODBO
          if( kfl_icodn > 0 )                         call postpr(kfl_codno,postp(1) % wopos(:,29),ittim,cutim,mcono,'ORIGI') ! CODNO
       end if

       !-----------------------------------------------------------------
       !
       ! Output FIELDS
       !
       !-----------------------------------------------------------------

       if( kfl_oumes(4) == 1  .or. mpio_not_light) then

          do ifiel = 1,nfiel

             if ( (kfl_field(6,ifiel) /= 1 ) .OR. (mpio_flag_geometry_export == PAR_MPIO_ON) ) then !not on demand or export

                if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
                   postp(1) % wopos(3,32) = 'NELEM'
                else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
                   postp(1) % wopos(3,32) = 'NPOIN'
                else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
                   postp(1) % wopos(3,32) = 'NBOUN'
                else
                   call runend('OUTDOM: UNDEFINED FIELD TYPE')
                end if
                kdime = kfl_field(1,ifiel)


                do istep = 1,kfl_field(4,ifiel)
                   if( kdime == 1 ) then
                      postp(1) % wopos(2,32) = 'SCALA'
                      if( associated(xfiel(ifiel) % a) ) then
                         dumm1 => xfiel(ifiel) % a(1,:,istep)
                         call postpr(dumm1,postp(1) % wopos(:,32),ittim,time_field(ifiel) % a(istep),'ORIGI',TAG1=ifiel,TAG2=istep)          ! XFIEL
                      else
                         call postpr(dumm1,postp(1) % wopos(:,32),ittim,cutim,'ORIGI',TAG1=ifiel,TAG2=istep)
                      end if
                   else
                      postp(1) % wopos(2,32) = 'VECTO'
                      if( associated(xfiel(ifiel) % a) ) then
                         dumm2 => xfiel(ifiel) % a(:,:,istep)
                         call postpr(dumm2,postp(1) % wopos(:,32),ittim,time_field(ifiel) % a(istep),kdime,'ORIGI',TAG1=ifiel,TAG2=istep)    ! XFIEL
                      else
                         call postpr(dumm2,postp(1) % wopos(:,32),ittim,cutim,kdime,'ORIGI',TAG1=ifiel,TAG2=istep)
                      end if
                   end if
                end do
             end if !not on demand
          end do
       end if
    end if

    !-----------------------------------------------------------------
    !
    ! Output boundary mesh in STL format
    !
    !-----------------------------------------------------------------

    if( .not. if_only_mesh .and. kfl_oustl /= 0 ) call output_stl()

1   format(i9)
2   format(10(1x,i9))

  end subroutine output_domain

end module mod_output
!> @}
