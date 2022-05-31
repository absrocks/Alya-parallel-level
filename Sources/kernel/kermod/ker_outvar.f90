subroutine ker_outvar(ivari)
  !------------------------------------------------------------------------
  !****f* Master/ker_output
  ! NAME
  !    ker_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    ker_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_elmtyp
  use def_kermod
  use def_domain
  use def_coupli
  use mod_ker_proper
  use mod_ker_vortex
  use mod_memory,         only : memory_size
  use mod_parall,         only : par_omp_num_colors
  use mod_parall,         only : par_omp_ia_colors
  use mod_parall,         only : par_omp_ja_colors
  use mod_parall,         only : par_omp_num_threads
  use mod_parall,         only : par_omp_nelem_chunk
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_ALLGATHER
  use mod_parall,         only : PAR_COMM_MY_CODE_ARRAY
  use mod_renumbering,    only : renumbering_elements
  use mod_filters,        only : filters_nodal
  use mod_graphs,         only : graphs_eleele
  use mod_graphs,         only : graphs_dealep
  use mod_gradie,         only : gradie
  use mod_projec,         only : projec_elements_to_nodes
  use mod_ker_nsw_visc2,  only : ker_nod_nsw_visc_0
  use mod_outvar,         only : outvar
  use mod_postpr
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: ibopo,ipoin,idime,ielpo,isubd,nenti,kdime
  integer(ip)             :: ielem,icont,iline,jpoin,iboun,inode,ithre
  integer(ip)             :: dummi,inodb,imate,incnt,kpoin,imesh,istack
  integer(ip)             :: nstack,izdom,icoup,icolo,kelem,ii,ifiel,ienti
  real(rp)                :: rutim,qmaxi,qmini,dummr,rmate
  character(5)            :: wopos(3)
  integer(ip), pointer    :: lstack(:)
  integer(ip), pointer    :: list_colors(:)
  real(rp),    pointer    :: xvalu(:)
  real(rp),    pointer    :: xcoo1(:)
  real(rp),    pointer    :: xcoo2(:)
  character(20)           :: wfiel

  if( ivari == 0 ) return
  !
  ! Define postprocess variable
  !
  rutim = cutim
  
  select case ( ivari )

  case ( 1_ip )
     !
     ! EXNOR
     !
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = exnor(idime,1,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 2_ip )
     !
     ! ???
     !
     return

  case ( 3_ip )
     !
     ! Geometrical local basis and type of point: LPOIN, SKCOS
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = skcos(idime,1,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 4_ip )
     !
     ! Geometrical local basis SKCOS(:,1,:)
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = skcos(idime,1,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 5_ip )
     !
     ! Geometrical local basis SKCOS(:,2,:)
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = skcos(idime,2,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 6_ip )
     !
     ! Geometrical local basis SKCOS(:,NDIME,:)
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = skcos(idime,ndime,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 7_ip )
     !
     ! Hanging nodes: LHANG
     !
     !if( nhang > 0 .and. INOTMASTER ) then
     !   call memgen(zero,npoin,zero)
     !   do ihang=1,nhang
     !      ipoin=lhang(1,ihang)
     !      if(ipoin<1.or.ipoin>npoin) then
     !         call runend('ERROR IN OUTDOM')
     !      end if
     !      gesca(ipoin)=real(lhang(0,ihang))
     !   end do
     !end if

  case ( 8_ip )
     !
     ! DISPM: Mesh displacement
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        if( associated(dispm) ) then
           do ipoin = 1,npoin
              do idime = 1,ndime
                 gevec(idime,ipoin) = dispm(idime,ipoin,1)
              end do
           end do
        else
           do ipoin = 1,npoin
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
        !if( nimbo > 0 ) then
        !   kpoin = npoin
        !   do iimbo = 1,nimbo
        !      do ipoib = 1,imbou(iimbo)%npoib
        !         kpoin = kpoin + 1
        !         do idime = 1,ndime
        !            gevec(idime,kpoin) = imbou(iimbo)%cooi2(idime,ipoib) - imbou(iimbo)%cooin(idime,ipoib)
        !         end do
        !      end do
        !   end do
        !end if
     end if
     !if( kfl_outib == 3 ) kfl_outib = 4

  case ( 9_ip )
     !
     ! Boundary points: LPOTY
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           if(lpoty(ipoin)/=0) gesca(ipoin)=1.0_rp
        end do
     end if

  case ( 10_ip )
     !
     ! PONUM: Point Numbering
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(lninv_loc(ipoin),rp)
        end do
     end if

  case ( 11_ip )
     !
     ! Boundary codes: KFL_CODNO
     !
     if( kfl_icodn > 0 ) then
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime = 1,ndime
                 gevec(idime,ipoin)=real(kfl_codno(idime,ipoin),rp)
              end do
           end do
        end if
     else
        return
     end if

  case ( 12_ip )
     !
     ! YWALP: Wall distance
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if( ibopo /= 0 ) gesca(ipoin)=ywalp(ibopo)
        end do
     end if

  case ( 13_ip )
     !
     ! LNTIB: node type
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin)=real(lntib(ipoin),rp)
        end do
     end if

  case ( 14_ip )
     !
     ! HOLES: with node type
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           if( lntib(ipoin) > 0 ) then
              gesca(ipoin) = 1.0_rp
           end if
        end do
     end if

  case ( 15_ip )
     !
     ! DENSI: Density
     !
     if( INOTEMPTY ) then
        call memgen(zero,npoin,zero)
        densi_ker % kfl_nedsm = 1
        call ker_proper('DENSI','NPOIN',dummi,dummi,gesca)
     end if

  case ( 16_ip )
     !
     ! VISCO: viscosity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('VISCO','NPOIN',dummi,dummi,gesca)
     end if

  case ( 17_ip )
     !
     !
     !
     return

  case ( 18_ip )
     !
     ! CONDU: Conductivity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('CONDU','NPOIN',dummi,dummi,gesca)
     end if

  case ( 19_ip )
     !
     ! SPECI: Specific Heat
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('SPHEA','NPOIN',dummi,dummi,gesca)
     end if

  case ( 20_ip )
     !
     ! DUMMY: Dummy variable
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('DUMMY','NPOIN',dummi,dummi,gesca)
     end if

  case ( 21_ip )
     !
     !
     !
     return

  case ( 22_ip )
     !
     ! LEVELS
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,0_ip)
        call meshin(-5_ip)
     end if

  case ( 23_ip )
     !
     ! LGROU_DOM
     !
     if( ngrou_dom <= 0 ) return
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(lgrou_dom(ipoin),rp)
        end do
     end if

  case ( 24_ip )
     !
     ! MASSM
     !
     gesca => vmass

  case ( 25_ip )
     !
     ! MASSC
     !
     gesca => vmasc

  case ( 26_ip )
     !
     ! KFL_GEONO
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if( ibopo /= 0 ) gesca(ipoin)=real(kfl_geono(ibopo),rp)
        end do
     end if

  case ( 27_ip )
     !
     ! Subdomains
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoi1
           gesca(ipoin) = real(kfl_paral,rp)
        end do
        do ipoin = npoi1+1,npoin
           gesca(ipoin) = -1.0_rp
        end do
        ! tuve que descomentariar estas lineas sino daba mal
        ! porque estan comentadas   - habra que agregar if iparall????
        do ipoin = npoi2,npoi3
           gesca(ipoin) = real(kfl_paral,rp)
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
     end if

  case ( 28_ip )
     !
     ! WALLD
     !
     if( kfl_walld == 0 ) return
     gesca => walld

  case ( 29_ip )
     !
     ! ROUGH
     !
     if( kfl_rough < 1 ) return
     gesca => rough

  case(30_ip)
     !
     ! LINEL: Linelets of preconditioner CG
     !
     if( INOTMASTER ) then
        icont = 0
        do ipoin = 1,npoin
           rhsid(ipoin) = 0.0_rp
        end do
        do iline = 1,solve(2) % nline
           icont = icont+1
           do ipoin = solve(2) % lline(iline),solve(2) % lline(iline+1)-1
              jpoin = solve(2) % lrenup(ipoin)
              rhsid(jpoin) = real(icont,rp)
           end do
        end do
        gesca => rhsid
     end if

  case ( 31_ip )
     !
     ! Boundary codes: CODBO
     !
     if( INOTMASTER ) then
        call memgen(zero,nboun,zero)
        do iboun = 1,nboun
           gesca(iboun)=real(kfl_codbo(iboun),rp)
        end do
     end if

  case ( 32_ip )
     !
     ! MATER: Material Numbering (Elemental)
     !
     if( INOTMASTER ) then
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = real(lmate(ielem),rp)
        end do
     end if

  case ( 33_ip )
     !
     ! Nodal material (maximum)
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do imate = 1,nmate
           rmate = real(imate,rp)
           do ii = 1,memory_size(lmatn(imate)%l)
              ipoin = lmatn(imate)%l(ii)
              gesca(ipoin) = max(rmate,gesca(ipoin))
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX')
     end if

  case ( 34_ip )
     !
     ! Node characteristic
     !
     !     if( INOTMASTER ) then
     !        call memgen(zero,npoin,zero)
     !        do ipoin = 1,npoin
     !           gesca(ipoin) = real(lnoch(ipoin),rp)
     !        end do
     !     end if
     return

  case ( 35_ip )
     !
     ! LELEV
     !
     if( INOTMASTER ) then
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = real(lelev(ielem),rp)
        end do
     end if

  case ( 36_ip )
     !
     ! Element quality
     !
     if( INOTMASTER ) then
        call memgen(zero,nelem,zero)
        call qualit(gesca,qmaxi,qmini)
     else
        call qualit(dummr,qmaxi,qmini)
     end if

  case ( 37_ip )
     !
     ! ELNUM: Global Element Number
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = real(leinv_loc(ielem),rp)
        end do
     end if

  case ( 38_ip )
     !
     ! Boundary codes: CODBB
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do iboun = 1,nboun
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              gesca(ipoin)=real(kfl_codbo(iboun),rp)
           end do
        end do
        call pararr('SMA',NPOIN_TYPE,npoin,gesca)
     end if

  case ( 39_ip )
     !
     ! LETIB
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              gesca(ipoin) = max(real(letib(ielem),rp),gesca(ipoin))
           end do
        end do
     end if

  case ( 40_ip )
     !
     ! LTYPE
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              gesca(ipoin) = real(ltype(ielem),rp)
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MIN','IN MY CODE')
     end if

  case ( 41_ip )
     !
     ! LELCH
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              dummr = real(lelch(ielem),rp)
              gesca(ipoin) = max(gesca(ipoin),dummr)
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
     end if

  case ( 42_ip )
     !
     ! Contact
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do incnt = 1,nncnt
           ipoin = lncnt(1,incnt)
           jpoin = lncnt(2,incnt)
           gesca(ipoin) = 1.0_rp
           gesca(jpoin) = 2.0_rp
        end do
     end if

  case ( 43_ip )
     !
     ! R_DOM
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(r_dom(ipoin+1)-r_dom(ipoin),rp)
        end do
     end if

  case ( 44_ip )
     !
     ! VORTX: extraction of vortex core
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
     else
        call memgen(zero,1_ip,1_ip)
     end if
     call ker_vortex(gevec)
     if( IMASTER ) call memgen(two,1_ip,1_ip)

  case ( 45_ip )
     !
     ! Element sets
     !
     if( neset < 1 ) return
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              gesca(ipoin) = max(gesca(ipoin),real(leset(ielem),rp))
           end do
        end do
        call pararr('SMA',NPOIN_TYPE,npoin,gesca)
     end if

  case ( 46_ip )
     !
     ! DISPL_KER
     !
     if( kfl_suppo == 0 ) return
     gevec => displ_ker

  case ( 47_ip )
     !
     ! VELOC
     !
     if( INOTMASTER ) then
        if( associated(veloc) ) then
           gevec => veloc(:,:,1)
        else if( associated(advec) ) then
           gevec => advec(:,:,1)
        end if
     end if

  case ( 48_ip )
     !
     ! LBSET: boundary sets
     !
     if( nbset < 1 ) return
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        gesca = 0.0_rp
        do iboun = 1,nboun
           do inodb = 1,nnode(abs(ltypb(iboun)))
              ipoin = lnodb(inodb,iboun)
              gesca(ipoin) = max(gesca(ipoin),real(lbset(iboun),rp))
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
     end if

  case ( 49_ip )
     !
     ! Connected meshes
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        allocate(lstack(npoin))
        do ipoin = 1,npoin
           lstack(ipoin) = 0
        end do

        imesh        = 0
        kpoin        = 0

        do while( kpoin /= npoin )

           imesh = imesh + 1
           ipoin = 1
           do while( gesca(ipoin) /= 0.0_rp )
              ipoin = ipoin + 1
           end do
           nstack       = 1
           lstack(1)    = ipoin
           gesca(ipoin) = real(imesh,rp)
           istack       = 0
           kpoin        = kpoin + 1

           do
              if( istack == nstack ) exit
              istack = istack + 1
              ipoin  = lstack(istack)
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( gesca(jpoin) == 0.0_rp ) then
                    gesca(jpoin)   = real(imesh,rp)
                    nstack         = nstack + 1
                    lstack(nstack) = jpoin
                    kpoin          = kpoin + 1
                 end if
              end do
           end do
        end do
        deallocate(lstack)

     end if

  case ( 50_ip )
     !
     ! TURBU: Turbulent viscosity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('TURBU','NPOIN',dummi,dummi,gesca)
     end if

  case ( 51_ip )
     !
     ! LNSUB: Element subdomain
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              gesca(ipoin) = max(gesca(ipoin),real(lesub(ielem),rp))
           end do
        end do
        call pararr('SMA',NPOIN_TYPE,npoin,gesca)
     end if

  case ( 52_ip )
     !
     ! Wet nodes
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do icoup = 1,mcoup
           do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
              ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
              if( ipoin /= 0 ) gesca(ipoin) = real(icoup,rp)
           end do
        end do
     end if

  case ( 53_ip )
     !
     ! LESUB: Element Subdomains
     !
     if( INOTMASTER ) then
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = max(gesca(ielem),real(lesub(ielem),rp))
        end do
     end if

  case ( 54_ip )
     !
     ! CANOPY HEIGHT
     !
     if( kfl_canhe < 1 ) return
     gesca => canhe

  case ( 55_ip )
     !
     ! HEIGHT OVER TERRAIN
     !
     if( kfl_heiov < 1 ) return
     gesca => heiov

  case ( 56_ip )
     !
     ! BEATR: Element subdomain + characteristic + material
     !
     if( INOTMASTER ) then
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           icont = 100 * lmate(ielem) + 10 * lesub(ielem) + lelch(ielem)
           gesca(ielem) = max(gesca(ielem),real(icont,rp))
        end do
     end if

  case ( 57_ip )
     !
     ! COLORING of the openmp stratgy
     !
     if( INOTMASTER ) then
        allocate(list_colors(nelem))
        call memgen(zero,nelem,zero)
        do icolo = 1,par_omp_num_colors
           do kelem = par_omp_ia_colors(icolo),par_omp_ia_colors(icolo+1)-1
              ielem = par_omp_ja_colors(kelem)
              list_colors(ielem) = icolo
           end do
        end do
        do ielem = 1,nelem
           gesca(ielem) = real(list_colors(ielem),rp)
        end do
        !do ipoin = 1,npoin
        !   gesca(ipoin) = huge(1.0_rp)
        !end do
        !do ielem = 1,nelem
        !   do inode = 1,lnnod(ielem)
        !      ipoin = lnods(inode,ielem)
        !      gesca(ipoin) = min(gesca(ipoin),real(list_colors(ielem),rp))
        !   end do
        !end do
        !call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MIN','IN MY CODE')
        deallocate(list_colors)
     end if

  case ( 58_ip )
     !
     ! WALLN
     !
     gevec => walln

  case ( 59_ip )
     !
     ! Local numbering
     ! what was here originally did not work well
     ! now I leave 2 option that both work well
     !
     if(1==1) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) = real(ipoin,rp)
           end do
        end if
     else if(1==2) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              ! esto da mal
              !gesca(ipoin) = real(ipoin,rp)
           end do
           do ipoin = 1,meshe(ndivi)%npoi1
              gesca(ipoin) = real(ipoin,rp)
           end do
           do ipoin = meshe(ndivi)%npoi1+1 , meshe(ndivi)%npoin
              gesca(ipoin) = real(PAR_COMM_MY_CODE_ARRAY(1) % node_number_in_owner(ipoin-meshe(ndivi)%npoi1),rp)
           end do
        end if
     else
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoi1
              gesca(ipoin) = real(ipoin,rp)
           end do
           do ipoin = npoi1+1,npoin
              gesca(ipoin) = -1.0_rp
           end do
           do ipoin = npoi2,npoi3
              gesca(ipoin) = real(ipoin,rp)
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
        end if
     end if

  case ( 60_ip )
     !
     ! Advection
     !
     gevec => advec(:,:,1)

  case ( 61_ip )
     !
     ! Interface renumbering
     !
!!$     if( INOTMASTER ) then
!!$       call memgen(zero,npoin,zero)
!!$       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_MY_CODE_WM4,PAR_CURRENT_RANK4,PAR_CURRENT_SIZE4)
!!$       allocate( num_interface_nodes(PAR_CURRENT_SIZE4) )
!!$       call PAR_ALLGATHER(meshe % npoi3-meshe % npoi2,num_interface_nodes,1_ip,'IN MY CODE WITHOUT MASTER')
!!$       ii = num_interface_nodes(1)
!!$       num_interface_nodes(1) = 0
!!$       do ipart = 2,PAR_CURRENT_SIZE4
!!$          num_interface_nodes(ipart) = num_interface_nodes(ipart) + ii
!!$          ii = num_interface_nodes(ipart)
!!$       end do
!!$       kpoin = num_interface_nodes(PAR_CURRENT_RANK4)
!!$       do ipoin = meshe % npoi2, meshe % npoi3
!!$          kpoin = kpoin + 1
!!$          gisca(ipoin) = kpoin
!!$       end do
!!$       call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX','IN MY CODE')
!!$       deallocate( interface_nodes )
!!$    end if

  case ( 62_ip )
     !
     ! RENEL: local element renumbering
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = real(ielem,rp)
        end do
        !call memgen(0_ip,nelem,zero)
        !allocate(permr(nelem))
        !call renumbering_elements(1_ip,1_ip,nelem,meshe(ndivi),permr)
        !gesca(1:nelem) = real(permr(1:nelem),rp)
        !deallocate(permr)
     end if

  case ( 63_ip )
     !
     ! RENPO: local node renumbering
     !
     if( INOTMASTER ) then
        call memgen(0_ip,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(ipoin,rp)
        end do
     end if

  case ( 64_ip )
     !
     ! RENPO: local node renumbering
     !
     if( INOTMASTER ) then
        call memgen(0_ip,npoin,zero)
        allocate(xvalu(npoin))
        allocate(xcoo1(npoin))
        allocate(xcoo2(npoin))

        do ipoin = 1,npoin
           !xvalu(ipoin) = cos(2.0_rp*pi*coord(1,ipoin)/0.4_rp) + 0.2_rp * rand()
        end do
        call filters_nodal(xvalu,gesca)

        do ipoin = 1001,1,-1
           write(97,*) coord(1,ipoin),xvalu(ipoin),gesca(ipoin)
        end do
        stop
        !do ipoin = 1,npoin
        !   !xvalu(ipoin) = cos(2.0_rp*pi*coord(1,ipoin)/0.4_rp)
        !   xvalu(ipoin) = 3.0_rp*(coord(1,ipoin)-0.5_rp)**2
        !end do
        !call memgen(0_ip,ndime,npoin)
        !call gradie(xvalu,gevec)
        !!call filters_nodal(xvalu,gesca)
        !deallocate(xvalu)

     end if

  case ( 65_ip )
     !
     ! OMPSS element subdomains
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,zero)
        do isubd = 1, size(ompss_domains)
           do kelem = 1,size(ompss_domains(isubd)%elements)
              ielem = ompss_domains(isubd)%elements(kelem)
              gesca(ielem) = real(isubd,rp)
           end do
        end do
        !
        ! Mark separator
        !
        ielpo = 0
!!$           if( .not. associated(lelel) ) then
!!$              call graphs_eleele(&
!!$                   nelem,npoin,mnode,mepoi,lnods,lnnod,&
!!$                   pelpo,lelpo,nedg1,medg1,pelel,lelel)
!!$           end if
!!$           do ielem = 1,nelem
!!$              do ielpo = pelel(ielem),pelel(ielem+1)-1
!!$                 kelem = lelel(ielpo)
!!$                 if( gesca(ielem) /= gesca(kelem) .and. gesca(ielem) > 0.0_rp .and. gesca(kelem) > 0.0_rp ) then
!!$                    if( gesca(ielem) < gesca(kelem) ) then
!!$                       gesca(ielem) = -abs(gesca(ielem))
!!$                    else
!!$                       gesca(kelem) = -abs(gesca(kelem))
!!$                    end if
!!$                 end if
!!$              end do
!!$           end do
!!$
!!$           ipoin = 1
!!$           dummi = 0
!!$           do while( ipoin > 0 )
!!$              dummi = dummi + 1
!!$              print*,'passs=',dummi
!!$              ipoin = 0
!!$              do ielem = 1,nelem
!!$                 if( gesca(ielem) < 0.0_rp ) then
!!$                    inode = 0
!!$                    do ielpo = pelel(ielem),pelel(ielem+1)-1
!!$                       kelem = lelel(ielpo)
!!$                       if( abs(gesca(ielem)) /= gesca(kelem) .or. gesca(kelem) > 0.0_rp ) then
!!$                          inode = 1
!!$                       end if
!!$                    end do
!!$                    if( inode == 0 ) then
!!$                       ipoin = ipoin + 1
!!$                       gesca(ielem) = abs(gesca(ielem))
!!$                    end if
!!$                 end if
!!$              end do
!!$           end do
!!$
!!$           gesca = max(gesca,0.0_rp)

        if( ielpo /= 0 ) then
           call graphs_dealep(pelel,lelel)
        end if

     end if

  case ( 66_ip )
     !
     ! OPENMP element chunks
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,zero)
        ielem = 0
        element_loop: do
           do ithre = 1,par_omp_num_threads
              do kelem = 1,par_omp_nelem_chunk
                 ielem = ielem + 1
                 if( ielem > nelem ) then
                    exit element_loop
                 else
                    gesca(ielem) = real(ithre,rp)
                 end if
              end do
           end do
        end do element_loop
     end if

  case ( 67_ip )
     !
     ! LGAUS: number of Gauss points
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = real(lgaus(ielem),rp)
        end do
     end if

  case ( 68_ip )
     !
     ! WALLD
     !
     if( kfl_walld == 2 .or. kfl_walld == 3) then
        do ipoin = 1,npoin
           gesca(ipoin) = real(wallo(ipoin),rp)
        end do
     endif

  case ( 69_ip )
     !
     ! VISCA: air viscosity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('MIXIN','NPOIN',dummi,dummi,gesca)
     end if

  case ( 70_ip )
     !
     ! FIELD
     !
     ! only if no prealoading is enabled
     do ifiel = 1,nfiel

        if (kfl_field(6,ifiel) /= 1 ) then !not on demand


            wfiel =  intost(ifiel)
            if(ifiel<10) then
               wopos(1)=postp(1)%wopos(1,70)(1:3)//'0'//trim(wfiel(1:2))
            else
               wopos(1)=postp(1)%wopos(1,70)(1:3)//trim(wfiel(1:2))
            end if
            if( kfl_field(1,ifiel) > 0 ) then
               if(      kfl_field(2,ifiel) == NPOIN_TYPE ) then
                  nenti    = npoin
                  wopos(3) = 'NPOIN'
               else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
                  nenti    = nboun
                  wopos(3) = 'NBOUN'
               else if( kfl_field(2,ifiel) == NELEM_TYPE ) then
                  nenti    = nelem
                  wopos(3) = 'NELEM'
               end if
               kdime = kfl_field(1,ifiel)
               if(      kdime == 1     ) then
                  wopos(2) = 'SCALA'
               else if( kdime == ndime ) then
                  wopos(2) = 'VECTO'
               else
                  kdime = 0
               end if
               
               if( kdime == 1 ) then
                  if( INOTMASTER ) then
                     call memgen(zero,nenti,0_ip)
                     do ienti = 1,nenti
                        gesca(ienti) = xfiel(ifiel) % a(1,ienti,1)
                     end do
                  end if
                  call postpr(gesca,wopos,ittim,rutim)
                  if( INOTMASTER ) call memgen(two,nenti,0_ip)
                  
               else if( kdime == ndime ) then
                  
                  if( INOTMASTER ) then
                     call memgen(zero,ndime,nenti,0_ip)
                     do ienti = 1,nenti
                        gevec(1:ndime,ienti) = xfiel(ifiel) % a(1:ndime,ienti,1)
                     end do
                  end if
                  call postpr(gevec,wopos,ittim,rutim)
                  if( INOTMASTER ) call memgen(two,ndime,nenti)
                  
               end if
               
            end if
        end if !not on demand
     end do
     
     return

  case ( 71_ip )
     !
     ! VELOM
     !
     if( INOTMASTER ) then 
        call memgen(zero,ndime,npoin)
        if( associated(velom) ) then
           gevec(1:ndime,1:npoin) = velom(1:ndime,1:npoin)
        else
           gevec = 0.0_rp
        end if
     end if
     
  case ( 72_ip )
     !
     ! OMPSS boundary subdomains
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nboun,zero)
        do isubd = 1, size(ompss_boundaries)
           do kelem = 1,size(ompss_boundaries(isubd)%elements)
              ielem = ompss_boundaries(isubd)%elements(kelem)
              gesca(ielem) = real(isubd,rp)
           end do
        end do
     end if
     
  case ( 73_ip )
     !
     ! NSWVI: nodal projection of no slip wall  viscosity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        if( kfl_noslw_ker /= 0 ) then
           call ker_nod_nsw_visc_0
           gesca = 0.0_rp
           call projec_elements_to_nodes(el_nsw_visc,gesca)
        end if
     end if

  case ( 74_ip )
     !
     ! NUMBERING
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(ipoin,rp)
        end do
     end if
     
  case (75_ip) 
        ! HEIGHT OVER TERRAIN
     !
     if( kfl_canla < 1 ) return
     gesca => canla
 
  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1) % wopos(:,ivari))

end subroutine ker_outvar
