subroutine par_boufac()
  !-----------------------------------------------------------------------
  !****f* Parall/par_boufac
  ! NAME
  !    par_boufac
  ! DESCRIPTION
  !    This subroutine determines neighbors that share the same edge 
  !    IEDGG (IEDGB)
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall, only : commd
  use mod_parall, only : PAR_COMM_MY_CODE4,PAR_INTEGER
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
  integer(4)           :: status(MPI_STATUS_SIZE)
#endif
  integer(4)           :: bsizer,bsizes,istat,bsize4
  integer(ip)          :: ipoin,ii,jj,bsize,dom_i,iowns,iothe
  integer(ip)          :: ifacg,ineig,kk,node1,node2,node3,node4
  integer(ip)          :: jpoin,kpoin,mpoin,ipoi3,ipoi4,t,i,j,kowne
  integer(ip)          :: ipoi1,ipoi2,lsize,knode,ll
  integer(ip), pointer :: nneis(:),nneir(:),giscl(:),lfnei(:)
  type(i1p),   pointer :: lneis(:),lneir(:),lfass(:)

  if( INOTMASTER ) then
     !
     ! Candidate boundary faces
     !
     allocate( lfnei(nfacg) , stat = istat )
     do ifacg = 1,nfacg
        lfnei(ifacg) = 0
        ipoin = lfacg(1,ifacg)
        jpoin = lfacg(2,ifacg)
        kpoin = lfacg(3,ifacg)
        mpoin = lfacg(4,ifacg)
        if( ipoin > npoi1 .and. jpoin > npoi1 .and. kpoin > npoi1 .and. mpoin > npoi1 ) then
           lfacg(6,ifacg) = 1
        end if
     end do
     !
     ! Count my edges with my neighbors   
     ! My own candidate edges: LNEIS
     ! My neighbor's candidate edges: LNEIR
     !
     allocate( nneis(nneig) , stat = istat )
     allocate( nneir(nneig) , stat = istat )
     allocate( lneis(nneig) , stat = istat )
     allocate( lneir(nneig) , stat = istat )
     allocate( lfass(nneig) , stat = istat )

     do ineig = 1,nneig
        nneir(ineig) = 0
        nneis(ineig) = 0
     end do

     do ifacg = 1,nfacg
        if( lfacg(4,ifacg) > 0 ) then
           ipoin = lfacg(1,ifacg)
           jpoin = lfacg(2,ifacg)
           kpoin = lfacg(3,ifacg)
           mpoin = lfacg(4,ifacg)
           do ineig = 1,nneig
              knode = 0
              jj    = commd % bound_size(ineig)
              dom_i = commd % neights(ineig)
              do while( jj <= commd%bound_size(ineig+1)-1 .and. knode < 4 )
                 if(       commd % bound_perm(jj) == ipoin &
                      .or. commd % bound_perm(jj) == jpoin &
                      .or. commd % bound_perm(jj) == kpoin &
                      .or. commd % bound_perm(jj) == mpoin ) then
                    knode = knode + 1
                 end if
                 jj = jj + 1
              end do
              if( knode == 4 ) then
                 nneis(ineig) = nneis(ineig) + 1
              end if
           end do
        end if
     end do
     !
     ! List of candidate common faces with my neighbors
     !    
     do ineig = 1,nneig
        allocate( lneis(ineig) % l(max(1_ip,4*nneis(ineig))) , stat = istat )
        allocate( lfass(ineig) % l(max(1_ip,  nneis(ineig))) , stat = istat )
        nneis(ineig) = 0
        nneir(ineig) = 0
     end do

     do ifacg = 1,nfacg
        if( lfacg(4,ifacg) > 0 ) then
           ipoin = lfacg(1,ifacg)
           jpoin = lfacg(2,ifacg)
           kpoin = lfacg(3,ifacg)
           mpoin = lfacg(4,ifacg)
           node1 = lninv_loc(ipoin)
           node2 = lninv_loc(jpoin)
           node3 = lninv_loc(kpoin)
           node4 = lninv_loc(mpoin)
           do ineig = 1,nneig
              knode = 0
              jj    = commd % bound_size(ineig)
              dom_i = commd % neights(ineig)
              do while( jj <= commd%bound_size(ineig+1)-1 .and. knode < 4 )
                 if(         commd % bound_perm(jj) == ipoin &
                      & .or. commd % bound_perm(jj) == jpoin &
                      & .or. commd % bound_perm(jj) == kpoin &
                      & .or. commd % bound_perm(jj) == mpoin ) then
                    knode = knode + 1
                 end if
                 jj = jj + 1
              end do
              if( knode == 4 ) then
                 nneis(ineig) = nneis(ineig) + 1
                 lfass(ineig) % l(nneis(ineig)) = ifacg
                 lneis(ineig) % l( (nneis(ineig)-1)*4+1 ) = node1
                 lneis(ineig) % l( (nneis(ineig)-1)*4+2 ) = node2         
                 lneis(ineig) % l( (nneis(ineig)-1)*4+3 ) = node3
                 lneis(ineig) % l( (nneis(ineig)-1)*4+4 ) = node4         

                 ii = (nneis(ineig)-1)*4+1
                 call sortin(4_ip,lneis(ineig) % l(ii))

              end if
           end do
        end if
     end do

#ifdef MPI_OFF
#else
     !
     ! Send/receive # my common edges to/from my neighbors
     !
     do ineig = 1, nneig

        dom_i  = commd%neights(ineig)
        bsize4 = int(1_ip,4)

        call MPI_Sendrecv( &
             nneis(ineig:), bsize4,        &
             PAR_INTEGER,  dom_i, 0_4,     &
             nneir(ineig:), bsize4,        &
             PAR_INTEGER,  dom_i, 0_4,     &
             PAR_COMM_MY_CODE4, status, istat )
     end do
#endif

     do ineig = 1,nneig
        lsize = max(1_ip,nneir(ineig))
        allocate( lneir(ineig) % l(4*lsize) , stat = istat )
     end do

#ifdef MPI_OFF
#else
     !
     ! Send/receive list of candidate common faces to/from my neighbors
     !
     do ineig = 1, nneig

        dom_i  = commd%neights(ineig)
        bsizes = int(max(1_ip,nneis(ineig)*4),4)
        bsizer = int(max(1_ip,nneir(ineig)*4),4) 

        call MPI_Sendrecv( &
             lneis(ineig) % l(1:), bsizes, &
             PAR_INTEGER,  dom_i, 0_4,     &
             lneir(ineig) % l(1:), bsizer, &
             PAR_INTEGER,  dom_i, 0_4,     &
             PAR_COMM_MY_CODE4, status, istat )

     end do
#endif

     !
     ! Eliminate edges not shared by any other subdomain: can be an interior edge
     ! NFACB: Number of boundary faces
     !
     nfacb = 0
     do ineig = 1,nneig
        dom_i = commd%neights(ineig)
        do ii = 1,nneis(ineig)
           ifacg = lfass(ineig) % l(ii)
           if( lfacg(4,ifacg) > 0 ) then
              ipoin = lneis(ineig) % l( (ii-1)*4+1 )
              jpoin = lneis(ineig) % l( (ii-1)*4+2 )
              kpoin = lneis(ineig) % l( (ii-1)*4+3 )
              mpoin = lneis(ineig) % l( (ii-1)*4+4 )
              knode = 0
              kk    = 0
              do while( kk < nneir(ineig) )
                 kk    = kk + 1
                 node1 = lneir(ineig) % l( (kk-1)*4+1 )
                 node2 = lneir(ineig) % l( (kk-1)*4+2 )
                 node3 = lneir(ineig) % l( (kk-1)*4+3 )
                 node4 = lneir(ineig) % l( (kk-1)*4+4 )
                 if( node1 == ipoin .and. node2 == jpoin .and. node3 == kpoin .and. node4 == mpoin ) then
                    knode = 4
                    kk = nneir(ineig)
                 end if
              end do
              if( knode == 4 ) then
                 !
                 ! Neighbor INEIG has also edge IFACG
                 !
                 nfacb          =  nfacb + 1
                 lfacg(6,ifacg) = -nfacb
                 lfnei(ifacg)   =  ineig
              end if
           end if
        end do
     end do
     do ifacg = 1,nfacg
        if( lfacg(6,ifacg) > 0 ) then
           lfacg(6,ifacg) =  0              ! Candidate face is not good
        else
           lfacg(6,ifacg) = -lfacg(6,ifacg) ! Candidate face is shared by another subdomain
        end if
     end do
     !
     ! Allocate memory for boundary edges
     !
     allocate( lfacb(6+1,nfacb), stat=istat )
     call memchk(zero,istat,mem_servi(1:2,servi),'LFACB','par_submsh',lfacb) 
     !
     ! Order boundary nodes
     !
     bsize = 0
     do ineig = 1,nneig
        bsize = bsize + commd%bound_size(ineig+1)-commd%bound_size(ineig)
     end do
     call memgen(1_ip,bsize,0_ip)
     allocate(giscl(bsize))
     ii = 0
     kk = 1
     do ineig = 1, nneig
        do jj = commd%bound_size(ineig),commd%bound_size(ineig+1)-1
           ii = ii + 1
           gisca(ii) = commd % bound_perm(ii)
           giscl(ii) = ii
        end do
        bsize = commd%bound_size(ineig+1)-commd%bound_size(ineig)
        call heapsorti2(2_ip,bsize,gisca(kk),giscl(kk))
        kk = kk + bsize 
     end do
     !
     ! Fill in boundary face
     !
     nfacb = 0
     do ifacg = 1,nfacg
        if( lfacg(6,ifacg) > 0 ) then
           nfacb          = nfacb + 1
           ipoi1          = lfacg(1,ifacg)
           ipoi2          = lfacg(2,ifacg)
           ipoi3          = lfacg(3,ifacg)
           ipoi4          = lfacg(4,ifacg)
           lfacb(1,nfacb) = ifacg
           lfacb(2,nfacb) = ipoi1          ! Local face node
           lfacb(3,nfacb) = ipoi2          ! Local face node
           lfacb(4,nfacb) = ipoi3          ! Local face node
           lfacb(5,nfacb) = ipoi4          ! Local face node
           do i = 1,4-1
              do j = i+1,4
                 if( lfacb(i+1,nfacb) > lfacb(j+1,nfacb) ) then
                    t                = lfacb(i+1,nfacb)
                    lfacb(i+1,nfacb) = lfacb(j+1,nfacb)
                    lfacb(j+1,nfacb) = t
                 end if
              end do
           end do
           ipoi1 = lfacb(2,nfacb)
           ipoi2 = lfacb(3,nfacb)
           ipoi3 = lfacb(4,nfacb)
           ipoi4 = lfacb(5,nfacb)
           ineig = lfnei(ifacg)
           dom_i = commd % neights(ineig)
           
           iowns = 0
           iothe = 0
           
           ii    = commd%bound_size(ineig) - 1  
           kk    = 0

           do while( ii < commd%bound_size(ineig+1)-1 )
              ii = ii + 1
              if(      kk == 0 .and. ipoi1 == gisca(ii) ) then
                 kk    = kk + 1
                 ll    = giscl(ii)
                 iothe = iothe + lownr_par(ll)
                 iowns = iowns + lowns_par(ll)
              else if( kk == 1 .and. ipoi2 == gisca(ii) ) then
                 kk    = kk + 1
                 ll    = giscl(ii)
                 iothe = iothe + lownr_par(ll)
                 iowns = iowns + lowns_par(ll)
              else if( kk == 2 .and. ipoi3 == gisca(ii) ) then
                 kk    = kk + 1
                 ll    = giscl(ii)
                 iothe = iothe + lownr_par(ll)
                 iowns = iowns + lowns_par(ll)
              else if( kk == 3 .and. ipoi4 == gisca(ii) ) then
                 kk    = kk + 1
                 ll    = giscl(ii)
                 iothe = iothe + lownr_par(ll)
                 iowns = iowns + lowns_par(ll)
                 ii    = commd%bound_size(ineig+1)
              end if
           end do
           !
           ! Order boundary face nodes
           !
           lfacb(2,nfacb) = lninv_loc(ipoi1)          ! global face node
           lfacb(3,nfacb) = lninv_loc(ipoi2)          ! global face node
           lfacb(4,nfacb) = lninv_loc(ipoi3)          ! global face node
           lfacb(5,nfacb) = lninv_loc(ipoi4)          ! global face node
           do i = 1,4-1
              do j = i+1,4
                 if( lfacb(i+1,nfacb) > lfacb(j+1,nfacb) ) then
                    t                = lfacb(i+1,nfacb)
                    lfacb(i+1,nfacb) = lfacb(j+1,nfacb)
                    lfacb(j+1,nfacb) = t
                 end if
              end do
           end do
           lfacb(7,nfacb) = ineig
           !
           ! Decide who owns the node
           !
           if( iowns > iothe ) then
              kowne = -1
           else if( iothe > iowns ) then
              kowne =  1
           else
              if( dom_i < kfl_paral ) then
                 kowne = -1
              else
                 kowne =  1
              end if
           end if
           if( iowns + iothe > 4 ) then
              print*,'PAR_BOUFAC: OUPS=',ineig,iowns,iothe
              call runend('PAR_BOUFAC: OUPS')               
           end if
           lfacg(5,ifacg) = abs(lfacg(5,ifacg)) * kowne
        end if
     end do
     !
     ! Deallocate memory
     !
     deallocate(giscl)
     call memgen(3_ip,bsize,0_ip)

     do ineig = 1, nneig
        deallocate( lfass(ineig) % l , stat = istat )
        deallocate( lneir(ineig) % l , stat = istat )
        deallocate( lneis(ineig) % l , stat = istat )
     end do

     deallocate( lfass , stat = istat )
     deallocate( lneir , stat = istat )
     deallocate( lneis , stat = istat )
     deallocate( nneis , stat = istat )
     deallocate( nneir , stat = istat )
     deallocate( lfnei , stat = istat )

  end if


end subroutine par_boufac
