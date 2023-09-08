!-----------------------------------------------------------------------
!
!> @addtogroup Domain
!> @{
!> @name    ToolBox for periodicity
!> @file    mod_periodicity.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox periodicity
!> @details ToolBox periodicity
!> @{
!
!-----------------------------------------------------------------------

module mod_periodicity

  use def_kintyp,    only : ip,rp
  use def_parame
  use def_master
  use def_domain
  use mod_memory
  use mod_messages, only : livinf

  implicit none

  integer(ip), pointer :: lnods_sav(:,:) 
  integer(ip)          :: mnode_sav

  public :: periodicity_preprocess
  public :: periodicity_postprocess
  public :: periodicity_enhance_node_graph
  public :: periodicity_setup

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   This routine modifies LNODS to account for periodicity
  !> @date    6/04/2017
  !> @details Periodicity is implemented as follows:
  !>          \verbatim
  !>
  !>          Input:  lnods(1:mnode,1:nelem)
  !>          Output: lnods(1:2*mnode,nelem)
  !>
  !>          1. nperi= numero de master-slave couples
  !>             NB: a master can have various slaves
  !>          2. lperi(2,1:nperi) = list of slave-master couples
  !>              lperi(1,iperi) = master
  !>              lperi(2,iperi) = slave
  !>          3. If ielem is an element containing a slave node:
  !>             add the corresponding master to the element connectivity
  !>             of ielem: lnods(1:2*mnode:ielem)
  !>             Examples: lnods(1:3,ielem) = 23,27,31
  !>                       lperi(1,56) = 127
  !>                       lperi(2,56) = 23
  !>             => lnods(1:4,ielem) = 23,27,31,127
  !>          4. Construct node-node graph (R_DOM,C_DOM) and node-element 
  !>             graph (PELPO,LELPO) using this modified lnods
  !>          5. Recover original lnode
  !>          6. After calling METIS, each time a slave is in a subdomain,
  !>             then the master will. Note that the partitioning will reflect
  !>             the presence of periodicity only if the partition is carried
  !>             out using the common-node criterion (not the face one) as
  !>             it uses LELPO and PELPO
  !>          7. During the assembly of a matrix in parallel, put all
  !>             the coefficients related to a slave to the master position
  !>             (row and columns). Do the same for RHS.
  !>          8. The slave wil have a null row and column. At the end
  !>             of the solver put the slave (jpoin) value to its master 
  !>             counterpart (ipoin): xx(jpoin) = xx(ipoin)
  !>
  !>          \endverbatim
  !-----------------------------------------------------------------------

  subroutine periodicity_preprocess()

    integer(ip)          :: ipoin,ielem,inode,jpoin,iperi,kperi
    integer(ip)          :: knode,ii,jnode
    logical(lg)          :: ifoun
    !integer(ip), pointer :: list_master(:)

    if( nperi > 0 .and. IMASTER ) then
       !
       ! Save MNODE
       !
       mnode_sav = mnode       
       !
       ! Compute LNNOD, number of nodes per element required to construct the graphs
       !
       nullify(lnods_sav)
       if( associated(lnnod) ) call runend('PERIODICITY: LNNOD ALREADY ASSOCIATED!')
       call memory_alloca(memor_dom,'LNNOD','periodicity',lnnod,nelem)
       do ielem = 1,nelem
          lnnod(ielem) = nnode(abs(ltype(ielem)))
       end do
       !
       ! Increase MNODE to receive master nodes
       !
       do ielem = 1,nelem
          knode = lnnod(ielem)
          do inode = 1,lnnod(ielem)
             ipoin = lnods(inode,ielem)
             if( lmast(ipoin) > 0 ) knode = knode + 1
          end do
          mnode = max(mnode,knode)
       end do
       !
       ! Reallocate LNODS
       !       
       call memory_copy(  memor_dom,'LNODS','periodicity',lnods,lnods_sav)       
       call memory_deallo(memor_dom,'LNODS','periodicity',lnods)
       call memory_alloca(memor_dom,'LNODS','periodicity',lnods,mnode,nelem)

       do ielem = 1,nelem
          do inode = 1,lnnod(ielem)
             lnods(inode,ielem) = lnods_sav(inode,ielem)
          end do
       end do
       !
       ! Recompute LNODS
       !
       do ielem = 1,nelem
          knode = lnnod(ielem)
          do inode = 1,lnnod(ielem)
             ipoin = lnods_sav(inode,ielem)
             if( lmast(ipoin) > 0 ) then
                ifoun = .false.
                do jnode = 1,lnnod(ielem)
                   if( lnods_sav(jnode,ielem) == lmast(ipoin) ) ifoun = .true.
                end do
                if( .not. ifoun ) then
                   knode = knode + 1
                   lnods(knode,ielem) = lmast(ipoin)
                end if
             end if
          end do
          lnnod(ielem) = knode
       end do

    end if

  end subroutine periodicity_preprocess

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   This routine modifies LNODS to account for periodicity
  !> @date    6/04/2017
  !> @details Recover original arrays
  !>
  !-----------------------------------------------------------------------
  
  subroutine periodicity_postprocess()

    integer(ip) :: ielem

    if( nperi > 0 .and. IMASTER ) then
       !
       ! Recover MNODE and LNODS
       !
       mnode =  mnode_sav

       call memory_deallo(memor_dom,'LNODS'    ,'periodicity_postprocess',lnods)
       call memory_copy(  memor_dom,'LNODS'    ,'periodicity_postprocess',lnods_sav,lnods)
       call memory_deallo(memor_dom,'LNODS_SAV','periodicity_postprocess',lnods_sav)
       call memory_deallo(memor_dom,'LLNOD'    ,'periodicity_postprocess',lnnod)

    end if

  end subroutine periodicity_postprocess

  !-----------------------------------------------------------------------
  !>
  !> @author   Guillaume Houzeaux
  !> @date     18/09/2012
  !> @brief    Change graph to account for periodicity
  !> @keywords periodicity,graph
  !> @details  Changes the graphs to account for periodicity by adding
  !>           the neighbors of the slave to the master.
  !>           LPERI(1,1:NPERI) = list of masters
  !>           LPERI(2,1:NPERI) = list of slaves
  !>           A master can have up to 4 slaves in 2D and 8 slaves in 3D
  !>
  !>           \verbatim
  !>
  !>                  +------------+----------------+
  !>                  |      CPU1  |                |
  !>                  o---o        |            x---x
  !>                  |   |        |            |   |
  !>           Master M---o--------+   CPU3     x---S Slave
  !>                  |   |        |            |   | M Master
  !>                  o---o        |            x---x
  !>                  |      CPU2  |                |
  !>                  +------------+----------------+
  !>
  !>           \endverbatim
  !>
  !>           o are neighbors of Master M in respective CPUs \n
  !>           x are neighbors of Slave  S in respective CPUs \n
  !>           By construction, the Master is present in all the CPUs where 
  !>           the slave is. Therefore, the row of the Slave can always
  !>           be copied to that of the Master.
  !>
  !> @} 
  !-----------------------------------------------------------------------

  subroutine periodicity_enhance_node_graph(onwhat)
    
    use def_kintyp, only : ip,i1p
    use def_domain, only : npoin
    use def_domain, only : r_dom
    use def_domain, only : r_dom_own
    use def_domain, only : r_sol
    use def_domain, only : c_dom
    use def_domain, only : c_dom_own
    use def_domain, only : c_sol
    use def_domain, only : nzdom
    use def_domain, only : nzdom_own
    use def_domain, only : nzsol
    use def_domain, only : nzsym
    use def_domain, only : lmast
    use def_domain, only : meshe
    use def_domain, only : memor_dom
    use def_master, only : INOTMASTER
    use def_kermod, only : ndivi
    use mod_memory, only : memory_alloca
    use mod_memory, only : memory_deallo
    use mod_memory, only : memory_resize
    implicit none
    character(*), intent(in)            :: onwhat
    integer(ip)                         :: ipoin,izdom,jpoin,jzdom,kpoin,ifoun
    integer(ip)                         :: kzdom,madde,lsize,isize
    integer(ip), pointer                :: lcdom_size(:) 
    integer(ip), pointer                :: list_of_masters(:) 
    type(i1p),   pointer                :: lcdom(:)
    
    integer(ip)                         :: npoin_loc
    integer(ip)                         :: nzdom_loc
    integer(ip), pointer                :: r_dom_loc(:)
    integer(ip), pointer                :: c_dom_loc(:)
    
    nullify(lcdom)
    nullify(lcdom_size)
    nullify(list_of_masters)

    if( nperi > 0 ) then
       !
       ! On which graph should we impose periodicity
       !
       if( onwhat == 'LOCAL GRAPH' ) then
          npoin_loc =  npoin
          nzdom_loc =  nzdom
          r_dom_loc => r_dom
          c_dom_loc => c_dom
       else if( onwhat == 'FULL ROWS GRAPH' ) then
          npoin_loc =  npoin_own
          nzdom_loc =  nzdom_own
          r_dom_loc => r_dom_own
          c_dom_loc => c_dom_own         
       end if
       
       call livinf(0_ip,'RECOMPUTE GRAPH FOR PERIODICITY',0_ip)

       if( INOTMASTER ) then
          !
          ! Allocate memory for LCDOM: temporary graph
          !
          call memory_alloca(memor_dom,'LCDOM'     ,'cresla',lcdom,     npoin_loc)
          call memory_alloca(memor_dom,'LCDOM_SIZE','cresla',lcdom_size,npoin_loc)
          do ipoin = 1,npoin_loc
             lsize = r_dom_loc(ipoin+1)-r_dom_loc(ipoin)
             lcdom_size(ipoin) = lsize
             call memory_alloca(memor_dom,'LCDOM(IPOIN) % L','cresla',lcdom(ipoin) % l,lsize)
             lsize = 0
             do izdom = r_dom_loc(ipoin),r_dom_loc(ipoin+1)-1
                lsize = lsize + 1
                lcdom(ipoin) % l(lsize) = c_dom_loc(izdom)
             end do
          end do
          !
          ! Put all slaves' columns into master's column
          ! KPOIN is neither master nor slave
          ! JPOIN is slave
          ! IPOIN is JPOIN's master
          !          
          do kpoin = 1,npoin_loc
             if( lmast(kpoin) == 0 ) then
                madde = 0
                do kzdom = r_dom_loc(kpoin),r_dom_loc(kpoin+1)-1
                   jpoin = c_dom_loc(kzdom)
                   if( lmast(jpoin) > 0 ) then
                      ipoin = lmast(jpoin)
                      ifoun = 0
                      lsize = lcdom_size(kpoin)
                      izdom1: do isize = 1,lsize
                         if( lcdom(kpoin) % l(isize) == ipoin )  then
                            ifoun = 1 
                            exit izdom1
                         end if
                      end do izdom1
                      if( ifoun == 0 ) then
                         call memory_resize(memor_dom,'LCDOM(KPOIN) % L','cresla',lcdom(kpoin) % l,lsize+1_ip)
                         lcdom_size(kpoin)         = lsize+1_ip
                         lcdom(kpoin) % l(lsize+1) = ipoin
                      end if
                   end if
                end do
             end if
          end do
          !
          ! Put all slaves' rows into master's row
          ! Put slave row to zero
          ! JPOIN:  slave
          ! IPOIN:  master
          ! Slave:  coef. JZDOM for JPOIN-KPOIN
          ! Master: coef. IZOMD for IPOIN-KPOIN
          !
          do jpoin = 1,npoin_loc
             if( lmast(jpoin) > 0 ) then
                madde = 0
                ipoin = lmast(jpoin)
                do jzdom = r_dom_loc(jpoin),r_dom_loc(jpoin+1)-1
                   kpoin = c_dom_loc(jzdom)     
                   if( lmast(kpoin) > 0 ) kpoin = lmast(kpoin)
                   ifoun = 0
                   lsize = lcdom_size(ipoin) 
                   izdom2: do isize = 1,lsize
                      if( lcdom(ipoin) % l(isize) == kpoin )  then
                         ifoun = 1 
                         exit izdom2
                      end if
                   end do izdom2
                   if( ifoun == 0 ) then
                      call memory_resize(memor_dom,'LCDOM(IPOIN) % L','cresla',lcdom(ipoin) % l,lsize+1_ip)
                      lcdom_size(ipoin)         = lsize+1_ip
                      lcdom(ipoin) % l(lsize+1) = kpoin
                   end if
                end do
             end if
          end do
          !
          ! Reallocate C_DOM and recompute NZDOM, R_DOM
          ! Reorder graph C_DOM 
          !
          nzdom_loc = 0
          do ipoin = 1,npoin_loc
             nzdom_loc = nzdom_loc + lcdom_size(ipoin)
          end do
             
          if( onwhat == 'LOCAL GRAPH' ) then
             
             nzdom = nzdom_loc
             call memory_alloca(memor_dom,'C_DOM','cresla',c_dom,nzdom,'REALLOCATE')
             c_dom_loc => c_dom
             
          else if( onwhat == 'FULL ROWS GRAPH' ) then
             
             nzdom_own = nzdom_loc
             call memory_alloca(memor_dom,'C_DOM_OWN','cresla',c_dom_own,nzdom_own,'REALLOCATE')
             c_dom_loc => c_dom_own
             
          end if
          
          r_dom_loc(1) = 1
          do ipoin = 2,npoin_loc+1  
             r_dom_loc(ipoin) = r_dom_loc(ipoin-1) + lcdom_size(ipoin-1)
          end do
          do ipoin = 1,npoin_loc
             isize = 0
             do izdom = r_dom_loc(ipoin),r_dom_loc(ipoin+1)-1
                isize = isize + 1
                c_dom_loc(izdom) = lcdom(ipoin) % l(isize)
             end do
             if( isize > 0 ) call heapsorti1(2_ip,isize,c_dom_loc(r_dom_loc(ipoin)))
          end do
          !
          ! Deallocate memory 
          !
          call memory_deallo(memor_dom,'LCDOM'     ,'cresla',lcdom)
          call memory_deallo(memor_dom,'LCDOM_SIZE','cresla',lcdom_size)   

       end if
       
    end if

    !-------------------------------------------------------------------
    !
    ! Solver dimensions and graphs
    ! 
    !-------------------------------------------------------------------

    if( onwhat == 'LOCAL GRAPH' ) then
       
       c_sol => c_dom
       r_sol => r_dom
       nzsol =  nzdom
       nzsym =  (nzsol-npoin)/2 + npoin
       
       meshe(ndivi) % r_dom     => r_dom  
       meshe(ndivi) % c_dom     => c_dom  
       meshe(ndivi) % nzdom     =  nzdom
       
    else if( onwhat == 'FULL ROWS GRAPH' ) then
       
       meshe(ndivi) % r_dom_own => r_dom_own
       meshe(ndivi) % c_dom_own => c_dom_own
       meshe(ndivi) % nzdom_own =  nzdom_own
       
    end if
       
  end subroutine periodicity_enhance_node_graph

  subroutine periodicity_setup()

    use mod_communications, only : PAR_ALLGATHER
    use mod_communications, only : PAR_ALLGATHERV
    
    integer(ip)          :: ipoin,iperi,kpoin
    integer(ip)          :: my_num_periodic_nodes
    integer(ip)          :: tot_num_periodic_nodes
    integer(ip), pointer :: my_list_periodic_nodes(:)
    integer(ip), pointer :: list_periodic_nodes(:)
    integer(ip), pointer :: num_periodic_nodes(:)
    !
    ! My list of masters
    !
    nullify(my_list_periodic_nodes)
    nullify(num_periodic_nodes)
    my_num_periodic_nodes = 0
    do ipoin = 1,npoin
       if( lmast(ipoin) /= 0 ) my_num_periodic_nodes = my_num_periodic_nodes + 1
    end do
    call memory_alloca(memor_dom,'MY_LIST_PERIODIC_NODES','periodicity_setup',my_list_periodic_nodes,my_num_periodic_nodes)
    call memory_alloca(memor_dom,'NUM_PERIODIC_NODES'    ,'periodicity_setup',num_periodic_nodes,npart+1,'INITIALIZE',0_ip)
    my_num_periodic_nodes = 0
    do ipoin = 1,npoin
       if( lmast(ipoin) /= 0 ) then
          my_num_periodic_nodes = my_num_periodic_nodes + 1
          my_list_periodic_nodes(my_num_periodic_nodes) = lmast(ipoin)
       end if
    end do
    !
    ! Gather list of masters
    !
    call PAR_ALLGATHER(my_num_periodic_nodes,num_periodic_nodes,1_4)
    tot_num_periodic_nodes = sum(num_periodic_nodes)
    
    call memory_alloca(memor_dom,'LIST_PERIODIC_NODES','periodicity_setup',list_periodic_nodes,tot_num_periodic_nodes)
    
    !call PAR_ALLGATHERV(my_list_periodic_nodes,list_periodic_nodes,tot_num_periodic_nodes)
    !
    ! Count
    !
    ipoin = 0
    do iperi = 1,size(list_periodic_nodes,KIND=ip)
       ipoin = 1
       kpoin = list_periodic_nodes(iperi)
       ipoin_loop: do while( ipoin < npoin )
          if( lninv_loc(ipoin) == kpoin ) then

          end if
       end do ipoin_loop
    end do
       
  end subroutine periodicity_setup
  
end module mod_periodicity
