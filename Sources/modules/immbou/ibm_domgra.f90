subroutine ibm_domgra()
  !-----------------------------------------------------------------------
  !****f* Domain/domgra
  ! NAME
  !    domgra
  ! DESCRIPTION
  !    This routine allocates memory for the coefficients of 
  !    the mesh graph.
  ! OUTPUT
  !    NZDOM ... Number of nonzero coefficients of the graph
  !    R_DOM ... Pointer to the array of rows r_dom(npoin+1) (r_dom(ipoin) = 
  !              coefficient of the graph where row ipoin starts)
  !    C_DOM ... Pointer to the array of columns c_dom(nzdom) (c_dom (izdom)
  !              = column of the izdom coefficient of mesh graph)
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use mod_memchk
  implicit none
  integer(ip)              :: ipoin,ielem,jelem,icoef,nzdom_2
  integer(ip)              :: izdom,ncoef,nlelp_2,mtouc,lsize
  integer(ip)              :: inode,jpoin
  integer(4)               :: istat
  integer(ip), allocatable :: lista(:),nepoi_2(:)!,lelpo_2(:),pelpo_2(:)
  logical(lg), allocatable :: touch(:)

  if( ISLAVE ) then
     !
     ! Allocate memory for NEPOI_2 and compute it
     !
     allocate(nepoi_2(npoin-npoi1),stat=istat)
     call memchk(zero,istat,memor_dom,'NEPOI_2','connpo',nepoi_2)
     do ielem = nelem+1,nelem_2
        do inode = 1,lnnod(ielem) 
           ipoin = lnods(inode,ielem)
           if (ipoin >= npoi1+1 .and. ipoin <= npoin) then
              nepoi_2(ipoin-npoi1)=nepoi_2(ipoin-npoi1)+1
           end if 
        end do
     end do
     !
     ! Allocate memory for PELPO_2 and compute it
     !
     allocate(pelpo_2(npoin-npoi1+1),stat=istat)
     call memchk(zero,istat,memor_dom,'PELPO_2','connpo',pelpo_2)
     pelpo_2(1)=1
     do jpoin = npoi1+1,npoin
        ipoin = jpoin - npoi1
        pelpo_2(ipoin+1)=pelpo_2(ipoin) + nepoi_2(ipoin)
     end do
     !
     ! Allocate memory for LELPO_2 and construct the list
     !
     nlelp_2 = pelpo_2(npoin-npoi1+1)
     allocate(lelpo_2(nlelp_2),stat=istat)
     call memchk(zero,istat,memor_dom,'LELPO_2','connpo',lelpo_2)
     do ielem = nelem+1,nelem_2
        do inode = 1,lnnod(ielem)
           ipoin = lnods(inode,ielem)
           if (ipoin >= npoi1+1 .and. ipoin <= npoin) then
              lelpo_2(pelpo_2(ipoin-npoi1)) = ielem
              pelpo_2(ipoin-npoi1) = pelpo_2(ipoin-npoi1)+1
           end if
        end do
     end do
     !
     ! Recompute PELPO_2 and maximum number of element neighbors MEPOI
     !
     pelpo_2(1) = 1
     mepoi = -1
     do jpoin = npoi1+1,npoin
        ipoin = jpoin - npoi1
        pelpo_2(ipoin+1)=pelpo_2(ipoin) + nepoi_2(ipoin)
        mepoi = max(mepoi,nepoi_2(ipoin))
     end do
     !
     ! Deallocate memory for temporary node/element connectivity
     !
     call memchk(two,istat,memor_dom,'NEPOI_2','connpo',nepoi_2)
     deallocate(nepoi_2,stat=istat)
     if(istat/=0) call memerr(two,'NEPOI_2','connpo',0_ip)
     !
     ! Compute node-node graph: C_DOM_2 and R_DOM_2
     !
     allocate(r_dom_2(npoin-npoi1+1),stat=istat)
     call memchk(zero,istat,memor_dom,'R_DOM_2','domgra',r_dom_2)

     mtouc = 0
     do ipoin = npoi1+1,npoin
        mtouc = max(mtouc,(pelpo_2(ipoin-npoi1+1)-pelpo_2(ipoin-npoi1))*mnode)
     end do
     nzdom_2 = 0

     allocate(touch(mtouc),stat=istat)
     call memchk(zero,istat,memor_dom,'TOUCH','domgra',touch)
     do ipoin = npoi1+1,npoin
        nlelp_2          = pelpo_2(ipoin-npoi1+1) - pelpo_2(ipoin-npoi1)
        ncoef          = nlelp_2 * mnode
        do icoef = 1,ncoef          
           touch(icoef) = .false.
        end do
        call ibm_nzecof(&
             nlelp_2,ncoef,nzdom_2,lelpo_2(pelpo_2(ipoin-npoi1)),touch,ipoin,npoin)
     end do

     !
     ! Construct the array of indexes
     !
     allocate(c_dom_2(nzdom_2),stat=istat)
     call memchk(zero,istat,memor_dom,'C_DOM_2','domgra',c_dom_2)
     izdom = 1
     do ipoin = npoi1+1,npoin
        nlelp_2 = pelpo_2(ipoin-npoi1+1) - pelpo_2(ipoin-npoi1)
        ncoef = nlelp_2 * mnode
        do icoef = 1,ncoef          
           touch(icoef) = .false.
        end do
        call ibm_arrind(&
             nlelp_2,ncoef,lelpo_2(pelpo_2(ipoin-npoi1)),touch,izdom,ipoin,npoin,&
             r_dom_2,c_dom_2)
     end do
     r_dom_2(npoin+1-npoi1) = nzdom_2 + 1
     call memchk(two,istat,memor_dom,'TOUCH','domgra',touch)
     deallocate(touch,stat=istat)
     if(istat/=0) call memerr(two,'TOUCH','domgra',0_ip) 

  end if

end subroutine ibm_domgra

subroutine ibm_nzecof(&
     nlist,ncoef,nzdom,liste,touch,ipoin,npoin)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the number of non-zero coefficients of a
  ! mesh graph stored in compressed sparse row (CSR) format 
  !                 
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  nelem,lnods,mnode,lnnod,r_dom,c_dom
  implicit none
  integer(ip), intent(in)    :: nlist,ncoef
  integer(ip), intent(in)    :: liste(nlist)
  integer(ip), intent(in)    :: ipoin,npoin
  integer(ip), intent(inout) :: nzdom
  logical(lg), intent(inout) :: touch(ncoef)
  integer(ip)                :: jelem,jnode,jpoin,nnodj,jposi,jlist
  integer(ip)                :: kelem,knode,kpoin,nnodk,kposi,klist,izdom
  logical(lg)                :: ifoun

  do jlist=1,nlist                                      ! Loop over those elements 
     jelem=liste(jlist)                                 ! where the point is
     nnodj=lnnod(jelem)
     do jnode=1,nnodj
        jpoin=lnods(jnode,jelem)
        if( jpoin /= ipoin ) then
           !
           ! In parallel, some nodes are not connected with their own neighbor nodes inside a subdomain.
           ! That is, these neighbor nodes are not stored in the data structure "r_dom"
           ! For this reason, the data structure "r_dom_2" has to store these neighbor nodes.  
           !
           ! Example:
           ! 1:         First  subdomian
           ! 2:         Second subdomian
           ! O:         A node
           ! A,B,C,D,E: Neighbor nodes of the node "O"
           !
           !     --------C-------D
           !    / \     / \     / \
           !   / 1 \ 1 / 2 \ 1 / 1 \ 
           !  /     \ /     \ /     \
           ! --------B-------O-------E
           !  \     / \     / \     /
           !   \ 1 / 1 \ 2 / 2 \ 1 /
           !    \ /     \ /     \ /
           !     --------A-------F
           !
           ! In this example, the node "O" is not connected with the nodes "A" and "B" inside the subdomain "1".
           ! This means that the nodes "A" and "B" were not added to the data structure "r_dom".
           ! The solution is add the nodes "A" and "B" to the data structure "r_dom_2".
           !           
           ifoun = .true.
           izdom = r_dom(ipoin)
           do while( izdom <= r_dom(ipoin+1)-1 )
              if( c_dom(izdom) == jpoin ) then
                 ifoun = .false.
                 izdom = r_dom(ipoin+1)
              end if
              izdom = izdom + 1              
           end do

           jposi=(jlist-1)*mnode+jnode
           if(.not.touch(jposi).and.ifoun) then            ! Position not touched           
              do klist=1,nlist                             ! Search other elements 
                 kelem=liste(klist)                        ! where JPOIN is and 
                 nnodk=lnnod(kelem)
                 do knode=1,nnodk                          ! touch their position
                    kpoin=lnods(knode,kelem)
                    if(kpoin==jpoin) then
                       kposi=(klist-1)*mnode+knode
                       touch(kposi)=.true.
                    end if
                 end do
              end do
              nzdom = nzdom+1
           end if
        end if
     end do
  end do

end subroutine ibm_nzecof

subroutine ibm_arrind(&
     nlist,ncoef,liste,touch,nzdom,ipoin,npoin,r_dom_2,c_dom_2)
  !-----------------------------------------------------------------------
  !                 
  ! This routine constructs the arrays of indexes for a mesh graph.
  ! These are organized as follows (CSR format):
  !   
  ! R_DOM(IPOIN) = coefficient of the graph where IPOIN starts,
  ! C_DOM(NZDOM) = column of the NZDOM coefficient of the graph.
  !              
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  mnode,lnods,lnnod,r_dom,c_dom
  use def_master, only       :  npoi1
  implicit none
  integer(ip), intent(in)    :: nlist,ncoef,ipoin,npoin
  integer(ip), intent(in)    :: liste(nlist)
  integer(ip), intent(inout) :: nzdom
  integer(ip), intent(out)   :: r_dom_2(*)
  integer(ip), intent(inout) :: c_dom_2(*)
  logical(lg), intent(inout) :: touch(ncoef)
  integer(ip)                :: jelem,jnode,jpoin,nnodj,jposi,jlist
  integer(ip)                :: kelem,knode,kpoin,nnodk,kposi,klist,izdom
  logical(lg)                :: ifoun

  r_dom_2(ipoin-npoi1) = nzdom

  do jlist = 1,nlist
     jelem = liste(jlist)
     nnodj = lnnod(jelem)
     do jnode = 1,nnodj
        jpoin = lnods(jnode,jelem)
        if( jpoin /= ipoin ) then        
           !
           ! In parallel, some nodes are not connected with their own neighbor nodes inside a subdomain.
           ! That is, these neighbor nodes are not stored in the data structure "r_dom"
           ! For this reason, the data structure "r_dom_2" has to store these neighbor nodes.  
           !
           ! Example:
           ! 1:         First  subdomian
           ! 2:         Second subdomian
           ! O:         A node
           ! A,B,C,D,E: Neighbor nodes of the node "O"
           !
           !     --------C-------D
           !    / \     / \     / \
           !   / 1 \ 1 / 2 \ 1 / 1 \ 
           !  /     \ /     \ /     \
           ! --------B-------O-------E
           !  \     / \     / \     /
           !   \ 1 / 1 \ 2 / 2 \ 1 /
           !    \ /     \ /     \ /
           !     --------A-------F
           !
           ! In this example, the node "O" is not connected with the nodes "A" and "B" inside the subdomain "1".
           ! This means that the nodes "A" and "B" were not added to the data structure "r_dom".
           ! The solution is add the nodes "A" and "B" to the data structure "r_dom_2".
           !           
           ifoun = .true.
           izdom = r_dom(ipoin)
           do while( izdom <= r_dom(ipoin+1)-1 )
              if( c_dom(izdom) == jpoin ) then
                 ifoun = .false.
                 izdom = r_dom(ipoin+1)
              end if
              izdom = izdom + 1              
           end do

           jposi = (jlist-1)*mnode+jnode
           if(.not.touch(jposi).and.ifoun) then
              do klist = 1,nlist
                 kelem = liste(klist)
                 nnodk = lnnod(kelem)
                 do knode = 1,nnodk
                    kpoin = lnods(knode,kelem)
                    if(kpoin==jpoin) then
                       kposi = (klist-1)*mnode+knode
                       touch(kposi) = .true.
                    end if
                 end do
              end do
              c_dom_2(nzdom) = jpoin
              nzdom = nzdom+1
           end if
        end if
     end do
  end do

end subroutine ibm_arrind
