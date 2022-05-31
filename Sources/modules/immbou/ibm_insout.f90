subroutine ibm_insout()
  !-----------------------------------------------------------------------
  !****f* ibm_insout/ibm_insout
  ! NAME
  !    ibm_insout
  ! DESCRIPTION
  !    This routines determine the nodes of the fluid mesh that are inside 
  !    the rigid bodies
  ! USED BY
  !    nsi_bounib
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_elmtyp
  use def_domain
  use def_immbou
  use mod_memchk
  use mod_kdtree
  use mod_communications, only : PAR_GHOST_NODE_EXCHANGE
  !use mod_cutele
  use mod_messages, only : livinf

  implicit none
  integer(ip)           :: iimbo,idime,inode,jnode,posix,posiy,posiz,kboun
  integer(ip)           :: ipoin,ielem,ielty,izdom,jpoin,knode,lnode,pnode
  integer(ip)           :: nnoib,nfaib,temp,subox(3,2),nboxe(3)
  integer(ip)           :: ntrav,ielpo,kface,updat,ihair,limit
  integer(4)            :: istat
  integer(ip), pointer  :: nodib(:),lnodl(:),lnodg(:)
  real(rp)              :: delta(3),proje(ndime),dumma(ndime)
  real(rp)              :: dista,dummr,disno,disti,distj

  if( ittim /= 0 ) then
     if( kfl_colli_ibm == 0 ) then  
        call livinf(165_ip,'IO,',0_ip)
     else
        call livinf( 56_ip,' ',modul) 
        call livinf(160_ip,' ',1_ip)  
        call livinf(165_ip,'IO,',modul)        
     end if
  end if
  if( kfl_coibm == 0 ) return
  if( IMASTER )        return

  !----------------------------------------------------------------------
  !
  ! Initialize
  !
  !----------------------------------------------------------------------

  do ipoin = 1,npoin_2
     lnti2(ipoin) = lntib(ipoin)
     lntib(ipoin) = 0_ip
  end do
  do ielem = 1,nelem_2
     letib(ielem) = 0_ip     
  end do
  if (kfl_diric_ibm == 1) then
     do ipoin = 1,npoin
        lntra(ipoin)     =  0_ip
     end do
  end if


  nboxe(3) = 1
  do idime = 1,ndime
     nboxe(idime) = nubox_ibm
     delta(idime) = (xmima_ibm(2,idime)-xmima_ibm(1,idime)) / real(nboxe(idime),rp)
     delta(idime) = delta(idime) + delta(idime) / 1000.0_rp
  end do

  subox(3,2) = 1_ip
  subox(3,1) = 1_ip

  !----------------------------------------------------------------------
  !
  ! Determine in which particles the nodes are in IIMBO
  ! LNTIB(IPOIN) = IIMBO
  !
  !----------------------------------------------------------------------

  do iimbo = 1,nimbo   

     !if( imbou(iimbo) % kfl_coupl == 0 ) then
     !   
     ! Number of faces
     ! 
     nfaib = imbou(iimbo) % nboib   
     imbou(iimbo) % kfl_insid = 1
     !
     ! Determine the boxes that contains the fluid nodes that maybe are inside the body
     !
     do idime = 1,ndime

        if( imbou(iimbo)%bobox(idime,1) > xmima_ibm(2,idime) ) imbou(iimbo) % kfl_insid = 0
        if( imbou(iimbo)%bobox(idime,2) < xmima_ibm(1,idime) ) imbou(iimbo) % kfl_insid = 0

        subox(idime,1) = int( ( imbou(iimbo) % bobox(idime,1) - xmima_ibm(1,idime) ) / delta(idime) , ip ) + 1_ip
        subox(idime,2) = int( ( imbou(iimbo) % bobox(idime,2) - xmima_ibm(1,idime) ) / delta(idime) , ip ) + 1_ip

     end do

     subox(1,1) = max(1_ip,subox(1,1))
     subox(1,2) = min(size(boxes_ibm,1,kind=ip),subox(1,2))
     subox(2,1) = max(1_ip,subox(2,1))
     subox(2,2) = min(size(boxes_ibm,2,kind=ip),subox(2,2))
     subox(3,1) = max(1_ip,subox(3,1))
     subox(3,2) = min(size(boxes_ibm,3,kind=ip),subox(3,2))

     if( imbou(iimbo) % kfl_insid == 1 ) then
        !
        ! Loop over the boxes
        !    
        do posix = subox(1,1),subox(1,2)
           do posiy = subox(2,1),subox(2,2)
              do posiz = subox(3,1),subox(3,2)

                 nnoib =  boxes_ibm(posix,posiy,posiz) % nnode
                 nodib => boxes_ibm(posix,posiy,posiz) % nodes
                 !
                 ! Loop over the nodes inside the boxes
                 !
                 do inode = 1,nnoib

                    ipoin =  nodib(inode)

                    temp  =  0_ip
                    do idime = 1,ndime
                       if ( coord(idime,ipoin) >= imbou(iimbo) % bobox(idime,1) .and. &
                            coord(idime,ipoin) <= imbou(iimbo) % bobox(idime,2) ) then
                          temp = temp + 1_ip
                       end if
                    end do
                    ! 
                    ! If the current node is inside the body bounding box 
                    !
                    if ( temp == ndime ) then   
                       call dpopar(1_ip,coord(:,ipoin),&
                            imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                            imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                            dista,dumma,proje,kboun,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                            imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                            imbou(iimbo) % lnele)
                       ! 
                       ! If the current node is inside the body
                       !                        
                       if( dista <= 0.0_rp ) then
                          !
                          ! LNTIB: Inside nodes
                          !
                          if( imbou(iimbo) % kfl_coupl == 0 ) then
                             lntib( ipoin ) =  iimbo               ! Dirichlet
                          else
                             lntib( ipoin ) = -iimbo               ! Force
                          end if
                       end if
                    end if

                 end do
              end do
           end do
        end do
     end if
  end do
  !
  ! Identify the nodes that are outside the particle but have almost one neighbor node inside the particle.
  !     
  do ipoin = 1,npoin
     if( lntib(ipoin) > 0 ) then   
        iimbo = lntib(ipoin)                            
        if( imbou(iimbo) % kfl_typeb == 0 .and. imbou(iimbo) % kfl_coupl == 0 ) then
           call dpopar(1_ip,coord(:,ipoin),&
                imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                disti,dumma,proje,kboun,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                imbou(iimbo) % lnele)
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
              jpoin = c_dom(izdom)
              if ( lntib(jpoin) == 0 ) then                 
                 call dpopar(1_ip,coord(:,jpoin),&
                      imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                      imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                      distj,dumma,proje,kboun,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                      imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                      imbou(iimbo) % lnele)
                 if (abs(distj) < abs(disti)) lntib(jpoin) = -iimbo
                 !lntib(jpoin) = -iimbo
              end if
           end do
        end if
     end if
     !
     ! In parallel, a node could have a neighbor node inside the particle in other subdomains
     !     
     if (npoin_2 > npoin .and. ipoin > npoi1) then           
        do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
           jpoin = c_dom_2(izdom)
           if (lntib(ipoin) == 0 .and. lntib(jpoin) > 0) then             
              iimbo = lntib(jpoin)              
              if( imbou(iimbo) % kfl_typeb == 0 .and. imbou(iimbo) % kfl_coupl == 0 ) then
                 call dpopar(1_ip,coord(:,ipoin),&
                      imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                      imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                      disti,dumma,proje,kboun,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                      imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                      imbou(iimbo) % lnele)
                 
                 call dpopar(1_ip,coord(:,jpoin),&
                      imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                      imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                      distj,dumma,proje,kboun,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                      imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                      imbou(iimbo) % lnele)
                 if (abs(disti) < abs(distj)) lntib(ipoin) = -iimbo
                 !lntib(ipoin) = -iimbo
              end if
           end if
        end do
     end if
  end do
  do ipoin = 1,npoin_2
     if( lntib(ipoin) < 0 ) then
        lntib(ipoin) = abs(lntib(ipoin))
     end if
  end do
  !
  ! Interchange information in the fringe geometry.
  !
  !call parari('FRI',NPOIN_TYPE,npoin_2,lntib)
  call PAR_GHOST_NODE_EXCHANGE(lntib,'SUBSTITUTE','IN MY CODE')

  !----------------------------------------------------------------------
  !
  ! LNTIB: Identify fringe nodes
  !
  !----------------------------------------------------------------------

  do ipoin = 1,npoin
     if( lntib(ipoin) > 0 ) then
        iimbo = abs(lntib(ipoin))
        !
        ! Check neighbors
        !
        izdom = r_dom(ipoin)
        do while( izdom <= r_dom(ipoin+1)-1 )
           jpoin = c_dom(izdom)
           if( lntib(jpoin) == 0 ) then
              lntib(ipoin) = -lntib(ipoin)
              izdom = r_dom(ipoin+1)
           end if
           izdom = izdom + 1
        end do
        !
        ! Check neighbors in the fringe geometry too
        !
        if (npoin_2 > npoin .and. ipoin > npoi1) then        
           izdom = r_dom_2(ipoin-npoi1)
           do while( izdom <= r_dom_2(ipoin-npoi1+1)-1 )
              jpoin = c_dom_2(izdom)
              if( lntib(jpoin) == 0 ) then
                 lntib(ipoin) = -iimbo
                 izdom = r_dom_2(ipoin-npoi1+1)
              end if
              izdom = izdom + 1
           end do
        end if
     end if
  end do
  !
  ! Interchange information in the fringe geometry.
  !
  !call parari('FRI',NPOIN_TYPE,npoin_2,lntib)
  call PAR_GHOST_NODE_EXCHANGE(lntib,'SUBSTITUTE','IN MY CODE')

  !----------------------------------------------------------------------
  !
  ! LETIB: Identify hole elements (use to remove deformed fringe nodes)
  ! LTYPE: Put negativea value for holes
  ! LELCH: Define holes
  !
  !----------------------------------------------------------------------
  !
  ! Recover original value
  !
  do ielem = 1,nelem_2
     lelch(ielem) = lelch_ibm(ielem)
     ltype(ielem) = abs(ltype(ielem)) ! Don't delete this line. 
     if( lelch(ielem) == ELHOL ) ltype(ielem) = -abs(ltype(ielem))
  end do
  !
  ! Create holes (or not if the method is body force)
  !
  do ielem = 1,nelem_2
     pnode        = nnode(abs(ltype(ielem)))
     knode        = 1
     lnode        = 0
     letib(ielem) = 0
     iimbo        = 0
     do while( knode <= pnode )
        ipoin = lnods(knode,ielem)
        !
        ! Check if the lntib of ipoin is defined (beware of fringe geometry)  
        !
        iimbo = abs(lntib(ipoin))
        if ( lntib(ipoin) /= 0 ) then
           lnode = lnode + 1
           if( imbou(iimbo) % kfl_coupl > 0 ) then
              lnode = pnode
              knode = pnode
           end if
        end if
        knode = knode + 1
     end do
     if( lnode == pnode ) then
        letib(ielem) = iimbo
        if( imbou(iimbo) % kfl_typeb == 0 .and. imbou(iimbo) % kfl_coupl == 0 ) then
           ltype(ielem) = -abs(ltype(ielem))
           lelch(ielem) =  ELHOL
        end if
     end if
  end do

  !----------------------------------------------------------------------
  !
  ! LNTIB: Remove hairs
  !
  !----------------------------------------------------------------------
  do ipoin = 1,npoin
     if(lntib(ipoin) < 0) then
        ihair = 1_ip
        do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
           ielem = lelpo(ielpo)           
           if (letib(ielem) > 0) ihair = 0_ip
        end do
        if (npoin_2 > npoin .and. ipoin > npoi1) then                
           do ielpo = pelpo_2(ipoin-npoi1),pelpo_2(ipoin-npoi1+1)-1
              ielem = lelpo_2(ielpo)
              if (letib(ielem) > 0) ihair = 0_ip
           end do
        end if
        if (ihair == 1_ip) lntib(ipoin) = 0
     end if
  end do

  !----------------------------------------------------------------------
  !
  ! LETIB: Identify fringe elements
  !
  !----------------------------------------------------------------------

!!$  do ielem = 1,nelem
!!$     jnode = 0_ip
!!$     knode = 0_ip
!!$     if( ltype(ielem) > 0 ) then
!!$        ielty = ltype(ielem)
!!$        pnode = nnode(ielty)
!!$        do inode = 1,pnode
!!$           ipoin = lnods(inode,ielem)
!!$           if ( lntib(ipoin) < 0 ) then
!!$              jnode = 1_ip
!!$              iimbo = abs(lntib(ipoin))
!!$           elseif ( lntib(ipoin) == 0 ) then
!!$              knode = 1_ip
!!$           end if
!!$        end do
!!$        !
!!$        ! It is a fringe element
!!$        !
!!$        if ( jnode == 1 .and. knode == 1 ) then
!!$
!!$           !----------------------------------------------------------------------
!!$           !
!!$           ! CUT FRINGE ELEMENT 
!!$           !
!!$           !----------------------------------------------------------------------
!!$           if( imbou(iimbo) % kfl_typeb == 0 .and. imbou(iimbo) % kfl_coupl > 0 ) then
!!$              letib(ielem) =  iimbo
!!$              !if( lelch_ibm(ielem) /= -1000 ) then
!!$              !   lelch(ielem)     = ELCUT
!!$              !   lelch_ibm(ielem) = ELCUT
!!$              call runend('MOD_CUTELE is under construction')
!!$              !end if
!!$           end if
!!$        end if
!!$     end if
!!$  end do




!!$  if( kfl_diric_ibm == 1 ) then
!!$     !----------------------------------------------------------------------
!!$     !
!!$     ! LNTRA: Identify travesties nodes and nearest node if i am using a 
!!$     !        dirichlet condition 
!!$     !
!!$     !----------------------------------------------------------------------
!!$     ! Identify travesties: 
!!$     !     - I am fringe or free: LNTIB <= 0 
!!$     !     - I was hole:          LNTRA >  0
!!$     ! Put LNTRA = -1 for travesties
!!$     !     LNTRA =  0 otherwise
!!$     !
!!$     ntrav = 0
!!$     do ipoin = 1,npoin
!!$        if ( lntib(ipoin) <= 0 .and. lnti2(ipoin) > 0 ) then 
!!$           lntra(ipoin) = -1_ip
!!$           ntrav = ntrav + 1
!!$        end if
!!$     end do
!!$     !
!!$     ! LNTRA(IPOIN) = JPOIN. JPOIN is nearest node for travesties only
!!$     !
!!$     do ipoin = 1,npoin
!!$        if( lntra(ipoin) < 0_ip ) then
!!$           !
!!$           ! Allocate arrays
!!$           !
!!$           limit = 0
!!$           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
!!$              jpoin = c_dom(izdom)
!!$              if ( lnti2(jpoin) <= 0 ) then
!!$                 limit = limit + 1
!!$              end if
!!$           end do
!!$           if (npoin_2 > npoin .and. ipoin > npoi1) then
!!$              do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
!!$                 jpoin = c_dom_2(izdom)
!!$                 if ( lnti2(jpoin) <= 0 ) then
!!$                    limit = limit + 1
!!$                 end if
!!$              end do
!!$           end if                    
!!$           allocate( lnodl(limit) )
!!$           allocate( lnodg(limit) )
!!$           !
!!$           ! Find the neiborgs nodes which status were fringe or free before 
!!$           !                   
!!$           inode = 0
!!$           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
!!$              jpoin = c_dom(izdom)
!!$              if ( lnti2(jpoin) <= 0 ) then
!!$                 inode = inode + 1
!!$                 lnodl(inode) = jpoin
!!$                 lnodg(inode) = lninv_loc(jpoin)
!!$              end if
!!$           end do
!!$           !
!!$           ! Find the neiborgs nodes which were fringe or free before in the fringe geometry 
!!$           !                    
!!$           if (npoin_2 > npoin .and. ipoin > npoi1) then
!!$              do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
!!$                 jpoin = c_dom_2(izdom)           
!!$                 if ( lnti2(jpoin) <= 0 ) then
!!$                    inode = inode + 1
!!$                    lnodl(inode) = jpoin
!!$                    lnodg(inode) = lninv_loc(jpoin)
!!$                 end if
!!$              end do
!!$           end if
!!$           !
!!$           ! Sort          
!!$           !
!!$           call  heapsorti2(2_ip,limit,lnodg,lnodl)
!!$           ! 
!!$           ! Determine the closer neiborg node
!!$           !
!!$           dista = 1.0e10_rp
!!$           do inode = 1,limit
!!$              jpoin = lnodl(inode)
!!$              dummr = 0.0_rp
!!$              do idime = 1,ndime
!!$                 dummr = dummr + (coord(idime,ipoin) - coord(idime,jpoin))**2.0_rp
!!$              end do
!!$
!!$              if (dummr < dista) then
!!$                 dista        = dummr
!!$                 lntra(ipoin) = jpoin
!!$              end if
!!$           end do
!!$           !
!!$           ! Deallocate arrays
!!$           !
!!$           deallocate(lnodl)
!!$           deallocate(lnodg)
!!$        end if
!!$     end do
!!$  end if

  !----------------------------------------------------------------------
  !
  ! LNOCH: Update list of node characteristic
  !
  !----------------------------------------------------------------------
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        lnoch(ipoin) = lnoch_ibm(ipoin)
        if( lntib(ipoin) > 0 ) then
           iimbo = lntib(ipoin)
           if( imbou(iimbo) % kfl_coupl == 0 ) then
              lnoch(ipoin) = NOHOL
           end if
        else if( lnoch(ipoin) /= NOFEM .and. lnoch(ipoin) /= NOHOL ) then
           call runend('IBM_INSOUT: NOT CODED')
        end if
     end do
  end if

  do ipoin = 1,npoin
     if (lntib(ipoin) < 0_ip .and. lntib(ipoin) /= lnti2(ipoin) .and. ittim >= 1) then
        print *,"--| ALYA     THERE IS A NEW NODE: ",lninv_loc(ipoin)
     end if
  end do

end subroutine ibm_insout
