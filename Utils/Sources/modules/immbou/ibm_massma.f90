subroutine ibm_massma()
  !------------------------------------------------------------------------
  !****f* Nastin/ibm_massma
  ! NAME 
  !    ibm_massma
  ! DESCRIPTION
  !    Compute the following arrays:
  !    - MASSC(2:NDIME+1,NPOIN) ... Compute restriction array
  !    - MASSC(1,NPOIN) ........... Hole boundary mass matrix (diagonal)
  !
  !    Where the restriction is for mass conservation
  !             +-
  !    R^a_i =  |  Na * ni ds
  !            -+ p
  ! USES
  ! USED BY
  !    immbou
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use mod_gradie
  use mod_kdtree
  use mod_memchk

  implicit none
  integer(ip)             :: ielem,inode,ipoin,inodb,iboun
  integer(ip)             :: idime
  integer(ip)             :: pnodb,pblty
  integer(4)              :: istat
  integer(ip)             :: iimbo,ielpo,ielty,iface,jelem
  integer(ip)             :: nnodb,pboun,salir,count
  integer(ip)             :: minpo,minp2,ipara,eltes(nelem_2)
  integer(ip), pointer    :: lnofb(:,:),ltyfb(:),lelfb(:)

  real(rp)                :: bocod(ndime,mnodb),baloc(9),gpdet
  real(rp)                :: elmas(npoin)
  real(rp)                :: norfa(ndime)


  if( IMASTER ) return
  !-------------------------------------------------------------------
  !
  ! Initialize mass matrix
  !
  !-------------------------------------------------------------------

  do ipoin = 1,npoin
     do idime = 1,ndime+1
        massc(idime,ipoin) = 0.0_rp
     end do
  end do

  !-------------------------------------------------------------------
  !
  ! Determine the boundary of the particle 
  !
  !-------------------------------------------------------------------

  do iimbo = 1,nimbo

     ! Temporal array for element boundary build with fringe nodes in a particle
     allocate(lnofb(mnodb,nelem),stat=istat)
     call memchk(zero,istat,memor_dom,'LNODE','bounib',lnofb)

     ! Temporal array for the types of element boundary build with fringe nodes in a particle
     allocate(ltyfb(nelem),stat=istat)
     call memchk(zero,istat,memor_dom,'LNODE','bounib',ltyfb)

     ! Temporal array for the types of element boundary build with fringe nodes in a particle
     allocate(lelfb(nelem),stat=istat)
     call memchk(zero,istat,memor_dom,'LNODE','bounib',lelfb)

     iboun = 0
     ! Loop over elements
     do ielem = 1,nelem
        ! If ielem element is inside the particle
        if ( letib(ielem) == iimbo ) Then           
           ielty=abs(ltype(ielem))
           ! Loop over faces
           do iface=1,nface(ielty)
              nnodb = nnode ( ltypf(ielty)%l(iface) )
              inodb = 0
              salir = 0
              ipara = 1
              iboun = iboun + 1             
              ! Nodes in each faces
              do while (salir == 0 .and. inodb < nnodb)
                 inodb = inodb + 1
                 ! Local face node (inside the element)
                 inode = lface(ielty)%l(inodb,iface)
                 ! Global face node (into domain)                
                 ipoin  = lnods(inode,ielem)
                 ! Test if the node is fringe
                 if ( lntib(ipoin) < 0 ) then                    
                    lnofb(inodb,iboun) = ipoin
                 else
                    salir = 1
                 end if
                 !
                 ! Test if the face is shared with another subdomain
                 !
                 if (ipoin <= npoi1)  ipara = 0
              end do
              !
              ! Test if the face is really a boundary of the particle
              !
              if (salir == 0) then
                 salir = 1      
                 count = 0
                 !
                 ! 1. Find the pair of elements that share the current face
                 !
                 do inodb = 1,nnodb
                    do ielpo = pelpo(lnofb(inodb,iboun)),pelpo(lnofb(inodb,iboun)+1)-1
                       jelem = lelpo(ielpo)
                       eltes(jelem) = 0_ip                 
                    end do
                    if (ipara == 1) then                          
                       do ielpo = pelpo_2(lnofb(inodb,iboun)-npoi1),pelpo_2(lnofb(inodb,iboun)-npoi1+1)-1
                          jelem = lelpo_2(ielpo)
                          eltes(jelem) = 0_ip                       
                       end do
                    end if
                 end do
                 do inodb = 1,nnodb
                    do ielpo = pelpo(lnofb(inodb,iboun)),pelpo(lnofb(inodb,iboun)+1)-1
                       jelem = lelpo(ielpo)
                       eltes(jelem) = eltes(jelem) + 1
                    end do
                    if (ipara == 1) then                          
                       do ielpo = pelpo_2(lnofb(inodb,iboun)-npoi1),pelpo_2(lnofb(inodb,iboun)-npoi1+1)-1
                          jelem = lelpo_2(ielpo)
                          eltes(jelem) = eltes(jelem) + 1
                       end do
                    end if
                 end do
                 !
                 ! 2. Check if one of the found elements is not a hole
                 ! 
                 inodb = 0 
                 do while (inodb < nnodb .and. salir == 1)
                    inodb = inodb + 1
                    ielpo = pelpo(lnofb(inodb,iboun)) - 1
                    do while (ielpo < pelpo(lnofb(inodb,iboun)+1)-1 .and. salir == 1)
                       ielpo = ielpo + 1                      
                       jelem = lelpo(ielpo)
                       if (eltes(jelem) == nnodb .and. jelem /= ielem) count = count + 1
                       if (eltes(jelem) == nnodb .and. jelem /= ielem .and. letib(jelem) == 0)  salir = 0
                    end do
                    if (ipara == 1) then                                                    
                       ielpo = pelpo_2(lnofb(inodb,iboun)-npoi1)-1
                       do while (ielpo < pelpo_2(lnofb(inodb,iboun)-npoi1+1)-1 .and. salir == 1)
                          ielpo = ielpo + 1                      
                          jelem = lelpo_2(ielpo)
                          if (eltes(jelem) == nnodb .and. jelem /= ielem) count = count + 1
                          if (eltes(jelem) == nnodb .and. jelem /= ielem .and. letib(jelem) == 0  ) salir = 0
                       end do
                    end if
                 end do
                 if (count == 0) salir = 0
              end if
              !
              ! Test if a new boundary element was built
              !
              if (salir == 1) then
                 iboun = iboun - 1
              else                                 
                 ltyfb(iboun) = ltypf(ielty)%l(iface)
                 lelfb(iboun) = ielem
              end if

           end do

        end if
     end do
     pboun=iboun

     ! Particle arrays for element boundary build with fringe nodes
     allocate(imbou(iimbo)%lnofb(mnodb,pboun),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo)%lnofb','BOUNIB',imbou(iimbo)%lnofb)

     allocate(imbou(iimbo)%ltyfb(pboun),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo)%lnofb','BOUNIB',imbou(iimbo)%ltyfb)

     allocate(imbou(iimbo)%lelfb(pboun),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo)%lnofb','BOUNIB',imbou(iimbo)%lelfb)


     imbou(iimbo)%npofb = mnodb 
     imbou(iimbo)%nbofb = pboun
     ! Put the information in the particle
     do iboun=1,pboun
        imbou(iimbo)%ltyfb(iboun) = ltyfb(iboun)
        imbou(iimbo)%lelfb(iboun) = lelfb(iboun)
        nnodb = nnode( ltyfb(iboun) )
        do inode=1,nnodb
           imbou(iimbo)%lnofb(inode,iboun) = lnofb(inode,iboun)
        end do
     end do
  end do
  !-------------------------------------------------------------------
  !
  ! Compute the matrices for mass conservation of the particles
  !
  !-------------------------------------------------------------------
  do iimbo = 1,nimbo
     boundaries: do iboun=1,imbou(iimbo)%nbofb
        !
        ! Boundary properties and dimensions        
        !
        pblty    = imbou(iimbo)%ltyfb(iboun)
        pnodb    = nnode(pblty)
        do inode = 1,pnodb
           ipoin = imbou(iimbo)%lnofb(inode,iboun)
           bocod(1,inode) = coord(1,ipoin)
           bocod(2,inode) = coord(2,ipoin)
           bocod(ndime,inode) = coord(ndime,ipoin)
        end do
        !
        ! Compute lumped mass matrix for boundary built with fringe node
        !
        do inode = 1,pnodb 
           call bouder(&
                pnodb,ndime,ndimb,elmar(pblty)%deric(1,1,inode),&               ! Cartesian derivative
                bocod,baloc,gpdet)                                              ! and Jacobian
           elmas(inode) = elmar(pblty)%weigc(inode)*gpdet           
        end do
        do inode = 1,pnodb
           ipoin = imbou(iimbo)%lnofb(inode,iboun)
           massc(1,ipoin) = massc(1,ipoin) + elmas(inode)
        end do
        !
        ! norfa: Exterior normal for boundary built with fringe node           
        !
        call extbou(1,pnodb,imbou(iimbo)%lnofb(:,iboun),coord,norfa)
        !
        ! Compute restriction matrix for boundary built with fringe node
        !
        do inode = 1,pnodb
           ipoin = imbou(iimbo)%lnofb(inode,iboun)      
           do idime = 1,ndime
              massc(idime+1,ipoin) = massc(idime+1,ipoin) + elmas(inode) * norfa(idime)
           end do
        end do
     end do boundaries
  end do


  if( INOTMASTER ) then
     call pararr('SLX',NPOIN_TYPE,(ndime+1)*npoin,massc)
  end if

end subroutine ibm_massma
