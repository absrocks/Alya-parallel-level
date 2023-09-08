module mod_vol

  use mod_sort
  use mod_cart

contains

  subroutine mshvol(ndim,npoin,lface,nnofa,nface,rsize,nelem,nnode,coor,&
       elem,nnosi)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    use mod_genbou
    use mod_optim, only : optim3d
    implicit none
    integer(ip),intent(in)    :: ndim,nface,nnofa,nnode,nnosi
    integer(ip),intent(inout) :: lface(nnofa,nface)
    integer(ip),intent(inout) :: npoin,nelem
    real(rp),pointer          :: rsize(:)
    real(rp),pointer          :: coor(:,:)
    integer(ip),pointer       :: elem(:,:)
    integer(ip),pointer       :: lfath(:),ltet(:),lmark(:),lelem(:),lfathf(:)
    integer(ip),pointer       :: eltoel(:,:),ptoel1(:),ptoel2(:),lboup(:),lhole(:)
    integer(ip),pointer       :: lhashf1(:,:),lhashf2(:)
    integer(ip),pointer       :: lcart(:),lpsur(:),lptype(:,:),lpsid(:) 
    type(cell),pointer        :: lcell(:) 
    real(rp)                  :: bboxp(ndim,2),bbox(ndim,2),rgrowth,stopcri,hmax
    real(rp),pointer          :: coorf(:,:),rqual(:),rnopo(:,:),rnofa(:,:)
    integer(ip)               :: npoinf,ipoin,nboup,nhashf,ncell,npnew,npoinf0
    integer(ip)               :: iface,nhole 
    integer(4)                :: istat
    !
    !     This is the main sub for volume mesh generation
    !
    !
    !     On input:
    !
    !             - lface: the surface mesh
    !             - elem:  the pointer to volume mesh
    !             - rsize: the size of the surface points 
    !
    !
    !
    !
    !     Data structure:
    ! 
    !
    !          Point database: 
    ! 
    !               - ltet: one element containing ipoin   
    !
    !
    allocate(ltet(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LTET','mshvol',ltet)
    !
    !     First copy coor to coorf
    !
    allocate(coorf(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'COORF','mshvol',coorf)
    allocate(lhole(max(1_ip,npoin/100_ip)),stat=istat)
    call memchk(zero,istat,memor_msh,'LHOLE','mshvol',lhole)
    do ipoin=1,npoin
       coorf(1,ipoin)=coor(1,ipoin) 
       coorf(2,ipoin)=coor(2,ipoin) 
       coorf(3,ipoin)=coor(3,ipoin) 
    enddo
    !
    !     Remember npoinf
    !
    npoinf=npoin
    npoinf0=npoinf
    !
    !     Initialize the allowed growth
    !
    rgrowth=1.0d+01
    !
    !     Initialize the expected ratio
    !
    stopcri=5.0d+00
    !
    !     Compute bbox of the boundary points
    !
    call boxbin(coorf,ndim,npoinf,bbox)
    !
    !     Remember it
    !
    bboxp(1,1)=bbox(1,1)
    bboxp(2,1)=bbox(2,1)
    bboxp(3,1)=bbox(3,1)
    bboxp(1,2)=bbox(1,2)
    bboxp(2,2)=bbox(2,2)
    bboxp(3,2)=bbox(3,2)
    !
    !     Compute bbox size to be conforming with the size distribution
    !
    call boxsiz(npoinf,ndim,rsize,bbox,rgrowth,stopcri,hmax)
    !
    !     Compute initial elements
    !
    call inivol(coor,elem,eltoel,ndim,npoin,bbox,ltet,nelem,nnode)
    !
    !     Build outer mesh point
    !
    call inibou(lface,nface,nnofa,elem,nelem,nnode,eltoel,coorf,npoinf,&
         ndim,rgrowth,rsize,coor,npoin,bbox,lfathf,hmax,bboxp)
    !
    !     DBG
    !
    if(npoin/=8)then
       write(*,*)'Error in mod_vol, npoin=',npoin
       stop
    endif
    !
    !     Add the points from the surface and outer points to coor
    !
    npnew=npoin+npoinf
    !
    !     Insert the points of the outer mesh
    !
    call memrea(npnew,memor_msh,'LTET','mshvol',ltet)
    call memrea(npnew,memor_msh,'COOR','mshvol',coor)
    allocate(lfath(npnew),stat=istat)
    call memchk(zero,istat,memor_msh,'LFATH','mshvol',lfath)

    do ipoin=1,8
       lfath(ipoin)=ipoin
    enddo

    do ipoin=1,npoinf
       npoin=npoin+1_ip
       coor(1,npoin)=coorf(1,ipoin)
       coor(2,npoin)=coorf(2,ipoin)
       coor(3,npoin)=coorf(3,ipoin)
       lfath(npoin)=lfathf(ipoin)+8_ip
    enddo
    !
    !     Shift the points in lface
    !
    do iface=1,nface
       lface(1,iface)=lface(1,iface)+8_ip
       lface(2,iface)=lface(2,iface)+8_ip
       lface(3,iface)=lface(3,iface)+8_ip
    enddo
    !call outface(nnofa,nface,npoin,ndim,lface,coor)
    !
    !     Insert the outer points
    !
    call insert(coor,npoin,ndim,elem,nelem,eltoel,nnode,lfath,ltet,npoinf0)
    !
    !     Regenerate the boundary
    !
    call genbou(lface,nface,nnofa,nelem,elem,nnode,coor,ndim,npoin,lfath,nnosi,&
         ltet,eltoel,lelem,rqual)
    !
    !     Generate the mesh
    ! 
    call advf3d(lface,nnofa,nface,ndim,npoin,coor,rsize,rnopo,lcart,lpsur,lmark,&
         lptype,lpsid,lcell,ncell)
    !
    !     Remove outer elements
    ! 
    call rmv3d(nnode,nelem,elem,lhashf2,lhashf1,nhashf,eltoel)
    !
    !     Get boudary points
    !
    nboup=npoinf
    allocate(lboup(nboup),stat=istat)
    call memchk(zero,istat,memor_msh,'LBOUP','mshvol',lboup)
    do ipoin=1,npoinf
       lboup(ipoin)=ipoin
    enddo
    !
    !     Optimize the mesh
    !
    call optim3d(nelem,ndim,npoin,nnode,elem,coor,eltoel,nnosi,&
         rsize,nnofa,lboup,lhole,nhole)

    call memchk(2_ip,istat,memor_msh,'COORF','mshvol',coorf)
    deallocate(coorf,stat=istat)
    if(istat/=0) call memerr(2_ip,'COORF','mshvol',0_ip)
    call memchk(2_ip,istat,memor_msh,'LTET','mshvol',ltet)
    deallocate(ltet,stat=istat)
    if(istat/=0) call memerr(2_ip,'LTET','mshvol',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFATH','mshvol',lfath)
    deallocate(lfath,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFATH','mshvol',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFATHF','mshvol',lfathf)
    deallocate(lfathf,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFATHF','mshvol',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHOLE','mshvol',lhole)
    deallocate(lhole,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHOLE','mshvol',0_ip)

  end subroutine mshvol

  subroutine boxsiz(npoin,ndim,rsize,bbox,rgrowth,stopcri,hmax)
    use def_kintyp, only :  ip,rp,lg
    implicit none
    integer(ip),intent(in)       :: npoin,ndim
    real(rp),intent(in)          :: rsize(npoin),rgrowth,stopcri
    real(rp),intent(inout)       :: bbox(ndim,2),hmax
    real(rp)              :: c12,rsmin,dx,dy,dz,hx,hy,hz,distx,disty,distz
    real(rp)              :: ratiox,ratioy,ratioz,c(ndim),dl,growthc
    integer(ip)           :: ipoin
    !
    !     This sub computes the bbox of the outer mesh to avoid a bad conditioning
    !
    c12=1.0d+00/2.0d+00
    !
    !     First get the smallest side for estimation
    !    
    rsmin=rsize(1)
    do ipoin=2,npoin
       if(rsize(ipoin)<rsmin)rsmin=rsize(ipoin)
    enddo
    !
    !     Bounding box dimensions
    !
    dx=bbox(1,2)-bbox(1,1) 
    dy=bbox(2,2)-bbox(2,1) 
    dz=bbox(3,2)-bbox(3,1)
    !
    !     Initialize iterative process in x
    !  
    hx=rsmin
    !
    !     Length 
    !
    distx=dx*c12+hx
    ratiox=distx/hx

    do

       if(ratiox<stopcri)exit

       hx=hx*rgrowth
       distx=distx+hx
       ratiox=distx/hx

    enddo
    !
    !     Initialize iterative process in y
    !  
    hy=rsmin
    !
    !     Length 
    !
    disty=dy*c12+hy
    ratioy=disty/hy

    do

       if(ratioy<stopcri)exit

       hy=hy*rgrowth
       disty=disty+hy
       ratioy=disty/hy

    enddo
    !
    !     Initialize iterative process in z
    !  
    hz=rsmin
    !
    !     Length 
    !
    distz=dz*c12+hz
    ratioz=distz/hz

    do

       if(ratioz<stopcri)exit

       hz=hz*rgrowth
       distz=distz+hz
       ratioz=distz/hz

    enddo
    !
    !     From this, define the size of the final Bbox
    !
    if(distx>disty)then
       dl=distx
       hmax=hx
    else 
       dl=disty
       hmax=hy
    endif
    if(dl<distz)then
       dl=distz 
       hmax=hz
    endif

    c(1)=(bbox(1,1)+bbox(1,2))*c12
    c(2)=(bbox(2,1)+bbox(2,2))*c12
    c(3)=(bbox(3,1)+bbox(3,2))*c12

    bbox(1,1)=c(1)-hmax
    bbox(2,1)=c(2)-hmax
    bbox(3,1)=c(3)-hmax
    bbox(1,2)=c(1)+hmax
    bbox(2,2)=c(2)+hmax
    bbox(3,2)=c(3)+hmax

  end subroutine boxsiz

  subroutine inivol(coor,elem,eltoel,ndim,npoin,bbox,ltet,nelem,nnode)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)       :: ndim,nnode
    integer(ip),intent(inout)    :: npoin,nelem
    real(rp),intent(in)          :: bbox(ndim,2)
    integer(ip),pointer          :: elem(:,:),eltoel(:,:),ltet(:)
    real(rp),pointer             :: coor(:,:)
    integer(4)                ::  istat
    !
    !     This sub generates the first 5 elements of the new mesh
    !
    npoin=8_ip
    call memrea(npoin,memor_msh,'COOR','inivol',coor)
    call memrea(npoin,memor_msh,'LTET','inivol',ltet)

    coor(1,1)=bbox(1,1) 
    coor(2,1)=bbox(2,1) 
    coor(3,1)=bbox(3,1) 

    coor(1,2)=bbox(1,2) 
    coor(2,2)=bbox(2,1) 
    coor(3,2)=bbox(3,1) 

    coor(1,3)=bbox(1,2) 
    coor(2,3)=bbox(2,2) 
    coor(3,3)=bbox(3,1) 

    coor(1,4)=bbox(1,1) 
    coor(2,4)=bbox(2,2) 
    coor(3,4)=bbox(3,1) 

    coor(1,5)=bbox(1,1) 
    coor(2,5)=bbox(2,1) 
    coor(3,5)=bbox(3,2) 

    coor(1,6)=bbox(1,2) 
    coor(2,6)=bbox(2,1) 
    coor(3,6)=bbox(3,2) 

    coor(1,7)=bbox(1,2) 
    coor(2,7)=bbox(2,2) 
    coor(3,7)=bbox(3,2) 

    coor(1,8)=bbox(1,1) 
    coor(2,8)=bbox(2,2) 
    coor(3,8)=bbox(3,2)

    nelem=5_ip
    allocate(elem(nnode,nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'ELEM','inivol',elem)

    elem(1,1)=1_ip
    elem(2,1)=4_ip
    elem(3,1)=5_ip
    elem(4,1)=2_ip

    elem(1,2)=8_ip
    elem(2,2)=5_ip
    elem(3,2)=4_ip
    elem(4,2)=7_ip

    elem(1,3)=6_ip
    elem(2,3)=7_ip
    elem(3,3)=2_ip
    elem(4,3)=5_ip

    elem(1,4)=2_ip
    elem(2,4)=7_ip
    elem(3,4)=3_ip
    elem(4,4)=4_ip

    elem(1,5)=7_ip
    elem(2,5)=5_ip
    elem(3,5)=4_ip
    elem(4,5)=2_ip

    ltet(1)=1_ip
    ltet(4)=1_ip
    ltet(5)=1_ip
    ltet(2)=1_ip
    ltet(8)=2_ip
    ltet(7)=2_ip
    ltet(6)=3_ip
    ltet(3)=4_ip

    allocate(eltoel(nnode,nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'ELTOEL','inivol',eltoel)

    eltoel(1,1)=5_ip
    eltoel(2,1)=0_ip
    eltoel(3,1)=0_ip
    eltoel(4,1)=0_ip

    eltoel(1,2)=5_ip
    eltoel(2,2)=0_ip
    eltoel(3,2)=0_ip
    eltoel(4,2)=0_ip

    eltoel(1,3)=5_ip
    eltoel(2,3)=0_ip
    eltoel(3,3)=0_ip
    eltoel(4,3)=0_ip

    eltoel(1,4)=0_ip
    eltoel(2,4)=0_ip
    eltoel(3,4)=5_ip
    eltoel(4,4)=0_ip

    eltoel(1,5)=1_ip
    eltoel(2,5)=4_ip
    eltoel(3,5)=3_ip
    eltoel(4,5)=2_ip

    !c(1)=(bbox(1,1)+bbox(1,2))*c12
    !c(2)=(bbox(2,1)+bbox(2,2))*c12
    !c(3)=(bbox(3,1)+bbox(3,2))*c12

    !dx=bbox(1,2)-c(1)
    !dy=bbox(2,2)-c(2)
    !dz=bbox(3,2)-c(3)

    !rini=sqrt(dx*dx+dy*dy+dz*dz)

  end subroutine inivol

  subroutine inibou(lface,nface,nnofa,elem,nelem,nnode,eltoel,coorf,npoinf,&
       ndim,rgrowth,rsize,coor,npoin,bbox,lfathf,hmax,bboxp)
    use def_kintyp, only :  ip,rp,lg,cell
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)           :: nface,nnofa,ndim,nnode 
    integer(ip),intent(inout)        :: npoin,nelem,npoinf 
    integer(ip),intent(inout)        :: lface(nnofa,nface) 
    real(rp),intent(in)              :: rgrowth,rsize(npoinf)
    real(rp),intent(in)              :: bbox(ndim,2),bboxp(ndim,2),hmax 
    type(cell),pointer               :: lcellp(:)
    integer(ip), pointer             :: elem(:,:),eltoel(:,:),ptoel1(:),ptoel2(:) 
    real(rp), pointer                :: coor(:,:),rnopo(:,:),rnofa(:,:),coorf(:,:)
    integer(ip),pointer              :: lheap(:),lcart(:),level(:),lfathf(:)
    integer(ip),pointer              :: lheapt(:)
    integer(ip)                      :: nstackt,nstack,ipoin,iguess,ipadv,ihost,ncellp  
    integer(ip)                      :: ichk,ipclos,istack,npmax,ilevel,ifclos  
    integer(ip)                      :: nlevel,maxlevel,ilev,ipnew,iplace,mcellp  
    integer(ip)                      :: ipstack,nstack0,noutwards  
    integer(ip),pointer              :: lstack(:),lmarkc(:),lctop(:,:),lstackt(:),lstackc(:)
    real(rp),pointer                 :: rsizef(:),rsizeft(:) 
    real(rp)                         :: rsloc,nrmal1(ndim),nrmal2(ndim),c00,epsil,c10
    real(rp)                         :: rl,pnew(ndim,5),rtol
    integer(4)                       :: istat
    !
    !     This sub generates the points around the initial surface mesh to avoid
    !     a bad conditioning
    !
    c00=0.0d+00
    c10=1.0d+00
    epsil=1.0d-08
    rtol=1.0d-10
    !
    !     Initialize maxlevel for sorting later on
    !
    maxlevel=1_ip
    !
    !     Nullify
    !
    nullify(ptoel1,ptoel2)
    !
    !     Allocate temporary arrays
    ! 
    allocate(rnofa(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOFA','inibou',rnofa)
    allocate(rnopo(ndim,npoinf),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOPO','inibou',rnopo)
    allocate(lstack(npoinf),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','inibou',lstack)
    ! 
    !     Allocate lcart  && lfathf
    ! 
    allocate(lcart(npoinf),stat=istat)
    call memchk(zero,istat,memor_msh,'LCART','inibou',lcart)
    allocate(lfathf(npoinf),stat=istat)
    call memchk(zero,istat,memor_msh,'LFATHF','inibou',lfathf)
    allocate(level(npoinf),stat=istat)
    call memchk(zero,istat,memor_msh,'LEVEL','inibou',level)
    !
    !     Initialize lfathf
    !
    do ipoin=1,npoinf
       lfathf(ipoin)=ipoin
    enddo
    !
    !     Initialize level
    !
    do ipoin=1,npoinf
       level(ipoin)=1_ip
    enddo
    !
    !     Set the max number of points for each cell
    !
    npmax=8_ip
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoinf,nnofa,ptoel1,ptoel2)
    !
    !     Get the face normals
    !
    call gtfnrl(lface,nface,nnofa,ndim,coorf,npoinf,rnofa)
    !
    !     Get the point normals
    !
    call gtpnrl(nface,rnopo,npoinf,ndim,ptoel1,ptoel2,rnofa)
    !
    !     Create the cartesian mesh point with all the surface points
    !
    call cartp(npoinf,ndim,lface,nnofa,nface,npmax,lcart,ptoel1,ptoel2,coorf,lcellp,ncellp,rtol,1_ip,bbox) 
    !
    !     Filter the initial close points on the surface 
    !
    call filterp(nnofa,nface,lface,ndim,npoinf,coorf,lstack,nstack,&
         rgrowth,rsize,ptoel1,ptoel2,bboxp)
    !
    !     Allocate cartesian related arrays
    !
    allocate(lmarkc(ncellp),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKC','inibou',lmarkc)
    allocate(lstackc(ncellp),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACKC','inibou',lstackc)
    allocate(lctop(npmax+1,ncellp),stat=istat)
    call memchk(zero,istat,memor_msh,'LCTOP','inibou',lctop)
    mcellp=ncellp
    !
    !     Interpolate the points of the stack in the cartesian mesh
    !
    call intercart(coorf,npoinf,ndim,lface,nface,nnofa,lcellp,ncellp,lcart,ptoel1,ptoel2,rtol)
    !
    !     Fill the linked list
    !
    call ctopnt2(lcart,npoinf,ncellp,lctop,npmax+1_ip )
    !
    !     Allocate rsizef for the filtered points
    ! 
    allocate(rsizef(nstack),stat=istat)
    call memchk(zero,istat,memor_msh,'RSIZEF','inibou',rsizef)
    allocate(rsizeft(nstack),stat=istat)
    call memchk(zero,istat,memor_msh,'RSIZEFT','inibou',rsizeft)

    do istack=1,nstack
       rsizef(istack)=rsize(lstack(istack))
    enddo
    !
    !     Allocate the heap for the filtered points
    ! 
    allocate(lheap(nstack),stat=istat)
    call memchk(zero,istat,memor_msh,'LHEAP','inibou',lheap)
    allocate(lheapt(nstack),stat=istat)
    call memchk(zero,istat,memor_msh,'LHEAPT','inibou',lheapt)
    !
    !     Initialize the heap
    !
    call buildheap(lheap,nstack,rsizef)
    allocate(lstackt(nstack),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACKT','inibou',lstackt)
    !
    !     Copy to temporary
    !
    do istack=1,nstack 
       rsizeft(istack)=rsizef(istack)
       lstackt(istack)=lstack(istack)
       lheapt(istack)=lheap(istack)
    enddo
    nstackt=nstack
    nstack0=nstack
    !
    !     First go outwards
    !  
    call advpoin(nface,nnofa,ndim,nnode,npoinf,lface,&
         rgrowth,rsizeft,coorf,bbox,lcellp,&
         rnopo,lheapt,lcart,level,lmarkc,lctop,ncellp,&
         mcellp,npmax,lfathf,rtol,nstackt,maxlevel,lstackt,hmax,&
         ptoel1,ptoel2,lstackc)
    !
    !     Invert the normals 
    ! 
    do ipoin=1,npoinf
       rnopo(1,ipoin)=-rnopo(1,ipoin)
       rnopo(2,ipoin)=-rnopo(2,ipoin)
       rnopo(3,ipoin)=-rnopo(3,ipoin)
    enddo
    !
    !     Then go inwards
    !
    call advpoin(nface,nnofa,ndim,nnode,npoinf,lface,&
         rgrowth,rsizef,coorf,bbox,lcellp,&
         rnopo,lheap,lcart,level,lmarkc,lctop,ncellp,&
         mcellp,npmax,lfathf,rtol,nstack,maxlevel,lstack,hmax,&
         ptoel1,ptoel2,lstackc)
    !
    !     Renumber the points depending of their level
    !
    call renulev(coorf,ndim,npoinf,maxlevel,level,lfathf,lface,nface,nnofa)

    call memchk(2_ip,istat,memor_msh,'RSIZEFT','inibou',rsizeft)
    deallocate(rsizeft,stat=istat)
    if(istat/=0) call memerr(2_ip,'RSIZEFT','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'RSIZEF','inibou',rsizef)
    deallocate(rsizef,stat=istat)
    if(istat/=0) call memerr(2_ip,'RSIZEF','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHEAPT','inibou',lheapt)
    deallocate(lheapt,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHEAPT','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHEAP','inibou',lheap)
    deallocate(lheap,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHEAP','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCART','inibou',lcart)
    deallocate(lcart,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCART','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACKT','inibou',lstackt)
    deallocate(lstackt,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACKT','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','inibou',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOPO','inibou',rnopo)
    deallocate(rnopo,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOPO','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOFA','inibou',rnofa)
    deallocate(rnofa,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOFA','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','inibou',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','inibou',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARKC','inibou',lmarkc)
    deallocate(lmarkc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARKC','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCTOP','inibou',lctop)
    deallocate(lctop,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCTOP','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACKC','inibou',lstackc)
    deallocate(lstackc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACKC','inibou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEVEL','inibou',level)
    deallocate(level,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEVEL','inibou',0_ip)

  end subroutine inibou

  subroutine filterp(nnofa,nface,lface,ndim,npoin,coor,lstack,nstack,&
       rgrowth,rsize,ptoel1,ptoel2,bbox)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: nface,nnofa,ndim,npoin 
    integer(ip),intent(inout)   :: nstack,lstack(npoin)
    integer(ip),intent(in)      :: lface(nnofa,nface) 
    integer(ip),intent(in)      :: ptoel2(npoin+1),ptoel1(*) 
    real(rp),intent(in)         :: coor(ndim,npoin),rgrowth,rsize(npoin) 
    real(rp),intent(in)         :: bbox(ndim,2)
    integer(ip),pointer         :: lmark(:)
    integer(4)                  :: istat
    !
    !     This sub filters the close surface points to generate a good cloud
    !     of points
    !
    !
    !     Allocate lmark
    ! 
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','filsur',lmark)
    !
    !     Filter the points on the surface
    !
    call filsur(nnofa,nface,lface,ndim,npoin,coor,lstack,nstack,rgrowth,&
         rsize,ptoel1,ptoel2,lmark)
    !
    !     Filter the remaining points extensively
    !
    call filvol(bbox,nstack,lstack,coor,ndim,npoin,rsize,lmark,rgrowth)
    !
    !     Rebuild ptoel1 and ptoel2 for the surviving points
    !
    !call rebuil(npoin,lmark,ptoel1,ptoel2,lstack,nstack,lface,nnofa,nface)

  end subroutine filterp

  subroutine filsur(nnofa,nface,lface,ndim,npoin,coor,lstack,nstack,rgrowth,&
       rsize,ptoel1,ptoel2,lmark)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: nface,nnofa,ndim,npoin 
    integer(ip),intent(inout)   :: nstack,lstack(npoin),lmark(npoin)
    integer(ip),intent(in)      :: lface(nnofa,nface),ptoel2(npoin+1),ptoel1(*) 
    real(rp),intent(in)         :: coor(ndim,npoin),rgrowth,rsize(npoin) 
    integer(ip)                 :: ipoin,jpoin,istackp,nstackp,nheap,ismall 
    integer(ip)                 :: ie,iface,inofa,kpoin 
    integer(ip),pointer         :: lstackp(:),lheap(:) 
    real(rp)                    :: rlen,dx,dy,dz,rdist
    integer(4)                  :: istat
    !
    !     This sub filters approximately the points on the surface by taking advantage
    !     of the surface connection. It is not exact as the surface may go out and come 
    !     back in the local bbox     
    ! 
    !
    !     Allocate the heap
    ! 
    allocate(lheap(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LHEAP','filsur',lheap)
    !
    !     Allocate lstackp
    ! 
    allocate(lstackp(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACKP','filsur',lstackp)
    !
    !     Initialize the heap
    !
    nheap=npoin
    call buildheap(lheap,nheap,rsize)
    !
    !     Initialize nstack
    !
    nstack=0_ip
    !
    !     Get the point with the smallest size
    !  
    do
       !
       !     Do we still have points?
       !
       if(nheap==0)exit
       !
       !     Get the point with the smallest size
       !
       call gtkey(lheap,nheap,rsize,ipoin,npoin)
       !
       !     Has this point been marked before?
       !
       if(lmark(ipoin)/=0)cycle 
       lmark(ipoin)=ipoin
       !
       !     Keep this point by putting it in the stack
       !
       nstack=nstack+1_ip
       lstack(nstack)=ipoin
       !
       !     Compute the local size
       !
       rlen=rgrowth*rsize(ipoin)
       !
       !     Initialize the local stack of points
       !
       nstackp=1_ip  
       lstackp(1)=ipoin
       istackp=0_ip
       !
       !     Loop on neighbors
       !
       do

          if(istackp==nstackp)exit
          istackp=istackp+1_ip

          jpoin=lstackp(istackp)
          !
          !     Loop on elements surrounding point
          !
          do ie=ptoel2(jpoin),ptoel2(jpoin+1)-1 
             iface=ptoel1(ie)
             do inofa=1,nnofa
                kpoin=lface(inofa,iface)
                !
                !     Has this point been marked?
                !
                if(lmark(kpoin)==0)then
                   !
                   !     Compute distance
                   !
                   dx=coor(1,ipoin)-coor(1,jpoin)                  
                   dy=coor(2,ipoin)-coor(2,jpoin)                  
                   dz=coor(3,ipoin)-coor(3,jpoin)                  
                   rdist=sqrt(dx*dx+dy*dy+dz*dz)
                   !
                   !     Is the point too close?
                   !
                   if(rdist<rlen)then                   
                      !
                      !     Add the point to the stack of points
                      !
                      nstackp=nstackp+1
                      lstackp(nstackp)=kpoin
                      !
                      !     Mark kpoin with ipoin
                      !
                      lmark(kpoin)=ipoin
                   endif
                endif
             enddo
          enddo
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LSTACKP','filsur',lstackp)
    deallocate(lstackp,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACKP','filsur',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHEAP','filsur',lheap)
    deallocate(lheap,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHEAP','filsur',0_ip)

  end subroutine filsur

  subroutine filvol(bboxp,nstack,lstack,coor,ndim,npoin,rsize,lmark,rgrowth)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: ndim,npoin 
    integer(ip),intent(inout)   :: nstack
    integer(ip),intent(inout)   :: lstack(nstack),lmark(npoin)
    real(rp),intent(in)         :: bboxp(ndim,2),coor(ndim,npoin),rsize(npoin) 
    real(rp),intent(in)         :: rgrowth 
    integer(ip)                 :: nbin,nbinx,istack,imin,jmin,kmin,jpoin
    integer(ip)                 :: ibin,ipoin,nbinxy,i,ntot,iplace,l,nstack0
    integer(ip)                 :: imax,jmax,kmax,k,kbinx,j,jbinx,isto0,isto1
    integer(ip),pointer         :: btopo1(:),btopo2(:)
    real(rp)                    :: dxbin,dybin,dzbin,rlen
    real(rp)                    :: vmin(ndim),vmax(ndim),c13
    real(rp)                    :: dx,dy,dz,rdist 
    integer(4)                  :: istat
    !
    !    This sub filters exactly the points that are too close by a bin sort
    !
    c13=1.0d+00/3.0d+00
    !
    !     Allocate the bin to point pointer
    !
    nbin=nstack
    nbinx=floor(nbin**c13)+50_ip
    nbin=nbinx*nbinx*nbinx
    nbinxy=nbinx*nbinx
    allocate(btopo2(nbin+1),stat=istat)
    call memchk(zero,istat,memor_msh,'BTOPO2','cart',btopo2)

    dxbin=(bboxp(1,2)-bboxp(1,1))/nbinx
    dybin=(bboxp(2,2)-bboxp(2,1))/nbinx
    dzbin=(bboxp(3,2)-bboxp(3,1))/nbinx

    do istack=1,nstack 
       ipoin=lstack(istack)

       imin=floor((coor(1,ipoin)-bboxp(1,1))/dxbin)+1
       jmin=floor((coor(2,ipoin)-bboxp(2,1))/dybin)+1
       kmin=floor((coor(3,ipoin)-bboxp(3,1))/dzbin)+1

       ibin=imin+(jmin-1)*nbinx+(kmin-1)*nbinxy+1_ip
       btopo2(ibin)=btopo2(ibin)+1_ip
    enddo
    !
    !     Reshuffle
    !
    btopo2(1)=1_ip
    do i=2,nbin+1_ip
       btopo2(i)=btopo2(i)+btopo2(i-1_ip)
    enddo
    !
    !     Allocate
    !
    ntot=btopo2(nbin+1)-1_ip
    allocate(btopo1(ntot),stat=istat)
    call memchk(zero,istat,memor_msh,'BTOEL1','cart',btopo1)
    ! 
    !     Store
    ! 
    do istack=1,nstack 
       ipoin=lstack(istack)

       imin=floor((coor(1,ipoin)-bboxp(1,1))/dxbin)+1
       jmin=floor((coor(2,ipoin)-bboxp(2,1))/dybin)+1
       kmin=floor((coor(3,ipoin)-bboxp(3,1))/dzbin)+1

       ibin=imin+(jmin-1)*nbinx+(kmin-1)*nbinxy
       iplace=btopo2(ibin)
       btopo1(iplace)=ipoin
       btopo2(ibin)=iplace+1_ip
    enddo
    !
    !     Clean up
    !
    do i=nbin+1,2,-1
       btopo2(i)=btopo2(i-1)
    enddo

    btopo2(1)=1
    !
    !     Loop on the points of the stack. They have already been sorted
    !     out with respect to their size
    !
    do istack=1,nstack 
       ipoin=lstack(istack)
       !
       !     Has ipoin been deleted before?
       !
       if(lmark(ipoin)/=ipoin)cycle

       imin=floor((coor(1,ipoin)-bboxp(1,1))/dxbin)+1
       jmin=floor((coor(2,ipoin)-bboxp(2,1))/dybin)+1
       kmin=floor((coor(3,ipoin)-bboxp(3,1))/dzbin)+1
       !
       !     Create the local bbox
       !
       rlen=rgrowth*rsize(ipoin)
       vmin(1)=coor(1,ipoin)-rlen
       vmin(2)=coor(2,ipoin)-rlen
       vmin(3)=coor(3,ipoin)-rlen
       vmax(1)=coor(1,ipoin)+rlen
       vmax(2)=coor(2,ipoin)+rlen
       vmax(3)=coor(3,ipoin)+rlen

       imin=floor((vmin(1)-bboxp(1,1))/dxbin)+1_ip
       jmin=floor((vmin(2)-bboxp(2,1))/dybin)+1_ip
       kmin=floor((vmin(3)-bboxp(3,1))/dzbin)+1_ip

       imax=floor((vmax(1)-bboxp(1,1))/dxbin)+1_ip
       jmax=floor((vmax(2)-bboxp(2,1))/dybin)+1_ip
       kmax=floor((vmax(3)-bboxp(3,1))/dzbin)+1_ip
       !
       !     Limit the extend
       !
       if(imax<1)then
          cycle
       else if(imin<1)then
          imin=1
       endif

       if(jmax<1)then
          cycle
       else if(jmin<1)then
          jmin=1
       endif

       if(kmax<1)then
          cycle
       else if(kmin<1)then
          kmin=1
       endif

       if(imin>nbinx)then
          cycle
       else if(imax>nbinx)then
          imax=nbinx
       endif
       if(jmin>nbinx)then
          cycle
       else if(jmax>nbinx)then
          jmax=nbinx
       endif
       if(kmin>nbinx)then
          cycle
       else if(kmax>nbinx)then
          kmax=nbinx
       endif

       do k=kmin,kmax
          kbinx=(k-1)*nbinxy
          do j=jmin,jmax
             jbinx=(j-1)*nbinx+kbinx
             do i=imin,imax
                ibin=i+jbinx

                isto0=btopo2(ibin)
                isto1=btopo2(ibin+1)-1

                do l=isto0,isto1

                   jpoin=btopo1(l)
                   !
                   !     Do not take ipoin into account
                   !
                   if(jpoin==ipoin)cycle
                   !
                   !     Compute distance
                   !
                   dx=coor(1,ipoin)-coor(1,jpoin)
                   dy=coor(2,ipoin)-coor(2,jpoin)
                   dz=coor(3,ipoin)-coor(3,jpoin)
                   rdist=sqrt(dx*dx+dy*dy+dz*dz)
                   !
                   !     Is the point too close?
                   !  
                   if(rdist<rlen)then
                      lmark(jpoin)=ipoin
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !     Compact the points
    !
    nstack0=nstack
    nstack=0_ip
    do istack=1,nstack0
       ipoin=lstack(istack)
       if(lmark(ipoin)==ipoin)then
          nstack=nstack+1
          lstack(nstack)=ipoin
       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'BTOPO1','filvol',btopo1)
    deallocate(btopo1,stat=istat)
    if(istat/=0) call memerr(2_ip,'BTOPO1','filvol',0_ip)
    call memchk(2_ip,istat,memor_msh,'BTOPO2','filvol',btopo2)
    deallocate(btopo2,stat=istat)
    if(istat/=0) call memerr(2_ip,'BTOPO2','filvol',0_ip)

  end subroutine filvol

  subroutine rebuil(npoin,lmark,ptoel1,ptoel2,lstack,nstack,lface,nnofa,nface)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip), intent(in)    :: npoin,nstack,nnofa,nface
    integer(ip),intent(in)     :: lmark(npoin),lstack(nstack),lface(nnofa,nface)
    integer(ip), intent(inout) :: ptoel1(*),ptoel2(npoin+1)   
    integer(ip),pointer        :: ptoel1t(:),ptoel2t(:),lmarkp(:),lstackp(:),lstack2(:)
    integer(ip)                :: istack,ipoin,istackp,nstackp,ief,inofa,jpoin,ntot   
    integer(ip)                :: iface,i,iplace,kpoin,istack2,nstack2   
    integer(4)                 :: istat
    logical(lg)                :: ifound
    !
    !     This sub rebuild ptoel1 and ptoel2 for the fast interpolation
    !     for the cartesian mesh. As a lot of points have been deleted,
    !     there are no faces to connect them.
    !       
    !
    !     Allocate ptoel2t
    !
    allocate(ptoel2t(nstack+1),stat=istat)
    call memchk(zero,istat,memor_msh,'PTOEL2T','rebuil',ptoel2t)
    !
    !     Allocate lmarkp
    !
    allocate(lmarkp(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKP','rebuil',lmarkp)
    !
    !     Allocate lstackp
    !
    allocate(lstackp(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACKP','rebuil',lstackp)
    allocate(lstack2(nstack),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK2','rebuil',lstack2)
    !
    !     Loop on the keeped points
    !
    do istack=1,nstack 
       ipoin=lstack(istack)
       !
       !     Initialize local stack
       !
       istackp=0_ip
       nstackp=1_ip 
       lstackp(1)=ipoin
       lmarkp(ipoin)=ipoin
       nstack2=0_ip 

       do 

          if(istackp==nstackp)exit
          istackp=istackp+1_ip 
          jpoin=lstackp(nstackp) 

          do ief=ptoel2(jpoin),ptoel2(jpoin+1)-1
             iface=ptoel1(ief)
             do inofa=1,nnofa
                kpoin=lface(inofa,iface) 
                if(lmarkp(kpoin)==ipoin)cycle
                lmarkp(kpoin)=ipoin
                !
                !     Does the point belong to the group ipoin? 
                !
                if(lmark(kpoin)==ipoin)then
                   nstackp=nstackp+1
                   lstackp(nstackp)=kpoin
                else
                   !
                   !     kpoin has not been marked with ipoin
                   !     Do we have already taken lmark(kpoin) into account?
                   !
                   ifound=.false. 
                   do istack2=1,nstack2   
                      if(lstack2(istack2)==lmark(kpoin))then
                         ifound=.true.
                         exit
                      endif
                   enddo

                   if(ifound.eqv. .false.)then
                      ptoel2t(istack+1)=ptoel2t(istack+1)+1_ip
                      nstack2=nstack2+1_ip
                      lstack2(nstack2)=lmark(kpoin)
                   endif
                endif
             enddo

          enddo

       enddo

    enddo
    !
    !     Reshuffle
    !
    ptoel2t(1)=1_ip
    do i=2,nstack+1_ip
       ptoel2t(i)=ptoel2t(i)+ptoel2t(i-1_ip)
    enddo
    !
    !     Allocate
    !
    ntot=ptoel2t(nstack+1)-1_ip
    allocate(ptoel1t(ntot),stat=istat)
    call memchk(zero,istat,memor_msh,'PTOEL1T','cart',ptoel1t)
    !
    !     Reset lmarkp
    !
    do ipoin=1,npoin
       lmarkp(ipoin)=0_ip
    enddo
    ! 
    !     Store
    ! 
    do istack=1,nstack 
       ipoin=lstack(istack)
       !
       !     Initialize local stack
       !
       istackp=0_ip
       nstackp=1_ip 
       lstackp(1)=ipoin
       lmarkp(ipoin)=ipoin
       nstack2=0_ip 

       do 

          if(istackp==nstackp)exit
          istackp=istackp+1_ip 
          jpoin=lstackp(nstackp) 

          do ief=ptoel2(jpoin),ptoel2(jpoin+1)-1
             iface=ptoel1(ief)
             do inofa=1,nnofa
                kpoin=lface(inofa,iface) 
                if(lmarkp(kpoin)==ipoin)cycle
                lmarkp(kpoin)=ipoin
                !
                !     Does the point belong to the group ipoin? 
                !
                if(lmark(kpoin)==ipoin)then
                   nstackp=nstackp+1
                   lstackp(nstackp)=kpoin
                else
                   !
                   !     kpoin has not been marked with ipoin
                   !     Do we have already taken lmark(kpoin) into account?
                   !
                   ifound=.false. 
                   do istack2=1,nstack2   
                      if(lstack2(istack2)==lmark(kpoin))then
                         ifound=.true.
                         exit
                      endif
                   enddo

                   if(ifound.eqv. .false.)then
                      iplace=ptoel2t(istack)
                      ptoel1t(iplace)=lmark(kpoin)
                      ptoel2t(istack)=iplace+1_ip
                      nstack2=nstack2+1_ip
                      lstack2(nstack2)=lmark(kpoin)
                   endif
                endif
             enddo

          enddo

       enddo

    enddo
    !
    !     Clean up
    !
    do istack=nstack+1,2,-1
       ptoel2t(istack)=ptoel2t(istack-1)
    enddo
    ptoel2t(1)=1
    !
    !     Copy ptoelt to ptoel (could be only smaller)
    !
    do istack=1,nstack+1
       ptoel2(istack)=ptoel2t(istack)
    enddo

    do i=1,ntot
       ptoel1(i)=ptoel1t(i)
    enddo
    call memchk(2_ip,istat,memor_msh,'LMARKP','rebuil',lmarkp)
    deallocate(lmarkp,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARKP','rebuil',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK2','rebuil',lstack2)
    deallocate(lstack2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK2','rebuil',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACKP','rebuil',lstackp)
    deallocate(lstackp,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACKP','rebuil',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1T','rebuil',ptoel1t)
    deallocate(ptoel1t,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1T','rebuil',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2T','rebuil',ptoel2t)
    deallocate(ptoel2t,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2T','rebuil',0_ip)

  end subroutine rebuil

  subroutine advpoin(nface,nnofa,ndim,nnode,npoin,lface,rgrowth,rsizef,coor,&
       bbox,lcellp,rnopo,lheap,lcart,level,lmarkc,lctop,ncellp,&
       mcellp,npmax,lfathf,rtol,nstack,maxlevel,lstack,hmax,&
       ptoel1,ptoel2,lstackc)
    use def_kintyp, only :  ip,rp,lg,cell
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)           :: nface,nnofa,ndim,nnode,ncellp,npmax 
    integer(ip),intent(inout)        :: mcellp,nstack,maxlevel,npoin 
    integer(ip),intent(in)           :: lface(nnofa,nface),ptoel1(*),ptoel2(npoin+1) 
    real(rp),intent(in)              :: rgrowth 
    real(rp),intent(in)              :: bbox(ndim,2),rtol 
    type(cell)                       :: lcellp(ncellp)
    real(rp), pointer                :: coor(:,:),rnopo(:,:)
    integer(ip),pointer              :: lheap(:),lcart(:),level(:),lstack(:),lfathf(:)
    integer(ip),intent(inout)        :: lmarkc(ncellp),lstackc(ncellp) 
    integer(ip),pointer              :: lctop(:,:) 
    integer(ip)                      :: nheap,ipoin,iguess,ipadv,ihost,npadv
    integer(ip)                      :: ichk,ipclos,istack,ilevel,ifclos  
    integer(ip)                      :: nlevel,ilev,iplace,iguessloc  
    real(rp),pointer                 :: rsizef(:),coort(:,:) 
    real(rp)                         :: rsloc,hmax,nrmal1(ndim),nrmal2(ndim),c00,epsil,c10
    real(rp)                         :: rl,pnew(ndim,5),c12,rsloc2
    integer(4)                       :: istat
    integer(ip),parameter            :: mstackp=500
    integer(ip)                      :: lstackp(mstackp),nstackp,npmax1,ifound

    c00=0.0d+00
    c10=1.0d+00
    c12=1.0d+00/2.0d+00
    epsil=1.0d-08 
    !
    !     Initialize the stack
    !
    nheap=nstack
    npmax1=npmax+1_ip
    !
    !     Loop on heap 
    !  
    do
       !
       !     Stop condition
       !
       if(nheap==0)exit
       !
       !     Get the point with the smallest size
       !  
       call gtkey(lheap,nheap,rsizef,istack,nstack)
       !
       !     Output mesh
       !
       call outGidpnt(ndim,npoin,coor,bbox)
       !
       !     Get the local size
       !
       rsloc=rsizef(istack)*rgrowth
       rsloc2=rsloc*c12
       ipoin=lstack(istack)
       ilevel=level(ipoin) 
       !
       !     Limit the size 
       !
       if(rsloc>hmax)rsloc=hmax
       !
       !     Get the tangent plane
       !
       nrmal2(1)=c00
       nrmal2(2)=c00
       nrmal2(3)=c00

       if(abs(rnopo(1,ipoin))<epsil)then
          nrmal2(1)=c10
       else if(abs(rnopo(2,ipoin))<epsil)then
          nrmal2(2)=c10
       else
          nrmal2(3)=c10
       endif

       nrmal1(1)= rnopo(2,ipoin)*nrmal2(3)-rnopo(3,ipoin)*nrmal2(2)
       nrmal1(2)=-rnopo(1,ipoin)*nrmal2(3)+rnopo(3,ipoin)*nrmal2(1)
       nrmal1(3)= rnopo(1,ipoin)*nrmal2(2)-rnopo(2,ipoin)*nrmal2(1)

       rl=sqrt(nrmal1(1)*nrmal1(1)+nrmal1(2)*nrmal1(2)+nrmal1(3)*nrmal1(3))
       rl=c10/rl
       nrmal1(1)=rl*nrmal1(1)
       nrmal1(2)=rl*nrmal1(2)
       nrmal1(3)=rl*nrmal1(3)

       nrmal2(1)= rnopo(2,ipoin)*nrmal1(3)-rnopo(3,ipoin)*nrmal1(2)
       nrmal2(2)=-rnopo(1,ipoin)*nrmal1(3)+rnopo(3,ipoin)*nrmal1(1)
       nrmal2(3)= rnopo(1,ipoin)*nrmal1(2)-rnopo(2,ipoin)*nrmal1(1)

       rl=sqrt(nrmal2(1)*nrmal2(1)+nrmal2(2)*nrmal2(2)+nrmal2(3)*nrmal2(3))
       rl=c10/rl
       nrmal2(1)=rl*nrmal2(1)
       nrmal2(2)=rl*nrmal2(2)
       nrmal2(3)=rl*nrmal2(3)
       !
       !     Are we on the surface?
       !
       if(ilevel==1)then 
          !
          !     Fill the stencil
          !
          pnew(1,1)=coor(1,ipoin)+rsloc*rnopo(1,ipoin)
          pnew(2,1)=coor(2,ipoin)+rsloc*rnopo(2,ipoin)
          pnew(3,1)=coor(3,ipoin)+rsloc*rnopo(3,ipoin)
          npadv=0_ip       
          !
          !     Filter with bbox
          !
          if(pnew(1,1)>bbox(1,1) .and. pnew(1,1)<bbox(1,2))then
             if(pnew(2,1)>bbox(2,1) .and. pnew(2,1)<bbox(2,2))then
                if(pnew(3,1)>bbox(3,1) .and. pnew(3,1)<bbox(3,2))then
                   npadv=1_ip 
                endif
             endif
          endif

       else
          !
          !     Fill the stencil
          !
          pnew(1,1)=coor(1,ipoin)+rsloc*rnopo(1,ipoin)
          pnew(2,1)=coor(2,ipoin)+rsloc*rnopo(2,ipoin)
          pnew(3,1)=coor(3,ipoin)+rsloc*rnopo(3,ipoin)
          npadv=0_ip
          !
          !     Filter with bbox
          !
          if(pnew(1,1)>bbox(1,1) .and. pnew(1,1)<bbox(1,2))then
             if(pnew(2,1)>bbox(2,1) .and. pnew(2,1)<bbox(2,2))then
                if(pnew(3,1)>bbox(3,1) .and. pnew(3,1)<bbox(3,2))then
                   npadv=1_ip
                endif
             endif
          endif
          pnew(1,2)=coor(1,ipoin)+rsloc*nrmal1(1)
          pnew(2,2)=coor(2,ipoin)+rsloc*nrmal1(2)
          pnew(3,2)=coor(3,ipoin)+rsloc*nrmal1(3)
          !
          !     Filter with bbox
          !
          if(pnew(1,2)>bbox(1,1) .and. pnew(1,2)<bbox(1,2))then
             if(pnew(2,2)>bbox(2,1) .and. pnew(2,2)<bbox(2,2))then
                if(pnew(3,2)>bbox(3,1) .and. pnew(3,2)<bbox(3,2))then
                   npadv=npadv+1_ip
                   pnew(1,npadv)=pnew(1,2)
                   pnew(2,npadv)=pnew(2,2)
                   pnew(3,npadv)=pnew(3,2)
                endif
             endif
          endif

          pnew(1,3)=coor(1,ipoin)-rsloc*nrmal1(1)
          pnew(2,3)=coor(2,ipoin)-rsloc*nrmal1(2)
          pnew(3,3)=coor(3,ipoin)-rsloc*nrmal1(3)
          ! 
          !     Filter with bbox
          !
          if(pnew(1,3)>bbox(1,1) .and. pnew(1,3)<bbox(1,2))then
             if(pnew(2,3)>bbox(2,1) .and. pnew(2,3)<bbox(2,2))then
                if(pnew(3,3)>bbox(3,1) .and. pnew(3,3)<bbox(3,2))then
                   npadv=npadv+1_ip
                   pnew(1,npadv)=pnew(1,3)
                   pnew(2,npadv)=pnew(2,3)
                   pnew(3,npadv)=pnew(3,3)
                endif
             endif
          endif

          pnew(1,4)=coor(1,ipoin)+rsloc*nrmal2(1)
          pnew(2,4)=coor(2,ipoin)+rsloc*nrmal2(2)
          pnew(3,4)=coor(3,ipoin)+rsloc*nrmal2(3)
          ! 
          !     Filter with bbox
          !
          if(pnew(1,4)>bbox(1,1) .and. pnew(1,4)<bbox(1,2))then
             if(pnew(2,4)>bbox(2,1) .and. pnew(2,4)<bbox(2,2))then
                if(pnew(3,4)>bbox(3,1) .and. pnew(3,4)<bbox(3,2))then
                   npadv=npadv+1_ip
                   pnew(1,npadv)=pnew(1,4)
                   pnew(2,npadv)=pnew(2,4)
                   pnew(3,npadv)=pnew(3,4)
                endif
             endif
          endif

          pnew(1,5)=coor(1,ipoin)-rsloc*nrmal2(1)
          pnew(2,5)=coor(2,ipoin)-rsloc*nrmal2(2)
          pnew(3,5)=coor(3,ipoin)-rsloc*nrmal2(3)
          !
          !     Filter with bbox
          !
          if(pnew(1,5)>bbox(1,1) .and. pnew(1,5)<bbox(1,2))then
             if(pnew(2,5)>bbox(2,1) .and. pnew(2,5)<bbox(2,2))then
                if(pnew(3,5)>bbox(3,1) .and. pnew(3,5)<bbox(3,2))then
                   npadv=npadv+1_ip
                   pnew(1,npadv)=pnew(1,5)
                   pnew(2,npadv)=pnew(2,5)
                   pnew(3,npadv)=pnew(3,5)
                endif
             endif
          endif

       endif
       !
       !     Initialize the cartesian point cell
       !
       iguess=lcart(ipoin)
       !
       !     Loop on stencil
       !
       do ipadv=1,npadv
          !
          !     Get the cell containing the new point
          !
          iguessloc=iguess 
          call gtelem(ipadv,pnew,6_ip,ndim,lcellp,ncellp,iguessloc,rtol)
          !
          !     Check if there are some points too close
          !
          ifound=0_ip
          call getClosIn3(pnew(1,ipadv),iguessloc,ndim,coor,npoin,lcellp,ncellp,lmarkc,lstackc,&
               lctop,rsloc2,npmax1,mcellp,-1_ip,ifound)
          !
          !     Do we have a close point
          ! 
          if(ifound/=0)cycle
          !
          !     Check if there are some faces too close
          !
          !call closfa(ifclos)
          !
          !     Do we have a close face
          ! 
          !if(ifclos/=0)cycle
          !
          !     Store coordinates
          ! 
          npoin=npoin+1_ip
          call memrea(npoin,memor_msh,'COOR','inivol',coor)
          call memrea(npoin,memor_msh,'RNOPO','inivol',rnopo)
          call memrea(npoin,memor_msh,'LCART','inivol',lcart)
          call memrea(npoin,memor_msh,'LFATHF','inivol',lfathf)
          call memrea(npoin,memor_msh,'LEVEL','inivol',level)
          coor(1,npoin)=pnew(1,ipadv)
          coor(2,npoin)=pnew(2,ipadv)
          coor(3,npoin)=pnew(3,ipadv)
          !
          !     Inherit normal and point cell number
          !  
          rnopo(1,npoin)=rnopo(1,ipoin)
          rnopo(2,npoin)=rnopo(2,ipoin)
          rnopo(3,npoin)=rnopo(3,ipoin)
          !
          !     Update pointer to father
          !  
          lfathf(npoin)=ipoin  
          !
          !     Update levels
          !
          nlevel=ilevel+1_ip
          level(nstack)=nlevel 
          if(maxlevel<nlevel)maxlevel=nlevel
          !
          !      Add the point to the stack
          ! 
          nstack=nstack+1
          call memrea(nstack,memor_msh,'LSTACK','inivol',lstack)
          call memrea(nstack,memor_msh,'RSIZEF','inivol',rsizef)
          lstack(nstack)=npoin
          !
          !     Store size
          !
          rsizef(nstack)=rsloc
          !
          !     Add the point to the heap
          !
          call addkey(lheap,nheap,rsizef,nstack,nstack)   
          !
          !     Add the point in the cartesian mesh point and interpolate lcart
          !
          call addcellp(npoin,coor,ndim,npoin,lcellp,ncellp,lctop,mcellp,iguessloc,npmax1,lcart,rtol)

       enddo

    enddo

  end subroutine advpoin

  subroutine closfa(ichk)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: ichk 

  end subroutine closfa

  subroutine renulev(coor,ndim,npoin,maxlevel,level,lfath,lface,nface,nnofa)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: ndim,npoin,maxlevel,nnofa,nface
    integer(ip),intent(in)      :: level(npoin)    
    integer(ip),intent(inout)   :: lface(nnofa,nface)    
    real(rp),intent(inout)      :: coor(ndim,npoin)
    integer(ip),intent(inout)   :: lfath(npoin) 
    integer(ip),pointer         :: lplev(:),lrenu(:),lfatht(:) 
    real(rp),pointer            :: coort(:,:) 
    integer(ip)                 :: istack,ilevel,ilev,iplace,ipnew,ipoin,iface 
    integer(4)                  :: istat
    !
    !     Allocate lplev 
    !
    allocate(lplev(maxlevel+1_ip),stat=istat)
    call memchk(zero,istat,memor_msh,'LPLEV','renulev',lplev)
    !
    !     Reorder the points with respect to their levels in decreasing order
    !
    do ipoin=1,npoin
       ilevel=level(ipoin)+1_ip
       lplev(ilevel)=lplev(ilevel)+1_ip
    enddo
    !
    !     Sum up
    !
    lplev(1)=1_ip
    do ilev=2,maxlevel+1
       lplev(ilev)=lplev(ilev)+lplev(ilev-1)
    enddo
    !
    !     Allocate lrenu
    !
    allocate(lrenu(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','renulev',lrenu)
    !
    !     Allocate lfatht
    !
    allocate(lfatht(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LFATHT','renulev',lfatht)
    !
    !     Allocate coort
    !
    allocate(coort(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'COORT','renulev',coort)
    !
    !     Store
    ! 
    do ipoin=1,npoin
       ilevel=level(ipoin)
       iplace=lplev(ilevel)
       lrenu(ipoin)=iplace
       lplev(ilevel)=iplace+1_ip
    enddo
    !
    !     Reverse the ordering to have the points further from the surface
    !     as  the first points
    !
    do ipoin=1,npoin
       lrenu(ipoin)=npoin-lrenu(ipoin)+1_ip
    enddo
    !
    !     Renumber
    ! 
    do ipoin=1,npoin
       ipnew=lrenu(ipoin)
       lfatht(ipnew)=lfath(ipnew) 
       coort(1,ipnew)=coor(1,ipoin) 
       coort(2,ipnew)=coor(2,ipoin) 
       coort(3,ipnew)=coor(3,ipoin) 
    enddo

    do ipoin=1,npoin
       lfath(ipoin)=lfatht(ipoin)
       coor(1,ipoin)=coort(1,ipoin) 
       coor(2,ipoin)=coort(2,ipoin) 
       coor(3,ipoin)=coort(3,ipoin) 
    enddo
    !
    !     Renumber lface
    !
    do iface=1,nface
       lface(1,iface)=lrenu(lface(1,iface))
       lface(2,iface)=lrenu(lface(2,iface))
       lface(3,iface)=lrenu(lface(3,iface))
    enddo

    call memchk(2_ip,istat,memor_msh,'COORT','renulev',coort)
    deallocate(coort,stat=istat)
    if(istat/=0) call memerr(2_ip,'COORT','renulev',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFATHT','renulev',lfatht)
    deallocate(lfatht,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFATHT','renulev',0_ip)
    call memchk(2_ip,istat,memor_msh,'LRENU','renulev',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','renulev',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPLEV','renulev',lplev)
    deallocate(lplev,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPLEV','renulev',0_ip)

  end subroutine renulev

  subroutine insert(coor,npoin,ndim,elem,nelem,eltoel,nnode,lfath,ltet,npoinf)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_voltol, only :  gtvol,gtqual,findin,split3d,swaploc3d
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: ndim,npoin,nnode,npoinf
    integer(ip),intent(inout)   :: nelem
    integer(ip),pointer         :: elem(:,:),eltoel(:,:)
    real(rp),intent(in)         :: coor(ndim,npoin) 
    integer(ip),intent(inout)   :: lfath(npoin),ltet(npoin) 
    integer(ip)                 :: ipoin,ielem,isplit,idir,iguess,ipclos,npboun
    integer(ip),pointer         :: lelem(:)
    real(rp),pointer            :: rqual(:)
    real(rp)                    :: rvol,rqual1
    integer(4)                  :: istat
    !
    !     This sub inserts the points on the outer boundary
    ! 
    allocate(lelem(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEM','insert',lelem)
    allocate(rqual(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'RQUAL','insert',rqual)
    !
    !     Compute the quality of the initial elements
    !
    do ielem=1,nelem
       call gtvol(elem,ielem,coor,rvol,nnode,nelem,npoin,ndim)
       call gtqual(elem,ielem,coor,rvol,nnode,nelem,ndim,rqual1,npoin)
       rqual(ielem)=rqual1 
    enddo
    !
    !     Invert lfath
    !
    call invlfa(lfath,npoin) 
    !
    !     Insert the points in order
    !
    npboun=npoin-npoinf 
    do ipoin=9,npboun
       !
       !     Get the parent
       !
       ipclos=lfath(ipoin)
       !
       !     Find an initial element guess
       !
       iguess=ltet(ipclos)  
       !
       !     Find the element containing ipoin
       !
       call findin(iguess,elem,nelem,npoin,coor,coor(1,ipoin),ndim,ielem,&
            nnode,eltoel,isplit,idir) 
       !
       !     Insert the point in the element
       !
       call split3d(ipoin,coor,ndim,npoin,nelem,elem,nnode,ielem,isplit,idir,ltet,eltoel)
       !
       !     Optimize the elements
       !
       call swaploc3d(ielem,elem,nnode,nelem,coor,npoin,ndim,lelem,eltoel,rqual,ltet)

    enddo

    call memchk(2_ip,istat,memor_msh,'RQUAL','insert',rqual)
    deallocate(rqual,stat=istat)
    if(istat/=0) call memerr(2_ip,'RQUAL','insert',0_ip)
    call memchk(2_ip,istat,memor_msh,'LELEM','insert',lelem)
    deallocate(lelem,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEM','insert',0_ip)

  end subroutine insert

  subroutine invlfa(lfath,npoin)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    implicit none
    integer(ip), intent(in)      :: npoin
    integer(ip), intent(inout)   :: lfath(npoin)
    integer(ip),pointer          :: lfatht(:)
    integer(4)                   :: istat
    integer(ip)                  :: ipfath,ipoin 
    !
    !     This sub inverts the relation for the outer boundary points
    !     It is assumed that the points have been sorted by levels
    !     The first points are the points further from the surface
    !
    !
    !     Allocate lfatht
    !
    allocate(lfatht(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LFATHT','invlfa',lfatht)
    !
    !     Loop on all the points, give the child to the father
    !  
    do ipoin=1,npoin
       ipfath=lfath(ipoin) 
       !
       !     For the points on the surface, lfath(ipoin)=ipoin
       !     As these points are at the end of the array, they could
       !     delete the information given by their children 
       !
       if(ipfath/=ipoin)then 
          lfatht(ipfath)=ipoin
       endif
    enddo
    !
    !     Find the points of the outer layer. They are the last children
    !     so they have not been marked.
    !     Initialize them with the first point already in the mesh
    ! 
    do ipoin=1,npoin
       if(lfatht(ipoin)==0)then 
          lfatht(ipoin)=1_ip
       endif
    enddo
    !
    !     Copy lfatht to lfath
    !   
    do ipoin=1,npoin
       lfath(ipoin)=lfatht(ipoin)
    enddo

    call memchk(2_ip,istat,memor_msh,'LFATHT','invlfa',lfatht)
    deallocate(lfatht,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFATHT','invlfa',0_ip)

  end subroutine invlfa

  subroutine frfront(lface,nnofa,nface,rnofa,coor,npoin,ndim,ltet,elem,nelem,&
       nnode,rsize,lmark,lelem,eltoel,rqual)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    use mod_voltol, only :  findin,split3d,swaploc3d
    implicit none
    integer(ip),intent(in)    :: nnofa,nface,ndim,nnode
    integer(ip),intent(inout) :: npoin,nelem
    integer(ip),intent(in)    :: lface(nnofa,nface)
    integer(ip),pointer       :: ltet(:),elem(:,:),lmark(:)
    integer(ip),pointer       :: lelem(:),eltoel(:,:)
    real(rp),intent(in)       :: rnofa(ndim,nface)
    real(rp),pointer          :: coor(:,:),rsize(:),rqual(:)
    real(rp)                  :: pnew(ndim)
    integer(ip)               :: iface,ipa,ipb,ipc,ipclos,iguess,isplit,idir,ihost
    !
    !     Loop on the front
    !
    do iface=1,nface

       ipa=lface(1,iface)
       ipb=lface(2,iface)
       ipc=lface(3,iface)
       !
       !     Create a best point
       !
       call newpnt(lface,iface,nnofa,rnofa,ndim,coor,npoin,nface,pnew,rsize) 
       !
       !     Find an initial guess
       !
       iguess=ltet(ipa)  
       !
       !     Add the point
       !
       npoin=npoin+1_ip
       call memrea(npoin,memor_msh,'COOR','frfront',coor)
       call memrea(npoin,memor_msh,'RSIZE','frfront',rsize)
       call memrea(npoin,memor_msh,'LMARK','frfront',lmark)
       !
       !     Find the element containing npoin
       !
       call findin(iguess,elem,nelem,npoin,coor,coor(:,npoin),ndim,ihost,&
            nnode,eltoel,isplit,idir) 
       !
       !     Check if there are points too close
       !
       call gtclos3d(ihost,npoin,elem,nelem,nnode,coor,ndim,npoin,rsize,ipclos,lmark,&
            lelem,eltoel,ipa,ipb,ipc)
       if(ipclos/=0)then
          !
          !     Do we have an imposed point? 
          !  


          npoin=npoin-1
          cycle
       endif
       !
       !     Insert the point in the element
       !
       call split3d(npoin,coor,ndim,npoin,nelem,elem,nnode,ihost,isplit,idir,ltet,eltoel)
       !
       !     Optimize the elements
       !
       call swaploc3d(ihost,elem,nnode,nelem,coor,npoin,ndim,lelem,eltoel,rqual,ltet)

    enddo

  end subroutine frfront

  subroutine newpnt(lface,iface,nnofa,rnofa,ndim,coor,npoin,nface,pnew,rsize)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    implicit none
    integer(ip),intent(in) :: nnofa,nface,ndim,npoin,iface
    integer(ip),intent(in) :: lface(nnofa,nface)
    real(rp),intent(in)    :: rnofa(ndim,nface),coor(ndim,npoin),rsize(npoin)
    real(rp),intent(inout) :: pnew(ndim)
    real(rp)               :: fact,c13,rnew,pmid(ndim) 
    integer(ip)            :: ip1,ip2,ip3 

    fact=sqrt(2.0d+00/3.0d+00)/3.0d+00 
    c13=1.0d+00/3.0d+00

    ip1=lface(1,iface)
    ip2=lface(2,iface)
    ip3=lface(3,iface)
    !
    !     Get the size
    !
    rnew=(rsize(ip1)+rsize(ip2)+rsize(ip3))*fact
    !
    !     Get the mid point
    !
    pmid(1)=(coor(1,ip1)+coor(1,ip2)+coor(1,ip3))*c13
    pmid(2)=(coor(2,ip1)+coor(2,ip2)+coor(1,ip3))*c13
    pmid(3)=(coor(3,ip1)+coor(3,ip2)+coor(1,ip3))*c13
    !
    !     Generate the new point
    ! 
    pnew(1)=pmid(1)+rnew*rnofa(1,iface)
    pnew(2)=pmid(2)+rnew*rnofa(2,iface)
    pnew(3)=pmid(3)+rnew*rnofa(3,iface)

  end subroutine newpnt

  subroutine gtclos3d(ihost,ipnew,elem,nelem,nnode,coor,ndim,npoin,rsize,ipclos,lmark,&
       lelem,eltoel,ipa,ipb,ipc)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nelem,ndim,npoin,ihost,nnode,ipnew,ipa,ipb,ipc
    integer(ip),intent(in)          :: elem(nnode,nelem),eltoel(nnode,nelem)
    real(rp),intent(in)             :: coor(ndim,npoin),rsize(npoin)
    integer(ip),intent(inout)       :: ipclos,lmark(npoin),lelem(nelem) 
    integer(ip)                     :: lstack(100),nstack,istack,ielem
    integer(ip)                     :: lstackp(100),nstackp,ichk,istackp,nstackp0
    integer(ip)                     :: lstackp2(100),nstackp2
    integer(ip)                     :: imin,ip1,ip2,ip3,ip4,ineigh,ipoin 
    real(rp)                        :: rlen,rstackp(100),rmin,vmin(3),vmax(3)
    real(rp)                        :: xmin,ymin,zmin,xmax,ymax,zmax,tollen
    !
    !     This subroutine checks if some point is too close to ipnew 
    !     ipa, ipb and ipcare not taken into account as they belong to the generating edge
    !
    tollen=1.0d+00/sqrt(2.0d+00)

    istack=0_ip
    nstack=1_ip
    lstack(1)=ihost
    lelem(ihost)=1_ip
    nstackp=0_ip
    nstackp2=0_ip
    !
    !     Compute the local box
    !
    rlen=rsize(ipnew)
    vmin(1)=coor(1,ipnew)-rlen
    vmin(2)=coor(2,ipnew)-rlen
    vmin(3)=coor(3,ipnew)-rlen
    vmax(1)=coor(1,ipnew)+rlen
    vmax(2)=coor(2,ipnew)+rlen
    vmax(3)=coor(3,ipnew)+rlen
    !
    !     Find the points in the box
    !
    do                  
       if(istack==nstack)exit
       istack=istack+1
       ielem=lstack(istack) 
       !
       !     Get the box of the element
       ! 

       ip1=elem(1,ielem)
       ip2=elem(2,ielem)
       ip3=elem(3,ielem)
       ip4=elem(4,ielem)

       xmin=min(coor(1,ip1),coor(1,ip2),coor(1,ip3),coor(1,ip4))
       ymin=min(coor(2,ip1),coor(2,ip2),coor(2,ip3),coor(2,ip4))
       zmin=min(coor(3,ip1),coor(3,ip2),coor(3,ip3),coor(3,ip4))

       xmax=max(coor(1,ip1),coor(1,ip2),coor(1,ip3),coor(1,ip4))
       ymax=max(coor(2,ip1),coor(2,ip2),coor(2,ip3),coor(2,ip4))
       zmax=max(coor(3,ip1),coor(3,ip2),coor(3,ip3),coor(3,ip4))

       !
       !     Check against vmin,vmax
       !

       if(vmin(1)>xmax)cycle
       if(vmin(2)>ymax)cycle
       if(vmin(3)>zmax)cycle
       if(vmax(1)<xmin)cycle
       if(vmax(2)<ymin)cycle
       if(vmax(3)<zmin)cycle


       ineigh=eltoel(1,ielem)
       if(ineigh/=0)then
          if(lelem(ineigh)==0)then
             lelem(ineigh)=1_ip
             nstack=nstack+1
             lstack(nstack)=ineigh
          endif
       endif
       ineigh=eltoel(2,ielem)
       if(ineigh/=0)then
          if(lelem(ineigh)==0)then
             lelem(ineigh)=1_ip
             nstack=nstack+1
             lstack(nstack)=ineigh
          endif
       endif
       ineigh=eltoel(3,ielem)
       if(ineigh/=0)then
          if(lelem(ineigh)==0)then
             lelem(ineigh)=1_ip
             nstack=nstack+1
             lstack(nstack)=ineigh
          endif
       endif
       !
       !     Check the points of the element     
       !
       if(lmark(ip1)==0)then
          lmark(ip1)=1_ip
          nstackp2=nstackp2+1
          lstackp2(nstackp2)=ip1
          call length(ip1,ipnew,rsize,coor,rlen)
          if(rlen<tollen)then
             nstackp=nstackp+1
             lstackp(nstackp)=ip1
             rstackp(nstackp)=rlen
          endif
       endif

       if(lmark(ip2)==0)then
          lmark(ip2)=1_ip
          nstackp2=nstackp2+1
          lstackp2(nstackp2)=ip2
          call length(ip2,ipnew,rsize,coor,rlen)
          if(rlen<tollen)then
             nstackp=nstackp+1
             lstackp(nstackp)=ip2
             rstackp(nstackp)=rlen
          endif
       endif

       if(lmark(ip3)==0)then
          lmark(ip3)=1_ip
          nstackp2=nstackp2+1
          lstackp2(nstackp2)=ip3
          call length(ip3,ipnew,rsize,coor,rlen)
          if(rlen<tollen)then
             nstackp=nstackp+1
             lstackp(nstackp)=ip3
             rstackp(nstackp)=rlen
          endif
       endif

    enddo
    !
    !     Clean up lelem
    !
    do istack=1,nstack
       lelem(lstack(istack))=0_ip
    enddo
    !
    !     Clean up lmark
    !
    do istackp=1,nstackp2
       lmark(lstackp2(istackp))=0_ip
    enddo
    !
    !     Do we have close points
    !
    if(nstackp==0)then
       ipclos=0 
       return
    endif
    !
    !     Take out ipa,ipb and ipc
    ! 
    nstackp0=nstackp
    nstackp=0_ip
    do istackp=1,nstackp
       ipoin=lstackp(istackp)
       if(ipoin/=ipa .and. ipoin/=ipb .and. ipoin /=ipc)then
          nstackp=nstackp+1
          lstackp(nstackp)=ipoin
       endif
    enddo
    !
    !     Do we have still close points
    !
    if(nstackp==0)then
       ipclos=0 
       return
    endif
    !
    !     Take the minimum
    !
    rmin=rstackp(1)
    imin=1_ip

    do istackp=2,nstackp
       if(rstackp(istackp)<rmin)then
          rmin=rstackp(istackp)
          imin=istackp
       endif
    enddo

    ipclos=lstackp(imin)

  end subroutine gtclos3d

  subroutine advf3d(lface,nnofa,nface,ndim,npoin,coor,rsize,rnopo,lcart,lpsur,lmark,lptype,lpsid,lcell,ncell)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: nnofa,nface,ndim,ncell
    integer(ip),intent(in)    :: lface(nnofa,nface)
    integer(ip),intent(inout) :: npoin
    real(rp),pointer      :: coor(:,:),rsize(:)      
    real(rp),pointer      :: rfront(:),rnopo(:,:)
    integer(ip)           :: lcart(:),lpsur(:),lmark(:),lptype(:,:),lpsid(:)
    type(cell)            :: lcell(ncell)   
    integer(ip),pointer   :: lheap(:),lfahol(:),lfhole(:),lpofa(:),lfapo(:,:),lfront(:,:)
    integer(ip)           :: nfront,nfahol,nfapo,nheap,nfhole,iface,iter,ismall,ip1,ip2,ip3
    integer(ip)           :: ipnew,icart,ipcros 
    integer(4)            :: istat
    real(rp)              :: pnew(3)
    !
    !     This subroutine drives the volume advancing front like mesh generator
    !

    !
    !     Allocate front related quantities
    !
    nfront=nface 
    allocate(lfront(4,nfront),stat=istat)
    call memchk(zero,istat,memor_msh,'LFRONT','refsmo',lfront)
    allocate(rfront(nfront),stat=istat)
    call memchk(zero,istat,memor_msh,'RFRONT','refsmo',rfront)
    allocate(lheap(nfront),stat=istat)
    call memchk(zero,istat,memor_msh,'LHEAP','refsmo',lheap)
    allocate(lfapo(7,nfront),stat=istat)
    call memchk(zero,istat,memor_msh,'LFAPO','refsmo',lfapo)
    allocate(lfahol(max(1_ip,nfront/10_ip)),stat=istat)
    call memchk(zero,istat,memor_msh,'LFAHOL','refsmo',lfahol)
    allocate(lfhole(max(1_ip,nfront/10_ip)),stat=istat)
    call memchk(zero,istat,memor_msh,'LFHOLE','refsmo',lfhole)
    !
    !     Allocate local arrays
    !
    allocate(lpofa(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOFA','advf3d',lpofa)
    !
    !     Initialize nfapo, nfahol && nfhole
    !
    nfapo=0_ip
    nfahol=0_ip
    nfhole=0_ip
    !
    !     Compute the face size
    !
    do iface=1,nface 
       !
       !    Get size
       ! 
       call getsiz(lface,nnofa,nface,ndim,npoin,coor,iface,rfront)
       !
       !     Add it in the heap
       !
       call addfac3d(nfront,nfapo,nfahol,coor,ndim,npoin,nfhole,nheap,&
       lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole,iface)

    enddo
    !
    !    Set the pointer lfront->lheap
    !
    call initheap3d(lheap,nheap,lfront,nfront)
    !
    !     Loop on the front
    !
    iter=1_ip
    do

       if(nfront==0)exit
       !
       !     Get the smallest edge of the front
       !
       call gtsmall3d(lheap,nheap,rfront,ismall,lfront,nfront+nfhole)
       !
       !     This edge face is: 
       !
       ip1=lfront(1,ismall)
       ip2=lfront(2,ismall)
       ip3=lfront(3,ismall)
       !
       !     Compute new point
       ! 
       call gtnewp(nnofa,nfront,lfront,coor,ndim,pnew,rsize,ismall,npoin,rfront)
       !
       !     Add the point to the point array
       !
       ipnew=npoin+1
       !call resizep(ipnew,coor,rnopo,lcart,lpsur,rsize,lptri,lmark,lpofa,lptype,lpsid)
       npoin=ipnew
       !call gtelem(ipnew,coor,npoin,ndim,lcell,ncell,icart,rtol)
       lcart(ipnew)=icart
       !call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipnew,ndim,coor,rsuni)
       !
       !     Find the host element in the new mesh
       !
       !call cross(ip1,ip2,lface,nnofa,nface,lptri,coor,npoin,ndim,eltoel,rnopo,rnofa,&
       !     lfmark,ihostn,d1,d2,d3,lmark,lpofa,lfapo,nfapo,nfahol,lfront,nfront+nfhole,&
       !     ipcros,ipnew,rsize,lelem,ierr)
       !
       !     Do we have a crossed front?
       !
       if(ipcros/=0)then
          !
          !     Does the point belong to an active current edge
          !
         ! if(lpofa(ipcros)==0)then
             !
             !     Try to move ipcros to ipnew
             !
         !    call move(coor(:,ipnew),coor,ndim,npoin,lface,nnofa,nface,ipcros,rnopo,&
         !         coorold,nfold,npold,rnofaold,lfold,eltoelold,eltoel,lptri,lpsur,&
         !         rnofa,lelemold,lsurfold,isurf,rsize,lmark,lelem,lpofa,lfmark,lcart,&
         !         lcell,ncell,rtol,rsuni,ihostn)
         ! endif
          !
          !  This is the point to reach by regeneration
          !
          npoin=npoin-1
          ipnew=ipcros
       else
          !
          !     Insert the new point in the triangulation
          !
        !  call insert(nface,ihostn,nnofa,ndim,npoin,d1,d2,d3,ipnew,coor,&
        !       lface,eltoel,rnofa,lptri,lelem,lfmark,rnopo)
       endif
       !
       !     Regenerate the missing edges
       !
       !call recover(ip1,ip2,lface,nnofa,nface,lptri,coor,npoin,ndim,eltoel,rnopo,ipnew,&
       !     rnofa,lfmark,ienew,iview,ierr)
       !
       !     Update the front
       !
       !call newfront3d(nfapo,nfront,ismall,ipnew,npoin,nfahol,coor,ndim,nfhole,nheap,&
       !     lfront,lheap,lpofa,lfapo,rsize,rfront,lfahol,lfhole)
       !
       !     Swap ahead of the newly introduced element
       !
       !call swaploc(ienew,lface,nnofa,nface,coor,npoin,ndim,rnopo,lelem,lfmark,eltoel,lptri,rnofa)

       iter=iter+1_ip

    enddo

    call memchk(2_ip,istat,memor_msh,'LFHOLE','advf3d',lfhole)
    deallocate(lfhole,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFHOLE','advf3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFAHOL','advf3d',lfahol)
    deallocate(lfahol,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFAHOL','advf3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFAPO','advf3d',lfapo)
    deallocate(lfapo,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFAPO','advf3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHEAP','advf3d',lheap)
    deallocate(lheap,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHEAP','advf3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'RFRONT','advf3d',rfront)
    deallocate(rfront,stat=istat)
    if(istat/=0) call memerr(2_ip,'RFRONT','advf3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFRONT','advf3d',lfront)
    deallocate(lfront,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFRONT','advf3d',0_ip)

  end subroutine advf3d

  subroutine getsiz(lface,nnofa,nface,ndim,npoin,coor,iface,rheap)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: nnofa,nface,ndim,iface
    integer(ip),intent(in)    :: lface(nnofa,nface)
    integer(ip),intent(inout) :: npoin
    real(rp),intent(in)   :: coor(ndim,npoin)      
    real(rp),intent(inout):: rheap(nface)      
    integer(ip)           :: ip1,ip2,ip3   
    real(rp)              :: rx1,ry1,rz1,rx2,ry2,rz2,rx3,ry3,rz3
    real(rp)              :: rlx,rly,rlz,rnl1,rnl2,rnl3

    ip1=lface(1,iface)
    ip2=lface(2,iface)
    ip3=lface(3,iface)

    rx1=coor(1,ip1)
    ry1=coor(2,ip1)
    rz1=coor(3,ip1)
    rx2=coor(1,ip2)
    ry2=coor(2,ip2)
    rz2=coor(3,ip2)
    rx3=coor(1,ip3)
    ry3=coor(2,ip3)
    rz3=coor(3,ip3)
    ! 
    !     Get smallest side 
    !
    rlx=rx1-rx2
    rly=ry1-ry2
    rlz=rz1-rz2
    rnl1=sqrt(rlx*rlx+rly*rly+rlz*rlz)
    rlx=rx1-rx3
    rly=ry1-ry3
    rlz=rz1-rz3
    rnl2=sqrt(rlx*rlx+rly*rly+rlz*rlz)
    rlx=rx2-rx3
    rly=ry2-ry3
    rlz=rz2-rz3
    rnl3=sqrt(rlx*rlx+rly*rly+rlz*rlz)

    rheap(iface)=min(rnl1,rnl2,rnl3)

  end subroutine getsiz

  subroutine newfront3d(nfapo,nfront,ifront,ipnew,npoin,nfahol,coor,ndim,nfhole,nheap,&
       lfront,lheap,lpofa,lfapo,rsize,rfront,lfahol,lfhole)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)     :: ifront,ipnew,npoin,ndim
    real(rp), intent(in)        :: coor(ndim,npoin)
    integer(ip),pointer         :: lfront(:,:),lheap(:),lpofa(:),lfapo(:,:)
    integer(ip),pointer         :: lfahol(:),lfhole(:)
    real(rp),pointer            :: rsize(:),rfront(:)
    integer(ip), intent(inout)  :: nfapo,nfront,nfahol,nfhole,nheap
    integer(ip)                 :: ip1,ip2,ip1old,ip2old,ipos,iplace,idel,ifface
    integer(ip)                 :: nfaloc,lfaloc(3,100),nfface,ifaloc,iface,ipa,ipb 


    ip1old=lfront(1,ifront)
    ip2old=lfront(2,ifront)

    !
    !     Delete ifront in lfapo for ip1 and ip2, the front endpoints
    !
    call delfac3d(ifront,lfront,nfront,lfapo,nfapo,lpofa,npoin,nfahol,nfhole,&
         lfahol,lfhole)
    !
    !     Add new faces if appropriate
    !
    if(lpofa(ipnew)==0)then
       !
       !     First side
       !
       call addfac3d(nfront,nfapo,nfahol,coor,ndim,npoin,nfhole,nheap,&
       lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole,iface)
       !
       !     Second side
       !
       call addfac3d(nfront,nfapo,nfahol,coor,ndim,npoin,nfhole,nheap,&
       lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole,iface)

    else
       !
       !     Must delete doubly defined faces when adding (ip2old,ipnew) and (ipnew,ip1old)
       !     or add them if not found
       ! 
       ipos=lpofa(ipnew)
       if(ipos<=0)then
          write(*,*)'Error in newfront, ipos=',ipos 
          stop
       endif
       !
       !     Copy all the faces in a local array
       ! 
       nfaloc=0_ip

       do 

          nfface=lfapo(7,ipos)
          if(nfface>0)then

             do ifface=1,nfface
  
                nfaloc=nfaloc+1
                lfaloc(1,nfaloc)=lfapo(ifface,ipos)
                lfaloc(2,nfaloc)=ifface
                lfaloc(3,nfaloc)=ipos

             enddo
             exit
              
           else
         
              do ifface=1,6
  
                nfaloc=nfaloc+1
                lfaloc(1,nfaloc)=lfapo(ifface,ipos)
                lfaloc(2,nfaloc)=ifface
                lfaloc(3,nfaloc)=ipos

             enddo
             ipos=-nfface

          endif
          !
          !     DBG 
          !
          if(ipos<=0)then
             write(*,*)'Error in newfront, ipos:',ipos
             stop
          endif  

       enddo
       !
       !     First side
       !
       idel=0_ip
       do ifaloc=1,nfaloc
          iface=lfaloc(1,ifaloc)    
          ipa=lfront(1,iface)
          ipb=lfront(2,iface)

          if(ipa==ip2old .and. ipb==ipnew)then
             !
             !     Face (ipnew,ip2) found
             !
             call delheap3d(lfront(3,iface),lheap,nheap,lfront,nfront+nfhole,rfront)
             call delfac3d(iface,lfront,nfront,lfapo,nfapo,lpofa,npoin,nfahol,nfhole,&
                  lfahol,lfhole)
             idel=1_ip
          endif
       enddo
       !
       !     Do we have to add this face if not found?
       !
       if(idel==0)then

       call addfac3d(nfront,nfapo,nfahol,coor,ndim,npoin,nfhole,nheap,&
       lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole,iface)

       endif
       !
       !     Second side
       !
       idel=0_ip
       do ifaloc=1,nfaloc
          iface=lfaloc(1,ifaloc)    
          ipa=lfront(1,iface)
          ipb=lfront(2,iface)

          if(ipa==ipnew .and. ipb==ip1old) then

             call delheap3d(lfront(3,iface),lheap,nheap,lfront,nfront+nfhole,rfront)
             call delfac3d(iface,lfront,nfront,lfapo,nfapo,lpofa,npoin,nfahol,nfhole,&
                  lfahol,lfhole)
             idel=1_ip
          endif
       enddo
       !
       !     Do we have to add this face if not found?
       !
       if(idel==0)then

          call addfac3d(nfront,nfapo,nfahol,coor,ndim,npoin,nfhole,nheap,&
               lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole,iface)

       endif

    endif

  end subroutine newfront3d

  subroutine delfac3d(ifront,lfront,nfront,lfapo,nfapo,lpofa,npoin,nfahol,nfhole,&
       lfahol,lfhole)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)     ::  ifront,npoin 
    integer(ip), intent(inout)  ::  nfront,nfapo,nfahol,nfhole 
    integer(ip), intent(inout)  ::  lfront(3,nfront+nfhole),lfapo(3,nfapo+nfahol),lpofa(npoin) 
    integer(ip), pointer        ::  lfahol(:),lfhole(:)
    integer(ip)                 ::  ip1old,ip2old,ipos,nfface,ip1,ip2,jpos,jplace,irow,iposn
    !
    !     Delete ifront in lfapo for ip1 and ip2, the front endpoints
    !
    ip1old=lfront(1,ifront)
    ip2old=lfront(2,ifront)
    !
    !     Delete for ip1old
    !
    ipos=lpofa(ip1old)
    iposn=0_ip
    jpos=0_ip
    jplace=0_ip
    irow=0_ip

    do
       nfface=lfapo(3,ipos)
       if(nfface==2)then
          if(lfapo(1,ipos)==ifront)then
             lfapo(1,ipos)=lfapo(2,ipos)
          else if(lfapo(2,ipos)==ifront)then   
          else
             lfapo(jplace,jpos)=lfapo(2,ipos)
          endif
          lfapo(2,ipos)=0_ip
          lfapo(3,ipos)=1_ip
          exit
       else if(nfface==1)then 
          if(lfapo(1,ipos)==ifront)then 
          else
             lfapo(jplace,jpos)=lfapo(1,ipos)
          endif
          lfapo(1,ipos)=0_ip
          lfapo(3,ipos)=0_ip
          irow=ipos
          if(iposn/=0)then
             lfapo(3,iposn)=2_ip
          endif
          exit
       else
          if(lfapo(1,ipos)==ifront)then
             jpos=ipos
             jplace=1_ip
          else if(lfapo(2,ipos)==ifront)then
             jpos=ipos
             jplace=2_ip
          endif

          iposn=ipos 
          ipos=-nfface
       endif

    enddo
    !
    !     Do we have an empty row
    !
    if(irow/=0)then 

       nfahol=nfahol+1
       call memrea(nfahol,memor_msh,'LFHOLE','delfac',lfahol)
       lfahol(nfahol)=irow
       nfapo=nfapo-1 
       ipos=lpofa(ip1old)
       if(lfapo(3,ipos)==0)then
          lpofa(ip1old)=-1_ip
       endif

    endif
    !
    !     Delete for ip2old
    !
    ipos=lpofa(ip2old)
    iposn=0_ip
    jpos=0_ip
    jplace=0_ip
    irow=0_ip

    do
       nfface=lfapo(3,ipos)
       if(nfface==2)then
          if(lfapo(1,ipos)==ifront)then
             lfapo(1,ipos)=lfapo(2,ipos)
          else if(lfapo(2,ipos)==ifront)then   
          else
             lfapo(jplace,jpos)=lfapo(2,ipos)
          endif
          lfapo(2,ipos)=0_ip
          lfapo(3,ipos)=1_ip
          exit
       else if(nfface==1)then 
          if(lfapo(1,ipos)==ifront)then 
          else
             lfapo(jplace,jpos)=lfapo(1,ipos)
          endif
          lfapo(1,ipos)=0_ip
          lfapo(3,ipos)=0_ip
          irow=ipos
          if(iposn/=0)then
             lfapo(3,iposn)=2_ip
          endif
          exit
       else
          if(lfapo(1,ipos)==ifront)then
             jpos=ipos
             jplace=1_ip
          else if(lfapo(2,ipos)==ifront)then
             jpos=ipos
             jplace=2_ip
          endif

          iposn=ipos
          ipos=-nfface
       endif

    enddo
    !
    !     Do we have an empty row
    !
    if(irow/=0)then 

       nfahol=nfahol+1
       call memrea(nfahol,memor_msh,'LFAHOL','delfac',lfahol)
       lfahol(nfahol)=irow
       nfapo=nfapo-1
       ipos=lpofa(ip2old)
       if(lfapo(3,ipos)==0)then
          lpofa(ip2old)=-1_ip
       endif

    endif
    !
    !     Delete ifront in lfront 
    !
    nfhole=nfhole+1
    call memrea(nfhole,memor_msh,'LFHOLE','delfac',lfhole)
    lfhole(nfhole)=ifront

    nfront=nfront-1    

  end subroutine delfac3d

  subroutine addfac3d(nfront,nfapo,nfahol,coor,ndim,npoin,nfhole,nheap,&
       lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole,iface)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)     :: ndim,npoin,iface
    real(rp), intent(in)        :: coor(ndim,npoin)
    integer(ip), intent(inout)  :: nfapo,nfront,nfahol,nfhole,nheap
    integer(ip), pointer        :: lfront(:,:),lpofa(:),lfapo(:,:),lfahol(:)
    integer(ip), pointer        :: lheap(:),lfhole(:)
    real(rp), pointer           :: rfront(:)
    integer(ip)                 :: ipos,nfface,nfnew,iposn,iposf,ipa,ipb
    real(rp)                    :: rx,ry,rz,rlen
    !
    !     This subroutine adds a face to the front and to lfapo
    !
    !
    !     Do we have holes in the datastructure?
    ! 
    if(nfhole>0)then

       iposf=lfhole(nfhole)
       nfhole=nfhole-1

    else 

       iposf=nfront+1 
       call memrea(iposf,memor_msh,'LFRONT','addfac',lfront)
       call memrea(iposf,memor_msh,'RFRONT','addfac',rfront)
       call memrea(iposf,memor_msh,'LHEAP','addfac',lheap)

    endif
    !
    !     Introduce face in heap
    !
    call addheap3d(lheap,nheap,rfront,iposf,lfront,nfront+nfhole)
    !
    !     Update face to point datastructure
    !
    ipos=lpofa(ipa)
    !
    !     Do we have already faces attached to ipa?
    !
    if(ipos==0)then
       if(nfahol>0)then
          iposn=lfahol(nfahol)
          nfahol=nfahol-1
       else
          iposn=nfapo+1
          call memrea(iposn,memor_msh,'LFAPO','addfac',lfapo)
       endif

       nfapo=nfapo+1
       lpofa(ipa)=iposn
       lfapo(1,iposn)=iposf
       lfapo(3,iposn)=1_ip

    else

       do 
          nfface=lfapo(3,ipos)
          if(nfface==2)then
             if(nfahol>0)then
                iposn=lfahol(nfahol)
                nfahol=nfahol-1
             else
                iposn=nfapo+1
                call memrea(iposn,memor_msh,'LFAPO','addfac',lfapo)
             endif

             nfapo=nfapo+1
             lfapo(1,iposn)=iposf
             lfapo(3,iposn)=1_ip
             lfapo(3,ipos)=-iposn
             exit
          else if(nfface==1)then
             lfapo(2,ipos)=iposf
             lfapo(3,ipos)=2_ip
             exit
          else  
             ipos=-nfface
          endif
       enddo

    endif

    ipos=lpofa(ipb)

    if(ipos==0)then
       if(nfahol>0)then
          iposn=lfahol(nfahol)
          nfahol=nfahol-1
       else
          iposn=nfapo+1
          call memrea(iposn,memor_msh,'LFAPO','addfac',lfapo)
       endif

       nfapo=nfapo+1
       lpofa(ipb)=iposn
       lfapo(1,iposn)=iposf
       lfapo(3,iposn)=1_ip

    else

       do 
          nfface=lfapo(3,ipos)
          if(nfface==2)then
             if(nfahol>0)then
                iposn=lfahol(nfahol)
                nfahol=nfahol-1
             else
                iposn=nfapo+1
                call memrea(iposn,memor_msh,'LFAPO','addfac',lfapo)
             endif

             nfapo=nfapo+1
             lfapo(1,iposn)=iposf
             lfapo(3,iposn)=1_ip
             lfapo(3,ipos)=-iposn
             exit
          else if(nfface==1)then
             lfapo(2,ipos)=iposf
             lfapo(3,ipos)=2_ip
             exit
          else  
             ipos=-nfface
          endif
       enddo
    endif

  end subroutine addfac3d

  subroutine gtnewp(nnofa,nfront,lfront,coor,ndim,pnew,rsize,ifront,npoin,rheap)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: nnofa,ndim,ifront,npoin,nfront
    integer(ip),intent(in)    :: lfront(nnofa,nfront)
    real(rp),intent(in)   :: coor(ndim,npoin),rheap(nfront)     
    real(rp),intent(in)   :: rsize(npoin)     
    real(rp),intent(inout):: pnew(3)     
    integer(ip)           :: ip1,ip2,ip3  
    real(rp)              :: c13,rx,ry,rz,rsiz,c10
    real(rp)              :: dx12,dy12,dz12,dx13,dy13,dz13,rnx,rny,rnz,rnl
    !
    !     This subroutine computes the best point
    !   
    c13=1.0d+00/3.0d+00
    c10=1.0d+00
    ip1=lfront(1,ifront)
    ip2=lfront(2,ifront)
    ip3=lfront(3,ifront)
    !
    !     Get barycenter
    ! 
    rx=c13*(coor(1,ip1)+coor(1,ip2)+coor(1,ip3))
    ry=c13*(coor(2,ip1)+coor(2,ip2)+coor(2,ip3))
    rz=c13*(coor(3,ip1)+coor(3,ip2)+coor(3,ip3))
    !
    !     Get size
    !
    rsiz=c13*(rsize(ip1)+rsize(ip2)+rsize(ip3))
    !
    !     Get normal
    !
    dx12=coor(1,ip2)-coor(1,ip1)
    dy12=coor(2,ip2)-coor(2,ip1)
    dz12=coor(3,ip2)-coor(3,ip1)

    dx13=coor(1,ip3)-coor(1,ip1)
    dy13=coor(2,ip3)-coor(2,ip1)
    dz13=coor(3,ip3)-coor(3,ip1)

    rnx= dy12*dz13-dz12*dy13
    rny=-dx12*dz13+dz12*dx13
    rnz= dx12*dy13-dy12*dx13

    rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
    rnl=c10/rnl

    rnx=rnx*rnl
    rny=rny*rnl
    rnz=rnz*rnl
    !
    !     Extrude it with the normal
    !
    pnew(1)=rx+rsiz*rnx 
    pnew(2)=rx+rsiz*rnx 
    pnew(3)=rx+rsiz*rnx 

  end subroutine gtnewp


  subroutine rmv3d(nnode,nelem,elem,lhashf2,lhashf1,nhashf,eltoel)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: nnode,nhashf
    integer(ip),intent(inout) :: nelem
    integer(ip),intent(in)    :: lhashf2(nhashf+1),lhashf1(2,*),eltoel(nnode,nelem)
    integer(ip),intent(inout) :: elem(nnode,nelem)
    integer(ip),pointer       :: lstack(:)
    logical(lg),pointer       :: lmark(:)
    integer(ip)               :: ielem,inode,imin,imax,ineigh,istack,nstack   
    integer(ip)               :: ip1,ip2,ip3,isum,isto,nelem0   
    integer(4)                :: istat
    integer(ip)               :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/)) 
    !
    !     This sub removes the outer mesh
    !
    allocate(lstack(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','rmv3d',lstack)
    allocate(lmark(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','rmv3d',lmark)
    !
    !     Get an outer element
    !   
    do ielem=1,nelem
       do inode=1,nnode
          ineigh=eltoel(inode,ielem)
          if(ineigh==0)goto 10
       enddo
    enddo
    !
    !     Error!
    ! 
    write(*,*)'Did not find outer element'
    stop

10  continue   

    lstack(1)=ielem
    nstack=1_ip
    istack=0_ip
    lmark(ielem)=.true.

    do 
       !
       !     Increase istack
       !
       if(istack==nstack)exit 
       istack=istack+1_ip
       !
       !     Get element
       !
       ielem=lstack(istack)
       do inode=1,nnode
          ineigh=eltoel(inode,ielem)
          !
          !   Is there a neighbor?
          !
          if(ineigh==0)cycle 
          !
          !     Has this neighbor been marked?
          !   
          if(lmark(ineigh).eqv. .true.)cycle
          !
          !     Check in the hash table
          !
          ip1=elem(ltab(1,inode),ielem)
          ip2=elem(ltab(2,inode),ielem)
          ip3=elem(ltab(3,inode),ielem)

          isum=ip1+ip2+ip3
          imin=min(ip1,ip2,ip3)
          imax=max(ip1,ip2,ip3)

          do isto=lhashf2(isum),lhashf2(isum+1)-1_ip
             if(lhashf1(1,isto)==imin .and. lhashf1(2,isto)==imax)then    
                cycle
             endif
          enddo
          !
          !     This face is not a boundary face
          !     Mark ineigh and add it to the stack
          !
          nstack=nstack+1_ip
          lstack(nstack)=ineigh
          lmark(ineigh)=.true.

       enddo

    enddo
    !
    !     Now compress the mesh
    !
    nelem0=nelem
    nelem=0_ip
    do ielem=1,nelem0
       if(lmark(ielem).eqv. .false.)then
          nelem=nelem+1_ip
          elem(1,nelem)=elem(1,ielem)
          elem(2,nelem)=elem(2,ielem)
          elem(3,nelem)=elem(3,ielem)
          elem(4,nelem)=elem(4,ielem)
       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'LSTACK','rmv3d',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','rmv3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','rmv3d',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','rmv3d',0_ip)

  end subroutine rmv3d


  subroutine outGidpnt(ndim,npoin,coor,bbox)
    use def_kintyp, only       : ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: npoin,ndim
    real(rp),intent(in)         :: coor(ndim,npoin),bbox(ndim,2)
    real(rp)                    :: rx,ry,rz
    integer(ip)                 :: ipoin

    open(unit=50,file='advpGid.msh',status='unknown')
    rewind 50
    write(50,7)
    write(50,2)
    write(50,3)

    do  ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(50,100)ipoin,rx,ry,rz
    enddo

    write(50,4)
    write(50,5)
    do  ipoin=1,npoin
       write(50,200)ipoin,ipoin
    enddo
    write(50,6)

    write(50,1)
    write(50,2)
    write(50,3)

    write(50,100)ipoin,bbox(1,1),bbox(2,1),bbox(3,1)
    ipoin=ipoin+1
    write(50,100)ipoin,bbox(1,2),bbox(2,1),bbox(3,1)
    ipoin=ipoin+1
    write(50,100)ipoin,bbox(1,2),bbox(2,2),bbox(3,1)
    ipoin=ipoin+1
    write(50,100)ipoin,bbox(1,1),bbox(2,2),bbox(3,1)
    ipoin=ipoin+1
    write(50,100)ipoin,bbox(1,1),bbox(2,1),bbox(3,2)
    ipoin=ipoin+1
    write(50,100)ipoin,bbox(1,2),bbox(2,1),bbox(3,2)
    ipoin=ipoin+1
    write(50,100)ipoin,bbox(1,2),bbox(2,2),bbox(3,2)
    ipoin=ipoin+1
    write(50,100)ipoin,bbox(1,1),bbox(2,2),bbox(3,2)

    write(50,4)
    write(50,5)

    write(50,600)npoin+1,npoin+1,npoin+2,npoin+3,npoin+4,npoin+5,npoin+6,npoin+7,npoin+8

    write(50,6)

    close(50)

1   format('MESH dimension 3 ElemType  Hexahedra Nnode 8')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(10i10)
300 format(4i10)
400 format(5i10)
600 format(9i10)
500 format(i10,e20.10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
7   format('MESH dimension 3 ElemType Point Nnode 1')


  end subroutine outGidpnt

  subroutine outface(nnofa,nface,npoin,ndim,lface,coor)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim
    integer(ip), intent(in)      :: lface(nnofa,nface)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: i,icont
    !
    !     This sub displays a triangular surface mesh with its surface patch number 
    !

    open(unit=50,file='outface.msh',status='unknown')
    rewind 50

    write(50,1)
    write(50,2)
    write(50,3)

    do i=1,npoin
       write(50,100)i,coor(1,i),coor(2,i),coor(3,i)
    enddo

    write(50,4)
    write(50,5)

    icont=0_ip
    do i=1, nface
       icont=icont+1
       write(50,300)icont,lface(1,i),lface(2,i),lface(3,i)
    enddo

    write(50,6)
    close(50)


1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
300 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')


  end subroutine outface


end module mod_vol


