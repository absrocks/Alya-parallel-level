module mod_delau

contains

  subroutine delaumsh(ndim,npoin,coor,nnode,elem,nelem,nnofa,nnosi,coold,npold,rsold,elold,neold)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    use mod_cart, only  : boxbin,cartp,intercart,gtelem
    use mod_mshtol, only  : ptoelm,tetote
    use mod_interp, only  : ctopnt,interpnt
    use mod_bin3d, only  : inibin,insertini
    use mod_optim, only  : optim3d
    implicit none
    integer(ip),intent(in)    :: ndim,nnode,nnofa,nnosi,npold,neold
    integer(ip),intent(in)    :: elold(nnode,neold)
    integer(ip),intent(inout) :: npoin,nelem
    real(rp)                  :: coold(ndim,npold),rsold(npold)
    real(rp),pointer          :: rsize(:)
    real(rp),pointer          :: coor(:,:),relem(:,:),coort(:,:)
    integer(ip),pointer       :: elem(:,:)
    integer(ip),pointer       :: lptet(:),lptetold(:),lmark(:),lhole(:),lmarkold(:)
    integer(ip),pointer       :: lbin1(:),lbin2(:),lstack(:),eltoel(:,:),eltoelold(:,:)
    integer(ip),pointer       :: lcart(:),lctop1(:),lctop2(:),lcmark(:),lcstack(:)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),pt1old(:),pt2old(:),lmarkp(:)
    type(cell),pointer        :: lcell(:)
    real(rp)                  :: bbox(3,2),bboxbin(3,2),csearch,rtol
    real(rp)                  :: dx,dy,dz,center(3),hmax,c05,c20,dxbin,dybin,dzbin
    real(rp)                  :: d1,d2,d3,d4,realtimenewpnt,realtimeinsertp
    integer(ip)               :: npoinf,ipoin,nboup,nhashf,ncell,npnew,npoinf0,npold1
    integer(ip)               :: iface,nxbin,nbin,ielem,nelem0,nstack,inode,ichk,jnode
    integer(ip)               :: ineigh,ip1,ip2,ip3,ip4,nybin,nzbin,npmax,nitermax,npoi0
    integer(ip)               :: iguess,ielold,icart,nhole,ihole,nele0,ierr,iter
    integer(4)                :: istat
    integer(ip)               :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))
    integer(ip),parameter     :: mstackp=500
    integer(ip)               :: lstackp(mstackp)
    integer(ip)               :: ela_tim(2),ir
    !
    !     Nullify private pointers
    !
    nullify(eltoel,lptet,lmark,lbin1,lbin2,lstack,ptoel1,ptoel2,relem,&
         pt1old,pt2old,eltoelold,lptetold,lmarkold) 
    !
    !     This is the main subroutine for Delaunay mesh generation
    !
    c05=0.5d+00
    c20=2.0d+00
    !
    !     Max points per cell
    !
    npmax=8_ip
    ierr=0_ip
    !
    !     Max iter for the close points
    !
    nitermax=5
    !
    !     Increase for the search size
    !
    csearch=2.0d+00
    !
    !     Do we have some input points?
    !
    if(npoin==0)then
       write(*,*)'Error in  delaumsh, no input point given'
       stop
    endif
    !
    !     Get the elements surrounding the points for the old mesh
    !
    call ptoelm(elold,neold,npold,nnode,pt1old,pt2old)
    !
    !     Get the elements surrounding elements for the old mesh
    !
    call tetote(elold,nnode,neold,pt1old,pt2old,npold,eltoelold)
    !
    !     Build the cartesian mesh with coold and elold
    !
    allocate(lcart(npold),stat=istat)
    call memchk(zero,istat,memor_msh,'LCART','delaumsh',lcart)
    call cartp(npold,ndim,elold,nnode,neold,npmax,lcart,pt1old,pt2old,&
         coold,lcell,ncell,rtol,0_ip,bbox)
    !
    !     Build the cell to old point pointers
    !
    call ctopnt(lcart,npold,ncell,lctop1,lctop2 )
    !
    !     Allocate the arrays for the cartesian mesh
    !
    allocate(lcmark(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LCMARK','delaumsh',lcmark)
    allocate(lcstack(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LCSTACK','delaumsh',lcstack)
    allocate(lmarkold(neold),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKOLD','delaumsh',lmarkold)
    !
    !     Copy imposed points in coort
    !
    npoi0=npoin
    allocate(coort(ndim,npoi0),stat=istat)
    call memchk(zero,istat,memor_msh,'COORT','delaumsh',coort)
    do ipoin=1,npoi0
       coort(1,ipoin)=coor(1,ipoin)
       coort(2,ipoin)=coor(2,ipoin)
       coort(3,ipoin)=coor(3,ipoin)
    enddo
    !
    !     Compute bbox of the points from the imposed points
    !
    call boxbin(coor,ndim,npoin,bbox)
    !
    !     Create bin
    !
    call inibin(bbox,bboxbin,lbin1,lbin2,nbin,nxbin,nybin,nzbin,npoin,dxbin,dybin,dzbin)
    !
    !     Multiply it by two
    ! 
    center(1)=c05*(bbox(1,2)+bbox(1,1))
    center(2)=c05*(bbox(2,2)+bbox(2,1))
    center(3)=c05*(bbox(3,2)+bbox(3,1))

    hmax=max(bbox(1,2)-bbox(1,1),bbox(2,2)-bbox(2,1),bbox(3,2)-bbox(3,1))

    bbox(1,1)=center(1)-c20*hmax
    bbox(2,1)=center(2)-c20*hmax
    bbox(3,1)=center(3)-c20*hmax
    bbox(1,2)=center(1)+c20*hmax
    bbox(2,2)=center(2)+c20*hmax
    bbox(3,2)=center(3)+c20*hmax
    !
    !     Create initial tetrahedra
    !
    call inivol(coor,elem,eltoel,ndim,npoin,bbox,lptet,nelem,nnode,relem)
    !
    !    Allocate lmark
    !
    allocate(lmark(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','delaumsh',lmark)
    !
    !     Add the imposed points (without checking for size compatibility)
    !
    allocate(lstack(npoi0),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','delaumsh',lstack)
    call memrea(npoin+npoi0,memor_msh,'COOR','delaumsh',coor)
    call memrea(npoin+npoi0,memor_msh,'LPTET','delaumsh',lptet)
    allocate(rsize(npoin+npoi0),stat=istat)
    call memchk(zero,istat,memor_msh,'RSIZE','delaumsh',rsize)
    allocate(lptetold(npoin+npoi0),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTETOLD','delaumsh',lptetold)
    allocate(lmarkp(npoin+npoi0),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARP','delaumsh',lmarkp)

    nstack=0_ip
    do ipoin=1,npoi0
       npoin=npoin+1
       coor(1,npoin)=coort(1,ipoin) 
       coor(2,npoin)=coort(2,ipoin) 
       coor(3,npoin)=coort(3,ipoin)
       !
       !     Interpolate size
       !
       if(ipoin==1)then
          iguess=0_ip
          icart=max(1_ip,ncell/2_ip)
       else
          iguess=lptetold(npoin-1)
          icart=lcart(npoin-1) 
       endif
       !
       !     Interpolate in the cartesian mesh
       !  
       call gtelem(npoin,coor,npoin,ndim,lcell,ncell,icart,rtol)
       lcart(npoin)=icart
       !
       !     Interpolate in the background mesh
       ! 
       call interpnt(nnode,nelem,elem,ndim,neold,elold,coor,coold,npoin,npold,&
            npoin,eltoelold,pt1old,pt2old,lcart,lcell,ncell,lcmark,&
            lcstack,lctop1,lctop2,nitermax,lstackp,mstackp,csearch,lmarkold,&
            ielold,d1,d2,d3,d4,0_ip,ierr)
       if(ierr==1)then
          write(*,*)'Error in delaunmsh, npoin=',npoin
          write(*,*)'x y z:',coor(1,npoin),coor(2,npoin),coor(3,npoin)
          return
       endif
       !
       !     Interpolate size
       !
       ip1=elold(1,ielold)
       ip2=elold(2,ielold)
       ip3=elold(3,ielold)
       ip4=elold(4,ielold)
       rsize(npoin)=d1*rsold(ip1)+d2*rsold(ip2)+&
            d3*rsold(ip3)+d4*rsold(ip4)
       lptetold(npoin)=ielold
       !
       !     Add to the stack
       !
       nstack=nstack+1
       lstack(nstack)=npoin 
    enddo
    call memchk(2_ip,istat,memor_msh,'COORT','delaumsh',coort)
    deallocate(coort,stat=istat)
    if(istat/=0) call memerr(2_ip,'COORT','delaumsh',0_ip)
    !
    !     Resize lcart
    !
    call memrea(npoin,memor_msh,'LCART','delaumsh',lcart)
    !
    !     Initialize nhole
    !
    nhole=0_ip
    allocate(lhole(100),stat=istat)
    call memchk(zero,istat,memor_msh,'LHOLE','delaumsh',lhole)
    !
    !     Prepare the bin for the imposed points
    ! 
    call insertini(coor,npoin,ndim,lbin1,lbin2,nbin,dxbin,dybin,dzbin,&
         nxbin,nybin,nzbin,bboxbin)    
    !
    !     Insert imposed points in the mesh
    !     
    call insertp(lstack,npoin,coor,elem,eltoel,nnode,nelem,rsize,nstack,lptet,&
         lbin1,lbin2,nbin,ndim,relem,bboxbin,nxbin,nybin,nzbin,dxbin,dybin,dzbin,&
         nhole,lhole,lmark,lmarkp)
    !
    !     And loop on point creation/insertion
    !
    iter=1_ip
    write(*,*)'Iteration number:',iter,'npoin=',npoin,'nelem=',nelem
    do
       !
       !     Remember current point number
       !
       npold1=npoin
       !
       !     Create new points
       ! 
       call system_clock(count=ela_tim(1),count_rate=ir) 
       call newpnt(coor,elem,ndim,npoin,nelem,nnode,lbin1,lbin2,nbin,rsize,&
            dxbin,dybin,dzbin,nxbin,nybin,nzbin,bboxbin,neold,elold,coold,&
            eltoelold,pt1old,pt2old,lcart,lcell,ncell,lctop1,lctop2,nitermax,&
            csearch,npold,lcmark,lcstack,rsold,lptetold,lptet,lmarkold,rtol,&
            ierr,lmarkp)
       call system_clock(count=ela_tim(2),count_rate=ir) 
       realtimenewpnt=(REAL(ela_tim(2),kind=ip)-REAL(ela_tim(1),kind=ip))/ REAL(ir,kind=ip)
       if(ierr==1)then
          goto 9999
       endif
       !
       !     Do we have new points?
       !
       if(npold1==npoin)exit
       !
       !     Fill the stack
       !
       call memrea(npoin-npold1,memor_msh,'LSTACK','delaumsh',lstack)
       nstack=0_ip
       do ipoin=npold1+1,npoin   
          nstack=nstack+1_ip
          lstack(nstack)=ipoin
       enddo
       !
       !     Insert points in the mesh
       !     
       call system_clock(count=ela_tim(1),count_rate=ir) 
       call insertp(lstack,npoin,coor,elem,eltoel,nnode,nelem,rsize,nstack,lptet,&
            lbin1,lbin2,nbin,ndim,relem,bboxbin,nxbin,nybin,nzbin,dxbin,dybin,dzbin,&
            nhole,lhole,lmark,lmarkp)
       realtimeinsertp=(REAL(ela_tim(2),kind=ip)-REAL(ela_tim(1),kind=ip))/ REAL(ir,kind=ip)
       !
       !     Clean up lmarkp
       !
       !do ipoin=1,npoin
       !   lmarkp(ipoin)=0_ip
       !enddo
       !
       !     Mark the points in an element connected to a corner
       !
       !do ielem=1,nelem
       !   if(elem(1,ielem)==0)cycle
       !   do inode=1,nnode
       !      ipoin=elem(inode,ielem)
       !      if(ipoin<9)then
       !         do jnode=1,nnode
       !            lmarkp(elem(jnode,ielem))=1_ip
       !         enddo
       !         exit
       !      endif
       !   enddo
       !enddo
       !
       !     Optimize the mesh
       !
       !call optim3d(nelem,ndim,npoin,nnode,elem,coor,eltoel,nnosi,rsize,&
       !     nnofa,lmarkp,lhole,nhole)
       !
       !     Clean up lmarkp
       !
       !do ipoin=1,npoin
       !   lmarkp(ipoin)=0_ip
       !enddo
       !
       !     Update lptet
       !
       !do ielem=1,nelem
       !   if(elem(1,ielem)==0)cycle
       !   lptet(elem(1,ielem))=ielem
       !   lptet(elem(2,ielem))=ielem
       !   lptet(elem(3,ielem))=ielem
       !   lptet(elem(4,ielem))=ielem
       !enddo  
       !
       !     Resize element based datastructure
       !
       !call memrea(nelem,memor_msh,'RELEM','delaumsh',relem)
       !call memrea(nelem,memor_msh,'LMARK','delaumsh',lmark)
       !
       !     Recompute center and radius
       !
       !call radius(elem,nnode,nelem,relem,coor,ndim,npoin) 
       !
       !     DBG
       !
       !call chkconf(eltoel,nnode,elem,nelem,ndim,npoin)

       iter=iter+1_ip
       write(*,*)'Iteration number:',iter,'npoin=',npoin,'nelem=',nelem
       write(*,*)'CPU time newpnt:',realtimenewpnt,' insertp:',realtimeinsertp 
    enddo
    !
    !     DBG
    !
    call chkconf(eltoel,nnode,elem,nelem,ndim,npoin)
    !
    !     Compress the mesh
    !
    do ihole=1,nhole
       elem(1,lhole(ihole))=0_ip
    enddo
    nele0=nelem
    nelem=0_ip
    do ielem=1,nele0
       if(elem(1,ielem)/=0)then
          nelem=nelem+1_ip
          elem(1,nelem)=elem(1,ielem)   
          elem(2,nelem)=elem(2,ielem)   
          elem(3,nelem)=elem(3,ielem)   
          elem(4,nelem)=elem(4,ielem)   
       endif
    enddo
    nhole=0_ip
    !
    !     Remove the elements touching the first eight points
    !
    nelem0=nelem
    nelem=0_ip
    do ielem=1,nelem0
       ichk=0
       do inode=1,nnode
          if(elem(inode,ielem)<9)then
             ichk=1_ip
          endif
       enddo
       if(ichk==0)then
          nelem=nelem+1_ip
          elem(1,nelem)=elem(1,ielem) 
          elem(2,nelem)=elem(2,ielem) 
          elem(3,nelem)=elem(3,ielem) 
          elem(4,nelem)=elem(4,ielem) 
       endif
    enddo
    !
    !     Get the elements surrounding the points
    !
    call ptoelm(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !    Resize lmark
    !
    call memrea(npoin,memor_msh,'LMARK','delaumsh',lmark)
    !
    !     Get the elements surrounding elements
    !
    call tetote(elem,nnode,nelem,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Clean up lmark
    !
    do ipoin=1,npoin
       lmarkp(ipoin)=0_ip
    enddo
    !
    !     Mark the boundary points
    !
    do ielem=1,nelem
       do inode=1,nnode
          ineigh=eltoel(inode,ielem)
          if(ineigh==0)then
             ip1=elem(ltab(1,inode),ielem)
             ip2=elem(ltab(2,inode),ielem)
             ip3=elem(ltab(3,inode),ielem)
             lmarkp(ip1)=1_ip
             lmarkp(ip2)=1_ip
             lmarkp(ip3)=1_ip
          endif
       enddo
    enddo
    call outdelau(nnode,nelem,elem,ndim,npoin,coor,rsize)
    !
    !     Optimize the mesh
    !
    call optim3d(nelem,ndim,npoin,nnode,elem,coor,eltoel,nnosi,rsize,nnofa,lmarkp,lhole,nhole)
    !
    !     Output the mesh
    !  
9999 continue 
    call outdelau(nnode,nelem,elem,ndim,npoin,coor,rsize)

    call memchk(2_ip,istat,memor_msh,'LMARKP','delaumsh',lmarkp)
    deallocate(lmarkp,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARKP','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHOLE','delaumsh',lhole)
    deallocate(lhole,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHOLE','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCART','delaumsh',lcart)
    deallocate(lcart,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCART','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','delaumsh',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','delaumsh',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTET','delaumsh',lptet)
    deallocate(lptet,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTET','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','delaumsh',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARKOLD','delaumsh',lmarkold)
    deallocate(lmarkold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARKOLD','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LBIN1','delaumsh',lbin1)
    deallocate(lbin1,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBIN1','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LBIN2','delaumsh',lbin2)
    deallocate(lbin2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBIN2','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','delaumsh',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOELOLD','delaumsh',eltoelold)
    deallocate(eltoelold,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOELOLD','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','delaumsh',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCTOP1','delaumsh',lctop1)
    deallocate(lctop1,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCTOP1','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCTOP2','delaumsh',lctop2)
    deallocate(lctop2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCTOP2','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCMARK','delaumsh',lcmark)
    deallocate(lcmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCMARK','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCSTACK','delaumsh',lcstack)
    deallocate(lcstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCSTACK','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'PT1OLD','delaumsh',pt1old)
    deallocate(pt1old,stat=istat)
    if(istat/=0) call memerr(2_ip,'PT1OLD','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'PT2OLD','delaumsh',pt2old)
    deallocate(pt2old,stat=istat)
    if(istat/=0) call memerr(2_ip,'PT2OLD','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCELL','delaumsh',lcell)
    deallocate(lcell,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCELL','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'RSIZE','delaumsh',rsize)
    deallocate(rsize,stat=istat)
    if(istat/=0) call memerr(2_ip,'RSIZE','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTETOLD','delaumsh',lptetold)
    deallocate(lptetold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTETOLD','delaumsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'RELEM','delaumsh',relem)
    deallocate(relem,stat=istat)
    if(istat/=0) call memerr(2_ip,'RELEM','delaumsh',0_ip)

  end subroutine delaumsh

  subroutine inivol(coor,elem,eltoel,ndim,npoin,bbox,lptet,nelem,nnode,relem)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)       :: ndim,nnode
    integer(ip),intent(inout)    :: npoin,nelem
    real(rp),intent(in)          :: bbox(ndim,2)
    integer(ip),pointer          :: elem(:,:),eltoel(:,:),lptet(:),lbin1(:),lbin2(:)
    real(rp),pointer             :: coor(:,:),relem(:,:)
    real(rp)                     :: c(3),rini,dx,dy,dz,c05 
    integer(4)                   ::  istat
    !
    !     This sub generates the first 5 elements of the new mesh
    !
    c05=0.5d+00 
    npoin=8_ip

    if(.not.associated(coor))then
       allocate(coor(ndim,npoin),stat=istat)
       call memchk(zero,istat,memor_msh,'COOR','inivol',coor)
    else
       call memrea(npoin,memor_msh,'COOR','inivol',coor)
    endif
    if(.not.associated(lptet))then
       allocate(lptet(npoin),stat=istat)
       call memchk(zero,istat,memor_msh,'LTET','inivol',lptet)
    else
       call memrea(npoin,memor_msh,'LTET','inivol',lptet)
    endif

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
    if(.not.associated(elem))then
       allocate(elem(nnode,nelem),stat=istat)
       call memchk(zero,istat,memor_msh,'ELEM','inivol',elem)
    else
       call memrea(nelem,memor_msh,'ELEM','inivol',elem)
    endif
    if(.not.associated(relem))then
       allocate(relem(5,nelem),stat=istat)
       call memchk(zero,istat,memor_msh,'RELEM','inivol',relem)
    else
       call memrea(nelem,memor_msh,'RELEM','inivol',relem)
    endif

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

    lptet(1)=1_ip
    lptet(4)=1_ip
    lptet(5)=1_ip
    lptet(2)=1_ip
    lptet(8)=2_ip
    lptet(7)=2_ip
    lptet(6)=3_ip
    lptet(3)=4_ip

    if(.not.associated(eltoel))then
       allocate(eltoel(nnode,nelem),stat=istat)
       call memchk(zero,istat,memor_msh,'ELTOEL','inivol',eltoel)
    else
       call memrea(nelem,memor_msh,'ELTOEL','inivol',eltoel)
    endif

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

    c(1)=(bbox(1,1)+bbox(1,2))*c05
    c(2)=(bbox(2,1)+bbox(2,2))*c05
    c(3)=(bbox(3,1)+bbox(3,2))*c05

    dx=bbox(1,2)-c(1)
    dy=bbox(2,2)-c(2)
    dz=bbox(3,2)-c(3)

    !rini=sqrt(dx*dx+dy*dy+dz*dz)
    rini=dx*dx+dy*dy+dz*dz

    relem(1,1)=c(1)
    relem(2,1)=c(2)
    relem(3,1)=c(3)
    relem(5,1)=rini
    relem(1,2)=c(1)
    relem(2,2)=c(2)
    relem(3,2)=c(3)
    relem(5,2)=rini
    relem(1,3)=c(1)
    relem(2,3)=c(2)
    relem(3,3)=c(3)
    relem(5,3)=rini
    relem(1,4)=c(1)
    relem(2,4)=c(2)
    relem(3,4)=c(3)
    relem(5,4)=rini
    relem(1,5)=c(1)
    relem(2,5)=c(2)
    relem(3,5)=c(3)
    relem(5,5)=rini

  end subroutine inivol

  subroutine insertp(lstack,npoin,coor,elem,eltoel,nnode,nelem,rsize,nstack,lptet,&
       lbin1,lbin2,nbin,ndim,relem,bboxbin,nx,ny,nz,dx,dy,dz,nhole,lhole,lmark,lmarkp)
    use def_kintyp, only :  ip,rp,lg
    use mod_bin3d, only :  gtclos,insertbin
    use mod_memchk
    use def_meshin, only : memor_msh
    implicit none
    integer(ip),intent(in)    :: nnode,ndim
    integer(ip),intent(in)    :: nx,ny,nz
    real(rp),intent(in)       :: dx,dy,dz,bboxbin(3,2)
    integer(ip),intent(inout) :: npoin,nelem,nbin,nhole,nstack
    integer(ip),intent(inout) :: lstack(nstack)
    integer(ip),intent(inout) :: lbin2(nbin+1),lbin1(*)
    real(rp),pointer          :: rsize(:),coor(:,:)
    integer(ip), pointer      :: elem(:,:),eltoel(:,:),lptet(:),lhole(:),lmarkp(:)
    integer(ip), pointer      :: lface(:,:),lmark(:),lenew(:),lmarkedg(:,:),lstackf(:)
    real(rp),pointer          :: relem(:,:),rface(:,:)
    integer(ip)               :: istack,itetinit,icav,ncav,ipclos,ielem,jstack 
    integer(ip)               :: ipoin,mface,badface,nface,nenew,nstack0,isto,ibin
    real(rp)                  :: rx,ry,rz,csliver,crelax,c10
    integer(ip),parameter     :: mcav=500
    integer(ip)               :: lcav(mcav)
    integer(4)                :: istat
    logical(lg)               :: insert
    integer(ip),save          :: iter=0  
    !
    !     This subroutine inserts a cloud of points in a Delaunay mesh
    !
    c10=1.0d+00
    !
    !     Set relaxation rate
    !
    crelax=0.8d+00
    !
    !      relem carries all the necessary geometry for the elements
    !            relem(1:3,ielem)  :    center of element
    !            relem(4,ielem)    :    old radius of element
    !            relem(5,ielem)    :    new radius of element
    !
    !
    !      rface carries all the necessary geometry for the cavity faces
    !            rface(1:3,iface)  :    center of element
    !            rface(4,iface)    :    volume of element
    !            rface(5,iface)    :    new radius of element
    !            rface(6-8,iface)  :    normal to face
    !
    !
    !
    !
    !     Given npoin resize arrays related to nelem
    !
    nenew=7*npoin
    call memrea(nenew,memor_msh,'RELEM','insertp',relem) 
    call memrea(nenew,memor_msh,'ELEM','insertp',elem) 
    call memrea(nenew,memor_msh,'ELTOEL','insertp',eltoel) 
    call memrea(nenew,memor_msh,'LMARK','insertp',lmark)
    !  
    !     Allocate local database
    !
    allocate(lface(6,100),stat=istat)
    call memchk(zero,istat,memor_msh,'LFACE','insertp',lface)
    allocate(rface(8,100),stat=istat)
    call memchk(zero,istat,memor_msh,'RFACE','insertp',rface)
    allocate(lenew(100),stat=istat)
    call memchk(zero,istat,memor_msh,'LENEW','insertp',lenew)
    allocate(lmarkedg(4,100),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKEDG','insertp',lmarkedg)
    allocate(lstackf(100),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACKF','insertp',lstackf)
    !
    !     Generate random insertion
    !
    call randins(lstack,nstack)
    !
    !     Initialize csliver
    !
    csliver=c10 
    jstack=0_ip
    !
    !     Loop on iterations
    !
    do
       !
       !     Loop on stack
       !
       do istack=1,nstack 

          jstack=jstack+1
          ipoin=lstack(istack)
          !write(*,*)ipoin,istack
          !if(nelem>=160641)then
          !write(*,*)ipoin,istack,eltoel(1,160641),eltoel(2,160641),eltoel(3,160641),eltoel(4,160641)
          !endif

          !if(npoin>1747)then
          !   if(lptet(1747)>0)then
          !   write(*,*)lptet(1747)
          !   write(*,*)elem(1,lptet(1747)),elem(2,lptet(1747)),elem(3,lptet(1747)),elem(4,lptet(1747)) 
          !endif 
          !endif
          !
          !     Get a close point with the bin
          !
          call gtclos(ipoin,lbin1,lbin2,ipclos,nbin,coor,ndim,bboxbin,&
               dx,dy,dz,nx,ny,nz,npoin)
          !
          !     Get element associated with this close point already in the mesh
          !
          itetinit=lptet(ipclos)
          !
          !     DBG
          !
          if(itetinit<=0)then
             write(*,*)'Error in insertp,itetinit=',itetinit
             stop
          endif
          if(elem(1,itetinit)==0)then
             write(*,*)'Error in insertp,itetinit=',itetinit
             write(*,*)'elem(1,itetinit)=0'  
             stop
          endif
          !
          !     Find element(s) containing ipoin
          !
          call findelem(itetinit,elem,nelem,npoin,coor,coor(:,ipoin),ndim,&
               nnode,eltoel,ncav,mcav,lcav)
          !
          !     Those elements belong to the cavity, mark them
          !
          do icav=1,ncav
             ielem=lcav(icav)
             lmark(ielem)=1_ip
             !
             !     Compute the old radius for those elements
             !
             rx=coor(1,ipoin)-relem(1,ielem)
             ry=coor(2,ipoin)-relem(2,ielem)
             rz=coor(3,ipoin)-relem(3,ielem)
             relem(4,ielem)=rx*rx+ry*ry+rz*rz

          enddo
          !
          !     Grow the cavity 
          !
          call growcav(ipoin,npoin,nelem,lmark,lcav,ncav,eltoel,nnode,coor,ndim,relem)
          !
          !     Output the elements of cavity
          !
          !if(iter==13 .and. ipoin==26761)call outelem(nnode,nelem,elem,ndim,npoin,coor,lcav,ncav)
          !call outelem(nnode,nelem,elem,ndim,npoin,coor,lcav,ncav)
          !
          !     Get the faces of the current cavity
          !   
          call gtface(nelem,lmark,lcav,ncav,eltoel,nnode,nface,lface,mface,elem)
          !
          !     Resize rface && lenew
          !
          call memrea(nface,memor_msh,'RFACE','insertp',rface) 
          call memrea(nface,memor_msh,'LENEW','insertp',lenew) 
          call memrea(nface,memor_msh,'LMARKEDG','insertp',lmarkedg) 
          call memrea(nface,memor_msh,'LSTACKF','insertp',lstackf) 
          !
          !     Correct the cavity
          ! 
          call corcav(ipoin,nelem,lmark,lcav,ncav,eltoel,nnode,npoin,lface,rsize,&
               coor,rface,relem,badface,elem,ndim,nface,insert,csliver,lmarkp)
          !
          !     Is is possible to insert ipoin
          !
          if(insert .eqv. .false.)cycle
          !
          !     Mark lstack as inserted
          !
          lstack(istack)=0_ip 
          !
          !     Output the faces of cavity
          !
          !if(iter==13 .and. ipoin==26761)call outcav(ndim,npoin,coor,nface,lface)
          !call outcav(ndim,npoin,coor,nface,lface)
          !
          !     Insert the point in the bin
          !      
          call insertbin(ipoin,coor,npoin,ndim,lbin1,lbin2,nbin,dx,dy,dz,nx,ny,nz,bboxbin) 
          !
          !     Remesh the cavity
          ! 
          call newcav(nelem,lmark,eltoel,nnode,npoin,badface,lhole,nhole,nface,lenew,&
               lface,relem,rface,ipoin,elem,lptet,coor,jstack,ndim,ncav,lcav,lmarkedg,&
               lstackf,nenew)

          !call outdelau(nnode,nelem,elem,ndim,npoin,coor,rsize)
          
       enddo
       !
       !     DBG check the mesh
       ! 
       !call chkconf(eltoel,nnode,elem,nelem,ndim,npoin)

       nstack0=nstack
       nstack=0_ip
       do istack=1,nstack0     
          if(lstack(istack)/=0)then
             nstack=nstack+1
             lstack(nstack)=lstack(istack)
          endif
       enddo
       !
       !     Do we have still points to insert?
       !
       if(nstack==0)exit 
       !
       !     Relax csliver
       !
       csliver=csliver*crelax

    enddo
    !
    !     DBG
    !
    do ibin=1,nbin
       do isto=lbin2(ibin),lbin2(ibin+1)-1
          if(lbin1(isto)==0)then
             write(*,*)'Error lbin1=0',isto,ibin
             stop
          endif
       enddo
    enddo

    iter=iter+1_ip

    call memchk(2_ip,istat,memor_msh,'LFACE','insertp',lface)
    deallocate(lface,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFACE','insertp',0_ip)
    call memchk(2_ip,istat,memor_msh,'RFACE','insertp',rface)
    deallocate(rface,stat=istat)
    if(istat/=0) call memerr(2_ip,'RFACE','insertp',0_ip)
    call memchk(2_ip,istat,memor_msh,'LENEW','insertp',lenew)
    deallocate(lenew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LENEW','insertp',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARKEDG','insertp',lmarkedg)
    deallocate(lmarkedg,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARKEDG','insertp',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACKF','insertp',lstackf)
    deallocate(lstackf,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACKF','insertp',0_ip)

  end subroutine insertp

  subroutine randins(lpoin,npoin)
    use def_kintyp, only :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: npoin
    integer(ip),intent(inout) :: lpoin(npoin)
    integer(ip)               :: ipoin,l,itemp
    real(rp)                  :: rl

    do ipoin=1,npoin
       call random_number(rl)
       l=floor((npoin-1)*rl)+1_ip
       itemp=lpoin(ipoin)
       lpoin(ipoin)=lpoin(l)
       lpoin(l)=itemp    
    enddo

  end subroutine randins

  subroutine growcav(ipoin,npoin,nelem,lmark,lcav,ncav,eltoel,nnode,coor,ndim,&
       relem)
    use def_kintyp, only :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: npoin,nelem,nnode,ndim,ipoin
    integer(ip),intent(inout) :: ncav
    real(rp),intent(in)       :: coor(ndim,npoin)
    real(rp),intent(inout)    :: relem(5,nelem)
    integer(ip),intent(inout) :: lmark(nelem),lcav(nelem)
    integer(ip),intent(in)    :: eltoel(nnode,nelem)
    integer(ip)               :: icav,inode,ineigh,ielem,ieleout
    integer(ip)               :: leleout(1000),neleout
    real(rp)                  :: rx,ry,rz,rdist
    !
    !     This subroutine grows the cavity with the Delaunay criterion
    !
    !
    !     On output, the elements of the cavity are marked with 1 in lmark
    ! 
    icav=0_ip
    neleout=0_ip 
    !
    !     Loop on the elements of the cavity, initialized by the base
    !
    do 
       !
       !     Are we done yet?
       !
       if(icav==ncav)exit  
       icav=icav+1_ip
       !
       !     Get the element
       !
       ielem=lcav(icav)
       !
       !     Loop on neighbors
       ! 
       do inode=1,nnode

          ineigh=eltoel(inode,ielem) 
          !
          !     Is it a boundary element?
          ! 
          if(ineigh==0)cycle 
          !
          !     Has the neighbor been already tested?
          ! 
          if(lmark(ineigh)==1)cycle
          !
          !     Mark the neighbor 
          ! 
          lmark(ineigh)=1_ip
          !
          !     Get distance to center
          !  
          rx=coor(1,ipoin)-relem(1,ineigh)
          ry=coor(2,ipoin)-relem(2,ineigh)
          rz=coor(3,ipoin)-relem(3,ineigh)
          rdist=rx*rx+ry*ry+rz*rz
          !
          !     And the test with the radius squared
          !
          if(rdist>relem(5,ineigh))then
             !
             !     This element does not belong to the cavity.
             !     Remember it for deletion 
             !
             neleout=neleout+1_ip
             leleout(neleout)=ineigh
          else
             !
             !     This element does belong to the cavity.
             !
             ncav=ncav+1 
             lcav(ncav)=ineigh
             relem(4,ineigh)=rdist

          endif

       enddo

    enddo
    !
    !     Remark the outer elements
    !
    do ieleout=1,neleout
       ielem=leleout(ieleout) 
       lmark(ielem)=0_ip 
    enddo

  end subroutine growcav

  subroutine gtface(nelem,lmark,lcav,ncav,eltoel,nnode,nface,lface,mface,elem)
    use def_kintyp, only :  ip,rp,lg
    use mod_memchk
    use def_meshin, only: memor_msh
    implicit none
    integer(ip),intent(in)    :: nelem,nnode,ncav,mface
    integer(ip),intent(in)    :: lmark(nelem),lcav(ncav),elem(nnode,nelem)
    integer(ip),intent(inout) :: nface,eltoel(nnode,nelem)
    integer(ip),pointer       :: lface(:,:) 
    integer(ip)               :: icav,ielem,ineigh,inode 
    integer(ip)               :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))
    !
    !     This subroutine finds the faces of the cavity
    !
    nface=0_ip
    !
    !     Loop on elements in the cavity
    !
    do icav=1,ncav
       !
       !     Get element
       !
       ielem=lcav(icav)
       !
       !     Loop on neighbors
       !
       do inode=1,nnode
          ineigh=eltoel(inode,ielem)
          !
          !     Do we have a boundary element or 
          !     Is the neighbor outside the cavity?
          !
          if(ineigh==0)then
             !
             !     Store face database
             !
             nface=nface+1_ip   
             call memrea(nface,memor_msh,'LFACE','gtface',lface) 
             lface(1,nface)=elem(ltab(1,inode),ielem)
             lface(2,nface)=elem(ltab(2,inode),ielem)
             lface(3,nface)=elem(ltab(3,inode),ielem)
             lface(4,nface)=ielem
             lface(5,nface)=inode
             lface(6,nface)=ineigh
             !
             !     Set pointer element to face
             !
             eltoel(inode,ielem)=-nface
             !
             !     Is ineigh outside the cavity?
             !
          else if(lmark(ineigh)==0)then
             !
             !     Store face database
             !
             nface=nface+1_ip   
             call memrea(nface,memor_msh,'LFACE','gtface',lface) 
             lface(1,nface)=elem(ltab(1,inode),ielem)
             lface(2,nface)=elem(ltab(2,inode),ielem)
             lface(3,nface)=elem(ltab(3,inode),ielem)
             lface(4,nface)=ielem
             lface(5,nface)=inode
             lface(6,nface)=ineigh
             !
             !     Set pointer element to face
             !
             eltoel(inode,ielem)=-nface
          endif
       enddo
    enddo

  end subroutine gtface

  subroutine corcav(ipoin,nelem,lmark,lcav,ncav,eltoel,nnode,npoin,lface,rsize,coor,rface,&
       relem,badface,elem,ndim,nface,insert,csliver,lmarkp)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: npoin,nelem,nnode,ndim,ipoin
    integer(ip),intent(in)    :: elem(nnode,nelem)
    integer(ip),intent(inout) :: ncav,badface,nface,lmarkp(npoin)
    logical(lg),intent(inout) :: insert
    real(rp),intent(in)       :: rsize(npoin),coor(ndim,npoin),csliver
    real(rp),intent(inout)    :: relem(5,nelem)
    real(rp),pointer          :: rface(:,:) 
    integer(ip),pointer       :: lface(:,:) 
    integer(ip),intent(inout) :: lmark(nelem),lcav(nelem)
    integer(ip),intent(inout) :: eltoel(nnode,nelem)
    integer(ip)               :: icav,ip1,ip2,ip3,inode,ineigh,iface,ielem
    integer(ip)               :: jnode,jneigh,iview,nstackf,nstacke,istacke,jpoin
    real(rp)                  :: hloc,sliver,correc,rx1,ry1,rz1,rx2,ry2,rz2,rx,ry,rz
    integer(ip)               :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))
    integer(ip),parameter     :: mstack=500
    integer(ip)               :: lstackf(mstack),lstacke(mstack)  
    !
    !     This sub obtains a valid cavity
    !
    !call outcav(ndim,npoin,coor,nface,lface)
    !
    !     Set tolerance
    !
    correc=1.0d-03
    hloc=rsize(ipoin)
    sliver=hloc*hloc*hloc*correc*csliver
    !
    !     Set insert
    !
    insert=.true.
    iface=0_ip
    badface=0_ip
    !
    !     Global loop on faces and points
    !
    do 
       !
       !     Loop on faces of the cavity
       !
       do 

          if(iface==nface)exit
          iface=iface+1_ip
          !
          !     Is this face still valid?
          !
          if(lface(6,iface)/=-1)then 
             !
             !     Get face information
             !
             ip1=lface(1,iface)
             ip2=lface(2,iface)
             ip3=lface(3,iface)
             ielem=lface(4,iface)
             inode=lface(5,iface)  
             ineigh=lface(6,iface)  
             !
             !     Get normal
             !
             rx1=coor(1,ip2)-coor(1,ip1)   
             ry1=coor(2,ip2)-coor(2,ip1)   
             rz1=coor(3,ip2)-coor(3,ip1)   

             rx2=coor(1,ip3)-coor(1,ip1)   
             ry2=coor(2,ip3)-coor(2,ip1)   
             rz2=coor(3,ip3)-coor(3,ip1)   

             rface(6,iface)= ry1*rz2-rz1*ry2
             rface(7,iface)=-rx1*rz2+rz1*rx2
             rface(8,iface)= rx1*ry2-ry1*rx2
             !
             !     Get additional vector
             !
             rx=coor(1,ip1)-coor(1,ipoin)
             ry=coor(2,ip1)-coor(2,ipoin)
             rz=coor(3,ip1)-coor(3,ipoin)
             !
             !     Get volume of the new element
             !
             rface(4,iface)=rface(6,iface)*rx+rface(7,iface)*ry+rface(8,iface)*rz
             !
             !     Test for sliver
             !
             if(rface(4,iface)>sliver)cycle
             !
             !     Delete elements of the initial cavity without removing the one of the base
             !
             !
             !     Restore eltoel
             !
             eltoel(inode,ielem)=ineigh
             !
             !     Delete ielem from the cavity
             !
             lmark(ielem)=0
             !
             !     Mark iface as bad face
             !
             lface(6,iface)=-1
             !
             !     Increment bad face counter
             !
             badface=badface+1_ip 
             !
             !     Then check if we have to add faces to the cavity
             !
             do jnode=1,3
                jneigh=eltoel(ltab(jnode,inode),ielem)
                !
                !     Is jneigh outside the cavity?
                !
                if(jneigh<0)then
                   !
                   !     Yes, this is a face of the cavity that must be deleted
                   !
                   badface=badface+1_ip
                   eltoel(ltab(jnode,inode),ielem)=lface(6,-jneigh)        
                   lface(6,-jneigh)=-1

                else
                   !
                   !     No, add this new face to the cavity
                   ! 
                   if(eltoel(1,jneigh)==ielem)then
                      iview=1_ip  
                   else if(eltoel(2,jneigh)==ielem)then
                      iview=2_ip  
                   else if(eltoel(3,jneigh)==ielem)then
                      iview=3_ip  
                   else
                      iview=4_ip  
                   endif
                   !
                   !     Resize lface and rface
                   !
                   nface=nface+1_ip 
                   call memrea(nface,memor_msh,'LFACE','gtface',lface) 
                   call memrea(nface,memor_msh,'RFACE','gtface',rface) 

                   eltoel(iview,jneigh)=-nface
                   lface(1,nface)=elem(ltab(1,iview),jneigh)
                   lface(2,nface)=elem(ltab(2,iview),jneigh)
                   lface(3,nface)=elem(ltab(3,iview),jneigh)
                   lface(4,nface)=jneigh
                   lface(5,nface)=iview
                   lface(6,nface)=ielem

                endif

             enddo

             !call outcav(ndim,npoin,coor,nface,lface)

          endif

       enddo
       !
       !     Now check if some points have been disconnected
       !
       nstacke=0_ip    
       do icav=1,ncav 
          ielem=lcav(icav)
          if(lmark(ielem)==1)then
             do inode=1,nnode 
                jpoin=elem(inode,ielem)
                !
                !     Mark lmarkp with 1
                !   
                if(lmarkp(jpoin)/=1)then
                   lmarkp(jpoin)=1 
                   nstacke=nstacke+1
                   lstacke(nstacke)=jpoin
                endif
             enddo
          endif
       enddo

       nstackf=0_ip
       do iface=1,nface
          !
          !     Is the face still valid?
          !
          if(lface(6,iface)==-1)cycle
          !
          !     Mark lmarkp with 2
          !   
          ip1=lface(1,iface)
          if(lmarkp(ip1)/=2)then
             lmarkp(ip1)=2
             nstackf=nstackf+1
             lstackf(nstackf)=ip1
          endif
          ip1=lface(2,iface)
          if(lmarkp(ip1)/=2)then
             lmarkp(ip1)=2
             nstackf=nstackf+1
             lstackf(nstackf)=ip1
          endif
          ip1=lface(3,iface)
          if(lmarkp(ip1)/=2)then
             lmarkp(ip1)=2
             nstackf=nstackf+1
             lstackf(nstackf)=ip1
          endif
       enddo
       !
       !     Is everybody here?
       !
       if(nstackf==nstacke)then
          do istacke=1,nstacke
             lmarkp(lstacke(istacke))=0_ip
          enddo
          exit
       endif
       !
       !     Will this happen?
       !
       write(*,*)'Disconnected point found in corcav'
       stop
       !
       !     Find a point not in the face point
       !
       jpoin=0_ip 
       do istacke=1,nstacke
          jpoin=lstacke(istacke)
          if(lmarkp(jpoin)==1)then
             exit
          endif
       enddo
       !
       !     Any error?
       !
       if(jpoin==0)then
          write(*,*)'Error in corcav, jpoin==0'
          stop
       endif
       !
       !     Get an element from the end of lcav (the farthest?)
       !
       do icav=ncav,1,-1 
          ielem=lcav(icav)
          if(elem(1,ielem)==jpoin)exit
          if(elem(2,ielem)==jpoin)exit
          if(elem(3,ielem)==jpoin)exit
          if(elem(4,ielem)==jpoin)exit
       enddo
       !
       !     Any error?
       !
       if(icav==0)then
          write(*,*)'Error in corcav, icav==0'
          stop
       endif
       !
       !     Add all the faces of ielem not already included in nface
       ! 
       do inode=1,nnode
          jneigh=eltoel(inode,ielem)
          !
          !     Is jneigh outside the cavity?
          !
          if(jneigh<0)then
             !
             !     Yes, this is a face of the cavity that must be deleted
             !
             badface=badface+1_ip
             eltoel(inode,ielem)=lface(6,-jneigh)        
             lface(6,-jneigh)=-1

          else
             !
             !     No, add this new face to the cavity
             ! 
             if(eltoel(1,jneigh)==ielem)then
                iview=1_ip  
             else if(eltoel(2,jneigh)==ielem)then
                iview=2_ip  
             else if(eltoel(3,jneigh)==ielem)then
                iview=3_ip  
             else
                iview=4_ip  
             endif
             !
             !     Resize lface and rface
             !
             nface=nface+1_ip 
             call memrea(nface,memor_msh,'LFACE','gtface',lface) 
             call memrea(nface,memor_msh,'RFACE','gtface',rface) 

             eltoel(iview,jneigh)=-nface
             lface(1,nface)=elem(ltab(1,iview),jneigh)
             lface(2,nface)=elem(ltab(2,iview),jneigh)
             lface(3,nface)=elem(ltab(3,iview),jneigh)
             lface(4,nface)=jneigh
             lface(5,nface)=iview
             lface(6,nface)=ielem

          endif

       enddo
       !
       !     Clean up lmarkp
       !     
       do istacke=1,nstacke
          lmarkp(lstacke(istacke))=0_ip
       enddo
       !
       !     And go back
       !
    enddo
    !
    !     Do we have still some faces?
    !
    if(badface==nface)then
       insert=.false.
    endif

  end  subroutine corcav

  subroutine newcav(nelem,lmark,eltoel,nnode,npoin,badface,lhole,nhole,nface,lenew,&
       lface,relem,rface,ipoin,elem,lptet,coor,istack,ndim,ncav,lcav,lmarkedg,lstack,&
       nenew)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: npoin,nnode,badface,ipoin,ndim,istack,ncav,nenew
    integer(ip),intent(in)    :: lcav(ncav)
    integer(ip),intent(inout) :: nhole,nelem,nface
    integer(ip),intent(inout) :: lenew(nface)
    integer(ip),pointer       :: lhole(:)
    integer(ip),intent(inout) :: lface(6,nface),lstack(nface),lmarkedg(4,nface)
    integer(ip),intent(inout) :: eltoel(nnode,*),elem(nnode,*),lmark(*)
    real(rp),intent(inout)    :: relem(5,*)
    integer(ip),intent(inout) :: lptet(npoin)
    real(rp),intent(in)       :: coor(ndim,npoin)
    real(rp),intent(inout)    :: rface(8,nface)
    integer(ip)               :: ielnew,nfac0,kpnt,iface,ielem,ip1,ip2,ip3,nstack 
    integer(ip)               :: inode,ineigh,inew,ienew,jenew 
    integer(ip)               :: kelem,idir,iview,iview2,ieold,iedg,jdir,jelem,jpnt,jnode,kneigh 
    integer(ip)               :: icav,jface,jstack 
    real(rp)                  :: c(3),c20
    real(rp)                  :: rx,ry,rz,rt,c16,rl,cx,cy,cz,rvol,radi,rnx,rny,rnz
    integer(ip)               :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))
    integer(ip)               :: ltab3(4,4)=RESHAPE((/0,2,3,4,2,0,4,3,2,3,0,4,2,4,3,0/),(/4,4/))
    integer(ip)               :: ltab4(4,4)=RESHAPE((/0,1,2,3,1,0,3,2,1,2,0,3,1,3,2,0/),(/4,4/))
    integer(ip)               :: ltab2(4)=RESHAPE((/2,3,4,1/),(/4/))
    c20=2.0d+00
    c16=1.0d+00/6.0d+00
    !
    !     Compress the faces
    ! 
    nfac0=nface
    nface=0_ip
    kpnt=0_ip 
    do iface=1,nfac0
       kpnt=kpnt+1 
       if(lface(6,iface)/=-1)then
          nface=nface+1_ip 
          lface(1,nface)=lface(1,iface) 
          lface(2,nface)=lface(2,iface) 
          lface(3,nface)=lface(3,iface) 
          ielem=lface(4,iface) 
          inode=lface(5,iface) 
          ineigh=lface(6,iface) 
          lface(4,nface)=ielem
          lface(5,nface)=inode
          lface(6,nface)=ineigh

          eltoel(inode,ielem)=-nface
          !
          !     Get center of old element
          !
          rface(1,nface)=relem(1,ielem)
          rface(2,nface)=relem(2,ielem)
          rface(3,nface)=relem(3,ielem)
          !
          !     Get volume of new element
          !
          rface(4,nface)=rface(4,iface)
          !
          !     Get radius 
          !
          rface(5,nface)=relem(4,ielem)-relem(5,ielem)           
          !
          !     Get normal 
          !
          rface(6,nface)=rface(6,iface)
          rface(7,nface)=rface(7,iface)
          rface(8,nface)=rface(8,iface)
          !
          !     Get new element number
          !
          if(nhole>0)then
             lenew(nface)=lhole(nhole)
             nhole=nhole-1_ip
          else
             nelem=nelem+1_ip
             lenew(nface)=nelem   
             if(nelem>nenew)then
                write(*,*)'Error nelem>nenew'
                stop
             endif
          endif

       endif
    enddo
    !
    !     DBG
    !
    if(nface==0)then
       write(*,*)'Error in newcav, no new cavity faces'
       write(*,*)'Check the imposed size'
       stop
    endif
    !
    !     Initialize the stack of faces
    !
    nstack=1
    lstack(1)=1_ip 
    jstack=0_ip

    do
       if(jstack==nstack)exit
       jstack=jstack+1_ip 

       iface=lstack(jstack)
       !
       !     Get topology of the face
       !  
       ip1=lface(1,iface) 
       ip2=lface(2,iface) 
       ip3=lface(3,iface) 
       ielem=lface(4,iface)
       inode=lface(5,iface)
       ineigh=lface(6,iface)
       !
       !     Get new element
       !
       ienew=lenew(iface)  
       !
       !     Get geometry of the face
       !
       cx=rface(1,iface)
       cy=rface(2,iface)
       cz=rface(3,iface)
       rvol=rface(4,iface)       
       radi=rface(5,iface)       
       rnx=rface(6,iface)
       rny=rface(7,iface)
       rnz=rface(8,iface)
       !
       !     Mark outer face
       !
       lmarkedg(4,iface)=istack
       !
       !     Create new element
       !
       elem(1,ienew)=ipoin
       elem(2,ienew)=lface(1,iface)
       elem(3,ienew)=lface(2,iface)
       elem(4,ienew)=lface(3,iface)
       !
       !     Get new center
       !
       rt=radi/(c20*rvol) 
       relem(1,ienew)=cx-rt*rnx
       relem(2,ienew)=cy-rt*rny
       relem(3,ienew)=cz-rt*rnz
       !
       !     DBG
       !
       !do inode=1,nnode
       !   rx=coor(1,elem(inode,ienew))-relem(1,ienew)
       !   ry=coor(2,elem(inode,ienew))-relem(2,ienew)
       !   rz=coor(3,elem(inode,ienew))-relem(3,ienew)
       !   rl=sqrt(rx*rx+ry*ry+rz*rz)
       !enddo   
       !
       !     Get new radius
       !
       rx=coor(1,ipoin)-relem(1,ienew)
       ry=coor(2,ipoin)-relem(2,ienew)
       rz=coor(3,ipoin)-relem(3,ienew)
       rl=rx*rx+ry*ry+rz*rz

       relem(5,ienew)=rl
       lmark(ienew)=0_ip
       !
       !     Do we really need the volume of the element?
       !
       !relem(6,ienew)=rvol*c16 
       !
       !     Get external neighbor
       !
       if(ineigh==0)then
          eltoel(1,ienew)=0_ip
       else

          if(eltoel(1,ineigh)==ielem)then
             iview=1_ip
          else if(eltoel(2,ineigh)==ielem)then
             iview=2_ip
          else if(eltoel(3,ineigh)==ielem)then
             iview=3_ip
          else
             iview=4_ip
          endif
          eltoel(iview,ineigh)=ienew
          eltoel(1,ienew)=ineigh
       endif
       !
       !     Update point to element pointer
       !
       lptet(ipoin)=ienew
       lptet(elem(2,ienew))=ienew  
       lptet(elem(3,ienew))=ienew  
       lptet(elem(4,ienew))=ienew  
       !
       !     Get internal neighbors
       !
       do iedg=1,3
          !
          !     Has this side been taken into account?
          ! 
          if(lmarkedg(iedg,iface)/=istack)then 

             lmarkedg(iedg,iface)=istack 
             !
             !     Get direction in old element
             !
             jdir=ltab(iedg,inode)
             jelem=ielem
             jpnt=elem(inode,ielem)
             jnode=inode
             !
             !     Roll over the edges inside the cavity
             !
             do

                kneigh=eltoel(jdir,jelem)
                !
                !     Did we reach a cavity face?
                !
                if(kneigh<0)then
                   !
                   !     Get the face number in the stack
                   !     
                   jface=-kneigh
                   !
                   !     Mark the edge found in jface
                   !
                   lmarkedg(ltab4(jnode,jdir),jface)=istack  
                   !
                   !     Update eltoel between ienew and lenew(jface)
                   !
                   jenew=lenew(jface) 
                   eltoel(iedg+1,ienew)=jenew 
                   eltoel(ltab3(jnode,jdir),jenew)=ienew
                   !
                   !     Do we have to add jface to the stack?
                   !
                   if(lmarkedg(4,jface)/=istack)then
                      lmarkedg(4,jface)=istack
                      nstack=nstack+1_ip   
                      lstack(nstack)=jface
                   endif
                   !
                   !     And go home
                   !
                   exit

                else 

                   if(elem(1,kneigh)==jpnt)then
                      iview=1_ip  
                   else if(elem(2,kneigh)==jpnt)then
                      iview=2_ip  
                   else if(elem(3,kneigh)==jpnt)then
                      iview=3_ip  
                   else
                      iview=4_ip  
                   endif

                   if(eltoel(1,kneigh)==jelem)then
                      iview2=1_ip  
                   else if(eltoel(2,kneigh)==jelem)then
                      iview2=2_ip  
                   else if(eltoel(3,kneigh)==jelem)then
                      iview2=3_ip  
                   else
                      iview2=4_ip  
                   endif

                   jpnt=elem(iview2,kneigh)
                   jnode=iview2          
                   jelem=kneigh
                   jdir=iview 

                endif

             enddo

          endif

       enddo

    enddo
    !
    !     DBG
    !
    if(nstack/=nface)then
       write(*,*)'Error in newcav'
       write(*,*)'nstack=',nstack,'nface=',nface
       stop
    endif
    !
    !     Get old elements of the cavity
    !
    do icav=1,ncav
       ielem=lcav(icav)
       if(lmark(ielem)==1)then
          elem(1,ielem)=0_ip
          nhole=nhole+1
          call memrea(nhole,memor_msh,'LHOLE','newcav',lhole) 
          lhole(nhole)=ielem 
          lmark(ielem)=0_ip
       endif
    enddo

  end subroutine newcav

  subroutine findelem(iguess,elem,nelem,npoin,coor,pnew,ndim,nnode,eltoel,ncav,mcav,lcav)
    use def_kintyp, only       :  ip,rp,lg
    use mod_mshtol, only       : orient3D
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode,iguess,mcav
    integer(ip),intent(in)     :: elem(nnode,nelem),eltoel(nnode,nelem)
    real(rp),intent(in)        :: coor(ndim,npoin),pnew(ndim)
    integer(ip),intent(inout)  :: ncav
    integer(ip),intent(inout)  :: lcav(mcav)
    integer(ip)                :: iter,maxiter,ip1,ip2,ip3,ip4,l,ihost
    real(rp)                   :: c40,p1(ndim),p2(ndim),p3(ndim)
    real(rp)                   :: dtot,rtol,rl,d1,d2,d3,d4
    !
    !     This subroutine finds the host element by a random walk through
    !     the mesh
    ! 
    c40=4.0d+00  
    maxiter=500_ip
    rtol=1.0d-06

    if(iguess==0)then
       write(*,*)'Error in iguess in findelem'  
       stop
    endif
    if(iguess<0)then
       write(*,*)'Error in findin, iguess=',iguess
       stop
    endif

    ihost=iguess

    do iter=1,maxiter 

       if(ihost==0)then
          write(*,*)'Error in findelem, ihost=0' 
          stop
       endif

       if(ihost<0 .or. ihost>nelem)then
          write(*,*)'Error findelem, ihost=',ihost
          stop
       endif

       call random_number(rl)
       l=floor(c40*rl)+1_ip


       ip1=elem(1,ihost)
       ip2=elem(2,ihost)
       ip3=elem(3,ihost)
       ip4=elem(4,ihost)
       !
       !     Compute volume for tolerances
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)
       p3(1)=coor(1,ip4)-coor(1,ip1)
       p3(2)=coor(2,ip4)-coor(2,ip1)
       p3(3)=coor(3,ip4)-coor(3,ip1)

       call orient3D(p1,p2,p3,dtot,ndim)
       !
       !     First face
       !
       if(l==1)then

          p1(1)=coor(1,ip2)-pnew(1)
          p1(2)=coor(2,ip2)-pnew(2)
          p1(3)=coor(3,ip2)-pnew(3)
          p2(1)=coor(1,ip3)-pnew(1)
          p2(2)=coor(2,ip3)-pnew(2)
          p2(3)=coor(3,ip3)-pnew(3)
          p3(1)=coor(1,ip4)-pnew(1)
          p3(2)=coor(2,ip4)-pnew(2)
          p3(3)=coor(3,ip4)-pnew(3)

          call orient3D(p1,p2,p3,d1,ndim)
          d1=d1/dtot 

          if(d1>rtol)then

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip4)-pnew(1)
             p2(2)=coor(2,ip4)-pnew(2)
             p2(3)=coor(3,ip4)-pnew(3)
             p3(1)=coor(1,ip3)-pnew(1)
             p3(2)=coor(2,ip3)-pnew(2)
             p3(3)=coor(3,ip3)-pnew(3)

             call orient3D(p1,p2,p3,d2,ndim)
             d2=d2/dtot 

             if(d2>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip2)-pnew(1)
                p2(2)=coor(2,ip2)-pnew(2)
                p2(3)=coor(3,ip2)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      ncav=1
                      lcav(ncav)=ihost
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(4,ihost)
                      return
                   endif


                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(3,ihost)
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,3_ip,4_ip,nelem,nnode)
                      return
                   endif

                endif

             else if(d2<-rtol)then

                ihost=eltoel(2,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip2)-pnew(1)
                p2(2)=coor(2,ip2)-pnew(2)
                p2(3)=coor(3,ip2)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(2,ihost) 
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,2_ip,4_ip,nelem,nnode)
                      return
                   endif

                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,2_ip,3_ip,nelem,nnode)
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 1'
                      stop
                      return
                   endif

                endif

             endif

          else if(d1<-rtol)then

             ihost=eltoel(1,ihost)

          else      

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip4)-pnew(1)
             p2(2)=coor(2,ip4)-pnew(2)
             p2(3)=coor(3,ip4)-pnew(3)
             p3(1)=coor(1,ip3)-pnew(1)
             p3(2)=coor(2,ip3)-pnew(2)
             p3(3)=coor(3,ip3)-pnew(3)

             call orient3D(p1,p2,p3,d2,ndim)
             d2=d2/dtot 

             if(d2>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip2)-pnew(1)
                p2(2)=coor(2,ip2)-pnew(2)
                p2(3)=coor(3,ip2)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(1,ihost)
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,1_ip,4_ip,nelem,nnode)
                      return
                   endif


                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,1_ip,3_ip,nelem,nnode)
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 2'
                      stop
                      return
                   endif

                endif

             else if(d2<-rtol)then

                ihost=eltoel(2,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip2)-pnew(1)
                p2(2)=coor(2,ip2)-pnew(2)
                p2(3)=coor(3,ip2)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,1_ip,2_ip,nelem,nnode)
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 3'
                      stop
                      return
                   endif

                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      write(*,*)'Point too close from previously inserted point 4'
                      stop
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else

                      write(*,*)'Error in findin, degenerate case'
                      stop
                   endif

                endif

             endif

          endif
          !
          !     Second face
          !
       else if(l==2)then

          p1(1)=coor(1,ip1)-pnew(1)
          p1(2)=coor(2,ip1)-pnew(2)
          p1(3)=coor(3,ip1)-pnew(3)
          p2(1)=coor(1,ip4)-pnew(1)
          p2(2)=coor(2,ip4)-pnew(2)
          p2(3)=coor(3,ip4)-pnew(3)
          p3(1)=coor(1,ip3)-pnew(1)
          p3(2)=coor(2,ip3)-pnew(2)
          p3(3)=coor(3,ip3)-pnew(3)

          call orient3D(p1,p2,p3,d2,ndim)
          d2=d2/dtot 

          if(d2>rtol)then

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip2)-pnew(1)
             p2(2)=coor(2,ip2)-pnew(2)
             p2(3)=coor(3,ip2)-pnew(3)
             p3(1)=coor(1,ip4)-pnew(1)
             p3(2)=coor(2,ip4)-pnew(2)
             p3(3)=coor(3,ip4)-pnew(3)

             call orient3D(p1,p2,p3,d3,ndim)
             d3=d3/dtot 

             if(d3>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip2)-pnew(1)
                p3(2)=coor(2,ip2)-pnew(2)
                p3(3)=coor(3,ip2)-pnew(3)

                call orient3D(p1,p2,p3,d4,ndim)
                d4=d4/dtot 

                if(d4>rtol)then

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      ncav=1
                      lcav(ncav)=ihost
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(1,ihost)
                      return
                   endif


                else if(d4<-rtol)then

                   ihost=eltoel(4,ihost)

                else

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(4,ihost)
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,4_ip,1_ip,nelem,nnode)
                      return
                   endif

                endif

             else if(d3<-rtol)then

                ihost=eltoel(3,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip2)-pnew(1)
                p3(2)=coor(2,ip2)-pnew(2)
                p3(3)=coor(3,ip2)-pnew(3)

                call orient3D(p1,p2,p3,d4,ndim)
                d4=d4/dtot 

                if(d4>rtol)then

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(3,ihost)
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,3_ip,1_ip,nelem,nnode)  
                      return
                   endif

                else if(d4<-rtol)then

                   ihost=eltoel(4,ihost)

                else

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,3_ip,4_ip,nelem,nnode)
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 5'
                      stop
                      return
                   endif

                endif

             endif

          else if(d2<-rtol)then

             ihost=eltoel(2,ihost)

          else      

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip2)-pnew(1)
             p2(2)=coor(2,ip2)-pnew(2)
             p2(3)=coor(3,ip2)-pnew(3)
             p3(1)=coor(1,ip4)-pnew(1)
             p3(2)=coor(2,ip4)-pnew(2)
             p3(3)=coor(3,ip4)-pnew(3)

             call orient3D(p1,p2,p3,d3,ndim)
             d3=d3/dtot 

             if(d3>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip2)-pnew(1)
                p3(2)=coor(2,ip2)-pnew(2)
                p3(3)=coor(3,ip2)-pnew(3)

                call orient3D(p1,p2,p3,d4,ndim)
                d4=d4/dtot 

                if(d4>rtol)then

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(2,ihost)
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,2_ip,1_ip,nelem,nnode)
                      return
                   endif


                else if(d4<-rtol)then

                   ihost=eltoel(4,ihost)

                else

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,2_ip,4_ip,nelem,nnode)
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 6'
                      stop
                      return
                   endif

                endif

             else if(d3<-rtol)then

                ihost=eltoel(3,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip2)-pnew(1)
                p3(2)=coor(2,ip2)-pnew(2)
                p3(3)=coor(3,ip2)-pnew(3)

                call orient3D(p1,p2,p3,d4,ndim)
                d4=d4/dtot 

                if(d4>rtol)then

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,2_ip,3_ip,nelem,nnode) 
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 7'
                      stop
                      return
                   endif

                else if(d4<-rtol)then

                   ihost=eltoel(4,ihost)

                else

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      write(*,*)'Point too close from previously inserted point 8'
                      stop
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else

                      write(*,*)'Error in findin, degenerate case'
                      stop
                   endif

                endif

             endif

          endif
          !
          !     Third face
          !
       else if(l==3)then

          p1(1)=coor(1,ip1)-pnew(1)
          p1(2)=coor(2,ip1)-pnew(2)
          p1(3)=coor(3,ip1)-pnew(3)
          p2(1)=coor(1,ip2)-pnew(1)
          p2(2)=coor(2,ip2)-pnew(2)
          p2(3)=coor(3,ip2)-pnew(3)
          p3(1)=coor(1,ip4)-pnew(1)
          p3(2)=coor(2,ip4)-pnew(2)
          p3(3)=coor(3,ip4)-pnew(3)

          call orient3D(p1,p2,p3,d3,ndim)
          d3=d3/dtot 

          if(d3>rtol)then

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip3)-pnew(1)
             p2(2)=coor(2,ip3)-pnew(2)
             p2(3)=coor(3,ip3)-pnew(3)
             p3(1)=coor(1,ip2)-pnew(1)
             p3(2)=coor(2,ip2)-pnew(2)
             p3(3)=coor(3,ip2)-pnew(3)

             call orient3D(p1,p2,p3,d4,ndim)
             d4=d4/dtot 

             if(d4>rtol)then

                p1(1)=coor(1,ip2)-pnew(1)
                p1(2)=coor(2,ip2)-pnew(2)
                p1(3)=coor(3,ip2)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      ncav=1
                      lcav(ncav)=ihost
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(2,ihost)
                      return
                   endif


                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(1,ihost)
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,1_ip,2_ip,nelem,nnode)
                      return
                   endif

                endif

             else if(d4<-rtol)then

                ihost=eltoel(4,ihost)

             else

                p1(1)=coor(1,ip2)-pnew(1)
                p1(2)=coor(2,ip2)-pnew(2)
                p1(3)=coor(3,ip2)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(4,ihost)
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,4_ip,2_ip,nelem,nnode) 
                      return
                   endif

                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,4_ip,1_ip,nelem,nnode)
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 9'
                      stop
                      return
                   endif

                endif

             endif

          else if(d3<-rtol)then

             ihost=eltoel(3,ihost)

          else      

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip3)-pnew(1)
             p2(2)=coor(2,ip3)-pnew(2)
             p2(3)=coor(3,ip3)-pnew(3)
             p3(1)=coor(1,ip2)-pnew(1)
             p3(2)=coor(2,ip2)-pnew(2)
             p3(3)=coor(3,ip2)-pnew(3)

             call orient3D(p1,p2,p3,d4,ndim)
             d4=d4/dtot 

             if(d4>rtol)then

                p1(1)=coor(1,ip2)-pnew(1)
                p1(2)=coor(2,ip2)-pnew(2)
                p1(3)=coor(3,ip2)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      ncav=2
                      lcav(1)=ihost 
                      lcav(2)=eltoel(3,ihost) 
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,3_ip,2_ip,nelem,nnode)
                      return
                   endif


                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,3_ip,1_ip,nelem,nnode)
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 10'
                      stop
                      return
                   endif

                endif

             else if(d4<-rtol)then

                ihost=eltoel(4,ihost)

             else

                p1(1)=coor(1,ip2)-pnew(1)
                p1(2)=coor(2,ip2)-pnew(2)
                p1(3)=coor(3,ip2)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,3_ip,4_ip,nelem,nnode)
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 11'
                      stop
                      return
                   endif

                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      write(*,*)'Point too close from previously inserted point 12'
                      stop
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else

                      write(*,*)'Error in findin, degenerate case'
                      stop
                   endif

                endif

             endif

          endif
          !
          !     Fourth face
          !
       else

          p1(1)=coor(1,ip1)-pnew(1)
          p1(2)=coor(2,ip1)-pnew(2)
          p1(3)=coor(3,ip1)-pnew(3)
          p2(1)=coor(1,ip3)-pnew(1)
          p2(2)=coor(2,ip3)-pnew(2)
          p2(3)=coor(3,ip3)-pnew(3)
          p3(1)=coor(1,ip2)-pnew(1)
          p3(2)=coor(2,ip2)-pnew(2)
          p3(3)=coor(3,ip2)-pnew(3)

          call orient3D(p1,p2,p3,d4,ndim)
          d4=d4/dtot 

          if(d4>rtol)then

             p1(1)=coor(1,ip2)-pnew(1)
             p1(2)=coor(2,ip2)-pnew(2)
             p1(3)=coor(3,ip2)-pnew(3)
             p2(1)=coor(1,ip3)-pnew(1)
             p2(2)=coor(2,ip3)-pnew(2)
             p2(3)=coor(3,ip3)-pnew(3)
             p3(1)=coor(1,ip4)-pnew(1)
             p3(2)=coor(2,ip4)-pnew(2)
             p3(3)=coor(3,ip4)-pnew(3)

             call orient3D(p1,p2,p3,d1,ndim)
             d1=d1/dtot 

             if(d1>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip4)-pnew(1)
                p2(2)=coor(2,ip4)-pnew(2)
                p2(3)=coor(3,ip4)-pnew(3)
                p3(1)=coor(1,ip3)-pnew(1)
                p3(2)=coor(2,ip3)-pnew(2)
                p3(3)=coor(3,ip3)-pnew(3)

                call orient3D(p1,p2,p3,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      ncav=1
                      lcav(ncav)=ihost
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(3,ihost)
                      return
                   endif


                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(2,ihost)
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,2_ip,3_ip,nelem,nnode)
                      return
                   endif

                endif

             else if(d1<-rtol)then

                ihost=eltoel(1,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip4)-pnew(1)
                p2(2)=coor(2,ip4)-pnew(2)
                p2(3)=coor(3,ip4)-pnew(3)
                p3(1)=coor(1,ip3)-pnew(1)
                p3(2)=coor(2,ip3)-pnew(2)
                p3(3)=coor(3,ip3)-pnew(3)

                call orient3D(p1,p2,p3,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(1,ihost)
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,1_ip,3_ip,nelem,nnode)
                      return
                   endif

                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,1_ip,2_ip,nelem,nnode)
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 13'
                      stop
                      return
                   endif

                endif

             endif

          else if(d4<-rtol)then

             ihost=eltoel(4,ihost)

          else      

             p1(1)=coor(1,ip2)-pnew(1)
             p1(2)=coor(2,ip2)-pnew(2)
             p1(3)=coor(3,ip2)-pnew(3)
             p2(1)=coor(1,ip3)-pnew(1)
             p2(2)=coor(2,ip3)-pnew(2)
             p2(3)=coor(3,ip3)-pnew(3)
             p3(1)=coor(1,ip4)-pnew(1)
             p3(2)=coor(2,ip4)-pnew(2)
             p3(3)=coor(3,ip4)-pnew(3)

             call orient3D(p1,p2,p3,d1,ndim)
             d1=d1/dtot 

             if(d1>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip4)-pnew(1)
                p2(2)=coor(2,ip4)-pnew(2)
                p2(3)=coor(3,ip4)-pnew(3)
                p3(1)=coor(1,ip3)-pnew(1)
                p3(2)=coor(2,ip3)-pnew(2)
                p3(3)=coor(3,ip3)-pnew(3)

                call orient3D(p1,p2,p3,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      ncav=2
                      lcav(1)=ihost
                      lcav(2)=eltoel(4,ihost)
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,4_ip,3_ip,nelem,nnode)
                      return
                   endif


                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,4_ip,2_ip,nelem,nnode)
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 14'
                      stop
                      return
                   endif

                endif

             else if(d1<-rtol)then

                ihost=eltoel(1,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip4)-pnew(1)
                p2(2)=coor(2,ip4)-pnew(2)
                p2(3)=coor(3,ip4)-pnew(3)
                p3(1)=coor(1,ip3)-pnew(1)
                p3(2)=coor(2,ip3)-pnew(2)
                p3(3)=coor(3,ip3)-pnew(3)

                call orient3D(p1,p2,p3,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      call gtshell(lcav,ncav,mcav,eltoel,elem,ihost,4_ip,1_ip,nelem,nnode)
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 15'
                      stop
                      return
                   endif

                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      write(*,*)'Point too close from previously inserted point 16'
                      stop
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      write(*,*)'Error in findin, degenerate case'
                      stop
                   endif

                endif

             endif

          endif

       endif

    enddo

  end subroutine findelem

  subroutine gtshell(lcav,ncav,mcav,eltoel,elem,ielem,idir,inode,nelem,nnode)
    use def_kintyp, only       :  ip,rp,lg
    use mod_mshtol, only       : orient3D
    implicit none
    integer(ip),intent(in)     :: nelem,nnode,ielem,idir,inode,mcav
    integer(ip),intent(in)     :: eltoel(nnode,nelem),elem(nnode,nelem)
    integer(ip),intent(inout)  :: ncav,lcav(mcav)
    integer(ip)                :: jelem,jnode,jpnt,iview,iview2,kneigh,kelem,jdir

    jelem=ielem
    jdir=idir 
    jpnt=elem(inode,ielem)
    jnode=inode
    ncav=1_ip
    lcav(1)=ielem
    !
    !     Roll over the edges inside the cavity
    !
    do

       kneigh=eltoel(jdir,jelem)

       if(kneigh==ielem)exit

       ncav=ncav+1_ip
       lcav(ncav)=kneigh
       !
       !     Did we reach a cavity face
       !
       if(kneigh<=0)then
          write(*,*)'Error in gtshell'
          stop 
       endif

       if(elem(1,kneigh)==jpnt)then
          iview=1_ip  
       else if(elem(2,kneigh)==jpnt)then
          iview=2_ip  
       else if(elem(3,kneigh)==jpnt)then
          iview=3_ip  
       else
          iview=4_ip  
       endif

       if(eltoel(1,kneigh)==jelem)then
          iview2=1_ip  
       else if(eltoel(2,kneigh)==jelem)then
          iview2=2_ip  
       else if(eltoel(3,kneigh)==jelem)then
          iview2=3_ip  
       else
          iview2=4_ip  
       endif

       jpnt=elem(iview2,kneigh)
       jnode=iview2          
       jelem=kneigh
       jdir=iview

    enddo

  end subroutine gtshell

  subroutine newpnt(coor,elem,ndim,npoin,nelem,nnode,lbin1,lbin2,nbin,rsize,dxbin,dybin,dzbin,&
       nxbin,nybin,nzbin,bboxbin,neold,elold,coold,eltoelold,pt1old,pt2old,lcart,lcell,ncell,&
       lctop1,lctop2,nitermax,csearch,npold,lcmark,lcstack,rsold,lptetold,lptet,lmarkold,rtol,&
       ierr,lmarkp)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    use mod_bin3d , only : filter
    use mod_mshtol, only  : ptoelmb
    use mod_interp, only  : interpnt
    use mod_cart, only  : gtelem
    implicit none
    integer(ip),intent(in)       :: ndim,nnode,nelem,nxbin,nybin,nzbin,nbin
    integer(ip),intent(in)       :: neold,npold,ncell,nitermax
    integer(ip),intent(inout)    :: npoin,ierr
    integer(ip),intent(in)       :: elem(nnode,nelem)
    integer(ip),intent(in)       :: elold(nnode,neold),eltoelold(nnode,neold)
    integer(ip),intent(in)       :: pt1old(*),pt2old(npold+1)
    integer(ip),intent(in)       :: lctop1(*),lctop2(ncell+1)
    integer(ip),intent(inout)    :: lcmark(ncell),lcstack(ncell),lmarkold(neold)
    type(cell)                   :: lcell(ncell)
    real(rp),intent(in)          :: coold(ndim,npold),csearch,rsold(npold),rtol
    real(rp),intent(in)          :: dxbin,dybin,dzbin,bboxbin(3,2)
    real(rp),pointer             :: rsize(:)
    real(rp),pointer             :: coor(:,:)
    integer(ip),pointer          :: lbin1(:),lcart(:),lptetold(:),lptet(:),lmarkp(:) 
    integer(ip),intent(inout)    :: lbin2(nbin+1) 
    integer(ip),pointer          :: ptoel1(:),ptoel2(:),lmark(:),ledge(:,:),lptoldn(:),lcartn(:)  
    integer(ip)                  :: nedge,ipoin,isto,ielem,inode,jpoin,npnew1,npgood
    real(rp)                     :: rx,ry,rz,length,l1,pmid(3),delta,c05,c10,d1,d2,d3,d4
    real(rp),pointer             :: coorn(:,:),rsizn(:)
    integer(ip)                  :: ip1,ip2,iedge,npnew,ipnew,ip3,ip4,ielold,iguess,icart
    integer(4)                   :: istat
    integer(ip),parameter     :: mstackp=500
    integer(ip)               :: lstackp(mstackp)
    !
    !     Fix tolerance
    !
    delta=sqrt(2.0d+00)
    c05=0.5d+00
    c10=0.5d+00
    nullify(ptoel1,ptoel2)
    !
    !     This subroutine creates new points along edges
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','newpnt',lmark)
    allocate(coorn(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'COORN','newpnt',coorn)
    allocate(rsizn(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'RSIZN','newpnt',rsizn)
    allocate(lptoldn(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTOLDN','newpnt',lptoldn)
    allocate(lcartn(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LCARTN','newpnt',lcartn)
    ! 
    !     First get the elements surrouding elements
    !
    call ptoelmb(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !     Then the edges
    !
    nedge=0_ip 
    do ipoin=9,npoin
       do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(isto)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(jpoin>ipoin)then
                if(lmark(jpoin)/=ipoin)then
                   lmark(jpoin)=ipoin
                   nedge=nedge+1_ip
                endif
             endif
          enddo
       enddo
    enddo

    allocate(ledge(2,nedge),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDGE','newpnt',ledge)

    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo

    nedge=0_ip 
    do ipoin=9,npoin
       lmark(ipoin)=ipoin 
       do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(isto)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(jpoin>ipoin)then
                if(lmark(jpoin)/=ipoin)then
                   lmark(jpoin)=ipoin
                   nedge=nedge+1_ip
                   ledge(1,nedge)=ipoin
                   ledge(2,nedge)=jpoin
                endif
             endif
          enddo
       enddo
    enddo
    !
    !     Then create points
    !
    npnew=0_ip
    do iedge=1,nedge
       ip1=ledge(1,iedge)
       ip2=ledge(2,iedge)

       rx=coor(1,ip1)-coor(1,ip2)
       ry=coor(2,ip1)-coor(2,ip2)
       rz=coor(3,ip1)-coor(3,ip2)
       length=sqrt(rx*rx+ry*ry+rz*rz)
       l1=length*c05*(c10/rsize(ip1)+c10/rsize(ip2)) 

       if(l1>delta)then
          !
          !     Get mid point
          !
          pmid(1)=c05*(coor(1,ip1)+coor(1,ip2))
          pmid(2)=c05*(coor(2,ip1)+coor(2,ip2))
          pmid(3)=c05*(coor(3,ip1)+coor(3,ip2))
          !
          !     Add it to list of new points
          ! 
          npnew=npnew+1_ip
          call memrea(npnew,memor_msh,'COORN','newpnt',coorn) 
          call memrea(npnew,memor_msh,'RSIZN','newpnt',rsizn) 
          call memrea(npnew,memor_msh,'LPTOLDN','newpnt',lptoldn) 
          call memrea(npnew,memor_msh,'LCARTN','newpnt',lcartn) 
          coorn(1,npnew)=pmid(1)
          coorn(2,npnew)=pmid(2)
          coorn(3,npnew)=pmid(3)
          !
          !     Get cartesian cell
          !
          icart=lcart(ip1)
          call gtelem(npnew,coorn,npnew,ndim,lcell,ncell,icart,rtol)
          lcartn(npnew)=icart
          !
          !     Get the size of the mid point
          !
          !write(*,*)iedge
          call interpnt(nnode,nelem,elem,ndim,neold,elold,coorn,coold,npnew,npold,&
               npnew,eltoelold,pt1old,pt2old,lcartn,lcell,ncell,lcmark,&
               lcstack,lctop1,lctop2,nitermax,lstackp,mstackp,csearch,lmarkold,&
               ielold,d1,d2,d3,d4,lptetold(ip1),ierr)
          if(ierr==1)then
             write(*,*)'Error in newpnt, npnew=',npnew
             write(*,*)'x y z:',coor(1,npnew),coor(2,npnew),coor(3,npnew)
             return
          endif
          !
          !     Interpolate size
          !
          ip1=elold(1,ielold)
          ip2=elold(2,ielold)
          ip3=elold(3,ielold)
          ip4=elold(4,ielold)
          rsizn(npnew)=d1*rsold(ip1)+d2*rsold(ip2)+&
               d3*rsold(ip3)+d4*rsold(ip4)
          lptoldn(npnew)=ielold

       endif

    enddo
    !
    !     Clean up lmark
    !
    call memrea(npnew,memor_msh,'LMARK','newpnt',lmark) 
    do ipoin=1,npnew
       lmark(ipoin)=0_ip
    enddo
    !
    !     Filter the new points with the bin
    ! 
    call filter(lbin1,lbin2,nbin,npnew,ndim,coorn,dxbin,dybin,dzbin,nxbin,nybin,nzbin,rsize,&
         coor,npoin,bboxbin,lmark,rsizn)
    !
    !     Obtain the new points
    !
    npgood=0_ip
    do ipnew=1,npnew
       if(lmark(ipnew)==0)then
          npgood=npgood+1_ip
       endif
    enddo
    !
    !     Resize nodal database
    !
    npnew1=npoin+npgood
    call memrea(npnew1,memor_msh,'COOR','newpnt',coor)
    call memrea(npnew1,memor_msh,'RSIZE','newpnt',rsize)
    call memrea(npnew1,memor_msh,'LPTET','newpnt',lptet)
    call memrea(npnew1,memor_msh,'LPTETOLD','newpnt',lptetold)
    call memrea(npnew1,memor_msh,'LCART','newpnt',lcart)
    call memrea(npnew1,memor_msh,'LMARKP','newpnt',lmarkp)
    !
    !     And transfer
    !
    do ipnew=1,npnew
       if(lmark(ipnew)==0)then
          npoin=npoin+1
          coor(1,npoin)=coorn(1,ipnew)
          coor(2,npoin)=coorn(2,ipnew)
          coor(3,npoin)=coorn(3,ipnew)
          rsize(npoin)=rsizn(ipnew)
          lptetold(npoin)=lptoldn(ipnew)
          lcart(npoin)=lcartn(ipnew)
       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'LCARTN','newpnt',lcartn)
    deallocate(lcartn,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCARTN','newpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTOLDN','newpnt',lptoldn)
    deallocate(lptoldn,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTOLDN','newpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'RSIZN','newpnt',rsizn)
    deallocate(rsizn,stat=istat)
    if(istat/=0) call memerr(2_ip,'RSIZN','newpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'COORN','newpnt',coorn)
    deallocate(coorn,stat=istat)
    if(istat/=0) call memerr(2_ip,'COORN','newpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','newpnt',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','newpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','newpnt',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','newpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','newpnt',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','newpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGE','newpnt',ledge)
    deallocate(ledge,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGE','newpnt',0_ip)

  end subroutine newpnt

  subroutine gtsize(pin)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    implicit none
    real(rp),intent(in)  :: pin(3)
  end subroutine gtsize

  subroutine outdelau(nnode,nelem,elem,ndim,npoin,coor,rsize)
    use def_kintyp, only       : ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nnode,nelem,npoin,ndim
    integer(ip),intent(in)      :: elem(nnode,nelem)
    real(rp), intent(in)        :: coor(ndim,npoin),rsize(npoin)
    real(rp)                    :: rx,ry,rz
    integer(ip)                 :: ipoin,ielem,ncont 
    !
    !     Print mesh
    !
    open(unit=70,file='outdelau.msh',status='unknown')
    rewind 70

1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coord_x   coord_y  coord_z')
100 format(i10,3e20.10e3)
200 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(70,1)
    write(70,2)
    write(70,3)
    do  ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(70,100)ipoin,rx,ry,rz
    enddo

    write(70,4)
    write(70,5)
    ncont=0_ip
    do ielem=1,nelem
       if(elem(1,ielem)/=0)then
          ncont=ncont+1_ip
          write(70,200)ncont,elem(1,ielem),elem(2,ielem),elem(3,ielem),elem(4,ielem)
       endif
    enddo
    write(70,6)
    close(70)

    open(unit=60,file='outdelau.res',status='unknown')
    rewind 60
    write(60,10)
    write(60,16)
    write(60,13)
    do  ipoin=1,npoin
       write(60,300)ipoin,rsize(ipoin)
    enddo
    write(60,14)
    write(60,15)
    close(60)

10  format('GID Post Results File 1.0')
13  format('Values')
14  format('End Values')
15  format('   ')
16  format('Result "Size" "Analysis/time" ','0.0',' Scalar OnNodes')
300 format(i10,1e20.10)

  end subroutine outdelau

  subroutine outelem(nnode,nelem,elem,ndim,npoin,coor,lcav,ncav)
    use def_kintyp, only       : ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nnode,nelem,npoin,ndim,ncav
    integer(ip),intent(in)      :: elem(nnode,nelem),lcav(ncav)
    real(rp), intent(in)        :: coor(ndim,npoin)
    real(rp)                    :: rx,ry,rz
    integer(ip)                 :: ipoin,ielem,ncont,icav 
    !
    !     Print mesh
    !
    open(unit=70,file='outelem.msh',status='unknown')
    rewind 70

1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coord_x   coord_y  coord_z')
100 format(i10,3e20.10e3)
200 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(70,1)
    write(70,2)
    write(70,3)
    do  ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(70,100)ipoin,rx,ry,rz
    enddo

    write(70,4)
    write(70,5)
    ncont=0_ip
    do icav=1,ncav
       ielem=lcav(icav)
       ncont=ncont+1_ip
       write(70,200)ncont,elem(1,ielem),elem(2,ielem),elem(3,ielem),elem(4,ielem)
    enddo
    write(70,6)
    close(70)

  end subroutine outelem

  subroutine outcav(ndim,npoin,coor,nface,lface)
    use def_kintyp, only       : ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nface,npoin,ndim
    integer(ip),intent(in)      :: lface(6,nface)
    real(rp), intent(in)        :: coor(ndim,npoin)
    real(rp)                    :: rx,ry,rz
    integer(ip)                 :: ipoin,iface,ncont
    !
    !     Print mesh
    !
    open(unit=70,file='outcav.msh',status='unknown')
    rewind 70

1   format('MESH dimension 3 ElemType Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coord_x   coord_y  coord_z')
100 format(i10,3e20.10e3)
200 format(4i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(70,1)
    write(70,2)
    write(70,3)
    do  ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(70,100)ipoin,rx,ry,rz
    enddo

    write(70,4)
    write(70,5)
    ncont=0_ip
    do iface=1,nface
       if(lface(6,iface)/=-1)then
          ncont=ncont+1
          write(70,200)ncont,lface(1,iface),lface(2,iface),lface(3,iface)
       endif
    enddo
    write(70,6)
    close(70)

  end subroutine outcav

  subroutine chkconf(eltoel,nnode,elem,nelem,ndim,npoin)
    use def_kintyp, only       : ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    use mod_mshtol, only: ptoelmb
    implicit none
    integer(ip),intent(in)      :: ndim,nnode,nelem,npoin
    integer(ip),intent(in)      :: elem(nnode,nelem),eltoel(nnode,nelem)
    integer(ip),pointer         :: lmark(:),ptoel1(:),ptoel2(:)
    integer(ip)                 :: ielem,ineigh,inode,iview,ip1,ip2,ip3,ip4,ipa,ipb,ipc
    integer(ip)                 :: ipoin
    integer(4)                :: istat
    integer(ip)               :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))

    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','chkconf',lmark)

    do ielem=1,nelem
       if(elem(1,ielem)==0)cycle

       ip1=elem(1,ielem)
       if(ip1>npoin .or. ip1<=0)then
          write(*,*)'Error in chkconf,ielem=',ielem
          write(*,*)'ip1=',ip1
          stop
       endif
       ip2=elem(2,ielem)
       if(ip2>npoin .or. ip2<=0)then
          write(*,*)'Error in chkconf,ielem=',ielem
          write(*,*)'ip2=',ip2
          stop
       endif
       ip3=elem(3,ielem)
       if(ip3>npoin .or. ip3<=0)then
          write(*,*)'Error in chkconf,ielem=',ielem
          write(*,*)'ip3=',ip3
          stop
       endif
       ip4=elem(4,ielem)
       if(ip4>npoin .or. ip4<=0)then
          write(*,*)'Error in chkconf,ielem=',ielem
          write(*,*)'ip4=',ip4
          stop
       endif
       do inode=1,nnode
          ineigh=eltoel(inode,ielem)
          if(ineigh==0)cycle
          if(elem(1,ineigh)==0)then
             write(*,*)'Error in chkconf'
             write(*,*)'ielem=',ielem,'ineigh=',ineigh
             stop
          endif
       enddo

    enddo

    do ielem=1,nelem
       if(elem(1,ielem)==0)cycle

       do inode=1,nnode
          ineigh=eltoel(inode,ielem)

          if(ineigh==0)cycle
          ip1=elem(ltab(1,inode),ielem) 
          ip2=elem(ltab(2,inode),ielem) 
          ip3=elem(ltab(3,inode),ielem) 



          lmark(ip1)=1
          lmark(ip2)=1
          lmark(ip3)=1

          if(eltoel(1,ineigh)==ielem)then
             iview=1 
          else if(eltoel(2,ineigh)==ielem)then
             iview=2 
          else if(eltoel(3,ineigh)==ielem)then
             iview=3 
          else 
             iview=4 
          endif

          ipa=elem(ltab(1,iview),ineigh)
          if(ipa>npoin .or. ipa<=0)then
             write(*,*)'Error in chkconf,ineigh=',ineigh
             write(*,*)'ipa=',ipa
             stop
          endif
          ipb=elem(ltab(2,iview),ineigh)
          if(ipb>npoin .or. ipb<=0)then
             write(*,*)'Error in chkconf,ineigh=',ineigh
             write(*,*)'ipb=',ipb
             stop
          endif
          ipc=elem(ltab(3,iview),ineigh)
          if(ipc>npoin .or. ipc<=0)then
             write(*,*)'Error in chkconf,ineigh=',ineigh
             write(*,*)'ipc=',ipc
             stop
          endif


          if(lmark(ipa)/=1)then
             write(*,*)'Error conform 1'
             stop
          endif
          if(lmark(ipb)/=1)then
             write(*,*)'Error conform 2'
             stop
          endif
          if(lmark(ipc)/=1)then
             write(*,*)'Error conform 3'
             stop
          endif

          lmark(ip1)=0
          lmark(ip2)=0
          lmark(ip3)=0

       enddo

    enddo

    call ptoelmb(elem,nelem,npoin,nnode,ptoel1,ptoel2)

    do ipoin=1,npoin
       if(ptoel2(ipoin+1)-ptoel2(ipoin)==0)then
          write(*,*)'Error in chkconf, disconnected point:',ipoin
          stop
       endif
    enddo 

    call memchk(2_ip,istat,memor_msh,'PTOEL1','chkconf',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','chkconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','chkconf',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','chkconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','chkconf',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','chkconf',0_ip)

  end subroutine chkconf

  subroutine radius(elem,nnode,nelem,relem,coor,ndim,npoin)
    use def_kintyp, only       : ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: nnode,nelem,ndim,npoin
    integer(ip),intent(in)      :: elem(nnode,nelem)
    real(rp),intent(inout)      :: relem(5,nelem)
    real(rp),intent(in)         :: coor(ndim,npoin)
    real(rp)                    :: rx1,ry1,rz1,rx2,ry2,rz2
    real(rp)                    :: rx3,ry3,rz3,rx4,ry4,rz4
    real(rp)                    :: c1,c2,c3,nx1,ny1,nz1,nx2,ny2,nz2,c05
    real(rp)                    :: l1,l2,l3,nx3,ny3,nz3,rvol,denom
    integer(ip)                 :: ielem,ip1,ip2,ip3,ip4

    c05=0.5d+00

    do ielem=1,nelem

       if(elem(1,ielem)==0)cycle  
       ip1=elem(1,ielem)
       ip2=elem(2,ielem)
       ip3=elem(3,ielem)
       ip4=elem(4,ielem)

       rx1=coor(1,ip2)-coor(1,ip1)
       ry1=coor(2,ip2)-coor(2,ip1)
       rz1=coor(3,ip2)-coor(3,ip1)

       rx2=coor(1,ip3)-coor(1,ip1)
       ry2=coor(2,ip3)-coor(2,ip1)
       rz2=coor(3,ip3)-coor(3,ip1)

       rx3=coor(1,ip4)-coor(1,ip1)
       ry3=coor(2,ip4)-coor(2,ip1)
       rz3=coor(3,ip4)-coor(3,ip1)

       nx1= ry1*rz2-ry2*rz1  
       ny1=-rx1*rz2+rx2*rz1  
       nz1= rx1*ry2-rx2*ry1  

       nx2= ry2*rz3-ry3*rz2  
       ny2=-rx2*rz3+rx3*rz2  
       nz2= rx2*ry3-rx3*ry2 

       nx3= ry3*rz1-ry1*rz3  
       ny3=-rx3*rz1+rx1*rz3  
       nz3= rx3*ry1-rx1*ry3

       rvol=rx1*nx2+ry1*ny2+rz1*nz2  
       denom=c05/rvol 

       l1=rx1*rx1+ry1*ry1+rz1*rz1
       l2=rx2*rx2+ry2*ry2+rz2*rz2
       l3=rx3*rx3+ry3*ry3+rz3*rz3

       c1=(l1*nx2+l2*nx3+l3*nx1)*denom
       c2=(l1*ny2+l2*ny3+l3*ny1)*denom
       c3=(l1*nz2+l2*nz3+l3*nz1)*denom

       !write(*,*)relem(1,ielem),c1+coor(1,ip1)
       !write(*,*)relem(2,ielem),c2+coor(2,ip1)
       !write(*,*)relem(3,ielem),c3+coor(3,ip1)
       !write(*,*)relem(5,ielem),c1*c1+c2*c2+c3*c3


       relem(5,ielem)=c1*c1+c2*c2+c3*c3
       relem(1,ielem)=c1+coor(1,ip1)
       relem(2,ielem)=c2+coor(2,ip1)
       relem(3,ielem)=c3+coor(3,ip1)

    enddo

  end subroutine radius

  subroutine optimloc(elem,nnode,nelem,coor,npoin,ndim,lmark,eltoel,&
             ltet,relem,lenew,nenew,menew)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnode,npoin,ndim,menew,nenew
    integer(ip), intent(inout)   :: nelem,ltet(npoin),lenew(menew)
    integer(ip), pointer         :: elem(:,:),eltoel(:,:),lmark(:)
    real(rp), pointer            :: relem(:)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: ielem,jelem,iopt,j,ienew,ineigh,k
    integer(ip)                  :: maxlevel,ioptglo,ioptel
    integer(ip),parameter        :: mqual=500
    real(rp)                     :: rqual(mqual) 
    !
    !     This subroutine optimizes locally the mesh
    !
    !
    !     The stack is initialized with the newly formed elements
    !
    do ienew=1,nenew
       ielem=lenew(ienew)
       lmark(ielem)=1_ip
    enddo   
    !
    !     Swap them
    !
    do 

       ioptglo=0_ip
       ienew=0_ip

       do  

          if(ienew==nenew)exit
          ienew=ienew+1_ip

          ielem=lenew(ienew)
          if(lmark(ielem)==1)then
             ioptel=0_ip
             do j=1,nnode
                jelem=eltoel(j,ielem)
                if(jelem/=0)then
                   if(lmark(jelem)==1)then
                      iopt=0_ip
                      if(eltoel(1,jelem)==ielem)then
                         k=1_ip
                      else if(eltoel(2,jelem)==ielem)then
                         k=2_ip
                      else if(eltoel(3,jelem)==ielem)then
                         k=3_ip
                      else
                         k=4_ip
                      endif

                      call swap23(ielem,jelem,elem,nnode,nelem,eltoel,coor,ndim,&
                           npoin,iopt,j,k,lmark,rqual,ltet,mqual)
                      if(iopt==1)then
                         ioptglo=1_ip
                         ioptel=1_ip 
                         !
                         !     Mark jelem and the neighbors of ielem and jelem 
                         !     if they belong to the patch
                         !
                         ineigh=eltoel(1,ielem)
                         if(ineigh/=0)then
                            if(lmark(ineigh)>1)then
                               lmark(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(2,ielem)
                         if(ineigh/=0)then
                            if(lmark(ineigh)>1)then
                               lmark(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(3,ielem)
                         if(ineigh/=0)then
                            if(lmark(ineigh)>1)then
                               lmark(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(4,ielem)
                         if(ineigh/=0)then
                            if(lmark(ineigh)>1)then
                               lmark(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(1,jelem)
                         if(ineigh/=0)then
                            if(lmark(ineigh)>1)then
                               lmark(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(2,jelem)
                         if(ineigh/=0)then
                            if(lmark(ineigh)>1)then
                               lmark(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(3,jelem)
                         if(ineigh/=0)then
                            if(lmark(ineigh)>1)then
                               lmark(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(4,jelem)
                         if(ineigh/=0)then
                            if(lmark(ineigh)>1)then
                               lmark(ineigh)=1_ip
                            endif
                         endif
                      endif
                   endif
                endif
             enddo
             !
             !     Did we optimize something?
             !
             if(ioptel==0)then
                !
                !     Mark the element to be not checked again
                !
                lmark(ielem)=2_ip
             endif
          endif
       enddo

       if(ioptglo==0)exit
    enddo
    !
    !     Clean up
    ! 
    do ienew=1,nenew
       lmark(lenew(ienew))=0_ip 
    enddo

  end subroutine optimloc

  subroutine swap23(iel1,iel2,elem,nnode,nelem,eltoel,coor,ndim,npoin,iopt,&
       idir1,idir2,lmark,rqual,ltet,mqual)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    use mod_voltol, only : gtvol,gtqual
    implicit none
    integer(ip), intent(in)    :: nnode,npoin,ndim,iel1,iel2,idir1,idir2,mqual
    integer(ip),pointer        :: elem(:,:),eltoel(:,:),lmark(:)
    real(rp)                   :: rqual(mqual)
    real(rp), intent(in)       :: coor(ndim,npoin)
    integer(ip),intent(inout)  :: iopt,nelem,ltet(npoin)
    real(rp)                   :: c00,rvol1,rvol2,rvol3,rqual1,rqual2,rqual3  
    real(rp)                   :: qmin  
    integer(ip)                :: nelemn,ineigh1,ineigh2,ineigh3,ineigh4,ineigh5,ineigh6
    integer(ip)                :: ineigh4t,ineigh5t,ineigh6t
    integer(ip)                :: ip1,ip2,ip3,ip4,ip5,ielem1,ielem2,ielem3,ielem4,ielem5,ielem6 
    integer(ip)                :: ltab(3,4),elmloc(nnode,3)

    ltab(1,1)=2
    ltab(2,1)=3
    ltab(3,1)=4

    ltab(1,2)=1
    ltab(2,2)=4
    ltab(3,2)=3

    ltab(1,3)=1
    ltab(2,3)=2
    ltab(3,3)=4

    ltab(1,4)=1
    ltab(2,4)=3
    ltab(3,4)=2
    !
    !     This subroutine swaps two elements sharing a face
    !
    c00=0.0d+00 
    iopt=0_ip
    !
    !     Get the relevant points
    !
    ip1=elem(idir1,iel1)
    ip2=elem(ltab(1,idir1),iel1)
    ip3=elem(ltab(2,idir1),iel1)
    ip4=elem(ltab(3,idir1),iel1)
    ip5=elem(idir2,iel2)

    elmloc(1,1)=ip1
    elmloc(2,1)=ip5
    elmloc(3,1)=ip2
    elmloc(4,1)=ip3

    elmloc(1,2)=ip1
    elmloc(2,2)=ip5
    elmloc(3,2)=ip3
    elmloc(4,2)=ip4

    elmloc(1,3)=ip1
    elmloc(2,3)=ip5
    elmloc(3,3)=ip4
    elmloc(4,3)=ip2
    !
    !     Compute volume
    ! 
    call gtvol(elmloc,1_ip,coor,rvol1,nnode,3_ip,npoin,ndim)
    if(rvol1<c00)return
    call gtvol(elmloc,2_ip,coor,rvol2,nnode,3_ip,npoin,ndim)
    if(rvol2<c00)return
    call gtvol(elmloc,3_ip,coor,rvol3,nnode,3_ip,npoin,ndim)
    if(rvol3<c00)return
    !
    !     Get old quality
    !
    qmin=min(rqual(iel1),rqual(iel2)) 
    !
    !     Compute new quality
    !  
    call gtqual(elmloc,1_ip,coor,rvol1,nnode,3_ip,ndim,rqual1,npoin)
    if(rqual1>qmin)return
    call gtqual(elmloc,2_ip,coor,rvol2,nnode,3_ip,ndim,rqual2,npoin)
    if(rqual2>qmin)return
    call gtqual(elmloc,3_ip,coor,rvol3,nnode,3_ip,ndim,rqual3,npoin)
    if(rqual3>qmin)return
    !
    !     Optimization successfull
    !
    nelemn=nelem+1_ip    
    call memrea(nelemn,memor_msh,'ELEM','swap23',elem)
    call memrea(nelemn,memor_msh,'ELTOEL','swap23',eltoel)

    ineigh1=eltoel(ltab(1,idir1),iel1)
    ineigh2=eltoel(ltab(2,idir1),iel1)
    ineigh3=eltoel(ltab(3,idir1),iel1)
    ineigh4t=eltoel(ltab(1,idir2),iel2)
    ineigh5t=eltoel(ltab(2,idir2),iel2)
    ineigh6t=eltoel(ltab(3,idir2),iel2)

    if(ip2==elem(ltab(1,idir2),iel2))then
       ineigh4=ineigh4t 
       ineigh5=ineigh6t
       ineigh6=ineigh5t
    else if(ip2==elem(ltab(2,idir2),iel2))then
       ineigh4=ineigh5t 
       ineigh5=ineigh4t
       ineigh6=ineigh6t
    else
       ineigh4=ineigh6t 
       ineigh5=ineigh5t
       ineigh6=ineigh4t
    endif

    ielem1=iel1
    ielem2=iel2
    ielem3=nelem+1

    elem(1,ielem1)=elmloc(1,1)
    elem(2,ielem1)=elmloc(2,1)
    elem(3,ielem1)=elmloc(3,1)
    elem(4,ielem1)=elmloc(4,1)

    eltoel(1,ielem1)=ineigh6
    eltoel(2,ielem1)=ielem2
    eltoel(3,ielem1)=ielem3
    eltoel(4,ielem1)=ineigh3

    elem(1,ielem2)=elmloc(1,2)
    elem(2,ielem2)=elmloc(2,2)
    elem(3,ielem2)=elmloc(3,2)
    elem(4,ielem2)=elmloc(4,2)

    eltoel(1,ielem2)=ineigh4
    eltoel(2,ielem2)=ielem3
    eltoel(3,ielem2)=ielem1
    eltoel(4,ielem2)=ineigh1

    elem(1,ielem3)=elmloc(1,3)
    elem(2,ielem3)=elmloc(2,3)
    elem(3,ielem3)=elmloc(3,3)
    elem(4,ielem3)=elmloc(4,3)

    eltoel(1,ielem3)=ineigh5
    eltoel(2,ielem3)=ielem1
    eltoel(3,ielem3)=ielem2
    eltoel(4,ielem3)=ineigh2

    nelem=nelemn

    rqual(ielem1)=rqual1
    rqual(ielem2)=rqual2
    rqual(ielem3)=rqual3

    ltet(ip1)=ielem1  
    ltet(ip5)=ielem1  
    ltet(ip2)=ielem1  
    ltet(ip3)=ielem1  
    ltet(ip4)=ielem2


  end subroutine swap23


end module mod_delau






