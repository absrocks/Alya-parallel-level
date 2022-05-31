module mod_href

  use mod_mshtol
  !
  !     Use openmp  $OMP1
  ! 

contains

  subroutine adapt(nnode,nnofa,nsid,nelem,npoin,nface,lmark,nlay,cpi,R,cvi,pres0,&
       elem,ndim,coor,htree,var,ntree,nvar,lface,lsurf,hftree,nftree,nlaymax)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nnode,nnofa,nsid,nlay,ndim,ntree,nvar,nftree,nlaymax
    integer(ip),intent(inout) :: nelem,npoin,nface
    integer(ip),intent(inout) :: lmark(nelem)
    real(rp),pointer          :: coor(:,:),var(:,:)
    integer(ip), pointer      :: lface(:,:),lsurf(:),hftree(:,:)
    integer(ip),pointer       :: htree(:,:),elem(:,:)
    real(rp),intent(in)       :: cpi,R,cvi,pres0
    !
    !     Mark the elements with respect to the pseudo estimator in lmark
    !
    call estim(nelem,lmark,nnode,elem,cpi,R,cvi,pres0,coor,var)
    !
    !     Modify the mesh accordingly
    !
    call href(nnode,nnofa,nsid,nelem,npoin,nface,lmark,nlay,ndim,coor,&
         htree,elem,var,ntree,nvar,lface,lsurf,hftree,nftree,nlaymax)

  end subroutine adapt

  subroutine estim(nelem,lmark,nnode,elem,cpi,R,cvi,pres0,coor,var)
    use def_kintyp, only       : ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nelem,nnode,elem(nnode,nelem)
    real(rp),intent(in)       :: cpi,R,cvi,pres0   
    integer(ip),intent(inout) :: lmark(nelem)
    real(rp),pointer          :: coor(:,:),var(:,:) 
    integer(ip)               :: ielem,ipoin,ino
    real(rp)                  :: tmin,tmax,coef1,coef2,rmassl
    real(rp)                  :: temppot,rho
    real(rp)                  :: c00,c10,c14,cm13,smax,smin,c05,c13,c16
    integer(ip),save          :: iter=0

    tmin=280.0d+00
    tmax=300.0d+00

    smax=100.0d+00
    smin=1.0d+00


    c00=0.0d+00
    c10=1.0d+00
    c05=0.5d+00
    c14=1.0d+00/4.0d+00
    cm13=-c10/3.0d+00
    c13=c10/3.0d+00
    c16=c10/6.0d+00

    coef1=(smax-smin)/(tmax-tmin)
    coef2=smax-coef1*tmax 

    do ielem=1,nelem
       lmark(ielem)=0_ip
    enddo

    if(iter<5)then 

       iter=iter+1

       do ielem=1,nelem
          !
          !     Compute potential temperature
          !  
          temppot=c00
          rmassl=c00
          rho=c00
          do ino=1,nnode
             ipoin=elem(ino,ielem) 
             rho=rho+var(1,ipoin)
             !rhoi=c10/var(1,ipoin)
             !u=var(2,ipoin)*rhoi
             !v=var(3,ipoin)*rhoi
             !w=var(4,ipoin)*rhoi
             !e=var(5,ipoin)*rhoi
             !v2=c05*u*u+v*v+w*w
             !temp=cvi*(e-v2)
             !prespot=(var(6,ipoin)/pres0)**(R*cpi)
             !temppot=temppot+temp/prespot
          enddo

          !temppot=temppot*c14 
          rho=rho*c14
          !
          !     Compute element size
          !
          !ip1=elem(1,ielem)
          !ip2=elem(2,ielem)
          !ip3=elem(3,ielem)
          !ip4=elem(4,ielem)

          !x1  =coor(1,ip1)
          !y1  =coor(2,ip1)
          !z1  =coor(3,ip1)
          !dx21=coor(1,ip2)-x1
          !dy21=coor(2,ip2)-y1
          !dz21=coor(3,ip2)-z1
          !dx31=coor(1,ip3)-x1
          !dy31=coor(2,ip3)-y1
          !dz31=coor(3,ip3)-z1
          !dx41=coor(1,ip4)-x1
          !dy41=coor(2,ip4)-y1
          !dz41=coor(3,ip4)-z1
          !dx32=coor(1,ip3)-coor(1,ip2)
          !dy32=coor(2,ip3)-coor(2,ip2)
          !dz32=coor(3,ip3)-coor(3,ip2)
          !dx42=coor(1,ip4)-coor(1,ip2)
          !dy42=coor(2,ip4)-coor(2,ip2)
          !dz42=coor(3,ip4)-coor(3,ip2)
          !dx43=coor(1,ip4)-coor(1,ip3)
          !dy43=coor(2,ip4)-coor(2,ip3)
          !dz43=coor(3,ip4)-coor(3,ip3)

          !rselem=c16*(dx41*(dy21*dz31-dz21*dy31)+dy41*(dz21*dx31-dx21*dz31)+dz41*(dx21*dy31-dy21*dx31))
          !rselem=rselem**c13

          !dl21=abs(dx21)+abs(dy21)+abs(dz21)
          !dl31=abs(dx31)+abs(dy31)+abs(dz31)
          !dl41=abs(dx41)+abs(dy41)+abs(dz41)
          !rselem=max(dl21,dl31,dl41) 
          !dl21=sqrt(dx21*dx21+dy21*dy21+dz21*dz21)
          !dl31=sqrt(dx31*dx31+dy31*dy31+dz31*dz31)
          !dl41=sqrt(dx41*dx41+dy41*dy41+dz41*dz41)
          !dl32=sqrt(dx32*dx32+dy32*dy32+dz32*dz32)
          !dl42=sqrt(dx42*dx42+dy42*dy42+dz42*dz42)
          !dl43=sqrt(dx43*dx43+dy43*dy43+dz43*dz43)

          !rselem=min(dl21,dl31,dl41,dl32,dl42,dl43)
          !rselem=0.9d+00*rselem 

          !temppot=max(temppot,tmin)
          !temppot=min(temppot,tmax)
          !
          !     Get expected size
          ! 
          !rsize=coef1*temppot+coef2
          !
          !     Compare to element size
          !
          !if(rsize<rselem)then
          !   lmark(ielem)=1_ip
          !endif

          !
          !     If positive mark the element
          !
          if(rho>c00)lmark(ielem)=1_ip 

       enddo

    else

       lmark=-1

    endif

  end subroutine estim

  subroutine href(nnode,nnofa,nsid,nelem,npoin,nface,lmark,nlay,ndim,coor,&
       htree,elem,var,ntree,nvar,lface,lsurf,hftree,nftree,nlaymax)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        : memor_msh
    use mod_mshtol 
    implicit none
    integer(ip),intent(in)    :: nnode,nnofa,nsid,nlay,ndim,ntree,nvar,nftree,nlaymax
    integer(ip),intent(inout) :: nelem,npoin,nface
    integer(ip),intent(inout) :: lmark(nelem)
    integer(ip), pointer      :: lsmark(:),lmarkp(:),elem(:,:),htree(:,:)
    integer(ip), pointer      :: lface(:,:),lsurf(:),hftree(:,:)
    integer(ip)               :: nedge,nelem0,npoin0
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),ledge(:,:)!,eltoel(:,:)
    integer(ip),pointer       :: ledglm(:,:),ledgfa(:,:),ptoed2(:)
    logical(lg),pointer       :: lpoin(:)
    real(rp),pointer          :: coor(:,:),var(:,:) 
    integer(ip)               :: ichk,ielem
    integer(4)                :: istat
    !
    !     This subroutine drives the h-refinement method
    !     
    nullify(ptoel1,ptoel2,ptoed2,ledglm,ledgfa,ledge)
    !
    !     On input: -lmark tells if the element has been marked for refinement
    !                     lmark(ielem)= 1 -> refine
    !                     lmark(ielem)=-1 -> coarsen
    !               -htree the tree of the mesh 
    !
    !     Know your ennemy
    !
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(lmark,nelem) &
    !$omp& private(ielem,ichk)

    do ielem=1,nelem 
       ichk=0_ip
       if(lmark(ielem)==-1)ichk=1
       if(lmark(ielem)== 0)ichk=1
       if(lmark(ielem)== 1)ichk=1

       if(ichk/=1)then
          write(*,*)'Error input href, lmark(ielem)=',lmark(ielem)
          stop
       endif
    enddo   
    if(nnode/=4)then
       write(*,*)'nnode must be 4'
       stop
    endif  
    if(nnofa/=3)then
       write(*,*)'nnofa must be 3'
       stop
    endif  
    if(nsid/=6)then
       write(*,*)'nsid must be 6'
       stop
    endif  
    !
    !    Allocate necessary help arrays 
    !
    allocate(ledglm(nsid,nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDGLM','href',ledglm)
    allocate(ledgfa(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDGFA','href',ledgfa)
    allocate(lpoin(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOIN','href',lpoin)
    !
    !     If first time, allocate htree
    !     Else it should have been resized and updated
    !
    if(.not.associated(htree))then
       allocate(htree(ntree,nelem),stat=istat)
       call memchk(zero,istat,memor_msh,'HTREE','href',htree)
    endif
    if(.not.associated(hftree))then
       allocate(hftree(nftree,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'LFTREE','href',hftree)
    endif
    !
    !     Verify coherence
    !
    call verif(nelem,nnode,elem,htree,ntree,hftree,nftree,nface)
    !
    !     Verify face to volume coherence
    !
    call verifface(nface,nnofa,lface,nnode,nelem,elem,npoin)
    !
    !     First reorder the elements and lmark
    !
    call reorder(elem,nnode,nelem,htree,lmark,ntree)
    !
    !     Get the elements surrounding the points
    !
    call ptoelm(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !     Get  all the  edges  
    !
    call ptoedg2(elem,nelem,npoin,nnode,ptoel1,ptoel2,nedge,ptoed2,ledge)
    !
    !     Get the element to edge correspondance
    !  
    call edglem(elem,nnode,nelem,ledglm,nsid,ptoed2,ledge,nedge)
    !
    !     Get the face to edge correspondance
    !  
    call edglem2(lface,nnofa,nface,ledgfa,nnofa,ptoed2,ledge,nedge)
    !
    !     Get the elements surrounding elements
    !
    !call tetote(elem,nnode,nelem,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Allocate lsmark
    !
    allocate(lsmark(nedge),stat=istat)
    call memchk(zero,istat,memor_msh,'LSMARK','href',lsmark)
    allocate(lmarkp(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKP','href',lmarkp)
    !
    !     Limit the level of refinement 
    !
    call limitref(nelem,ntree,htree,lmark,nlaymax)
    !
    !     Add protective layers
    !
    call addlayer(nelem,elem,lmark,nnode,nlay,npoin)
    !
    !     Remember the element and point number for coarsening after
    !
    nelem0=nelem
    npoin0=npoin
    !
    !     Get a conforming subdvision
    !
    call divconf(lmark,lsmark,nedge,ledge,npoin,elem,nnode,nelem,ledglm,nsid,&
         nlay,htree,ntree)
    !
    !     Correct the 4 elements on the boundaries (Not necessary neither general) 
    !
    !call corbou(nelem,nnode,eltoel,nsid,ledglm,lsmark,nedge,htree)
    !
    !     Output the marked elements
    !
    call outhref(nnode,nelem,elem,ndim,npoin,coor,lmark)
    !xmin=1962.0d+00
    !ymin=249.0d+00 
    !zmin=-0.1d+00
    !xmax=2209.0d+00
    !ymax=501.0d+00 
    !zmax=0.1d+00

    !do iface=1,nface
    !   ichk=0_ip 
    !   do inofa=1,nnofa 
    !      ipoin=lface(inofa,iface)
    !      if(coor(1,ipoin)>xmax)ichk=1
    !      if(coor(2,ipoin)>ymax)ichk=1
    !      if(coor(3,ipoin)>zmax)ichk=1
    !      if(coor(1,ipoin)<xmin)ichk=1
    !      if(coor(2,ipoin)<ymin)ichk=1
    !      if(coor(3,ipoin)<zmin)ichk=1
    !   enddo
       
    !   if(ichk==0)then 
    !      write(*,*)'Trouve, iface=',iface
    !   endif 
    !enddo
    !
    !     Get a conforming coarsening
    !
    call coaconf(nelem,nnode,npoin,nedge,htree,lmark,elem,&
                 lsmark,ledglm,nsid,lmarkp,ntree,lpoin)
    !
    !     Get the new points
    ! 
    call newpoi(npoin,nelem,nsid,ledglm,lsmark,nedge,htree,&
         ledge,var,coor,ntree,nvar,lface,nface,nnofa,ledgfa,&
         lsurf,lmarkp,hftree)
    !
    !     Refine the mesh
    !
    !call verif(nelem,nnode,elem,htree,ntree,hftree,nftree,nface)
    call hrefi(ledglm,nelem,nnode,nsid,nedge,lmark,elem,htree,&
               ntree,coor,ndim,npoin)
    !call verif(nelem,nnode,elem,htree,ntree,hftree,nftree,nface)
    !
    !     DBG
    !
    !call ptoelm(elem,nelem,npoin,nnode)
    !call tetote(elem,nnode,nelem,ptoel1,ptoel2,npoin,eltoel)
    !call dbghref2(eltoel,nnode,nelem,elem,htree)
    !
    !     Coarsen the mesh
    !
    call hcoar(lmarkp,nelem0,nnode,htree,elem,lmark,npoin0,ntree,coor,ndim,lpoin)
    !
    !     Reorder the surviving elements
    !
    call reorderd(elem,nnode,nelem,htree,lmarkp,coor,ndim,npoin,npoin0,var,&
         ntree,nvar,nface,nnofa,lface)
    !
    !     Update htree 
    !
    call update(nelem,htree,elem,nnode,ntree,hftree,nftree,nface)
    !
    !     DBG 
    !
    !call outhref2(nnode,nelem,elem,ndim,npoin,coor)
    !xmin=1962.0d+00
    !ymin=249.0d+00 
    !zmin=-0.1d+00
    !xmax=2209.0d+00
    !ymax=501.0d+00 
    !zmax=0.1d+00

    !do iface=1,nface
    !   ichk=0_ip 
    !   do inofa=1,nnofa 
    !      ipoin=lface(inofa,iface)
    !      if(coor(1,ipoin)>xmax)ichk=1
    !      if(coor(2,ipoin)>ymax)ichk=1
    !      if(coor(3,ipoin)>zmax)ichk=1
    !      if(coor(1,ipoin)<xmin)ichk=1
    !      if(coor(2,ipoin)<ymin)ichk=1
    !      if(coor(3,ipoin)<zmin)ichk=1
    !   enddo
    !   
    !   if(ichk==0)then 
    !      write(*,*)'Trouve, iface=',iface
    !   endif 
    !enddo
 
    !
    !     Verify coherence
    !
    call verif(nelem,nnode,elem,htree,ntree,hftree,nftree,nface)
    !
    !     Verify face to volume coherence
    !
    call verifface(nface,nnofa,lface,nnode,nelem,elem,npoin)

    call memchk(2_ip,istat,memor_msh,'LMARKP','href',lmarkp)
    deallocate(lmarkp,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARKP','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSMARK','href',lsmark)
    deallocate(lsmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSMARK','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGLM','href',ledglm)
    deallocate(ledglm,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGLM','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGFA','href',ledgfa)
    deallocate(ledgfa,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGFA','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPOIN','href',lpoin)
    deallocate(lpoin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPOIN','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','href',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','href',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOED2','href',ptoed2)
    deallocate(ptoed2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOED2','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGE','href',ledge)
    deallocate(ledge,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGE','href',0_ip)

  end subroutine href

  subroutine addlayer(nelem,elem,lmark,nnode,nlay,npoin)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)    :: nnode,nelem,nlay,npoin
    integer(ip),intent(in)    :: elem(nnode,nelem)
    integer(ip),intent(inout) :: lmark(nelem)
    integer(ip),pointer       :: lmarkp(:)
    integer(ip)               :: ilay,ielem,ncont0,ncont1
    integer(4)                :: istat
    !
    !     This sub adds protective layer to the marked elements
    !
    allocate(lmarkp(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKP','addlayer',lmarkp)
    !
    !     Diagnostic purposes
    !
    ncont0=0
  
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(lmark,nelem) &
    !$omp& reduction(+:ncont0) &
    !$omp& private(ielem)
    
    do ielem=1,nelem
       if(lmark(ielem)==1)then
          ncont0=ncont0+1
       endif
    enddo
    !
    !     Loop on layer number, element-to-point-to-element version
    !
    do ilay=1,nlay 

       lmarkp=0_ip
       do ielem=1,nelem
          if(lmark(ielem)==1)then
             lmarkp(elem(1,ielem))=1
             lmarkp(elem(2,ielem))=1
             lmarkp(elem(3,ielem))=1
             lmarkp(elem(4,ielem))=1
          endif
       enddo

       do ielem=1,nelem
          if(lmarkp(elem(1,ielem))==1)then
             lmark(ielem)=1
             cycle
          endif
          if(lmarkp(elem(2,ielem))==1)then
             lmark(ielem)=1
             cycle
          endif
          if(lmarkp(elem(3,ielem))==1)then
             lmark(ielem)=1
             cycle
          endif
          if(lmarkp(elem(4,ielem))==1)then
             lmark(ielem)=1
             cycle
          endif
       enddo

       ncont1=0
       do ielem=1,nelem
          if(lmark(ielem)==1)then
             ncont1=ncont1+1
          endif
       enddo
       write(*,*)'Iteration:',ilay,'Marked:',ncont1

    enddo

    ncont1=0
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(lmark,nelem) &
    !$omp& reduction(+:ncont1) &
    !$omp& private(ielem)
   
    do ielem=1,nelem
       if(lmark(ielem)==1)then
          ncont1=ncont1+1
       endif
    enddo

    write(*,*)'Elements marked before added layers:',ncont0
    write(*,*)'Elements marked after added layers:',ncont1

    call memchk(2_ip,istat,memor_msh,'LMARKP','addlayer',lmarkp)
    deallocate(lmarkp,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARKP','addlayer',0_ip)

  end subroutine addlayer

  subroutine divconf(lmark,lsmark,nedge,ledge,npoin,elem,nnode,nelem,&
       ledglm,nsid,nlay,htree,ntree)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)    :: nelem,nedge,npoin,nnode,nlay,nsid,ntree     
    integer(ip),intent(in)    :: elem(nnode,nelem)     
    integer(ip),intent(in)    :: ledglm(nsid,nelem),ledge(2,nedge)
    integer(ip),intent(in)    :: htree(ntree,nelem)
    integer(ip),intent(inout) :: lsmark(nedge),lmark(nelem)     
    integer(ip)               :: iedg,iedge,ielem,ncoun,lnode(4)
    integer(ip)               :: lsid(2,6),ltab2(6,6,6),idg,iftyp,iter,icont
    integer(ip)               :: ilevel,jneigh,istack,ine,ineigh,neighlevel
    integer(ip)               :: nstack,jelem,kelem,ncont,nconte,jne,ledgl(6)
    integer(ip)               :: jlevel,klevel,iopt,ifath,icontn,ison1,ison2
    integer(ip),pointer       :: lstack(:),lemark(:)
    integer(4)                :: istat

    allocate(lstack(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','divconf',lstack)
    allocate(lemark(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LEMARK','divconf',lemark)

    lsid(1,1)=1
    lsid(2,1)=2
    lsid(1,2)=1
    lsid(2,2)=3
    lsid(1,3)=1
    lsid(2,3)=4
    lsid(1,4)=2
    lsid(2,4)=3
    lsid(1,5)=2
    lsid(2,5)=4
    lsid(1,6)=3
    lsid(2,6)=4
    !
    !     Obtain correct face
    !
    ltab2=0
    ltab2(1,2,4)=4
    ltab2(1,3,5)=3
    ltab2(2,3,6)=2  
    ltab2(4,5,6)=1
    !
    !     First loop on the marked elements to be divided
    !     to mark the edges
    !     Connection through edges can not be reached through neighbors
    !
    ncont=0_ip
    do ielem=1,nelem
       if(lmark(ielem)==1)then
          call mrkedg(htree,ielem,nelem,nsid,ledglm,lsmark,nedge,ntree)
          ncont=ncont+1_ip
       endif
    enddo
    icontn=0_ip
    do iedge=1,nedge
       if(lsmark(iedge)==1)then  
          icontn=icontn+1_ip
       endif
    enddo
    write(*,*)'Edges marked for refinement initially:',icontn,'out of:',nedge
    !
    !     Loop on the mesh
    !  
    iter=0

    do  
       iter=iter+1
       icont=0
       !
       !     Loop from the smallest elements
       !
       do ielem=nelem,1,-1
          !
          !     Has this element been marked for refinement?
          !
          !if(lmark(ielem)/=1)then
          !     
          !     Check sides
          !
          !write(*,*)'iter=',iter,'ielem=',ielem
          !if(nedge>=845034)then
          !   if(lsmark(845034)==1)then 
          !      write(*,*)'Edge 845034 ',lsmark(845034)
          !      write(*,*)'iter=',iter,'ielem=',ielem
          !   endif
          !endif
          nconte=0_ip
          do idg=1,nsid
             if(lsmark(ledglm(idg,ielem))==1)then
                nconte=nconte+1
                ledgl(nconte)=idg
             endif
          enddo
          !
          !     Do we have something to do?
          !   
          if(nconte==0)cycle
          !
          !     Check cases
          !
          if(htree(10,ielem)==0 .or. abs(htree(10,ielem))==8)then
             !
             !     Case 0 or 8
             !
             !
             !     If one side marked, will be divided in two
             !
             if(nconte==1)cycle 
             !
             !     If one face marked, will be divided in four 
             !     if the sides belong to the same face
             !
             if(nconte==3)then  
                iftyp=ltab2(ledgl(1),ledgl(2),ledgl(3))
                if(iftyp/=0)cycle
             endif

          else if(htree(10,ielem)==2)then 
             !
             !     If a 2 element with the second side marked, will be divided in four  
             !     else must go for division in eight             
             !
             if(nconte==1 )then
                if(lsmark(ledglm(2,ielem))==1)then
                   !
                   !     Get son 
                   !
                   ison1=htree(1,ielem)
                   !
                   !     Mark appropriate edge
                   ! 
                   iedge=ledglm(4,ison1)
                   lsmark(iedge)=1_ip
                   cycle

                else if(lsmark(ledglm(3,ielem))==1)then
                   !
                   !     Get son 
                   !
                   ison1=htree(1,ielem)
                   !
                   !     Mark appropriate edge
                   ! 
                   iedge=ledglm(5,ison1)
                   lsmark(iedge)=1_ip
                   cycle

                endif
             endif

          else if(htree(10,ielem)==-2)then
             !
             !     If a -2 element with the fourth side marked, will be divided in four
             !     else must go for division in 8 
             !             
             if(nconte==1)then
                if(lsmark(ledglm(4,ielem))==1)then
                   !
                   !     Get father 
                   !
                   ifath=htree(8,ielem)
                   !
                   !     Mark appropriate edge
                   ! 
                   iedge=ledglm(2,ifath)
                   lsmark(iedge)=1_ip
                   cycle

                else if(lsmark(ledglm(5,ielem))==1)then
                   !
                   !     Get father 
                   !
                   ifath=htree(8,ielem)
                   !
                   !     Mark appropriate edge
                   ! 
                   iedge=ledglm(3,ifath)
                   lsmark(iedge)=1_ip
                   cycle

                endif
             endif

          else  if(htree(10,ielem)==4)then 
             !
             !     If a 4 element and nconte /=0 --> must refine
             !     Check for 4:8+
             ! 
             if(lsmark(ledglm(1,ielem))==1 .and. lsmark(ledglm(2,ielem))==1)then
                lsmark(ledglm(4,ielem))=1_ip
             endif

          else if(htree(10,ielem)==-4)then
             !
             !     Third son case
             !
             if(htree(9,ielem)==3)then
                !
                !     Must avoid only two marked sides
                !     Count inner marked sides
                !             
                nconte=lsmark(ledglm(1,ielem))+lsmark(ledglm(2,ielem))+&
                     lsmark(ledglm(4,ielem)) 
             
                if(nconte==2)then
                   !
                   !     Get father and sons as the mesh may have been renumbered
                   !
                   ifath=htree(8,ielem)

                   if(lsmark(ledglm(1,ielem))==0)then
                      lsmark(ledglm(1,ielem))=1_ip
                      ison1=htree(1,ifath)
                      lsmark(ledglm(1,ison1))=1_ip
                      lsmark(ledglm(4,ison1))=1_ip
                   else if(lsmark(ledglm(4,ielem))==0)then
                      lsmark(ledglm(4,ielem))=1_ip
                      ison2=htree(2,ifath)
                      lsmark(ledglm(1,ison2))=1_ip
                      lsmark(ledglm(4,ison2))=1_ip
                   else 
                      lsmark(ledglm(2,ielem))=1_ip
                      lsmark(ledglm(1,ifath))=1_ip
                      lsmark(ledglm(2,ifath))=1_ip
                   endif
                endif
             else 
                !
                !     If a -4 element and nconte /=0 --> must refine
                !     Check for 4:8+
                ! 
                if(lsmark(ledglm(1,ielem))==1 .and. lsmark(ledglm(4,ielem))==1)then
                   lsmark(ledglm(2,ielem))=1_ip
                endif
             endif

          else
             write(*,*)'Error divconf case unknown'
             stop 
          endif
          !
          !     Mark the element
          !     
          !icont=icont+1
          lmark(ielem)=1_ip
          call mrkedg(htree,ielem,nelem,nsid,ledglm,lsmark,nedge,ntree)

          !endif
       enddo
       !
       !     Now count the marked sides
       !
       icont=0_ip
       do iedge=1,nedge
          if(lsmark(iedge)==1)then  
             icont=icont+1_ip
          endif
       enddo
       write(*,*)'Refine iter:',iter,'Marked sides:',icont
       if(icont==icontn)exit
       icontn=icont  

    enddo
    !
    !     Finally put lmark to zero for elements with at least one side
    !     to be coherent with coarsening
    !
    do ielem=nelem,1,-1
       if(lmark(ielem)/=1)then
          do idg=1,nsid
             if(lsmark(ledglm(idg,ielem))==1)then
                lmark(ielem)=0_ip
             endif
          enddo
       endif
    enddo   

    call memchk(2_ip,istat,memor_msh,'LEMARK','divconf',lemark)
    deallocate(lemark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEMARK','divconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','divconf',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','divconf',0_ip)

  end subroutine divconf

  subroutine mrkedg(htree,ielem,nelem,nsid,ledglm,lsmark,nedge,ntree)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nelem,ielem,nsid,nedge,ntree
    integer(ip),intent(in)      :: htree(ntree,nelem),ledglm(nsid,nelem)
    integer(ip), intent(inout)  :: lsmark(nedge)
    integer(ip)                 :: isid,iedge,ifath,ison1,ison2,ison3
    !
    !     This subroutine marks the edges depending of the element type
    !
    if(abs(htree(10,ielem))==8 .or. htree(10,ielem)==0)then
       !
       !     First case 8 or 0 --> mark all the edges
       ! 
       do isid=1,nsid
          iedge=ledglm(isid,ielem)
          lsmark(iedge)=1
       enddo

    else if(htree(10,ielem)==2)then
       !
       !     Second case 2 --> mark the edges of the highest level
       !
       iedge=ledglm(2,ielem)
       lsmark(iedge)=1

       iedge=ledglm(3,ielem)
       lsmark(iedge)=1

       iedge=ledglm(6,ielem)
       lsmark(iedge)=1
       !
       !     The son must refine
       !
       ison1=htree(1,ielem)

       iedge=ledglm(4,ison1)
       lsmark(iedge)=1

       iedge=ledglm(5,ison1)
       lsmark(iedge)=1

    else if(htree(10,ielem)==-2)then

       iedge=ledglm(4,ielem)
       lsmark(iedge)=1

       iedge=ledglm(5,ielem)
       lsmark(iedge)=1

       iedge=ledglm(6,ielem)
       lsmark(iedge)=1
       !
       !     The father must refine
       !
       ifath=htree(8,ielem)
       iedge=ledglm(2,ifath)
       lsmark(iedge)=1

       iedge=ledglm(3,ifath)
       lsmark(iedge)=1

    else if(htree(10,ielem)==4)then
       !
       !     Third case  4 --> mark the highest edges
       !
       iedge=ledglm(3,ielem)
       lsmark(iedge)=1
       !
       !     The sons must refine
       !
       ison1=htree(1,ielem)
       ison2=htree(2,ielem)

       iedge=ledglm(5,ison1)
       lsmark(iedge)=1
       iedge=ledglm(5,ison2)
       lsmark(iedge)=1

    else if(htree(10,ielem)==-4)then 
       !
       !     Check the position of the son (The third son is completely surrounded)
       !
       if(htree(9,ielem)/=3)then
          !
          !     The father and sons must refine
          !
          ifath=htree(8,ielem)

          iedge=ledglm(3,ifath)
          lsmark(iedge)=1

          ison1=htree(1,ifath)
          ison2=htree(2,ifath)

          iedge=ledglm(5,ison1)
          lsmark(iedge)=1
          iedge=ledglm(5,ison2)
          lsmark(iedge)=1

       endif

    endif

  end subroutine mrkedg

  subroutine newpoi(npoin,nelem,nsid,ledglm,lsmark,nedge,htree,&
       ledge,var,coor,ntree,nvar,lface,nface,nnofa,ledgfa,&
       lsurf,lmarkp,hftree)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)    :: nelem,nedge,nsid,ntree,nvar,nnofa
    integer(ip),intent(inout) :: npoin,nface
    integer(ip),intent(in)    :: htree(ntree,nelem),ledge(2,nedge)
    integer(ip),intent(in)    :: ledgfa(nnofa,nface),lmarkp(npoin)
    integer(ip),intent(inout) :: ledglm(nsid,nelem),lsmark(nedge)
    real(rp),pointer          :: var(:,:),coor(:,:)
    integer(ip),pointer       :: lface(:,:),lsurf(:),hftree(:,:)
    integer(ip),pointer       :: eltoel(:,:),ptoel1(:),ptoel2(:),lrenu(:)
    logical(lg),pointer       :: lmark(:)
    integer(ip)               :: ielem,npoi0,isid,iedge,ivar,lsid(3,3)  
    integer(ip)               :: ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8,ncont  
    integer(ip)               :: iface1,iface2,iface3,iface4,iface5,iface6,nfac0,iface,inofa 
    integer(ip)               :: jp1,jp2,jp3,nffound,lffound(100),jface,jneigh
    integer(ip)               :: ledgl(2,3),iview,isto,ipoin,ncont2,jnofa,jedge
    integer(ip)               :: ifacet,iedge1,iedge2,iview1,iview2,nfnew 
    integer(ip)               :: ineigh,iffound,icoun,lcoun(100),ipkeep 
    real(rp)                  :: c05,c10,rx1,ry1,rz1,rx2,ry2,rz2,rnl,cscal1,cscal2
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    integer(4)                :: istat

    c05=0.5d+00
    c10=1.0d+00
    lsid(1,1)=2
    lsid(2,1)=3
    lsid(3,1)=1
    lsid(1,2)=3
    lsid(2,2)=1
    lsid(3,2)=2
    lsid(1,3)=1
    lsid(2,3)=2
    lsid(3,3)=3
    !
    !     This sub creates the new points for conformity
    !     We assume that lsmark has been marked to 1 for the new edges
    !     lsmark is destructed and ledglm modified
    !
    nullify(ptoel1,ptoel2,eltoel)
    allocate(lmark(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','href',lmark)
    !
    npoi0=npoin
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Loop on all the elements 
    !
    do ielem=1,nelem
       !
       !     Loop on sides of the element
       !  
       do isid=1,nsid
          iedge=ledglm(isid,ielem)
          if(lsmark(iedge)==1)then
             !
             !     Has the side been marked for refinement?
             !
             npoin=npoin+1
             lsmark(iedge)=-npoin
             ledglm(isid,ielem)=npoin
             !
             !     Reallocate point arrays
             !
             call memrea(npoin,memor_msh,'COOR','newpoi',coor)
             ip1=ledge(1,iedge)
             ip2=ledge(2,iedge)
             coor(1,npoin)=c05*(coor(1,ip1)+coor(1,ip2))  
             coor(2,npoin)=c05*(coor(2,ip1)+coor(2,ip2))  
             coor(3,npoin)=c05*(coor(3,ip1)+coor(3,ip2))  
             !
             !     Reallocate variable
             !
             if(nvar/=0)then
                call memrea(npoin,memor_msh,'VAR','newpoi',var)
             endif

             do ivar=1,nvar
                var(ivar,npoin)=c05*(var(ivar,ip1)+var(ivar,ip2))
             enddo

          else if(lsmark(iedge)<0)then
             !
             !     Has the side already been assigned a point number?
             !     Pick up the new point number
             !
             ledglm(isid,ielem)=-lsmark(iedge)

          else

             ledglm(isid,ielem)=0_ip 

          endif

       enddo
    enddo
    !
    !     Loop on the faces for refinement
    !
    nfac0=nface
    do iface=1,nfac0
       !
       !     Did we already mark the face to avoid it?
       !
       if(lmark(iface) .eqv. .true.)cycle
       !
       !     DBG
       !
       !if(nface>=756)then
       !write(*,*)'iface=',iface,lface(1,756)        
       !endif
       !write(*,*)'iface=',iface        
       !
       !     If a 2-son, do not consider as it will be considered by the father
       !
       if(hftree(6,iface)==-2)cycle

       ncont=0_ip
       do inofa=1,nnofa
          !
          !     Get side number
          !
          iedge=ledgfa(inofa,iface)
          if(lsmark(iedge)<0)then
             ncont=ncont+1
             ledgl(1,ncont)=inofa
             ledgl(2,ncont)=iedge
          endif
       enddo
       !
       !     Do we have something to do?
       !
       if(ncont==0)cycle
       !
       !     One side marked?
       !  
       if(ncont==1)then
          !
          !     Two possibilities:
          !     - It may come from a 1:2 refinement      
          !     - It may come from a 4:8 refinement
          !     --> Must find the other pair if it exists
          !
          isid=ledgl(1,1)
          iedge=ledgl(2,1)
          ip1=lsid(1,isid)
          ip2=lsid(2,isid)
          ip3=lsid(3,isid)  
          ip1=lface(ip1,iface)       
          ip2=lface(ip2,iface)       
          ip3=lface(ip3,iface)       
          ip4=-lsmark(iedge)
          !
          !     Do we have a neighbor associated with iface 
          !
          if(hftree(6,iface)/=2)then 
             !
             !     No other face is associated with iface, split in two 1:2
             !
             iface1=iface
             iface2=nface+1
             lmark(iface)=.true.

             nface=iface2
             call memrea(nface,memor_msh,'LFACE','newpoi',lface)
             call memrea(nface,memor_msh,'LSURF','newpoi',lsurf)
             call memrea(nface,memor_msh,'LFTREE','newpoi',hftree)
             lface(1,iface1)=ip1
             lface(2,iface1)=ip4
             lface(3,iface1)=ip3
             hftree(1,iface1)=iface2
             hftree(2,iface1)=0_ip
             hftree(3,iface1)=0_ip
             hftree(6,iface1)=2_ip
             hftree(7,iface1)=hftree(7,iface)+1_ip
             hftree(8,iface1)=hftree(8,iface)+1_ip
             lsurf(iface1)=lsurf(iface)

             lface(1,iface2)=ip4
             lface(2,iface2)=ip2
             lface(3,iface2)=ip3
             hftree(1,iface2)=0_ip
             hftree(2,iface2)=0_ip
             hftree(3,iface2)=0_ip
             hftree(4,iface2)=iface1
             hftree(5,iface2)=1_ip
             hftree(6,iface2)=-2_ip
             hftree(7,iface2)=0_ip
             hftree(8,iface2)=hftree(8,iface1)
             lsurf(iface2)=lsurf(iface)

          else if(hftree(6,iface)==2)then
             !
             !     We have a face associated with iface  
             !     Split in four and maybe more 2:4+ 
             !
             jface=hftree(1,iface)  
             lmark(jface)=.true.
             lmark(iface)=.true.
             !
             !     Get the points that matter
             !
             ip1=lface(1,iface)
             ip4=lface(2,iface)
             ip3=lface(3,iface)
             ip2=lface(2,jface)
             ip5=-lsmark(ledgfa(1,jface))
             ip6=-lsmark(ledgfa(2,iface))
             !
             !     Resize with the worst case
             !
             nfnew=nface+3 
             call memrea(nfnew,memor_msh,'LFACE','newpoi',lface)
             call memrea(nfnew,memor_msh,'LSURF','newpoi',lsurf)
             call memrea(nface,memor_msh,'LFTREE','newpoi',hftree)

             iface1=iface
             iface2=jface
             iface3=nface+1
             iface4=iface3+1
             iface5=iface4+1

             lface(1,iface1)=ip6
             lface(2,iface1)=ip4
             lface(3,iface1)=ip5
             lsurf(iface1)=lsurf(iface)
             hftree(1,iface1)=iface2
             hftree(2,iface1)=iface3
             hftree(3,iface1)=iface4
             hftree(6,iface1)=4_ip             
             hftree(7,iface1)=hftree(7,iface)
             hftree(8,iface1)=hftree(8,iface)

             lface(1,iface2)=ip5
             lface(2,iface2)=ip3
             lface(3,iface2)=ip6
             lsurf(iface2)=lsurf(iface)
             hftree(1,iface2)=0_ip
             hftree(2,iface2)=0_ip
             hftree(3,iface2)=0_ip
             hftree(4,iface2)=iface1
             hftree(5,iface2)=1_ip
             hftree(6,iface2)=-4_ip
             hftree(7,iface2)=0_ip
             hftree(8,iface2)=hftree(8,iface1)

             nface=nface+1
             lface(1,iface3)=ip1
             lface(2,iface3)=ip4
             lface(3,iface3)=ip6
             lsurf(iface3)=lsurf(iface)
             hftree(1,iface3)=0_ip
             hftree(2,iface3)=0_ip
             hftree(3,iface3)=0_ip
             hftree(4,iface3)=iface1
             hftree(5,iface3)=2_ip
             hftree(6,iface3)=-4_ip
             hftree(7,iface3)=0_ip
             hftree(8,iface3)=hftree(8,iface1)
             !
             !     Now check the sub edge (only the one of jface)
             ! 
             iedge2=ledgfa(3,jface)
             if(lsmark(iedge2)<0)then

                ip8=-lsmark(iedge2) 
                nface=nface+1
                lface(1,iface4)=ip4
                lface(2,iface4)=ip8
                lface(3,iface4)=ip5
                lsurf(iface4)=lsurf(iface)
                hftree(1,iface4)=iface5
                hftree(2,iface4)=0_ip
                hftree(3,iface4)=0_ip
                hftree(4,iface4)=iface1
                hftree(5,iface4)=3_ip
                hftree(6,iface4)=2_ip
                hftree(7,iface4)=1_ip
                hftree(8,iface4)=hftree(8,iface1)+1_ip

                nface=nface+1
                lface(1,iface5)=ip8
                lface(2,iface5)=ip2
                lface(3,iface5)=ip5
                lsurf(iface5)=lsurf(iface)
                hftree(1,iface5)=0_ip
                hftree(2,iface5)=0_ip
                hftree(3,iface5)=0_ip
                hftree(4,iface5)=iface4
                hftree(5,iface5)=1_ip
                hftree(6,iface5)=-2_ip
                hftree(7,iface5)=0_ip
                hftree(8,iface5)=hftree(8,iface4)

             else

                nface=nface+1
                lface(1,iface4)=ip4
                lface(2,iface4)=ip2
                lface(3,iface4)=ip5
                lsurf(iface4)=lsurf(iface)
                hftree(1,iface4)=0_ip
                hftree(2,iface4)=0_ip
                hftree(3,iface4)=0_ip
                hftree(4,iface4)=iface1
                hftree(5,iface4)=3_ip
                hftree(6,iface4)=-4_ip
                hftree(7,iface4)=0_ip
                hftree(8,iface4)=hftree(8,iface1)

             endif

          else

             write(*,*)'Error in hftree'
             stop
          endif

       else if(ncont==2)then
          !
          !     Should be a face of a 2 element marked for refinement 2:4+
          !
          if(hftree(6,iface)/=2)then
             write(*,*)'Error in 2 face'
             stop 
          endif

          jface=hftree(1,iface)
          lmark(jface)=.true.
          lmark(iface)=.true.
          !
          !     Get the points that matter
          !
          ip4=lface(2,iface)
          ip3=lface(3,iface)
          ip1=lface(1,iface)
          ip2=lface(2,jface)
          ip5=-lsmark(ledgfa(1,jface))
          ip6=-lsmark(ledgfa(2,iface))
          !
          !     Resize with the worst case
          !
          nfnew=nface+4 
          call memrea(nfnew,memor_msh,'LFACE','newpoi',lface)
          call memrea(nfnew,memor_msh,'LSURF','newpoi',lsurf)
          call memrea(nface,memor_msh,'LFTREE','newpoi',hftree)

          iface1=iface
          iface2=jface
          iface3=nface+1
          iface4=iface3+1
          iface5=iface4+1
          iface6=iface5+1

          lface(1,iface1)=ip6
          lface(2,iface1)=ip4
          lface(3,iface1)=ip5
          lsurf(iface1)=lsurf(iface)
          hftree(1,iface1)=iface2
          hftree(2,iface1)=iface3
          hftree(3,iface1)=iface4
          hftree(4,iface1)=hftree(4,iface)
          hftree(5,iface1)=hftree(5,iface)
          hftree(6,iface1)=4_ip             
          hftree(7,iface1)=hftree(7,iface)
          hftree(8,iface1)=hftree(8,iface)

          lface(1,iface2)=ip5
          lface(2,iface2)=ip3
          lface(3,iface2)=ip6
          lsurf(iface2)=lsurf(iface)
          hftree(1,iface2)=0_ip
          hftree(2,iface2)=0_ip
          hftree(3,iface2)=0_ip
          hftree(4,iface2)=iface1
          hftree(5,iface2)=1_ip
          hftree(6,iface2)=-4_ip
          hftree(7,iface2)=0_ip
          hftree(8,iface2)=hftree(8,iface1)
          !
          !     Now check the sub edges
          ! 
          iedge1=ledgfa(3,iface)
          if(lsmark(iedge1)<0)then

             ip7=-lsmark(iedge1) 
             nface=nface+1
             lface(1,iface3)=ip1
             lface(2,iface3)=ip7
             lface(3,iface3)=ip6
             lsurf(iface3)=lsurf(iface)
             hftree(1,iface3)=iface5
             hftree(2,iface3)=0_ip
             hftree(3,iface3)=0_ip
             hftree(4,iface3)=iface1
             hftree(5,iface3)=2_ip
             hftree(6,iface3)=2_ip
             hftree(7,iface3)=1_ip
             hftree(8,iface3)=hftree(8,iface1)+1_ip

             nface=nface+1
             lface(1,iface5)=ip7
             lface(2,iface5)=ip4
             lface(3,iface5)=ip6
             lsurf(iface5)=lsurf(iface)
             hftree(1,iface5)=0_ip
             hftree(2,iface5)=0_ip
             hftree(3,iface5)=0_ip
             hftree(4,iface5)=iface3
             hftree(5,iface5)=1_ip
             hftree(6,iface5)=-2_ip
             hftree(7,iface5)=0_ip
             hftree(8,iface5)=hftree(8,iface3)

          else

             nface=nface+1
             lface(1,iface3)=ip1
             lface(2,iface3)=ip4
             lface(3,iface3)=ip6
             lsurf(iface3)=lsurf(iface)
             hftree(1,iface3)=0_ip
             hftree(2,iface3)=0_ip
             hftree(3,iface3)=0_ip
             hftree(4,iface3)=iface1
             hftree(5,iface3)=2_ip
             hftree(6,iface3)=-4_ip
             hftree(7,iface3)=0_ip
             hftree(8,iface3)=hftree(8,iface1)

          endif

          iedge2=ledgfa(3,jface)
          if(lsmark(iedge2)<0)then

             ip8=-lsmark(iedge2) 
             nface=nface+1
             lface(1,iface4)=ip4
             lface(2,iface4)=ip8
             lface(3,iface4)=ip5
             lsurf(iface4)=lsurf(iface)
             hftree(1,iface4)=iface6
             hftree(2,iface4)=0_ip
             hftree(3,iface4)=0_ip
             hftree(4,iface4)=iface1
             hftree(5,iface4)=3_ip
             hftree(6,iface4)=2_ip
             hftree(7,iface4)=1_ip
             hftree(8,iface4)=hftree(8,iface1)+1_ip

             nface=nface+1
             lface(1,iface6)=ip8
             lface(2,iface6)=ip2
             lface(3,iface6)=ip5
             lsurf(iface6)=lsurf(iface)
             hftree(1,iface6)=0_ip
             hftree(2,iface6)=0_ip
             hftree(3,iface6)=0_ip
             hftree(4,iface6)=iface4
             hftree(5,iface6)=1_ip
             hftree(6,iface6)=-2_ip
             hftree(7,iface6)=0_ip
             hftree(8,iface6)=hftree(8,iface4)

          else

             nface=nface+1
             lface(1,iface4)=ip4
             lface(2,iface4)=ip2
             lface(3,iface4)=ip5
             lsurf(iface4)=lsurf(iface)
             hftree(1,iface4)=0_ip
             hftree(2,iface4)=0_ip
             hftree(3,iface4)=0_ip
             hftree(4,iface4)=iface1
             hftree(5,iface4)=3_ip
             hftree(6,iface4)=-4_ip
             hftree(7,iface4)=0_ip
             hftree(8,iface4)=hftree(8,iface1)

          endif

       else if(ncont==3)then
          !
          !     1:4 refinement
          !
          iface1=iface
          iface2=nface+1
          iface3=iface2+1 
          iface4=iface3+1 
          lmark(iface)=.true.

          nface=iface4
          call memrea(nface,memor_msh,'LFACE','newpoi',lface)
          call memrea(nface,memor_msh,'LSURF','newpoi',lsurf)
          call memrea(nface,memor_msh,'LFTREE','newpoi',hftree)
          ip1=lface(1,iface)       
          ip2=lface(2,iface)       
          ip3=lface(3,iface)       
          ip4=-lsmark(ledgl(2,1))
          ip5=-lsmark(ledgl(2,2))
          ip6=-lsmark(ledgl(2,3))

          lface(1,iface1)=ip4
          lface(2,iface1)=ip5
          lface(3,iface1)=ip6
          lsurf(iface1)=lsurf(iface)  
          hftree(1,iface1)=iface2
          hftree(2,iface1)=iface3
          hftree(3,iface1)=iface4
          hftree(4,iface1)=hftree(4,iface)
          hftree(5,iface1)=hftree(5,iface)
          hftree(6,iface1)=4_ip
          hftree(7,iface1)=hftree(7,iface)+1_ip
          hftree(8,iface1)=hftree(8,iface)+1_ip

          lface(1,iface2)=ip4
          lface(2,iface2)=ip3
          lface(3,iface2)=ip5
          lsurf(iface2)=lsurf(iface)  
          hftree(1,iface2)=0_ip
          hftree(2,iface2)=0_ip
          hftree(3,iface2)=0_ip
          hftree(4,iface2)=iface1
          hftree(5,iface2)=1_ip
          hftree(6,iface2)=-4_ip
          hftree(7,iface2)=0_ip
          hftree(8,iface2)=hftree(8,iface)

          lface(1,iface3)=ip1
          lface(2,iface3)=ip6
          lface(3,iface3)=ip5
          lsurf(iface3)=lsurf(iface)  
          hftree(1,iface3)=0_ip
          hftree(2,iface3)=0_ip
          hftree(3,iface3)=0_ip
          hftree(4,iface3)=iface1
          hftree(5,iface3)=2_ip
          hftree(6,iface3)=-4_ip
          hftree(7,iface3)=0_ip
          hftree(8,iface3)=hftree(8,iface)

          lface(1,iface4)=ip6
          lface(2,iface4)=ip2
          lface(3,iface4)=ip4
          lsurf(iface4)=lsurf(iface)  
          hftree(1,iface4)=0_ip
          hftree(2,iface4)=0_ip
          hftree(3,iface4)=0_ip
          hftree(4,iface4)=iface1
          hftree(5,iface4)=3_ip
          hftree(6,iface4)=-4_ip
          hftree(7,iface4)=0_ip
          hftree(8,iface4)=hftree(8,iface)

       else

          write(*,*)'Error in newpoi, iface=',iface,'ncont=',ncont
          stop
       endif

    enddo
    !
    !     Now take care of coarsened faces
    !
    do iface=1,nfac0
       !
       !     DBG
       !
       !if(nface>=756)then
       !write(*,*)'iface=',iface,lface(1,756)        
       !endif
       !write(*,*)'iface=',iface        
       !
       !     Only consider parent
       !
       if(hftree(6,iface)<=0)cycle 
       !
       !     Has this face been already marked by refinement?
       !     No needed in coarsening as new faces have hftree(10,iface)=0
       !  
       if(lmark(iface) .eqv. .true.)cycle
       !
       !     Check the 4 elements
       !
       if(hftree(6,iface)==4)then  
          !
          !     Get the remaining points
          !
          ip1=lface(1,iface)
          ip2=lface(2,iface)
          ip3=lface(3,iface)

          ncont=lmarkp(ip1)+lmarkp(ip2)+lmarkp(ip3)
          !
          !     Will this face be deleted?
          !
          if(ncont==0)cycle
          !
          !     Collapse 4:1
          !
          if(ncont==3)then  
             !
             !     Get the three sons
             !
             iface2=hftree(1,iface)
             iface3=hftree(2,iface)
             iface4=hftree(3,iface)
             !
             !     Find the three remaining points
             !
             if(eltoel(1,iface2)==iface)then
                iview=1_ip
             else if(eltoel(2,iface2)==iface)then
                iview=2_ip
             else
                iview=3_ip
             endif
             ip1=lface(iview,iface2)  
             if(eltoel(1,iface3)==iface)then
                iview=1_ip
             else if(eltoel(2,iface3)==iface)then
                iview=2_ip
             else
                iview=3_ip
             endif
             ip2=lface(iview,iface3)  
             if(eltoel(1,iface4)==iface)then
                iview=1_ip
             else if(eltoel(2,iface4)==iface)then
                iview=2_ip
             else
                iview=3_ip
             endif
             ip3=lface(iview,iface4)  
             !
             !     DBG
             !
             if(lmarkp(ip1)+lmarkp(ip2)+lmarkp(ip3)/=0)then
                write(*,*)'Error in newpoi coarsen b' 
                stop
             endif

             lface(1,iface)=ip1
             lface(2,iface)=ip2
             lface(3,iface)=ip3
             lsurf(iface)=lsurf(iface)  
             hftree(1,iface)=0_ip
             hftree(2,iface)=0_ip
             hftree(3,iface)=0_ip
             hftree(6,iface)=0_ip
             hftree(7,iface)=hftree(7,iface)-1_ip
             hftree(8,iface)=hftree(8,iface)-1_ip

             lface(1,iface2)=0_ip
             hftree(1,iface2)=0_ip
             hftree(2,iface2)=0_ip
             hftree(3,iface2)=0_ip
             hftree(4,iface2)=0_ip
             hftree(5,iface2)=0_ip
             hftree(6,iface2)=0_ip
             hftree(7,iface2)=0_ip
             hftree(8,iface2)=0_ip

             lface(1,iface3)=0_ip
             hftree(1,iface3)=0_ip
             hftree(2,iface3)=0_ip
             hftree(3,iface3)=0_ip
             hftree(4,iface3)=0_ip
             hftree(5,iface3)=0_ip
             hftree(6,iface3)=0_ip
             hftree(7,iface3)=0_ip
             hftree(8,iface3)=0_ip

             lface(1,iface4)=0_ip
             hftree(1,iface4)=0_ip
             hftree(2,iface4)=0_ip
             hftree(3,iface4)=0_ip
             hftree(4,iface4)=0_ip
             hftree(5,iface4)=0_ip
             hftree(6,iface4)=0_ip
             hftree(7,iface4)=0_ip
             hftree(8,iface4)=0_ip

          else if(ncont==2)then 
             !
             !     Collapse 4:2
             !
             !
             !     Get non marked point
             !
             if(lmarkp(lface(1,iface))/=1)then
                ipkeep=1_ip
             else if(lmarkp(lface(2,iface))/=1)then
                ipkeep=2_ip
             else if(lmarkp(lface(3,iface))/=1)then
                ipkeep=3_ip
             else 
                write(*,*)'Error in newpoi, ipkeep not found'
                stop
             endif
             !
             !     Get the three neighboring faces
             !
             iface2=eltoel(ltab(1,ipkeep),iface)
             iface3=eltoel(ltab(2,ipkeep),iface)
             iface4=eltoel(ipkeep,iface)
             !
             !     Find the three remaining points
             !
             if(eltoel(1,iface2)==iface)then
                iview=1_ip
             else if(eltoel(2,iface2)==iface)then
                iview=2_ip
             else
                iview=3_ip
             endif
             ip1=lface(iview,iface2)  
             if(eltoel(1,iface3)==iface)then
                iview=1_ip
             else if(eltoel(2,iface3)==iface)then
                iview=2_ip
             else
                iview=3_ip
             endif
             ip2=lface(iview,iface3)  
             if(eltoel(1,iface4)==iface)then
                iview=1_ip
             else if(eltoel(2,iface4)==iface)then
                iview=2_ip
             else
                iview=3_ip
             endif
             ip3=lface(iview,iface4)  
             ipkeep=lface(ipkeep,iface) 
             !
             !     DBG
             !
             if(lmarkp(ip1)+lmarkp(ip2)+lmarkp(ip3)/=0)then
                write(*,*)'Error in newpoi coarsen a' 
                stop
             endif

             lface(1,iface)=ip1
             lface(2,iface)=ipkeep
             lface(3,iface)=ip3
             lsurf(iface)=lsurf(iface)  
             hftree(1,iface)=iface2
             hftree(2,iface)=0_ip
             hftree(3,iface)=0_ip
             hftree(6,iface)=2_ip
             hftree(7,iface)=hftree(7,iface)
             hftree(8,iface)=hftree(8,iface)

             lface(1,iface2)=ipkeep
             lface(2,iface2)=ip2
             lface(3,iface2)=ip3
             lsurf(iface2)=lsurf(iface)  
             hftree(1,iface2)=0_ip
             hftree(2,iface2)=0_ip
             hftree(3,iface2)=0_ip
             hftree(4,iface2)=iface 
             hftree(5,iface2)=1_ip 
             hftree(6,iface2)=-2_ip
             hftree(7,iface2)=0_ip
             hftree(8,iface2)=hftree(8,iface)

             lface(1,iface3)=0_ip
             hftree(1,iface3)=0_ip
             hftree(2,iface3)=0_ip
             hftree(3,iface3)=0_ip
             hftree(4,iface3)=0_ip
             hftree(5,iface3)=0_ip
             hftree(6,iface3)=0_ip
             hftree(7,iface3)=0_ip
             hftree(8,iface3)=0_ip

             lface(1,iface4)=0_ip
             hftree(1,iface4)=0_ip
             hftree(2,iface4)=0_ip
             hftree(3,iface4)=0_ip
             hftree(4,iface4)=0_ip
             hftree(5,iface4)=0_ip
             hftree(6,iface4)=0_ip
             hftree(7,iface4)=0_ip
             hftree(8,iface4)=0_ip

          else
             write(*,*)'Strange collapse'
             stop
          endif

       else if(hftree(6,iface)==2)then
          !
          !     Collapse 2:1
          !
          ip1=lface(2,iface)
          !
          !     Do we have something to do?
          !
          if(lmarkp(ip1)/=1)cycle          
          !
          !     Get son
          !
          jface=hftree(1,iface)
          !
          !     Get the points that matter
          !   
          ip1=lface(1,iface)
          ip4=lface(2,iface)
          ip3=lface(3,iface)
          ip2=lface(2,jface)

          lface(1,iface)=ip1
          lface(2,iface)=ip2
          lface(3,iface)=ip3
          lsurf(iface)=lsurf(iface) 
          hftree(1,iface)=0_ip
          hftree(2,iface)=0_ip
          hftree(3,iface)=0_ip
          hftree(6,iface)=0_ip
          hftree(7,iface)=hftree(7,iface)-1_ip
          hftree(8,iface)=hftree(8,iface)-1_ip

          lface(1,jface)=0_ip
          hftree(1,jface)=0_ip
          hftree(2,jface)=0_ip
          hftree(3,jface)=0_ip
          hftree(4,jface)=0_ip
          hftree(5,jface)=0_ip
          hftree(6,jface)=0_ip
          hftree(7,jface)=0_ip
          hftree(8,jface)=0_ip

       endif

    enddo
    !
    !     Allocate lrenu
    !
    allocate(lrenu(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','newpoi',lrenu)
    !
    !     Compact the faces
    ! 
    nfac0=nface
    nface=0_ip
    do iface=1,nfac0
       if(lface(1,iface)/=0)then
          nface=nface+1_ip
          lrenu(iface)=nface
          lface(1,nface)=lface(1,iface)
          lface(2,nface)=lface(2,iface)
          lface(3,nface)=lface(3,iface)     
          lsurf(nface)=lsurf(iface) 
          hftree(1,nface)=hftree(1,iface) 
          hftree(2,nface)=hftree(2,iface) 
          hftree(3,nface)=hftree(3,iface) 
          hftree(4,nface)=hftree(4,iface) 
          hftree(5,nface)=hftree(5,iface) 
          hftree(6,nface)=hftree(6,iface) 
          hftree(7,nface)=hftree(7,iface) 
          hftree(8,nface)=hftree(8,iface) 
       endif
    enddo
    !
    !     Update hftree
    !     Reorder sons if necessary 
    !     and father   
    ! 
    do iface=1,nface
       ineigh=hftree(1,iface)
       if(ineigh>0)then
          hftree(1,iface)=lrenu(ineigh)
       endif
       ineigh=hftree(2,iface)
       if(ineigh>0)then
          hftree(2,iface)=lrenu(ineigh)
       endif
       ineigh=hftree(3,iface)
       if(ineigh>0)then
          hftree(3,iface)=lrenu(ineigh)
       endif
       ineigh=hftree(4,iface)
       if(ineigh>0)then
          hftree(4,iface)=lrenu(ineigh)
       endif
    enddo




    do iface=1,nface
       if(hftree(7,iface)>hftree(8,iface))then
          write(*,*)'Error in hftree'
          stop
       endif
       do inofa=1,nnofa 
          if(lface(inofa,iface)==0)then
             write(*,*)'Error in lface'
             stop
          endif
          ipoin=lface(inofa,iface)
          if(ipoin<=npoi0)then
             if(lmarkp(ipoin)==1)then
                write(*,*)'Error in lface lmarkp'
                stop
             endif
          endif
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','href',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'LRENU','href',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','href',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','href',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','href',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','href',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','href',0_ip)

  end subroutine newpoi

  subroutine coaconf(nelem,nnode,npoin,nedge,htree,lmark,elem,lsmark,&
             ledglm,nsid,lmarkp,ntree,lpoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nelem,nedge,npoin,nnode,ntree     
    integer(ip),intent(in)    :: elem(nnode,nelem)
    integer(ip),intent(in)    :: nsid,ledglm(nsid,nelem),lsmark(nedge)     
    integer(ip),intent(inout) :: lmark(nelem),htree(ntree,nelem),lmarkp(npoin)     
    logical(lg),intent(inout) :: lpoin(npoin)
    integer(ip)               :: ielem,ncoun,ncou,ltab2(6,6,6)
    integer(ip)               :: nchil,je,jelem,ifath,nploc,lploc(10)
    integer(ip)               :: iftyp,ilevel,iter,ipos,ipoin,ncont,iel,lpout(2,4),npout
    integer(ip)               :: nelem1,nelem2,nelem3,nelem4,nelem5,nelem6,nelem7
    integer(ip)               :: ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8,ip9,ip10
    !
    !     This sub gives conforming coarsening
    !
    !
    !     Templates for elements:
    !
    !
    !     8-->1: elements  1,2,3,8:          1 5 6 7                   
    !                                        5 2 8 9
    !                                        8 3 6 10
    !                                        7 9 10 4 
    !
    !            Then depends on the diagonal ...                  
    !
    !
    !
    !
    !     Obtain correct face
    !
    ltab2=0
    ltab2(1,2,4)=4
    ltab2(1,3,5)=3
    ltab2(2,3,6)=2  
    ltab2(4,5,6)=1
    !
    !     Count elements marked for coarsening
    !
    ncont=0_ip
    do ielem=1,nelem
       if(lmark(ielem)==-1)then
          ncont=ncont+1_ip
       endif
    enddo
    write(*,*)'Elements marked for coarsening initially:',ncont,'out of:',nelem
    !
    !     First filter out non compatible father/son relations and original elements
    !
    do ielem=1,nelem
       !
       !     Has the element been marked for coarsening?
       !
       if(lmark(ielem)==-1)then
          !
          !     First consider parents
          !      
          if(htree(10,ielem)>0)then
             !
             !     Level of refinement
             !
             ilevel=htree(12,ielem)
             !
             !     Number of sons
             !
             nchil=htree(10,ielem)-1  
             !
             !     Examine the children (marked + same level)
             !     for 8 as may have been divided further
             ! 
             do je=1,nchil
                jelem=htree(je,ielem)
                !
                !     Do we have more than one level?
                !
                if(htree(12,jelem)/=ilevel)exit 
                !
                !     Have the sons been marked?
                !
                if(lmark(jelem)/=-1)exit

             enddo
             !
             !     Are the children compatible?
             !     If not, reset lmark
             !
             if(je/=nchil+1)then

                lmark(ielem)=0_ip  
                !
                !     Reset coarsening for sons of same level
                !
                do je=1,nchil
                   jelem=htree(je,ielem)
                   !
                   !     Do we have more than one level?
                   !
                   if(htree(12,jelem)==ilevel)then
                      !
                      !     Has this element been marked for coarsening? 
                      !
                      if(lmark(jelem)==-1)then

                         lmark(jelem)=0_ip

                      endif
                   endif
                enddo

             endif

          else if(htree(10,ielem)==0)then
             !
             !     Do not coarsen initial elements 
             !
             lmark(ielem)=0

          else
             !
             !     Check father/son relationship 
             !
             ifath=htree(8,ielem)
             !
             !     If the father has not the same type, mark the son
             ! 
             if(htree(10,ifath)/=-htree(10,ielem))then
                lmark(ielem)=0_ip
             else
                !
                !     If the father has not the son, mark the son 
                !     May occur with 8/-8 relationship
                !
                nchil=htree(10,ifath)-1
                do je=1,nchil
                   if(htree(je,ifath)==ielem)exit
                enddo

                if(je==nchil+1)then 
                   lmark(ielem)=0_ip
                else
                   !
                   !     If the father has not been marked for deletion, mark the son
                   !
                   if(lmark(ifath)/=-1)then
                      lmark(ielem)=0_ip
                   endif 
                endif
             endif
          endif
       endif
    enddo
    !
    !     Monitoring purposes
    !
    ncont=0_ip
    do ielem=1,nelem
       if(lmark(ielem)==-1)then
          ncont=ncont+1_ip
       endif
    enddo
    write(*,*)'Elements marked for coarsening after first filtering:',ncont,'out of:',nelem
    !
    !     Obtain the list of "total deletion points"
    !     by marking the points to be deleted for each type of element
    !
    lmarkp=0_ip

    do ielem=1,nelem
       !
       !     Has the element been marked for deletion?
       !
       if(lmark(ielem)==-1 .and. htree(10,ielem)>0)then
          !
          !     Which type of division do we have?
          !
          if(htree(10,ielem)==8)then
             !
             !     Case 8. These elements have all their sons with -8.
             !     The points are not ordered  
             ! 
             !
             !     Get four inner sons
             !
             jelem=htree(4,ielem)
             lmarkp(elem(1,jelem))=1_ip
             lmarkp(elem(2,jelem))=1_ip
             lmarkp(elem(3,jelem))=1_ip
             lmarkp(elem(4,jelem))=1_ip
             jelem=htree(5,ielem)
             lmarkp(elem(1,jelem))=1_ip
             lmarkp(elem(2,jelem))=1_ip
             lmarkp(elem(3,jelem))=1_ip
             lmarkp(elem(4,jelem))=1_ip
             jelem=htree(6,ielem)
             lmarkp(elem(1,jelem))=1_ip
             lmarkp(elem(2,jelem))=1_ip
             lmarkp(elem(3,jelem))=1_ip
             lmarkp(elem(4,jelem))=1_ip
             jelem=htree(7,ielem)
             lmarkp(elem(1,jelem))=1_ip
             lmarkp(elem(2,jelem))=1_ip
             lmarkp(elem(3,jelem))=1_ip
             lmarkp(elem(4,jelem))=1_ip

          else if(htree(10,ielem)==4)then
             !
             !     Case 4 
             !     The points are ordered
             !
             jelem=htree(3,ielem)
             lmarkp(elem(1,jelem))=1_ip
             lmarkp(elem(2,jelem))=1_ip
             lmarkp(elem(3,jelem))=1_ip

          else if(htree(10,ielem)==2)then
             !
             !     Case 2
             !     The points are ordered
             !
             lmarkp(elem(2,ielem))=1_ip  

          endif
       endif
    enddo
    !
    !     Then mark the points of the untouched elements
    !
    do ielem=1,nelem
       if(lmark(ielem)/=-1)then
          lmarkp(elem(1,ielem))=0_ip   
          lmarkp(elem(2,ielem))=0_ip   
          lmarkp(elem(3,ielem))=0_ip   
          lmarkp(elem(4,ielem))=0_ip   
       endif
    enddo
    !
    !     And remark the points to be kept in the elements to be coarsened
    ! 
    do ielem=1,nelem
       !
       !     Has the element been marked for deletion?
       !
       if(lmark(ielem)==-1 .and. htree(10,ielem)>0)then
          !
          !     Which type of division do we have?
          !
          if(htree(10,ielem)==8)then
             !
             !     Case 8. These elements have all their sons with -8.
             !     Not ordered, must mark the outer points as kept 
             !

             !
             !     First mark the inner points in lpoin
             !
             do ipos=4,7
                iel=htree(ipos,ielem)
                lpoin(elem(1,iel))=.true.
                lpoin(elem(2,iel))=.true.
                lpoin(elem(3,iel))=.true.
                lpoin(elem(4,iel))=.true.
             enddo
             !
             !     Then identify the outer points
             ! 
             if(lpoin(elem(1,ielem)).eqv. .false.)then
                lmarkp(elem(1,ielem))=0_ip
             else if(lpoin(elem(2,ielem)).eqv. .false.)then
                lmarkp(elem(2,ielem))=0_ip
             else if(lpoin(elem(3,ielem)).eqv. .false.)then
                lmarkp(elem(3,ielem))=0_ip
             else
                lmarkp(elem(4,ielem))=0_ip
             endif

             do ipos=1,3
                iel=htree(ipos,ielem)  
                if(lpoin(elem(1,iel)).eqv. .false.)then
                   lmarkp(elem(1,iel))=0_ip
                else if(lpoin(elem(2,iel)).eqv. .false.)then
                   lmarkp(elem(2,iel))=0_ip
                else if(lpoin(elem(3,iel)).eqv. .false.)then
                   lmarkp(elem(3,iel))=0_ip
                else
                   lmarkp(elem(4,iel))=0_ip
                endif
             enddo
             !
             !     Clean up
             !  
             do ipos=4,7
                iel=htree(ipos,ielem)
                lpoin(elem(1,iel))=.false.
                lpoin(elem(2,iel))=.false.
                lpoin(elem(3,iel))=.false.
                lpoin(elem(4,iel))=.false.
             enddo

          else if(htree(10,ielem)==4)then
             !
             !     Case 4
             !
             lmarkp(elem(1,ielem))=0_ip
             lmarkp(elem(4,ielem))=0_ip
             jelem=htree(1,ielem)
             lmarkp(elem(2,jelem))=0_ip
             jelem=htree(2,ielem)
             lmarkp(elem(2,jelem))=0_ip

          else if(htree(10,ielem)==2)then
             !
             !     Case 2
             !
             lmarkp(elem(1,ielem))=0_ip  
             lmarkp(elem(3,ielem))=0_ip  
             lmarkp(elem(4,ielem))=0_ip  
             jelem=htree(1,ielem)
             lmarkp(elem(2,jelem))=0_ip  

          endif
       endif
    enddo
    !
    !     Diagnostic
    !
    ncoun=0_ip
    do ipoin=1,npoin
       if(lmarkp(ipoin)==1)then
          ncoun=ncoun+1
       endif
    enddo
    write(*,*)'Points marked before iteration:',ncoun
    !
    !     Loop until conforming
    !
    iter=0_ip
    do 
       !
       !     Initialize counter
       !
       ncoun=0_ip
       !
       !     Loop on the elements
       !
       do ielem=1,nelem
          !
          !     Has this element been marked for coarsening ?
          !
          if(lmark(ielem)==-1 .and. htree(10,ielem)>0)then
             !
             !     Check element type, only consider parents
             !   
             if(htree(10,ielem)==8)then
                !
                !     Case 8   
                !     The points are not ordered due to possible rotations 
                !  

                !
                !     First mark the inner points in lpoin
                !
                do ipos=4,7
                   iel=htree(ipos,ielem)
                   lpoin(elem(1,iel))=.true.
                   lpoin(elem(2,iel))=.true.
                   lpoin(elem(3,iel))=.true.
                   lpoin(elem(4,iel))=.true.
                enddo
                !
                !     Then identify the outer points
                ! 
                npout=0_ip
                if(lpoin(elem(1,ielem)).eqv. .false.)then
                   npout=npout+1_ip
                   lpout(1,npout)=1
                   lpout(2,npout)=elem(1,ielem)
                else if(lpoin(elem(2,ielem)).eqv. .false.)then
                   npout=npout+1_ip
                   lpout(1,npout)=2
                   lpout(2,npout)=elem(2,ielem)
                else if(lpoin(elem(3,ielem)).eqv. .false.)then
                   npout=npout+1_ip
                   lpout(1,npout)=3
                   lpout(2,npout)=elem(3,ielem)
                else
                   npout=npout+1_ip
                   lpout(1,npout)=4
                   lpout(2,npout)=elem(4,ielem)
                endif

                do ipos=1,3
                   iel=htree(ipos,ielem)  
                   if(lpoin(elem(1,iel)).eqv. .false.)then
                      npout=npout+1_ip
                      lpout(1,npout)=1
                      lpout(2,npout)=elem(1,iel)
                   else if(lpoin(elem(2,iel)).eqv. .false.)then
                      npout=npout+1_ip
                      lpout(1,npout)=2
                      lpout(2,npout)=elem(2,iel)
                   else if(lpoin(elem(3,iel)).eqv. .false.)then
                      npout=npout+1_ip
                      lpout(1,npout)=3
                      lpout(2,npout)=elem(3,iel)
                   else
                      npout=npout+1_ip
                      lpout(1,npout)=4
                      lpout(2,npout)=elem(4,iel)
                   endif
                enddo
                !
                !     Clean up
                !  
                do ipos=4,7
                   iel=htree(ipos,ielem)
                   lpoin(elem(1,iel))=.false.
                   lpoin(elem(2,iel))=.false.
                   lpoin(elem(3,iel))=.false.
                   lpoin(elem(4,iel))=.false.
                enddo
                !
                !     Get the sons
                ! 
                nelem1=htree(1,ielem)
                nelem2=htree(2,ielem)
                nelem3=htree(3,ielem)
                nelem4=htree(4,ielem)
                nelem5=htree(5,ielem)
                nelem6=htree(6,ielem)
                nelem7=htree(7,ielem)
                !
                !     Get the points
                ! 
                ip1=lpout(2,1) 
                ip2=lpout(2,2) 
                ip3=lpout(2,3) 
                ip4=lpout(2,4) 
                !
                !     Mark points of the father
                !
                lpoin(elem(1,ielem))=.true.
                lpoin(elem(2,ielem))=.true.
                lpoin(elem(3,ielem))=.true.
                lpoin(elem(4,ielem))=.true.
                !
                !   Get ip5
                !
                if(lpoin(elem(1,nelem1)).eqv. .true.)then
                   ip5=elem(1,nelem1)
                else if(lpoin(elem(2,nelem1)).eqv. .true.)then
                   ip5=elem(2,nelem1)
                else if(lpoin(elem(3,nelem1)).eqv. .true.)then
                   ip5=elem(3,nelem1)
                else
                   ip5=elem(4,nelem1)
                endif
                !
                !   Get ip6
                !
                if(lpoin(elem(1,nelem2)).eqv. .true.)then
                   ip6=elem(1,nelem2)
                else if(lpoin(elem(2,nelem2)).eqv. .true.)then
                   ip6=elem(2,nelem2)
                else if(lpoin(elem(3,nelem2)).eqv. .true.)then
                   ip6=elem(3,nelem2)
                else
                   ip6=elem(4,nelem2)
                endif
                !
                !   Get ip7
                !
                if(lpoin(elem(1,nelem3)).eqv. .true.)then
                   ip7=elem(1,nelem3)
                else if(lpoin(elem(2,nelem3)).eqv. .true.)then
                   ip7=elem(2,nelem3)
                else if(lpoin(elem(3,nelem3)).eqv. .true.)then
                   ip7=elem(3,nelem3)
                else
                   ip7=elem(4,nelem3)
                endif
                !
                !     Clean up
                !
                lpoin(elem(1,ielem))=.false.
                lpoin(elem(2,ielem))=.false.
                lpoin(elem(3,ielem))=.false.
                lpoin(elem(4,ielem))=.false.
                !
                !     Mark points of the nelem1
                !
                lpoin(elem(1,nelem1))=.true.
                lpoin(elem(2,nelem1))=.true.
                lpoin(elem(3,nelem1))=.true.
                lpoin(elem(4,nelem1))=.true.
                !
                !   Get ip8
                !
                if(lpoin(elem(1,nelem2)).eqv. .true.)then
                   ip8=elem(1,nelem2)
                else if(lpoin(elem(2,nelem2)).eqv. .true.)then
                   ip8=elem(2,nelem2)
                else if(lpoin(elem(3,nelem2)).eqv. .true.)then
                   ip8=elem(3,nelem2)
                else
                   ip8=elem(4,nelem2)
                endif
                !
                !   Get ip9
                !
                if(lpoin(elem(1,nelem3)).eqv. .true.)then
                   ip9=elem(1,nelem3)
                else if(lpoin(elem(2,nelem3)).eqv. .true.)then
                   ip9=elem(2,nelem3)
                else if(lpoin(elem(3,nelem3)).eqv. .true.)then
                   ip9=elem(3,nelem3)
                else
                   ip9=elem(4,nelem3)
                endif
                !
                !     Clean up
                !
                lpoin(elem(1,nelem1))=.false.
                lpoin(elem(2,nelem1))=.false.
                lpoin(elem(3,nelem1))=.false.
                lpoin(elem(4,nelem1))=.false.
                !
                !     Mark points of the nelem2
                !
                lpoin(elem(1,nelem2))=.true.
                lpoin(elem(2,nelem2))=.true.
                lpoin(elem(3,nelem2))=.true.
                lpoin(elem(4,nelem2))=.true.
                !
                !   Get ip10
                !
                if(lpoin(elem(1,nelem3)).eqv. .true.)then
                   ip10=elem(1,nelem3)
                else if(lpoin(elem(2,nelem3)).eqv. .true.)then
                   ip10=elem(2,nelem3)
                else if(lpoin(elem(3,nelem3)).eqv. .true.)then
                   ip10=elem(3,nelem3)
                else
                   ip10=elem(4,nelem3)
                endif
                !
                !     Clean up
                !
                lpoin(elem(1,nelem2))=.false.
                lpoin(elem(2,nelem2))=.false.
                lpoin(elem(3,nelem2))=.false.
                lpoin(elem(4,nelem2))=.false.
                !
                !     Count number of deleted points
                !
                ncou=lmarkp(ip5)+lmarkp(ip6)+lmarkp(ip7)+lmarkp(ip8)+lmarkp(ip9)+lmarkp(ip10)
                !
                !     Cases 8:1 and 8:2
                !          
                if(ncou==6 .or. ncou==5)cycle           
                ! 
                !     Case 8:4 more demanding case
                !
                if(ncou==3)then
                   !
                   !     Check if point belong to the same face
                   !
                   !
                   !     Find the non marked points    
                   !
                   nploc=0
                   if(lmarkp(ip5)==0)then
                      nploc=nploc+1
                      lploc(nploc)=1_ip
                   endif
                   if(lmarkp(ip6)==0)then
                      nploc=nploc+1
                      lploc(nploc)=2_ip
                   endif
                   if(lmarkp(ip7)==0)then
                      nploc=nploc+1
                      lploc(nploc)=3_ip
                   endif
                   if(lmarkp(ip8)==0)then
                      nploc=nploc+1
                      lploc(nploc)=4_ip
                   endif
                   if(lmarkp(ip9)==0)then
                      nploc=nploc+1
                      lploc(nploc)=5_ip
                   endif
                   if(lmarkp(ip10)==0)then
                      nploc=nploc+1
                      lploc(nploc)=6_ip
                   endif

                   iftyp=ltab2(lploc(1),lploc(2),lploc(3))
                   !
                   !     Do the three non marked points belong to the same face?
                   !     If yes, go home 
                   !
                   if(iftyp/=0)cycle

                endif
                !
                !     Test not passed, the points do not belong to the same face, 
                !     Take out elements and points
                !
                ncoun=ncoun+1
                lmark(ielem)=0_ip
                lmark(nelem1)=0_ip
                lmark(nelem2)=0_ip
                lmark(nelem3)=0_ip
                lmark(nelem4)=0_ip
                lmark(nelem5)=0_ip
                lmark(nelem6)=0_ip
                lmark(nelem7)=0_ip
                lmarkp(ip5)=0_ip
                lmarkp(ip6)=0_ip
                lmarkp(ip7)=0_ip
                lmarkp(ip8)=0_ip
                lmarkp(ip9)=0_ip
                lmarkp(ip10)=0_ip

             else if(htree(10,ielem)==4)then
                !
                !     Case 4
                !     The points are ordered 
                !
                nelem1=htree(1,ielem)
                nelem2=htree(2,ielem)
                nelem3=htree(3,ielem)

                ip1=elem(1,ielem) 
                ip2=elem(2,nelem1)
                ip3=elem(2,nelem2)
                ip4=elem(4,nelem3)

                ip5=elem(2,ielem)
                ip6=elem(3,nelem1)
                ip7=elem(3,ielem)

                ncou=lmarkp(ip5)+lmarkp(ip6)+lmarkp(ip7)
                !
                !     Case 4:1 and 4:2
                !
                if(ncou==3 .or. ncou==2)cycle
                !
                !     Test not passed, take out elements and points
                !
                ncoun=ncoun+1
                lmark(ielem)=0_ip
                lmark(nelem1)=0_ip
                lmark(nelem2)=0_ip
                lmark(nelem3)=0_ip
                lmarkp(ip5)=0_ip
                lmarkp(ip6)=0_ip
                lmarkp(ip7)=0_ip

             else if(htree(10,ielem)==2)then
                !
                !     Case 2:1
                !     The points are ordered 
                !
                nelem1=htree(1,ielem)

                ip1=elem(1,ielem) 
                ip3=elem(3,ielem)
                ip4=elem(4,ielem)
                ip2=elem(2,nelem1)       
                ip5=elem(2,ielem)

                if(lmarkp(ip5)==1)cycle
                !
                !     Test not passed,  
                !     Take out elements and point
                !
                ncoun=ncoun+1
                lmark(ielem)=0_ip
                lmark(nelem1)=0_ip
                lmarkp(ip5)=0_ip

             endif

          endif

       enddo

       iter=iter+1
       write(*,*)'Coarse iter:',iter,'Deleted elements:',ncoun
       if(ncoun==0)exit

    enddo

    !if(iter>1)then
    !    do ipoin=1,npoin
    !       write(*,*)ipoin,lmarkp(ipoin)
    !    enddo 
    !  endif 


  end subroutine coaconf

  subroutine hrefi(ledglm,nelem,nnode,nsid,nedge,lmark,elem,htree,ntree,coor,ndim,npoin)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)    :: nnode,nsid,nedge,ntree,ndim,npoin
    integer(ip),intent(inout) :: nelem
    integer(ip),intent(in)    :: ledglm(nsid,nelem)
    real(rp),intent(in)       :: coor(ndim,npoin)
    integer(ip),intent(inout) :: lmark(nelem)
    integer(ip)               :: ncont,ielem,ltab(2,4,4),lsid(2,6),ltri(2,4),itri
    integer(ip)               :: ltab2(6,6,6),ltab3(7,4),ledgl(6),nelem0
    integer(ip)               :: idg,iedg,ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8,ip9,ip10,inode,ipoin
    integer(ip)               :: elem0,nelem1,nelem2,nelem3,nelem4,nelem5 
    integer(ip)               :: nelem6,nelem7,nelem8,ifath 
    integer(ip)               :: iedg1,iedg2,iedg3,iftyp,ison1,ison2,ison3
    integer(ip)               :: ipedg1,ipedg2,jelem,newelem
    integer(ip)               :: lptab(10),letab(8),lptab2(10)
    integer(ip),pointer       :: elem(:,:),htree(:,:) 
    !
    !     Point to side 
    !
    lsid(1,1)=1
    lsid(2,1)=2
    lsid(1,2)=1
    lsid(2,2)=3
    lsid(1,3)=1
    lsid(2,3)=4
    lsid(1,4)=2
    lsid(2,4)=3
    lsid(1,5)=2
    lsid(2,5)=4
    lsid(1,6)=3
    lsid(2,6)=4

    ltab=0_ip
    !
    !     Obtain correct edge end points for 1:2
    !     Given two points of a tetra, obtain the two oriented others     
    !
    ltab(1,1,2)=3
    ltab(2,1,2)=4
    ltab(1,1,3)=4
    ltab(2,1,3)=2
    ltab(1,1,4)=2
    ltab(2,1,4)=3
    ltab(1,2,1)=4
    ltab(2,2,1)=3
    ltab(1,2,3)=1
    ltab(2,2,3)=4
    ltab(1,2,4)=3
    ltab(2,2,4)=1
    ltab(1,3,1)=2
    ltab(2,3,1)=4
    ltab(1,3,2)=4
    ltab(2,3,2)=1
    ltab(1,3,4)=1
    ltab(2,3,4)=2
    ltab(1,4,1)=3
    ltab(2,4,1)=2
    ltab(1,4,2)=1
    ltab(2,4,2)=3
    ltab(1,4,3)=2
    ltab(2,4,3)=1
    !
    !     Obtain correct face for 1:4
    !
    ltab2=0
    ltab2(1,2,4)=4
    ltab2(1,3,5)=3
    ltab2(2,3,6)=2  
    ltab2(4,5,6)=1
    !
    !     Points and edges of the face for 1:4
    !
    ltab3(1,1)=2
    ltab3(2,1)=4
    ltab3(3,1)=3
    ltab3(4,1)=1
    ltab3(5,1)=5
    ltab3(6,1)=6
    ltab3(7,1)=4

    ltab3(1,2)=3
    ltab3(2,2)=4
    ltab3(3,2)=1
    ltab3(4,2)=2
    ltab3(5,2)=6
    ltab3(6,2)=3
    ltab3(7,2)=2

    ltab3(1,3)=1
    ltab3(2,3)=4
    ltab3(3,3)=2
    ltab3(4,3)=3
    ltab3(5,3)=3
    ltab3(6,3)=5
    ltab3(7,3)=1

    ltab3(1,4)=1
    ltab3(2,4)=2
    ltab3(3,4)=3
    ltab3(4,4)=4
    ltab3(5,4)=1
    ltab3(6,4)=4
    ltab3(7,4)=2
    !
    !     DBG 
    !
    !call dbghref(htree,elem,nelem,nnode,ntree)

    !
    !    DBG
    !
    !do ielem=1,nelem
    !   if(htree(8,ielem)==1)then
    !      write(*,*)ielem
    !   endif
    !enddo  
    !
    !     DBG: verif that if the sons of 2 and 4 have been marked,
    !          the parents have been marked too
    !
    do ielem=1,nelem
       if(lmark(ielem)==1)then
          if(htree(10,ielem)==-4 .or. htree(10,ielem)==-2)then
             ifath=htree(8,ielem) 
             if(lmark(ifath)/=1)then
                write(*,*)'Error in hrefi, father not marked' 
                stop
             endif
          endif
       endif
    enddo


    nelem0=nelem

    do ielem=1,nelem0
       !
       !     How many sides are marked for refinement
       !
       ncont=0_ip
       do idg=1,nsid
          if(ledglm(idg,ielem)>0)then
             ncont=ncont+1
             ledgl(ncont)=idg
          endif
       enddo
       !
       !     DBG
       ! 
       !write(*,*)'ielem=',ielem
       !
       !     Only consider parents for 2 and 4
       !            
       if(htree(10,ielem)==-2 .or. htree(10,ielem)==-4)cycle
       !
       !     Has the element been marked for refinement?
       !
       if(lmark(ielem)==0)then
          !
          !     One possibility: 
          !                         - we are in a new element, 
          !     but may be divided in two or four
          !
          !     Do we have something to do?
          !  
          if(ncont==0)cycle
          !
          !     Find pattern
          !
          if(ncont==1)then
             !
             !     Check error
             !
             if(abs(htree(10,ielem))==4)then
                write(*,*)'Error in hrefi configuration not allowed'
                stop
             endif
             !
             !     2:4 configuration
             !  
             if(htree(10,ielem)==2)then
                !
                !     Get son
                !
                ison1=htree(1,ielem)
                !
                !     Verify that the son has also been marked
                !  
                if(lmark(ison1)==1)then
                   write(*,*)'Error in 2-4 configuration a'
                   stop
                endif

                nelem=nelem+2
                call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
                call memrea(nelem,memor_msh,'HTREE','hrefi',htree)
                !
                !     Two cases
                ! 
                if(ledglm(2,ielem)>0)then 

                   lptab(1)=elem(1,ielem) 
                   lptab(5)=elem(2,ielem) 
                   lptab(3)=elem(3,ielem) 
                   lptab(4)=elem(4,ielem) 
                   lptab(2)=elem(2,ison1)
                   lptab(7)=ledglm(2,ielem)
                   lptab(6)=ledglm(4,ison1)
                   letab(1)=ielem 
                   letab(2)=ison1
                   letab(3)=nelem-1 
                   letab(4)=nelem

                else if(ledglm(3,ielem)>0)then

                   lptab(1)=elem(1,ielem) 
                   lptab(5)=ledglm(3,ielem)
                   lptab(3)=elem(2,ison1) 
                   lptab(4)=elem(3,ielem) 
                   lptab(2)=elem(4,ielem)
                   lptab(7)=elem(2,ielem)
                   lptab(6)=ledglm(5,ison1)
                   letab(1)=ielem 
                   letab(2)=ison1
                   letab(3)=nelem-1 
                   letab(4)=nelem

                else

                   write(*,*)'Error case 2:4 a'
                   stop 

                endif
                !
                !     Correct level
                !
                htree(11,ielem)=htree(11,ielem)-1_ip
                htree(12,ielem)=htree(12,ielem)-1_ip
                !
                !     Mark the son to avoid it later on 
                !
                lmark(ison1)=2_ip

                call ref4(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                     letab,lptab)
                !
                !     1:2 configuration
                !  
             else if(abs(htree(10,ielem))==8 .or. htree(10,ielem)==0)then


                iedg=ledgl(1)      

                ip1=lsid(1,iedg)
                ip2=lsid(2,iedg)

                ip3=ltab(1,ip1,ip2) 
                ip4=ltab(2,ip1,ip2) 

                nelem=nelem+1
                !call dbghref(htree,elem,nelem,nnode,ntree)
                call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
                call memrea(nelem,memor_msh,'HTREE','hrefi',htree)
                !call dbghref(htree,elem,nelem,nnode,ntree)
                lptab(1)=elem(ip1,ielem) 
                lptab(2)=elem(ip2,ielem) 
                lptab(3)=elem(ip3,ielem) 
                lptab(4)=elem(ip4,ielem) 
                lptab(5)=ledglm(iedg,ielem)
                letab(1)=ielem 
                letab(2)=nelem 

                call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                     letab,lptab)

             else

                write(*,*)'Error with not marked element, strange configuration'
                stop

             endif

          else if(ncont==3)then  
             !
             !   1:4 configuration
             !
             !
             !     Check error
             !
             if(abs(htree(10,ielem))==2 .or. abs(htree(10,ielem))==4)then
                write(*,*)'Error in hrefi configuration not allowed'
                stop
             endif
             !
             !  
             !    Which face has been marked for refinement?
             ! 
             iedg1=ledgl(1)
             iedg2=ledgl(2)
             iedg3=ledgl(3)

             iftyp=ltab2(iedg1,iedg2,iedg3)
             if(iftyp==0)then
                write(*,*)'Error in hrefi 2' 
                stop  
             endif
             !
             !     Obtain the points of this face
             !
             nelem=nelem+3
             call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
             call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

             lptab(1)=elem(ltab3(1,iftyp),ielem) 
             lptab(2)=elem(ltab3(2,iftyp),ielem) 
             lptab(3)=elem(ltab3(3,iftyp),ielem) 
             lptab(4)=elem(ltab3(4,iftyp),ielem)
             lptab(5)=ledglm(ltab3(5,iftyp),ielem)
             lptab(6)=ledglm(ltab3(6,iftyp),ielem)
             lptab(7)=ledglm(ltab3(7,iftyp),ielem)
             letab(1)=ielem 
             letab(2)=nelem-2 
             letab(3)=nelem-1 
             letab(4)=nelem 

             call ref4(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                  letab,lptab)

          else
             !
             !     Error
             !
             write(*,*)'Error in hrefi for non marked element'
             stop
          endif

       else if(lmark(ielem)==1)then 
          !
          !     This element has been marked for refinement
          !     Only consider parents for 2 and 4 cases as 
          !     children will be considered at the same time
          !
          !     What kind of refinement do we have?
          ! 
          !
          !     First remember the sons before deleting them
          !
          ison1=htree(1,ielem)
          ison2=htree(2,ielem)
          ison3=htree(3,ielem)
          !
          !     Three possibilities: 
          !                         - we are in a new element, 
          !                         - we are in an element divided by 2 
          !                         - we are in an element divided by 4
          !
          if(abs(htree(10,ielem))==8 .or. htree(10,ielem)==0)then
             !
             !     Check error
             !
             if(ncont/=6)then 
                write(*,*)'Error 1:8, ncont=',ncont  
                stop
             endif
             !
             !     1:8 configuration
             !
             nelem1=ielem 
             nelem2=nelem+1
             nelem3=nelem2+1
             nelem4=nelem3+1
             nelem5=nelem4+1
             nelem6=nelem5+1
             nelem7=nelem6+1
             nelem8=nelem7+1
             !
             !     Reallocate elem and htree
             !
             nelem=nelem8 
             call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
             call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

             lptab(1)=elem(1,ielem) 
             lptab(2)=elem(2,ielem) 
             lptab(3)=elem(3,ielem) 
             lptab(4)=elem(4,ielem)
             lptab(5)=ledglm(1,ielem)
             lptab(6)=ledglm(2,ielem)
             lptab(7)=ledglm(3,ielem)
             lptab(8)=ledglm(4,ielem)
             lptab(9)=ledglm(5,ielem)
             lptab(10)=ledglm(6,ielem) 

             letab(1)=nelem1 
             letab(2)=nelem2 
             letab(3)=nelem3 
             letab(4)=nelem4
             letab(5)=nelem5 
             letab(6)=nelem6 
             letab(7)=nelem7 
             letab(8)=nelem8 

             call ref8(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                  letab,lptab)

          else if(htree(10,ielem)==2)then
             !
             !     Case 2
             !  
             if(ncont==1 .and. ledglm(2,ielem)>0)then 
                !
                !     2:4, only side 2 should have been marked
                ! 
                !
                !     Verify that the son has also been marked
                !  
                if(lmark(ison1)/=1)then
                   write(*,*)'Error in 2-4 configuration b'
                   stop
                endif

                nelem=nelem+2
                call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
                call memrea(nelem,memor_msh,'HTREE','hrefi',htree)
                !
                !     Two cases
                ! 
                if(ledglm(2,ielem)>0)then 

                   lptab(1)=elem(1,ielem) 
                   lptab(5)=elem(2,ielem) 
                   lptab(3)=elem(3,ielem) 
                   lptab(4)=elem(4,ielem) 
                   lptab(2)=elem(2,ison1)
                   lptab(7)=ledglm(2,ielem)
                   lptab(6)=ledglm(4,ison1)
                   letab(1)=ielem 
                   letab(2)=ison1
                   letab(3)=nelem-1 
                   letab(4)=nelem

                else if(ledglm(3,ielem)>0)then

                   lptab(1)=elem(1,ielem) 
                   lptab(5)=ledglm(3,ielem)
                   lptab(3)=elem(2,ison1) 
                   lptab(4)=elem(3,ielem) 
                   lptab(2)=elem(4,ielem)
                   lptab(7)=elem(2,ielem)
                   lptab(6)=ledglm(5,ison1)
                   letab(1)=ielem 
                   letab(2)=ison1
                   letab(3)=nelem-1 
                   letab(4)=nelem

                else

                   write(*,*)'Error case 2:4 a'
                   stop 

                endif
                !
                !     Mark the son to avoid it later on 
                !
                lmark(ison1)=2_ip
                !
                !     Correct level
                !
                htree(11,ielem)=htree(11,ielem)-1_ip
                htree(12,ielem)=htree(12,ielem)-1_ip

                call ref4(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                     letab,lptab)

             else if(ncont==3 .or. ncont==4)then   
                !
                !     2:8 configuration
                !     Only consider parents as sons MUST have been marked  
                !     ncont==3  --> divide in 8
                !     ncont==4  --> divide in 8+2
                !
                nelem1=ielem
                nelem2=ison1 
                nelem3=nelem+1
                nelem4=nelem3+1
                nelem5=nelem4+1
                nelem6=nelem5+1
                nelem7=nelem6+1
                nelem8=nelem7+1
                !
                !     Verify that the son has also been marked
                !  
                if(lmark(ison1)/=1)then
                   write(*,*)'Error in 2-8 configuration'
                   stop
                endif
                !
                !     Mark the son to avoid it later on 
                !
                lmark(nelem2)=2_ip
                !
                !     Reallocate elem and htree
                !
                nelem=nelem8 
                call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
                call memrea(nelem,memor_msh,'HTREE','hrefi',htree)
                !
                !     Verify that the deleted edges have not been marked 
                !     for refinement
                !
                if(ledglm(4,ielem)>0 .or. ledglm(5,ielem)>0)then
                   write(*,*)'Error in 2-8 refinement'
                   write(*,*)'Deleted edge marked for refinement'
                   stop 
                endif

                lptab(1)=elem(1,ielem)
                lptab(5)=elem(2,ielem) 
                lptab(3)=elem(3,ielem) 
                lptab(4)=elem(4,ielem) 
                lptab(2)=elem(2,ison1) 
                lptab(6)=ledglm(2,ielem)
                lptab(7)=ledglm(3,ielem)
                lptab(8)=ledglm(4,ison1)
                lptab(9)=ledglm(5,ison1)
                lptab(10)=ledglm(6,ison1)

                letab(1)=nelem1 
                letab(2)=nelem2 
                letab(3)=nelem3 
                letab(4)=nelem4
                letab(5)=nelem5 
                letab(6)=nelem6 
                letab(7)=nelem7 
                letab(8)=nelem8 
                !
                !     Correct level
                !
                htree(11,nelem1)=htree(11,nelem1)-1_ip
                htree(12,nelem1)=htree(12,nelem1)-1_ip
                !
                !     Divide in 8
                !
                call ref8(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                     letab,lptab)
                !
                !     Do we have to further subdivide the elements?
                !     Check the side that may have been marked in each of the 2 original elements 
                !     ielem and ison1 which appear in the new elements nelem1 and nelem2
                !  
                if(ledglm(1,ielem)/=0)then

                   ipedg1=ledglm(1,ielem)

                   nelem=nelem+1
                   call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
                   call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

                   lptab2(1)=lptab(1) 
                   lptab2(2)=lptab(5) 
                   lptab2(3)=lptab(6)
                   lptab2(4)=lptab(7) 
                   lptab2(5)=ipedg1
                   letab(1)=ielem 
                   letab(2)=nelem 

                   call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                        letab,lptab2)

                endif

                if(ledglm(1,ison1)/=0)then

                   ipedg2=ledglm(1,ison1)

                   nelem=nelem+1
                   call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
                   call memrea(nelem,memor_msh,'HTREE','hrefi',htree)
                   lptab2(1)=lptab(5) 
                   lptab2(2)=lptab(2) 
                   lptab2(3)=lptab(8)
                   lptab2(4)=lptab(9) 
                   lptab2(5)=ipedg2
                   letab(1)=nelem2
                   letab(2)=nelem 

                   call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                        letab,lptab2)

                endif
             else 
                write(*,*)'Error case 2 ncont=',ncont
                stop
             endif

          else if(htree(10,ielem)==4)then
             !
             !     Check error
             !
             if(ncont==5 .or. ncont==6)then 
                write(*,*)'Error 4:8, ncont=',ncont  
                stop
             endif
             !
             !     Verify that the sons have also been marked
             !  
             if(lmark(ison1)/=1)then
                write(*,*)'Error in 4-8 configuration 1'
                stop
             endif
             if(lmark(ison2)/=1)then
                write(*,*)'Error in 4-8 configuration 2'
                stop
             endif
             !
             !     ison3 is inside and may not have been marked
             !
             !if(lmark(ison3)/=1)then
             !   write(*,*)'Error in 4-8 configuration 3'
             !   stop
             !endif
             !
             !     Verify that the deleted edges have not been marked 
             !     for refinement
             !
             if(ledglm(3,ison3)>0 .or. ledglm(5,ison3)>0 .or. ledglm(6,ison3)>0)then
                write(*,*)'Error in 4-8 refinement'
                write(*,*)'Deleted edge marked for refinement'
                stop 
             endif
             !
             !     4:8 configuration
             !     Only consider parents as sons MUST have been marked  
             !
             nelem1=ielem
             nelem2=ison1 
             nelem3=ison2
             nelem4=ison3
             nelem5=nelem+1
             nelem6=nelem5+1
             nelem7=nelem6+1
             nelem8=nelem7+1
             !
             !     Mark the three sons to avoid them later on
             !
             lmark(nelem2)=2_ip
             lmark(nelem3)=2_ip
             lmark(nelem4)=2_ip
             !
             !     Reallocate elem and htree
             !
             nelem=nelem8 
             call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
             call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

             lptab(1)=elem(1,ielem)
             lptab(5)=elem(2,ielem) 
             lptab(6)=elem(3,ielem) 
             lptab(4)=elem(4,ielem) 
             lptab(2)=elem(2,ison1) 
             lptab(8)=elem(3,ison1)
             lptab(3)=elem(2,ison2)
             lptab(7)=ledglm(3,ielem)
             lptab(9)=ledglm(5,ison1)
             lptab(10)=ledglm(5,ison2)

             letab(1)=nelem1 
             letab(2)=nelem2 
             letab(3)=nelem3 
             letab(4)=nelem4
             letab(5)=nelem5 
             letab(6)=nelem6 
             letab(7)=nelem7 
             letab(8)=nelem8 
             !
             !     Correct level
             !
             htree(11,nelem1)=htree(11,nelem1)-1_ip
             htree(12,nelem1)=htree(12,nelem1)-1_ip
             !
             !     Divide in 8
             !
             call ref8(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                  letab,lptab)
             !
             !     Do we have to further subdivide the elements?
             !     Check the side that may have been marked in each of the 4 original elements 
             !     ielem, ison1, ison2 and ison3 in sides 1,2 and 4
             !     which appear in the new four elements nelem1, nelem2, nelem3, nelem8
             ! 
             call ref84(nnode,nelem,elem,npoin,ielem,ledglm,nsid,lptab,ntree,htree,coor,ndim,&
                  ison1,ison2,ison3,nelem1,nelem2,nelem3,nelem8)
             !
             !     Two inner sides with two inner elements to be divided in 1:2
             !     if marked
             !     Necessary for conformity inside tetra 
             !
             if(ledglm(4,ison3)/=0)then
                !
                !     Get the edge 4 (8-6) of ison3 now nelem5   
                !
                nelem=nelem+1
                call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
                call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

                lptab2(1)=lptab(6)    
                lptab2(2)=lptab(8) 
                lptab2(3)=lptab(10) 
                lptab2(4)=lptab(9) 
                lptab2(5)=ledglm(4,ison3)
                letab(1)=nelem5
                letab(2)=nelem 

                call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                     letab,lptab2)

             endif

             if(ledglm(2,ison3)/=0)then
                !
                !     Get the edge 2 (6-5) of ielem now nelem7   
                !
                nelem=nelem+1
                call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
                call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

                lptab2(1)=lptab(5) 
                lptab2(2)=lptab(6) 
                lptab2(3)=lptab(7) 
                lptab2(4)=lptab(9) 
                lptab2(5)=ledglm(2,ison3)
                letab(1)=nelem7
                letab(2)=nelem 

                call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
                     letab,lptab2)

             endif

          endif

       endif

    enddo
    !
    !     DBG
    !
    do ielem=1,nelem
       do inode=1,nnode
          ipoin=elem(inode,ielem)
          if(ipoin==0)then
             write(*,*)'Error at points, ielem=',ielem
             stop
          endif
       enddo
    enddo
    !
    !     DBG 
    !
    !call dbghref(htree,elem,nelem,nnode,ntree)

  end subroutine hrefi

  subroutine ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
       lelem,lpoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nelem,nnode,ntree,ndim,npoin
    integer(ip),intent(in)    :: lelem(2),lpoin(5)
    real(rp),intent(in)       :: coor(ndim,npoin)
    integer(ip),intent(inout) :: htree(ntree,nelem),elem(nnode,nelem)
    real(rp)                  :: rvol,c00
    integer(ip)               :: nelem1,nelem2
    c00=0.0d+00 
    !
    !     1:2 configuration
    !  
    !
    !         1 5 3 4 
    !         5 2 3 4
    !
    nelem1=lelem(1)
    nelem2=lelem(2)
    elem(1,nelem1)=lpoin(1) 
    elem(2,nelem1)=lpoin(5)
    elem(3,nelem1)=lpoin(3) 
    elem(4,nelem1)=lpoin(4) 
    call gtvol(elem,nelem1,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref2 a'
       stop  
    endif

    elem(1,nelem2)=lpoin(5) 
    elem(2,nelem2)=lpoin(2)
    elem(3,nelem2)=lpoin(3) 
    elem(4,nelem2)=lpoin(4) 
    call gtvol(elem,nelem2,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref2 b'
       stop  
    endif
    !
    !     Update htree 
    !
    htree(1,nelem1)=nelem2
    htree(2,nelem1)=0
    htree(3,nelem1)=0
    htree(4,nelem1)=0
    htree(5,nelem1)=0
    htree(6,nelem1)=0
    htree(7,nelem1)=0
    htree(10,nelem1)=2_ip
    htree(11,nelem1)=1_ip+htree(11,nelem1)
    htree(12,nelem1)=1_ip+htree(12,nelem1)

    htree(1,nelem2)=0
    htree(2,nelem2)=0
    htree(3,nelem2)=0
    htree(4,nelem2)=0
    htree(5,nelem2)=0
    htree(6,nelem2)=0
    htree(7,nelem2)=0
    htree(8,nelem2)=nelem1
    htree(9,nelem2)=1_ip
    htree(10,nelem2)=-2_ip
    htree(11,nelem2)=0_ip 
    htree(12,nelem2)=htree(12,nelem1) 

  end subroutine ref2

  subroutine ref4(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
       lelem,lpoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nelem,nnode,ntree
    integer(ip),intent(in)    :: ndim,npoin    
    integer(ip),intent(in)    :: lelem(4),lpoin(7)    
    real(rp),intent(in)       :: coor(ndim,npoin)
    integer(ip),intent(inout) :: htree(ntree,nelem),elem(nnode,nelem)
    real(rp)                  :: rvol,c00
    integer(ip)               :: nelem1,nelem2,nelem3,nelem4
    c00=0.0d+00 
    !
    !   1:4 configuration
    !
    !
    !   1 5 7 4
    !   5 2 6 4 
    !   6 3 7 4 
    !   5 6 7 4
    nelem1=lelem(1)
    nelem2=lelem(2)
    nelem3=lelem(3)
    nelem4=lelem(4)

    elem(1,nelem1)=lpoin(1)
    elem(2,nelem1)=lpoin(5)
    elem(3,nelem1)=lpoin(7)
    elem(4,nelem1)=lpoin(4)
    call gtvol(elem,nelem1,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref4 a'
       stop  
    endif

    elem(1,nelem2)=lpoin(5)
    elem(2,nelem2)=lpoin(2)
    elem(3,nelem2)=lpoin(6)
    elem(4,nelem2)=lpoin(4)
    call gtvol(elem,nelem2,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref4 b'
       stop  
    endif

    elem(1,nelem3)=lpoin(6)
    elem(2,nelem3)=lpoin(3)
    elem(3,nelem3)=lpoin(7)
    elem(4,nelem3)=lpoin(4)
    call gtvol(elem,nelem3,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref4 c'
       stop  
    endif

    elem(1,nelem4)=lpoin(5)
    elem(2,nelem4)=lpoin(6)
    elem(3,nelem4)=lpoin(7)
    elem(4,nelem4)=lpoin(4)
    call gtvol(elem,nelem4,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref4 d'
       stop  
    endif
    !
    !     Update htree
    !
    htree(1,nelem1)=nelem2
    htree(2,nelem1)=nelem3
    htree(3,nelem1)=nelem4
    htree(4,nelem1)=0
    htree(5,nelem1)=0
    htree(6,nelem1)=0
    htree(7,nelem1)=0
    htree(10,nelem1)=4_ip
    htree(11,nelem1)=1_ip+htree(11,nelem1)
    htree(12,nelem1)=1_ip+htree(12,nelem1)

    htree(1,nelem2)=0
    htree(2,nelem2)=0
    htree(3,nelem2)=0
    htree(4,nelem2)=0
    htree(5,nelem2)=0
    htree(6,nelem2)=0
    htree(7,nelem2)=0
    htree(8,nelem2)=nelem1
    htree(9,nelem2)=1_ip
    htree(10,nelem2)=-4_ip
    htree(11,nelem2)=0_ip 
    htree(12,nelem2)=htree(12,nelem1) 

    htree(1,nelem3)=0
    htree(2,nelem3)=0
    htree(3,nelem3)=0
    htree(4,nelem3)=0
    htree(5,nelem3)=0
    htree(6,nelem3)=0
    htree(7,nelem3)=0
    htree(8,nelem3)=nelem1
    htree(9,nelem3)=2_ip
    htree(10,nelem3)=-4_ip
    htree(11,nelem3)=0_ip 
    htree(12,nelem3)=htree(12,nelem1) 

    htree(1,nelem4)=0
    htree(2,nelem4)=0
    htree(3,nelem4)=0
    htree(4,nelem4)=0
    htree(5,nelem4)=0
    htree(6,nelem4)=0
    htree(7,nelem4)=0
    htree(8,nelem4)=nelem1
    htree(9,nelem4)=3_ip
    htree(10,nelem4)=-4_ip
    htree(11,nelem4)=0_ip 
    htree(12,nelem4)=htree(12,nelem1) 

  end subroutine ref4

  subroutine ref8(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
       lelem,lpoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nnode,ntree
    integer(ip),intent(in)    :: nelem,ndim,npoin
    integer(ip),intent(in)    :: lelem(8),lpoin(10)
    real(rp),intent(in)       :: coor(ndim,npoin)
    integer(ip),intent(inout) :: htree(ntree,nelem),elem(nnode,nelem)
    real(rp)                  :: rvol,c00
    integer(ip)               :: nelem4,nelem5,nelem6,nelem7,nelem8
    integer(ip)               :: nelem1,nelem2,nelem3

    c00=0.0d+00 
    !
    !     1:8 configuration
    !
    !   1 5 6 7 
    !   5 2 8 9 
    !   8 3 6 10
    !   7 9 10 4 
    !   8 10 6 9 
    !   10 7 6 9 
    !   7 5 6 9 
    !   5 8 6 9 
    ! 
    !   This choice has consequences in 4:8+ refinements
    !
    nelem1=lelem(1)
    nelem2=lelem(2)
    nelem3=lelem(3)
    nelem4=lelem(4)
    nelem5=lelem(5)
    nelem6=lelem(6)
    nelem7=lelem(7)
    nelem8=lelem(8)

    elem(1,nelem1)=lpoin(1)
    elem(2,nelem1)=lpoin(5)
    elem(3,nelem1)=lpoin(6)
    elem(4,nelem1)=lpoin(7)
    call gtvol(elem,nelem1,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref8 a'
       stop  
    endif

    elem(1,nelem2)=lpoin(5)
    elem(2,nelem2)=lpoin(2)
    elem(3,nelem2)=lpoin(8)
    elem(4,nelem2)=lpoin(9)
    call gtvol(elem,nelem2,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref8 b'
       stop  
    endif

    elem(1,nelem3)=lpoin(8)
    elem(2,nelem3)=lpoin(3)
    elem(3,nelem3)=lpoin(6)
    elem(4,nelem3)=lpoin(10)
    call gtvol(elem,nelem3,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref8 c'
       stop  
    endif

    elem(1,nelem4)=lpoin(7)
    elem(2,nelem4)=lpoin(9)
    elem(3,nelem4)=lpoin(10)
    elem(4,nelem4)=lpoin(4)
    call gtvol(elem,nelem4,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref8 d'
       stop  
    endif

    elem(1,nelem5)=lpoin(8)
    elem(2,nelem5)=lpoin(10)
    elem(3,nelem5)=lpoin(6)
    elem(4,nelem5)=lpoin(9)
    call gtvol(elem,nelem5,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref8 e'
       stop  
    endif

    elem(1,nelem6)=lpoin(10)
    elem(2,nelem6)=lpoin(7)
    elem(3,nelem6)=lpoin(6)
    elem(4,nelem6)=lpoin(9)
    call gtvol(elem,nelem6,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref8 f'
       stop  
    endif

    elem(1,nelem7)=lpoin(7)
    elem(2,nelem7)=lpoin(5)
    elem(3,nelem7)=lpoin(6)
    elem(4,nelem7)=lpoin(9)
    call gtvol(elem,nelem7,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref8 g'
       stop  
    endif

    elem(1,nelem8)=lpoin(5)
    elem(2,nelem8)=lpoin(8)
    elem(3,nelem8)=lpoin(6)
    elem(4,nelem8)=lpoin(9)
    call gtvol(elem,nelem8,coor,rvol,nnode,nelem,npoin,ndim)
    if(rvol<=c00)then
       write(*,*)'Negative volume ref8 h'
       stop  
    endif !
    !     Update htree
    !
    htree(1,nelem1)=nelem2
    htree(2,nelem1)=nelem3
    htree(3,nelem1)=nelem4
    htree(4,nelem1)=nelem5
    htree(5,nelem1)=nelem6
    htree(6,nelem1)=nelem7
    htree(7,nelem1)=nelem8
    htree(10,nelem1)=8_ip
    htree(11,nelem1)=1_ip+htree(11,nelem1)
    htree(12,nelem1)=1_ip+htree(12,nelem1)

    htree(1,nelem2)=0
    htree(2,nelem2)=0
    htree(3,nelem2)=0
    htree(4,nelem2)=0
    htree(5,nelem2)=0
    htree(6,nelem2)=0
    htree(7,nelem2)=0
    htree(8,nelem2)=nelem1
    htree(9,nelem2)=1_ip
    htree(10,nelem2)=-8_ip
    htree(11,nelem2)=0_ip 
    htree(12,nelem2)=htree(12,nelem1) 

    htree(1,nelem3)=0
    htree(2,nelem3)=0
    htree(3,nelem3)=0
    htree(4,nelem3)=0
    htree(5,nelem3)=0
    htree(6,nelem3)=0
    htree(7,nelem3)=0
    htree(8,nelem3)=nelem1
    htree(9,nelem3)=2_ip
    htree(10,nelem3)=-8_ip
    htree(11,nelem3)=0_ip 
    htree(12,nelem3)=htree(12,nelem1) 

    htree(1,nelem4)=0
    htree(2,nelem4)=0
    htree(3,nelem4)=0
    htree(4,nelem4)=0
    htree(5,nelem4)=0
    htree(6,nelem4)=0
    htree(7,nelem4)=0
    htree(8,nelem4)=nelem1
    htree(9,nelem4)=3_ip
    htree(10,nelem4)=-8_ip
    htree(11,nelem4)=0_ip 
    htree(12,nelem4)=htree(12,nelem1)

    htree(1,nelem5)=0
    htree(2,nelem5)=0
    htree(3,nelem5)=0
    htree(4,nelem5)=0
    htree(5,nelem5)=0
    htree(6,nelem5)=0
    htree(7,nelem5)=0
    htree(8,nelem5)=nelem1
    htree(9,nelem5)=4_ip
    htree(10,nelem5)=-8_ip
    htree(11,nelem5)=0_ip 
    htree(12,nelem5)=htree(12,nelem1) 

    htree(1,nelem6)=0
    htree(2,nelem6)=0
    htree(3,nelem6)=0
    htree(4,nelem6)=0
    htree(5,nelem6)=0
    htree(6,nelem6)=0
    htree(7,nelem6)=0
    htree(8,nelem6)=nelem1
    htree(9,nelem6)=5_ip
    htree(10,nelem6)=-8_ip
    htree(11,nelem6)=0_ip 
    htree(12,nelem6)=htree(12,nelem1) 

    htree(1,nelem7)=0
    htree(2,nelem7)=0
    htree(3,nelem7)=0
    htree(4,nelem7)=0
    htree(5,nelem7)=0
    htree(6,nelem7)=0
    htree(7,nelem7)=0
    htree(8,nelem7)=nelem1
    htree(9,nelem7)=6_ip
    htree(10,nelem7)=-8_ip
    htree(11,nelem7)=0_ip 
    htree(12,nelem7)=htree(12,nelem1)

    htree(1,nelem8)=0
    htree(2,nelem8)=0
    htree(3,nelem8)=0
    htree(4,nelem8)=0
    htree(5,nelem8)=0
    htree(6,nelem8)=0
    htree(7,nelem8)=0
    htree(8,nelem8)=nelem1
    htree(9,nelem8)=7_ip
    htree(10,nelem8)=-8_ip
    htree(11,nelem8)=0_ip 
    htree(12,nelem8)=htree(12,nelem1)

  end subroutine ref8

  subroutine hcoar(lmarkp,nelem,nnode,htree,elem,lmark,npoin,ntree,coor,ndim,lpoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nelem,npoin,nnode,ntree,ndim
    integer(ip),intent(in)    :: lmarkp(npoin)
    real(rp),intent(in)       :: coor(ndim,npoin)
    integer(ip),intent(inout) :: htree(ntree,nelem),elem(nnode,nelem),lmark(nelem)
    logical(lg),intent(inout) :: lpoin(npoin)
    integer(ip)               :: ncont,ielem,ltab(7,4),ltab3(4,6),ltab4(4,3)
    integer(ip)               :: nelem1,nelem2,nelem3,nelem4,nelem5,nelem6,nelem7
    integer(ip)               :: ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8,ip9,ip10
    integer(ip)               :: ltab2(6,6,6),nploc,lploc(10),lploc2(10),ncoun,nploc2,ipmark,iftyp
    integer(ip)               :: ip11,ip12,ip13,ip14,ip21,ip22,ip23,ip24
    integer(ip)               :: ip31,ip32,ip33,ip34,ip41,ip42,ip43,ip44
    integer(ip)               :: ip51,ip52,ip53,ip54,ip61,ip62,ip63,ip64
    integer(ip)               :: ip71,ip72,ip73,ip74,ip81,ip82,ip83,ip84
    integer(ip)               :: ipos,iel,npout,lpout(2,4),ipoin
    real(rp)                  :: c00,rvol
    c00=0.0d+00
    !
    !     Get the points of the face
    !
    ltab(1,1)=3
    ltab(2,1)=2
    ltab(3,1)=4
    ltab(4,1)=1
    ltab(5,1)=5
    ltab(6,1)=6
    ltab(7,1)=7

    ltab(1,2)=1
    ltab(2,2)=3
    ltab(3,2)=4
    ltab(4,2)=2
    ltab(5,2)=5
    ltab(6,2)=7
    ltab(7,2)=6

    ltab(1,3)=2
    ltab(2,3)=1
    ltab(3,3)=4
    ltab(4,3)=3
    ltab(5,3)=5
    ltab(6,3)=6
    ltab(7,3)=7

    ltab(1,4)=1
    ltab(2,4)=2
    ltab(3,4)=3
    ltab(4,4)=4
    ltab(5,4)=5
    ltab(6,4)=7
    ltab(7,4)=6
    !
    !     Obtain correct face
    !
    ltab2=0
    ltab2(1,2,4)=4
    ltab2(1,3,5)=3
    ltab2(2,3,6)=2  
    ltab2(4,5,6)=1
    !
    !     Obtain correct edge for 8:2 configuration
    !
    ltab3(1,1)=1
    ltab3(2,1)=2
    ltab3(3,1)=3
    ltab3(4,1)=4

    ltab3(1,2)=1
    ltab3(2,2)=3
    ltab3(3,2)=4
    ltab3(4,2)=2

    ltab3(1,3)=1
    ltab3(2,3)=4
    ltab3(3,3)=2
    ltab3(4,3)=3

    ltab3(1,4)=2
    ltab3(2,4)=3
    ltab3(3,4)=1
    ltab3(4,4)=4

    ltab3(1,5)=2
    ltab3(2,5)=4
    ltab3(3,5)=3
    ltab3(4,5)=1

    ltab3(1,6)=3
    ltab3(2,6)=4
    ltab3(3,6)=1
    ltab3(4,6)=2
    !
    !     Obtain correct edge for 4:2 configuration
    !
    ltab4(1,1)=1
    ltab4(2,1)=2
    ltab4(3,1)=3
    ltab4(4,1)=4

    ltab4(1,2)=2
    ltab4(2,2)=3
    ltab4(3,2)=1
    ltab4(4,2)=4

    ltab4(1,3)=3
    ltab4(2,3)=1
    ltab4(3,3)=2
    ltab4(4,3)=4
    !
    !     DBG
    !
    do ipoin=1,npoin
       if(lpoin(ipoin).eqv. .true.)then
          write(*,*)'Error, lpoin not clean hcoar'
          stop
       endif
    enddo
    !
    !     Loop on elements
    !
    do ielem=1,nelem
       !
       !     Has this element been marked to be coarsened?
       !  
       if(lmark(ielem)/=-1)cycle
       !
       !     Only consider parents
       !
       if(htree(10,ielem)>0)then
          !
          !     Which type do we have?
          !
          if(htree(10,ielem)==8)then
             !
             !     Case 8
             !  
             !
             !     First mark the inner points in lpoin
             !
             do ipos=4,7
                iel=htree(ipos,ielem)
                lpoin(elem(1,iel))=.true.
                lpoin(elem(2,iel))=.true.
                lpoin(elem(3,iel))=.true.
                lpoin(elem(4,iel))=.true.
             enddo
             !
             !     Then identify the outer points
             ! 
             npout=0_ip
             if(lpoin(elem(1,ielem)).eqv. .false.)then
                npout=npout+1_ip
                lpout(1,npout)=1
                lpout(2,npout)=elem(1,ielem)
             else if(lpoin(elem(2,ielem)).eqv. .false.)then
                npout=npout+1_ip
                lpout(1,npout)=2
                lpout(2,npout)=elem(2,ielem)
             else if(lpoin(elem(3,ielem)).eqv. .false.)then
                npout=npout+1_ip
                lpout(1,npout)=3
                lpout(2,npout)=elem(3,ielem)
             else
                npout=npout+1_ip
                lpout(1,npout)=4
                lpout(2,npout)=elem(4,ielem)
             endif

             do ipos=1,3
                iel=htree(ipos,ielem)  
                if(lpoin(elem(1,iel)).eqv. .false.)then
                   npout=npout+1_ip
                   lpout(1,npout)=1
                   lpout(2,npout)=elem(1,iel)
                else if(lpoin(elem(2,iel)).eqv. .false.)then
                   npout=npout+1_ip
                   lpout(1,npout)=2
                   lpout(2,npout)=elem(2,iel)
                else if(lpoin(elem(3,iel)).eqv. .false.)then
                   npout=npout+1_ip
                   lpout(1,npout)=3
                   lpout(2,npout)=elem(3,iel)
                else
                   npout=npout+1_ip
                   lpout(1,npout)=4
                   lpout(2,npout)=elem(4,iel)
                endif
             enddo
             !
             !     Clean up
             !  
             do ipos=4,7
                iel=htree(ipos,ielem)
                lpoin(elem(1,iel))=.false.
                lpoin(elem(2,iel))=.false.
                lpoin(elem(3,iel))=.false.
                lpoin(elem(4,iel))=.false.
             enddo
             !
             !     Get the sons
             ! 
             nelem1=htree(1,ielem)
             nelem2=htree(2,ielem)
             nelem3=htree(3,ielem)
             nelem4=htree(4,ielem)
             nelem5=htree(5,ielem)
             nelem6=htree(6,ielem)
             nelem7=htree(7,ielem)
             !
             !     Get the points
             ! 
             ip1=lpout(2,1) 
             ip2=lpout(2,2) 
             ip3=lpout(2,3) 
             ip4=lpout(2,4) 
             !
             !     Mark points of the father
             !
             lpoin(elem(1,ielem))=.true.
             lpoin(elem(2,ielem))=.true.
             lpoin(elem(3,ielem))=.true.
             lpoin(elem(4,ielem))=.true.
             !
             !   Get ip5
             !
             if(lpoin(elem(1,nelem1)).eqv. .true.)then
                ip5=elem(1,nelem1)
             else if(lpoin(elem(2,nelem1)).eqv. .true.)then
                ip5=elem(2,nelem1)
             else if(lpoin(elem(3,nelem1)).eqv. .true.)then
                ip5=elem(3,nelem1)
             else
                ip5=elem(4,nelem1)
             endif
             !
             !   Get ip6
             !
             if(lpoin(elem(1,nelem2)).eqv. .true.)then
                ip6=elem(1,nelem2)
             else if(lpoin(elem(2,nelem2)).eqv. .true.)then
                ip6=elem(2,nelem2)
             else if(lpoin(elem(3,nelem2)).eqv. .true.)then
                ip6=elem(3,nelem2)
             else
                ip6=elem(4,nelem2)
             endif
             !
             !   Get ip7
             !
             if(lpoin(elem(1,nelem3)).eqv. .true.)then
                ip7=elem(1,nelem3)
             else if(lpoin(elem(2,nelem3)).eqv. .true.)then
                ip7=elem(2,nelem3)
             else if(lpoin(elem(3,nelem3)).eqv. .true.)then
                ip7=elem(3,nelem3)
             else
                ip7=elem(4,nelem3)
             endif
             !
             !     Clean up
             !
             lpoin(elem(1,ielem))=.false.
             lpoin(elem(2,ielem))=.false.
             lpoin(elem(3,ielem))=.false.
             lpoin(elem(4,ielem))=.false.
             !
             !     Mark points of the nelem1
             !
             lpoin(elem(1,nelem1))=.true.
             lpoin(elem(2,nelem1))=.true.
             lpoin(elem(3,nelem1))=.true.
             lpoin(elem(4,nelem1))=.true.
             !
             !   Get ip8
             !
             if(lpoin(elem(1,nelem2)).eqv. .true.)then
                ip8=elem(1,nelem2)
             else if(lpoin(elem(2,nelem2)).eqv. .true.)then
                ip8=elem(2,nelem2)
             else if(lpoin(elem(3,nelem2)).eqv. .true.)then
                ip8=elem(3,nelem2)
             else
                ip8=elem(4,nelem2)
             endif
             !
             !   Get ip9
             !
             if(lpoin(elem(1,nelem3)).eqv. .true.)then
                ip9=elem(1,nelem3)
             else if(lpoin(elem(2,nelem3)).eqv. .true.)then
                ip9=elem(2,nelem3)
             else if(lpoin(elem(3,nelem3)).eqv. .true.)then
                ip9=elem(3,nelem3)
             else
                ip9=elem(4,nelem3)
             endif
             !
             !     Clean up
             !
             lpoin(elem(1,nelem1))=.false.
             lpoin(elem(2,nelem1))=.false.
             lpoin(elem(3,nelem1))=.false.
             lpoin(elem(4,nelem1))=.false.
             !
             !     Mark points of the nelem2
             !
             lpoin(elem(1,nelem2))=.true.
             lpoin(elem(2,nelem2))=.true.
             lpoin(elem(3,nelem2))=.true.
             lpoin(elem(4,nelem2))=.true.
             !
             !   Get ip10
             !
             if(lpoin(elem(1,nelem3)).eqv. .true.)then
                ip10=elem(1,nelem3)
             else if(lpoin(elem(2,nelem3)).eqv. .true.)then
                ip10=elem(2,nelem3)
             else if(lpoin(elem(3,nelem3)).eqv. .true.)then
                ip10=elem(3,nelem3)
             else
                ip10=elem(4,nelem3)
             endif
             !
             !     Clean up
             !
             lpoin(elem(1,nelem2))=.false.
             lpoin(elem(2,nelem2))=.false.
             lpoin(elem(3,nelem2))=.false.
             lpoin(elem(4,nelem2))=.false.
             !
             !     Count number of deleted points
             !
             ncoun=lmarkp(ip5)+lmarkp(ip6)+lmarkp(ip7)+lmarkp(ip8)+lmarkp(ip9)+lmarkp(ip10)

             if(ncoun==6)then
                !
                !     Case 8:1
                !
                elem(1,ielem)=ip1
                elem(2,ielem)=ip2
                elem(3,ielem)=ip3
                elem(4,ielem)=ip4

                call gtvol(elem,ielem,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 8:1'
                   stop  
                endif
                !
                !     Update htree
                !
                htree(1,ielem)=0_ip
                htree(2,ielem)=0_ip
                htree(3,ielem)=0_ip
                htree(4,ielem)=0_ip
                htree(5,ielem)=0_ip
                htree(6,ielem)=0_ip
                htree(7,ielem)=0_ip
                htree(10,ielem)=0_ip
                htree(11,ielem)=htree(11,ielem)-1
                htree(12,ielem)=htree(12,ielem)-1

                lmark(nelem1)=0_ip
                elem(1,nelem1)=0_ip
                htree(1,nelem1)=0_ip
                htree(2,nelem1)=0_ip
                htree(3,nelem1)=0_ip
                htree(4,nelem1)=0_ip
                htree(5,nelem1)=0_ip
                htree(6,nelem1)=0_ip
                htree(7,nelem1)=0_ip
                htree(8,nelem1)=0_ip
                htree(9,nelem1)=0_ip
                htree(10,nelem1)=0_ip
                htree(11,nelem1)=0_ip
                htree(12,nelem1)=0_ip

                lmark(nelem2)=0_ip
                elem(1,nelem2)=0_ip
                htree(1,nelem2)=0_ip
                htree(2,nelem2)=0_ip
                htree(3,nelem2)=0_ip
                htree(4,nelem2)=0_ip
                htree(5,nelem2)=0_ip
                htree(6,nelem2)=0_ip
                htree(7,nelem2)=0_ip
                htree(8,nelem2)=0_ip
                htree(9,nelem2)=0_ip
                htree(10,nelem2)=0_ip
                htree(11,nelem2)=0_ip
                htree(12,nelem2)=0_ip

                lmark(nelem3)=0_ip
                elem(1,nelem3)=0_ip
                htree(1,nelem3)=0_ip
                htree(2,nelem3)=0_ip
                htree(3,nelem3)=0_ip
                htree(4,nelem3)=0_ip
                htree(5,nelem3)=0_ip
                htree(6,nelem3)=0_ip
                htree(7,nelem3)=0_ip
                htree(8,nelem3)=0_ip
                htree(9,nelem3)=0_ip
                htree(10,nelem3)=0_ip
                htree(11,nelem3)=0_ip
                htree(12,nelem3)=0_ip

                lmark(nelem4)=0_ip
                elem(1,nelem4)=0_ip
                htree(1,nelem4)=0_ip
                htree(2,nelem4)=0_ip
                htree(3,nelem4)=0_ip
                htree(4,nelem4)=0_ip
                htree(5,nelem4)=0_ip
                htree(6,nelem4)=0_ip
                htree(7,nelem4)=0_ip
                htree(8,nelem4)=0_ip
                htree(9,nelem4)=0_ip
                htree(10,nelem4)=0_ip
                htree(11,nelem4)=0_ip
                htree(12,nelem4)=0_ip

                lmark(nelem5)=0_ip
                elem(1,nelem5)=0_ip
                htree(1,nelem5)=0_ip
                htree(2,nelem5)=0_ip
                htree(3,nelem5)=0_ip
                htree(4,nelem5)=0_ip
                htree(5,nelem5)=0_ip
                htree(6,nelem5)=0_ip
                htree(7,nelem5)=0_ip
                htree(8,nelem5)=0_ip
                htree(9,nelem5)=0_ip
                htree(10,nelem5)=0_ip
                htree(11,nelem5)=0_ip
                htree(12,nelem5)=0_ip

                lmark(nelem6)=0_ip
                elem(1,nelem6)=0_ip
                htree(1,nelem6)=0_ip
                htree(2,nelem6)=0_ip
                htree(3,nelem6)=0_ip
                htree(4,nelem6)=0_ip
                htree(5,nelem6)=0_ip
                htree(6,nelem6)=0_ip
                htree(7,nelem6)=0_ip
                htree(8,nelem6)=0_ip
                htree(9,nelem6)=0_ip
                htree(10,nelem6)=0_ip
                htree(11,nelem6)=0_ip
                htree(12,nelem6)=0_ip

                lmark(nelem7)=0_ip
                elem(1,nelem7)=0_ip
                htree(1,nelem7)=0_ip
                htree(2,nelem7)=0_ip
                htree(3,nelem7)=0_ip
                htree(4,nelem7)=0_ip
                htree(5,nelem7)=0_ip
                htree(6,nelem7)=0_ip
                htree(7,nelem7)=0_ip
                htree(8,nelem7)=0_ip
                htree(9,nelem7)=0_ip
                htree(10,nelem7)=0_ip
                htree(11,nelem7)=0_ip
                htree(12,nelem7)=0_ip

             else if(ncoun==5)then
                !
                !     Case 8:2
                !
                lploc(1)=ip1
                lploc(2)=ip2
                lploc(3)=ip3
                lploc(4)=ip4
                !
                !     Find the non marked point    
                !
                if(lmarkp(ip5)/=1)then
                   ipmark=1_ip
                   lploc(5)=ip5
                endif
                if(lmarkp(ip6)/=1)then
                   ipmark=2_ip
                   lploc(5)=ip6
                endif
                if(lmarkp(ip7)/=1)then
                   ipmark=3_ip
                   lploc(5)=ip7
                endif
                if(lmarkp(ip8)/=1)then
                   ipmark=4_ip
                   lploc(5)=ip8
                endif
                if(lmarkp(ip9)/=1)then
                   ipmark=5_ip
                   lploc(5)=ip9
                endif
                if(lmarkp(ip10)/=1)then
                   ipmark=6_ip
                   lploc(5)=ip10
                endif

                ip1=lploc(ltab3(1,ipmark))
                ip2=lploc(ltab3(2,ipmark))
                ip3=lploc(ltab3(3,ipmark))
                ip4=lploc(ltab3(4,ipmark))
                ip5=lploc(5)

                elem(1,ielem)=ip1
                elem(2,ielem)=ip5
                elem(3,ielem)=ip3
                elem(4,ielem)=ip4
                call gtvol(elem,ielem,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 8:2 a'
                   stop  
                endif

                elem(1,nelem1)=ip5
                elem(2,nelem1)=ip2
                elem(3,nelem1)=ip3
                elem(4,nelem1)=ip4
                call gtvol(elem,nelem1,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 8:2 b'
                   stop  
                endif
                !
                !     Update htree
                !
                htree(1,ielem)=nelem1
                htree(2,ielem)=0_ip
                htree(3,ielem)=0_ip
                htree(4,ielem)=0_ip
                htree(5,ielem)=0_ip
                htree(6,ielem)=0_ip
                htree(7,ielem)=0_ip
                htree(10,ielem)=2_ip
                htree(11,ielem)=htree(11,ielem)
                htree(12,ielem)=htree(12,ielem)

                lmark(nelem1)=0_ip
                htree(1,nelem1)=0_ip
                htree(2,nelem1)=0_ip
                htree(3,nelem1)=0_ip
                htree(4,nelem1)=0_ip
                htree(5,nelem1)=0_ip
                htree(6,nelem1)=0_ip
                htree(7,nelem1)=0_ip
                htree(8,nelem1)=ielem
                htree(9,nelem1)=1_ip
                htree(10,nelem1)=-2_ip
                htree(11,nelem1)=0_ip
                htree(12,nelem1)=htree(12,ielem)

                lmark(nelem2)=0_ip
                elem(1,nelem2)=0_ip
                htree(1,nelem2)=0_ip
                htree(2,nelem2)=0_ip
                htree(3,nelem2)=0_ip
                htree(4,nelem2)=0_ip
                htree(5,nelem2)=0_ip
                htree(6,nelem2)=0_ip
                htree(7,nelem2)=0_ip
                htree(8,nelem2)=0_ip
                htree(9,nelem2)=0_ip
                htree(10,nelem2)=0_ip
                htree(11,nelem2)=0_ip
                htree(12,nelem2)=0_ip

                lmark(nelem3)=0_ip
                elem(1,nelem3)=0_ip
                htree(1,nelem3)=0_ip
                htree(2,nelem3)=0_ip
                htree(3,nelem3)=0_ip
                htree(4,nelem3)=0_ip
                htree(5,nelem3)=0_ip
                htree(6,nelem3)=0_ip
                htree(7,nelem3)=0_ip
                htree(8,nelem3)=0_ip
                htree(9,nelem3)=0_ip
                htree(10,nelem3)=0_ip
                htree(11,nelem3)=0_ip
                htree(12,nelem3)=0_ip

                lmark(nelem4)=0_ip
                elem(1,nelem4)=0_ip
                htree(1,nelem4)=0_ip
                htree(2,nelem4)=0_ip
                htree(3,nelem4)=0_ip
                htree(4,nelem4)=0_ip
                htree(5,nelem4)=0_ip
                htree(6,nelem4)=0_ip
                htree(7,nelem4)=0_ip
                htree(8,nelem4)=0_ip
                htree(9,nelem4)=0_ip
                htree(10,nelem4)=0_ip
                htree(11,nelem4)=0_ip
                htree(12,nelem4)=0_ip

                lmark(nelem5)=0_ip
                elem(1,nelem5)=0_ip
                htree(1,nelem5)=0_ip
                htree(2,nelem5)=0_ip
                htree(3,nelem5)=0_ip
                htree(4,nelem5)=0_ip
                htree(5,nelem5)=0_ip
                htree(6,nelem5)=0_ip
                htree(7,nelem5)=0_ip
                htree(8,nelem5)=0_ip
                htree(9,nelem5)=0_ip
                htree(10,nelem5)=0_ip
                htree(11,nelem5)=0_ip
                htree(12,nelem5)=0_ip

                lmark(nelem6)=0_ip
                elem(1,nelem6)=0_ip
                htree(1,nelem6)=0_ip
                htree(2,nelem6)=0_ip
                htree(3,nelem6)=0_ip
                htree(4,nelem6)=0_ip
                htree(5,nelem6)=0_ip
                htree(6,nelem6)=0_ip
                htree(7,nelem6)=0_ip
                htree(8,nelem6)=0_ip
                htree(9,nelem6)=0_ip
                htree(10,nelem6)=0_ip
                htree(11,nelem6)=0_ip
                htree(12,nelem6)=0_ip

                lmark(nelem7)=0_ip
                elem(1,nelem7)=0_ip
                htree(1,nelem7)=0_ip
                htree(2,nelem7)=0_ip
                htree(3,nelem7)=0_ip
                htree(4,nelem7)=0_ip
                htree(5,nelem7)=0_ip
                htree(6,nelem7)=0_ip
                htree(7,nelem7)=0_ip
                htree(8,nelem7)=0_ip
                htree(9,nelem7)=0_ip
                htree(10,nelem7)=0_ip
                htree(11,nelem7)=0_ip
                htree(12,nelem7)=0_ip

             else if(ncoun==3)then
                !
                !     Case 8:4
                !
                lploc(1)=ip1
                lploc(2)=ip2
                lploc(3)=ip3
                lploc(4)=ip4
                nploc=4 
                !
                !     Find the non marked points    
                !
                nploc2=0
                if(lmarkp(ip5)==0)then
                   nploc=nploc+1
                   lploc(nploc)=ip5
                   nploc2=nploc2+1
                   lploc2(nploc2)=1_ip
                endif
                if(lmarkp(ip6)==0)then
                   nploc=nploc+1
                   lploc(nploc)=ip6
                   nploc2=nploc2+1
                   lploc2(nploc2)=2_ip
                endif
                if(lmarkp(ip7)==0)then
                   nploc=nploc+1
                   lploc(nploc)=ip7
                   nploc2=nploc2+1
                   lploc2(nploc2)=3_ip
                endif
                if(lmarkp(ip8)==0)then
                   nploc=nploc+1
                   lploc(nploc)=ip8
                   nploc2=nploc2+1
                   lploc2(nploc2)=4_ip
                endif
                if(lmarkp(ip9)==0)then
                   nploc=nploc+1
                   lploc(nploc)=ip9
                   nploc2=nploc2+1
                   lploc2(nploc2)=5_ip
                endif
                if(lmarkp(ip10)==0)then
                   nploc=nploc+1
                   lploc(nploc)=ip10
                   nploc2=nploc2+1
                   lploc2(nploc2)=6_ip
                endif

                iftyp=ltab2(lploc2(1),lploc2(2),lploc2(3))

                if(iftyp==0)then
                   write(*,*)'Error in hcoar iftyp' 
                   stop  
                endif

                ip1=lploc(ltab(1,iftyp))
                ip2=lploc(ltab(2,iftyp))
                ip3=lploc(ltab(3,iftyp))
                ip4=lploc(ltab(4,iftyp))
                ip5=lploc(ltab(5,iftyp))
                ip6=lploc(ltab(6,iftyp))
                ip7=lploc(ltab(7,iftyp))

                elem(1,ielem)=ip1
                elem(2,ielem)=ip5
                elem(3,ielem)=ip7
                elem(4,ielem)=ip4
                call gtvol(elem,ielem,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 8:4 a'
                   stop  
                endif

                elem(1,nelem1)=ip5
                elem(2,nelem1)=ip2
                elem(3,nelem1)=ip6
                elem(4,nelem1)=ip4
                call gtvol(elem,nelem1,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 8:4 b'
                   stop  
                endif

                elem(1,nelem2)=ip6
                elem(2,nelem2)=ip3
                elem(3,nelem2)=ip7
                elem(4,nelem2)=ip4
                call gtvol(elem,nelem2,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 8:4 c'
                   stop  
                endif

                elem(1,nelem3)=ip5
                elem(2,nelem3)=ip6
                elem(3,nelem3)=ip7
                elem(4,nelem3)=ip4
                call gtvol(elem,nelem3,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 8:4 c'
                   stop  
                endif
                !
                !     Update htree
                !
                htree(1,ielem)=nelem1
                htree(2,ielem)=nelem2
                htree(3,ielem)=nelem3
                htree(4,ielem)=0_ip
                htree(5,ielem)=0_ip
                htree(6,ielem)=0_ip
                htree(7,ielem)=0_ip
                htree(10,ielem)=4_ip
                htree(11,ielem)=htree(11,ielem)
                htree(12,ielem)=htree(12,ielem)

                lmark(nelem1)=0_ip
                htree(1,nelem1)=0_ip
                htree(2,nelem1)=0_ip
                htree(3,nelem1)=0_ip
                htree(4,nelem1)=0_ip
                htree(5,nelem1)=0_ip
                htree(6,nelem1)=0_ip
                htree(7,nelem1)=0_ip
                htree(8,nelem1)=ielem
                htree(9,nelem1)=1_ip
                htree(10,nelem1)=-4_ip
                htree(11,nelem1)=0_ip
                htree(12,nelem1)=htree(12,ielem)

                lmark(nelem2)=0_ip
                htree(1,nelem2)=0_ip
                htree(2,nelem2)=0_ip
                htree(3,nelem2)=0_ip
                htree(4,nelem2)=0_ip
                htree(5,nelem2)=0_ip
                htree(6,nelem2)=0_ip
                htree(7,nelem2)=0_ip
                htree(8,nelem2)=ielem
                htree(9,nelem2)=2_ip
                htree(10,nelem2)=-4_ip
                htree(11,nelem2)=0_ip
                htree(12,nelem2)=htree(12,ielem)

                lmark(nelem3)=0_ip
                htree(1,nelem3)=0_ip
                htree(2,nelem3)=0_ip
                htree(3,nelem3)=0_ip
                htree(4,nelem3)=0_ip
                htree(5,nelem3)=0_ip
                htree(6,nelem3)=0_ip
                htree(7,nelem3)=0_ip
                htree(8,nelem3)=ielem
                htree(9,nelem3)=3_ip
                htree(10,nelem3)=-4_ip
                htree(11,nelem3)=0_ip
                htree(12,nelem3)=htree(12,ielem)

                lmark(nelem4)=0_ip
                elem(1,nelem4)=0_ip
                htree(1,nelem4)=0_ip
                htree(2,nelem4)=0_ip
                htree(3,nelem4)=0_ip
                htree(4,nelem4)=0_ip
                htree(5,nelem4)=0_ip
                htree(6,nelem4)=0_ip
                htree(7,nelem4)=0_ip
                htree(8,nelem4)=0_ip
                htree(9,nelem4)=0_ip
                htree(10,nelem4)=0_ip
                htree(11,nelem4)=0_ip
                htree(12,nelem4)=0_ip

                lmark(nelem5)=0_ip
                elem(1,nelem5)=0_ip
                htree(1,nelem5)=0_ip
                htree(2,nelem5)=0_ip
                htree(3,nelem5)=0_ip
                htree(4,nelem5)=0_ip
                htree(5,nelem5)=0_ip
                htree(6,nelem5)=0_ip
                htree(7,nelem5)=0_ip
                htree(8,nelem5)=0_ip
                htree(9,nelem5)=0_ip
                htree(10,nelem5)=0_ip
                htree(11,nelem5)=0_ip
                htree(12,nelem5)=0_ip

                lmark(nelem6)=0_ip
                elem(1,nelem6)=0_ip
                htree(1,nelem6)=0_ip
                htree(2,nelem6)=0_ip
                htree(3,nelem6)=0_ip
                htree(4,nelem6)=0_ip
                htree(5,nelem6)=0_ip
                htree(6,nelem6)=0_ip
                htree(7,nelem6)=0_ip
                htree(8,nelem6)=0_ip
                htree(9,nelem6)=0_ip
                htree(10,nelem6)=0_ip
                htree(11,nelem6)=0_ip
                htree(12,nelem6)=0_ip

                lmark(nelem7)=0_ip
                elem(1,nelem7)=0_ip
                htree(1,nelem7)=0_ip
                htree(2,nelem7)=0_ip
                htree(3,nelem7)=0_ip
                htree(4,nelem7)=0_ip
                htree(5,nelem7)=0_ip
                htree(6,nelem7)=0_ip
                htree(7,nelem7)=0_ip
                htree(8,nelem7)=0_ip
                htree(9,nelem7)=0_ip
                htree(10,nelem7)=0_ip
                htree(11,nelem7)=0_ip
                htree(12,nelem7)=0_ip

             else 

                write(*,*)'Error in hcoar, ncoun=',ncoun
                stop

             endif

          else if(htree(10,ielem)==4)then
             !
             !    Case 4
             !
             nelem1=htree(1,ielem)
             nelem2=htree(2,ielem)
             nelem3=htree(3,ielem)

             ip1=elem(1,ielem) 
             ip2=elem(2,nelem1)
             ip3=elem(2,nelem2)
             ip4=elem(4,nelem3)

             ip5=elem(2,ielem)
             ip6=elem(3,nelem1)
             ip7=elem(3,ielem)

             ncoun=lmarkp(ip5)+lmarkp(ip6)+lmarkp(ip7)

             if(ncoun==3)then
                !
                !    Case 4:1
                !
                elem(1,ielem)=ip1
                elem(2,ielem)=ip2
                elem(3,ielem)=ip3
                elem(4,ielem)=ip4
                call gtvol(elem,ielem,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 4:1'
                   stop  
                endif
                !
                !     Update htree
                !
                htree(1,ielem)=0_ip
                htree(2,ielem)=0_ip
                htree(3,ielem)=0_ip
                htree(4,ielem)=0_ip
                htree(5,ielem)=0_ip
                htree(6,ielem)=0_ip
                htree(7,ielem)=0_ip
                htree(10,ielem)=0_ip
                htree(11,ielem)=htree(11,ielem)-1
                htree(12,ielem)=htree(12,ielem)-1

                lmark(nelem1)=0_ip
                elem(1,nelem1)=0_ip
                htree(1,nelem1)=0_ip
                htree(2,nelem1)=0_ip
                htree(3,nelem1)=0_ip
                htree(4,nelem1)=0_ip
                htree(5,nelem1)=0_ip
                htree(6,nelem1)=0_ip
                htree(7,nelem1)=0_ip
                htree(8,nelem1)=0_ip
                htree(9,nelem1)=0_ip
                htree(10,nelem1)=0_ip
                htree(11,nelem1)=0_ip
                htree(12,nelem1)=0_ip

                lmark(nelem2)=0_ip
                elem(1,nelem2)=0_ip
                htree(1,nelem2)=0_ip
                htree(2,nelem2)=0_ip
                htree(3,nelem2)=0_ip
                htree(4,nelem2)=0_ip
                htree(5,nelem2)=0_ip
                htree(6,nelem2)=0_ip
                htree(7,nelem2)=0_ip
                htree(8,nelem2)=0_ip
                htree(9,nelem2)=0_ip
                htree(10,nelem2)=0_ip
                htree(11,nelem2)=0_ip
                htree(12,nelem2)=0_ip

                lmark(nelem3)=0_ip
                elem(1,nelem3)=0_ip
                htree(1,nelem3)=0_ip
                htree(2,nelem3)=0_ip
                htree(3,nelem3)=0_ip
                htree(4,nelem3)=0_ip
                htree(5,nelem3)=0_ip
                htree(6,nelem3)=0_ip
                htree(7,nelem3)=0_ip
                htree(8,nelem3)=0_ip
                htree(9,nelem3)=0_ip
                htree(10,nelem3)=0_ip
                htree(11,nelem3)=0_ip
                htree(12,nelem3)=0_ip

             else if(ncoun==2)then
                !
                !    Case 4:2
                !
                lploc(1)=ip1
                lploc(2)=ip2
                lploc(3)=ip3
                lploc(4)=ip4
                !
                !     Find the non marked point    
                !
                if(lmarkp(ip5)==0)then
                   ipmark=1_ip
                   lploc(5)=ip5
                endif
                if(lmarkp(ip6)==0)then
                   ipmark=2_ip
                   lploc(5)=ip6
                endif
                if(lmarkp(ip7)==0)then
                   ipmark=3_ip
                   lploc(5)=ip7
                endif

                ip1=lploc(ltab4(1,ipmark))
                ip2=lploc(ltab4(2,ipmark))
                ip3=lploc(ltab4(3,ipmark))
                ip4=lploc(ltab4(4,ipmark))
                ip5=lploc(5)

                elem(1,ielem)=ip1
                elem(2,ielem)=ip5
                elem(3,ielem)=ip3
                elem(4,ielem)=ip4
                call gtvol(elem,ielem,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 4:2 a'
                   stop  
                endif

                elem(1,nelem1)=ip5
                elem(2,nelem1)=ip2
                elem(3,nelem1)=ip3
                elem(4,nelem1)=ip4
                call gtvol(elem,nelem1,coor,rvol,nnode,nelem,npoin,ndim)
                if(rvol<=c00)then
                   write(*,*)'Negative volume 4:2 b'
                   stop  
                endif
                !
                !     Update htree
                !
                htree(1,ielem)=nelem1
                htree(2,ielem)=0_ip
                htree(3,ielem)=0_ip
                htree(4,ielem)=0_ip
                htree(5,ielem)=0_ip
                htree(6,ielem)=0_ip
                htree(7,ielem)=0_ip
                htree(10,ielem)=2_ip
                htree(11,ielem)=htree(11,ielem)
                htree(12,ielem)=htree(12,ielem)

                lmark(nelem1)=0_ip
                htree(1,nelem1)=0_ip
                htree(2,nelem1)=0_ip
                htree(3,nelem1)=0_ip
                htree(4,nelem1)=0_ip
                htree(5,nelem1)=0_ip
                htree(6,nelem1)=0_ip
                htree(7,nelem1)=0_ip
                htree(8,nelem1)=ielem
                htree(9,nelem1)=1_ip
                htree(10,nelem1)=-2_ip
                htree(11,nelem1)=0_ip
                htree(12,nelem1)=htree(12,ielem)

                lmark(nelem2)=0_ip
                elem(1,nelem2)=0_ip
                htree(1,nelem2)=0_ip
                htree(2,nelem2)=0_ip
                htree(3,nelem2)=0_ip
                htree(4,nelem2)=0_ip
                htree(5,nelem2)=0_ip
                htree(6,nelem2)=0_ip
                htree(7,nelem2)=0_ip
                htree(8,nelem2)=0_ip
                htree(9,nelem2)=0_ip
                htree(10,nelem2)=0_ip
                htree(11,nelem2)=0_ip
                htree(12,nelem2)=0_ip

                lmark(nelem3)=0_ip
                elem(1,nelem3)=0_ip
                htree(1,nelem3)=0_ip
                htree(2,nelem3)=0_ip
                htree(3,nelem3)=0_ip
                htree(4,nelem3)=0_ip
                htree(5,nelem3)=0_ip
                htree(6,nelem3)=0_ip
                htree(7,nelem3)=0_ip
                htree(8,nelem3)=0_ip
                htree(9,nelem3)=0_ip
                htree(10,nelem3)=0_ip
                htree(11,nelem3)=0_ip
                htree(12,nelem3)=0_ip

             else

                write(*,*)'Error coarsening 4:1'  
                stop

             endif

          else if(htree(10,ielem)==2)then
             !
             !     Case 2:1
             !
             nelem1=htree(1,ielem)

             ip5=elem(2,ielem)

             if(lmarkp(ip5)/=1)then
                write(*,*)'Error coarsening 2:1'
                stop
             endif

             elem(2,ielem)=elem(2,nelem1)
             call gtvol(elem,ielem,coor,rvol,nnode,nelem,npoin,ndim)
             if(rvol<=c00)then
                write(*,*)'Negative volume 2:1'
                stop  
             endif
             !
             !     Update htree
             !
             htree(1,ielem)=0_ip
             htree(2,ielem)=0_ip
             htree(3,ielem)=0_ip
             htree(4,ielem)=0_ip
             htree(5,ielem)=0_ip
             htree(6,ielem)=0_ip
             htree(7,ielem)=0_ip
             htree(10,ielem)=0_ip
             htree(11,ielem)=htree(11,ielem)-1
             htree(12,ielem)=htree(12,ielem)-1

             lmark(nelem1)=0_ip
             elem(1,nelem1)=0_ip
             htree(1,nelem1)=0_ip
             htree(2,nelem1)=0_ip
             htree(3,nelem1)=0_ip
             htree(4,nelem1)=0_ip
             htree(5,nelem1)=0_ip
             htree(6,nelem1)=0_ip
             htree(7,nelem1)=0_ip
             htree(8,nelem1)=0_ip
             htree(9,nelem1)=0_ip
             htree(10,nelem1)=0_ip
             htree(11,nelem1)=0_ip
             htree(12,nelem1)=0_ip

          endif

       endif

    enddo

  end subroutine hcoar

  subroutine reorder(elem,nnode,nelem,htree,lmark,ntree)
    use mod_memchk
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only        : memor_msh
    implicit none
    integer(ip),intent(in)    :: nnode,nelem,ntree
    integer(ip),intent(inout) :: elem(nnode,nelem),htree(ntree,nelem),lmark(nelem)
    integer(ip)               :: ielem,level(500),ilevel,ienew,ilevelmax
    integer(ip), pointer      :: elemt(:,:),renum(:)
    integer(4)                :: istat
    !
    !     This subroutine reorders the elements with respect to their level in the refinement
    !
    !
    !     Allocate help arrays
    !  
    allocate(elemt(nnode,nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'ELEMT','reorder',elemt)
    allocate(renum(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'RENUM','reorder',renum)
    !
    !     First count the elements in each level
    !
    level=0
    ilevelmax=-1
    do ielem=1,nelem
       ilevel=htree(12,ielem)+2
       if(ilevel>ilevelmax)ilevelmax=ilevel
       level(ilevel)=level(ilevel)+1
    enddo
    !
    !     Sum up
    !
    level(1)=1_ip
    do ilevel=2,ilevelmax
       level(ilevel)=level(ilevel)+level(ilevel-1)
    enddo
    !
    !     Get the reordering
    !
    do ielem=1,nelem
       ilevel=htree(12,ielem)+1
       renum(ielem)=level(ilevel)
       level(ilevel)=level(ilevel)+1
    enddo
    !
    !     Reorder everything
    !
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,elem,renum) &
    !$omp& private(ielem,ienew)

    do ielem=1,nelem
       ienew=renum(ielem) 
       if(ienew==0)then
          write(*,*)'Error reordering ielem=',ielem
          stop
       endif
       elemt(1,ienew)=elem(1,ielem)
       elemt(2,ienew)=elem(2,ielem)
       elemt(3,ienew)=elem(3,ielem)
       elemt(4,ienew)=elem(4,ielem)
    enddo
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,elem) &
    !$omp& private(ielem,ienew)

    do ielem=1,nelem
       elem(1,ielem)=elemt(1,ielem)
       elem(2,ielem)=elemt(2,ielem)
       elem(3,ielem)=elemt(3,ielem)
       elem(4,ielem)=elemt(4,ielem)
    enddo
    !
    !     Renum htree
    !
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,htree,renum) &
    !$omp& private(ielem,ienew)

    do ielem=1,nelem
       ienew=renum(ielem)
       if(htree(1,ielem)/=0)then 
          elemt(1,ienew)=renum(htree(1,ielem))
       else
          elemt(1,ienew)=0_ip
       endif
       if(htree(2,ielem)/=0)then 
          elemt(2,ienew)=renum(htree(2,ielem))
       else
          elemt(2,ienew)=0_ip
       endif
       if(htree(3,ielem)/=0)then 
          elemt(3,ienew)=renum(htree(3,ielem))
       else
          elemt(3,ienew)=0_ip
       endif
       if(htree(4,ielem)/=0)then 
          elemt(4,ienew)=renum(htree(4,ielem))
       else
          elemt(4,ienew)=0_ip
       endif
    enddo
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,htree) &
    !$omp& private(ielem)

    do ielem=1,nelem
       htree(1,ielem)=elemt(1,ielem)
       htree(2,ielem)=elemt(2,ielem)
       htree(3,ielem)=elemt(3,ielem)
       htree(4,ielem)=elemt(4,ielem)
    enddo
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,htree,renum) &
    !$omp& private(ielem,ienew)

    do ielem=1,nelem
       ienew=renum(ielem) 
       if(htree(5,ielem)/=0)then 
          elemt(1,ienew)=renum(htree(5,ielem))
       else
          elemt(1,ienew)=0_ip
       endif
       if(htree(6,ielem)/=0)then 
          elemt(2,ienew)=renum(htree(6,ielem))
       else
          elemt(2,ienew)=0_ip
       endif
       if(htree(7,ielem)/=0)then 
          elemt(3,ienew)=renum(htree(7,ielem))
       else
          elemt(3,ienew)=0_ip
       endif
       if(htree(8,ielem)/=0)then 
          elemt(4,ienew)=renum(htree(8,ielem))
       else
          elemt(4,ienew)=0_ip
       endif
    enddo
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,htree) &
    !$omp& private(ielem)

    do ielem=1,nelem
       htree(5,ielem)=elemt(1,ielem)
       htree(6,ielem)=elemt(2,ielem)
       htree(7,ielem)=elemt(3,ielem)
       htree(8,ielem)=elemt(4,ielem)
    enddo
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,htree,renum) &
    !$omp& private(ielem,ienew)

    do ielem=1,nelem
       ienew=renum(ielem) 
       elemt(1,ienew)=htree(9,ielem)
       elemt(2,ienew)=htree(10,ielem)
       elemt(3,ienew)=htree(11,ielem)
       elemt(4,ienew)=htree(12,ielem)
    enddo
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,htree) &
    !$omp& private(ielem)

    do ielem=1,nelem
       htree(9,ielem)=elemt(1,ielem)
       htree(10,ielem)=elemt(2,ielem)
       htree(11,ielem)=elemt(3,ielem)
       htree(12,ielem)=elemt(4,ielem)
    enddo
    !
    !     Renum lmark
    !
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,lmark,renum) &
    !$omp& private(ielem,ienew)

    do ielem=1,nelem
       ienew=renum(ielem) 
       elemt(1,ienew)=lmark(ielem)
    enddo
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elemt,lmark) &
    !$omp& private(ielem)

    do ielem=1,nelem
       lmark(ielem)=elemt(1,ielem)
    enddo

    call memchk(2_ip,istat,memor_msh,'RENUM','reorder',renum)
    deallocate(renum,stat=istat)
    if(istat/=0) call memerr(2_ip,'RENUM','reorder',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELEMT','reorder',elemt)
    deallocate(elemt,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELEMT','reorder',0_ip)

  end subroutine reorder

  subroutine reorderd(elem,nnode,nelem,htree,lmarkp,coor,ndim,npoin,npoin0,&
       var,ntree,nvar,nface,nnofa,lface)
    use mod_memchk
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)    :: nnode,ndim,npoin0,ntree,nvar,nnofa,nface
    integer(ip),intent(inout) :: nelem,npoin
    integer(ip),intent(in)    :: lmarkp(npoin)
    integer(ip),intent(inout) :: lface(nnofa,nface)
    integer(ip),intent(inout) :: elem(nnode,nelem),htree(ntree,nelem)
    real(rp),intent(inout)    :: coor(ndim,npoin)
    real(rp),pointer          :: var(:,:)
    integer(ip)               :: ielem,ienew,ilevelmax,nelem0,npoin1,ipoin,ivar
    integer(ip)               :: iface
    real(rp)                  :: c10  
    integer(ip),pointer       :: renum(:) 
    integer(4)                :: istat

    allocate(renum(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'RENUM','reorder',renum)
    c10=1.0d+00
    !
    !     Get the reordering
    !
    nelem0=nelem 
    nelem=0_ip
    do ielem=1,nelem0
       if(elem(1,ielem)/=0)then
          nelem=nelem+1
          renum(ielem)=nelem
          elem(1,nelem)=elem(1,ielem)
          elem(2,nelem)=elem(2,ielem)
          elem(3,nelem)=elem(3,ielem)
          elem(4,nelem)=elem(4,ielem)
          htree(1,nelem)=htree(1,ielem) 
          htree(2,nelem)=htree(2,ielem) 
          htree(3,nelem)=htree(3,ielem) 
          htree(4,nelem)=htree(4,ielem) 
          htree(5,nelem)=htree(5,ielem) 
          htree(6,nelem)=htree(6,ielem) 
          htree(7,nelem)=htree(7,ielem) 
          htree(8,nelem)=htree(8,ielem) 
          htree(9,nelem)=htree(9,ielem) 
          htree(10,nelem)=htree(10,ielem) 
          htree(11,nelem)=htree(11,ielem) 
          htree(12,nelem)=htree(12,ielem) 
       endif
    enddo
    !
    !    Update htree
    !    Renumber sons if necessary and father
    !
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,htree,renum) &
    !$omp& private(ielem)

    do ielem=1,nelem
       if(htree(1,ielem)/=0)then 
          htree(1,ielem)=renum(htree(1,ielem))
       endif
       if(htree(2,ielem)/=0)then 
          htree(2,ielem)=renum(htree(2,ielem))
       endif
       if(htree(3,ielem)/=0)then 
          htree(3,ielem)=renum(htree(3,ielem))
       endif
       if(htree(4,ielem)/=0)then 
          htree(4,ielem)=renum(htree(4,ielem))
       endif
       if(htree(5,ielem)/=0)then 
          htree(5,ielem)=renum(htree(5,ielem))
       endif
       if(htree(6,ielem)/=0)then 
          htree(6,ielem)=renum(htree(6,ielem))
       endif
       if(htree(7,ielem)/=0)then 
          htree(7,ielem)=renum(htree(7,ielem))
       endif
       if(htree(8,ielem)/=0)then 
          htree(8,ielem)=renum(htree(8,ielem))
       endif
    enddo
    !
    !     Reorder the coordinates and the unknowns until newly created points
    !
    renum=0
    npoin1=npoin 
    npoin=0
    do ipoin=1,npoin0
       if(lmarkp(ipoin)==0)then 
          npoin=npoin+1
          renum(ipoin)=npoin
          coor(1,npoin)=coor(1,ipoin)
          coor(2,npoin)=coor(2,ipoin)
          coor(3,npoin)=coor(3,ipoin)
          do ivar=1,nvar
             var(ivar,npoin)=var(ivar,ipoin)
          enddo
       endif
    enddo
    !
    !     Reorder the coordinates and the unknowns remaining
    !
    do ipoin=npoin0+1,npoin1
       npoin=npoin+1
       renum(ipoin)=npoin
       coor(1,npoin)=coor(1,ipoin)
       coor(2,npoin)=coor(2,ipoin)
       coor(3,npoin)=coor(3,ipoin)
       do ivar=1,nvar
          var(ivar,npoin)=var(ivar,ipoin)
       enddo
    enddo
    !
    !     Reorder the points in the elements
    ! 
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,elem,renum) &
    !$omp& private(ielem)

    do ielem=1,nelem
       elem(1,ielem)=renum(elem(1,ielem))
       elem(2,ielem)=renum(elem(2,ielem))
       elem(3,ielem)=renum(elem(3,ielem))
       elem(4,ielem)=renum(elem(4,ielem))
    enddo
    !
    !     Renumber the faces 
    !
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nface,lface,renum) &
    !$omp& private(iface)

    do iface=1,nface
       lface(1,iface)=renum(lface(1,iface))
       lface(2,iface)=renum(lface(2,iface))
       lface(3,iface)=renum(lface(3,iface))
    enddo

    call memchk(2_ip,istat,memor_msh,'RENUM','reorder',renum)
    deallocate(renum,stat=istat)
    if(istat/=0) call memerr(2_ip,'RENUM','reorder',0_ip)

  end subroutine reorderd

  subroutine outhref(nnode,nelem,elem,ndim,npoin,coor,lmark)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nnode,nelem,npoin,ndim
    integer(ip),intent(in)      :: elem(nnode,nelem),lmark(nelem)
    real(rp), intent(in)        :: coor(ndim,npoin)
    integer(ip)                 :: ielem,ipoin       
    real(rp)                    :: rx,ry,rz 
    !
    !     Print mesh
    !
    open(unit=50,file='outhref.msh',status='unknown')
    rewind 50

1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coord_x   coord_y  coord_z')
100 format(i10,3e20.10)
200 format(6i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(50,1)
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
    do  ielem=1,nelem
       write(50,200)ielem,elem(1,ielem),elem(2,ielem),elem(3,ielem),elem(4,ielem),lmark(ielem)+2
    enddo
    write(50,6)
    close(50)

  end subroutine outhref

  subroutine outhref2(nnode,nelem,elem,ndim,npoin,coor)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nnode,nelem,npoin,ndim
    integer(ip),intent(in)      :: elem(nnode,nelem)
    real(rp), intent(in)        :: coor(ndim,npoin)
    integer(ip)                 :: ielem,ipoin       
    real(rp)                    :: rx,ry,rz 
    !
    !     Print mesh
    !
    open(unit=50,file='outhref2.msh',status='unknown')
    rewind 50

1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coord_x   coord_y  coord_z')
100 format(i10,3e20.10)
200 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(50,1)
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
    do  ielem=1,nelem
       write(50,200)ielem,elem(1,ielem),elem(2,ielem),elem(3,ielem),elem(4,ielem)
    enddo
    write(50,6)
    close(50)

  end subroutine outhref2

  subroutine dbghref(htree,elem,nelem,nnode,ntree)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nnode,nelem,ntree
    integer(ip),intent(in)      :: elem(nnode,nelem),htree(ntree,nelem)

  end subroutine dbghref

  subroutine dbghref2(eltoel,nnode,nelem,elem,htree,ntree)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nnode,nelem,ntree
    integer(ip),intent(in)      :: elem(nnode,nelem),htree(ntree,nelem),eltoel(nnode,nelem)

  end subroutine dbghref2

  subroutine corbou(nelem,nnode,eltoel,nsid,ledglm,lsmark,nedge,htree,ntree)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nelem,nnode,nsid,nedge,ntree 
    integer(ip),intent(in)    :: eltoel(nnode,nelem)
    integer(ip),intent(in)    :: ledglm(nsid,nelem),htree(ntree,nelem)  
    integer(ip),intent(inout) :: lsmark(nedge)
    integer(ip)               :: ielem 
    !
    !     This subroutine corrects the 8:4 that have been marked with only two sides  
    !     due to the boundary 
    ! 
    do ielem=1,nelem
       if(htree(10,ielem)==4)then  
          !
          !     Must be on the boundary
          !
          if(eltoel(4,ielem)==0)then
             !
             !     Must have two sides marked   
             !
             if(lsmark(ledglm(1,ielem))==1 .and. lsmark(ledglm(2,ielem))==1)then
                lsmark(ledglm(4,ielem))=1
             endif
          endif

       else if(htree(10,ielem)==-4)then  
          !
          !     Must be on the boundary
          !
          if(eltoel(4,ielem)==0)then
             !
             !     Must not be the third son
             ! 
             if(htree(9,ielem)==3)cycle
             !
             !     Must have two sides marked   
             !
             if(lsmark(ledglm(1,ielem))==1 .and. lsmark(ledglm(4,ielem))==1)then
                lsmark(ledglm(2,ielem))=1
             endif
          endif

       endif

    enddo

  end subroutine corbou

  subroutine update(nelem,htree,elem,nnode,ntree,hftree,nftree,nface)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nelem,nnode,ntree,nftree,nface
    integer(ip),intent(in)      :: elem(nnode,nelem)
    integer(ip),intent(inout)   :: htree(ntree,nelem),hftree(nftree,nface)
    integer(ip)                 :: nchil,ichil,ityp,ielem,ilevel,jlevel,ipos,ifath,jelem
    integer(ip)                 :: ison,iface,jface,je
    !
    !     This sub reconstructs htree deleted during the coarsening
    !

    !
    !     First send the information from the sons to the fathers
    !
    do ielem=nelem,1,-1
       !
       !     Find the coarsened elements that have lost their parents/sons
       !
       !
       !     The sons must have been at the very least divided
       ! 
       if(htree(12,ielem)==0)cycle
       !
       !     Get global level
       !
       ilevel=htree(12,ielem)
       !
       !     Get father
       ! 
       ifath=htree(8,ielem)
       !
       !     Do not consider coarsest elements
       !
       if(ifath==0)cycle
       !
       !     The global level between the son and the father must be compatible
       !     The father can not have a deeper level than the son
       !
       if(htree(12,ifath)>ilevel)cycle      
       !
       !     The parent must have been coarsened 
       !
       if(htree(10,ifath)==0)then
          ipos=htree(9,ielem)
          if(ipos>7)then
             write(*,*)'Error ipos'
             stop
          endif

          if(htree(11,ifath)==0)then
             write(*,*)'Strange case update'
             stop
          endif
          !
          !     Check if already written
          !
          if(htree(ipos,ifath)==0)then
             htree(ipos,ifath)=ielem
          else
             !
             !     Take the one with the lowest local level
             !
             jelem=htree(ipos,ifath)
             if(htree(11,ielem)<htree(11,jelem))then
                htree(ipos,ifath)=ielem
             else if(htree(11,ielem)==htree(11,jelem))then
                if(htree(12,ielem)>htree(12,jelem))then
                   htree(ipos,ifath)=ielem
                endif 
             endif
          endif
       endif
    enddo
    !
    !     Then decide the type
    !
    do ielem=1,nelem
       !
       !     The element should have been coarsened
       !     and not be at the coarsest level  
       !
       if(htree(10,ielem)==0 .and. htree(12,ielem)>0)then
          do ichil=1,7
             if(htree(ichil,ielem)==0)exit
          enddo
          !
          !     Test the child number
          !
          if(ichil==1)then
             !
             !     No children, It must be a -8
             ! 
             htree(10,ielem)=-8  
             cycle

          else if(ichil==2)then 
             !
             !     Two children, end of the tree
             !
             htree(10,ielem)=2_ip
             htree(10,htree(1,ielem))=-2_ip

          else if(ichil==4)then
             !
             !     Four children, end of the tree
             !
             htree(10,ielem)=4_ip
             htree(10,htree(1,ielem))=-4_ip
             htree(10,htree(2,ielem))=-4_ip
             htree(10,htree(3,ielem))=-4_ip

          else if(ichil==8)then
             !
             !     Eight children, must check for children
             !     of children if refined 
             !
             htree(10,ielem)=8_ip
             ison=htree(1,ielem)
             if(htree(1,ison)==0)then
                htree(10,ison)=-8_ip
             endif
             ison=htree(2,ielem)
             if(htree(1,ison)==0)then
                htree(10,ison)=-8_ip
             endif
             ison=htree(3,ielem)
             if(htree(1,ison)==0)then
                htree(10,ison)=-8_ip
             endif
             ison=htree(4,ielem)
             if(htree(1,ison)==0)then
                htree(10,ison)=-8_ip
             endif
             ison=htree(5,ielem)
             if(htree(1,ison)==0)then
                htree(10,ison)=-8_ip
             endif
             ison=htree(6,ielem)
             if(htree(1,ison)==0)then
                htree(10,ison)=-8_ip
             endif
             ison=htree(7,ielem)
             if(htree(1,ison)==0)then
                htree(10,ison)=-8_ip
             endif
          else
             write(*,*)'Error in update htree'
             write(*,*)'ichil=',ichil
             stop 
          endif
       !else if(htree(10,ielem)>0)then
          !
          !     If the sons have coarsened, remark them
          !   
       !   nchil=htree(10,ielem)-1
       !   do je=1,nchil
       !      ichil=htree(je,ielem) 
       !      if(htree(10,ichil)==0)then
       !         htree(10,ichil)=-htree(10,ielem)
       !      endif
       !   enddo

       endif
    enddo
    !
    !     Now update faces
    !
    !
    !     First send the information from the sons to the fathers
    !
    do iface=nface,1,-1
       !
       !     Find the coarsened elements that have lost their parents/sons
       !
       !write(*,*)'In update iface=',iface
       !
       !     The sons must have been at the very least divided
       ! 
       if(hftree(8,iface)==0)cycle
       !
       !     Get global level
       !
       ilevel=hftree(8,iface)
       !
       !     Get father
       ! 
       ifath=hftree(4,iface)
       !
       !     Do not consider coarsest elements
       !
       if(ifath==0)cycle
       !
       !     The global level between the son and the father must be compatible
       !     The father can not have a deeper level than the son
       !
       if(hftree(8,ifath)>ilevel)cycle      
       !
       !     The parent must have been coarsened 
       !
       if(hftree(6,ifath)==0)then
          ipos=hftree(5,iface)
          if(ipos>3)then
             write(*,*)'Error ipos'
             stop
          endif

          if(hftree(7,ifath)==0)then
             write(*,*)'Strange case update'
             stop
          endif
          !
          !     Check if already written
          !
          if(hftree(ipos,ifath)==0)then
             hftree(ipos,ifath)=iface
          else
             !
             !     Take the one with the lowest local level
             !
             jface=hftree(ipos,ifath)
             if(hftree(7,iface)<hftree(7,jface))then
                hftree(ipos,ifath)=iface
             else if(hftree(7,iface)==hftree(7,jface))then
                if(hftree(8,iface)>hftree(8,jface))then
                   hftree(ipos,ifath)=iface
                endif 
             endif
          endif
       endif
    enddo
    !
    !     Then decide the type
    !
    do iface=1,nface
       !write(*,*)iface,hftree(6,874)
       !
       !     The element should have been coarsened
       !     and not be at the coarsest level  
       !
       if(hftree(6,iface)==0 .and. hftree(8,iface)>0)then
          do ichil=1,3
             if(hftree(ichil,iface)==0)exit
          enddo
          !
          !     Test the child number
          !
          if(ichil==1)then
             !
             !     No children, it must be a -4
             ! 
             hftree(6,iface)=-4
             cycle

          else if(ichil==2)then 
             !
             !     Two children, end of the tree
             !
             hftree(6,iface)=2_ip
             hftree(6,hftree(1,iface))=-2_ip

          else if(ichil==4)then
             !
             !     Four children, end of the tree must check for children
             !     of children if refined 
             !
             !
             hftree(6,iface)=4_ip
             ison=hftree(1,iface)
             if(hftree(1,ison)==0)then
                hftree(6,ison)=-4_ip
             endif
             ison=hftree(2,iface)
             if(hftree(1,ison)==0)then
                hftree(6,ison)=-4_ip
             endif
             ison=hftree(3,iface)
             if(hftree(1,ison)==0)then
                hftree(6,ison)=-4_ip
             endif

          else
             write(*,*)'Error in update hftree'
             write(*,*)'ichil=',ichil
             stop 
          endif

       !else if(hftree(6,iface)>0)then
          !
          !     If the sons have coarsened, remark them
          !   
       !   nchil=hftree(6,iface)-1
       !   do je=1,nchil
       !      ichil=hftree(je,iface) 
       !      if(hftree(6,ichil)==0)then
       !         hftree(6,ichil)=-hftree(6,iface)
       !      endif
       !   enddo

       endif
    enddo

  end subroutine update

  subroutine verif(nelem,nnode,elem,htree,ntree,hftree,nftree,nface)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)  :: nelem,nnode,ntree,nftree,nface
    integer(ip), intent(in)  :: elem(nnode,nelem),htree(ntree,nelem)
    integer(ip), intent(in)  :: hftree(nftree,nface)
    integer(ip)              :: ielem,ifath,j,k,ison
    integer(ip)               :: ip11,ip12,ip13,ip14,ip21,ip22,ip23,ip24
    integer(ip)               :: ip31,ip32,ip33,ip34,ip41,ip42,ip43,ip44
    integer(ip)               :: ip51,ip52,ip53,ip54,ip61,ip62,ip63,ip64
    integer(ip)               :: ip71,ip72,ip73,ip74,ip81,ip82,ip83,ip84
    integer(ip)               :: nelem1,nelem2,nelem3,nelem4,nelem5,nelem6,nelem7
    integer(ip)               :: ichil,ichild,ncont,iface 
    !
    !     This sub checks the coherence of the tree and the element
    !
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nelem,htree) &
    !$omp& private(ielem,ifath,j,ison)

    do ielem=1,nelem
       !
       !     Check sons
       !
       if(htree(10,ielem)<0)then
          !
          !     Get father
          !
          ifath=htree(8,ielem)
          if(ifath==0)then
             write(*,*)'Error verif ifath==0'
             stop 
          endif
          !
          !     Must be the same type at the same level
          !
          if(htree(10,ifath)==-htree(10,ielem) .and.  htree(12,ifath)==htree(12,ielem))then
             do j=1,8
                if(htree(j,ifath)==ielem)exit
             enddo

             if(j==9)then
                write(*,*)'Error in verif'
                stop
             endif
          endif

          if(htree(10,ielem)/=-8)then
             do j=1,7
                if(htree(j,ielem)/=0)then
                   write(*,*)'Error in verif sons'
                   stop
                endif
             enddo
          endif
       else if(htree(10,ielem)>0)then
          !
          !     Check father
          !
          do j=1,7 
             ison=htree(j,ielem)
             if(ison==0)exit
             if(htree(8,ison)/=ielem)then 
                write(*,*)'Error in verif 2'
                stop
             endif
             if(htree(9,ison)/=j)then 
                write(*,*)'Error in verif 3'
                stop
             endif
          enddo
          if(j<htree(10,ielem))then
             write(*,*)'Error, sons not found'
             stop
          endif
          do j=max(1_ip,htree(10,ielem)),7
             if(htree(j,ielem)/=0)then
                write(*,*)'Error in verif 4'
                stop
             endif
          enddo
       else if(htree(10,ielem)==0)then
          if(htree(12,ielem)/=0)then
             write(*,*)'Error conformity htree(10,ielem)=0'
             stop
          endif
       endif
    enddo


    !
    !     The rest is false, points may change and rotate dues to succesive 2 and 4 refinement and coarsening.
    !


    !do ielem=1,nelem
    !   if(htree(10,ielem)==8)then
    !
    !     Check that no son has been refined
    !
    !      ncont=0_ip
    !      do ichil=1,7
    !         ichild=htree(ichil,ielem)
    !         if(htree(1,ichild)==0)then
    !            ncont=ncont+1
    !         endif
    !      enddo   
    !      if(ncont/=7)cycle
    !
    !     Get sons
    ! 
    !      nelem1=htree(1,ielem)
    !      nelem2=htree(2,ielem)
    !      nelem3=htree(3,ielem)
    !      nelem4=htree(4,ielem)
    !      nelem5=htree(5,ielem)
    !      nelem6=htree(6,ielem)
    !      nelem7=htree(7,ielem)
    !
    !     DBG 
    !
    !     ip11=elem(1,ielem)
    !     ip12=elem(2,ielem)
    !     ip13=elem(3,ielem)
    !     ip14=elem(4,ielem)
    !     ip21=elem(1,nelem1)
    !     ip22=elem(2,nelem1)
    !     ip23=elem(3,nelem1)
    !     ip24=elem(4,nelem1)
    !     ip31=elem(1,nelem2)
    !     ip32=elem(2,nelem2)
    !     ip33=elem(3,nelem2)
    !     ip34=elem(4,nelem2)
    !     ip41=elem(1,nelem3)
    !     ip42=elem(2,nelem3)
    !     ip43=elem(3,nelem3)
    !     ip44=elem(4,nelem3)
    !     ip51=elem(1,nelem4)
    !     ip52=elem(2,nelem4)
    !     ip53=elem(3,nelem4)
    !     ip54=elem(4,nelem4)
    !     ip61=elem(1,nelem5)
    !     ip62=elem(2,nelem5)
    !     ip63=elem(3,nelem5)
    !     ip64=elem(4,nelem5)
    !     ip71=elem(1,nelem6)
    !     ip72=elem(2,nelem6)
    !     ip73=elem(3,nelem6)
    !     ip74=elem(4,nelem6)
    !     ip81=elem(1,nelem7)
    !     ip82=elem(2,nelem7)
    !     ip83=elem(3,nelem7)
    !     ip84=elem(4,nelem7)
    !   1 5 6 7 
    !   5 2 8 9 
    !   8 3 6 10 
    !   5 8 6 9   should be 5 9 8 6 but easier for 8:4+
    !   7 9 5 6 
    !   10 9 7 6 
    !   8 9 10 6 
    !   7 9 10 4 
    !     if(ip12/=ip21 .or. ip21/=ip41 .or. ip41/=ip53)then
    !        write(*,*)'Error verif 8 ip5'
    !        stop
    !     endif
    !     if(ip13/=ip33 .or. ip33/=ip43 .or. ip43/=ip54 .or.ip54 /=ip64 .or.ip64/=ip74)then
    !        write(*,*)'Error verif 8 ip6'
    !        stop
    !     endif
    !     if(ip14/=ip51 .or. ip51/=ip63 .or. ip63/=ip81)then
    !        write(*,*)'Error verif 8 ip7'
    !        stop
    !     endif
    !     if(ip23/=ip31 .or. ip31/=ip42 .or. ip42/=ip71)then
    !        write(*,*)'Error verif 8 ip8'
    !        stop
    !     endif
    !     if(ip24/=ip44 .or. ip44/=ip52 .or. ip52/=ip62 .or. ip62/=ip72 .or. ip72/=ip82)then
    !        write(*,*)'Error verif 8 ip9'
    !        stop
    !     endif
    !     if(ip34/=ip61 .or. ip61/=ip73 .or. ip73/=ip83)then
    !        write(*,*)'Error verif 8 ip10'
    !        stop
    !     endif
    !       else if(htree(10,ielem)==4)then
    !
    !     Get sons
    ! 
    !      nelem1=htree(1,ielem)
    !      nelem2=htree(2,ielem)
    !      nelem3=htree(3,ielem)
    !      !
    !     DBG 
    !
    !      ip11=elem(1,ielem)
    !      ip12=elem(2,ielem)
    !      ip13=elem(3,ielem)
    !      ip14=elem(4,ielem)
    !      ip21=elem(1,nelem1)
    !      ip22=elem(2,nelem1)
    !      ip23=elem(3,nelem1)
    !      ip24=elem(4,nelem1)
    !      ip31=elem(1,nelem2)
    !      ip32=elem(2,nelem2)
    !      ip33=elem(3,nelem2)
    !      ip34=elem(4,nelem2)
    !      ip41=elem(1,nelem3)
    !      ip42=elem(2,nelem3)
    !      ip43=elem(3,nelem3)
    !      ip44=elem(4,nelem3)
    !   1 5 7 4
    !   5 2 6 4 
    !   6 3 7 4 
    !   5 6 7 4 

    !      if(ip14/=ip24 .or. ip24/=ip34 .or. ip34/=ip44)then
    !         write(*,*)'Error verif 4 ip4'
    !         stop
    !      endif
    !      if(ip12/=ip21 .or. ip21/=ip41 )then
    !         write(*,*)'Error verif 4 ip1'
    !         stop
    !      endif
    !      if(ip13/=ip33 .or. ip33/=ip43 )then
    !         write(*,*)'Error verif 4 ip7'
    !         stop
    !      endif
    !      if(ip23/=ip31 .or. ip31/=ip42 )then
    !         write(*,*)'Error verif 4 ip6'
    !         stop
    !      endif

    !   else if(htree(10,ielem)==2)then
    !
    !     Get sons
    ! 
    !      nelem1=htree(1,ielem)
    !
    !     DBG 
    !
    !      ip11=elem(1,ielem)
    !      ip12=elem(2,ielem)
    !      ip13=elem(3,ielem)
    !      ip14=elem(4,ielem)
    !      ip21=elem(1,nelem1)
    !      ip22=elem(2,nelem1)
    !      ip23=elem(3,nelem1)
    !      ip24=elem(4,nelem1)
    !      !         1 5 3 4 
    !         5 2 3 4

    !      if(ip12/=ip21 )then
    !         write(*,*)'Error verif 2 ip5'
    !         stop
    !      endif
    !      if(ip13/=ip23  )then
    !         write(*,*)'Error verif 4 ip3'
    !         stop
    !      endif
    !      if(ip14/=ip24 )then
    !         write(*,*)'Error verif 4 ip4'
    !         stop
    !      endif
    !

    !      endif


    !   enddo

    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nface,hftree) &
    !$omp& private(iface,ifath,j,ison)

    do iface=1,nface
       !write(*,*)iface
       !
       !     Check sons
       !
       if(hftree(6,iface)<0)then
          !
          !     Get father
          !
          ifath=hftree(4,iface)
          if(ifath==0)then
             write(*,*)'Error verif ifath==0'
             stop 
          endif
          ! 
          !     Must be the same type at the same level
          !
          if(hftree(6,ifath)==-hftree(6,iface) .and.  hftree(8,ifath)==hftree(8,iface))then
             do j=1,3
                if(hftree(j,ifath)==iface)exit
             enddo

             if(j==4)then
                write(*,*)'Error in verif face'
                stop
             endif
          endif
          if(hftree(6,iface)/=-4)then
             do j=1,3
                if(hftree(j,iface)/=0)then
                   write(*,*)'Error in verif sons face'
                   stop
                endif
             enddo
          endif
       else if(hftree(6,iface)>0)then
          do j=1,3 
             ison=hftree(j,iface)
             if(ison==0)exit
             if(hftree(4,ison)/=iface)then 
                write(*,*)'Error in verif 2 face'
                stop
             endif
             if(hftree(5,ison)/=j)then 
                write(*,*)'Error in verif 3 face'
                stop
             endif
          enddo

          if(j<hftree(6,iface))then
             write(*,*)'Error, sons not found'
             stop 
          endif

          do j=max(1_ip,hftree(6,iface)),3
             if(hftree(j,iface)/=0)then
                write(*,*)'Error in verif 4 face'
                stop
             endif
          enddo
        else if(hftree(6,iface)==0)then
          if(hftree(8,iface)/=0)then
             write(*,*)'Error conformity hftree(6,iface)=0'
             stop
          endif
       endif
    enddo


  end subroutine verif

  subroutine gtvol(elem,ielem,coor,rvol,nnode,nelem,npoin,ndim)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnode,nelem,npoin,ndim,ielem
    integer(ip), intent(in)      :: elem(nnode,nelem)
    real(rp),intent(in)          :: coor(ndim,npoin)
    real(rp),intent(inout)       :: rvol
    integer(ip)                  :: ip1,ip2,ip3,ip4
    real(rp)                     :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,c16  
    real(rp)                     :: x21,y21,z21,x31,y31,z31,x41,y41,z41  

    c16=1.0d+00/6.0d+00
    ip1=elem(1,ielem)
    ip2=elem(2,ielem)
    ip3=elem(3,ielem)
    ip4=elem(4,ielem)

    x1=coor(1,ip1)
    y1=coor(2,ip1)
    z1=coor(3,ip1)
    x2=coor(1,ip2)
    y2=coor(2,ip2)
    z2=coor(3,ip2)
    x3=coor(1,ip3)
    y3=coor(2,ip3)
    z3=coor(3,ip3)
    x4=coor(1,ip4)
    y4=coor(2,ip4)
    z4=coor(3,ip4)

    x21=x2-x1
    y21=y2-y1
    z21=z2-z1

    x31=x3-x1
    y31=y3-y1
    z31=z3-z1

    x41=x4-x1
    y41=y4-y1
    z41=z4-z1

    rvol=c16*(x21*(y31*z41-z31*y41) + &
         x31*(z21*y41-y21*z41) + &
         x41*(y21*z31-z21*y31))

  end subroutine gtvol

  subroutine verifface(nface,nnofa,lface,nnode,nelem,elem,npoin)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        : memor_msh
    implicit none
    integer(ip), intent(in)  :: nelem,nnode,nface,nnofa,npoin
    integer(ip), intent(in)  :: elem(nnode,nelem),lface(nnofa,nface)
    integer(ip)              :: ielem,iface,inode,ineigh,ncont,isto
    integer(ip)              :: ip1,ip2,ip3
    integer(ip),pointer      :: ptoel1(:),ptoel2(:),ptofa1(:),ptofa2(:) 
    integer(ip),pointer      :: eltoel(:,:),lmark(:) 
    integer(4)                :: istat
    integer(ip)               :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))
    !
    !     This sub checks coherence between faces and volumes
    !
    nullify(ptoel1,ptoel2,ptofa1,ptofa2,eltoel)
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','verifface',lmark)
    !
    !     Get the elem surrounding the points
    !
    call ptoelm(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptofa1,ptofa2)
    !
    !     Get the faces surrounding faces
    !
    call tetote(elem,nnode,nelem,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Loop on faces
    !
    do ielem=1,nelem
       do inode=1,nnode
          ineigh=eltoel(inode,ielem) 
          if(ineigh==0)then
             ip1=elem(ltab(1,inode),ielem) 
             ip2=elem(ltab(2,inode),ielem) 
             ip3=elem(ltab(3,inode),ielem) 

             lmark(ip1)=1_ip 
             lmark(ip2)=1_ip 
             lmark(ip3)=1_ip 

             do isto=ptofa2(ip1),ptofa2(ip1+1)-1_ip
                iface=ptofa1(isto)
                ncont=0_ip 
                ncont=ncont+lmark(lface(1,iface))
                ncont=ncont+lmark(lface(2,iface))
                ncont=ncont+lmark(lface(3,iface))
                if(ncont==3)exit

             enddo

             if(isto==ptofa2(ip1+1))then
                write(*,*)'Error in verifface'
                stop
             endif

             lmark(ip1)=1_ip 
             lmark(ip2)=1_ip 
             lmark(ip3)=1_ip 

          endif
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'ELTOEL','verifface',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','verifface',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOFA1','verifface',ptofa1)
    deallocate(ptofa1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOFA1','verifface',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOFA2','verifface',ptofa2)
    deallocate(ptofa2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOFA2','verifface',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','verifface',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','verifface',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','verifface',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','verifface',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','verifface',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','verifface',0_ip)

  end subroutine verifface

  subroutine ref84(nnode,nelem,elem,npoin,ielem,ledglm,nsid,lptab,ntree,htree,coor,ndim,&
       ison1,ison2,ison3,nelem1,nelem2,nelem3,nelem8)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        : memor_msh
    implicit none
    integer(ip),intent(in)    :: nnode,npoin,ielem,nsid,ntree,ndim
    integer(ip),intent(in)    :: ison1,ison2,ison3,nelem1,nelem2,nelem3,nelem8
    integer(ip),intent(inout) :: nelem
    integer(ip),intent(in)    :: ledglm(nsid,nelem),lptab(10)
    real(rp),intent(in)       :: coor(ndim,npoin)
    integer(ip),pointer       :: elem(:,:),htree(:,:)
    integer(ip)               :: ncont,ledgl(3),lptab2(10),letab(8)
    integer(ip)               :: ipedg1
    !
    !     Count first triangle
    ! 
    ncont=0
    if(ledglm(1,ielem)/=0)then
       ncont=ncont+1
       ledgl(ncont)=1
    endif
    if(ledglm(2,ielem)/=0)then
       ncont=ncont+1
       ledgl(ncont)=2
    endif
    if(ledglm(4,ielem)/=0)then
       ncont=ncont+1
       ledgl(ncont)=4
    endif
    !
    !     How many edges are marked for refinement?
    !
    if(ncont==1)then

       if(ledglm(1,ielem)/=0)then

          ipedg1=ledglm(1,ielem)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(1) 
          lptab2(2)=lptab(5) 
          lptab2(3)=lptab(6)
          lptab2(4)=lptab(7) 
          lptab2(5)=ipedg1
          letab(1)=nelem1
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)


       else if(ledglm(2,ielem)/=0)then

          ipedg1=ledglm(2,ielem)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(6) 
          lptab2(2)=lptab(1) 
          lptab2(3)=lptab(5)
          lptab2(4)=lptab(7) 
          lptab2(5)=ipedg1
          letab(1)=nelem1
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else if(ledglm(4,ielem)/=0)then

          ipedg1=ledglm(4,ielem)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(5) 
          lptab2(2)=lptab(6) 
          lptab2(3)=lptab(1)
          lptab2(4)=lptab(7) 
          lptab2(5)=ipedg1
          letab(1)=nelem1
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else

          write(*,*)'Error in ref84 a'
          stop
       endif

    else if(ncont==3)then

       nelem=nelem+3
       call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
       call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

       lptab2(1)=lptab(1) 
       lptab2(2)=lptab(5) 
       lptab2(3)=lptab(6)
       lptab2(4)=lptab(7)
       lptab2(5)=ledglm(1,ielem)
       lptab2(6)=ledglm(4,ielem)
       lptab2(7)=ledglm(2,ielem)

       letab(1)=nelem1
       letab(2)=nelem-2 
       letab(3)=nelem-1 
       letab(4)=nelem 

       call ref4(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
            letab,lptab2)

    endif
    !
    !     Count second triangle
    ! 
    ncont=0
    if(ledglm(1,ison1)/=0)then
       ncont=ncont+1
       ledgl(ncont)=1
    endif
    if(ledglm(2,ison1)/=0)then
       ncont=ncont+1
       ledgl(ncont)=2
    endif
    if(ledglm(4,ison1)/=0)then
       ncont=ncont+1
       ledgl(ncont)=4
    endif
    !
    !     How many edges are marked for refinement?
    !
    if(ncont==1)then

       if(ledglm(1,ison1)/=0)then

          ipedg1=ledglm(1,ison1)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(5) 
          lptab2(2)=lptab(2) 
          lptab2(3)=lptab(8)
          lptab2(4)=lptab(9) 
          lptab2(5)=ipedg1
          letab(1)=nelem2 
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else if(ledglm(2,ison1)/=0)then

          ipedg1=ledglm(2,ison1)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(8) 
          lptab2(2)=lptab(5) 
          lptab2(3)=lptab(2)
          lptab2(4)=lptab(9) 
          lptab2(5)=ipedg1
          letab(1)=nelem2 
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else if(ledglm(4,ison1)/=0)then

          ipedg1=ledglm(4,ison1)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(2) 
          lptab2(2)=lptab(8) 
          lptab2(3)=lptab(5)
          lptab2(4)=lptab(9) 
          lptab2(5)=ipedg1
          letab(1)=nelem2 
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else

          write(*,*)'Error in ref84 b'
          stop

       endif

    else if(ncont==3)then

       nelem=nelem+3
       call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
       call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

       lptab2(1)=lptab(5) 
       lptab2(2)=lptab(2) 
       lptab2(3)=lptab(8)
       lptab2(4)=lptab(9)
       lptab2(5)=ledglm(1,ison1)
       lptab2(6)=ledglm(4,ison1)
       lptab2(7)=ledglm(2,ison1)

       letab(1)=nelem2 
       letab(2)=nelem-2 
       letab(3)=nelem-1 
       letab(4)=nelem 

       call ref4(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
            letab,lptab2)

    endif
    !
    !     Count third triangle
    ! 
    ncont=0
    if(ledglm(1,ison2)/=0)then
       ncont=ncont+1
       ledgl(ncont)=1
    endif
    if(ledglm(2,ison2)/=0)then
       ncont=ncont+1
       ledgl(ncont)=2
    endif
    if(ledglm(4,ison2)/=0)then
       ncont=ncont+1
       ledgl(ncont)=4
    endif
    !
    !     How many edges are marked for refinement?
    !
    if(ncont==1)then

       if(ledglm(1,ison2)/=0)then

          ipedg1=ledglm(1,ison2)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(8) 
          lptab2(2)=lptab(3) 
          lptab2(3)=lptab(6)
          lptab2(4)=lptab(10) 
          lptab2(5)=ipedg1
          letab(1)=nelem3 
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else if(ledglm(2,ison2)/=0)then

          ipedg1=ledglm(2,ison2)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(6) 
          lptab2(2)=lptab(8) 
          lptab2(3)=lptab(3)
          lptab2(4)=lptab(10) 
          lptab2(5)=ipedg1
          letab(1)=nelem3 
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else if(ledglm(4,ison2)/=0)then

          ipedg1=ledglm(4,ison2)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(3) 
          lptab2(2)=lptab(6) 
          lptab2(3)=lptab(8)
          lptab2(4)=lptab(10) 
          lptab2(5)=ipedg1
          letab(1)=nelem3 
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else

          write(*,*)'Error in ref84 b'
          stop

       endif

    else if(ncont==3)then

       nelem=nelem+3
       call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
       call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

       lptab2(1)=lptab(8) 
       lptab2(2)=lptab(3) 
       lptab2(3)=lptab(6)
       lptab2(4)=lptab(10)
       lptab2(5)=ledglm(1,ison2)
       lptab2(6)=ledglm(4,ison2)
       lptab2(7)=ledglm(2,ison2)

       letab(1)=nelem3 
       letab(2)=nelem-2 
       letab(3)=nelem-1 
       letab(4)=nelem 

       call ref4(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
            letab,lptab2)

    endif
    !
    !     Count fourth triangle
    ! 
    ncont=0
    if(ledglm(1,ison3)/=0)then
       ncont=ncont+1
       ledgl(ncont)=1
    endif
    if(ledglm(2,ison3)/=0)then
       ncont=ncont+1
       ledgl(ncont)=2
    endif
    if(ledglm(4,ison3)/=0)then
       ncont=ncont+1
       ledgl(ncont)=4
    endif
    !
    !     How many edges are marked for refinement?
    !
    if(ncont==0)then
       !
       !     Nothing to do
       !
    else if(ncont==1)then

       if(ledglm(1,ison3)/=0)then

          ipedg1=ledglm(1,ison3)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(5) 
          lptab2(2)=lptab(8) 
          lptab2(3)=lptab(6)
          lptab2(4)=lptab(9) 
          lptab2(5)=ipedg1
          letab(1)=nelem8 
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else if(ledglm(4,ison3)/=0)then

          ipedg1=ledglm(4,ison3)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(8) 
          lptab2(2)=lptab(6) 
          lptab2(3)=lptab(5)
          lptab2(4)=lptab(9) 
          lptab2(5)=ipedg1
          letab(1)=nelem8 
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else if(ledglm(2,ison3)/=0)then

          ipedg1=ledglm(2,ison3)

          nelem=nelem+1
          call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
          call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

          lptab2(1)=lptab(6) 
          lptab2(2)=lptab(5) 
          lptab2(3)=lptab(8)
          lptab2(4)=lptab(9) 
          lptab2(5)=ipedg1
          letab(1)=nelem8 
          letab(2)=nelem 

          call ref2(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
               letab,lptab2)

       else

          write(*,*)'Error in ref84 aa'
          stop

       endif

    else if(ncont==3)then

       nelem=nelem+3
       call memrea(nelem,memor_msh,'ELEM','hrefi',elem)
       call memrea(nelem,memor_msh,'HTREE','hrefi',htree)

       lptab2(1)=lptab(5) 
       lptab2(2)=lptab(8) 
       lptab2(3)=lptab(6)
       lptab2(4)=lptab(9)
       lptab2(5)=ledglm(1,ison3)
       lptab2(6)=ledglm(4,ison3)
       lptab2(7)=ledglm(2,ison3)

       letab(1)=nelem8
       letab(2)=nelem-2 
       letab(3)=nelem-1 
       letab(4)=nelem 

       call ref4(nelem,elem,htree,nnode,ntree,coor,ndim,npoin,&
            letab,lptab2)

    else

       write(*,*)'Error in ref84'
       stop

    endif

  end subroutine ref84

  subroutine limitref(nelem,ntree,htree,lmark,nlaymax)
   use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        : memor_msh
    implicit none
    integer(ip),intent(in)    :: nelem,ntree,nlaymax
    integer(ip),intent(in)    :: htree(ntree,nelem) 
    integer(ip),intent(inout) :: lmark(nelem) 
    integer(ip)               :: ielem  
    
    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(htree,nelem,lmark,nlaymax) &
    !$omp& private(ielem)

    do ielem=1,nelem
       if(lmark(ielem)==1)then
          if(htree(12,ielem)>=nlaymax)then
             lmark(ielem)=0
          endif
       endif  
    enddo

  end subroutine limitref
  


end module mod_href









