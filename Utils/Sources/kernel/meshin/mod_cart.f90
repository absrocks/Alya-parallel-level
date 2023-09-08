module mod_cart
  !
  !     Use openmp  $OMP1
  ! 
  use mod_mshtol

contains

  subroutine carmsh(npoin,ndim,nface,nnofa,nsour,rsuni,rscal,ifbox,boxin,rtol,&
       lcell,ncell,coor,lface,rsour,rsgeo,ismoo,isizcrit,iinter,icartout,irefsurf,&
       isizuni)
    use def_meshin, only          : memor_msh
    use def_meshin, only          : 
    use def_kintyp, only          : ip,rp,lg,cell
    use mod_memchk
    implicit none
    integer(ip),intent(in)    ::  ndim,npoin,nface,nnofa,nsour,ifbox,ismoo
    integer(ip),intent(in)    ::  isizcrit,iinter,icartout,irefsurf,isizuni
    real(rp),intent(in)       ::  rsuni,rscal,boxin(ndim,2),coor(ndim,npoin)
    integer(ip),intent(in)    ::  lface(nnofa,nface)
    real(rp),intent(inout)    ::  rtol
    integer(ip),intent(inout) ::  ncell
    type(cell), pointer       ::  lcell(:)
    real(rp),pointer          ::  rsour(:,:),rsgeo(:,:,:)
    real(rp)                  ::  bboxbin(ndim,2),dxbin,dybin,dzbin
    real(rp)                  ::  dx,dy,dz,bbox(ndim,2),c13,cx,cy,cz,rsmin,c05 
    integer(ip)               ::  nbin,nbinx,nvoxx,nvoxy,nvoxz,niter,iface,icell
    integer(ip)               ::  nfacet,nbint 
    integer(ip),pointer       ::  limp(:),eltoel(:,:),btoel1(:),btoel2(:)  
    integer(ip),pointer       ::  ptoel1(:),ptoel2(:)  
    real(rp),pointer          ::  trisize(:),tribox(:,:),rnofa(:,:)
    integer(4)                ::  istat
    !
    !     The %marked is:            -2 at the creation of the cell
    !                                -1 if do not intersect face
    !                                -3 if intersect but agrees with the size
    !                                 1 if marked  
    !
    c13=1.0d+00/3.0d+00
    c05=0.5d+00
    !
    !     Do we have faces or just sources
    !
    if(nface/=0_ip)then

       allocate(rnofa(ndim,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'RNOFA','cart',rnofa)
       !
       !     Get the faces surrounding the points
       !
       call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
       !
       !     Get the face normals
       !
       call gtfnrl(lface,nface,nnofa,ndim,coor,npoin,rnofa)
       !
       !     Get the faces surrounding faces
       !
       call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
       !
       allocate(trisize(nface),stat=istat)
       call memchk(zero,istat,memor_msh,'TRISIZE','cart',trisize)
       allocate(tribox(6,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'TRIBOX','cart',tribox)
       !
       !     Allocate the bin to face pointer
       !
       nbin=nface
       nbinx=floor(nbin**c13)+50_ip
       nbin=nbinx*nbinx*nbinx
       allocate(btoel2(nbin+1),stat=istat)
       call memchk(zero,istat,memor_msh,'BTOEL2','cart',btoel2)
       !
       !     Compute bbox of the triangulation
       !
       call boxbin(coor,ndim,npoin,bboxbin)
       !
       !     Compute the bbox of each triangle and their size
       !
       call bboxtri(coor,lface,nface,nnofa,tribox,trisize,ndim,npoin)
       !
       !     Get the minimum size of the faces
       !
       rsmin=trisize(1)
       do iface=2,nface
          if(rsmin<trisize(iface))rsmin=trisize(iface)
       enddo
       !
       !     Compute extent of the bin
       !
       dxbin=(bboxbin(1,2)-bboxbin(1,1))*c05
       dybin=(bboxbin(2,2)-bboxbin(2,1))*c05
       dzbin=(bboxbin(3,2)-bboxbin(3,1))*c05
       !
       !     Do we have a degenerate case?
       !
       if(dxbin<rsmin)dxbin=rsmin
       if(dybin<rsmin)dybin=rsmin
       if(dzbin<rsmin)dzbin=rsmin
       !
       !     Recompute bboxbin
       !
       cx=(bboxbin(1,2)+bboxbin(1,1))*c05
       cy=(bboxbin(2,2)+bboxbin(2,1))*c05
       cz=(bboxbin(3,2)+bboxbin(3,1))*c05

       bboxbin(1,1)=cx-dxbin  
       bboxbin(2,1)=cy-dybin  
       bboxbin(3,1)=cz-dzbin  
       bboxbin(1,2)=cx+dxbin  
       bboxbin(2,2)=cy+dybin  
       bboxbin(3,2)=cz+dzbin  
       !
       !     Compute dxbin,dybin,dzbin of the triangulation
       !
       dxbin=(bboxbin(1,2)-bboxbin(1,1))/nbinx
       dybin=(bboxbin(2,2)-bboxbin(2,1))/nbinx
       dzbin=(bboxbin(3,2)-bboxbin(3,1))/nbinx
       !
       !     Fill the bin
       !
       call tri2bin(tribox,nface,bboxbin,dxbin,dybin,dzbin,nbin,nbinx,&
            trisize,ndim,btoel2,btoel1)

    else
       !
       !     Allocate dummy variables
       !  
       nfacet=1_ip
       nbint=1_ip
       allocate(rnofa(ndim,nfacet),stat=istat)
       call memchk(zero,istat,memor_msh,'RNOFA','cart',rnofa)
       allocate(trisize(nfacet),stat=istat)
       call memchk(zero,istat,memor_msh,'TRISIZE','cart',trisize)
       allocate(tribox(6,nfacet),stat=istat)
       call memchk(zero,istat,memor_msh,'TRIBOX','cart',tribox)
       allocate(btoel2(nbint+1),stat=istat)
       call memchk(zero,istat,memor_msh,'BTOEL2','cart',btoel2)
       allocate(btoel1(nbint),stat=istat)
       call memchk(zero,istat,memor_msh,'BTOEL1','cart',btoel1)
       allocate(ptoel1(1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOEL1','cart',ptoel1)
       allocate(ptoel2(1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOEL2','cart',ptoel2)
       allocate(eltoel(1,1),stat=istat)
       call memchk(zero,istat,memor_msh,'ELTOEL','cart',eltoel)

    endif
    !
    !     Compute dx,dy,dz from Bbox, the bounding box of the cartesian mesh 
    !
    if(ifbox==0) then
       call boxcart(ndim,ncell,nvoxx,nvoxy,nvoxz,bboxbin,bbox,dx,dy,dz,rsuni)
    else
       call boxcart(ndim,ncell,nvoxx,nvoxy,nvoxz,boxin,  bbox,dx,dy,dz,rsuni)
    end if
    !
    !     Is everything ok?
    !
    if(ncell<=0)then
       write(*,*)'Error in cartesian mesh, ncell=',ncell
       write(*,*)'Check uniform size:',rsuni
       write(*,*)'Box dimension',bbox(1,1)-bbox(1,2),bbox(2,1)-bbox(2,2),bbox(3,1)-bbox(3,2)
       stop
    endif
    !
    !     Get tolerance
    !
    rtol=min(dx,dy,dz)*1.0d-08
    rtol=0.0d+00
    !
    !     Build the initial mesh
    !
    call buildIni(nvoxx,nvoxy,nvoxz,bbox,dx,dy,dz,ndim,lcell,ncell)
    !
    !     Subdivide the mesh
    ! 
    call mesh(btoel1,btoel2,tribox,bboxbin,dxbin,dybin,dzbin,nface,nbin,nbinx,&
         trisize,lface,coor,npoin,ndim,eltoel,rnofa,nnofa,nsour,rsour,rsgeo,niter,&
         lcell,ncell,rtol,rsuni,iinter,irefsurf)
    !
    !     The mesh has been adapted with a given criterion
    !     What is the size criterion?
    !    
    !
    !     Do we have to apply a uniform size first?
    !
    if(isizuni==1)then
       !
       !     Apply uniform size and delete the previous size distribution
       !
       call unisiz(ncell,lcell,rsuni)

    endif
    !
    !     Now decide depending on isizcrit
    !
    if(isizcrit==1)then
       !
       !     Do nothing
       !
    else if(isizcrit==2)then
       !
       !     Apply minimum between cell size and uniform size
       !
       call unisizmin(ncell,lcell,rsuni)

    else if(isizcrit==3)then 
       !
       !     Verify that we have sources
       !
       if(nsour==0)then
          write(*,*)'Error in input, size given by sources but there are no sources'
          stop  
       endif
       !
       !     Apply source size
       !
       call soursiz(ncell,lcell,ndim,nsour,rsour,rsgeo,rtol,rsuni)

    else
       write(*,*)'isizcrit not defined'
       stop

    endif
    !
    !     Scale the size
    !
    call scalsize(ncell,lcell,rscal)
    !
    !     Do we want to smooth the size
    !
    if(ismoo==1)then
       !
       !     Allocate limp, the array of imposed size
       !
       allocate(limp(ncell),stat=istat)
       call memchk(zero,istat,memor_msh,'LIMP','cart',limp)
       !
       !     Find the minima 
       !
       call findMin(ncell,lcell,limp)
       !
       !     Smooth out the size
       !
       call smoothsize(lcell,ncell,niter,limp)
       !
       !     Deallocate limp
       !
       call memchk(2_ip,istat,memor_msh,'LIMP','cart',limp)
       deallocate(limp,stat=istat)
       if(istat/=0) call memerr(2_ip,'LIMP','cart',0_ip)

    endif
    !
    !     Renumber the cells
    !
    !call renucart(lcell,ncell)

    if(icartout==1)then
       write(*,*)'Cartesian mesh done, writing output'
       call writeGidconf(ncell,lcell,lface,nface,coor,npoin,ndim,nnofa)
    endif

    call memchk(2_ip,istat,memor_msh,'BTOEL1','cart',btoel1)
    deallocate(btoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'BTOEL1','cart',0_ip)
    call memchk(2_ip,istat,memor_msh,'BTOEL2','cart',btoel2)
    deallocate(btoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'BTOEL2','cart',0_ip)
    call memchk(2_ip,istat,memor_msh,'TRIBOX','cart',tribox)
    deallocate(tribox,stat=istat)
    if(istat/=0) call memerr(2_ip,'TRIBOX','cart',0_ip)
    call memchk(2_ip,istat,memor_msh,'TRISIZE','cart',trisize)
    deallocate(trisize,stat=istat)
    if(istat/=0) call memerr(2_ip,'TRISIZE','cart',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOFA','cart',rnofa)
    deallocate(rnofa,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOFA','cart',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','cart',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','cart',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','cart',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','cart',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','cart',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','cart',0_ip)

  end subroutine carmsh


  subroutine boxbin(coor,ndim,npoin,bboxbin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in) ::  ndim,npoin
    real(rp),intent(in)    ::  coor(ndim,npoin)
    real(rp),intent(inout) ::  bboxbin(ndim,2)
    real(rp)               ::  rmin1,rmin2,rmin3,rmax1,rmax2,rmax3,dif(ndim),c(ndim)
    integer(ip)            ::  ipoin
    real(rp)               ::   c05,coef1
    !
    !   This sub computes the bbox of the triangulation
    !

    c05=0.5d+00
    coef1=1.07d+00  

    rmin1=coor(1,1)
    rmin2=coor(2,1)
    rmin3=coor(3,1)

    rmax1=coor(1,1)
    rmax2=coor(2,1)
    rmax3=coor(3,1)

    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(coor,npoin) &
    !$omp& reduction(min:rmin1)&
    !$omp& reduction(min:rmin2)&
    !$omp& reduction(min:rmin3)&
    !$omp& reduction(max:rmax1)&
    !$omp& reduction(max:rmax2)&
    !$omp& reduction(max:rmax3)&
    !$omp& private(ipoin)

    do ipoin=2,npoin

       if(coor(1,ipoin)<rmin1)rmin1=coor(1,ipoin)
       if(coor(2,ipoin)<rmin2)rmin2=coor(2,ipoin)
       if(coor(3,ipoin)<rmin3)rmin3=coor(3,ipoin)

       if(coor(1,ipoin)>rmax1)rmax1=coor(1,ipoin)
       if(coor(2,ipoin)>rmax2)rmax2=coor(2,ipoin)
       if(coor(3,ipoin)>rmax3)rmax3=coor(3,ipoin)

    enddo
    !
    !      Make it a litter bigger
    !
    c(1)=(rmin1+rmax1)*c05
    c(2)=(rmin2+rmax2)*c05
    c(3)=(rmin3+rmax3)*c05

    dif(1)=rmax1-c(1) 
    dif(2)=rmax2-c(2) 
    dif(3)=rmax3-c(3) 

    rmin1=c(1)-coef1*dif(1)
    rmin2=c(2)-coef1*dif(2)
    rmin3=c(3)-coef1*dif(3)

    rmax1=c(1)+coef1*dif(1)
    rmax2=c(2)+coef1*dif(2)
    rmax3=c(3)+coef1*dif(3)

    bboxbin(1,1)=rmin1
    bboxbin(2,1)=rmin2
    bboxbin(3,1)=rmin3

    bboxbin(1,2)=rmax1
    bboxbin(2,2)=rmax2
    bboxbin(3,2)=rmax3

  end subroutine boxbin

  subroutine bboxtri(coor,lface,nface,nnofa,tribox,trisize,ndim,npoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in) ::  npoin,nface,ndim,nnofa
    integer(ip),intent(in) ::  lface(nnofa,nface)
    real(rp),intent(in)    ::  coor(ndim,npoin)
    real(rp),intent(inout) ::  tribox(6,nface),trisize(nface)
    integer(ip)            ::  iface,p1,p2,p3
    real(rp)               ::  dx,dy,dz,dl,xmin,ymin,zmin,xmax,ymax,zmax,dl1,dl2,dl3
    !
    !     This sub builds the bbox of the triangulation
    !

    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(lface,nface,coor,tribox,trisize) &
    !$omp& private(iface,xmin,ymin,zmin,xmax,ymax,zmax)&
    !$omp& private(p1,p2,p3,dx,dy,dz,dl1,dl2,dl3,dl)

    do iface=1,nface

       p1=lface(1,iface)
       p2=lface(2,iface)
       p3=lface(3,iface)
       !
       !    Initialize min,max
       !
       xmin=coor(1,p1)
       ymin=coor(2,p1)
       zmin=coor(3,p1)  
       xmax=xmin
       ymax=ymin
       zmax=zmin
       !
       !     Compare with p2
       !
       if(coor(1,p2)<xmin)xmin=coor(1,p2)
       if(coor(1,p2)>xmax)xmax=coor(1,p2)

       if(coor(2,p2)<ymin)ymin=coor(2,p2)
       if(coor(2,p2)>ymax)ymax=coor(2,p2)

       if(coor(3,p2)<zmin)zmin=coor(3,p2)
       if(coor(3,p2)>zmax)zmax=coor(3,p2)
       !
       !     Compare with p3
       !
       if(coor(1,p3)<xmin)xmin=coor(1,p3)
       if(coor(1,p3)>xmax)xmax=coor(1,p3)

       if(coor(2,p3)<ymin)ymin=coor(2,p3)
       if(coor(2,p3)>ymax)ymax=coor(2,p3)

       if(coor(3,p3)<zmin)zmin=coor(3,p3)
       if(coor(3,p3)>zmax)zmax=coor(3,p3)
       !
       !     Remember bbox and size
       !
       tribox(1,iface)=xmin
       tribox(2,iface)=ymin
       tribox(3,iface)=zmin  
       tribox(4,iface)=xmax
       tribox(5,iface)=ymax
       tribox(6,iface)=zmax  

       !dx=xmax-xmin 
       !dy=ymax-ymin 
       !dz=zmax-zmin
       !dl=dx
       !if(dy>dl)dl=dy
       !if(dz>dl)dl=dz

       !
       !     Compute shortest edge
       !
       dx=coor(1,p1)-coor(1,p2) 
       dy=coor(2,p1)-coor(2,p2) 
       dz=coor(3,p1)-coor(3,p2) 
       dl1=sqrt(dx*dx+dy*dy+dz*dz)
       dx=coor(1,p1)-coor(1,p3) 
       dy=coor(2,p1)-coor(2,p3) 
       dz=coor(3,p1)-coor(3,p3) 
       dl2=sqrt(dx*dx+dy*dy+dz*dz)
       dx=coor(1,p2)-coor(1,p3) 
       dy=coor(2,p2)-coor(2,p3) 
       dz=coor(3,p2)-coor(3,p3) 
       dl3=sqrt(dx*dx+dy*dy+dz*dz)
       dl=min(dl1,dl2,dl3)

       trisize(iface)=dl

    enddo

  end subroutine bboxtri


  subroutine tri2bin(tribox,nface,bboxbin,dxbin,dybin,dzbin,nbin,nbinx,trisize,ndim,btoel2,btoel1)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       :  memor_msh 
    use mod_memchk
    implicit none
    integer(ip),intent(in)        :: nbinx,nface,ndim,nbin
    integer(ip),intent(inout)     :: btoel2(nbin+1)
    integer(ip),pointer           :: btoel1(:)
    real(rp),intent(in)           :: trisize(nface),tribox(6,nface),bboxbin(ndim,2),dxbin,dybin,dzbin
    integer(ip)                   :: i,itri,j,k,iplace,nbinxy
    integer(ip)                   :: imin,jmin,kmin,imax,jmax,kmax,kbinx,jbinx,ibin,ntot
    real(rp)                      :: xmin,ymin,zmin,xmax,ymax,zmax,tol 
    integer(4)                 :: istat
    !
    !     Fill bins with triangulation 
    !

    nbinxy=nbinx*nbinx
    tol=1.0d-12


    !     Loop on the triangular faces to give the faces to the bins 

    !first count

    do itri=1,nface

       if(trisize(itri)>tol)then

          xmin=tribox(1,itri)
          ymin=tribox(2,itri)
          zmin=tribox(3,itri)  
          xmax=tribox(4,itri)
          ymax=tribox(5,itri)
          zmax=tribox(6,itri)  

          !Find the bin containing the bbox of the triangle

          imin=floor((xmin-bboxbin(1,1))/dxbin)+1
          jmin=floor((ymin-bboxbin(2,1))/dybin)+1
          kmin=floor((zmin-bboxbin(3,1))/dzbin)+1

          imax=floor((xmax-bboxbin(1,1))/dxbin)+1
          jmax=floor((ymax-bboxbin(2,1))/dybin)+1
          kmax=floor((zmax-bboxbin(3,1))/dzbin)+1


          !loop on the bins

          do k=kmin,kmax
             kbinx=(k-1)*nbinxy
             do j=jmin,jmax
                jbinx=(j-1)*nbinx+kbinx
                do i=imin,imax
                   ibin=i+jbinx+1_ip
                   btoel2(ibin)=btoel2(ibin)+1_ip
                enddo
             enddo
          enddo
       endif
    enddo

    !reshuffle

    btoel2(1)=1_ip

    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nbin,btoel2) &
    !$omp& private(i)

    do i=2,nbin+1_ip
       btoel2(i)=btoel2(i)+btoel2(i-1_ip)
    enddo

    !Allocate   

    ntot=btoel2(nbin+1)-1_ip
    allocate(btoel1(ntot),stat=istat)
    call memchk(zero,istat,memor_msh,'BTOEL1','cart',btoel1)

    !Store

    do itri=1,nface

       if(trisize(itri)>tol)then

          xmin=tribox(1,itri)
          ymin=tribox(2,itri)
          zmin=tribox(3,itri)  
          xmax=tribox(4,itri)
          ymax=tribox(5,itri)
          zmax=tribox(6,itri)  

          !Find the bin containing the bbox of the triangle

          imin=floor((xmin-bboxbin(1,1))/dxbin)+1_ip
          jmin=floor((ymin-bboxbin(2,1))/dybin)+1_ip
          kmin=floor((zmin-bboxbin(3,1))/dzbin)+1_ip

          imax=floor((xmax-bboxbin(1,1))/dxbin)+1_ip
          jmax=floor((ymax-bboxbin(2,1))/dybin)+1_ip
          kmax=floor((zmax-bboxbin(3,1))/dzbin)+1_ip

          !Loop on the bins

          do k=kmin,kmax
             kbinx=(k-1)*nbinxy
             do j=jmin,jmax
                jbinx=(j-1)*nbinx+kbinx
                do i=imin,imax
                   ibin=i+jbinx
                   iplace=btoel2(ibin)
                   btoel1(iplace)=itri
                   btoel2(ibin)=iplace+1_ip
                enddo
             enddo
          enddo
       endif
    enddo

    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(nbin,btoel2) &
    !$omp& private(i)

    do i=nbin+1,2,-1
       btoel2(i)=btoel2(i-1)
    enddo

    btoel2(1)=1_ip

  end subroutine tri2bin

  subroutine buildIni(nvoxx,nvoxy,nvoxz,bbox,dx,dy,dz,ndim,lcell,ncell)
    use def_kintyp, only       :  ip,rp,lg,cell
    use def_meshin, only          : memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)        :: nvoxx,nvoxy,nvoxz,ndim,ncell           
    real(rp),intent(in)           :: bbox(ndim,2),dx,dy,dz
    type(cell), pointer           :: lcell(:)
    integer(ip)                   :: i,j,k,np,kvoxx,jvoxx,ivox,ineigh
    integer(ip)                   :: nvoxxy
    integer(4)                 :: istat
    !
    !     This subroutine  build the lcell/face relation
    !
    nvoxxy=nvoxx*nvoxy
    allocate(lcell(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LCELL','buildIni',lcell)

    !$omp parallel do &
    !$omp& schedule(static) &
    !$omp& default(none) &
    !$omp& shared(ncell,lcell) &
    !$omp& private(ivox)

    do ivox=1,ncell
       lcell(ivox)%neigh(1)=0_ip
       lcell(ivox)%neigh(2)=0_ip
       lcell(ivox)%neigh(3)=0_ip
       lcell(ivox)%neigh(4)=0_ip
       lcell(ivox)%neigh(5)=0_ip
       lcell(ivox)%neigh(6)=0_ip
    enddo
    !
    !     Building lcell relationship, create lcells on the fly
    !
    do k=1,nvoxz
       kvoxx=nvoxxy*(k-1)
       do j=1,nvoxy
          jvoxx=kvoxx+(j-1)*nvoxx
          do i=1,nvoxx
             !
             !     Current voxel
             ! 
             ivox=jvoxx+i
             !write(*,*)ivox
             !
             !     Create lcells
             !
             lcell(ivox)%marked=-2_ip
             lcell(ivox)%level=1_ip
             lcell(ivox)%coor(1,1)=bbox(1,1)+(i-1)*dx
             lcell(ivox)%coor(2,1)=bbox(2,1)+(j-1)*dy
             lcell(ivox)%coor(3,1)=bbox(3,1)+(k-1)*dz
             lcell(ivox)%coor(1,2)=lcell(ivox)%coor(1,1)+dx
             lcell(ivox)%coor(2,2)=lcell(ivox)%coor(2,1)+dy
             lcell(ivox)%coor(3,2)=lcell(ivox)%coor(3,1)+dz
             call gtcartsize(lcell,ivox,ncell)
             !
             !     Loop on the 6 faces
             !
             !
             !      On x
             !
             if(lcell(ivox)%neigh(1)==0)then

                if(i>1)then           
                   ineigh=ivox-1_ip
                   lcell(ivox)%neigh(1)=ineigh
                   lcell(ineigh)%neigh(2)=ivox
                endif
             endif


             if(lcell(ivox)%neigh(2)==0)then

                if(i<nvoxx)then           
                   ineigh=ivox+1_ip
                   lcell(ivox)%neigh(2)=ineigh
                   lcell(ineigh)%neigh(1)=ivox
                endif
             endif
             !
             !       On y
             !
             if(lcell(ivox)%neigh(3)==0)then

                if(j>1)then           
                   ineigh=kvoxx+(j-2)*nvoxx+i
                   lcell(ivox)%neigh(3)=ineigh
                   lcell(ineigh)%neigh(4)=ivox
                endif
             endif

             if(lcell(ivox)%neigh(4)==0)then

                if(j<nvoxy)then           
                   ineigh=kvoxx+j*nvoxx+i
                   lcell(ivox)%neigh(4)=ineigh
                   lcell(ineigh)%neigh(3)=ivox
                endif

             endif
             !
             !      On z
             !
             if(lcell(ivox)%neigh(5)==0)then

                if(k>1)then           
                   ineigh=(k-2)*nvoxxy+(j-1)*nvoxx+i
                   lcell(ivox)%neigh(5)=ineigh
                   lcell(ineigh)%neigh(6)=ivox
                endif

             endif

             if(lcell(ivox)%neigh(6)==0)then

                if(k<nvoxz)then           
                   ineigh=k*nvoxxy+(j-1)*nvoxx+i
                   lcell(ivox)%neigh(6)=ineigh
                   lcell(ineigh)%neigh(5)=ivox
                endif

             endif

          enddo
       enddo
    enddo

  end subroutine buildIni

  subroutine mesh(btoel1,btoel2,tribox,bboxbin,dxbin,dybin,dzbin,nface,nbin,&
       nbinx,trisize,lface,coor,npoin,ndim,eltoel,rnofa,nnofa,nsour,rsour,rsgeo,niter,&
       lcell,ncell,rtol,rsuni,iinter,irefsurf)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ndim,nbin,nbinx,npoin,nnofa,nsour,iinter
    integer(ip),intent(in)        :: btoel1(*),btoel2(*)
    integer(ip),intent(in)        :: nface,irefsurf 
    integer(ip),intent(inout)     :: niter,ncell 
    integer(ip),intent(in)        :: lface(ndim,nface)
    type(cell), pointer           :: lcell(:) 
    real(rp),intent(in)           :: dxbin,dybin,dzbin,trisize(nface),coor(ndim,npoin)
    real(rp),intent(in)           :: tribox(6,nface),bboxbin(ndim,2),rsuni
    real(rp),intent(in)           :: rnofa(ndim,nface),rtol
    real(rp),pointer              :: rsour(:,:),rsgeo(:,:,:)
    integer(ip),intent(in)        :: eltoel(nnofa,nface)
    integer(ip)                   :: imin,jmin,kmin,imax,jmax
    integer(ip)                   :: ncell0,iter,itermax,marked,maxlevel,level,maxlev2,icell
    real(rp)                      :: xmin,ymin,zmin,xmax,ymax,zmax
    iter=1_ip
    itermax=50_ip

    maxlevel=8_ip
    maxlev2=int(maxlevel/2.0d+00)
    !
    !     Loop on refinement
    !
    level=0_ip

    do

       level=level+1
       !
       !     Mark the lcells to be refined
       ! 
       marked=0_ip
       ncell0=ncell 
       !
       !     if(level<maxlev2)then
       !     Mark the cells comparing the surface size to the cell size
       !


       !
       !     DBG
       !
       !xmin=1189.77
       !ymin=-6.3
       !zmin=56.24
       !xmax=1191.4
       !ymax=-4.73
       !zmax=58.0

       !do icell=1,ncell
       !   if(lcell(icell)%coor(1,2)<xmax .and. lcell(icell)%coor(2,2)<ymax .and. lcell(icell)%coor(3,2)<zmax .and.   &  
       !   lcell(icell)%coor(1,1)>xmin .and. lcell(icell)%coor(2,1)>ymin .and. lcell(icell)%coor(3,1)>zmin      )then
       !    write(*,*)'trouve 1, icell=',icell
       !   endif 
       !
       !enddo
       !
       !     Do we want to refine depending on the surface mesh ?
       !
       if(irefsurf==1)then

          call markCell(btoel1,btoel2,nbin,nbinx,ncell0,lcell,dxbin,dybin,dzbin,&
               bboxbin,trisize,marked,tribox,ndim,nface)

       endif

       !     else

       !        if(level==maxlev2)then
       !           do icell=1,ncell
       !              lcell(icell)%marked=-2_ip 
       !           enddo  
       !        endif

       !call markCellGeo(btoel1,btoel2,nbin,nbinx,ncell,lcell,dxbin,dybin,dzbin,bboxbin,marked,tribox,ndim,nface,eltoel,rnofa,nnofa,maxlevel)
       !
       !     Do we want to refine regarding some sources ?
       !
       if(nsour>0)then
          !
          !     Mark the cells comparing the source size to the cell size
          !
          call markCellSou(ncell,lcell,marked,ndim,nsour,rsour,rsgeo,iter,rtol,rsuni)

       endif

       !call markCellUni(ncell,lcell,marked,rsuni)

       !     endif

       write(*,*) marked,' lcells marked for refinement on',ncell0, 'total lcells for iter ',iter 
       !
       !     Do we have to do something?
       !
       if(marked==0)exit
       !
       !     DBG
       !
       do icell=1,ncell
          if(lcell(icell)%marked==1)then
             if(lcell(icell)%level/=iter)then
                write(*,*)'Error, elements marked not at the last level'
                stop          
             endif
          endif
       enddo
       !
       !     Sweep on the levels of refinement
       !
       call smooth(ncell,iter,ndim,lcell)
       ncell0=ncell 
       !
       !     debug
       !
       call dbgconform(lcell,ncell)
       !
       !     Refine the lcells 
       !
       call divideCell(ncell,ndim,lcell)
       !write(*,*)ncell
       !
       !     Refine the faces
       !
       call divideFace(ncell0,lcell)
       !
       !     debug
       !
       call dbgconform(lcell,ncell)
       !
       !     Compact arrays
       !
       call compact(ncell,lcell)
       !write(*,*)ncell
       !
       !     Debug
       !
       if(iter==itermax)then
          exit
       endif
       iter=iter+1

    enddo

    niter=iter
    !
    !     Do we want to intersect the cartesian mesh?
    ! 
    if(iinter==1)then
       !
       !     Now , run the full intersection test on marked cells
       !
       call realCut(btoel1, btoel2,nbin,nbinx,ncell,lcell,dxbin,dybin,dzbin,bboxbin,trisize,marked,tribox,&
            ndim,nface,lface,nnofa,coor,npoin)

    endif
    !Cut lcells

    ! call  cutCell(btoel1,btoel2,nbin,nbinx,ncell0,lcell,dxbin,dybin,dzbin,bboxbin,trisize,marked,tribox,lface,coor,npo,mpo,nface)

  end subroutine mesh

  subroutine divideCell(ncell,ndim,lcell)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only       :  memor_msh 
    implicit none
    integer(ip),intent(inout)     :: ncell
    integer(ip),intent(in)        :: ndim
    type(cell), pointer           :: lcell(:)           
    integer(ip)                   :: i,ncell1,ncell0,ivox,ineigh,level,ncell01
    integer(ip)                   :: imin,jmin,kmin,imax,jmax,kmax,kbinx,jbinx,ibin
    real(rp)                      :: lmin(ndim),lmax(ndim),c(ndim),c05,c20  


    c05=0.5d+00
    c20=2.0d+00

    ncell01=ncell
    do i=1,ncell01

       !Has the cell been marked for refinement?


       if(lcell(i)%marked==1)then

          !Give the pointer to the children  

          ncell0=ncell+1
          lcell(i)%marked=ncell0

          level=lcell(i)%level+1      

          lmin(1)=lcell(i)%coor(1,1)
          lmin(2)=lcell(i)%coor(2,1)
          lmin(3)=lcell(i)%coor(3,1)


          lmax(1)=lcell(i)%coor(1,2)
          lmax(2)=lcell(i)%coor(2,2)
          lmax(3)=lcell(i)%coor(3,2)

          c(1)=(lmin(1)+lmax(1))/c20 
          c(2)=(lmin(2)+lmax(2))/c20 
          c(3)=(lmin(3)+lmax(3))/c20 

          !Resize

          ncell1=ncell+8
          call memrea(ncell1,memor_msh,'LCELL','divideCell',lcell)

          !Create the  lcell children

          ncell=ncell+1
          lcell(ncell)%marked=-2
          lcell(ncell)%level=level
          lcell(ncell)%coor(1,1)=lmin(1)
          lcell(ncell)%coor(2,1)=lmin(2)
          lcell(ncell)%coor(3,1)=lmin(3)
          lcell(ncell)%coor(1,2)=c(1)
          lcell(ncell)%coor(2,2)=c(2)
          lcell(ncell)%coor(3,2)=c(3)
          call gtcartsize(lcell,ncell,ncell)
          ncell=ncell+1
          lcell(ncell)%marked=-2
          lcell(ncell)%level=level
          lcell(ncell)%coor(1,1)=lmin(1)
          lcell(ncell)%coor(2,1)=lmin(2)
          lcell(ncell)%coor(3,1)=c(3)
          lcell(ncell)%coor(1,2)=c(1)
          lcell(ncell)%coor(2,2)=c(2)
          lcell(ncell)%coor(3,2)=lmax(3)
          call gtcartsize(lcell,ncell,ncell)


          ncell=ncell+1
          lcell(ncell)%marked=-2
          lcell(ncell)%level=level
          lcell(ncell)%coor(1,1)=lmin(1)
          lcell(ncell)%coor(2,1)=c(2)
          lcell(ncell)%coor(3,1)=lmin(3)
          lcell(ncell)%coor(1,2)=c(1)
          lcell(ncell)%coor(2,2)=lmax(2)
          lcell(ncell)%coor(3,2)=c(3)
          call gtcartsize(lcell,ncell,ncell)
          ncell=ncell+1
          lcell(ncell)%marked=-2
          lcell(ncell)%level=level
          lcell(ncell)%coor(1,1)=lmin(1)
          lcell(ncell)%coor(2,1)=c(2)
          lcell(ncell)%coor(3,1)=c(3)
          lcell(ncell)%coor(1,2)=c(1)
          lcell(ncell)%coor(2,2)=lmax(2)
          lcell(ncell)%coor(3,2)=lmax(3)
          call gtcartsize(lcell,ncell,ncell)


          ncell=ncell+1
          lcell(ncell)%marked=-2
          lcell(ncell)%level=level
          lcell(ncell)%coor(1,1)=c(1)
          lcell(ncell)%coor(2,1)=lmin(2)
          lcell(ncell)%coor(3,1)=lmin(3)
          lcell(ncell)%coor(1,2)=lmax(1)
          lcell(ncell)%coor(2,2)=c(2)
          lcell(ncell)%coor(3,2)=c(3)
          call gtcartsize(lcell,ncell,ncell)
          ncell=ncell+1
          lcell(ncell)%marked=-2
          lcell(ncell)%level=level
          lcell(ncell)%coor(1,1)=c(1)
          lcell(ncell)%coor(2,1)=lmin(2)
          lcell(ncell)%coor(3,1)=c(3)
          lcell(ncell)%coor(1,2)=lmax(1)
          lcell(ncell)%coor(2,2)=c(2)
          lcell(ncell)%coor(3,2)=lmax(3)
          call gtcartsize(lcell,ncell,ncell)


          ncell=ncell+1
          lcell(ncell)%marked=-2
          lcell(ncell)%level=level
          lcell(ncell)%coor(1,1)=c(1)
          lcell(ncell)%coor(2,1)=c(2)
          lcell(ncell)%coor(3,1)=lmin(3)
          lcell(ncell)%coor(1,2)=lmax(1)
          lcell(ncell)%coor(2,2)=lmax(2)
          lcell(ncell)%coor(3,2)=c(3)
          call gtcartsize(lcell,ncell,ncell)
          ncell=ncell+1
          lcell(ncell)%marked=-2
          lcell(ncell)%level=level
          lcell(ncell)%coor(1,1)=c(1)
          lcell(ncell)%coor(2,1)=c(2)
          lcell(ncell)%coor(3,1)=c(3)
          lcell(ncell)%coor(1,2)=lmax(1)
          lcell(ncell)%coor(2,2)=lmax(2)
          lcell(ncell)%coor(3,2)=lmax(3)
          call gtcartsize(lcell,ncell,ncell)

          !Create the inner faces

          !In x

          ivox=ncell0
          ineigh=ncell0+4 
          lcell(ivox)%neigh(2)=ineigh
          lcell(ineigh)%neigh(1)=ivox

          ivox=ncell0+1
          ineigh=ncell0+5 
          lcell(ivox)%neigh(2)=ineigh
          lcell(ineigh)%neigh(1)=ivox

          ivox=ncell0+2
          ineigh=ncell0+6 
          lcell(ivox)%neigh(2)=ineigh
          lcell(ineigh)%neigh(1)=ivox

          ivox=ncell0+3
          ineigh=ncell0+7 
          lcell(ivox)%neigh(2)=ineigh
          lcell(ineigh)%neigh(1)=ivox

          !In y

          ivox=ncell0
          ineigh=ncell0+2 
          lcell(ivox)%neigh(4)=ineigh
          lcell(ineigh)%neigh(3)=ivox

          ivox=ncell0+4
          ineigh=ncell0+6 
          lcell(ivox)%neigh(4)=ineigh
          lcell(ineigh)%neigh(3)=ivox

          ivox=ncell0+1
          ineigh=ncell0+3 
          lcell(ivox)%neigh(4)=ineigh
          lcell(ineigh)%neigh(3)=ivox

          ivox=ncell0+5
          ineigh=ncell0+7 
          lcell(ivox)%neigh(4)=ineigh
          lcell(ineigh)%neigh(3)=ivox


          !In z

          ivox=ncell0
          ineigh=ncell0+1 
          lcell(ivox)%neigh(6)=ineigh
          lcell(ineigh)%neigh(5)=ivox

          ivox=ncell0+4
          ineigh=ncell0+5 
          lcell(ivox)%neigh(6)=ineigh
          lcell(ineigh)%neigh(5)=ivox

          ivox=ncell0+2
          ineigh=ncell0+3 
          lcell(ivox)%neigh(6)=ineigh
          lcell(ineigh)%neigh(5)=ivox

          ivox=ncell0+6
          ineigh=ncell0+7 
          lcell(ivox)%neigh(6)=ineigh
          lcell(ineigh)%neigh(5)=ivox

       endif
    enddo

  end subroutine divideCell

  subroutine divideCellSmoo(ncell,lcell,remem,ndim)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    type(cell)                    :: lcell(*)
    integer(ip),intent(in)        :: ndim,ncell 
    integer(ip),intent(in)        :: remem(ncell) 
    integer(ip)                   :: i,ncell1,ncell0,ivox,ineigh,level,iplace
    real(rp)                      :: lmin(ndim),lmax(ndim),c(ndim),c05

    c05=0.5d+00

    do i=1,ncell

       !Has the lcell been marked for refinement?

       if(remem(i)/=0)then

          !Give the pointer to the children  

          ncell0=remem(i)
          lcell(i)%marked=ncell0

          level=lcell(i)%level+1      

          lmin(1)=lcell(i)%coor(1,1)
          lmin(2)=lcell(i)%coor(2,1)
          lmin(3)=lcell(i)%coor(3,1)

          lmax(1)=lcell(i)%coor(1,2)
          lmax(2)=lcell(i)%coor(2,2)
          lmax(3)=lcell(i)%coor(3,2)

          c(1)=(lmin(1)+lmax(1))*c05
          c(2)=(lmin(2)+lmax(2))*c05
          c(3)=(lmin(3)+lmax(3))*c05

          !Create the  lcell children

          iplace=ncell0

          !Be carefull: if the marking criterion is not only verified 
          !by the highest level of cell, it will loose neighbor conformity                
          lcell(iplace)%marked=-2  
          !lcell(iplace)%marked=-1
          lcell(iplace)%level=level
          lcell(iplace)%coor(1,1)=lmin(1)
          lcell(iplace)%coor(2,1)=lmin(2)
          lcell(iplace)%coor(3,1)=lmin(3)
          lcell(iplace)%coor(1,2)=c(1)
          lcell(iplace)%coor(2,2)=c(2)
          lcell(iplace)%coor(3,2)=c(3)
          call gtcartsize(lcell,iplace,ncell)
          iplace=iplace+1
          lcell(iplace)%marked=-2
          !lcell(iplace)%marked=-1
          lcell(iplace)%level=level
          lcell(iplace)%coor(1,1)=lmin(1)
          lcell(iplace)%coor(2,1)=lmin(2)
          lcell(iplace)%coor(3,1)=c(3)
          lcell(iplace)%coor(1,2)=c(1)
          lcell(iplace)%coor(2,2)=c(2)
          lcell(iplace)%coor(3,2)=lmax(3)
          call gtcartsize(lcell,iplace,ncell)
          iplace=iplace+1


          lcell(iplace)%marked=-2
          !lcell(iplace)%marked=-1
          lcell(iplace)%level=level
          lcell(iplace)%coor(1,1)=lmin(1)
          lcell(iplace)%coor(2,1)=c(2)
          lcell(iplace)%coor(3,1)=lmin(3)
          lcell(iplace)%coor(1,2)=c(1)
          lcell(iplace)%coor(2,2)=lmax(2)
          lcell(iplace)%coor(3,2)=c(3)
          call gtcartsize(lcell,iplace,ncell)
          iplace=iplace+1
          lcell(iplace)%marked=-2
          !lcell(iplace)%marked=-1
          lcell(iplace)%level=level
          lcell(iplace)%coor(1,1)=lmin(1)
          lcell(iplace)%coor(2,1)=c(2)
          lcell(iplace)%coor(3,1)=c(3)
          lcell(iplace)%coor(1,2)=c(1)
          lcell(iplace)%coor(2,2)=lmax(2)
          lcell(iplace)%coor(3,2)=lmax(3)
          call gtcartsize(lcell,iplace,ncell)
          iplace=iplace+1


          lcell(iplace)%marked=-2
          !lcell(iplace)%marked=-1
          lcell(iplace)%level=level
          lcell(iplace)%coor(1,1)=c(1)
          lcell(iplace)%coor(2,1)=lmin(2)
          lcell(iplace)%coor(3,1)=lmin(3)
          lcell(iplace)%coor(1,2)=lmax(1)
          lcell(iplace)%coor(2,2)=c(2)
          lcell(iplace)%coor(3,2)=c(3)
          call gtcartsize(lcell,iplace,ncell)
          iplace=iplace+1
          lcell(iplace)%marked=-2
          !lcell(iplace)%marked=-1
          lcell(iplace)%level=level
          lcell(iplace)%coor(1,1)=c(1)
          lcell(iplace)%coor(2,1)=lmin(2)
          lcell(iplace)%coor(3,1)=c(3)
          lcell(iplace)%coor(1,2)=lmax(1)
          lcell(iplace)%coor(2,2)=c(2)
          lcell(iplace)%coor(3,2)=lmax(3)
          call gtcartsize(lcell,iplace,ncell)
          iplace=iplace+1


          lcell(iplace)%marked=-2
          !lcell(iplace)%marked=-1
          lcell(iplace)%level=level
          lcell(iplace)%coor(1,1)=c(1)
          lcell(iplace)%coor(2,1)=c(2)
          lcell(iplace)%coor(3,1)=lmin(3)
          lcell(iplace)%coor(1,2)=lmax(1)
          lcell(iplace)%coor(2,2)=lmax(2)
          lcell(iplace)%coor(3,2)=c(3)
          call gtcartsize(lcell,iplace,ncell)
          iplace=iplace+1
          lcell(iplace)%marked=-2
          !lcell(iplace)%marked=-1
          lcell(iplace)%level=level
          lcell(iplace)%coor(1,1)=c(1)
          lcell(iplace)%coor(2,1)=c(2)
          lcell(iplace)%coor(3,1)=c(3)
          lcell(iplace)%coor(1,2)=lmax(1)
          lcell(iplace)%coor(2,2)=lmax(2)
          lcell(iplace)%coor(3,2)=lmax(3)
          call gtcartsize(lcell,iplace,ncell)
          iplace=iplace+1

          !Create the inner faces

          !In x

          ivox=ncell0
          ineigh=ncell0+4 
          lcell(ivox)%neigh(2)=ineigh
          lcell(ineigh)%neigh(1)=ivox

          ivox=ncell0+1
          ineigh=ncell0+5 
          lcell(ivox)%neigh(2)=ineigh
          lcell(ineigh)%neigh(1)=ivox

          ivox=ncell0+2
          ineigh=ncell0+6 
          lcell(ivox)%neigh(2)=ineigh
          lcell(ineigh)%neigh(1)=ivox

          ivox=ncell0+3
          ineigh=ncell0+7 
          lcell(ivox)%neigh(2)=ineigh
          lcell(ineigh)%neigh(1)=ivox

          !In y

          ivox=ncell0
          ineigh=ncell0+2 
          lcell(ivox)%neigh(4)=ineigh
          lcell(ineigh)%neigh(3)=ivox

          ivox=ncell0+4
          ineigh=ncell0+6 
          lcell(ivox)%neigh(4)=ineigh
          lcell(ineigh)%neigh(3)=ivox

          ivox=ncell0+1
          ineigh=ncell0+3 
          lcell(ivox)%neigh(4)=ineigh
          lcell(ineigh)%neigh(3)=ivox

          ivox=ncell0+5
          ineigh=ncell0+7 
          lcell(ivox)%neigh(4)=ineigh
          lcell(ineigh)%neigh(3)=ivox


          !In z

          ivox=ncell0
          ineigh=ncell0+1 
          lcell(ivox)%neigh(6)=ineigh
          lcell(ineigh)%neigh(5)=ivox

          ivox=ncell0+4
          ineigh=ncell0+5 
          lcell(ivox)%neigh(6)=ineigh
          lcell(ineigh)%neigh(5)=ivox

          ivox=ncell0+2
          ineigh=ncell0+3 
          lcell(ivox)%neigh(6)=ineigh
          lcell(ineigh)%neigh(5)=ivox

          ivox=ncell0+6
          ineigh=ncell0+7 
          lcell(ivox)%neigh(6)=ineigh
          lcell(ineigh)%neigh(5)=ivox

       endif
    enddo

  end subroutine divideCellSmoo



  subroutine divideFace(ncell,lcell)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)            :: ncell
    type(cell)                    :: lcell(*)
    integer(ip)                   :: i,iorient,imark1,imark2,ineigh,inew,ineighold
    integer(ip)                   :: level,neighlevel,imarka,jpos,jcont,j,tab(3)

    do i=1,ncell

       !Has the cell been marked?

       if(lcell(i)%marked>0)then


          imark1=lcell(i)%marked
          level=lcell(i)%level

          !On x

          !Has this face been marked before?

          ineighold=lcell(i)%neigh(1)
          if(ineighold/=-2)then       

             lcell(i)%neigh(1)=-2

             !Is there a neighbor?

             if(ineighold==0)then
                lcell(imark1)%neigh(1)=0
                lcell(imark1+1)%neigh(1)=0
                lcell(imark1+2)%neigh(1)=0
                lcell(imark1+3)%neigh(1)=0
             else

                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Are the cells on the same level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      lcell(imark1)%neigh(1)=ineighold
                      lcell(imark1+1)%neigh(1)=ineighold
                      lcell(imark1+2)%neigh(1)=ineighold
                      lcell(imark1+3)%neigh(1)=ineighold
                      lcell(ineighold)%neigh(2)=imark1

                   else 

                      inew=imark1
                      ineigh=imark2+4
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      inew=imark1+1
                      ineigh=imark2+5
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      inew=imark1+2
                      ineigh=imark2+6
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      inew=imark1+3
                      ineigh=imark2+7
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew

                      lcell(ineighold)%neigh(2)=-2 

                   endif

                else
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      jpos=ineighold+1_ip
                      jcont=1_ip
                      do j=1,8
                         if(lcell(jpos)%neigh(2)==i)then
                            tab(jcont)=jpos
                            jcont=jcont+1_ip
                            if(jcont==4)exit 
                         endif
                         jpos=jpos+1  
                      enddo

                      inew=imark1
                      ineigh=ineighold
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      inew=imark1+1
                      ineigh=tab(1)
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      inew=imark1+2
                      ineigh=tab(2)
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      inew=imark1+3
                      ineigh=tab(3)
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew

                   else

                      if(level/=neighlevel)then
                         write(*,*)'Different levels 1'
                         stop
                      endif
                   endif
                endif
             endif
          endif

          ineighold=lcell(i)%neigh(2)
          if(ineighold/=-2)then

             lcell(i)%neigh(2)=-2

             if(ineighold==0)then
                lcell(imark1+4)%neigh(2)=0
                lcell(imark1+5)%neigh(2)=0
                lcell(imark1+6)%neigh(2)=0
                lcell(imark1+7)%neigh(2)=0
             else

                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Are the cells on the same level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      lcell(imark1+4)%neigh(2)=ineighold
                      lcell(imark1+5)%neigh(2)=ineighold
                      lcell(imark1+6)%neigh(2)=ineighold
                      lcell(imark1+7)%neigh(2)=ineighold
                      lcell(ineighold)%neigh(1)=imark1+4    

                   else 

                      inew=imark1+4
                      ineigh=imark2
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      inew=imark1+5
                      ineigh=imark2+1
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      inew=imark1+6
                      ineigh=imark2+2
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      inew=imark1+7
                      ineigh=imark2+3
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew

                      lcell(ineighold)%neigh(1)=-2    

                   endif

                else
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      jpos=ineighold+1_ip
                      jcont=1_ip
                      do j=1,8
                         if(lcell(jpos)%neigh(1)==i)then
                            tab(jcont)=jpos
                            jcont=jcont+1_ip
                            if(jcont==4)exit 
                         endif
                         jpos=jpos+1  
                      enddo

                      inew=imark1+4
                      ineigh=ineighold
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      inew=imark1+5
                      ineigh=tab(1)
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      inew=imark1+6
                      ineigh=tab(2)
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      inew=imark1+7
                      ineigh=tab(3)
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew

                   else

                      if(level/=neighlevel)then
                         write(*,*)'Different levels 2'
                         !   stop
                      endif

                   endif

                endif
             endif
          endif

          !On y

          ineighold=lcell(i)%neigh(3)

          if(ineighold/=-2)then

             lcell(i)%neigh(3)=-2

             !Only one side marked
             if(ineighold==0)then
                lcell(imark1)%neigh(3)=0
                lcell(imark1+1)%neigh(3)=0
                lcell(imark1+4)%neigh(3)=0
                lcell(imark1+5)%neigh(3)=0

             else 

                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Are the cells on the same level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      lcell(imark1)%neigh(3)=ineighold
                      lcell(imark1+1)%neigh(3)=ineighold
                      lcell(imark1+4)%neigh(3)=ineighold
                      lcell(imark1+5)%neigh(3)=ineighold
                      lcell(ineighold)%neigh(4)=imark1    

                   else  

                      inew=imark1
                      ineigh=imark2+2
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      inew=imark1+1
                      ineigh=imark2+3
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      inew=imark1+4
                      ineigh=imark2+6
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      inew=imark1+5
                      ineigh=imark2+7
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew

                      lcell(ineighold)%neigh(4)=-2

                   endif

                else
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      jpos=ineighold+1_ip
                      jcont=1_ip
                      do j=1,8
                         if(lcell(jpos)%neigh(4)==i)then
                            tab(jcont)=jpos
                            jcont=jcont+1_ip
                            if(jcont==4)exit 
                         endif
                         jpos=jpos+1  
                      enddo

                      inew=imark1
                      ineigh=ineighold
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      inew=imark1+1
                      ineigh=tab(1)
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      inew=imark1+4
                      ineigh=tab(2)
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      inew=imark1+5
                      ineigh=tab(3)
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew

                   else

                      if(level/=neighlevel)then
                         write(*,*)'Different levels 3'
                         !   stop
                      endif
                   endif
                endif
             endif
          endif

          ineighold=lcell(i)%neigh(4)
          if(ineighold/=-2)then

             lcell(i)%neigh(4)=-2

             !Only one side marked
             if(ineighold==0)then
                lcell(imark1+2)%neigh(4)=0
                lcell(imark1+3)%neigh(4)=0
                lcell(imark1+6)%neigh(4)=0
                lcell(imark1+7)%neigh(4)=0

             else              

                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Are the cells on the same level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then


                      lcell(imark1+2)%neigh(4)=ineighold
                      lcell(imark1+3)%neigh(4)=ineighold
                      lcell(imark1+6)%neigh(4)=ineighold
                      lcell(imark1+7)%neigh(4)=ineighold
                      lcell(ineighold)%neigh(3)=imark1+2    

                   else  

                      inew=imark1+2
                      ineigh=imark2
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      inew=imark1+3
                      ineigh=imark2+1
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      inew=imark1+6
                      ineigh=imark2+4
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      inew=imark1+7
                      ineigh=imark2+5
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew

                      lcell(ineighold)%neigh(3)=-2    

                   endif

                else
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      jpos=ineighold+1_ip
                      jcont=1_ip
                      do j=1,8
                         if(lcell(jpos)%neigh(3)==i)then
                            tab(jcont)=jpos
                            jcont=jcont+1_ip
                            if(jcont==4)exit 
                         endif
                         jpos=jpos+1  
                      enddo

                      inew=imark1+2
                      ineigh=ineighold
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      inew=imark1+3
                      ineigh=tab(1)
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      inew=imark1+6
                      ineigh=tab(2)
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      inew=imark1+7
                      ineigh=tab(3)
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew

                   else

                      if(level/=neighlevel)then
                         write(*,*)'Different levels 4'
                         !   stop
                      endif
                   endif


                endif
             endif
          endif

          !On z

          ineighold=lcell(i)%neigh(5)
          if(ineighold/=-2)then

             lcell(i)%neigh(5)=-2

             !Only one side marked

             if(ineighold==0)then
                lcell(imark1)%neigh(5)=0
                lcell(imark1+2)%neigh(5)=0
                lcell(imark1+4)%neigh(5)=0
                lcell(imark1+6)%neigh(5)=0

             else      

                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Are the cells on the same level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      lcell(imark1)%neigh(5)=ineighold
                      lcell(imark1+2)%neigh(5)=ineighold
                      lcell(imark1+4)%neigh(5)=ineighold
                      lcell(imark1+6)%neigh(5)=ineighold
                      lcell(ineighold)%neigh(6)=imark1    

                   else  

                      inew=imark1
                      ineigh=imark2+1
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      inew=imark1+2
                      ineigh=imark2+3
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      inew=imark1+4
                      ineigh=imark2+5
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      inew=imark1+6
                      ineigh=imark2+7
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew

                      lcell(ineighold)%neigh(6)=-2    

                   endif

                else

                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      jpos=ineighold+1_ip
                      jcont=1_ip
                      do j=1,8
                         if(lcell(jpos)%neigh(6)==i)then
                            tab(jcont)=jpos
                            jcont=jcont+1_ip
                            if(jcont==4)exit 
                         endif
                         jpos=jpos+1  
                      enddo
                      inew=imark1
                      ineigh=ineighold
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      inew=imark1+2
                      ineigh=tab(1)
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      inew=imark1+4
                      ineigh=tab(2)
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      inew=imark1+6
                      ineigh=tab(3)
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew

                   else

                      if(level/=neighlevel)then
                         write(*,*)'Different levels 5'
                         !   stop
                      endif
                   endif
                endif
             endif
          endif

          ineighold=lcell(i)%neigh(6)
          if(ineighold/=-2)then

             lcell(i)%neigh(6)=-2

             !Only one side marked
             if(ineighold==0)then
                lcell(imark1+1)%neigh(6)=0
                lcell(imark1+3)%neigh(6)=0
                lcell(imark1+5)%neigh(6)=0
                lcell(imark1+7)%neigh(6)=0

             else 

                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Are the cells on the same level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      lcell(imark1+1)%neigh(6)=ineighold
                      lcell(imark1+3)%neigh(6)=ineighold
                      lcell(imark1+5)%neigh(6)=ineighold
                      lcell(imark1+7)%neigh(6)=ineighold
                      lcell(ineighold)%neigh(5)=imark1+1    

                   else  

                      inew=imark1+1
                      ineigh=imark2
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      inew=imark1+3
                      ineigh=imark2+2
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      inew=imark1+5
                      ineigh=imark2+4
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      inew=imark1+7
                      ineigh=imark2+6
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew

                      lcell(ineighold)%neigh(5)=-2

                   endif

                else
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(imark2<0)then

                      jpos=ineighold+1_ip
                      jcont=1_ip
                      do j=1,8
                         if(lcell(jpos)%neigh(5)==i)then
                            tab(jcont)=jpos
                            jcont=jcont+1_ip
                            if(jcont==4)exit 
                         endif
                         jpos=jpos+1  
                      enddo

                      inew=imark1+1
                      ineigh=ineighold
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      inew=imark1+3
                      ineigh=tab(1)
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      inew=imark1+5
                      ineigh=tab(2)
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      inew=imark1+7
                      ineigh=tab(3)
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew

                   else               

                      if(level/=neighlevel)then
                         write(*,*)'Different levels 6'
                         !   stop
                      endif
                   endif
                endif
             endif
          endif
       endif
    enddo

  end subroutine divideFace

  subroutine divideFaceSmoo(ncell,lcell,remem)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ncell
    integer(ip)                   :: remem(ncell) 
    type(cell)                    :: lcell(*)
    integer(ip)                   :: i,iorient,imark1,imark2,ineigh,inew,ineighold,irem,j
    integer(ip)                   :: level,neighlevel,imarka,jpos,jcont,tab(3)

    !
    !     This sub takes care of the neigbor update.
    !     It handles the cases where neighbor cells may have to be refined
    !     at the same time.
    !     Cells are processed from the lowest level to the highest, so 
    !     the neighbor pointer is the one from the lowest cell.
    !


    do i=1,ncell
       !
       !     Has the lcell been marked?
       !
       if(remem(i)/=0)then


          imark1=lcell(i)%marked
          level=lcell(i)%level
          !
          !     On x
          !
          !     Has this face been marked before?
          !
          ineighold=lcell(i)%neigh(1)
          if(ineighold/=-2)then    

             lcell(i)%neigh(1)=-2_ip
             !
             !     Is there a neighbor?
             !
             if(ineighold==0)then

                lcell(imark1)%neigh(1)=0_ip
                lcell(imark1+1)%neigh(1)=0_ip
                lcell(imark1+2)%neigh(1)=0_ip
                lcell(imark1+3)%neigh(1)=0_ip

             else 

                irem=remem(ineighold)
                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Is there a difference of level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(irem==0)then

                      lcell(imark1)%neigh(1)=ineighold
                      lcell(imark1+1)%neigh(1)=ineighold
                      lcell(imark1+2)%neigh(1)=ineighold
                      lcell(imark1+3)%neigh(1)=ineighold
                      lcell(ineighold)%neigh(2)=imark1

                   else 

                      inew=imark1
                      ineigh=imark2+4_ip
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      inew=imark1+1_ip
                      ineigh=imark2+5_ip
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      inew=imark1+2_ip
                      ineigh=imark2+6_ip
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      inew=imark1+3_ip
                      ineigh=imark2+7_ip
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew

                      lcell(ineighold)%neigh(2)=-2_ip 

                   endif

                else
                   !
                   !     There is a level difference
                   !     Cell i has one level lower than cell neighlevel
                   ! 

                   jpos=ineighold+1_ip
                   jcont=1_ip
                   do j=1,8
                      if(lcell(jpos)%neigh(2)==i)then
                         tab(jcont)=jpos
                         jcont=jcont+1_ip
                         if(jcont==4)exit 
                      endif
                      jpos=jpos+1  
                   enddo
                   !
                   !     Must check if each neighbor will be refined
                   !
                   inew=imark1
                   if(irem==0)then
                      ineigh=ineighold
                      lcell(inew)%neigh(1_ip)=ineigh
                      lcell(ineigh)%neigh(2_ip)=inew
                   else 
                      ineigh=imark2+4_ip
                      lcell(ineighold)%neigh(2)=-2_ip
                      lcell(inew)%neigh(1_ip)=ineigh
                      lcell(ineigh)%neigh(2_ip)=inew
                      ineigh=imark2+5_ip
                      lcell(ineigh)%neigh(2_ip)=inew
                      ineigh=imark2+6_ip
                      lcell(ineigh)%neigh(2_ip)=inew
                      ineigh=imark2+7_ip
                      lcell(ineigh)%neigh(2_ip)=inew
                   endif

                   irem=remem(tab(1))
                   imark2=lcell(tab(1))%marked
                   inew=imark1+1_ip
                   if(irem==0)then
                      ineigh=tab(1)                            
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                   else 
                      ineigh=imark2+4_ip                          
                      lcell(tab(1))%neigh(2)=-2_ip
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      ineigh=imark2+5_ip                          
                      lcell(ineigh)%neigh(2)=inew
                      ineigh=imark2+6_ip                          
                      lcell(ineigh)%neigh(2)=inew
                      ineigh=imark2+7_ip                          
                      lcell(ineigh)%neigh(2)=inew
                   endif

                   irem=remem(tab(2))
                   imark2=lcell(tab(2))%marked
                   inew=imark1+2_ip
                   if(irem==0)then
                      ineigh=tab(2)                            
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                   else 
                      ineigh=imark2+4_ip                          
                      lcell(tab(2))%neigh(2)=-2_ip
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      ineigh=imark2+5_ip                          
                      lcell(ineigh)%neigh(2)=inew
                      ineigh=imark2+6_ip                          
                      lcell(ineigh)%neigh(2)=inew
                      ineigh=imark2+7_ip                          
                      lcell(ineigh)%neigh(2)=inew
                   endif

                   irem=remem(tab(3))
                   imark2=lcell(tab(3))%marked
                   inew=imark1+3_ip
                   if(irem==0)then
                      ineigh=tab(3)                            
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                   else 
                      ineigh=imark2+4_ip                          
                      lcell(tab(3))%neigh(2)=-2_ip
                      lcell(inew)%neigh(1)=ineigh
                      lcell(ineigh)%neigh(2)=inew
                      ineigh=imark2+5_ip                          
                      lcell(ineigh)%neigh(2)=inew
                      ineigh=imark2+6_ip                          
                      lcell(ineigh)%neigh(2)=inew
                      ineigh=imark2+7_ip                          
                      lcell(ineigh)%neigh(2)=inew
                   endif
                endif
             endif
          endif

          ineighold=lcell(i)%neigh(2)
          if(ineighold/=-2)then

             lcell(i)%neigh(2)=-2_ip

             if(ineighold==0)then

                lcell(imark1+4)%neigh(2)=0_ip
                lcell(imark1+5)%neigh(2)=0_ip
                lcell(imark1+6)%neigh(2)=0_ip
                lcell(imark1+7)%neigh(2)=0_ip

             else 

                irem=remem(ineighold)
                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Is there a difference of level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(irem==0)then

                      lcell(imark1+4)%neigh(2)=ineighold
                      lcell(imark1+5)%neigh(2)=ineighold
                      lcell(imark1+6)%neigh(2)=ineighold
                      lcell(imark1+7)%neigh(2)=ineighold
                      lcell(ineighold)%neigh(1)=imark1+4_ip    

                   else 

                      inew=imark1+4_ip
                      ineigh=imark2
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      inew=imark1+5_ip
                      ineigh=imark2+1_ip
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      inew=imark1+6_ip
                      ineigh=imark2+2_ip
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      inew=imark1+7_ip
                      ineigh=imark2+3_ip
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew

                      lcell(ineighold)%neigh(1)=-2_ip  

                   endif

                else
                   !
                   !     There is a level difference
                   !     Cell i has one level lower than cell neighlevel
                   ! 

                   jpos=ineighold+1_ip
                   jcont=1_ip
                   do j=1,8
                      if(lcell(jpos)%neigh(1)==i)then
                         tab(jcont)=jpos
                         jcont=jcont+1_ip
                         if(jcont==4)exit 
                      endif
                      jpos=jpos+1_ip  
                   enddo
                   !
                   !     Must check if each neighbor will be refined
                   !
                   inew=imark1+4_ip
                   if(irem==0)then
                      ineigh=ineighold
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                   else
                      ineigh=imark2
                      lcell(ineighold)%neigh(1)=-2_ip
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+1_ip
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+2_ip
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+3_ip
                      lcell(ineigh)%neigh(1)=inew
                   endif

                   irem=remem(tab(1))
                   imark2=lcell(tab(1))%marked
                   inew=imark1+5_ip
                   if(irem==0)then
                      ineigh=tab(1)   
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                   else
                      ineigh=imark2
                      lcell(tab(1))%neigh(1)=-2_ip
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+1_ip
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+2_ip
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+3_ip
                      lcell(ineigh)%neigh(1)=inew

                   endif

                   irem=remem(tab(2))
                   imark2=lcell(tab(2))%marked
                   inew=imark1+6_ip
                   if(irem==0)then
                      ineigh=tab(2)   
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                   else
                      ineigh=imark2
                      lcell(tab(2))%neigh(1)=-2_ip
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+1_ip
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+2_ip
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+3_ip
                      lcell(ineigh)%neigh(1)=inew

                   endif

                   irem=remem(tab(3))
                   imark2=lcell(tab(3))%marked
                   inew=imark1+7_ip
                   if(irem==0)then
                      ineigh=tab(3)   
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                   else
                      ineigh=imark2
                      lcell(tab(3))%neigh(1)=-2_ip
                      lcell(inew)%neigh(2)=ineigh
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+1_ip
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+2_ip
                      lcell(ineigh)%neigh(1)=inew
                      ineigh=imark2+3_ip
                      lcell(ineigh)%neigh(1)=inew

                   endif
                endif
             endif
          endif

          !On y

          ineighold=lcell(i)%neigh(3)

          if(ineighold/=-2)then

             lcell(i)%neigh(3)=-2_ip

             !Only one side marked
             if(ineighold==0)then
                lcell(imark1)%neigh(3)=0_ip
                lcell(imark1+1)%neigh(3)=0_ip
                lcell(imark1+4)%neigh(3)=0_ip
                lcell(imark1+5)%neigh(3)=0_ip

             else 

                irem=remem(ineighold)
                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Is there a difference of level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(irem==0)then

                      lcell(imark1)%neigh(3)=ineighold
                      lcell(imark1+1)%neigh(3)=ineighold
                      lcell(imark1+4)%neigh(3)=ineighold
                      lcell(imark1+5)%neigh(3)=ineighold
                      lcell(ineighold)%neigh(4)=imark1

                   else 

                      inew=imark1
                      ineigh=imark2+2_ip
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      inew=imark1+1_ip
                      ineigh=imark2+3_ip
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      inew=imark1+4_ip
                      ineigh=imark2+6_ip
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      inew=imark1+5_ip
                      ineigh=imark2+7_ip
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew

                      lcell(ineighold)%neigh(4)=-2_ip

                   endif

                else
                   !
                   !     There is a level difference
                   !     Cell i has one level lower than cell neighlevel
                   ! 

                   jpos=ineighold+1_ip
                   jcont=1_ip
                   do j=1,8
                      if(lcell(jpos)%neigh(4)==i)then
                         tab(jcont)=jpos
                         jcont=jcont+1_ip
                         if(jcont==4)exit 
                      endif
                      jpos=jpos+1_ip  
                   enddo
                   !
                   !     Must check if each neighbor will be refined
                   !
                   inew=imark1
                   if(irem==0)then
                      ineigh=ineighold                            
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                   else
                      ineigh=imark2+2_ip 
                      lcell(ineighold)%neigh(4)=-2_ip
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+3_ip 
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+6_ip 
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+7_ip 
                      lcell(ineigh)%neigh(4)=inew
                   endif

                   irem=remem(tab(1))
                   imark2=lcell(tab(1))%marked
                   inew=imark1+1_ip
                   if(irem==0)then
                      ineigh=tab(1)                            
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                   else
                      ineigh=imark2+2_ip 
                      lcell(tab(1))%neigh(4)=-2_ip
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+3_ip 
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+6_ip 
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+7_ip 
                      lcell(ineigh)%neigh(4)=inew

                   endif

                   irem=remem(tab(2))
                   imark2=lcell(tab(2))%marked
                   inew=imark1+4_ip
                   if(irem==0)then
                      ineigh=tab(2)                            
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                   else
                      ineigh=imark2+2_ip 
                      lcell(tab(2))%neigh(4)=-2_ip
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+3_ip 
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+6_ip 
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+7_ip 
                      lcell(ineigh)%neigh(4)=inew

                   endif

                   irem=remem(tab(3))
                   imark2=lcell(tab(3))%marked
                   inew=imark1+5_ip
                   if(irem==0)then
                      ineigh=tab(3)                            
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                   else
                      ineigh=imark2+2_ip 
                      lcell(tab(3))%neigh(4)=-2_ip
                      lcell(inew)%neigh(3)=ineigh
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+3_ip 
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+6_ip 
                      lcell(ineigh)%neigh(4)=inew
                      ineigh=imark2+7_ip 
                      lcell(ineigh)%neigh(4)=inew

                   endif
                endif
             endif
          endif

          ineighold=lcell(i)%neigh(4)
          if(ineighold/=-2)then

             lcell(i)%neigh(4)=-2_ip

             !Only one side marked
             if(ineighold==0)then
                lcell(imark1+2)%neigh(4)=0_ip
                lcell(imark1+3)%neigh(4)=0_ip
                lcell(imark1+6)%neigh(4)=0_ip
                lcell(imark1+7)%neigh(4)=0_ip

             else              

                irem=remem(ineighold)
                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Is there a difference of level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(irem==0)then

                      lcell(imark1+2)%neigh(4)=ineighold
                      lcell(imark1+3)%neigh(4)=ineighold
                      lcell(imark1+6)%neigh(4)=ineighold
                      lcell(imark1+7)%neigh(4)=ineighold
                      lcell(ineighold)%neigh(3)=imark1+2_ip    

                   else 

                      inew=imark1+2_ip
                      ineigh=imark2
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      inew=imark1+3_ip
                      ineigh=imark2+1_ip
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      inew=imark1+6_ip
                      ineigh=imark2+4_ip
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      inew=imark1+7_ip
                      ineigh=imark2+5_ip
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew

                      lcell(ineighold)%neigh(3)=-2_ip

                   endif
                else
                   !
                   !     There is a level difference
                   !     Cell i has one level lower than cell neighlevel
                   ! 

                   jpos=ineighold+1_ip
                   jcont=1_ip
                   do j=1,8
                      if(lcell(jpos)%neigh(3)==i)then
                         tab(jcont)=jpos
                         jcont=jcont+1
                         if(jcont==4)exit 
                      endif
                      jpos=jpos+1_ip  
                   enddo
                   !
                   !     Must check if each neighbor will be refined
                   !
                   inew=imark1+2_ip
                   if(irem==0)then           
                      ineigh=ineighold    
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                   else
                      ineigh=imark2  
                      lcell(ineighold)%neigh(3)=-2_ip
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+1_ip  
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+4_ip  
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+5_ip  
                      lcell(ineigh)%neigh(3)=inew
                   endif

                   irem=remem(tab(1))
                   imark2=lcell(tab(1))%marked
                   inew=imark1+3_ip
                   if(irem==0)then           
                      ineigh=tab(1)    
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                   else
                      ineigh=imark2  
                      lcell(tab(1))%neigh(3)=-2_ip
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+1_ip  
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+4_ip  
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+5_ip  
                      lcell(ineigh)%neigh(3)=inew

                   endif

                   irem=remem(tab(2))
                   imark2=lcell(tab(2))%marked
                   inew=imark1+6_ip
                   if(irem==0)then           
                      ineigh=tab(2)    
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                   else
                      ineigh=imark2  
                      lcell(tab(2))%neigh(3)=-2_ip
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+1_ip  
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+4_ip  
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+5_ip  
                      lcell(ineigh)%neigh(3)=inew

                   endif

                   irem=remem(tab(3))
                   imark2=lcell(tab(3))%marked
                   inew=imark1+7_ip
                   if(irem==0)then           
                      ineigh=tab(3)    
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                   else
                      ineigh=imark2  
                      lcell(tab(3))%neigh(3)=-2_ip
                      lcell(inew)%neigh(4)=ineigh
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+1_ip  
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+4_ip  
                      lcell(ineigh)%neigh(3)=inew
                      ineigh=imark2+5_ip  
                      lcell(ineigh)%neigh(3)=inew

                   endif
                endif
             endif
          endif

          ! On z


          ineighold=lcell(i)%neigh(5)
          if(ineighold/=-2)then

             lcell(i)%neigh(5)=-2_ip

             !Only one side marked

             if(ineighold==0)then
                lcell(imark1)%neigh(5)=0_ip
                lcell(imark1+2)%neigh(5)=0_ip
                lcell(imark1+4)%neigh(5)=0_ip
                lcell(imark1+6)%neigh(5)=0_ip

             else      

                irem=remem(ineighold)
                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level 
                !
                !     Is there a difference of level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(irem==0)then

                      lcell(imark1)%neigh(5)=ineighold
                      lcell(imark1+2)%neigh(5)=ineighold
                      lcell(imark1+4)%neigh(5)=ineighold
                      lcell(imark1+6)%neigh(5)=ineighold
                      lcell(ineighold)%neigh(6)=imark1    

                   else

                      inew=imark1
                      ineigh=imark2+1_ip
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      inew=imark1+2_ip
                      ineigh=imark2+3_ip
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      inew=imark1+4_ip
                      ineigh=imark2+5_ip
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      inew=imark1+6_ip
                      ineigh=imark2+7_ip
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew

                      lcell(ineighold)%neigh(6)=-2_ip    

                   endif
                else
                   !
                   !     There is a level difference
                   !     Cell i has one level lower than cell neighlevel
                   ! 

                   jpos=ineighold+1_ip
                   jcont=1_ip
                   do j=1,8
                      if(lcell(jpos)%neigh(6)==i)then
                         tab(jcont)=jpos
                         jcont=jcont+1
                         if(jcont==4)exit 
                      endif
                      jpos=jpos+1_ip  
                   enddo
                   !
                   !     Must check if each neighbor will be refined
                   !
                   inew=imark1
                   if(irem==0)then
                      ineigh=ineighold
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                   else
                      ineigh=imark2+1_ip                    
                      lcell(ineighold)%neigh(6)=-2_ip
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+3_ip                    
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+5_ip                    
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+7_ip                    
                      lcell(ineigh)%neigh(6)=inew
                   endif

                   irem=remem(tab(1))
                   imark2=lcell(tab(1))%marked
                   inew=imark1+2_ip
                   if(irem==0)then
                      ineigh=tab(1)
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                   else
                      ineigh=imark2+1_ip                    
                      lcell(tab(1))%neigh(6)=-2_ip
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+3_ip                    
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+5_ip                    
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+7_ip                    
                      lcell(ineigh)%neigh(6)=inew
                   endif

                   irem=remem(tab(2))
                   imark2=lcell(tab(2))%marked
                   inew=imark1+4_ip
                   if(irem==0)then
                      ineigh=tab(2)
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                   else
                      ineigh=imark2+1_ip                    
                      lcell(tab(2))%neigh(6)=-2_ip
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+3_ip                    
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+5_ip                    
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+7_ip                    
                      lcell(ineigh)%neigh(6)=inew
                   endif

                   irem=remem(tab(3))
                   imark2=lcell(tab(3))%marked
                   inew=imark1+6_ip
                   if(irem==0)then
                      ineigh=tab(3)
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                   else
                      ineigh=imark2+1_ip                    
                      lcell(tab(3))%neigh(6)=-2_ip
                      lcell(inew)%neigh(5)=ineigh
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+3_ip                    
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+5_ip                    
                      lcell(ineigh)%neigh(6)=inew
                      ineigh=imark2+7_ip                    
                      lcell(ineigh)%neigh(6)=inew
                   endif
                endif
             endif
          endif

          ineighold=lcell(i)%neigh(6)
          if(ineighold/=-2)then

             lcell(i)%neigh(6)=-2_ip

             !Only one side marked
             if(ineighold==0)then
                lcell(imark1+1)%neigh(6)=0_ip
                lcell(imark1+3)%neigh(6)=0_ip
                lcell(imark1+5)%neigh(6)=0_ip
                lcell(imark1+7)%neigh(6)=0_ip

             else 

                irem=remem(ineighold)
                imark2=lcell(ineighold)%marked
                neighlevel=lcell(ineighold)%level
                !
                !     Is there a difference of level
                !
                if(level==neighlevel)then
                   !
                   !     Has the neighbor been marked for refinement?           
                   !
                   if(irem==0)then

                      lcell(imark1+1)%neigh(6)=ineighold
                      lcell(imark1+3)%neigh(6)=ineighold
                      lcell(imark1+5)%neigh(6)=ineighold
                      lcell(imark1+7)%neigh(6)=ineighold
                      lcell(ineighold)%neigh(5)=imark1+1_ip    

                   else

                      inew=imark1+1_ip
                      ineigh=imark2
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      inew=imark1+3_ip
                      ineigh=imark2+2_ip
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      inew=imark1+5_ip
                      ineigh=imark2+4_ip
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      inew=imark1+7_ip
                      ineigh=imark2+6_ip
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew

                      lcell(ineighold)%neigh(5)=-2_ip

                   endif
                else
                   !
                   !     There is a level difference
                   !     Cell i has one level lower than cell neighlevel
                   ! 

                   jpos=ineighold+1_ip
                   jcont=1_ip
                   do j=1,8
                      if(lcell(jpos)%neigh(5)==i)then
                         tab(jcont)=jpos
                         jcont=jcont+1
                         if(jcont==4)exit 
                      endif
                      jpos=jpos+1_ip  
                   enddo
                   !
                   !     Must check if each neighbor will be refined
                   ! 
                   inew=imark1+1_ip
                   if(irem==0)then
                      ineigh=ineighold
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                   else
                      ineigh=imark2
                      lcell(ineighold)%neigh(5)=-2_ip
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+2_ip
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+4_ip
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+6_ip
                      lcell(ineigh)%neigh(5)=inew
                   endif

                   irem=remem(tab(1))
                   imark2=lcell(tab(1))%marked
                   inew=imark1+3_ip
                   if(irem==0)then
                      ineigh=tab(1)
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                   else
                      ineigh=imark2
                      lcell(tab(1))%neigh(5)=-2_ip
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+2_ip
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+4_ip
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+6_ip
                      lcell(ineigh)%neigh(5)=inew
                   endif

                   irem=remem(tab(2))
                   imark2=lcell(tab(2))%marked
                   inew=imark1+5_ip
                   if(irem==0)then
                      ineigh=tab(2)
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                   else
                      ineigh=imark2
                      lcell(tab(2))%neigh(5)=-2_ip
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+2_ip
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+4_ip
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+6_ip
                      lcell(ineigh)%neigh(5)=inew
                   endif

                   irem=remem(tab(3))
                   imark2=lcell(tab(3))%marked
                   inew=imark1+7_ip
                   if(irem==0)then
                      ineigh=tab(3)
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                   else
                      ineigh=imark2
                      lcell(tab(3))%neigh(5)=-2_ip
                      lcell(inew)%neigh(6)=ineigh
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+2_ip
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+4_ip
                      lcell(ineigh)%neigh(5)=inew
                      ineigh=imark2+6_ip
                      lcell(ineigh)%neigh(5)=inew
                   endif
                endif
             endif
          endif
       endif
    enddo
  end subroutine divideFaceSmoo

  subroutine compact(ncell, lcell)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(inout)     :: ncell
    type(cell),intent(inout)      :: lcell(*)
    integer(ip)                   :: i,ipoint,icell,ineigh
    integer(ip),pointer           :: renum(:)
    integer(4)                 :: istat

    allocate(renum(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'RENUM','compact',renum)

    !
    !    First loop, count the non marked cells
    !
    ipoint=0
    do i=1,ncell
       if(lcell(i)%marked<0)then
          ipoint=ipoint+1
          renum(i)=ipoint
       endif
    enddo

    !
    !     Second loop, renumber the non marked cells
    !     The second loop is needed for the neighbors
    ! 

    do i=1,ncell
       if(lcell(i)%marked<0)then
          icell=renum(i)
          lcell(icell)%marked=lcell(i)%marked
          lcell(icell)%level=lcell(i)%level
          lcell(icell)%coor(1,1)=lcell(i)%coor(1,1)
          lcell(icell)%coor(2,1)=lcell(i)%coor(2,1)
          lcell(icell)%coor(3,1)=lcell(i)%coor(3,1)
          lcell(icell)%coor(1,2)=lcell(i)%coor(1,2)
          lcell(icell)%coor(2,2)=lcell(i)%coor(2,2)
          lcell(icell)%coor(3,2)=lcell(i)%coor(3,2)
          lcell(icell)%rsize=lcell(i)%rsize

          ineigh=lcell(i)%neigh(1)
          if(ineigh==0)then
             lcell(icell)%neigh(1)=0 
          else 
             lcell(icell)%neigh(1)=renum(ineigh) 
          endif

          ineigh=lcell(i)%neigh(2)
          if(ineigh==0)then
             lcell(icell)%neigh(2)=0 
          else 
             lcell(icell)%neigh(2)=renum(ineigh)
          endif

          ineigh=lcell(i)%neigh(3)
          if(ineigh==0)then
             lcell(icell)%neigh(3)=0 
          else 
             lcell(icell)%neigh(3)=renum(ineigh)
          endif

          ineigh=lcell(i)%neigh(4)
          if(ineigh==0)then
             lcell(icell)%neigh(4)=0 
          else 
             lcell(icell)%neigh(4)=renum(ineigh)
          endif

          ineigh=lcell(i)%neigh(5)
          if(ineigh==0)then
             lcell(icell)%neigh(5)=0 
          else 
             lcell(icell)%neigh(5)=renum(ineigh)
          endif

          ineigh=lcell(i)%neigh(6)
          if(ineigh==0)then
             lcell(icell)%neigh(6)=0 
          else 
             lcell(icell)%neigh(6)=renum(ineigh)
          endif

       endif
    enddo


    ncell=ipoint

    call memchk(2_ip,istat,memor_msh,'RENUM','compact',renum)
    deallocate(renum,stat=istat)
    if(istat/=0) call memerr(2_ip,'RENUM','compact',0_ip)

  end subroutine compact

  subroutine compactSmoo(ncell, lcell, remem)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(inout)     :: ncell 
    integer(ip),intent(inout)     :: remem(ncell) 
    type(cell)                    :: lcell(*)
    integer(ip)                   :: i,ipoint,icell,ineigh


    !First loop, counting

    ipoint=0
    do i=1,ncell
       if(remem(i)==0)then
          ipoint=ipoint+1
          remem(i)=ipoint
       else
          remem(i)=0 
       endif
    enddo

    do i=1,ncell
       if(remem(i)>0)then
          icell=remem(i)
          lcell(icell)%marked=lcell(i)%marked
          lcell(icell)%level=lcell(i)%level
          lcell(icell)%coor(1,1)=lcell(i)%coor(1,1)
          lcell(icell)%coor(2,1)=lcell(i)%coor(2,1)
          lcell(icell)%coor(3,1)=lcell(i)%coor(3,1)
          lcell(icell)%coor(1,2)=lcell(i)%coor(1,2)
          lcell(icell)%coor(2,2)=lcell(i)%coor(2,2)
          lcell(icell)%coor(3,2)=lcell(i)%coor(3,2)
          lcell(icell)%rsize=lcell(i)%rsize

          ineigh=lcell(i)%neigh(1)
          if(ineigh==0)then
             lcell(icell)%neigh(1)=0 
          else 
             lcell(icell)%neigh(1)=remem(ineigh)
          endif

          ineigh=lcell(i)%neigh(2)
          if(ineigh==0)then
             lcell(icell)%neigh(2)=0 
          else 
             lcell(icell)%neigh(2)=remem(ineigh)  
          endif

          ineigh=lcell(i)%neigh(3)
          if(ineigh==0)then
             lcell(icell)%neigh(3)=0 
          else 
             lcell(icell)%neigh(3)=remem(ineigh)
          endif

          ineigh=lcell(i)%neigh(4)
          if(ineigh==0)then
             lcell(icell)%neigh(4)=0 
          else 
             lcell(icell)%neigh(4)=remem(ineigh)
          endif

          ineigh=lcell(i)%neigh(5)
          if(ineigh==0)then
             lcell(icell)%neigh(5)=0 
          else 
             lcell(icell)%neigh(5)=remem(ineigh)
          endif

          ineigh=lcell(i)%neigh(6)
          if(ineigh==0)then
             lcell(icell)%neigh(6)=0 
          else 
             lcell(icell)%neigh(6)=remem(ineigh)
          endif

       endif
    enddo

    ncell=ipoint

  end subroutine compactSmoo

  subroutine markCell(btoel1, btoel2,nbin,nbinx,ncell,lcell,dxbin,dybin,dzbin,bboxbin,&
       trisize,marked,tribox,ndim,nface)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: nbinx,ndim,ncell,nbin,nface
    integer(ip),intent(in)        :: btoel1(*),btoel2(*)
    real(rp),intent(in)           :: dxbin,dybin,dzbin,trisize(nface),tribox(6,nface),bboxbin(ndim,2)
    integer(ip),intent(inout)     :: marked
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: imin,jmin,kmin,imax,jmax,kmax,kbinx,jbinx,ibin,isto0,isto1
    integer(ip)                   :: icell,i,j,k,l,iface,nbinxy
    logical(lg)                   :: icheck
    real(rp)                      :: xmin,ymin,zmin,xmax,ymax,zmax,rsize,rcellsize
    real(rp)                      :: xdist,ydist,zdist,xtmin,ytmin,ztmin
    real(rp)                      :: xtmax,ytmax,ztmax,xc,yc,zc,coefmax,c05,c11,c01,c14,c00
    !
    !     This sub marked the cells with respect to the size of the (pseudo) intersected triangles.
    !     The intersection is box-box so it includes more intersections than it should be.
    !     Only the newly created cells (marked with -2) are checked so there can not be any
    !     conforming problems as only the cells of the last level (created by divideCell) 
    !     are checked. 
    !     On output, the cells of the last level are marked with: 
    !          -  1 for intersected and refinement needed     
    !          - -1 for not intersected 
    !          - -3 for intersected and no refinement needed
    ! 
    !
    nbinxy=nbinx*nbinx
    coefmax=2.0d+00
    c00=0.0d+00
    c05=0.5d+00
    c14=1.0d+00/4.0d+00
    c01=0.1d+00
    c11=1.1d+00
    !
    !     Loop on the lcells
    !
    do icell=1,ncell
       !
       !   Check only the new cells
       !
       if(lcell(icell)%marked==-2)then 

          xmin=lcell(icell)%coor(1,1)
          ymin=lcell(icell)%coor(2,1)
          zmin=lcell(icell)%coor(3,1)  
          xmax=lcell(icell)%coor(1,2)
          ymax=lcell(icell)%coor(2,2)
          zmax=lcell(icell)%coor(3,2)  


          !xdist=xmax-xmin 
          !ydist=ymax-ymin 
          !zdist=zmax-zmin 

          !rcellsize=xdist+ydist+zdist
          rcellsize=lcell(icell)%rsize
          !
          !     Tolerance
          ! 
          rcellsize=rcellsize*c05 
          !rcellsize=rcellsize*c14
          ! 
          !     Expand the lcell a bit
          !
          xc=(xmax+xmin)*c05
          yc=(ymax+ymin)*c05
          zc=(zmax+zmin)*c05

          xmin=-xc*c01+xmin*c11
          ymin=-yc*c01+ymin*c11
          zmin=-zc*c01+zmin*c11

          xmax=-xc*c01+xmax*c11
          ymax=-yc*c01+ymax*c11
          zmax=-zc*c01+zmax*c11
          ! 
          !     Find the bin containing the lcell
          !
          imin=floor((xmin-bboxbin(1,1))/dxbin)+1_ip
          jmin=floor((ymin-bboxbin(2,1))/dybin)+1_ip
          kmin=floor((zmin-bboxbin(3,1))/dzbin)+1_ip

          imax=floor((xmax-bboxbin(1,1))/dxbin)+1_ip
          jmax=floor((ymax-bboxbin(2,1))/dybin)+1_ip
          kmax=floor((zmax-bboxbin(3,1))/dzbin)+1_ip
          !
          !     Limit the extend
          !
          if(imax<1)then 
             lcell(icell)%marked=-1_ip
             cycle
          else if(imin<1)then
             imin=1
          endif

          if(jmax<1)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(jmin<1)then
             jmin=1
          endif

          if(kmax<1)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(kmin<1)then
             kmin=1
          endif


          if(imin>nbinx)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(imax>nbinx)then
             imax=nbinx
          endif
          if(jmin>nbinx)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(jmax>nbinx)then
             jmax=nbinx
          endif
          if(kmin>nbinx)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(kmax>nbinx)then
             kmax=nbinx  
          endif

          icheck=.false.
          !
          !     loop on the voxels
          !
          do k=kmin,kmax
             kbinx=(k-1)*nbinxy
             do j=jmin,jmax
                jbinx=(j-1)*nbinx+kbinx
                do i=imin,imax
                   ibin=i+jbinx

                   isto0=btoel2(ibin)
                   isto1=btoel2(ibin+1)-1             

                   do l=isto0,isto1
                      !
                      !     Real test
                      !
                      iface=btoel1(l)    
                      xtmin=tribox(1,iface)
                      ytmin=tribox(2,iface)
                      ztmin=tribox(3,iface)
                      xtmax=tribox(4,iface)
                      ytmax=tribox(5,iface)
                      ztmax=tribox(6,iface)

                      if(xtmin>xmax)cycle
                      if(ytmin>ymax)cycle
                      if(ztmin>zmax)cycle
                      if(xtmax<xmin)cycle
                      if(ytmax<ymin)cycle
                      if(ztmax<zmin)cycle

                      rsize=trisize(iface)
                      !
                      !     Is the lcell too big ?
                      !
                      if(rsize<rcellsize)then

                         lcell(icell)%marked=1_ip
                         marked=marked+1_ip
                         goto 10  

                      else 
                         !
                         !     No refinement needed, mark the lcell 
                         !     for further passes as intersected but ok with size
                         !
                         icheck=.true.   
                         !
                         !     Transfer the size
                         ! 
                         lcell(icell)%rsize=min(rsize,lcell(icell)%rsize)
                         if(lcell(icell)%rsize==c00)then
                            write(*,*)'Error size c00'
                            stop 
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo

          if(icheck.eqv. .true.)then 
             lcell(icell)%marked=-3_ip
          else 
             lcell(icell)%marked=-1_ip
          endif

10        continue

       endif
    enddo
  end subroutine markCell
  subroutine realCut(btoel1, btoel2,nbin,nbinx,ncell,lcell,dxbin,dybin,dzbin,bboxbin,&
       trisize,marked,tribox,ndim,nface,lface,nnofa,coor,npoin)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: nbinx,ndim,ncell,nbin,nface,nnofa,npoin
    integer(ip),intent(in)        :: btoel1(*),btoel2(*),lface(nnofa,nface)
    real(rp),intent(in)           :: dxbin,dybin,dzbin,trisize(nface),coor(ndim,npoin)
    real(rp),intent(in)           :: tribox(6,nface),bboxbin(ndim,2)
    integer(ip),intent(inout)     :: marked
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: imin,jmin,kmin,imax,jmax,kmax,kbinx,jbinx,ibin,isto0,isto1
    integer(ip)                   :: icell,i,j,k,l,iface,nbinxy,ip1,ip2,ip3,icont,ifac,jp1,jp2,jp3
    logical(lg)                   :: icheck
    real(rp)                      :: xmin,ymin,zmin,xmax,ymax,zmax,rsize,rcellsize
    real(rp)                      :: xdist,ydist,zdist,xtmin,ytmin,ztmin,x,y,z,coorb(ndim,8)
    real(rp)                      :: xtmax,ytmax,ztmax,xc,yc,zc,coefmax,c05,c11,c01,c14,c00
    integer(ip)                   :: lfac(3,12)=RESHAPE(&
         (/5,1,2,5,2,6,6,2,3,6,3,7,7,3,4,7,4,8,8,4,1,8,1,5,5,6,7,5,7,8,2,1,4,2,4,3/),(/3,12/))
    !
    !     This sub marked the cells with respect to the size of the (pseudo) intersected triangles.
    !     The intersection is box-box so it includes more intersections than it should be.
    !     Only the newly created cells (marked with -2) are checked so there can not be any
    !     conforming problems as only the cells of the last level (created by divideCell) 
    !     are checked. 
    !     On output, the cells of the last level are marked with: 
    !          -  1 for intersected and refinement needed     
    !          - -1 for not intersected 
    !          - -3 for intersected and no refinement needed
    ! 
    !
    nbinxy=nbinx*nbinx
    coefmax=2.0d+00
    c00=0.0d+00
    c05=0.5d+00
    c14=1.0d+00/4.0d+00
    c01=0.1d+00
    c11=1.1d+00
    !
    !     Loop on the lcells
    !
    do icell=1,ncell
       !
       !   Check only the marked faces
       !
       if(lcell(icell)%marked==-3)then 

          xmin=lcell(icell)%coor(1,1)
          ymin=lcell(icell)%coor(2,1)
          zmin=lcell(icell)%coor(3,1)  
          xmax=lcell(icell)%coor(1,2)
          ymax=lcell(icell)%coor(2,2)
          zmax=lcell(icell)%coor(3,2)  


          !xdist=xmax-xmin 
          !ydist=ymax-ymin 
          !zdist=zmax-zmin 

          !rcellsize=xdist+ydist+zdist
          rcellsize=lcell(icell)%rsize
          !
          !     Tolerance
          ! 
          rcellsize=rcellsize*c05 
          !rcellsize=rcellsize*c14
          ! 
          !     Expand the lcell a bit
          !
          xc=(xmax+xmin)*c05
          yc=(ymax+ymin)*c05
          zc=(zmax+zmin)*c05

          xmin=-xc*c01+xmin*c11
          ymin=-yc*c01+ymin*c11
          zmin=-zc*c01+zmin*c11

          xmax=-xc*c01+xmax*c11
          ymax=-yc*c01+ymax*c11
          zmax=-zc*c01+zmax*c11
          ! 
          !     Find the bin containing the lcell
          !
          imin=floor((xmin-bboxbin(1,1))/dxbin)+1_ip
          jmin=floor((ymin-bboxbin(2,1))/dybin)+1_ip
          kmin=floor((zmin-bboxbin(3,1))/dzbin)+1_ip

          imax=floor((xmax-bboxbin(1,1))/dxbin)+1_ip
          jmax=floor((ymax-bboxbin(2,1))/dybin)+1_ip
          kmax=floor((zmax-bboxbin(3,1))/dzbin)+1_ip
          !
          !     Limit the extend
          !
          if(imax<1)then 
             lcell(icell)%marked=-1_ip
             cycle
          else if(imin<1)then
             imin=1
          endif

          if(jmax<1)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(jmin<1)then
             jmin=1
          endif

          if(kmax<1)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(kmin<1)then
             kmin=1
          endif


          if(imin>nbinx)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(imax>nbinx)then
             imax=nbinx
          endif
          if(jmin>nbinx)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(jmax>nbinx)then
             jmax=nbinx
          endif
          if(kmin>nbinx)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(kmax>nbinx)then
             kmax=nbinx  
          endif

          icheck=.false.
          !
          !     loop on the voxels
          !
          do k=kmin,kmax
             kbinx=(k-1)*nbinxy
             do j=jmin,jmax
                jbinx=(j-1)*nbinx+kbinx
                do i=imin,imax
                   ibin=i+jbinx

                   isto0=btoel2(ibin)
                   isto1=btoel2(ibin+1)-1             

                   do l=isto0,isto1
                      !
                      !     First test with bounding box
                      !
                      iface=btoel1(l)    
                      xtmin=tribox(1,iface)
                      ytmin=tribox(2,iface)
                      ztmin=tribox(3,iface)
                      xtmax=tribox(4,iface)
                      ytmax=tribox(5,iface)
                      ztmax=tribox(6,iface)

                      if(xtmin>xmax)cycle
                      if(ytmin>ymax)cycle
                      if(ztmin>zmax)cycle
                      if(xtmax<xmin)cycle
                      if(ytmax<ymin)cycle
                      if(ztmax<zmin)cycle
                      !
                      !     Second test if one point of the triangle is inside the box
                      !
                      ip1=lface(1,iface) 
                      ip2=lface(2,iface) 
                      ip3=lface(3,iface) 

                      x=coor(1,ip1) 
                      y=coor(2,ip1) 
                      z=coor(3,ip1) 
                      if(x<xmax .and. x>xmin .and. y<ymax .and. y>ymin .and. z<zmax .and. z>zmin)then
                         icont=1_ip   
                      else
                         icont=-1_ip   
                      endif
                      x=coor(1,ip2) 
                      y=coor(2,ip2) 
                      z=coor(3,ip2) 
                      if(x<xmax .and. x>xmin .and. y<ymax .and. y>ymin .and. z<zmax .and. z>zmin)then
                         icont=icont+1_ip   
                      else
                         icont=icont-1_ip   
                      endif
                      x=coor(1,ip3) 
                      y=coor(2,ip3) 
                      z=coor(3,ip3) 
                      if(x<xmax .and. x>xmin .and. y<ymax .and. y>ymin .and. z<zmax .and. z>zmin)then
                         icont=icont+1_ip   
                      else
                         icont=icont-1_ip   
                      endif
                      !
                      !   If the triangle is not completely inside or outside, it intersects the box
                      !
                      if(abs(icont)/= 3)goto 10 
                      !
                      !     Copy the box coordinates
                      !
                      xtmin=lcell(icell)%coor(1,1)
                      ytmin=lcell(icell)%coor(2,1)
                      ztmin=lcell(icell)%coor(3,1)
                      xtmax=lcell(icell)%coor(1,2)
                      ytmax=lcell(icell)%coor(2,2)
                      ztmax=lcell(icell)%coor(3,2)

                      coorb(1,1)=xmin
                      coorb(2,1)=ymin
                      coorb(3,1)=zmin
                      coorb(1,2)=xmax
                      coorb(2,2)=ymin
                      coorb(3,2)=zmin
                      coorb(1,3)=xmax
                      coorb(2,3)=ymax
                      coorb(3,3)=zmin
                      coorb(1,4)=xmin
                      coorb(2,4)=ymax
                      coorb(3,4)=zmin
                      coorb(1,5)=xmin
                      coorb(2,5)=ymin
                      coorb(3,5)=zmax
                      coorb(1,6)=xmax
                      coorb(2,6)=ymin
                      coorb(3,6)=zmax
                      coorb(1,7)=xmax
                      coorb(2,7)=ymax
                      coorb(3,7)=zmax
                      coorb(1,8)=xmin
                      coorb(2,8)=ymax
                      coorb(3,8)=zmax
                      !
                      !     Last test, full intersection
                      !
                      do ifac=1,12
                         jp1=lfac(1,ifac) 
                         jp2=lfac(2,ifac) 
                         jp3=lfac(3,ifac) 
                         icheck=.false.
                         call intersect(ip1,ip2,ip3,coor,jp1,jp2,jp3,coorb,npoin,ndim,icheck) 
                         if(icheck.eqv. .true.)goto 10
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
          !     As we reached this part, mark the cell as not intersected
          ! 
          lcell(icell)%marked=-1_ip

10        continue

       endif
    enddo

  end subroutine realCut

  subroutine intersect(ip1,ip2,ip3,coor,jp1,jp2,jp3,coorb,npoin,ndim,icheck) 
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ndim,npoin
    integer(ip),intent(in)        :: ip1,ip2,ip3,jp1,jp2,jp3
    logical(lg),intent(inout)     :: icheck
    real(rp),intent(in)           :: coor(ndim,npoin),coorb(ndim,12)
    integer(ip)                   :: ltab2(3),ip11,ip12,jp11,jp12,iedg
    integer(ip)                   :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    real(rp)                      :: rtx1,rty1,rtz1,rnx1,rny1,rnz1,rtl,c10 
    real(rp)                      :: rnx2,rny2,rnz2,rnl2,rnl22,rnl1,rnl11 
    real(rp)                      :: rtx2,rty2,rtz2,epsil,rv1,rv2,rv3,rnl 

    c10=1.0d+00
    epsil=1.0d-12

    ltab2(1)=ip1
    ltab2(2)=ip2
    ltab2(3)=ip3
    !
    !     Get normal to triangle (jp1,jp2,jp3)
    !
    rtx1=coorb(1,jp2)-coorb(1,jp1)
    rty1=coorb(2,jp2)-coorb(2,jp1)
    rtz1=coorb(3,jp2)-coorb(3,jp1)

    rtx2=coorb(1,jp3)-coorb(1,jp1)
    rty2=coorb(2,jp3)-coorb(2,jp1)
    rtz2=coorb(3,jp3)-coorb(3,jp1)

    rnx1= rty1*rtz2-rtz1*rty2
    rny1=-rtx1*rtz2+rtz1*rtx2 
    rnz1= rtx1*rty2-rty1*rtx2
    rnl1=sqrt(rnx1*rnx1+rny1*rny1+rnz1*rnz1)
    rnl11=rnl1
    rnl1=c10/rnl1
    rnx1=rnx1*rnl1
    rny1=rny1*rnl1
    rnz1=rnz1*rnl1
    !
    !     Get normal to triangle (ip1,ip2,ip3)
    !
    rtx1=coor(1,ip2)-coor(1,ip1)
    rty1=coor(2,ip2)-coor(2,ip1)
    rtz1=coor(3,ip2)-coor(3,ip1)

    rtx2=coor(1,ip3)-coor(1,ip1)
    rty2=coor(2,ip3)-coor(2,ip1)
    rtz2=coor(3,ip3)-coor(3,ip1)

    rnx2= rty1*rtz2-rtz1*rty2
    rny2=-rtx1*rtz2+rtz1*rtx2 
    rnz2= rtx1*rty2-rty1*rtx2
    rnl2=sqrt(rnx2*rnx2+rny2*rny2+rnz2*rnz2)
    rnl22=rnl2
    rnl2=c10/rnl2
    rnx2=rnx2*rnl2
    rny2=rny2*rnl2
    rnz2=rnz2*rnl2
    !
    !     Get meaningful epsilon 
    ! 
    rnl=min(rnl11,rnl22)
    epsil=sqrt(rnl)*epsil
    !
    !     First part, intersection triangle (jp1,jp2,jp3) against edges (ip1,ip2),(ip2,ip3),(ip3,ip1)
    !  
    do iedg=1,3

       ip11=ltab2(ltab(1,iedg))  
       ip12=ltab2(ltab(2,iedg))  

       rtx1=coor(1,ip11)-coorb(1,jp1) 
       rty1=coor(2,ip11)-coorb(2,jp1) 
       rtz1=coor(3,ip11)-coorb(3,jp1) 
       rtl=sqrt(rtx1*rtx1+rty1*rty1+rtz1*rtz1)
       rtl=c10/rtl
       rtx1=rtx1*rtl
       rty1=rty1*rtl
       rtz1=rtz1*rtl

       rtx2=coor(1,ip12)-coorb(1,jp1) 
       rty2=coor(2,ip12)-coorb(2,jp1) 
       rtz2=coor(3,ip12)-coorb(3,jp1) 
       rtl=sqrt(rtx2*rtx2+rty2*rty2+rtz2*rtz2)
       rtl=c10/rtl
       rtx2=rtx2*rtl
       rty2=rty2*rtl
       rtz2=rtz2*rtl

       rv1=rtx1*rnx1+rty1*rny1+rtz1*rnz1
       rv2=rtx2*rnx1+rty2*rny1+rtz2*rnz1
       !
       !     At least (jp1,jp2,jp3) must separe ip11 and ip12
       !
       if(rv1*rv2<epsil)then

          call gtvolc(ip11,ip12,coor,jp1,jp2,coorb,ndim,npoin,rv1) 
          call gtvolc(ip11,ip12,coor,jp2,jp3,coorb,ndim,npoin,rv2) 
          call gtvolc(ip11,ip12,coor,jp3,jp1,coorb,ndim,npoin,rv3) 

          if(rv1>epsil .and. rv2>epsil .and. rv3>epsil)then
             icheck=.true.
             return
          endif
          if(rv1<-epsil .and. rv2<-epsil .and. rv3<-epsil)then
             icheck=.true.
             return
          endif

       endif

    enddo

    ltab2(1)=jp1
    ltab2(2)=jp2
    ltab2(3)=jp3
    !
    !     Second part, intersection triangle (ip1,ip2,ip3) against edges (jp1,jp2),(jp2,jp3),(jp3,jp1)
    !  
    do iedg=1,3

       jp11=ltab2(ltab(1,iedg))  
       jp12=ltab2(ltab(2,iedg))  

       rtx1=coorb(1,jp11)-coor(1,ip1) 
       rty1=coorb(2,jp11)-coor(2,ip1) 
       rtz1=coorb(3,jp11)-coor(3,ip1) 
       rtl=sqrt(rtx1*rtx1+rty1*rty1+rtz1*rtz1)
       rtl=c10/rtl
       rtx1=rtx1*rtl
       rty1=rty1*rtl
       rtz1=rtz1*rtl

       rtx2=coorb(1,jp12)-coor(1,ip1) 
       rty2=coorb(2,jp12)-coor(2,ip1) 
       rtz2=coorb(3,jp12)-coor(3,ip1) 
       rtl=sqrt(rtx2*rtx2+rty2*rty2+rtz2*rtz2)
       rtl=c10/rtl
       rtx2=rtx2*rtl
       rty2=rty2*rtl
       rtz2=rtz2*rtl

       rv1=rtx1*rnx2+rty1*rny2+rtz1*rnz2
       rv2=rtx2*rnx2+rty2*rny2+rtz2*rnz2
       !
       !     At least (ip1,ip2,ip3) must separe jp11 and jp12
       !
       if(rv1*rv2<epsil)then

          call gtvolc(ip1,ip2,coor,jp11,jp12,coorb,ndim,npoin,rv1) 
          call gtvolc(ip2,ip3,coor,jp11,jp12,coorb,ndim,npoin,rv2) 
          call gtvolc(ip3,ip1,coor,jp11,jp12,coorb,ndim,npoin,rv3) 

          if(rv1>epsil .and. rv2>epsil .and. rv3>epsil)then
             icheck=.true.
             return
          endif
          if(rv1<-epsil .and. rv2<-epsil .and. rv3<-epsil)then
             icheck=.true.
             return
          endif

       endif

    enddo

  end subroutine intersect

  subroutine gtvolc(ip1,ip2,coor,jp1,jp2,coorb,ndim,npoin,rvol)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)     :: ip1,ip2,jp1,jp2,ndim,npoin
    real(rp),intent(in)        :: coor(ndim,npoin),coorb(ndim,8)
    real(rp),intent(inout)     :: rvol
    real(rp)                   :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    real(rp)                   :: x21,y21,z21,x31,y31,z31,x41,y41,z41,c16

    c16=1.0d+00/6.0d+00

    x1=coor(1,ip1)
    y1=coor(2,ip1)
    z1=coor(3,ip1)
    x2=coor(1,ip2)
    y2=coor(2,ip2)
    z2=coor(3,ip2)
    x3=coorb(1,jp1)
    y3=coorb(2,jp1)
    z3=coorb(3,jp1)
    x4=coorb(1,jp2)
    y4=coorb(2,jp2)
    z4=coorb(3,jp2)

    x21=x2-x1
    y21=y2-y1
    z21=z2-z1

    x31=x3-x1
    y31=y3-y1
    z31=z3-z1

    x41=x4-x1
    y41=y4-y1
    z41=z4-z1

    rvol=c16*(x21*(y31*z41-z31*y41) + x31*(z21*y41-y21*z41) + x41*(y21*z31-z21*y31))

  end subroutine gtvolc


  subroutine markCellGeo(&
       btoel1,btoel2,nbin,nbinx,ncell,lcell,dxbin,dybin,dzbin,bboxbin,marked,&
       tribox,ndim,nface,eltoel,rnofa,nnofa,maxlevel)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(4)                :: istat
    integer(ip),intent(in)        :: nbinx,ndim,ncell,nbin,nface,nnofa,maxlevel
    integer(ip),intent(in)        :: btoel1(*),btoel2(*),eltoel(nnofa,nface)
    real(rp),intent(in)           :: dxbin,dybin,dzbin,tribox(6,nface),bboxbin(ndim,2)
    real(rp),intent(in)           :: rnofa(ndim,nface)
    integer(ip),intent(inout)     :: marked
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: imin,jmin,kmin,imax,jmax,kmax,kbinx,jbinx,ibin,isto0,isto1
    integer(ip)                   :: icell,i,j,k,l,iface,icheck,nbinxy,nstack,istack,ineigh
    real(rp)                      :: xmin,ymin,zmin,xmax,ymax,zmax,rsize,xdist,ydist,zdist,xtmin,ytmin,ztmin
    real(rp)                      :: xtmax,ytmax,ztmax,xc,yc,zc,coefmax,c05,c11,c01,cscal,tolscal
    integer(ip),pointer           :: lmark(:),lstack(:)


    nbinxy=nbinx*nbinx
    coefmax=2.0d+00
    c05=0.5d+00
    c01=0.1d+00
    c11=1.1d+00
    tolscal=0.3d+00  

    allocate(lmark(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','markCellGeo',lmark)
    allocate(lstack(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','markCellGeo',lstack)
    !
    !     Loop on the lcells
    !
    do icell=1,ncell
       !
       !     Check only the new faces
       !
       if(lcell(icell)%marked==-2)then 
          !
          !     Did we reach the maximum
          ! 
          if(lcell(icell)%level==maxlevel)then
             lcell(icell)%marked=-1
             cycle
          endif


          xmin=lcell(icell)%coor(1,1)
          ymin=lcell(icell)%coor(2,1)
          zmin=lcell(icell)%coor(3,1)  
          xmax=lcell(icell)%coor(1,2)
          ymax=lcell(icell)%coor(2,2)
          zmax=lcell(icell)%coor(3,2)  


          xdist=xmax-xmin 
          ydist=ymax-ymin 
          zdist=zmax-zmin 
          !
          !     Expand the lcell a bit
          !
          xc=(xmax+xmin)*c05
          yc=(ymax+ymin)*c05
          zc=(zmax+zmin)*c05

          xmin=-xc*c01+xmin*c11
          ymin=-yc*c01+ymin*c11
          zmin=-zc*c01+zmin*c11

          xmax=-xc*c01+xmax*c11
          ymax=-yc*c01+ymax*c11
          zmax=-zc*c01+zmax*c11
          !
          !     Find the bin containing the lcell
          !
          imin=floor((xmin-bboxbin(1,1))/dxbin)+1_ip
          jmin=floor((ymin-bboxbin(2,1))/dybin)+1_ip
          kmin=floor((zmin-bboxbin(3,1))/dzbin)+1_ip

          imax=floor((xmax-bboxbin(1,1))/dxbin)+1_ip
          jmax=floor((ymax-bboxbin(2,1))/dybin)+1_ip
          kmax=floor((zmax-bboxbin(3,1))/dzbin)+1_ip
          !
          !     Limit the extend
          !
          if(imax<1)then 
             lcell(icell)%marked=-1_ip
             cycle
          else if(imin<1)then
             imin=1
          endif

          if(jmax<1)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(jmin<1)then
             jmin=1
          endif

          if(kmax<1)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(kmin<1)then
             kmin=1
          endif


          if(imin>nbinx)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(imax>nbinx)then
             imax=nbinx
          endif
          if(jmin>nbinx)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(jmax>nbinx)then
             jmax=nbinx
          endif
          if(kmin>nbinx)then
             lcell(icell)%marked=-1_ip
             cycle
          else if(kmax>nbinx)then
             kmax=nbinx  
          endif

          nstack=0_ip
          !
          !     loop on the voxels
          ! 
          do k=kmin,kmax
             kbinx=(k-1)*nbinxy
             do j=jmin,jmax
                jbinx=(j-1)*nbinx+kbinx
                do i=imin,imax
                   ibin=i+jbinx

                   isto0=btoel2(ibin)
                   isto1=btoel2(ibin+1)-1             

                   do l=isto0,isto1
                      !
                      !     Mark the faces belonging to this voxel
                      !
                      iface=btoel1(l)    
                      xtmin=tribox(1,iface)
                      ytmin=tribox(2,iface)
                      ztmin=tribox(3,iface)
                      xtmax=tribox(4,iface)
                      ytmax=tribox(5,iface)
                      ztmax=tribox(6,iface)

                      if(xtmin>xmax)cycle
                      if(ytmin>ymax)cycle
                      if(ztmin>zmax)cycle
                      if(xtmax<xmin)cycle
                      if(ytmax<ymin)cycle
                      if(ztmax<zmin)cycle
                      ! 
                      !     Add to the stack
                      !
                      if(lmark(iface)/=icell)then 
                         lmark(iface)=icell
                         nstack=nstack+1_ip
                         lstack(nstack)=iface
                      endif
                   enddo
                enddo
             enddo
          enddo

          if(nstack==0)then
             lcell(icell)%marked=-1_ip
             cycle
          endif
          !
          !     Compare scalar products
          !
          do istack=1,nstack 
             iface=lstack(istack)
             !
             !     Mark the face to avoid to compute twice the scalar product
             !  
             lmark(iface)=0_ip
             ineigh=eltoel(1,iface)
             if(ineigh/=0)then 
                if(lmark(ineigh)==icell)then
                   cscal=rnofa(1,iface)*rnofa(1,ineigh)+rnofa(2,iface)*rnofa(2,ineigh)+rnofa(3,iface)*rnofa(3,ineigh)
                   if(cscal<tolscal)then
                      lcell(icell)%marked=1_ip
                      marked=marked+1_ip
                      goto 10  
                   endif
                endif
             endif
             ineigh=eltoel(2,iface)
             if(ineigh/=0)then 
                if(lmark(ineigh)==icell)then
                   cscal=rnofa(1,iface)*rnofa(1,ineigh)+rnofa(2,iface)*rnofa(2,ineigh)+rnofa(3,iface)*rnofa(3,ineigh)
                   if(cscal<tolscal)then
                      lcell(icell)%marked=1_ip
                      marked=marked+1_ip
                      goto 10  
                   endif
                endif
             endif
             ineigh=eltoel(3,iface)
             if(ineigh/=0)then 
                if(lmark(ineigh)==icell)then
                   cscal=rnofa(1,iface)*rnofa(1,ineigh)+rnofa(2,iface)*rnofa(2,ineigh)+rnofa(3,iface)*rnofa(3,ineigh)
                   if(cscal<tolscal)then
                      lcell(icell)%marked=1_ip
                      marked=marked+1_ip
                      goto 10  
                   endif
                endif
             endif
          enddo

          lcell(icell)%marked=-1_ip

10        continue 

          if(lcell(icell)%marked==1)then
             !write(*,*)'cscal=',cscal 
             !write(*,*)lface(1,iface),lface(2,iface),lface(3,iface)
             !write(*,*)lface(1,ineigh),lface(2,ineigh),lface(3,ineigh)
             write(*,*)lcell(icell)%level
          endif

       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','markCellGeo',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','markCellGeo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','markCellGeo',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','markCellGeo',0_ip)

  end subroutine markCellGeo

  subroutine markCellSou(ncell,lcell,marked,ndim,nsour,rsour,rsgeo,iter,rtol,rsuni)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(4)                :: istat
    integer(ip),intent(in)        :: ndim,ncell,nsour,iter
    real(rp),intent(in)           :: rsgeo(ndim,ndim,nsour),rtol,rsuni
    real(rp),pointer              :: rsour(:,:)
    integer(ip),intent(inout)     :: marked
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: icell,isour,l,neigh,icart,nstack,istack
    real(rp)                      :: xc,yc,zc,c10,csca,cscal,rpoin(ndim),rcellsize
    real(rp)                      :: rvec1x,rvec1y,rvec1z,rvec2x,rvec2y,rvec2z,cbig
    real(rp)                      :: rnx,rny,rnz,rface(ndim),rnl,p1(ndim),p2(ndim)
    real(rp)                      :: dtot,d1,d2,d3,rdist,epsil,epsil1,epsil0,c13,c05,rsize
    real(rp)                      :: rdirx,rdiry,rdirz,xpro,ypro,zpro,rl,rdist1,rdist2,rdist3,rfact
    integer(ip),pointer           :: lmark(:),lstack(:)

    c10=1.0d+00
    c05=0.5d+00
    c13=1.0d+00/3.0d+00
    epsil=1.0d-06
    epsil1=1.0d+00+1.0d-04
    epsil0=1.0d-04
    cbig=1.0d+12*rsuni 
    !
    !     Ratio between cell size and source size
    ! 
    rfact=0.25d+00
    rfact=1.0d+00

    !allocate(lmark(ncell),stat=istat)
    !call memchk(zero,istat,memor_msh,'LMARK','markCellSou',lmark)
    !allocate(lstack(ncell),stat=istat)
    !call memchk(zero,istat,memor_msh,'LSTACK','markCellSou',lstack)
    !
    !     Initialize the cell size to a huge value
    !
    !do icell=1,ncell
    !   lcell(icell)%rsize=cbig
    !enddo
    !
    !     Loop on the sources
    !
    do isour=1,nsour
       !
       !     Compute the center of the triangle
       !
       rpoin(1)=(rsgeo(1,1,isour)+rsgeo(1,2,isour)+rsgeo(1,3,isour))*c13
       rpoin(2)=(rsgeo(2,1,isour)+rsgeo(2,2,isour)+rsgeo(2,3,isour))*c13
       rpoin(3)=(rsgeo(3,1,isour)+rsgeo(3,2,isour)+rsgeo(3,3,isour))*c13
       !
       !     Compute source side vectors
       !
       rvec1x=rsgeo(1,2,isour)-rsgeo(1,1,isour)
       rvec1y=rsgeo(2,2,isour)-rsgeo(2,1,isour)
       rvec1z=rsgeo(3,2,isour)-rsgeo(3,1,isour)

       rvec2x=rsgeo(1,3,isour)-rsgeo(1,1,isour)
       rvec2y=rsgeo(2,3,isour)-rsgeo(2,1,isour)
       rvec2z=rsgeo(3,3,isour)-rsgeo(3,1,isour)
       !
       !     Compute normal to source
       !
       rnx= rvec1y*rvec2z-rvec1z*rvec2y 
       rny=-rvec1x*rvec2z+rvec1z*rvec2x 
       rnz= rvec1x*rvec2y-rvec1y*rvec2x 

       rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
       !
       !     Remember face area (*2)
       !
       dtot=rnl
       rnl=c10/rnl

       rface(1)=rnl*rnx
       rface(2)=rnl*rny
       rface(3)=rnl*rnz
       !
       !     Find the cell containing this point
       !
       !icart=1_ip
       !call gtelem(1,rpoin,1,ndim,lcell,ncell,icart,rtol)
       !if(icart==0)then
       !   write(*,*)'Error in markCellSou, point out of domain'
       !   stop 
       !endif
       !
       !     Initialize the stack with icart
       !
       !lstack(1)=icart
       !lmark(icart)=isour
       !nstack=1_ip
       !istack=0_ip
       !
       !     Propagate the source size
       !
       do icell=1,ncell
          !
          !   Check only the new cells
          !
          if(lcell(icell)%marked==-2)then 
             !do  

             !   if(istack==nstack)exit
             !   istack=istack+1 
             !
             !     Get the cell number
             !  
             !   icell=lstack(istack)
             !
             !     Get cell center
             !
             xc=(lcell(icell)%coor(1,1)+lcell(icell)%coor(1,2))*c05
             yc=(lcell(icell)%coor(2,1)+lcell(icell)%coor(2,2))*c05
             zc=(lcell(icell)%coor(3,1)+lcell(icell)%coor(3,2))*c05

             rcellsize=(lcell(icell)%coor(1,2)-xc)+(lcell(icell)%coor(2,2)-yc)+(lcell(icell)%coor(3,2)-zc)*rfact
             !
             !     Project the point on the plane of the source
             !
             rdirx=xc-rsgeo(1,1,isour)
             rdiry=yc-rsgeo(2,1,isour)
             rdirz=zc-rsgeo(3,1,isour)

             csca=rface(1)*rdirx+rface(2)*rdiry+rface(3)*rdirz

             xpro=xc-csca*rface(1)
             ypro=yc-csca*rface(2)
             zpro=zc-csca*rface(3)
             !
             !     Side (ip2,ip3)
             !
             p1(1)=rsgeo(1,2,isour)-xpro
             p1(2)=rsgeo(2,2,isour)-ypro
             p1(3)=rsgeo(3,2,isour)-zpro
             p2(1)=rsgeo(1,3,isour)-xpro
             p2(2)=rsgeo(2,3,isour)-ypro
             p2(3)=rsgeo(3,3,isour)-zpro
             call orient3D(p1,p2,rface,d1,ndim)
             d1=d1/dtot
             !
             !     Side (ip3,ip1)
             !
             p1(1)=rsgeo(1,3,isour)-xpro
             p1(2)=rsgeo(2,3,isour)-ypro
             p1(3)=rsgeo(3,3,isour)-zpro
             p2(1)=rsgeo(1,1,isour)-xpro
             p2(2)=rsgeo(2,1,isour)-ypro
             p2(3)=rsgeo(3,1,isour)-zpro
             call orient3D(p1,p2,rface,d2,ndim)
             d2=d2/dtot
             !
             !     Side (ip1,ip2)
             !
             p1(1)=rsgeo(1,1,isour)-xpro
             p1(2)=rsgeo(2,1,isour)-ypro
             p1(3)=rsgeo(3,1,isour)-zpro
             p2(1)=rsgeo(1,2,isour)-xpro
             p2(2)=rsgeo(2,2,isour)-ypro
             p2(3)=rsgeo(3,2,isour)-zpro
             call orient3D(p1,p2,rface,d3,ndim)
             d3=d3/dtot
             !
             !     Is the projected point inside the source triangle?
             !
             if(d1>-epsil .and. d2>-epsil .and. d3>-epsil)then  
                rdist=abs(csca)
             else
                !
                !     Evaluate distance to lines
                !

                !
                !     Side (ip2,ip3)
                !
                p1(1)=rsgeo(1,3,isour)-rsgeo(1,2,isour)
                p1(2)=rsgeo(2,3,isour)-rsgeo(2,2,isour)
                p1(3)=rsgeo(3,3,isour)-rsgeo(3,2,isour)
                rl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
                rl=c10/rl
                p1(1)=rl*p1(1) 
                p1(2)=rl*p1(2) 
                p1(3)=rl*p1(3) 

                p2(1)=xc-rsgeo(1,2,isour)
                p2(2)=yc-rsgeo(2,2,isour)
                p2(3)=zc-rsgeo(3,2,isour)

                csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
                cscal=csca*rl 

                if(cscal<-epsil0)then

                   rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

                else if(cscal>epsil1)then

                   p2(1)=xc-rsgeo(1,3,isour)
                   p2(2)=yc-rsgeo(2,3,isour)
                   p2(3)=zc-rsgeo(3,3,isour)
                   rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

                else

                   p2(1)=rsgeo(1,2,isour)+csca*p1(1)-xc
                   p2(2)=rsgeo(2,2,isour)+csca*p1(2)-yc
                   p2(3)=rsgeo(3,2,isour)+csca*p1(3)-zc
                   rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 

                endif
                !
                !     Side (ip3,ip1)
                !
                p1(1)=rsgeo(1,1,isour)-rsgeo(1,3,isour)
                p1(2)=rsgeo(2,1,isour)-rsgeo(2,3,isour)
                p1(3)=rsgeo(3,1,isour)-rsgeo(3,3,isour)
                rl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
                rl=c10/rl
                p1(1)=rl*p1(1) 
                p1(2)=rl*p1(2) 
                p1(3)=rl*p1(3) 

                p2(1)=xc-rsgeo(1,3,isour)
                p2(2)=yc-rsgeo(2,3,isour)
                p2(3)=zc-rsgeo(3,3,isour)

                csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
                cscal=csca*rl 

                if(cscal<-epsil0)then

                   rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

                else if(cscal>epsil1)then

                   p2(1)=xc-rsgeo(1,1,isour)
                   p2(2)=yc-rsgeo(2,1,isour)
                   p2(3)=zc-rsgeo(3,1,isour)
                   rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

                else

                   p2(1)=rsgeo(1,3,isour)+csca*p1(1)-xc
                   p2(2)=rsgeo(2,3,isour)+csca*p1(2)-yc
                   p2(3)=rsgeo(3,3,isour)+csca*p1(3)-zc
                   rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 

                endif
                !
                !     Side (ip1,ip2)
                !
                p1(1)=rsgeo(1,2,isour)-rsgeo(1,1,isour)
                p1(2)=rsgeo(2,2,isour)-rsgeo(2,1,isour)
                p1(3)=rsgeo(3,2,isour)-rsgeo(3,1,isour)
                rl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
                rl=c10/rl
                p1(1)=rl*p1(1) 
                p1(2)=rl*p1(2) 
                p1(3)=rl*p1(3) 

                p2(1)=xc-rsgeo(1,1,isour)
                p2(2)=yc-rsgeo(2,1,isour)
                p2(3)=zc-rsgeo(3,1,isour)

                csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
                cscal=csca*rl 

                if(cscal<-epsil0)then

                   rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

                else if(cscal>epsil1)then

                   p2(1)=xc-rsgeo(1,2,isour)
                   p2(2)=yc-rsgeo(2,2,isour)
                   p2(3)=zc-rsgeo(3,2,isour)
                   rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

                else

                   p2(1)=rsgeo(1,1,isour)+csca*p1(1)-xc
                   p2(2)=rsgeo(2,1,isour)+csca*p1(2)-yc
                   p2(3)=rsgeo(3,1,isour)+csca*p1(3)-zc
                   rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 

                endif
                !
                !     Take the min over the side distances
                ! 
                rdist=min(rdist1,rdist2,rdist3)

             endif
             !
             !     The distance has been computed, now compute the desired size
             !
             !
             !     Are we inside the sphere of constant size
             !
             if(rdist<rsour(2,isour))then
                rsize=rsour(1,isour)
             else
                rsize=rsour(3,isour)*(rdist-rsour(2,isour))+rsour(1,isour)
             endif
             !
             !     Now compare to the cell size
             !
             !if(rsize<lcell(icell)%rsize)then 
             !   lcell(icell)%rsize=rsize
             !endif
             !
             !     As a refinement criterion, compare the cell size
             !     with its dimensions
             !    
             if(rsize<rcellsize)then
                !
                !     Was the cell already marked before?
                !   
                if(lcell(icell)%marked/=1)then
                   !
                   !     Is the cell at the last iteration?
                   ! 
                   !if(lcell(icell)%level==iter)then
                   lcell(icell)%marked=1_ip
                   marked=marked+1_ip
                   !endif
                endif
             endif
             !
             !    Add cells to the stack
             !    For the moment, add all the cells for each source...
             !
             !do l=1,6
             !   neigh=lcell(icell)%neigh(l)
             !   if(neigh==0)cycle
             !   if(lmark(neigh)/=isour)then
             !      nstack=nstack+1
             !      lstack(nstack)=neigh
             !      lmark(neigh)=isour
             !   endif
             !enddo
          endif
       enddo
    enddo

    !call memchk(2_ip,istat,memor_msh,'LMARK','markCellSou',lmark)
    !deallocate(lmark,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LMARK','markCellSou',0_ip)
    !!call memchk(2_ip,istat,memor_msh,'LSTACK','markCellSou',lstack)
    !deallocate(lstack,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LSTACK','markCellSou',0_ip)

  end subroutine markCellSou

  subroutine markCellUni(ncell,lcell,marked,rsuni)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)        :: ncell
    real(rp),intent(in)           :: rsuni
    integer(ip),intent(inout)     :: marked
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: icell
    real(rp)                      :: xc,yc,zc,rcellsize,c05

    !
    !     This sub refines the cartesian mesh given the uniform size distribution
    !     Only the cells inside the volume should be refined (Not ready for the moment...)
    !

    c05=0.5d+00

    do icell=1,ncell

       !xc=(lcell(icell)%coor(1,1)+lcell(icell)%coor(1,2))*c05
       !yc=(lcell(icell)%coor(2,1)+lcell(icell)%coor(2,2))*c05
       !zc=(lcell(icell)%coor(3,1)+lcell(icell)%coor(3,2))*c05

       !rcellsize=(lcell(icell)%coor(1,2)-xc)+(lcell(icell)%coor(2,2)-yc)+(lcell(icell)%coor(3,2)-zc)
       rcellsize=lcell(icell)%coor(1,2)-lcell(icell)%coor(1,1)
       if(rcellsize>rsuni)then
          marked=marked+1_ip
          lcell(icell)%marked=1_ip
       else
          lcell(icell)%rsize=rsuni
       endif

    enddo

  end subroutine markCellUni

  subroutine smooth(ncell,iter,ndim,lcell)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)        :: ndim,iter
    integer(ip),intent(inout)     :: ncell
    type(cell), pointer           :: lcell(:)  
    integer(ip)                   :: i,j,l,nstack,level,neigh,neighlevel,icell,ilay,ibegin,ipt,ineigh
    integer(ip)                   :: ipoint,mstack,nnew,ncont
    integer(ip),pointer           :: stack(:),remem(:),lay(:),renum(:)
    integer(4)                 :: istat
    !
    !     This sub takes care of jumps between the marked elements 
    !     It also takes care of the case when two neighboring cells 
    !     must be refined at the same time in divideFaceSmoo,
    !     which does not do divideFace, supposing everything has been previously
    !     done here.  
    !     At this point, the marked cells are marked with %marked=1
    !
    mstack=512
    !
    !     Check iter for malloc
    !
    if(iter==1)then
       return 
    endif
    !
    !    Avoid a too large jump
    !
    allocate(stack(mstack),stat=istat)
    call memchk(zero,istat,memor_msh,'STACK','smooth',stack)
    allocate(remem(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'REMEM','smooth',remem)
    allocate(lay(iter+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LAY','smooth',lay)
    !
    !     Initialize ncont
    !
    ncont=0
    !
    !     Loop on the lcells in decreasing order 
    !
    do i=ncell,1,-1
       !
       !     Has the cell been marked for refinement?
       !
       if(lcell(i)%marked>0)then
          !
          !     Initialize the stack
          !
          nstack=1
          stack(1)=i
          !
          !     Loop on the stack
          !
          do j=1,nstack

             icell=stack(j)
             !
             !     Loop on the neighbors
             !
             do l=1,6

                neigh=lcell(icell)%neigh(l)
                level=lcell(icell)%level

                if(neigh/=0)then           
                   !
                   !     Has the neighbor been marked for refinement here with 2
                   !
                   if(lcell(neigh)%marked/=2)then

                      neighlevel=lcell(neigh)%level
                      !
                      !     Is the difference level larger than 1
                      !
                      if((level-neighlevel)==1)then
                         !
                         !     Mark this cell to be refined with 2, as it considers
                         !     as much the marked cells than the non marked cells
                         !
                         remem(neigh)=1
                         lay(neighlevel)=lay(neighlevel)+1  
                         lcell(neigh)%marked=2_ip

                         call memrea(nstack+1_ip,memor_msh,'STACK','smooth',stack)
                         stack(nstack)=neigh
                         nstack=nstack+1
                         ncont=ncont+1

                      endif
                   endif
                endif
             enddo
          enddo
       endif
    enddo
    !
    !     Clean up the 2 mark  (do not clean up as these cells will be deleted)
    !
    !do icell=1,ncell
    !   if(lcell(icell)%marked==2)lcell(icell)%marked=1_ip
    !enddo  

    !printf("In smoothing ncont=%d\n",ncont) 
    !
    !     Expands lcell to add the newly marked lcells
    !
    !     Resize
    !
    nnew=(ncell)+(8*ncont) 
    call memrea(nnew,memor_msh,'LCELL','smooth',lcell)
    !
    !     Resize remem
    !
    call memrea(nnew,memor_msh,'REMEM','smooth',remem)
    !
    !     Renumber
    !
    allocate(renum(nnew),stat=istat)
    call memchk(zero,istat,memor_msh,'RENUM','smooth',renum)
    !
    !     DBG
    !  
    !do i=iter+1,1,-1
    !   write(*,*)i,lay(i)
    !enddo
    !
    !     Initialize to expand
    ! 
    ipt=iter+1
    ipoint=nnew

    do i=ncell,1,-1
       !
       !     Reset ipt at each level change
       !     
       !
       if(lcell(i)%level/=ipt)then
          !write(*,*)ipt,lcell(i)%level
          ipt=ipt-1
          !
          !     Get the new element number at this level
          !
          ilay=lay(ipt)
          ibegin=ipoint
          ipoint=ipoint-8*ilay
          !
          !     All those new elements do not have sons
          !
          do j=ibegin,ipoint+1,-1
             remem(j)=0
          enddo

       endif
       !
       !     If marked, store the son position at the old i position
       !
       if(remem(i)/=0)then
          ibegin=ibegin-8
          remem(i)=ibegin+1
       endif
       !
       !    Transfer from the old i position to the new ipoint position
       !
       remem(ipoint)=remem(i)
       !
       !     Remember the renumbering
       !
       renum(i)=ipoint
       !
       !     Decrease new pointer
       !
       ipoint=ipoint-1

    enddo

    call dbgsub(lcell,ncell,remem,renum)

    do i=ncell,1,-1
       icell=renum(i)
       lcell(icell)%marked=lcell(i)%marked
       lcell(icell)%level=lcell(i)%level
       lcell(icell)%coor(1,1)=lcell(i)%coor(1,1)
       lcell(icell)%coor(2,1)=lcell(i)%coor(2,1)
       lcell(icell)%coor(3,1)=lcell(i)%coor(3,1)
       lcell(icell)%coor(1,2)=lcell(i)%coor(1,2)
       lcell(icell)%coor(2,2)=lcell(i)%coor(2,2)
       lcell(icell)%coor(3,2)=lcell(i)%coor(3,2)
       lcell(icell)%rsize=lcell(i)%rsize

       ineigh=lcell(i)%neigh(1)
       if(ineigh==0)then
          lcell(icell)%neigh(1)=0 
       else 
          lcell(icell)%neigh(1)=renum(ineigh)
       endif

       ineigh=lcell(i)%neigh(2)
       if(ineigh==0)then
          lcell(icell)%neigh(2)=0 
       else 
          lcell(icell)%neigh(2)=renum(ineigh)
       endif

       ineigh=lcell(i)%neigh(3)
       if(ineigh==0)then
          lcell(icell)%neigh(3)=0  
       else 
          lcell(icell)%neigh(3)=renum(ineigh) 
       endif

       ineigh=lcell(i)%neigh(4)
       if(ineigh==0)then
          lcell(icell)%neigh(4)=0 
       else 
          lcell(icell)%neigh(4)=renum(ineigh)
       endif

       ineigh=lcell(i)%neigh(5)
       if(ineigh==0)then
          lcell(icell)%neigh(5)=0 
       else 
          lcell(icell)%neigh(5)=renum(ineigh)
       endif

       ineigh=lcell(i)%neigh(6)
       if(ineigh==0)then
          lcell(icell)%neigh(6)=0 
       else 
          lcell(icell)%neigh(6)=renum(ineigh)
       endif

    enddo

    call dbgsub(lcell,ncell,remem,renum)

    ncell=nnew

    call divideCellSmoo(ncell,lcell,remem,ndim) 

    call divideFaceSmoo(ncell,lcell,remem)

    call compactSmoo(ncell,lcell,remem)

    call dbgsub(lcell,ncell,remem,renum)

    call memchk(2_ip,istat,memor_msh,'RENUM','smooth',renum)
    deallocate(renum,stat=istat)
    if(istat/=0) call memerr(2_ip,'RENUM','smooth',0_ip)
    call memchk(2_ip,istat,memor_msh,'LAY','smooth',lay)
    deallocate(lay,stat=istat)
    if(istat/=0) call memerr(2_ip,'LAY','smooth',0_ip)
    call memchk(2_ip,istat,memor_msh,'STACK','smooth',stack)
    deallocate(stack,stat=istat)
    if(istat/=0) call memerr(2_ip,'STACK','smooth',0_ip)
    call memchk(2_ip,istat,memor_msh,'REMEM','smooth',remem)
    deallocate(remem,stat=istat)
    if(istat/=0) call memerr(2_ip,'REMEM','smooth',0_ip)

  end subroutine smooth

  subroutine dbgconform(lcell,ncell)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ncell                    
    type(cell)                    :: lcell(*)
    integer(ip)                   :: i,j,neigh,level,neighlevel,k,kpos,kcont,tab(6),tab3(6),ilev,icell

    tab(1)=2
    tab(2)=1
    tab(3)=4
    tab(4)=3
    tab(5)=6
    tab(6)=5

    tab3(1)=2
    tab3(2)=1
    tab3(3)=4
    tab3(4)=3
    tab3(5)=6
    tab3(6)=5

    do i=1,ncell

       if(lcell(i)%marked<0)then

          level=lcell(i)%level

          do j=1,6      

             neigh=lcell(i)%neigh(j)

             if(neigh/=0)then

                neighlevel=lcell(neigh)%level

                if(neighlevel==level)then
                   !
                   !     Check coherence in neighbors
                   !
                   if(lcell(neigh)%neigh(tab(j))/=i)then
                      write(*,*)'error conform 1'   
                      stop
                   endif

                else if((neighlevel-level)==1)then
                   !
                   !     Check coherence from the coarsest
                   !
                   if(lcell(neigh)%neigh(tab3(j))/=i)then
                      write(*,*)'error conform 2'   
                      stop
                   endif
                   !
                   !     There should be four next cells that should see i as neigh is the smallest
                   !
                   kcont=1
                   kpos=neigh+1            
                   do k=1,8
                      if(lcell(kpos)%neigh(tab3(j))==i)then
                         kcont=kcont+1
                         if(kcont==4)exit
                      endif
                      kpos=kpos+1
                   enddo
                   if(k==9) then
                      write(*,*)'error conform 3'   
                      stop
                   endif
                endif
             endif
          enddo
       endif
    enddo


    !
    !     Check increasing level 
    !
    ilev=lcell(1)%level
    do icell=1,ncell
       if(lcell(icell)%level/=ilev)then
          ilev=ilev+1_ip
          if(ilev/=lcell(icell)%level)then
             write(*,*)'Error conform level'
             stop
          endif
       endif
    enddo

  end subroutine dbgconform

  subroutine dbgsub(lcell,ncell,remem,renum)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ncell                    
    type(cell)                    :: lcell(*)
    integer(ip)                   :: remem(*),renum(*),a12

    a12=1

  end subroutine dbgsub

  subroutine writeGidconf(ncell,lcell,lface,nface,coor,npoin,ndim,nnofa)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)        :: npoin,nface,ndim,nnofa
    integer(ip),intent(inout)     :: ncell
    integer(ip),intent(in)        :: lface(nnofa,nface)
    real(rp),intent(in)           :: coor(ndim,npoin)                    
    integer(ip)                   :: icont,iconte,icont1,icont2,icont3,icont4,icont5,icont6,icont7,icont8,i,icontt,index1
    integer(ip)                   :: neigh,level,pt1,pt2,j,k,ptnew,icontpp,iter,ncell8,jpos,jcont,tab4(3),npnew,ncont,ncell0
    integer(ip)                   :: tab(4,6),tab2(4,6),icell,ilevel,ineigh,jlevel,jneigh,ipoin,ipnew,npold,iplace,jplace
    integer(ip)                   :: tab3(6),tab5(4,6),tab6(4,6),ipo
    integer(ip),pointer           :: renum(:),renum2(:),lmark(:),lhang(:,:),lsize(:),lmarkc(:)
    real(rp), pointer             :: rsize(:),coorc(:,:)                   
    integer(4)                 :: istat
    type(cell)                    :: lcell(*)


    tab(1,1)=1
    tab(2,1)=5
    tab(3,1)=8
    tab(4,1)=4
    tab(1,2)=2
    tab(2,2)=3
    tab(3,2)=7
    tab(4,2)=6
    tab(1,3)=1
    tab(2,3)=2
    tab(3,3)=6
    tab(4,3)=5
    tab(1,4)=3
    tab(2,4)=4
    tab(3,4)=8
    tab(4,4)=7
    tab(1,5)=4
    tab(2,5)=3
    tab(3,5)=2
    tab(4,5)=1
    tab(1,6)=5
    tab(2,6)=6
    tab(3,6)=7
    tab(4,6)=8

    tab2(1,1)=2
    tab2(2,1)=6
    tab2(3,1)=7
    tab2(4,1)=3
    tab2(1,2)=1
    tab2(2,2)=4
    tab2(3,2)=8
    tab2(4,2)=5
    tab2(1,3)=4
    tab2(2,3)=3
    tab2(3,3)=7
    tab2(4,3)=8
    tab2(1,4)=2
    tab2(2,4)=1
    tab2(3,4)=5
    tab2(4,4)=6
    tab2(1,5)=8
    tab2(2,5)=7
    tab2(3,5)=6
    tab2(4,5)=5
    tab2(1,6)=1
    tab2(2,6)=2
    tab2(3,6)=3
    tab2(4,6)=4

    tab3(1)=2
    tab3(2)=1
    tab3(3)=4
    tab3(4)=3
    tab3(5)=6
    tab3(6)=5

    tab5(1,1)=1
    tab5(2,1)=5
    tab5(3,1)=4
    tab5(4,1)=8
    tab5(1,2)=2
    tab5(2,2)=6
    tab5(3,2)=3
    tab5(4,2)=7
    tab5(1,3)=1
    tab5(2,3)=5
    tab5(3,3)=2
    tab5(4,3)=6
    tab5(1,4)=4
    tab5(2,4)=8
    tab5(3,4)=3
    tab5(4,4)=7
    tab5(1,5)=1
    tab5(2,5)=4
    tab5(3,5)=2
    tab5(4,5)=3
    tab5(1,6)=5
    tab5(2,6)=8
    tab5(3,6)=6
    tab5(4,6)=7

    tab6(1,1)=2
    tab6(2,1)=6
    tab6(3,1)=3
    tab6(4,1)=7
    tab6(1,2)=1
    tab6(2,2)=5
    tab6(3,2)=4
    tab6(4,2)=8
    tab6(1,3)=4
    tab6(2,3)=8
    tab6(3,3)=3
    tab6(4,3)=7
    tab6(1,4)=1
    tab6(2,4)=5
    tab6(3,4)=2
    tab6(4,4)=6
    tab6(1,5)=5
    tab6(2,5)=8
    tab6(3,5)=6
    tab6(4,5)=7
    tab6(1,6)=1
    tab6(2,6)=4
    tab6(3,6)=2
    tab6(4,6)=3

    !
    !     Allocate renum and renum2
    !
    ncell8=ncell*8 
    allocate(renum(ncell8),stat=istat)
    call memchk(zero,istat,memor_msh,'RENUM','writeGidconf',renum)
    allocate(renum2(ncell8),stat=istat)
    call memchk(zero,istat,memor_msh,'RENUM2','writeGidconf',renum2)
    !
    !     Initialize renum
    !
    do i=1,ncell8
       renum(i)=i
    enddo
    !
    !     Conformize the cells
    !
    do
       !
       !     Initialize the "something done" check 
       !
       iter=0
       !
       !     Loop on the cells
       !
       do icell=1,ncell

          ilevel=lcell(icell)%level
          !
          !     Loop on the neighbors
          !
          do i=1,6

             ineigh=lcell(icell)%neigh(i)

             if(ineigh==0)cycle

             jlevel=lcell(ineigh)%level

             if(jlevel==ilevel)then
                !
                !     Loop on the points of the faces
                !
                do k=1,4

                   pt1=(icell-1)*8+tab(k,i) 
                   pt2=(ineigh-1)*8+tab2(k,i)
                   if(renum(pt1)<renum(pt2))then
                      renum(pt2)=renum(pt1)
                      iter=1
                   else if(renum(pt2)<renum(pt1))then   
                      renum(pt1)=renum(pt2)
                      iter=1
                   endif

                enddo

             else if(ilevel<jlevel)then

                jpos=ineigh+1_ip
                jcont=1_ip
                do j=1,8
                   if(lcell(jpos)%neigh(tab3(i))==icell)then
                      tab4(jcont)=jpos
                      jcont=jcont+1_ip
                      if(jcont==4)exit 
                   endif
                   jpos=jpos+1  
                enddo

                pt1=(icell-1)*8+tab5(1,i)
                pt2=(ineigh-1)*8+tab6(1,i)
                if(renum(pt1)<renum(pt2))then
                   renum(pt2)=renum(pt1)
                   iter=1
                else if(renum(pt2)<renum(pt1))then   
                   renum(pt1)=renum(pt2)
                   iter=1
                endif

                pt1=(icell-1)*8+tab5(2,i)
                pt2=(tab4(1)-1)*8+tab6(2,i)
                if(renum(pt1)<renum(pt2))then
                   renum(pt2)=renum(pt1)
                   iter=1
                else if(renum(pt2)<renum(pt1))then   
                   renum(pt1)=renum(pt2)
                   iter=1
                endif

                pt1=(icell-1)*8+tab5(3,i)
                pt2=(tab4(2)-1)*8+tab6(3,i)
                if(renum(pt1)<renum(pt2))then
                   renum(pt2)=renum(pt1)
                   iter=1
                else if(renum(pt2)<renum(pt1))then   
                   renum(pt1)=renum(pt2)
                   iter=1
                endif

                pt1=(icell-1)*8+tab5(4,i)
                pt2=(tab4(3)-1)*8+tab6(4,i)
                if(renum(pt1)<renum(pt2))then
                   renum(pt2)=renum(pt1)
                   iter=1
                else if(renum(pt2)<renum(pt1))then   
                   renum(pt1)=renum(pt2)
                   iter=1
                endif

             endif
          enddo
       enddo

       do i=1,ncell8
          renum(i)=renum(renum(i))
       enddo

       if(iter==0)exit
       exit  

    enddo
    !
    !     Get the number and the new renumbering of the conformed points
    !
    icont=0
    do i=1,ncell8
       if(renum(i)==i)then
          icont=icont+1
          renum2(i)=icont
       endif
    enddo
    !
    !     Remember the conforming point number
    ! 
    npnew=icont
    !
    !     Put the new renumbering in renum
    !
    do i=1,ncell8
       renum(i)=renum2(renum(i))
    enddo
    !
    !     Allocate lmark for the hanging nodes
    !
    allocate(lmark(npnew),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','writeGidconf',lmark)
    allocate(lhang(4,npnew),stat=istat)
    call memchk(zero,istat,memor_msh,'LHANG','writeGidconf',lhang)
    !
    !     Get the hanging nodes 
    !
    do icell=1,ncell
       ilevel=lcell(icell)%level
       !
       !     In X -
       !  
       ineigh=lcell(icell)%neigh(1)
       if(ineigh/=0)then
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(2)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 1'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(2)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo
             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+7_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
                !lhang(1,ipoin)=renum(8*(ineigh-1)+2_ip) 
                !lhang(2,ipoin)=renum(8*(tab4(1)-1)+6_ip) 
                !lhang(3,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
                !lhang(4,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+6_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+2_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+6_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+3_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+2_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+6_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
          endif
       endif
       !
       !     In X +
       !  
       ineigh=lcell(icell)%neigh(2)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(1)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 2'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(1)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo
             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+8_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+4_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+8_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+5_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+4_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+4_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+8_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+8_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+8_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+8_ip) 
             endif
          endif
       endif
       !
       !     In Y -
       !  
       ineigh=lcell(icell)%neigh(3)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(4)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 3'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(4)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo
             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+7_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+8_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+3_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
          endif
       endif
       !
       !     In Y +
       !  
       ineigh=lcell(icell)%neigh(4)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(3)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 4'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(3)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo
             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+6_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+6_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+5_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+2_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+6_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+6_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+6_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+6_ip) 
             endif
          endif
       endif
       !
       !     In Z -
       !  
       ineigh=lcell(icell)%neigh(5)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(6)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 5'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(6)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo
             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+7_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+6_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+8_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+6_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+6_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+6_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
          endif
       endif
       !
       !     In Z +
       !  
       ineigh=lcell(icell)%neigh(6)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(5)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 6'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(5)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo
             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+3_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+4_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+3_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+4_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+4_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+2_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+3_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+3_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+3_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+3_ip) 
             endif
          endif
       endif
    enddo
    !
    !     Get the inner cells 
    !
    allocate(lmarkc(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKC','takeout',lmarkc)
    !call carout(ncell,lcell,lmarkc)
    do icell=1,ncell
       lmarkc(icell)=1_ip
    enddo
    !
    !     Mark the exterior points
    !
    do icell=1,ncell 
       if(lmarkc(icell)==1)then
          iplace=8*(icell-1)+1 
          ipoin=renum(iplace)
          renum2(ipoin)=1
          iplace=iplace+1 
          ipoin=renum(iplace)
          renum2(ipoin)=1
          iplace=iplace+1 
          ipoin=renum(iplace)
          renum2(ipoin)=1
          iplace=iplace+1 
          ipoin=renum(iplace)
          renum2(ipoin)=1
          iplace=iplace+1 
          ipoin=renum(iplace)
          renum2(ipoin)=1
          iplace=iplace+1 
          ipoin=renum(iplace)
          renum2(ipoin)=1
          iplace=iplace+1 
          ipoin=renum(iplace)
          renum2(ipoin)=1
          iplace=iplace+1 
          ipoin=renum(iplace)
          renum2(ipoin)=1
       endif
    enddo
    !
    !     Mark the interior points already marked
    !
    do icell=1,ncell 
       if(lmarkc(icell)/=1)then
          iplace=8*(icell-1)+1 
          ipoin=renum(iplace)
          if(renum2(ipoin)==1)then
             renum2(ipoin)=2
          endif
          iplace=iplace+1 
          ipoin=renum(iplace)
          if(renum2(ipoin)==1)then
             renum2(ipoin)=2
          endif
          ipoin=renum(iplace)
          if(renum2(ipoin)==1)then
             renum2(ipoin)=2
          endif
          ipoin=renum(iplace)
          if(renum2(ipoin)==1)then
             renum2(ipoin)=2
          endif
          ipoin=renum(iplace)
          if(renum2(ipoin)==1)then
             renum2(ipoin)=2
          endif
          ipoin=renum(iplace)
          if(renum2(ipoin)==1)then
             renum2(ipoin)=2
          endif
          ipoin=renum(iplace)
          if(renum2(ipoin)==1)then
             renum2(ipoin)=2
          endif
          ipoin=renum(iplace)
          if(renum2(ipoin)==1)then
             renum2(ipoin)=2
          endif
       endif
    enddo
    !
    !     Get the inner points to be deleted
    !
    do icell=1,ncell 
       if(lmarkc(icell)/=1)then     
          do j=1,8
             iplace=8*(icell-1)+j
             ipoin=renum(iplace)
             !if(renum2(ipoin)/=2)then   !NOT NECESSARY AS MARKED BY OUTER  
             renum(iplace)=-ipoin     ! ELEMENTS
             !endif
          enddo
       endif
    enddo
    !
    !     Clean up renum2     
    ! 
    do ipoin=1,npnew
       renum2(ipoin)=0_ip
    enddo
    !
    !     Transfer the valid points from cells to points 
    !
    do ipoin=1,ncell8 
       ipnew=renum(ipoin)
       if(ipnew>0)then
          renum2(ipnew)=1_ip
       endif
    enddo
    !
    !     Renumber the points
    !
    icont=0_ip
    do ipoin=1,npnew
       if(renum2(ipoin)==1)then
          icont=icont+1_ip 
          renum2(ipoin)=icont
       endif
    enddo
    npold=npnew
    npnew=icont
    !
    !     Renumber the hanging nodes
    !
    icont=0_ip
    do ipoin=1,npold
       if(renum2(ipoin)>0)then   
          icont=icont+1_ip
          if(lmark(ipoin)==1)then

             lmark(icont)=1_ip
            ! HN lhang(1,icont)=renum2(lhang(1,ipoin))
            ! HN lhang(2,icont)=renum2(lhang(2,ipoin))
             if(lhang(3,ipoin)/=0)then
               ! HN lhang(3,icont)=renum2(lhang(3,ipoin))
             else 
               ! HN lhang(3,icont)=0_ip
             endif
             if(lhang(4,ipoin)/=0)then
               ! HN lhang(4,icont)=renum2(lhang(4,ipoin))
             else 
               ! HN lhang(4,icont)=0_ip
             endif
          else
             lmark(icont)=0_ip
            ! HN lhang(1,icont)=0_ip
            ! HN lhang(2,icont)=0_ip
            ! HN lhang(3,icont)=0_ip
            ! HN lhang(4,icont)=0_ip
          endif
       endif
    enddo
    !     
    !     Renumber the cells
    !
    icont=0_ip
    do icell=1,ncell
       if(lmarkc(icell)==1)then
          icont=icont+1_ip
          lmarkc(icell)=icont
       endif
    enddo
    !
    !     Compact the lcell array
    !
    ncont=0_ip
    ncell0=ncell
    do icell=1,ncell0
       if(lmarkc(icell)>0)then
          ncont=ncont+1
          lcell(ncont)%marked    = lcell(icell)%marked
          lcell(ncont)%level     = lcell(icell)%level
          lcell(ncont)%coor(1,1) = lcell(icell)%coor(1,1)
          lcell(ncont)%coor(2,1) = lcell(icell)%coor(2,1)
          lcell(ncont)%coor(3,1) = lcell(icell)%coor(3,1)
          lcell(ncont)%coor(1,2) = lcell(icell)%coor(1,2)
          lcell(ncont)%coor(2,2) = lcell(icell)%coor(2,2)
          lcell(ncont)%coor(3,2) = lcell(icell)%coor(3,2)
          lcell(ncont)%rsize     = lcell(icell)%rsize
          !
          !     Renumber the neighbors
          !
          ineigh=lcell(icell)%neigh(1)
          if(ineigh==0)then
             lcell(ncont)%neigh(1)=0 
          else 
             lcell(ncont)%neigh(1)=lmarkc(ineigh) 
          endif

          ineigh=lcell(icell)%neigh(2)
          if(ineigh==0)then
             lcell(ncont)%neigh(2)=0 
          else 
             lcell(ncont)%neigh(2)=lmarkc(ineigh)
          endif

          ineigh=lcell(icell)%neigh(3)
          if(ineigh==0)then
             lcell(ncont)%neigh(3)=0 
          else 
             lcell(ncont)%neigh(3)=lmarkc(ineigh)
          endif

          ineigh=lcell(icell)%neigh(4)
          if(ineigh==0)then
             lcell(ncont)%neigh(4)=0 
          else 
             lcell(ncont)%neigh(4)=lmarkc(ineigh)
          endif

          ineigh=lcell(icell)%neigh(5)
          if(ineigh==0)then
             lcell(ncont)%neigh(5)=0 
          else 
             lcell(ncont)%neigh(5)=lmarkc(ineigh)
          endif

          ineigh=lcell(icell)%neigh(6)
          if(ineigh==0)then
             lcell(ncont)%neigh(6)=0 
          else 
             lcell(ncont)%neigh(6)=lmarkc(ineigh)
          endif
          !
          !     Compact renum 
          !
          iplace=8*(ncont-1_ip)+1_ip
          jplace=8*(icell-1_ip)+1_ip
          renum(iplace)=renum2(renum(jplace))
          iplace=iplace+1_ip
          jplace=jplace+1_ip 
          renum(iplace)=renum2(renum(jplace))
          iplace=iplace+1_ip
          jplace=jplace+1_ip 
          renum(iplace)=renum2(renum(jplace))
          iplace=iplace+1_ip
          jplace=jplace+1_ip 
          renum(iplace)=renum2(renum(jplace))
          iplace=iplace+1_ip
          jplace=jplace+1_ip 
          renum(iplace)=renum2(renum(jplace))
          iplace=iplace+1_ip
          jplace=jplace+1_ip 
          renum(iplace)=renum2(renum(jplace))
          iplace=iplace+1_ip
          jplace=jplace+1_ip 
          renum(iplace)=renum2(renum(jplace))
          iplace=iplace+1_ip
          jplace=jplace+1_ip 
          renum(iplace)=renum2(renum(jplace))
       endif
    enddo

    ncell=ncont  
    !
    !     Get the coor array
    !   
    allocate(coorc(ndim,npnew),stat=istat)
    call memchk(zero,istat,memor_msh,'COORC','writeGidconf',coorc)

    icont=0_ip
    do i=1,ncell

       icont=icont+1_ip
       ipoin=renum(icont)
       coorc(1,ipoin)=lcell(i)%coor(1,1) 
       coorc(2,ipoin)=lcell(i)%coor(2,1) 
       coorc(3,ipoin)=lcell(i)%coor(3,1) 

       icont=icont+1_ip
       ipoin=renum(icont)
       coorc(1,ipoin)=lcell(i)%coor(1,2) 
       coorc(2,ipoin)=lcell(i)%coor(2,1) 
       coorc(3,ipoin)=lcell(i)%coor(3,1) 

       icont=icont+1_ip
       ipoin=renum(icont)
       coorc(1,ipoin)=lcell(i)%coor(1,2) 
       coorc(2,ipoin)=lcell(i)%coor(2,2) 
       coorc(3,ipoin)=lcell(i)%coor(3,1) 

       icont=icont+1_ip
       ipoin=renum(icont)
       coorc(1,ipoin)=lcell(i)%coor(1,1) 
       coorc(2,ipoin)=lcell(i)%coor(2,2) 
       coorc(3,ipoin)=lcell(i)%coor(3,1) 

       icont=icont+1_ip
       ipoin=renum(icont)
       coorc(1,ipoin)=lcell(i)%coor(1,1) 
       coorc(2,ipoin)=lcell(i)%coor(2,1) 
       coorc(3,ipoin)=lcell(i)%coor(3,2) 

       icont=icont+1_ip
       ipoin=renum(icont)
       coorc(1,ipoin)=lcell(i)%coor(1,2) 
       coorc(2,ipoin)=lcell(i)%coor(2,1) 
       coorc(3,ipoin)=lcell(i)%coor(3,2) 

       icont=icont+1_ip
       ipoin=renum(icont)
       coorc(1,ipoin)=lcell(i)%coor(1,2) 
       coorc(2,ipoin)=lcell(i)%coor(2,2) 
       coorc(3,ipoin)=lcell(i)%coor(3,2) 

       icont=icont+1_ip
       ipoin=renum(icont)
       coorc(1,ipoin)=lcell(i)%coor(1,1) 
       coorc(2,ipoin)=lcell(i)%coor(2,2) 
       coorc(3,ipoin)=lcell(i)%coor(3,2) 

    enddo
    !
    !     Write the cells
    !
    open(unit=50,file='cartGid.msh',status='unknown')
    rewind 50
    write(50,1)
    write(50,2)
    write(50,3)

    do ipoin=1,npnew
       write(50,100)ipoin,coorc(1,ipoin),coorc(2,ipoin),coorc(3,ipoin)
    enddo

    write(50,4)
    write(50,5)

    call memchk(2_ip,istat,memor_msh,'COORC','writeGidconf',coorc)
    deallocate(coorc,stat=istat)
    if(istat/=0) call memerr(2_ip,'COORC','writeGidconf',0_ip)
    !
    !     Output connectivity
    !
    icont=0
    icontpp=0
    do i=1,ncell

       icontpp=icontpp+1 
       icont1=renum(icontpp)
       icontpp=icontpp+1 
       icont2=renum(icontpp)  
       icontpp=icontpp+1 
       icont3=renum(icontpp)  
       icontpp=icontpp+1 
       icont4=renum(icontpp)  
       icontpp=icontpp+1 
       icont5=renum(icontpp)  
       icontpp=icontpp+1 
       icont6=renum(icontpp)  
       icontpp=icontpp+1 
       icont7=renum(icontpp)  
       icontpp=icontpp+1 
       icont8=renum(icontpp)  
       !index1=-lcell(i)%marked
       index1=lcell(i)%level 
       !if(index1==3)then
       icont=icont+1
       write(50,200) icont,icont1,icont2,icont3,icont4,icont5,icont6,icont7,icont8,index1
       !write(50,200) icont,icont1,icont2,icont3,icont4,icont5,icont6,icont7,icont8,lmarkc(i)+1
       !endif
    enddo
    write(50,6)

    icontt=icont
    !
    !     Output the faces of the triangulation
    !
    write(50,7)
    write(50,2)
    write(50,3)

    icont=npnew
    do i=1,npoin
       icont=icont+1
       write(50,100)icont,coor(1,i),coor(2,i),coor(3,i)
    enddo
    write(50,4) 
    write(50,5) 

    icont=icontt
    do i=1, nface
       icont=icont+1
       write(50,300)icont,lface(1,i)+npnew,lface(2,i)+npnew,lface(3,i)+npnew     
    enddo
    write(50,6)      
    close(50)
    !
    !     Write the marked points
    !
    open(unit=60,file='carthan.msh',status='unknown')
    rewind 60

    do ipoin=1,npnew
       if(lmark(ipoin)==1)then
          write(60,400)ipoin,lhang(1,ipoin),lhang(2,ipoin),lhang(3,ipoin),lhang(4,ipoin)
       endif
    enddo

    close(60)
    !
    !     Output the size
    !
    open(unit=70,file='cartGid.res',status='unknown')
    rewind 70

    allocate(rsize(npnew),stat=istat)
    call memchk(zero,istat,memor_msh,'RSIZE','writeGidconf',rsize)
    allocate(lsize(npnew),stat=istat)
    call memchk(zero,istat,memor_msh,'LSIZE','writeGidconf',lsize)

    do icell=1,ncell
       do ipo=1,8
          ipoin=renum(8*(icell-1)+ipo)
          rsize(ipoin)=rsize(ipoin)+lcell(icell)%rsize
          lsize(ipoin)=lsize(ipoin)+1
       enddo
    enddo

    do  ipoin=1,npnew
       if(lsize(ipoin)==0)then
          write(*,*)'Error in writeGidconf, lsize=0 at ipoin:',ipoin
          stop 
       endif
       rsize(ipoin)=rsize(ipoin)/real(lsize(ipoin))
    enddo


    write(70,10)
    write(70,11)
    write(70,13)
    do  ipoin=1,npnew
       write(70,500)ipoin,rsize(ipoin)
    enddo
    write(70,14)
    write(70,15)
    close(70)

    call memchk(2_ip,istat,memor_msh,'LSIZE','writeGidconf',lsize)
    deallocate(lsize,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSIZE','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RSIZE','writeGidconf',rsize)
    deallocate(rsize,stat=istat)
    if(istat/=0) call memerr(2_ip,'RSIZE','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARKC','writeGidconf',lmarkc)
    deallocate(lmarkc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARKC','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHANG','writeGidconf',lhang)
    deallocate(lhang,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHANG','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','writeGidconf',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RENUM','writeGidconf',renum)
    deallocate(renum,stat=istat)
    if(istat/=0) call memerr(2_ip,'RENUM','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RENUM2','writeGidconf',renum2)
    deallocate(renum2,stat=istat)
    if(istat/=0) call memerr(2_ip,'RENUM2','writeGidconf',0_ip)


1   format('MESH dimension 3 ElemType  Hexahedra Nnode 8')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(10i10)
300 format(4i10)
400 format(5i10)
500 format(i10,e20.10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
7   format('MESH dimension 3 ElemType Triangle Nnode 3')


10  format('GID Post Results File 1.0')
11  format('Result "Size" "Analysis/time" 1 Scalar OnNodes')
13  format('Values')
14  format('End Values')
15  format('   ')

  end subroutine writeGidconf

  subroutine intercart(coor,npoin,ndim,lface,nface,nnofa,lcell,ncell,lcart,ptoel1,ptoel2,rtol)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)        :: ncell,npoin,nface,ndim,nnofa
    integer(ip),intent(in)        :: lface(nnofa,nface),ptoel1(*),ptoel2(*)
    real(rp),intent(in)           :: coor(ndim,npoin),rtol                   
    integer(ip),intent(inout)     :: lcart(npoin) 
    integer(ip),pointer           :: lstack(:) 
    logical(lg),pointer           :: lmark(:)
    integer(4)                    :: istat
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: nstack,ip1,ip2,j,ie,ipoin,ielem,icart,ifirst 
    !
    !     This sub finds the cartesian cell containing the points of lface for each point
    !     and store it in lcart. The point to element graph is given in ptoel
    !
    allocate(lstack(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'lstack','intercart',lstack)
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'lstack','intercart',lstack)
    do ipoin=1,npoin
       lmark(ipoin)=.false.
    end do
    !
    !     Initialize ifirst
    !
    ifirst=max(ncell/2_ip,1_ip)
    !
    !     Initialize the stack
    !
    lstack(1)=1_ip
    nstack=1_ip
    lcart(1)=ifirst
    lmark(1)=.true.
    !
    !     Main loop on stack
    !
    do ip1=1,npoin

       ipoin=lstack(ip1)
       !
       !     Do we have various connex components
       !
       if(ipoin==0)then
          do j=1,npoin
             if(.not.lmark(j)) then
                exit
             endif
          enddo

          if(j==npoin+1)then
             write(*,*)'Error intercart, new point  not found'
             stop
          endif
          ipoin=j
          nstack=nstack+1
          lstack(nstack)=ipoin
          lcart(ipoin)=ifirst
          lmark(ipoin)=.true.
       endif

       icart=lcart(ipoin)
       call gtelem(ipoin,coor,npoin,ndim,lcell,ncell,icart,rtol)
       if(icart==0)then
          write(*,*)'Error in intercart, point out of domain'
          stop 
       endif
       !
       !     Remember the element containing the point
       !
       lcart(ipoin)=icart
       !
       !     Add the neighbors to the stack and initialize the initial guess
       !
       do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(ie)
          do j=1,nnofa
             ip2=lface(j,ielem)
             if(.not.lmark(ip2)) then
                lmark(ip2)=.true.      
                nstack=nstack+1
                lstack(nstack)=ip2
                lcart(ip2)=icart
             endif
          enddo
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','intercart',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','intercart',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','intercart',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','intercart',0_ip)

  end subroutine intercart

  subroutine gtelem(ipoin,coor,npoin,ndim,lcell,ncell,icart,rtol)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)        :: ncell,npoin,ndim,ipoin
    real(rp),intent(in)           :: coor(ndim,npoin),rtol     
    integer(ip), intent(inout)    :: icart               
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: icell
    real(rp)                      :: rx,ry,rz,rxmax,rymax,rzmax

    rx=coor(1,ipoin)  
    ry=coor(2,ipoin)  
    rz=coor(3,ipoin)  

    do 

       if(icart==0 .or. icart>ncell)then
          write(*,*)'Error gtelem, point out of domain :',ipoin
          write(*,*)coor(1,ipoin),coor(2,ipoin),coor(3,ipoin)

          rxmax=lcell(1)%coor(1,2) 
          rymax=lcell(1)%coor(2,2) 
          rzmax=lcell(1)%coor(3,2) 
          do icell=2,ncell
             if(rxmax<lcell(icell)%coor(1,2))rxmax=lcell(icell)%coor(1,2)
             if(rymax<lcell(icell)%coor(2,2))rymax=lcell(icell)%coor(2,2)
             if(rzmax<lcell(icell)%coor(3,2))rzmax=lcell(icell)%coor(3,2)
          enddo
          write(*,*)'Extent of cartesian mesh  xmin,ymin,zmin  xmax,ymax,zmax:'
          write(*,*)lcell(1)%coor(1,1),lcell(1)%coor(2,1),lcell(1)%coor(3,1)
          write(*,*)rxmax,rymax,rzmax
          icart=0_ip
          return  
          !stop
       endif

       if(rx-lcell(icart)%coor(1,1)<-rtol)then
          icart=lcell(icart)%neigh(1)    
          cycle
       endif
       if(ry-lcell(icart)%coor(2,1)<-rtol)then
          icart=lcell(icart)%neigh(3)    
          cycle
       endif
       if(rz-lcell(icart)%coor(3,1)<-rtol)then
          icart=lcell(icart)%neigh(5)    
          cycle
       endif
       if(rx-lcell(icart)%coor(1,2)>=rtol)then
          icart=lcell(icart)%neigh(2)    
          cycle
       endif
       if(ry-lcell(icart)%coor(2,2)>=rtol)then
          icart=lcell(icart)%neigh(4)    
          cycle
       endif
       if(rz-lcell(icart)%coor(3,2)>=rtol)then
          icart=lcell(icart)%neigh(6)    
          cycle
       endif

       exit 

    enddo

  end subroutine gtelem

  subroutine gtsize(ncell,lcell,npoin,rsize,lcart,ndim,coor)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only              : memor_msh
    implicit none
    integer(ip),intent(in)        :: ncell,npoin,ndim
    integer(ip),intent(in)        :: lcart(npoin)
    type(cell),intent(in)         :: lcell(ncell)
    real(rp),intent(inout)        :: rsize(npoin),coor(ndim,npoin)
    integer(ip)                   :: ipoin,icell,ilevel,ineigh,nx,ny,nz
    real(rp)                      :: xc,yc,zc,c05,rcsize,c20,c00,c10,c23,c43 
    real(rp)                      :: hsize,dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3 
    real(rp)                      :: dxs,dys,dzs,dx,dy,dz
    real(rp)                      :: cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,ct
    !
    !     This sub transfers the size of the bg grid to the point
    !
    c05=0.5d+00
    c20=2.0d+00
    c00=0.0d+00
    c10=1.0d+00
    c23=2.0d+00/3.0d+00
    c43=4.0d+00/3.0d+00

    do ipoin=1,npoin
       !
       !     Get the cell number
       !
       icell=lcart(ipoin)
       ilevel=lcell(icell)%level
       !
       !     Get the cell center
       ! 
       xc=(lcell(icell)%coor(1,1)+lcell(icell)%coor(1,2))*c05
       yc=(lcell(icell)%coor(2,1)+lcell(icell)%coor(2,2))*c05
       zc=(lcell(icell)%coor(3,1)+lcell(icell)%coor(3,2))*c05
       !
       !     Get the cell size
       !
       !rcsize=(lcell(icell)%coor(1,2)-xc)+(lcell(icell)%coor(2,2)-yc)+(lcell(icell)%coor(3,2)-zc)
       !rcsize=lcell(icell)%coor(1,2)-lcell(icell)%coor(1,1)
       hsize=c10/(lcell(icell)%coor(1,2)-lcell(icell)%coor(1,1))
       rcsize=lcell(icell)%rsize
       !
       !     Get the neighbors
       !
       nx=0_ip
       ineigh=lcell(icell)%neigh(1)
       if(ineigh/=0)then 
          if(lcell(ineigh)%level<ilevel)then
             dx1=(rcsize-lcell(ineigh)%rsize)*c23
          else if(lcell(ineigh)%level>ilevel)then
             dx1=(rcsize-lcell(ineigh)%rsize)*c43
          else 
             dx1=rcsize-lcell(ineigh)%rsize
          endif
          nx=nx+1_ip
       else
          dx1=c00
       endif

       ineigh=lcell(icell)%neigh(2) 
       if(ineigh/=0)then 
          if(lcell(ineigh)%level<ilevel)then
             dx2=(lcell(ineigh)%rsize-rcsize)*c23
          else if(lcell(ineigh)%level>ilevel)then
             dx2=(lcell(ineigh)%rsize-rcsize)*c43
          else 
             dx2=lcell(ineigh)%rsize-rcsize
          endif
          nx=nx+1_ip
       else
          dx2=c00
       endif

       ny=0_ip
       ineigh=lcell(icell)%neigh(3)
       if(ineigh/=0)then 
          if(lcell(ineigh)%level<ilevel)then
             dy1=(rcsize-lcell(ineigh)%rsize)*c23
          else if(lcell(ineigh)%level>ilevel)then
             dy1=(rcsize-lcell(ineigh)%rsize)*c43
          else 
             dy1=rcsize-lcell(ineigh)%rsize
          endif
          ny=ny+1_ip 
       else
          dy1=c00
       endif

       ineigh=lcell(icell)%neigh(4) 
       if(ineigh/=0)then 
          if(lcell(ineigh)%level<ilevel)then
             dy2=(lcell(ineigh)%rsize-rcsize)*c23
          else if(lcell(ineigh)%level>ilevel)then
             dy2=(lcell(ineigh)%rsize-rcsize)*c43
          else 
             dy2=lcell(ineigh)%rsize-rcsize
          endif
          ny=ny+1_ip
       else
          dy2=c00
       endif

       nz=0_ip
       ineigh=lcell(icell)%neigh(5)
       if(ineigh/=0)then 
          if(lcell(ineigh)%level<ilevel)then
             dz1=(rcsize-lcell(ineigh)%rsize)*c23
          else if(lcell(ineigh)%level>ilevel)then
             dz1=(rcsize-lcell(ineigh)%rsize)*c43
          else 
             dz1=rcsize-lcell(ineigh)%rsize
          endif
          nz=nz+1_ip
       else
          dz1=c00
       endif

       ineigh=lcell(icell)%neigh(6) 
       if(ineigh/=0)then 
          if(lcell(ineigh)%level<ilevel)then
             dz2=(lcell(ineigh)%rsize-rcsize)*c23
          else if(lcell(ineigh)%level>ilevel)then
             dz2=(lcell(ineigh)%rsize-rcsize)*c43
          else 
             dz2=lcell(ineigh)%rsize-rcsize
          endif
          nz=nz+1_ip
       else
          dz2=c00
       endif
       !
       !     Compute gradient of the size   
       !
       if(nx>0)then  
          dxs=(dx2+dx1)*hsize/real(nx)
       else
          dxs=c00
       endif
       if(ny>0)then  
          dys=(dy2+dy1)*hsize/real(ny)
       else
          dys=c00
       endif
       if(nz>0)then  
          dzs=(dz2+dz1)*hsize/real(nz)
       else
          dzs=c00
       endif
       !
       !     Compute dx
       !
       dx=coor(1,ipoin)-xc 
       dy=coor(2,ipoin)-yc 
       dz=coor(3,ipoin)-zc 
       !
       !     Get the final size
       !
       !rsize(ipoin)=rcsize+dxs*dx+dys*dy+dzs*dz
       rsize(ipoin)=rcsize
       !
       !     DBG
       !
       if(rsize(ipoin)<c00)then
          write(*,*)'Error in gtsize, negative size at ipoin:',ipoin
          stop
       endif

       !rsize(ipoin)=1.0d-04
    enddo

  end subroutine gtsize

  subroutine gtsiz2(ncell,lcell,npoin,rsize,lcart,ipoin,ndim,coor,rsuni)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)        :: ncell,npoin,ndim,ipoin
    real(rp),intent(in)           :: rsuni
    integer(ip),intent(in)        :: lcart(npoin)
    type(cell),intent(in)         :: lcell(ncell)
    real(rp),intent(inout)        :: rsize(npoin)
    real(rp),intent(in)           :: coor(ndim,npoin)
    integer(ip)                   :: icell,ilevel,ineigh,nx,ny,nz
    real(rp)                      :: xc,yc,zc,c05,rcsize,c20,c10,c23,c43,c00 
    real(rp)                      :: dxs,dys,dzs,dx,dy,dz,dx1,dx2,dy1,dy2,dz1,dz2,hsize
    !
    !     This sub transfers the size of the bg grid to the point
    !
    c05=0.5d+00
    c20=2.0d+00
    c00=0.0d+00
    c10=1.0d+00
    c23=2.0d+00/3.0d+00
    c43=4.0d+00/3.0d+00
    !
    !     Get the cell number
    !
    icell=lcart(ipoin)
    !
    !     Do we have a point outside the volume
    !
    if(icell==0)then
       rsize(ipoin)=rsuni
       return
    endif


    ilevel=lcell(icell)%level
    !
    !     Get the cell center
    ! 
    xc=(lcell(icell)%coor(1,1)+lcell(icell)%coor(1,2))*c05
    yc=(lcell(icell)%coor(2,1)+lcell(icell)%coor(2,2))*c05
    zc=(lcell(icell)%coor(3,1)+lcell(icell)%coor(3,2))*c05
    !
    !     Get the cell size
    !
    !rcsize=(lcell(icell)%coor(1,2)-xc)+(lcell(icell)%coor(2,2)-yc)+(lcell(icell)%coor(3,2)-zc)
    hsize=c10/(lcell(icell)%coor(1,2)-lcell(icell)%coor(1,1))
    rcsize=lcell(icell)%rsize
    !
    !     Get the neighbors
    !
    nx=0_ip
    ineigh=lcell(icell)%neigh(1)
    if(ineigh/=0)then 
       if(lcell(ineigh)%level<ilevel)then
          dx1=(rcsize-lcell(ineigh)%rsize)*c23
       else if(lcell(ineigh)%level>ilevel)then
          dx1=(rcsize-lcell(ineigh)%rsize)*c43
       else 
          dx1=rcsize-lcell(ineigh)%rsize
       endif
       nx=nx+1_ip
    else
       dx1=c00
    endif

    ineigh=lcell(icell)%neigh(2) 
    if(ineigh/=0)then 
       if(lcell(ineigh)%level<ilevel)then
          dx2=(lcell(ineigh)%rsize-rcsize)*c23
       else if(lcell(ineigh)%level>ilevel)then
          dx2=(lcell(ineigh)%rsize-rcsize)*c43
       else 
          dx2=lcell(ineigh)%rsize-rcsize
       endif
       nx=nx+1_ip
    else
       dx2=c00
    endif

    ny=0_ip
    ineigh=lcell(icell)%neigh(3)
    if(ineigh/=0)then 
       if(lcell(ineigh)%level<ilevel)then
          dy1=(rcsize-lcell(ineigh)%rsize)*c23
       else if(lcell(ineigh)%level>ilevel)then
          dy1=(rcsize-lcell(ineigh)%rsize)*c43
       else 
          dy1=rcsize-lcell(ineigh)%rsize
       endif
       ny=ny+1_ip 
    else
       dy1=c00
    endif

    ineigh=lcell(icell)%neigh(4) 
    if(ineigh/=0)then 
       if(lcell(ineigh)%level<ilevel)then
          dy2=(lcell(ineigh)%rsize-rcsize)*c23
       else if(lcell(ineigh)%level>ilevel)then
          dy2=(lcell(ineigh)%rsize-rcsize)*c43
       else 
          dy2=lcell(ineigh)%rsize-rcsize
       endif
       ny=ny+1_ip
    else
       dy2=c00
    endif

    nz=0_ip
    ineigh=lcell(icell)%neigh(5)
    if(ineigh/=0)then 
       if(lcell(ineigh)%level<ilevel)then
          dz1=(rcsize-lcell(ineigh)%rsize)*c23
       else if(lcell(ineigh)%level>ilevel)then
          dz1=(rcsize-lcell(ineigh)%rsize)*c43
       else 
          dz1=rcsize-lcell(ineigh)%rsize
       endif
       nz=nz+1_ip
    else
       dz1=c00
    endif

    ineigh=lcell(icell)%neigh(6) 
    if(ineigh/=0)then 
       if(lcell(ineigh)%level<ilevel)then
          dz2=(lcell(ineigh)%rsize-rcsize)*c23
       else if(lcell(ineigh)%level>ilevel)then
          dz2=(lcell(ineigh)%rsize-rcsize)*c43
       else 
          dz2=lcell(ineigh)%rsize-rcsize
       endif
       nz=nz+1_ip
    else
       dz2=c00
    endif
    !
    !     Compute gradient of the size   
    ! 
    dxs=(dx2+dx1)*hsize/real(nx)
    dys=(dy2+dy1)*hsize/real(ny)
    dzs=(dz2+dz1)*hsize/real(nz)
    !
    !     Compute dx
    !
    dx=coor(1,ipoin)-xc 
    dy=coor(2,ipoin)-yc 
    dz=coor(3,ipoin)-zc 
    !
    !     Get the final size
    !
    !rsize(ipoin)=rcsize+dxs*dx+dys*dy+dzs*dz
    rsize(ipoin)=rcsize
    if(rsize(ipoin)<c00)then
       write(*,*)'Error in gtsiz2, negative size at ipoin:',ipoin
       stop
    endif

    !rsize(ipoin)=1.0d-04

  end subroutine gtsiz2

  subroutine carout(ncell,lcell,lmark)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)    :: ncell
    integer(ip),intent(inout) :: lmark(ncell)
    type(cell)                :: lcell(ncell)
    integer(ip)               :: icell,j,istack,ineigh,nstack,icont,ncell0,ncont
    integer(ip)               :: ilevel,jlevel,icont1,icont2
    integer(ip),pointer       :: lstack(:)
    integer(4)                :: istat

    allocate(lstack(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','takeout',lstack)
    !
    !     Find an external cell  
    !
    do icell=1,ncell
       do j=1,6
          if(lcell(icell)%neigh(j)==0 .and. lcell(icell)%marked/=-3)goto 100  
       enddo
    enddo
    !
    !     Error, no external cell found
    ! 
    write(*,*)'No external cell found'
    do icell=1,ncell
       lmark(icell)=1_ip
    enddo
    return

100 continue
    !
    !     Now mark all the cells without crossing a cell marked -3
    ! 
    lstack(1)=icell
    lmark(icell)=1_ip
    nstack=1_ip
    istack=0_ip
    !
    !     Loop on the stack
    !
    do

       if(istack==nstack)exit
       istack=istack+1_ip

       icell=lstack(istack)
       !
       !    Loop on the neighbors
       !    
       do j=1,6
          ineigh=lcell(icell)%neigh(j) 
          if(ineigh==0)cycle
          if(lcell(ineigh)%marked==-3)cycle
          if(lmark(ineigh)==1)cycle
          !
          !     Add to the stack
          !   
          nstack=nstack+1_ip
          lstack(nstack)=ineigh
          lmark(ineigh)=1_ip

       enddo

    enddo
    !
    !     Add the -3 to the stack
    !
    do icell=1,ncell
       if(lcell(icell)%marked==-3)then
          lmark(icell)=1_ip
       endif
    enddo
    !
    !     DBG
    !
    !icont1=0_ip
    !icont2=0_ip
    !do icell=1,ncell
    !   if(lmark(icell)==0)then
    !      write(*,*)'icell=',icell
    !      icont1=icont1+1
    !      else
    !      icont2=icont2+1
    !  endif
    !enddo
    !write(*,*)'icont1=',icont1
    !write(*,*)'icont2=',icont2
    !
    !     For the hanging nodes, force the small deleted cells
    !     beside large one to be kept
    ! 
    !do icell=1,ncell
    !   if(lmark(icell)==0)then
    !      ilevel=lcell(icell)%level 
    !      do j=1,6
    !         ineigh=lcell(icell)%neigh(j)
    !         if(ineigh==0)cycle
    !         jlevel=lcell(ineigh)%level 
    !
    !     Only consider larger cells
    !
    !         if(jlevel<ilevel)then    
    !            if(lmark(ineigh)==1)then
    !               lmark(icell)=1_ip
    !               exit
    !            endif  
    !         endif
    !      enddo
    !   endif
    !enddo

    !
    !     Put the neighbors to zero if not marked
    !
    !do icell=1,ncell
    !   if(lcell(icell)%marked==-3)then
    !      do j=1,6
    !         ineigh=lcell(icell)%neigh(j)
    !         if(ineigh==0)cycle
    !         if(lmark(ineigh)==0)then
    !            lcell(icell)%neigh(j)=0_ip
    !         endif
    !      enddo
    !   endif
    !enddo
    !
    !     Renumber the cells
    !
    !icont=0_ip
    !do icell=1,ncell
    !   if(lmark(icell)==1)then
    !      icont=icont+1_ip
    !      lmark(icell)=icont
    !   endif
    !enddo
    !
    !     Compact the lcell array
    !
    !ncont=0_ip
    !ncell0=ncell
    !do icell=1,ncell0
    !   if(lmark(icell)>0)then
    !      ncont=ncont+1
    !      lcell(ncont)%marked    = lcell(icell)%marked
    !      lcell(ncont)%level     = lcell(icell)%level
    !      lcell(ncont)%coor(1,1) = lcell(icell)%coor(1,1)
    !      lcell(ncont)%coor(2,1) = lcell(icell)%coor(2,1)
    !      lcell(ncont)%coor(3,1) = lcell(icell)%coor(3,1)
    !      lcell(ncont)%coor(1,2) = lcell(icell)%coor(1,2)
    !      lcell(ncont)%coor(2,2) = lcell(icell)%coor(2,2)
    !      lcell(ncont)%coor(3,2) = lcell(icell)%coor(3,2)
    !      lcell(ncont)%rsize     = lcell(icell)%rsize
    !
    !     Renumber the neighbors
    !
    !      ineigh=lcell(icell)%neigh(1)
    !      if(ineigh==0)then
    !         lcell(ncont)%neigh(1)=0 
    !      else 
    !         lcell(ncont)%neigh(1)=lmark(ineigh) 
    !      endif

    !      ineigh=lcell(icell)%neigh(2)
    !      if(ineigh==0)then
    !         lcell(ncont)%neigh(2)=0 
    !      else 
    !         lcell(ncont)%neigh(2)=lmark(ineigh)
    !      endif

    !      ineigh=lcell(icell)%neigh(3)
    !      if(ineigh==0)then
    !         lcell(ncont)%neigh(3)=0 
    !      else 
    !         lcell(ncont)%neigh(3)=lmark(ineigh)
    !      endif

    !      ineigh=lcell(icell)%neigh(4)
    !      if(ineigh==0)then
    !         lcell(ncont)%neigh(4)=0 
    !      else 
    !         lcell(ncont)%neigh(4)=lmark(ineigh)
    !      endif

    !      ineigh=lcell(icell)%neigh(5)
    !      if(ineigh==0)then
    !         lcell(ncont)%neigh(5)=0 
    !      else 
    !         lcell(ncont)%neigh(5)=lmark(ineigh)
    !      endif

    !      ineigh=lcell(icell)%neigh(6)
    !      if(ineigh==0)then
    !         lcell(ncont)%neigh(6)=0 
    !      else 
    !         lcell(ncont)%neigh(6)=lmark(ineigh)
    !      endif
    !   endif
    !enddo

    !ncell=ncont     


    !call memchk(2_ip,istat,memor_msh,'LMARK','takeout',lmark)
    !deallocate(lmark,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LMARK','takeout',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','takeout',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','takeout',0_ip)

  end subroutine carout

  subroutine gtcartsize(lcell,icell,ncell)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in) :: ncell,icell
    type(cell)             :: lcell(ncell)
    real(rp)               :: dx,dy,dz,rsize,c14

    c14=1.0d+00/4.0d+00 

    dx=lcell(icell)%coor(1,2)-lcell(icell)%coor(1,1)
    dy=lcell(icell)%coor(2,2)-lcell(icell)%coor(2,1)
    dz=lcell(icell)%coor(3,2)-lcell(icell)%coor(3,1)
    !rsize=(dx+dy+dz)*c14 
    !rsize=max(dx,dy,dz)  more stringent if non isotropic
    rsize=min(dx,dy,dz)
    !rsize=min(rsize,rsuni)
    lcell(icell)%rsize=rsize

  end subroutine gtcartsize

  subroutine scalsize(ncell,lcell,rscal)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(inout) :: ncell
    type(cell)             :: lcell(ncell)
    real(rp),intent(in)       :: rscal
    real(rp)                  :: c00
    integer(ip)               :: icell
    !
    !     Check input
    !
    c00=0.0d+00
    if(rscal<c00)then
       write(*,*)'Error in scalsize, rcoef negative'
       stop
    endif

    do icell=1,ncell
       lcell(icell)%rsize=lcell(icell)%rsize*rscal
    enddo

  end subroutine scalsize

  subroutine smoothsize(lcell,ncell,niter,limp)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in) :: ncell,niter
    integer(ip),intent(in) :: limp(ncell)
    type(cell)             :: lcell(ncell)
    real(rp)               :: c00,rsizenew,dsize,epsil,drange,c23,c43,c10
    integer(ip)            :: icell,iter,j,ineigh,imax,ilevel
    real(rp)               :: rcsize,dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3 
    real(rp)               :: cx1,cx2,cx3,cy1,cy2,cy3,cz1,cz2,cz3,ct,c05 
    !
    !     This sub smoothes out the size 
    !
    c00=0.0d+00
    c10=1.0d+00
    epsil=1.0d-12
    c23=2.0d+00/3.0d+00
    c43=4.0d+00/3.0d+00
    c05=0.5d+00

    do iter=1,niter

       do icell=ncell,1,-1
          !
          !     Do not modify the imposed cells
          !
          if(limp(icell)==1)cycle
          !
          !     Get the cell level, size, and size distribution
          !
          ilevel=lcell(icell)%level 
          rcsize=lcell(icell)%rsize
          !
          !     Get the neighbors
          !
          ineigh=lcell(icell)%neigh(1)
          if(ineigh/=0)then 
             if(lcell(ineigh)%level<ilevel)then
                dx1=(rcsize-lcell(ineigh)%rsize)*c23
                cx1=c23
             else if(lcell(ineigh)%level>ilevel)then
                dx1=(rcsize-lcell(ineigh)%rsize)*c43
                cx1=c43
             else 
                dx1=rcsize-lcell(ineigh)%rsize
                cx1=c10
             endif
          else
             dx1=c00
             cx1=c00 
          endif

          ineigh=lcell(icell)%neigh(2) 
          if(ineigh/=0)then 
             if(lcell(ineigh)%level<ilevel)then
                dx2=(lcell(ineigh)%rsize-rcsize)*c23
                cx2=c23
             else if(lcell(ineigh)%level>ilevel)then
                dx2=(lcell(ineigh)%rsize-rcsize)*c43
                cx2=c43
             else 
                dx2=lcell(ineigh)%rsize-rcsize
                cx2=c10
             endif
          else
             dx2=c00
             cx2=c00 
          endif

          ineigh=lcell(icell)%neigh(3)
          if(ineigh/=0)then 
             if(lcell(ineigh)%level<ilevel)then
                dy1=(rcsize-lcell(ineigh)%rsize)*c23
                cy1=c23
             else if(lcell(ineigh)%level>ilevel)then
                dy1=(rcsize-lcell(ineigh)%rsize)*c43
                cy1=c43
             else 
                dy1=rcsize-lcell(ineigh)%rsize
                cy1=c10
             endif
          else
             dy1=c00
             cy1=c00 
          endif

          ineigh=lcell(icell)%neigh(4) 
          if(ineigh/=0)then 
             if(lcell(ineigh)%level<ilevel)then
                dy2=(lcell(ineigh)%rsize-rcsize)*c23
                cy2=c23
             else if(lcell(ineigh)%level>ilevel)then
                dy2=(lcell(ineigh)%rsize-rcsize)*c43
                cy2=c43
             else 
                dy2=lcell(ineigh)%rsize-rcsize
                cy2=c10
             endif
          else
             dy2=c00
             cy2=c00 
          endif

          ineigh=lcell(icell)%neigh(5)
          if(ineigh/=0)then 
             if(lcell(ineigh)%level<ilevel)then
                dz1=(rcsize-lcell(ineigh)%rsize)*c23
                cz1=c23
             else if(lcell(ineigh)%level>ilevel)then
                dz1=(rcsize-lcell(ineigh)%rsize)*c43
                cz1=c43
             else 
                dz1=rcsize-lcell(ineigh)%rsize
                cz1=c10
             endif
          else
             dz1=c00
             cz1=c00 
          endif

          ineigh=lcell(icell)%neigh(6) 
          if(ineigh/=0)then 
             if(lcell(ineigh)%level<ilevel)then
                dz2=(lcell(ineigh)%rsize-rcsize)*c23
                cz2=c23
             else if(lcell(ineigh)%level>ilevel)then
                dz2=(lcell(ineigh)%rsize-rcsize)*c43
                cz2=c43
             else 
                dz2=lcell(ineigh)%rsize-rcsize
                cz2=c10
             endif
          else
             dz2=c00
             cz2=c00 
          endif
          !
          !     Compute constants for the Laplacian
          !
          ct=cx1+cx2+cy1+cy2+cz1+cz2
          !
          !     Compute Laplacian of the size   
          ! 
          dsize=((dx2-dx1)+(dy2-dy1)+(dz2-dz1))/ct
          !
          !     Get the final smoothed size 
          !
          rsizenew=rcsize+dsize
          if(rsizenew<c00)then
             write(*,*)'Error in smoothsize, negative size'
             stop
          endif
          !
          !     Final size
          ! 
          call gtcartsize(lcell,icell,ncell)
          lcell(icell)%rsize=max(min(rsizenew,rcsize),lcell(icell)%rsize*c05)
          !lcell(icell)%rsize=rsizenew
       enddo

    enddo

  end subroutine smoothsize

  subroutine unisizmin(ncell,lcell,rsuni)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in) :: ncell
    real(rp),intent(in)    :: rsuni
    type(cell)             :: lcell(ncell)
    integer(ip)            :: icell 

    do icell=1,ncell
       lcell(icell)%rsize=min(lcell(icell)%rsize,rsuni)
    enddo

  end subroutine unisizmin

  subroutine unisiz(ncell,lcell,rsuni)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in) :: ncell
    real(rp),intent(in)    :: rsuni
    type(cell)             :: lcell(ncell)
    integer(ip)            :: icell 

    do icell=1,ncell
       lcell(icell)%rsize=rsuni
    enddo

  end subroutine unisiz

  subroutine soursiz(ncell,lcell,ndim,nsour,rsour,rsgeo,rtol,rsuni)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(4)                :: istat
    integer(ip),intent(in)        :: ndim,ncell,nsour
    real(rp),intent(in)           :: rsgeo(ndim,ndim,nsour),rtol,rsuni
    real(rp),pointer              :: rsour(:,:)
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: icell,isour,l,neigh,icart,nstack,istack
    real(rp)                      :: xc,yc,zc,c10,csca,cscal,rpoin(ndim),rcellsize
    real(rp)                      :: rvec1x,rvec1y,rvec1z,rvec2x,rvec2y,rvec2z,cbig
    real(rp)                      :: rnx,rny,rnz,rface(ndim),rnl,p1(ndim),p2(ndim)
    real(rp)                      :: dtot,d1,d2,d3,rdist,epsil,epsil1,epsil0,c13,c05,rsize
    real(rp)                      :: rdirx,rdiry,rdirz,xpro,ypro,zpro,rl,rdist1,rdist2,rdist3
    integer(ip),pointer           :: lmark(:),lstack(:)

    c10=1.0d+00
    c05=0.5d+00
    c13=1.0d+00/3.0d+00
    epsil=1.0d-06
    epsil1=1.0d+00+1.0d-04
    epsil0=1.0d-04
    cbig=1.0d+12*rsuni 

    !allocate(lmark(ncell),stat=istat)
    !call memchk(zero,istat,memor_msh,'LMARK','soursiz',lmark)
    !allocate(lstack(ncell),stat=istat)
    !call memchk(zero,istat,memor_msh,'LSTACK','soursiz',lstack)
    !
    !     Initialize the cell size to a huge value
    !
    !do icell=1,ncell
    !   lcell(icell)%rsize=cbig
    !enddo
    !
    !     Loop on the sources
    !
    do isour=1,nsour
       !
       !     Compute the center of the triangle
       !
       rpoin(1)=(rsgeo(1,1,isour)+rsgeo(1,2,isour)+rsgeo(1,3,isour))*c13
       rpoin(2)=(rsgeo(2,1,isour)+rsgeo(2,2,isour)+rsgeo(2,3,isour))*c13
       rpoin(3)=(rsgeo(3,1,isour)+rsgeo(3,2,isour)+rsgeo(3,3,isour))*c13
       !
       !     Compute source side vectors
       !
       rvec1x=rsgeo(1,2,isour)-rsgeo(1,1,isour)
       rvec1y=rsgeo(2,2,isour)-rsgeo(2,1,isour)
       rvec1z=rsgeo(3,2,isour)-rsgeo(3,1,isour)

       rvec2x=rsgeo(1,3,isour)-rsgeo(1,1,isour)
       rvec2y=rsgeo(2,3,isour)-rsgeo(2,1,isour)
       rvec2z=rsgeo(3,3,isour)-rsgeo(3,1,isour)
       !
       !     Compute normal to source
       !
       rnx= rvec1y*rvec2z-rvec1z*rvec2y 
       rny=-rvec1x*rvec2z+rvec1z*rvec2x 
       rnz= rvec1x*rvec2y-rvec1y*rvec2x 

       rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
       !
       !     Remember face area (*2)
       !
       dtot=rnl
       rnl=c10/rnl

       rface(1)=rnl*rnx
       rface(2)=rnl*rny
       rface(3)=rnl*rnz
       !
       !     Find the cell containing this point
       !
       !icart=1_ip
       !call gtelem(1,rpoin,1,ndim,lcell,ncell,icart,rtol)
       !if(icart==0)then
       !   write(*,*)'Error in soursiz, point out of domain'
       !   stop 
       !endif
       !
       !     Initialize the stack with icart
       !
       !lstack(1)=icart
       !lmark(icart)=isour
       !nstack=1_ip
       !istack=0_ip
       !
       !     Propagate the source size
       !
       !do  
       do icell=1,ncell

          !if(istack==nstack)exit
          !istack=istack+1 
          !
          !     Get the cell number
          !  
          !icell=lstack(istack)
          !
          !     Get cell center
          !
          xc=(lcell(icell)%coor(1,1)+lcell(icell)%coor(1,2))*c05
          yc=(lcell(icell)%coor(2,1)+lcell(icell)%coor(2,2))*c05
          zc=(lcell(icell)%coor(3,1)+lcell(icell)%coor(3,2))*c05

          rcellsize=(lcell(icell)%coor(1,2)-xc)+(lcell(icell)%coor(2,2)-yc)+(lcell(icell)%coor(3,2)-zc)
          !
          !     Project the point on the plane of the source
          !
          rdirx=xc-rsgeo(1,1,isour)
          rdiry=yc-rsgeo(2,1,isour)
          rdirz=zc-rsgeo(3,1,isour)

          csca=rface(1)*rdirx+rface(2)*rdiry+rface(3)*rdirz

          xpro=xc-csca*rface(1)
          ypro=yc-csca*rface(2)
          zpro=zc-csca*rface(3)
          !
          !     Side (ip2,ip3)
          !
          p1(1)=rsgeo(1,2,isour)-xpro
          p1(2)=rsgeo(2,2,isour)-ypro
          p1(3)=rsgeo(3,2,isour)-zpro
          p2(1)=rsgeo(1,3,isour)-xpro
          p2(2)=rsgeo(2,3,isour)-ypro
          p2(3)=rsgeo(3,3,isour)-zpro
          call orient3D(p1,p2,rface,d1,ndim)
          d1=d1/dtot
          !
          !     Side (ip3,ip1)
          !
          p1(1)=rsgeo(1,3,isour)-xpro
          p1(2)=rsgeo(2,3,isour)-ypro
          p1(3)=rsgeo(3,3,isour)-zpro
          p2(1)=rsgeo(1,1,isour)-xpro
          p2(2)=rsgeo(2,1,isour)-ypro
          p2(3)=rsgeo(3,1,isour)-zpro
          call orient3D(p1,p2,rface,d2,ndim)
          d2=d2/dtot
          !
          !     Side (ip1,ip2)
          !
          p1(1)=rsgeo(1,1,isour)-xpro
          p1(2)=rsgeo(2,1,isour)-ypro
          p1(3)=rsgeo(3,1,isour)-zpro
          p2(1)=rsgeo(1,2,isour)-xpro
          p2(2)=rsgeo(2,2,isour)-ypro
          p2(3)=rsgeo(3,2,isour)-zpro
          call orient3D(p1,p2,rface,d3,ndim)
          d3=d3/dtot
          !
          !     Is the projected point inside the source triangle?
          !
          if(d1>-epsil .and. d2>-epsil .and. d3>-epsil)then  
             rdist=abs(csca)
          else
             !
             !     Evaluate distance to lines
             !

             !
             !     Side (ip2,ip3)
             !
             p1(1)=rsgeo(1,3,isour)-rsgeo(1,2,isour)
             p1(2)=rsgeo(2,3,isour)-rsgeo(2,2,isour)
             p1(3)=rsgeo(3,3,isour)-rsgeo(3,2,isour)
             rl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
             rl=c10/rl
             p1(1)=rl*p1(1) 
             p1(2)=rl*p1(2) 
             p1(3)=rl*p1(3) 

             p2(1)=xc-rsgeo(1,2,isour)
             p2(2)=yc-rsgeo(2,2,isour)
             p2(3)=zc-rsgeo(3,2,isour)

             csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
             cscal=csca*rl 

             if(cscal<-epsil0)then

                rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

             else if(cscal>epsil1)then

                p2(1)=xc-rsgeo(1,3,isour)
                p2(2)=yc-rsgeo(2,3,isour)
                p2(3)=zc-rsgeo(3,3,isour)
                rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

             else

                p2(1)=rsgeo(1,2,isour)+csca*p1(1)-xc
                p2(2)=rsgeo(2,2,isour)+csca*p1(2)-yc
                p2(3)=rsgeo(3,2,isour)+csca*p1(3)-zc
                rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 

             endif
             !
             !     Side (ip3,ip1)
             !
             p1(1)=rsgeo(1,1,isour)-rsgeo(1,3,isour)
             p1(2)=rsgeo(2,1,isour)-rsgeo(2,3,isour)
             p1(3)=rsgeo(3,1,isour)-rsgeo(3,3,isour)
             rl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
             rl=c10/rl
             p1(1)=rl*p1(1) 
             p1(2)=rl*p1(2) 
             p1(3)=rl*p1(3) 

             p2(1)=xc-rsgeo(1,3,isour)
             p2(2)=yc-rsgeo(2,3,isour)
             p2(3)=zc-rsgeo(3,3,isour)

             csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
             cscal=csca*rl 

             if(cscal<-epsil0)then

                rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

             else if(cscal>epsil1)then

                p2(1)=xc-rsgeo(1,1,isour)
                p2(2)=yc-rsgeo(2,1,isour)
                p2(3)=zc-rsgeo(3,1,isour)
                rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

             else

                p2(1)=rsgeo(1,3,isour)+csca*p1(1)-xc
                p2(2)=rsgeo(2,3,isour)+csca*p1(2)-yc
                p2(3)=rsgeo(3,3,isour)+csca*p1(3)-zc
                rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 

             endif
             !
             !     Side (ip1,ip2)
             !
             p1(1)=rsgeo(1,2,isour)-rsgeo(1,1,isour)
             p1(2)=rsgeo(2,2,isour)-rsgeo(2,1,isour)
             p1(3)=rsgeo(3,2,isour)-rsgeo(3,1,isour)
             rl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
             rl=c10/rl
             p1(1)=rl*p1(1) 
             p1(2)=rl*p1(2) 
             p1(3)=rl*p1(3) 

             p2(1)=xc-rsgeo(1,1,isour)
             p2(2)=yc-rsgeo(2,1,isour)
             p2(3)=zc-rsgeo(3,1,isour)

             csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
             cscal=csca*rl 

             if(cscal<-epsil0)then

                rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

             else if(cscal>epsil1)then

                p2(1)=xc-rsgeo(1,2,isour)
                p2(2)=yc-rsgeo(2,2,isour)
                p2(3)=zc-rsgeo(3,2,isour)
                rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  

             else

                p2(1)=rsgeo(1,1,isour)+csca*p1(1)-xc
                p2(2)=rsgeo(2,1,isour)+csca*p1(2)-yc
                p2(3)=rsgeo(3,1,isour)+csca*p1(3)-zc
                rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 

             endif
             !
             !     Take the min over the side distances
             ! 
             rdist=min(rdist1,rdist2,rdist3)

          endif
          !
          !     The distance has been computed, now compute the desired size
          !
          !
          !     Are we inside the sphere of constant size
          !
          if(rdist<rsour(2,isour))then
             rsize=rsour(1,isour)
          else
             rsize=rsour(3,isour)*(rdist-rsour(2,isour))+rsour(1,isour)
          endif
          !
          !     Now compare to the cell size
          !
          lcell(icell)%rsize=min(rsize,lcell(icell)%rsize)
          !
          !    Add cells to the stack
          !    For the moment, add all the cells for each source...
          !
          !do l=1,6
          !   neigh=lcell(icell)%neigh(l)
          !   if(neigh==0)cycle
          !   if(lmark(neigh)/=isour)then
          !      nstack=nstack+1
          !      lstack(nstack)=neigh
          !      lmark(neigh)=isour
          !   endif
          !enddo
       enddo
    enddo

    !call memchk(2_ip,istat,memor_msh,'LMARK','soursiz',lmark)
    !deallocate(lmark,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LMARK','soursiz',0_ip)
    !call memchk(2_ip,istat,memor_msh,'LSTACK','soursiz',lstack)
    !deallocate(lstack,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LSTACK','soursiz',0_ip)

  end subroutine soursiz

  subroutine renucart(lcell,ncell)
    use def_kintyp, only       :  ip,rp,lg,cell
    use def_meshin, only          :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: ncell
    type(cell)                  :: lcell(ncell)
    integer(ip),pointer         :: lstack(:),lrenu(:)
    real(rp),pointer            :: rtemp(:) 
    integer(ip)                 :: icell,icont,j,istack,nstack,ineigh,k,icellnew
    integer(4)                  :: istat
    !
    !     Allocate lstack
    ! 
    allocate(lrenu(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','renucart',lrenu)
    allocate(lstack(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','renucart',lstack)
    !
    !     First, find a corner of the cartesian mesh
    !
    do icell=1,ncell
       icont=0_ip
       do j=1,6
          ineigh=lcell(icell)%neigh(j)
          if(ineigh==0)icont=icont+1_ip
       enddo
       !
       !     Do we have a corner?
       !
       if(icont==3)goto 10

    enddo
    !
    !     Error, we did not find a corner
    ! 
    write(*,*)'Error in renucart, no corner found'
    stop

10  continue 
    !
    !     Get the renumbering
    !
    lstack(1)=icell
    nstack=1_ip
    lrenu(icell)=1_ip 
    istack=0_ip 
    !
    !     Loop on the stack
    !
    do  
       write(*,*)istack,nstack 
       if(istack==nstack)exit
       istack=istack+1_ip
       icell=lstack(istack) 
       write(*,*)icell 
       do j=1,6
          ineigh=lcell(icell)%neigh(j)
          if(ineigh==0)cycle
          if(lrenu(ineigh)/=0)cycle
          nstack=nstack+1_ip
          lstack(nstack)=ineigh
          lrenu(ineigh)=nstack
       enddo
    enddo
    !
    !     Is everything ok?
    !
    if(nstack/=ncell)then
       do icell=1,ncell
          if(lrenu(icell)==0)then
             write(*,*)'Error in  renucart, lrenu=0 for icell:',icell
             stop
          endif
       enddo
    endif
    !
    !     Renumber the cells, first the integer part
    !
    do j=1,6
       do icell=1,ncell
          ineigh=lcell(icell)%neigh(j)
          if(ineigh==0)then
             lstack(icell)=0_ip
          else 
             lstack(icell)=lrenu(ineigh)
          endif
       enddo

       do icell=1,ncell
          lcell(icell)%neigh(j)=lstack(icell)
       enddo
    enddo

    do icell=1,ncell
       icellnew=lrenu(icell)
       lstack(icellnew)=lcell(icell)%level
    enddo

    do icell=1,ncell
       lcell(icell)%level=lstack(icell)
    enddo

    do icell=1,ncell
       icellnew=lrenu(icell)
       lstack(icellnew)=lcell(icell)%marked
    enddo

    do icell=1,ncell
       lcell(icell)%marked=lstack(icell)
    enddo
    !
    !     Deallocate lstack
    !
    call memchk(2_ip,istat,memor_msh,'LSTACK','renucart',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','renucart',0_ip)
    !
    !     Allocate rtemp
    !
    allocate(rtemp(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'RTEMP','renucart',rtemp)

    do icell=1,ncell
       icellnew=lrenu(icell)
       rtemp(icellnew)=lcell(icell)%rsize
    enddo

    do icell=1,ncell
       lcell(icell)%rsize=rtemp(icell)
    enddo

    do j=1,3
       do k=1,2

          do icell=1,ncell
             icellnew=lrenu(icell)
             rtemp(icellnew)=lcell(icell)%coor(j,k)
          enddo

          do icell=1,ncell
             lcell(icell)%coor(j,k)=rtemp(icell)
          enddo

       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'RTEMP','renucart',rtemp)
    deallocate(rtemp,stat=istat)
    if(istat/=0) call memerr(2_ip,'RTEMP','renucart',0_ip)
    call memchk(2_ip,istat,memor_msh,'LRENU','renucart',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','renucart',0_ip)


  end subroutine renucart


  subroutine findMin(ncell,lcell,limp)
    use def_kintyp, only       :  ip,rp,lg,cell
    use def_meshin, only          :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: ncell
    type(cell)                  :: lcell(ncell)
    integer(ip),intent(inout)   :: limp(ncell)
    integer(ip),pointer         :: lmark(:),lstack(:)
    integer(ip)                 :: icell,j,istack,nstack,ineigh,ilevel,jlevel,istackold,i
    integer(4)                  :: istat

    !
    !     Allocate lmark
    !
    allocate(lmark(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','renucart',lmark)
    allocate(lstack(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','renucart',lstack)

    !
    !     We must find the isolated minima.
    !     Begin from the highest level and mark all the superior levels
    !
    nstack=0_ip

    do i=ncell,1,-1
       if(lmark(i)==0)then
          !
          !     We have found a minimum
          !     Remember it in limp 
          !
          !write(*,*)'Seed is ',i
          !write(*,*)'nstack is ',nstack
          limp(i)=1_ip
          !
          !     Put it in the stack
          !
          nstack=nstack+1_ip
          lstack(nstack)=i
          lmark(i)=1_ip
          !
          !     Mark the neighbors of the same level
          !
          istack=nstack-1_ip
          istackold=istack
          do
             !write(*,*)istack,nstack 
             if(istack==nstack)exit
             istack=istack+1_ip
             icell=lstack(istack)
             ilevel=lcell(icell)%level
             do j=1,6 
                ineigh=lcell(icell)%neigh(j)
                if(ineigh==0)cycle
                if(lmark(ineigh)/=0)cycle 
                jlevel=lcell(ineigh)%level          
                if(jlevel/=ilevel)cycle
                nstack=nstack+1_ip
                lstack(nstack)=ineigh
                lmark(ineigh)=1_ip
                limp(ineigh)=1_ip
             enddo
          enddo
          !
          !     Mark the neighbors of superior levels
          !
          !write(*,*)'nstack=',nstack
          istack=istackold
          do
             !write(*,*)istack,nstack
             if(istack==nstack)exit
             istack=istack+1_ip
             icell=lstack(istack)
             ilevel=lcell(icell)%level
             do j=1,6 
                ineigh=lcell(icell)%neigh(j)
                if(ineigh==0)cycle
                if(lmark(ineigh)/=0)cycle 
                jlevel=lcell(ineigh)%level          
                if(jlevel>ilevel)cycle
                nstack=nstack+1_ip
                lstack(nstack)=ineigh
                lmark(ineigh)=1_ip
             enddo
          enddo
       endif
    enddo

    !
    !     Find the external cells
    !

    do icell=1,ncell
       do j=1,6
          ineigh=lcell(icell)%neigh(j)
          if(ineigh==0)then
             limp(icell)=1_ip 
             exit
          endif
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LSTACK','renucart',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','renucart',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','renucart',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','renucart',0_ip)


  end subroutine findMin

  subroutine boxcart(ndim,ncell,nvoxx,nvoxy,nvoxz,bboxbin,bbox,dx,dy,dz,rsuni)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)      :: ndim
    integer(ip),intent(inout)   :: ncell,nvoxx,nvoxy,nvoxz
    real(rp),intent(in)         :: bboxbin(ndim,2),rsuni
    real(rp),intent(inout)      :: bbox(ndim,2),dx,dy,dz
    real(rp)                    :: c(ndim),rx,ry,rz,c05

    c05=0.500034567d+00  
    !
    !     Take dx=dy=dz=rsuni
    !
    dx=rsuni
    dy=rsuni
    dz=rsuni
    !
    !     Compute center of the scene
    ! 
    c(1)=(bboxbin(1,1)+bboxbin(1,2))*c05
    c(2)=(bboxbin(2,1)+bboxbin(2,2))*c05
    c(3)=(bboxbin(3,1)+bboxbin(3,2))*c05
    !
    !     Compute number of cells in each direction
    !
    nvoxx=int((bboxbin(1,2)-bboxbin(1,1))/dx)+1
    nvoxy=int((bboxbin(2,2)-bboxbin(2,1))/dy)+1
    nvoxz=int((bboxbin(3,2)-bboxbin(3,1))/dz)+1
    ncell=nvoxx*nvoxy*nvoxz
    !
    !     Recompute extent
    !
    rx=real(nvoxx)*dx*c05 
    ry=real(nvoxy)*dy*c05 
    rz=real(nvoxz)*dz*c05 
    !
    !     Center the scene
    !
    bbox(1,1)=c(1)-rx 
    bbox(2,1)=c(2)-ry 
    bbox(3,1)=c(3)-rz 

    bbox(1,2)=c(1)+rx 
    bbox(2,2)=c(2)+ry 
    bbox(3,2)=c(3)+rz 

  end subroutine boxcart

  subroutine outbou(ncell,nnofa,npoin,nface,lcell,lmark,lface,renum,renum2)
    use def_kintyp, only          :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          :  memor_msh
    implicit none
    integer(ip),intent(in)    :: ncell,nnofa,npoin
    type(cell)                :: lcell(ncell)
    integer(ip),intent(in)    :: lmark(ncell),renum(8*ncell),renum2(npoin)
    integer(ip),intent(inout) :: nface
    integer(ip),pointer       :: lface(:,:)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:)
    integer(ip)               :: iface,inofa,icell,j,ineigh,tab(4,6)
    integer(ip)               :: ifa,jface,ip1,ichk,tab7(6)
    integer(4)                :: istat

    tab7(1)=1
    tab7(2)=2
    tab7(3)=3
    tab7(4)=4
    tab7(5)=5
    tab7(6)=6


    tab(1,1)=1
    tab(2,1)=5
    tab(3,1)=8
    tab(4,1)=4
    tab(1,2)=2
    tab(2,2)=3
    tab(3,2)=7
    tab(4,2)=6
    tab(1,3)=1
    tab(2,3)=2
    tab(3,3)=6
    tab(4,3)=5
    tab(1,4)=3
    tab(2,4)=4
    tab(3,4)=8
    tab(4,4)=7
    tab(1,5)=4
    tab(2,5)=3
    tab(3,5)=2
    tab(4,5)=1
    tab(1,6)=5
    tab(2,6)=6
    tab(3,6)=7
    tab(4,6)=8

    !
    !     Check that no outer face has been deleted
    ! 
    do iface=1,nface
       do inofa=1,nnofa
          ip1=lface(inofa,iface)
          if(renum2(ip1)==0)then
             write(*,*)'Outer face destroyed, must stop'
             stop
          endif
       enddo
    enddo
    !
    !     Get the face surrounding points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get the new faces
    !
    do icell=1,ncell
       do j=1,6
          ineigh=lcell(icell)%neigh(j)
          if(ineigh==0)then
             ip1=renum(8*(icell-1)+tab(1,j))
             !
             !     Find this face
             !        
             ichk=0_ip
             do ifa=ptoel2(ip1),ptoel2(ip1+1)-1
                jface=ptoel1(ifa)
                if(jface==iface)then
                   goto 100
                endif
             enddo

             nface=nface+1
             call memrea(nface,memor_msh,'LFACE','addfac',lface) 
             lface(1,nface)=renum(8*(icell-1)+tab(1,j))
             lface(2,nface)=renum(8*(icell-1)+tab(2,j))
             lface(3,nface)=renum(8*(icell-1)+tab(3,j))
             lface(4,nface)=renum(8*(icell-1)+tab(4,j))
             lface(5,nface)=tab(1,j)
             lface(6,nface)=tab(2,j)
             lface(7,nface)=tab(3,j)
             lface(8,nface)=tab(4,j)
             lface(9,nface)=icell
             lface(10,nface)=tab7(j)

100          continue
          endif
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'PTOEL1','addfac',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','addfac',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','addfac',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','addfac',0_ip)

  end subroutine outbou


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !     This section deals with cartesian meshes of points
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cartp(npoin,ndim,elem,nnode,nelem,npmax,lcart,ptoel1,ptoel2,&
       coor,lcell,ncell,rtol,ifbox,bboxin)
    use def_meshin, only       :  memor_msh
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    implicit none
    integer(ip),intent(in)    ::  ndim,npoin,nnode,nelem,npmax,ifbox
    integer(ip),intent(in)    ::  elem(nnode,nelem)
    integer(ip),intent(inout) ::  lcart(npoin),ncell
    integer(ip),intent(in)    ::  ptoel1(*),ptoel2(npoin+1)
    real(rp),intent(in)       ::  coor(ndim,npoin)
    real(rp),intent(inout)    ::  rtol
    type(cell),pointer        ::  lcell(:)
    logical(lg),pointer       ::  lpoin(:)
    real(rp)                  ::  bboxbin(ndim,2),bboxin(ndim,2)
    real(rp)                  ::  dx,dy,dz,bbox(ndim,2)
    integer(ip)               ::  nvoxx,nvoxy,nvoxz,niter
    integer(ip)               ::  ielem,inode,nploc,ipoin
    integer(4)                ::  istat
    !
    !     The %marked is:            -2 at the creation of the cell
    !                                -1 if do not intersect face
    !                                -3 if intersect but agrees with the size
    !                                 1 if marked  
    allocate(lpoin(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOIN','cartp',lpoin)
    !
    !     Mark the points belonging to elem
    !
    do ielem=1,nelem
       do inode=1,nnode
          lpoin(elem(inode,ielem))=.true.
       enddo
    enddo
    !
    !     Count how many effective points we have 
    ! 
    nploc=0_ip
    do ipoin=1,npoin
       if(lpoin(ipoin).eqv. .true.)nploc=nploc+1_ip
    enddo
    !
    !     Compute bbox of the cloud of points
    !
    call boxbinp(coor,ndim,npoin,bboxbin,lpoin)
    !
    !     Compute dx,dy,dz from Bbox, the bounding box of the cartesian mesh 
    !
    if(ifbox==0) then
       call boxcartp(ndim,ncell,nvoxx,nvoxy,nvoxz,bboxbin,bbox,dx,dy,dz,nploc)
    else
       call boxcartp(ndim,ncell,nvoxx,nvoxy,nvoxz,bboxin, bbox,dx,dy,dz,nploc)
    end if
    !
    !     Get tolerance
    !
    rtol=min(dx,dy,dz)*1.0d-08
    rtol=0.0d+00
    !
    !     Build the initial mesh
    !
    call buildIni(nvoxx,nvoxy,nvoxz,bbox,dx,dy,dz,ndim,lcell,ncell)
    !
    !     Subdivide the mesh
    ! 
    call meshp(coor,npoin,ndim,niter,nelem,nnode,elem,ptoel1,ptoel2,npmax,&
         lcart,lcell,ncell,rtol,lpoin,nploc)
    !
    !     Renumber the cells
    !
    !call renucart(lcell,ncell)

    !  write(*,*)'Cartesian mesh done, writing output'

    call writeGidconfp(ncell,lcell,coor,npoin,ndim)

    call memchk(2_ip,istat,memor_msh,'LPOIN','cartp',lpoin)
    deallocate(lpoin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPOIN','cartp',0_ip)

  end subroutine cartp

  subroutine boxbinp(coor,ndim,npoin,bboxbin,lpoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in) ::  ndim,npoin
    logical(lg),intent(in) ::  lpoin(npoin)
    real(rp),intent(in)    ::  coor(ndim,npoin)
    real(rp),intent(inout) ::  bboxbin(ndim,2)
    real(rp)               ::  rmin(ndim),rmax(ndim),dif(ndim),c(ndim)
    integer(ip)            ::  i,j
    real(rp)               ::   c05,coef1
    !
    !   This sub computes the bbox of the triangulation
    !

    c05=0.5d+00
    coef1=1.07 
    !
    !     Get first point
    !
    do j=1,npoin
       if(lpoin(j).eqv. .true.)then
          exit
       endif
    enddo

    rmin(1)=coor(1,j)
    rmin(2)=coor(2,j)
    rmin(3)=coor(3,j)

    rmax(1)=coor(1,j)
    rmax(2)=coor(2,j)
    rmax(3)=coor(3,j)

    do i=j+1,npoin

       if(lpoin(i).eqv. .true.)then

          if(coor(1,i)<rmin(1))rmin(1)=coor(1,i)
          if(coor(2,i)<rmin(2))rmin(2)=coor(2,i)
          if(coor(3,i)<rmin(3))rmin(3)=coor(3,i)

          if(coor(1,i)>rmax(1))rmax(1)=coor(1,i)
          if(coor(2,i)>rmax(2))rmax(2)=coor(2,i)
          if(coor(3,i)>rmax(3))rmax(3)=coor(3,i)

       endif
    enddo
    !
    !      Make it a litter bigger
    !
    c(1)=(rmin(1)+rmax(1))*c05
    c(2)=(rmin(2)+rmax(2))*c05
    c(3)=(rmin(3)+rmax(3))*c05

    dif(1)=rmax(1)-c(1) 
    dif(2)=rmax(2)-c(2) 
    dif(3)=rmax(3)-c(3) 

    rmin(1)=c(1)-coef1*dif(1)
    rmin(2)=c(2)-coef1*dif(2)
    rmin(3)=c(3)-coef1*dif(3)

    rmax(1)=c(1)+coef1*dif(1)
    rmax(2)=c(2)+coef1*dif(2)
    rmax(3)=c(3)+coef1*dif(3)

    bboxbin(1,1)=rmin(1)
    bboxbin(2,1)=rmin(2)
    bboxbin(3,1)=rmin(3)

    bboxbin(1,2)=rmax(1)
    bboxbin(2,2)=rmax(2)
    bboxbin(3,2)=rmax(3)

  end subroutine boxbinp


  subroutine boxcartp(ndim,ncell,nvoxx,nvoxy,nvoxz,bboxbin,bbox,dx,dy,dz,npoin)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)      :: ndim,npoin
    integer(ip),intent(inout)   :: ncell,nvoxx,nvoxy,nvoxz
    integer(ip)                 :: npdim
    real(rp),intent(in)         :: bboxbin(ndim,2)
    real(rp),intent(inout)      :: bbox(ndim,2),dx,dy,dz
    real(rp)                    :: c(ndim),rx,ry,rz,c05,c13,rmin,rsuni,coef1

    c05=0.500034567d+00  
    c13=1.0d+00/3.0d+00
    coef1=0.1d+00
    !
    !     Estimate the size of the first coarse mesh
    !  
    rx=bboxbin(1,2)-bboxbin(1,1)
    ry=bboxbin(2,2)-bboxbin(2,1)
    rz=bboxbin(3,2)-bboxbin(3,1)
    rmin=max(rx,ry,rz)       ! Safer if surfaces aligned with axis
    npdim=floor(npoin**c13)+10_ip
    rsuni=rmin/(npdim*coef1) 
    !
    !     Take dx=dy=dz=rsuni
    !
    dx=rsuni
    dy=rsuni
    dz=rsuni
    !
    !     Compute center of the scene
    ! 
    c(1)=(bboxbin(1,1)+bboxbin(1,2))*c05
    c(2)=(bboxbin(2,1)+bboxbin(2,2))*c05
    c(3)=(bboxbin(3,1)+bboxbin(3,2))*c05
    !
    !     Compute number of cells in each direction
    !
    nvoxx=int((bboxbin(1,2)-bboxbin(1,1))/dx)+1
    nvoxy=int((bboxbin(2,2)-bboxbin(2,1))/dy)+1
    nvoxz=int((bboxbin(3,2)-bboxbin(3,1))/dz)+1
    ncell=nvoxx*nvoxy*nvoxz
    !
    !     Recompute extent
    !
    rx=real(nvoxx)*dx*c05 
    ry=real(nvoxy)*dy*c05 
    rz=real(nvoxz)*dz*c05 
    !
    !     Center the scene
    !
    bbox(1,1)=c(1)-rx 
    bbox(2,1)=c(2)-ry 
    bbox(3,1)=c(3)-rz 

    bbox(1,2)=c(1)+rx 
    bbox(2,2)=c(2)+ry 
    bbox(3,2)=c(3)+rz 

  end subroutine boxcartp

  subroutine meshp(coor,npoin,ndim,niter,nelem,nnode,elem,ptoel1,ptoel2,npmax,lcart,&
       lcell,ncell,rtol,lpoin,nploc)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ndim,npoin,nnode,nelem,npmax,nploc
    integer(ip),intent(inout)     :: niter,ncell
    logical(lg),intent(in)        :: lpoin(npoin)
    type(cell), pointer           :: lcell(:) 
    integer(ip),intent(inout)     :: lcart(npoin)
    integer(ip),intent(in)        :: elem(nnode,nelem)
    integer(ip),intent(in)        :: ptoel1(*),ptoel2(npoin+1)
    real(rp),intent(in)           :: coor(ndim,npoin),rtol
    integer(ip)                   :: ncell0,iter,itermax,marked
    integer(ip)                   :: maxlevel,level,maxlev2,icell

    iter=1_ip
    itermax=50_ip

    maxlevel=8_ip
    maxlev2=int(maxlevel/2.0d+00)

    !Loop on refinement
    level=0_ip

    do

       level=level+1

       !Mark the lcells to be refined

       marked=0_ip
       ncell0=ncell 


       call markCellPoint(ndim,ptoel1,ptoel2,ncell,elem,nelem,nnode,coor,npoin,npmax,&
            lcart,lcell,marked,rtol,lpoin,nploc)

       write(*,*) marked,' lcells marked for refinement on',ncell0, 'total lcells for iter ',iter 
       !
       !     Do we have to do something?
       !
       if(marked==0)exit
       !
       !     DBG
       !
       do icell=1,ncell
          if(lcell(icell)%marked==1)then
             if(lcell(icell)%level/=iter)then
                write(*,*)'Error, elements marked not at the last level'
                stop          
             endif
          endif
       enddo
       !
       !     Sweep on the levels of refinement
       !
       call smooth(ncell,iter,ndim,lcell)
       ncell0=ncell 
       !
       !     debug
       !
       call dbgconform(lcell,ncell)
       !
       !     Refine the lcells 
       !
       call divideCell(ncell,ndim,lcell)
       !write(*,*)ncell
       !
       !     Refine the faces
       !
       call divideFace(ncell0,lcell)
       !
       !     debug
       !
       call dbgconform(lcell,ncell)
       !
       !     Compact arrays
       !
       call compact(ncell,lcell)
       !write(*,*)ncell
       !
       !     Debug
       !
       if(iter==itermax)then
          exit
       endif
       iter=iter+1

    enddo


    niter=iter


    !Full intersection test

    !realCut(btoel1,btoel2,nbin,nbinx,ncell0,lcell,dxbin,dybin,dzbin,bboxbin,trisize,marked,tribox,lface,*coor)

    !Cut lcells

    ! call  cutCell(btoel1,btoel2,nbin,nbinx,ncell0,lcell,dxbin,dybin,dzbin,bboxbin,trisize,marked,tribox,lface,coor,npo,mpo,nface)

  end subroutine meshp

  subroutine markCellPoint(ndim,ptoel1,ptoel2,ncell,elem,nelem,nnode,coor,npoin,npmax,&
       lcart,lcell,marked,rtol,lpoin,nploc)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)        :: ndim,npoin,ncell,nelem,nnode,npmax,nploc
    integer(ip),intent(in)        :: elem(nnode,nelem)
    logical(lg),intent(in)        :: lpoin(npoin)
    type(cell), intent(inout)     :: lcell(ncell)
    integer(ip),intent(inout)     :: lcart(npoin),marked
    real(rp),intent(in)           :: coor(ndim,npoin),rtol
    integer(ip),intent(in)        :: ptoel1(*),ptoel2(npoin+1) 
    integer(ip),pointer           :: lcmark(:)
    integer(ip)                   :: ipoin,icell 
    integer(4)                ::  istat
    !
    !     This sub mark the cells for refinement if the number of points in the cell
    !     is greater than npmax
    !
    !
    !     Allocate lcmark
    !
    allocate(lcmark(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LCMARK','markCellPoint',lcmark)
    !
    !     Clean up lcart (not very usefull)
    !
    do ipoin=1,npoin
       lcart(ipoin)=0_ip
    enddo
    !
    !     Get the elements of the cart grid
    !
    call intercartp(coor,npoin,ndim,elem,nelem,nnode,lcell,ncell,lcart,ptoel1,ptoel2,rtol,&
         lpoin,nploc)
    !
    !     Count how many points we have in each cell
    !
    do ipoin=1,npoin
       if(lpoin(ipoin).eqv. .true.)then
          icell=lcart(ipoin)
          lcmark(icell)=lcmark(icell)+1_ip
       endif
    enddo
    !
    !     If we have more than npmax points in the cell --> refine 
    !  
    do icell=1,ncell
       if(lcmark(icell)>npmax)then
          lcell(icell)%marked=1_ip
          marked=marked+1_ip
       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'LCMARK','renucart',lcmark)
    deallocate(lcmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCMARK','renucart',0_ip)

  end subroutine markCellPoint

  subroutine intercartp(coor,npoin,ndim,lface,nface,nnofa,lcell,ncell,lcart,ptoel1,ptoel2,rtol,&
       lpoin,nploc)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)        :: ncell,npoin,nface,ndim,nnofa,nploc
    integer(ip),intent(in)        :: lface(nnofa,nface),ptoel1(*),ptoel2(*)
    logical(lg),intent(in)        :: lpoin(npoin)
    real(rp),intent(in)           :: coor(ndim,npoin),rtol                   
    integer(ip),intent(inout)     :: lcart(npoin) 
    integer(ip),pointer           :: lstack(:) 
    logical(lg),pointer           :: lmark(:)
    integer(4)                    :: istat
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: nstack,ip1,ip2,j,ie,ipoin,ielem,icart
    integer(ip)                   :: istack,ifirst,jp1 
    !
    !     This sub finds the cartesian cell containing the points of lface for each point
    !     and store it in lcart. The point to element graph is given in ptoel
    !
    allocate(lstack(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'lstack','intercart',lstack)
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'lstack','intercart',lstack)

    do ipoin=1,npoin
       if(lpoin(ipoin).eqv. .true.)then
          lmark(ipoin)=.false.
       else
          lmark(ipoin)=.true.
       endif
    end do
    !
    !     Initialize ifirst
    !
    ifirst=max(ncell/2_ip,1_ip)
    !
    !     Initialize the stack
    !
    do jp1=1,npoin
       if(lmark(jp1).eqv. .false.)then
          exit
       endif
    enddo

    lstack(1)=jp1
    nstack=1_ip
    lcart(jp1)=ifirst
    lmark(jp1)=.true.
    !
    !     Main loop on stack
    !
    istack=0_ip

    do 

       if(istack==nstack)then
          !
          !     Do we have various connex components
          !
          do j=1,npoin
             if(.not.lmark(j)) then
                exit
             endif
          enddo
          !
          !     Did we find some unmarked points?
          !
          if(j==npoin+1)then
             !
             !     Did we mark all the points to be marked?
             !
             if(nstack==nploc)then
                exit 
             else
                write(*,*)'Error intercartp, new point not found'
                stop
             endif
          endif
          ipoin=j
          nstack=nstack+1
          lstack(nstack)=ipoin
          lcart(ipoin)=ifirst
          lmark(ipoin)=.true.
       endif

       istack=istack+1_ip
       ipoin=lstack(istack)

       icart=lcart(ipoin)
       call gtelem(ipoin,coor,npoin,ndim,lcell,ncell,icart,rtol)
       if(icart==0)then
          write(*,*)'Error in intercartp, point out of domain'
          stop 
       endif
       !
       !     Remember the element containing the point
       !
       lcart(ipoin)=icart
       !
       !     Add the neighbors to the stack and initialize the initial guess
       !
       do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(ie)
          do j=1,nnofa
             ip2=lface(j,ielem)
             if(.not.lmark(ip2)) then
                lmark(ip2)=.true.      
                nstack=nstack+1
                lstack(nstack)=ip2
                lcart(ip2)=icart
             endif
          enddo
       enddo
    enddo


    call memchk(2_ip,istat,memor_msh,'LMARK','intercart',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','intercart',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','intercart',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','intercart',0_ip)

  end subroutine intercartp

  subroutine ctopnt( lcart, npoin , ncell, lctop1,lctop2 )
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: npoin,ncell
    integer(ip), intent(in)            :: lcart(npoin)
    integer(ip), pointer               :: lctop1(:),lctop2(:)
    integer(ip)                        :: iplace,ipoin,nctopnt,icell
    integer(4)                 :: istat

    if(.not.associated(lctop2))then
       allocate(lctop2(ncell+1),stat=istat)
       call memchk(zero,istat,memor_msh,'LCTOP2','ctopnt',lctop2)
    else
       call memrea(ncell+1_ip,memor_msh,'LCTOP2','ctopnt',lctop2)
       lctop2=0_ip
    endif
    !
    !    Store ahead
    !
    do ipoin=1,npoin
       icell=lcart(ipoin)+1
       lctop2(icell)=lctop2(icell)+1
    enddo
    !
    !     Sum up
    !
    lctop2(1)=1
    do icell=2,ncell+1
       lctop2(icell)=lctop2(icell)+lctop2(icell-1)
    enddo
    !
    !      Allocate lctop1
    !
    nctopnt=lctop2(ncell+1)-1


    if(.not.associated(lctop1))then
       allocate(lctop1(nctopnt),stat=istat)
       call memchk(zero,istat,memor_msh,'LCTOP1','ctopnt',lctop1)
    else
       call memrea(nctopnt,memor_msh,'LCTOP1','ctopnt',lctop1)
    endif

    !
    !     Store in lctop1
    !
    do ipoin=1,npoin
       icell=lcart(ipoin)
       iplace=lctop2(icell)
       lctop1(iplace)=ipoin
       lctop2(icell)=iplace+1
    enddo
    !
    !     Clean up
    !
    do icell=ncell+1,2,-1
       lctop2(icell)=lctop2(icell-1)
    enddo

    lctop2(1)=1

  end subroutine ctopnt

  subroutine ctopnt2(lcart,npoin,ncell,lctop,npmax )
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: npoin,ncell,npmax
    integer(ip), intent(in)            :: lcart(npoin)
    integer(ip)                        :: lctop(npmax,ncell)
    integer(ip)                        :: iplace,ipoin,nctopnt,icell
    integer(4)                 :: istat
    !
    !     This subroutine assumes that the points in one cell will not be
    !     bigger than npmax-1
    !
    !
    !     Store in lctop
    !
    do ipoin=1,npoin
       icell=lcart(ipoin)
       if(icell==0)cycle
       iplace=lctop(npmax,icell)+1_ip
       if(iplace==npmax)then
          write(*,*)'Error in ctopnt2'
          stop
       endif
       lctop(iplace,icell)=ipoin
       lctop(npmax,icell)=iplace   
    enddo

  end subroutine ctopnt2

  subroutine getClose(pnew,icellinit,ndim,coor,npoin,lcell,ncell,lmark,lstack,&
       lctop1,lctop2,ipclos)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ndim,npoin,ncell,icellinit
    integer(ip),intent(in)        :: lctop1(*),lctop2(ncell+1)
    real(rp),intent(in)           :: coor(ndim,npoin),pnew(ndim)
    integer(ip),intent(inout)     :: lmark(ncell),lstack(ncell),ipclos
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: icell,jpo,jpoin,jmin,nstack,istack,i,ineigh
    real(rp)                      :: rx,ry,rz,rdistmin,vmin(ndim),vmax(ndim),rdist
    !
    !     This sub finds the points closest to pnew.
    !     Be carefull, pnew must not be in coor, otherwise it will be its own closest point
    !     For speed, the assumption is that lcart ALREADY CONTAINS the cartesian cell  
    !     containing pnew
    !

    if(icellinit<=0)then
       write(*,*)'Error in getClose, icellinit=',icellinit 
       stop
    endif
    !
    !     DBG
    ! 
    !do icell=1,ncell
    !   if(lmark(icell)/=0)then
    !      write(*,*)'Error in getClose, lmark not clean'
    !      stop
    !   endif
    !enddo
    !
    !     Initialize the stack
    !
    lmark(icellinit)=1_ip 
    lstack(1)=icellinit
    nstack=1_ip
    istack=0_ip
    !
    !     First extend until finding a point from the initial cell 
    !
    do

       if(istack==nstack)exit
       istack=istack+1_ip
       icell=lstack(istack)

       do jpo=lctop2(icell),lctop2(icell+1)-1
          if(jpo==0)then
             write(*,*)'Error in GetClose'
             stop
          endif
          jpoin=lctop1(jpo)
          !if(jpoin/=ipoin) exit
          goto 100
       enddo

       do i=1,6
          ineigh=lcell(icell)%neigh(i)
          if(ineigh==0)cycle
          if(lmark(ineigh)==1)cycle
          nstack=nstack+1
          lstack(nstack)=ineigh
          lmark(ineigh)=1_ip
       enddo

    enddo
100 continue     

    !
    !     Compute reference distance
    !
    rx=coor(1,jpoin)-pnew(1)
    ry=coor(2,jpoin)-pnew(2)
    rz=coor(3,jpoin)-pnew(3)
    rdistmin=sqrt(rx*rx+ry*ry+rz*rz)
    jmin=jpoin
    !
    !     Compute the local box
    !
    vmin(1)=pnew(1)-rdistmin
    vmin(2)=pnew(2)-rdistmin
    vmin(3)=pnew(3)-rdistmin
    vmax(1)=pnew(1)+rdistmin
    vmax(2)=pnew(2)+rdistmin
    vmax(3)=pnew(3)+rdistmin
    !
    !     Reset istack
    !  
    istack=istack-1_ip
    !
    !     Now that we have found jpoin, only add if distance is less
    !
    do

       if(istack==nstack)exit
       istack=istack+1_ip
       icell=lstack(istack)
       !
       !     Check against vmin,vmax
       !
       if(vmin(1)>lcell(icell)%coor(1,2))cycle
       if(vmin(2)>lcell(icell)%coor(2,2))cycle
       if(vmin(3)>lcell(icell)%coor(3,2))cycle
       if(vmax(1)<lcell(icell)%coor(1,1))cycle
       if(vmax(2)<lcell(icell)%coor(2,1))cycle
       if(vmax(3)<lcell(icell)%coor(3,1))cycle
       !
       !     Brute force inside icell
       ! 
       do jpo=lctop2(icell),lctop2(icell+1)-1 
          if(jpo==0)then
             write(*,*)'Error in GetClose 2'
             stop
          endif
          jpoin=lctop1(jpo)
          !if(ipoin==jpoin)cycle
          rx=coor(1,jpoin)-pnew(1)
          ry=coor(2,jpoin)-pnew(2)
          rz=coor(3,jpoin)-pnew(3)
          rdist=sqrt(rx*rx+ry*ry+rz*rz)
          if(rdist<rdistmin)then
             jmin=jpoin
             rdistmin=rdist
          endif
       enddo
       !
       !     Push neighbors into the stack
       !
       do i=1,6
          ineigh=lcell(icell)%neigh(i)
          if(ineigh==0)cycle
          if(lmark(ineigh)==1)cycle
          nstack=nstack+1
          lstack(nstack)=ineigh
          lmark(ineigh)=1_ip
       enddo

    enddo

    !
    !     Clean up lmark
    !
    do istack=1,nstack
       lmark(lstack(istack))=0_ip
    enddo

    ipclos=jmin


  end subroutine getClose

  subroutine getClosIn(pnew,icellinit,ndim,coor,npoin,lcell,ncell,lmark,lstack,&
       lctop1,lctop2,rdist,lstackp,nstackp,mstackp)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ndim,npoin,ncell,mstackp,icellinit
    integer(ip),intent(inout)     :: nstackp
    integer(ip),intent(in)        :: lctop1(*),lctop2(ncell+1)
    real(rp),intent(in)           :: coor(ndim,npoin),pnew(ndim)
    integer(ip),intent(inout)     :: lmark(ncell),lstack(ncell),lstackp(mstackp)
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: icell,jpo,jpoin,jmin,nstack,istack,i,ineigh
    real(rp)                      :: rx,ry,rz,vmin(ndim),vmax(ndim),rdist
    !
    !     This subroutine gives the points whose cartesian cell intersects
    !     the search area
    !

    !
    !     DBG
    ! 
    !do icell=1,ncell
    !   if(lmark(icell)/=0)then
    !      write(*,*)'Error in getCloseIn, lmark not clean'
    !      stop
    !   endif
    !enddo

    !
    !     Initialize the stack
    !
    lmark(icellinit)=1_ip 
    lstack(1)=icellinit
    nstack=1_ip
    istack=0_ip
    !
    !     Initialize nstackp
    !
    nstackp=0_ip
    !
    !     Compute the local box
    !
    vmin(1)=pnew(1)-rdist
    vmin(2)=pnew(2)-rdist
    vmin(3)=pnew(3)-rdist
    vmax(1)=pnew(1)+rdist
    vmax(2)=pnew(2)+rdist
    vmax(3)=pnew(3)+rdist
    !
    !     Now get the points in the search region
    !
    do

       if(istack==nstack)exit
       istack=istack+1_ip
       icell=lstack(istack)
       !
       !     Check against vmin,vmax
       !
       if(vmin(1)>lcell(icell)%coor(1,2))cycle
       if(vmin(2)>lcell(icell)%coor(2,2))cycle
       if(vmin(3)>lcell(icell)%coor(3,2))cycle
       if(vmax(1)<lcell(icell)%coor(1,1))cycle
       if(vmax(2)<lcell(icell)%coor(2,1))cycle
       if(vmax(3)<lcell(icell)%coor(3,1))cycle
       !
       !     Brute force inside icell
       ! 
       do jpo=lctop2(icell),lctop2(icell+1)-1 
          if(jpo==0)then
             write(*,*)'Error in GetClosin'
             stop
          endif
          jpoin=lctop1(jpo)
          !if(ipoin==jpoin)cycle
          if(nstackp==mstackp)goto 100
          nstackp=nstackp+1_ip
          lstackp(nstackp)=jpoin
       enddo
       !
       !     Push neighbors into the stack
       !
       do i=1,6
          ineigh=lcell(icell)%neigh(i)
          if(ineigh==0)cycle
          if(lmark(ineigh)==1)cycle
          nstack=nstack+1
          lstack(nstack)=ineigh
          lmark(ineigh)=1_ip
       enddo

    enddo
100 continue
    !
    !     Clean up lmark
    !
    do istack=1,nstack
       lmark(lstack(istack))=0_ip
    enddo

  end subroutine getClosIn

  subroutine getClosIn2(pnew,icellinit,ndim,coor,npoin,lcell,ncell,lmark,lstack,&
       lctop,rdist,lstackp,nstackp,mstackp,npmax,mcell,ipoin)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ndim,npoin,ncell,mstackp,icellinit,npmax,mcell
    integer(ip),intent(inout)     :: nstackp
    integer(ip),intent(in)        :: lctop(npmax,mcell),ipoin
    real(rp),intent(in)           :: coor(ndim,npoin),pnew(ndim),rdist
    integer(ip),intent(inout)     :: lmark(ncell),lstack(ncell),lstackp(mstackp)
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: icell,jcell,jpo,jpoin,jmin,nstack,istack,i
    integer(ip)                   :: ineigh,npcell,nptop,nstackp0,istackp
    real(rp)                      :: dx,dy,dz,vmin(ndim),vmax(ndim),rnl
    !
    !     This subroutine gives the points in the sphere around pnew of radio rdist 
    ! 

    !
    !     DBG
    ! 
    do icell=1,ncell
       if(lmark(icell)/=0)then
          write(*,*)'Error in getCloseIn, lmark not clean icell'
          stop
       endif
    enddo
    nptop=npmax-1_ip
    !
    !     Initialize the stack
    !
    lmark(icellinit)=1_ip 
    lstack(1)=icellinit
    nstack=1_ip
    istack=0_ip
    !
    !     Initialize nstackp
    !
    nstackp=0_ip
    !
    !     Compute the local box
    !
    vmin(1)=pnew(1)-rdist
    vmin(2)=pnew(2)-rdist
    vmin(3)=pnew(3)-rdist
    vmax(1)=pnew(1)+rdist
    vmax(2)=pnew(2)+rdist
    vmax(3)=pnew(3)+rdist
    !
    !     Now get the points in the search region
    !
    do

       if(istack==nstack)exit
       istack=istack+1_ip
       icell=lstack(istack)
       !
       !     Check against vmin,vmax
       !
       if(vmin(1)>lcell(icell)%coor(1,2))cycle
       if(vmin(2)>lcell(icell)%coor(2,2))cycle
       if(vmin(3)>lcell(icell)%coor(3,2))cycle
       if(vmax(1)<lcell(icell)%coor(1,1))cycle
       if(vmax(2)<lcell(icell)%coor(2,1))cycle
       if(vmax(3)<lcell(icell)%coor(3,1))cycle
       !
       !     Brute force inside icell
       !
       jcell=icell

       do 

          npcell=lctop(npmax,jcell)
          if(npcell>=0)then

             do jpo=1,npcell
                jpoin=lctop(jpo,jcell)
                if(nstackp==mstackp)then
                   write(*,*)'Error getClosIn2, nstackp=mstackp'
                   stop
                endif
                nstackp=nstackp+1_ip
                lstackp(nstackp)=jpoin
             enddo
             exit

          else

             do jpo=1,nptop
                jpoin=lctop(jpo,jcell)
                if(nstackp==mstackp)then
                   write(*,*)'Error getClosIn2, nstackp=mstackp'
                   stop
                endif
                nstackp=nstackp+1_ip
                lstackp(nstackp)=jpoin
             enddo

             jcell=-lctop(npmax,jcell)

          endif
       enddo
       !
       !     Push neighbors into the stack
       !
       do i=1,6
          ineigh=lcell(icell)%neigh(i)
          if(ineigh==0)cycle
          if(lmark(ineigh)==1)cycle
          nstack=nstack+1
          lstack(nstack)=ineigh
          lmark(ineigh)=1_ip
       enddo

    enddo
    !
    !     Clean up lmark
    !
    do istack=1,nstack
       lmark(lstack(istack))=0_ip
    enddo
    !
    !     Compress the stack
    !
    nstackp0=nstackp
    nstackp=0_ip
    do istackp=1,nstackp0
       jpoin=lstackp(istackp)
       if(jpoin/=ipoin)then
          dx=coor(1,jpoin)-pnew(1)
          dy=coor(2,jpoin)-pnew(2)
          dz=coor(3,jpoin)-pnew(3)
          rnl=sqrt(dx*dx+dy*dy+dz*dz)
          if(rnl<rdist)then
             nstackp=nstackp+1_ip
             lstackp(nstackp)=jpoin
          endif
       endif
    enddo

  end subroutine getClosIn2

  subroutine getClosIn3(pnew,icellinit,ndim,coor,npoin,lcell,ncell,lmark,lstack,&
       lctop,rdist,npmax,mcell,ipoin,ifound)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)        :: ndim,npoin,ncell,icellinit,npmax,mcell
    integer(ip),intent(in)        :: lctop(npmax,mcell),ipoin
    real(rp),intent(in)           :: coor(ndim,npoin),pnew(ndim),rdist
    integer(ip),intent(inout)     :: lmark(ncell),lstack(ncell),ifound
    type(cell)                    :: lcell(ncell)
    integer(ip)                   :: icell,jcell,jpo,jpoin,jmin,nstack,istack,i
    integer(ip)                   :: ineigh,npcell,nptop,nstackp0,istackp
    real(rp)                      :: dx,dy,dz,vmin(ndim),vmax(ndim),rnl
    !
    !     This subroutine returns as soon as a point has
    !     been found in the sphere around pnew of radio rdist 
    ! 

    !
    !     DBG
    ! 
    do icell=1,ncell
       if(lmark(icell)/=0)then
          write(*,*)'Error in getCloseIn, lmark not clean icell'
          stop
       endif
    enddo
    nptop=npmax-1_ip
    !
    !     Initialize the stack
    !
    lmark(icellinit)=1_ip 
    lstack(1)=icellinit
    nstack=1_ip
    istack=0_ip
    !
    !     Compute the local box
    !
    vmin(1)=pnew(1)-rdist
    vmin(2)=pnew(2)-rdist
    vmin(3)=pnew(3)-rdist
    vmax(1)=pnew(1)+rdist
    vmax(2)=pnew(2)+rdist
    vmax(3)=pnew(3)+rdist
    !
    !     Now get the points in the search region
    !
    do

       if(istack==nstack)exit
       istack=istack+1_ip
       icell=lstack(istack)
       !
       !     Check against vmin,vmax
       !
       if(vmin(1)>lcell(icell)%coor(1,2))cycle
       if(vmin(2)>lcell(icell)%coor(2,2))cycle
       if(vmin(3)>lcell(icell)%coor(3,2))cycle
       if(vmax(1)<lcell(icell)%coor(1,1))cycle
       if(vmax(2)<lcell(icell)%coor(2,1))cycle
       if(vmax(3)<lcell(icell)%coor(3,1))cycle
       !
       !     Brute force inside icell
       !
       jcell=icell

       do 

          npcell=lctop(npmax,jcell)
          if(npcell>=0)then

             do jpo=1,npcell
                jpoin=lctop(jpo,jcell)
                if(jpoin/=ipoin)then
                   dx=coor(1,jpoin)-pnew(1)
                   dy=coor(2,jpoin)-pnew(2)
                   dz=coor(3,jpoin)-pnew(3)
                   rnl=sqrt(dx*dx+dy*dy+dz*dz)
                   if(rnl<rdist)then
                      ifound=1_ip
                      goto 10
                   endif
                endif
             enddo
             exit

          else

             do jpo=1,nptop
                jpoin=lctop(jpo,jcell)
                if(jpoin/=ipoin)then
                   dx=coor(1,jpoin)-pnew(1)
                   dy=coor(2,jpoin)-pnew(2)
                   dz=coor(3,jpoin)-pnew(3)
                   rnl=sqrt(dx*dx+dy*dy+dz*dz)
                   if(rnl<rdist)then
                      ifound=1_ip
                      goto 10
                   endif
                endif
             enddo

             jcell=-lctop(npmax,jcell)

          endif

       enddo
       !
       !     Push neighbors into the stack
       !
       do i=1,6
          ineigh=lcell(icell)%neigh(i)
          if(ineigh==0)cycle
          if(lmark(ineigh)==1)cycle
          nstack=nstack+1
          lstack(nstack)=ineigh
          lmark(ineigh)=1_ip
       enddo

    enddo

10  continue

    !
    !     Clean up lmark
    !
    do istack=1,nstack
       lmark(lstack(istack))=0_ip
    enddo

  end subroutine getClosIn3

  subroutine rmvpnt(ipoin,lctop,ncell,mcell,icell,npmax)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)     :: ipoin,icell,npmax,ncell,mcell   
    integer(ip),intent(inout)  :: lctop(npmax,mcell)
    integer(ip)                :: nptop,iplace,jpo,npcell,jpoin
    !
    !     This sub takes out a point from its containing cartesian cell
    ! 
    !
    !     Get the correct place
    !
    nptop=npmax-1_ip
    !
    !     Get correct location
    !
    if(lctop(npmax,icell)<0)then

       write(*,*)'dynamic rmvpnt not yet ready'
       iplace=-lctop(npmax,icell)

       do
          if(lctop(npmax,iplace)>0)exit
          iplace=-lctop(npmax,iplace)
       enddo


    else
       !
       !     Get the point number in this cell
       !
       npcell=lctop(npmax,icell)
       !
       !     DBG 
       !
       if(npcell==0)then
          write(*,*)'Error in rmvpnt, npcell=0'
          stop
       endif
       !
       !     Find ipoin in lctop
       !
       do jpo=1,npcell
          jpoin=lctop(jpo,icell)
          if(jpoin==ipoin)exit
       enddo
       !
       !     Any error?
       ! 
       if(jpo==npcell+1)then
          write(*,*)'Error in rmvpnt, point not found' 
          stop
       endif
       !
       !     Update  lctop
       !
       lctop(jpo,icell)=lctop(npcell,icell)
       lctop(npmax,icell)=lctop(npmax,icell)-1_ip

    endif

  end subroutine rmvpnt

  subroutine addcellp(ipnew,coor,ndim,npoin,lcell,ncell,lctop,mcell,icellinit,&
       npmax,lcart,rtol)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip),intent(in)     :: ipnew,ndim,npoin,ncell,icellinit,npmax   
    integer(ip),intent(inout)  :: mcell,lcart(npoin)   
    type(cell)                 :: lcell(ncell)
    integer(ip),pointer        :: lctop(:,:)
    real(rp)                   :: coor(ndim,npoin),rtol
    integer(ip)                :: icart,nptop,jcart,npcel

    nptop=npmax-1_ip
    !
    !     Get the cell containing ipnew
    !
    icart=icellinit  
    call gtelem(ipnew,coor,npoin,ndim,lcell,ncell,icart,rtol)
    lcart(ipnew)=icart
    !
    !     Add it to lctop
    !    
    jcart=icart 

    do 

       npcel=lctop(npmax,jcart)

       if(npcel>=0)then

          if(npcel/=nptop)then 
             npcel=npcel+1_ip
             lctop(npcel,icart)=ipnew
             lctop(npmax,icart)=npcel
             return
          else
             mcell=mcell+1   
             call memrea(mcell,memor_msh,'LCTOP','addcellp',lctop)
             lctop(npmax,icart)=-mcell 
             lctop(1,mcell)=ipnew
             lctop(npmax,mcell)=1_ip
             return
          endif

       else   

          jcart=-npcel  

       endif

    enddo

  end subroutine addcellp

  subroutine writeGidconfp(ncell,lcell,coor,npoin,ndim)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)  :: ncell,npoin,ndim
    real(rp),intent(in)     :: coor(ndim,npoin)                    
    integer(ip)             :: icont,iconte,icont1,icont2,icont3,icont4
    integer(ip)             :: icont5,icont6,icont7,icont8,i,icontp,icontt,index1
    integer(ip)             :: neigh,level,pt1,pt2,j,k,ptnew,icontpp,iter,ncell8,jpos,jcont,tab4(3)
    integer(ip)             :: tab(4,6),tab2(4,6),icell,ilevel,ineigh,jlevel,jneigh,ipoin
    integer(ip)             :: tab3(6),tab5(4,6),tab6(4,6),ipo
    integer(ip),pointer     :: renum(:),renum2(:),lmark(:),lhang(:,:),lsize(:)
    real(rp), pointer       :: rsize(:)                   
    integer(4)              :: istat
    type(cell)              :: lcell(*)



    open(unit=50,file='cartpGid.msh',status='unknown')
    rewind 50

    ncell8=ncell*8 
    allocate(renum(ncell8),stat=istat)
    call memchk(zero,istat,memor_msh,'RENUM','writeGidconf',renum)
    allocate(renum2(ncell8),stat=istat)
    call memchk(zero,istat,memor_msh,'RENUM2','writeGidconf',renum2)

    do i=1,ncell8
       renum(i)=i
    enddo

    tab(1,1)=1
    tab(2,1)=5
    tab(3,1)=8
    tab(4,1)=4
    tab(1,2)=2
    tab(2,2)=3
    tab(3,2)=7
    tab(4,2)=6
    tab(1,3)=1
    tab(2,3)=2
    tab(3,3)=6
    tab(4,3)=5
    tab(1,4)=3
    tab(2,4)=4
    tab(3,4)=8
    tab(4,4)=7
    tab(1,5)=4
    tab(2,5)=3
    tab(3,5)=2
    tab(4,5)=1
    tab(1,6)=5
    tab(2,6)=6
    tab(3,6)=7
    tab(4,6)=8


    tab2(1,1)=2
    tab2(2,1)=6
    tab2(3,1)=7
    tab2(4,1)=3
    tab2(1,2)=1
    tab2(2,2)=4
    tab2(3,2)=8
    tab2(4,2)=5
    tab2(1,3)=4
    tab2(2,3)=3
    tab2(3,3)=7
    tab2(4,3)=8
    tab2(1,4)=2
    tab2(2,4)=1
    tab2(3,4)=5
    tab2(4,4)=6
    tab2(1,5)=8
    tab2(2,5)=7
    tab2(3,5)=6
    tab2(4,5)=5
    tab2(1,6)=1
    tab2(2,6)=2
    tab2(3,6)=3
    tab2(4,6)=4

    tab3(1)=2
    tab3(2)=1
    tab3(3)=4
    tab3(4)=3
    tab3(5)=6
    tab3(6)=5

    tab5(1,1)=1
    tab5(2,1)=5
    tab5(3,1)=4
    tab5(4,1)=8
    tab5(1,2)=2
    tab5(2,2)=6
    tab5(3,2)=3
    tab5(4,2)=7
    tab5(1,3)=1
    tab5(2,3)=5
    tab5(3,3)=2
    tab5(4,3)=6
    tab5(1,4)=4
    tab5(2,4)=8
    tab5(3,4)=3
    tab5(4,4)=7
    tab5(1,5)=1
    tab5(2,5)=4
    tab5(3,5)=2
    tab5(4,5)=3
    tab5(1,6)=5
    tab5(2,6)=8
    tab5(3,6)=6
    tab5(4,6)=7

    tab6(1,1)=2
    tab6(2,1)=6
    tab6(3,1)=3
    tab6(4,1)=7
    tab6(1,2)=1
    tab6(2,2)=5
    tab6(3,2)=4
    tab6(4,2)=8
    tab6(1,3)=4
    tab6(2,3)=8
    tab6(3,3)=3
    tab6(4,3)=7
    tab6(1,4)=1
    tab6(2,4)=5
    tab6(3,4)=2
    tab6(4,4)=6
    tab6(1,5)=5
    tab6(2,5)=8
    tab6(3,5)=6
    tab6(4,5)=7
    tab6(1,6)=1
    tab6(2,6)=4
    tab6(3,6)=2
    tab6(4,6)=3
    !
    !     Conformize the cells
    !
    do
       !
       !     Initialize the "something done" check 
       !
       iter=0
       !
       !     Loop on the cells
       !
       do icell=1,ncell

          ilevel=lcell(icell)%level
          !
          !     Loop on the neighbors
          !
          do i=1,6

             ineigh=lcell(icell)%neigh(i)

             if(ineigh==0)cycle

             jlevel=lcell(ineigh)%level

             if(jlevel==ilevel)then
                !
                !     Loop on the points of the faces
                !
                do k=1,4

                   pt1=(icell-1)*8+tab(k,i) 
                   pt2=(ineigh-1)*8+tab2(k,i)
                   if(renum(pt1)<renum(pt2))then
                      renum(pt2)=renum(pt1)
                      iter=1
                   else if(renum(pt2)<renum(pt1))then   
                      renum(pt1)=renum(pt2)
                      iter=1
                   endif

                enddo

             else if(ilevel<jlevel)then

                jpos=ineigh+1_ip
                jcont=1_ip
                do j=1,8
                   if(lcell(jpos)%neigh(tab3(i))==icell)then
                      tab4(jcont)=jpos
                      jcont=jcont+1_ip
                      if(jcont==4)exit 
                   endif
                   jpos=jpos+1  
                enddo

                pt1=(icell-1)*8+tab5(1,i)
                pt2=(ineigh-1)*8+tab6(1,i)
                if(renum(pt1)<renum(pt2))then
                   renum(pt2)=renum(pt1)
                   iter=1
                else if(renum(pt2)<renum(pt1))then   
                   renum(pt1)=renum(pt2)
                   iter=1
                endif

                pt1=(icell-1)*8+tab5(2,i)
                pt2=(tab4(1)-1)*8+tab6(2,i)
                if(renum(pt1)<renum(pt2))then
                   renum(pt2)=renum(pt1)
                   iter=1
                else if(renum(pt2)<renum(pt1))then   
                   renum(pt1)=renum(pt2)
                   iter=1
                endif

                pt1=(icell-1)*8+tab5(3,i)
                pt2=(tab4(2)-1)*8+tab6(3,i)
                if(renum(pt1)<renum(pt2))then
                   renum(pt2)=renum(pt1)
                   iter=1
                else if(renum(pt2)<renum(pt1))then   
                   renum(pt1)=renum(pt2)
                   iter=1
                endif

                pt1=(icell-1)*8+tab5(4,i)
                pt2=(tab4(3)-1)*8+tab6(4,i)
                if(renum(pt1)<renum(pt2))then
                   renum(pt2)=renum(pt1)
                   iter=1
                else if(renum(pt2)<renum(pt1))then   
                   renum(pt1)=renum(pt2)
                   iter=1
                endif

             endif
          enddo
       enddo

       do i=1,ncell8
          renum(i)=renum(renum(i))
       enddo

       if(iter==0)exit
       exit  

    enddo

    icont=0
    do i=1,ncell8
       if(renum(i)==i)then
          icont=icont+1
          renum2(i)=icont
       endif
    enddo

    write(50,1)
    write(50,2)
    write(50,3)

    icont=0
    icontpp=0 

    do i=1,ncell

       icontpp=icontpp+1  
       if(renum(icontpp)==icontpp)then
          icont=icont+1
          write(50,100)icont,lcell(i)%coor(1,1),lcell(i)%coor(2,1),lcell(i)%coor(3,1)
       endif

       icontpp=icontpp+1  
       if(renum(icontpp)==icontpp)then
          icont=icont+1
          write(50,100)icont,lcell(i)%coor(1,2),lcell(i)%coor(2,1),lcell(i)%coor(3,1)
       endif

       icontpp=icontpp+1  
       if(renum(icontpp)==icontpp)then
          icont=icont+1
          write(50,100)icont,lcell(i)%coor(1,2),lcell(i)%coor(2,2),lcell(i)%coor(3,1)
       endif

       icontpp=icontpp+1  
       if(renum(icontpp)==icontpp)then
          icont=icont+1
          write(50,100)icont,lcell(i)%coor(1,1),lcell(i)%coor(2,2),lcell(i)%coor(3,1)
       endif

       icontpp=icontpp+1  
       if(renum(icontpp)==icontpp)then
          icont=icont+1
          write(50,100)icont,lcell(i)%coor(1,1),lcell(i)%coor(2,1),lcell(i)%coor(3,2)
       endif

       icontpp=icontpp+1  
       if(renum(icontpp)==icontpp)then
          icont=icont+1
          write(50,100)icont,lcell(i)%coor(1,2),lcell(i)%coor(2,1),lcell(i)%coor(3,2)
       endif

       icontpp=icontpp+1  
       if(renum(icontpp)==icontpp)then
          icont=icont+1
          write(50,100)icont,lcell(i)%coor(1,2),lcell(i)%coor(2,2),lcell(i)%coor(3,2)
       endif

       icontpp=icontpp+1  
       if(renum(icontpp)==icontpp)then
          icont=icont+1
          write(50,100)icont,lcell(i)%coor(1,1),lcell(i)%coor(2,2),lcell(i)%coor(3,2)
       endif

    enddo
    write(50,4)


    icontp=icont 

    do i=1,ncell8
       renum(i)=renum2(renum(i))
    enddo
    write(50,5)

    icont=0
    icontpp=0
    do i=1,ncell 
       icont=icont+1
       icontpp=icontpp+1 
       icont1=renum(icontpp)
       icontpp=icontpp+1 
       icont2=renum(icontpp)  
       icontpp=icontpp+1 
       icont3=renum(icontpp)  
       icontpp=icontpp+1 
       icont4=renum(icontpp)  
       icontpp=icontpp+1 
       icont5=renum(icontpp)  
       icontpp=icontpp+1 
       icont6=renum(icontpp)  
       icontpp=icontpp+1 
       icont7=renum(icontpp)  
       icontpp=icontpp+1 
       icont8=renum(icontpp)  
       !index1=-lcell(i)%marked
       index1=lcell(i)%level 
       write(50,200) icont,icont1,icont2,icont3,icont4,icont5,icont6,icont7,icont8,index1
    enddo
    write(50,6)

    icontt=icont

    !     Faces of the triangulation

    write(50,7)
    write(50,2)
    write(50,3)


    icont=icontp
    do i=1,npoin
       icont=icont+1
       write(50,100)icont,coor(1,i),coor(2,i),coor(3,i)
    enddo
    write(50,4) 
    write(50,5) 

    icont=icontp
    do i=1,npoin
       icont=icont+1
       write(50,300)icont,icont     
    enddo
    write(50,6)      


    close(50)

    !
    !     Allocate lmark for the hanging nodes
    !
    allocate(lmark(icontp),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','writeGidconf',lmark)
    allocate(lhang(4,icontp),stat=istat)
    call memchk(zero,istat,memor_msh,'LHANG','writeGidconf',lhang)

    !
    !     Get the hanging nodes 
    !
    open(unit=60,file='carthan.msh',status='unknown')
    rewind 60

    do icell=1,ncell
       ilevel=lcell(icell)%level
       !
       !     In X -
       !  
       ineigh=lcell(icell)%neigh(1)
       if(ineigh/=0)then
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(2)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 1'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(2)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo
             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+7_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+2_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+6_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+6_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+2_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+6_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+3_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+2_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+6_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
          endif
       endif
       !
       !     In X +
       !  

       ineigh=lcell(icell)%neigh(2)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(1)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 2'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(1)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo

             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+8_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+4_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+8_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+5_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+4_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+4_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+8_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+8_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+8_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+8_ip) 
             endif
          endif
       endif
       !
       !     In Y -
       !  

       ineigh=lcell(icell)%neigh(3)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(4)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 3'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(4)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo

             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+7_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+8_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+3_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
          endif
       endif
       !
       !     In Y +
       !  
       ineigh=lcell(icell)%neigh(4)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(3)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 4'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(3)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo

             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+6_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+6_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+5_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+2_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+6_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+6_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+6_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+6_ip) 
             endif
          endif
       endif
       !
       !     In Z -
       !  
       ineigh=lcell(icell)%neigh(5)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(6)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 5'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(6)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo

             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+7_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+6_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+8_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+6_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+5_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+6_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+7_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+6_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
             endif
          endif
       endif
       !
       !     In Z +
       !  
       ineigh=lcell(icell)%neigh(6)
       if(ineigh/=0)then  
          jlevel=lcell(ineigh)%level
          if(ilevel<jlevel)then
             jneigh=lcell(ineigh)%neigh(5)
             if(icell/=jneigh)then 
                write(*,*)'Error conformity in hanging nodes 6'
                stop
             endif

             jpos=ineigh+1_ip
             jcont=1_ip
             do j=1,8
                if(lcell(jpos)%neigh(5)==icell)then
                   tab4(jcont)=jpos
                   jcont=jcont+1_ip
                   if(jcont==4)exit 
                endif
                jpos=jpos+1  
             enddo

             !
             !     Node at the center of the face
             !     
             ipoin=renum(8*(ineigh-1)+3_ip) 
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+4_ip) 
               ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
               ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+3_ip) 
             endif
             !
             !     Nodes at the edges
             !
             ipoin=renum(8*(ineigh-1)+4_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+4_ip) 
             endif
             ipoin=renum(8*(ineigh-1)+2_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
             endif
             ipoin=renum(8*(tab4(1)-1)+3_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+4_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+3_ip) 
             endif
             ipoin=renum(8*(tab4(2)-1)+3_ip)
             if(lmark(ipoin)==0)then
                lmark(ipoin)=1_ip
               ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
               ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+3_ip) 
             endif
          endif
       endif
    enddo
    !
    !     Write the marked points
    !
    do ipoin=1,icontp
       if(lmark(ipoin)==1)then
          write(60,400)ipoin,lhang(1,ipoin),lhang(2,ipoin),lhang(3,ipoin),lhang(4,ipoin)
       endif
    enddo

    close(60)


    !
    !     Output the size
    !

    !open(unit=70,file='cartGid.res',status='unknown')
    !rewind 70

    !allocate(rsize(icontp),stat=istat)
    !call memchk(zero,istat,memor_msh,'RSIZE','writeGidconf',rsize)
    !allocate(lsize(icontp),stat=istat)
    !call memchk(zero,istat,memor_msh,'LSIZE','writeGidconf',lsize)

    !do icell=1,ncell
    !   do ipo=1,8
    !      ipoin=renum(8*(icell-1)+ipo)
    !      rsize(ipoin)=rsize(ipoin)+lcell(icell)%rsize
    !      lsize(ipoin)=lsize(ipoin)+1
    !   enddo
    !enddo

    !do  ipoin=1,icontp
    !   if(lsize(ipoin)==0)then
    !      write(*,*)'Error in writeGidconf, lsize=0 at ipoin:',ipoin
    !      stop 
    !   endif
    !   rsize(ipoin)=rsize(ipoin)/real(lsize(ipoin))
    !enddo


    !write(70,10)
    !write(70,11)
    !write(70,13)
    !do  ipoin=1,icontp
    !   write(70,500)ipoin,rsize(ipoin)
    !enddo
    !write(70,14)
    !write(70,15)
    !close(70)

    !call memchk(2_ip,istat,memor_msh,'LSIZE','writeGidconf',lsize)
    !deallocate(lsize,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LSIZE','writeGidconf',0_ip)
    !call memchk(2_ip,istat,memor_msh,'RSIZE','writeGidconf',rsize)
    !deallocate(rsize,stat=istat)
    !if(istat/=0) call memerr(2_ip,'RSIZE','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHANG','writeGidconf',lhang)
    deallocate(lhang,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHANG','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','writeGidconf',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RENUM','writeGidconf',renum)
    deallocate(renum,stat=istat)
    if(istat/=0) call memerr(2_ip,'RENUM','writeGidconf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RENUM2','writeGidconf',renum2)
    deallocate(renum2,stat=istat)
    if(istat/=0) call memerr(2_ip,'RENUM2','writeGidconf',0_ip)



1   format('MESH dimension 3 ElemType  Hexahedra Nnode 8')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(10i10)
300 format(4i10)
400 format(5i10)
500 format(i10,e20.10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
7   format('MESH dimension 3 ElemType Point Nnode 1')


10  format('GID Post Results File 1.0')
11  format('Result "Size" "Analysis/time" 1 Scalar OnNodes')
13  format('Values')
14  format('End Values')
15  format('   ')

  end subroutine writeGidconfp

end module mod_cart
