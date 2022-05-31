module mod_interp
  use mod_cart

contains

  subroutine inter(elem,elold,nelem,nnode,neold,npoin,npold,coor,coold,ndim,lboup,&
       nboup,lelem,rint)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only    : memor_msh
    implicit none
    integer(ip),intent(in)    :: nelem,neold,npoin,npold,nnode,ndim,nboup
    integer(ip),intent(in)    :: elem(nnode,nelem),elold(nnode,neold)
    integer(ip),intent(in)    :: lboup(nboup)
    integer(ip),intent(inout) :: lelem(npoin)
    real(rp),intent(in)       :: coor(ndim,npoin),coold(ndim,npold)
    real(rp),intent(inout)    :: rint(nnode,nelem)
    integer(ip),pointer       :: lstack(:),pt1old(:),pt2old(:),eltoel(:,:)
    integer(ip),pointer       :: lcart(:),lctop1(:),lctop2(:),lcmark(:),lcstack(:)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),lbpoin(:),lmark(:) 
    type(cell),pointer        :: lcell(:)
    integer(ip)               :: inode,ie,ielem,jpoin,ihost,iguess,ncell 
    integer(ip)               :: istack,nstack,ipoin,ipclos,ipold,ivar,iboup 
    integer(ip)               :: ip1,ip2,ip3,ip4,npmax,nitermax,nstackp,iter,ierr
    real(rp)                  :: d1,d2,d3,d4,rx,ry,rz,rinit,csearch,rtol,bbox(ndim,2)
    integer(4)                :: istat
    integer(ip),parameter     :: mstackp=500
    integer(ip)               :: lstackp(mstackp)               
    !
    !     This subroutine interpolates the mesh elem from the mesh elold
    !
    !
    !     Allocate arrays for the new points
    !
    allocate(lstack(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','inter',lstack)
    allocate(lbpoin(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LBPOIN','inter',lbpoin)
    !
    !     Allocate arrays for the old points
    !
    allocate(lcart(npold),stat=istat)
    call memchk(zero,istat,memor_msh,'LCART','inter',lcart)
    allocate(lmark(neold),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','inter',lmark)
    !
    !     Max points per cell
    !
    npmax=8_ip
    ierr=0
    !
    !     Max iter for the close points
    !
    nitermax=5
    !
    !     Increase for the search size
    !
    csearch=2.0d+00
    !
    !     Mark the boundary points
    !
    do iboup=1,nboup
       lbpoin(lboup(iboup))=1_ip
    enddo
    !
    !     Get the elements surrounding the points for the old mesh
    !
    call ptoelm(elold,neold,npold,nnode,pt1old,pt2old)
    !
    !     Get the elements surrounding elements for the old mesh
    !
    call tetote(elold,nnode,neold,pt1old,pt2old,npold,eltoel)
    !
    !     Get the elements surrounding the points for the new mesh
    !
    call ptoelm(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !     Build the cartesian mesh with coold and elold 
    !
    call cartp(npold,ndim,elold,nnode,neold,npmax,lcart,pt1old,pt2old,&
         coold,lcell,ncell,rtol,0_ip,bbox)
    !
    !     Build the cell to old point pointers
    !
    call ctopnt(lcart,npold,ncell,lctop1,lctop2 )
    !
    !     Interpolate the new points in the cartesian mesh
    !
    call memrea(npoin,memor_msh,'LCART','inter',lcart)
    call intercart(coor,npoin,ndim,elem,nelem,nnode,lcell,ncell,lcart,ptoel1,ptoel2,rtol)
    !
    !     Allocate the arrays for the cartesian mesh
    !
    allocate(lcmark(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LCMARK','inter',lcmark)
    allocate(lcstack(ncell),stat=istat)
    call memchk(zero,istat,memor_msh,'LCSTACK','inter',lcstack)
    !
    !     Initialize the stack
    !   
    do ipoin=1,npoin
       if(lbpoin(ipoin)==0)exit
    enddo
    nstack=1_ip
    lstack(nstack)=ipoin
    !
    !     First interpolate the non boundary points
    !
    call interp(nnode,nelem,elem,ndim,neold,elold,coor,coold,npoin,npold,ptoel1,ptoel2,&
         lstack,nstack,eltoel,lelem,pt1old,pt2old,rint,lcart,lcell,ncell,lcmark,&
         lcstack,lctop1,lctop2,nitermax,lstackp,mstackp,csearch,lbpoin,lmark,ierr)
         if(ierr==1)then
            write(*,*)'Error in interpolating interior points'
            stop
         endif 

    !
    !     Initialize the new stack
    ! 
    do iboup=1,nboup 
       lstack(iboup)=lboup(iboup)
    enddo
    nstack=nboup
    !
    !     Then interpolate the boundary points
    !
    call interp(nnode,nelem,elem,ndim,neold,elold,coor,coold,npoin,npold,ptoel1,ptoel2,&
         lstack,nstack,eltoel,lelem,pt1old,pt2old,rint,lcart,lcell,ncell,lcmark,&
         lcstack,lctop1,lctop2,nitermax,lstackp,mstackp,csearch,lbpoin,lmark,ierr)
         if(ierr==1)then
            write(*,*)'Error in interpolating boundary points'
            stop
         endif 
    !
    !     If some points still have not been interpolated, mark them as boundary points
    ! 
    nstack=0_ip
    do ipoin=1,npoin
       if(lelem(ipoin)==0)then
          nstack=nstack+1_ip
          lstack(nstack)=1_ip
          lbpoin(ipoin)=1_ip 
       endif
    enddo
    !
    !      Interpolate this last serie
    !
    call interp(nnode,nelem,elem,ndim,neold,elold,coor,coold,npoin,npold,ptoel1,ptoel2,&
         lstack,nstack,eltoel,lelem,pt1old,pt2old,rint,lcart,lcell,ncell,lcmark,&
         lcstack,lctop1,lctop2,nitermax,lstackp,mstackp,csearch,lbpoin,lmark,ierr)
         if(ierr==1)then
            write(*,*)'Error in interpolating last points'
            stop
         endif 

    call memchk(2_ip,istat,memor_msh,'LBPOIN','inter',lbpoin)
    deallocate(lbpoin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBPOIN','inter',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCSTACK','inter',lcstack)
    deallocate(lcstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCSTACK','inter',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCMARK','inter',lcmark)
    deallocate(lcmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCMARK','inter',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','inter',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','inter',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCART','inter',lcart)
    deallocate(lcart,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCART','inter',0_ip)
    call memchk(2_ip,istat,memor_msh,'PT1OLD','inter',pt1old)
    deallocate(pt1old,stat=istat)
    if(istat/=0) call memerr(2_ip,'PT1OLD','inter',0_ip)
    call memchk(2_ip,istat,memor_msh,'PT2OLD','inter',pt2old)
    deallocate(pt2old,stat=istat)
    if(istat/=0) call memerr(2_ip,'PT2OLD','inter',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','inter',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','inter',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','inter',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','inter',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','inter',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','inter',0_ip)

  end subroutine inter

  subroutine brutint(elem,nelem,nnode,pnew,ndim,npoin,ihost,coor,d1,d2,d3,d4)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode
    integer(ip),intent(in)     :: elem(nnode,nelem)
    integer(ip),intent(inout)  :: ihost
    real(rp),intent(in)        :: coor(ndim,npoin),pnew(ndim)
    real(rp),intent(inout)     :: d1,d2,d3,d4 
    integer(ip)                :: ielem,ip1,ip2,ip3,ip4
    real(rp)                   :: p1(ndim),p2(ndim),p3(ndim)
    real(rp)                   :: rtol
    !
    !     This subroutine finds the host element by checking all the elements 
    !
    rtol=-1.0d-04

    do ielem=1,nelem

       call gtshape(elem,nelem,ielem,nnode,pnew,ndim,npoin,coor,d1,d2,d3,d4)

       if(d1>rtol .and. d2>rtol .and. d3>rtol .and. d4>rtol)then
          ihost=ielem
          return
       endif

    enddo

    !
    !     We still did not reach the element with the brute force...
    !
    ihost=0_ip


  end subroutine brutint

  subroutine chkSurElem(elem,nelem,nnode,pnew,ndim,npoin,coor,ipclos,ihost,ptoel1,ptoel2,d1,d2,d3,d4)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode,ipclos
    integer(ip),intent(in)     :: elem(nnode,nelem)
    integer(ip),intent(in)     :: ptoel1(*),ptoel2(npoin+1)
    integer(ip),intent(inout)  :: ihost
    real(rp),intent(in)        :: coor(ndim,npoin),pnew(ndim)
    real(rp),intent(inout)     :: d1,d2,d3,d4
    integer(ip)                :: ie,ielem
    real(rp)                   :: rtol

    rtol=-1.0d-04
    !
    !     Check the elements surrouding the closest point
    !  
    do ie=ptoel2(ipclos),ptoel2(ipclos+1)-1
       ielem=ptoel1(ie)

       call gtshape(elem,nelem,ielem,nnode,pnew,ndim,npoin,coor,d1,d2,d3,d4)

       if(d1>rtol .and. d2>rtol .and. d3>rtol .and. d4>rtol)then
          ihost=ielem
          return
       endif

    enddo

    !
    !     We still did not reach the element with the brute force...
    !
    ihost=0_ip

  end subroutine chkSurElem

  subroutine gtshape(elem,nelem,ielem,nnode,pnew,ndim,npoin,coor,d1,d2,d3,d4)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode
    integer(ip),intent(in)     :: elem(nnode,nelem)
    real(rp),intent(in)        :: coor(ndim,npoin),pnew(ndim)
    real(rp),intent(inout)     :: d1,d2,d3,d4
    integer(ip)                :: ielem,ip1,ip2,ip3,ip4
    real(rp)                   :: p1(ndim),p2(ndim),p3(ndim)
    real(rp)                   :: dtot,c10

    c10=1.0d+00
    !
    !     This subroutine compute the shape functions 
    !
    ip1=elem(1,ielem)
    ip2=elem(2,ielem)
    ip3=elem(3,ielem)
    ip4=elem(4,ielem)
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
    dtot=c10/dtot

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
    d1=d1*dtot 

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
    d2=d2*dtot 

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
    d3=d3*dtot 

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
    d4=d4*dtot 

  end subroutine gtshape

  subroutine chkSurElem2(elem,nelem,nnode,pnew,ndim,npoin,coor,ihost,ptoel1,ptoel2,&
       nstackp,lstackp,d1,d2,d3,d4)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode,nstackp
    integer(ip),intent(in)     :: elem(nnode,nelem),lstackp(nstackp)
    integer(ip),intent(in)     :: ptoel1(*),ptoel2(npoin+1)
    integer(ip),intent(inout)  :: ihost
    real(rp),intent(in)        :: coor(ndim,npoin),pnew(ndim)
    real(rp),intent(inout)     :: d1,d2,d3,d4
    integer(ip)                :: ie,ielem,istackp,ipoin
    real(rp)                   :: rtol

    rtol=-1.0d-04
    !
    !     Loop on surrounding points
    !
    do istackp=1,nstackp

       ipoin=lstackp(istackp)
       !
       !     Check the elements surrouding the closest point
       !  
       do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(ie)

          call gtshape(elem,nelem,ielem,nnode,pnew,ndim,npoin,coor,d1,d2,d3,d4)

          if(d1>rtol .and. d2>rtol .and. d3>rtol .and. d4>rtol)then
             ihost=ielem
             return
          endif

       enddo

    enddo
    !
    !     We still did not reach the element with the brute force...
    !
    ihost=0_ip

  end subroutine chkSurElem2

  subroutine interp(nnode,nelem,elem,ndim,neold,elold,coor,coold,npoin,npold,ptoel1,ptoel2,&
       lstack,nstack,eltoel,lelem,pt1old,pt2old,rint,lcart,lcell,ncell,lcmark,&
       lcstack,lctop1,lctop2,nitermax,lstackp,mstackp,csearch,lbpoin,lmark,ierr)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode,npold,neold,ncell,nitermax,mstackp
    integer(ip),intent(in)     :: elem(nnode,nelem),elold(nnode,neold)
    integer(ip),intent(in)     :: ptoel1(*),ptoel2(npoin+1),eltoel(nnode,neold)
    integer(ip),intent(in)     :: pt1old(*),pt2old(npold+1),lbpoin(npoin)
    integer(ip),intent(in)     :: lctop1(*),lctop2(ncell+1)
    type(cell)                 :: lcell(ncell) 
    integer(ip),intent(inout)  :: lcart(npoin),lcmark(ncell),lcstack(ncell),ierr 
    integer(ip),intent(inout)  :: nstack,lstack(npoin),lelem(npoin),lstackp(mstackp) 
    integer(ip),intent(inout)  :: lmark(neold) 
    real(rp),intent(in)        :: coor(ndim,npoin),coold(ndim,npold),csearch
    real(rp),intent(inout)     :: rint(nnode,npoin)
    integer(ip)                :: istack,ipoin,iguess,ihost,ipclos,iter
    integer(ip)                :: inode,jpoin,ielem,ie,nstackp
    real(rp)                   :: rtol,d1,d2,d3,d4,rx,ry,rz,rinit

    !
    !     Loop on the points to be interpolated
    !
    istack=0_ip

    do 

       if(istack==nstack)exit
       istack=istack+1

       ipoin=lstack(istack)
       iguess=lelem(ipoin)  

       call findelem(iguess,elold,neold,npold,coold,coor(1,ipoin),ndim,ihost,&
            nnode,eltoel,d1,d2,d3,d4)
       !
       !     Do we have a valid host?
       !
       if(ihost==0)then
          !
          !     Must go deeper, get the closest point  
          !
          call getClose(coor(1,ipoin),lcart(ipoin),ndim,coold,npold,lcell,ncell,&
               lcmark,lcstack,lctop1,lctop2,ipclos)
          !
          !     First check the elements surrounding the close point and its neighbors
          ! 
          call chkSurElem(elold,neold,nnode,coor(1,ipoin),ndim,npold,coold,ipclos,ihost,&
               pt1old,pt2old,d1,d2,d3,d4)
          !
          !     Do we have a valid host?
          !
          if(ihost==0)then
             !
             !     Compute initial search region
             !
             rx=coold(1,ipclos)-coor(1,ipoin)
             ry=coold(2,ipclos)-coor(2,ipoin)
             rz=coold(3,ipclos)-coor(3,ipoin)
             rinit=sqrt(rx*rx+ry*ry+rz*rz)
             !
             !     Iteratively increase the search region, testing the elements surrounding the close points
             !
             do iter=1,nitermax

                call getClosIn(coor(1,ipoin),lcart(ipoin),ndim,coold,npold,lcell,ncell,lcmark,&
                     lcstack,lctop1,lctop2,rinit,lstackp,nstackp,mstackp)
                !
                !     Did we reach mpstackp?
                !    
                if(nstackp==mstackp)then
                   ihost=0_ip
                   exit
                endif
                call chkSurElem2(elold,neold,nnode,coor(1,ipoin),ndim,npold,coold,ihost,pt1old,&
                     pt2old,nstackp,lstackp,d1,d2,d3,d4)

                if(ihost/=0)exit
                !
                !     Increase the search area
                !
                rinit=rinit*csearch

             enddo
             !
             !     Do we have a valid host?
             !
             if(ihost==0)then
                !
                !     Brute force
                !
                !call brutint(elold,neold,nnode,coor(1,ipoin),ndim,npold,ihost,coold,d1,d2,d3,d4)
                call distint(ndim,neold,npold,nnode,elold,eltoel,coold,coor(1,ipoin),iguess,&
                     ihost,lmark,d1,d2,d3,d4,ierr)
                if(ierr==1)then
                   write(*,*)'Error in interp point:',ipoin 
                   return
                endif
                !
                !     Do we finally have a valid host?
                !
                if(ihost==0)then
                   !
                   !     Outside the domain, interpolate with the nearest neighbor
                   !
                   lelem(ipoin)=-ipclos
                   ihost=-ipclos
                endif
             endif
          endif
       endif
       !
       !     Do we have a valid host element or we interpolate with the nearest neighbor?
       !
       if(ihost<0)cycle 
       !
       !     Store host element
       !
       lelem(ipoin)=ihost
       !
       !     Store interpolation coefficients
       ! 
       rint(1,ipoin)=d1
       rint(2,ipoin)=d2
       rint(3,ipoin)=d3
       rint(4,ipoin)=d4
       !
       !     The element has been found, mark the neighbor points with it
       !
       do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(ie)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             !
             !     Has jpoin been marked?
             !
             if(lelem(jpoin)==0)then
                lelem(jpoin)=ihost
                !
                !     Is jpoin a boundary point?
                !
                if(lbpoin(jpoin)==0)then
                   nstack=nstack+1
                   lstack(nstack)=jpoin
                endif
             endif
          enddo
       enddo

    enddo

  end subroutine interp

  subroutine distTet(elem,nelem,nnode,pnew,ndim,npoin,ihost,coor,d1,d2,d3,d4,icasemin,rdistmin,ipmin1,ipmin2,ielem)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode,ielem
    integer(ip),intent(in)     :: elem(nnode,nelem)
    integer(ip),intent(inout)  :: ihost,icasemin,ipmin1,ipmin2
    real(rp),intent(in)        :: coor(ndim,npoin),pnew(ndim)
    real(rp),intent(inout)     :: d1,d2,d3,d4,rdistmin 
    integer(ip)                :: ip1,ip2,ip3,ip4
    integer(ip)                :: icase(nnode),ipoin1(nnode),ipoin2(nnode)
    real(rp)                   :: p1(ndim),p2(ndim),p3(ndim),d(nnode)
    real(rp)                   :: rtol,rx,ry,rz,rnx,rny,rnz,rnl,c00
    real(rp)                   :: rface(ndim),pproj(ndim),dtot,rscal,rdistf
    real(rp)                   :: rdistt(3),df(4,nnode),rl,c10,de(4,3)
    real(rp)                   :: rdist(nnode),csca,cscal,epsil0,epsil1
    integer(ip)                :: ipoin1t(3),ipoin2t(3),icaset(3) 
    integer(ip)                :: ipos1,ipos2,ipos3,ipos4 
    integer(ip)                :: ltab(3,4),inode,icedg(3) 

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
    !     This subroutine gives the distance from pnew to ielem 
    !     On output - ihost gives the element number containing pnew if not 0,
    !               - icasemin is 0 if the point is inside the element  
    !                             1 if the smallest distance is for a point inside a face
    !                             2 if the smallest distance is for a point inside a side
    !                             3 if the smallest distance is for a point of the element
    !               - ipmin1 gives the point of smallest distance if icasemin=3 
    !               - ipmin1 and ipmin2 give the points of the edge of smallest distance
    !                 if icasemin=2 
    ! 
    rtol=-1.0d-08
    c00=0.0d+00
    c10=1.0d+00
    epsil1=1.0d+00+1.0d-08
    epsil0=1.0d-08
    !
    !     Initialize ihost, icase, ip1, ip2
    !
    ihost=0_ip
    icasemin=0_ip 
    ipmin1=0_ip
    ipmin2=0_ip
    !
    !     Compute shape function
    !
    call gtshape(elem,nelem,ielem,nnode,pnew,ndim,npoin,coor,d1,d2,d3,d4)
    !
    !     Are we inside the element?
    ! 
    if(d1>rtol .and. d2>rtol .and. d3>rtol .and. d4>rtol)then
       ihost=ielem
       return
    endif
    !
    !     As we are outside, compute the distance to faces
    ! 
    !
    !     Compute the projected point on the faces
    !  
    do inode=1,nnode
       !
       !     Get the points of the face
       !
       ipos1=ltab(1,inode)
       ipos2=ltab(2,inode)
       ipos3=ltab(3,inode)

       ip1=elem(ipos1,ielem)
       ip2=elem(ipos2,ielem)
       ip3=elem(ipos3,ielem)

       d(1)=c00
       d(2)=c00
       d(3)=c00
       d(4)=c00
       !
       !     Compute face normal
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)

       rnx= p1(2)*p2(3)-p1(3)*p2(2)
       rny=-p1(1)*p2(3)+p1(3)*p2(1)
       rnz= p1(1)*p2(2)-p1(2)*p2(1)
       rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
       rnl=c10/rnl
       rface(1)=rnx*rnl
       rface(2)=rny*rnl
       rface(3)=rnz*rnl
       !
       !     Compute face area to normalize
       !
       dtot=rnl
       !
       !     Get projected point and distance to face
       !
       rx=coor(1,ip1)-pnew(1) 
       ry=coor(2,ip1)-pnew(2) 
       rz=coor(3,ip1)-pnew(3) 
       rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
       rdistf=abs(rscal)         

       pproj(1)= pnew(1)-rscal*rface(1)
       pproj(2)= pnew(2)-rscal*rface(2)
       pproj(3)= pnew(3)-rscal*rface(3)
       !
       !     Get barycentric coordinates
       !
       !
       !     Side (ip2,ip3)
       !
       p1(1)=coor(1,ip2)-pproj(1)
       p1(2)=coor(2,ip2)-pproj(2)
       p1(3)=coor(3,ip2)-pproj(3)
       p2(1)=coor(1,ip3)-pproj(1)
       p2(2)=coor(2,ip3)-pproj(2)
       p2(3)=coor(3,ip3)-pproj(3)

       call orient3D(p1,p2,rface,d(ipos1),ndim)
       d(ipos1)=d(ipos1)*dtot
       !
       !     Side (ip3,ip1)
       !
       p1(1)=coor(1,ip3)-pproj(1)
       p1(2)=coor(2,ip3)-pproj(2)
       p1(3)=coor(3,ip3)-pproj(3)
       p2(1)=coor(1,ip1)-pproj(1)
       p2(2)=coor(2,ip1)-pproj(2)
       p2(3)=coor(3,ip1)-pproj(3)

       call orient3D(p1,p2,rface,d(ipos2),ndim)
       d(ipos2)=d(ipos2)*dtot
       !
       !     Side (ip1,ip2)
       !
       p1(1)=coor(1,ip1)-pproj(1)
       p1(2)=coor(2,ip1)-pproj(2)
       p1(3)=coor(3,ip1)-pproj(3)
       p2(1)=coor(1,ip2)-pproj(1)
       p2(2)=coor(2,ip2)-pproj(2)
       p2(3)=coor(3,ip2)-pproj(3)

       call orient3D(p1,p2,rface,d(ipos3),ndim)
       d(ipos3)=d(ipos3)*dtot
       !
       !     Is the projected point inside the face?
       !
       if(d(ipos1)>rtol .and. d(ipos2)>rtol .and. d(ipos3)>rtol )then

          icase(inode)=1_ip 
          rdist(inode)=rdistf
          df(1,inode)=d(1) 
          df(2,inode)=d(2) 
          df(3,inode)=d(3) 
          df(4,inode)=d(4)

       else

          !
          !     Evaluate distance to sides
          !
          de(1,1)=c00
          de(2,1)=c00
          de(3,1)=c00
          de(4,1)=c00
          de(1,2)=c00
          de(2,2)=c00
          de(3,2)=c00
          de(4,2)=c00
          de(1,3)=c00
          de(2,3)=c00
          de(3,3)=c00
          de(4,3)=c00
          !
          !     Side (ip2,ip3)
          !
          p1(1)=coor(1,ip3)-coor(1,ip2)
          p1(2)=coor(2,ip3)-coor(2,ip2)
          p1(3)=coor(3,ip3)-coor(3,ip2)
          rl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
          rl=c10/rl
          p1(1)=rl*p1(1) 
          p1(2)=rl*p1(2) 
          p1(3)=rl*p1(3) 

          p2(1)=pnew(1)-coor(1,ip2)
          p2(2)=pnew(2)-coor(2,ip2)
          p2(3)=pnew(3)-coor(3,ip2)

          csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
          cscal=csca*rl 

          if(cscal<-epsil0)then

             rdistt(1)=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
             ipoin1t(1)=ip2
             icedg(1)=1_ip
             de(ipos2,1)=c10 

          else if(cscal>epsil1)then

             p2(1)=pnew(1)-coor(1,ip3)
             p2(2)=pnew(2)-coor(2,ip3)
             p2(3)=pnew(3)-coor(3,ip3)
             rdistt(1)=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
             ipoin1t(1)=ip3
             icedg(1)=1_ip
             de(ipos3,1)=c10

          else

             p2(1)=coor(1,ip2)+csca*p1(1)-pnew(1)
             p2(2)=coor(2,ip2)+csca*p1(2)-pnew(2)
             p2(3)=coor(3,ip2)+csca*p1(3)-pnew(3)
             rdistt(1)=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 
             ipoin1t(1)=ip2
             ipoin2t(1)=ip3
             icedg(1)=2_ip  
             de(ipos2,1)=c10-cscal
             de(ipos3,1)=cscal 

          endif
          !
          !     Side (ip3,ip1)
          !
          p1(1)=coor(1,ip1)-coor(1,ip3)
          p1(2)=coor(2,ip1)-coor(2,ip3)
          p1(3)=coor(3,ip1)-coor(3,ip3)
          rl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
          rl=c10/rl
          p1(1)=rl*p1(1) 
          p1(2)=rl*p1(2) 
          p1(3)=rl*p1(3) 

          p2(1)=pnew(1)-coor(1,ip3)
          p2(2)=pnew(2)-coor(2,ip3)
          p2(3)=pnew(3)-coor(3,ip3)

          csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
          cscal=csca*rl 

          if(cscal<-epsil0)then

             rdistt(2)=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
             ipoin1t(2)=ip3
             icedg(2)=1_ip
             de(ipos3,2)=c10 

          else if(cscal>epsil1)then

             p2(1)=pnew(1)-coor(1,ip1)
             p2(2)=pnew(2)-coor(2,ip1)
             p2(3)=pnew(3)-coor(3,ip1)
             rdistt(2)=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
             ipoin1t(2)=ip1
             icedg(2)=1_ip
             de(ipos1,2)=c10 

          else

             p2(1)=coor(1,ip3)+csca*p1(1)-pnew(1)
             p2(2)=coor(2,ip3)+csca*p1(2)-pnew(2)
             p2(3)=coor(3,ip3)+csca*p1(3)-pnew(3)
             rdistt(2)=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 
             ipoin1t(2)=ip3
             ipoin2t(2)=ip1
             icedg(2)=2_ip
             de(ipos3,2)=c10-cscal 
             de(ipos1,2)=cscal

          endif
          !
          !     Side (ip1,ip2)
          !
          p1(1)=coor(1,ip2)-coor(1,ip1)
          p1(2)=coor(2,ip2)-coor(2,ip1)
          p1(3)=coor(3,ip2)-coor(3,ip1)
          rl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
          rl=c10/rl
          p1(1)=rl*p1(1) 
          p1(2)=rl*p1(2) 
          p1(3)=rl*p1(3) 

          p2(1)=pnew(1)-coor(1,ip1)
          p2(2)=pnew(2)-coor(2,ip1)
          p2(3)=pnew(3)-coor(3,ip1)

          csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
          cscal=csca*rl 

          if(cscal<-epsil0)then

             rdistt(3)=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
             ipoin1t(3)=ip1
             icedg(3)=1_ip
             de(ipos1,3)=c10 

          else if(cscal>epsil1)then

             p2(1)=pnew(1)-coor(1,ip2)
             p2(2)=pnew(2)-coor(2,ip2)
             p2(3)=pnew(3)-coor(3,ip2)
             rdistt(3)=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
             ipoin1t(3)=ip2
             icedg(3)=1_ip
             de(ipos2,3)=c10 

          else

             p2(1)=coor(1,ip1)+csca*p1(1)-pnew(1)
             p2(2)=coor(2,ip1)+csca*p1(2)-pnew(2)
             p2(3)=coor(3,ip1)+csca*p1(3)-pnew(3)
             rdistt(3)=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 
             ipoin1t(3)=ip1
             ipoin2t(3)=ip2
             icedg(3)=2_ip
             de(ipos1,3)=c10-cscal 
             de(ipos2,3)=cscal

          endif
          !
          !     Take the min over the side distances
          !
          if(rdistt(1)<rdistt(2))then 
             rdist(inode)=rdistt(1)
             if(icedg(1)==1)then
                icase(inode)=3_ip
                ipoin1(inode)=ipoin1t(1)
             else 
                icase(inode)=2_ip 
                ipoin1(inode)=ipoin1t(1)
                ipoin2(inode)=ipoin2t(1) 
             endif

             df(1,inode)=de(1,1)
             df(2,inode)=de(2,1)
             df(3,inode)=de(3,1)
             df(4,inode)=de(4,1)

          else 
             rdist(inode)=rdistt(2)
             if(icedg(2)==1)then
                icase(inode)=3_ip
                ipoin1(inode)=ipoin1t(2)
             else 
                icase(inode)=2_ip
                ipoin1(inode)=ipoin1t(2)
                ipoin2(inode)=ipoin2t(2) 
             endif

             df(1,inode)=de(1,2)
             df(2,inode)=de(2,2)
             df(3,inode)=de(3,2)
             df(4,inode)=de(4,2)

          endif

          if(rdistt(3)<rdist(inode))then 
             rdist(inode)=rdistt(3)
             if(icedg(3)==1)then
                icase(inode)=3_ip
                ipoin1(inode)=ipoin1t(3)
             else
                icase(inode)=2_ip 
                ipoin1(inode)=ipoin1t(3)
                ipoin2(inode)=ipoin2t(3)
             endif

             df(1,inode)=de(1,3)
             df(2,inode)=de(2,3)
             df(3,inode)=de(3,3)
             df(4,inode)=de(4,3)

          endif

       endif

    enddo

    !
    !     Find the min over the faces
    !
    if(rdist(1)<rdist(2))then

       rdistmin=rdist(1)
       icasemin=icase(1)
       if(icasemin==2)then
          ipmin1=ipoin1(1)
          ipmin2=ipoin2(1)
       else if(icasemin==3)then
          ipmin1=ipoin1(1)
       endif

       d1=df(1,1)
       d2=df(2,1)
       d3=df(3,1)
       d4=df(4,1)

    else

       rdistmin=rdist(2)
       icasemin=icase(2)
       if(icasemin==2)then
          ipmin1=ipoin1(2)
          ipmin2=ipoin2(2)
       else if(icasemin==3)then
          ipmin1=ipoin1(2)
       endif

       d1=df(1,2)
       d2=df(2,2)
       d3=df(3,2)
       d4=df(4,2)

    endif

    if(rdist(3)<rdistmin)then

       rdistmin=rdist(3)
       icasemin=icase(3)
       if(icasemin==2)then
          ipmin1=ipoin1(3)
          ipmin2=ipoin2(3)
       else if(icasemin==3)then
          ipmin1=ipoin1(3)
       endif

       d1=df(1,3)
       d2=df(2,3)
       d3=df(3,3)
       d4=df(4,3)

    endif

    if(rdist(4)<rdistmin)then
       rdistmin=rdist(4)
       icasemin=icase(4)
       if(icasemin==2)then
          ipmin1=ipoin1(4)
          ipmin2=ipoin2(4)
       else if(icasemin==3)then
          ipmin1=ipoin1(4)
       endif

       d1=df(1,4)
       d2=df(2,4)
       d3=df(3,4)
       d4=df(4,4)

    endif

  end subroutine distTet

  subroutine distint(ndim,nelem,npoin,nnode,elem,eltoel,coor,pnew,iguess,&
       ihost,lmark,d1,d2,d3,d4,ierr)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode,iguess
    integer(ip),intent(in)     :: elem(nnode,nelem),eltoel(nnode,nelem)
    integer(ip),intent(inout)  :: ihost,lmark(nelem),ierr
    real(rp),intent(inout)     :: d1,d2,d3,d4
    real(rp),intent(in)        :: coor(ndim,npoin),pnew(ndim)
    integer(ip),parameter      :: mstack=500
    integer(ip)                :: lstack(mstack),lstack2(mstack),nstack,nstack2
    integer(ip)                :: istack,jelem,ielem,ineigh 
    integer(ip)                :: icase,ip1,ip2,icaset,ip1t,ip2t,inode 
    real(rp)                   :: rdist,rdistt 
    real(rp)                   :: d1t,d2t,d3t,d4t
    !
    !     DBG
    !
    do ielem=1,nelem
       if(lmark(ielem)/=0)then
          write(*,*)'Error in distint'
          ierr=1
          return
       endif
    enddo
    !
    !     Initialize the stack
    !
    lstack(1)=iguess
    lstack2(1)=iguess
    nstack=1_ip
    nstack2=1_ip
    istack=0_ip
    lmark(iguess)=1_ip
    !
    !     Initialize the distance to the element iguess
    !
    call distTet(elem,nelem,nnode,pnew,ndim,npoin,ihost,coor,d1,d2,d3,d4,&
         icase,rdist,ip1,ip2,iguess)
    !
    !     Loop on the stack
    ! 
    do 

       if(istack==nstack)exit
       if(istack==mstack)then
          write(*,*)'Error in distint, istack==mstack'
          ierr=1
          return
       endif


       istack=istack+1_ip
       ielem=lstack(istack)
       !
       !     Loop on the neighbors
       !
       do inode=1,nnode
          ineigh=eltoel(inode,ielem)
          !
          !     Do we have a neighbor?
          !
          if(ineigh==0)cycle
          !
          !     Has the element been already checked ?
          ! 
          if(lmark(ineigh)==1)cycle
          lmark(ineigh)=1_ip
          nstack2=nstack2+1_ip
          if(nstack2==mstack)then
             write(*,*)'Error in distint, nstack2=mstack'
             ierr=1
             return
          endif
          lstack2(nstack2)=ineigh

          call distTet(elem,nelem,nnode,pnew,ndim,npoin,ihost,coor,d1t,d2t,d3t,d4t,&
               icaset,rdistt,ip1t,ip2t,ineigh)
          !
          !     Did we find the host element?
          !
          if(ihost/=0)then
             goto 100
          endif
          !
          !     Do we have the same configuration than the current minimum distance?
          !
          if(icaset==icase)then
             !
             !     Do we roll over an edge and is this the same edge? 
             !
             if((icaset==2).and.((ip1==ip1 .and. ip2==ip2t) .or. (ip1==ip2t .and.ip2==ip1t)))then
                ! 
                !     Add to the stack
                !    
                nstack=nstack+1_ip
                lstack(nstack)=ineigh   
                !
                !     Do we roll over a point and is this the same point?
                !
             else if(icase==3 .and. ip1==ip1t)then
                ! 
                !     Add to the stack
                !    
                nstack=nstack+1_ip
                lstack(nstack)=ineigh   

             else if(rdistt<rdist)then 
                !
                !     Rely on distance
                !
                nstack=nstack+1_ip
                lstack(nstack)=ineigh   
                rdist=rdistt
                icase=icaset 
                ip1=ip1t
                ip2=ip2t
                d1=d1t         ! Not usefull
                d2=d2t         !
                d3=d3t         !
                d4=d4t         !

             endif
          else if(rdistt<rdist)then
             !
             !     Rely on distance
             !
             nstack=nstack+1_ip
             lstack(nstack)=ineigh   
             rdist=rdistt
             icase=icaset 
             ip1=ip1t
             ip2=ip2t 
             d1=d1t      ! Not usefull
             d2=d2t      !
             d3=d3t      !
             d4=d4t      !

          endif

       enddo

    enddo

100 continue

    !
    !     Clean up lmark
    !
    do istack=1,nstack2
       lmark(lstack2(istack))=0_ip
    enddo

    !
    !     This is ihost
    !
    ihost=lstack(nstack) 


  end subroutine distint

  subroutine interpnt(nnode,nelem,elem,ndim,neold,elold,coor,coold,npoin,npold,&
       ipoin,eltoel,pt1old,pt2old,lcart,lcell,ncell,lcmark,&
       lcstack,lctop1,lctop2,nitermax,lstackp,mstackp,csearch,lmark,&
       ielem,d1,d2,d3,d4,iguess,ierr)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode,npold,neold,ncell,nitermax,mstackp
    integer(ip),intent(in)     :: elem(nnode,nelem),elold(nnode,neold),iguess
    integer(ip),intent(in)     :: eltoel(nnode,neold)
    integer(ip),intent(in)     :: pt1old(*),pt2old(npold+1)
    integer(ip),intent(in)     :: lctop1(*),lctop2(ncell+1)
    type(cell)                 :: lcell(ncell) 
    integer(ip),intent(inout)  :: lcart(npoin),lcmark(ncell),lcstack(ncell),ierr 
    integer(ip),intent(inout)  :: lstackp(mstackp),ielem 
    integer(ip),intent(inout)  :: lmark(neold) 
    real(rp),intent(in)        :: coor(ndim,npoin),coold(ndim,npold),csearch
    real(rp),intent(inout)     :: d1,d2,d3,d4
    integer(ip)                :: istack,ipoin,ihost,ipclos,iter
    integer(ip)                :: inode,jpoin,ie,nstackp
    real(rp)                   :: rtol,rx,ry,rz,rinit
    !
    !     Interpolated ipoin
    !
    call findelem(iguess,elold,neold,npold,coold,coor(1,ipoin),ndim,ihost,&
         nnode,eltoel,d1,d2,d3,d4)
    !
    !     Do we have a valid host?
    !
    if(ihost==0)then
       !
       !     Must go deeper, get the closest point  
       !
       call getClose(coor(1,ipoin),lcart(ipoin),ndim,coold,npold,lcell,ncell,&
            lcmark,lcstack,lctop1,lctop2,ipclos)
       !
       !     First check the elements surrounding the close point and its neighbors
       ! 
       call chkSurElem(elold,neold,nnode,coor(1,ipoin),ndim,npold,coold,ipclos,ihost,&
            pt1old,pt2old,d1,d2,d3,d4)
       !
       !     Do we have a valid host?
       !
       if(ihost==0)then
          !
          !     Compute initial search region
          !
          rx=coold(1,ipclos)-coor(1,ipoin)
          ry=coold(2,ipclos)-coor(2,ipoin)
          rz=coold(3,ipclos)-coor(3,ipoin)
          rinit=sqrt(rx*rx+ry*ry+rz*rz)
          !
          !     Iteratively increase the search region, testing the elements surrounding the close points
          !
          do iter=1,nitermax

             call getClosIn(coor(1,ipoin),lcart(ipoin),ndim,coold,npold,lcell,ncell,lcmark,&
                  lcstack,lctop1,lctop2,rinit,lstackp,nstackp,mstackp)
             !
             !     Did we reach mpstackp?
             !    
             if(nstackp==mstackp)then
                ihost=0_ip
                exit
             endif
             call chkSurElem2(elold,neold,nnode,coor(1,ipoin),ndim,npold,coold,ihost,pt1old,&
                  pt2old,nstackp,lstackp,d1,d2,d3,d4)

             if(ihost/=0)exit
             !
             !     Increase the search area
             !
             rinit=rinit*csearch

          enddo
          !
          !     Do we have a valid host?
          !
          if(ihost==0)then
             !
             !     Brute force
             !
             !call brutint(elold,neold,nnode,coor(1,ipoin),ndim,npold,ihost,coold,d1,d2,d3,d4)
             call distint(ndim,neold,npold,nnode,elold,eltoel,coold,coor(1,ipoin),iguess,&
                  ihost,lmark,d1,d2,d3,d4,ierr)
             if(ierr==1)then
                write(*,*)'Error in interpnt point:',ipoin 
                
                return
             endif
             !
             !     Do we finally have a valid host?
             !
             if(ihost==0)then
                !
                !     Outside the domain, interpolate with the nearest neighbor
                !
                ihost=-ipclos
             endif
          endif
       endif
    endif
    !
    !     Do we have a valid host element or we interpolate with the nearest neighbor?
    !
    if(ihost<0)then
       write(*,*)'Error in interpnt, element not found'
       ierr=1_ip
       return
    endif
    !
    !     Store host element
    !
    ielem=ihost

  end subroutine interpnt

end module mod_interp



