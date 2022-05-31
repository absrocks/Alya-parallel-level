module mod_genbou

contains

  subroutine genbou(lface,nface,nnofa,nelem,elem,nnode,coor,ndim,npoin,lfath,nnosi,&
       ltet,eltoel,lelem,rqual)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    use mod_mshtol, only : ptoelm,ptoedg
    use mod_voltol, only : hashface,hashedge
    implicit none
    integer(ip),intent(in)    :: ndim,nface,nnofa,nnode,nnosi
    integer(ip),intent(inout) :: lface(nnofa,nface)
    integer(ip),intent(inout) :: npoin,nelem
    integer(ip),pointer       :: elem(:,:),lfath(:),ltet(:),eltoel(:,:),lelem(:)   
    real(rp),pointer          :: coor(:,:),rqual(:)   
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),ledge(:,:)
    integer(ip),pointer       :: lhashf1(:,:),lhashf2(:)
    integer(ip),pointer       :: lhashe1(:),lhashe2(:)
    integer(ip),pointer       :: lstacke(:),lmarke(:)
    integer(ip)               :: nedge,nhashf,nhashe,iface
    integer(4)                :: istat 
    !
    !     This subroutine regenerates the boundary by inserting the points
    !     with an advancing front approach 
    !

    !
    !     Get the faces surrounding the points 
    !     The surface points are from npoin-npoinf to npoin 
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get  all the  edges
    !
    call ptoedg(lface,nface,npoin,nnode,ptoel1,ptoel2,nedge,ledge)
    !
    !     Hash the faces
    !
    call hashface(nface,nnofa,lface,lhashf1,lhashf2,nhashf)
    !
    !     Hash the edges
    !
    call hashedge(nedge,nnosi,ledge,lhashe1,lhashe2,nhashe)
    !
    !     Insert first triangle
    ! 
    call insertface(lface,nface,nnofa,nelem,elem,nnode,coor,ndim,npoin,lfath,1_ip,&
         ltet,lelem,eltoel,rqual,lstacke,lmarke)
    !
    !     Now regenerate all surface 
    !
    call advgenbou(ndim,nface,nnofa,nnode,nnosi,lface,npoin,nelem,elem,ltet,eltoel,&
       lelem,coor,rqual,ptoel1,ptoel2,1_ip,lstacke,lmarke)

    call memchk(2_ip,istat,memor_msh,'LEDGE','genbou',ledge)
    deallocate(ledge,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGE','genbou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHASHE1','genbou',lhashe1)
    deallocate(lhashe1,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHASHE1','genbou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHASHE2','genbou',lhashe2)
    deallocate(lhashe2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHASHE2','genbou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHASHF1','genbou',lhashf1)
    deallocate(lhashf1,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHASHF1','genbou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHASHF2','genbou',lhashf2)
    deallocate(lhashf2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHASHF2','genbou',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','genbou',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','genbou',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','genbou',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','genbou',0_ip)


  end subroutine genbou

  subroutine insertface(lface,nface,nnofa,nelem,elem,nnode,coor,ndim,npoin,lfath,iface,&
       ltet,lelem,eltoel,rqual,lstack,lmark)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    use mod_voltol, only : findin,split3d,swaploc3dc
    implicit none
    integer(ip),intent(in)    :: ndim,nface,nnofa,nnode,iface
    integer(ip),intent(inout) :: lface(nnofa,nface)
    integer(ip),intent(inout) :: npoin,nelem
    integer(ip),pointer       :: elem(:,:),lfath(:),ltet(:),eltoel(:,:),lelem(:)   
    integer(ip),pointer       :: lstack(:),lmark(:)   
    real(rp),pointer          :: coor(:,:),rqual(:)   
    integer(ip)               :: ip1,ip2,ip3,ipfath,iguess,ielem,idir,isplit,idone
    !
    !     This subroutine inserts the first triangle and regenerate it in 
    !     the tet mesh
    ! 
    ip1=lface(1,iface)
    ip2=lface(2,iface)
    ip3=lface(3,iface)
    !
    !     Get father  
    !
    ipfath=lfath(ip1)
    !
    !     Get initial guess
    !
    iguess=ltet(ipfath)
    !
    !     Find the element containing ip1
    !
    call findin(iguess,elem,nelem,npoin,coor,coor(:,ip1),ndim,ielem,&
         nnode,eltoel,isplit,idir)
    !
    !     Insert the point in the element
    !
    call split3d(ip1,coor,ndim,npoin,nelem,elem,nnode,ielem,isplit,idir,ltet,eltoel)
    !
    !     Optimize the elements
    !
    call swaploc3dc(ielem,elem,nnode,nelem,coor,npoin,ndim,lelem,eltoel,rqual,ltet)
    !
    !     Insert ip2
    !  
    iguess=ltet(ip1)
    !
    !     Find the element containing ip2
    !
    call findin(iguess,elem,nelem,npoin,coor,coor(:,ip2),ndim,ielem,&
         nnode,eltoel,isplit,idir)
    !
    !     Insert the point in the element
    !
    call split3d(ip2,coor,ndim,npoin,nelem,elem,nnode,ielem,isplit,idir,ltet,eltoel)
    !
    !     Optimize the elements
    !
    call swaploc3dc(ielem,elem,nnode,nelem,coor,npoin,ndim,lelem,eltoel,rqual,ltet)
    !
    !     Regenerate (ip1,ip2)
    !
    call genedg(ip1,ip2,elem,nnode,nelem,ltet,eltoel,coor,ndim,npoin,&
                    lstack,lmark,idone)

  end  subroutine insertface

  subroutine genedg(ip1,ip2,elem,nnode,nelem,ltet,eltoel,coor,ndim,npoin,&
                    lstack,lmark,idone)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_mshtol, only :  orient3D
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: nelem,ndim,ip1,ip2,nnode 
    integer(ip),intent(inout) :: npoin,idone
    integer(ip),pointer       :: elem(:,:),eltoel(:,:) 
    integer(ip),pointer       :: lstack(:),lmark(:),ltet(:) 
    real(rp),pointer          :: coor(:,:)  
    integer(ip)               :: nstack,istack,ielem,iview,ineigh 
    integer(ip)               :: ipa,ipb,ipc
    real(rp)                  :: d1,d2,d3,p1(3),p2(3),p3(3),c00,dtot
    integer(ip)           :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))
    !
    !     This subroutine regenerates edge (ip1,ip2)
    !
    c00=0.0d+00 
    !
    !     First get the ball of ip1
    !
    ielem=ltet(ip1)
    istack=0_ip
    nstack=1_ip 
    lmark(ielem)=1_ip
    lstack(1)=ielem

    do 
       if(istack==nstack)exit
       istack=istack+1_ip
       ielem=lstack(istack)
       !
       !     Find ip1
       !
       if(elem(1,ielem)==ip1)then 
          iview=1_ip 
       else if(elem(2,ielem)==ip1)then 
          iview=2_ip 
       else if(elem(3,ielem)==ip1)then 
          iview=3_ip 
       else
          iview=4_ip 
       endif
       !
       !     Get opposite face
       !
       ipa=elem(ltab(1,iview),ielem)  
       ipb=elem(ltab(2,iview),ielem)  
       ipc=elem(ltab(3,iview),ielem)
       !
       !     Does the edge already exist?
       !
       if(ipa==ip2)then
          idone=1_ip
          return
       endif 
       if(ipb==ip2)then
          idone=1_ip
          return
       endif 
       if(ipc==ip2)then
          idone=1_ip
          return
       endif 
       !
       !     Is ip2 on the other side of face (ipa,ipb,ipc)?
       !
       p1(1)=coor(1,ipa)-coor(1,ip2)
       p1(2)=coor(2,ipa)-coor(2,ip2)
       p1(3)=coor(3,ipa)-coor(3,ip2)
       p2(1)=coor(1,ipb)-coor(1,ip2)
       p2(2)=coor(2,ipb)-coor(2,ip2)
       p2(3)=coor(3,ipb)-coor(3,ip2)
       p3(1)=coor(1,ipc)-coor(1,ip2)
       p3(2)=coor(2,ipc)-coor(2,ip2)
       p3(3)=coor(3,ipc)-coor(3,ip2)

       call orient3D(p1,p2,p3,dtot,ndim) 

       if(dtot<c00)then
          !
          !     Does edge (ip1,ip2) intersect face (ipa,ipb,ipc)?
          !
          p1(1)=coor(1,ip2)-coor(1,ip1)
          p1(2)=coor(2,ip2)-coor(2,ip1)
          p1(3)=coor(3,ip2)-coor(3,ip1)
          p2(1)=coor(1,ipa)-coor(1,ip1)
          p2(2)=coor(2,ipa)-coor(2,ip1)
          p2(3)=coor(3,ipa)-coor(3,ip1)
          p3(1)=coor(1,ipb)-coor(1,ip1)
          p3(2)=coor(2,ipb)-coor(2,ip1)
          p3(3)=coor(3,ipb)-coor(3,ip1)
          call orient3D(p1,p2,p3,d1,ndim) 
          p1(1)=coor(1,ip2)-coor(1,ip1)
          p1(2)=coor(2,ip2)-coor(2,ip1)
          p1(3)=coor(3,ip2)-coor(3,ip1)
          p2(1)=coor(1,ipb)-coor(1,ip1)
          p2(2)=coor(2,ipb)-coor(2,ip1)
          p2(3)=coor(3,ipb)-coor(3,ip1)
          p3(1)=coor(1,ipc)-coor(1,ip1)
          p3(2)=coor(2,ipc)-coor(2,ip1)
          p3(3)=coor(3,ipc)-coor(3,ip1)
          call orient3D(p1,p2,p3,d2,ndim)
          p1(1)=coor(1,ip2)-coor(1,ip1)
          p1(2)=coor(2,ip2)-coor(2,ip1)
          p1(3)=coor(3,ip2)-coor(3,ip1)
          p2(1)=coor(1,ipc)-coor(1,ip1)
          p2(2)=coor(2,ipc)-coor(2,ip1)
          p2(3)=coor(3,ipc)-coor(3,ip1)
          p3(1)=coor(1,ipa)-coor(1,ip1)
          p3(2)=coor(2,ipa)-coor(2,ip1)
          p3(3)=coor(3,ipa)-coor(3,ip1)
          call orient3D(p1,p2,p3,d3,ndim) 
          if(d1>c00 .and. d2>c00 .and.d3>c00)then
             exit
          endif  
 
       endif
       !
       !     Add the neigbors to the stack if not marked
       !
       ineigh=eltoel(ltab(1,iview),ielem)
       if(lmark(ineigh)==0)then
          lmark(ineigh)=1_ip
          nstack=nstack+1_ip
          lstack(nstack)=ineigh
       endif 
       ineigh=eltoel(ltab(2,iview),ielem)
       if(lmark(ineigh)==0)then
          lmark(ineigh)=1_ip
          nstack=nstack+1_ip
          lstack(nstack)=ineigh
       endif 
       ineigh=eltoel(ltab(3,iview),ielem)
       if(lmark(ineigh)==0)then
          lmark(ineigh)=1_ip
          nstack=nstack+1_ip
          lstack(nstack)=ineigh
       endif 


    enddo  







  end subroutine genedg

  subroutine genfac(ip1,ip2,ip3)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)  :: ip1,ip2,ip3 


  end subroutine genfac

  subroutine advgenbou(ndim,nface,nnofa,nnode,nnosi,lface,npoin,nelem,elem,ltet,eltoel,&
       lelem,coor,rqual,ptoel1,ptoel2,jface,lstacke,lmarke)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    use mod_mshtol, only : trtotr
    use mod_voltol, only : findin,split3d,swaploc3dc
    implicit none
    integer(ip),intent(in)    :: ndim,nface,nnofa,nnode,nnosi,jface
    integer(ip),intent(inout) :: lface(nnofa,nface)
    integer(ip),intent(inout) :: npoin,nelem
    integer(ip),pointer       :: elem(:,:),lfath(:),ltet(:),eltoel(:,:),lelem(:)   
    real(rp),pointer          :: coor(:,:),rqual(:)   
    integer(ip),pointer       :: ptoel1(:),ptoel2(:)
    integer(ip),pointer       :: lmarke(:),lstacke(:)
    integer(ip),pointer       :: eltoelf(:,:),lmark(:),lstack(:,:),lmarkp(:)
    integer(ip)               :: nstack,istack,inofa,ineigh,iview,ipos,ip1,ip2,ipnew
    integer(ip)               :: iface,ielem,isplit,idir,iguess,idone
    integer(4)                :: istat 
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This subroutine regenerates the boundary with an advancing front approach
    ! 
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Allocate help arrays
    !
    allocate(eltoelf(1,1),stat=istat)
    call runend('advgenbou: USELESS')
    allocate(lstack(2,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','advgenbou',lstack)
    allocate(lmark(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','advgenbou',lmark)
    allocate(lmarkp(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKP','advgenbou',lmarkp)
    !
    !     Iface has been regenerated, march through the surface
    !
    !
    !     Initialize the front
    !
    nstack=1_ip
    istack=1_ip
    lstack(1,1)=iface
    lmark(iface)=1_ip 
    !
    !     Find neighbors
    !
    do inofa=1,nnofa
       ineigh=eltoelf(inofa,iface)
       if(ineigh/=0)then
          if(lmark(ineigh)==0)then
             !
             !     Find viewer
             !
             if(eltoelf(1,ineigh)==iface)then
                iview=1_ip 
             else if(eltoelf(2,ineigh)==iface)then
                iview=2_ip 
             else
                iview=3_ip
             endif
             nstack=nstack+1_ip
             lstack(1,istack)=ineigh
             lstack(2,istack)=iview
          endif
       endif
    enddo
    !
    !     Loop on stack
    !
    do 
       if(istack==nstack)exit 
       istack=istack+1_ip
       iface=lstack(1,istack)
       lmark(iface)=1_ip
       ipos=lstack(2,istack)
       !
       !     Get edge end points
       !
       ip1=lface(ltab(1,ipos),iface)
       ip2=lface(ltab(2,ipos),iface)
       !
       !     Get the new point
       ! 
       ipnew=lface(ipos,iface)
       !
       !     Do we have to insert ipnew
       !
       if(lmarkp(ipnew)==0)then
          !
          !     Get initial guess
          !
          iguess=ltet(ip1)
          !
          !     Find the element containing ip1
          !
          call findin(iguess,elem,nelem,npoin,coor,coor(:,ip1),ndim,ielem,&
               nnode,eltoel,isplit,idir)
          !
          !     Insert the point in the element
          !
          call split3d(ipnew,coor,ndim,npoin,nelem,elem,nnode,ielem,isplit,idir,ltet,eltoel)
          !
          !     Optimize the elements
          !
          call swaploc3dc(ielem,elem,nnode,nelem,coor,npoin,ndim,lelem,eltoel,rqual,ltet)

          lmarkp(ipnew)=1_ip

       endif
       !
       !     Regenerate edge (ip1,ipnew)        
       !
       call genedg(ip1,ipnew,elem,nnode,nelem,ltet,eltoel,coor,ndim,npoin,&
                    lstacke,lmark,idone) 
       !
       !     Regenerate edge (ip2,ipnew)        
       !
       call genedg(ip2,ipnew,elem,nnode,nelem,ltet,eltoel,coor,ndim,npoin,&
                    lstacke,lmark,idone) 
       !
       !     Regenerate face (ip2,ipnew)        
       !
       call genfac(ip1,ip2,ipnew) 
       !
       !     Update the stack
       !
       ineigh=eltoelf(ltab(1,ipos),iface)
       if(ineigh/=0)then
          if(lmark(ineigh)==0)then
             !
             !     Find viewer
             !
             if(eltoelf(1,ineigh)==iface)then
                iview=1_ip 
             else if(eltoelf(2,ineigh)==iface)then
                iview=2_ip 
             else
                iview=3_ip
             endif
             nstack=nstack+1_ip
             lstack(1,istack)=ineigh
             lstack(2,istack)=iview
          endif
       endif

       ineigh=eltoelf(ltab(2,ipos),iface)
       if(ineigh/=0)then
          if(lmark(ineigh)==0)then
             !
             !     Find viewer
             !
             if(eltoelf(1,ineigh)==iface)then
                iview=1_ip 
             else if(eltoelf(2,ineigh)==iface)then
                iview=2_ip 
             else
                iview=3_ip
             endif
             nstack=nstack+1_ip
             lstack(1,istack)=ineigh
             lstack(2,istack)=iview
          endif
       endif

    enddo

    call memchk(2_ip,istat,memor_msh,'ELTOEL','advgenbou',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','advgenbou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','advgenbou',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','advgenbou',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','advgenbou',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','advgenbou',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','advgenbou',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','advgenbou',0_ip)

  end subroutine advgenbou

end module mod_genbou

