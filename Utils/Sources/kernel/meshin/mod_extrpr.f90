module mod_extrpr
  use mod_mshtol

contains


  subroutine blmesh(nface,nnofa,nnode,ndim,nelem,npoin,nblay,npsur,&
       thick,reason,rsize,lface,coor,elem,rblay,lsurf)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh 
    implicit none
    integer(ip), intent(in)   :: nnofa,nnode,ndim,nblay,npsur
    integer(ip), intent(inout):: nface 
    integer(ip), intent(inout):: npoin,nelem
    real(rp), intent(in)      :: rsize,rblay(nblay) 
    real(rp),intent(in)       :: thick,reason 
    real(rp),pointer          :: coor(:,:)
    integer(ip),pointer       :: elem(:,:),lface(:,:),lsurf(:)
    integer(ip),pointer       :: lmark(:),ptoel1(:),ptoel2(:),eltoel(:,:)
    integer(ip),pointer       :: ldiag(:),ledge(:,:),lftoed(:,:),lplay(:)
    integer(ip),pointer       :: ltest(:),lefac(:)
    real(rp),pointer          :: rnopo(:,:),rnofa(:,:)
    integer(ip)               :: nedge,iface,ielem,ipoin
    integer(4)                :: istat
    !
    !     Allocate
    !
    allocate(rnopo(ndim,npsur),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOPO','blmesh',rnopo)
    allocate(rnofa(ndim,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOFA','blmesh',rnofa)
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface, nface , npsur, nnofa,ptoel1,ptoel2 )
    !
    !     Get the face normals
    !
    call gtfnrl(lface,nface,nnofa,ndim,coor,npsur,rnofa)
    !
    !     Get the point normals
    !
    call gtpnrl(nface,rnopo,npsur,ndim,ptoel1,ptoel2,rnofa) 
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npsur,eltoel)
    !
    !     Get the faces surrounding edges and edges
    !
    call fatoed(nedge,nface,lface,nnofa,eltoel,ledge,lftoed)    
    !
    !     Get the orientation of the diagonal
    !
    allocate(ldiag(nedge),stat=istat)
    call memchk(zero,istat,memor_msh,'LDIAG','blmesh',ldiag)

    call diagdi(ldiag,npsur,lface,nface,nnofa,lftoed,nedge,ptoel1,ptoel2,ledge)
    !
    !     Extrude the prisms
    !
    call extrpr(lface,nface,nnofa,nelem,npoin,nblay,nnode,lftoed,nedge,&
         thick,reason,rnopo,npsur,ndim,elem,coor,rblay,lplay,lefac,eltoel)
    !
    !     Allocate lmark
    !
    allocate(lmark(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','blmesh',lmark)
    allocate(ltest(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LTEST','blmesh',ltest)
    !
    !     Set lmark to 1
    !
    do ielem=1,nelem
       lmark(ielem)=1_ip 
    enddo
    !
    !     Delete bad elements
    !
    !      call dlblmsh(coor,ndim,npoin,elem,nnode,nelem,rsize )
    !
    !     Take out one more level
    !
    call onemore(nnode,nelem,elem,lmark,lplay,npoin)
    !
    !     Output the mesh
    !   
    !call outaniface3(nnofa,nfnew,npoin,ndim,lfnew,coor,lfmark)
    !
    !     Allow for one level jump only
    !
    call rmvjmp(elem,nelem,npoin,nnode,lplay,lmark,nface,nblay)
    !
    !     Mark the bottom faces to 1 and invert them
    !
    do iface=1,nface
       lsurf(iface)=1_ip
       ipoin=lface(1,iface)
       lface(1,iface)=lface(2,iface)
       lface(2,iface)=ipoin
    enddo
    !
    !     Add elements to make a smooth transition
    !
    call addcover(elem,nelem,nnode,lmark,lplay,npoin,nnofa,nface,lface,&
         ltest,ndim,coor,eltoel,lsurf,lefac)
    !
    !     Compact elements and points
    !
    call compac(nelem,elem,coor,nnode,npoin,lmark,ndim,lface,nnofa,nface)
    !
    !     DBG
    !
    call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)
    call chkface(nelem,elem,nnode,nface,nnofa,lface,npoin)

    call memchk(2_ip,istat,memor_msh,'LTEST','blmesh',ltest)
    deallocate(ltest,stat=istat)
    if(istat/=0) call memerr(2_ip,'LTEST','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPLAY','blmesh',lplay)
    deallocate(lplay,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPLAY','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','blmesh',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','blmesh',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','blmesh',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','blmesh',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGE','blmesh',ledge)
    deallocate(ledge,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGE','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFTOED','blmesh',lftoed)
    deallocate(lftoed,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFTOED','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LDIAG','blmesh',ldiag)
    deallocate(ldiag,stat=istat)
    if(istat/=0) call memerr(2_ip,'LDIAG','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOFA','blmesh',rnofa)
    deallocate(rnofa,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOFA','blmesh',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOPO','blmesh',rnopo)
    deallocate(rnopo,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOPO','blmesh',0_ip)

  end subroutine blmesh

  subroutine diagdi(ldiag,npoin,lface,nface,nnofa,lftoed,nedge,ptoel1,ptoel2,ledge)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)         :: npoin,nface,nnofa,nedge 
    integer(ip),intent(inout)      :: ldiag(nedge),lftoed(nnofa,nface)
    integer(ip),intent(in)         :: lface(nnofa,nface),ledge(2,nedge)
    integer(ip),intent(in)         :: ptoel1(*),ptoel2(npoin+1) 
    integer(ip), pointer           :: lpoin(:)                    
    integer(ip)                :: ipoin,ifa,iface,ip1,iedge,inode,icont  
    integer(4)                 :: istat
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    allocate(lpoin(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOIN','diagdi',lpoin)

    do ipoin=1,npoin
       if(lpoin(ipoin)==0_ip)then

          lpoin(ipoin)=1_ip
          !
          !     Loop on surrounding faces
          !    
          do ifa=ptoel2(ipoin),ptoel2(ipoin+1)-1

             iface=ptoel1(ifa)
             do inode=1,nnofa
                ip1=lface(inode,iface)
                if(ip1==ipoin)then
                   !
                   !     First edge to visit in this element
                   !
                   iedge=lftoed(ltab(2,inode),iface) 
                   !
                   !     Has the edge been marked before
                   !
                   if(ldiag(iedge)==0)then
                      !
                      !     Mark the edge
                      !
                      if(ledge(1,iedge)==ipoin)then
                         ldiag(iedge)=-1
                      else
                         ldiag(iedge)=1
                      endif
                   endif
                   !
                   !     Second edge to visit in this element
                   !
                   iedge=lftoed(ltab(1,inode),iface) 
                   !
                   !     Has the edge been marked before
                   !
                   if(ldiag(iedge)==0)then
                      !
                      !     Mark the edge
                      !
                      if(ledge(1,iedge)==ipoin)then
                         ldiag(iedge)=-1
                      else
                         ldiag(iedge)=1
                      endif
                   endif
                endif
             enddo
          enddo
       endif
    enddo
    !
    !     Transfer to lftoed
    !
    do iface=1,nface
       do inode=1,nnofa
          ipoin=lface(ltab(1,inode),iface)
          iedge=lftoed(inode,iface) 
          if(ipoin==ledge(1,iedge))then
             lftoed(inode,iface)=ldiag(iedge)
          else
             lftoed(inode,iface)=-ldiag(iedge)
          endif
       enddo
    enddo
    !
    !     DBG
    !
    do iface=1,nface
       icont=lftoed(1,iface)+lftoed(2,iface)+lftoed(3,iface)
       if(abs(icont)==3)then
          write(*,*)'Error in diagdi, diagonal not correct'
          stop
       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'LPOIN','diagdi',lpoin)
    deallocate(lpoin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPOIN','diagdi',0_ip)

  end subroutine diagdi

  subroutine extrpr(lface,nface,nnofa,nelem,npoin,nblay,nnode,lftoed,nedge,thick,reason,&
       rnopo,npsur,ndim,elem,coor,rblay,lplay,lefac,eltoel)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)    :: nface,nnofa,nnode,nedge,nblay,npsur,ndim
    integer(ip), intent(inout) :: nelem,npoin  
    integer(ip), intent(inout) :: eltoel(nnofa,nface)  
    integer(ip), intent(in)    :: lface(nnofa,nface),lftoed(nnofa,nedge)
    integer(ip)                :: iface,inode,npnew,isid1,isid2,iba,ibb,ibc,ita,itb,itc,ilay,ipoin
    integer(ip),pointer        :: lpoin(:)
    integer(ip),pointer        :: elem(:,:),lplay(:),lefac(:)      
    real(rp),pointer           :: coor(:,:) 
    real(rp),intent(in)        :: rblay(nblay)  
    integer(4)                 :: istat
    integer(ip)                :: nenew,ineig1,ineig2,ineig3 
    integer(ip)                :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    real(rp)                   :: thick,reason,thickness,coef,x0,y0,z0,reasont 
    real(rp),intent(in)        :: rnopo(ndim,npsur)  
    !
    !     Allocate the elements
    !
    nenew=nface*3*nblay
    allocate(elem(nnode,nenew),stat=istat)
    call memchk(zero,istat,memor_msh,'ELEM','extrpr',elem)
    allocate(lefac(nenew),stat=istat)
    call memchk(zero,istat,memor_msh,'ELEM','extrpr',lefac)
    !
    !     Allocate help arrays
    !
    allocate(lpoin(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOIN','extrpr',lpoin)
    !
    !     Count the number of new points
    !
    npnew=npoin
    do iface=1,nface
       do inode=1,nnofa
          ipoin=lface(inode,iface)
          if(lpoin(ipoin)==0)then
             npnew=npnew+nblay
             lpoin(ipoin)=1
          endif
       enddo
    enddo
    !
    !     Reallocate lpoin && coor
    ! 
    call memrea(npnew,memor_msh,'LPOIN','extrpr',lpoin)
    do ipoin=1,npnew
       lpoin(ipoin)=0_ip  
    enddo

    call memrea(npnew,memor_msh,'COOR','extrpr',coor)
    npoin=npsur
    !
    !     Allocate lplay
    !
    allocate(lplay(npnew),stat=istat)
    call memchk(zero,istat,memor_msh,'LPLAY','extrpr',lplay)
    !
    !     Generate the new points
    !
    !     Loop on the faces
    !
    do iface=1,nface

       do inode=1,nnofa

          ipoin=lface(inode,iface)

          if(lpoin(ipoin)==0)then
             !
             !     Build the point position
             !                  
             x0=coor(1,ipoin)
             y0=coor(2,ipoin)
             z0=coor(3,ipoin)

             lpoin(ipoin)=npoin+1
             lplay(ipoin)=0_ip
             reasont=1.0d+00

             !                  do ilay=1,nblay

             !                      coef=reasont
             !                      thickness=coef*thick
             !                      npoin=npoin+1  
             !                      coor(1,npoin)=x0+rnopo(1,ipoin)*thickness             
             !                      coor(2,npoin)=y0+rnopo(2,ipoin)*thickness             
             !                      coor(3,npoin)=z0+rnopo(3,ipoin)*thickness             

             !                      x0=coor(1,npoin)
             !                      y0=coor(2,npoin)
             !                      z0=coor(3,npoin)
             !                      reasont=reasont*reason 
             !                  enddo

             do ilay=1,nblay

                thickness=rblay(ilay)
                npoin=npoin+1  
                coor(1,npoin)=x0+rnopo(1,ipoin)*thickness             
                coor(2,npoin)=y0+rnopo(2,ipoin)*thickness             
                coor(3,npoin)=z0+rnopo(3,ipoin)*thickness             
                lplay(npoin)=ilay
                lpoin(npoin)=ipoin 

                x0=coor(1,npoin)
                y0=coor(2,npoin)
                z0=coor(3,npoin)

             enddo


          endif
       enddo
    enddo
    !
    !     Generate the new elements
    !
    nelem=0_ip
    do iface=1,nface
       !
       !     Find the alternate diagonal configuration
       !
       loop_side:do inode=1,nnofa
          isid1=inode
          isid2=ltab(1,inode)  
          if(lftoed(isid1,iface)==-1 .and. lftoed(isid2,iface)==1)then
             exit loop_side             
          endif
       enddo loop_side

       iba=lface(ltab(1,inode),iface)
       ibb=lface(ltab(2,inode),iface)
       ibc=lface(inode,iface)

       ineig1=eltoel(ltab(1,inode),iface)
       ineig2=eltoel(ltab(2,inode),iface)
       ineig3=eltoel(inode,iface)
       !
       !     Modify eltoel
       !
       eltoel(1,iface)=ineig1
       eltoel(2,iface)=ineig2
       eltoel(3,iface)=ineig3

       ita=lpoin(iba)
       itb=lpoin(ibb)
       itc=lpoin(ibc)
       !
       !     Build the two kinds of prisms
       !
       if(lftoed(ltab(2,inode),iface)==1)then
          !
          !     -++ prism
          ! 
          do ilay=1,nblay

             nelem=nelem+1
             elem(1,nelem)=iba
             elem(2,nelem)=ibb
             elem(3,nelem)=ibc
             elem(4,nelem)=ita
             lefac(nelem)=iface 

             nelem=nelem+1
             elem(1,nelem)=ita
             elem(2,nelem)=ibb
             elem(3,nelem)=ibc
             elem(4,nelem)=itc
             lefac(nelem)=iface 

             nelem=nelem+1
             elem(1,nelem)=ita
             elem(2,nelem)=ibb
             elem(3,nelem)=itc
             elem(4,nelem)=itb
             lefac(nelem)=iface 

             iba=ita
             ibb=itb
             ibc=itc

             ita=ita+1
             itb=itb+1
             itc=itc+1

          enddo

       else
          !
          !     -+- prism
          ! 
          do ilay=1,nblay

             nelem=nelem+1
             elem(1,nelem)=iba
             elem(2,nelem)=ibb
             elem(3,nelem)=ibc
             elem(4,nelem)=itc
             lefac(nelem)=iface 

             nelem=nelem+1
             elem(1,nelem)=iba
             elem(2,nelem)=ibb
             elem(3,nelem)=itc
             elem(4,nelem)=ita
             lefac(nelem)=iface 

             nelem=nelem+1
             elem(1,nelem)=ita
             elem(2,nelem)=ibb
             elem(3,nelem)=itc
             elem(4,nelem)=itb
             lefac(nelem)=iface 

             iba=ita
             ibb=itb
             ibc=itc

             ita=ita+1_ip
             itb=itb+1_ip
             itc=itc+1_ip

          enddo

       endif

    enddo
    !
    !     End loop on faces
    !
    call memchk(2_ip,istat,memor_msh,'LPOIN','extrpr',lpoin)
    deallocate(lpoin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPOIN','extrpr',0_ip)

  end subroutine extrpr

  subroutine dlblmsh(coor,ndim,npoin,elem,nnode,nelem,rsize )
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnode,ndim
    real(rp),intent(in)        :: rsize
    integer(ip),intent(inout)  :: nelem,npoin,elem(nnode,nelem)
    integer(ip),pointer        :: lelem(:)
    real(rp),pointer           :: rvol(:) 
    real(rp)                   :: coor(ndim,npoin) 
    integer(4)                 :: istat

    allocate(lelem(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEM','dlblmsh',lelem)
    allocate(rvol(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'RVOL','dlblmsh',rvol)

    !
    !     Check volume of elements
    !
    call chkvol(nelem,elem,coor,nnode,npoin,lelem,ndim,rvol)
    !
    !     Check aspect ratio of element
    !
    call chkasp(nelem,elem,coor,nnode,npoin,lelem,ndim)
    !
    !     Check size of element
    !
    call chksiz(nelem,elem,coor,nnode,npoin,lelem,rvol,ndim,rsize)
    !
    !     Check intersection
    !


    call memchk(2_ip,istat,memor_msh,'LELEM','dlblmsh',lelem)
    deallocate(lelem,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEM','dlblmsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'RVOL','dlblmsh',rvol)
    deallocate(rvol,stat=istat)
    if(istat/=0) call memerr(2_ip,'RVOL','dlblmsh',0_ip)

  end subroutine dlblmsh

  subroutine chkvol(nelem,elem,coor,nnode,npoin,lelem,ndim,rvol)
    use mod_memchk
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)  :: nelem,npoin,nnode,ndim
    integer(ip),intent(in)   :: elem(nnode,nelem)
    real(rp),intent(in)      :: coor(ndim,npoin)
    integer(ip),intent(inout):: lelem(nelem)
    real(rp),intent(inout)   :: rvol(nelem)
    integer(ip)              :: ielem,ip1,ip2,ip3,ip4
    real(rp)                 :: x1,y1,z1,dx21,dy21,dz21,dx31,dy31,dz31,dx41,dy41,dz41
    !
    !     Compute volume
    !  
    do ielem=1,nelem
       ip1=elem(1,ielem) 
       ip2=elem(2,ielem) 
       ip3=elem(3,ielem) 
       ip4=elem(4,ielem) 

       x1  =coor(1,ip1)
       y1  =coor(2,ip1)
       z1  =coor(3,ip1)
       dx21=coor(1,ip2)-x1
       dy21=coor(2,ip2)-y1
       dz21=coor(3,ip2)-z1
       dx31=coor(1,ip3)-x1
       dy31=coor(2,ip3)-y1
       dz31=coor(3,ip3)-z1
       dx41=coor(1,ip4)-x1
       dy41=coor(2,ip4)-y1
       dz41=coor(3,ip4)-z1

       rvol(ielem)=dx41*(dy21*dz31-dz21*dy31) &
            +dy41*(dz21*dx31-dx21*dz31)   &
            +dz41*(dx21*dy31-dy21*dx31)

    enddo
    !
    !     Check if positive
    !
    do ielem=1,nelem
       if(rvol(ielem)<zero)then
          lelem(ielem)=0
       endif
    enddo

  end subroutine chkvol

  subroutine chkasp(nelem,elem,coor,nnode,npoin,lelem,ndim)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)  :: nelem,npoin,nnode,ndim
    integer(ip),intent(in)   :: elem(nnode,nelem)
    integer(ip),intent(inout):: lelem(nelem)
    integer(ip)              :: ielem,ip1,ip2,ip3,ip4,jelem
    integer(ip)              :: iba,ibb,ibc,ita,itb,itc
    real(rp)                 :: dx,dy,dz,dl1,dl2,dl3,dzmax,dh1,dh2,dh3,dhmin,TOL 
    real(rp),intent(in)      :: coor(ndim,npoin)
    !
    !     Set tolerance
    !
    TOL=1.0d+01

    do ielem=1,nelem,3

       jelem=ielem+2

       iba=elem(1,ielem)
       ibb=elem(2,ielem)
       ibc=elem(3,ielem)

       ita=elem(1,jelem)
       itb=elem(2,jelem)
       itc=elem(3,jelem)


       dx=coor(1,iba)-coor(1,ita) 
       dy=coor(2,iba)-coor(2,ita) 
       dz=coor(3,iba)-coor(3,ita) 
       dl1=sqrt(dx*dx+dy*dy+dz*dz)

       dx=coor(1,ibb)-coor(1,itb) 
       dy=coor(2,ibb)-coor(2,itb) 
       dz=coor(3,ibb)-coor(3,itb) 
       dl2=sqrt(dx*dx+dy*dy+dz*dz)

       dx=coor(1,ibc)-coor(1,itc) 
       dy=coor(2,ibc)-coor(2,itc) 
       dz=coor(3,ibc)-coor(3,itc) 
       dl3=sqrt(dx*dx+dy*dy+dz*dz)

       dzmax=max(dl1,dl2,dl3)

       dx=coor(1,iba)-coor(1,ibc) 
       dy=coor(2,iba)-coor(2,ibc) 
       dz=coor(3,iba)-coor(3,ibc) 
       dh1=sqrt(dx*dx+dy*dy+dz*dz)

       dx=coor(1,iba)-coor(1,ibb) 
       dy=coor(2,iba)-coor(2,ibb) 
       dz=coor(3,iba)-coor(3,ibb) 
       dh2=sqrt(dx*dx+dy*dy+dz*dz)

       dh3=sqrt(dh1*dh1+dh2*dh2)

       dhmin=min(dh1,dh2,dh3)


       if(dhmin/dzmax<TOL)then
          lelem(ielem)=0
          lelem(ielem+1)=0
          lelem(ielem+2)=0
       endif

    enddo


  end subroutine chkasp

  subroutine chksiz(nelem,elem,coor,nnode,npoin,lelem,rvol,ndim,rsize)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)  :: nelem,npoin,nnode,ndim
    integer(ip),intent(in)   :: elem(nnode,nelem)
    integer(ip),intent(inout):: lelem(nelem)
    real(rp),intent(in)      :: coor(ndim,npoin),rvol(nelem),rsize
    integer(ip)              :: ielem,ip1,ip2,ip3,ip4
    integer(4)               :: istat
    real(rp)                 :: rexp,rsiz 
    real(rp),pointer         :: rpexp(:)
    !
    !     Allocate array for point size
    !
    allocate(rpexp(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'RPEXP','chksiz',rpexp)
    !
    !     Find the expected size for each element 
    !  
    call expsiz(rpexp,npoin,rsize)
    !
    !     Compare size to expected size
    !
    do ielem=1,nelem

       ip1=elem(1,ielem)
       ip2=elem(2,ielem)
       ip3=elem(3,ielem)
       ip4=elem(4,ielem)

       rexp=(rpexp(ip1)+rpexp(ip2)+rpexp(ip3)+rpexp(ip4))*0.25d+00
       rsiz=rvol(ielem)**0.33333
       if(rsiz>rexp)then
          lelem(ielem)=0
       endif
    enddo


    call memchk(2_ip,istat,memor_msh,'RPEXP','chksiz',rpexp)
    deallocate(rpexp,stat=istat)
    if(istat/=0) call memerr(2_ip,'RPEXP','chksiz',0_ip)

  end subroutine chksiz


  subroutine expsiz(rpexp,npoin,rsize)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: npoin
    real(rp),intent(inout)     :: rpexp(npoin)
    integer(ip)                :: ipoin
    real(rp),intent(in)        :: rsize
    !
    !     For the moment, very simple routine
    !
    do ipoin=1,npoin
       rpexp(ipoin)=rsize
    enddo

  end subroutine expsiz

  subroutine compac(nelem,elem,coor,nnode,npoin,lelem,ndim,lface,nnofa,nface)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnode,ndim,nnofa,nface
    integer(ip),intent(inout)  :: npoin,nelem
    real(rp),intent(inout)     :: coor(ndim,npoin)
    integer(ip),intent(inout)  :: elem(nnode,nelem),lelem(nelem),lface(nnofa,nface)
    integer(ip)                :: ielem,ipoin,nele0,npoi0,iface
    integer(ip),pointer        :: lmark(:)
    integer(4)                 :: istat

    nele0=nelem
    nelem=0
    !
    !     Compact the elements
    !
    do ielem=1,nele0

       if(lelem(ielem)==1)then
          nelem=nelem+1 
          elem(1,nelem)=elem(1,ielem)          
          elem(2,nelem)=elem(2,ielem)          
          elem(3,nelem)=elem(3,ielem)          
          elem(4,nelem)=elem(4,ielem)          
       endif

    enddo
    !
    !     Allocate lmark
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','seltoel',lmark)

    do ielem=1,nelem
       lmark(elem(1,ielem))=1_ip
       lmark(elem(2,ielem))=1_ip
       lmark(elem(3,ielem))=1_ip
       lmark(elem(4,ielem))=1_ip
    enddo
    !
    !     Compress the points
    !
    npoi0=npoin
    npoin=0

    do ipoin=1,npoi0 
       if(lmark(ipoin)==1_ip)then
          npoin=npoin+1
          coor(1,npoin)=coor(1,ipoin)
          coor(2,npoin)=coor(2,ipoin)
          coor(3,npoin)=coor(3,ipoin)
          lmark(ipoin)=npoin
       endif
    enddo
    !
    !     Renumber elements
    !
    do ielem=1,nelem
       elem(1,ielem)=lmark(elem(1,ielem))
       elem(2,ielem)=lmark(elem(2,ielem))
       elem(3,ielem)=lmark(elem(3,ielem))
       elem(4,ielem)=lmark(elem(4,ielem))
    enddo
    !
    !     Renumber faces
    !
    do iface=1,nface
       lface(1,iface)=lmark(lface(1,iface))
       lface(2,iface)=lmark(lface(2,iface))
       lface(3,iface)=lmark(lface(3,iface))
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','trtotr',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','trtotr',0_ip)

  end subroutine compac

  subroutine readbl(nblay,rblay)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(inout)  :: nblay
    integer(ip)                :: iblay,iostatus
    real(rp)                   :: rval
    real(rp),pointer           :: rblay(:)
    integer(4)                 :: istat

    open(unit=77,file='nblay.dat',status='old',iostat=iostatus)
    if(iostatus>0)then
       write(*,*)'Boundary file not found in readbl'
       write(*,*)'Set nblay to zero'
       nblay=0
       allocate(rblay(1),stat=istat)
       call memchk(zero,istat,memor_msh,'RBLAY','readbl',rblay)
       return
    endif
    write(*,*)'Reading nblay.dat'
    read(77,*)nblay
    write(*,*)'nblay=',nblay
    !
    !      Allocate rblay
    !
    allocate(rblay(nblay),stat=istat)
    call memchk(zero,istat,memor_msh,'RBLAY','readbl',rblay)


    do iblay=1,nblay
       read(77,*)rval
       rblay(iblay)=rval
    enddo

    nblay=nblay-1

    !
    !     Compute the delta
    !
    do iblay=1,nblay
       rblay(iblay)=rblay(iblay+1)-rblay(iblay)
       !write(*,*)iblay,rblay(iblay)
    enddo

    close(77)
    write(*,*)'End reading nblay.dat'

  end subroutine readbl

  subroutine addcover(elem,nelem,nnode,lmark,lplay,npoin,nnofa,nface,lface,&
       ltest,ndim,coor,eltoel,lsurf,lefac)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(inout) :: nelem,nface
    integer(ip),intent(in)    :: npoin,nnode,nnofa,ndim
    integer(ip),pointer       :: elem(:,:),lmark(:),lface(:,:),ltest(:)   
    integer(ip),pointer       :: lpoin(:),ptoel1(:),ptoel2(:),lsurf(:) 
    integer(ip)               :: eltoel(nnofa,nface),lefac(nelem)  
    real(rp)                  :: coor(ndim,npoin) 
    integer(ip)               :: lplay(npoin)
    integer(ip)               :: ip1,ip2,ip3,ip4,ielem,nele0,jelem,ielay    
    integer(ip)               :: ip1t,ip2t,ip3t,ip1h,ip2h,ip3h,ip1b,jp1b,jp2b,jp3b,jelay,ldiag(3)   
    integer(ip)               :: icont,iph,jp1t,jp2t,jp3t,ip2b,ip3b,ipmax,iface1,jface1
    integer(ip)               :: kelem,imin,imax,lebas(3),letop(3),ie2,ie3
    integer(ip)               :: ielem2,ielem3,icase,ip1a,ip2a,ip3a,ip4a,ip4b,nfnew
    integer(ip)               :: iface
    real(rp)                  :: x1,y1,z1,dx21,dy21,dz21,dx31,dy31,dz31,dx41,dy41,dz41,v1,v2,c00 
    integer(4)                :: istat
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    integer(ip),parameter     :: ncase=16
    integer(ip)               :: lside(3)=(/2,3,1/),ldiagno(ncase)
    c00=0.0d+00
    !
    !     This sub adds elements to avoid bad aspect faces and update the face list
    !     The initial faces have lsurf(iface)=1
    !     The lateral faces have lsurf(iface)=2
    !     The top faces have lsurf(iface)=3
    !    
    !  
    allocate(lpoin(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOIN','addcover',lpoin) 
    !
    !     Set up the known diagonals
    !
    ldiag(1)=-1_ip
    ldiag(3)= 1_ip
    !
    !     Mark all the points touched by non deleted elements
    !
    do ielem=1,nelem
       if(lmark(ielem)/=0)then
          ip1=elem(1,ielem) 
          ip2=elem(2,ielem) 
          ip3=elem(3,ielem) 
          ip4=elem(4,ielem) 
          lpoin(ip1)=1_ip 
          lpoin(ip2)=1_ip 
          lpoin(ip3)=1_ip 
          lpoin(ip4)=1_ip 
       endif
    enddo
    !
    !     Loop on elements
    !
    nele0=nelem 
    ielem=-2_ip

    do 
       !
       !     Find the next element marked for deletion
       !
       ielem=ielem+3_ip
       !
       !     Did we exhaust the list?
       !
       if(ielem>nele0)exit
       !
       !     Bottom face
       !  
       ip1b=elem(1,ielem)
       ip2b=elem(2,ielem)
       ip3b=elem(3,ielem)
       ielay=lplay(ip1b)
       !
       !     Find the points of the top face
       !
       kelem=ielem+2
       ip1t=elem(1,kelem)
       ip2t=elem(4,kelem)
       ip3t=elem(3,kelem)
       !
       !     Remember the type of the prism
       ! 
       if(elem(4,ielem).eq.ip1t)then
          ldiag(2)=1_ip
       else
          ldiag(2)=-1_ip
       endif
       !
       !     Is the element still active?
       !
       if(lmark(ielem)==1)then
          !
          !    Get the face the element comes from
          !
          iface=lefac(ielem) 
          !
          !     Check if we are on the edge of the surface
          !
          if(eltoel(1,iface)==0)then

             nfnew=nface+2_ip
             call memrea(nfnew,memor_msh,'LFACE','addcover',lface)
             call memrea(nfnew,memor_msh,'LSURF','addcover',lsurf)

             nface=nface+1_ip
             lface(1,nface)=ip2b            
             lface(2,nface)=ip3b            
             lface(3,nface)=ip3t
             lsurf(nface)=2_ip            
             nface=nface+1_ip
             lface(1,nface)=ip2b            
             lface(2,nface)=ip3t            
             lface(3,nface)=ip2t            
             lsurf(nface)=2_ip           

          endif

          if(eltoel(3,iface)==0)then

             nfnew=nface+2_ip
             call memrea(nfnew,memor_msh,'LFACE','addcover',lface)
             call memrea(nfnew,memor_msh,'LSURF','addcover',lsurf)

             nface=nface+1_ip
             lface(1,nface)=ip1b            
             lface(2,nface)=ip2b            
             lface(3,nface)=ip1t
             lsurf(nface)=2_ip            
             nface=nface+1_ip
             lface(1,nface)=ip1t            
             lface(2,nface)=ip2b            
             lface(3,nface)=ip2t            
             lsurf(nface)=2_ip           

          endif

          if(eltoel(2,iface)==0)then

             nfnew=nface+2_ip
             call memrea(nfnew,memor_msh,'LFACE','addcover',lface)
             call memrea(nfnew,memor_msh,'LSURF','addcover',lsurf)

             if(ldiag(2)==1)then 

                nface=nface+1_ip
                lface(1,nface)=ip3b            
                lface(2,nface)=ip1b            
                lface(3,nface)=ip1t
                lsurf(nface)=2_ip            
                nface=nface+1_ip
                lface(1,nface)=ip3b            
                lface(2,nface)=ip1t           
                lface(3,nface)=ip3t            
                lsurf(nface)=2_ip

             else

                nface=nface+1_ip
                lface(1,nface)=ip3b            
                lface(2,nface)=ip1b            
                lface(3,nface)=ip3t
                lsurf(nface)=2_ip            
                nface=nface+1_ip
                lface(1,nface)=ip3t            
                lface(2,nface)=ip1b            
                lface(3,nface)=ip1t            
                lsurf(nface)=2_ip

             endif

          endif
          !
          !     Did we reach the top of a column
          !
          kelem=ielem+3_ip
          if(kelem<=nele0)then
             !
             !     Get level
             !          
             jp1b=elem(1,kelem)
             jelay=lplay(jp1b)  
             !
             !     Compare levels in the column
             ! 
             if(jelay<ielay)then

                nfnew=nface+1_ip
                call memrea(nfnew,memor_msh,'LFACE','addcover',lface)
                call memrea(nfnew,memor_msh,'LSURF','addcover',lsurf)

                kelem=ielem+2
                nface=nface+1_ip
                lface(1,nface)=elem(1,kelem)            
                lface(2,nface)=elem(4,kelem)            
                lface(3,nface)=elem(3,kelem)            
                lsurf(nface)=3_ip    
 
             endif

          else
                !
                !     We are in the last prism
                !
                nfnew=nface+1_ip
                call memrea(nfnew,memor_msh,'LFACE','addcover',lface)
                call memrea(nfnew,memor_msh,'LSURF','addcover',lsurf)

                kelem=ielem+2
                nface=nface+1_ip
                lface(1,nface)=elem(1,kelem)            
                lface(2,nface)=elem(4,kelem)            
                lface(3,nface)=elem(3,kelem)            
                lsurf(nface)=3_ip 

          endif

          cycle
       endif
       !
       !     Initialize the levels
       ! 
       ip1h=0_ip
       ip2h=0_ip
       ip3h=0_ip
       !
       !     Are these points still active?
       !
       if(lpoin(ip1t)==1)ip1h=ip1h+1_ip
       if(lpoin(ip2t)==1)ip2h=ip2h+1_ip
       if(lpoin(ip3t)==1)ip3h=ip3h+1_ip
       !
       !     Go upwards in the column
       !
       jelem=ielem

       do 
          !
          !     Update counter
          !
          jelem=jelem+3 
          !
          !     Did we exhaust the list?
          !
          if(jelem>nele0)exit
          !
          !     Bottom face
          !
          jp1b=elem(1,jelem)
          jp2b=elem(2,jelem)
          jp3b=elem(3,jelem)
          jelay=lplay(jp1b)
          !
          !     Did we exhaust the column?
          !
          if(jelay<=ielay)exit
          !
          !     Top face
          !
          kelem=jelem+2
          jp1t=elem(1,kelem)
          jp2t=elem(4,kelem)
          jp3t=elem(3,kelem)
          !
          !     Are these points still active
          !
          if(lpoin(jp1t)==1)ip1h=ip1h+1_ip
          if(lpoin(jp2t)==1)ip2h=ip2h+1_ip
          if(lpoin(jp3t)==1)ip3h=ip3h+1_ip

       enddo
       !
       !     Get the lowest point
       !
       imin=1_ip
       iph=ip1h 
       if(ip2h<iph)then
          imin=2_ip
          iph=ip2h
       endif
       if(ip3h<iph)then
          imin=3_ip
          iph=ip3h
       endif
       !
       !     Get the highest point
       !
       imax=1_ip
       iph=ip1h
       if(ip2h>iph)then
          imax=2_ip
          iph=ip2h
       endif
       if(ip3h>iph)then
          imax=3_ip
          iph=ip3h
       endif
       !
       !     Check that the jump is not higher than 1
       !
       if(iph>1)then
          write(*,*)'Jump higher than 1  --> stop '
          stop
       endif
       ! 
       !    Get the points of the top and bottom faces
       ! 
       lebas(1)=ip1b
       lebas(2)=ip2b
       lebas(3)=ip3b
       letop(1)=ip1t
       letop(2)=ip2t
       letop(3)=ip3t
       !
       !    Reorder && create new elements
       !
       icont=ip1h+ip2h+ip3h
       !
       !    Is this already at the same level (low)?
       !
       if(icont.eq.0)then
          ielem=jelem
          cycle
       endif
       !
       !    If it is at the same high level check why it has been deleted
       !
       if(icont.eq.3)then
          ie2=ielem+1
          ie3=ielem+2
          !
          !    If the volumes are not negatives, reintroduce these elements 
          !
          if(ltest(ielem).ne.1 .and. ltest(ielem2).ne.1 .and.  ltest(ielem3).ne.1)then

             lmark(ielem)=1_ip
             lmark(ielem2)=1_ip
             lmark(ielem3)=1_ip

             icase=ltest(ielem)
             if(icase.ne.0)then
                ldiagno(icase)=ldiagno(icase)+1
             endif
             icase=ltest(ielem2)
             if(icase.ne.0)then
                ldiagno(icase)=ldiagno(icase)+1
             endif
             icase=ltest(ielem3)
             if(icase.ne.0)then
                ldiagno(icase)=ldiagno(icase)+1
             endif
             !
             !     Add top face
             !   
             nfnew=nface+1_ip
             call memrea(nfnew,memor_msh,'LFACE','addcover',lface)
             call memrea(nfnew,memor_msh,'LSURF','addcover',lsurf)

             nface=nface+1_ip
             lface(1,nface)=ip1t            
             lface(2,nface)=ip2t            
             lface(3,nface)=ip3t            
             lsurf(nface)=3_ip

             ltest(ielem)=0
             ltest(ielem2)=0
             ltest(ielem3)=0

          endif

          ielem=jelem
          cycle

       endif
       !
       !    Has the face only one point at a higher level?
       !
       if(icont.eq.1)then

          ip1=lebas(imax)
          ip2=lebas(ltab(1,imax))
          ip3=lebas(ltab(2,imax))
          ip4=letop(imax)
          !
          !     Test volume
          !
          x1  =coor(1,ip1)
          y1  =coor(2,ip1)
          z1  =coor(3,ip1)
          dx21=coor(1,ip2)-x1
          dy21=coor(2,ip2)-y1
          dz21=coor(3,ip2)-z1
          dx31=coor(1,ip3)-x1
          dy31=coor(2,ip3)-y1
          dz31=coor(3,ip3)-z1
          dx41=coor(1,ip4)-x1
          dy41=coor(2,ip4)-y1
          dz41=coor(3,ip4)-z1

          v1=dx41*(dy21*dz31-dz21*dy31)        &
               +dy41*(dz21*dx31-dx21*dz31) &
               +dz41*(dx21*dy31-dy21*dx31)

          if(v1.lt.c00)then
             jelem=jelem
             cycle
          endif

          nelem=nelem+1
          elem(1,nelem)=ip1
          elem(2,nelem)=ip2
          elem(3,nelem)=ip3
          elem(4,nelem)=ip4
          lmark(  nelem)=1
          ltest(  nelem)=33
          !
          !     Add top face
          !   
          nfnew=nface+1_ip
          call memrea(nfnew,memor_msh,'LFACE','addcover',lface)
          call memrea(nfnew,memor_msh,'LSURF','addcover',lsurf)

          nface=nface+1_ip
          lface(1,nface)=ip1t            
          lface(2,nface)=ip2t            
          lface(3,nface)=ip3t            
          lsurf(nface)=3_ip

       else

          if(ldiag(lside(imin)).eq.1)then

             ip1a=lebas(imin)
             ip2a=lebas(ltab(1,imin))
             ip3a=lebas(ltab(2,imin))
             ip4a=letop(ltab(2,imin))

             ip1b=lebas(imin)
             ip2b=letop(ltab(2,imin))
             ip3b=letop(ltab(1,imin))
             ip4b=lebas(ltab(1,imin))

          else

             ip1a=lebas(imin)
             ip2a=lebas(ltab(1,imin))
             ip3a=lebas(ltab(2,imin))
             ip4a=letop(ltab(1,imin))

             ip1b=lebas(imin)
             ip2b=letop(ltab(2,imin))
             ip3b=letop(ltab(1,imin))
             ip4b=lebas(ltab(2,imin))

          endif
          !
          !    Check volumes
          !
          x1  =coor(1,ip1a)
          y1  =coor(2,ip1a)
          z1  =coor(3,ip1a)
          dx21=coor(1,ip2a)-x1
          dy21=coor(2,ip2a)-y1
          dz21=coor(3,ip2a)-z1
          dx31=coor(1,ip3a)-x1
          dy31=coor(2,ip3a)-y1
          dz31=coor(3,ip3a)-z1
          dx41=coor(1,ip4a)-x1
          dy41=coor(2,ip4a)-y1
          dz41=coor(3,ip4a)-z1

          v1=dx41*(dy21*dz31-dz21*dy31)     & 
               +dy41*(dz21*dx31-dx21*dz31) &
               +dz41*(dx21*dy31-dy21*dx31)

          if(v1.lt.c00)then
             ielem=jelem
             cycle
          endif

          x1  =coor(1,ip1b)
          y1  =coor(2,ip1b)
          z1  =coor(3,ip1b)
          dx21=coor(1,ip2b)-x1
          dy21=coor(2,ip2b)-y1
          dz21=coor(3,ip2b)-z1
          dx31=coor(1,ip3b)-x1
          dy31=coor(2,ip3b)-y1
          dz31=coor(3,ip3b)-z1
          dx41=coor(1,ip4b)-x1
          dy41=coor(2,ip4b)-y1
          dz41=coor(3,ip4b)-z1

          v2=dx41*(dy21*dz31-dz21*dy31)       &
               +dy41*(dz21*dx31-dx21*dz31) &
               +dz41*(dx21*dy31-dy21*dx31)

          if(v2.lt.c00)then
             ielem=jelem
             cycle
          endif

          nelem=nelem+1
          elem(1,nelem)=ip1a
          elem(2,nelem)=ip2a
          elem(3,nelem)=ip3a
          elem(4,nelem)=ip4a
          lmark(  nelem)=1
          ltest(  nelem)=33

          nelem=nelem+1
          elem(1,nelem)=ip1b
          elem(2,nelem)=ip2b
          elem(3,nelem)=ip3b
          elem(4,nelem)=ip4b
          lmark(  nelem)=1
          ltest(  nelem)=33
          !
          !     Add top face
          !   
          nfnew=nface+1_ip
          call memrea(nfnew,memor_msh,'LFACE','addcover',lface)
          call memrea(nfnew,memor_msh,'LSURF','addcover',lsurf)

          nface=nface+1_ip
          lface(1,nface)=ip1t            
          lface(2,nface)=ip2t            
          lface(3,nface)=ip3t            
          lsurf(nface)=3_ip

       endif

    enddo
    !
    !     Resize lmark
    !
    call memrea(nelem,memor_msh,'LMARK','addcover',lmark)
    call memrea(nelem,memor_msh,'LTEST','addcover',ltest)
    !
    !     Mark the new elements as kept
    !
    do ielem=nele0+1,nelem
       lmark(ielem)=1_ip
    enddo

    call memchk(2_ip,istat,memor_msh,'LPOIN','addcover',lpoin)
    deallocate(lpoin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPOIN','addcover',0_ip)

  end subroutine addcover

  subroutine rmvjmp(elem,nelem,npoin,nnofa,lplay,lmark,nface,nblay)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(in)    :: nnofa,npoin,nface,nblay
    integer(ip),intent(inout) :: nelem
    integer(ip),intent(in)    :: elem(nnofa,nface)
    integer(ip),pointer       :: ptoel2(:),ptoel1(:)
    integer(ip),intent(in)    :: lplay(npoin)  
    integer(ip),intent(inout) :: lmark(nelem)  
    integer(ip)               :: ielem,jelem,kelem,isto,ip1,ip2,ip3,ip4,kelay,jelay
    integer(ip)               :: inofa,ipoin,lplist(3),iface2,iflay,iplist,ielay
    logical(lg)               :: imarked
    integer(4)                :: istat
    !
    !     This sub guarantees that the jump between neighbors is not higher than 1
    !
    !
    !     Get the elements surrounding the points
    !
    call ptoelm(elem,nelem,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Loop on iteration
    !
    do 
       !
       !     Loop on elements
       !  
       ielem=-1_ip
       !
       !     Set flag for exit condition
       ! 
       imarked=.false.

       do 
          !
          !     Update counter
          !
          ielem=ielem+3 
          !
          !     Did we exhaust the list?
          !
          if(ielem>nelem)exit
          !
          !     Is the element marked for deletion?
          !
          if(lmark(ielem)==0)then
             !
             !     Get the level
             !
             ielay=lplay(elem(1,ielem))
             !
             !     Create the list of points
             !
             jelem=ielem+2
             lplist(1)=elem(1,kelem)
             lplist(2)=elem(4,kelem)
             lplist(3)=elem(3,kelem)
             !
             !     Mark the elements that are too high in the column
             !
             do iplist=1,3
                ipoin=lplist(iplist) 
                do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1
                   kelem=ptoel1(isto)  
                   !
                   !     Has this face already been deleted before
                   !
                   if(lmark(kelem)==1)then
                      ip1=elem(1,kelem)
                      ip2=elem(2,kelem)
                      ip3=elem(3,kelem)
                      ip4=elem(4,kelem)
                      kelay=min(lplay(ip1),lplay(ip2),lplay(ip3),lplay(ip4))
                      if(kelay.gt.ielay) then
                         lmark(kelem)=0
                         imarked=.true.
                      endif
                   endif
                enddo
             enddo
          endif
       enddo
       !
       !     Exit condition
       !
       if(imarked.eqv. .false.)exit
       !
       !     Mark the elements not marked in the column
       ! 
       call mrkcol(nface,nblay,lmark,nelem)

    enddo

    call memchk(2_ip,istat,memor_msh,'PTOEL1','rmvjmp',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','rmvjmp',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','rmvjmp',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','rmvjmp',0_ip)

  end subroutine rmvjmp

  subroutine onemore(nnode,nelem,elem,lmark,lplay,npoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nelem,nnode,npoin
    integer(ip),intent(inout) :: lmark(nelem)
    integer(ip),intent(in)    :: elem(nnode,nelem),lplay(npoin)
    integer(ip)               :: ielem,jelem,ip1b,jp1b,jelay,ielay 
    !
    !     This subroutine marks all the elements for which the 
    !     next element have been already marked
    !

    !
    !    Loop on elements
    !  
    ielem=-2_ip

    do 
       !
       !     Find the next element marked for deletion
       !
       ielem=ielem+3_ip
       !
       !     Did we exhaust the list?
       !
       if(ielem>nelem)exit
       !
       !     Is the element still active?
       !
       if(lmark(ielem)==1)cycle 
       !
       !     We have found a marked element
       !     Get ielay
       !
       ip1b=elem(1,ielem)
       ielay=lplay(ip1b)
       !
       !     Now go three steps before 
       !
       jelem=ielem-3_ip
       !
       !     Are we still in the range?
       !
       if(jelem<1)cycle
       !
       !     Has jelem been deleted?
       !
       if(lmark(jelem)==0)cycle
       !
       !     Get jelay
       !
       jp1b=elem(1,jelem)
       jelay=lplay(jp1b)
       !
       !     Is ielem above jelem?
       !
       if(jelay>=ielay)cycle
       !
       !     Mark jelem
       !
       lmark(jelem)=0_ip
       lmark(jelem+1)=0_ip
       lmark(jelem+2)=0_ip

    enddo

  end subroutine onemore

  subroutine mrkcol(nface,nblay,lmark,nelem)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(in)    :: nface,nelem,nblay
    integer(ip),intent(inout) :: lmark(nelem)
    integer(ip)               :: nele,iface,iblay,jblay
    !
    !     Initialize the pointer
    !
    nele=0_ip
    !
    !     Loop on sides to check elements in column no marked
    !   
    do iface=1,nface
       !
       !     Loop on layers
       !  
       do iblay=1,nblay

          nele=nele+1_ip
          !
          !     Has this first element been marked for deletion
          ! 
          if(lmark(nele)==0)then
             !
             !     Mark next elements in the column
             !
             nele=nele+1_ip 
             lmark(nele)=0_ip
             nele=nele+1_ip 
             lmark(nele)=0_ip

             do jblay=iblay+1,nblay
                nele=nele+1_ip
                lmark(nele)=0_ip
                nele=nele+1_ip
                lmark(nele)=0_ip 
                nele=nele+1_ip
                lmark(nele)=0_ip 
             enddo
             !
             !     And go to the next face
             ! 
             exit

          endif

          nele=nele+1_ip
          !
          !     Has this second element been marked for deletion
          ! 
          if(lmark(nele)==0)then
             !
             !     Mark next elements in the column
             !
             lmark(nele-1)=0_ip
             nele=nele+1_ip 
             lmark(nele)=0_ip

             do jblay=iblay+1,nblay
                nele=nele+1_ip
                lmark(nele)=0_ip
                nele=nele+1_ip
                lmark(nele)=0_ip
                nele=nele+1_ip
                lmark(nele)=0_ip 
             enddo
             !
             !     And go to the next face
             ! 
             exit

          endif

          nele=nele+1_ip
          !
          !     Has this third element been marked for deletion
          ! 
          if(lmark(nele)==0)then
             !
             !     Mark next elements in the column
             !
             lmark(nele-2)=0_ip
             lmark(nele-1)=0_ip

             do jblay=iblay+1,nblay
                nele=nele+1_ip
                lmark(nele)=0_ip
                nele=nele+1_ip
                lmark(nele)=0_ip 
                nele=nele+1_ip
                lmark(nele)=0_ip 
             enddo
             !
             !     And go to the next face
             ! 
             exit

          endif

       enddo

    enddo

  end subroutine mrkcol

  subroutine chkface(nelem,elem,nnode,nface,nnofa,lface,npoin)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(in)       :: nnode,nnofa,nface,npoin,nelem
    integer(ip),intent(in)       :: lface(nnofa,nface),elem(nnode,nelem)                 
    integer(ip),pointer          :: ptoel1(:),ptoel2(:),ptofa1(:),ptofa2(:)  
    integer(ip),pointer          :: lmark(:),eltoel(:,:)  
    integer(ip)                  :: ielem,iface,inode,isto,ip1,ip2,ip3,ipa,ipb,ipc  
    logical(lg)                  :: ifound
    integer(4)                   :: istat
    integer(ip)                  :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/)) 
    !
    !     This sub checks that the faces are consistent with the volume
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','chkface',lmark)
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptofa1,ptofa2 )
    !
    !     Get the elements surrounding the points
    !
    call ptoelm(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !     Get the elements surrounding elements
    !
    call tetote(elem,nnode,nelem,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Loop on elements
    !  
    do ielem=1,nelem
         !
         !     Loop on nodes
         !
         do inode=1,nnode 
            !
            !     Do we have a boundary face?
            ! 
            if(eltoel(inode,ielem)==0)then
               !
               !     Get the points
               !
               ip1=elem(ltab(1,inode),ielem) 
               ip2=elem(ltab(2,inode),ielem) 
               ip3=elem(ltab(3,inode),ielem) 
               !
               !     Mark the points
               !
               lmark(ip1)=1_ip 
               lmark(ip2)=1_ip 
               lmark(ip3)=1_ip 
               !
               !     Check if we have the surface
               ! 
               ifound=.false.
               do isto=ptofa2(ip1),ptofa2(ip1+1)-1
                  iface=ptofa1(isto)
                  ipa=lface(1,iface)
                  ipb=lface(2,iface)
                  ipc=lface(3,iface)
                  if(lmark(ipa)+lmark(ipb)+lmark(ipc)==3)then
                     ifound=.true.
                     exit
                  endif
               enddo
               if(ifound.eqv. .false.)then
                  write(*,*)'Error in chkface, face not found, iface=',iface
                  stop
               endif
               !
               !     Clean up lmark
               !
               lmark(ip1)=0_ip 
               lmark(ip2)=0_ip 
               lmark(ip3)=0_ip 

            endif
         enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'PTOFA1','chkface',ptofa1)
    deallocate(ptofa1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOFA1','chkface',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOFA2','chkface',ptofa2)
    deallocate(ptofa2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOFA2','chkface',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','chkface',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','chkface',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','chkface',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','chkface',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','chkface',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','chkface',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','chkface',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','chkface',0_ip)

  end subroutine chkface

 subroutine outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim
    integer(ip), intent(in)      :: lface(nnofa,nface),lsurf(nface)
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
       write(50,300)icont,lface(1,i),lface(2,i),lface(3,i),lsurf(i)
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








end module mod_extrpr




