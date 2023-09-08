module mod_clean

contains 

  subroutine colpnt(lface,nnofa,nface,coor,ndim,npoin,nnosi,nside,lside,lsurf,nsurf,lline,nline)
    use def_kintyp, only :  ip,rp,lg
    use mod_memchk
    use def_meshin, only    : memor_msh
    use def_meshin, only    :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    use mod_cart, only    : boxbin,cartp 
    use mod_mshtol, only    : ptoelm 
    implicit none
    integer(ip),intent(in)        :: ndim,nnofa,nnosi,nsurf,nline
    integer(ip),intent(inout)     :: npoin,nface,nside
    integer(ip),intent(inout)     :: lsurf(nface),lline(nside)
    integer(ip),intent(inout)     :: lface(nnofa,nface),lside(nnosi,nside)
    real(rp),intent(inout)        :: coor(ndim,npoin)
    integer(ip),pointer           :: ptoel1(:),ptoel2(:),lcart(:),lrenu(:)
    real(rp),pointer              :: rsize(:)
    real(rp)                      :: c00,rtol,redge,rx,ry,rz,bbox(2,3) 
    integer(ip)                   :: ipoin,nedge,isto,iface,inofa,jpoin,ichk
    integer(ip)                   :: iside,npoi0,nfac0,nsid0
    integer(4)                    :: istat
    !
    !     This sub eliminates the points which are too close
    !
    nullify(ptoel1,ptoel2)
    c00=0.0d+00
    !
    !     Allocate section
    !
    allocate(rsize(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'RSIZE','colpnt',rsize)
    allocate(lcart(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LCART','colpnt',lcart)
    allocate(lrenu(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','colpnt',lrenu)
    !
    !     Get the faces surrounding the points
    ! 
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2) 
    !
    !     Compute bbox of the boundary points
    !
    call boxbin(coor,ndim,npoin,bbox)
    !
    !     Give the size
    !  
    do ipoin=1,npoin
       lrenu(ipoin)=ipoin
       redge=c00
       nedge=0_ip
       do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(isto)
          do inofa=1,nnofa
             jpoin=lface(inofa,iface)
             if(lrenu(jpoin)/=ipoin)then
                lrenu(jpoin)=ipoin
                rx=coor(1,ipoin)-coor(1,jpoin)
                ry=coor(2,ipoin)-coor(2,jpoin)
                rz=coor(3,ipoin)-coor(3,jpoin)
                redge=redge+sqrt(rx*rx+ry*ry+rz*rz)
                nedge=nedge+1_ip
             endif
          enddo
       enddo
       rsize(ipoin)=redge/real(nedge)
    enddo
    !
    !     Filter the points
    !
    call filvol(bbox,coor,ndim,npoin,rsize,lrenu)
    !
    !     Did we change something?
    !
    ichk=0_ip
    npoi0=npoin
    npoin=0_ip
    do ipoin=1,npoi0
       if(lrenu(ipoin)==ipoin)then
          npoin=npoin+1
          lrenu(ipoin)=-npoin
       else
          ichk=1
       endif
    enddo
    
    if(ichk==1)then
       !
       !     Give new numbering to deleted points
       !
       do ipoin=1,npoi0
          if(lrenu(ipoin)>0)then
             lrenu(ipoin)=abs(lrenu(lrenu(ipoin)))
          endif 
       enddo  
       !
       !     Get renumbering
       ! 
       npoin=0_ip
       do ipoin=1,npoi0
          if(lrenu(ipoin)<0)then
             npoin=npoin+1
             lrenu(ipoin)=npoin
             coor(1,npoin)=coor(1,ipoin) 
             coor(2,npoin)=coor(2,ipoin) 
             coor(3,npoin)=coor(3,ipoin) 
          endif
       enddo
       !
       !     Renumber the faces
       !
       do iface=1,nface
          do inofa=1,nnofa
             lface(inofa,iface)=lrenu(lface(inofa,iface))
          enddo
       enddo
       !
       !     Do we have repeated points ?
       !
       nfac0=nface
       nface=0_ip
       do iface=1,nfac0
          if(lface(1,iface)/=lface(2,iface) .and. lface(2,iface)/=lface(3,iface) &
             .and. lface(3,iface)/=lface(1,iface))then
             nface=nface+1
             lface(1,nface)=lface(1,iface)
             lface(2,nface)=lface(2,iface)
             lface(3,nface)=lface(3,iface)
             lsurf(nface)=lsurf(iface)
          endif 
       enddo
       !
       !     Renumber the sides
       !
       do iside=1,nside
          lside(1,iside)=lrenu(lside(1,iside))
          lside(2,iside)=lrenu(lside(2,iside))
       enddo
       !
       !     Do we have repeated points ?
       !
       nsid0=nside
       nside=0_ip
       do iside=1,nsid0
          if(lside(1,iside)/=lside(2,iside))then
             nside=nside+1
             lface(1,nside)=lface(1,iside)
             lface(2,nside)=lface(2,iside)
             lline(nside)=lline(iside)
          endif 
       enddo
       write(*,*)'Degenerate points found -->  collapsed '
       write(*,*)'npoin=',npoin
       write(*,*)'nface=',nface
       write(*,*)'nside=',nside
       !
       !     Write new file
       !
       write(*,*)'Writing new file dekhna_new.geo'
       open(unit=50,file='dekhna_new.geo',status='unknown')
       rewind 50
       write(50,*)'npoin=',npoin
       do ipoin=1,npoin
          write(50,*)ipoin,coor(1,ipoin),coor(2,ipoin),coor(3,ipoin)
       enddo
       write(50,*)'nline=',nline
       write(50,*)'nridge=',nside
       do iside=1,nside
          write(50,*)iside,lside(1,iside),lside(2,iside),lline(iside)
       enddo
       write(50,*)'nsurf=',nsurf
       write(50,*)'nface=',nface
       do iface=1,nface
          write(50,*)iface,lface(1,iface),lface(2,iface),lface(3,iface),lsurf(iface)
       enddo

       close(50)
    else
   
       write(*,*)'No close points found' 

    endif

    call memchk(2_ip,istat,memor_msh,'RSIZE','colpnt',rsize)
    deallocate(rsize,stat=istat)
    if(istat/=0) call memerr(2_ip,'RSIZE','colpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCART','colpnt',lcart)
    deallocate(lcart,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCART','colpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'LRENU','colpnt',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','colpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','colpnt',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','colpnt',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','colpnt',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','colpnt',0_ip)

  end subroutine colpnt

  subroutine filvol(bboxp,coor,ndim,npoin,rsize,lrenu)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: ndim,npoin 
    integer(ip),intent(inout)   :: lrenu(npoin)
    real(rp),intent(in)         :: bboxp(ndim,2),coor(ndim,npoin),rsize(npoin) 
    integer(ip)                 :: nbin,nbinx,istack,imin,jmin,kmin,jpoin
    integer(ip)                 :: ibin,ipoin,nbinxy,i,ntot,iplace,l,ichk
    integer(ip)                 :: imax,jmax,kmax,k,kbinx,j,jbinx,isto0,isto1
    integer(ip),pointer         :: btopo1(:),btopo2(:)
    real(rp)                    :: dxbin,dybin,dzbin,rlen
    real(rp)                    :: vmin(ndim),vmax(ndim),c13
    real(rp)                    :: dx,dy,dz,rdist,rtol 
    integer(4)                  :: istat
    !
    !    This sub filters exactly the points that are too close by a bin sort
    !
    c13=1.0d+00/3.0d+00
    rtol=1.0d-6
    !
    !     Allocate the bin to point pointer
    !
    nbin=npoin
    nbinx=floor(nbin**c13)+50_ip
    nbin=nbinx*nbinx*nbinx
    nbinxy=nbinx*nbinx
    allocate(btopo2(nbin+1),stat=istat)
    call memchk(zero,istat,memor_msh,'BTOPO2','cart',btopo2)

    dxbin=(bboxp(1,2)-bboxp(1,1))/nbinx
    dybin=(bboxp(2,2)-bboxp(2,1))/nbinx
    dzbin=(bboxp(3,2)-bboxp(3,1))/nbinx

    do ipoin=1,npoin 

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
    do ipoin=1,npoin 

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
    !     Initialize lrenu
    !
    do ipoin=1,npoin
       lrenu(ipoin)=ipoin
    enddo
    !
    !     Loop on iterations
    ! 
    do
       !
       !     Initialize the 'something done' flag
       !
       ichk=0_ip
       !
       !     Loop on the points
       !
       do ipoin=1,npoin 
          !
          !     Create the local bbox
          !
          rlen=rsize(ipoin)
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
                      if(lrenu(jpoin)==lrenu(ipoin))cycle
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
                      if(rdist<rlen*rtol)then
                         if(lrenu(ipoin)<lrenu(jpoin))then
                            lrenu(jpoin)=lrenu(ipoin)
                         else if(lrenu(jpoin)<lrenu(ipoin))then
                            lrenu(ipoin)=lrenu(jpoin)
                         endif
                         ichk=1_ip 
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo

       do ipoin=1,npoin
          lrenu(ipoin)=lrenu(lrenu(ipoin))
       enddo

       if(ichk==0)exit

    enddo

    call memchk(2_ip,istat,memor_msh,'BTOPO1','filvol',btopo1)
    deallocate(btopo1,stat=istat)
    if(istat/=0) call memerr(2_ip,'BTOPO1','filvol',0_ip)
    call memchk(2_ip,istat,memor_msh,'BTOPO2','filvol',btopo2)
    deallocate(btopo2,stat=istat)
    if(istat/=0) call memerr(2_ip,'BTOPO2','filvol',0_ip)

  end subroutine filvol


end module mod_clean
