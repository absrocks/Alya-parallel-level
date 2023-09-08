module mod_bin3d

contains

  subroutine inibin(bbox,bboxbin,lbin1,lbin2,nbin,nbinx,nbiny,nbinz,&
       npoin,dxbin,dybin,dzbin)
    use def_kintyp, only :  ip,rp,lg
    use mod_memchk
    use def_meshin, only :memor_msh  
    implicit none
    integer(ip),intent(in)    :: npoin
    integer(ip),intent(inout) :: nbin,nbinx,nbiny,nbinz
    real(rp),intent(in)       :: bbox(3,2)
    real(rp),intent(inout)    :: bboxbin(3,2),dxbin,dybin,dzbin
    integer(ip), pointer      :: lbin1(:),lbin2(:)
    real(rp)                  :: lx,ly,lz,c13,rscal,rcenter(3),c05
    integer(4)                :: istat 

    c13=1.0d+00/3.0d+00
    c05=0.5d+00
    rscal=1.1d+00  

    nbinx=npoin**c13
    nbinx=max(nbinx,50_ip)
    nbin=nbinx*nbinx*nbinx
    nbiny=nbinx
    nbinz=nbinx  

    rcenter(1)=c05*(bbox(1,1)+bbox(1,2)) 
    rcenter(2)=c05*(bbox(2,1)+bbox(2,2)) 
    rcenter(3)=c05*(bbox(3,1)+bbox(3,2)) 

    bboxbin(1,1)=rcenter(1)+rscal*(bbox(1,1)-rcenter(1))
    bboxbin(2,1)=rcenter(2)+rscal*(bbox(2,1)-rcenter(2))
    bboxbin(3,1)=rcenter(3)+rscal*(bbox(3,1)-rcenter(3))

    bboxbin(1,2)=rcenter(1)+rscal*(bbox(1,2)-rcenter(1))
    bboxbin(2,2)=rcenter(2)+rscal*(bbox(2,2)-rcenter(2))
    bboxbin(3,2)=rcenter(3)+rscal*(bbox(3,2)-rcenter(3))

    lx=bboxbin(1,2)-bboxbin(1,1)
    ly=bboxbin(2,2)-bboxbin(2,1)
    lz=bboxbin(3,2)-bboxbin(3,1)

    dxbin=lx/nbinx 
    dybin=ly/nbinx 
    dzbin=lz/nbinx 
    allocate(lbin2(nbin+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LBIN2','inibin',lbin2)
    allocate(lbin1(nbin),stat=istat)
    call memchk(zero,istat,memor_msh,'LBIN1','inibin',lbin1)

  end subroutine inibin

  subroutine gtclos(ipoin,lbin1,lbin2,ipclos,nbin,coor,ndim,bboxbin,&
       dx,dy,dz,nx,ny,nz,npoin)
    use def_kintyp, only :  ip,rp,lg
    implicit none
    integer(ip),intent(in)       ::  ipoin,nbin,ndim,nx,ny,nz,npoin
    integer(ip),intent(in)       ::  lbin2(nbin+1),lbin1(*)
    real(rp),intent(in)          ::  coor(ndim,npoin),bboxbin(3,2),dx,dy,dz
    integer(ip),intent(inout)    ::  ipclos
    integer(ip)                  :: ix,iy,iz,ibin,nlay,nyz
    integer(ip)                  :: jx,jy,jz,ilay,isto,kx,ky,kz
    integer(ip),save             :: npin=0_ip      

    if(npin<50)then
       npin=npin+1_ip
       ipclos=1_ip
       return
    endif

    nyz=ny*nz
    ix=floor((coor(1,ipoin)-bboxbin(1,1))/dx)+1_ip 
    iy=floor((coor(2,ipoin)-bboxbin(2,1))/dy)+1_ip 
    iz=floor((coor(3,ipoin)-bboxbin(3,1))/dz)+1_ip 

    if(ix>=1 .and.ix<=nx)then
       if(iy>=1 .and.iy<=ny)then
          if(iz>=1 .and.iz<=nz)then
             ibin=(ix-1)*nyz+(iy-1)*nz+iz
             !
             !     Try ibin first
             !
             do isto=lbin2(ibin),lbin2(ibin+1)-1
                ipclos=lbin1(isto) 
                if(ipclos==0)cycle 
                return 
             enddo
          endif
       endif
    endif
    !
    !     Loop until point found
    ! 
    nlay=0_ip
    do 
       nlay=nlay+1
       !
       !     Top layer center
       !
       kx=ix
       ky=iy
       kz=iz+nlay
       if(kx>=1 .and.kx<=nx)then
          if(ky>=1 .and.ky<=ny)then
             if(kz>=1 .and.kz<=nz)then
                ibin=(kx-1)*nyz+(ky-1)*nz+kz
                do isto=lbin2(ibin),lbin2(ibin+1)-1
                   ipclos=lbin1(isto)
                   if(ipclos==0)cycle 
                   return 
                enddo

             endif
          endif
       endif

       do ilay=1,nlay 
          !
          !     Top layer 
          !
          do jy=-ilay+1,ilay
             kx=ix+ilay
             ky=iy+jy
             kz=iz+nlay
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

          do jx=-ilay+1,ilay
             kx=ix-jx
             ky=iy+ilay
             kz=iz+nlay
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

          do jy=-ilay+1,ilay
             kx=ix-ilay
             ky=iy-jy
             kz=iz+nlay
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

          do jx=-ilay+1,ilay
             kx=ix+jx
             ky=iy-ilay
             kz=iz+nlay
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo
       enddo
       !
       !     Center
       !
       do jz=nlay-1,-nlay+1,-1

          do jy=-nlay+1,nlay
             kx=ix+nlay
             ky=iy+jy
             kz=iz+jz
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

          do jx=-nlay+1,nlay
             kx=ix-jx
             ky=iy+nlay
             kz=iz+jz
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

          do jy=-nlay+1,nlay
             kx=ix-nlay
             ky=iy-jy
             kz=iz+jz
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

          do jx=-nlay+1,nlay
             kx=ix+jx
             ky=iy-nlay
             kz=iz+jz
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

       enddo
       !
       !    Lower part center
       !
       kx=ix
       ky=iy
       kz=iz-nlay
       if(kx>=1 .and.kx<=nx)then
          if(ky>=1 .and.ky<=ny)then
             if(kz>=1 .and.kz<=nz)then
                ibin=(kx-1)*nyz+(ky-1)*nz+kz
                do isto=lbin2(ibin),lbin2(ibin+1)-1
                   ipclos=lbin1(isto)
                   if(ipclos==0)cycle 
                   return 
                enddo
             endif
          endif
       endif

       do ilay=1,nlay 
          !
          !     Bottom layer 
          !
          do jy=-ilay+1,ilay
             kx=ix+ilay
             ky=iy+jy
             kz=iz-nlay
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

          do jx=-ilay+1,ilay
             kx=ix-jx
             ky=iy+ilay
             kz=iz-nlay
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

          do jy=-ilay+1,ilay
             kx=ix-ilay
             ky=iy-jy
             kz=iz-nlay
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo

          do jx=-ilay+1,ilay
             kx=ix+jx
             ky=iy-ilay
             kz=iz-nlay
             if(kx>=1 .and.kx<=nx)then
                if(ky>=1 .and.ky<=ny)then
                   if(kz>=1 .and.kz<=nz)then
                      ibin=(kx-1)*nyz+(ky-1)*nz+kz
                      do isto=lbin2(ibin),lbin2(ibin+1)-1
                         ipclos=lbin1(isto)
                         if(ipclos==0)cycle 
                         return 
                      enddo
                   endif
                endif
             endif
          enddo
       enddo
    enddo

  end subroutine gtclos

  subroutine insertbin(ipoin,coor,npoin,ndim,lbin1,lbin2,nbin,&
       dx,dy,dz,nx,ny,nz,bboxbin)
    use def_kintyp, only :  ip,rp,lg
    implicit none
    integer(ip),intent(in)   :: npoin,ipoin,ndim,nbin,nx,ny,nz
    integer(ip),intent(in)   :: lbin2(nbin+1)
    integer(ip),intent(inout)   :: lbin1(*)
    real(rp),intent(in)      :: coor(ndim,npoin),dx,dy,dz,bboxbin(3,2)
    integer(ip)              :: ix,iy,iz,ibin,nyz,isto
    nyz=ny*nz


    ix=floor((coor(1,ipoin)-bboxbin(1,1))/dx)+1_ip 
    iy=floor((coor(2,ipoin)-bboxbin(2,1))/dy)+1_ip 
    iz=floor((coor(3,ipoin)-bboxbin(3,1))/dz)+1_ip 
    ibin=(ix-1)*nyz+(iy-1)*nz+iz

    if(ibin>nbin)then
       write(*,*)'Error in insertbin:',ibin,nbin
       stop
    endif

    do isto=lbin2(ibin),lbin2(ibin+1)-1
       if(lbin1(isto)==0)then
          lbin1(isto)=ipoin
          return
       endif
    enddo
    !
    !     Error
    !
    write(*,*)'Error in insertbin, place not found'
    stop

  end subroutine insertbin

  subroutine filter(lbin1,lbin2,nbin,npoin,ndim,coor,dx,dy,dz,&
       nx,ny,nz,rsizn,coorold,npold,bboxbin,lmark,rsize)
    use def_kintyp, only :  ip,rp,lg
    use mod_memchk
    use def_meshin, only :memor_msh  
    implicit none
    integer(ip),intent(in)    ::  nbin,npoin,nx,ny,nz,ndim,npold
    integer(ip),intent(inout) ::  lbin2(nbin+1),lmark(npoin)
    real(rp),intent(in)       ::  rsize(npoin),coor(ndim,npoin),coorold(ndim,npold)
    real(rp),intent(in)       ::  dx,dy,dz,bboxbin(3,2),rsizn(npold)
    integer(ip),pointer       ::  lbin1(:)
    integer(ip),pointer       ::  lbint(:),lbin(:)
    integer(ip)               ::  nbnew,ix,iy,iz,ibin,jx,jy,jz,jbin,ipoin,nyz,iplace,jpoin
    integer(ip)               ::  isto,ipnt,ipnt0,nxmax,nymax,nzmax,nxmin,nymin,nzmin,jsto,lbinmin
    real(rp)                  ::  rscal,rscal2,rx,ry,rz,rl2
    integer(4)                :: istat
    !
    !     This subroutine filters the points given in coor
    !     It is assumed that the arrays lbin1,lbin2 are compact, ie there is no
    !     blank in the grid.
    !     On output, there will be empty spaces in the bin for the points to be inserted later on.
    !
    !
    !     Allocate help array
    !
    allocate(lbin(nbin),stat=istat)
    call memchk(zero,istat,memor_msh,'LBIN','filter',lbin)
    nyz=ny*nz
    !
    !     Insert the new points in the bin
    !
    nbnew=0_ip
    do ipoin=1,npoin   
       ix=floor((coor(1,ipoin)-bboxbin(1,1))/dx)+1_ip 
       iy=floor((coor(2,ipoin)-bboxbin(2,1))/dy)+1_ip 
       iz=floor((coor(3,ipoin)-bboxbin(3,1))/dz)+1_ip 
       ibin=(ix-1)*nyz+(iy-1)*nz+iz
       !
       !     Accumulate in lbin
       !
       lbin(ibin)=lbin(ibin)+1_ip
       nbnew=nbnew+1_ip
    enddo
    !
    !     Resize the bin
    !
    ipnt=lbin2(nbin+1)-1+nbnew
    allocate(lbint(ipnt),stat=istat)
    call memchk(zero,istat,memor_msh,'LBINT','filter',lbint)
    !
    !     Shift
    !
    do ibin=nbin,1,-1
       !
       !     Get some new space
       !
       ipnt=ipnt-lbin(ibin)
       ipnt0=ipnt
       !
       !     Write the old values
       !
       do isto=lbin2(ibin+1)-1,lbin2(ibin),-1
          lbint(ipnt)=lbin1(isto)
          ipnt=ipnt-1
       enddo
       !
       !     Rewrite lbin2
       !
       lbin2(ibin+1)=ipnt0+1

    enddo

    call memchk(2_ip,istat,memor_msh,'LBIN1','filter',lbin1)
    deallocate(lbin1,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBIN1','filter',0_ip)

    lbin1=>lbint
    !
    !     Write the new points in the bin
    !
    do ipoin=1,npoin   
       ix=floor((coor(1,ipoin)-bboxbin(1,1))/dx)+1_ip 
       iy=floor((coor(2,ipoin)-bboxbin(2,1))/dy)+1_ip 
       iz=floor((coor(3,ipoin)-bboxbin(3,1))/dz)+1_ip 
       ibin=(ix-1)*nyz+(iy-1)*nz+iz
       iplace=lbin2(ibin+1)
       lbin2(ibin+1)=lbin2(ibin+1)+1
       lbin1(iplace)=ipoin
    enddo
    !
    !     Then filter the points 
    !
    do ix=1,nx
       do iy=1,ny
          do iz=1,nz
             !
             !     Loop on the inner points
             !
             ibin=(ix-1)*nyz+(iy-1)*nz+iz 
             !  
             !     First the loop on the old points
             !
             do isto=lbin2(ibin),lbin2(ibin+1)-lbin(ibin)-1
                ipoin=lbin1(isto)
                !
                !     Get the scale
                !
                rscal=rsizn(ipoin)
                rscal2=rscal*rscal 
                !
                !     Get the neighboring region
                !
                nxmin=floor((coorold(1,ipoin)-rscal-bboxbin(1,1))/dx)+1_ip 
                nymin=floor((coorold(2,ipoin)-rscal-bboxbin(2,1))/dy)+1_ip 
                nzmin=floor((coorold(3,ipoin)-rscal-bboxbin(3,1))/dz)+1_ip 
                nxmax=floor((coorold(1,ipoin)+rscal-bboxbin(1,1))/dx)+1_ip 
                nymax=floor((coorold(2,ipoin)+rscal-bboxbin(2,1))/dy)+1_ip 
                nzmax=floor((coorold(3,ipoin)+rscal-bboxbin(3,1))/dz)+1_ip 
                !
                !     Cut off
                ! 
                nxmin=max(1_ip,nxmin)
                nymin=max(1_ip,nymin)
                nzmin=max(1_ip,nzmin)
                nxmax=min(nxmax,nx)
                nymax=min(nymax,ny)
                nzmax=min(nzmax,nz)
                !
                !     Loop on local cells
                !
                do jx=nxmin,nxmax
                   do jy=nymin,nymax
                      do jz=nzmin,nzmax
                         jbin=(jx-1)*nyz+(jy-1)+jz 
                         !
                         !     Loop on local points 
                         !
                         do jsto=lbin2(jbin+1)-lbin(jbin),lbin2(jbin+1)-1
                            jpoin=lbin1(jsto)
                            if(jpoin==ipoin)cycle
                            if(jpoin<0)cycle
                            !
                            !     Compute distance
                            ! 
                            rx=coor(1,ipoin)-coor(1,jpoin)
                            ry=coor(2,ipoin)-coor(2,jpoin)
                            rz=coor(3,ipoin)-coor(3,jpoin)
                            rl2=rx*rx+ry*ry+rz*rz
                            !
                            !     Is jpoin too close from ipoin?
                            !  
                            if(rl2<rscal2)then
                               lmark(jpoin)=1_ip
                               lbin1(jsto)=-lbin1(jsto)
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    !  
    !     Then the new points
    !
    do ix=1,nx
       do iy=1,ny
          do iz=1,nz
             !
             !     Loop on the inner points
             !
             ibin=(ix-1)*nyz+(iy-1)*nz+iz 
 
             do isto=lbin2(ibin+1)-lbin(ibin),lbin2(ibin+1)-1
                ipoin=lbin1(isto)
                !
                !   Has the point been marked before?
                !
                if(ipoin<0)cycle
                !
                !     Get the scale
                !
                rscal=rsize(ipoin)
                rscal2=rscal*rscal 
                !
                !     Get the neighboring region
                !
                nxmin=floor((coor(1,ipoin)-rscal-bboxbin(1,1))/dx)+1_ip 
                nymin=floor((coor(2,ipoin)-rscal-bboxbin(2,1))/dy)+1_ip 
                nzmin=floor((coor(3,ipoin)-rscal-bboxbin(3,1))/dz)+1_ip 
                nxmax=floor((coor(1,ipoin)+rscal-bboxbin(1,1))/dx)+1_ip 
                nymax=floor((coor(2,ipoin)+rscal-bboxbin(2,1))/dy)+1_ip 
                nzmax=floor((coor(3,ipoin)+rscal-bboxbin(3,1))/dz)+1_ip 
                !
                !     Cut off
                ! 
                nxmin=max(1_ip,nxmin)
                nymin=max(1_ip,nymin)
                nzmin=max(1_ip,nzmin)
                nxmax=min(nxmax,nx)
                nymax=min(nymax,ny)
                nzmax=min(nzmax,nz)
                !
                !     Loop on local cells
                !
                do jx=nxmin,nxmax
                   do jy=nymin,nymax
                      do jz=nzmin,nzmax
                         jbin=(jx-1)*nyz+(jy-1)*nz+jz 
                         !
                         !     Loop on local points 
                         !
                         do jsto=lbin2(jbin+1)-lbin(jbin),lbin2(jbin+1)-1
                            jpoin=lbin1(jsto)
                            if(jpoin==ipoin)cycle
                            if(jpoin<0)cycle
                            !
                            !     Compute distance
                            ! 
                            rx=coor(1,ipoin)-coor(1,jpoin)
                            ry=coor(2,ipoin)-coor(2,jpoin)
                            rz=coor(3,ipoin)-coor(3,jpoin)
                            rl2=rx*rx+ry*ry+rz*rz
                            !
                            !     Is jpoin too close from ipoin?
                            !  
                            if(rl2<rscal2)then
                               lmark(jpoin)=1_ip
                               lbin1(jsto)=-lbin1(jsto) 
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !    Take out marked points 
    !  
    lbinmin=1
    ipnt=0_ip
    do ibin=1,nbin
       do isto=lbinmin,lbin2(ibin+1)-1
          ipoin=lbin1(isto)
          if(ipoin>0)then
             ipnt=ipnt+1_ip
             lbin1(ipnt)=ipoin
          else
             lbin(ibin)=lbin(ibin)-1  
          endif
       enddo
       lbinmin=lbin2(ibin+1)
       lbin2(ibin+1)=ipnt+1

    enddo
    !
    !     Do not compress lbin1, leave blanks
    !
    do ibin=1,nbin
       do isto=lbin2(ibin+1)-lbin(ibin),lbin2(ibin+1)-1
          lbin1(isto)=0
       enddo
    enddo

    !call compresbin(lbin2,lbin1,nbin)

    call memchk(2_ip,istat,memor_msh,'LBIN','filter',lbin)
    deallocate(lbin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBIN','filter',0_ip)

  end subroutine filter

  subroutine compresbin(lbin2,lbin1,nbin)
    use def_kintyp, only :  ip,rp,lg
    implicit none
    integer(ip),intent(in)      ::  nbin
    integer(ip),intent(inout)   ::  lbin2(nbin+1),lbin1(*)
    integer(ip)                 :: ipnt,lbin2m,ipoin,isto
    integer(ip)               ::  ibin

    ipnt=0_ip
    lbin2m=lbin2(ibin+1)-1
    do ibin=1,nbin
       do isto=lbin2(ibin),lbin2m
          ipoin=lbin1(isto)
          if(ipoin/=0)then
             ipnt=ipnt+1
             lbin1(ipnt)=ipoin
          endif
       enddo
       lbin2m=lbin2(ibin+1)-1
       lbin2(ibin+1)=ipnt+1
    enddo

  end subroutine compresbin

  subroutine insertini(coor,npoin,ndim,lbin1,lbin2,nbin,dx,dy,dz,nx,ny,nz,bboxbin)
    use def_kintyp, only :  ip,rp,lg
    use mod_memchk
    use def_meshin, only : memor_msh
    implicit none
    integer(ip),intent(in)    :: npoin,ndim,nbin,nx,ny,nz
    integer(ip),intent(inout) :: lbin2(nbin+1)
    integer(ip),pointer       :: lbin1(:)
    real(rp),intent(in)       :: coor(ndim,npoin),dx,dy,dz,bboxbin(3,2)
    integer(ip)               :: ix,iy,iz,ibin,nyz,isto,ipoin,iplace
    nyz=ny*nz

    do ibin=1,nbin+1
       lbin2(ibin)=0_ip
    enddo

    do ipoin=9,npoin
       ix=floor((coor(1,ipoin)-bboxbin(1,1))/dx)+1_ip 
       iy=floor((coor(2,ipoin)-bboxbin(2,1))/dy)+1_ip 
       iz=floor((coor(3,ipoin)-bboxbin(3,1))/dz)+1_ip 
       ibin=(ix-1)*nyz+(iy-1)*nz+iz+1
       lbin2(ibin)=lbin2(ibin)+1
    enddo

    lbin2(1)=1_ip
    do ibin=2,nbin+1
       lbin2(ibin)=lbin2(ibin)+lbin2(ibin-1)
    enddo

    call memrea(lbin2(nbin+1)-1_ip,memor_msh,'LBIN1','insertini',lbin1) 

    do ipoin=9,npoin
       ix=floor((coor(1,ipoin)-bboxbin(1,1))/dx)+1_ip 
       iy=floor((coor(2,ipoin)-bboxbin(2,1))/dy)+1_ip 
       iz=floor((coor(3,ipoin)-bboxbin(3,1))/dz)+1_ip 
       ibin=(ix-1)*nyz+(iy-1)*nz+iz
       iplace=lbin2(ibin)
       lbin1(iplace)=ipoin
       lbin2(ibin)=iplace+1
    enddo

    do ibin=nbin+1,2,-1
       lbin2(ibin)=lbin2(ibin-1)
    enddo
    lbin2(1)=1_ip
    !
    !     Clean up lbin1
    !
    do ibin=1,nbin
       do isto=lbin2(ibin),lbin2(ibin+1)-1
          lbin1(isto)=0_ip
       enddo
    enddo

  end subroutine insertini

end module mod_bin3d

