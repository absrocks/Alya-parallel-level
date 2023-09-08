module mod_voltol

contains

 subroutine hashface(nface,nnofa,lface,lhashf1,lhashf2,nhashf)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: nface,nnofa
    integer(ip),intent(inout) :: nhashf
    integer(ip),intent(in)    :: lface(nnofa,nface)         
    integer(ip),pointer       :: lhashf2(:),lhashf1(:,:)
    integer(ip)               :: ipmax,iface,inofa,ip1,isum,iplace,nsto,ihashf         
    integer(4)                :: istat
    !
    !     This subroutine hashes the boundary faces
    !

    !
    !     Get range
    !
    ipmax=0_ip
    do iface=1,nface
       do inofa=1,nnofa
          ip1=lface(inofa,iface)  
          if(ip1>ipmax)ipmax=ip1
       enddo
    enddo
    nhashf=3*ipmax
    !
    !     Allocate lhashf2
    ! 
    allocate(lhashf2(nhashf+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LHASHF2','hashface',lhashf2)
    do iface=1,nface
       isum=lface(1,iface)+lface(2,iface)+lface(3,iface)+1_ip
       lhashf2(isum)=lhashf2(isum)+1_ip
    enddo
    !
    !     Sum up
    !
    lhashf2(1)=1_ip
    do ihashf=2,nhashf+1
       lhashf2(ihashf)=lhashf2(ihashf)+lhashf2(ihashf-1)
    enddo
    !
    !     Allocate lhashf1
    !
    nsto=lhashf2(nhashf)  
    allocate(lhashf1(2,nsto),stat=istat)
    call memchk(zero,istat,memor_msh,'LHASHF1','hashface',lhashf1)
    !
    !     Now store
    !
    do iface=1,nface
       isum=lface(1,iface)+lface(2,iface)+lface(3,iface)
       iplace=lhashf2(isum)
       lhashf1(1,iplace)=min(lface(1,iface),lface(2,iface),lface(3,iface))
       lhashf1(2,iplace)=max(lface(1,iface),lface(2,iface),lface(3,iface))
       lhashf2(isum)=iplace+1_ip
    enddo

  end subroutine hashface

 subroutine hashedge(nedge,nnosi,ledge,lhashe1,lhashe2,nhashe)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh 
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: nedge,nnosi
    integer(ip),intent(inout) :: nhashe
    integer(ip),intent(in)    :: ledge(nnosi,nedge)         
    integer(ip),pointer       :: lhashe2(:),lhashe1(:)
    integer(ip)               :: ipmax,iedge,inosi,ip1,isum,iplace,nsto,ihashe         
    integer(4)                :: istat
    !
    !     This subroutine hashes the boundary edges
    !

    !
    !     Get range
    !
    ipmax=0_ip
    do iedge=1,nedge
       do inosi=1,nnosi
          ip1=ledge(inosi,iedge)  
          if(ip1>ipmax)ipmax=ip1
       enddo
    enddo
    nhashe=2*ipmax
    !
    !     Allocate lhashe2
    ! 
    allocate(lhashe2(nhashe+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LHASHF2','hashedge',lhashe2)
    do iedge=1,nedge
       isum=ledge(1,iedge)+ledge(2,iedge)+1_ip
       lhashe2(isum)=lhashe2(isum)+1_ip
    enddo
    !
    !     Sum up
    !
    lhashe2(1)=1_ip
    do ihashe=2,nhashe+1
       lhashe2(ihashe)=lhashe2(ihashe)+lhashe2(ihashe-1)
    enddo
    !
    !     Allocate lhashe1
    !
    nsto=lhashe2(nhashe)  
    allocate(lhashe1(nsto),stat=istat)
    call memchk(zero,istat,memor_msh,'LHASHF1','hashedge',lhashe1)
    !
    !     Now store
    !
    do iedge=1,nedge
       isum=ledge(1,iedge)+ledge(2,iedge)
       iplace=lhashe2(isum)
       lhashe1(iplace)=min(ledge(1,iedge),ledge(2,iedge))
       lhashe2(isum)=iplace+1_ip
    enddo

  end subroutine hashedge

  subroutine findin(iguess,elem,nelem,npoin,coor,pnew,ndim,ihost,nnode,eltoel,isplit,idir)
    use def_kintyp, only       :  ip,rp,lg
    use mod_mshtol, only       : orient3D
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode,iguess
    integer(ip),intent(in)     :: elem(nnode,nelem),eltoel(nnode,nelem)
    real(rp),intent(in)        :: coor(ndim,npoin),pnew(ndim)
    integer(ip),intent(inout)  :: ihost,isplit,idir
    integer(ip)                :: iter,maxiter,ip1,ip2,ip3,ip4,l
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
       ihost=0_ip 
       return
    endif
    if(iguess<0)then
       write(*,*)'Error in findin, iguess=',iguess
       stop
    endif

    ihost=iguess
    isplit=1_ip
    idir=0_ip

    do iter=1,maxiter 

       if(ihost==0)return 

       if(ihost<0 .or. ihost>nelem)then
          write(*,*)'Error findin, ihost=',ihost
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
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      isplit=2
                      idir=4_ip
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
                      idir=3_ip
                      isplit=2_ip
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      idir=1_ip
                      isplit=3_ip 
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
                      isplit=2_ip
                      idir=2_ip
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      isplit=3_ip
                      idir=2_ip
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
                      isplit=3_ip
                      idir=3_ip
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 11'
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
                      isplit=2_ip
                      idir=1_ip
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      isplit=3_ip
                      idir=4_ip 
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
                      isplit=2_ip
                      idir=5_ip 
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 12'
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
                      isplit=3_ip
                      idir=6_ip
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 13'
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
                      write(*,*)'Point too close from previously inserted point 14'
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
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      isplit=1_ip
                      idir=1_ip
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
                      isplit=1_ip
                      idir=4_ip 
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      isplit=2_ip
                      idir=4_ip
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
                      isplit=1_ip
                      idir=3_ip
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      isplit=2_ip
                      idir=5_ip 
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
                      isplit=2_ip
                      idir=1_ip
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 14'
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
                      isplit=1_ip
                      idir=2_ip
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      isplit=2_ip
                      idir=6_ip 
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
                      isplit=2_ip
                      idir=2_ip
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 14'
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
                      isplit=1_ip
                      idir=3_ip
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      write(*,*)'Point too close from previously inserted point 14'
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
                      write(*,*)'Point too close from previously inserted point 14'
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
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
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
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
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
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
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
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
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
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
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
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
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
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
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
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
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
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
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
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
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
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
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
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
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
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
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
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
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
    !
    !     Host element not found
    !
    ihost=0_ip

  end subroutine findin

  subroutine swaploc3d(ifirst,elem,nnode,nelem,coor,npoin,ndim,lelem,eltoel,rqual,ltet)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnode,npoin,ndim,ifirst
    integer(ip), intent(inout)   :: nelem,ltet(npoin)
    integer(ip), pointer         :: elem(:,:),eltoel(:,:),lelem(:)
    real(rp), pointer            :: rqual(:)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: ielem,jelem,iopt,j,istack,lstack(500),jstack,nstack,ineigh,k
    integer(ip)                  :: maxlevel,level(100),ilevel,ioptglo,ioptel
    !
    !     This subroutine optimizes locally the mesh
    !

    !
    !     Set the level of neighboring elements to swap
    !     In confined places, 4 is too low 
    !
    maxlevel=4
    !
    !     Build the stack of swapable elements
    !  
    lelem(ifirst)=1_ip
    !
    !     Initialize with the newly created element
    !
    nstack=1_ip
    lstack(1)=ifirst
    level(1)=1_ip
    istack=0_ip

    do 

       if(istack==nstack)exit
       istack=istack+1
       ielem=lstack(istack)
       ilevel=level(istack)

       if(ilevel==maxlevel)exit

       do j=1,nnode
          ineigh=eltoel(j,ielem)
          !
          !     Is there a neighbor?
          !
          if(ineigh/=0)then
             !
             !     Is the element already in the stack?
             !
             if(lelem(ineigh)==0)then 
                !
                !     Add to the stack and mark in lelem with 1     
                !     
                nstack=nstack+1
                lstack(nstack)=ineigh 
                lelem(ineigh)=1_ip
                level(nstack)=ilevel+1 
             endif
          endif
       enddo
    enddo
    !
    !     Swap them
    !
    do 

       ioptglo=0_ip

       do istack=1,nstack

          ielem=lstack(istack)
          if(lelem(ielem)==1)then
             ioptel=0_ip
             do j=1,nnode
                jelem=eltoel(j,ielem)
                if(jelem/=0)then
                   if(lelem(jelem)==1)then
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
                           npoin,iopt,j,k,lelem,rqual,ltet)
                      if(iopt==1)then
                         ioptglo=1_ip
                         ioptel=1_ip 
                         !
                         !     Mark jelem and the neighbors of ielem and jelem 
                         !     if they belong to the patch
                         !
                         ineigh=eltoel(1,ielem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(2,ielem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(3,ielem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(4,ielem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(1,jelem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(2,jelem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(3,jelem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(4,jelem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
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
                lelem(ielem)=2_ip
             endif
          endif
       enddo

       if(ioptglo==0)exit
    enddo
    !
    !     Clean up
    ! 
    do istack=1,nstack
       lelem(lstack(istack))=0_ip 
    enddo

  end subroutine swaploc3d

  subroutine swap23(iel1,iel2,elem,nnode,nelem,eltoel,coor,ndim,npoin,iopt,&
       idir1,idir2,lelem,rqual,ltet)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip), intent(in)    :: nnode,npoin,ndim,iel1,iel2,idir1,idir2
    integer(ip),pointer        :: elem(:,:),eltoel(:,:),lelem(:)
    real(rp),pointer           :: rqual(:)
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
    call memrea(nelemn,memor_msh,'RQUAL','swap23',rqual)

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

  subroutine swaploc3dc(ifirst,elem,nnode,nelem,coor,npoin,ndim,lelem,eltoel,rqual,ltet)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnode,npoin,ndim,ifirst
    integer(ip), intent(inout)   :: nelem,ltet(npoin)
    integer(ip), pointer         :: elem(:,:),eltoel(:,:),lelem(:)
    real(rp), pointer            :: rqual(:)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: ielem,jelem,iopt,j,istack,lstack(500),jstack,nstack,ineigh,k
    integer(ip)                  :: maxlevel,level(100),ilevel,ioptglo,ioptel
    !
    !     This subroutine optimizes locally the mesh
    !     This is the constrained version
    !

    !
    !     Set the level of neighboring elements to swap
    !     In confined places, 4 is too low 
    !
    maxlevel=4
    !
    !     Build the stack of swapable elements
    !  
    lelem(ifirst)=1_ip
    !
    !     Initialize with the newly created element
    !
    nstack=1_ip
    lstack(1)=ifirst
    level(1)=1_ip
    istack=0_ip

    do 

       if(istack==nstack)exit
       istack=istack+1
       ielem=lstack(istack)
       ilevel=level(istack)

       if(ilevel==maxlevel)exit

       do j=1,nnode
          ineigh=eltoel(j,ielem)
          !
          !     Is there a neighbor?
          !
          if(ineigh/=0)then
             !
             !     Is the element already in the stack?
             !
             if(lelem(ineigh)==0)then 
                !
                !     Add to the stack and mark in lelem with 1     
                !     
                nstack=nstack+1
                lstack(nstack)=ineigh 
                lelem(ineigh)=1_ip
                level(nstack)=ilevel+1 
             endif
          endif
       enddo
    enddo
    !
    !     Swap them
    !
    do 

       ioptglo=0_ip

       do istack=1,nstack

          ielem=lstack(istack)
          if(lelem(ielem)==1)then
             ioptel=0_ip
             do j=1,nnode
                jelem=eltoel(j,ielem)
                if(jelem/=0)then
                   if(lelem(jelem)==1)then
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
                           npoin,iopt,j,k,lelem,rqual,ltet)
                      if(iopt==1)then
                         ioptglo=1_ip
                         ioptel=1_ip 
                         !
                         !     Mark jelem and the neighbors of ielem and jelem 
                         !     if they belong to the patch
                         !
                         ineigh=eltoel(1,ielem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(2,ielem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(3,ielem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(4,ielem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(1,jelem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(2,jelem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(3,jelem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(4,jelem)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
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
                lelem(ielem)=2_ip
             endif
          endif
       enddo

       if(ioptglo==0)exit
    enddo
    !
    !     Clean up
    ! 
    do istack=1,nstack
       lelem(lstack(istack))=0_ip 
    enddo

  end subroutine swaploc3dc


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

  subroutine  gtqual(elem,ielem,coor,rvol,nnode,nelem,ndim,rqual,npoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnode,nelem,npoin,ndim,ielem
    integer(ip), intent(in)      :: elem(nnode,nelem)
    real(rp),intent(in)          :: coor(ndim,npoin),rvol
    real(rp),intent(inout)       :: rqual
    integer(ip)          :: ip1,ip2,ip3,ip4
    real(rp)             :: rx12,ry12,rz12,rx13,ry13,rz13,rx14,ry14,rz14
    real(rp)             :: rx23,ry23,rz23,rx24,ry24,rz24,rx34,ry34,rz34
    real(rp)             :: rsur,rsur1,rsur2,rsur3,rsur4,rl1,rl2,rl3,rl4,rl5,rl6
    real(rp)             :: rx,ry,rz,c60,c12,c05,c13,alpha,hmax

    c60=6.0d+00
    c12=12.0d+00 
    c05=0.5d+00
    c13=1.0d+00/3.0d+00
    alpha=sqrt(c60)/c12
    !
    !     Get the points
    !
    ip1=elem(1,ielem)
    ip2=elem(2,ielem)
    ip3=elem(3,ielem)
    ip4=elem(4,ielem)
    !
    !     Compute length
    !
    rx12=coor(1,ip2)-coor(1,ip1)
    ry12=coor(2,ip2)-coor(2,ip1)
    rz12=coor(3,ip2)-coor(3,ip1)
    rl1=sqrt(rx12*rx12+ry12*ry12+rz12*rz12)

    rx13=coor(1,ip3)-coor(1,ip1)
    ry13=coor(2,ip3)-coor(2,ip1)
    rz13=coor(3,ip3)-coor(3,ip1)
    rl2=sqrt(rx13*rx13+ry13*ry13+rz13*rz13)

    rx14=coor(1,ip4)-coor(1,ip1)
    ry14=coor(2,ip4)-coor(2,ip1)
    rz14=coor(3,ip4)-coor(3,ip1)
    rl3=sqrt(rx14*rx14+ry14*ry14+rz14*rz14)

    rx23=coor(1,ip3)-coor(1,ip2)
    ry23=coor(2,ip3)-coor(2,ip2)
    rz23=coor(3,ip3)-coor(3,ip2)
    rl4=sqrt(rx23*rx23+ry23*ry23+rz23*rz23)

    rx24=coor(1,ip4)-coor(1,ip2)
    ry24=coor(2,ip4)-coor(2,ip2)
    rz24=coor(3,ip4)-coor(3,ip2)
    rl5=sqrt(rx24*rx24+ry24*ry24+rz24*rz24)

    rx34=coor(1,ip4)-coor(1,ip3)
    ry34=coor(2,ip4)-coor(2,ip3)
    rz34=coor(3,ip4)-coor(3,ip3)
    rl6=sqrt(rx34*rx34+ry34*ry34+rz34*rz34)
    !
    !     Compute area
    !
    rx= ry23*rz24-rz23*ry24  
    ry=-rx23*rz24+rz23*rx24  
    rz= rx23*ry24-ry23*rx24  
    rsur1=sqrt(rx*rx+ry*ry+rz*rz)

    rx= ry14*rz13-rz14*ry13  
    ry=-rx14*rz13+rz14*rx13  
    rz= rx14*ry13-ry14*rx13  
    rsur2=sqrt(rx*rx+ry*ry+rz*rz)

    rx= ry12*rz14-rz12*ry14  
    ry=-rx12*rz14+rz12*rx14  
    rz= rx12*ry14-ry12*rx14  
    rsur3=sqrt(rx*rx+ry*ry+rz*rz)

    rx= ry13*rz12-rz13*ry12  
    ry=-rx13*rz12+rz13*rx12  
    rz= rx13*ry12-ry13*rx12  
    rsur4=sqrt(rx*rx+ry*ry+rz*rz)

    rsur=(rsur1+rsur2+rsur3+rsur4)*c05

    hmax=max(rl1,rl2,rl3,rl4,rl5,rl6)
    rqual=alpha*hmax*rsur/rvol*c13

  end subroutine gtqual

subroutine split3d(ipnew,coor,ndim,npoin,nelem,elem,nnode,ielem,isplit,idir,ltet,eltoel)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)        :: nnode,ndim,ipnew,ielem,isplit,idir,npoin
    integer(ip),intent(inout)     :: nelem
    integer(ip),pointer           :: elem(:,:),eltoel(:,:)       
    integer(ip)                   :: ielem2,idir2
    real(rp),intent(in)           :: coor(ndim,npoin)
    integer(ip),intent(inout)     :: ltet(npoin) 
    !
    !     This subroutine splits directly one, two or a ring of elements
    !

    if(isplit==1)then 

       call splvol(ipnew,coor,ndim,npoin,nelem,elem,nnode,ielem,ltet,eltoel)

    else if(isplit==2)then

       ielem2=eltoel(idir,ielem) 
       if(eltoel(1,ielem2)==ielem)then
          idir2=1_ip  
       else if(eltoel(2,ielem2)==ielem)then
          idir2=2_ip  
       else if(eltoel(3,ielem2)==ielem)then
          idir2=3_ip  
       else
          idir2=4_ip  
       endif

       call splfac(ielem,ielem2,ipnew,coor,ndim,npoin,nelem,elem,nnode,idir,idir2,ltet,eltoel)

    else

       call splring(ipnew,coor,ndim,npoin,nelem,elem,nnode,ielem,idir,ltet,eltoel)

    endif

  end subroutine split3d

  subroutine splvol(ipnew,coor,ndim,npoin,nelem,elem,nnode,ielem,ltet,eltoel)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)        :: nnode,ndim,ipnew,ielem,npoin
    integer(ip),intent(inout)     :: nelem,ltet(npoin)
    integer(ip),pointer           :: elem(:,:),eltoel(:,:)       
    integer(ip)                   :: nelemn,ineigh1,ineigh2,ineigh3,ineigh4
    integer(ip)                   :: ip1,ip2,ip3,ip4,ielem1,ielem2,ielem3,ielem4 
    real(rp)                      :: coor(ndim,npoin)

    nelemn=nelem+3_ip    
    call memrea(nelemn,memor_msh,'ELEM','split3d',elem)
    call memrea(nelemn,memor_msh,'ELTOEL','split3d',eltoel)

    ip1=elem(1,ielem)
    ip2=elem(2,ielem)
    ip3=elem(3,ielem)
    ip4=elem(4,ielem)

    ineigh1=eltoel(1,ielem)
    ineigh2=eltoel(2,ielem)
    ineigh3=eltoel(3,ielem)
    ineigh4=eltoel(4,ielem)

    ielem1=ielem
    ielem2=nelem+1_ip
    ielem3=ielem2+1_ip
    ielem4=ielem3+1_ip

    elem(1,ielem1)=ip1
    elem(2,ielem1)=ip2
    elem(3,ielem1)=ipnew
    elem(4,ielem1)=ip4

    eltoel(1,ielem1)=ielem2
    eltoel(2,ielem1)=ielem3
    eltoel(3,ielem1)=ineigh3
    eltoel(4,ielem1)=ielem4

    elem(1,ielem2)=ip2
    elem(2,ielem2)=ip3
    elem(3,ielem2)=ipnew
    elem(4,ielem2)=ip4

    eltoel(1,ielem2)=ielem3
    eltoel(2,ielem2)=ielem1
    eltoel(3,ielem2)=ineigh1
    eltoel(4,ielem2)=ielem4

    elem(1,ielem3)=ip3
    elem(2,ielem3)=ip1
    elem(3,ielem3)=ipnew
    elem(4,ielem3)=ip4

    eltoel(1,ielem3)=ielem1
    eltoel(2,ielem3)=ielem2
    eltoel(3,ielem3)=ineigh2
    eltoel(4,ielem3)=ielem4

    elem(1,ielem4)=ip1
    elem(2,ielem4)=ip2
    elem(3,ielem4)=ip3
    elem(4,ielem4)=ipnew

    eltoel(1,ielem4)=ielem2
    eltoel(2,ielem4)=ielem3
    eltoel(3,ielem4)=ielem1
    eltoel(4,ielem4)=ineigh4

    nelem=nelemn

    ltet(ip1)=ielem1
    ltet(ip2)=ielem1
    ltet(ipnew)=ielem1
    ltet(ip4)=ielem1
    ltet(ip3)=ielem2

  end subroutine splvol

  subroutine splfac(iel1,iel2,ipnew,coor,ndim,npoin,nelem,elem,nnode,idir1,idir2,ltet,eltoel)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)        :: nnode,ndim,ipnew,iel1,iel2,idir1,idir2,npoin
    integer(ip),intent(inout)     :: nelem,ltet(npoin)
    integer(ip),pointer           :: elem(:,:),eltoel(:,:)  
    real(rp),intent(in)           :: coor(ndim,npoin)     
    integer(ip)                   :: nelemn,ineigh1,ineigh2,ineigh3,ineigh4,ineigh5,ineigh6
    integer(ip)                   :: ineigh4t,ineigh5t,ineigh6t
    integer(ip)                   :: ip1,ip2,ip3,ip4,ip5 
    integer(ip)                   :: ielem1,ielem2,ielem3,ielem4,ielem5,ielem6 
    integer(ip)                   :: ltab(3,4)

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

    nelemn=nelem+4_ip    
    call memrea(nelemn,memor_msh,'ELEM','split3d',elem)
    call memrea(nelemn,memor_msh,'ELTOEL','split3d',eltoel)

    ip1=elem(idir1,iel1)
    ip2=elem(ltab(1,idir1),iel1)
    ip3=elem(ltab(2,idir1),iel1)
    ip4=elem(ltab(3,idir1),iel1)
    ip5=elem(idir2,iel2)

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
    ielem4=ielem3+1_ip
    ielem5=ielem4+1_ip
    ielem6=ielem5+1_ip

    elem(1,ielem1)=ip1
    elem(2,ielem1)=ip4
    elem(3,ielem1)=ip2
    elem(4,ielem1)=ipnew

    eltoel(1,ielem1)=ielem4
    eltoel(2,ielem1)=ielem2
    eltoel(3,ielem1)=ielem3
    eltoel(4,ielem1)=ineigh2

    elem(1,ielem2)=ip1
    elem(2,ielem2)=ip2
    elem(3,ielem2)=ip3
    elem(4,ielem2)=ipnew

    eltoel(1,ielem2)=ielem5
    eltoel(2,ielem2)=ielem3
    eltoel(3,ielem2)=ielem1
    eltoel(4,ielem2)=ineigh3

    elem(1,ielem3)=ip1
    elem(2,ielem3)=ip3
    elem(3,ielem3)=ip4
    elem(4,ielem3)=ipnew

    eltoel(1,ielem3)=ielem6
    eltoel(2,ielem3)=ielem1
    eltoel(3,ielem3)=ielem2
    eltoel(4,ielem3)=ineigh1

    elem(1,ielem4)=ip5
    elem(2,ielem4)=ip2
    elem(3,ielem4)=ip4
    elem(4,ielem4)=ipnew

    eltoel(1,ielem4)=ielem1
    eltoel(2,ielem4)=ielem6
    eltoel(3,ielem4)=ielem5
    eltoel(4,ielem4)=ineigh5

    elem(1,ielem5)=ip5
    elem(2,ielem5)=ip3
    elem(3,ielem5)=ip2
    elem(4,ielem5)=ipnew

    eltoel(1,ielem5)=ielem2
    eltoel(2,ielem5)=ielem4
    eltoel(3,ielem5)=ielem6
    eltoel(4,ielem5)=ineigh6

    elem(1,ielem5)=ip5
    elem(2,ielem5)=ip4
    elem(3,ielem5)=ip3
    elem(4,ielem5)=ipnew

    eltoel(1,ielem5)=ielem3
    eltoel(2,ielem5)=ielem5
    eltoel(3,ielem5)=ielem4
    eltoel(4,ielem5)=ineigh4

    nelem=nelemn

    ltet(ip1)=ielem1
    ltet(ip2)=ielem1
    ltet(ip3)=ielem2
    ltet(ip4)=ielem1
    ltet(ipnew)=ielem1

  end subroutine splfac

  subroutine splring(ipnew,coor,ndim,npoin,nelem,elem,nnode,ielem,iedg,ltet,eltoel)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)        :: nnode,ndim,ipnew,ielem,iedg,npoin
    integer(ip),intent(inout)     :: nelem,ltet(npoin)
    integer(ip),pointer           :: elem(:,:),eltoel(:,:)  
    real(rp),intent(in)           :: coor(ndim,npoin)     
    integer(ip),parameter         :: mstack=100_ip
    integer(ip)                   :: lstack(mstack),nstack,nelemn,ip1,ip2,ip3,ip4
    integer(ip)                   :: lneigh(2,mstack),lring(2,mstack)
    integer(ip)                   :: iel,jel,ienext,ielast,ltab(2,6),istack 

    ltab(1,1)=1_ip
    ltab(2,1)=2_ip
    ltab(1,2)=1_ip
    ltab(2,2)=3_ip
    ltab(1,3)=1_ip
    ltab(2,3)=4_ip
    ltab(1,4)=2_ip
    ltab(2,4)=3_ip
    ltab(1,5)=2_ip
    ltab(2,5)=4_ip
    ltab(1,6)=3_ip
    ltab(2,6)=4_ip
    !
    !     Get the points of the edge
    ! 
    ip1=elem(ltab(1,iedg),ielem)
    ip2=elem(ltab(2,iedg),ielem)
    !
    !     Find the ring of elements surrounding edge iedg
    !
    call findring(lstack,nstack,eltoel,nelem,nnode,iedg,elem,ielem,mstack,lneigh,lring) 
    !
    !     Allocate new elements     
    !
    nelemn=nelem+nstack
    call memrea(nelemn,memor_msh,'ELEM','split3d',elem)
    call memrea(nelemn,memor_msh,'ELTOEL','split3d',eltoel)
    !
    !     Split the ring
    !
    do istack=1,nstack
       iel=lstack(istack)
       jel=nelem+istack        
       ienext=mod(istack,nstack)+1_ip
       ielast=mod(istack-2_ip,nstack)+1_ip

       elem(1,iel)=ipnew
       elem(2,iel)=ip2
       elem(3,iel)=lring(1,istack)
       elem(4,iel)=lring(2,istack)

       ltet(lring(1,istack))=iel
       ltet(lring(2,istack))=iel

       eltoel(1,iel)=lneigh(1,istack)
       eltoel(2,iel)=jel
       eltoel(3,iel)=lstack(ienext)
       eltoel(4,iel)=lstack(ielast)

       elem(1,jel)=ip1
       elem(2,jel)=ipnew
       elem(3,jel)=lring(1,istack)
       elem(4,jel)=lring(2,istack)

       eltoel(1,jel)=iel
       eltoel(2,jel)=lneigh(2,istack)
       eltoel(3,jel)=lstack(ienext+nelem)
       eltoel(4,jel)=lstack(ielast+nelem)

    enddo

    ltet(ip1)=lstack(1)
    ltet(ip2)=nelem+1_ip

    nelem=nelemn

  end subroutine splring

  subroutine findring(lstack,nstack,eltoel,nelem,nnode,iedg,elem,ielem,mstack,lneigh,lring)
    use def_kintyp, only :  ip,rp,lg
    use def_meshin, only :  memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)        :: nnode,ielem,iedg,nelem,mstack
    integer(ip),intent(in)        :: elem(nnode,nelem),eltoel(nnode,nelem)
    integer(ip),intent(inout)     :: nstack,lneigh(2,mstack),lring(2,mstack)
    integer(ip),intent(inout)     :: lstack(mstack)
    integer(ip)                   :: ipoin,idir,iel,ienext,ltab(2,6),ipnt
    integer(ip)                   :: ltab2(2,4,4)
    !
    !     This subroutine find the elements surrounding edge iedg
    !     The edge is supposed to be INTERNAL.
    !
    ltab(1,1)=3_ip
    ltab(2,1)=4_ip
    ltab(1,2)=4_ip
    ltab(2,2)=2_ip
    ltab(1,3)=2_ip
    ltab(2,3)=3_ip
    ltab(1,4)=1_ip
    ltab(2,4)=4_ip
    ltab(1,5)=3_ip
    ltab(2,5)=1_ip
    ltab(1,6)=1_ip
    ltab(2,6)=2_ip

    ltab2(1,1,2)=3_ip
    ltab2(2,1,2)=4_ip
    ltab2(1,1,3)=4_ip
    ltab2(2,1,3)=2_ip
    ltab2(1,1,4)=2_ip
    ltab2(2,1,4)=3_ip
    ltab2(1,2,3)=1_ip
    ltab2(2,2,3)=4_ip
    ltab2(1,2,4)=3_ip
    ltab2(2,2,4)=1_ip
    ltab2(1,3,4)=1_ip
    ltab2(2,3,4)=2_ip


    lstack(1)=ielem
    nstack=1_ip
    idir=ltab(1,iedg)
    ipnt=ltab(2,iedg)
    ipoin=elem(ipnt,ielem)
    ienext=eltoel(idir,ielem)

    lring(1,1)=elem(idir,ielem)
    lring(2,1)=ipoin

    lneigh(1,1)=eltoel(ltab2(1,idir,ipnt),ielem)
    lneigh(2,1)=eltoel(ltab2(2,idir,ipnt),ielem)

    iel=ielem

    do

       if(ienext==ielem)exit

       nstack=nstack+1_ip
       lstack(nstack)=ienext

       if(elem(1,ienext)==ipoin)then
          idir=1_ip 
       else if(elem(2,ienext)==ipoin)then
          idir=2_ip 
       else if(elem(3,ienext)==ipoin)then
          idir=3_ip 
       else 
          idir=4_ip 
       endif

       if(eltoel(1,ienext)==iel)then
          ipnt=1_ip 
       else if(eltoel(2,ienext)==iel)then
          ipnt=2_ip 
       else if(eltoel(3,ienext)==iel)then
          ipnt=3_ip 
       else 
          ipnt=4_ip 
       endif

       lring(1,nstack)=elem(idir,ienext)
       lring(2,nstack)=elem(ipnt,ienext)

       lneigh(1,nstack)=eltoel(ltab2(1,idir,ipnt),ienext)
       lneigh(2,nstack)=eltoel(ltab2(2,idir,ipnt),ienext)

       iel=ienext
       ipoin=elem(ipnt,iel)
       ienext=eltoel(idir,iel)

    enddo

  end subroutine findring








end module mod_voltol
