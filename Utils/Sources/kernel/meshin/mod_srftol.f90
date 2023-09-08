module mod_srftol

  use mod_mshtol 

contains

  subroutine gthost(pnew,iguess,ihost,rnofa,ndim,nface,npoin,lface,nnofa,coor,d1,d2,d3,eltoel,&
       lelem,lsurf,isurf,ierr,rsiz)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: ndim,npoin,nface,iguess,nnofa,isurf
    integer(ip),intent(in)    :: lface(nnofa,nface)
    integer(ip),intent(in)    :: eltoel(nnofa,nface),lsurf(nface)
    integer(ip),intent(inout) :: ihost,lelem(nface),ierr
    real(rp),intent(in)       :: pnew(ndim),rnofa(ndim,nface),coor(ndim,npoin),rsiz 
    real(rp),intent(inout)    :: d1,d2,d3 
    real(rp)                  :: pproj(ndim),rface(ndim),rscal,rx,ry,rz,rscaln
    real(rp)                  :: p1(ndim),p2(ndim),p3(ndim),rl,rtol,d1t,d2t,d3t,c10,c00 
    integer(ip)               :: ip1,ip2,ip3,l,nstack,istack,ielem,j,nstack2,nstack3
    real(rp)                  :: rtl1,dtot,tolrel,rdist,rdistt,rtol2 
    real(rp)                  :: pint(3),pintt(3) 
    integer(ip)               :: ineigh,maxiter,iface,iter,ipoin,ipoint,ipointt
    integer(ip)               :: ichk,ilast,jstack,inofa,ihostf
    integer(ip), parameter    :: mstack=1000_ip
    integer(ip)               :: lstack(mstack),lstack2(mstack),lstack3(mstack)
    !
    !     This subroutine finds the host element of pnew on the old surface 
    !     Begins with iguess. Layered approach. 
    !
    tolrel=-0.1d+00
    rtol=-1.0d-02
    rtol2=-0.1d+00
    c00=0.0d+00
    c10=1.0d+00
    !
    !     First, be sure iguess is on the surface
    !
    ihostf=iguess
    ihost=iguess

    if(lsurf(ihost)/=isurf)then 

       nstack=1_ip
       lstack(1)=ihost  
       lelem(ihost)=1_ip
       istack=0_ip 

       do 
          if(istack==nstack)then
             write(*,*)'No face belonging to isurf found'
             ierr=1_ip
             return
          endif
          istack=istack+1_ip
          iface=lstack(istack) 

          do inofa=1,nnofa
             ineigh=eltoel(inofa,iface)
             if(ineigh==0)cycle
             if(lelem(ineigh)==0)then
                if(lsurf(ineigh)==isurf)then
                   ihostf=ineigh
                   ihost=ineigh
                   goto  2
                else
                   nstack=nstack+1_ip
                   lstack(nstack)=ineigh
                   lelem(ineigh)=1_ip 
                endif
             endif
          enddo

       enddo

2      continue
       do istack=1,nstack
          lelem(lstack(istack))=0_ip
       enddo

    endif
    !
    !     First layer, fastest, move in the triangulation
    !
    nstack=0_ip
    maxiter=500_ip

    do iter=1,maxiter
       !
       !     Did we stay on the surface?
       !
       if(ihost==0)then
          write(*,*)'Element out of the surface'
          stop     
       endif
       ! 
       !     Are we still on the same patch ? 
       !
       if(lsurf(ihost)/=isurf)goto 5
       !
       !     Did we come back where we were?
       !
       if(lelem(ihost)==1)goto 5
       lelem(ihost)=1_ip
       !
       !     Store ihost
       !
       nstack=nstack+1
       lstack(nstack)=ihost
       !
       !     Approximation element
       !
       ip1=lface(1,ihost)
       ip2=lface(2,ihost)
       ip3=lface(3,ihost)
       !
       !     Host face normal
       !
       rface(1)=rnofa(1,ihost)
       rface(2)=rnofa(2,ihost)
       rface(3)=rnofa(3,ihost)
       !
       !     Compute area to normalize
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)
       call orient3D(p1,p2,rface,dtot,ndim)
       !
       !     Project on the surface
       !
       rx=pnew(1)-coor(1,ip1)
       ry=pnew(2)-coor(2,ip1)
       rz=pnew(3)-coor(3,ip1)

       rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
       pproj(1)= pnew(1)-rscal*rface(1)
       pproj(2)= pnew(2)-rscal*rface(2)
       pproj(3)= pnew(3)-rscal*rface(3)
       rscal=abs(rscal)

       call random_number(rl)
       l=floor(3*rl)+1_ip

       if(l==1)then

          !
          !     First side
          !
          p1(1)=coor(1,ip2)-pproj(1)
          p1(2)=coor(2,ip2)-pproj(2)
          p1(3)=coor(3,ip2)-pproj(3)
          p2(1)=coor(1,ip3)-pproj(1)
          p2(2)=coor(2,ip3)-pproj(2)
          p2(3)=coor(3,ip3)-pproj(3)

          call orient3D(p1,p2,rface,d1,ndim)
          d1=d1/dtot 

          if(d1>rtol)then

             p1(1)=coor(1,ip3)-pproj(1)
             p1(2)=coor(2,ip3)-pproj(2)
             p1(3)=coor(3,ip3)-pproj(3)
             p2(1)=coor(1,ip1)-pproj(1)
             p2(2)=coor(2,ip1)-pproj(2)
             p2(3)=coor(3,ip1)-pproj(3)

             call orient3D(p1,p2,rface,d2,ndim)
             d2=d2/dtot 

             if(d2>rtol)then

                p1(1)=coor(1,ip1)-pproj(1)
                p1(2)=coor(2,ip1)-pproj(2)
                p1(3)=coor(3,ip1)-pproj(3)
                p2(1)=coor(1,ip2)-pproj(1)
                p2(2)=coor(2,ip2)-pproj(2)
                p2(3)=coor(3,ip2)-pproj(3)

                call orient3D(p1,p2,rface,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then
                   goto 10
                else if (d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else
                   iface=eltoel(3,ihost) 
                   call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(3,ihost)
                   else
                      goto 10
                   endif

                endif

             else if(d2<-rtol)then

                ihost=eltoel(2,ihost)

             else
                iface=eltoel(2,ihost)
                call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                if(rscal>rscaln)then 
                   ihost=eltoel(2,ihost)
                else

                   p1(1)=coor(1,ip1)-pproj(1)
                   p1(2)=coor(2,ip1)-pproj(2)
                   p1(3)=coor(3,ip1)-pproj(3)
                   p2(1)=coor(1,ip2)-pproj(1)
                   p2(2)=coor(2,ip2)-pproj(2)
                   p2(3)=coor(3,ip2)-pproj(3)

                   call orient3D(p1,p2,rface,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      goto 10
                   else if (d3<-rtol)then

                      ihost=eltoel(3,ihost)

                   else
                      iface=eltoel(3,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(3,ihost)
                      else
                         goto 10
                      endif

                   endif
                endif

             endif

          else if(d1<-rtol)then

             ihost=eltoel(1,ihost)

          else
             iface=eltoel(1,ihost)
             call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
             if(rscal>rscaln)then 
                ihost=eltoel(1,ihost)
             else

                p1(1)=coor(1,ip3)-pproj(1)
                p1(2)=coor(2,ip3)-pproj(2)
                p1(3)=coor(3,ip3)-pproj(3)
                p2(1)=coor(1,ip1)-pproj(1)
                p2(2)=coor(2,ip1)-pproj(2)
                p2(3)=coor(3,ip1)-pproj(3)
                call orient3D(p1,p2,rface,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pproj(1)
                   p1(2)=coor(2,ip1)-pproj(2)
                   p1(3)=coor(3,ip1)-pproj(3)
                   p2(1)=coor(1,ip2)-pproj(1)
                   p2(2)=coor(2,ip2)-pproj(2)
                   p2(3)=coor(3,ip2)-pproj(3)

                   call orient3D(p1,p2,rface,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      goto 10
                   else if (d3<-rtol)then

                      ihost=eltoel(3,ihost)

                   else
                      iface=eltoel(3,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(3,ihost)
                      else
                         goto 10
                      endif

                   endif

                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else
                   iface=eltoel(2,ihost)
                   call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(2,ihost)
                   else

                      p1(1)=coor(1,ip1)-pproj(1)
                      p1(2)=coor(2,ip1)-pproj(2)
                      p1(3)=coor(3,ip1)-pproj(3)
                      p2(1)=coor(1,ip2)-pproj(1)
                      p2(2)=coor(2,ip2)-pproj(2)
                      p2(3)=coor(3,ip2)-pproj(3)

                      call orient3D(p1,p2,rface,d3,ndim)
                      d3=d3/dtot 

                      if(d3>rtol)then
                         goto 10
                      else if (d3<-rtol)then

                         ihost=eltoel(3,ihost)

                      else

                         write(*,*)'Degenerate case in gthost 1'
                         stop

                      endif
                   endif

                endif
             endif

          endif

       else if(l==2)then
          !
          !     Second side
          !
          p1(1)=coor(1,ip3)-pproj(1)
          p1(2)=coor(2,ip3)-pproj(2)
          p1(3)=coor(3,ip3)-pproj(3)
          p2(1)=coor(1,ip1)-pproj(1)
          p2(2)=coor(2,ip1)-pproj(2)
          p2(3)=coor(3,ip1)-pproj(3)

          call orient3D(p1,p2,rface,d2,ndim)
          d2=d2/dtot 

          if(d2>rtol)then

             p1(1)=coor(1,ip1)-pproj(1)
             p1(2)=coor(2,ip1)-pproj(2)
             p1(3)=coor(3,ip1)-pproj(3)
             p2(1)=coor(1,ip2)-pproj(1)
             p2(2)=coor(2,ip2)-pproj(2)
             p2(3)=coor(3,ip2)-pproj(3)

             call orient3D(p1,p2,rface,d3,ndim)
             d3=d3/dtot 

             if(d3>rtol)then

                p1(1)=coor(1,ip2)-pproj(1)
                p1(2)=coor(2,ip2)-pproj(2)
                p1(3)=coor(3,ip2)-pproj(3)
                p2(1)=coor(1,ip3)-pproj(1)
                p2(2)=coor(2,ip3)-pproj(2)
                p2(3)=coor(3,ip3)-pproj(3)

                call orient3D(p1,p2,rface,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then
                   goto 10
                else if (d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else
                   iface=eltoel(1,ihost) 
                   call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(1,ihost)
                   else
                      goto 10
                   endif

                endif

             else if(d3<-rtol)then

                ihost=eltoel(3,ihost)

             else
                iface=eltoel(3,ihost)
                call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                if(rscal>rscaln)then
                   ihost=eltoel(3,ihost)
                else

                   p1(1)=coor(1,ip2)-pproj(1)
                   p1(2)=coor(2,ip2)-pproj(2)
                   p1(3)=coor(3,ip2)-pproj(3)
                   p2(1)=coor(1,ip3)-pproj(1)
                   p2(2)=coor(2,ip3)-pproj(2)
                   p2(3)=coor(3,ip3)-pproj(3)

                   call orient3D(p1,p2,rface,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      goto 10
                   else if (d1<-rtol)then

                      ihost=eltoel(1,ihost)

                   else
                      iface=eltoel(1,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(1,ihost)
                      else
                         goto 10
                      endif

                   endif

                endif

             endif

          else if(d2<-rtol)then

             ihost=eltoel(2,ihost)

          else
             iface=eltoel(2,ihost)
             call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
             if(rscal>rscaln)then
                ihost=eltoel(2,ihost)
             else

                p1(1)=coor(1,ip1)-pproj(1)
                p1(2)=coor(2,ip1)-pproj(2)
                p1(3)=coor(3,ip1)-pproj(3)
                p2(1)=coor(1,ip2)-pproj(1)
                p2(2)=coor(2,ip2)-pproj(2)
                p2(3)=coor(3,ip2)-pproj(3)

                call orient3D(p1,p2,rface,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip2)-pproj(1)
                   p1(2)=coor(2,ip2)-pproj(2)
                   p1(3)=coor(3,ip2)-pproj(3)
                   p2(1)=coor(1,ip3)-pproj(1)
                   p2(2)=coor(2,ip3)-pproj(2)
                   p2(3)=coor(3,ip3)-pproj(3)

                   call orient3D(p1,p2,rface,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      goto 10
                   else if (d1<-rtol)then

                      ihost=eltoel(1,ihost)

                   else
                      iface=eltoel(1,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(1,ihost)
                      else
                         goto 10
                      endif

                   endif

                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else
                   iface=eltoel(3,ihost)
                   call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then
                      ihost=eltoel(3,ihost)
                   else

                      p1(1)=coor(1,ip2)-pproj(1)
                      p1(2)=coor(2,ip2)-pproj(2)
                      p1(3)=coor(3,ip2)-pproj(3)
                      p2(1)=coor(1,ip3)-pproj(1)
                      p2(2)=coor(2,ip3)-pproj(2)
                      p2(3)=coor(3,ip3)-pproj(3)

                      call orient3D(p1,p2,rface,d1,ndim)
                      d1=d1/dtot 

                      if(d1>rtol)then
                         goto 10
                      else if (d1<-rtol)then

                         ihost=eltoel(1,ihost)

                      else

                         write(*,*)'Degenerate case in gthost 2'
                         stop

                      endif

                   endif

                endif

             endif

          endif

       else 
          !
          !     Third side
          !
          p1(1)=coor(1,ip1)-pproj(1)
          p1(2)=coor(2,ip1)-pproj(2)
          p1(3)=coor(3,ip1)-pproj(3)
          p2(1)=coor(1,ip2)-pproj(1)
          p2(2)=coor(2,ip2)-pproj(2)
          p2(3)=coor(3,ip2)-pproj(3)

          call orient3D(p1,p2,rface,d3,ndim)
          d3=d3/dtot 

          if(d3>rtol)then

             p1(1)=coor(1,ip2)-pproj(1)
             p1(2)=coor(2,ip2)-pproj(2)
             p1(3)=coor(3,ip2)-pproj(3)
             p2(1)=coor(1,ip3)-pproj(1)
             p2(2)=coor(2,ip3)-pproj(2)
             p2(3)=coor(3,ip3)-pproj(3)

             call orient3D(p1,p2,rface,d1,ndim)
             d1=d1/dtot 

             if(d1>rtol)then

                p1(1)=coor(1,ip3)-pproj(1)
                p1(2)=coor(2,ip3)-pproj(2)
                p1(3)=coor(3,ip3)-pproj(3)
                p2(1)=coor(1,ip1)-pproj(1)
                p2(2)=coor(2,ip1)-pproj(2)
                p2(3)=coor(3,ip1)-pproj(3)

                call orient3D(p1,p2,rface,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then
                   goto 10
                else if (d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else
                   iface=eltoel(2,ihost)
                   call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(2,ihost)
                   else
                      goto 10
                   endif

                endif

             else if(d1<-rtol)then

                ihost=eltoel(1,ihost)

             else
                iface=eltoel(1,ihost) 
                call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                if(rscal>rscaln)then 
                   ihost=eltoel(1,ihost)
                else

                   p1(1)=coor(1,ip3)-pproj(1)
                   p1(2)=coor(2,ip3)-pproj(2)
                   p1(3)=coor(3,ip3)-pproj(3)
                   p2(1)=coor(1,ip1)-pproj(1)
                   p2(2)=coor(2,ip1)-pproj(2)
                   p2(3)=coor(3,ip1)-pproj(3)

                   call orient3D(p1,p2,rface,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      goto 10
                   else if (d2<-rtol)then

                      ihost=eltoel(2,ihost)

                   else
                      iface=eltoel(2,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(2,ihost)
                      else
                         goto 10
                      endif

                   endif

                endif

             endif

          else if(d3<-rtol)then

             ihost=eltoel(3,ihost)

          else
             iface=eltoel(3,ihost)
             call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
             if(rscal>rscaln)then 
                ihost=eltoel(3,ihost)
             else

                p1(1)=coor(1,ip2)-pproj(1)
                p1(2)=coor(2,ip2)-pproj(2)
                p1(3)=coor(3,ip2)-pproj(3)
                p2(1)=coor(1,ip3)-pproj(1)
                p2(2)=coor(2,ip3)-pproj(2)
                p2(3)=coor(3,ip3)-pproj(3)

                call orient3D(p1,p2,rface,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip3)-pproj(1)
                   p1(2)=coor(2,ip3)-pproj(2)
                   p1(3)=coor(3,ip3)-pproj(3)
                   p2(1)=coor(1,ip1)-pproj(1)
                   p2(2)=coor(2,ip1)-pproj(2)
                   p2(3)=coor(3,ip1)-pproj(3)

                   call orient3D(p1,p2,rface,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      goto 10
                   else if (d2<-rtol)then

                      ihost=eltoel(2,ihost)

                   else
                      iface=eltoel(2,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(2,ihost)
                      else
                         goto 10
                      endif

                   endif

                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else
                   iface=eltoel(1,ihost)
                   call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(1,ihost)
                   else

                      p1(1)=coor(1,ip3)-pproj(1)
                      p1(2)=coor(2,ip3)-pproj(2)
                      p1(3)=coor(3,ip3)-pproj(3)
                      p2(1)=coor(1,ip1)-pproj(1)
                      p2(2)=coor(2,ip1)-pproj(2)
                      p2(3)=coor(3,ip1)-pproj(3)

                      call orient3D(p1,p2,rface,d2,ndim)
                      d2=d2/dtot 

                      if(d2>rtol)then
                         goto 10
                      else if (d2<-rtol)then

                         ihost=eltoel(2,ihost)

                      else

                         write(*,*)'Degenerate case in gthost 3'
                         stop

                      endif

                   endif

                endif

             endif

          endif

       endif

    enddo
5   continue
    !
    !     Did we reach maxiter
    !
    if(iter==maxiter)then
       write(*,*)'Error in gthost, max iter number reached'
       stop
    endif

    !
    !     Second layer if needed, brute force from ihostf
    !

    !
    !     Clean up
    !
    do istack=1,nstack
       lelem(lstack(istack))=0_ip
    enddo
    !  
    !     We could not reach the host element
    !     Brute force search beginning from ihostf and controlled by the distance to the face
    !
    istack=0_ip
    nstack=1_ip
    nstack2=1_ip
    lstack(1)=ihostf
    lstack2(1)=ihostf
    lelem(ihostf)=1_ip 
    !
    !     Compute first pnew
    ! 
    !
    !     Initialize the distance to the face ihostf
    ! 
    call distTri(lface,nface,nnofa,pnew,coor,ndim,npoin,rnofa,rdist,d1,d2,d3,ipoin,ihostf,pint)
    lstack3(1)=ipoin

    do

       if(istack==nstack)exit
       istack=istack+1

       ielem=lstack(istack)
       do j=1,nnofa
          ineigh=eltoel(j,ielem)

          if(ineigh==0)cycle
          if(lsurf(ineigh)/=isurf)cycle
          if(lelem(ineigh)/=0)cycle
          lelem(ineigh)=1_ip
          if(nstack2==mstack)then
             exit
          endif
          nstack2=nstack2+1_ip
          lstack2(nstack2)=ineigh
          !
          !     Get the distance to the face
          ! 
          call distTri(lface,nface,nnofa,pnew,coor,ndim,npoin,rnofa,rdistt,d1t,d2t,d3t,ipoint,ineigh,pintt)
          !
          !     Did we find the host element?
          !
          if(d1t>rtol2 .and. d2t>rtol2 .and. d3t>rtol2)then
             ihost=ineigh
             d1=d1t
             d2=d2t
             d3=d3t
             rdist=rdistt
             pint(1)=pintt(1)
             pint(2)=pintt(2)
             pint(3)=pintt(3)
             goto 7
          endif
          !
          !     Do we have the same configuration than the current minimum distance?
          !
          if(ipoint/=0 .and. ipoin==ipoint)then

             if(nstack==mstack)then
                write(*,*)'Overflow 2 in stack in gthost for iguess:',iguess
                ierr=1_ip
                return
             endif
             !
             !     Add to the stack
             ! 
             nstack=nstack+1
             lstack(nstack)=ineigh
             lstack3(nstack)=ipoint

          else if(rdistt<rdist)then

             if(nstack==mstack)then
                exit
             endif
             !
             !     Add to the stack
             ! 
             nstack=nstack+1
             lstack(nstack)=ineigh
             lstack3(nstack)=ipoint
             rdist=rdistt        
             ipoin=ipoint         
             ihost=ineigh         !Just for debug as not usefull
             d1=d1t               !
             d2=d2t               !
             d3=d3t               !
             pint(1)=pintt(1)
             pint(2)=pintt(2)
             pint(3)=pintt(3)

          endif

       enddo

    enddo
    !
    !     Brute force
    ! 
    call brutefface(lface,nface,nnofa,pnew,coor,ndim,npoin,rnofa,rdist,d1,d2,d3,ipoin,&
         pint,ihost,isurf,lsurf,rsiz)
    goto 7
    !
    !     The host element has still no been reached, we have the closest element---> diagnostic
    !     Check that the point is out of the surface
    !
    ichk=0_ip
    !
    !     Initialize ipoin, as in the case for which the elements are around the same point,
    !     the last element may not be the limit of the surface mesh for deciding if the point is outside 
    !     the mesh
    !    
    ipoin=lstack3(nstack)

    do istack=nstack,1,-1

       ilast=lstack(istack)

       if(lstack3(istack)==ipoin)then

          if(d1<rtol)then
             ineigh=eltoel(1,ilast)
             if(ineigh==0)then
                ichk=1_ip
             else
                if(lsurf(ineigh)/=isurf)then
                   ichk=1_ip
                endif
             endif
          endif
          if(d2<rtol)then
             ineigh=eltoel(2,ilast)
             if(ineigh==0)then
                ichk=1_ip
             else
                if(lsurf(ineigh)/=isurf)then
                   ichk=1_ip
                endif
             endif
          endif
          if(d3<rtol)then
             ineigh=eltoel(3,ilast)
             if(ineigh==0)then
                ichk=1_ip
             else
                if(lsurf(ineigh)/=isurf)then
                   ichk=1_ip
                endif
             endif
          endif

          if(ichk==1)then

             ihost=0_ip
             !
             !     Clean up lelem
             !
             do jstack=1,nstack2
                lelem(lstack2(jstack))=0_ip
             enddo
             !
             !     Go home
             !
             return

          endif

       endif

    enddo
    !
    !     We did not find an element for which the way to go to pnew is outside the surface
    !
    if(ichk==0)then
       write(*,*)'Strange configuration in gthost'
       ierr=1_ip
       !
       !     Clean up lelem
       !
       do jstack=1,nstack2
          lelem(lstack2(jstack))=0_ip
       enddo
       return
    endif

7   continue
    !
    !     Clean up lelem
    !
    do istack=1,nstack2
       lelem(lstack2(istack))=0_ip
    enddo
    !
    !     Set nstack to 0 for lelem clean up 
    !
    nstack=0_ip

10  continue
    !
    !     Last part: we were successfull
    !
    !                - on input:  a good ihost
    !                - on output: a better ihost
    !
    !     Clean up
    !
    do istack=1,nstack
       lelem(lstack(istack))=0_ip
    enddo
    !
    !     Try to find better
    !
    istack=0_ip
    nstack=1_ip
    lstack(1)=ihost
    nstack2=1_ip
    lstack2(1)=ihost 
    lelem(ihost)=1_ip
    !
    !     Initialize the distance to the face ihost
    ! 
    call distTri(lface,nface,nnofa,pnew,coor,ndim,npoin,rnofa,rdist,d1,d2,d3,ipoin,ihost,pint)

    do         
       if(istack==nstack)exit
       istack=istack+1

       ielem=lstack(istack)
       do j=1,nnofa
          ineigh=eltoel(j,ielem)
          if(ineigh==0)cycle
          if(lsurf(ineigh)/=isurf)cycle
          if(lelem(ineigh)/=0)cycle
          lelem(ineigh)=1_ip
          if(nstack2==mstack)then
             write(*,*)'Overflow 4 in stack in gthost for iguess:',iguess
             ierr=1_ip
             return
          endif
          nstack2=nstack2+1_ip
          lstack2(nstack2)=ineigh
          !
          !     Get the distance to the face
          ! 
          call distTri(lface,nface,nnofa,pnew,coor,ndim,npoin,rnofa,rdistt,d1t,d2t,d3t,ipoint,ineigh,pintt)

          if(ipoint/=0 .and. ipoin==ipoint)then

             if(nstack==mstack)then
                write(*,*)'Overflow 5 in stack in gthost for iguess:',iguess
                ierr=1_ip
                return
             endif
             !
             !     Add to the stack
             ! 
             nstack=nstack+1
             lstack(nstack)=ineigh

          else if(rdistt<rdist)then
             !
             !     Check the shape function of the projected point
             !
             if(d1t>rtol .and. d2t>rtol .and. d3t>rtol)then
                ihost=ineigh
                d1=d1t
                d2=d2t
                d3=d3t
                rdist=rdistt
                pint(1)=pintt(1)
                pint(2)=pintt(2)
                pint(3)=pintt(3)
             endif

             if(nstack==mstack)then
                write(*,*)'Overflow 6 in stack in gthost for iguess:',iguess
                ierr=1_ip
                return
             endif
             !
             !     Add to the stack
             ! 
             nstack=nstack+1
             lstack(nstack)=ineigh
             rdist=rdistt        
             ipoin=ipoint         
             ihost=ineigh         !Just for debug as not usefull
             d1=d1t               !
             d2=d2t               !
             d3=d3t               !
             pint(1)=pintt(1)
             pint(2)=pintt(2)
             pint(3)=pintt(3)

          endif
       enddo
    enddo
    !
    !     Clean up lelem
    !
    do istack=1,nstack2
       lelem(lstack2(istack))=0_ip
    enddo
    !
    !     Recompute shape functions for points outside
    !
    if(d1<rtol .or. d2<rtol .or. d3<rtol)then

       ip1=lface(1,ihost)
       ip2=lface(2,ihost)
       ip3=lface(3,ihost)
       !
       !     Compute area to normalize
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)
       call orient3D(p1,p2,rface,dtot,ndim)
       !
       !     Get shape functions
       !
       p1(1)=coor(1,ip2)-pint(1)
       p1(2)=coor(2,ip2)-pint(2)
       p1(3)=coor(3,ip2)-pint(3)
       p2(1)=coor(1,ip3)-pint(1)
       p2(2)=coor(2,ip3)-pint(2)
       p2(3)=coor(3,ip3)-pint(3)

       call orient3D(p1,p2,rface,d1,ndim)
       d1=d1/dtot 

       p1(1)=coor(1,ip3)-pint(1)
       p1(2)=coor(2,ip3)-pint(2)
       p1(3)=coor(3,ip3)-pint(3)
       p2(1)=coor(1,ip1)-pint(1)
       p2(2)=coor(2,ip1)-pint(2)
       p2(3)=coor(3,ip1)-pint(3)

       call orient3D(p1,p2,rface,d2,ndim)
       d2=d2/dtot 

       p1(1)=coor(1,ip1)-pint(1)
       p1(2)=coor(2,ip1)-pint(2)
       p1(3)=coor(3,ip1)-pint(3)
       p2(1)=coor(1,ip2)-pint(1)
       p2(2)=coor(2,ip2)-pint(2)
       p2(3)=coor(3,ip2)-pint(3)

       call orient3D(p1,p2,rface,d3,ndim)
       d3=d3/dtot 

    endif


  end subroutine gthost

  subroutine gthost2(ipa,ipb,pnew,iguess,ihost,rnofa,ndim,nface,npoin,lface,nnofa,&
       coor,d1,d2,d3,eltoel,rsize,coornew,npnew,lsurf,isurf,lelem,ierr,lptri,&
       lfnew,nfnew,eltoelnew)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)              :: ndim,npoin,nface,iguess,nnofa,npnew,isurf,nfnew
    integer(ip),intent(in)              :: lface(nnofa,nface),lptri(npoin),lfnew(nnofa,nfnew)
    integer(ip),intent(in)              :: eltoel(nnofa,nface),ipa,ipb,lsurf(nface),eltoelnew(nnofa,nfnew)
    integer(ip),intent(inout)           :: ihost,lelem(nface),ierr
    real(rp), intent(in)                :: rsize(npnew)
    real(rp),intent(in)                 :: rnofa(ndim,nface),coor(ndim,npoin),coornew(ndim,npnew) 
    real(rp),intent(inout)              :: pnew(ndim) 
    real(rp),intent(inout)              :: d1,d2,d3
    real(rp)                            :: pproj(ndim),rface(ndim),rscal,rx,ry,rz,rscaln,rdistt,rdist
    real(rp)                            :: p1(ndim),p2(ndim),p3(ndim),rl,rtol,d1t,d2t,d3t,c05,c10 
    integer(ip)                         :: ip1,ip2,ip3,l,nstack,nstack2,istack,jstack,ielem,j,rtol2
    integer(ip)                         :: ineigh,iface,ilast,ichk,niter,ihostf,inofa
    real(rp)                            :: c1,c2,det,d1a,d2a,pmid(ndim),rtx,rty,rtz,epsil,rtl,a,b,c,rtl1
    real(rp)                            :: cscal,x1,y1,z1,rnew,ux,uy,uz,rkl,ul,t1,fact,dtot,tolrel,c00
    real(rp)                            :: tolratio,ratio,pint(3),pintt(3)
    integer(ip),parameter               :: mstack=1000
    integer(ip)                         :: lstack(mstack),lstack2(mstack),lstack3(mstack),ipoin,ipoint
    integer(ip)                         :: iturn,iter,nitermax,iview1,ienext1,ienext2,ipnew,ie
    real(rp)                            :: pnewt(ndim)
    real(rp)                            :: vmin(3),vmax(3),rlen,c20
    real(rp)                            :: xmax,ymax,zmax,xmin,ymin,zmin,coef 
    integer(ip)                         :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))


    c00=0.0d+00
    c05=0.5d+00
    c10=1.0d+00
    c20=2.0d+00
    epsil=1.0d-12
    fact=sqrt(3.0d+00)/2.0d+00
    rtol=-1.0d-02
    rtol2=-0.1d+00
    tolratio=1.0d+01
    nitermax=100
    coef=10.0d+00
    !tolrel=-0.1d+00 not used
    !
    !     DBG
    !  
    !do iface=1,nface
    !   if(lelem(iface)/=0)then
    !      write(*,*)'lelem not clean in gthost2'
    !   endif
    !enddo
    !
    !     Get the side vector
    !
    rtx=coornew(1,ipa)-coornew(1,ipb)
    rty=coornew(2,ipa)-coornew(2,ipb)
    rtz=coornew(3,ipa)-coornew(3,ipb)
    rtl1=sqrt(rtx*rtx+rty*rty+rtz*rtz)
    rtl=c10/rtl1
    rtx=rtx*rtl
    rty=rty*rtl
    rtz=rtz*rtl
    !
    !     Get the size
    !
    !rnew=c05*(c05*(rsize(ipa)+rsize(ipb))+rtl1)*fact
    rnew=(c05*(rsize(ipa)+rsize(ipb)))*fact
    !
    !     Compare size with side vector
    !
    ratio=rnew*rtl
    if(max(ratio,c10/ratio)>tolratio)then
       write(*,*)'Warning: bad size distribution ipa:',ipa,'ipb:',ipb
       !rnew=c20*min(rsize(rnew),rtl1)

       ie=lptri(ipa)
       iturn=1_ip
       if(ie<=0 .or. ie>nface)then
          write(*,*)'Error in gthost2 1'
          write(*,*)ie,nface
          stop
       endif

       iter=0
       do

          iter=iter+1
          if(iter==nitermax)then
             write(*,*)'Error in gthost2 finding ipa,ipb'
             ierr=1
             return
          endif
          !
          !     Find element containing (ipa,ipb)
          !
          if(lfnew(1,ie)==ipa)then
             ienext1=eltoelnew(2,ie)
             ienext2=eltoelnew(3,ie)
             iview1=1_ip
          else if(lfnew(2,ie)==ipa)then
             ienext1=eltoelnew(3,ie)
             ienext2=eltoelnew(1,ie)
             iview1=2_ip
          else
             ienext1=eltoelnew(1,ie)
             ienext2=eltoelnew(2,ie)
             iview1=3_ip
          endif

          if(ipb==lfnew(ltab(1,iview1),ie))exit

          if(iturn==1)then

             if(ienext1/=0)then
                ie=ienext1
             else
                ie=ienext2
                iturn=0_ip
             endif
          else

             ie=ienext2

          endif

       enddo
       !
       !     Get the new point as the third point
       !
       if(lfnew(1,ie)/=ipa .and. lfnew(1,ie)/=ipb)then

          ipnew=lfnew(1,ie)
          d1=c10
          d2=c00
          d3=c00

       else if(lfnew(2,ie)/=ipa .and. lfnew(2,ie)/=ipb)then

          ipnew=lfnew(1,ie)
          d1=c00
          d2=c10
          d3=c00

       else

          ipnew=lfnew(3,ie)
          d1=c00
          d2=c00
          d3=c10

       endif

       pnew(1)=coor(1,ipnew)
       pnew(2)=coor(2,ipnew)
       pnew(3)=coor(3,ipnew)
       ihost=ie

    endif
    !
    !     Get the mid point
    !
    pmid(1)=(coornew(1,ipa)+coornew(1,ipb))*c05
    pmid(2)=(coornew(2,ipa)+coornew(2,ipb))*c05
    pmid(3)=(coornew(3,ipa)+coornew(3,ipb))*c05
    !write(*,*)pmid(1),pmid(2),pmid(3)
    !
    !     First  part: use a random walk from iguess
    !                  on input:  iguess
    !                  on output: ihost 
    !                             goto 10 if successfull 
    !                             goto 5 if not
    ! 
    ! 
    !     Find host on the old surface
    !
    ihost=iguess
    ihostf=iguess
    !
    !     First, be sure iguess is on the surface
    !
    if(lsurf(ihost)/=isurf)then 

       nstack=1_ip
       lstack(1)=ihost  
       lelem(ihost)=1_ip
       istack=0_ip 

       do 
          if(istack==nstack)then
             write(*,*)'No face belonging to isurf found'
             ierr=1_ip
             return
          endif
          istack=istack+1_ip
          iface=lstack(istack) 

          do inofa=1,nnofa
             ineigh=eltoel(inofa,iface)
             if(ineigh==0)cycle
             if(lelem(ineigh)==0)then
                if(lsurf(ineigh)==isurf)then
                   ihostf=ineigh
                   ihost=ineigh
                   goto  2
                else
                   nstack=nstack+1
                   lstack(nstack)=ineigh
                   lelem(ineigh)=1
                endif
             endif
          enddo

       enddo

2      continue

       do istack=1,nstack
          lelem(lstack(istack))=0_ip
       enddo

    endif

    nstack=0_ip

    do
       !
       !     Did we stay on the surface ?
       !
       if(ihost==0)goto 5
       ! 
       !     Are we still on the same patch ? 
       !
       if(lsurf(ihost)/=isurf)goto 5
       ! 
       !     Did we come back where we were?
       !
       if(lelem(ihost)==1)goto 5
       lelem(ihost)=1_ip
       !
       !     Store ihost
       !
       nstack=nstack+1
       lstack(nstack)=ihost
       !
       !     Compute the new point position
       !
       call newpnt(nface,npoin,lface,ihost,nnofa,ndim,coor,rnofa,pnew,&
            pmid,rtx,rty,rtz,ichk,rnew,ierr)
       if(ierr==1)return

       if(ichk==-1)goto 5 
       !
       !     Get the points of the face
       !  
       ip1=lface(1,ihost)
       ip2=lface(2,ihost)
       ip3=lface(3,ihost)
       !
       !     Get the face normal
       !  
       rface(1)=rnofa(1,ihost)
       rface(2)=rnofa(2,ihost)
       rface(3)=rnofa(3,ihost)
       !
       !     Project on the surface
       !
       rx=pnew(1)-coor(1,ip1)
       ry=pnew(2)-coor(2,ip1)
       rz=pnew(3)-coor(3,ip1)

       rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
       pproj(1)= pnew(1)-rscal*rface(1)
       pproj(2)= pnew(2)-rscal*rface(2)
       pproj(3)= pnew(3)-rscal*rface(3)
       rscal=abs(rscal)
       !
       !     Compute area to normalize
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)
       call orient3D(p1,p2,rface,dtot,ndim)

       call random_number(rl)
       l=floor(3*rl)+1_ip

       if(l==1)then

          !
          !     First side
          !
          p1(1)=coor(1,ip2)-pproj(1)
          p1(2)=coor(2,ip2)-pproj(2)
          p1(3)=coor(3,ip2)-pproj(3)
          p2(1)=coor(1,ip3)-pproj(1)
          p2(2)=coor(2,ip3)-pproj(2)
          p2(3)=coor(3,ip3)-pproj(3)

          call orient3D(p1,p2,rface,d1,ndim)
          d1=d1/dtot
          if(d1>rtol)then

             p1(1)=coor(1,ip3)-pproj(1)
             p1(2)=coor(2,ip3)-pproj(2)
             p1(3)=coor(3,ip3)-pproj(3)
             p2(1)=coor(1,ip1)-pproj(1)
             p2(2)=coor(2,ip1)-pproj(2)
             p2(3)=coor(3,ip1)-pproj(3)

             call orient3D(p1,p2,rface,d2,ndim)
             d2=d2/dtot

             if(d2>rtol)then

                p1(1)=coor(1,ip1)-pproj(1)
                p1(2)=coor(2,ip1)-pproj(2)
                p1(3)=coor(3,ip1)-pproj(3)
                p2(1)=coor(1,ip2)-pproj(1)
                p2(2)=coor(2,ip2)-pproj(2)
                p2(3)=coor(3,ip2)-pproj(3)

                call orient3D(p1,p2,rface,d3,ndim)
                d3=d3/dtot

                if(d3>rtol)then
                   goto 10
                else if (d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else
                   iface=eltoel(3,ihost)
                   call gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(3,ihost)
                   else
                      goto 10
                   endif

                endif

             else if(d2<-rtol)then

                ihost=eltoel(2,ihost)

             else

                iface=eltoel(2,ihost)
                call gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                if(rscal>rscaln)then
                   ihost=eltoel(2,ihost)
                else

                   p1(1)=coor(1,ip1)-pproj(1)
                   p1(2)=coor(2,ip1)-pproj(2)
                   p1(3)=coor(3,ip1)-pproj(3)
                   p2(1)=coor(1,ip2)-pproj(1)
                   p2(2)=coor(2,ip2)-pproj(2)
                   p2(3)=coor(3,ip2)-pproj(3)

                   call orient3D(p1,p2,rface,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      goto 10
                   else if (d3<-rtol)then
                      ihost=eltoel(3,ihost)
                   else
                      iface=eltoel(3,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(3,ihost)
                      else
                         goto 10
                      endif

                   endif
                endif

             endif

          else if(d1<-rtol)then

             ihost=eltoel(1,ihost)

          else
             iface=eltoel(1,ihost)
             call gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
             if(rscal>rscaln)then 
                ihost=eltoel(1,ihost)
             else 
                p1(1)=coor(1,ip3)-pproj(1)
                p1(2)=coor(2,ip3)-pproj(2)
                p1(3)=coor(3,ip3)-pproj(3)
                p2(1)=coor(1,ip1)-pproj(1)
                p2(2)=coor(2,ip1)-pproj(2)
                p2(3)=coor(3,ip1)-pproj(3)
                call orient3D(p1,p2,rface,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pproj(1)
                   p1(2)=coor(2,ip1)-pproj(2)
                   p1(3)=coor(3,ip1)-pproj(3)
                   p2(1)=coor(1,ip2)-pproj(1)
                   p2(2)=coor(2,ip2)-pproj(2)
                   p2(3)=coor(3,ip2)-pproj(3)

                   call orient3D(p1,p2,rface,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      goto 10
                   else if (d3<-rtol)then

                      ihost=eltoel(3,ihost)

                   else
                      iface=eltoel(3,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(3,ihost)
                      else
                         goto 10
                      endif

                   endif

                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else
                   iface=eltoel(2,ihost)
                   call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(2,ihost)
                   else

                      p1(1)=coor(1,ip1)-pproj(1)
                      p1(2)=coor(2,ip1)-pproj(2)
                      p1(3)=coor(3,ip1)-pproj(3)
                      p2(1)=coor(1,ip2)-pproj(1)
                      p2(2)=coor(2,ip2)-pproj(2)
                      p2(3)=coor(3,ip2)-pproj(3)

                      call orient3D(p1,p2,rface,d3,ndim)
                      d3=d3/dtot 

                      if(d3>rtol)then
                         goto 10
                      else if (d3<-rtol)then

                         ihost=eltoel(3,ihost)

                      else

                         write(*,*)'Degenerate case in gthost2 1'
                         stop

                      endif
                   endif

                endif
             endif

          endif

       else if(l==2)then

          !
          !     Second side
          !
          p1(1)=coor(1,ip3)-pproj(1)
          p1(2)=coor(2,ip3)-pproj(2)
          p1(3)=coor(3,ip3)-pproj(3)
          p2(1)=coor(1,ip1)-pproj(1)
          p2(2)=coor(2,ip1)-pproj(2)
          p2(3)=coor(3,ip1)-pproj(3)

          call orient3D(p1,p2,rface,d2,ndim)
          d2=d2/dtot

          if(d2>rtol)then

             p1(1)=coor(1,ip1)-pproj(1)
             p1(2)=coor(2,ip1)-pproj(2)
             p1(3)=coor(3,ip1)-pproj(3)
             p2(1)=coor(1,ip2)-pproj(1)
             p2(2)=coor(2,ip2)-pproj(2)
             p2(3)=coor(3,ip2)-pproj(3)

             call orient3D(p1,p2,rface,d3,ndim)
             d3=d3/dtot

             if(d3>rtol)then

                p1(1)=coor(1,ip2)-pproj(1)
                p1(2)=coor(2,ip2)-pproj(2)
                p1(3)=coor(3,ip2)-pproj(3)
                p2(1)=coor(1,ip3)-pproj(1)
                p2(2)=coor(2,ip3)-pproj(2)
                p2(3)=coor(3,ip3)-pproj(3)

                call orient3D(p1,p2,rface,d1,ndim)
                d1=d1/dtot

                if(d1>rtol)then
                   goto 10
                else if (d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else
                   iface=eltoel(1,ihost)
                   call gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(1,ihost)
                   else
                      goto 10
                   endif

                endif

             else if(d3<-rtol)then

                ihost=eltoel(3,ihost)

             else
                iface=eltoel(3,ihost)
                call gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                if(rscal>rscaln)then 
                   ihost=eltoel(3,ihost)
                else
                   p1(1)=coor(1,ip2)-pproj(1)
                   p1(2)=coor(2,ip2)-pproj(2)
                   p1(3)=coor(3,ip2)-pproj(3)
                   p2(1)=coor(1,ip3)-pproj(1)
                   p2(2)=coor(2,ip3)-pproj(2)
                   p2(3)=coor(3,ip3)-pproj(3)

                   call orient3D(p1,p2,rface,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      goto 10
                   else if (d1<-rtol)then

                      ihost=eltoel(1,ihost)

                   else
                      iface=eltoel(1,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(1,ihost)
                      else
                         goto 10
                      endif

                   endif
                endif

             endif

          else if(d2<-rtol)then

             ihost=eltoel(2,ihost)

          else
             iface=eltoel(2,ihost)
             call gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
             if(rscal>rscaln)then 
                ihost=eltoel(2,ihost)
             else

                p1(1)=coor(1,ip1)-pproj(1)
                p1(2)=coor(2,ip1)-pproj(2)
                p1(3)=coor(3,ip1)-pproj(3)
                p2(1)=coor(1,ip2)-pproj(1)
                p2(2)=coor(2,ip2)-pproj(2)
                p2(3)=coor(3,ip2)-pproj(3)

                call orient3D(p1,p2,rface,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip2)-pproj(1)
                   p1(2)=coor(2,ip2)-pproj(2)
                   p1(3)=coor(3,ip2)-pproj(3)
                   p2(1)=coor(1,ip3)-pproj(1)
                   p2(2)=coor(2,ip3)-pproj(2)
                   p2(3)=coor(3,ip3)-pproj(3)

                   call orient3D(p1,p2,rface,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      goto 10
                   else if (d1<-rtol)then

                      ihost=eltoel(1,ihost)

                   else
                      iface=eltoel(1,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(1,ihost)
                      else
                         goto 10
                      endif

                   endif

                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else
                   iface=eltoel(3,ihost)
                   call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then
                      ihost=eltoel(3,ihost)
                   else
                      p1(1)=coor(1,ip2)-pproj(1)
                      p1(2)=coor(2,ip2)-pproj(2)
                      p1(3)=coor(3,ip2)-pproj(3)
                      p2(1)=coor(1,ip3)-pproj(1)
                      p2(2)=coor(2,ip3)-pproj(2)
                      p2(3)=coor(3,ip3)-pproj(3)

                      call orient3D(p1,p2,rface,d1,ndim)
                      d1=d1/dtot 

                      if(d1>rtol)then
                         goto 10
                      else if (d1<-rtol)then

                         ihost=eltoel(1,ihost)

                      else

                         write(*,*)'Degenerate case in gthost2 2'
                         stop

                      endif

                   endif

                endif
             endif

          endif

       else 
          !
          !     Third side
          !
          p1(1)=coor(1,ip1)-pproj(1)
          p1(2)=coor(2,ip1)-pproj(2)
          p1(3)=coor(3,ip1)-pproj(3)
          p2(1)=coor(1,ip2)-pproj(1)
          p2(2)=coor(2,ip2)-pproj(2)
          p2(3)=coor(3,ip2)-pproj(3)

          call orient3D(p1,p2,rface,d3,ndim)
          d3=d3/dtot

          if(d3>rtol)then

             p1(1)=coor(1,ip2)-pproj(1)
             p1(2)=coor(2,ip2)-pproj(2)
             p1(3)=coor(3,ip2)-pproj(3)
             p2(1)=coor(1,ip3)-pproj(1)
             p2(2)=coor(2,ip3)-pproj(2)
             p2(3)=coor(3,ip3)-pproj(3)

             call orient3D(p1,p2,rface,d1,ndim)
             d1=d1/dtot

             if(d1>rtol)then

                p1(1)=coor(1,ip3)-pproj(1)
                p1(2)=coor(2,ip3)-pproj(2)
                p1(3)=coor(3,ip3)-pproj(3)
                p2(1)=coor(1,ip1)-pproj(1)
                p2(2)=coor(2,ip1)-pproj(2)
                p2(3)=coor(3,ip1)-pproj(3)

                call orient3D(p1,p2,rface,d2,ndim)
                d2=d2/dtot

                if(d2>rtol)then
                   goto 10
                else if (d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else
                   iface=eltoel(2,ihost)
                   call gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(2,ihost)
                   else
                      goto 10
                   endif

                endif

             else if(d1<-rtol)then

                ihost=eltoel(1,ihost)

             else
                iface=eltoel(1,ihost)
                call gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                if(rscal>rscaln)then 
                   ihost=eltoel(1,ihost)
                else

                   p1(1)=coor(1,ip3)-pproj(1)
                   p1(2)=coor(2,ip3)-pproj(2)
                   p1(3)=coor(3,ip3)-pproj(3)
                   p2(1)=coor(1,ip1)-pproj(1)
                   p2(2)=coor(2,ip1)-pproj(2)
                   p2(3)=coor(3,ip1)-pproj(3)

                   call orient3D(p1,p2,rface,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      goto 10
                   else if (d2<-rtol)then

                      ihost=eltoel(2,ihost)

                   else
                      iface=eltoel(2,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(2,ihost)
                      else
                         goto 10
                      endif

                   endif
                endif

             endif

          else if(d3<-rtol)then

             ihost=eltoel(3,ihost)

          else
             iface=eltoel(3,ihost)
             call gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
             if(rscal>rscaln)then 
                ihost=eltoel(3,ihost)
             else

                p1(1)=coor(1,ip2)-pproj(1)
                p1(2)=coor(2,ip2)-pproj(2)
                p1(3)=coor(3,ip2)-pproj(3)
                p2(1)=coor(1,ip3)-pproj(1)
                p2(2)=coor(2,ip3)-pproj(2)
                p2(3)=coor(3,ip3)-pproj(3)

                call orient3D(p1,p2,rface,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip3)-pproj(1)
                   p1(2)=coor(2,ip3)-pproj(2)
                   p1(3)=coor(3,ip3)-pproj(3)
                   p2(1)=coor(1,ip1)-pproj(1)
                   p2(2)=coor(2,ip1)-pproj(2)
                   p2(3)=coor(3,ip1)-pproj(3)

                   call orient3D(p1,p2,rface,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      goto 10
                   else if (d2<-rtol)then

                      ihost=eltoel(2,ihost)

                   else
                      iface=eltoel(2,ihost)
                      call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                      if(rscal>rscaln)then 
                         ihost=eltoel(2,ihost)
                      else
                         goto 10
                      endif

                   endif

                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else
                   iface=eltoel(1,ihost)
                   call  gtbscal(lface,iface,pnew,coor,rscaln,ndim,npoin,rnofa,nface,nnofa)
                   if(rscal>rscaln)then 
                      ihost=eltoel(1,ihost)
                   else

                      p1(1)=coor(1,ip3)-pproj(1)
                      p1(2)=coor(2,ip3)-pproj(2)
                      p1(3)=coor(3,ip3)-pproj(3)
                      p2(1)=coor(1,ip1)-pproj(1)
                      p2(2)=coor(2,ip1)-pproj(2)
                      p2(3)=coor(3,ip1)-pproj(3)

                      call orient3D(p1,p2,rface,d2,ndim)
                      d2=d2/dtot 
                      if(d2>rtol)then
                         goto 10
                      else if (d2<-rtol)then

                         ihost=eltoel(2,ihost)

                      else

                         write(*,*)'Degenerate case in gthost2 3'
                         stop

                      endif

                   endif

                endif
             endif

          endif

       endif

    enddo

5   continue
    !
    !     Clean up
    !
    do istack=1,nstack
       lelem(lstack(istack))=0_ip
    enddo
    !  
    !     We could not reach the host element
    !     Brute force search beginning from ihostf and controlled by the distance to the face
    !
    istack=0_ip
    nstack=1_ip
    nstack2=1_ip
    lstack(1)=ihostf
    lstack2(1)=ihostf
    lelem(ihostf)=1_ip 
    !
    !     Compute first pnew
    ! 
    !
    !     Compute the new point position
    !
    call newpnt(nface,npoin,lface,ihostf,nnofa,ndim,coor,rnofa,pnew,pmid,&
         rtx,rty,rtz,ichk,rnew,ierr)
    if(ierr==1)return
    if(ichk==-1)then

       ipoin=-1   
       rdist=1.0d+32 

    else
       !
       !     Initialize the distance to the face ihostf
       ! 
       call distTri(lface,nface,nnofa,pnew,coor,ndim,npoin,rnofa,rdist,d1,d2,d3,ipoin,ihostf,pint)

    endif

    lstack3(1)=ipoin

    do

       if(istack==nstack)exit
       istack=istack+1

       ielem=lstack(istack)
       do j=1,nnofa
          ineigh=eltoel(j,ielem)

          if(ineigh==0)cycle
          if(lsurf(ineigh)/=isurf)cycle
          if(lelem(ineigh)/=0)cycle
          lelem(ineigh)=1_ip
          if(nstack2==mstack)then
             write(*,*)'Overflow in stack in gthost2 for iguess:',iguess
             ierr=1_ip
             return
          endif
          nstack2=nstack2+1_ip
          lstack2(nstack2)=ineigh
          !
          !     Compute the new point position
          !
          call newpnt(nface,npoin,lface,ineigh,nnofa,ndim,coor,rnofa,pnewt,&
               pmid,rtx,rty,rtz,ichk,rnew,ierr)
          if(ierr==1)return
          if(ichk==-1)cycle
          !
          !     Get the distance to the face
          ! 
          call distTri(lface,nface,nnofa,pnewt,coor,ndim,npoin,rnofa,rdistt,d1t,d2t,d3t,ipoint,ineigh,pintt)
          !
          !     Did we find the host element?
          !
          if(d1t>rtol2 .and. d2t>rtol2 .and. d3t>rtol2)then
             ihost=ineigh
             d1=d1t
             d2=d2t
             d3=d3t
             rdist=rdistt
             pnew(1)=pnewt(1)
             pnew(2)=pnewt(2)
             pnew(3)=pnewt(3)
             pint(1)=pintt(1)
             pint(2)=pintt(2)
             pint(3)=pintt(3)
             goto 7
          endif
          !
          !     Do we have the same configuration than the current minimum distance?
          !
          if(ipoint/=0 .and. ipoin==ipoint)then

             if(nstack==mstack)then
                write(*,*)'Overflow in stack in gthost2 for iguess:',iguess
                ierr=1_ip
                return
             endif
             !
             !     Add to the stack
             ! 
             nstack=nstack+1
             lstack(nstack)=ineigh
             lstack3(nstack)=ipoint

          else if(rdistt<rdist)then

             if(nstack==mstack)then
                write(*,*)'Overflow in stack in gthost2 for iguess:',iguess
                ierr=1_ip
                return
             endif
             !
             !     Add to the stack
             ! 
             nstack=nstack+1
             lstack(nstack)=ineigh
             lstack3(nstack)=ipoint
             rdist=rdistt        
             ipoin=ipoint         
             ihost=ineigh         !Just for debug as not usefull
             d1=d1t               !
             d2=d2t               !
             d3=d3t               !
             pnew(1)=pnewt(1)
             pnew(2)=pnewt(2)
             pnew(3)=pnewt(3)
             pint(1)=pintt(1)
             pint(2)=pintt(2)
             pint(3)=pintt(3)

          endif

       enddo

    enddo
    !
    !     Brute force
    !
    rdist=1.0d+12
    rlen=coef*rnew
    vmin(1)=pmid(1)-rlen
    vmin(2)=pmid(2)-rlen
    vmin(3)=pmid(3)-rlen
    vmax(1)=pmid(1)+rlen
    vmax(2)=pmid(2)+rlen
    vmax(3)=pmid(3)+rlen

    do iface=1,nface

       if(lsurf(iface)==isurf)then 
          !
          !     Get the points of the face
          !
          ip1=lface(1,iface)
          ip2=lface(2,iface)
          ip3=lface(3,iface)
          !
          !     Host face normal
          !
          rface(1)=rnofa(1,iface)
          rface(2)=rnofa(2,iface)
          rface(3)=rnofa(3,iface)
          !
          !     Filter with bbox
          !
          xmin=coor(1,ip1)
          xmax=xmin
          ymin=coor(2,ip1)
          ymax=ymin
          zmin=coor(3,ip1)
          zmax=zmin

          if(coor(1,ip2)<xmin)xmin=coor(1,ip2)
          if(coor(2,ip2)<ymin)ymin=coor(2,ip2)
          if(coor(3,ip2)<zmin)zmin=coor(3,ip2)
          if(coor(1,ip2)>xmax)xmax=coor(1,ip2)
          if(coor(2,ip2)>ymax)ymax=coor(2,ip2)
          if(coor(3,ip2)>zmax)zmax=coor(3,ip2)

          if(coor(1,ip3)<xmin)xmin=coor(1,ip3)
          if(coor(2,ip3)<ymin)ymin=coor(2,ip3)
          if(coor(3,ip3)<zmin)zmin=coor(3,ip3)
          if(coor(1,ip3)>xmax)xmax=coor(1,ip3)
          if(coor(2,ip3)>ymax)ymax=coor(2,ip3)
          if(coor(3,ip3)>zmax)zmax=coor(3,ip3)
          !
          !     Check against vmin,vmax
          !
          if(vmin(1)>xmax)cycle
          if(vmin(2)>ymax)cycle
          if(vmin(3)>zmax)cycle
          if(vmax(1)<xmin)cycle
          if(vmax(2)<ymin)cycle
          if(vmax(3)<zmin)cycle
          !
          !     Get new point position
          !
          call newpnt(nface,npoin,lface,iface,nnofa,ndim,coor,rnofa,pnewt,&
               pmid,rtx,rty,rtz,ichk,rnew,ierr)
          !
          !     Get distance  
          !  
          call distTri(lface,nface,nnofa,pnewt,coor,ndim,npoin,rnofa,rdistt,d1t,d2t,d3t,ipoint,iface,pintt)
          if(rdistt<rdist)then    
             rdist=rdistt        
             ipoin=ipoint         
             ihost=iface          !Just for debug as not usefull
             d1=d1t               !
             d2=d2t               !
             d3=d3t               !
             pnew(1)=pnewt(1)
             pnew(2)=pnewt(2)
             pnew(3)=pnewt(3)
             pint(1)=pintt(1)
             pint(2)=pintt(2)
             pint(3)=pintt(3)
          endif
       endif
    enddo
    goto 7
    !
    !     The host element has still not been reached, we have the closest element---> diagnostic
    !     Check that the point is out of the surface
    !
    ichk=0_ip
    !
    !     Initialize ipoin, as in the case for which the elements are around the same point,
    !     the last element may not be the limit of the surface mesh for deciding if the point
    !     is outside the mesh
    !    
    ipoin=lstack3(nstack)

    do istack=nstack,1,-1

       ilast=lstack(istack)

       if(lstack3(istack)==ipoin)then

          if(d1<rtol)then
             ineigh=eltoel(1,ilast)
             if(ineigh==0)then
                ichk=1_ip
             else
                if(lsurf(ineigh)/=isurf)then
                   ichk=1_ip
                endif
             endif
          endif
          if(d2<rtol)then
             ineigh=eltoel(2,ilast)
             if(ineigh==0)then
                ichk=1_ip
             else
                if(lsurf(ineigh)/=isurf)then
                   ichk=1_ip
                endif
             endif
          endif
          if(d3<rtol)then
             ineigh=eltoel(3,ilast)
             if(ineigh==0)then
                ichk=1_ip
             else
                if(lsurf(ineigh)/=isurf)then
                   ichk=1_ip
                endif
             endif
          endif

          if(ichk==1)then

             ihost=0_ip
             !
             !     Clean up lelem
             !
             do jstack=1,nstack2
                lelem(lstack2(jstack))=0_ip
             enddo
             !
             !     Go home
             !
             return

          endif

       endif

    enddo
    !
    !     We did not find an element for which the way to go to pnew is outside the surface
    !
    if(ichk==0)then
       write(*,*)'Strange configuration in gthost2'
       ierr=1_ip
       !
       !     Clean up lelem
       !
       do jstack=1,nstack2
          lelem(lstack2(jstack))=0_ip
       enddo
       return
    endif

7   continue
    !
    !     Clean up lelem
    !
    do istack=1,nstack2
       lelem(lstack2(istack))=0_ip
    enddo
    !
    !     Set nstack to 0 for lelem clean up 
    !
    nstack=0_ip

10  continue
    !
    !     Last part: we were successfull
    !
    !                - on input:  a good ihost
    !                - on output: a better ihost
    !
    !     Clean up
    !
    do istack=1,nstack
       lelem(lstack(istack))=0_ip
    enddo
    if(ihost==0)then
       write(*,*)'Error, ihost=0'
       ierr=1
       return
    endif
    !
    !     Try to find better
    !
    istack=0_ip
    nstack=1_ip
    lstack(1)=ihost
    nstack2=1_ip
    lstack2(1)=ihost 
    lelem(ihost)=1_ip
    !
    !     Initialize the distance to the face ihost
    ! 
    call distTri(lface,nface,nnofa,pnew,coor,ndim,npoin,rnofa,rdist,d1,d2,d3,ipoin,ihost,pint)

    do         
       if(istack==nstack)exit
       istack=istack+1

       ielem=lstack(istack)
       do j=1,nnofa
          ineigh=eltoel(j,ielem)
          if(ineigh==0)cycle
          if(lsurf(ineigh)/=isurf)cycle
          if(lelem(ineigh)/=0)cycle
          lelem(ineigh)=1_ip
          if(nstack2==mstack)then
             write(*,*)'Overflow in stack in gthost2 for iguess:',iguess
             ierr=1_ip
             return
          endif
          nstack2=nstack2+1_ip
          lstack2(nstack2)=ineigh
          !
          !     Compute the new point position
          !
          call newpnt(nface,npoin,lface,ineigh,nnofa,ndim,coor,rnofa,pnewt,pmid,&
               rtx,rty,rtz,ichk,rnew,ierr)
          if(ierr==1)return
          if(ichk==-1)cycle
          !
          !     Get the distance to the face
          ! 
          call distTri(lface,nface,nnofa,pnewt,coor,ndim,npoin,rnofa,rdistt,d1t,d2t,d3t,ipoint,ineigh,pintt)

          if(ipoint/=0 .and. ipoin==ipoint)then

             if(nstack==mstack)then
                write(*,*)'Overflow in stack in gthost2 for iguess:',iguess
                ierr=1_ip
                return
             endif
             !
             !     Add to the stack
             ! 
             nstack=nstack+1
             lstack(nstack)=ineigh

          else if(rdistt<rdist)then
             !
             !     Check the shape function of the projected point
             !
             if(d1t>rtol .and. d2t>rtol .and. d3t>rtol)then
                ihost=ineigh
                d1=d1t
                d2=d2t
                d3=d3t
                rdist=rdistt
                pint(1)=pintt(1)
                pint(2)=pintt(2)
                pint(3)=pintt(3)
                pnew(1)=pnewt(1)
                pnew(2)=pnewt(2)
                pnew(3)=pnewt(3)

             endif

             if(nstack==mstack)then
                write(*,*)'Overflow in stack in gthost2 for iguess:',iguess
                ierr=1_ip
                return
             endif
             !
             !     Add to the stack
             ! 
             nstack=nstack+1
             lstack(nstack)=ineigh
             rdist=rdistt        
             ipoin=ipoint         
             ihost=ineigh         !Just for debug as not usefull
             d1=d1t               !
             d2=d2t               !
             d3=d3t               !
             pint(1)=pintt(1)
             pint(2)=pintt(2)
             pint(3)=pintt(3)
             pnew(1)=pnewt(1)
             pnew(2)=pnewt(2)
             pnew(3)=pnewt(3)

          endif
       enddo
    enddo
    !
    !     Clean up lelem
    !
    do istack=1,nstack2
       lelem(lstack2(istack))=0_ip
    enddo
    !
    !     Recompute shape functions for points outside
    !
    if(d1<rtol .or. d2<rtol .or. d3<rtol)then

       ip1=lface(1,ihost)
       ip2=lface(2,ihost)
       ip3=lface(3,ihost)
       !
       !     Compute area to normalize
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)
       call orient3D(p1,p2,rface,dtot,ndim)
       !
       !     Get shape functions
       !
       p1(1)=coor(1,ip2)-pint(1)
       p1(2)=coor(2,ip2)-pint(2)
       p1(3)=coor(3,ip2)-pint(3)
       p2(1)=coor(1,ip3)-pint(1)
       p2(2)=coor(2,ip3)-pint(2)
       p2(3)=coor(3,ip3)-pint(3)

       call orient3D(p1,p2,rface,d1,ndim)
       d1=d1/dtot 

       p1(1)=coor(1,ip3)-pint(1)
       p1(2)=coor(2,ip3)-pint(2)
       p1(3)=coor(3,ip3)-pint(3)
       p2(1)=coor(1,ip1)-pint(1)
       p2(2)=coor(2,ip1)-pint(2)
       p2(3)=coor(3,ip1)-pint(3)

       call orient3D(p1,p2,rface,d2,ndim)
       d2=d2/dtot 

       p1(1)=coor(1,ip1)-pint(1)
       p1(2)=coor(2,ip1)-pint(2)
       p1(3)=coor(3,ip1)-pint(3)
       p2(1)=coor(1,ip2)-pint(1)
       p2(2)=coor(2,ip2)-pint(2)
       p2(3)=coor(3,ip2)-pint(3)

       call orient3D(p1,p2,rface,d3,ndim)
       d3=d3/dtot 

    endif

  end subroutine gthost2

  subroutine gtbscal(lface,iface,pnew,coor,rscal,ndim,npoin,rnofa,nface,nnofa)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nface,ndim,npoin,iface,nnofa
    integer(ip),intent(in)          :: lface(nnofa,nface)
    real(rp),intent(in)             :: coor(ndim,npoin),pnew(ndim),rnofa(ndim,nface)
    real(rp),intent(inout)          :: rscal
    integer(ip)                     :: ip1
    real(rp)                        :: rx,ry,rz


    ip1=lface(1,iface)

    rx=pnew(1)-coor(1,ip1)
    ry=pnew(2)-coor(2,ip1)
    rz=pnew(3)-coor(3,ip1)

    rscal=rnofa(1,iface)*rx+rnofa(2,iface)*ry+rnofa(3,iface)*rz
    rscal=abs(rscal)

  end subroutine gtbscal

  subroutine newpnt(nface,npoin,lface,iface,nnofa,ndim,coor,rnofa,pnew,&
       pmid,rtx,rty,rtz,ichk,rnew,ierr)
    use def_kintyp, only       :  ip,rp,lg
    implicit none 
    integer(ip),intent(in)    :: npoin,nface,nnofa,ndim,iface
    integer(ip),intent(in)    :: lface(nnofa,nface)
    integer(ip),intent(inout) :: ichk,ierr 
    real(rp),intent(in)       :: coor(ndim,npoin),rnofa(ndim,nface),pmid(ndim),rtx,rty,rtz,rnew
    real(rp),intent(inout)    :: pnew(ndim)
    integer(ip)               :: ip1,ip2,ip3 
    real(rp)                  :: ux,uy,uz,rkl,ul,d1a,d2a,cscal,det,c1,c2,x1,y1,z1,a,b,c,rface(ndim),c10,t1,c00

    c10=1.0d+00
    c00=0.0d+00
    !
    !     Initialize state flag 
    !
    ichk=0_ip 
    !
    !     Compute the vector carrying the line intersection of the plane normal to the current side 
    !     and normal to the current face normal
    !
    ip1=lface(1,iface)
    ip2=lface(2,iface)
    ip3=lface(3,iface)

    rface(1)=rnofa(1,iface)
    rface(2)=rnofa(2,iface)
    rface(3)=rnofa(3,iface)

    ux= rty*rface(3)-rtz*rface(2)     
    uy=-rtx*rface(3)+rtz*rface(1)
    uz= rtx*rface(2)-rty*rface(1)
    rkl=sqrt(ux*ux+uy*uy+uz*uz)
    ul=c10/rkl
    ux=ux*ul
    uy=uy*ul
    uz=uz*ul
    !
    !     Define first plane through pmid with rt
    !
    d1a=pmid(1)*rtx+pmid(2)*rty+pmid(3)*rtz
    !
    !     Define second plane through ip1 with rface
    !
    d2a=coor(1,ip1)*rface(1)+coor(2,ip1)*rface(2)+coor(3,ip1)*rface(3)
    !
    !     The equation of the line is:
    !     
    !      p = c1 N1 + c2 N2 + u N1 * N2
    !
    !     Solve for c1 and c2
    !
    cscal=rface(1)*rtx+rface(2)*rty+rface(3)*rtz
    det=c10-cscal*cscal 
    c1=(d1a-d2a*cscal)/det
    c2=(d2a-d1a*cscal)/det
    !
    !     Compute one point on the line  
    !
    x1=c1*rtx+c2*rface(1) 
    y1=c1*rty+c2*rface(2) 
    z1=c1*rtz+c2*rface(3) 
    !
    !     Compute the intersection of the sphere of radius rnew 
    !
    a=c10 
    b=2.0d+00*((ux*(x1-pmid(1)))+(uy*(y1-pmid(2)))+(uz*(z1-pmid(3))))
    c=pmid(1)*pmid(1)+pmid(2)*pmid(2)+pmid(3)*pmid(3)+x1*x1+y1*y1+z1*z1&
         -2.0d+00*(pmid(1)*x1+pmid(2)*y1+pmid(3)*z1)-rnew*rnew
    !
    !     Check for solutions
    ! 
    det=b*b-4.0d+00*a*c
    if(det<c00)then
       !
       !     We have gone too far
       !
       ichk=-1_ip
       return    
    else if(det>c00)then
       det=sqrt(det)
       t1=(-b+det)/(2.0d+00*a)
       !t1=(-b-det)/(2.0d+00*a)
       pnew(1)=x1+t1*ux
       pnew(2)=y1+t1*uy
       pnew(3)=z1+t1*uz

       !write(*,*)'radius',sqrt((pnew(1)-pmid(1))**2+(pnew(2)-pmid(2))**2+(pnew(3)-pmid(3))**2)

    else
       write(*,*)'error det 2 c'
       !ihost=-1
       !return
       ierr=1 
       return
    endif

  end subroutine newpnt

  subroutine distTri(lface,nface,nnofa,pnew,coor,ndim,npoin,rnofa,rdist,d1,d2,d3,ipoin,iface,pint)
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    use def_kintyp, only       :  ip,rp,lg
    implicit none 
    integer(ip),intent(in)    :: npoin,nface,nnofa,ndim,iface
    integer(ip),intent(in)    :: lface(nnofa,nface)
    real(rp),intent(in)       :: coor(ndim,npoin),rnofa(ndim,nface),pnew(ndim)
    real(rp),intent(inout)    :: rdist,d1,d2,d3,pint(3)
    integer(ip),intent(inout) :: ipoin
    integer(ip)               :: ip1,ip2,ip3,ipoin1,ipoin2,ipoin3
    real(rp)                  :: rdirx,rdiry,rdirz,csca,cscal,xpro,ypro,zpro,epsil,rl,epsil0,epsil1
    real(rp)                  :: rdist1,rdist2,rdist3,dtot,p1(ndim),p2(ndim),c10,rface(ndim)
    real(rp)                  :: rtx1,rty1,rtz1,rtx2,rty2,rtz2,rtl,pint1(3),pint2(3),pint3(3)

    c10=1.0d+00 
    epsil=1.0d-06
    epsil1=1.0d+00+1.0d-04
    epsil0=1.0d-04
    !
    !     Set ipoin to zero
    ! 
    ipoin=0_ip
    ipoin1=0_ip
    ipoin2=0_ip
    ipoin3=0_ip
    !
    !     Get the points of the face
    !
    ip1=lface(1,iface)
    ip2=lface(2,iface)
    ip3=lface(3,iface)
    !
    !     Host face normal
    !
    rface(1)=rnofa(1,iface)
    rface(2)=rnofa(2,iface)
    rface(3)=rnofa(3,iface)
    !
    !     Project the point on the plane of the source
    !
    rdirx=pnew(1)-coor(1,ip1)
    rdiry=pnew(2)-coor(2,ip1)
    rdirz=pnew(3)-coor(3,ip1)

    csca=rface(1)*rdirx+rface(2)*rdiry+rface(3)*rdirz

    xpro=pnew(1)-csca*rface(1)
    ypro=pnew(2)-csca*rface(2)
    zpro=pnew(3)-csca*rface(3)
    !
    !     Compute area to normalize
    !
    p1(1)=coor(1,ip2)-coor(1,ip1)
    p1(2)=coor(2,ip2)-coor(2,ip1)
    p1(3)=coor(3,ip2)-coor(3,ip1)
    p2(1)=coor(1,ip3)-coor(1,ip1)
    p2(2)=coor(2,ip3)-coor(2,ip1)
    p2(3)=coor(3,ip3)-coor(3,ip1)
    call orient3D(p1,p2,rface,dtot,ndim)
    !
    !     Side (ip2,ip3)
    !
    p1(1)=coor(1,ip2)-xpro
    p1(2)=coor(2,ip2)-ypro
    p1(3)=coor(3,ip2)-zpro
    p2(1)=coor(1,ip3)-xpro
    p2(2)=coor(2,ip3)-ypro
    p2(3)=coor(3,ip3)-zpro
    call orient3D(p1,p2,rface,d1,ndim)
    d1=d1/dtot
    !
    !     Side (ip3,ip1)
    !
    p1(1)=coor(1,ip3)-xpro
    p1(2)=coor(2,ip3)-ypro
    p1(3)=coor(3,ip3)-zpro
    p2(1)=coor(1,ip1)-xpro
    p2(2)=coor(2,ip1)-ypro
    p2(3)=coor(3,ip1)-zpro
    call orient3D(p1,p2,rface,d2,ndim)
    d2=d2/dtot
    !
    !     Side (ip1,ip2)
    !
    p1(1)=coor(1,ip1)-xpro
    p1(2)=coor(2,ip1)-ypro
    p1(3)=coor(3,ip1)-zpro
    p2(1)=coor(1,ip2)-xpro
    p2(2)=coor(2,ip2)-ypro
    p2(3)=coor(3,ip2)-zpro
    call orient3D(p1,p2,rface,d3,ndim)
    d3=d3/dtot
    !
    !     Is the projected point inside the source triangle?
    !
    if(d1>-epsil .and. d2>-epsil .and. d3>-epsil)then  
       rdist=abs(csca)
       pint(1)=xpro
       pint(2)=ypro
       pint(3)=zpro
    else
       !
       !     Evaluate distance to lines
       !

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

          rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
          ipoin1=ip2
          pint1(1)=coor(1,ip2)
          pint1(2)=coor(2,ip2)
          pint1(3)=coor(3,ip2)

       else if(cscal>epsil1)then

          p2(1)=pnew(1)-coor(1,ip3)
          p2(2)=pnew(2)-coor(2,ip3)
          p2(3)=pnew(3)-coor(3,ip3)
          rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
          ipoin1=ip3
          pint1(1)=coor(1,ip3)
          pint1(2)=coor(2,ip3)
          pint1(3)=coor(3,ip3)

       else

          p2(1)=coor(1,ip2)+csca*p1(1)-pnew(1)
          p2(2)=coor(2,ip2)+csca*p1(2)-pnew(2)
          p2(3)=coor(3,ip2)+csca*p1(3)-pnew(3)
          rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 
          pint1(1)=coor(1,ip2)+csca*p1(1)
          pint1(2)=coor(2,ip2)+csca*p1(2)
          pint1(3)=coor(3,ip2)+csca*p1(3)

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

          rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
          ipoin2=ip3
          pint2(1)=coor(1,ip3)
          pint2(2)=coor(2,ip3)
          pint2(3)=coor(3,ip3)

       else if(cscal>epsil1)then

          p2(1)=pnew(1)-coor(1,ip1)
          p2(2)=pnew(2)-coor(2,ip1)
          p2(3)=pnew(3)-coor(3,ip1)
          rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
          ipoin2=ip1
          pint2(1)=coor(1,ip1)
          pint2(2)=coor(2,ip1)
          pint2(3)=coor(3,ip1)

       else

          p2(1)=coor(1,ip3)+csca*p1(1)-pnew(1)
          p2(2)=coor(2,ip3)+csca*p1(2)-pnew(2)
          p2(3)=coor(3,ip3)+csca*p1(3)-pnew(3)
          rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 
          pint2(1)=coor(1,ip3)+csca*p1(1)
          pint2(2)=coor(2,ip3)+csca*p1(2)
          pint2(3)=coor(3,ip3)+csca*p1(3)

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

          rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
          ipoin3=ip1
          pint3(1)=coor(1,ip1)
          pint3(2)=coor(2,ip1)
          pint3(3)=coor(3,ip1)

       else if(cscal>epsil1)then

          p2(1)=pnew(1)-coor(1,ip2)
          p2(2)=pnew(2)-coor(2,ip2)
          p2(3)=pnew(3)-coor(3,ip2)
          rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
          ipoin3=ip2
          pint3(1)=coor(1,ip2)
          pint3(2)=coor(2,ip2)
          pint3(3)=coor(3,ip2)

       else

          p2(1)=coor(1,ip1)+csca*p1(1)-pnew(1)
          p2(2)=coor(2,ip1)+csca*p1(2)-pnew(2)
          p2(3)=coor(3,ip1)+csca*p1(3)-pnew(3)
          rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 
          pint3(1)=coor(1,ip1)+csca*p1(1)
          pint3(2)=coor(2,ip1)+csca*p1(2)
          pint3(3)=coor(3,ip1)+csca*p1(3)

       endif
       !
       !     Take the min over the side distances
       !
       if(rdist1<rdist2)then 
          rdist=rdist1
          ipoin=ipoin1
          pint(1)=pint1(1)
          pint(2)=pint1(2)
          pint(3)=pint1(3)
       else 
          rdist=rdist2
          ipoin=ipoin2
          pint(1)=pint2(1)
          pint(2)=pint2(2)
          pint(3)=pint2(3)
       endif

       if(rdist3<rdist)then 
          rdist=rdist3
          ipoin=ipoin3
          pint(1)=pint3(1)
          pint(2)=pint3(2)
          pint(3)=pint3(3)
       endif

    endif

  end subroutine distTri

  subroutine chkgeo(lface,nface,nnofa,iface,minscal,rface,rnopo,npoin,ndim)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)     :: nface,nnofa,iface,npoin,ndim
    integer(ip),intent(in)     :: lface(nnofa,nface)
    real(rp),intent(in)        :: rface(ndim),rnopo(ndim,npoin)
    real(rp),intent(inout)     :: minscal
    integer(ip)                :: ip1,ip2,ip3
    real(rp)                   :: rscal              
    !
    !     Compare the normal of the face to the normal of its points
    !


    ip1=lface(1,iface)
    ip2=lface(2,iface)
    ip3=lface(3,iface)

    minscal=rnopo(1,ip1)*rface(1)+rnopo(2,ip1)*rface(2)+rnopo(3,ip1)*rface(3)
    rscal=rnopo(1,ip2)*rface(1)+rnopo(2,ip2)*rface(2)+rnopo(3,ip2)*rface(3)
    if(rscal<minscal)minscal=rscal
    rscal=rnopo(1,ip3)*rface(1)+rnopo(2,ip3)*rface(2)+rnopo(3,ip3)*rface(3)
    if(rscal<minscal)minscal=rscal

  end subroutine chkgeo

  subroutine quality(lface,nface,iface,nnofa,coor,ndim,npoin,q,vt)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nface,ndim,npoin,nnofa,iface
    integer(ip),intent(in)          :: lface(nnofa,nface) 
    real(rp),intent(in)             :: coor(ndim,npoin),vt
    real(rp),intent(out)            :: q
    integer(ip)                     :: p1,p2,p3
    real(rp)                        :: rx,ry,rz,rl1,rl2,rl3,c05,hmax,perim,alpha 
    real(rp)                        :: dx12,dy12,dz12,dx13,dy13,dz13,rnl,rnx,rny,rnz
    !
    !     This sub computes  the quality of a tetrahedron with the Gamma quality
    !
    c05=0.5d+00
    alpha=sqrt(3.0d+00)/6.0d+00
    p1=lface(1,iface)
    p2=lface(2,iface)
    p3=lface(3,iface)

    rx=coor(1,p2)-coor(1,p1)
    ry=coor(2,p2)-coor(2,p1)
    rz=coor(3,p2)-coor(3,p1)
    rl1=sqrt(rx*rx+ry*ry+rz*rz)

    rx=coor(1,p3)-coor(1,p2)
    ry=coor(2,p3)-coor(2,p2)
    rz=coor(3,p3)-coor(3,p2)
    rl2=sqrt(rx*rx+ry*ry+rz*rz)

    rx=coor(1,p1)-coor(1,p3)
    ry=coor(2,p1)-coor(2,p3)
    rz=coor(3,p1)-coor(3,p3)
    rl3=sqrt(rx*rx+ry*ry+rz*rz)

    hmax=max(rl1,rl2,rl3)
    perim=(rl1+rl2+rl3)*c05
    q=alpha*hmax*perim/vt

    !if(q<1.0d+00)then

    !   dx12=coor(1,p2)-coor(1,p1)
    !   dy12=coor(2,p2)-coor(2,p1)
    !   dz12=coor(3,p2)-coor(3,p1)

    !   dx13=coor(1,p3)-coor(1,p1)
    !   dy13=coor(2,p3)-coor(2,p1)
    !   dz13=coor(3,p3)-coor(3,p1)

    !   rnx= dy12*dz13-dz12*dy13
    !   rny=-dx12*dz13+dz12*dx13
    !   rnz= dx12*dy13-dy12*dx13

    !   rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
    !   rnl=rnl/2.0d+00

    !   write(*,*)rnl,vt
    !   write(*,*)'error quality'
    !   stop
    !endif

  end subroutine quality

  subroutine insurf(nface,npoin,nnofa,ndim,nnosi,nside,nsurf,nline,&
       lface,coor,lsurf,lside,lline)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(inout)   :: nface,npoin,nsurf,nside,nline
    integer(ip), intent(in)      :: nnofa,ndim,nnosi
    integer(ip),pointer          :: lface(:,:),lsurf(:),lside(:,:),lline(:)
    real(rp),pointer             :: coor(:,:)
    integer(ip)                  :: iface,ipoin,jface,i,inofa,icont
    character*80              :: ntext
    integer(ip)                  :: iostatus
    integer(4)                   :: istat

    open(unit=77,file='dekhna.geo',status='old',iostat=iostatus)
    if(iostatus>0)then
       write(*,*)'Surface mesh file not found in insurf'
       stop
    endif

    write(*,*)'Reading surface mesh file'
    !
    !     Read points
    !
    read(77,*)ntext,npoin
    write(*,*)'npoin=',npoin
    allocate(coor(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'COOR','insurf',coor)

    do i=1,npoin
       read(77,*) ipoin,coor(1,i),coor(2,i),coor(3,i)
    end do
    !
    !     Read ridges
    !
    read(77,*)ntext,nline
    write(*,*)'nline=',nline

    read(77,*)ntext,nside
    write(*,*)'nside=',nside

    allocate(lside(nnosi,nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LSIDE','insurf',lside)
    allocate(lline(nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LLINE','insurf',lline) 
    do i=1,nside
       read(77,*) ipoin,lside(1,i),lside(2,i),lline(i)
    end do
    !
    !     Read faces
    !
    read(77,*)ntext,nsurf
    write(*,*)'nsurf=',nsurf

    read(77,*)ntext,nface
    write(*,*)'nface=',nface
    allocate(lface(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFACE','insurf',lface) 
    allocate(lsurf(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSURF','insurf',lsurf) 
    do iface=1,nface
       !read(77,*) jface,lface(1,iface),lface(2,iface),lface(3,iface),lsurf(iface)
       read(77,*) jface,lface(1,iface),lface(2,iface),lface(3,iface),lsurf(iface)
       !lsurf(iface)=1_ip  
    enddo

    close(77)
    write(*,*)'End reading surface mesh file'

    !
    !     Compact the points 
    !
    !allocate(lmark(npoin),stat=istat)
    !call memchk(zero,istat,memor_msh,'LMARK','insurf',lmark) 
    !do iface=1,nface
    !   do inofa=1,nnofa
    !      lmark(lface(inofa,iface))=1_ip
    !   enddo
    !enddo   

    !icont=0_ip
    !do ipoin=1,npoin
    !   if(lmark(ipoin)==1)then
    !      icont=icont+1_ip
    !      lmark(ipoin)=icont
    !      coor(1,icont)=coor(1,ipoin)
    !      coor(2,icont)=coor(2,ipoin)
    !      coor(3,icont)=coor(3,ipoin)
    !   endif
    !enddo    
    !npoin=icont

    !do iface=1,nface
    !  lface(1,iface)=lmark(lface(1,iface))
    !  lface(2,iface)=lmark(lface(2,iface))
    !  lface(3,iface)=lmark(lface(3,iface))
    !enddo

    !    call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)
    !   stop

  end subroutine insurf

  subroutine swapf(ielem1,ielem2,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,&
       iopt,iqual,lptri,rnofa)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nface,ndim,npoin,nnofa,ielem1,ielem2,iqual
    integer(ip),intent(inout) :: lface(nnofa,nface),eltoel(nnofa,nface),iopt,lptri(npoin) 
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin)
    real(rp),intent(inout)    :: rnofa(nnofa,nface)
    integer(ip)               :: lfacn(nnofa,2),iview1,iview2,p1,p2,p3,p4   
    integer(ip)               :: ineigh11,ineigh12,ineigh21,ineigh22
    real(rp)                  :: Qgeo1,Qgeo2,Qgeotol,rfloc(ndim,2),modul(2),c05,q1,q2,qold
    real(rp)                  :: modulold(2),rflocold(ndim,2),q1old,q2old,c00
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))

    c05=0.5d+00
    c00=0.0d+00
    Qgeotol=0.0d+00
    Qgeotol=0.6d+00
    !
    !     This subroutine swaps two elements ielem1 and ielem2 if it is possible
    !     If iqual==1 swaping is done if quality is improved 
    !     Compared to swap, this subroutine updates rnofa & lptri
    !     which are both needed during refsmo
    !
    !     Find viewer
    !
    if(eltoel(1,ielem1)==ielem2)then
       iview1=1
    else if(eltoel(2,ielem1)==ielem2)then
       iview1=2
    else
       iview1=3
    endif

    if(eltoel(1,ielem2)==ielem1)then
       iview2=1
    else if(eltoel(2,ielem2)==ielem1)then
       iview2=2
    else
       iview2=3
    endif

    p1=lface(iview1,ielem1)
    p2=lface(ltab(1,iview1),ielem1)
    p3=lface(ltab(2,iview1),ielem1)
    p4=lface(iview2,ielem2)

    lfacn(1,1)=p1 
    lfacn(2,1)=p2 
    lfacn(3,1)=p4 

    lfacn(1,2)=p1 
    lfacn(2,2)=p4 
    lfacn(3,2)=p3 


    call gtfnr2(lfacn,2_ip,nnofa,ndim,coor,npoin,rfloc(1,1),1_ip,modul(1))
    if(modul(1)==c00)return
    call chkgeo(lfacn,2_ip,nnofa,1_ip,Qgeo1,rfloc(1,1),rnopo,npoin,ndim)
    if(Qgeo1<Qgeotol)return

    call gtfnr2(lfacn,2_ip,nnofa,ndim,coor,npoin,rfloc(1,2),2_ip,modul(2))
    if(modul(2)==c00)return
    call chkgeo(lfacn,2_ip,nnofa,2_ip,Qgeo2,rfloc(1,2),rnopo,npoin,ndim)
    if(Qgeo2<Qgeotol)return
    !
    !     Do we care about increasing the quality?
    !
    if(iqual==1)then
       !
       !     Compute surface of old elements
       !
       call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rflocold(1,1),ielem1,modulold(1))
       call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rflocold(1,2),ielem2,modulold(2))
       !
       !     Divide the area by 2
       !
       modulold(1)=modulold(1)*c05
       modulold(2)=modulold(2)*c05
       !
       !     Compute quality of old elements
       !
       call quality(lface,nface,ielem1,nnofa,coor,ndim,npoin,q1old,modulold(1))
       call quality(lface,nface,ielem2,nnofa,coor,ndim,npoin,q2old,modulold(2))
       qold=max(q1old,q2old)
       !
       !     Divide the area by 2 of the new elements
       !
       modul(1)=modul(1)*c05
       modul(2)=modul(2)*c05
       !
       !  Compute shape quality of the new elements
       !
       call quality(lfacn,2_ip,1_ip,nnofa,coor,ndim,npoin,q1,modul(1))
       if(q1>=qold)return
       call quality(lfacn,2_ip,2_ip,nnofa,coor,ndim,npoin,q2,modul(2))
       if(q2>=qold)return

    endif
    !
    !     Swap successfull
    !
    iopt=1_ip
    ineigh11=eltoel(ltab(1,iview1),ielem1)
    ineigh12=eltoel(ltab(2,iview1),ielem1)
    ineigh21=eltoel(ltab(1,iview2),ielem2)
    ineigh22=eltoel(ltab(2,iview2),ielem2)
    !
    !    Update lface && rnofa
    !
    lface(1,ielem1)=lfacn(1,1)
    lface(2,ielem1)=lfacn(2,1)
    lface(3,ielem1)=lfacn(3,1)
    rnofa(1,ielem1)=rfloc(1,1)  
    rnofa(2,ielem1)=rfloc(2,1)  
    rnofa(3,ielem1)=rfloc(3,1)  

    lface(1,ielem2)=lfacn(1,2)
    lface(2,ielem2)=lfacn(2,2)
    lface(3,ielem2)=lfacn(3,2)
    rnofa(1,ielem2)=rfloc(1,2)  
    rnofa(2,ielem2)=rfloc(2,2)  
    rnofa(3,ielem2)=rfloc(3,2)  
    !
    !     Update eltoel
    !
    eltoel(1,ielem1)=ineigh21
    eltoel(2,ielem1)=ielem2
    eltoel(3,ielem1)=ineigh12

    eltoel(1,ielem2)=ineigh22
    eltoel(2,ielem2)=ineigh11
    eltoel(3,ielem2)=ielem1
    !
    !     Update outside
    !
    if(ineigh21/=0)then
       if(eltoel(1,ineigh21)==ielem2)then
          eltoel(1,ineigh21)=ielem1
       else if(eltoel(2,ineigh21)==ielem2)then
          eltoel(2,ineigh21)=ielem1
       else
          eltoel(3,ineigh21)=ielem1
       endif
    endif

    if(ineigh11/=0)then
       if(eltoel(1,ineigh11)==ielem1)then
          eltoel(1,ineigh11)=ielem2
       else if(eltoel(2,ineigh11)==ielem1)then
          eltoel(2,ineigh11)=ielem2
       else
          eltoel(3,ineigh11)=ielem2
       endif
    endif

    lptri(p2)=ielem1
    lptri(p3)=ielem2

    !call chkrnofa(lface,nface,nnofa,ndim,coor,npoin,rnofa)

  end subroutine swapf

  subroutine swaprecover(ielem1,ielem2,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,&
       iopt,lptri,rnofa)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nface,ndim,npoin,nnofa,ielem1,ielem2
    integer(ip),intent(inout) :: lface(nnofa,nface),eltoel(nnofa,nface),iopt,lptri(npoin) 
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin)
    real(rp),intent(inout)    :: rnofa(nnofa,nface)
    integer(ip)               :: lfacn(nnofa,2),iview1,iview2,p1,p2,p3,p4   
    integer(ip)               :: ineigh11,ineigh12,ineigh21,ineigh22
    real(rp)                  :: Qgeo1,Qgeo2,Qgeotol,rfloc(ndim,2),modul(2),c05,q1,q2,qold
    real(rp)                  :: modulold(2),rflocold(ndim,2),q1old,q2old,c00
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))

    c05=0.5d+00
    c00=0.0d+00
    Qgeotol=0.0d+00
    Qgeotol=0.0001d+00
    !
    !     This subroutine swaps two elements ielem1 and ielem2 if it is possible
    !     Swaping is not done regarding quality  
    !     Compared to swap, this subroutine updates rnofa & lptri
    !     which are both needed during refsmo
    !     Compared to swapf, Qgeo has been lowered drastically
    !
    !     Find viewer
    !
    if(eltoel(1,ielem1)==ielem2)then
       iview1=1
    else if(eltoel(2,ielem1)==ielem2)then
       iview1=2
    else
       iview1=3
    endif

    if(eltoel(1,ielem2)==ielem1)then
       iview2=1
    else if(eltoel(2,ielem2)==ielem1)then
       iview2=2
    else
       iview2=3
    endif

    p1=lface(iview1,ielem1)
    p2=lface(ltab(1,iview1),ielem1)
    p3=lface(ltab(2,iview1),ielem1)
    p4=lface(iview2,ielem2)

    lfacn(1,1)=p1 
    lfacn(2,1)=p2 
    lfacn(3,1)=p4 

    lfacn(1,2)=p1 
    lfacn(2,2)=p4 
    lfacn(3,2)=p3 


    call gtfnr2(lfacn,2_ip,nnofa,ndim,coor,npoin,rfloc(1,1),1_ip,modul(1))
    if(modul(1)==c00)return
    call chkgeo(lfacn,2_ip,nnofa,1_ip,Qgeo1,rfloc(1,1),rnopo,npoin,ndim)
    if(Qgeo1<Qgeotol)return

    call gtfnr2(lfacn,2_ip,nnofa,ndim,coor,npoin,rfloc(1,2),2_ip,modul(2))
    if(modul(2)==c00)return
    call chkgeo(lfacn,2_ip,nnofa,2_ip,Qgeo2,rfloc(1,2),rnopo,npoin,ndim)
    if(Qgeo2<Qgeotol)return
    !
    !     Swap successfull
    !
    iopt=1_ip
    ineigh11=eltoel(ltab(1,iview1),ielem1)
    ineigh12=eltoel(ltab(2,iview1),ielem1)
    ineigh21=eltoel(ltab(1,iview2),ielem2)
    ineigh22=eltoel(ltab(2,iview2),ielem2)
    !
    !    Update lface && rnofa
    !
    lface(1,ielem1)=lfacn(1,1)
    lface(2,ielem1)=lfacn(2,1)
    lface(3,ielem1)=lfacn(3,1)
    rnofa(1,ielem1)=rfloc(1,1)  
    rnofa(2,ielem1)=rfloc(2,1)  
    rnofa(3,ielem1)=rfloc(3,1)  

    lface(1,ielem2)=lfacn(1,2)
    lface(2,ielem2)=lfacn(2,2)
    lface(3,ielem2)=lfacn(3,2)
    rnofa(1,ielem2)=rfloc(1,2)  
    rnofa(2,ielem2)=rfloc(2,2)  
    rnofa(3,ielem2)=rfloc(3,2)  
    !
    !     Update eltoel
    !
    eltoel(1,ielem1)=ineigh21
    eltoel(2,ielem1)=ielem2
    eltoel(3,ielem1)=ineigh12

    eltoel(1,ielem2)=ineigh22
    eltoel(2,ielem2)=ineigh11
    eltoel(3,ielem2)=ielem1
    !
    !     Update outside
    !
    if(ineigh21/=0)then
       if(eltoel(1,ineigh21)==ielem2)then
          eltoel(1,ineigh21)=ielem1
       else if(eltoel(2,ineigh21)==ielem2)then
          eltoel(2,ineigh21)=ielem1
       else
          eltoel(3,ineigh21)=ielem1
       endif
    endif

    if(ineigh11/=0)then
       if(eltoel(1,ineigh11)==ielem1)then
          eltoel(1,ineigh11)=ielem2
       else if(eltoel(2,ineigh11)==ielem1)then
          eltoel(2,ineigh11)=ielem2
       else
          eltoel(3,ineigh11)=ielem2
       endif
    endif

    lptri(p2)=ielem1
    lptri(p3)=ielem2

    !call chkrnofa(lface,nface,nnofa,ndim,coor,npoin,rnofa)

  end subroutine swaprecover

  subroutine colsmo(lface,nnofa,nface,lptype,lmark,npoin,lelem,coor,ndim,ptoel1,ptoel2,&
       rnopo,tollen,rsize,eltoel,rmax,tolref,ledge,redge,ledg2,lstack,tolgeo)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    use mod_sort, only        :  sort3
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,nface,npoin
    integer(ip),intent(in)    :: lptype(2,npoin)
    integer(ip),intent(inout) :: lface(nnofa,nface),ledge(2,3*nface),ledg2(3*nface)
    real(rp),intent(inout)    :: redge(3*nface)
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin),tollen,rsize(npoin),tolref,tolgeo    
    integer(ip),intent(inout) :: lelem(nface)
    integer(ip),intent(in)    :: eltoel(nnofa,nface)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:)
    integer(ip),intent(inout) :: lmark(npoin),lstack(npoin)
    real(rp),intent(inout)    :: rmax
    integer(ip)               :: iface,ip1,ip2,nedge,iedge,jedge,icol,icollapse,ichk,iedg
    integer(ip)               :: ipoin,isto,jpoin,inofa
    real(rp)                  :: rlen 
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub collapses the short edges for each patch
    !
    !
    !     Initialize rmax the largest edge
    !
    rmax=0.0d+00
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Two loops, if collapse is controled or aggresive 
    !
    do icollapse=1,2
       !
       !     Initialize edges found to be collapsed
       !
       nedge=0_ip
       !
       !     Clean up lstack
       ! 
       do ipoin=1,npoin
          lstack(ipoin)=0_ip
       enddo
       ! 
       !     Loop on points
       !
       do ipoin=1,npoin
          !
          !     Loop on faces surrounding points
          !
          do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1

             iface=ptoel1(isto)  
             !
             !     Has this face been already marked ?
             !
             if(lelem(iface)/=0)cycle

             do inofa=1,nnofa
                jpoin=lface(inofa,iface)
                if(jpoin>ipoin)then
                   if(lstack(jpoin)/=ipoin)then
                      lstack(jpoin)=ipoin   
                      !
                      !     Get the edges to collapse
                      !  
                      ichk=0_ip
                      if(icollapse==1)then
                         if(lmark(ipoin)==1 .and. lmark(jpoin)==-1)then
                            !
                            !     jpoin must be a smooth point
                            !
                            if(lptype(1,jpoin)==ID_SMOOTH)then  
                               call length(ipoin,jpoin,rsize,coor,rlen)
                               if(rlen>rmax)rmax=rlen
                               if(rlen<tollen)then
                                  nedge=nedge+1
                                  ledge(1,nedge)=ipoin
                                  ledge(2,nedge)=jpoin
                                  ledg2(nedge)=nedge
                                  redge(nedge)=rlen
                               endif
                            endif
                         else if(lmark(ipoin)==-1 .and. lmark(jpoin)==1)then
                            !
                            !     ipoin must be a smooth point
                            !
                            if(lptype(1,ipoin)==ID_SMOOTH)then  
                               call length(ipoin,jpoin,rsize,coor,rlen)
                               if(rlen>rmax)rmax=rlen
                               if(rlen<tollen)then
                                  nedge=nedge+1
                                  ledge(1,nedge)=jpoin
                                  ledge(2,nedge)=ipoin
                                  ledg2(nedge)=nedge
                                  redge(nedge)=rlen
                               endif
                            endif
                         endif
                      else 
                         !
                         !     jpoin must be a smooth point
                         !
                         if(lptype(1,jpoin)==ID_SMOOTH)then  
                            call length(ipoin,jpoin,rsize,coor,rlen)
                            if(rlen>rmax)rmax=rlen
                            if(rlen<tollen)then
                               nedge=nedge+1
                               ledge(1,nedge)=ipoin
                               ledge(2,nedge)=jpoin
                               ledg2(nedge)=nedge
                               redge(nedge)=rlen
                            endif
                            !
                            !     ipoin must be a smooth point
                            !
                         else if(lptype(1,ipoin)==ID_SMOOTH)then  
                            call length(ipoin,jpoin,rsize,coor,rlen)
                            if(rlen>rmax)rmax=rlen
                            if(rlen<tollen)then
                               nedge=nedge+1
                               ledge(1,nedge)=jpoin
                               ledge(2,nedge)=ipoin
                               ledg2(nedge)=nedge
                               redge(nedge)=rlen
                            endif
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
       !
       !     Sort the edges
       !
       call  sort3(ledg2,redge,nedge)
       !
       !     Collapse in order
       !
       do iedg=1,nedge
          iedge=ledg2(iedg) 
          ip1=ledge(1,iedge)     
          ip2=ledge(2,iedge)
          !icol=0
          call collapsedg2(ip1,ip2,lface,coor,ndim,nnofa,nface,npoin,lelem,ptoel1,ptoel2,&
               lmark,rnopo,icol,eltoel,rsize,tolref,tolgeo)
       enddo

    enddo

  end subroutine colsmo

  subroutine collapsedg2(ip1,ip2,lface,coor,ndim,nnofa,nface,npoin,lelem,ptoel1,&
       ptoel2,lmark,rnopo,icol,eltoel,rsize,tollen,tolgeo)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only     :  memor_msh
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,ip1,ip2
    integer(ip),intent(in)    :: ptoel1(*),ptoel2(*)
    integer(ip),intent(in)    :: nface,npoin
    integer(ip),intent(inout) :: icol
    integer(ip),intent(inout) :: lelem(nface)
    integer(ip),intent(inout) :: lface(nnofa,nface),lmark(npoin)
    integer(ip),intent(in)    :: eltoel(nnofa,nface)
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin),rsize(npoin),tollen,tolgeo
    integer(ip)               :: lfacn(nnofa,100),lballn(100),lold(2),nold
    integer(ip)               :: p1,p2,ielem,j,ienew,i,ie,iview1,iview2,ielem1,ielem2
    integer(ip)               :: ineigh,iplace,ielem0,nballn,icount,ip3
    integer(ip)               :: iter,niter,iview,ipa,ipb
    real(rp)                  :: rfloc(ndim,100),modul(100),Qgeo,c00,rlen
    real(rp)                  :: tolvol,rx,ry,rz
    integer(4)                :: istat
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This subroutine collapses an edge given by (ip1,ip2)
    !     This is the version used by colsmo to collapse smooth edges
    !     - eltoel is not updated as the collapsing is performed on the whole patch 
    !       each time, and eltoel is reconstructed each time
    !     - No optimization is performed. It is done on the global level
    !     Compared to collapsed: 
    !            - the sides are not updated (no need)
    !
    !
    !     Give geometric tolerance
    !
    !write(*,*)'ip1=',ip1,'ip2=',ip2
    c00=0.0d+00
    niter=10_ip
    !
    !     Try to collapse ip2 on ip1, replacing ip2 by ip1
    !
    !
    !     Get characteristic length
    ! 
    rx=coor(1,ip1)-coor(1,ip2)
    ry=coor(2,ip1)-coor(2,ip2)
    rz=coor(3,ip1)-coor(3,ip2)
    rlen=rx*rx+ry*ry+rz*rz
    tolvol=rlen*1.0d-02
    !
    !     Find the new elements = ball of ip2 -(ball of ip1 inter ball of pb )
    !
    nballn=0_ip
    nold=0_ip
    do ie=ptoel2(ip2),ptoel2(ip2+1)-1
       ielem=ptoel1(ie)
       !
       !     Go home if an element has been touched 
       !
       if(lelem(ielem)/=0)return
       !
       !     Count if ip1 and ip2 belong to ielem
       !
       icount=0_ip
       do j=1,nnofa
          ip3=lface(j,ielem)
          if(ip3==ip1)then
             icount=icount+1
          else if(ip3==ip2)then
             icount=icount+1
          endif
       enddo
       if(icount/=2)then
          nballn=nballn+1
          lballn(nballn)=ielem
          do j=1,nnofa 
             p1=lface(j,ielem)
             if(p1==ip2)then
                lfacn(j,nballn)=ip1
             else
                lfacn(j,nballn)=p1
             endif
          enddo
       else
          if(nold==2)then 
             write(*,*)'Error in collapsedg2, nold=2'
             stop
          endif
          nold=nold+1
          lold(nold)=ielem
       endif
    enddo
    !
    !     Compute normals and test the geometric quality of the new elements
    !
    do ie=1,nballn
       call gtfnr2(lfacn,nballn,nnofa,ndim,coor,npoin,rfloc(1,ie),ie,modul(ie))
       if(modul(ie)<tolvol)return
       call chkgeo(lfacn,nballn,nnofa,ie,Qgeo,rfloc(1,ie),rnopo,npoin,ndim)
       if(Qgeo<tolgeo)return
    enddo
    !
    !     Check the length of the new edges
    !
    !     THIS IS DUMB
    !     We need to collapse this edge anyway
    !
    ! 
    !do ie=1,nballn 
    !   if(lfacn(1,ie)==ip1)then 
    !      iview=1_ip
    !   else if(lfacn(2,ie)==ip1)then 
    !      iview=2_ip
    !   else 
    !      iview=3_ip
    !   endif

    !   ipa=lfacn(ltab(1,iview),ie)
    !   ipb=lfacn(ltab(2,iview),ie)
    !   call length(ip1,ipa,rsize,coor,rlen)
    !   if(rlen>tollen)return
    !   call length(ip1,ipb,rsize,coor,rlen)
    !   if(rlen>tollen)return
    !enddo
    !
    !     Finally collapse
    !
    do ie=1,nballn
       ielem=lballn(ie)
       lface(1,ielem)=lfacn(1,ie)
       lface(2,ielem)=lfacn(2,ie)
       lface(3,ielem)=lfacn(3,ie)
    enddo
    !
    !    Delete elements
    !
    ielem1=lold(1)
    ielem2=lold(2)
    lelem(ielem1)=-1_ip 
    lelem(ielem2)=-1_ip 
    !
    !     Mark the elements of ip2 as touched
    !  
    do ie=ptoel2(ip2),ptoel2(ip2+1)-1
       ielem=ptoel1(ie)
       if(lelem(ielem)==0)then
          lelem(ielem)=1_ip
       endif
    enddo
    !
    !     Delete ip2
    !
    lmark(ip2)=-2
    !
    !     Remember collapse successful
    !
    icol=1_ip
    !
    !   Debug
    !
    !write(*,*)'from collapsedg2 3'
    !call chkcons(lface,nnofa,nface,eltoel,lelem,coor,npoin,ndim,lmark)
    !call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)

  end subroutine collapsedg2

  subroutine mark(lstack,nstack,ptoel1,ptoel2,nnofa,lface,nface,npoin,lptype,&
       lmark,ncol,coor,ndim,rsize,nline,nside,nnosi,nsurf,lsurf,ptosi1,ptosi2,lline,lside)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only     :  memor_msh
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)        :: nnofa,ndim,nline,nside,nnosi,nsurf,npoin,nface
    real(rp),intent(in)           :: rsize(npoin),coor(ndim,npoin)          
    integer(ip),intent(in)        :: lface(nnofa,nface),ptoel1(*),ptoel2(*),lsurf(nface)
    integer(ip),pointer           :: ptosi1(:),ptosi2(:),lline(:),lside(:,:)
    integer(ip),intent(inout)     :: lstack(npoin),lptype(2,npoin),lmark(npoin),nstack,ncol       
    integer(ip)                   :: iside,iline,i,iface,isurf,ijump,inosi,ipa,inofa,jsurf,ichk,imark1,imarkm1 
    integer(ip)                   :: ipoin,ip1,ip2,ip3,j,it,itype,istack,ie,ielem,nstackold,ibegin,istackold,ntop
    integer(ip),pointer           :: lmline(:),lmsurf(:),lelem(:)
    integer(4)                    :: istat
    !
    !     This sub mark the points with lmark depending on the BREP geometry:
    !           -for an untouched point:       lmark(ip)=0            
    !           -for a point to be kept:       lmark(ip)=1           
    !           -for a point to be deleted:    lmark(ip)=-1 
    !
    !     Go from corners/cusps to ridges to smooth points     
    ! 

    allocate(lmline(nline),stat=istat)
    call memchk(zero,istat,memor_msh,'LMLINE','mark',lmline) 
    allocate(lmsurf(nsurf),stat=istat)
    call memchk(zero,istat,memor_msh,'LMSURF','mark',lmsurf) 
    allocate(lelem(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEM','mark',lelem) 
    !
    !     Initialization
    !
    nstack=0
    !
    !     Get the corner and cusp points to keep them always 
    !
    do ipoin=1,npoin
       if(lptype(1,ipoin)==ID_CORNER)then
          nstack=nstack+1
          lstack(nstack)=ipoin
          lmark(ipoin)=1_ip
       else if(lptype(1,ipoin)==ID_CUSP)then
          nstack=nstack+1
          lstack(nstack)=ipoin
          lmark(ipoin)=1_ip
       endif
    enddo
    !
    !     Mark each line touched by the corner points
    !
    do istack=1,nstack
       ipoin=lstack(istack)
       do i=ptosi2(ipoin),ptosi2(ipoin+1)-1
          iside=ptosi1(i)
          iline=lline(iside)
          lmline(iline)=1_ip
       enddo
    enddo
    !
    !     Check if some lines have been untouched (there are no corner or 
    !     cusps points on these lines)
    !
    do iside=1,nside
       iline=lline(iside)      
       if(lmline(iline)==0)then
          ip1=lside(1,iside)
          if(lmark(ip1)==0)then
             lmark(ip1)=1_ip
             lmline(iline)=1_ip
             nstack=nstack+1
             lstack(nstack)=ip1
             cycle
          endif
          ip2=lside(2,iside)
          if(lmark(ip2)==0)then
             lmark(ip2)=1_ip
             lmline(iline)=1_ip
             nstack=nstack+1
             lstack(nstack)=ip2
          endif
       endif
    enddo
    !
    !     Loop on point lines 
    !
    istack=0_ip

    do 
       if(istack==nstack)exit
       istack=istack+1
       ipoin=lstack(istack)
       !
       !     If ipoin should be kept, mark the neighbors on the line to be deleted
       !                     
       if(lmark(ipoin)==1)then
          !
          !     Loop on the neighboring sides, lelem only used for speed
          !         
          do ie=ptosi2(ipoin),ptosi2(ipoin+1)-1
             ielem=ptosi1(ie)
             if(lelem(ielem)==0)then
                lelem(ielem)=1_ip   
                do j=1,nnosi
                   ip1=lside(j,ielem)
                   if(lmark(ip1)==0)then 
                      lmark(ip1)=-1_ip
                      !
                      !     Add ip1 to the stack
                      ! 
                      nstack=nstack+1
                      lstack(nstack)=ip1
                      exit    
                   endif
                enddo
             endif
          enddo
       else
          !
          !     ipoin should be deleted, mark the first free neighbor to 
          !     be kept and jump to this neighbor
          !
          do ie=ptosi2(ipoin),ptosi2(ipoin+1)-1
             ielem=ptosi1(ie)
             if(lelem(ielem)==0)then
                lelem(ielem)=1_ip   
                do j=1,nnosi
                   ip1=lside(j,ielem)
                   if(lmark(ip1)==0)then
                      lmark(ip1)=1_ip
                      nstack=nstack+1
                      lstack(nstack)=ip1 
                      exit  
                   endif
                enddo
             endif
          enddo
       endif
    enddo
    !
    !     Mark each surface
    !
    do istack=1,nstack
       ipoin=lstack(istack)
       if(lmark(ipoin)==1)then
          do i=ptoel2(ipoin),ptoel2(ipoin+1)-1
             iface=ptoel1(i)
             isurf=lsurf(iface)
             if(lmsurf(isurf)==0)then
                lmsurf(isurf)=ipoin
             endif
          enddo
       endif
    enddo
    !
    !     Clean up lelem
    ! 
    do iside=1,nside
       lelem(iside)=0_ip
    enddo
    !
    !     Check that all the surfaces have been marked
    !
    do iface=1,nface
       isurf=lsurf(iface)
       if(lmsurf(isurf)==0)then
          ip1=lface(1,iface) 
          if(lmark(ip1)==0)then
             lmark(ip1)=1_ip
             lmsurf(isurf)=ip1
             nstack=nstack+1
             lstack(nstack)=ip1
             cycle
          endif
          ip2=lface(2,iface) 
          if(lmark(ip2)==0)then
             lmark(ip2)=1_ip
             lmsurf(isurf)=ip2
             nstack=nstack+1
             lstack(nstack)=ip2
             cycle
          endif
          ip3=lface(3,iface) 
          if(lmark(ip3)==0)then
             lmark(ip3)=1_ip
             lmsurf(isurf)=ip3
             nstack=nstack+1
             lstack(nstack)=ip3
             cycle
          endif
       endif
    enddo
    !
    !     Loop on the surfaces
    !
    do isurf=1,nsurf

       ipoin=lmsurf(isurf)
       !
       !     DBG
       !
       if(ipoin<=0)then

          cycle
          !
          !     Do we really have deleted a surface?
          !
          do iface=1,nface
             if(lsurf(iface)==isurf)exit
          enddo

          if(iface==nface+1)then
             cycle
          else
             write(*,*)'Error in mark, ipoin<=0'
             stop
          endif
       endif
       lstack(1)=ipoin
       nstack=1_ip
       istackold=1_ip
       ijump=0_ip
       !
       !     Initialize the stack for that surface
       !     Slightly different from after as this point is at the boundary
       !     of the surface and we do not want to check for the surface number
       !     each time. The points must then be interior.
       !
       do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(ie)
          if(lsurf(ielem)==isurf)then
             lelem(ielem)=1_ip   
             do j=1,nnofa
                ip1=lface(j,ielem)
                if(lmark(ip1)==0)then 
                   lmark(ip1)=-1_ip
                   !
                   !     Add ip1 to the stack
                   ! 
                   nstack=nstack+1
                   lstack(nstack)=ip1
                endif
             enddo
          endif
       enddo
       !
       !     Check that some points have been added to the stack
       !     As we add points from a corner, this corner may be overconstrained for that surface
       !
       if(nstack==1)then
          !
          !     Brute force on the ridges
          !     
          do iside=1,nside
             do inosi=1,nnosi 
                ip1=lside(inosi,iside)
                !
                !     Must be a point to be kept
                !
                if(lmark(ip1)==1)then
                   do ie=ptoel2(ip1),ptoel2(ip1+1)-1
                      ielem=ptoel1(ie)
                      if(lsurf(ielem)==isurf)then
                         lelem(ielem)=1_ip   
                         do j=1,nnofa
                            ipa=lface(j,ielem)
                            if(lmark(ipa)==0)then 
                               lmark(ipa)=-1_ip
                               !
                               !     Add ip1 to the stack
                               ! 
                               nstack=nstack+1
                               lstack(nstack)=ipa
                            endif
                         enddo
                      endif
                   enddo
                   !
                   !     Did the stack increase?
                   ! 
                   if(nstack/=1)goto 100

                else if(lmark(ip1)==-1)then
                   do ie=ptoel2(ip1),ptoel2(ip1+1)-1
                      ielem=ptoel1(ie)
                      if(lsurf(ielem)==isurf)then
                         lelem(ielem)=1_ip   
                         do j=1,nnofa
                            ipa=lface(j,ielem)
                            if(lmark(ipa)==0)then 
                               lmark(ipa)=1_ip
                               !
                               !     Add ip1 to the stack
                               ! 
                               nstack=nstack+1
                               lstack(nstack)=ipa
                               goto 100
                            endif
                         enddo
                      endif
                   enddo
                endif
             enddo
          enddo
          !
          !     Did not find other points to mark for that surface.
          !     There are no inner smooth points for that surface ---> go home
          !
          !write(*,*)'No inner points for surface:',isurf
          !
          !     Verifies that all the points of this surface have been marked
          !
          do iface=1,nface
             jsurf=lsurf(iface)
             if(jsurf==isurf)then
                do inofa=1,nnofa
                   ipa=lface(inofa,iface)    
                   if(lmark(ipa)==0)then
                      write(*,*)'Error in mark point not marked',ipa,'isurf=',isurf
                      stop
                   endif
                enddo
             endif
          enddo

          cycle

100       continue

       endif
       !
       !     Loop on the seeds 
       !
200    continue

       !
       !     Loop on the points marked for deletion
       !
       ntop=1_ip
       do
          !write(*,*)'Inner loop on progression'
          !
          !     Loop on the current expansion
          ! 
          istack=istackold
          ! 
          do 

             !write(*,*)'Inner loop : istack=',istack,'nstack=',nstack
             !
             !     Did we reach the end of the stack
             !
             if(istack==nstack)exit
             istack=istack+1
             ipoin=lstack(istack)
             !
             !     If ipoin should be kept, mark the neighbors to be deleted
             !                     
             if(lmark(ipoin)==1)then
                !
                !     Loop on the neighboring elements
                !         
                do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
                   ielem=ptoel1(ie)
                   if(lelem(ielem)==0)then
                      lelem(ielem)=1_ip   
                      do j=1,nnofa
                         ip1=lface(j,ielem)
                         if(lmark(ip1)==0)then 
                            lmark(ip1)=-1_ip
                            !
                            !     Add ip1 to the stack
                            ! 
                            nstack=nstack+1
                            lstack(nstack)=ip1
                         endif
                      enddo
                   endif
                enddo
             else
                !
                !     Remember this position
                !
                ntop=istack 
                !
                !     ipoin should be deleted, mark the first free neighbor to 
                !     be kept and jump to this neighbor
                !
                loop:do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
                   ielem=ptoel1(ie)
                   do j=1,nnofa
                      ip1=lface(j,ielem)
                      if(lmark(ip1)==0)then
                         lmark(ip1)=1_ip
                         istack=nstack 
                         nstack=nstack+1
                         lstack(nstack)=ip1 
                         ijump=1
                         exit loop 
                      endif
                   enddo
                enddo loop
                !
                !     Remember this point as visited if no jump commited
                !
                if(ijump==0)istackold=istack  
             endif
          enddo
          !
          !     We have reached the end of this expansion because we can not add more points 
          !     Once an expansion has ended, the only way to add new points is to go back 
          !     to a deleted point
          !     Did we reach the highest deleted point?  
          !
          if(istackold==ntop)exit 
          !
          !     Reset the jumps
          !     
          ijump=0_ip
       enddo
       !
       !     Expand the points of this patch
       !
       lstack=0_ip
       do iface=1,nface
          jsurf=lsurf(iface)
          if(jsurf==isurf)then
             do inofa=1,nnofa
                ipa=lface(inofa,iface)    
                lstack(ipa)=1_ip
             enddo
          endif
       enddo
       !
       !     See if we have still non marked points 
       !      
       do ipoin=1,npoin
          if(lstack(ipoin)==1)then
             if(lmark(ipoin)==0)then
                imark1=0_ip
                imarkm1=0_ip
                do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
                   ielem=ptoel1(ie)
                   if(lsurf(ielem)==isurf)then
                      do j=1,nnofa
                         ip1=lface(j,ielem)
                         if(ip1/=ipoin)then
                            if(lmark(ip1)==1)then 
                               imark1=imark1+1_ip
                            else if(lmark(ip1)==-1)then
                               imarkm1=imarkm1+1 
                            endif
                         endif
                      enddo
                   endif
                enddo
                if(imark1>imarkm1)then
                   lmark(ipoin)=-1_ip
                else
                   lmark(ipoin)=1_ip
                endif
                lstack(1)=ipoin
                nstack=1_ip
                lstack(2)=ipoin
                nstack=2_ip
                istackold=1_ip
                ijump=0_ip
                goto 200
             endif
          endif
       enddo
    enddo
    !
    !     DBG
    !
    !do ipoin=1,npoin
    !   if(lmark(ipoin)==0)then
    !      write(*,*)'Error in mark, not all the points have been marked, ipoin=',ipoin
    !      stop   
    !   endif
    !enddo

    call memchk(2_ip,istat,memor_msh,'LELEM','mark',lelem)
    deallocate(lelem,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEM','mark',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMLINE','mark',lmline)
    deallocate(lmline,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMLINE','mark',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMSURF','mark',lmsurf)
    deallocate(lmsurf,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMSURF','mark',0_ip)

  end subroutine mark

  subroutine swapsurf(lface,nnofa,nface,coor,npoin,ndim,rnopo,eltoel,niter,ptosi1,&
       ptosi2,lside,nnosi,nside)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       : memor_msh
    use mod_memchk
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim,niter,nnosi,nside
    integer(ip), intent(in)      :: ptosi2(npoin+1),ptosi1(*)
    integer(ip), intent(inout)   :: lface(nnofa,nface),eltoel(nnofa,nface)
    integer(ip), intent(in)      :: lside(nnosi,nside)
    real(rp), intent(in)         :: coor(ndim,npoin),rnopo(ndim,npoin)
    integer(ip)                  :: iface,jface,iopt,j,iter,icont
    integer(ip)                  :: isto,inofa,iside,ip1,ip2,ipc,ineigh,icheck
    integer(ip),pointer          :: lmark(:)
    integer(4)                   :: istat 
    integer(ip)                  :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !
    !     This sub tries to optimize the hole surface triangulation in lface
    !     No update is made for rface as swap recomputes rface
    !
    !
    !     If we want to iterate, allocate an element marked
    !
    if(niter==1)then  

       do iface=1,nface
          do j=1,nnofa
             jface=eltoel(j,iface)
             if(jface/=0)then
                ip1=lface(ltab(1,j),iface)
                ip2=lface(ltab(2,j),iface)
                icheck=1_ip
                do isto=ptosi2(ip1),ptosi2(ip1+1)-1
                   iside=ptosi1(isto)
                   if(lside(1,iside)==ip1)then
                      ipc=lside(2,iside)
                   else
                      ipc=lside(1,iside)
                   endif
                   if(ipc==ip2)then
                      icheck=0_ip               
                      exit
                   endif
                enddo

                if(icheck==1)then

                   call swap(iface,jface,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,iopt)
                endif
             endif
          enddo
       enddo

    else

       allocate(lmark(nface),stat=istat)
       call memchk(zero,istat,memor_msh,'LMARK','swapsurf',lmark) 

       do iter=1,niter
          do iface=1,nface
             if(lmark(iface)==0)then
                icont=0_ip 
                do j=1,nnofa
                   jface=eltoel(j,iface)
                   if(jface/=0)then

                      ip1=lface(ltab(1,j),iface)
                      ip2=lface(ltab(2,j),iface)
                      icheck=1_ip
                      do isto=ptosi2(ip1),ptosi2(ip1+1)-1
                         iside=ptosi1(isto)
                         if(lside(1,iside)==ip1)then
                            ipc=lside(2,iside)
                         else
                            ipc=lside(1,iside)
                         endif
                         if(ipc==ip2)then
                            icheck=0_ip               
                            exit
                         endif
                      enddo
                      if(icheck==1)then

                         iopt=0_ip
                         call swap(iface,jface,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,iopt)
                         !
                         !     Is the swap successfull?
                         ! 
                         if(iopt==1)then
                            icont=icont+1_ip 
                            do inofa=1,nnofa
                               ineigh=eltoel(inofa,iface)
                               if(ineigh/=0)then
                                  lmark(ineigh)=0_ip
                               endif
                            enddo
                            do inofa=1,nnofa
                               ineigh=eltoel(inofa,jface)
                               if(ineigh/=0)then
                                  lmark(ineigh)=0_ip
                               endif
                            enddo

                         endif

                      endif

                   endif
                enddo
                if(icont==0)then
                   lmark(iface)=1_ip
                endif
             endif
          enddo
       enddo

       call memchk(2_ip,istat,memor_msh,'LMARK','swapsurf',lmark)
       deallocate(lmark,stat=istat)
       if(istat/=0) call memerr(2_ip,'LMARK','swapsurf',0_ip)

    endif

  end subroutine swapsurf

  subroutine swap(ielem1,ielem2,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,iopt)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nface,ndim,npoin,nnofa,ielem1,ielem2
    integer(ip),intent(inout) :: lface(nnofa,nface),eltoel(nnofa,nface),iopt 
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin)
    integer(ip)               :: lfacn(nnofa,2),iview1,iview2,p1,p2,p3,p4   
    integer(ip)               :: ineigh11,ineigh12,ineigh21,ineigh22
    real(rp)                  :: Qgeo1,Qgeo2,Qgeo,Qgeotol,rfloc(ndim,2),modul(2),c05,q1,q2,qold
    real(rp)                  :: Qgeo1old,Qgeo2old,Qgeoold,rtolani,rlmax,rlmin,rx,ry,rz,rl1,rl2,rl3,rl4,rl5
    real(rp)                  :: modulold(2),rflocold(ndim,2),q1old,q2old,c00,qtolmax
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This subroutine swaps two elements ielem1 and ielem2 if it is possible
    !     If iqual==1 swaping is done if quality is improved (SHOULD ALWAYS BE CONTROLED BY QUALITY) 
    !     Normals and quality are recomputed for the old faces
    !     It is called by colsmo through swapsurf which do not need updates of rnofa
    ! 
    c05=0.5d+00
    c00=0.0d+00
    Qgeotol=0.95d+00
    Qgeotol=0.88d+00
    Qgeotol=0.85d+00
    qtolmax=1.2d+00
    rtolani=4.0d+00
    !  write(*,*)lface(1,14),lface(2,14),lface(3,14)
    !  write(*,*)lmark(lface(3,14))
    !  write(*,*)'ielem1=',ielem1,'ielem2',ielem2
    !
    !     Find viewer
    !
    if(eltoel(1,ielem1)==ielem2)then
       iview1=1
    else if(eltoel(2,ielem1)==ielem2)then
       iview1=2
    else
       iview1=3
    endif

    if(eltoel(1,ielem2)==ielem1)then
       iview2=1
    else if(eltoel(2,ielem2)==ielem1)then
       iview2=2
    else
       iview2=3
    endif

    p1=lface(iview1,ielem1)
    p2=lface(ltab(1,iview1),ielem1)
    p3=lface(ltab(2,iview1),ielem1)
    p4=lface(iview2,ielem2)

    lfacn(1,1)=p1 
    lfacn(2,1)=p2 
    lfacn(3,1)=p4 

    lfacn(1,2)=p1 
    lfacn(2,2)=p4 
    lfacn(3,2)=p3 
    !
    !     Compute old geometric quality
    !
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rflocold(1,1),ielem1,modulold(1))
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rflocold(1,2),ielem2,modulold(2))
    call chkgeo(lface,nface,nnofa,ielem1,Qgeo1old,rflocold(1,1),rnopo,npoin,ndim)
    call chkgeo(lface,nface,nnofa,ielem2,Qgeo2old,rflocold(1,2),rnopo,npoin,ndim)
    Qgeoold=min(Qgeo1old,Qgeo2old)
    !
    !     Compute new geometric quality
    !
    call gtfnr2(lfacn,2_ip,nnofa,ndim,coor,npoin,rfloc(1,1),1_ip,modul(1))
    if(modul(1)==c00)return
    call chkgeo(lfacn,2_ip,nnofa,1_ip,Qgeo1,rfloc(1,1),rnopo,npoin,ndim)
    if(Qgeo1<Qgeotol)return

    call gtfnr2(lfacn,2_ip,nnofa,ndim,coor,npoin,rfloc(1,2),2_ip,modul(2))
    if(modul(2)==c00)return
    call chkgeo(lfacn,2_ip,nnofa,2_ip,Qgeo2,rfloc(1,2),rnopo,npoin,ndim)
    if(Qgeo2<Qgeotol)return

    !if(Qgeo<Qgeoold)return 
    !
    !     If the context is anisotropic, the isotropic quality has no sense --> go home
    !
    !rx=coor(1,p2)-coor(1,p1)
    !ry=coor(2,p2)-coor(2,p1)
    !rz=coor(3,p2)-coor(3,p1)
    !rl1=sqrt(rx*rx+ry*ry+rz*rz) 
    !rx=coor(1,p3)-coor(1,p1)
    !ry=coor(2,p3)-coor(2,p1)
    !rz=coor(3,p3)-coor(3,p1)
    !rl2=sqrt(rx*rx+ry*ry+rz*rz) 
    !rx=coor(1,p3)-coor(1,p2)
    !ry=coor(2,p3)-coor(2,p2)
    !rz=coor(3,p3)-coor(3,p2)
    !rl3=sqrt(rx*rx+ry*ry+rz*rz) 
    !rx=coor(1,p4)-coor(1,p2)
    !ry=coor(2,p4)-coor(2,p2)
    !rz=coor(3,p4)-coor(3,p2)
    !rl4=sqrt(rx*rx+ry*ry+rz*rz) 
    !rx=coor(1,p4)-coor(1,p3)
    !ry=coor(2,p4)-coor(2,p3)
    !rz=coor(3,p4)-coor(3,p3)
    !rl5=sqrt(rx*rx+ry*ry+rz*rz) 
    !rlmax=max(rl1,rl2,rl3,rl4,rl5) 
    !rlmin=min(rl1,rl2,rl3,rl4,rl5) 
    !if(rlmax/rlmin>rtolani)return
    !
    !     Do we care about increasing the quality?
    !
    !
    !     Compute surface of old elements
    !
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rflocold(1,1),ielem1,modulold(1))
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rflocold(1,2),ielem2,modulold(2))
    !
    !     Divide the area by 2
    !
    modulold(1)=modulold(1)*c05
    modulold(2)=modulold(2)*c05
    !
    !     Compute quality of old elements
    !
    call quality(lface,nface,ielem1,nnofa,coor,ndim,npoin,q1old,modulold(1))
    call quality(lface,nface,ielem2,nnofa,coor,ndim,npoin,q2old,modulold(2))
    qold=max(q1old,q2old)
    !
    !   Is it worth following
    !
    if(qold<qtolmax)return
    !
    !     Divide the area by 2 of the new elements
    !
    modul(1)=modul(1)*c05
    modul(2)=modul(2)*c05
    !
    !  Compute shape quality of the new elements
    !
    call quality(lfacn,2_ip,1_ip,nnofa,coor,ndim,npoin,q1,modul(1))
    if(q1>qold)return
    call quality(lfacn,2_ip,2_ip,nnofa,coor,ndim,npoin,q2,modul(2))
    if(q2>qold)return

    !write(*,*)'swap done'
    !  call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
    !
    !     Swap successfull
    !
    iopt=1_ip
    ineigh11=eltoel(ltab(1,iview1),ielem1)
    ineigh12=eltoel(ltab(2,iview1),ielem1)
    ineigh21=eltoel(ltab(1,iview2),ielem2)
    ineigh22=eltoel(ltab(2,iview2),ielem2)
    !
    !    Update lface
    !
    lface(1,ielem1)=lfacn(1,1)
    lface(2,ielem1)=lfacn(2,1)
    lface(3,ielem1)=lfacn(3,1)
    lface(1,ielem2)=lfacn(1,2)
    lface(2,ielem2)=lfacn(2,2)
    lface(3,ielem2)=lfacn(3,2)
    !
    !     Update eltoel
    !
    eltoel(1,ielem1)=ineigh21
    eltoel(2,ielem1)=ielem2
    eltoel(3,ielem1)=ineigh12

    eltoel(1,ielem2)=ineigh22
    eltoel(2,ielem2)=ineigh11
    eltoel(3,ielem2)=ielem1
    !
    !     Update outside
    !
    if(ineigh21/=0)then
       if(eltoel(1,ineigh21)==ielem2)then
          eltoel(1,ineigh21)=ielem1
       else if(eltoel(2,ineigh21)==ielem2)then
          eltoel(2,ineigh21)=ielem1
       else
          eltoel(3,ineigh21)=ielem1
       endif
    endif

    if(ineigh11/=0)then
       if(eltoel(1,ineigh11)==ielem1)then
          eltoel(1,ineigh11)=ielem2
       else if(eltoel(2,ineigh11)==ielem1)then
          eltoel(2,ineigh11)=ielem2
       else
          eltoel(3,ineigh11)=ielem2
       endif
    endif
    !write(*,*)'from swap' 
    !call chkcons(lface,nnofa,nface,eltoel,lelem,coor,npoin,ndim,lmark)

  end subroutine swap

  subroutine brutefface(lface,nface,nnofa,pnew,coor,ndim,npoin,rnofa,rdist,&
       d1,d2,d3,ipoin,pint,ihost,isurf,lsurf,rsize)
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    use def_kintyp, only       :  ip,rp,lg
    implicit none 
    integer(ip),intent(in)    :: npoin,nface,nnofa,ndim,isurf
    integer(ip),intent(in)    :: lface(nnofa,nface),lsurf(nface)
    integer(ip),intent(inout) :: ihost
    real(rp),intent(in)       :: coor(ndim,npoin),rnofa(ndim,nface),pnew(ndim),rsize
    real(rp),intent(inout)    :: rdist,d1,d2,d3,pint(3)
    integer(ip),intent(inout) :: ipoin
    integer(ip)               :: ip1,ip2,ip3,ipoin1,ipoin2,ipoin3,iface
    real(rp)                  :: rdirx,rdiry,rdirz,csca,cscal,xpro,ypro,zpro,epsil,rl,epsil0,epsil1
    real(rp)                  :: rdist1,rdist2,rdist3,dtot,p1(ndim),p2(ndim),c10,rface(ndim)
    real(rp)                  :: rtx1,rty1,rtz1,rtx2,rty2,rtz2,rtl,pint1(3),pint2(3),pint3(3)
    real(rp)                  :: d1t,d2t,d3t,pintt(3),rdistt,vmin(3),vmax(3),rlen
    real(rp)                  :: xmax,ymax,zmax,xmin,ymin,zmin,coef 

    c10=1.0d+00 
    epsil=1.0d-06
    epsil1=1.0d+00+1.0d-04
    epsil0=1.0d-04
    coef=10.0d+00
    rdist=1.0d+12 

    rlen=coef*rsize
    vmin(1)=pnew(1)-rlen
    vmin(2)=pnew(2)-rlen
    vmin(3)=pnew(3)-rlen
    vmax(1)=pnew(1)+rlen
    vmax(2)=pnew(2)+rlen
    vmax(3)=pnew(3)+rlen

    do iface=1,nface

       if(lsurf(iface)==isurf)then
          !
          !     Get the points of the face
          !
          ip1=lface(1,iface)
          ip2=lface(2,iface)
          ip3=lface(3,iface)
          !
          !     Host face normal
          !
          rface(1)=rnofa(1,iface)
          rface(2)=rnofa(2,iface)
          rface(3)=rnofa(3,iface)
          !
          !     Filter with bbox
          !
          xmin=coor(1,ip1)
          xmax=xmin
          ymin=coor(2,ip1)
          ymax=ymin
          zmin=coor(3,ip1)
          zmax=zmin

          if(coor(1,ip2)<xmin)xmin=coor(1,ip2)
          if(coor(2,ip2)<ymin)ymin=coor(2,ip2)
          if(coor(3,ip2)<zmin)zmin=coor(3,ip2)
          if(coor(1,ip2)>xmax)xmax=coor(1,ip2)
          if(coor(2,ip2)>ymax)ymax=coor(2,ip2)
          if(coor(3,ip2)>zmax)zmax=coor(3,ip2)

          if(coor(1,ip3)<xmin)xmin=coor(1,ip3)
          if(coor(2,ip3)<ymin)ymin=coor(2,ip3)
          if(coor(3,ip3)<zmin)zmin=coor(3,ip3)
          if(coor(1,ip3)>xmax)xmax=coor(1,ip3)
          if(coor(2,ip3)>ymax)ymax=coor(2,ip3)
          if(coor(3,ip3)>zmax)zmax=coor(3,ip3)
          !
          !     Check against vmin,vmax
          !
          if(vmin(1)>xmax)cycle
          if(vmin(2)>ymax)cycle
          if(vmin(3)>zmax)cycle
          if(vmax(1)<xmin)cycle
          if(vmax(2)<ymin)cycle
          if(vmax(3)<zmin)cycle
          !
          !     Set ipoin to zero
          ! 
          ipoin=0_ip
          ipoin1=0_ip
          ipoin2=0_ip
          ipoin3=0_ip
          !
          !     Project the point on the plane of the source
          !
          rdirx=pnew(1)-coor(1,ip1)
          rdiry=pnew(2)-coor(2,ip1)
          rdirz=pnew(3)-coor(3,ip1)

          csca=rface(1)*rdirx+rface(2)*rdiry+rface(3)*rdirz

          xpro=pnew(1)-csca*rface(1)
          ypro=pnew(2)-csca*rface(2)
          zpro=pnew(3)-csca*rface(3)
          !
          !     Compute area to normalize
          !
          p1(1)=coor(1,ip2)-coor(1,ip1)
          p1(2)=coor(2,ip2)-coor(2,ip1)
          p1(3)=coor(3,ip2)-coor(3,ip1)
          p2(1)=coor(1,ip3)-coor(1,ip1)
          p2(2)=coor(2,ip3)-coor(2,ip1)
          p2(3)=coor(3,ip3)-coor(3,ip1)
          call orient3D(p1,p2,rface,dtot,ndim)
          !
          !     Side (ip2,ip3)
          !
          p1(1)=coor(1,ip2)-xpro
          p1(2)=coor(2,ip2)-ypro
          p1(3)=coor(3,ip2)-zpro
          p2(1)=coor(1,ip3)-xpro
          p2(2)=coor(2,ip3)-ypro
          p2(3)=coor(3,ip3)-zpro
          call orient3D(p1,p2,rface,d1t,ndim)
          d1t=d1t/dtot
          !
          !     Side (ip3,ip1)
          !
          p1(1)=coor(1,ip3)-xpro
          p1(2)=coor(2,ip3)-ypro
          p1(3)=coor(3,ip3)-zpro
          p2(1)=coor(1,ip1)-xpro
          p2(2)=coor(2,ip1)-ypro
          p2(3)=coor(3,ip1)-zpro
          call orient3D(p1,p2,rface,d2t,ndim)
          d2t=d2t/dtot
          !
          !     Side (ip1,ip2)
          !
          p1(1)=coor(1,ip1)-xpro
          p1(2)=coor(2,ip1)-ypro
          p1(3)=coor(3,ip1)-zpro
          p2(1)=coor(1,ip2)-xpro
          p2(2)=coor(2,ip2)-ypro
          p2(3)=coor(3,ip2)-zpro
          call orient3D(p1,p2,rface,d3t,ndim)
          d3t=d3t/dtot
          !
          !     Is the projected point inside the source triangle?
          !
          if(d1t>-epsil .and. d2t>-epsil .and. d3t>-epsil)then 
 
             rdistt=abs(csca)
             pintt(1)=xpro
             pintt(2)=ypro
             pintt(3)=zpro

          else
             !
             !     Evaluate distance to lines
             !

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

                rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
                ipoin1=ip2
                pint1(1)=coor(1,ip2)
                pint1(2)=coor(2,ip2)
                pint1(3)=coor(3,ip2)

             else if(cscal>epsil1)then

                p2(1)=pnew(1)-coor(1,ip3)
                p2(2)=pnew(2)-coor(2,ip3)
                p2(3)=pnew(3)-coor(3,ip3)
                rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
                ipoin1=ip3
                pint1(1)=coor(1,ip3)
                pint1(2)=coor(2,ip3)
                pint1(3)=coor(3,ip3)

             else

                p2(1)=coor(1,ip2)+csca*p1(1)-pnew(1)
                p2(2)=coor(2,ip2)+csca*p1(2)-pnew(2)
                p2(3)=coor(3,ip2)+csca*p1(3)-pnew(3)
                rdist1=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 
                pint1(1)=coor(1,ip2)+csca*p1(1)
                pint1(2)=coor(2,ip2)+csca*p1(2)
                pint1(3)=coor(3,ip2)+csca*p1(3)

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

                rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
                ipoin2=ip3
                pint2(1)=coor(1,ip3)
                pint2(2)=coor(2,ip3)
                pint2(3)=coor(3,ip3)

             else if(cscal>epsil1)then

                p2(1)=pnew(1)-coor(1,ip1)
                p2(2)=pnew(2)-coor(2,ip1)
                p2(3)=pnew(3)-coor(3,ip1)
                rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
                ipoin2=ip1
                pint2(1)=coor(1,ip1)
                pint2(2)=coor(2,ip1)
                pint2(3)=coor(3,ip1)

             else

                p2(1)=coor(1,ip3)+csca*p1(1)-pnew(1)
                p2(2)=coor(2,ip3)+csca*p1(2)-pnew(2)
                p2(3)=coor(3,ip3)+csca*p1(3)-pnew(3)
                rdist2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 
                pint2(1)=coor(1,ip3)+csca*p1(1)
                pint2(2)=coor(2,ip3)+csca*p1(2)
                pint2(3)=coor(3,ip3)+csca*p1(3)

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

                rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
                ipoin3=ip1
                pint3(1)=coor(1,ip1)
                pint3(2)=coor(2,ip1)
                pint3(3)=coor(3,ip1)

             else if(cscal>epsil1)then

                p2(1)=pnew(1)-coor(1,ip2)
                p2(2)=pnew(2)-coor(2,ip2)
                p2(3)=pnew(3)-coor(3,ip2)
                rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))  
                ipoin3=ip2
                pint3(1)=coor(1,ip2)
                pint3(2)=coor(2,ip2)
                pint3(3)=coor(3,ip2)

             else

                p2(1)=coor(1,ip1)+csca*p1(1)-pnew(1)
                p2(2)=coor(2,ip1)+csca*p1(2)-pnew(2)
                p2(3)=coor(3,ip1)+csca*p1(3)-pnew(3)
                rdist3=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3)) 
                pint3(1)=coor(1,ip1)+csca*p1(1)
                pint3(2)=coor(2,ip1)+csca*p1(2)
                pint3(3)=coor(3,ip1)+csca*p1(3)

             endif
             !
             !     Take the min over the side distances
             !
             if(rdist1<rdist2)then 
                rdistt=rdist1
                ipoin=ipoin1
                pint(1)=pint1(1)
                pint(2)=pint1(2)
                pint(3)=pint1(3)
             else 
                rdistt=rdist2
                ipoin=ipoin2
                pint(1)=pint2(1)
                pint(2)=pint2(2)
                pint(3)=pint2(3)
             endif

             if(rdist3<rdist)then 
                rdistt=rdist3
                ipoin=ipoin3
                pint(1)=pint3(1)
                pint(2)=pint3(2)
                pint(3)=pint3(3)
             endif

          endif

          if(rdistt<rdist)then 

             ihost=iface
             pint(1)=pintt(1)
             pint(2)=pintt(2)
             pint(3)=pintt(3)
             d1=d1t
             d2=d2t
             d3=d3t

          endif

       endif

    enddo


  end subroutine brutefface


end module mod_srftol
