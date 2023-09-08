module mod_surf

  use mod_cart
  use mod_sort
  use mod_srftol
  use mod_clean

contains

  subroutine mshsrf(ndim,nnofa,nface,npoin,nboup,nsurf,nline,nnosi,nside,&
       rsuni,nblay,rblay,tolscal,lsmark,rtol,lcell,ncell,lface,rsize,coor,&
       lsurf,lside,lline)
    use def_kintyp, only :  ip,rp,lg
    use mod_memchk
    use def_meshin, only    : memor_msh
    use def_meshin, only    :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)        :: ndim,nnofa,nnosi,nblay,ncell
    real(rp),intent(in)           :: rsuni,rblay(nblay),tolscal,rtol
    integer(ip),intent(inout)     :: nface,npoin,nsurf,nline,nside,nboup
    integer(ip),pointer           :: lfold(:,:),eltoelold(:,:),lptypeold(:,:),lsurfold(:)
    integer(ip),pointer           :: lsmark(:),lstotr2(:),lblmsh(:),lstotr2old(:),lsold(:,:)
    real(rp),pointer              :: coorold(:,:),rnofaold(:,:),rmax(:)
    real(rp)                      :: tolcol,tolcolr,tolref,tolrefr
    type(cell)                    :: lcell(ncell)
    integer(ip),pointer           :: lface(:,:),lsurf(:),lside(:,:),lline(:)
    integer(ip),pointer           :: ptoel1(:),ptoel2(:),eltoel(:,:),llinold(:)
    integer(ip),pointer           :: ptosi1(:),ptosi2(:),sitosiold(:,:)
    integer(ip),pointer           :: ptosi1old(:),ptosi2old(:)
    integer(ip),pointer           :: ptoel1old(:),ptoel2old(:)
    integer(ip),pointer           :: lcart(:),lpsur(:),lpsid(:),lptype(:,:)
    integer(ip),pointer           :: lstof(:),lftoed(:)
    real(rp),pointer              :: rnopo(:,:),rnofa(:,:),rsize(:)
    real(rp),pointer              :: coor(:,:)
    integer(ip)                   :: npold,nfold,nsold
    integer(ip)                   :: iface,ipoin,ipo,iside,isurf,iline,nsid
    integer(4)                    :: istat
    !
    !     Basic database used in first part collapse
    !
    !
    !
    !     Database for the points:
    ! 
    !     - rnopo          normal at point
    !     - lcart          cartesian element containing the point
    !     - coor           coordinates
    !     - rsize          size at the point
    !     - lpsur          triangle of the old surface mesh
    !     - lpsid          sides of the old surface mesh
    !
    !
    !     Database for the triangles:
    !
    !     - lsurf          surface number
    !     - eltoel         neighboring element
    !     - rnofa          face normal
    !     - lface          triangle connectivity
    !
    !
    !
    !    
    !
    !     Extended database used for frontal refinement
    !
    !
    !       For the points:
    !      
    !       lpofa            pointer to active faces in front  
    !           >0           the beginning of the faces
    !           =0           untouched
    !           =-1          front has already gone through this point          
    !       lmark            help array
    !       
    !
    !       For the front:
    !
    !
    !       lfront           the edges of the front 
    !       lfhole           small help array for deleted front edges 
    !       lfapo            the adjacent active edges of the front 
    !       lfahol           small help array for deleted adjacent front edges 
    !
    !
    nullify(ptoel1,ptoel2,eltoel,ptosi1,ptosi2,sitosiold,lstof)
    !
    !     First check data
    !
    call chkorient(nface,nnofa,npoin,lface)
    !
    !     Check for degenerate points
    !
    call colpnt(lface,nnofa,nface,coor,ndim,npoin,nnosi,nside,lside,&
         lsurf,nsurf,lline,nline)
    !
    !     Then reorder faces with respect to surfaces
    !  
    allocate(lstotr2(nsurf+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTOTR2','mshsrf',lstotr2) 
    call renuface(nface,nnofa,nsurf,lface,lsurf,lstotr2)
    !
    !     Allocate the database (containing element in cartesian mesh, size at point)
    ! 
    allocate(lcart(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LCART','mshsrf',lcart) 
    allocate(rsize(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'RSIZE','mshsrf',rsize) 
    allocate(rnopo(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOPO','mshsrf',rnopo) 
    allocate(rnofa(ndim,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOFA','mshsrf',rnofa) 
    allocate(lpsur(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPSUR','mshsrf',lpsur) 
    allocate(lpsid(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPSID','mshsrf',lpsid) 
    !
    !     Allocate temporal arrays
    !
    allocate(lptype(2,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LFOLD','mshsrf',lptype) 
    allocate(rmax(nsurf),stat=istat)
    call memchk(zero,istat,memor_msh,'RMAX','mshsrf',rmax) 
    allocate(lblmsh(nsurf),stat=istat)
    call memchk(zero,istat,memor_msh,'LBLMSH','mshsrf',lblmsh) 
    !
    !     Do we have sides?
    !
    if(nside>0)then
       !
       !     Count and renumber the boundary points
       !
       call renucorner(nface,ndim,npoin,nnofa,nnosi,nside,lface,coor,lside,nboup,&
            lptype,lline,lsurf,tolscal)

    else
       !
       !     We do not have sides, there is only one surface
       !
       do ipoin=1,npoin
          lptype(1,ipoin)=ID_SMOOTH 
       enddo
       !
       !     Set boundary point number
       !
       nboup=0_ip

    endif
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get the face normals
    !
    call gtfnrl(lface,nface,nnofa,ndim,coor,npoin,rnofa)
    !
    !     Get the point normals
    !
    call gtpnrl(nface,rnopo,npoin,ndim,ptoel1,ptoel2,rnofa)
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Do we have sides?
    !
    if(nside>0)then
       !
       !     Get the sides surrounding points with the new renumerotation
       !
       call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)
       !
       !     Get the sides surrounding sides
       !
       call sitosa(lside,nside,npoin,nnosi,ptosi1,ptosi2,sitosiold)
       !
       !     Get the side to face pointer
       !
       call sitofa(lside,nside,nnosi,lface,nnofa,nface,npoin,ptoel1,ptoel2,lstof)
       !
       !     Copy ptosi1 && ptosi2
       !
       allocate(ptosi2old(npoin+1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOSI2OLD','mshsrf',ptosi2old) 
       allocate(ptosi1old(ptosi2(npoin+1)-1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOSI1OLD','mshsrf',ptosi1old) 
       ptosi2old=ptosi2 
       ptosi1old=ptosi1 
       !
       !     Check line data around corner
       !
       !call chkline(npoin,nside,nnosi,lptype,ptosi1,ptosi2,lline,lside)
    else    

       allocate(ptosi2(npoin+1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOSI2','mshsrf',ptosi2) 
       allocate(ptosi1(1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOSI1','mshsrf',ptosi1)
       allocate(sitosiold(1,1),stat=istat)
       call memchk(zero,istat,memor_msh,'SITOSIOLD','mshsrf',sitosiold) 
       allocate(lstof(1),stat=istat)
       call memchk(zero,istat,memor_msh,'LSTOF','mshsrf',lstof) 
       allocate(ptosi2old(npoin+1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOSI2OLD','mshsrf',ptosi2old) 
       allocate(ptosi1old(1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOSI1OLD','mshsrf',ptosi1old) 

    endif
    !
    !     Do we need to smooth out the initial surface
    ! 
    write(*,*)'Smoothing the initial surface' 
    !call presmoo(eltoel,lface,nface,nnofa,ndim,coor,npoin,rnofa,ptoel1,ptoel2,&
    !     nboup,lsurf,lptype,rnopo,nsurf,lside,nside,lline,nline)
    !
    !     Give a surface number to smooth points, a line number to ridge and cusp points
    !
    do ipoin=1,npoin
       if(lptype(1,ipoin)==ID_SMOOTH)then
          isurf=lsurf(ptoel1(ptoel2(ipoin)))
          lptype(2,ipoin)=isurf
       else if(lptype(1,ipoin)==ID_CUSP)then    
          iline=lline(ptosi1(ptosi2(ipoin)))
          lptype(2,ipoin)=iline
       else if(lptype(1,ipoin)==ID_RIDGE)then
          iline=lline(ptosi1(ptosi2(ipoin)))
          lptype(2,ipoin)=iline
       endif
    enddo
    !
    !     Then copy the old mesh (faces,coor,face normals, neigh to neigh connection, ...)
    !
    allocate(lfold(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFOLD','mshsrf',lfold) 
    allocate(coorold(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'COOROLD','mshsrf',coorold) 
    allocate(rnofaold(ndim,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOFAOLD','mshsrf',rnofaold) 
    allocate(eltoelold(ndim,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'ELTOELOLD','mshsrf',eltoelold) 
    allocate(lptypeold(2,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTYPEOLD','mshsrf',lptypeold) 
    allocate(lsurfold(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSURFOLD','mshsrf',lsurfold) 
    allocate(lstotr2old(nsurf+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTOTR2OLD','mshsrf',lstotr2old) 
    allocate(ptoel2old(npoin+1),stat=istat)
    call memchk(zero,istat,memor_msh,'PTOEL2OLD','mshsrf',ptoel2old) 
    allocate(ptoel1old(ptoel2(npoin+1)-1),stat=istat)
    call memchk(zero,istat,memor_msh,'PTOEL1OLD','mshsrf',ptoel1old) 

    nfold=nface
    npold=npoin
    nsold=nside

    do iface=1,nface
       lfold(1,iface)=lface(1,iface)
       lfold(2,iface)=lface(2,iface)
       lfold(3,iface)=lface(3,iface)
       lsurfold(iface)=lsurf(iface)
       rnofaold(1,iface)=rnofa(1,iface)
       rnofaold(2,iface)=rnofa(2,iface)
       rnofaold(3,iface)=rnofa(3,iface)
       eltoelold(1,iface)=eltoel(1,iface)
       eltoelold(2,iface)=eltoel(2,iface)
       eltoelold(3,iface)=eltoel(3,iface)
    enddo

    do ipoin=1,npoin
       coorold(1,ipoin)=coor(1,ipoin)
       coorold(2,ipoin)=coor(2,ipoin)
       coorold(3,ipoin)=coor(3,ipoin)
       lptypeold(1,ipoin)=lptype(1,ipoin)
       lptypeold(2,ipoin)=lptype(2,ipoin)
    enddo

    do ipoin=1,npoin+1
       ptoel2old(ipoin)=ptoel2(ipoin)
    enddo

    do ipoin=1,ptoel2(npoin+1)-1
       ptoel1old(ipoin)=ptoel1(ipoin)
    enddo

    do isurf=1,nsurf+1
       lstotr2old(isurf)=lstotr2(isurf) 
    enddo
    !
    !     Do we have sides?
    !
    if(nside>0)then

       allocate(lsold(nnosi,nside),stat=istat)
       call memchk(zero,istat,memor_msh,'LSIDEOLD','mshsrf',lsold) 
       allocate(llinold(nside),stat=istat)
       call memchk(zero,istat,memor_msh,'LLINOLD','mshsrf',llinold) 
       do iside=1,nside
          lsold(1,iside)=lside(1,iside)
          lsold(2,iside)=lside(2,iside)
          llinold(iside)=lline(iside)
       enddo
       !
       !     Initialize lpsid
       ! 
       do ipoin=1,npoin
          if(lptype(1,ipoin)==ID_CUSP .or. lptype(1,ipoin)==ID_RIDGE)then
             lpsid(ipoin)=ptosi1(ptosi2(ipoin))
          endif
       enddo

    else

       allocate(lsold(1,1),stat=istat)
       call memchk(zero,istat,memor_msh,'LSIDEOLD','mshsrf',lsold) 
       allocate(llinold(1),stat=istat)
       call memchk(zero,istat,memor_msh,'LLINOLD','mshsrf',llinold) 

    endif
    !
    !     Get the elements of the cart grid
    !
    write(*,*)'Interpolating the cartesian mesh'
    call intercart(coor,npoin,ndim,lface,nface,nnofa,lcell,ncell,lcart,ptoel1,ptoel2,rtol)
    !
    !     Get the size of the points
    ! 
    write(*,*)'Getting the point size'
    call gtsize(ncell,lcell,npoin,rsize,lcart,ndim,coor)
    !
    !     Initialize lpsur
    ! 
    do ipoin=1,npoin
       lpsur(ipoin)=ptoel1(ptoel2(ipoin))
    enddo
    !
    !     Define the tolerance for the edge length in the normalized space
    !     The classical measure is [1/sqrt(2),sqrt(2)] 
    !     We want the final mesh to be twice as big before optimization
    !     to leave some room to the optimization
    !     However, the ridges SHOULD be consistent with the size distribution
    !
    !     Tolerance on length for collapsing for smooth edges 
    !
    !tollen=1.0d+00/sqrt(2.0d+00)   too permissive, risk run over
    !tolcol=2.0d+00/sqrt(2.0d+00)
    tolcol=1.0d+00
    !
    !     Tolerance on length for collapsing for ridges, as desired 
    !
    tolcolr=1.0d+00/sqrt(2.0d+00)            
    !
    !     Tolerance on length for refining for smooth edges, twice as big as desired 
    !
    tolref=2.0d+00*sqrt(2.0d+00)
    !
    !     Tolerance on length for refining for ridges, as desired 
    !
    tolrefr=sqrt(2.0d+00)
    !
    !     Here begins the action
    !
    !
    !     Collapse the edges first
    ! 
    call colglo(coor,ndim,npoin,lface,nface,nnofa,lptype,ptoel1,ptoel2,nboup,rsize,&
         lcart,rnopo,nsurf,lsurf,eltoel,lpsur,nline,nside,nnosi,rmax,tolcol,tolcolr,&
         tolref,tolrefr,lpsid,lside,lline,lsmark,lstotr2,rnofa,ptosi1,ptosi2,eltoelold,&
         lsurfold,rnofaold,npold,lfold,lsold,nsold,lstof,nfold)
    !
    !     Write intermediate file
    !  
    !call outinter(nnofa,nface,npoin,ndim,lface,coor,lsurf,lside,lline,nside,nnosi,nline)
    !
    !     Refine the mesh
    !
    call refglo(nface,nnofa,npoin,ndim,lfold,nfold,coorold,npold,lcell,ncell,eltoelold,&
         rnofaold,nsurf,nboup,nside,nnosi,lptypeold,lsurfold,rmax,tolref,tolrefr,tolcol,&
         nsold,rsuni,ptosi1old,ptosi2old,lline,lside,nblay,rblay,tolscal,lsmark,rtol,tolcolr,&
         lstotr2,lblmsh,lface,coor,rnopo,lpsur,lptype,lpsid,lcart,rsize,lsurf,lsold,lstof,&
         sitosiold,ptosi1,ptosi2,nline,llinold,ptoel1old,ptoel2old)
    !
    !     Optimize the mesh
    !  
    call optglo(lface,nnofa,nface,lptype,npoin,coor,ndim,rnopo,eltoelold,npold,&
         nfold,coorold,rnofaold,lpsur,rsize,lfold,nsurf,lsurf,lsurfold,ptoel1,ptoel2,&
         eltoel,lblmsh,lstotr2,ptosi1,ptosi2,lside,nnosi,nside)
    !
    !     Reposition the points with a higher order recovery
    !
    call horder(nnofa,npoin,lfold,nfold,coor,coorold,npold,rnofaold,ndim,&
         lpsur,eltoelold,lstotr2old,lsurfold,nsurf)
    !
    !     Output surface mesh
    !    
    write(*,*)'Dumping surface mesh after optimization in outface.msh'
    call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)
    !
    call memchk(2_ip,istat,memor_msh,'LLINOLD','mshsrf',llinold)
    deallocate(llinold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LLINOLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSOLD','mshsrf',lsold)
    deallocate(lsold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSOLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI1','mshsrf',ptosi1)
    deallocate(ptosi1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI1','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI2','mshsrf',ptosi2)
    deallocate(ptosi2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI2','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI1OLD','mshsrf',ptosi1old)
    deallocate(ptosi1old,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI1OLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI2OLD','mshsrf',ptosi2old)
    deallocate(ptosi2old,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI2OLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'SITOSIOLD','mshsrf',sitosiold)
    deallocate(sitosiold,stat=istat)
    if(istat/=0) call memerr(2_ip,'SITOSIOLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTOF','mshsrf',lstof)
    deallocate(lstof,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTOF','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1OLD','mshsrf',ptoel1old)
    deallocate(ptoel1old,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1OLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2OLD','mshsrf',ptoel2old)
    deallocate(ptoel2old,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2OLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','mshsrf',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','mshsrf',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','mshsrf',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPSID','mshsrf',lpsid)
    deallocate(lpsid,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPSID','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSURFOLD','mshsrf',lsurfold)
    deallocate(lsurfold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSURFOLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTYPEOLD','mshsrf',lptypeold)
    deallocate(lptypeold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTYPEOLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOELOLD','mshsrf',eltoelold)
    deallocate(eltoelold,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOELOLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOFAOLD','mshsrf',rnofaold)
    deallocate(rnofaold,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOFAOLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'COODROLD','mshsrf',coorold)
    deallocate(coorold,stat=istat)
    if(istat/=0) call memerr(2_ip,'COOROLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFOLD','mshsrf',lfold)
    deallocate(lfold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFOLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTYPE','mshsrf',lptype)
    deallocate(lptype,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTYPE','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RMAX','mshsrf',rmax)
    deallocate(rmax,stat=istat)
    if(istat/=0) call memerr(2_ip,'RMAX','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTOTR2OLD','mshsrf',lstotr2old)
    deallocate(lstotr2old,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTOTR2OLD','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTOTR2','mshsrf',lstotr2)
    deallocate(lstotr2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTOTR2','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LBLMSH','mshsrf',lblmsh)
    deallocate(lblmsh,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBLMSH','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOPO','mshsrf',rnopo)
    deallocate(rnopo,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOPO','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCART','mshsrf',lcart)
    deallocate(lcart,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCART','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPSUR','mshsrf',lpsur)
    deallocate(lpsur,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPSUR','mshsrf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOFA','mshsrf',rnofa)
    deallocate(rnofa,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOFA','mshsrf',0_ip)

  end subroutine mshsrf

  subroutine colglo(coor,ndim,npoin,lface,nface,nnofa,lptype,ptoel1,ptoel2,&
       nboup,rsize,lcart,rnopo,nsurf,lsurf,eltoel,lpsur,nline,nside,&
       nnosi,rmax,tolcol,tolcolr,tolref,tolrefr,lpsid,lside,lline,lsmark,lstotr2,&
       rnofa,ptosi1,ptosi2,eltoelold,lsurfold,rnofaold,npold,lfold,lsold,&
       nsold,lstof,nfold)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only     :  memor_msh
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)        :: ndim,nnofa,nline,nnosi,nsurf
    integer(ip),intent(in)        :: npold,nsold,nfold
    integer(ip),intent(in)        :: lfold(nnofa,nfold),lsurfold(nfold)
    integer(ip),intent(in)        :: lsold(nnosi,nsold)
    integer(ip),intent(in)        :: eltoelold(nnofa,nfold),lstof(nsold)
    real(rp),intent(in)           :: rnofaold(nnofa,nfold)
    integer(ip),intent(inout)     :: lstotr2(nsurf+1)
    integer(ip),intent(inout)     :: nface,npoin,nboup
    integer(ip),intent(inout)     :: lptype(2,npoin),lsurf(nface)
    integer(ip),intent(inout)     :: lface(nnofa,nface),lcart(npoin)
    integer(ip),pointer           :: ptoel1(:),ptoel2(:),eltoel(:,:),lelemold(:)
    integer(ip),pointer           :: ptosi1(:),ptosi2(:),ledge(:,:),ledg2(:)
    integer(ip),pointer           :: lside(:,:),lline(:),lsmark(:)
    real(rp),pointer              :: redge(:)
    integer(ip),intent(inout)     :: nside
    integer(ip),intent(inout)     :: lpsur(npoin),lpsid(npoin)
    real(rp),   intent(in)        :: tolcol,tolcolr,tolref,tolrefr
    real(rp),   intent(inout)     :: rsize(npoin),rnofa(ndim,nface)
    real(rp),   intent(inout)     :: coor(ndim,npoin),rnopo(ndim,npoin),rmax(nsurf)
    integer(ip),pointer           :: lstack(:),lmark(:),lelem(:),lelemloc(:)
    integer(ip),pointer           :: lfloc(:,:),lfnew(:,:),lstotr2new(:)
    integer(ip)                   :: nstack,ncol,ipoin,iface,isurf,nfloc,ipnt 
    integer(ip)                   :: iboun,ie,ielem,nfoldt,iside,j,iter,niter,nfnew 
    integer(ip)                   :: isid1,isid2,ipmin1,ipmin2,ipmax1,ipmax2,nsid0 
    integer(4)                    :: istat
    real(rp)                      :: rmaxloc,tolgeo 
    !
    !     This subroutine drives the collapsing by marking and collapsing the surface mesh 
    !
    !
    !     Allocate local arrays
    !
    allocate(lstack(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','colglo',lstack) 
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','colglo',lmark) 
    allocate(lelem(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEM','colglo',lelem) 
    allocate(lelemloc(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEMLOC','colglo',lelemloc) 
    allocate(lfloc(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFLOC','colglo',lfloc) 
    allocate(lfnew(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFNEW','colglo',lfnew) 
    allocate(lstotr2new(nsurf+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTOTR2NEW','colglo',lstotr2new) 
    allocate(ledge(2,3*nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDGE','reglo',ledge) 
    allocate(redge(3*nface),stat=istat)
    call memchk(zero,istat,memor_msh,'REDGE','reglo',redge) 
    allocate(ledg2(3*nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDG2','reglo',ledg2) 
    allocate(lelemold(max(nfold,nsold)),stat=istat)           ! Take the max of faces and sides as used by gthostsid also
    call memchk(zero,istat,memor_msh,'LELEMOLD','refglo',lelemold) 
    !
    !     DBG
    !    
    !call outerrort(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,&
    !        lline,nline)
    !
    !     Initialize exit condition
    !   
    nfoldt=nface
    !
    !     Initialize tolgeo
    !
    tolgeo=0.8d+00
    !
    !     Initialize lstotr2new
    !
    lstotr2new(1)=1_ip
    !
    !     Iter on collapsing iterations
    !     Each iteration is meant to divide the point mesh size by 2
    !     in a multilevel approach
    !
    iter=1_ip

    do 
       write(*,*)'Coarsening Iteration',iter,'nface=',nface
       !
       !     Initialize new face number for this round
       ! 
       nfnew=0_ip
       !
       !     Get the faces surrounding the points for the whole triangulation
       !
       call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
       !
       !     Do we have sides?
       !
       if(nside>0)then
          !
          !     Get the sides surrounding points
          ! 
          call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)

       endif
       !
       !     Mark the points to be collapsed on the whole triangulation with lmark
       !
       call mark(lstack,nstack,ptoel1,ptoel2,nnofa,lface,nface,npoin,&
            lptype,lmark,ncol,coor,ndim,rsize,nline,nside,nnosi,nsurf,lsurf,&
            ptosi1,ptosi2,lline,lside)
       !call outmark(lface,nnofa,nface,lelem,coor,npoin,ndim,lmark)
       !
       !     Collapse edges on ridges only 
       !
       if(nboup/=0)then
          !
          !     Collapse the ridges
          ! 
          call colridg(lface,nnofa,nface,lptype,lmark,npoin,lelem,coor,ndim,ptoel1,&
               ptoel2,rnopo,tolcolr,rsize,ptosi1,ptosi2,lside,nnosi,nside,lline,&
               tolrefr,tolgeo)

       endif
       !
       !     Loop on surface patches
       !
       do isurf=1,nsurf
          !
          !     Copy global to local faces
          !
          nfloc=0_ip
          do iface=lstotr2(isurf),lstotr2(isurf+1)-1
             nfloc=nfloc+1_ip
             lfloc(1,nfloc)=lface(1,iface) 
             lfloc(2,nfloc)=lface(2,iface) 
             lfloc(3,nfloc)=lface(3,iface)
             lelemloc(nfloc)=lelem(iface)
          enddo
          !
          !     Correct geometry for the boundary points to be compatible 
          !     with the current patch
          !  
          call corgeo(nfloc,nnofa,lfloc,npoin,eltoelold,lsurfold,isurf,nfold,lside,&
               nnosi,nside,rnopo,rnofaold,ndim,lpsur,lptype,npold,lfold,lsurf,lelemold,&
               lsold,nsold,lstof,lpsid)
          !
          !     Collapse edges inside each patch
          !
          call colsmo(lfloc,nnofa,nfloc,lptype,lmark,npoin,lelemloc,coor,ndim,ptoel1,ptoel2,&
               rnopo,tolcol,rsize,eltoel,rmaxloc,tolref,ledge,redge,ledg2,lstack,tolgeo)
          !
          !     Remember the largest edge for possible refinement
          !
          rmax(isurf)=rmaxloc
          ! 
          !     Copy local to global   (The mesh can only be smaller)
          !
          do iface=1,nfloc
             if(lelemloc(iface)/=-1)then
                nfnew=nfnew+1
                lfnew(1,nfnew)=lfloc(1,iface) 
                lfnew(2,nfnew)=lfloc(2,iface) 
                lfnew(3,nfnew)=lfloc(3,iface)
                if(lmark(lfnew(1,nfnew))==-2)then
                   write(*,*)'error colglo 1'
                   write(*,*)lelemloc(iface)
                   stop
                endif
                if(lmark(lfnew(2,nfnew))==-2)then
                   write(*,*)'error colglo 2'
                   write(*,*)lelemloc(iface)
                   stop
                endif
                if(lmark(lfnew(3,nfnew))==-2)then
                   write(*,*)'error colglo 3'
                   write(*,*)lelemloc(iface)
                   stop
                endif

                lelem(nfnew)=0_ip
                lsurf(nfnew)=isurf
             endif
             lelemloc(iface)=0_ip
          enddo
          !
          !     Update lstotr2new
          !
          lstotr2new(isurf+1)=nfnew+1_ip 

          !call outerror(lfnew,nnofa,nfnew,lelem,coor,npoin,ndim)

       enddo
       !
       !     Copy lfnew in lface
       !
       nface=nfnew
       do iface=1,nfnew
          lface(1,iface)=lfnew(1,iface)
          lface(2,iface)=lfnew(2,iface)
          lface(3,iface)=lfnew(3,iface)
       enddo
       !
       !     Copy lstotr2new
       !
       do isurf=2,nsurf+1        
          lstotr2(isurf)=lstotr2new(isurf)
       enddo
       !
       !     Compact the points (Should be done at global level as points are global)
       !
       call compact2(nface,npoin,lface,nnofa,coor,ndim,lmark,lcart,rsize,lptype,&
            rnopo,lpsur,nnosi,nside,lpsid,lline,lside,lsmark,ptosi2,ptosi1)
       !
       !     Did we do something?
       !
       !     call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
       if(nface==nfoldt)exit
       nfoldt=nface
       iter=iter+1_ip
       !
       !     Optimize the whole mesh by swapping
       !     Need to call swapglo after compacting the points
       !
       !call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
       niter=2_ip
       !call swapglo(lface,nnofa,nface,npoin,coor,ndim,nsurf,lsurf,niter,&
       !     ptoel1,ptoel2,eltoel,rnofa,lside,nnosi,nside,ptosi1,ptosi2)
       !call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
       call outerrort(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,&
            lline,nline)
       !
       !     Relax geometrical tolerance
       ! 
       tolgeo=tolgeo*0.8d+00 


    enddo
    !
    !     Check that ridge points still have ridges
    !     and  corner points still have at least two ridges
    do ipoin=1,npoin 
       if(lptype(1,ipoin)==ID_RIDGE)then
          !
          !     How many ridges do we have?
          !  
          if(ptosi2(ipoin+1)==ptosi2(ipoin))then
             !
             !     We do not have ridges anymore --> change the type
             !
             lptype(1,ipoin)=ID_SMOOTH

          else if(ptosi2(ipoin+1)==ptosi2(ipoin)+1)then
             !
             !     We only have one ridge left --> change type
             !
             lptype(1,ipoin)=ID_CUSP

          else if(ptosi2(ipoin+1)==ptosi2(ipoin)+2)then
             !
             !     Everything ok
             !
          else
             write(*,*)'strange ridge point in colglo'
             stop
          endif

       else if(lptype(1,ipoin)==ID_CORNER)then
          !
          !     How many ridges do we have?
          !  
          if(ptosi2(ipoin+1)==ptosi2(ipoin)+1)then
             !
             !     We only have one ridge left --> change the type
             !
             lptype(1,ipoin)=ID_CUSP
          endif
       endif
    enddo

    nsid0=nside
    nside=0_ip
    do iside=1,nsid0
       if(lside(1,iside)/=-1)then
          nside=nside+1
          lside(1,nside)=lside(1,iside)
          lside(2,nside)=lside(2,iside)
       endif
    enddo

    write(*,*)'Dumping outerrort.msh after collapsing'
    !call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
    call outerrort(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,&
         lline,nline)
    !
    !     Free memory
    !
    call memchk(2_ip,istat,memor_msh,'LELEMOLD','colglo',lelemold)
    deallocate(lelemold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEMOLD','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTOTR2NEW','colglo',lstotr2new)
    deallocate(lstotr2new,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTOTR2NEW','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFLOC','colglo',lfloc)
    deallocate(lfloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFLOC','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFNEW','colglo',lfnew)
    deallocate(lfnew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFNEW','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LELEMLOC','colglo',lelemloc)
    deallocate(lelemloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEMLOC','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LELEM','colglo',lelem)
    deallocate(lelem,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEM','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','colglo',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','colglo',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','colglo',0_ip)

    if(nside>0)then
       call memchk(2_ip,istat,memor_msh,'PTOSI1','colglo',ptosi1)
       deallocate(ptosi1,stat=istat)
       if(istat/=0) call memerr(2_ip,'PTOSI1','colglo',0_ip)
       call memchk(2_ip,istat,memor_msh,'PTOSI2','colglo',ptosi2)
       deallocate(ptosi2,stat=istat)
       if(istat/=0) call memerr(2_ip,'PTOSI2','colglo',0_ip)
    endif
    call memchk(2_ip,istat,memor_msh,'REDGE','refglo',redge)
    deallocate(redge,stat=istat)
    if(istat/=0) call memerr(2_ip,'REDGE','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDG2','refglo',ledg2)
    deallocate(ledg2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDG2','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGE','refglo',ledge)
    deallocate(ledge,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGE','refglo',0_ip)

  end subroutine colglo

  subroutine colridg(lface,nnofa,nface,lptype,lmark,npoin,lelem,coor,ndim,ptoel1,&
       ptoel2,rnopo,tollen,rsize,ptosi1,ptosi2,lside,nnosi,nside,lline,tolrefr,tolgeo)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,nface,npoin
    integer(ip),intent(in)    :: nnosi,nside
    integer(ip),intent(inout) :: lface(nnofa,nface)
    integer(ip),intent(in)    :: lptype(2,npoin)
    integer(ip),pointer       :: lside(:,:),ptosi1(:),ptosi2(:),lline(:)
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin),tollen 
    real(rp),intent(in)       :: rsize(npoin),tolrefr,tolgeo  
    integer(ip),intent(inout) :: lelem(nface)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:) 
    integer(ip),intent(inout) :: lmark(npoin)
    integer(ip)               :: j,ip1,ip2,iedg,jedg,nedg,ledg(2,100),icol,ipoin
    integer(ip)               :: icollapse,ichk,is,iside,iline,ip3,ip4,nsid1,nsid2
    integer(ip)               :: ipa,ipb,ipmax,ipmin,ipmina,ipmaxa,ifound
    real(rp)                  :: rlen,redg(100) 
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub collapses the short ridges of the whole surface triangulation
    !
    !
    !     Two loops, if collapse is controled or aggresive 
    !
    do icollapse=1,2

       do iside=1,nside
          !
          !     Get the edges to collapse 
          !  
          ip1=lside(1,iside)  
          ip2=lside(2,iside) 
          !
          !     Has this edge been already collapsed?
          !
          if(ip1==-1)cycle  
          ichk=0_ip
          if(icollapse==1)then
             if(lmark(ip1)==1 .and. lmark(ip2)==-1)then      
                ichk=1_ip
             else if(lmark(ip1)==-1 .and. lmark(ip2)==1)then      
                ichk=-1_ip
             endif
          else
             ichk=1_ip 
          endif
          !
          !     Check that the edge must be collapsed
          !
          if(ichk==1)then
             !
             !     Will collapse ip2 if it is a ridge point
             !
             if(lptype(1,ip2)==ID_RIDGE)then
                !
                !     Check length
                !      
                call length(ip1,ip2,rsize,coor,rlen)

                if(rlen<tollen)then

                   call collapsedg(ip1,ip2,lface,coor,ndim,nnofa,nface,npoin,lelem,ptoel1,&
                        ptoel2,lmark,rnopo,icol,ptosi1,ptosi2,lside,nnosi,nside,rsize,&
                        tolrefr,tolgeo)

                endif
             endif

          else if(ichk==-1)then  
             !
             !     Will collapse ip1 if it is a ridge point
             !
             if(lptype(1,ip1)==ID_RIDGE)then 
                !
                !     Check length
                !      
                call length(ip1,ip2,rsize,coor,rlen)

                if(rlen<tollen)then

                   call collapsedg(ip2,ip1,lface,coor,ndim,nnofa,nface,npoin,lelem,ptoel1,&
                        ptoel2,lmark,rnopo,icol,ptosi1,ptosi2,lside,nnosi,nside,rsize,&
                        tolrefr,tolgeo)

                endif
             endif
          endif

       enddo

    enddo

  end subroutine colridg

  subroutine collapsedg(ip1,ip2,lface,coor,ndim,nnofa,nface,npoin,lelem,ptoel1,&
       ptoel2,lmark,rnopo,icol,ptosi1,ptosi2,lside,nnosi,nside,rsize,tollen,tolgeo)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only     :  memor_msh
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,ip1,ip2
    integer(ip),intent(in)    :: nnosi,nside
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),ptosi1(:),ptosi2(:)
    integer(ip),intent(in)    :: nface,npoin
    integer(ip),intent(inout) :: lelem(nface),icol
    integer(ip),intent(inout) :: lface(nnofa,nface),lmark(npoin),lside(nnosi,nside)
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin),rsize(npoin),tollen,tolgeo
    integer(ip)               :: lfacn(nnofa,100),lballn(100),lold(2),nold
    integer(ip)               :: p1,p2,ielem,j,ienew,i,ie,ielem1,ielem2,iview
    integer(ip)               :: ineigh,iplace,ielem0,nballn,icount,ip3,is,iside,ipa,ipb
    real(rp)                  :: rfloc(ndim,100),modul(100),Qgeo,c00,rlen
    real(rp)                  :: rx,ry,rz,tolvol
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    c00=0.0d+00 
    !
    !     This subroutine collapses an edge given by (ip1,ip2)
    !     This is the version used by colridg to collapse ridges 
    !     - eltoel is not updated as the collapsing is performed on the whole patch 
    !       each time, and eltoel is reconstructed each time
    !     - No optimization is performed. It is done on the global level
    !     The sides are updated
    !
    !
    !     Try to collapse ip2 on ip1, replacing ip2 by ip1, deleting ip2
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
    !     Find the new elements = ball of ip2 -(ball of ip1 inter ball of ip2 )
    !
    nballn=0_ip
    nold=0_ip
    do ie=ptoel2(ip2),ptoel2(ip2+1)-1
       ielem=ptoel1(ie)
       !
       !     If an element has already been marked, go home
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
             write(*,*)'Error in collapsedg, nold=2:',ip1,ip2
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
       lelem(ielem)=1_ip
       lface(1,ielem)=lfacn(1,ie)
       lface(2,ielem)=lfacn(2,ie)
       lface(3,ielem)=lfacn(3,ie)
    enddo
    !
    !     Update the sides
    !
    do is=ptosi2(ip2),ptosi2(ip2+1)-1
       iside=ptosi1(is)
       ipa=lside(1,iside)
       ipb=lside(2,iside)
       if((ipa==ip2 .and. ipb==ip1) .or. (ipa==ip1 .and. ipb==ip2))then
          lside(1,iside)=-1
       else
          if(ipa==ip2)then
             lside(1,iside)=ip1
          else
             lside(2,iside)=ip1
          endif
       endif
    enddo
    
    !
    !    Delete elements 
    !
    ielem1=lold(1)
    lelem(ielem1)=-1_ip
    !
    !     Are we collapsing on an isolated ridge
    !
    if(nold==2)then 
       ielem2=lold(2)
       lelem(ielem2)=-1_ip 
    endif
    !
    !     Delete ip2
    !
    lmark(ip2)=-2
    !
    !     Remember collapse successful
    !
    icol=1_ip

  end subroutine collapsedg

  subroutine compact2(nface,npoin,lface,nnofa,coor,ndim,lmark,lcart,rsize,lptype,&
       rnopo,lpsur,nnosi,nside,lpsid,lline,lside,lsmark,ptosi2,ptosi1)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only     :  memor_msh
    implicit none
    integer(ip),intent(inout)           :: nface,npoin
    integer(ip),intent(in)              :: ndim,nnofa,nnosi
    integer(ip),intent(in)              :: ptosi2(npoin+1),ptosi1(*) 
    integer(ip),intent(inout)           :: lface(nnofa,nface),lptype(2,npoin),lcart(npoin)
    integer(ip),intent(inout)           :: lmark(npoin),lpsur(npoin),lpsid(npoin)
    integer(ip),intent(inout)           :: nside
    real(rp),intent(inout)              :: coor(ndim,npoin),rsize(npoin),rnopo(ndim,npoin)
    integer(ip),pointer                 :: lline(:),lside(:,:),lsmark(:)
    integer(ip)                         :: nfold,npold,ifpnt,ipnt,iface,ipoin,iboun  
    integer(ip)                         :: iside,ip1,ip2,kside,inosi,is  
    integer(ip)                         :: ipmin,ipmax,jpmin,jpmax,js,jside
    integer(ip),pointer                 :: lrenum(:)
    integer(4)                          :: istat 

    allocate(lrenum(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENUM','compact2',lrenum) 
    !
    !     Check that the sides have not been repeated
    !     It may happen for very distorded geometries
    !
    do kside=1,nside
       if(lside(1,kside)==-1)cycle
       do inosi=1,nnosi
          ipoin=lside(inosi,kside)
          do is=ptosi2(ipoin),ptosi2(ipoin+1)-2
             iside=ptosi1(is)
             if(lside(1,iside)==-1)cycle
             ipmin=min(lside(1,iside),lside(2,iside))
             ipmax=max(lside(1,iside),lside(2,iside))
             do js=is+1,ptosi2(ipoin+1)-1 
                jside=ptosi1(js)
                if(lside(1,jside)==-1)cycle 
                jpmin=min(lside(1,jside),lside(2,jside))
                jpmax=max(lside(1,jside),lside(2,jside))
                if(ipmin==jpmin .and. ipmax==jpmax)then
                   lside(1,jside)=-1
                endif
             enddo
          enddo
       enddo
    enddo
    !
    !     Compress the side array
    !
    ipnt=0_ip
    do iside=1,nside
       if(lside(1,iside)/=-1)then
          ipnt=ipnt+1
          lside(1,ipnt)=lside(1,iside) 
          lside(2,ipnt)=lside(2,iside) 
          lline(ipnt)=lline(iside)
          lsmark(ipnt)=lsmark(iside)
       endif
    enddo
    nside=ipnt 
    !
    !     Compress the coor array
    !
    npold=npoin
    ipnt=0_ip
    do ipoin=1,npold
       if(lmark(ipoin)/=-2)then 
          ipnt=ipnt+1
          coor(1,ipnt)=coor(1,ipoin)
          coor(2,ipnt)=coor(2,ipoin)
          coor(3,ipnt)=coor(3,ipoin)
          rnopo(1,ipnt)=rnopo(1,ipoin)
          rnopo(2,ipnt)=rnopo(2,ipoin)
          rnopo(3,ipnt)=rnopo(3,ipoin)
          lcart(ipnt)=lcart(ipoin)
          lpsur(ipnt)=lpsur(ipoin)
          lpsid(ipnt)=lpsid(ipoin)
          rsize(ipnt)=rsize(ipoin)
          lptype(1,ipnt)=lptype(1,ipoin)
          lptype(2,ipnt)=lptype(2,ipoin)
          lmark(ipnt)=0_ip
          lrenum(ipoin)=ipnt
       endif
    enddo
    npoin=ipnt
    !
    !     Renumber the lface array
    !
    do iface=1,nface

       if(lrenum(lface(1,iface))==0)then
          write(*,*)'error compact2 1'
          stop
       endif
       if(lrenum(lface(2,iface))==0)then
          write(*,*)'error compact2 2'
          stop
       endif
       if(lrenum(lface(3,iface))==0)then
          write(*,*)'error compact2 3'
          stop
       endif

       lface(1,iface)=lrenum(lface(1,iface))
       lface(2,iface)=lrenum(lface(2,iface))
       lface(3,iface)=lrenum(lface(3,iface))

    enddo
    !
    !     Renumber the lside array
    !
    do iside=1,nside

       ip1=lside(1,iside)
       ip2=lside(2,iside)

       if(lrenum(ip1)==0)then
          write(*,*)'error compact2 1'
          stop
       endif
       if(lrenum(ip2)==0)then
          write(*,*)'error compact2 2'
          stop
       endif

       lside(1,iside)=lrenum(ip1)
       lside(2,iside)=lrenum(ip2)

    enddo

    call memchk(2_ip,istat,memor_msh,'LRENUM','compact2',lrenum)
    deallocate(lrenum,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENUM','compact2',0_ip)

  end subroutine compact2

  subroutine refglo(nface,nnofa,npoin,ndim,lfold,nfold,coorold,npold,lcell,ncell,&
       eltoelold,rnofaold,nsurf,nboup,nside,nnosi,lptypeold,lsurfold,rmax,&
       tolref,tolrefr,tolcol,nsold,rsuni,ptosi1old,ptosi2old,lline,lside,&
       nblay,rblay,tolscal,lsmark,rtol,tolcolr,lstotr2,lblmsh,lface,coor,&
       rnopo,lpsur,lptype,lpsid,lcart,rsize,lsurf,lsold,lstof,sitosiold,&
       ptosi1,ptosi2,nline,llinold,ptoel1old,ptoel2old)
    use def_kintyp, only  :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only  :  memor_msh
    use mod_anisurf, only      :  anisurfl
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,npold,nfold,nnosi,nsold
    integer(ip),intent(in)    :: ncell,nsurf,nboup,nblay,nline
    integer(ip),intent(inout) :: nface,npoin,nside
    integer(ip),intent(in)    :: lfold(nnofa,nfold)
    integer(ip),intent(in)    :: lptypeold(2,npold),lsurfold(nfold)
    integer(ip),intent(in)    :: eltoelold(nnofa,nfold),llinold(nsold)
    integer(ip),intent(in)    :: ptoel2old(npold+1),ptoel1old(*)
    integer(ip),intent(in)    :: sitosiold(:,:),lsold(:,:)
    integer(ip),pointer       :: ptosi1old(:),ptosi2old(:) 
    integer(ip),pointer       :: lside(:,:),lline(:),lsmark(:),lstotr2new(:) 
    integer(ip),pointer       :: lfnew(:,:),lblmshloc(:,:),lelem(:),lfmark(:) 
    integer(ip),pointer       :: lsurfloc(:),eltoel(:,:) 
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),lfloc(:,:),lptri(:) 
    integer(ip),pointer       :: lpofa(:),lmark(:),lcart(:),lface(:,:) 
    integer(ip),pointer       :: lpsur(:),lptype(:,:),lpsid(:),lsurf(:) 
    integer(ip),pointer       :: lstof(:),ledge(:,:),lptoedg2(:),ledg2(:) 
    real(rp),pointer          :: rnopo(:,:),rsize(:),rnofa(:,:),coor(:,:) 
    real(rp),pointer          :: redge(:)
    integer(ip),intent(inout) :: lstotr2(nsurf+1),lblmsh(nsurf)
    type(cell)                :: lcell(ncell)
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold)
    real(rp),intent(in)       :: tolref,tolrefr,tolcol,rsuni,tolcolr
    real(rp),intent(in)       :: rblay(nblay),tolscal,rtol
    real(rp),intent(inout)    :: rmax(nsurf)
    integer(ip)               :: iface,isurf,nfloc,ipnt,nfnew,nblmsh,ierr
    integer(ip)               :: nfloct,iside,ie,j,ipoin,nfnew2,npcusp
    integer(ip)               :: ielem,ipos,idone,npcorner,npoin0,iter,nsloc,ip1,ip2
    integer(4)                :: istat 
    integer(ip), pointer      :: lelemold(:),ldone(:),lsloc(:,:),lcolor(:),lsmloc(:)
    integer(ip), pointer      :: ptosi1(:),ptosi2(:)
    real(rp)                  :: rmaxloc  
    !
    !     This subroutine drives the refinement and optimization of the surface meshing 
    !
    nullify(ptoel1,ptoel2,eltoel,ptosi1,ptosi2)
    !
    !     If there are no sides, dummy allocate
    !
    if(nboup==0)then   
       allocate(ptosi1(1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOSI1','refglo',ptosi1)
       allocate(ptosi2(npoin+1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOSI2','refglo',ptosi2)
    endif
    !
    !     Initialize error flag
    !
    ierr=0_ip
    !
    !     Allocate local arrays for faces
    !
    allocate(lblmshloc(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LBLMSH','refglo',lblmshloc)
    allocate(lfnew(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFNEW','refglo',lfnew)
    allocate(lfloc(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFLOC','refglo',lfloc)
    allocate(lsurfloc(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSURFLOC','refglo',lsurfloc)
    allocate(lsloc(nnosi,nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LSLOC','refglo',lsloc)
    allocate(lsmloc(nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LSMLOC','refglo',lsmloc)
    allocate(lptri(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTRI','refglo',lptri)
    allocate(lelemold(max(nfold,nsold)),stat=istat)           ! Take the max of faces and sides as used by gthostsid also
    call memchk(zero,istat,memor_msh,'LELEMOLD','refglo',lelemold) 
    allocate(lelem(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEM','refglo',lelem) 
    allocate(rnofa(ndim,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOFA','refglo',rnofa) 
    allocate(lfmark(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFMARK','refglo',lfmark) 
    allocate(ldone(nsurf),stat=istat)
    call memchk(zero,istat,memor_msh,'LDONE','refglo',ldone) 
    allocate(lstotr2new(nsurf+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTOTR2NEW','refglo',lstotr2new) 
    allocate(ledge(3,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDGE','reglo',ledge) 
    allocate(redge(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'REDGE','reglo',redge) 
    allocate(ledg2(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDG2','reglo',ledg2) 
    !
    !     Initialize lstotr2new
    !
    lstotr2new(1)=1_ip
    !
    !     Allocate for points
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','refglo',lmark)
    allocate(lpofa(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOFA','refglo',lpofa)
    allocate(lcolor(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LCOLOR','refglo',lcolor)
    allocate(lptoedg2(npoin+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTOEDG2','refglo',lptoedg2)
    !
    !     Main loop on the surfaces:
    !           - first split them with refsplit until reasonable size distribution 
    !           - then smooth the ridges 
    !           - then call advancing front with refsmo 
    !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Loop on refinement by splitting
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    iter=1_ip

    do 
       write(*,*)'Splitting Iteration',iter,'nface=',nface
       !
       !     Initialize new face number for this round
       ! 
       nfnew=0_ip
       !
       !     Initialize "something has been done" flag
       !      
       idone=0_ip
       !
       !     Get the faces surrounding the points for the whole triangulation
       !
       call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
       !
       !     Refine the ridges first
       !
       if(nboup/=0)then

          nfloc=0_ip
          do iside=1,nside
             do j=1,2
                ipoin=lside(j,iside)
                do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
                   ielem=ptoel1(ie) 
                   if(lelem(ielem)==0)then
                      nfloc=nfloc+1
                      lelem(ielem)=nfloc
                      lfloc(1,nfloc)=lface(1,ielem)
                      lfloc(2,nfloc)=lface(2,ielem)
                      lfloc(3,nfloc)=lface(3,ielem)
                      lsurfloc(nfloc)=lsurf(ielem)
                   endif
                enddo
             enddo
          enddo

          nfloct=nfloc

          call refridg(ndim,nnofa,nfold,npold,coorold,rnofaold,nfloc,&
               npoin,tolrefr,lcell,ncell,nside,lelemold,sitosiold,nsold,lsold,nnosi,lstof,rsuni,&
               ptosi1old,ptosi2old,lline,lside,lpsid,lfloc,eltoel,lcart,lsurfloc,rnopo,&
               rsize,ptoel1,ptoel2,lptype,coor,lpsur,lptri,lmark,lpofa,rtol,lsmark,llinold)
          !
          !     Copy lfloc to fill the wholes
          !
          do iface=1,nface
             ielem=lelem(iface)
             if(ielem/=0)then
                lface(1,iface)=lfloc(1,ielem)
                lface(2,iface)=lfloc(2,ielem)
                lface(3,iface)=lfloc(3,ielem)
                lsurf(iface)=lsurfloc(ielem)
                lelem(iface)=0_ip
             endif
          enddo
          !
          !     Reallocate lface (The mesh can only be bigger)
          !  
          ipos=nface
          nface=nface+(nfloc-nfloct)
          call memrea(nface,memor_msh,'LFACE','refglo',lface) 
          call memrea(nface,memor_msh,'LSURF','refglo',lsurf) 
          call memrea(nface,memor_msh,'LELEM','refglo',lelem)
          !
          !     Now copy the new faces
          !
          do iface=nfloct+1,nfloc
             ipos=ipos+1
             lface(1,ipos)=lfloc(1,iface)
             lface(2,ipos)=lfloc(2,iface)
             lface(3,ipos)=lfloc(3,iface)
             lsurf(ipos)=lsurfloc(iface)
          enddo
          !
          !     Renumber the faces with respect to surfaces
          !
          call renuface(nface,nnofa,nsurf,lface,lsurf,lstotr2)
          !
          !     Reallocate local arrays
          !
          call memrea(nface,memor_msh,'LFLOC','refglo',lfloc) 
          call memrea(nface,memor_msh,'LSURFLOC','refglo',lsurfloc) 
          call memrea(nside,memor_msh,'LSLOC','refglo',lsloc) 

       endif
       !
       !     Resize lcolor
       !         
       call memrea(npoin,memor_msh,'LCOLOR','refglo',lcolor) 
       !call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)
       !
       !     Remember the number of points before refining surfaces
       !
       npcorner=npoin+1
       if(nboup/=0)then
          !
          !     Get the sides surrounding points
          !
          call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)

       endif
       !
       !     Resize lptoedg2
       !
       call memrea(npoin+1_ip,memor_msh,'LPTOEDG2','refglo',lptoedg2) 
       !
       !    Loop on the surface patches
       !
       do isurf=1,nsurf
          !
          !     Do we have to do something ?
          ! 
          if(ldone(isurf)/=1)then
             !
             !     Set "something has been done" flag for all the surfaces  
             !
             idone=1_ip
             !
             !     Copy global to local faces
             !
             nfloc=0_ip
             do iface=lstotr2(isurf),lstotr2(isurf+1)-1
                nfloc=nfloc+1_ip
                lfloc(1,nfloc)=lface(1,iface) 
                lfloc(2,nfloc)=lface(2,iface) 
                lfloc(3,nfloc)=lface(3,iface)
                lptri(lfloc(1,nfloc))=nfloc
                lptri(lfloc(2,nfloc))=nfloc
                lptri(lfloc(3,nfloc))=nfloc
                lcolor(lfloc(1,nfloc))=isurf
                lcolor(lfloc(2,nfloc))=isurf
                lcolor(lfloc(3,nfloc))=isurf
             enddo
             !
             !     Get local sides
             !
             nsloc=0_ip
             do iside=1,nside
                ip1=lside(1,iside)
                ip2=lside(2,iside)
                if(lcolor(ip1)==isurf .and. lcolor(ip2)==isurf)then
                   nsloc=nsloc+1
                   lsloc(1,nsloc)=lside(1,iside) 
                   lsloc(2,nsloc)=lside(2,iside) 
                endif
             enddo
             !
             !     Correct geometry for the boundary points to be compatible 
             !     with the current patch
             !  
             call corgeo(nfloc,nnofa,lfloc,npoin,eltoelold,lsurfold,isurf,nfold,lsloc,&
                  nnosi,nsloc,rnopo,rnofaold,ndim,lpsur,lptype,npold,lfold,lsurf,lelemold,&
                  lsold,nsold,lstof,lpsid)
             !
             !     DBG  
             !
             if(ldone(isurf)/=0)then
                write(*,*)'perdu'
                stop
             endif
             !
             !     Do we have to split this patch 
             ! 
             if(ldone(isurf)==0)then
                !
                !     Remember point number
                !
                npoin0=npoin
                !
                !     Split the edges
                !
                call refsplit(ndim,nnofa,nfold,npold,eltoelold,lfold,coorold,rnofaold,nfloc,&
                     npoin,tolref,lcell,ncell,lelemold,isurf,rmaxloc,npcorner,tolcol,rsuni,&
                     lsurfold,lfloc,eltoel,lcart,rnopo,rsize,ptoel1,ptoel2,lptype,&
                     coor,lpsur,lptri,lmark,lpofa,rtol,lpsid,lelem,lptoedg2,ledge,&
                     redge,ledg2,ptosi1,ptosi2,nboup,lside,nnosi,nside,ptoel2old,ptoel1old)
                !
                !     Remember maximal edge length 
                !
                rmax(isurf)=rmaxloc 
                !
                !     Did we do something or can we go to refsmo?
                !       
                if(npoin0==npoin)then
                   ldone(isurf)=1_ip
                endif

             endif
             !
             !     Reallocate lfnew  (The mesh can only be bigger)
             ! 
             ipos=nfnew 
             nfnew=nfnew+nfloc
             call memrea(nfnew,memor_msh,'LFNEW','refglo',lfnew) 
             call memrea(nfnew,memor_msh,'LSURF','refglo',lsurf) 
             call memrea(nfnew,memor_msh,'LELEM','refglo',lelem) 
             call memrea(npoin,memor_msh,'LCOLOR','refglo',lcolor) 
             ! 
             !     Copy local to global   
             !
             do iface=1,nfloc
                ipos=ipos+1
                lfnew(1,ipos)=lfloc(1,iface) 
                lfnew(2,ipos)=lfloc(2,iface) 
                lfnew(3,ipos)=lfloc(3,iface)
                lsurf(ipos)=isurf
                lelem(ipos)=0_ip
             enddo
             !
             !     Update lstotr2new
             !
             lstotr2new(isurf+1)=nfnew+1_ip 
             !
             !     Reallocate local arrays
             !
             call memrea(nfnew,memor_msh,'LFLOC','refglo',lfloc) 

          else
             !
             !     Reallocate lfnew  (The mesh can only be bigger)
             ! 
             ipos=nfnew 
             nfnew=nfnew+lstotr2(isurf+1)-lstotr2(isurf)
             call memrea(nfnew,memor_msh,'LFNEW','refglo',lfnew) 
             call memrea(nfnew,memor_msh,'LSURF','refglo',lsurf) 
             call memrea(nfnew,memor_msh,'LELEM','refglo',lelem) 
             call memrea(npoin,memor_msh,'LCOLOR','refglo',lcolor) 
             ! 
             !     Copy local to global   
             !
             do iface=lstotr2(isurf),lstotr2(isurf+1)-1
                ipos=ipos+1
                lfnew(1,ipos)=lface(1,iface) 
                lfnew(2,ipos)=lface(2,iface) 
                lfnew(3,ipos)=lface(3,iface)
                lsurf(ipos)=isurf
                lelem(ipos)=0_ip
             enddo
             !
             !     Update lstotr2new
             !
             lstotr2new(isurf+1)=nfnew+1_ip 
             !
             !     Reallocate local arrays
             !
             call memrea(nfnew,memor_msh,'LFLOC','refglo',lfloc) 

          endif

       enddo
       !
       !     Copy lfnew in lface
       !
       nface=nfnew
       call memrea(nface,memor_msh,'LFACE','refglo',lface) 
       do iface=1,nfnew
          lface(1,iface)=lfnew(1,iface)
          lface(2,iface)=lfnew(2,iface)
          lface(3,iface)=lfnew(3,iface)
       enddo
       !
       !     Copy lstotr2new
       !
       do isurf=2,nsurf+1        
          lstotr2(isurf)=lstotr2new(isurf)
       enddo
       !call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)
       !
       !     Do we still have something to do?
       !
       if(idone==0)exit

       write(*,*)'Dumping the whole surface in outerrort.msh for iter',iter
       !call outerror4(lface,nnofa,nface,coor,npoin,ndim)
       call outerrort(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,&
            lline,nline)
       !
       !    Resize lptoedg2
       !   
       call memrea(npoin,memor_msh,'LPTOEDG2','refglo',lptoedg2) 

       iter=iter+1

    enddo
    !
    !     Remember the point number at this point for possible cusps
    !
    npcusp=npoin
    !
    !     Do we have boundary points?
    !
    if(nboup/=0)then
       !
       !     Get the sides surrounding points
       !
       call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)
       !
       !     Smooth the ridges
       !
       call optridg(npoin,lptype,coor,ndim,nside,lside,nnosi,coorold,&
            nfold,npold,rnofaold,lcell,ncell,lelemold,&
            sitosiold,lsold,nsold,lstof,ptosi2old,ptosi1old,rsuni,&
            rsize,lpsid,rnopo,lcart,lpsur,rtol,ptosi1,ptosi2,ierr,lline,llinold)
       if(ierr==1)then 
          call outerrort(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,&
               lline,nline)
          stop
       endif


    endif

    write(*,*)'Dumping the whole surface in outerrort.msh after side smoothing'
    !call outerror4(lface,nnofa,nface,coor,npoin,ndim)
    call outerrort(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,&
         lline,nline)
    !
    !     Clean up lcolor
    !
    do ipoin=1,npoin
       lcolor(ipoin)=0_ip
    enddo
    !
    !     Resize lsmloc
    !         
    call memrea(nside,memor_msh,'LSMLOC','refglo',lsmloc) 
    !
    !     Those arrays are not usefull anymore
    !
    call memchk(2_ip,istat,memor_msh,'REDGE','refglo',redge)
    deallocate(redge,stat=istat)
    if(istat/=0) call memerr(2_ip,'REDGE','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDG2','refglo',ledg2)
    deallocate(ledg2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDG2','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTOEDG2','refglo',lptoedg2)
    deallocate(lptoedg2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTOEDG2','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGE','refglo',ledge)
    deallocate(ledge,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGE','refglo',0_ip)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Loop on refinement by advancing front
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    write(*,*)'Remeshing with advancing front, nface=',nface
    !
    !     Clean up lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Initialize new face number for this round
    ! 
    nfnew=0_ip
    !
    !     Initialize "something has been done" flag
    !      
    idone=0_ip
    !
    !    Loop on the surface patches
    !
    do isurf=1,nsurf
       !
       !     Set boundary layer mesh
       ! 
       nblmsh=0_ip 
       !
       !     Do we have to do something ?
       ! 
       if(ldone(isurf)/=2)then
          !
          !     Set "something has been done" flag   
          !
          idone=1_ip
          !
          !     Copy global to local faces
          !
          nfloc=0_ip
          do iface=lstotr2(isurf),lstotr2(isurf+1)-1
             nfloc=nfloc+1_ip
             lfloc(1,nfloc)=lface(1,iface) 
             lfloc(2,nfloc)=lface(2,iface) 
             lfloc(3,nfloc)=lface(3,iface)
             lptri(lfloc(1,nfloc))=nfloc
             lptri(lfloc(2,nfloc))=nfloc
             lptri(lfloc(3,nfloc))=nfloc
             lcolor(lfloc(1,nfloc))=isurf
             lcolor(lfloc(2,nfloc))=isurf
             lcolor(lfloc(3,nfloc))=isurf
          enddo
          !
          !     Get local sides
          !
          nsloc=0_ip
          do iside=1,nside
             ip1=lside(1,iside)
             ip2=lside(2,iside)
             if(lcolor(ip1)==isurf .and. lcolor(ip2)==isurf)then
                nsloc=nsloc+1
                lsloc(1,nsloc)=lside(1,iside) 
                lsloc(2,nsloc)=lside(2,iside) 
                lsmloc(nsloc)=lsmark(iside)
             endif
          enddo
          !
          !     Correct geometry for the boundary points to be compatible 
          !     with the current patch
          !  
          call corgeo(nfloc,nnofa,lfloc,npoin,eltoelold,lsurfold,isurf,nfold,lsloc,&
               nnosi,nsloc,rnopo,rnofaold,ndim,lpsur,lptype,npold,lfold,lsurf,lelemold,&
               lsold,nsold,lstof,lpsid)
          !
          !     DBG  
          !
          if(ldone(isurf)/=1)then
             write(*,*)'perdu 2 '
             stop
          endif

          if(ldone(isurf)==1)then
             !
             !     Do we want to generate a boundary layer mesh?
             ! 
             if(nblay/=0)then
                !
                !     Generate anisotropic mesh from lfloc
                !     Write it in lblmshloc and modify lfloc 
                !            
                call anisurfl(lfloc,nnofa,nfloc,ndim,nblay,coor,npoin,rblay,&
                     lsloc,nsloc,lsmloc,tolscal,nnosi,coorold,eltoelold,lfold,nfold,&
                     npold,rnofaold,lelemold,lsurfold,isurf,rsize,rsuni,rnopo,&
                     lcell,ncell,lcart,lpsur,lptri,lblmshloc,nblmsh,lptype,&
                     tolcol,tolcolr,tolref,tolrefr,lline,nboup)
                !
                !     Copy the new anisotropic mesh in lfnew
                ! 
                ipnt=nfnew 
                nfnew=nblmsh+nfnew
                call memrea(nfnew,memor_msh,'LFACE','refglo',lfnew) 
                call memrea(nfnew,memor_msh,'LSURF','refglo',lsurf) 
                call memrea(nfnew,memor_msh,'LELEM','refglo',lelem) 
                call resizep(npoin,coor,rnopo,lcart,lpsur,rsize,lptri,lmark,lpofa,lptype,lpsid)
                call memrea(npoin,memor_msh,'LCOLOR','refglo',lcolor) 
                ! 
                !     Copy anisotropic mesh to new mesh   
                !
                do iface=1,nblmsh
                   ipnt=ipnt+1
                   lfnew(1,ipnt)=lblmshloc(1,iface) 
                   lfnew(2,ipnt)=lblmshloc(2,iface) 
                   lfnew(3,ipnt)=lblmshloc(3,iface)
                   lsurf(ipnt)=isurf
                   lelem(ipnt)=0_ip
                enddo

             endif
             !
             !     Refine inside each patch
             !
             write(*,*)'Remeshing surface:',isurf
             call refsmo(nfloc,nnofa,npoin,ndim,lfold,nfold,coorold,npold,lcell,ncell,&
                  eltoelold,rnofaold,isurf,lelemold,lsurfold,rsuni,rnopo,lcart,lpsur,&
                  rsize,lfloc,eltoel,ptoel1,ptoel2,lptri,lfmark,rnofa,&
                  lmark,lelem,lpofa,lptype,coor,rtol,lpsid,ptosi1,ptosi2,npcusp,lside,&
                  nside,nnosi)
             !
             !     Mark the surface as done
             !
             ldone(isurf)=2_ip 
             !
          endif
          !
          !     Reallocate lfnew  (The mesh can only be bigger)
          !
          ipnt=nfnew
          nfnew=nfnew+nfloc
          call memrea(nfnew,memor_msh,'LFNEW','refglo',lfnew) 
          call memrea(nfnew,memor_msh,'LSURF','refglo',lsurf) 
          call memrea(nfnew,memor_msh,'LELEM','refglo',lelem) 
          call memrea(npoin,memor_msh,'LCOLOR','refglo',lcolor) 
          ! 
          !     Copy local to global   
          !
          do iface=1,nfloc
             ipnt=ipnt+1
             lfnew(1,ipnt)=lfloc(1,iface) 
             lfnew(2,ipnt)=lfloc(2,iface) 
             lfnew(3,ipnt)=lfloc(3,iface)
             lsurf(ipnt)=isurf
             lelem(ipnt)=0_ip
          enddo
          !
          !     Update lstotr2new
          !
          lstotr2new(isurf+1)=nfnew+1_ip 
          !
          !     Remember boundary layer
          !
          lblmsh(isurf)=nblmsh
          !
          !     Reallocate local arrays
          !
          call memrea(nfnew,memor_msh,'LFLOC','refglo',lfloc) 

       else
          !
          !     Reallocate lfnew  (The mesh can only be bigger)
          ! 
          ipos=nfnew 
          nfnew=nfnew+lstotr2(isurf+1)-lstotr2(isurf)
          call memrea(nfnew,memor_msh,'LFNEW','refglo',lfnew) 
          call memrea(nfnew,memor_msh,'LSURF','refglo',lsurf) 
          call memrea(nfnew,memor_msh,'LELEM','refglo',lelem) 
          call memrea(npoin,memor_msh,'LCOLOR','refglo',lcolor) 
          ! 
          !     Copy local to global   
          !
          do iface=lstotr2(isurf),lstotr2(isurf+1)-1
             ipos=ipos+1
             lfnew(1,ipos)=lface(1,iface) 
             lfnew(2,ipos)=lface(2,iface) 
             lfnew(3,ipos)=lface(3,iface)
             lsurf(ipos)=isurf
             lelem(ipos)=0_ip
          enddo
          !
          !     Update lstotr2new
          !
          lstotr2new(isurf+1)=nfnew+1_ip 
          !
          !     Remember boundary layer
          !
          lblmsh(isurf)=nblmsh
          !
          !     Reallocate local arrays
          !
          call memrea(nfnew,memor_msh,'LFLOC','refglo',lfloc) 

       endif

    enddo

    !call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)
    !
    !     Copy lfnew in lface
    !
    nface=nfnew
    call memrea(nface,memor_msh,'LFACE','refglo',lface) 
    do iface=1,nfnew
       lface(1,iface)=lfnew(1,iface)
       lface(2,iface)=lfnew(2,iface)
       lface(3,iface)=lfnew(3,iface)
    enddo
    !
    !     Copy lstotr2new
    !
    do isurf=2,nsurf+1        
       lstotr2(isurf)=lstotr2new(isurf)
    enddo

    write(*,*)'Dumping surface mesh after global refining in outface.msh'
    call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)

    call memchk(2_ip,istat,memor_msh,'LCOLOR','refglo',lcolor)
    deallocate(lcolor,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCOLOR','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPOFA','refglo',lpofa)
    deallocate(lpofa,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPOFA','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','refglo',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTOTR2NEW','refglo',lstotr2new)
    deallocate(lstotr2new,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTOTR2NEW','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LDONE','refglo',ldone)
    deallocate(ldone,stat=istat)
    if(istat/=0) call memerr(2_ip,'LDONE','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LELEM','refglo',lelem)
    deallocate(lelem,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEM','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LELEMOLD','refglo',lelemold)
    deallocate(lelemold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEMOLD','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTRI','refglo',lptri)
    deallocate(lptri,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTRI','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSMLOC','refglo',lsmloc)
    deallocate(lsmloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSMLOC','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSLOC','refglo',lsloc)
    deallocate(lsloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSLOC','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSURFLOC','refglo',lsurfloc)
    deallocate(lsurfloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSURFLOC','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFLOC','refglo',lfloc)
    deallocate(lfloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFLOC','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFNEW','refglo',lfnew)
    deallocate(lfnew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFNEW','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LBLMSH','refglo',lblmshloc)
    deallocate(lblmshloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBLMSH','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOFA','refglo',rnofa)
    deallocate(rnofa,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOFA','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','refglo',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','refglo',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','refglo',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFMARK','refglo',lfmark)
    deallocate(lfmark,stat=istat)

  end subroutine refglo

  subroutine refridg(ndim,nnofa,nfold,npold,coorold,rnofaold,nface,npoin,&
       tollen,lcell,ncell,nside,lelemold,sitosiold,nsold,lsold,nnosi,lstof,rsuni,ptosi1old,&
       ptosi2old,lline,lside,lpsid,lface,eltoel,lcart,lsurf,rnopo,rsize,ptoel1,ptoel2,&
       lptype,coor,lpsur,lptri,lmark,lpofa,rtol,lsmark,llinold)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only        :  memor_msh
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,nfold,npold,ncell,nsold,nnosi
    integer(ip),intent(in)    :: sitosiold(nnosi,nsold),lsold(nnosi,nsold),lstof(nsold)
    integer(ip),intent(in)    :: llinold(nsold)
    type(cell),intent(in)     :: lcell(ncell)
    integer(ip),pointer       :: ptosi1old(:),ptosi2old(:),lside(:,:),lline(:)
    integer(ip),pointer       :: lpsid(:),lface(:,:),eltoel(:,:),lcart(:),lsurf(:)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),lptype(:,:),lpsur(:) 
    integer(ip),pointer       :: lptri(:),lmark(:),lpofa(:),lsmark(:) 
    real(rp),pointer          :: rnopo(:,:),rsize(:),coor(:,:) 
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold),tollen,rsuni,rtol 
    integer(ip),intent(inout) :: nface,npoin,nside,lelemold(nfold)
    integer(ip)               :: ipo,ipoin,ielem,j,ip1,ip2,ie,nsidnew,ierr
    integer(ip)               :: ipa,ipb,iside,iline,nsid0,ichk,ismark
    real(rp)                  :: rlen,tolsplit
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub splits the ridges 
    !     Split without considering the length of the newly created edges
    ! 
    ichk=0_ip
    tolsplit=1_ip   !not used
    ierr=0_ip
    !call outerror4(lface,nnofa,nface,coor,npoin,ndim)
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface, nface , npoin, nnofa ,ptoel1,ptoel2)
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)

    nsid0=nside
    do iside=1,nsid0
       ipa=lside(1,iside)
       ipb=lside(2,iside)
       !write(*,*)iside,ipa,ipb
       iline=lline(iside)
       ismark=lsmark(iside)
       !
       !     Loop on surrounding elements finding the elements with two boundary points
       !
       do ie=ptoel2(ipa),ptoel2(ipa+1)-1
          ielem=ptoel1(ie)
          do j=1,nnofa
             ip1=lface(ltab(1,j),ielem)
             ip2=lface(ltab(2,j),ielem)
             if((ip1==ipa .and. ip2==ipb) .or.  (ip1==ipb .and. ip2==ipa))then
                !
                !     Check length
                ! 
                call length(ip1,ip2,rsize,coor,rlen)
                !
                !     Split if too big
                !
                if(rlen>tollen)then

                   call splitsid(ielem,j,npoin,ndim,coorold,nface,&
                        nfold,npold,rnofaold,lcell,ncell,lelemold,tolsplit,ichk,sitosiold,&
                        lsold,nsold,nnosi,lstof,ptosi2old,ptosi1old,rsuni,ptoel1,ptoel2,lpsid, &
                        lface,eltoel,lcart,rnopo,lpsur,rsize,&
                        lptri,lmark,lpofa,lsurf,coor,rtol,lptype,ierr,iline,llinold)
                   if(ierr==1)then
                      write(*,*)'Error in refridg, dumping outerror4'
                      call outerror4(lface,nface,nnofa,coor,npoin,ndim)
                      stop
                   endif
                   !
                   !     Update lside
                   !
                   nsidnew=nside+1
                   call memrea(nsidnew,memor_msh,'LSIDE','refridg',lside)
                   call memrea(nsidnew,memor_msh,'LLINE','refridg',lline)
                   call memrea(nsidnew,memor_msh,'LSMARK','refridg',lsmark)
                   lside(1,iside)=ip1 
                   lside(2,iside)=npoin
                   lline(iside)=iline
                   lsmark(iside)=ismark
                   nside=nside+1 
                   lside(1,nside)=npoin
                   lside(2,nside)=ip2
                   lline(nside)=iline
                   lsmark(nside)=ismark
                   !
                   !     Update point type
                   !
                   lptype(1,npoin)=ID_RIDGE
                   lptype(2,npoin)=iline
                   !write(*,*)iside,nside
                endif
             endif
          enddo
       enddo
    enddo
    !
    !     Should smooth out the mesh here
    !

  end subroutine refridg

  subroutine refsplit(ndim,nnofa,nfold,npold,eltoelold,lfold,coorold,rnofaold,nface, &
       npoin,tollen,lcell,ncell,lelemold,isurf,rmax,npcorner,tolcol,rsuni,lsurfold, &
       lface,eltoel,lcart,rnopo,rsize,ptoel1,ptoel2,lptype, &
       coor,lpsur,lptri,lmark,lpofa,rtol,lpsid,lelem,lptoedg2,&
       ledge,redge,ledg2,ptosi1,ptosi2,nboup,lside,nnosi,nside,&
       ptoel2old,ptoel1old)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only        :  memor_msh
    use mod_sort, only          :  sort3 
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,nfold,npold,ncell,npcorner,isurf,nboup
    integer(ip),intent(in)    :: nnosi,nside
    integer(ip),intent(in)    :: eltoelold(nnofa,nfold),lfold(nnofa,nfold)
    integer(ip),intent(in)    :: lsurfold(nfold),lside(nnosi,nside)
    integer(ip),intent(in)    :: ptoel2old(npold+1),ptoel1old(*)
    type(cell),intent(in)     :: lcell(ncell)
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold) 
    real(rp),intent(in)       :: tollen,tolcol,rsuni,rtol 
    integer(ip),intent(inout) :: nface,npoin,lelemold(nfold),lptoedg2(npoin+1)
    integer(ip),pointer       :: lface(:,:),eltoel(:,:),lcart(:),lelem(:)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),lptype(:,:),lpsid(:),ledg2(:)
    integer(ip),pointer       :: lpsur(:),lptri(:),lmark(:),lpofa(:),ledge(:,:)
    integer(ip),pointer       :: ptosi1(:),ptosi2(:)
    real(rp),pointer          :: coor(:,:),rnopo(:,:),rsize(:),redge(:)
    real(rp), intent(inout)   :: rmax
    integer(ip)               :: ipo,ipoin,ielem,j,ip1,ip2,ie,ielem2,nsidnew,ierr
    integer(ip)               :: ipa,ipb,nfold1,iter,ichk,ledg(3),jnofa,jpmin
    integer(ip)               :: iedg,jedg,nedg,isto,iface,inofa,jpoin,iedge
    integer(ip)               :: jp1,jp2,jpmax,nedge,jedge,npcornerm1,nedg0
    integer(ip)               :: jsto,ipmin,ipmax,ifound,iside
    real(rp)                  :: rlen,c05,redg(3)
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub splits the edges taking into account the length of the newly 
    !     created edges
    !
    ichk=1_ip
    rmax=0.0d+00
    c05=0.5d+00
    ierr=0_ip 
    npcornerm1=npcorner-1
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Clean up lmark
    ! 
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Get the edges
    !
    lptoedg2(1)=1_ip
    nedge=0_ip
    do ipoin=1,npcornerm1
       do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(isto)    
          do inofa=1,nnofa
             jpoin=lface(inofa,iface)   
             if(jpoin>ipoin)then
                if(lmark(jpoin)/=ipoin)then
                   lmark(jpoin)=ipoin
                   !
                   !     Is this edge a ridge?
                   !
                   ifound=0_ip
                   do jsto=ptosi2(ipoin),ptosi2(ipoin+1)-1
                      iside=ptosi1(jsto)
                      ipmin=min(lside(1,iside),lside(2,iside))
                      ipmax=max(lside(1,iside),lside(2,iside))
                      if(ipmin==ipoin .and. ipmax==jpoin)then
                         ifound=1
                         exit
                      endif
                   enddo

                   if(ifound==0)then
                      nedge=nedge+1           
                      call memrea(nedge,memor_msh,'LEDGE','refsplit',ledge)
                      ledge(1,nedge)=ipoin
                      ledge(2,nedge)=jpoin
                      ledge(3,nedge)=iface
                   endif
                endif
             endif
          enddo
       enddo
       lptoedg2(ipoin+1)=nedge+1
    enddo
    !
    !     Reallocate redge and ledg2
    ! 
    call memrea(nedge,memor_msh,'REDGE','refsplit',redge)
    call memrea(nedge,memor_msh,'LEDG2','refsplit',ledg2)
    !
    !     Compute length
    !
    nedg0=0_ip
    do iedge=1,nedge
       ip1=ledge(1,iedge)
       ip2=ledge(2,iedge)
       call length(ip1,ip2,rsize,coor,rlen)
       if(rlen>rmax)rmax=rlen
       !
       !     Is the edge too large?
       !
       if(rlen>tollen)then
          nedg0=nedg0+1
          redge(nedg0)=rlen
          ledg2(nedg0)=iedge
       endif
    enddo
    !
    !     Sort the edges
    !
    call  sort3(ledg2,redge,nedg0)            
    !
    !     Loop on edges
    !    
    do jedge=nedg0,1,-1

       iedge=ledg2(jedge)

       ip1=ledge(1,iedge)
       ip2=ledge(2,iedge)
       iface=ledge(3,iedge)

       !write(*,*)jedge,ledge(1,9574),ledge(2,9574),lface(1,ledge(3,9574)),lface(2,ledge(3,9574)),lface(3,ledge(3,9574))
       do inofa=1,nnofa
          jp1=lface(ltab(1,inofa),iface)
          jp2=lface(ltab(2,inofa),iface)
          jpmax=max(jp1,jp2)
          jpmin=min(jp1,jp2)
          if(jpmin==ip1 .and. jpmax==ip2)then
             exit
          endif
       enddo
       !
       !     Any error?
       ! 
       if(inofa==nnofa+1)then
          write(*,*)'Error in refsplit, edge not found'
          stop
       endif

       ichk=1_ip 
       call split(iface,inofa,nnofa,nface,npoin,ndim,coorold,eltoelold,lfold,nfold,&
            npold,rnofaold,lcell,ncell,lelemold,tolcol,ichk,rsuni,lsurfold,isurf,lface,&
            eltoel,lcart,rnopo,lpsur,rsize,lptri,lmark,lpofa,coor,rtol,lptype,lpsid,&
            lelem,lptoedg2,ledge,nedge,npcorner,ierr,ptoel2old,ptoel1old)
       if(ierr==1)then
          write(*,*)'Error in refsplit'
          call outerror4(lface,nnofa,nface,coor,npoin,ndim)
          stop
       endif
       !
       !     Update point type
       !
       !
       !     Was the point successfully created so that all the new edges
       !     agree with the desired length?
       !
       if(ichk==1)then
          lptype(1,npoin)=ID_SMOOTH
          lptype(2,npoin)=isurf
       endif

    enddo
    !
    !     Iteration number
    !  
    iter=2_ip
    !
    !     Do we have sides ?
    !
    if(nboup>0)then
       !
       !     Get the sides surrounding the points
       !
       call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)
    else 
       call memrea(npoin+1_ip,memor_msh,'PTOSI2','refsplit',ptosi2)
       do ipoin=1,npoin+1
          ptosi2(ipoin)=0_ip
       enddo

    endif
    !
    !    Optimize the patch only by swapping 
    !
    !call swapsurf(lface,nnofa,nface,coor,npoin,ndim,rnopo,eltoel,iter,ptosi1,&
    !     ptosi2,lside,nnosi,nside)

    !write(*,*)'Dumping surface in outerror4.msh after splitting for isurf=',isurf
    !call outerror4(lface,nnofa,nface,coor,npoin,ndim)


  end subroutine refsplit

  subroutine split(ielem1,iview1,nnofa,nface,npoin,ndim,coorold,eltoelold,lfold,nfold,&
       npold,rnofaold,lcell,ncell,lelemold,tolsplit,ichk,rsuni,lsurfold,isurf,lface,&
       eltoel,lcart,rnopo,lpsur,rsize,lptri,lmark,lpofa,coor,rtol,lptype,lpsid,lelem,&
       lptoedg2,ledge,nedge,npcorner,ierr,ptoel2old,ptoel1old)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only     :  memor_msh
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,ielem1,iview1,nfold,npold,ncell
    integer(ip),intent(inout) :: nface,npoin,lelemold(nfold),ichk
    integer(ip),intent(in)    :: npcorner,isurf,nedge
    integer(ip),intent(in)    :: eltoelold(nnofa,nfold),lfold(nnofa,nfold),lsurfold(nfold)
    integer(ip),intent(in)    :: lptoedg2(npoin+1),ptoel2old(npold+1),ptoel1old(*)
    integer(ip),intent(inout) :: ledge(3,nedge),ierr
    type(cell),intent(in)     :: lcell(ncell) 
    integer(ip),pointer       :: lface(:,:),eltoel(:,:),lcart(:),lpsur(:),lptri(:)
    integer(ip),pointer       :: lmark(:),lpofa(:),lptype(:,:),lpsid(:),lelem(:)
    real(rp),pointer          :: rnopo(:,:),rsize(:),coor(:,:)
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold),tolsplit,rsuni,rtol 
    integer(ip)               :: ip1,ip2,ip3,ip4,iview2,ip1o,ip2o,ip3o,ipclos,iguess
    integer(ip)               :: ineigh11,ineigh12,ineigh21,ineigh22,ipnew,icart
    integer(ip)               :: isurf1,isurf2,nfnew,ielem3,ielem4,npnew,ihost,ielem2
    integer(ip)               :: ipmin,ipmax,iedge,ifaold,isurfold,isto,lfloc(3,2)
    real(rp)                  :: d1,d2,d3,dtot,c05,rlen,c00,csca,cscal,rnl3,rnofloc(3,2)
    real(rp)                  :: rx,ry,rz,rx1,ry1,rz1,rnl,rnli,rx2,ry2,rz2,c10,modul
    real(rp)                  :: rnx,rny,rnz,cscat,rstot,rs2,rs3,rsiz
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub splits the edge separing elements ielem1 and ielem2
    !     by inserting a point at the middle of the edge
    !     Be careful: lptype is updated outside
    !                 eltoel is updated
    !
    c00=0.0d+00
    c05=0.5d+00
    c10=1.0d+00
    !if(eltoel(1,ielem1)==ielem2)then
    !   iview1=1_ip
    !else if(eltoel(2,ielem1)==ielem2)then
    !   iview1=2_ip
    !else if(eltoel(3,ielem1)==ielem2)then
    !   iview1=3_ip
    !else
    !   write(*,*)'error in split 1'
    !   stop  
    !endif

    ineigh11=eltoel(ltab(1,iview1),ielem1)
    ineigh12=eltoel(ltab(2,iview1),ielem1)

    ip1=lface(iview1,ielem1)
    ip2=lface(ltab(1,iview1),ielem1)
    ip3=lface(ltab(2,iview1),ielem1)

    ielem2=eltoel(iview1,ielem1)

    if(ielem2==0)then
       write(*,*)'Strange ielem2=0 found' 
       stop
    endif

    if(ielem2/=0)then

       if(eltoel(1,ielem2)==ielem1)then
          iview2=1_ip
       else if(eltoel(2,ielem2)==ielem1)then
          iview2=2_ip
       else if(eltoel(3,ielem2)==ielem1)then
          iview2=3_ip
       else
          write(*,*)'error in split 2'
          stop  
       endif

       ip4=lface(iview2,ielem2)
       ineigh21=eltoel(ltab(1,iview2),ielem2)
       ineigh22=eltoel(ltab(2,iview2),ielem2)

    endif
    !
    !     Allocate data base for ipnew
    !
    npnew=npoin+1
    ipnew=npnew
    call resizep(npnew,coor,rnopo,lcart,lpsur,rsize,lptri,lmark,lpofa,lptype,lpsid)
    npoin=npnew
    !
    !     Update the point database for ipnew
    !
    rstot=rsize(ip2)+rsize(ip3)
    rs2=rsize(ip3)/rstot
    rs3=rsize(ip2)/rstot
    coor(1,ipnew)=rs2*coor(1,ip2)+rs3*coor(1,ip3)
    coor(2,ipnew)=rs2*coor(2,ip2)+rs3*coor(2,ip3)
    coor(3,ipnew)=rs2*coor(3,ip2)+rs3*coor(3,ip3)
    !
    !     Interpolated from a smooth point
    !
    if(lptype(1,ip1)==ID_SMOOTH)then
       iguess=lpsur(ip1)
    else if(lptype(1,ip2)==ID_SMOOTH)then
       iguess=lpsur(ip2)
       !
       !   Can not find a smooth point, try with cusp and ridge point
       !
    else if(lptype(1,ip1)/=ID_CORNER)then 
       iguess=lpsur(ip1)
    else if(lptype(1,ip2)/=ID_CORNER)then 
       iguess=lpsur(ip2)
    else 
       !
       !     Find the face from isurf mostly aligned with both normals 
       !
       lfloc(1,1)=ip1
       lfloc(2,1)=ip2
       lfloc(3,1)=ip3
       call gtfnr2(lfloc,2_ip,nnofa,ndim,coor,npoin,rnofloc(1,1),1_ip,modul)
       lfloc(1,2)=ip4
       lfloc(2,2)=ip3
       lfloc(3,2)=ip2
       call gtfnr2(lfloc,2_ip,nnofa,ndim,coor,npoin,rnofloc(1,2),2_ip,modul)
       rnx=rnofloc(1,1)+rnofloc(1,2)
       rny=rnofloc(2,1)+rnofloc(2,2)
       rnz=rnofloc(3,1)+rnofloc(3,2)
       rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
       rnl=c10/rnl
       rnx=rnx*rnl
       rny=rny*rnl
       rnz=rnz*rnl
       csca=-1.01d+00
       do isto=ptoel2old(ip1),ptoel2old(ip1+1)-1
          ifaold=ptoel1old(isto)
          isurfold=lsurfold(ifaold)
          if(isurfold==isurf)then
             cscat=rnx*rnofaold(1,ifaold)+rny*rnofaold(2,ifaold)+rnz*rnofaold(3,ifaold)
             if(cscat>csca)then
                csca=cscat
                iguess=ifaold
             endif
          endif
       enddo
    endif

    rsiz=c05*(rsize(ip2)+rsize(ip3)) 

    !write(*,*)'ipnew=',ipnew,ip2,ip3
    call gthost(coor(:,ipnew),iguess,ihost,rnofaold,ndim,nfold,npold,lfold,&
         nnofa,coorold,d1,d2,d3,eltoelold,lelemold,lsurfold,isurf,ierr,rsiz)
    if(ierr==1)then
       write(*,*)'Error in split',ip2,ip3 
       return
    endif
    !
    !     Get the projected point
    !
    !dtot=d1+d2+d3
    !d1=d1/dtot  
    !d2=d2/dtot  
    !d3=d3/dtot  
    ip1o=lfold(1,ihost)
    ip2o=lfold(2,ihost)
    ip3o=lfold(3,ihost)

    !write(*,*)'split'
    !write(*,*)d1,d2,d3
    !write(*,*)ihost,ip1o,ip2o,ip3o

    coor(1,ipnew)=d1*coorold(1,ip1o)+d2*coorold(1,ip2o)+d3*coorold(1,ip3o)
    coor(2,ipnew)=d1*coorold(2,ip1o)+d2*coorold(2,ip2o)+d3*coorold(2,ip3o)
    coor(3,ipnew)=d1*coorold(3,ip1o)+d2*coorold(3,ip2o)+d3*coorold(3,ip3o)
    rnopo(1,ipnew)=rnofaold(1,ihost)
    rnopo(2,ipnew)=rnofaold(2,ihost)
    rnopo(3,ipnew)=rnofaold(3,ihost)
    icart=lcart(ip1)

    !write(*,*)coor(1,ipnew),coor(2,ipnew),coor(3,ipnew)
    !write(*,*)ihost,ip1o,ip2o,ip3o

    call gtelem(ipnew,coor,npoin,ndim,lcell,ncell,icart,rtol)
    lcart(ipnew)=icart
    call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipnew,ndim,coor,rsuni)
    lpsur(ipnew)=ihost 
    lptri(ipnew)=ielem1
    lmark(ipnew)=0_ip
    lpofa(ipnew)=0_ip
    !
    !     Do we want to check the length of the newly created edges?
    !
    if(ichk==1)then
       call length(ip1,ipnew,rsize,coor,rlen)
       if(rlen<tolsplit)then
          npoin=npoin-1
          ichk=0_ip
          return
       endif

       if(ielem2/=0)then
          call length(ip4,ipnew,rsize,coor,rlen)
          if(rlen<tolsplit)then
             npoin=npoin-1
             ichk=0_ip
             return
          endif
       endif
       !
       !     Check for close points
       !
       call gtclos3(ielem1,ipnew,lface,nface,nnofa,coor,ndim,npoin,rsize,&
            lelem,eltoel,ipclos,tolsplit)
       if(ipclos/=0)then 
          npoin=npoin-1
          ichk=0_ip
          return
       endif
       !
       !     Check for close sides
       !
       if(ineigh11==0)then
          !
          !   Get distance to side
          !
          rx1=coor(1,ip3)-coor(1,ip1)
          ry1=coor(2,ip3)-coor(2,ip1)
          rz1=coor(3,ip3)-coor(3,ip1)
          rnl=sqrt(rx1*rx1+ry1*ry1+rz1*rz1) 
          rnli=c10/rnl
          rx1=rx1*rnli
          ry1=ry1*rnli
          rz1=rz1*rnli
          rx2=coor(1,ipnew)-coor(1,ip1)
          ry2=coor(2,ipnew)-coor(2,ip1)
          rz2=coor(3,ipnew)-coor(3,ip1)
          csca=rx1*rx2+ry1*ry2+rz1*rz2
          cscal=csca*rnli
          if(cscal<c00)then
             rx=coor(1,ipnew)-coor(1,ip1)
             ry=coor(2,ipnew)-coor(2,ip1)
             rz=coor(3,ipnew)-coor(3,ip1)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          else if(cscal>c10)then 
             rx=coor(1,ipnew)-coor(1,ip3)
             ry=coor(2,ipnew)-coor(2,ip3)
             rz=coor(3,ipnew)-coor(3,ip3)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          else 
             rx=coor(1,ip1)+csca*rx1-coor(1,ipnew)
             ry=coor(2,ip1)+csca*ry1-coor(2,ipnew)
             rz=coor(3,ip1)+csca*rz1-coor(3,ipnew)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          endif
          if(rnl3<tolsplit*rsize(ipnew))then
             npoin=npoin-1
             ichk=0_ip
             return
          endif
       endif

       if(ineigh12==0)then
          rx1=coor(1,ip2)-coor(1,ip1)
          ry1=coor(2,ip2)-coor(2,ip1)
          rz1=coor(3,ip2)-coor(3,ip1)
          rnl=sqrt(rx1*rx1+ry1*ry1+rz1*rz1) 
          rnli=c10/rnl
          rx1=rx1*rnli
          ry1=ry1*rnli
          rz1=rz1*rnli
          rx2=coor(1,ipnew)-coor(1,ip1)
          ry2=coor(2,ipnew)-coor(2,ip1)
          rz2=coor(3,ipnew)-coor(3,ip1)
          csca=rx1*rx2+ry1*ry2+rz1*rz2
          cscal=csca*rnli
          if(cscal<c00)then
             rx=coor(1,ipnew)-coor(1,ip1)
             ry=coor(2,ipnew)-coor(2,ip1)
             rz=coor(3,ipnew)-coor(3,ip1)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          else if(cscal>c10)then 
             rx=coor(1,ipnew)-coor(1,ip2)
             ry=coor(2,ipnew)-coor(2,ip2)
             rz=coor(3,ipnew)-coor(3,ip2)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          else 
             rx=coor(1,ip1)+csca*rx1-coor(1,ipnew)
             ry=coor(2,ip1)+csca*ry1-coor(2,ipnew)
             rz=coor(3,ip1)+csca*rz1-coor(3,ipnew)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          endif
          if(rnl3<tolsplit*rsize(ipnew))then
             npoin=npoin-1
             ichk=0_ip
             return
          endif
       endif

       if(ineigh21==0)then
          !
          !   Get distance to side
          !
          rx1=coor(1,ip2)-coor(1,ip4)
          ry1=coor(2,ip2)-coor(2,ip4)
          rz1=coor(3,ip2)-coor(3,ip4)
          rnl=sqrt(rx1*rx1+ry1*ry1+rz1*rz1) 
          rnli=c10/rnl
          rx1=rx1*rnli
          ry1=ry1*rnli
          rz1=rz1*rnli
          rx2=coor(1,ipnew)-coor(1,ip4)
          ry2=coor(2,ipnew)-coor(2,ip4)
          rz2=coor(3,ipnew)-coor(3,ip4)
          csca=rx1*rx2+ry1*ry2+rz1*rz2
          cscal=csca*rnli
          if(cscal<c00)then
             rx=coor(1,ipnew)-coor(1,ip4)
             ry=coor(2,ipnew)-coor(2,ip4)
             rz=coor(3,ipnew)-coor(3,ip4)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          else if(cscal>c10)then 
             rx=coor(1,ipnew)-coor(1,ip2)
             ry=coor(2,ipnew)-coor(2,ip2)
             rz=coor(3,ipnew)-coor(3,ip2)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          else 
             rx=coor(1,ip4)+csca*rx1-coor(1,ipnew)
             ry=coor(2,ip4)+csca*ry1-coor(2,ipnew)
             rz=coor(3,ip4)+csca*rz1-coor(3,ipnew)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          endif
          if(rnl3<tolsplit*rsize(ipnew))then
             npoin=npoin-1
             ichk=0_ip
             return
          endif
       endif

       if(ineigh22==0)then
          rx1=coor(1,ip3)-coor(1,ip4)
          ry1=coor(2,ip3)-coor(2,ip4)
          rz1=coor(3,ip3)-coor(3,ip4)
          rnl=sqrt(rx1*rx1+ry1*ry1+rz1*rz1) 
          rnli=c10/rnl
          rx1=rx1*rnli
          ry1=ry1*rnli
          rz1=rz1*rnli
          rx2=coor(1,ipnew)-coor(1,ip4)
          ry2=coor(2,ipnew)-coor(2,ip4)
          rz2=coor(3,ipnew)-coor(3,ip4)
          csca=rx1*rx2+ry1*ry2+rz1*rz2
          cscal=csca*rnli
          if(cscal<c00)then
             rx=coor(1,ipnew)-coor(1,ip4)
             ry=coor(2,ipnew)-coor(2,ip4)
             rz=coor(3,ipnew)-coor(3,ip4)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          else if(cscal>c10)then 
             rx=coor(1,ipnew)-coor(1,ip3)
             ry=coor(2,ipnew)-coor(2,ip3)
             rz=coor(3,ipnew)-coor(3,ip3)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          else 
             rx=coor(1,ip4)+csca*rx1-coor(1,ipnew)
             ry=coor(2,ip4)+csca*ry1-coor(2,ipnew)
             rz=coor(3,ip4)+csca*rz1-coor(3,ipnew)
             rnl3=sqrt(rx*rx+ry*ry+rz*rz)
          endif
          if(rnl3<tolsplit*rsize(ipnew))then
             npoin=npoin-1
             ichk=0_ip
             return
          endif
       endif

    endif

    if(ielem2/=0)then
       !
       !     Allocate the 2 new elements
       !
       nfnew=nface+2
       call memrea(nfnew,memor_msh,'LFACE','split',lface)
       call memrea(nfnew,memor_msh,'ELTOEL','split',eltoel)
       call memrea(nfnew,memor_msh,'LELEM','split',lelem)
       nface=nface+1
       ielem3=nface
       nface=nface+1
       ielem4=nface
       !
       !     Create new elements
       ! 
       lface(1,ielem1)=ip1  
       lface(2,ielem1)=ip2  
       lface(3,ielem1)=ipnew
       eltoel(1,ielem1)=ielem2
       eltoel(2,ielem1)=ielem3
       eltoel(3,ielem1)=ineigh12
       lelem(ielem1)=0_ip

       lface(1,ielem2)=ip2  
       lface(2,ielem2)=ip4  
       lface(3,ielem2)=ipnew
       eltoel(1,ielem2)=ielem4
       eltoel(2,ielem2)=ielem1
       eltoel(3,ielem2)=ineigh21
       lelem(ielem2)=0_ip

       lface(1,ielem4)=ip4  
       lface(2,ielem4)=ip3  
       lface(3,ielem4)=ipnew
       eltoel(1,ielem4)=ielem3
       eltoel(2,ielem4)=ielem2
       eltoel(3,ielem4)=ineigh22
       lelem(ielem4)=0_ip

       lface(1,ielem3)=ip3  
       lface(2,ielem3)=ip1  
       lface(3,ielem3)=ipnew
       eltoel(1,ielem3)=ielem1
       eltoel(2,ielem3)=ielem4
       eltoel(3,ielem3)=ineigh11
       lelem(ielem3)=0_ip


       !lptri(ip3)=ielem3

       !
       !     Update outside
       !
       if(ineigh22/=0)then
          if(eltoel(1,ineigh22)==ielem2)then
             eltoel(1,ineigh22)=ielem4
          else if(eltoel(2,ineigh22)==ielem2)then
             eltoel(2,ineigh22)=ielem4
          else
             eltoel(3,ineigh22)=ielem4
          endif
          !
          !     Update edge
          !  
          if(ip3<npcorner .and. ip4<npcorner)then
             if(ip3<ip4)then
                ipmin=ip3
                ipmax=ip4
             else
                ipmin=ip4
                ipmax=ip3
             endif
             !
             !     Look for ipmax in ipmin
             ! 
             do iedge=lptoedg2(ipmin),lptoedg2(ipmin+1)-1
                if(ledge(2,iedge)==ipmax)then
                   ledge(3,iedge)=ielem4
                   exit
                endif
             enddo
             !
             !     Any error? This is a ridge
             !  
             !if(iedge==lptoedg2(ipmin+1))then
             !   write(*,*)'Error in split, edge not found 1'
             !   stop
             !endif
          endif

       endif

       if(ineigh11/=0)then 
          if(eltoel(1,ineigh11)==ielem1)then
             eltoel(1,ineigh11)=ielem3
          else if(eltoel(2,ineigh11)==ielem1)then
             eltoel(2,ineigh11)=ielem3
          else
             eltoel(3,ineigh11)=ielem3
          endif
          !
          !     Update edge
          ! 
          if(ip1<npcorner .and. ip3<npcorner)then
             if(ip1<ip3)then
                ipmin=ip1
                ipmax=ip3
             else
                ipmin=ip3
                ipmax=ip1
             endif
             !
             !     Look for ipmax in ipmin
             ! 
             do iedge=lptoedg2(ipmin),lptoedg2(ipmin+1)-1
                if(ledge(2,iedge)==ipmax)then
                   ledge(3,iedge)=ielem3
                   exit
                endif
             enddo
             !
             !     Any error? this is a ridge
             !  
             !if(iedge==lptoedg2(ipmin+1))then
             !   write(*,*)'Error in split, edge not found 2'
             !   stop
             !endif
          endif
       endif

    else
       !
       !     Allocate 1 new element
       !
       nfnew=nface+1
       call memrea(nfnew,memor_msh,'Lface','split',lface)
       call memrea(nfnew,memor_msh,'ELTOEL','split',eltoel)
       call memrea(nfnew,memor_msh,'ELTOEL','split',lelem)
       nface=nface+1
       ielem3=nface
       !
       !     Create new elements
       ! 
       lface(1,ielem1)=ip1  
       lface(2,ielem1)=ip2  
       lface(3,ielem1)=ipnew
       eltoel(1,ielem1)=0_ip
       eltoel(2,ielem1)=ielem3
       eltoel(3,ielem1)=ineigh12
       lelem(ielem1)=0_ip

       lface(1,ielem3)=ip3  
       lface(2,ielem3)=ip1  
       lface(3,ielem3)=ipnew
       eltoel(1,ielem3)=ielem1
       eltoel(2,ielem3)=0_ip
       eltoel(3,ielem3)=ineigh11
       lelem(ielem3)=0_ip
       !
       !     Update outside
       !

       if(ineigh11/=0)then 
          if(eltoel(1,ineigh11)==ielem1)then
             eltoel(1,ineigh11)=ielem3
          else if(eltoel(2,ineigh11)==ielem1)then
             eltoel(2,ineigh11)=ielem3
          else
             eltoel(3,ineigh11)=ielem3
          endif
       endif

    endif

  end subroutine split

  subroutine splitsid(ielem1,iview1,npoin,ndim,coorold,nface,&
       nfold,npold,rnofaold,lcell,ncell,lelemold,tolsplit,ichk,sitosiold,lsold,nsold,&
       nnosi,lstof,ptosi2old,ptosi1old,rsuni,ptoel1,ptoel2,lpsid, &
       lface,eltoel,lcart,rnopo,lpsur,rsize,&
       lptri,lmark,lpofa,lsurf,coor,rtol,lptype,ierr,iline,llinold )
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(in)    :: ndim,ielem1,iview1,nfold,npold,ncell,nsold,nnosi,iline
    integer(ip),intent(inout) :: npoin,ichk,nface,ierr
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),lpsid(:),lface(:,:),eltoel(:,:)
    integer(ip),pointer       :: lcart(:),lpsur(:),lptri(:),lmark(:) 
    integer(ip),pointer       :: lsurf(:),lpofa(:),lptype(:,:)
    real(rp),pointer          :: rnopo(:,:),rsize(:),coor(:,:)
    integer(ip),intent(in)    :: lstof(nsold)
    integer(ip),intent(in)    :: sitosiold(nnosi,nsold),lsold(nnosi,nsold),ptosi2old(npold+1)
    integer(ip),intent(in)    :: ptosi1old(*),llinold(nsold)
    integer(ip),intent(inout) :: lelemold(nfold)
    type(cell),intent(in)     :: lcell(ncell) 
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold),tolsplit,rsuni,rtol 
    integer(ip)               :: ip1,ip2,ip3,ip4,iview2,ip1o,ip2o,ip3o
    integer(ip)               :: ipa,ipb,ipc,iside,is,ie,ielem
    integer(ip)               :: ineigh11,ineigh12,ineigh21,ineigh22,ipnew,icart,iguess
    integer(ip)               :: isurf1,isurf2,nfnew,ielem3,ielem4,npnew,ihost,ielem2,ihostf
    real(rp)                  :: d1,d2,c05,rlen,rx,ry,rz,rx1,ry1,rz1,rl,scamax,c10,sca,rstot,rs2,rs3
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub splits the ridge separing elements ielem1 and ielem2
    !     by inserting a point at the middle of the edge
    !     Be careful: lptype is updated outside
    !                 eltoel is updated
    !
    c05=0.5d+00
    c10=1.0d+00
    !if(eltoel(1,ielem1)==ielem2)then
    !   iview1=1_ip
    !else if(eltoel(2,ielem1)==ielem2)then
    !   iview1=2_ip
    !else if(eltoel(3,ielem1)==ielem2)then
    !   iview1=3_ip
    !else
    !   write(*,*)'error in split 1'
    !   stop  
    !endif

    ineigh11=eltoel(ltab(1,iview1),ielem1)
    ineigh12=eltoel(ltab(2,iview1),ielem1)

    ip1=lface(iview1,ielem1)
    ip2=lface(ltab(1,iview1),ielem1)
    ip3=lface(ltab(2,iview1),ielem1)

    ielem2=eltoel(iview1,ielem1)

    if(ielem2/=0)then

       if(eltoel(1,ielem2)==ielem1)then
          iview2=1_ip
       else if(eltoel(2,ielem2)==ielem1)then
          iview2=2_ip
       else if(eltoel(3,ielem2)==ielem1)then
          iview2=3_ip
       else
          write(*,*)'error in split 2'
          stop  
       endif

       ip4=lface(iview2,ielem2)
       ineigh21=eltoel(ltab(1,iview2),ielem2)
       ineigh22=eltoel(ltab(2,iview2),ielem2)

    endif
    !
    !     Allocate data base for ipnew
    !
    npnew=npoin+1
    ipnew=npnew
    call resizep(npnew,coor,rnopo,lcart,lpsur,rsize,lptri,lmark,lpofa,lptype,lpsid)
    npoin=npnew
    !
    !     Update the point database for ipnew
    !
    rstot=rsize(ip2)+rsize(ip3)
    rs2=rsize(ip3)/rstot
    rs3=rsize(ip2)/rstot
    coor(1,ipnew)=rs2*coor(1,ip2)+rs3*coor(1,ip3)
    coor(2,ipnew)=rs2*coor(2,ip2)+rs3*coor(2,ip3)
    coor(3,ipnew)=rs2*coor(3,ip2)+rs3*coor(3,ip3)

    !write(*,*)'ipnew=',ipnew,ip2,ip3
    if(lpsid(ip2)/=0)then
       iguess=lpsid(ip2) 
    else if(lpsid(ip3)/=0)then
       iguess=lpsid(ip3) 
    else
       !
       !     It should be a side with two corners
       !
       do is=ptosi2old(ip2),ptosi2old(ip2+1)-1
          iside=ptosi1old(is)
          ipa=lsold(1,iside) 
          ipb=lsold(2,iside)
          if(ipa==ip2)then
             ipc=ipb
          else
             ipc=ipa
          endif
          if(ipc==ip3)then 
             iguess=iside
             goto 10
          endif
       enddo
       !
       !     It still could be a collapsed side with two corners
       !     Take the side mostly aligned with (ip2,ip3)
       ! 
       rx=coor(1,ip3)-coor(1,ip2)
       ry=coor(2,ip3)-coor(2,ip2)
       rz=coor(3,ip3)-coor(3,ip2)
       rl=sqrt(rx*rx+ry*ry+rz*rz)
       rl=c10/rl
       rx=rx*rl
       ry=ry*rl
       rz=rz*rl
       scamax=-1.01d+00

       do is=ptosi2old(ip2),ptosi2old(ip2+1)-1
          iside=ptosi1old(is)
          ipa=lsold(1,iside) 
          ipb=lsold(2,iside)
          if(ipa==ip2)then
             ipc=ipb
          else
             ipc=ipa
          endif

          rx1=coorold(1,ipc)-coorold(1,ip2)
          ry1=coorold(2,ipc)-coorold(2,ip2)
          rz1=coorold(3,ipc)-coorold(3,ip2)
          rl=sqrt(rx1*rx1+ry1*ry1+rz1*rz1)
          rl=c10/rl
          rx1=rx1*rl
          ry1=ry1*rl
          rz1=rz1*rl

          sca=rx*rx1+ry*ry1+rz*rz1 

          if(sca>scamax)then
             scamax=sca
             iguess=iside
          endif
       enddo

10     continue

    endif
    call gthostsid(coor(:,ipnew),nsold,nnosi,npold,iguess,ihost,ndim,lsold,sitosiold,&
         coorold,d1,d2,ptosi2old,ptosi1old,npoin,lelemold,ierr,llinold,iline)
    if(ierr==1)then
       write(*,*)'Error in splitsid',ip2,ip3
       return 
    endif
    !
    !     Get a face attached to this side
    !
    ihostf=lstof(ihost)
    !
    !     Get the projected point
    !
    !dtot=d1+d2+d3
    !d1=d1/dtot  
    !d2=d2/dtot  
    !d3=d3/dtot  
    ip1o=lsold(1,ihost)
    ip2o=lsold(2,ihost)

    !write(*,*)'split'
    !write(*,*)d1,d2,d3
    !write(*,*)ihost,ip1o,ip2o,ip3o

    coor(1,ipnew)=d1*coorold(1,ip1o)+d2*coorold(1,ip2o)
    coor(2,ipnew)=d1*coorold(2,ip1o)+d2*coorold(2,ip2o)
    coor(3,ipnew)=d1*coorold(3,ip1o)+d2*coorold(3,ip2o)
    rnopo(1,ipnew)=rnofaold(1,ihostf)
    rnopo(2,ipnew)=rnofaold(2,ihostf)
    rnopo(3,ipnew)=rnofaold(3,ihostf)
    icart=lcart(ip2)
    !icart=lcart(ip1)

    !write(*,*)coor(1,ipnew),coor(2,ipnew),coor(3,ipnew)
    !write(*,*)ihost,ip1o,ip2o,ip3o

    call gtelem(ipnew,coor,npoin,ndim,lcell,ncell,icart,rtol)
    lcart(ipnew)=icart
    call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipnew,ndim,coor,rsuni)
    lpsur(ipnew)=ihostf
    lpsid(ipnew)=ihost 
    lptri(ipnew)=ielem1
    lmark(ipnew)=0_ip
    lpofa(ipnew)=0_ip
    !
    !     Do we want to check the length of the newly created edges?
    !
    if(ichk==1)then
       call length(ip1,ipnew,rsize,coor,rlen)
       if(rlen<tolsplit)then
          npoin=npoin-1
          ichk=0_ip 
          return
       endif

       if(ielem2/=0)then
          call length(ip4,ipnew,rsize,coor,rlen)
          if(rlen<tolsplit)then
             npoin=npoin-1
             ichk=0_ip 
             return
          endif
       endif
    endif

    if(ielem2/=0)then
       !
       !     Allocate the 2 new elements
       !
       nfnew=nface+2
       call memrea(nfnew,memor_msh,'LFLOC','split',lface)
       call memrea(nfnew,memor_msh,'LSURF','split',lsurf)
       call memrea(nfnew,memor_msh,'ELTOEL','split',eltoel)
       nface=nface+1
       ielem3=nface
       nface=nface+1
       ielem4=nface

       isurf1=lsurf(ielem1)
       isurf2=lsurf(ielem2)

       !
       !     Create new elements
       ! 
       lface(1,ielem1)=ip1  
       lface(2,ielem1)=ip2  
       lface(3,ielem1)=ipnew
       eltoel(1,ielem1)=ielem2
       eltoel(2,ielem1)=ielem3
       eltoel(3,ielem1)=ineigh12
       lsurf(ielem1)=isurf1

       lface(1,ielem2)=ip2  
       lface(2,ielem2)=ip4  
       lface(3,ielem2)=ipnew
       eltoel(1,ielem2)=ielem4
       eltoel(2,ielem2)=ielem1
       eltoel(3,ielem2)=ineigh21
       lsurf(ielem2)=isurf2

       lface(1,ielem4)=ip4  
       lface(2,ielem4)=ip3  
       lface(3,ielem4)=ipnew
       eltoel(1,ielem4)=ielem3
       eltoel(2,ielem4)=ielem2
       eltoel(3,ielem4)=ineigh22
       lsurf(ielem4)=isurf2

       lface(1,ielem3)=ip3  
       lface(2,ielem3)=ip1  
       lface(3,ielem3)=ipnew
       eltoel(1,ielem3)=ielem1
       eltoel(2,ielem3)=ielem4
       eltoel(3,ielem3)=ineigh11
       lsurf(ielem3)=isurf1


       !lptri(ip3)=ielem3

       !
       !     Update outside
       !
       if(ineigh22/=0)then
          if(eltoel(1,ineigh22)==ielem2)then
             eltoel(1,ineigh22)=ielem4
          else if(eltoel(2,ineigh22)==ielem2)then
             eltoel(2,ineigh22)=ielem4
          else
             eltoel(3,ineigh22)=ielem4
          endif
       endif

       if(ineigh11/=0)then 
          if(eltoel(1,ineigh11)==ielem1)then
             eltoel(1,ineigh11)=ielem3
          else if(eltoel(2,ineigh11)==ielem1)then
             eltoel(2,ineigh11)=ielem3
          else
             eltoel(3,ineigh11)=ielem3
          endif
       endif
       !
       !     Must update ptoel1 if more than one side should be refined in the same 
       !     triangle
       !
       do ie=ptoel2(ip3),ptoel2(ip3+1)-1
          ielem=ptoel1(ie)
          if(ielem==ielem1)then
             ptoel1(ie)=ielem3
             exit
          endif
       enddo

       do ie=ptoel2(ip3),ptoel2(ip3+1)-1
          ielem=ptoel1(ie)
          if(ielem==ielem2)then
             ptoel1(ie)=ielem4
             exit
          endif
       enddo

    else
       !
       !     Allocate 1 new element
       !
       nfnew=nface+1
       call memrea(nfnew,memor_msh,'LFLOC','split',lface)
       call memrea(nfnew,memor_msh,'LSURF','split',lsurf)
       call memrea(nfnew,memor_msh,'ELTOEL','split',eltoel)
       nface=nface+1
       ielem3=nface

       isurf1=lsurf(ielem1)

       !
       !     Create new elements
       ! 
       lface(1,ielem1)=ip1  
       lface(2,ielem1)=ip2  
       lface(3,ielem1)=ipnew
       eltoel(1,ielem1)=0_ip
       eltoel(2,ielem1)=ielem3
       eltoel(3,ielem1)=ineigh12
       lsurf(ielem1)=isurf1

       lface(1,ielem3)=ip3  
       lface(2,ielem3)=ip1  
       lface(3,ielem3)=ipnew
       eltoel(1,ielem3)=ielem1
       eltoel(2,ielem3)=0_ip
       eltoel(3,ielem3)=ineigh11
       lsurf(ielem3)=isurf1
       !
       !     Update outside
       !
       if(ineigh11/=0)then 
          if(eltoel(1,ineigh11)==ielem1)then
             eltoel(1,ineigh11)=ielem3
          else if(eltoel(2,ineigh11)==ielem1)then
             eltoel(2,ineigh11)=ielem3
          else
             eltoel(3,ineigh11)=ielem3
          endif
       endif
       !
       !     Must update ptoel1 if more than one side should be refined in the same 
       !     triangle
       !
       do ie=ptoel2(ip3),ptoel2(ip3+1)-1
          ielem=ptoel1(ie)
          if(ielem==ielem1)then
             ptoel1(ie)=ielem3
             exit
          endif
       enddo

    endif

  end subroutine splitsid

  subroutine split2(ielem1,iview1,nnofa,nface,npoin,ndim,ipnew,lptri,&
       lface,eltoel,lfmark,rnofa,lelem,coor,rnopo,ierr,nehole,lehole )
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only     :  memor_msh
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,ielem1,iview1,ipnew,npoin
    integer(ip),intent(inout) :: nface,ierr,nehole
    integer(ip),intent(inout) :: lptri(npoin),lehole(nehole)
    integer(ip),pointer       :: lface(:,:),eltoel(:,:),lfmark(:),lelem(:)
    real(rp),pointer          :: rnofa(:,:)
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin)
    integer(ip)               :: ip1,ip2,ip3,ip4,iview2,ielem2
    integer(ip)               :: ineigh11,ineigh12,ineigh21,ineigh22
    integer(ip)               :: nfnew,ielem3,ielem4
    real(rp)                  :: modul,qgeo2,qgeotol 
    integer(ip)               :: lfloc(3,4)
    real(rp)                  :: rnofloc(3,4)
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub splits the edge separing elements ielem1 and ielem2
    !     by inserting a point near the middle of the edge
    !     Compared to split, point ipnew and its database have already been computed
    !
    !
    qgeotol=0.1d+00
    ielem2=eltoel(iview1,ielem1)  

    ip1=lface(iview1,ielem1)
    ip2=lface(ltab(1,iview1),ielem1)
    ip3=lface(ltab(2,iview1),ielem1)

    if(eltoel(1,ielem2)==ielem1)then
       iview2=1_ip
    else if(eltoel(2,ielem2)==ielem1)then
       iview2=2_ip
    else
       iview2=3_ip
    endif

    ip4=lface(iview2,ielem2)
    !
    !     Create new elements
    !
    lfloc(1,1)=ip1
    lfloc(2,1)=ip2
    lfloc(3,1)=ipnew
    call gtfnr2(lfloc,4_ip,nnofa,ndim,coor,npoin,rnofloc(1,1),1_ip,modul)
    call chkgeo(lfloc,4_ip,nnofa,1_ip,qgeo2,rnofloc(1,1),rnopo,npoin,ndim)
    if(qgeo2<qgeotol)then
       ierr=1_ip
       return
    endif

    lfloc(1,2)=ip2
    lfloc(2,2)=ip4
    lfloc(3,2)=ipnew
    call gtfnr2(lfloc,4_ip,nnofa,ndim,coor,npoin,rnofloc(1,2),2_ip,modul)
    call chkgeo(lfloc,4_ip,nnofa,2_ip,qgeo2,rnofloc(1,2),rnopo,npoin,ndim)
    if(qgeo2<qgeotol)then
       ierr=1_ip
       return
    endif

    lfloc(1,3)=ip4
    lfloc(2,3)=ip3
    lfloc(3,3)=ipnew
    call gtfnr2(lfloc,4_ip,nnofa,ndim,coor,npoin,rnofloc(1,3),3_ip,modul)
    call chkgeo(lfloc,4_ip,nnofa,3_ip,qgeo2,rnofloc(1,3),rnopo,npoin,ndim)
    if(qgeo2<qgeotol)then
       ierr=1_ip
       return
    endif

    lfloc(1,4)=ip3
    lfloc(2,4)=ip1
    lfloc(3,4)=ipnew
    call gtfnr2(lfloc,4_ip,nnofa,ndim,coor,npoin,rnofloc(1,4),4_ip,modul)
    call chkgeo(lfloc,4_ip,nnofa,4_ip,qgeo2,rnofloc(1,4),rnopo,npoin,ndim)
    if(qgeo2<qgeotol)then
       ierr=1_ip
       return
    endif
    !
    !     We can split
    !
    ineigh11=eltoel(ltab(1,iview1),ielem1)
    ineigh12=eltoel(ltab(2,iview1),ielem1)
    ineigh21=eltoel(ltab(1,iview2),ielem2)
    ineigh22=eltoel(ltab(2,iview2),ielem2)
    !
    !     Allocate the 2 new elements
    !
    if(nehole>0)then
       ielem3=lehole(nehole) 
       nehole=nehole-1
    else 
       ielem3=nface+1
       call memrea(ielem3,memor_msh,'LFACE','insert',lface)
       call memrea(ielem3,memor_msh,'RNOFA','insert',rnofa)
       call memrea(ielem3,memor_msh,'ELTOEL','insert',eltoel)
       call memrea(ielem3,memor_msh,'LELEM','insert',lelem) 
       call memrea(ielem3,memor_msh,'LFMARK','insert',lfmark) 
       nface=ielem3
    endif
    
    if(nehole>0)then
       ielem4=lehole(nehole) 
       nehole=nehole-1
    else 
       ielem4=nface+1
       call memrea(ielem4,memor_msh,'LFACE','insert',lface)
       call memrea(ielem4,memor_msh,'RNOFA','insert',rnofa)
       call memrea(ielem4,memor_msh,'ELTOEL','insert',eltoel)
       call memrea(ielem4,memor_msh,'LELEM','insert',lelem) 
       call memrea(ielem4,memor_msh,'LFMARK','insert',lfmark) 
       nface=ielem4
    endif

    !isurf1=lsurf(ielem1)
    !isurf2=lsurf(ielem2)
    !
    !     Create new elements
    !
    lface(1,ielem1)=lfloc(1,1)
    lface(2,ielem1)=lfloc(2,1)
    lface(3,ielem1)=lfloc(3,1) 
    eltoel(1,ielem1)=ielem2
    eltoel(2,ielem1)=ielem3
    eltoel(3,ielem1)=ineigh12
    !lsurf(ielem1)=isurf1
    lelem(ielem1)=0_ip
    rnofa(1,ielem1)=rnofloc(1,1)
    rnofa(2,ielem1)=rnofloc(2,1)
    rnofa(3,ielem1)=rnofloc(3,1)
    lfmark(ielem1)=0_ip

    lface(1,ielem2)=lfloc(1,2)
    lface(2,ielem2)=lfloc(2,2)
    lface(3,ielem2)=lfloc(3,2)
    eltoel(1,ielem2)=ielem4
    eltoel(2,ielem2)=ielem1
    eltoel(3,ielem2)=ineigh21
    !lsurf(ielem2)=isurf2
    lelem(ielem2)=0_ip
    rnofa(1,ielem2)=rnofloc(1,2)
    rnofa(2,ielem2)=rnofloc(2,2)
    rnofa(3,ielem2)=rnofloc(3,2)
    lfmark(ielem2)=0_ip

    lface(1,ielem4)=lfloc(1,3)
    lface(2,ielem4)=lfloc(2,3)
    lface(3,ielem4)=lfloc(3,3)
    eltoel(1,ielem4)=ielem3
    eltoel(2,ielem4)=ielem2
    eltoel(3,ielem4)=ineigh22
    !lsurf(ielem4)=isurf2
    lelem(ielem4)=0_ip
    rnofa(1,ielem4)=rnofloc(1,3)
    rnofa(2,ielem4)=rnofloc(2,3)
    rnofa(3,ielem4)=rnofloc(3,3)
    lfmark(ielem4)=0_ip

    lface(1,ielem3)=lfloc(1,4)
    lface(2,ielem3)=lfloc(2,4)
    lface(3,ielem3)=lfloc(3,4)
    eltoel(1,ielem3)=ielem1
    eltoel(2,ielem3)=ielem4
    eltoel(3,ielem3)=ineigh11
    !lsurf(ielem3)=isurf1
    lelem(ielem3)=0_ip
    rnofa(1,ielem3)=rnofloc(1,4)
    rnofa(2,ielem3)=rnofloc(2,4)
    rnofa(3,ielem3)=rnofloc(3,4)
    lfmark(ielem3)=0_ip

    lptri(ip3)=ielem3
    lptri(ipnew)=ielem1
    !
    !     Update outside
    !
    if(ineigh22/=0)then
       if(eltoel(1,ineigh22)==ielem2)then
          eltoel(1,ineigh22)=ielem4
       else if(eltoel(2,ineigh22)==ielem2)then
          eltoel(2,ineigh22)=ielem4
       else
          eltoel(3,ineigh22)=ielem4
       endif
    endif

    if(ineigh11/=0)then 
       if(eltoel(1,ineigh11)==ielem1)then
          eltoel(1,ineigh11)=ielem3
       else if(eltoel(2,ineigh11)==ielem1)then
          eltoel(2,ineigh11)=ielem3
       else
          eltoel(3,ineigh11)=ielem3
       endif
    endif

  end subroutine split2

  subroutine split3(ielem1,ielem2,nnofa,nface,npoin,ndim,coorold,eltoelold,lfold,nfold,&
       npold,rnofaold,lcell,ncell,lelemold,rsuni,lsurfold,isurf,&
       lface,eltoel,lcart,rnopo,rnofa,lpsur,rsize,lptri,lmark,&
       lelem,lpofa,lsurf,coor,rtol,lptype,lpsid)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only     :  memor_msh
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,ielem1,ielem2,nfold,npold,ncell,isurf
    integer(ip),intent(in)    :: eltoelold(nnofa,nfold),lfold(nnofa,nfold),lsurfold(nfold)
    type(cell),intent(in)     :: lcell(ncell) 
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold),rsuni,rtol 
    integer(ip),intent(inout) :: nface,npoin,lelemold(nfold)
    integer(ip),pointer       :: lface(:,:),eltoel(:,:),lcart(:),lpsur(:),lptri(:) 
    integer(ip),pointer       :: lmark(:),lelem(:),lpofa(:),lsurf(:),lptype(:,:),lpsid(:)
    real(rp),pointer          :: rnopo(:,:),rnofa(:,:),rsize(:),coor(:,:) 
    integer(ip)               :: ip1,ip2,ip3,ip4,iview1,iview2,ip1o,ip2o,ip3o
    integer(ip)               :: ineigh11,ineigh12,ineigh21,ineigh22,ipnew,icart
    integer(ip)               :: isurf1,isurf2,nfnew,ielem3,ielem4,npnew,ihost
    integer(ip)               :: aneigh11,aneigh12,aneigh21,aneigh22,ierr
    real(rp)                  :: d1,d2,d3,dtot,c05,rsiz
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub splits the edge separing elements ielem1 and ielem2
    !     by inserting a point at the middle of the edge
    !     be careful: lptype is updated outside
    !     same as split but handle negative entries in eltoel
    !
    !     VOLUMES ARE NOT CHECKED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     THIS SUB IS NOT USED AND NOT CALLED FROM ANYWHERE
    !
    !  
    c05=0.5d+00

    if(eltoel(1,ielem1)==ielem2)then
       iview1=1_ip
    else if(eltoel(2,ielem1)==ielem2)then
       iview1=2_ip
    else if(eltoel(3,ielem1)==ielem2)then
       iview1=3_ip
    else
       write(*,*)'error in split3 1'
       stop  
    endif

    ip1=lface(iview1,ielem1)
    ip2=lface(ltab(1,iview1),ielem1)
    ip3=lface(ltab(2,iview1),ielem1)

    if(eltoel(1,ielem2)==ielem1)then
       iview2=1_ip
    else if(eltoel(2,ielem2)==ielem1)then
       iview2=2_ip
    else if(eltoel(3,ielem2)==ielem1)then
       iview2=3_ip
    else
       write(*,*)'error in split3 2'
       write(*,*)eltoel(1,ielem1),eltoel(2,ielem1),eltoel(3,ielem1)
       write(*,*)eltoel(1,ielem2),eltoel(2,ielem2),eltoel(3,ielem2)
       stop  
    endif

    ip4=lface(iview2,ielem2)

    aneigh11=eltoel(ltab(1,iview1),ielem1)
    aneigh12=eltoel(ltab(2,iview1),ielem1)
    aneigh21=eltoel(ltab(1,iview2),ielem2)
    aneigh22=eltoel(ltab(2,iview2),ielem2)

    ineigh11=abs(aneigh11)
    ineigh12=abs(aneigh12)
    ineigh21=abs(aneigh21)
    ineigh22=abs(aneigh22)
    !
    !     Allocate data base for ipnew
    !
    npnew=npoin+1
    ipnew=npnew
    call resizep(npnew,coor,rnopo,lcart,lpsur,rsize,lptri,lmark,lpofa,lptype,lpsid)
    npoin=npnew
    !
    !     Update the point database for ipnew
    !
    coor(1,ipnew)=(coor(1,ip2)+coor(1,ip3))*c05
    coor(2,ipnew)=(coor(2,ip2)+coor(2,ip3))*c05
    coor(3,ipnew)=(coor(3,ip2)+coor(3,ip3))*c05

    !write(*,*)'ipnew=',ipnew,ip2,ip3
    call gthost(coor(:,ipnew),lpsur(ip2),ihost,rnofaold,ndim,nfold,npold,lfold,&
         nnofa,coorold,d1,d2,d3,eltoelold,lelemold,lsurfold,isurf,ierr,rsiz)
    !
    !     Get the projected point
    !
    !dtot=d1+d2+d3
    !d1=d1/dtot  
    !d2=d2/dtot  
    !d3=d3/dtot  
    ip1o=lfold(1,ihost)
    ip2o=lfold(2,ihost)
    ip3o=lfold(3,ihost)

    !write(*,*)'split'
    !write(*,*)d1,d2,d3
    !write(*,*)ihost,ip1o,ip2o,ip3o

    coor(1,ipnew)=d1*coorold(1,ip1o)+d2*coorold(1,ip2o)+d3*coorold(1,ip3o)
    coor(2,ipnew)=d1*coorold(2,ip1o)+d2*coorold(2,ip2o)+d3*coorold(2,ip3o)
    coor(3,ipnew)=d1*coorold(3,ip1o)+d2*coorold(3,ip2o)+d3*coorold(3,ip3o)
    rnopo(1,ipnew)=rnofaold(1,ihost)
    rnopo(2,ipnew)=rnofaold(2,ihost)
    rnopo(3,ipnew)=rnofaold(3,ihost)
    icart=lcart(ip1)

    !write(*,*)coor(1,ipnew),coor(2,ipnew),coor(3,ipnew)
    !write(*,*)ihost,ip1o,ip2o,ip3o

    call gtelem(ipnew,coor,npoin,ndim,lcell,ncell,icart,rtol)
    lcart(ipnew)=icart
    call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipnew,ndim,coor,rsuni)
    lpsur(ipnew)=ihost 
    lptri(ipnew)=ielem1
    lmark(ipnew)=0_ip
    lpofa(ipnew)=0_ip
    !
    !     Allocate the 2 new elements
    !
    nfnew=nface+2
    call memrea(nfnew,memor_msh,'LFACE','split',lface)
    call memrea(nfnew,memor_msh,'LSURF','split',lsurf)
    call memrea(nfnew,memor_msh,'ELTOEL','split',eltoel)
    nface=nface+1
    ielem3=nface
    nface=nface+1
    ielem4=nface

    isurf1=lsurf(ielem1)
    isurf2=lsurf(ielem2)
    !
    !     Create new elements && mark the new sides negatively
    ! 
    lface(1,ielem1)=ip1  
    lface(2,ielem1)=ip2  
    lface(3,ielem1)=ipnew
    eltoel(1,ielem1)=-ielem2
    eltoel(2,ielem1)=-ielem3
    if(aneigh12>0)then
       eltoel(3,ielem1)=ineigh12
    else
       eltoel(3,ielem1)=-ineigh12
    endif
    lsurf(ielem1)=isurf1

    lface(1,ielem2)=ip2  
    lface(2,ielem2)=ip4  
    lface(3,ielem2)=ipnew
    eltoel(1,ielem2)=-ielem4
    eltoel(2,ielem2)=-ielem1
    if(aneigh21>0)then
       eltoel(3,ielem2)=ineigh21
    else
       eltoel(3,ielem2)=-ineigh21
    endif
    lsurf(ielem2)=isurf2

    lface(1,ielem4)=ip4  
    lface(2,ielem4)=ip3  
    lface(3,ielem4)=ipnew
    eltoel(1,ielem4)=-ielem3
    eltoel(2,ielem4)=-ielem2
    eltoel(3,ielem4)=ineigh22
    lsurf(ielem4)=isurf2

    lface(1,ielem3)=ip3  
    lface(2,ielem3)=ip1  
    lface(3,ielem3)=ipnew
    eltoel(1,ielem3)=-ielem1
    eltoel(2,ielem3)=-ielem4
    eltoel(3,ielem3)=ineigh11
    lsurf(ielem3)=isurf1

    !lptri(ip3)=ielem3

    !
    !     Update outside
    !
    if(ineigh22/=0)then
       if(abs(eltoel(1,ineigh22))==ielem2)then
          iview1=1
       else if(abs(eltoel(2,ineigh22))==ielem2)then
          iview1=2
       else
          iview1=3
       endif

       if(aneigh22>0)then
          eltoel(iview1,ineigh22)=ielem4
       else 
          eltoel(iview1,ineigh22)=-ielem4
          eltoel(3,ielem4)=-ineigh22
       endif
    endif

    if(ineigh11/=0)then 
       if(abs(eltoel(1,ineigh11))==ielem1)then
          iview1=1
       else if(abs(eltoel(2,ineigh11))==ielem1)then
          iview1=2
       else
          iview1=3
       endif

       if(aneigh11>0)then
          eltoel(iview1,ineigh11)=ielem3
       else   
          eltoel(iview1,ineigh11)=-ielem3
          eltoel(3,ielem3)=-ineigh11
       endif
    endif

  end subroutine split3

  subroutine refsmo(nface,nnofa,npoin,ndim,lfold,nfold,coorold,npold,lcell,ncell,eltoelold,&
       rnofaold,isurf,lelemold,lsurfold,rsuni,rnopo,lcart,lpsur,rsize,lface,eltoel,&
       ptoel1,ptoel2,lptri,lfmark,rnofa,lmark,lelem,lpofa,lptype,coor,rtol,lpsid,&
       ptosi1,ptosi2,npcusp,lside,nside,nnosi)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    use def_meshin, only     :  memor_msh
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,npold,nfold,ncell,isurf,nside,nnosi
    integer(ip),intent(inout) :: nface,npoin
    integer(ip),intent(in)    :: lfold(nnofa,nfold),lsurfold(nfold)
    integer(ip),intent(in)    :: eltoelold(nnofa,nfold)
    integer(ip),intent(inout) :: lelemold(nfold)
    type(cell)                :: lcell(ncell)
    integer(ip),pointer       :: lcart(:),lpsur(:),lfront(:,:),lface(:,:),eltoel(:,:)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),lptri(:),lfmark(:),lheap(:),lside(:,:) 
    integer(ip),pointer       :: lmark(:),lelem(:),lfapo(:,:),lfahol(:),lfhole(:) 
    integer(ip),pointer       :: lptype(:,:),lpofa(:),lpsid(:),ptosi1(:),ptosi2(:) 
    integer(ip),pointer       :: lphole(:),lehole(:) 
    real(rp),pointer          :: rnopo(:,:),rsize(:),rnofa(:,:),rfront(:),coor(:,:) 
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold),rsuni,rtol
    integer(ip)               :: iface,ipoin,ifront,nfront,p1,p2,j,icart,ifhole,iehole
    integer(ip)               :: ihost,ihostn,ip1,ip2,ip3,ipclos,nfnew,ipa,ipb,ipc,isto 
    integer(ip)               :: jpmin,jpmax,jnofa,ipmin,ipmax,iside,niterprint,nprint
    integer(ip)               :: nphole,nehole,iphole,npoi0,inofa,nfac0
    integer(ip)               :: ifron,ismall,ifnew,ipnew,ienew,iview,ierr,npcusp,ismoo 
    integer(ip)               :: ip1o,ip2o,ip3o,i,ipos,nfapo,nfahol,nfhole,iter,nheap,ipcros
    real(rp)                  :: pnew(ndim),d1,d2,d3,dtot,rpoin(ndim),rx,ry,rz,rlen
    integer(4)                :: istat 
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub meshes a discrete patch with an advancing front approach
    !
    !
    !     On input:   
    !                 - lface: the current mesh
    !                 - lptri: a face number of the current mesh for each point 
    !                 - lpsur: a face number of the old mesh for each point 
    !                 - lfold: the old mesh
    !
    !     On output: 
    !                 - lface: the new mesh
    !
    !
    !
    !     Do we have some faces?
    !     Could have been completely collapsed
    !
    if(nface==0)return
    !
    !    Allocate database arrays
    !
    if(.not.associated(lfmark))then
       allocate(lfmark(nface),stat=istat)
       call memchk(zero,istat,memor_msh,'LFMARK','refsmo',lfmark)
    else
       call memrea(nface,memor_msh,'LFMARK','refsmo',lfmark)
       do iface=1,nface
          lfmark(iface)=0_ip
       enddo
    endif
    if(.not.associated(lelem))then
       allocate(lelem(nface),stat=istat)
       call memchk(zero,istat,memor_msh,'LELEM','refsmo',lelem)
    else
       call memrea(nface,memor_msh,'LELEM','refsmo',lelem)
       do iface=1,nface
          lelem(iface)=0_ip
       enddo
    endif
    if(.not.associated(rnofa))then
       allocate(rnofa(ndim,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'RNOFA','refsmo',rnofa)
    else
       call memrea(nface,memor_msh,'RNOFA','refsmo',rnofa)
       do iface=1,nface  
          rnofa(1,iface)=0_rp
          rnofa(2,iface)=0_rp
          rnofa(3,iface)=0_rp
       enddo
    endif
    !
    !     Set error flag
    !
    ierr=0_ip
    !
    !     Initialize nphole
    !
    nphole=0_ip
    allocate(lphole(10),stat=istat)
    call memchk(zero,istat,memor_msh,'LPHOLE','refsmo',lphole)
    nehole=0_ip
    allocate(lehole(10),stat=istat)
    call memchk(zero,istat,memor_msh,'LEHOLE','refsmo',lehole)
    !
    !     Clean up lpofa
    !
    do ipoin=1,npoin
       lpofa(ipoin)=0_ip
    enddo
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Get the face normals
    !
    call gtfnrl(lface,nface,nnofa,ndim,coor,npoin,rnofa)
    !
    !     Get the boundary edges
    !
    ifront=0_ip
    do iface=1,nface
       do j=1,nnofa
          if(eltoel(j,iface)==0)then
             ifront=ifront+1
          endif
       enddo
    enddo
    !
    !     Count inner sides
    !
    do iface=1,nface
       do j=1,nnofa
          ip1=lface(ltab(1,j),iface)
          ip2=lface(ltab(2,j),iface)
          !
          !     Are ip1 or ip2 a potential cusp point?
          !
          if(ip1<=npcusp .or. ip2<=npcusp)then 
             !
             !     Look for side (ip2,ip1)
             !     Loop on sides surrounding ip1
             !
             do isto=ptosi2(ip1),ptosi2(ip1+1)-1
                iside=ptosi1(isto)
                ipa=lside(1,iside)
                ipb=lside(2,iside) 
                if(ipa==ip1)then
                   ipc=ipb
                else
                   ipc=ipa
                endif
                if(ipc==ip2)then
                   !
                   !     Do we have a boundary?
                   !
                   if(eltoel(j,iface)/=0)then
                      !
                      !     This is an internal side
                      ! 
                      ifront=ifront+1_ip
                   endif
                endif
             enddo
          endif
       enddo
    enddo

    nfront=ifront
    !
    !     Initialize nfapo, nfahol && nfhole 
    ! 
    nfapo=0_ip
    nfahol=0_ip  
    nfhole=0_ip  
    !
    !     Is the front empty?
    !
    if(nfront==0)then
       !
       !     Try to find an edge large enough to split it
       !  
       call chkcons2(lface,nnofa,nface,eltoel,lptri,npoin)
       call frstfr(nface,nnofa,ndim,npoin,nfront,coorold,eltoelold,lfold,nfold,&
            npold,rnofaold,lcell,ncell,nfapo,nheap,lelemold,isurf,rsuni,lsurfold,&
            lface,eltoel,lcart,&
            rnopo,rnofa,lpsur,rsize,lfront,lptri,lfmark,lmark,rfront,lheap,lelem,lfapo,&
            lpofa,lfahol,lfhole,lptype,coor,rtol,lpsid)
       call chkcons2(lface,nnofa,nface,eltoel,lptri,npoin)

    else
       !
       !     Allocate the front array
       !
       allocate(lfront(3,nfront),stat=istat)
       call memchk(zero,istat,memor_msh,'LFRONT','refsmo',lfront) 
       allocate(rfront(nfront),stat=istat)
       call memchk(zero,istat,memor_msh,'RFRONT','refsmo',rfront) 
       allocate(lheap(nfront),stat=istat)
       call memchk(zero,istat,memor_msh,'LHEAP','refsmo',lheap) 
       allocate(lfapo(3,nfront),stat=istat)
       call memchk(zero,istat,memor_msh,'LFAPO','refsmo',lfapo) 
       allocate(lfahol(max(1_ip,nfront/10_ip)),stat=istat)
       call memchk(zero,istat,memor_msh,'LFAHOL','refsmo',lfahol) 
       allocate(lfhole(max(1_ip,nfront/10_ip)),stat=istat)
       call memchk(zero,istat,memor_msh,'LFHOLE','refsmo',lfhole) 
       !
       !     Store the front
       !
       nfront=0_ip
       nheap=0_ip
       do iface=1,nface
          do j=1,nnofa
             if(eltoel(j,iface)==0)then
                p1=lface(ltab(1,j),iface)
                p2=lface(ltab(2,j),iface)
                call addfac(nfront,nfapo,p1,p2,nfahol,coor,ndim,npoin,nfhole,nheap,&
                     lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole)
             endif
          enddo
       enddo
    endif
    !
    !     Add inner sides
    !
    do iface=1,nface
       do j=1,nnofa
          ip1=lface(ltab(1,j),iface)
          ip2=lface(ltab(2,j),iface)
          !
          !     Are ip1 or ip2 a potential cusp point?
          !
          if(ip1<=npcusp .or. ip2<=npcusp)then 
             !
             !     Loop on sides surrounding ip1
             !
             do isto=ptosi2(ip1),ptosi2(ip1+1)-1
                iside=ptosi1(isto)
                ipa=lside(1,iside)
                ipb=lside(2,iside) 
                if(ipa==ip1)then
                   ipc=ipb
                else
                   ipc=ipa
                endif
                if(ipc==ip2)then
                   !
                   !     Do we have a boundary?
                   !
                   if(eltoel(j,iface)/=0)then
                      !
                      !     This is an internal side
                      ! 
                      call addfac(nfront,nfapo,ip1,ip2,nfahol,coor,ndim,npoin,nfhole,nheap,&
                           lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole)
                   endif
                endif
             enddo
          endif
       enddo
    enddo
    !
    !     Output new mesh
    !
    !call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor)
    !
    !    Set the pointer lfront->lheap
    !
    call initheap(lheap,nheap,lfront,nfront)
    !
    !     Loop on the front
    ! 
    iter=1_ip
    niterprint=10000
    nprint=niterprint

    do 

       if(nfront==0)exit

       if(iter==niterprint)then

          write(*,*)'Front iteration:',iter,'nfront=',nfront
          if(niterprint==10*nprint)then
             nprint=10*nprint
             niterprint=2*nprint
          else
             niterprint=niterprint+nprint
          endif
       endif


       !call chkarray(lpsur,npoin,rnofa,rnopo,nnofa,nface,ndim,rnofaold,nfold)
       !
       !     Get the smallest edge of the front
       !
       call gtsmall(lheap,nheap,rfront,ismall,lfront,nfront+nfhole)
       !
       !     DBG
       !
       !do ifhole=1,nfhole
       !   if(lfhole(ifhole)==ismall)then
       !      write(*,*)'error ismall'
       !      stop
       !   endif
       !enddo
       !
       !     This edge is
       !
       ip1=lfront(1,ismall)
       ip2=lfront(2,ismall)
       !write(*,*)'ip1=',ip1,'ip2=',ip2
       !write(*,*)'Front iteration:',iter,ip1,ip2
       !
       !     Get the best point
       !
       !call bestpt(ip1,ip2,rsize,pnew,rnopo,ndim,npoin,coor)
       !
       !     Find the host face in the old mesh
       !
       !call gthost(pnew,lpsur(ip1),ihost,rnofaold,ndim,nfold,npoin,lfold,nnofa,coorold,d1,d2,d3,eltoelold)
       call gthost2(ip1,ip2,pnew,lpsur(ip1),ihost,rnofaold,ndim,nfold,npold,lfold,&
            nnofa,coorold,d1,d2,d3,eltoelold,rsize,coor,npoin,lsurfold,isurf,lelemold,&
            ierr,lptri,lface,nface,eltoel)
       if(ihost==-1 .or. ierr==1)then
          write(*,*)'iter==',iter 
          write(*,*)'Front: ip1=',ip1,'ip2=',ip2 
          call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor,nehole,lehole)
          stop
       endif
       !
       !     Check if the surface is smooth locally
       !
       !call chksmo()
       !
       !     Add the point to the point array
       !
       if(nphole>0)then
          ipnew=lphole(nphole)
          nphole=nphole-1
       else  
          ipnew=npoin+1
          call resizep(ipnew,coor,rnopo,lcart,lpsur,rsize,lptri,lmark,lpofa,lptype,lpsid)
          npoin=ipnew
       endif
       !
       !     Does the newly created point belong to the surface?
       !
       if(ihost/=0)then
          !
          !     Get the projected point
          !
          !dtot=d1+d2+d3
          !d1=d1/dtot  
          !d2=d2/dtot  
          !d3=d3/dtot  
          ip1o=lfold(1,ihost)
          ip2o=lfold(2,ihost)
          ip3o=lfold(3,ihost)

          pnew(1)=d1*coorold(1,ip1o)+d2*coorold(1,ip2o)+d3*coorold(1,ip3o)
          pnew(2)=d1*coorold(2,ip1o)+d2*coorold(2,ip2o)+d3*coorold(2,ip3o)
          pnew(3)=d1*coorold(3,ip1o)+d2*coorold(3,ip2o)+d3*coorold(3,ip3o)
          rpoin(1)=rnofaold(1,ihost)
          rpoin(2)=rnofaold(2,ihost)
          rpoin(3)=rnofaold(3,ihost)
          !
          !     Update the point database
          !     
          coor(1,ipnew)=pnew(1)
          coor(2,ipnew)=pnew(2)
          coor(3,ipnew)=pnew(3)
          rnopo(1,ipnew)=rpoin(1)
          rnopo(2,ipnew)=rpoin(2)
          rnopo(3,ipnew)=rpoin(3)
          icart=lcart(ip1)
          call gtelem(ipnew,coor,npoin,ndim,lcell,ncell,icart,rtol)
          lcart(ipnew)=icart
          call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipnew,ndim,coor,rsuni)
          lpsur(ipnew)=ihost 
          !lptri(ipnew)=nface+1                          updated in insert()
          lmark(ipnew)=0_ip
          lpofa(ipnew)=0_ip
          lptype(1,ipnew)=ID_SMOOTH
          lptype(2,ipnew)=isurf

       else
          !
          !     Update the point database, only coordinate and size required by cross
          !     
          coor(1,ipnew)=pnew(1)
          coor(2,ipnew)=pnew(2)
          coor(3,ipnew)=pnew(3)
          !
          !     As a size, give distance to ip2
          !   
          rx=coor(1,ipnew)-coor(1,ip2)   
          ry=coor(2,ipnew)-coor(2,ip2)   
          rz=coor(3,ipnew)-coor(3,ip2)
          rsize(ipnew)=sqrt(rx*rx+ry*ry+rz*rz)

       endif
       !
       !     Find the host face in the new mesh
       !
       call cross(ip1,ip2,lface,nnofa,nface,lptri,coor,npoin,ndim,eltoel,rnopo,rnofa,&
            lfmark,ihostn,d1,d2,d3,lmark,lpofa,lfapo,nfapo,nfahol,lfront,nfront+nfhole,&
            ipcros,ipnew,rsize,lelem,ierr,lside,nnosi,nside,npcusp,ptosi1,ptosi2)
       if(ierr==1)then
          write(*,*)'iter==',iter 
          write(*,*)'Front: ip1=',ip1,'ip2=',ip2 
          call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor,nehole,lehole)
          stop 
       endif
       !
       !     Do we have a crossed front?
       !    
       if(ipcros/=0)then 
          !
          !     Does the point belong to an active current edge   
          !      
          if(lpofa(ipcros)==0)then
             !
             !     Try to move ipcros to ipnew
             !
             call move(coor(:,ipnew),coor,ndim,npoin,lface,nnofa,nface,ipcros,rnopo,&
                  coorold,nfold,npold,rnofaold,lfold,eltoelold,eltoel,lptri,lpsur,&
                  rnofa,lelemold,lsurfold,isurf,rsize,lmark,lelem,lpofa,lfmark,lcart,&
                  lcell,ncell,rtol,rsuni,ihostn)
          endif
          !
          !  This is the point to reach by regeneration   
          !
          nphole=nphole+1
          call memrea(nphole,memor_msh,'LPHOLE','refsmo',lphole)
          lphole(nphole)=ipnew
          ipnew=ipcros

       else
          !
          !     Insert the new point in the triangulation
          !
          call insert(nface,ihostn,nnofa,ndim,npoin,d1,d2,d3,ipnew,coor,&
               lface,eltoel,rnofa,lptri,lelem,lfmark,rnopo,lehole,nehole)
       endif
       !
       !     Regenerate the missing edges    
       !
       !call chkrnofa(lface,nface,nnofa,ndim,coor,npoin,rnofa)
       !if(ip1== .or. ip2==)call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor)
       call recover(ip1,ip2,lface,nnofa,nface,lptri,coor,npoin,ndim,eltoel,rnopo,ipnew,&
            rnofa,lfmark,ienew,iview,ierr,lside,nnosi,nside,ptosi1,ptosi2,npcusp,&
            lphole,nphole,lehole,nehole,lelem)
       if(ierr==1)then
          write(*,*)'iter==',iter 
          write(*,*)'Front: ip1=',ip1,'ip2=',ip2 
          call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor,nehole,lehole)
          stop 
       endif

       !call chkrnofa(lface,nface,nnofa,ndim,coor,npoin,rnofa)
       !
       !     Update the front
       !
       call newfront(nfapo,nfront,ismall,ipnew,npoin,nfahol,coor,ndim,nfhole,nheap,&
            lfront,lheap,lpofa,lfapo,rsize,rfront,lfahol,lfhole)
       !
       !     Swap ahead of the newly introduced element
       !
       !call swpahd(ienew,iview,lface,eltoel,ndim,npoin,coor,nnofa,nface,rnopo,lfmark,lptri,rnofa)
       !if(isurf==38)call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor)
       !call chkrnofa(lface,nface,nnofa,ndim,coor,npoin,rnofa)

       call swaploc(ienew,lface,nnofa,nface,coor,npoin,ndim,rnopo,lelem,lfmark,&
            eltoel,lptri,rnofa,lside,nnosi,nside,ptosi1,ptosi2,npcusp)
       !
       !     Output new mesh
       !
       !call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor)
       !
       !     DBG
       !

       !call chkcons2(lface,nnofa,nface,eltoel,lptri,npoin)
       !if(isurf==2)then
       !if(isurf==2)call dbgfront(nface,lface,nnofa,nfront,lfront,eltoel,lfmark,lfhole,nfhole,nfahol,&
       !     lptri,npoin,ierr,ptosi1,ptosi2,lside,nnosi,nside,nehole,lehole)
       !if(ierr==1)then
       !   write(*,*)'iter==',iter 
       !   write(*,*)'Front: ip1=',ip1,'ip2=',ip2 
       !   call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor)
       !   stop 
       !endif
       !endif
       !call chkrnofa(lface,nface,nnofa,ndim,coor,npoin,rnofa)
       iter=iter+1 
       !if(ip1==406 .or. ip2==406)call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor)
          !if(isurf==2)call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor,nehole,lehole)

    enddo
    !
    !     Compress the points
    !
    if(nphole>0)then 

       do iphole=1,nphole
          lmark(lphole(iphole))=-1       
       enddo 
      
       do iehole=1,nehole
          lfmark(lehole(iehole))=-1
       enddo     
 
       npoi0=npoin
       npoin=0   
       do ipoin=1,npoi0
          if(lmark(ipoin)==0)then
             npoin=npoin+1_ip 
             lmark(ipoin)=npoin
          endif
       enddo  
      
       do iface=1,nface
          if(lfmark(iface)/=-1)then
             do inofa=1,nnofa 
                lface(inofa,iface)=lmark(lface(inofa,iface))
             enddo
          endif  
       enddo 
       
       do ipoin=1,npoi0
          lmark(ipoin)=0
       enddo

    endif
    !
    !     Compress the faces
    !
    if(nehole>0)then
 
       do iehole=1,nehole
          lfmark(lehole(iehole))=-1
       enddo     
    
       nfac0=nface
       nface=0
     
       do iface=1,nface
          if(lfmark(iface)/=-1)then
             nface=nface+1_ip    
             lface(1,nface)=lface(1,iface)   
             lface(2,nface)=lface(2,iface)   
             lface(3,nface)=lface(3,iface)   
          endif
          lfmark(iface)=0_ip
       enddo

    endif
    !call outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lface,nnofa,coor)
    !
    !     Deallocate local arrays
    !
    call memchk(2_ip,istat,memor_msh,'LFHOLE','refsmo',lfhole)
    deallocate(lfhole,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFHOLE','refsmo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFAHOL','refsmo',lfahol)
    deallocate(lfahol,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFAHOL','refsmo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFAPO','refsmo',lfapo)
    deallocate(lfapo,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFAPO','refsmo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHEAP','refsmo',lheap)
    deallocate(lheap,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHEAP','refsmo',0_ip)
    call memchk(2_ip,istat,memor_msh,'RFRONT','refsmo',rfront)
    deallocate(rfront,stat=istat)
    if(istat/=0) call memerr(2_ip,'RFRONT','refsmo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFRONT','refsmo',lfront)
    deallocate(lfront,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFRONT','refsmo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPHOLE','refsmo',lphole)
    deallocate(lphole,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPHOLE','refsmo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEHOLE','refsmo',lehole)
    deallocate(lehole,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEHOLE','refsmo',0_ip)

  end subroutine refsmo

  subroutine frstfr(nface,nnofa,ndim,npoin,nfront,coorold,eltoelold,lfold,nfold,npold,&
       rnofaold,lcell,ncell,nfapo,nheap,lelemold,isurf,rsuni,lsurfold,lface,eltoel,lcart,&
       rnopo,rnofa,lpsur,rsize,lfront,lptri,lfmark,lmark,rfront,lheap,lelem,lfapo,&
       lpofa,lfahol,lfhole,lptype,coor,rtol,lpsid )
    use def_kintyp, only       : ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    use def_meshin, only       : memor_msh 
    implicit none

    integer(ip),intent(in)    :: nnofa,ndim,nfold,npold,ncell,isurf
    integer(ip),intent(inout) :: nface,npoin,nfront,nfapo,nheap
    integer(ip),intent(inout) :: lelemold(nfold)
    type(cell),intent(in)     :: lcell(ncell)
    integer(ip),intent(in)    :: eltoelold(nnofa,nfold),lfold(nnofa,nfold),lsurfold(nfold) 
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold),rsuni,rtol 
    integer(ip),pointer       :: lface(:,:),eltoel(:,:),lcart(:),lpsur(:),lfront(:,:)
    integer(ip),pointer       :: lptri(:),lfmark(:),lheap(:),lelem(:),lfapo(:,:),lpsid(:)
    integer(ip),pointer       :: lpofa(:),lfahol(:),lfhole(:),lptype(:,:),lmark(:)
    real(rp),pointer          :: rnopo(:,:),rnofa(:,:),rsize(:),rfront(:),coor(:,:)
    integer(ip)               :: iface,j,ip1ref,ip2ref
    integer(ip)               :: ielem1,ielem2,ineigh11,ineigh12,ineigh21,ineigh22 
    integer(ip)               :: ip1,ip2,ip3,ip4,iview1,iview2,npnew,ipnew,nfnew,ielem3,ielem4,ip1o,ip2o,ip3o,ip4o 
    integer(ip)               :: nfrontn,ihost,icart,ierr
    real(rp)                  :: rlen,tollen,c10,rlenref,c05,vect(3),rl,pnew(3),d1,d2,d3,dtot,rx,ry,rz,rface(3),modul
    real(rp)                  :: rsiz
    integer(4)                :: istat 
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))

    c10=1.0d+00
    c05=0.5d+00
    tollen=2.0d+00*sqrt(2.0d+00)
    rlenref=1.0d+10
    !
    !     Loop on faces
    !
    do iface=1,nface
       do j=1,nnofa
          ip1=lface(ltab(1,j),iface)
          ip2=lface(ltab(2,j),iface)
          !
          !     Check length
          !
          call length(ip1,ip2,rsize,coor,rlen)
          !
          !     Remember it for later
          !
          !write(*,*)iface,j,rlen
          if(abs(rlen-c10)<abs(rlenref-c10))then
             ip1ref=ip1
             ip2ref=ip2
          endif
          !
          !     We want an edge at least twice as big as the tolerance
          !
          if(rlen>tollen)then
             !
             !     Edge found, split it
             !      
             iview1=j
             ielem1=iface
             ielem2=eltoel(j,ielem1)

             ip1=lface(iview1,ielem1)
             ip2=lface(ltab(1,iview1),ielem1)
             ip3=lface(ltab(2,iview1),ielem1)

             if(eltoel(1,ielem2)==ielem1)then
                iview2=1_ip
             else if(eltoel(2,ielem2)==ielem1)then
                iview2=2_ip
             else
                iview2=3_ip
             endif

             ip4=lface(iview2,ielem2)

             ineigh11=eltoel(ltab(1,iview1),ielem1)
             ineigh12=eltoel(ltab(2,iview1),ielem1)
             ineigh21=eltoel(ltab(1,iview2),ielem2)
             ineigh22=eltoel(ltab(2,iview2),ielem2)
             !
             !     Allocate data base for ipnew
             !
             npnew=npoin+1
             ipnew=npnew
             call resizep(ipnew,coor,rnopo,lcart,lpsur,rsize,lptri,lmark,lpofa,lptype,lpsid)
             npoin=npnew
             !
             !     Update point database for ipnew 
             !
             vect(1)=coor(1,ip3)-coor(1,ip2)
             vect(2)=coor(2,ip3)-coor(2,ip2)
             vect(3)=coor(3,ip3)-coor(3,ip2)
             rl=sqrt(vect(1)*vect(1)+vect(2)*vect(2)+vect(3)*vect(3))
             rl=c10/rl
             vect(1)=rl*vect(1)
             vect(2)=rl*vect(2)
             vect(3)=rl*vect(3)

             pnew(1)=coor(1,ip2)+rsize(ip2)*vect(1)
             pnew(2)=coor(2,ip2)+rsize(ip2)*vect(2)
             pnew(3)=coor(3,ip2)+rsize(ip2)*vect(3)

             rsiz=c05*(rsize(ip2)+rsize(ip3))
             !
             !     Find the host face in the old mesh
             !
             call gthost(pnew,lpsur(ip2),ihost,rnofaold,ndim,nfold,npold,lfold,nnofa,&
                  coorold,d1,d2,d3,eltoelold,lelemold,lsurfold,isurf,ierr,rsiz)
             !
             !     Get the projected point
             !
             !dtot=d1+d2+d3
             !d1=d1/dtot  
             !d2=d2/dtot  
             !d3=d3/dtot  
             ip1o=lfold(1,ihost)
             ip2o=lfold(2,ihost)
             ip3o=lfold(3,ihost)

             coor(1,ipnew)=d1*coorold(1,ip1o)+d2*coorold(1,ip2o)+d3*coorold(1,ip3o)
             coor(2,ipnew)=d1*coorold(2,ip1o)+d2*coorold(2,ip2o)+d3*coorold(2,ip3o)
             coor(3,ipnew)=d1*coorold(3,ip1o)+d2*coorold(3,ip2o)+d3*coorold(3,ip3o)
             rnopo(1,ipnew)=rnofaold(1,ihost)
             rnopo(2,ipnew)=rnofaold(2,ihost)
             rnopo(3,ipnew)=rnofaold(3,ihost)
             icart=lcart(ip1)
             call gtelem(ipnew,coor,npoin,ndim,lcell,ncell,icart,rtol)
             lcart(ipnew)=icart
             call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipnew,ndim,coor,rsuni)
             lpsur(ipnew)=ihost 
             lptri(ipnew)=ielem1
             lmark(ipnew)=0_ip
             lpofa(ipnew)=0_ip
             lptype(1,ipnew)=ID_SMOOTH
             lptype(2,ipnew)=isurf
             !
             !     Allocate the 2 new elements
             !
             nfnew=nface+2
             call memrea(nfnew,memor_msh,'LFLOC','frstfr',lface)
             call memrea(nfnew,memor_msh,'RNOFA','frstfr',rnofa)
             !call memrea(nfnew,memor_msh,'LSURF','frstfr',lsurf)
             call memrea(nfnew,memor_msh,'ELTOEL','frstfr',eltoel)
             call memrea(nfnew,memor_msh,'LELEM','frstfr',lelem) 
             call memrea(nfnew,memor_msh,'LFMARK','frstfr',lfmark)

             nface=nface+1
             ielem3=nface
             nface=nface+1
             ielem4=nface

             !isurf1=lsurf(ielem1)
             !isurf2=lsurf(ielem2)

             !
             !     Create new elements
             ! 
             lface(1,ielem1)=ip1  
             lface(2,ielem1)=ip2  
             lface(3,ielem1)=ipnew
             eltoel(1,ielem1)=ielem2
             eltoel(2,ielem1)=ielem3
             eltoel(3,ielem1)=ineigh12
             !lsurf(ielem1)=isurf1
             call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,ielem1),ielem1,modul)

             lface(1,ielem2)=ip2  
             lface(2,ielem2)=ip4  
             lface(3,ielem2)=ipnew
             eltoel(1,ielem2)=ielem4
             eltoel(2,ielem2)=ielem1
             eltoel(3,ielem2)=ineigh21
             !lsurf(ielem2)=isurf2
             call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,ielem2),ielem2,modul)

             lface(1,ielem4)=ip4  
             lface(2,ielem4)=ip3  
             lface(3,ielem4)=ipnew
             eltoel(1,ielem4)=ielem3
             eltoel(2,ielem4)=ielem2
             eltoel(3,ielem4)=ineigh22
             !lsurf(ielem4)=isurf2
             lfmark(ielem4)=0_ip
             lelem(ielem4)=0_ip
             call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,ielem4),ielem4,modul)

             lface(1,ielem3)=ip3  
             lface(2,ielem3)=ip1  
             lface(3,ielem3)=ipnew
             eltoel(1,ielem3)=ielem1
             eltoel(2,ielem3)=ielem4
             eltoel(3,ielem3)=ineigh11
             !lsurf(ielem3)=isurf1
             lfmark(ielem3)=0_ip
             lelem(ielem3)=0_ip
             call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,ielem3),ielem3,modul)

             lptri(ip3)=ielem3


             !
             !     Outside
             !
             if(ineigh22/=0)then
                if(eltoel(1,ineigh22)==ielem2)then
                   eltoel(1,ineigh22)=ielem4
                else if(eltoel(2,ineigh22)==ielem2)then
                   eltoel(2,ineigh22)=ielem4
                else
                   eltoel(3,ineigh22)=ielem4
                endif
             endif

             if(ineigh11/=0)then 
                if(eltoel(1,ineigh11)==ielem1)then
                   eltoel(1,ineigh11)=ielem3
                else if(eltoel(2,ineigh11)==ielem1)then
                   eltoel(2,ineigh11)=ielem3
                else
                   eltoel(3,ineigh11)=ielem3
                endif
             endif
             !
             !     Resize the front
             !
             nfrontn=max(2,floor(sqrt(real(nface))))
             allocate(lfront(3,nfrontn),stat=istat)
             call memchk(zero,istat,memor_msh,'LFRONT','frstfr',lfront) 
             allocate(rfront(nfrontn),stat=istat)
             call memchk(zero,istat,memor_msh,'RFRONT','frstfr',rfront) 
             allocate(lheap(nfrontn),stat=istat)
             call memchk(zero,istat,memor_msh,'LHEAP','frstfr',lheap) 
             allocate(lfapo(3,nfrontn),stat=istat)
             call memchk(zero,istat,memor_msh,'LFAPO','frstfr',lfapo)
             allocate(lfahol(max(1_ip,nfrontn/10_ip)),stat=istat)
             call memchk(zero,istat,memor_msh,'LFAHOL','frstfr',lfahol) 
             allocate(lfhole(max(1_ip,nfrontn/10_ip)),stat=istat)
             call memchk(zero,istat,memor_msh,'LFHOLE','frstfr',lfhole) 
             nfront=2_ip
             lfront(1,1)=ip2
             lfront(2,1)=ipnew
             lfront(1,2)=ipnew
             lfront(2,2)=ip2
             lpofa(ip2)=1_ip
             lpofa(ipnew)=2_ip
             lfapo(1,1)=1_ip
             lfapo(2,1)=2_ip
             lfapo(3,1)=2_ip
             lfapo(1,2)=1_ip
             lfapo(2,2)=2_ip
             lfapo(3,2)=2_ip
             nfapo=2_ip
             nheap=2_ip 
             !
             !     Compute real length
             !
             rx=coor(1,ip2)-coor(1,ipnew)
             ry=coor(2,ip2)-coor(2,ipnew)
             rz=coor(3,ip2)-coor(3,ipnew)
             rlen=sqrt(rx*rx+ry*ry+rz*rz)
             rfront(1)=rlen
             rfront(2)=rlen
             lheap(1)=1_ip
             lheap(2)=2_ip

             return
          endif

       enddo
    enddo
    !
    !     Did not find an edge twice as big as the tolerance. 
    !     Select the nearest to 1 and moves ip1 to the optimal value
    !
    ip1=ip1ref
    ip2=ip2ref

    vect(1)=coor(1,ip1)-coor(1,ip2)
    vect(2)=coor(2,ip1)-coor(2,ip2)
    vect(3)=coor(3,ip1)-coor(3,ip2)
    rl=sqrt(vect(1)*vect(1)+vect(2)*vect(2)+vect(3)*vect(3))
    rl=c10/rl
    vect(1)=rl*vect(1)
    vect(2)=rl*vect(2)
    vect(3)=rl*vect(3)
    !
    !     Store the point at the end of the array
    !
    npnew=npoin+1
    call memrea(npnew,memor_msh,'COOR','frstfr',coor)
    coor(1,npnew)=coor(1,ip2)+rsize(ip2)*vect(1)
    coor(2,npnew)=coor(2,ip2)+rsize(ip2)*vect(2)
    coor(3,npnew)=coor(3,ip2)+rsize(ip2)*vect(3)

    call move(coor(:,npnew),coor,ndim,npoin,lface,nnofa,nface,ip1,rnopo,coorold,nfold,&
         npold,rnofaold,lfold,eltoelold,eltoel,lptri,lpsur,rnofa,lelemold,lsurfold,isurf,&
         rsize,lmark,lelem,lpofa,lfmark,lcart,lcell,ncell,rtol,rsuni,lptri(ip1))
    !
    !     Resize the front
    !
    nfrontn=max(2,floor(sqrt(real(nface))))
    allocate(lfront(3,nfrontn),stat=istat)
    call memchk(zero,istat,memor_msh,'LFRONT','frstfr',lfront) 
    allocate(rfront(nfrontn),stat=istat)
    call memchk(zero,istat,memor_msh,'RFRONT','frstfr',rfront) 
    allocate(lheap(nfrontn),stat=istat)
    call memchk(zero,istat,memor_msh,'LHEAP','frstfr',lheap) 
    call memchk(zero,istat,memor_msh,'LHEAP','frstfr',lheap) 
    allocate(lfapo(3,nfrontn),stat=istat)
    call memchk(zero,istat,memor_msh,'LFAPO','frstfr',lfapo)
    allocate(lfahol(max(1_ip,nfrontn/10_ip)),stat=istat)
    call memchk(zero,istat,memor_msh,'LFAHOL','frstfr',lfahol) 
    allocate(lfhole(max(1_ip,nfrontn/10_ip)),stat=istat)
    call memchk(zero,istat,memor_msh,'LFHOLE','frstfr',lfhole)   
    nfront=2_ip
    lfront(1,1)=ip2
    lfront(2,1)=ip1
    lfront(1,2)=ip1
    lfront(2,2)=ip2
    lpofa(ip2)=1_ip
    lpofa(ip1)=2_ip
    lfapo(1,1)=1_ip
    lfapo(2,1)=2_ip
    lfapo(3,1)=2_ip
    lfapo(1,2)=1_ip
    lfapo(2,2)=2_ip
    lfapo(3,2)=2_ip
    nfapo=2_ip
    nheap=2_ip 
    !
    !     Compute real length
    !
    rx=coor(1,ip2)-coor(1,ip1)
    ry=coor(2,ip2)-coor(2,ip1)
    rz=coor(3,ip2)-coor(3,ip1)
    rlen=sqrt(rx*rx+ry*ry+rz*rz)
    rfront(1)=rlen
    rfront(2)=rlen
    lheap(1)=1_ip
    lheap(2)=2_ip


  end subroutine frstfr

  subroutine bestpt(ip1,ip2,rsize,pnew,rnopo,ndim,npoin,coor)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)              :: ip1,ip2,ndim,npoin
    real(rp),intent(inout)              :: pnew(3) 
    real(rp),intent(in)                 :: rnopo(ndim,npoin),rsize(npoin),coor(ndim,npoin) 
    real(rp)                            :: rnx,rny,rnz,rtx,rty,rtz,rnl,rtl,rkx,rky,rkz,rkl,c10,c05,rnew
    real(rp)                            :: pmid(3) 

    c05=0.5d+00
    c10=1.0d+00
    !
    !     Interpolate the normal at the mid point
    !
    rnx=rnopo(1,ip1)+rnopo(1,ip2) 
    rny=rnopo(2,ip1)+rnopo(2,ip2) 
    rnz=rnopo(3,ip1)+rnopo(3,ip2) 
    rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
    rnl=c10/rnl
    rnx=rnx*rnl
    rny=rny*rnl
    rnz=rnz*rnl
    !
    !     Get the side vector
    !
    rtx=coor(1,ip2)-coor(1,ip1)
    rty=coor(2,ip2)-coor(2,ip1)
    rtz=coor(3,ip2)-coor(3,ip1)
    rtl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
    rtl=c10/rtl
    rtx=rtx*rtl
    rty=rty*rtl
    rtz=rtz*rtl
    !
    !     Get the vectorial product with the sides
    !
    rkx= rny*rtz-rnz*rty     
    rky=-rnx*rtz+rnz*rtx
    rkz= rnx*rty-rny*rtx
    rkl=sqrt(rkx*rkx+rky*rky+rkz*rkz)
    rkl=c10/rkl
    rkx=rkx*rkl
    rky=rky*rkl
    rkz=rkz*rkl 
    !
    !     Get the size
    !
    rnew=(rsize(ip1)+rsize(ip2))*c05
    !
    !     Get the mid point
    !
    pmid(1)=(coor(1,ip1)+coor(1,ip2))*c05
    pmid(2)=(coor(2,ip1)+coor(2,ip2))*c05
    pmid(3)=(coor(3,ip1)+coor(3,ip2))*c05

    !
    !     Create the new point
    !  
    pnew(1)= pmid(1)+rnew*rkx
    pnew(2)= pmid(2)+rnew*rky
    pnew(3)= pmid(3)+rnew*rkz

  end subroutine bestpt

  subroutine gthostsid(pnew,nside,nnosi,npoin,iguess,ihost,ndim,lside,sitosi,coor,d1,d2,&
       ptosi2old,ptosi1old,npnew,lelem,ierr,lline,iline)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nside,nnosi,npoin,iguess,ndim,npnew,iline
    integer(ip),intent(inout) :: ihost,lelem(nside),ierr
    real(rp),intent(inout)    :: d1,d2
    integer(ip),intent(in)    :: lside(nnosi,nside),sitosi(nnosi,nside)
    integer(ip),intent(in)    :: ptosi2old(npoin),ptosi1old(*),lline(nside)
    real(rp),intent(in)       :: coor(ndim,npoin),pnew(ndim)
    integer(ip)               :: maxiter,iter,ip1,ip2,istack,nstack,nstack2,iside,ipn,iclos
    integer(ip)               :: ineigh,inosi,ihostn,idir,ipoin,nsid,is,ipoint,js,jside,ihostf
    real(rp)                  :: rtx,rty,rtz,rnl,c10,rx,ry,rz,d1t,d2t
    real(rp)                  :: rtol,csca,cscal,cscat,cscalt
    real(rp)                  :: rdist,rdistt,rpx,rpy,rpz
    integer(ip), parameter    :: mstack=100_ip
    integer(ip)               :: lstack(mstack),lstack2(mstack)

    c10=1.0d+00
    rtol=-1.0d-05
    !
    !     Do we have a valid initial guess
    ! 
    if(iguess==0)then
       write(*,*)'Error in gthostsid, iguess=0'
       stop 
    endif

    ihostf=iguess 
    !
    !     Make sure iguess belongs to iline
    !
    if(lline(ihostf)/=iline)then

       nstack=1_ip
       istack=0_ip
       lstack(1)=ihostf
       lelem(ihostf)=1_ip 

       do 

          if(istack==nstack)then
             !
             !    Error
             !
             write(*,*)'Error in gthostsid, iline not found'
             stop
          endif
          istack=istack+1
          iside=lstack(istack) 

          do inosi=1,nnosi
             ip1=lside(inosi,iside)
             do is=ptosi2old(ip1),ptosi2old(ip1+1)-1
                iside=ptosi1old(is)
                if(lline(iside)==iline)then
                   ihostf=iside
                   goto 2
                endif
                if(lelem(iside)/=1)then
                   lelem(iside)=1
                   nstack=nstack+1
                   lstack(nstack)=iside
                endif

             enddo

          enddo
       enddo

2      continue

       do istack=1,nstack
          lelem(lstack(istack))=0_ip
       enddo

    endif
    !
    !     This subroutine interpolates a point created on a ridge
    !
    ihost=ihostf
    maxiter=500_ip
    do iter=1,maxiter
       !
       !     Do we have a valid host?
       !
       if(ihost==0)then
          ipn=lside(idir,ihostn)
          !
          !     Test the ridges around ipn
          !
          do is=ptosi2old(ipn),ptosi2old(ipn+1)-1 
             iside=ptosi1old(is)
             if(lline(iside)/=iline)cycle
             !
             !     Host points 
             !
             ip1=lside(1,iside) 
             ip2=lside(2,iside) 
             !
             !     Host tangent
             !
             rtx=coor(1,ip2)-coor(1,ip1)
             rty=coor(2,ip2)-coor(2,ip1)
             rtz=coor(3,ip2)-coor(3,ip1)
             rnl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
             rnl=c10/rnl 
             rtx=rnl*rtx  
             rty=rnl*rty  
             rtz=rnl*rtz
             !
             !     Project
             !
             rx=pnew(1)-coor(1,ip1) 
             ry=pnew(2)-coor(2,ip1) 
             rz=pnew(3)-coor(3,ip1) 
             csca=rtx*rx+rty*ry+rtz*rz
             cscal=csca*rnl
             d2=cscal
             d1=c10-cscal
             if(d1>rtol .and. d2>rtol)then
                ihost=iside
                goto 10
             endif
          enddo
          !
          !     Brute force 
          !
          call brutefside(lside,nnosi,nside,pnew,ihost,coor,ndim,npoin,lline,iline)
          goto 10

       else
          !
          !     Did we stay on the same line
          !
          if(lline(ihost)/=iline)exit
          !
          !     Host points 
          !
          ip1=lside(1,ihost) 
          ip2=lside(2,ihost) 
          !
          !     Host tangent
          !
          rtx=coor(1,ip2)-coor(1,ip1)
          rty=coor(2,ip2)-coor(2,ip1)
          rtz=coor(3,ip2)-coor(3,ip1)
          rnl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
          rnl=c10/rnl 
          rtx=rnl*rtx  
          rty=rnl*rty  
          rtz=rnl*rtz
          !
          !     Project
          !
          rx=pnew(1)-coor(1,ip1) 
          ry=pnew(2)-coor(2,ip1) 
          rz=pnew(3)-coor(3,ip1) 
          csca=rtx*rx+rty*ry+rtz*rz
          cscal=csca*rnl 
          d2=cscal
          d1=c10-cscal
          if(d2>rtol)then
             if(d1>rtol)then
                goto 10
             else
                ihostn=ihost
                ihost=sitosi(2,ihost)
                idir=2_ip 
             endif
          else
             ihostn=ihost
             ihost=sitosi(1,ihost)
             idir=1_ip 
          endif
       endif

    enddo
    !
    !     Element not found, test on distance
    !
    nstack=1_ip
    istack=0_ip
    lstack(1)=ihostf
    lelem(ihostf)=1_ip 
    nstack2=1_ip
    lstack2(1)=ihostf
    !
    !     Initialize distance from side ihostf to ipnew
    !     
    call distEdg(lside,nside,nnosi,ihostf,pnew,coor,ndim,npoin,rdist,d1,d2,ipoin)
    !
    !     Loop on stack
    !
    do 

       if(istack==nstack)exit
       istack=istack+1

       iside=lstack(istack)
       !
       !     Loop on neighbors
       ! 
       do inosi=1,nnosi
          ineigh=sitosi(inosi,iside) 
          !
          !     Do we have a valid host?
          !
          if(ineigh==0)then
             ipn=lside(inosi,iside)
             nsid=ptosi2old(ipn+1)-ptosi2old(ipn)
             !
             !     Test the ridges around corner ipn
             !
             do js=ptosi2old(ipn),ptosi2old(ipn+1)-1

                jside=ptosi1old(js)
                !
                !     Must pick the same line
                !
                if(lline(jside)/=iline)cycle
                !
                !     Has the side been already marked?
                ! 
                if(lelem(jside)==1)cycle
                lelem(jside)=1_ip
                if(nstack2==mstack)then
                   write(*,*)'Overflow in stack in gthostsid 1a for iside:',iside
                   stop
                endif
                nstack2=nstack2+1
                lstack2(nstack2)=jside
                !
                !     Compute distance
                !
                call distEdg(lside,nside,nnosi,jside,pnew,coor,ndim,npoin,rdistt,d1t,d2t,ipoint)

                if(ipoint/=0 .and. ipoint==ipoin)then

                   if(nstack==mstack)then
                      write(*,*)'Overflow in stack in gthostsid 1aa for iside:',iside
                      stop
                   endif
                   !
                   !     Add to the stack
                   ! 
                   nstack=nstack+1
                   lstack(nstack)=jside

                else if(rdistt<rdist)then
                   !
                   !     Check the shape function of the projected point
                   !
                   if(d1t>rtol .and. d2t>rtol)then
                      ihost=ineigh
                      d1=d1t
                      d2=d2t
                      ipoin=ipoint
                      goto 7
                   else
                      if(nstack==mstack)then
                         write(*,*)'Overflow in stack in gthostsid 1b for iside:',iside
                         stop
                      endif
                      nstack=nstack+1
                      lstack(nstack)=jside  
                      rdist=rdistt
                      ihost=ineigh
                      d2=d2t
                      d1=d1t
                      ipoin=ipoint
                   endif
                endif
             enddo

          else
             !
             !     Did we stay on the same line?
             !
             if(lline(ineigh)/=iline)cycle
             !
             !     Has the side been already marked?
             ! 
             if(lelem(ineigh)==1)cycle

             lelem(ineigh)=1_ip
             if(nstack2==mstack)then
                write(*,*)'Overflow in stack in gthostsid 1c for iside:',iside
                stop
             endif
             nstack2=nstack2+1_ip
             lstack2(nstack2)=ineigh
             !
             !     Compute distance
             !
             call distEdg(lside,nside,nnosi,ineigh,pnew,coor,ndim,npoin,rdistt,d1t,d2t,ipoint)

             if(ipoint/=0 .and. ipoint==ipoin)then

                if(nstack==mstack)then
                   write(*,*)'Overflow in stack in gthostsid 1cc for iside:',iside
                   stop
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
                if(d1t>rtol .and. d2t>rtol)then
                   ihost=ineigh
                   ipoin=ipoint
                   d2=d2t
                   d1=d1t
                   goto 7
                else
                   if(nstack==mstack)then
                      write(*,*)'Overflow in stack in gthostsid 1d for iside:',iside
                      stop
                   endif
                   nstack=nstack+1
                   lstack(nstack)=ineigh  
                   ihost=ineigh
                   rdist=rdistt
                   d2=d2t
                   d1=d1t
                   ipoin=ipoint
                endif
             endif
          endif
       enddo
    enddo
    !
    !     Brute force 
    !
    call brutefside(lside,nnosi,nside,pnew,ihost,coor,ndim,npoin,lline,iline)

7   continue
    !
    !     Clean up lelem
    !
    do istack=1,nstack2
       lelem(lstack2(istack))=0_ip
    enddo

10  continue  
    !
    !     Try to find better
    !
    istack=0_ip
    nstack=1_ip
    lstack(1)=ihost
    lstack2(1)=ihost
    nstack2=1_ip
    lelem(ihost)=1_ip 
    iclos=0_ip
    !
    !     Initialize distance from side ihost to ipnew
    !     
    call distEdg(lside,nside,nnosi,ihost,pnew,coor,ndim,npoin,rdist,d1,d2,ipoin)
    !
    !     Loop on the stack
    !
    do 
       if(istack==nstack)exit
       istack=istack+1

       iside=lstack(istack)
       do inosi=1,nnosi
          ineigh=sitosi(inosi,iside)
          !
          !     We may have ineigh==0 for cusp and corner points
          !
          if(ineigh==0)then

             ipn=lside(inosi,iside)
             !
             !     Test the ridges around corner ipn
             !
             do js=ptosi2old(ipn),ptosi2old(ipn+1)-1 
                jside=ptosi1old(js)
                !
                !     Must pick the same line
                !
                if(lline(jside)/=iline)cycle
                !
                !     Has the side been already marked?
                ! 
                if(lelem(jside)==1)cycle
                lelem(jside)=1_ip
                if(nstack2==mstack)then
                   write(*,*)'Overflow in stack in gthostsid 2a for iside:',iside
                   stop
                endif
                nstack2=nstack2+1
                lstack2(nstack2)=jside
                !
                !     Compute distance
                !
                call distEdg(lside,nside,nnosi,jside,pnew,coor,ndim,npoin,rdistt,d1t,d2t,ipoint)

                if(ipoint/=0 .and. ipoint==ipoin)then

                   if(nstack==mstack)then
                      write(*,*)'Overflow in stack in gthostsid 2aa for iside:',iside
                      stop
                   endif
                   !
                   !     Add to the stack
                   ! 
                   nstack=nstack+1
                   lstack(nstack)=jside

                else if(rdistt<rdist)then
                   !
                   !     Check the shape function of the projected point
                   !
                   if(d1t>rtol .and. d2t>rtol)then
                      ihost=ineigh
                      d1=d1t
                      d2=d2t
                      ipoin=ipoint
                      rdist=rdistt
                   else
                      if(nstack==mstack)then
                         write(*,*)'Overflow in stack in gthostsid 2b for iside:',iside
                         stop
                      endif
                      nstack=nstack+1
                      lstack(nstack)=jside  
                      rdist=rdistt
                      ihost=ineigh
                      d2=d2t
                      d1=d1t
                      ipoin=ipoint
                   endif
                endif

             enddo

          else
             !
             !     Did we stay on the same line?
             !        
             if(lline(ineigh)/=iline)cycle
             !
             !     Has the side been already marked?
             ! 
             if(lelem(ineigh)==1)cycle
             lelem(ineigh)=1_ip
             if(nstack2==mstack)then
                write(*,*)'Overflow in stack in gthostsid 2c for iside:',iside
                stop
             endif
             nstack2=nstack2+1_ip
             lstack2(nstack2)=ineigh
             !
             !     Compute distance
             !
             call distEdg(lside,nside,nnosi,ineigh,pnew,coor,ndim,npoin,rdistt,d1t,d2t,ipoint)

             if(ipoint/=0 .and. ipoint==ipoin)then

                if(nstack==mstack)then
                   write(*,*)'Overflow in stack in gthostsid for 2cc iside:',iside
                   stop
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
                if(d1>rtol .and. d2>rtol)then
                   ihost=ineigh
                   rdist=rdistt
                   d2=d2t
                   d1=d1t
                   ipoin=ipoint
                   iclos=0_ip

                else
                   if(nstack==mstack)then
                      write(*,*)'Overflow in stack in gthostsid 2d for iside:',iside
                      stop
                   endif
                   nstack=nstack+1
                   lstack(nstack)=ineigh  
                   ihost=ineigh
                   rdist=rdistt
                   d2=d2t
                   d1=d1t
                   ipoin=ipoint
                   iclos=1
                endif
             endif
          endif
       enddo
    enddo

20  continue
    !
    !     Clean up lelem
    !
    do istack=1,nstack2
       lelem(lstack2(istack))=0_ip
    enddo
    !
    !     Did we find a closer side not containing pnew?
    !
    if(iclos==1)then 
       call brutefside(lside,nnosi,nside,pnew,ihost,coor,ndim,npoin,lline,iline)
    endif

  end subroutine gthostsid

  subroutine gtclos(ihost,ipnew,lface,nface,nnofa,coor,ndim,npoin,rsize,ipclos,lmark,&
       lelem,eltoel,lpofa,lfmark,ipa,ipb)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nface,ndim,npoin,ihost,nnofa,ipnew,ipa,ipb
    integer(ip),intent(in)          :: lface(nnofa,nface),eltoel(nnofa,nface),lpofa(npoin),lfmark(nface)
    real(rp),intent(in)             :: coor(ndim,npoin),rsize(npoin)
    integer(ip),intent(inout)       :: ipclos,lmark(npoin),lelem(nface) 
    integer(ip)                     :: lstack(100),nstack,istack,ielem,lstackp(100),nstackp,ichk,istackp,nstackp0
    integer(ip)                     :: lstackp2(100),nstackp2
    integer(ip)                     :: imin,ip1,ip2,ip3,ineigh,iface,ipoin 
    real(rp)                        :: rlen,rstackp(100),rmin,vmin(3),vmax(3)
    real(rp)                        :: xmin,ymin,zmin,xmax,ymax,zmax,tollen
    !
    !     This subroutine checks if some point is too close to ipnew 
    !     ipa and ipb are not taken into account as they belong to the generating edge
    !
    tollen=1.0d+00/sqrt(2.0d+00)

    !
    !     DBG
    ! 
    !do ipoin=1,npoin
    !   if(lmark(ipoin)/=0)then
    !      write(*,*)'error in gtclos, lmark not clean at ipoin=',ipoin
    !      stop
    !   endif
    !enddo

    !do iface=1,nface
    !   if(lelem(iface)/=0)then
    !      write(*,*)'error in gtclos, lelem not clean at iface=',iface
    !      stop
    !   endif
    !enddo

    istack=0_ip
    nstack=1_ip
    lstack(1)=ihost
    lelem(ihost)=1_ip
    nstackp=0_ip
    nstackp2=0_ip
    !
    !     Compute the local box
    !
    rlen=rsize(ipnew)
    vmin(1)=coor(1,ipnew)-rlen
    vmin(2)=coor(2,ipnew)-rlen
    vmin(3)=coor(3,ipnew)-rlen
    vmax(1)=coor(1,ipnew)+rlen
    vmax(2)=coor(2,ipnew)+rlen
    vmax(3)=coor(3,ipnew)+rlen
    !write(*,*)vmin(1),vmin(2),vmin(3)
    !write(*,*)vmax(1),vmax(2),vmax(3)
    !
    !     Find the points in the box
    !
    do                  
       if(istack==nstack)exit
       istack=istack+1
       ielem=lstack(istack) 
       !
       !     Get the box of the element
       ! 

       ip1=lface(1,ielem)
       ip2=lface(2,ielem)
       ip3=lface(3,ielem)

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
       !write(*,*)xmin,ymin,zmin
       !write(*,*)xmax,ymax,zmax

       if(vmin(1)>xmax)cycle
       if(vmin(2)>ymax)cycle
       if(vmin(3)>zmax)cycle
       if(vmax(1)<xmin)cycle
       if(vmax(2)<ymin)cycle
       if(vmax(3)<zmin)cycle


       ineigh=eltoel(1,ielem)
       if(ineigh/=0)then
          if(lfmark(ineigh)/=1)then
             if(lelem(ineigh)==0)then
                lelem(ineigh)=1_ip
                nstack=nstack+1
                lstack(nstack)=ineigh
             endif
          endif
       endif
       ineigh=eltoel(2,ielem)
       if(ineigh/=0)then
          if(lfmark(ineigh)/=1)then
             if(lelem(ineigh)==0)then
                lelem(ineigh)=1_ip
                nstack=nstack+1
                lstack(nstack)=ineigh
             endif
          endif
       endif
       ineigh=eltoel(3,ielem)
       if(ineigh/=0)then
          if(lfmark(ineigh)/=1)then
             if(lelem(ineigh)==0)then
                lelem(ineigh)=1_ip
                nstack=nstack+1
                lstack(nstack)=ineigh
             endif
          endif
       endif
       !
       !     Check the points of the element     
       !
       if(lmark(ip1)==0)then
          lmark(ip1)=1_ip
          nstackp2=nstackp2+1
          lstackp2(nstackp2)=ip1
          if(lpofa(ip1)/=-1)then 
             call length(ip1,ipnew,rsize,coor,rlen)
             if(rlen<tollen)then
                nstackp=nstackp+1
                lstackp(nstackp)=ip1
                rstackp(nstackp)=rlen
             endif
          endif
       endif

       if(lmark(ip2)==0)then
          lmark(ip2)=1_ip
          nstackp2=nstackp2+1
          lstackp2(nstackp2)=ip2
          if(lpofa(ip2)/=-1)then 
             call length(ip2,ipnew,rsize,coor,rlen)
             if(rlen<tollen)then
                nstackp=nstackp+1
                lstackp(nstackp)=ip2
                rstackp(nstackp)=rlen
             endif
          endif
       endif

       if(lmark(ip3)==0)then
          lmark(ip3)=1_ip
          nstackp2=nstackp2+1
          lstackp2(nstackp2)=ip3
          if(lpofa(ip3)/=-1)then 
             call length(ip3,ipnew,rsize,coor,rlen)
             if(rlen<tollen)then
                nstackp=nstackp+1
                lstackp(nstackp)=ip3
                rstackp(nstackp)=rlen
             endif
          endif
       endif

    enddo
    !
    !     Clean up lelem
    !
    do istack=1,nstack
       lelem(lstack(istack))=0_ip
    enddo
    !
    !     Clean up lmark
    !
    do istackp=1,nstackp2
       lmark(lstackp2(istackp))=0_ip
    enddo
    !
    !     Do we have close points?
    !
    if(nstackp==0)then
       ipclos=0 
       return
    endif
    !
    !     Take out ipa and ipb
    ! 
    nstackp0=nstackp
    nstackp=0_ip
    do istackp=1,nstackp
       ipoin=lstackp(istackp)
       if(ipoin/=ipa .and. ipoin/=ipb)then
          nstackp=nstackp+1
          lstackp(nstackp)=ipoin
       endif
    enddo
    !
    !     Do we still have close points?
    !
    if(nstackp==0)then
       ipclos=0 
       return
    endif
    !
    !     Take the minimum
    !
    rmin=rstackp(1)
    imin=1_ip

    do istackp=2,nstackp
       if(rstackp(istackp)<rmin)then
          rmin=rstackp(istackp)
          imin=istackp
       endif
    enddo

    ipclos=lstackp(imin)

    !
    !     DBG
    ! 
    !do ipoin=1,npoin
    !   if(lmark(ipoin)/=0)then
    !      write(*,*)'error in gtclos 2, lmark not clean at ipoin=',ipoin
    !      stop
    !   endif
    !enddo

    !do iface=1,nface
    !   if(lelem(iface)/=0)then
    !      write(*,*)'error in gtclos , lelem not clean at iface=',iface
    !      stop
    !   endif
    !enddo

  end subroutine gtclos

  subroutine gtclos2(ihost,ipnew,lface,nface,nnofa,coor,ndim,npoin,rsize,lmark,&
       lelem,eltoel,lpofa,lfmark,lstackp,nstackp,mstackp)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nface,ndim,npoin,ihost,nnofa,ipnew,mstackp
    integer(ip),intent(in)          :: lface(nnofa,nface),eltoel(nnofa,nface),lpofa(npoin),lfmark(nface)
    real(rp),intent(in)             :: coor(ndim,npoin),rsize(npoin)
    integer(ip),intent(inout)       :: lmark(npoin),lelem(nface),lstackp(mstackp),nstackp 
    integer(ip)                     :: lstack(100),nstack,istack,ielem,ichk,istackp,nstackp0
    integer(ip)                     :: lstackp2(100),nstackp2
    integer(ip)                     :: imin,ip1,ip2,ip3,ineigh,iface,ipoin 
    real(rp)                        :: rlen,rstackp(100),rmin,vmin(3),vmax(3)
    real(rp)                        :: xmin,ymin,zmin,xmax,ymax,zmax,tollen,c20
    !
    !     This subroutine checks if some point is too close to ipnew 
    !
    tollen=1.0d+00/sqrt(2.0d+00)
    c20=2.0d+00
    !
    !     DBG
    ! 
    !do ipoin=1,npoin
    !   if(lmark(ipoin)/=0)then
    !      write(*,*)'error in gtclos, lmark not clean at ipoin=',ipoin
    !      stop
    !   endif
    !enddo

    !do iface=1,nface
    !   if(lelem(iface)/=0)then
    !      write(*,*)'error in gtclos, lelem not clean at iface=',iface
    !      stop
    !   endif
    !enddo

    istack=0_ip
    nstack=1_ip
    lstack(1)=ihost
    lelem(ihost)=1_ip
    nstackp=0_ip
    nstackp2=0_ip
    !
    !     Compute the local box
    !
    rlen=c20*rsize(ipnew)
    vmin(1)=coor(1,ipnew)-rlen
    vmin(2)=coor(2,ipnew)-rlen
    vmin(3)=coor(3,ipnew)-rlen
    vmax(1)=coor(1,ipnew)+rlen
    vmax(2)=coor(2,ipnew)+rlen
    vmax(3)=coor(3,ipnew)+rlen
    !write(*,*)vmin(1),vmin(2),vmin(3)
    !write(*,*)vmax(1),vmax(2),vmax(3)
    !
    !     Find the points in the box
    !
    do                  
       if(istack==nstack)exit
       istack=istack+1
       ielem=lstack(istack) 
       !
       !     Get the box of the element
       ! 

       ip1=lface(1,ielem)
       ip2=lface(2,ielem)
       ip3=lface(3,ielem)

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
       !write(*,*)xmin,ymin,zmin
       !write(*,*)xmax,ymax,zmax

       if(vmin(1)>xmax)cycle
       if(vmin(2)>ymax)cycle
       if(vmin(3)>zmax)cycle
       if(vmax(1)<xmin)cycle
       if(vmax(2)<ymin)cycle
       if(vmax(3)<zmin)cycle


       ineigh=eltoel(1,ielem)
       if(ineigh/=0)then
          if(lfmark(ineigh)/=1)then
             if(lelem(ineigh)==0)then
                lelem(ineigh)=1_ip
                nstack=nstack+1
                lstack(nstack)=ineigh
             endif
          endif
       endif
       ineigh=eltoel(2,ielem)
       if(ineigh/=0)then
          if(lfmark(ineigh)/=1)then
             if(lelem(ineigh)==0)then
                lelem(ineigh)=1_ip
                nstack=nstack+1
                lstack(nstack)=ineigh
             endif
          endif
       endif
       ineigh=eltoel(3,ielem)
       if(ineigh/=0)then
          if(lfmark(ineigh)/=1)then
             if(lelem(ineigh)==0)then
                lelem(ineigh)=1_ip
                nstack=nstack+1
                lstack(nstack)=ineigh
             endif
          endif
       endif
       !
       !     Check the points of the element     
       !
       if(lmark(ip1)==0)then
          lmark(ip1)=1_ip
          nstackp2=nstackp2+1
          lstackp2(nstackp2)=ip1
          if(lpofa(ip1)>0)then 
             nstackp=nstackp+1
             lstackp(nstackp)=ip1
          endif
       endif

       if(lmark(ip2)==0)then
          lmark(ip2)=1_ip
          nstackp2=nstackp2+1
          lstackp2(nstackp2)=ip2
          if(lpofa(ip2)>0)then 
             nstackp=nstackp+1
             lstackp(nstackp)=ip2
          endif
       endif

       if(lmark(ip3)==0)then
          lmark(ip3)=1_ip
          nstackp2=nstackp2+1
          lstackp2(nstackp2)=ip3
          if(lpofa(ip3)>0)then 
             nstackp=nstackp+1
             lstackp(nstackp)=ip3
          endif
       endif

    enddo
    !
    !     Clean up lelem
    !
    do istack=1,nstack
       lelem(lstack(istack))=0_ip
    enddo
    !
    !     Clean up lmark
    !
    do istackp=1,nstackp2
       lmark(lstackp2(istackp))=0_ip
    enddo
    !
    !     Remove ipnew
    !
    nstackp0=nstackp
    nstackp=0_ip
    do istackp=1,nstackp0 
       ipoin=lstackp(istackp)
       if(ipoin/=ipnew)then
          nstackp=nstackp+1
          lstackp(nstackp)=ipoin
       endif
    enddo
    !
    !     DBG
    ! 
    !do ipoin=1,npoin
    !   if(lmark(ipoin)/=0)then
    !      write(*,*)'error in gtclos 2, lmark not clean at ipoin=',ipoin
    !      stop
    !   endif
    !enddo

    !do iface=1,nface
    !   if(lelem(iface)/=0)then
    !      write(*,*)'error in gtclos , lelem not clean at iface=',iface
    !      stop
    !   endif
    !enddo

  end subroutine gtclos2

  subroutine gtclos3(ihost,ipnew,lface,nface,nnofa,coor,ndim,npoin,rsize,lelem,eltoel,ipclos,tolsplit)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nface,ndim,npoin,ihost,nnofa,ipnew
    integer(ip),intent(in)          :: lface(nnofa,nface),eltoel(nnofa,nface)
    real(rp),intent(in)             :: coor(ndim,npoin),rsize(npoin),tolsplit
    integer(ip),intent(inout)       :: lelem(nface),ipclos
    integer(ip),parameter           :: mstack=1000
    integer(ip)                     :: lstack(mstack),nstack,istack,ielem,ichk,istackp,nstackp0
    integer(ip)                     :: imin,ip1,ip2,ip3,ineigh,iface,ipoin 
    real(rp)                        :: rlen,rmin,vmin(3),vmax(3),rnl
    real(rp)                        :: xmin,ymin,zmin,xmax,ymax,zmax,tollen,c20
    !
    !     This subroutine checks if some point is too close to ipnew 
    !
    tollen=1.0d+00/sqrt(2.0d+00)
    c20=2.0d+00

    ipclos=0_ip 
    istack=0_ip
    nstack=1_ip
    lstack(1)=ihost
    lelem(ihost)=1_ip
    !
    !     Compute the local box
    !
    rlen=c20*rsize(ipnew)
    vmin(1)=coor(1,ipnew)-rlen
    vmin(2)=coor(2,ipnew)-rlen
    vmin(3)=coor(3,ipnew)-rlen
    vmax(1)=coor(1,ipnew)+rlen
    vmax(2)=coor(2,ipnew)+rlen
    vmax(3)=coor(3,ipnew)+rlen
    !write(*,*)vmin(1),vmin(2),vmin(3)
    !write(*,*)vmax(1),vmax(2),vmax(3)
    !
    !     Find the points in the box
    !
    do                  
       if(istack==nstack)exit
       istack=istack+1
       ielem=lstack(istack) 
       !
       !     Get the box of the element
       ! 

       ip1=lface(1,ielem)
       ip2=lface(2,ielem)
       ip3=lface(3,ielem)

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


       ineigh=eltoel(1,ielem)
       if(ineigh/=0)then
          if(lelem(ineigh)==0)then
             lelem(ineigh)=1_ip
             if(nstack==mstack)then
                write(*,*)'Error in gtclos3, nstack==mstack 1'
                stop
             endif
             nstack=nstack+1
             lstack(nstack)=ineigh
          endif
       endif
       ineigh=eltoel(2,ielem)
       if(ineigh/=0)then
          if(lelem(ineigh)==0)then
             lelem(ineigh)=1_ip
             if(nstack==mstack)then
                write(*,*)'Error in gtclos3, nstack==mstack 2'
                stop
             endif
             nstack=nstack+1
             lstack(nstack)=ineigh
          endif
       endif
       ineigh=eltoel(3,ielem)
       if(ineigh/=0)then
          if(lelem(ineigh)==0)then
             lelem(ineigh)=1_ip
             if(nstack==mstack)then
                write(*,*)'Error in gtclos3, nstack==mstack 3'
                stop
             endif
             nstack=nstack+1
             lstack(nstack)=ineigh
          endif
       endif
       !
       !     Check the points of the element     
       !
       if(ip1/=ipnew)then
          call length(ip1,ipnew,rsize,coor,rnl)
          if(rnl<tolsplit)then
             ipclos=ip1
             exit
          endif
       endif
       if(ip2/=ipnew)then
          call length(ip2,ipnew,rsize,coor,rnl)
          if(rnl<tolsplit)then
             ipclos=ip2
             exit
          endif
       endif
       if(ip3/=ipnew)then
          call length(ip3,ipnew,rsize,coor,rnl)
          if(rnl<tolsplit)then
             ipclos=ip3
             exit
          endif
       endif

    enddo
    !
    !     Clean up lelem
    !
    do istack=1,nstack
       lelem(lstack(istack))=0_ip
    enddo

  end subroutine gtclos3

  subroutine insert(nface,ihost,nnofa,ndim,npoin,d1,d2,d3,ipnew,coor,&
       lface,eltoel,rnofa,lptri,lelem,lfmark,rnopo,lehole,nehole)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(inout)       :: nface,nehole
    integer(ip),intent(inout)       :: lehole(nehole)
    integer(ip),intent(in)          :: ndim,npoin,ihost,nnofa,ipnew
    real(rp),intent(in)             :: d1,d2,d3
    real(rp),intent(in)             :: coor(ndim,npoin),rnopo(ndim,npoin)
    integer(ip),pointer             :: lface(:,:),eltoel(:,:),lptri(:)
    integer(ip),pointer             :: lelem(:),lfmark(:)
    real(rp),pointer                :: rnofa(:,:)
    real(rp)                        :: epsil,modul
    integer(ip)                     :: ip1,ip2,ip3,neigh1,neigh2,nfac1,nfnew
    integer(ip)                     :: ineigh1,ineigh2,ineigh3,ierr,ielem2,ielem3
    epsil=1.0d-01
    ierr=0_ip 
    !
    !     Did we cross the front
    !
    if(lfmark(ihost)==1)then
       write(*,*)'Error insert, face already marked',ihost 
       stop
    endif
    !
    !     Get host face
    !      
    ip1=lface(1,ihost)
    ip2=lface(2,ihost)
    ip3=lface(3,ihost)
    !
    !     Get the neighbors
    !
    ineigh1=eltoel(1,ihost) 
    ineigh2=eltoel(2,ihost) 
    ineigh3=eltoel(3,ihost) 
    !
    !     Check the simplicial coordinates
    !
    if(d1<epsil)then
       if(ineigh1/=0)then
          if(lfmark(ineigh1)==0)then
             !
             !     Split on first side
             !  
             call split2(ihost,1_ip,nnofa,nface,npoin,ndim,ipnew,lptri,&
                  lface,eltoel,lfmark,rnofa,lelem,coor,rnopo,ierr,nehole,lehole)
             !
             !     The point  was inserted --> go home
             !
             if(ierr==0)then
                return 
             endif
          endif
       endif

    else if(d2<epsil)then
       if(ineigh2/=0)then
          if(lfmark(ineigh2)==0)then 
             !
             !     Split on second side
             !  
             call split2(ihost,2_ip,nnofa,nface,npoin,ndim,ipnew,lptri,&
                  lface,eltoel,lfmark,rnofa,lelem,coor,rnopo,ierr,nehole,lehole)
             !
             !     The point  was inserted --> go home
             !
             if(ierr==0)then
                return 
             endif
          endif
       endif

    else if(d3<epsil)then
       if(ineigh3/=0)then
          if(lfmark(ineigh3)==0)then 
             !
             !     Split on third side
             !  
             call split2(ihost,3_ip,nnofa,nface,npoin,ndim,ipnew,lptri,&
                  lface,eltoel,lfmark,rnofa,lelem,coor,rnopo,ierr,nehole,lehole)
             !
             !     The point  was inserted --> go home
             !
             if(ierr==0)then
                return 
             endif
          endif
       endif
    endif
    !
    !     It is possible that some of the shape function be below epsil,
    !     but the point is in front of a good element because the old element
    !     in front is very large. Therefore, the element in front must be split
    !     in three elements.  
    !

    !
    !     Insert npoin in the surface mesh 
    !   
    lface(3,ihost)=ipnew
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,ihost),ihost,modul)
    !
    !     Realloc lface 
    !
    if(nehole>0)then
       ielem2=lehole(nehole) 
       nehole=nehole-1
    else 
       ielem2=nface+1
       call memrea(ielem2,memor_msh,'LFACE','insert',lface)
       call memrea(ielem2,memor_msh,'RNOFA','insert',rnofa)
       call memrea(ielem2,memor_msh,'ELTOEL','insert',eltoel)
       call memrea(ielem2,memor_msh,'LELEM','insert',lelem) 
       call memrea(ielem2,memor_msh,'LFMARK','insert',lfmark) 
       nface=ielem2
    endif

     if(nehole>0)then
       ielem3=lehole(nehole) 
       nehole=nehole-1
    else 
       ielem3=nface+1
       call memrea(ielem3,memor_msh,'LFACE','insert',lface)
       call memrea(ielem3,memor_msh,'RNOFA','insert',rnofa)
       call memrea(ielem3,memor_msh,'ELTOEL','insert',eltoel)
       call memrea(ielem3,memor_msh,'LELEM','insert',lelem) 
       call memrea(ielem3,memor_msh,'LFMARK','insert',lfmark) 
       nface=ielem3
    endif

    lface(1,ielem2)=ip2
    lface(2,ielem2)=ip3
    lface(3,ielem2)=ipnew
    lfmark(ielem2)=0_ip
    lelem(ielem2)=0_ip
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,ielem2),ielem2,modul)

    lface(1,ielem3)=ip3
    lface(2,ielem3)=ip1
    lface(3,ielem3)=ipnew
    lfmark(ielem3)=0_ip
    lelem(ielem3)=0_ip
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,ielem3),ielem3,modul)
    !
    !     Update eltoel
    !
    neigh1=eltoel(1,ihost)
    neigh2=eltoel(2,ihost)
    eltoel(1,ihost)=ielem2
    eltoel(2,ihost)=ielem3

    eltoel(1,ielem2)=ielem3
    eltoel(2,ielem2)=ihost
    eltoel(3,ielem2)=neigh1

    eltoel(1,ielem3)=ihost
    eltoel(2,ielem3)=ielem2
    eltoel(3,ielem3)=neigh2
    !
    !     Update eltoel outside
    !
    if(neigh1/=0)then
       if(eltoel(1,neigh1)==ihost)then
          eltoel(1,neigh1)=ielem2
       else if(eltoel(2,neigh1)==ihost)then
          eltoel(2,neigh1)=ielem2
       else
          eltoel(3,neigh1)=ielem2
       endif
    endif

    if(neigh2/=0)then
       if(eltoel(1,neigh2)==ihost)then
          eltoel(1,neigh2)=ielem3
       else if(eltoel(2,neigh2)==ihost)then
          eltoel(2,neigh2)=ielem3
       else
          eltoel(3,neigh2)=ielem3
       endif
    endif

    lptri(ipnew)=ihost
    lptri(ip3)=ielem2

  end subroutine insert

  subroutine move(ptarget,coor,ndim,npoin,lface,nnofa,nface,ipclos,rnopo,coorold,nfold,&
       npold,rnofaold,lfold,eltoelold,eltoel,lptri,lpsur,rnofa,lelemold,lsurfold,isurf,rsize,&
       lmark,lelem,lpofa,lfmark,lcart,lcell,ncell,rtol,rsuni,ihostn)
    use def_kintyp, only       :  ip,rp,lg,cell
    implicit none
    integer(ip),intent(in)          :: nface,ndim,npoin,nnofa,ipclos,isurf
    integer(ip),intent(in)          :: npold,nfold,ncell,ihostn
    integer(ip),intent(in)          :: lface(nnofa,nface),lptri(npoin),lpofa(npoin)
    integer(ip),intent(in)          :: lfold(nnofa,nfold),eltoel(nnofa,nface),lfmark(nface)
    integer(ip),intent(in)          :: eltoelold(nnofa,nfold)
    integer(ip),intent(inout)       :: lpsur(npoin),lelemold(nfold),lmark(npoin),lelem(nface),lcart(npoin)
    integer(ip),intent(in)          :: lsurfold(nfold)
    type(cell),intent(in)           :: lcell(ncell)
    real(rp),intent(inout)          :: coor(ndim,npoin),rnopo(ndim,npoin),rnofa(ndim,nface)
    real(rp),intent(inout)          :: rsize(npoin)
    real(rp),intent(in)             :: coorold(ndim,npold),ptarget(ndim)
    real(rp),intent(in)             :: rnofaold(ndim,nfold),rtol,rsuni
    integer(ip)                     :: ip1,ip2,ip3,iefirst,ienext,ipsur,ierr,lpclos(500),npclos
    integer(ip)                     :: ie,iter,maxiter,lball(100),nball,iface,lfacn(nnofa,100),ielem,ihost 
    integer(ip)                     :: icart,ipo 
    integer(ip)                     :: mpclos=500 
    real(rp)                        :: pold(ndim),rnopoold(ndim),c10,wei1,wei2,Qgeo,Qgeotol,rface(ndim),pnew(ndim) 
    real(rp)                        :: d1,d2,d3,dtot,rpoin(ndim),rfloc(ndim,100),modul(100),c00,rlen,tollen,rsiz
    !
    !     This subroutine tries to move iteratively ipclos to ptarget 
    !
    c10=1.0d+00
    c00=0.0d+00
    Qgeotol=0.0d+00
    tollen=1.0d+00/sqrt(2.0d+00)
    !
    !     At this point, we have a good position but an already existing point is too close 
    !

    !
    !     Get the list of close points
    !  
    call gtclos2(ihostn,ipclos,lface,nface,nnofa,coor,ndim,npoin,rsize,lmark,&
         lelem,eltoel,lpofa,lfmark,lpclos,npclos,mpclos)
    !
    !     Fix local size
    !
    rsiz=rsize(ipclos)
    !
    !     Get the ball of the point
    !
    ie=lptri(ipclos)
    lfacn(1,1)=lface(1,ie)
    lfacn(2,1)=lface(2,ie)
    lfacn(3,1)=lface(3,ie)
    iefirst=ie
    nball=1_ip
    lball(1)=ie


    do 

       if(lface(1,ie)==ipclos)then
          ienext=eltoel(2,ie)
       else if(lface(2,ie)==ipclos)then
          ienext=eltoel(3,ie)
       else
          ienext=eltoel(1,ie)
       endif

       if(ienext==iefirst)exit

       nball=nball+1
       lball(nball)=ienext
       lfacn(1,nball)=lface(1,ienext)
       lfacn(2,nball)=lface(2,ienext)
       lfacn(3,nball)=lface(3,ienext)
       ie=ienext

    enddo
    !
    !     Remember the position, normal and lpsur of ipclos
    !
    pold(1)=coor(1,ipclos)
    pold(2)=coor(2,ipclos)
    pold(3)=coor(3,ipclos)
    rnopoold(1)=rnopo(1,ipclos)
    rnopoold(2)=rnopo(2,ipclos)
    rnopoold(3)=rnopo(3,ipclos)
    ipsur=lpsur(ipclos)
    !
    !     Loop iteratively
    !   
    maxiter=5_ip

    do iter=1,maxiter
       !
       !     Compute new position
       !
       wei1=(iter-1)/maxiter
       wei2=c10-wei1

       pnew(1)=wei1*pold(1)+wei2*ptarget(1)
       pnew(2)=wei1*pold(2)+wei2*ptarget(2)
       pnew(3)=wei1*pold(3)+wei2*ptarget(3)
       !
       !     Find the host face in the old mesh
       !
       call gthost(pnew,lpsur(ipclos),ihost,rnofaold,ndim,nfold,npold,lfold,nnofa,&
            coorold,d1,d2,d3,eltoelold,lelemold,lsurfold,isurf,ierr,rsiz)
       lpsur(ipclos)=ihost
       !
       !     Get the projected point
       !
       !dtot=d1+d2+d3
       !d1=d1/dtot  
       !d2=d2/dtot  
       !d3=d3/dtot  
       ip1=lfold(1,ihost)
       ip2=lfold(2,ihost)
       ip3=lfold(3,ihost)

       pnew(1)=d1*coorold(1,ip1)+d2*coorold(1,ip2)+d3*coorold(1,ip3)
       pnew(2)=d1*coorold(2,ip1)+d2*coorold(2,ip2)+d3*coorold(2,ip3)
       pnew(3)=d1*coorold(3,ip1)+d2*coorold(3,ip2)+d3*coorold(3,ip3)
       rpoin(1)=rnofaold(1,ihost)
       rpoin(2)=rnofaold(2,ihost)
       rpoin(3)=rnofaold(3,ihost)
       !
       !     Transfer data to ipclos
       !
       icart=lcart(ipclos)
       call gtelem(ipclos,coor,npoin,ndim,lcell,ncell,icart,rtol)
       lcart(ipclos)=icart
       call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipclos,ndim,coor,rsuni)
       coor(1,ipclos)=pnew(1)
       coor(2,ipclos)=pnew(2)
       coor(3,ipclos)=pnew(3)
       rnopo(1,ipclos)=rpoin(1)
       rnopo(2,ipclos)=rpoin(2)
       rnopo(3,ipclos)=rpoin(3)
       !
       !     Try to move ipopt to ptarget as much as possible
       !
       do ie=1,nball
          call gtfnr2(lfacn,nball,nnofa,ndim,coor,npoin,rfloc(1,ie),ie,modul(ie))
          if(modul(ie)==c00)exit
          call chkgeo(lfacn,nball,nnofa,ie,Qgeo,rfloc(1,ie),rnopo,npoin,ndim)
          if(Qgeo<Qgeotol)then 
             exit
          endif
       enddo
       !
       !     Is the quality check ok?
       !
       if(ie>nball)then
          !
          !     Check that the new position is not too close from other points
          !
          do ipo=1,npclos
             ip1=lpclos(ipo)
             call length(ip1,ipclos,rsize,coor,rlen)
             if(rlen<tollen)then
                ie=nball-1_ip
                exit 
             endif
          enddo
       endif

       if(ie>nball)then
          !
          !     Successful, transfer normals and go home
          !
          do ie=1,nball
             ielem=lball(ie)   
             rnofa(1,ielem)=rfloc(1,ie) 
             rnofa(2,ielem)=rfloc(2,ie) 
             rnofa(3,ielem)=rfloc(3,ie) 
          enddo

          return

       else
          !
          !     Not successfull
          ! 

          cycle
       endif

    enddo
    !
    !     Optimization not successfull, restore rnopo and coor
    !
    coor(1,ipclos)=pold(1)
    coor(2,ipclos)=pold(2)
    coor(3,ipclos)=pold(3)
    rnopo(1,ipclos)=rnopoold(1)
    rnopo(2,ipclos)=rnopoold(2)
    rnopo(3,ipclos)=rnopoold(3)
    lpsur(ipclos)=ipsur

  end subroutine move


  subroutine chkcons(lface,nnofa,nface,eltoel,lelem,coor,npoin,ndim,lmark)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)   :: nnofa,nface,ndim,npoin
    integer(ip), intent(in)   :: lface(nnofa,nface)
    integer(ip), intent(in)   :: eltoel(nnofa,nface),lelem(nface),lmark(npoin)
    real(rp), intent(in)      :: coor(ndim,npoin)
    integer(ip)               :: iface,neigh,j,iview2,ip1,ip2,ip3,ip4
    real(rp)                  :: rface(ndim),rface1(nnofa,ndim),modul,cscal  
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))

    do iface=1,nface

       if(lelem(iface)==-1)cycle
       !
       !     Get normal
       !  

       call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rface,iface,modul)

       !
       !     Check coherence
       !
       if(lmark(lface(1,iface))==-2)then
          write(*,*)'error conform 6a'
          write(*,*)'iface=',iface
          stop
       endif
       if(lmark(lface(2,iface))==-2)then
          write(*,*)'error conform 6b'
          write(*,*)'iface=',iface
          stop
       endif
       if(lmark(lface(3,iface))==-2)then
          write(*,*)'error conform 6c'
          write(*,*)'iface=',iface
          stop
       endif



       do j=1,nnofa

          neigh=eltoel(j,iface)


          if(neigh/=0)then

             if(lelem(neigh)==-1)then
                write(*,*)'error conform 3'
                write(*,*)'iface=',iface,'neigh=',neigh
                write(*,*)lelem(iface),lelem(neigh)
                call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
                stop
             endif

             call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rface1(1,j),iface,modul)

             cscal=rface(1)*rface1(1,j)+rface(2)*rface1(2,j)+rface(3)*rface1(3,j)

             if(cscal<0.0d+00)then
                write(*,*)'error conform 4'
                write(*,*)'iface=',iface,'neigh=',neigh
                write(*,*)lelem(iface),lelem(neigh)
                call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
                stop
             endif




             if(eltoel(1,neigh)==iface)then
                iview2=1_ip
             else if(eltoel(2,neigh)==iface)then
                iview2=2_ip
             else if(eltoel(3,neigh)==iface)then
                iview2=3_ip
             else
                write(*,*)'error conform 1'
                write(*,*)'iface=',iface,'neigh=',neigh
                call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
                stop
             endif

             ip1=lface(ltab(1,j),iface)
             ip2=lface(ltab(2,j),iface)

             ip3=lface(ltab(1,iview2),neigh)
             ip4=lface(ltab(2,iview2),neigh)

             if(ip1/=ip4  .or. ip2/=ip3)then
                write(*,*)'error conform 2'
                write(*,*)'iface=',iface,'neigh=',neigh
                write(*,*)ip1,ip2,ip3,ip4  
                call outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
                stop
             endif



          endif
       enddo
    enddo

  end subroutine chkcons


  subroutine chkcons2(lface,nnofa,nface,eltoel,lptri,npoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)   :: nnofa,nface,npoin
    integer(ip), intent(in)   :: lface(nnofa,nface)
    integer(ip), intent(in)   :: eltoel(nnofa,nface),lptri(npoin)
    integer(ip)               :: iface,neigh,j,iview2,ip1,ip2,ip3,ip4,ichk,ipoin,jface
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))

    do iface=1,nface

       do j=1,nnofa

          neigh=eltoel(j,iface)


          if(neigh/=0)then

             if(eltoel(1,neigh)==iface)then
                iview2=1_ip
             else if(eltoel(2,neigh)==iface)then
                iview2=2_ip
             else if(eltoel(3,neigh)==iface)then
                iview2=3_ip
             else
                write(*,*)'error conform 1'
                write(*,*)'iface=',iface,'neigh=',neigh
                stop
             endif

             ip1=lface(ltab(1,j),iface)
             ip2=lface(ltab(2,j),iface)

             ip3=lface(ltab(1,iview2),neigh)
             ip4=lface(ltab(2,iview2),neigh)

             if(ip1/=ip4  .or. ip2/=ip3)then
                write(*,*)'error conform 2'
                write(*,*)'iface=',iface,'neigh=',neigh
                write(*,*)ip1,ip2,ip3,ip4  
                stop
             endif
          endif
       enddo
    enddo

    do iface=1,nface
       do j=1,nnofa
          ipoin=lface(j,iface)
          jface=lptri(ipoin)
          ichk=0_ip 
          if(lface(1,jface)==ipoin)then
             ichk=1_ip
          else if(lface(2,jface)==ipoin)then
             ichk=1_ip
          else if(lface(3,jface)==ipoin)then
             ichk=1_ip
          endif

          if(ichk==0)then
             write(*,*)'error point/triangle for ipoin:',ipoin
             stop
          endif
       enddo
    enddo

    do iface=1,nface
       ip1=lface(1,iface)
       ip2=lface(2,iface)
       ip3=lface(3,iface)

       if(ip1==ip2)then
          write(*,*)'error chkcons2 aa'
          stop
       endif
       if(ip1==ip3)then
          write(*,*)'error chkcons2 ab'
          stop
       endif
       if(ip2==ip3)then
          write(*,*)'error chkcons2 ac'
          stop
       endif
    enddo
  end subroutine chkcons2


  subroutine outerror(lface,nnofa,nface,lelem,coor,npoin,ndim)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim
    integer(ip), intent(in)      :: lface(nnofa,nface)
    integer(ip), intent(in)      :: lelem(nface)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: i,icont
    !
    !     This sub displays a triangular surface mesh
    !     lelem may be displayed for marked elements to collapse
    !
    open(unit=50,file='outerr.msh',status='unknown')
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
       if(lelem(i)/=-1)then
          icont=icont+1
          write(50,300) icont,lface(1,i),lface(2,i),lface(3,i)
       endif
    enddo
    write(50,6)

    close(50)


1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(10i10)
300 format(4i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')

  end subroutine outerror

  subroutine outerror4(lface,nnofa,nface,coor,npoin,ndim)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim
    integer(ip), intent(in)      :: lface(nnofa,nface)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: i,icont
    !
    !     This sub displays only a triangular surface mesh 
    !
    open(unit=50,file='outerr4.msh',status='unknown')
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
       write(50,300)i,lface(1,i),lface(2,i),lface(3,i)
    enddo
    write(50,6)

    close(50)


1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(10i10)
300 format(4i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')


  end subroutine outerror4

  subroutine outerror2(nface,npoin,ndim,lfmark,nfront,nfhole,lfront,lfhole,rnopo,lfloc,&
             nnofa,coor,nehole,lehole)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)      :: nface,npoin,ndim,nfhole,nfront,nnofa,nehole
    integer(ip), intent(in)      :: lfmark(nface),lfront(3,nfhole+nfront),lehole(nehole)
    integer(ip), intent(in)      :: lfhole(nfhole),lfloc(nnofa,nface)
    real(rp), intent(in)         :: rnopo(ndim,npoin),coor(ndim,npoin)
    integer(ip)                  :: i,icont,ifront,ifhole,ipoin,iehole
    integer(ip), pointer         :: lmark(:),lemark(:)
    integer(4)                   :: istat 
    !
    !     This sub displays a triangular surface mesh, with the front edges 
    !
    allocate(lmark(nfront+nfhole),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','outerror2',lmark) 
    allocate(lemark(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LEMARK','outerror2',lemark) 
    !
    !     Mark the deleted sides of the front
    ! 
    do ifhole=1,nfhole
       lmark(lfhole(ifhole))=1_ip
    enddo
    !
    !     Mark the deleted faces
    ! 
    do iehole=1,nehole
       lemark(lehole(iehole))=1_ip
    enddo

    open(unit=50,file='outerr2.msh',status='unknown')
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
    do i=1,nface
       if(lemark(i)/=1)then
          icont=icont+1
          write(50,300)icont,lfloc(1,i),lfloc(2,i),lfloc(3,i),lfmark(i)
       endif
    enddo
    write(50,6)
    !
    !     Edges
    !
    write(50,7)
    write(50,2)
    write(50,3)
    do i=1,npoin
       write(50,100)npoin+i,coor(1,i),coor(2,i),coor(3,i)
    enddo
    write(50,4)
    write(50,5)


    icont=0_ip
    do ifront=1,nfront+nfhole
       if(lmark(ifront)==0)then
          icont=icont+1
          write(50,400)icont,lfront(1,ifront)+npoin,lfront(2,ifront)+npoin
       endif
    enddo
    write(50,6)

    close(50)

    open(unit=60,file='outerr2.res',status='unknown')
    rewind 60

10  format('GID Post Results File 1.0')
11  format('Result "normals" "Analysis/time" 1 Vector OnNodes')
12  format('ComponentNames "Vx", "Vy" , "Vz" ')
13  format('Values')
14  format('End Values')
15  format('   ')
    write(60,10)
    write(60,11)
    write(60,12)
    write(60,13)
    do  ipoin=1,npoin
       write(60,100)ipoin,rnopo(1,ipoin),rnopo(2,ipoin),rnopo(3,ipoin)
    enddo
    write(60,14)
    write(60,15)
    close(60)

    call memchk(2_ip,istat,memor_msh,'LMARK','outerror2',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','outerror2',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEMARK','outerror2',lemark)
    deallocate(lemark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEMARK','outerror2',0_ip)

1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(10i10)
300 format(5i10)
400 format(3i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
7   format('MESH dimension 3 ElemType Linear Nnode 2')


  end subroutine outerror2

  subroutine outerrort(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,&
       lline,nline)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH,memor_msh
    use mod_memchk
    implicit none 
    integer(ip),intent(in)    :: nface,nnofa,ndim,npoin,nsurf,nside,nline 
    integer(ip),intent(in)    :: lface(nnofa,nface),lside(2,nside),lline(nside)
    integer(ip),intent(in)    :: lsurf(nface)
    real(rp),intent(in)       :: coor(ndim,npoin)
    integer(ip)               :: i,icont

    open(unit=50,file='outerrort.msh',status='unknown')
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
       !write(50,300)icont,lface(1,i),lface(2,i),lface(3,i),lsurf(i)
       write(50,300)icont,lface(1,i),lface(2,i),lface(3,i),1
    enddo

    write(50,6)

    write(50,2)
    write(50,3)
    write(50,4)
    write(50,11)
    write(50,5)

    icont=0_ip
    do i=1, nside
       icont=icont+1
       !write(50,400)icont,lside(1,i),lside(2,i),lline(i)
       write(50,400)icont,lside(1,i),lside(2,i),2
    enddo

    write(50,6)

    close(50)

1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
11  format('MESH dimension 3 ElemType  Linear Nnode 2')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
300 format(5i10)
400 format(4i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')


  end subroutine outerrort

  subroutine outmark(lface,nnofa,nface,lelem,coor,npoin,ndim,lmark)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim
    integer(ip), intent(in)      :: lface(nnofa,nface)
    integer(ip), intent(in)      :: lelem(nface),lmark(npoin)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: i,icont,ipoin
    !
    !     This sub displays a triangular surface mesh with the points marked for collapsing 
    !
    open(unit=50,file='outmark.msh',status='unknown')
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
       if(lelem(i)/=-1)then
          icont=icont+1
          write(50,300)icont,lface(1,i),lface(2,i),lface(3,i)
       endif
    enddo
    write(50,6)


    close(50)


1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(2i10)
300 format(4i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')


    open(unit=60,file='outmark.res',status='unknown')
    rewind 60

10  format('GID Post Results File 1.0')
11  format('Result "lmark" "Analysis/time" 1 Scalar OnNodes')
13  format('Values')
14  format('End Values')
15  format('   ')
    write(60,10)
    write(60,11)
    write(60,13)
    do  ipoin=1,npoin
       write(60,200)ipoin,lmark(ipoin)
    enddo
    write(60,14)
    write(60,15)
    close(60)


  end subroutine outmark

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

  subroutine optsmo(lface,nnofa,nface,coor,npoin,ndim,rnopo,eltoel,niter,coorold,&
       nfold,npold,rnofaold,lfold,eltoelold,lptri,lpsur,lelemold,rsize,lptype,lsurfold,isurf,&
       ptosi1,ptosi2,lside,nnosi,nside)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       : memor_msh
    use mod_memchk
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim,niter,isurf,nfold,npold,nnosi,nside
    integer(ip), intent(in)      :: lptype(2,npoin),lside(nnosi,nside)
    integer(ip), intent(inout)   :: lface(nnofa,nface),eltoel(nnofa,nface)
    real(rp), intent(inout)      :: rnopo(ndim,npoin),coor(ndim,npoin)
    integer(ip),intent(in)       :: eltoelold(nnofa,nfold),lfold(nnofa,nfold),lsurfold(nfold)
    integer(ip),intent(inout)    :: lpsur(npoin),lelemold(nfold),lptri(npoin)
    real(rp),intent(in)          :: coorold(ndim,npold)
    real(rp),intent(in)          :: rnofaold(ndim,nfold),rsize(npoin)
    integer(ip)                  :: iface,jface,iopt,j,iter,icont,ipoin,ipo,ipnt
    integer(ip)                  :: isto,inofa,iside,ip1,ip2,ipc,ineigh,icheck
    integer(ip),pointer          :: lmark(:),ptosi1(:),ptosi2(:)
    integer(4)                   :: istat 
    integer(ip)                  :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub tries to optimize the hole surface triangulation in lface
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','optsmo',lmark) 
    !
    !     Get the points of this surface
    !
    ipnt=0_ip
    do ipoin=1,npoin
       if(lptri(ipoin)/=0)then
          ipnt=ipnt+1
          lmark(ipnt)=ipoin
       endif
    enddo
    !
    !     Loop on optimization
    !
    do iter=1,niter
       !
       !     Swap the faces
       !
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

                   iopt=0_ip
                   call swap(iface,jface,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,iopt)
                   !
                   !     Is the swap sucessfull?
                   ! 
                   if(iopt==1)then
                      !
                      !     Update lptri outside swap
                      !
                      lptri(lface(1,iface))=iface
                      lptri(lface(2,iface))=iface
                      lptri(lface(3,iface))=iface
                      lptri(lface(1,jface))=jface
                      lptri(lface(2,jface))=jface
                      lptri(lface(3,jface))=jface
                   endif
                endif
             endif
          enddo
       enddo
       !call outerror4(lface,nnofa,nface,coor,npoin,ndim)
       !
       !     Optimize the point placement
       !
       !write(*,*)'Iter=',iter 

       do ipo=1,ipnt
          ipoin=lmark(ipo)
          !write(*,*)'ipoin=',ipoin
          !if(isurf==2)call outerror4(lface,nnofa,nface,coor,npoin,ndim)

          call moveopt(ipoin,coor,ndim,npoin,lface,nnofa,nface,rnopo,coorold,nfold,&
               npold,rnofaold,lfold,eltoelold,eltoel,lptri,lpsur,lelemold,rsize,lptype,&
               lsurfold,isurf)
       enddo
       !call outerror4(lface,nnofa,nface,coor,npoin,ndim)

    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','optsmo',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','optsmo',0_ip)

  end subroutine optsmo

  subroutine recover(ipa,ipb,lface,nnofa,nface,lptri,coor,npoin,ndim,eltoel,rnopo,&
       ipnew,rnofa,lfmark,ienew,iview,ierr,lside,nnosi,nside,ptosi1,ptosi2,npcusp,&
       lphole,nphole,lehole,nehole,lelem)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       : memor_msh
    use mod_memchk
    implicit none
    integer(ip), intent(in)   :: ipa,ipb,nnofa,nface,npoin,ndim,ipnew,npcusp,nnosi,nside
    integer(ip), intent(inout):: nphole,nehole
    integer(ip), intent(inout):: lface(nnofa,nface),eltoel(nnofa,nface),lelem(nface)
    integer(ip),pointer       :: lehole(:),lphole(:)
    integer(ip), intent(inout):: lptri(npoin),lfmark(nface),ienew,iview,ierr
    integer(ip), intent(in)   :: lside(nnosi,nside),ptosi2(npcusp+1),ptosi1(*)
    real(rp), intent(in)      :: coor(ndim,npoin),rnopo(ndim,npoin)
    real(rp), intent(inout)   :: rnofa(nnofa,nface)
    integer(ip)               :: ie,j,ip1,ip2,ip3,ienext,iopt,ielem2,iview1,nitermax,iter
    integer(ip)               :: npipe,ieold,ipipe,ienext1,ienext2,iturn,ifirst
    integer(ip)               :: ipaa,ipbb,isto,ipmin,ipmax,ipcusp,ifound,iside,inofa,iface 
    real(rp)                  :: rface(ndim),rx,ry,rz,rscal,pproj(ndim),pnew(ndim)
    real(rp)                  :: p1(ndim),p2(ndim),csca,c10
    real(rp)                  :: d1,d2,d3,epsil,dtot,rtx,rty,rtz,rtl,rnx,rny,rnz,rnl  
    real(rp)                  :: rnopx,rnopy,rnopz  
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    integer(ip),parameter     :: mstack=500,mpipe=100
    integer(ip)               :: lstack(mstack),lstackp(mstack),lpipe(mpipe)
    integer(ip)               :: nstack,nstackp,ineigh1,ineigh2,ineigh3,istack,istackp   
    integer(ip)               :: is11,is12,is21,is22,is31,is32,ineigh,ipoin   
    !
    !     This sub recovers the sides (ipa,ipnew) and (ipb,ipnew)
    !
    !
    !     On output:  ienew, the new element 
    !                 iview, the viewer in ienew
    !
    !  
    c10=1.0d+00
    epsil=1.0d-09
    nitermax=10000_ip
    pnew(1)=coor(1,ipnew)
    pnew(2)=coor(2,ipnew)
    pnew(3)=coor(3,ipnew)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Find the  element containing edge (ipa,ipb)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ie=lptri(ipa)
    iturn=1_ip

    do 

       if(lface(1,ie)==ipa)then
          ienext1=eltoel(2,ie)
          ienext2=eltoel(3,ie)
          iview1=1_ip
       else if(lface(2,ie)==ipa)then
          ienext1=eltoel(3,ie)
          ienext2=eltoel(1,ie)
          iview1=2_ip
       else
          ienext1=eltoel(1,ie)
          ienext2=eltoel(2,ie)
          iview1=3_ip
       endif

       if(ipb==lface(ltab(1,iview1),ie))exit
       !
       !     Did we reach the boundary of the patch?
       !
       if(iturn==1)then

          if(ienext1/=0)then
             ie=ienext1
          else
             ie=ienext2
             if(ie==0)exit
             iturn=0_ip
          endif
       else
          ie=ienext2
          if(ie==0)exit
       endif

    enddo
    !
    !     Did we find ie?
    !
    if(ie==0)then
       !
       !     Surfaces may be disconnected
       !
       do iface=1,nface           
          do inofa=1,nnofa 
             if(lface(ltab(1,inofa),iface)==ipa .and. lface(ltab(2,inofa),iface)==ipb)then
                ie=iface
                goto 33
             endif
          enddo
       enddo
       write(*,*)'Error in recover finding ipa,ipb'
       ierr=1
       return

    endif
33  continue
    !
    !     Remember it
    !
    ifirst=ie
    !
    !    Loop until the edges have been regenerated
    !

    !
    !     First regenerate edge (ipa,ipnew)
    !
5   continue 
    !
    !     Find the element in the ball of ipa cut by edge (ipa,ipnew)
    !
    ie=ifirst
    iter=0_ip

    rtx=pnew(1)-coor(1,ipa)
    rty=pnew(2)-coor(2,ipa)
    rtz=pnew(3)-coor(3,ipa)
    rtl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
    rtl=c10/rtl
    rtx=rtx*rtl
    rty=rty*rtl
    rtz=rtz*rtl
    !
    !     Initialize the normal to evaluate intersection
    !
    rnopx=rnopo(1,ipa)+rnopo(1,ipnew)
    rnopy=rnopo(2,ipa)+rnopo(2,ipnew)
    rnopz=rnopo(3,ipa)+rnopo(3,ipnew)
    rnl=sqrt(rnopx*rnopx+rnopy*rnopy+rnopz*rnopz)
    rnl=c10/rnl
    rnopx=rnopx*rnl 
    rnopy=rnopy*rnl 
    rnopz=rnopz*rnl 
    !
    !     Get normal to intersection plane
    !
    rnx= rnopy*rtz-rnopz*rty
    rny=-rnopx*rtz+rnopz*rtx
    rnz= rnopx*rty-rnopy*rtx
    rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
    rnl=c10/rnl
    rnx=rnx*rnl
    rny=rny*rnl
    rnz=rnz*rnl

    do
       iter=iter+1_ip
       if(iter==nitermax)then
          write(*,*)'Recover 1 iter==nitermax' 
          ierr=1
          return
       endif

       if(lface(1,ie)==ipa)then
          iview1=1_ip
       else if(lface(2,ie)==ipa)then
          iview1=2_ip
       else
          iview1=3_ip
       endif

       ip1=ipa
       ip2=lface(ltab(1,iview1),ie)
       ip3=lface(ltab(2,iview1),ie)

       if(ip2==ipnew)goto 10
       if(ip3==ipnew)goto 10
       !
       !     Find which side is crossed
       ! 
       rface(1)=rnofa(1,ie)
       rface(2)=rnofa(2,ie)
       rface(3)=rnofa(3,ie)
       !
       !     Side (ip1,ip3)
       !
       rx=coor(1,ip3)-coor(1,ipa)
       ry=coor(2,ip3)-coor(2,ipa)
       rz=coor(3,ip3)-coor(3,ipa)

       csca=rnx*rx+rny*ry+rnz*rz

       if(csca<-epsil)then
          ie=eltoel(ltab(1,iview1),ie)
          cycle
       endif
       !
       !     Side (ip1,ip2)
       !
       rx=coor(1,ip2)-coor(1,ipa)
       ry=coor(2,ip2)-coor(2,ipa)
       rz=coor(3,ip2)-coor(3,ipa)

       csca=rnx*rx+rny*ry+rnz*rz

       if(csca>epsil)then
          ie=eltoel(ltab(2,iview1),ie)
          cycle
       endif
       !
       !     Element found swap it
       !
       !
       !     Has ie been imposed before? 
       !
       if(lfmark(ie)==1)then
          write(*,*)'First element already marked'
          ierr=1 
          return
       endif
       ielem2=eltoel(iview1,ie)
       !
       !     Has ielem2 been imposed before? 
       !
       if(lfmark(ielem2)==1)then
          write(*,*)'Second element already marked'
          ierr=1
          return
       endif
       !
       !     Do we have inner sides?
       !
       ip1=min(lface(ltab(1,iview1),ie),lface(ltab(2,iview1),ie))
       ip2=max(lface(ltab(1,iview1),ie),lface(ltab(2,iview1),ie))

       if(ip1<=npcusp .or. ip2<=npcusp)then

          if(ip1<=npcusp)then
             ipcusp=ip1
          else
             ipcusp=ip2
          endif

          ifound=0
          do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
             iside=ptosi1(isto)
             ipaa=min(lside(1,iside),lside(2,iside)) 
             ipbb=max(lside(1,iside),lside(2,iside)) 
             if(ip1==ipaa .and. ip2==ipbb)then
                ifound=1
                exit
             endif
          enddo

          if(ifound==1)then
             write(*,*)'Error in recover, must swap an inner side'
             ierr=1
             return
          endif

       endif

       iopt=0_ip
       call swaprecover(ie,ielem2,lface,nnofa,nface,eltoel,coor,rnopo,&
            ndim,npoin,iopt,lptri,rnofa)
       if(iopt==0)then
          exit
       endif


    enddo
    !
    !     At this point, we could not swap the first element
    !     Find the pipe
    !
    lpipe(1)=ie
    npipe=1_ip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Go through the pipe
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ieold=ie
    ie=ielem2
    iter=0_ip

    do

       iter=iter+1_ip
       if(iter==nitermax)then
          write(*,*)'Recover 2 iter==nitermax' 
          ierr=1
          return
       endif

       npipe=npipe+1
       lpipe(npipe)=ie

       if(eltoel(1,ie)==ieold)then
          iview1=1_ip
       else if(eltoel(2,ie)==ieold)then
          iview1=2_ip
       else
          iview1=3_ip
       endif

       ip1=lface(iview1,ie)
       ip2=lface(ltab(1,iview1),ie)
       ip3=lface(ltab(2,iview1),ie)
       !
       !     Did we cross the front
       !
       if(lfmark(ie)==1)then
          write(*,*)'Element already marked'
          ierr=1_ip
          return
       endif
       !
       !     Does ipnew belong to this element
       !
       if(ipnew==ip1)exit
       !
       !     Do we have inner sides?
       !
       ipmin=min(ip2,ip3)
       ipmax=max(ip2,ip3)

       if(ipmin<=npcusp .or. ipmax<=npcusp)then

          if(ipmin<=npcusp)then
             ipcusp=ipmin
          else
             ipcusp=ipmax
          endif

          ifound=0
          do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
             iside=ptosi1(isto)
             ipaa=min(lside(1,iside),lside(2,iside)) 
             ipbb=max(lside(1,iside),lside(2,iside)) 
             if(ipmin==ipaa .and. ipmax==ipbb)then
                ifound=1
                exit
             endif
          enddo

          if(ifound==1)then
             write(*,*)'Error in recover, must swap an inner side'
             ierr=1
             return
          endif

       endif
       !
       !     Find which side is crossed
       ! 
       rface(1)=rnofa(1,ie)
       rface(2)=rnofa(2,ie)
       rface(3)=rnofa(3,ie)
       !
       !     Which side between (ip1,ip2) and (ip1,ip3) is being crossed ?
       !
       rx=coor(1,ip1)-coor(1,ipa)
       ry=coor(2,ip1)-coor(2,ipa)
       rz=coor(3,ip1)-coor(3,ipa)

       csca=rx*rnx+ry*rny+rz*rnz

       if(csca>epsil)then
          ieold=ie
          ie=eltoel(ltab(1,iview1),ie)
       else 
          ieold=ie
          ie=eltoel(ltab(2,iview1),ie)
       endif

    enddo
    !
    !     Swap the hole pipe except the first element
    !
    if(npipe==2)then
       write(*,*)'Elements not swappable 1'
       ierr=1
       return
    endif
    !
    !     Loop on iterations
    !
    iter=0_ip
    do  
       iter=iter+1_ip
       if(iter==nitermax)then
          write(*,*)'Recover 3 iter==nitermax' 
          ierr=1
          return
       endif

       ipipe=1_ip
       !
       !      Loop on the pipe
       !
       do    

          ipipe=ipipe+1
          if(ipipe>=npipe)exit

          ie=lpipe(ipipe)
          ielem2=lpipe(ipipe+1)

          iopt=0_ip
          call swaprecover(ie,ielem2,lface,nnofa,nface,eltoel,coor,rnopo,&
               ndim,npoin,iopt,lptri,rnofa)
          if(iopt==1)then
             call updatepipe(ipipe,ipipe+2_ip,lpipe,npipe,ipa,ipnew,rnofa,nnofa,&
                  nface,coor,ndim,npoin,lfmark,pnew,eltoel,lface,ierr,&
                  rnx,rny,rnz,ipa,mpipe)
             if(ierr==1)then
                write(*,*)'Error in update pipe 1' 
                return
             endif
             !
             !     We are almost done
             !
             if(npipe==2)goto 5

          endif

       enddo

    enddo

10  continue
    !
    !     Then regenerate edge (ipb,ipnew)
    !
    !
    !     Find the element in the ball of ipb cut by edge (ipb,ipnew)
    !
    ie=ifirst
    iter=0_ip

    rtx=pnew(1)-coor(1,ipb)
    rty=pnew(2)-coor(2,ipb)
    rtz=pnew(3)-coor(3,ipb)
    rtl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
    rtl=c10/rtl
    rtx=rtx*rtl
    rty=rty*rtl
    rtz=rtz*rtl
    !
    !     Initialize the normal to evaluate intersection
    !
    rnopx=rnopo(1,ipb)+rnopo(1,ipnew)
    rnopy=rnopo(2,ipb)+rnopo(2,ipnew)
    rnopz=rnopo(3,ipb)+rnopo(3,ipnew)
    rnl=sqrt(rnopx*rnopx+rnopy*rnopy+rnopz*rnopz)
    rnl=c10/rnl
    rnopx=rnopx*rnl 
    rnopy=rnopy*rnl 
    rnopz=rnopz*rnl 
    !
    !     Compute normal to plane of intersection
    !
    rnx= rnopy*rtz-rnopz*rty
    rny=-rnopx*rtz+rnopz*rtx
    rnz= rnopx*rty-rnopy*rtx
    rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
    rnl=c10/rnl
    rnx=rnx*rnl
    rny=rny*rnl
    rnz=rnz*rnl

    do
       iter=iter+1_ip
       if(iter==nitermax)then
          write(*,*)'Recover 4 iter==nitermax' 
          ierr=1
          return
       endif

       if(lface(1,ie)==ipb)then
          iview1=1_ip
       else if(lface(2,ie)==ipb)then
          iview1=2_ip
       else
          iview1=3_ip
       endif

       ip1=ipb
       ip2=lface(ltab(1,iview1),ie)
       ip3=lface(ltab(2,iview1),ie)

       if(ip2==ipnew)goto 20
       if(ip3==ipnew)goto 20
       !
       !     Find which side is crossed
       ! 
       rface(1)=rnofa(1,ie)
       rface(2)=rnofa(2,ie)
       rface(3)=rnofa(3,ie)
       !
       !     Side (ip1,ip2)
       !
       rx=coor(1,ip2)-coor(1,ipb)
       ry=coor(2,ip2)-coor(2,ipb)
       rz=coor(3,ip2)-coor(3,ipb)

       csca=rnx*rx+rny*ry+rnz*rz

       if(csca>epsil)then
          ie=eltoel(ltab(2,iview1),ie)
          cycle
       endif
       !
       !     Side (ip1,ip3)
       !
       rx=coor(1,ip3)-coor(1,ipb)
       ry=coor(2,ip3)-coor(2,ipb)
       rz=coor(3,ip3)-coor(3,ipb)

       csca=rnx*rx+rny*ry+rnz*rz

       if(csca<-epsil)then
          ie=eltoel(ltab(1,iview1),ie)
          cycle
       endif
       !
       !     Element found swap it
       !
       !
       !     Has ie been imposed before ? 
       !
       if(lfmark(ie)==1)then
          write(*,*)'First element already marked a'
          ierr=1
          return
       endif
       ielem2=eltoel(iview1,ie)
       !
       !     Has ielem2 been imposed before ? 
       !
       if(lfmark(ielem2)==1)then
          write(*,*)'Second element already marked b'
          ierr=1
          return
       endif
       !
       !     Do we have inner sides?
       !
       ip1=min(lface(ltab(1,iview1),ie),lface(ltab(2,iview1),ie))
       ip2=max(lface(ltab(1,iview1),ie),lface(ltab(2,iview1),ie))

       if(ip1<=npcusp .or. ip2<=npcusp)then

          if(ip1<=npcusp)then
             ipcusp=ip1
          else
             ipcusp=ip2
          endif

          ifound=0
          do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
             iside=ptosi1(isto)
             ipaa=min(lside(1,iside),lside(2,iside)) 
             ipbb=max(lside(1,iside),lside(2,iside)) 
             if(ip1==ipaa .and. ip2==ipbb)then
                ifound=1
                exit
             endif
          enddo

          if(ifound==1)then
             write(*,*)'Error in recover, must swap an inner side'
             ierr=1
             return
          endif

       endif

       iopt=0_ip
       call swaprecover(ie,ielem2,lface,nnofa,nface,eltoel,coor,rnopo,&
            ndim,npoin,iopt,lptri,rnofa)
       if(iopt==0)then
          exit
       endif

    enddo
    !
    !     At this point, we could not swap the first element
    !     Find the pipe
    !
    lpipe(1)=ie
    npipe=1_ip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Go through the pipe
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ieold=ie
    ie=ielem2
    iter=0_ip

    do
       iter=iter+1_ip
       if(iter==nitermax)then
          write(*,*)'Recover 5 iter==nitermax' 
          ierr=1
          return
       endif

       npipe=npipe+1
       lpipe(npipe)=ie

       if(eltoel(1,ie)==ieold)then
          iview1=1_ip
       else if(eltoel(2,ie)==ieold)then
          iview1=2_ip
       else
          iview1=3_ip
       endif

       ip1=lface(iview1,ie)
       ip2=lface(ltab(1,iview1),ie)
       ip3=lface(ltab(2,iview1),ie)
       !
       !     Did we cross the front
       !
       if(lfmark(ie)==1)then
          write(*,*)'Element already marked'
          ierr=1_ip
          return
       endif
       !
       !     Does ipnew belong to this element
       !
       if(ipnew==ip1)exit
       !
       !     Do we have inner sides?
       !
       ipmin=min(ip2,ip3)
       ipmax=max(ip2,ip3)

       if(ipmin<=npcusp .or. ipmax<=npcusp)then

          if(ipmin<=npcusp)then
             ipcusp=ipmin
          else
             ipcusp=ipmax
          endif

          ifound=0
          do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
             iside=ptosi1(isto)
             ipaa=min(lside(1,iside),lside(2,iside)) 
             ipbb=max(lside(1,iside),lside(2,iside)) 
             if(ipmin==ipaa .and. ipmax==ipbb)then
                ifound=1
                exit
             endif
          enddo

          if(ifound==1)then
             write(*,*)'Error in recover, must swap an inner side'
             ierr=1
             return
          endif

       endif
       !
       !     Find which side is crossed
       ! 
       rface(1)=rnofa(1,ie)
       rface(2)=rnofa(2,ie)
       rface(3)=rnofa(3,ie)
       !
       !     Which side between (ip1,ip2) and (ip1,ip3) is being crossed ?
       !
       rx=coor(1,ip1)-coor(1,ipb)
       ry=coor(2,ip1)-coor(2,ipb)
       rz=coor(3,ip1)-coor(3,ipb)

       csca=rx*rnx+ry*rny+rz*rnz      

       if(csca>epsil)then
          ieold=ie
          ie=eltoel(ltab(1,iview1),ie)
       else
          ieold=ie
          ie=eltoel(ltab(2,iview1),ie)
       endif

    enddo
    !
    !     Swap the hole pipe except the first element
    !
    if(npipe==2)then
       write(*,*)'Elements not swappable 2'
       ierr=1
       return
    endif
    !
    !     Loop on iterations
    !
    iter=0_ip
    do  
       iter=iter+1_ip
       if(iter==nitermax)then
          write(*,*)'Recover 6 iter==nitermax' 
          ierr=1
          return
       endif

       ipipe=1_ip
       !
       !      Loop on the pipe
       !
       do    

          ipipe=ipipe+1
          if(ipipe>=npipe)exit

          ie=lpipe(ipipe)
          ielem2=lpipe(ipipe+1)

          iopt=0_ip
          call swaprecover(ie,ielem2,lface,nnofa,nface,eltoel,coor,&
               rnopo,ndim,npoin,iopt,lptri,rnofa)
          if(iopt==1)then
             call updatepipe(ipipe,ipipe+2_ip,lpipe,npipe,ipb,ipnew,rnofa,&
                  nnofa,nface,coor,ndim,npoin,lfmark,pnew,eltoel,lface,ierr,&
                  rnx,rny,rnz,ipa,mpipe)
             if(ierr==1)then
                write(*,*)'Error in update pipe 2' 
                return
             endif
             !
             !     We are almost done
             !
             if(npipe==2)goto 10
          endif

       enddo

    enddo

20  continue
    !
    !     Find the element limited by the front
    !
    if(lface(1,ie)==ipb)then 
       iview1=1_ip 
    else if(lface(2,ie)==ipb)then 
       iview1=2_ip 
    else 
       iview1=3_ip 
    endif
    !
    !     We may have reached ipnew from the neighbor element not containing ipa,ipb
    !
    if(lface(ltab(1,iview1),ie)/=ipnew)then

       ie=eltoel(ltab(1,iview1),ie)

       if(lface(1,ie)==ipb)then 
          iview1=1_ip 
       else if(lface(2,ie)==ipb)then 
          iview1=2_ip 
       else 
          iview1=3_ip 
       endif

    endif

    ienew=ie
    iview=ltab(1,iview1) 
    lfmark(ie)=1_ip
    !
    !     Did we go over some points, so that ipnew does not belong to ienew 
    !
    if(lface(ltab(1,iview1),ie)/=ipnew)then
       !
       !     Add all the elements contained in (ipa,ipb,ipnew)
       !
       nstack=1
       lstack(1)=ie
       lelem(ie)=1
       istack=0
       is11=min(ipa,ipb)
       is12=max(ipa,ipb)  
       is21=min(ipb,ipnew)
       is22=max(ipb,ipnew)  
       is31=min(ipnew,ipa)
       is32=max(ipnew,ipa)  

       do 

          if(istack==nstack)exit
          istack=istack+1  
          iface=lstack(istack)
          do inofa=1,nnofa
             ineigh=eltoel(inofa,iface)
             if(ineigh==0)cycle
             if(lelem(ineigh)/=0)cycle   
             lelem(ineigh)=1      

             ip1=lface(ltab(1,inofa),iface)  
             ip2=lface(ltab(2,inofa),iface)
             ipmin=min(ip1,ip2)
             ipmax=max(ip1,ip2)  
             !
             !     Can not cross sides 
             !
             if(ipmin==is11 .and. ipmax==is12)then
                ineigh1=ineigh  
                cycle
             endif
             if(ipmin==is21 .and. ipmax==is22)then
                ineigh2=ineigh  
                cycle
             endif
             if(ipmin==is31 .and. ipmax==is32)then
                ineigh3=ineigh 
                cycle
             endif
             !
             !     Add to the stack
             !
             nstack=nstack+1
             lstack(nstack)=ineigh

          enddo

       enddo
       !
       !     Now create the new element
       ! 
       lface(1,ie)=ipa 
       lface(2,ie)=ipb 
       lface(3,ie)=ipnew 
       iview=3_ip

       write(*,*)'New element:',ipa,ipb,ipnew

       eltoel(1,ie)=ineigh2 
       eltoel(2,ie)=ineigh3
       eltoel(3,ie)=ineigh1
    
       lptri(ipa)=ie
       lptri(ipb)=ie
       lptri(ipnew)=ie
       !
       !     Add points to lphole
       ! 
       nstackp=0_ip
       do istack=1,nstack
          iface=lstack(istack)
          do inofa=1,nnofa
             ipoin=lface(inofa,iface)
             if(ipoin/=ipa .and. ipoin/=ipb .and. ipoin/=ipnew)then
                do istackp=1,nstackp
                   if(ipoin==lstackp(istackp))exit
                enddo
                if(istackp==nstackp+1)then
                   nstackp=nstackp+1
                   lstackp(nstackp)=ipoin
                endif
             endif
          enddo
       enddo

       call memrea(nstackp+nphole,memor_msh,'LPHOLE','recover',lphole) 

       do istackp=1,nstackp
          write(*,*)'nphole=',nphole,'lstackp(istackp)=',lstackp(istackp)
          nphole=nphole+1
          lphole(nphole)=lstackp(istackp)
       enddo
       !
       !     Add the elements to lehole
       !
       lelem(ie)=0_ip

       call memrea(nstack-1_ip+nehole,memor_msh,'LEHOLE','recover',lehole) 

       do istack=2,nstack
          nehole=nehole+1
          write(*,*)'nehole=',nehole,'lstack(istack)=',lstack(istack)
          lehole(nehole)=lstack(istack)
          lelem(lstack(istack))=0
       enddo

    endif

  end subroutine recover

  subroutine newfront(nfapo,nfront,ifront,ipnew,npoin,nfahol,coor,ndim,nfhole,nheap,&
       lfront,lheap,lpofa,lfapo,rsize,rfront,lfahol,lfhole)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)     :: ifront,ipnew,npoin,ndim
    real(rp), intent(in)        :: coor(ndim,npoin)
    integer(ip),pointer         :: lfront(:,:),lheap(:),lpofa(:),lfapo(:,:)
    integer(ip),pointer         :: lfahol(:),lfhole(:)
    real(rp),pointer            :: rsize(:),rfront(:)
    integer(ip), intent(inout)  :: nfapo,nfront,nfahol,nfhole,nheap
    integer(ip)                 :: ip1,ip2,ip1old,ip2old,ipos,iplace,idel
    integer(ip)                 :: nfaloc,lfaloc(3,100),nfface,ifaloc,iface,ipa,ipb 


    ip1old=lfront(1,ifront)
    ip2old=lfront(2,ifront)

    !
    !     Delete ifront in lfapo for ip1 and ip2, the front endpoints
    !
    call delfac(ifront,lfront,nfront,lfapo,nfapo,lpofa,npoin,nfahol,nfhole,&
         lfahol,lfhole)
    !
    !     Add new faces if appropriate
    !
    if(lpofa(ipnew)==0)then
       !
       !     First side
       !
       call addfac(nfront,nfapo,ipnew,ip2old,nfahol,coor,ndim,npoin,nfhole,nheap,&
            lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole)
       !
       !     Second side
       !
       call addfac(nfront,nfapo,ip1old,ipnew,nfahol,coor,ndim,npoin,nfhole,nheap,&
            lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole)

    else
       !
       !     Must delete doubly defined faces when adding (ip2old,ipnew) and (ipnew,ip1old)
       !     or add them if not found
       ! 
       ipos=lpofa(ipnew)
       if(ipos<=0)then
          write(*,*)'Error in newfront, ipos=',ipos 
          stop
       endif
       !
       !     Copy all the faces in a local array
       ! 
       nfaloc=0_ip

       do 

          nfface=lfapo(3,ipos)
          if(nfface==2)then
             nfaloc=nfaloc+1
             lfaloc(1,nfaloc)=lfapo(1,ipos)
             lfaloc(2,nfaloc)=1
             lfaloc(3,nfaloc)=ipos

             nfaloc=nfaloc+1
             lfaloc(1,nfaloc)=lfapo(2,ipos)
             lfaloc(2,nfaloc)=2
             lfaloc(3,nfaloc)=ipos
             exit
          else if(nfface==1)then
             nfaloc=nfaloc+1
             lfaloc(1,nfaloc)=lfapo(1,ipos)
             lfaloc(2,nfaloc)=1
             lfaloc(3,nfaloc)=ipos
             exit
          else 
             nfaloc=nfaloc+1
             lfaloc(1,nfaloc)=lfapo(1,ipos)
             lfaloc(2,nfaloc)=1
             lfaloc(3,nfaloc)=ipos
             nfaloc=nfaloc+1
             lfaloc(1,nfaloc)=lfapo(2,ipos)
             lfaloc(2,nfaloc)=2
             lfaloc(3,nfaloc)=ipos
          endif
          ipos=-nfface
          !
          !     DBG 
          !
          if(ipos<=0)then
             write(*,*)'Error in newfront, ipos:',ipos
             stop
          endif

       enddo
       !
       !     First side
       !
       idel=0_ip
       do ifaloc=1,nfaloc
          iface=lfaloc(1,ifaloc)    
          ipa=lfront(1,iface)
          ipb=lfront(2,iface)

          if(ipa==ip2old .and. ipb==ipnew)then
             !
             !     Face (ipnew,ip2) found
             !
             call delheap(lfront(3,iface),lheap,nheap,lfront,nfront+nfhole,rfront)
             call delfac(iface,lfront,nfront,lfapo,nfapo,lpofa,npoin,nfahol,nfhole,&
                  lfahol,lfhole)
             idel=1_ip
          endif
       enddo
       !
       !     Do we have to add this face if not found?
       !
       if(idel==0)then

          call addfac(nfront,nfapo,ipnew,ip2old,nfahol,coor,ndim,npoin,nfhole,nheap,&
               lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole)

       endif
       !
       !     Second side
       !
       idel=0_ip
       do ifaloc=1,nfaloc
          iface=lfaloc(1,ifaloc)    
          ipa=lfront(1,iface)
          ipb=lfront(2,iface)

          if(ipa==ipnew .and. ipb==ip1old) then

             call delheap(lfront(3,iface),lheap,nheap,lfront,nfront+nfhole,rfront)
             call delfac(iface,lfront,nfront,lfapo,nfapo,lpofa,npoin,nfahol,nfhole,&
                  lfahol,lfhole)
             idel=1_ip
          endif
       enddo
       !
       !     Do we have to add this face if not found?
       !
       if(idel==0)then

          call addfac(nfront,nfapo,ip1old,ipnew,nfahol,coor,ndim,npoin,nfhole,nheap,&
               lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole)

       endif

    endif

  end subroutine newfront

  subroutine delfac(ifront,lfront,nfront,lfapo,nfapo,lpofa,npoin,nfahol,nfhole,&
       lfahol,lfhole)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)     ::  ifront,npoin 
    integer(ip), intent(inout)  ::  nfront,nfapo,nfahol,nfhole 
    integer(ip), intent(inout)  ::  lfront(3,nfront+nfhole),lfapo(3,nfapo+nfahol),lpofa(npoin) 
    integer(ip), pointer        ::  lfahol(:),lfhole(:)
    integer(ip)                 ::  ip1old,ip2old,ipos,nfface,ip1,ip2,jpos,jplace,irow,iposn
    !
    !     Delete ifront in lfapo for ip1 and ip2, the front endpoints
    !
    ip1old=lfront(1,ifront)
    ip2old=lfront(2,ifront)
    !
    !     Delete for ip1old
    !
    ipos=lpofa(ip1old)
    iposn=0_ip
    jpos=0_ip
    jplace=0_ip
    irow=0_ip

    do
       nfface=lfapo(3,ipos)
       if(nfface==2)then
          if(lfapo(1,ipos)==ifront)then
             lfapo(1,ipos)=lfapo(2,ipos)
          else if(lfapo(2,ipos)==ifront)then   
          else
             lfapo(jplace,jpos)=lfapo(2,ipos)
          endif
          lfapo(2,ipos)=0_ip
          lfapo(3,ipos)=1_ip
          exit
       else if(nfface==1)then 
          if(lfapo(1,ipos)==ifront)then 
          else
             lfapo(jplace,jpos)=lfapo(1,ipos)
          endif
          lfapo(1,ipos)=0_ip
          lfapo(3,ipos)=0_ip
          irow=ipos
          if(iposn/=0)then
             lfapo(3,iposn)=2_ip
          endif
          exit
       else
          if(lfapo(1,ipos)==ifront)then
             jpos=ipos
             jplace=1_ip
          else if(lfapo(2,ipos)==ifront)then
             jpos=ipos
             jplace=2_ip
          endif

          iposn=ipos 
          ipos=-nfface
       endif

    enddo
    !
    !     Do we have an empty row
    !
    if(irow/=0)then 

       nfahol=nfahol+1
       call memrea(nfahol,memor_msh,'LFHOLE','delfac',lfahol)
       lfahol(nfahol)=irow
       nfapo=nfapo-1 
       ipos=lpofa(ip1old)
       if(lfapo(3,ipos)==0)then
          lpofa(ip1old)=-1_ip
       endif

    endif
    !
    !     Delete ip2old
    !
    ipos=lpofa(ip2old)
    iposn=0_ip
    jpos=0_ip
    jplace=0_ip
    irow=0_ip

    do
       nfface=lfapo(3,ipos)
       if(nfface==2)then
          if(lfapo(1,ipos)==ifront)then
             lfapo(1,ipos)=lfapo(2,ipos)
          else if(lfapo(2,ipos)==ifront)then   
          else
             lfapo(jplace,jpos)=lfapo(2,ipos)
          endif
          lfapo(2,ipos)=0_ip
          lfapo(3,ipos)=1_ip
          exit
       else if(nfface==1)then 
          if(lfapo(1,ipos)==ifront)then 
          else
             lfapo(jplace,jpos)=lfapo(1,ipos)
          endif
          lfapo(1,ipos)=0_ip
          lfapo(3,ipos)=0_ip
          irow=ipos
          if(iposn/=0)then
             lfapo(3,iposn)=2_ip
          endif
          exit
       else
          if(lfapo(1,ipos)==ifront)then
             jpos=ipos
             jplace=1_ip
          else if(lfapo(2,ipos)==ifront)then
             jpos=ipos
             jplace=2_ip
          endif

          iposn=ipos
          ipos=-nfface
       endif

    enddo
    !
    !     Do we have an empty row
    !
    if(irow/=0)then 

       nfahol=nfahol+1
       call memrea(nfahol,memor_msh,'LFAHOL','delfac',lfahol)
       lfahol(nfahol)=irow
       nfapo=nfapo-1
       ipos=lpofa(ip2old)
       if(lfapo(3,ipos)==0)then
          lpofa(ip2old)=-1_ip
       endif

    endif
    !
    !     Delete ifront in lfront 
    !
    nfhole=nfhole+1
    call memrea(nfhole,memor_msh,'LFHOLE','delfac',lfhole)
    lfhole(nfhole)=ifront

    nfront=nfront-1    

  end subroutine delfac

  subroutine addfac(nfront,nfapo,ipa,ipb,nfahol,coor,ndim,npoin,nfhole,nheap,&
       lfront,lpofa,lfapo,rfront,lheap,lfahol,lfhole)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)     :: ipa,ipb,ndim,npoin
    real(rp), intent(in)        :: coor(ndim,npoin)
    integer(ip), intent(inout)  :: nfapo,nfront,nfahol,nfhole,nheap
    integer(ip), pointer        :: lfront(:,:),lpofa(:),lfapo(:,:),lfahol(:)
    integer(ip), pointer        :: lheap(:),lfhole(:)
    real(rp), pointer           :: rfront(:)
    integer(ip)                 :: ipos,nfface,nfnew,iposn,iposf
    real(rp)                    :: rx,ry,rz,rlen
    !
    !     This subroutine adds a face to the front and to lfapo
    !
    !
    !     Do we have holes in the datastructure?
    ! 
    if(nfhole>0)then

       iposf=lfhole(nfhole)
       nfhole=nfhole-1

    else 

       iposf=nfront+1 
       call memrea(iposf,memor_msh,'LFRONT','addfac',lfront)
       call memrea(iposf,memor_msh,'RFRONT','addfac',rfront)
       call memrea(iposf,memor_msh,'LHEAP','addfac',lheap)

    endif

    nfront=nfront+1
    lfront(1,iposf)=ipa
    lfront(2,iposf)=ipb
    !
    !     Compute face length
    !
    rx=coor(1,ipb)-coor(1,ipa)
    ry=coor(2,ipb)-coor(2,ipa)
    rz=coor(3,ipb)-coor(3,ipa)
    rlen=sqrt(rx*rx+ry*ry+rz*rz)
    rfront(iposf)=rlen
    !
    !     Introduce face in heap
    !
    call addheap(lheap,nheap,rfront,iposf,lfront,nfront+nfhole)
    !
    !     Update face to point datastructure
    !
    ipos=lpofa(ipa)
    !
    !     Do we have already faces attached to ipa?
    !
    if(ipos==0)then
       if(nfahol>0)then
          iposn=lfahol(nfahol)
          nfahol=nfahol-1
       else
          iposn=nfapo+1
          call memrea(iposn,memor_msh,'LFAPO','addfac',lfapo)
       endif

       nfapo=nfapo+1
       lpofa(ipa)=iposn
       lfapo(1,iposn)=iposf
       lfapo(3,iposn)=1_ip

    else

       do 
          nfface=lfapo(3,ipos)
          if(nfface==2)then
             if(nfahol>0)then
                iposn=lfahol(nfahol)
                nfahol=nfahol-1
             else
                iposn=nfapo+1
                call memrea(iposn,memor_msh,'LFAPO','addfac',lfapo)
             endif

             nfapo=nfapo+1
             lfapo(1,iposn)=iposf
             lfapo(3,iposn)=1_ip
             lfapo(3,ipos)=-iposn
             exit
          else if(nfface==1)then
             lfapo(2,ipos)=iposf
             lfapo(3,ipos)=2_ip
             exit
          else  
             ipos=-nfface
          endif
       enddo

    endif

    ipos=lpofa(ipb)

    if(ipos==0)then
       if(nfahol>0)then
          iposn=lfahol(nfahol)
          nfahol=nfahol-1
       else
          iposn=nfapo+1
          call memrea(iposn,memor_msh,'LFAPO','addfac',lfapo)
       endif

       nfapo=nfapo+1
       lpofa(ipb)=iposn
       lfapo(1,iposn)=iposf
       lfapo(3,iposn)=1_ip

    else

       do 
          nfface=lfapo(3,ipos)
          if(nfface==2)then
             if(nfahol>0)then
                iposn=lfahol(nfahol)
                nfahol=nfahol-1
             else
                iposn=nfapo+1
                call memrea(iposn,memor_msh,'LFAPO','addfac',lfapo)
             endif

             nfapo=nfapo+1
             lfapo(1,iposn)=iposf
             lfapo(3,iposn)=1_ip
             lfapo(3,ipos)=-iposn
             exit
          else if(nfface==1)then
             lfapo(2,ipos)=iposf
             lfapo(3,ipos)=2_ip
             exit
          else  
             ipos=-nfface
          endif
       enddo
    endif

  end subroutine addfac

  subroutine swpahd(ienew,iview,lface,eltoel,ndim,npoin,coor,nnofa,nface,rnopo,lfmark,lptri,rnofa)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: ienew,ndim,npoin,iview,nface,nnofa
    real(rp), intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin)
    real(rp), intent(inout)    :: rnofa(nnofa,nface)
    integer(ip), intent(in)    :: lfmark(nface) 
    integer(ip), intent(inout) :: lface(nnofa,nface),eltoel(nnofa,nface),lptri(npoin)
    integer(ip)                :: ineigh,iopt,ielem,iview1 
    integer(ip)                :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))

    ineigh=eltoel(ltab(1,iview),ienew)
    if(ineigh/=0)then
       if(lfmark(ineigh)==0)then
          ielem=ineigh
          if(eltoel(1,ielem)==ienew)then
             iview1=1_ip
          else if(eltoel(2,ielem)==ienew)then
             iview1=2_ip
          else    
             iview1=3_ip
          endif
          ineigh=eltoel(ltab(1,iview1),ielem)
          if(lfmark(ineigh)==0)then
             !
             !     Swap is controlled by quality --> iqual=1
             !
             call swapf(ielem,ineigh,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,iopt,1_ip,lptri,rnofa)
          endif
       endif
    endif

    ineigh=eltoel(ltab(2,iview),ienew)
    if(ineigh/=0)then
       if(lfmark(ineigh)==0)then
          ielem=ineigh
          if(eltoel(1,ielem)==ienew)then
             iview1=1_ip
          else if(eltoel(2,ielem)==ienew)then
             iview1=2_ip
          else    
             iview1=3_ip
          endif
          ineigh=eltoel(ltab(2,iview1),ielem)
          if(lfmark(ineigh)==0)then
             !
             !     Swap is controlled by quality --> iqual=1
             !
             call swapf(ielem,ineigh,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,iopt,1_ip,lptri,rnofa)
          endif
       endif
    endif

  end subroutine swpahd

  subroutine chkarray(lpsur,npoin,rnofa,rnopo,nnofa,nface,ndim,rnofaold,nfold)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)     :: npoin,nnofa,nface,ndim,nfold
    integer(ip),intent(in)     :: lpsur(npoin)
    real(rp),intent(in)        :: rnofa(nnofa,nface),rnopo(ndim,npoin)
    real(rp),intent(in)        :: rnofaold(nnofa,nfold)

  end subroutine chkarray


  subroutine dbgfront(nface,lface,nnofa,nfront,lfront,eltoel,lfmark,lfhole,nfhole,nfahol,&
       lptri,npoin,ierr,ptosi1,ptosi2,lside,nnosi,nside,nehole,lehole)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nface,nnofa,nfront,nfahol,npoin,nfhole,nside,nnosi,nehole
    integer(ip),intent(in)     :: lface(nnofa,nface),lfront(3,nfront+nfhole),lehole(nehole)
    integer(ip),intent(in)     :: ptosi2(npoin+1),ptosi1(*),lside(nnosi,nside)
    integer(ip),intent(in)     :: eltoel(nnofa,nface),lfmark(nface),lfhole(nfhole+nfahol),lptri(npoin) 
    integer(ip),intent(inout)  :: ierr 
    integer(ip)                :: iface,ifacen,j,ip1,ip2,ifront,ifhole,ichk,ineigh,ifirst,iview,iturn
    integer(ip)                :: ipa,ipb,iside,isto,ifound,iexist,iehole
    integer(ip),pointer        :: lmark(:),lemark(:) 
    integer(4)                 :: istat 
    integer(ip)                :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))

    allocate(lmark(nfront+nfhole),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','dbgfront',lmark) 
    allocate(lemark(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LEMARK','dbgfront',lemark)
    !
    !     Mark the deleted sides of the front
    ! 
    do ifhole=1,nfhole
       lmark(lfhole(ifhole))=1_ip
    enddo
    !
    !     Mark the deleted faces 
    ! 
    do iehole=1,nehole
       lemark(lehole(ifhole))=1_ip
    enddo
    !
    !     Check that a front side is at the edge of a marked/non marked face 
    !
    do iface=1,nface
       if(lemark(iface)==1)cycle
       if(lfmark(iface)==1)then
          do j=1,nnofa
             ineigh=eltoel(j,iface)
             if(ineigh/=0)then
                if(lfmark(ineigh)==0)then
                   ip1=lface(ltab(1,j),iface)
                   ip2=lface(ltab(2,j),iface)

                   ichk=0_ip
                   do ifront=1,nfront+nfhole
                      if(lmark(ifront)==0)then
                         if(ip2==lfront(1,ifront) .and. ip1==lfront(2,ifront))then
                            ichk=1_ip
                            exit
                         endif
                      endif
                   enddo

                   if(ichk==0)then
                      write(*,*)'Error in dbgfront'
                      write(*,*)'iface=',iface
                      write(*,*)'ip1=',ip1,'ip2=',ip2
                      ierr=1
                      return
                   endif
                endif
             endif
          enddo
       endif
    enddo
    !
    !     Check that the front is well oriented
    !
    do ifront=1,nfront+nfhole
       if(lmark(ifront)==0)then
          ip1=lfront(1,ifront)
          ip2=lfront(2,ifront)

          iexist=0

          iface=lptri(ip1)
          ifirst=iface 
          iturn=1
          do

             if(lface(1,iface)==ip1)then
                iview=1_ip 
             else if(lface(2,iface)==ip1)then
                iview=2_ip 
             else       
                iview=3_ip 
             endif


             if(lface(ltab(1,iview),iface)==ip2)then
                !
                !     The front side exists in the mesh
                !
                iexist=1
                !
                !     iface is in front of the front, lfmark should not 
                !     be marked
                !
                if(lfmark(iface)==1)then
                   write(*,*)'Error with lfmark 1'
                   write(*,*)'ip1=',ip1,'ip2=',ip2
                   ierr=1
                   return
                endif
                !
                !     iface is behind the front, lfmark should  
                !     be marked
                !
                ineigh=eltoel(ltab(2,iview),iface)
                if(ineigh/=0)then
                   if(lfmark(ineigh)/=1)then
                      !
                      !     It may be a ridge inside the face
                      !
                      ifound=0_ip
                      do isto=ptosi2(ip1),ptosi2(ip1+1)-1
                         iside=ptosi1(isto)
                         ipa=min(lside(1,iside),lside(2,iside))
                         ipb=max(lside(1,iside),lside(2,iside))
                         if(ipa==min(ip1,ip2) .and. ipb==max(ip1,ip2))then 
                            ifound=1_ip
                            exit
                         endif
                      enddo

                      if(ifound==0)then
                         write(*,*)'Error with lfmark 2'
                         ierr=1
                         return
                      endif
                   endif
                endif
             endif

             if(iturn==1)then
                ifacen=eltoel(ltab(1,iview),iface)
                if(ifacen/=0)then
                   iface=ifacen
                else
                   iface=eltoel(ltab(2,iview),iface) 
                   if(iface==0)exit
                   iturn=0
                endif
             else
                ifacen=eltoel(ltab(2,iview),iface)
                if(ifacen/=0)then
                   iface=ifacen
                else
                   exit
                endif
             endif

             if(iface==ifirst .and. iturn==1)exit

          enddo
          !
          !     Check that the front was found
          !
          if(iexist==0)then
             write(*,*)'Front side not found, surface swapped?'
             ierr=1
             return
          endif

       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'LEMARK','dbgfront',lemark)
    deallocate(lemark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEMARK','dbgfront',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','dbgfront',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','dbgfront',0_ip)

  end subroutine dbgfront

  subroutine swaploc(ifirst,lface,nnofa,nface,coor,npoin,ndim,rnopo,lelem,lfmark,eltoel,lptri,rnofa,&
       lside,nnosi,nside,ptosi1,ptosi2,npcusp)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim,ifirst,nnosi,nside,npcusp
    integer(ip), intent(in)      :: lfmark(nface)
    integer(ip), intent(in)      :: lside(nnosi,nside),ptosi2(npcusp+1),ptosi1(*)
    integer(ip), intent(inout)   :: lface(nnofa,nface),eltoel(nnofa,nface),lptri(npoin),lelem(nface)
    real(rp), intent(in)         :: coor(ndim,npoin),rnopo(ndim,npoin)
    real(rp), intent(inout)      :: rnofa(ndim,nface)
    integer(ip)                  :: iface,jface,iopt,j,istack,lstack(100),jstack,nstack,ineigh,isto,iside
    integer(ip)                  :: maxlevel,level(100),ilevel,ioptglo,ioptel,ifound,ipa,ipb,ip1,ip2,ipcusp
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     Set the level of neighboring elements to swap
    !     In confined places, 4 is too low 
    !
    maxlevel=6
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
       iface=lstack(istack)
       ilevel=level(istack)

       if(ilevel==maxlevel)exit

       do j=1,nnofa
          ineigh=eltoel(j,iface)
          !
          !     Is there a neighbor?
          !
          if(ineigh/=0)then
             !
             !     Are we crossing the front?
             ! 
             if(lfmark(ineigh)==0)then
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
          endif
       enddo
    enddo
    !
    !     Delete the first element from the stack
    !
    lelem(ifirst)=0_ip
    !
    !     Swap them
    !
    do 
       ioptglo=0_ip
       !
       !     Do not consider the first element
       !
       do istack=2,nstack

          iface=lstack(istack)
          if(lelem(iface)==1)then
             ioptel=0_ip
             do j=1,nnofa
                !
                !     Do we have inner sides?
                !
                ip1=min(lface(ltab(1,j),iface),lface(ltab(2,j),iface))
                ip2=max(lface(ltab(1,j),iface),lface(ltab(2,j),iface))

                if(ip1<=npcusp .or. ip2<=npcusp)then

                   if(ip1<=npcusp)then
                      ipcusp=ip1
                   else
                      ipcusp=ip2
                   endif

                   ifound=0
                   do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
                      iside=ptosi1(isto)
                      ipa=min(lside(1,iside),lside(2,iside)) 
                      ipb=max(lside(1,iside),lside(2,iside)) 
                      if(ip1==ipa .and. ip2==ipb)then
                         ifound=1
                         exit
                      endif
                   enddo

                   if(ifound==1)cycle

                endif

                jface=eltoel(j,iface)
                if(jface/=0)then
                   if(lelem(jface)==1)then
                      iopt=0_ip
                      !
                      !     Swap is controlled by quality --> iqual=1
                      !
                      call swapf(iface,jface,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,iopt,1_ip,lptri,rnofa)
                      !call outerror4(lface,nnofa,nface,coor,npoin,ndim)
                      if(iopt==1)then
                         ioptglo=1_ip
                         ioptel=1_ip 
                         !
                         !     Mark jface and the neighbors of iface and jface if they belong to the patch
                         !
                         ineigh=eltoel(1,iface)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(2,iface)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(3,iface)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(1,jface)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(2,jface)
                         if(ineigh/=0)then
                            if(lelem(ineigh)>1)then
                               lelem(ineigh)=1_ip
                            endif
                         endif
                         ineigh=eltoel(3,jface)
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
                lelem(iface)=2_ip
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

  end subroutine swaploc

  subroutine chkrnofa(lface,nface,nnofa,ndim,coor,npoin,rnofa)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim
    integer(ip), intent(in)      :: lface(nnofa,nface)
    real(rp), intent(in)         :: coor(ndim,npoin),rnofa(ndim,nface)
    integer(ip)                  :: iface
    real(rp)                     :: rface(3),modul,rx,ry,rz,rl

    do iface=1,nface

       call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rface,iface,modul)

       rx=rnofa(1,iface)-rface(1)
       ry=rnofa(2,iface)-rface(2)
       rz=rnofa(3,iface)-rface(3)
       rl=sqrt(rx*rx+ry*ry+rz*rz)

       if(rl>1e-9)then
          write(*,*)'Error chkrnofa, iface=',iface
          write(*,*)rface(1),rface(2),rface(3)
          write(*,*)rnofa(1,iface),rnofa(2,iface),rnofa(3,iface)
          stop

       endif

    enddo

  endsubroutine chkrnofa

  subroutine cross(ipa,ipb,lface,nnofa,nface,lptri,coor,npoin,ndim,eltoel,rnopo,rnofa,lfmark,&
       ienew,d1,d2,d3,lmark,lpofa,lfapo,nfapo,nfahol,lfront,nfront,ipcros,ipfirst,rsize,lelem,&
       ierr,lside,nnosi,nside,npcusp,ptosi1,ptosi2)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: ipa,ipb,nnofa,nface,npoin,ndim,nfront,ipfirst,npcusp,nnosi,nside
    integer(ip), intent(inout) :: nfapo,nfahol,ienew,ipcros 
    integer(ip), intent(in)    :: lface(nnofa,nface),eltoel(nnofa,nface),lptri(npoin),lfmark(nface)
    integer(ip), intent(in)    :: lside(nnosi,nside),ptosi2(npcusp+1),ptosi1(*)
    integer(ip), intent(inout) :: lmark(npoin),lpofa(npoin),lfapo(nfapo+nfahol) 
    integer(ip), intent(inout) :: lfront(3,nfront),lelem(nface),ierr 
    real(rp), intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin),rnofa(nnofa,nface)
    real(rp), intent(inout)    :: d1,d2,d3
    real(rp), intent(in)       :: rsize(npoin)
    integer(ip)                :: ie,j,ip1,ip2,ip3,ienext,ielem2,iview1,ifirst,npf,lpf(100),ieold
    real(rp)                   :: rface(ndim),rx,ry,rz,rscal,pproj(ndim),p1(ndim),p2(ndim),p3(ndim)
    real(rp)                   :: epsil,rx1,ry1,rz1,rx2,ry2,rz2,rx3,ry3,rz3,rlen1,rlen2,rlen3,tollen
    real(rp)                   :: tolgeo,q1,qmax,tollen2,tolside,rnli,rtx,rty,rtz,rtl,rnx,rny,rnz  
    real(rp)                   :: rnopx,rnopy,rnopz  
    real(rp)                   :: modul,c05,qgeo,rqual(100),rswap,dtot,rnl,rnl2,c10,csca,pnew(ndim),c20,qbig,c00,cscal
    integer(ip)                :: lpstack1(100),lpstack2(100),mstack1,mstack2,npstack1,npstack2,ipstack
    integer(ip)                :: ipoin,lfacn(3),ireach,ipnew,jpstack,jpoin,lpopt(100),npopt,nitermax,iter
    integer(ip)                :: ienext1,ienext2,iturn,ipclos,ineigh,iback,nestack,lestack(100),iestack,ichk
    integer(ip)                :: ipmin,ipmax,iside,isto,ipcusp,ipaa,ipbb,iface,inofa
    integer(ip)                :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This sub is responsible for expanding the front
    !
    !     We try to insert the optimal point  
    !
    !              - If we reach this point:
    !                       without crossing the front, 
    !                       without crossing an inner side, 
    !                       without being too close from another point, 
    !                       and without being too close from a side of the front                  
    !                  then:
    !                       the point is accepted
    !                  else:
    !                       we build a list of close point and choose 
    !                       the one with the best quality and go back searching
    !
    !              - If we find a point too close, we stop the progression
    !              - If we cross the front, we stop the progression
    !              - If we cross an inner side, we stop the progression
    !              - If the point is too close from a side of the front, we stop the progression
    !
    !
    !     We mark lmark with: 1 if it is in lpstack2       
    !                         2 if it is in lpstack1
    !                         3 if it is in lpstack1 and it has stopped a progression
    !
    epsil=1.0d-09
    tollen=sqrt(2.0d+00)
    tollen2=1.0d+00/sqrt(2.0d+00)
    tolside=0.3d+00 
    tolgeo=0.0d+00
    c05=0.5d+00
    c00=0.0d+00
    c10=1.0d+00
    c20=2.0d+00
    npf=0_ip
    qbig=1.0d+12
    nitermax=10000

    mstack1=100
    mstack2=100
    !
    !     Mark the points of the current edge of the front
    !
    lmark(ipa)=1_ip
    lmark(ipb)=1_ip
    !
    !     Initialize the stacks
    !     lpstack2 will contain the points marked in lmark
    !     lpstack1 will contain the close points candidates for best point  
    !     lestack will contain the faces visited
    !
    npstack1=0_ip
    npstack2=0_ip
    nestack=0_ip
    !
    !     Initialize the point to check by ipfirst
    !
    ipnew=ipfirst
    !
    !     Initialize the number of optimal points found
    !
    npopt=0_ip
    !
    !    Check if some edges of the front will be crossed moving 
    !    from ipa or ipb to ipnew
    !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Find the  element containing edge (ipa,ipb)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ie=lptri(ipa)
    iturn=1_ip
    if(ie<=0 .or. ie>nface)then
       write(*,*)'Error in cross 1'
       write(*,*)ie,nface
       stop
    endif

    iter=0
    do 

       iter=iter+1
       if(iter==nitermax)then
          !
          !     Surfaces may be disconnected
          !
          do iface=1,nface           
             do inofa=1,nnofa 
                if(lface(ltab(1,inofa),iface)==ipa .and. lface(ltab(2,inofa),iface)==ipb)then
                   ie=iface
                   goto 33
                endif
             enddo
          enddo
          write(*,*)'Error in cross finding ipa,ipb'
          ierr=1
          return
       endif

       if(lface(1,ie)==ipa)then
          ienext1=eltoel(2,ie)
          ienext2=eltoel(3,ie)
          iview1=1_ip
       else if(lface(2,ie)==ipa)then
          ienext1=eltoel(3,ie)
          ienext2=eltoel(1,ie)
          iview1=2_ip
       else
          ienext1=eltoel(1,ie)
          ienext2=eltoel(2,ie)
          iview1=3_ip
       endif

       if(ipb==lface(ltab(1,iview1),ie))goto 33

       if(iturn==1)then

          if(ienext1/=0)then
             ie=ienext1
          else
             ie=ienext2
             if(ie==0)exit
             iturn=0_ip
          endif
       else

          ie=ienext2
          if(ie==0)exit

       endif


    enddo
    !
    !     Did we find ie?
    !
    if(ie==0)then
       !
       !     Surfaces may be disconnected
       !
       do iface=1,nface           
          do inofa=1,nnofa 
             if(lface(ltab(1,inofa),iface)==ipa .and. lface(ltab(2,inofa),iface)==ipb)then
                ie=iface
                goto 33
             endif
          enddo
       enddo
       write(*,*)'Error in cross finding ipa,ipb'
       ierr=1
       return

    endif
    !
    !     Remember it
    !
33  continue
    ifirst=ie

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Loop tentatively to add a new point
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do
       !
       !     Clean up lmark
       !
       do ipstack=1,npstack2
          lmark(lpstack2(ipstack))=0_ip
       enddo
       lmark(ipa)=2_ip
       lmark(ipb)=2_ip
       npstack1=0_ip
       npstack2=0_ip
       !
       !     Clean up lelem
       !
       do iestack=1,nestack
          lelem(lestack(iestack))=0_ip
       enddo
       nestack=0_ip
       !
       !     Initialize pnew
       !
       pnew(1)=coor(1,ipnew)
       pnew(2)=coor(2,ipnew)
       pnew(3)=coor(3,ipnew)
       !
       !     Initialize the first element with ifirst
       !
       ie=ifirst
       !
       !     Initialize the reach marker to sucessfull
       !
       ireach=1_ip
       !
       !     Initialize the back marker to sucessfull
       !
       iback=0_ip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !     Find the element of the ball of ipa crossed by edge (ipa,ipnew)
       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       rtx=pnew(1)-coor(1,ipa)
       rty=pnew(2)-coor(2,ipa)
       rtz=pnew(3)-coor(3,ipa)
       rtl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
       rtl=c10/rtl
       rtx=rtx*rtl
       rty=rty*rtl
       rtz=rtz*rtl
       !
       !     Initialize the normal to evaluate intersection
       !
       rnopx=rnopo(1,ipa)+rnopo(1,ipnew)
       rnopy=rnopo(2,ipa)+rnopo(2,ipnew)
       rnopz=rnopo(3,ipa)+rnopo(3,ipnew)
       rnl=sqrt(rnopx*rnopx+rnopy*rnopy+rnopz*rnopz)
       rnl=c10/rnl
       rnopx=rnopx*rnl 
       rnopy=rnopy*rnl 
       rnopz=rnopz*rnl 

       do
          !
          !    Are we still on the surface ?
          !
          if(ie==0)then
             ireach=0_ip
             goto 10
          endif
          !
          !     Are we going back ?
          !
          if(lelem(ie)==1)then
             goto 10  
          endif
          lelem(ie)=1_ip
          nestack=nestack+1
          lestack(nestack)=ie

          if(lface(1,ie)==ipa)then
             iview1=1_ip
          else if(lface(2,ie)==ipa)then
             iview1=2_ip
          else
             iview1=3_ip
          endif

          ip1=ipa
          ip2=lface(ltab(1,iview1),ie)
          ip3=lface(ltab(2,iview1),ie)
          !
          !     Does ipnew belong to this element ?
          !
          if(ipnew==ip2)goto 10
          if(ipnew==ip3)goto 10
          !
          !     Has this point been already checked ?    
          !
          if(lmark(ip2)==0)then 
             lmark(ip2)=1_ip
             npstack2=npstack2+1
             lpstack2(npstack2)=ip2
             !
             !     Check distance to ip2
             !
             call length(ip2,ipnew,rsize,coor,rlen2)
             !
             !     Is the point too close ?
             !
             if(rlen2<tollen2)then
                ireach=0_ip
                lmark(ip2)=3_ip
                npstack1=npstack1+1
                lpstack1(npstack1)=ip2
                goto 10
             endif
             !
             !     Is the point in the local box ?
             !
             if(rlen2<tollen)then
                !
                !     Is this point a point of the current front ?
                !     If yes, add all the points of the current front inside the local box
                !
                if(lpofa(ip2)>0)then
                   call chkboun(ip2,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                        nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                endif
             endif
          endif
          !
          !     Has this point been already checked ?    
          !
          if(lmark(ip3)==0)then 
             lmark(ip3)=1_ip
             npstack2=npstack2+1
             lpstack2(npstack2)=ip3
             !
             !     Check distance to ip3
             !
             call length(ip3,ipnew,rsize,coor,rlen3)
             !
             !     Is the point too close ?
             !
             if(rlen3<tollen2)then
                ireach=0_ip
                lmark(ip3)=3_ip
                npstack1=npstack1+1
                lpstack1(npstack1)=ip3
                goto 10
             endif
             !
             !     Is the point in the local box ?
             !
             if(rlen3<tollen)then
                !
                !     Is this point a point of the current front ?
                !     If yes, add all the points of the current front inside the local box
                !
                if(lpofa(ip3)>0)then
                   call chkboun(ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                        nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                endif
             endif
          endif
          !
          !     Find which side is crossed
          ! 
          rface(1)=rnofa(1,ie)
          rface(2)=rnofa(2,ie)
          rface(3)=rnofa(3,ie)

          rx=pnew(1)-coor(1,ip1)
          ry=pnew(2)-coor(2,ip1)
          rz=pnew(3)-coor(3,ip1)

          rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
          pproj(1)= pnew(1)-rscal*rface(1)
          pproj(2)= pnew(2)-rscal*rface(2)
          pproj(3)= pnew(3)-rscal*rface(3)
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
          !     Side (ip1,ip3)
          !
          p1(1)=coor(1,ip3)-pproj(1)
          p1(2)=coor(2,ip3)-pproj(2)
          p1(3)=coor(3,ip3)-pproj(3)
          p2(1)=coor(1,ip1)-pproj(1)
          p2(2)=coor(2,ip1)-pproj(2)
          p2(3)=coor(3,ip1)-pproj(3)

          call orient3D(p1,p2,rface,d1,ndim)
          d1=d1/dtot
          !
          !     Compute the visibility by checking if edge (ip1,ip3) crosses the plane (ip1,pnew,rface(iface))
          ! 
          !     Compute normal to plane  
          !
          !rnx= rface(2)*rtz-rface(3)*rty
          !rny=-rface(1)*rtz+rface(3)*rtx
          !rnz= rface(1)*rty-rface(2)*rtx
          rnx= rnopy*rtz-rnopz*rty
          rny=-rnopx*rtz+rnopz*rtx
          rnz= rnopx*rty-rnopy*rtx
          rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
          rnl=c10/rnl
          rnx=rnx*rnl
          rny=rny*rnl
          rnz=rnz*rnl

          rx=coor(1,ip3)-coor(1,ipa)
          ry=coor(2,ip3)-coor(2,ipa)
          rz=coor(3,ip3)-coor(3,ipa)

          csca=rnx*rx+rny*ry+rnz*rz
          !
          !    And the test
          ! 
          !if(d1<-epsil)then
          if(csca<-epsil)then
             ie=eltoel(ltab(1,iview1),ie)
             !
             !     Did we cross the border ?
             !  
             if(ie==0)then
                ireach=0_ip
                call chkboun2(ip1,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                     nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                goto 10 
             endif
             !
             !     Did we cross the front ?
             !
             if(lfmark(ie)==1)then
                ireach=0_ip
                call chkboun2(ip1,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                     nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                goto 10 
             endif
             !
             !     Did we cross an inner side ?
             !
             ipmin=min(ip1,ip3)
             ipmax=max(ip1,ip3)

             if(ipmin<=npcusp .or. ipmax<=npcusp)then

                if(ipmin<=npcusp)then
                   ipcusp=ipmin
                else
                   ipcusp=ipmax
                endif

                do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
                   iside=ptosi1(isto)
                   ipaa=min(lside(1,iside),lside(2,iside)) 
                   ipbb=max(lside(1,iside),lside(2,iside)) 
                   if(ipmin==ipaa .and. ipmax==ipbb)then
                      ireach=0_ip
                      call chkboun2(ip1,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                           nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                      goto 10 
                   endif
                enddo
             endif

             cycle
          endif
          !
          !     Side (ip1,ip2)
          !
          p1(1)=coor(1,ip2)-pproj(1)
          p1(2)=coor(2,ip2)-pproj(2)
          p1(3)=coor(3,ip2)-pproj(3)
          p2(1)=coor(1,ip1)-pproj(1)
          p2(2)=coor(2,ip1)-pproj(2)
          p2(3)=coor(3,ip1)-pproj(3)

          call orient3D(p1,p2,rface,d2,ndim)
          d2=d2/dtot
          ! 
          !     Compute normal to plane  
          !
          rx=coor(1,ip2)-coor(1,ipa)
          ry=coor(2,ip2)-coor(2,ipa)
          rz=coor(3,ip2)-coor(3,ipa)

          csca=rnx*rx+rny*ry+rnz*rz
          !
          !     And the test
          !
          !if(d2>epsil)then
          if(csca>epsil)then
             ie=eltoel(ltab(2,iview1),ie)
             !
             !     Did we cross the border ?
             !
             if(ie==0)then
                ireach=0_ip
                call chkboun2(ip1,ip2,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                     nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                goto 10 
             endif
             !
             !     Did we cross the front ?
             !
             if(lfmark(ie)==1)then
                ireach=0_ip
                call chkboun2(ip1,ip2,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                     nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                goto 10 
             endif
             !
             !     Did we cross an inner side ?
             !
             ipmin=min(ip1,ip2)
             ipmax=max(ip1,ip2)

             if(ipmin<=npcusp .or. ipmax<=npcusp)then

                if(ipmin<=npcusp)then
                   ipcusp=ipmin
                else
                   ipcusp=ipmax
                endif

                do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
                   iside=ptosi1(isto)
                   ipaa=min(lside(1,iside),lside(2,iside)) 
                   ipbb=max(lside(1,iside),lside(2,iside)) 
                   if(ipmin==ipaa .and. ipmax==ipbb)then
                      ireach=0_ip
                      call chkboun2(ip1,ip2,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                           nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                      goto 10 
                   endif
                enddo
             endif

             cycle
          endif
          !
          !     Element found 
          !
          exit
       enddo
       !
       !     Are we already in the final element?
       !
       p1(1)=coor(1,ip2)-pproj(1)
       p1(2)=coor(2,ip2)-pproj(2)
       p1(3)=coor(3,ip2)-pproj(3)
       p2(1)=coor(1,ip3)-pproj(1)
       p2(2)=coor(2,ip3)-pproj(2)
       p2(3)=coor(3,ip3)-pproj(3)

       call orient3D(p1,p2,rface,d3,ndim)
       d3=d3/dtot

       if(d3>epsil)then
          goto 10
       endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !     Go through the pipe
       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ieold=ie
       ie=eltoel(iview1,ie)

       do
          !
          !    Are we still on the surface
          !
          if(ie==0)then
             ireach=0_ip
             goto 10
          endif

          if(eltoel(1,ie)==ieold)then
             iview1=1_ip
          else if(eltoel(2,ie)==ieold)then
             iview1=2_ip
          else
             iview1=3_ip
          endif

          ip1=lface(iview1,ie)
          ip2=lface(ltab(1,iview1),ie)
          ip3=lface(ltab(2,iview1),ie)
          !
          !     Did we cross the front ?
          !
          if(lfmark(ie)==1)then
             ireach=0_ip
             call chkboun2(ip2,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,&
                  nfapo,nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
             goto 10 
          endif
          !
          !     Did we cross an inner side ?
          !
          ipmin=min(ip2,ip3)
          ipmax=max(ip2,ip3)

          if(ipmin<=npcusp .or. ipmax<=npcusp)then

             if(ipmin<=npcusp)then
                ipcusp=ipmin
             else
                ipcusp=ipmax
             endif

             do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
                iside=ptosi1(isto)
                ipaa=min(lside(1,iside),lside(2,iside)) 
                ipbb=max(lside(1,iside),lside(2,iside)) 
                if(ipmin==ipaa .and. ipmax==ipbb)then
                   ireach=0_ip
                   call chkboun2(ip2,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                        nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                   goto 10 
                endif
             enddo
          endif
          !
          !     Are we going back ?
          !
          if(lelem(ie)==1)then
             goto 10  
          endif
          lelem(ie)=1_ip
          nestack=nestack+1
          lestack(nestack)=ie
          !
          !     Does ipnew belong to this element
          !
          if(ipnew==ip1)goto 10
          !
          !     Check if the point has already been marked, if it belongs to the front, and the distances
          !
          if(lmark(ip1)==0)then
             lmark(ip1)=1_ip
             npstack2=npstack2+1
             lpstack2(npstack2)=ip1
             call length(ip1,ipnew,rsize,coor,rlen1)
             !
             !     Is the point too close ?
             !
             if(rlen1<tollen2)then
                ireach=0_ip
                lmark(ip1)=3_ip
                npstack1=npstack1+1
                lpstack1(npstack1)=ip1
                goto 10
             endif
             !
             !     Is the point in the local box ?
             !
             if(rlen1<tollen)then
                if(lpofa(ip1)>0)then
                   call chkboun(ip1,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                        nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                endif
             endif
          endif
          !
          !     Find which side is crossed
          ! 
          rface(1)=rnofa(1,ie)
          rface(2)=rnofa(2,ie)
          rface(3)=rnofa(3,ie)

          rx=pnew(1)-coor(1,ip1)
          ry=pnew(2)-coor(2,ip1)
          rz=pnew(3)-coor(3,ip1)

          rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
          pproj(1)= pnew(1)-rscal*rface(1)
          pproj(2)= pnew(2)-rscal*rface(2)
          pproj(3)= pnew(3)-rscal*rface(3)
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
          !     Side (ip1,ip2)
          !
          p1(1)=coor(1,ip2)-pproj(1)
          p1(2)=coor(2,ip2)-pproj(2)
          p1(3)=coor(3,ip2)-pproj(3)
          p2(1)=coor(1,ip1)-pproj(1)
          p2(2)=coor(2,ip1)-pproj(2)
          p2(3)=coor(3,ip1)-pproj(3)

          call orient3D(p1,p2,rface,d2,ndim)
          d2=d2/dtot
          !
          !     Side (ip1,ip3)
          !
          p1(1)=coor(1,ip1)-pproj(1)
          p1(2)=coor(2,ip1)-pproj(2)
          p1(3)=coor(3,ip1)-pproj(3)
          p2(1)=coor(1,ip3)-pproj(1)
          p2(2)=coor(2,ip3)-pproj(2)
          p2(3)=coor(3,ip3)-pproj(3)

          call orient3D(p1,p2,rface,d3,ndim)
          d3=d3/dtot
          !
          !     Did we reach the last element
          !
          if(d2<-epsil .and. d3<-epsil)then
             !
             !     This is the last element
             !
             exit

          else
             !
             !     Which side between (ip1,ip2) and (ip1,ip3) is intersected by (ipa,pnew)? 
             !
             p1(1)=coor(1,ip1)-coor(1,ipa)
             p1(2)=coor(2,ip1)-coor(2,ipa)
             p1(3)=coor(3,ip1)-coor(3,ipa)
             p2(1)=pproj(1)-coor(1,ipa)
             p2(2)=pproj(2)-coor(2,ipa)
             p2(3)=pproj(3)-coor(3,ipa)

             call orient3D(p1,p2,rface,d2,ndim)
             d2=d2/dtot
             ! 
             !     Compute normal to plane  
             !
             !rnx= rface(2)*rtz-rface(3)*rty
             !rny=-rface(1)*rtz+rface(3)*rtx
             !rnz= rface(1)*rty-rface(2)*rtx
             rnx= rnopy*rtz-rnopz*rty
             rny=-rnopx*rtz+rnopz*rtx
             rnz= rnopx*rty-rnopy*rtx
             rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
             rnl=c10/rnl
             rnx=rnx*rnl
             rny=rny*rnl
             rnz=rnz*rnl

             rx=coor(1,ip1)-coor(1,ipa)
             ry=coor(2,ip1)-coor(2,ipa)
             rz=coor(3,ip1)-coor(3,ipa)

             csca=rx*rnx+ry*rny+rz*rnz
             !
             !     And the test
             !
             !if(d2>epsil)then
             if(csca>epsil)then
                ieold=ie
                ie=eltoel(ltab(1,iview1),ie)
             else
                ieold=ie
                ie=eltoel(ltab(2,iview1),ie)
             endif
          endif

       enddo

10     continue
       !
       !     Initialize with ifirst
       !
       ie=ifirst
       !
       !     Clean up lelem
       !
       do iestack=1,nestack
          lelem(lestack(iestack))=0_ip
       enddo
       nestack=0_ip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !     Find the element of the ball of ipb crossed by edge (ipb,ipnew)
       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       rtx=pnew(1)-coor(1,ipb)
       rty=pnew(2)-coor(2,ipb)
       rtz=pnew(3)-coor(3,ipb)
       rtl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
       rtl=c10/rtl
       rtx=rtx*rtl
       rty=rty*rtl
       rtz=rtz*rtl
       !
       !     Initialize the normal to evaluate intersection
       !
       rnopx=rnopo(1,ipb)+rnopo(1,ipnew)
       rnopy=rnopo(2,ipb)+rnopo(2,ipnew)
       rnopz=rnopo(3,ipb)+rnopo(3,ipnew)
       rnl=sqrt(rnopx*rnopx+rnopy*rnopy+rnopz*rnopz)
       rnl=c10/rnl
       rnopx=rnopx*rnl 
       rnopy=rnopy*rnl 
       rnopz=rnopz*rnl 

       do
          !
          !    Are we still on the surface
          !
          if(ie==0)then
             ireach=0_ip
             goto 20
          endif
          !
          !     Are we going back ?
          !
          if(lelem(ie)==1)then
             iback=1_ip
             call rescueside(ie,rnofa,nnofa,nface,d1,d2,d3,pnew,ndim,coor,npoin,lface,ierr)
             if(ierr==1)return   
             goto 20  
          endif
          lelem(ie)=1_ip
          nestack=nestack+1
          lestack(nestack)=ie

          if(lface(1,ie)==ipb)then
             iview1=1_ip
          else if(lface(2,ie)==ipb)then
             iview1=2_ip
          else
             iview1=3_ip
          endif

          ip1=ipb
          ip2=lface(ltab(1,iview1),ie)
          ip3=lface(ltab(2,iview1),ie)
          !
          !     Does ipnew belong to this element
          !
          if(ipnew==ip2)goto 20
          if(ipnew==ip3)goto 20
          !
          !     Has this point been already checked ?    
          !
          if(lmark(ip2)==0)then
             lmark(ip2)=1_ip
             npstack2=npstack2+1
             lpstack2(npstack2)=ip2
             !
             !     Check distance to ip2
             !
             call length(ip2,ipnew,rsize,coor,rlen2)
             !
             !     Is the point too close ?
             !
             if(rlen2<tollen2)then
                ireach=0_ip
                lmark(ip2)=3_ip
                npstack1=npstack1+1
                lpstack1(npstack1)=ip2
                goto 20
             endif
             !
             !     Is the point in the local box ?
             !
             if(rlen2<tollen)then
                if(lpofa(ip2)>0)then
                   !
                   !     Is this point a point of the current front ?
                   !     If yes, add all the points of the current front inside the local box
                   !
                   call chkboun(ip2,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                        nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                endif
             endif
             !
             !     Has this point been reached already in the first part as too close
             !
          else if(lmark(ip2)==3)then
             goto 20
          endif
          !
          !     Has this point been already checked ?    
          !
          if(lmark(ip3)==0)then
             lmark(ip3)=1_ip
             npstack2=npstack2+1
             lpstack2(npstack2)=ip3
             !
             !     Check distance to ip3
             !
             call length(ip3,ipnew,rsize,coor,rlen3)
             !
             !     Is the point too close ?
             !
             if(rlen3<tollen2)then
                ireach=0_ip
                lmark(ip3)=3_ip
                npstack1=npstack1+1
                lpstack1(npstack1)=ip3
                goto 20
             endif
             !
             !     Is the point in the local box ?
             !
             if(rlen3<tollen)then
                if(lpofa(ip3)>0)then
                   !
                   !     Is this point a point of the current front ?
                   !     If yes, add all the points of the current front inside the local box
                   !
                   call chkboun(ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,&
                        nfapo,nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                endif
             endif
             !
             !     Has this point been reached already in the first part as too close
             !
          else if(lmark(ip3)==3)then
             goto 20
          endif
          !
          !     Find which side is crossed
          ! 
          rface(1)=rnofa(1,ie)
          rface(2)=rnofa(2,ie)
          rface(3)=rnofa(3,ie)

          rx=pnew(1)-coor(1,ip1)
          ry=pnew(2)-coor(2,ip1)
          rz=pnew(3)-coor(3,ip1)

          rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
          pproj(1)= pnew(1)-rscal*rface(1)
          pproj(2)= pnew(2)-rscal*rface(2)
          pproj(3)= pnew(3)-rscal*rface(3)
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
          !     Side (ip1,ip2)
          !
          p1(1)=coor(1,ip2)-pproj(1)
          p1(2)=coor(2,ip2)-pproj(2)
          p1(3)=coor(3,ip2)-pproj(3)
          p2(1)=coor(1,ip1)-pproj(1)
          p2(2)=coor(2,ip1)-pproj(2)
          p2(3)=coor(3,ip1)-pproj(3)

          call orient3D(p1,p2,rface,d1,ndim)
          d1=d1/dtot
          !
          !     Compute the visibility by checking if edge (ip1,ip2) crosses the plane (ip1,pnew,rface(iface))
          ! 
          !     Compute normal to plane  
          !
          !rnx= rface(2)*rtz-rface(3)*rty
          !rny=-rface(1)*rtz+rface(3)*rtx
          !rnz= rface(1)*rty-rface(2)*rtx
          rnx= rnopy*rtz-rnopz*rty
          rny=-rnopx*rtz+rnopz*rtx
          rnz= rnopx*rty-rnopy*rtx
          rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
          rnl=c10/rnl
          rnx=rnx*rnl
          rny=rny*rnl
          rnz=rnz*rnl

          rx=coor(1,ip2)-coor(1,ipb)
          ry=coor(2,ip2)-coor(2,ipb)
          rz=coor(3,ip2)-coor(3,ipb)

          csca=rnx*rx+rny*ry+rnz*rz
          !
          !    And the test
          ! 
          !if(d1>epsil)then
          if(csca>epsil)then
             ie=eltoel(ltab(2,iview1),ie)
             !
             !     Did we cross the border ?
             !
             if(ie==0)then
                ireach=0_ip
                call chkboun2(ip1,ip2,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,&
                     nfapo,nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                goto 20 
             endif
             !
             !     Did we cross the front ?
             !
             if(lfmark(ie)==1)then
                ireach=0_ip
                call chkboun2(ip1,ip2,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                     nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                goto 20 
             endif
             !
             !     Did we cross an inner side ?
             !
             ipmin=min(ip1,ip2)
             ipmax=max(ip1,ip2)

             if(ipmin<=npcusp .or. ipmax<=npcusp)then

                if(ipmin<=npcusp)then
                   ipcusp=ipmin
                else
                   ipcusp=ipmax
                endif

                do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
                   iside=ptosi1(isto)
                   ipaa=min(lside(1,iside),lside(2,iside)) 
                   ipbb=max(lside(1,iside),lside(2,iside)) 
                   if(ipmin==ipaa .and. ipmax==ipbb)then
                      ireach=0_ip
                      call chkboun2(ip1,ip2,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                           nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                      goto 20 
                   endif
                enddo
             endif

             cycle
          endif
          !
          !     Side (ip1,ip3)
          !
          p1(1)=coor(1,ip3)-pproj(1)
          p1(2)=coor(2,ip3)-pproj(2)
          p1(3)=coor(3,ip3)-pproj(3)
          p2(1)=coor(1,ip1)-pproj(1)
          p2(2)=coor(2,ip1)-pproj(2)
          p2(3)=coor(3,ip1)-pproj(3)

          call orient3D(p1,p2,rface,d2,ndim)
          d2=d2/dtot
          ! 
          !     Compute normal to plane  
          !
          rx=coor(1,ip3)-coor(1,ipb)
          ry=coor(2,ip3)-coor(2,ipb)
          rz=coor(3,ip3)-coor(3,ipb)

          csca=rnx*rx+rny*ry+rnz*rz
          !
          !     And the test
          !
          !if(d2<-epsil)then
          if(csca<-epsil)then
             ie=eltoel(ltab(1,iview1),ie)
             !
             !     Did we cross the border ?
             !
             if(ie==0)then
                ireach=0_ip
                call chkboun2(ip1,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,&
                     nfapo,nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                goto 20 
             endif
             !
             !     Did we cross the front ?
             !
             if(lfmark(ie)==1)then
                ireach=0_ip
                call chkboun2(ip1,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                     nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                goto 20 
             endif
             !
             !     Did we cross an inner side ?
             !
             ipmin=min(ip1,ip3)
             ipmax=max(ip1,ip3)

             if(ipmin<=npcusp .or. ipmax<=npcusp)then

                if(ipmin<=npcusp)then
                   ipcusp=ipmin
                else
                   ipcusp=ipmax
                endif

                do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
                   iside=ptosi1(isto)
                   ipaa=min(lside(1,iside),lside(2,iside)) 
                   ipbb=max(lside(1,iside),lside(2,iside)) 
                   if(ipmin==ipaa .and. ipmax==ipbb)then
                      ireach=0_ip
                      call chkboun2(ip1,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                           nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                      goto 20 
                   endif
                enddo
             endif
             cycle
          endif
          !
          !     Element found 
          !
          exit
       enddo
       !
       !     Are we already in the final element
       !
       p1(1)=coor(1,ip2)-pproj(1)
       p1(2)=coor(2,ip2)-pproj(2)
       p1(3)=coor(3,ip2)-pproj(3)
       p2(1)=coor(1,ip3)-pproj(1)
       p2(2)=coor(2,ip3)-pproj(2)
       p2(3)=coor(3,ip3)-pproj(3)

       call orient3D(p1,p2,rface,d3,ndim)
       d3=d3/dtot

       if(d3>epsil)then
          goto 20
       endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !     Go through the pipe
       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ieold=ie
       ie=eltoel(iview1,ie)

       do
          !
          !    Are we still on the surface
          !
          if(ie==0)then
             ireach=0_ip
             goto 20
          endif

          if(eltoel(1,ie)==ieold)then
             iview1=1_ip
          else if(eltoel(2,ie)==ieold)then
             iview1=2_ip
          else
             iview1=3_ip
          endif

          ip1=lface(iview1,ie)
          ip2=lface(ltab(1,iview1),ie)
          ip3=lface(ltab(2,iview1),ie)
          !
          !     Did we cross the front ? 
          !
          if(lfmark(ie)==1)then
             ireach=0_ip
             call chkboun2(ip2,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                  nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
             goto 20 
          endif
          !
          !     Did we cross an inner side ?
          !
          ipmin=min(ip2,ip3)
          ipmax=max(ip2,ip3)

          if(ipmin<=npcusp .or. ipmax<=npcusp)then

             if(ipmin<=npcusp)then
                ipcusp=ipmin
             else
                ipcusp=ipmax
             endif

             do isto=ptosi2(ipcusp),ptosi2(ipcusp+1)-1
                iside=ptosi1(isto)
                ipaa=min(lside(1,iside),lside(2,iside)) 
                ipbb=max(lside(1,iside),lside(2,iside)) 
                if(ipmin==ipaa .and. ipmax==ipbb)then
                   ireach=0_ip
                   call chkboun2(ip2,ip3,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                        nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                   goto 20 
                endif
             enddo
          endif
          !
          !     Are we going back ?
          !
          if(lelem(ie)==1)then
             iback=1_ip
             call rescueside(ie,rnofa,nnofa,nface,d1,d2,d3,pnew,ndim,coor,npoin,lface,ierr)
             if(ierr==1)return   
             goto 20  
          endif
          lelem(ie)=1_ip
          nestack=nestack+1
          lestack(nestack)=ie
          !
          !     Does ipnew belong to this element
          !
          if(ipnew==ip1)goto 20
          !
          !     Check if the point has already been marked, if it belongs to the front, and the distances
          !
          if(lmark(ip1)==0)then
             lmark(ip1)=1_ip
             npstack2=npstack2+1
             lpstack2(npstack2)=ip1
             call length(ip1,ipnew,rsize,coor,rlen1)
             !
             !     Is the point too close ?
             !
             if(rlen1<tollen2)then
                ireach=0_ip
                lmark(ip1)=3_ip
                npstack1=npstack1+1
                lpstack1(npstack1)=ip1
                goto 20
             endif
             !
             !     Is the point in the local box ?
             !
             if(rlen1<tollen)then
                if(lpofa(ip1)>0)then
                   call chkboun(ip1,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
                        nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
                endif
             endif
             !
             !     Has this point been reached already in the first part as too close
             !
          else if(lmark(ip1)==3)then
             goto 20
          endif
          !
          !     Find which side is crossed
          ! 
          rface(1)=rnofa(1,ie)
          rface(2)=rnofa(2,ie)
          rface(3)=rnofa(3,ie)

          rx=pnew(1)-coor(1,ip1)
          ry=pnew(2)-coor(2,ip1)
          rz=pnew(3)-coor(3,ip1)

          rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
          pproj(1)= pnew(1)-rscal*rface(1)
          pproj(2)= pnew(2)-rscal*rface(2)
          pproj(3)= pnew(3)-rscal*rface(3)
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
          !     Side (ip1,ip2)
          !
          p1(1)=coor(1,ip2)-pproj(1)
          p1(2)=coor(2,ip2)-pproj(2)
          p1(3)=coor(3,ip2)-pproj(3)
          p2(1)=coor(1,ip1)-pproj(1)
          p2(2)=coor(2,ip1)-pproj(2)
          p2(3)=coor(3,ip1)-pproj(3)

          call orient3D(p1,p2,rface,d2,ndim)
          d2=d2/dtot
          !
          !     Side (ip1,ip3)
          !
          p1(1)=coor(1,ip1)-pproj(1)
          p1(2)=coor(2,ip1)-pproj(2)
          p1(3)=coor(3,ip1)-pproj(3)
          p2(1)=coor(1,ip3)-pproj(1)
          p2(2)=coor(2,ip3)-pproj(2)
          p2(3)=coor(3,ip3)-pproj(3)

          call orient3D(p1,p2,rface,d3,ndim)
          d3=d3/dtot

          if(d2<-epsil .and. d3<-epsil)then
             !
             !     This is the last element
             !
             exit
          else
             ! 
             !     Which side between (ip1,ip2) and (ip1,ip3) is intersected by (ipb,pnew)? 
             !
             p1(1)=coor(1,ip1)-coor(1,ipb)
             p1(2)=coor(2,ip1)-coor(2,ipb)
             p1(3)=coor(3,ip1)-coor(3,ipb)
             p2(1)=pproj(1)-coor(1,ipb)
             p2(2)=pproj(2)-coor(2,ipb)
             p2(3)=pproj(3)-coor(3,ipb)

             call orient3D(p1,p2,rface,d2,ndim)
             d2=d2/dtot
             ! 
             !     Compute normal to plane  
             !
             !rnx= rface(2)*rtz-rface(3)*rty
             !rny=-rface(1)*rtz+rface(3)*rtx
             !rnz= rface(1)*rty-rface(2)*rtx
             rnx= rnopy*rtz-rnopz*rty
             rny=-rnopx*rtz+rnopz*rtx
             rnz= rnopx*rty-rnopy*rtx
             rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
             rnl=c10/rnl
             rnx=rnx*rnl
             rny=rny*rnl
             rnz=rnz*rnl

             rx=coor(1,ip1)-coor(1,ipb)
             ry=coor(2,ip1)-coor(2,ipb)
             rz=coor(3,ip1)-coor(3,ipb)

             csca=rx*rnx+ry*rny+rz*rnz
             !
             !     And the test
             ! 
             !if(d2>epsil)then
             if(csca>epsil)then
                ieold=ie
                ie=eltoel(ltab(1,iview1),ie)
             else 
                ieold=ie
                ie=eltoel(ltab(2,iview1),ie)
             endif
          endif

       enddo

20     continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !           End of the march through ipnew
       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       !
       !     DBG
       !
       !do ipstack=1,npstack1
       !   if(lmark(lpstack1(ipstack))<2)then
       !      write(*,*)'Error lpstack1'
       !      stop
       !   endif
       !enddo



       !
       !   Did we cross the front?
       ! 
       if(ireach==1)then 
          !
          !     If this point is not ipfirst, it belongs already to the mesh
          !     -->   do not test for close points or close sides
          !   
          if(ipnew/=ipfirst)goto 50
          !
          !     Clean up lmark
          ! 
          do ipstack=1,npstack2
             lmark(lpstack2(ipstack))=0_ip
          enddo
          lmark(ipa)=0_ip
          lmark(ipb)=0_ip
          !
          !     Clean up lelem
          !
          do iestack=1,nestack
             lelem(lestack(iestack))=0_ip
          enddo
          nestack=0_ip
          !
          !     We did not cross the front, check for close points 
          !
          call gtclos(ie,ipnew,lface,nface,nnofa,coor,ndim,npoin,rsize,ipclos,lmark,lelem,&
               eltoel,lpofa,lfmark,ipa,ipb)
          ! 
          !     Reset lmark
          !
          do ipstack=1,npstack2
             lmark(lpstack2(ipstack))=1_ip
          enddo
          do ipstack=1,npstack1
             lmark(lpstack1(ipstack))=2_ip
          enddo
          lmark(ipa)=2_ip
          lmark(ipb)=2_ip

          if(ipclos==0)then
             !
             !     There are no close points, check for close faces
             !
             ip1=lface(1,ie)
             ip2=lface(2,ie)
             ip3=lface(3,ie)

             rface(1)=rnofa(1,ie)
             rface(2)=rnofa(2,ie)
             rface(3)=rnofa(3,ie)

             rx=pnew(1)-coor(1,ip1)
             ry=pnew(2)-coor(2,ip1)
             rz=pnew(3)-coor(3,ip1)

             rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
             pproj(1)= pnew(1)-rscal*rface(1)
             pproj(2)= pnew(2)-rscal*rface(2)
             pproj(3)= pnew(3)-rscal*rface(3)
             !
             !     Side (ip2,ip3)
             !
             ineigh=eltoel(1,ie)
             ichk=0_ip
             if(ineigh==0)then
                ichk=1_ip
             else if(lfmark(ineigh)==1)then 
                ichk=1_ip
             endif
             if(ichk==1)then 
                p1(1)=coor(1,ip3)-coor(1,ip2)
                p1(2)=coor(2,ip3)-coor(2,ip2)
                p1(3)=coor(3,ip3)-coor(3,ip2)
                rnl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
                rnli=c10/rnl
                p1(1)=p1(1)*rnli
                p1(2)=p1(2)*rnli
                p1(3)=p1(3)*rnli
                p2(1)=pproj(1)-coor(1,ip2)
                p2(2)=pproj(2)-coor(2,ip2)
                p2(3)=pproj(3)-coor(3,ip2)
                csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
                cscal=csca*rnli
                if(cscal<c00)then
                   p3(1)=pproj(1)-coor(1,ip2)
                   p3(2)=pproj(2)-coor(2,ip2)
                   p3(3)=pproj(3)-coor(3,ip2)
                else if(cscal>c10)then
                   p3(1)=pproj(1)-coor(1,ip3)
                   p3(2)=pproj(2)-coor(2,ip3)
                   p3(3)=pproj(3)-coor(3,ip3)
                else
                   p3(1)=coor(1,ip2)+csca*p1(1)-pproj(1)
                   p3(2)=coor(2,ip2)+csca*p1(2)-pproj(2)
                   p3(3)=coor(3,ip2)+csca*p1(3)-pproj(3)
                endif

                rnl2=sqrt(p3(1)*p3(1)+p3(2)*p3(2)+p3(3)*p3(3))
                if(rnl2*rnli<tolside)then
                   !
                   !     Add both points to lpstack1 if necessary
                   !  
                   if(lmark(ip2)<2)then
                      npstack1=npstack1+1
                      lpstack1(npstack1)=ip2
                      lmark(ip2)=2_ip
                   endif
                   if(lmark(ip3)<2)then
                      npstack1=npstack1+1
                      lpstack1(npstack1)=ip3
                      lmark(ip3)=2_ip
                   endif
                   goto 40
                endif
             endif
             !
             !     Side (ip3,ip1)
             !
             ineigh=eltoel(2,ie)
             ichk=0_ip
             if(ineigh==0)then
                ichk=1_ip
             else if(lfmark(ineigh)==1)then 
                ichk=1_ip
             endif
             if(ichk==1)then 
                p1(1)=coor(1,ip1)-coor(1,ip3)
                p1(2)=coor(2,ip1)-coor(2,ip3)
                p1(3)=coor(3,ip1)-coor(3,ip3)
                rnl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
                rnli=c10/rnl
                p1(1)=p1(1)*rnli
                p1(2)=p1(2)*rnli
                p1(3)=p1(3)*rnli
                p2(1)=pproj(1)-coor(1,ip3)
                p2(2)=pproj(2)-coor(2,ip3)
                p2(3)=pproj(3)-coor(3,ip3)
                csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
                cscal=csca*rnli
                if(cscal<c00)then
                   p3(1)=pproj(1)-coor(1,ip3)
                   p3(2)=pproj(2)-coor(2,ip3)
                   p3(3)=pproj(3)-coor(3,ip3)
                else if(cscal>c10)then
                   p3(1)=pproj(1)-coor(1,ip1)
                   p3(2)=pproj(2)-coor(2,ip1)
                   p3(3)=pproj(3)-coor(3,ip1)
                else
                   p3(1)=coor(1,ip3)+csca*p1(1)-pproj(1)
                   p3(2)=coor(2,ip3)+csca*p1(2)-pproj(2)
                   p3(3)=coor(3,ip3)+csca*p1(3)-pproj(3)
                endif

                rnl2=sqrt(p3(1)*p3(1)+p3(2)*p3(2)+p3(3)*p3(3))
                if(rnl2*rnli<tolside)then
                   !
                   !     Add both points to lpstack1 if necessary
                   !  
                   if(lmark(ip1)<2)then
                      npstack1=npstack1+1
                      lpstack1(npstack1)=ip1
                      lmark(ip1)=2_ip
                   endif
                   if(lmark(ip3)<2)then
                      npstack1=npstack1+1
                      lpstack1(npstack1)=ip3
                      lmark(ip3)=2_ip
                   endif
                   goto 40
                endif
             endif
             !
             !     Side (ip1,ip2)
             !
             ineigh=eltoel(3,ie)
             ichk=0_ip
             if(ineigh==0)then
                ichk=1_ip
             else if(lfmark(ineigh)==1)then 
                ichk=1_ip
             endif
             if(ichk==1)then 
                p1(1)=coor(1,ip2)-coor(1,ip1)
                p1(2)=coor(2,ip2)-coor(2,ip1)
                p1(3)=coor(3,ip2)-coor(3,ip1)
                rnl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
                rnli=c10/rnl
                p1(1)=p1(1)*rnli
                p1(2)=p1(2)*rnli
                p1(3)=p1(3)*rnli
                p2(1)=pproj(1)-coor(1,ip1)
                p2(2)=pproj(2)-coor(2,ip1)
                p2(3)=pproj(3)-coor(3,ip1)
                csca=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
                cscal=csca*rnli
                if(cscal<c00)then
                   p3(1)=pproj(1)-coor(1,ip1)
                   p3(2)=pproj(2)-coor(2,ip1)
                   p3(3)=pproj(3)-coor(3,ip1)
                else if(cscal>c10)then
                   p3(1)=pproj(1)-coor(1,ip2)
                   p3(2)=pproj(2)-coor(2,ip2)
                   p3(3)=pproj(3)-coor(3,ip2)
                else
                   p3(1)=coor(1,ip1)+csca*p1(1)-pproj(1)
                   p3(2)=coor(2,ip1)+csca*p1(2)-pproj(2)
                   p3(3)=coor(3,ip1)+csca*p1(3)-pproj(3)
                endif

                rnl2=sqrt(p3(1)*p3(1)+p3(2)*p3(2)+p3(3)*p3(3))
                if(rnl2*rnli<tolside)then
                   !
                   !     Add both points to lpstack1 if necessary
                   !  
                   if(lmark(ip2)<2)then
                      npstack1=npstack1+1
                      lpstack1(npstack1)=ip2
                      lmark(ip2)=2_ip
                   endif
                   if(lmark(ip1)<2)then
                      npstack1=npstack1+1
                      lpstack1(npstack1)=ip1
                      lmark(ip1)=2_ip
                   endif
                   goto 40
                endif
             endif
             !
             !     As we have passed all the tests succefully, go home
             !
             goto 50 

          else
             !
             !     ipclos is too close to ipnew
             !     Add ipclos to lpstack1
             ! 
             if(lmark(ipclos)==0)then
                lmark(ipclos)=2_ip
                npstack2=npstack2+1
                lpstack2(npstack2)=ipclos         
                npstack1=npstack1+1
                lpstack1(npstack1)=ipclos         
             else if(lmark(ipclos)<2)then
                lmark(ipclos)=2_ip
                npstack1=npstack1+1
                lpstack1(npstack1)=ipclos         
             endif

          endif

       endif

40     continue
       !
       !     We could not reach ipnew 
       !     Compute best point
       !
       do ipstack=1,npstack1
          ipoin=lpstack1(ipstack)
          lfacn(1)=ipa
          lfacn(2)=ipb
          lfacn(3)=ipoin
          call gtfnr2(lfacn,1_ip,nnofa,ndim,coor,npoin,rface,1_ip,modul)
          if(modul==c00)then
             rqual(ipstack)=qbig
             cycle
          endif
          call chkgeo(lfacn,1_ip,nnofa,1_ip,Qgeo,rface,rnopo,npoin,ndim)
          if(Qgeo<tolgeo)then
             rqual(ipstack)=qbig 
             cycle     
          endif
          modul=modul*c05
          call quality(lfacn,1_ip,1_ip,nnofa,coor,ndim,npoin,q1,modul)
          rqual(ipstack)=q1 
          !write(*,*)ipoin,q1
       enddo
       !
       !     Order the points with respect to their quality
       !        
       do ipstack=1,npstack1
          do jpstack=ipstack+1,npstack1
             ipoin=lpstack1(ipstack)
             jpoin=lpstack1(jpstack)
             if(rqual(ipstack)>rqual(jpstack))then
                rswap=rqual(jpstack)
                lpstack1(jpstack)=ipoin
                rqual(jpstack)=rqual(ipstack)
                lpstack1(ipstack)=jpoin
                rqual(ipstack)=rswap  
             endif
          enddo
       enddo

       !
       !     Assign to ipnew the new point to reach if not tried before
       !
       do ipstack=1,npstack1

          ipnew=lpstack1(ipstack)
          do jpstack=1,npopt
             if(lpopt(jpstack)==ipnew)goto 60
          enddo

          npopt=npopt+1
          lpopt(npopt)=ipnew
          goto 70

60        continue

       enddo
       !
       !     There are no more optimal points
       !
       !
       !     Enlarge tollen, and try again with ipfirst
       !
       !write(*,*)'Enlarging local bbox'
       tollen=tollen*c20
       if(tollen>100*rsize(ipa))then
          write(*,*)'Error in cross, point not found'
          ierr=1
          return
       endif
       ipnew=ipfirst

70     continue


    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     End of the choice of the optimal point
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
50  continue 

    !
    !     Did we reach ipfirst
    !
    if(ipnew==ipfirst)then
       ipcros=0_ip
       ienew=ie
       !
       !     Did we really reach the element ?
       !
       if(iback==0)then 
          !
          !     Transcribe element containing the new point
          !
          ip1=lface(1,ie)
          ip2=lface(2,ie)
          ip3=lface(3,ie)

          rface(1)=rnofa(1,ie)
          rface(2)=rnofa(2,ie)
          rface(3)=rnofa(3,ie)

          rx=pnew(1)-coor(1,ip1)
          ry=pnew(2)-coor(2,ip1)
          rz=pnew(3)-coor(3,ip1)

          rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
          pproj(1)= pnew(1)-rscal*rface(1)
          pproj(2)= pnew(2)-rscal*rface(2)
          pproj(3)= pnew(3)-rscal*rface(3)
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

          p1(1)=coor(1,ip2)-pproj(1)
          p1(2)=coor(2,ip2)-pproj(2)
          p1(3)=coor(3,ip2)-pproj(3)
          p2(1)=coor(1,ip3)-pproj(1)
          p2(2)=coor(2,ip3)-pproj(2)
          p2(3)=coor(3,ip3)-pproj(3)

          call orient3D(p1,p2,rface,d1,ndim)
          d1=d1/dtot
          !
          !     Side (ip3,ip1)
          !
          p1(1)=coor(1,ip3)-pproj(1)
          p1(2)=coor(2,ip3)-pproj(2)
          p1(3)=coor(3,ip3)-pproj(3)
          p2(1)=coor(1,ip1)-pproj(1)
          p2(2)=coor(2,ip1)-pproj(2)
          p2(3)=coor(3,ip1)-pproj(3)

          call orient3D(p1,p2,rface,d2,ndim)
          d2=d2/dtot
          !
          !     Side (ip1,ip2)
          !
          p1(1)=coor(1,ip1)-pproj(1)
          p1(2)=coor(2,ip1)-pproj(2)
          p1(3)=coor(3,ip1)-pproj(3)
          p2(1)=coor(1,ip2)-pproj(1)
          p2(2)=coor(2,ip2)-pproj(2)
          p2(3)=coor(3,ip2)-pproj(3)

          call orient3D(p1,p2,rface,d3,ndim)
          d3=d3/dtot
       endif

    else
       ipcros=ipnew
       ienew=ifirst
    endif
    !
    !     Clean up lmark
    !
    do ipstack=1,npstack2
       lmark(lpstack2(ipstack))=0_ip
    enddo
    lmark(ipa)=0_ip
    lmark(ipb)=0_ip
    !
    !     Clean up lelem
    !
    do iestack=1,nestack
       lelem(lestack(iestack))=0_ip
    enddo

  end subroutine cross

  subroutine chkboun(ipini,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,nfahol,&
       tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)     :: ipini,npoin,nfapo,nfahol,nfront,ipnew
    integer(ip), intent(in)     :: lpofa(npoin),lfapo(3,nfapo+nfahol)
    integer(ip), intent(in)     :: mstack1,mstack2,ndim,lfront(3,nfront) 
    integer(ip), intent(inout)  :: lmark(npoin),npstack1,npstack2,lpstack1(mstack1),lpstack2(mstack2)
    real(rp), intent(in)        :: tollen,coor(ndim,npoin),pnew(3),rsize(npoin)
    integer(ip)                 :: ipos,nfface,i,j,iface,ipstack,ip1,ipoin,npstackold  
    real(rp)                    :: rl 
    !
    !     This sub adds close points belonging to the current front in lpstack1
    !     The search begins with ipini
    !
    if(lmark(ipini)<2)then
       ipstack=npstack1
       npstack1=npstack1+1
       lpstack1(npstack1)=ipini
       lmark(ipini)=2_ip
    endif

    do

       if(ipstack==npstack1)exit
       ipstack=ipstack+1
       !
       !     Retrieve point from stack
       !
       ipoin=lpstack1(ipstack)
       !
       !     Loop on faces of the front attached to ipoin
       ! 
       ipos=lpofa(ipoin)

       do
          !
          !     Get face number attached to point
          ! 
          nfface=lfapo(3,ipos)
          !
          !     If we have two faces, know them explicitely
          !
          if(nfface==2)then
             do i=1,2
                iface=lfapo(i,ipos)
                do j=1,2
                   ip1=lfront(j,iface) 
                   if(lmark(ip1)==0)then
                      lmark(ip1)=1_ip
                      npstack2=npstack2+1
                      lpstack2(npstack2)=ip1
                      !
                      !     Is ip1 in the area of interest?
                      !
                      call length(ipnew,ip1,rsize,coor,rl)
                      if(rl<tollen)then
                         lmark(ip1)=2_ip
                         npstack1=npstack1+1
                         lpstack1(npstack1)=ip1
                      endif
                   endif
                enddo
             enddo
             exit
             !
             !     If we have one face, know it explicitely
             !
          else if(nfface==1)then
             iface=lfapo(1,ipos)
             do j=1,2
                ip1=lfront(j,iface) 
                if(lmark(ip1)==0)then
                   lmark(ip1)=1_ip
                   npstack2=npstack2+1
                   lpstack2(npstack2)=ip1
                   call length(ipnew,ip1,rsize,coor,rl)
                   !
                   !     Is ip1 in the area of interest?
                   !
                   if(rl<tollen)then
                      lmark(ip1)=2_ip
                      npstack1=npstack1+1
                      lpstack1(npstack1)=ip1
                   endif
                endif
             enddo
             exit
          else 
             !
             !     Should iterate through face database
             ! 
             do i=1,2
                iface=lfapo(i,ipos)
                do j=1,2
                   ip1=lfront(j,iface) 
                   if(lmark(ip1)==0)then
                      lmark(ip1)=1_ip
                      npstack2=npstack2+1
                      lpstack2(npstack2)=ip1
                      call length(ipnew,ip1,rsize,coor,rl)
                      if(rl<tollen)then
                         lmark(ip1)=2_ip
                         npstack1=npstack1+1
                         lpstack1(npstack1)=ip1
                      endif
                   endif
                enddo
             enddo
             ipos=-nfface
          endif
       enddo
    enddo

  end subroutine chkboun

  subroutine chkboun2(ipa,ipb,lmark,npoin,lpstack1,mstack1,lpstack2,mstack2,lfapo,nfapo,&
       nfahol,tollen,coor,ndim,npstack1,npstack2,lpofa,lfront,nfront,pnew,ipnew,rsize)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)     :: ipa,ipb,npoin,nfapo,nfahol,nfront,ipnew
    integer(ip), intent(in)     :: lpofa(npoin),lfapo(3,nfapo+nfahol)
    integer(ip), intent(in)     :: mstack1,mstack2,ndim,lfront(3,nfront) 
    integer(ip), intent(inout)  :: lmark(npoin),npstack1,npstack2,lpstack1(mstack1),lpstack2(mstack2)
    real(rp), intent(in)        :: tollen,coor(ndim,npoin),pnew(3),rsize(npoin)
    integer(ip)                 :: ipos,nfface,i,j,iface,ipstack,ip1,ipoin,npstackold  
    real(rp)                    :: rl
    !
    !     Contrarily to chkboun, this sub adds points if the front is crossed,
    !     not relying on a distance argument for the first two points
    !
    !     The search is initialized by the crossed edge with end points (ipa,ipb)
    !
    !
    ipstack=npstack1

    if(lmark(ipa)==0)then
       npstack2=npstack2+1   
       lpstack2(npstack2)=ipa
       npstack1=npstack1+1
       lpstack1(npstack1)=ipa
       lmark(ipa)=2_ip 
    else if(lmark(ipa)==1)then
       npstack1=npstack1+1
       lpstack1(npstack1)=ipa
       lmark(ipa)=2_ip 
    endif

    if(lmark(ipb)==0)then
       npstack2=npstack2+1   
       lpstack2(npstack2)=ipb
       npstack1=npstack1+1
       lpstack1(npstack1)=ipb
       lmark(ipb)=2_ip 
    else if(lmark(ipb)==1)then
       npstack1=npstack1+1
       lpstack1(npstack1)=ipb
       lmark(ipb)=2_ip 
    endif

    do

       if(ipstack==npstack1)exit
       ipstack=ipstack+1
       !
       !     Retrieve point from stack
       !
       ipoin=lpstack1(ipstack)
       !
       !     Loop on faces of the front attached to ipoin
       ! 

       ipos=lpofa(ipoin)

       do 
          !
          !     Get face number attached to point
          ! 
          nfface=lfapo(3,ipos)
          if(nfface==2)then
             !
             !     If we have two faces, know them explicitely
             !     and go home
             !
             do i=1,2
                iface=lfapo(i,ipos)
                do j=1,2
                   ip1=lfront(j,iface) 
                   if(lmark(ip1)==0)then
                      lmark(ip1)=1_ip
                      npstack2=npstack2+1
                      lpstack2(npstack2)=ip1
                      call length(ipnew,ip1,rsize,coor,rl)
                      if(rl<tollen)then
                         lmark(ip1)=2_ip
                         npstack1=npstack1+1
                         lpstack1(npstack1)=ip1
                      endif
                   endif
                enddo
             enddo
             exit
          else if(nfface==1)then
             !
             !     If we have one face, know it explicitely
             !     and go home
             !
             iface=lfapo(1,ipos)
             do j=1,2
                ip1=lfront(j,iface) 
                if(lmark(ip1)==0)then
                   lmark(ip1)=1_ip
                   npstack2=npstack2+1
                   lpstack2(npstack2)=ip1
                   call length(ipnew,ip1,rsize,coor,rl)
                   if(rl<tollen)then
                      lmark(ip1)=2_ip
                      npstack1=npstack1+1
                      lpstack1(npstack1)=ip1
                   endif
                endif
             enddo
             exit
          else  
             !
             !     Should iterate through face database
             ! 
             do i=1,2
                iface=lfapo(i,ipos)
                do j=1,2
                   ip1=lfront(j,iface) 
                   if(lmark(ip1)==0)then
                      lmark(ip1)=1_ip
                      npstack2=npstack2+1
                      lpstack2(npstack2)=ip1
                      call length(ipnew,ip1,rsize,coor,rl)
                      if(rl<tollen)then
                         lmark(ip1)=2_ip
                         npstack1=npstack1+1
                         lpstack1(npstack1)=ip1
                      endif
                   endif
                enddo
             enddo
             ipos=-nfface
          endif
       enddo

    enddo

  end subroutine chkboun2

  subroutine corgeo(nface,nnofa,lface,npoin,eltoelold,lsurfold,isurf,nfold,lside,nnosi,&
       nside,rnopo,rnofaold,ndim,lpsur,lptype,npold,lfold,lsurf,lelemold,lsold,nsold,lstof,&
       lpsid)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip), intent(in)     :: nface,nnofa,npoin,nside,isurf,npold
    integer(ip), intent(in)     :: nfold,ndim,nnosi,nsold
    integer(ip), intent(in)     :: lface(nnofa,nface),lside(nnosi,nside) 
    integer(ip), intent(in)     :: lfold(nnofa,nfold),lsurf(nface) 
    integer(ip), intent(in)     :: lsold(nnosi,nsold),lstof(nsold),lpsid(npoin)
    integer(ip), intent(in)     :: eltoelold(nnofa,nfold),lsurfold(nfold)
    integer(ip), intent(in)     :: lptype(2,npoin) 
    real(rp), intent(in)        :: rnofaold(ndim,nfold)
    real(rp), intent(inout)     :: rnopo(ndim,npoin)
    integer(ip),intent(inout)   :: lpsur(npoin),lelemold(nfold)
    integer(ip)                 :: iface,j,iside,is,ipoin,ifacn,ipold,iview,ineigh
    integer(ip)                 :: lstack(500),nstack,istack,iturn,ncont,inofa,ip1,ip2 
    integer(ip)                 :: lfloc(500),nfloc,ifloc,ifirst,ifnew,isidold 
    integer(4)                  :: istat
    integer(ip)                 :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    real(rp)                    :: rnx,rny,rnz,rnl,c00,c10
    !
    !     This sub corrects the point normals for corners and ridges 
    !     inside each patch
    !     We assume that lface contains only faces of isurf
    !
    c00=0.0d+00
    c10=1.0d+00
    !
    !     Loop on sides
    !
    do iside=1,nside
       do is=1,nnosi
          ipoin=lside(is,iside)
          !
          !     Find normal for ridge points
          !
          if(lptype(1,ipoin)==ID_RIDGE)then

             isidold=lpsid(ipoin) 
             ip1=lsold(1,isidold)
             ip2=lsold(2,isidold)
             iface=lstof(isidold)
             !  
             !     We are looking for a face containing ip1 and ip2 
             !
             lstack(1)=iface
             lelemold(iface)=1_ip
             nstack=1_ip
             istack=0_ip

             if(lsurfold(iface)==isurf)then
                nfloc=1_ip
                lfloc(1)=iface
             else
                nfloc=0_ip
             endif

             do 
                if(istack==nstack)exit
                istack=istack+1_ip
                iface=lstack(istack)
                do inofa=1,nnofa
                   ineigh=eltoelold(inofa,iface)  
                   if(ineigh==0)cycle
                   if(lelemold(ineigh)==1)cycle
                   ncont=0_ip
                   if(lfold(1,ineigh)==ip1)then
                      ncont=ncont+1
                   else if(lfold(2,ineigh)==ip1)then
                      ncont=ncont+1
                   else if(lfold(3,ineigh)==ip1)then
                      ncont=ncont+1
                   endif
                   if(lfold(1,ineigh)==ip2)then
                      ncont=ncont+1
                   else if(lfold(2,ineigh)==ip2)then
                      ncont=ncont+1
                   else if(lfold(3,ineigh)==ip2)then
                      ncont=ncont+1
                   endif

                   if(ncont>=1)then
                      lelemold(ineigh)=1
                      nstack=nstack+1
                      lstack(nstack)=ineigh
                      !
                      !     We only want isurf
                      ! 
                      if(lsurfold(ineigh)==isurf)then
                         nfloc=nfloc+1
                         lfloc(nfloc)=ineigh
                      endif
                   endif

                enddo
             enddo
             !
             !     Clean up lelemold
             !
             do istack=1,nstack
                lelemold(lstack(istack))=0_ip
             enddo
             !
             !     Compute normal
             !
             rnx=c00   
             rny=c00   
             rnz=c00   
             do ifloc=1,nfloc
                iface=lfloc(ifloc)
                rnx=rnx+rnofaold(1,iface)
                rny=rny+rnofaold(2,iface)
                rnz=rnz+rnofaold(3,iface)
             enddo
             rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
             rnl=c10/rnl
             rnopo(1,ipoin)=rnx*rnl
             rnopo(2,ipoin)=rny*rnl
             rnopo(3,ipoin)=rnz*rnl
             lpsur(ipoin)=lfloc(1)   

          else if(lptype(1,ipoin)==ID_CORNER)then 
             !
             !     The corner point can NOT be deleted and are at the beginning of the point array 
             !     so they keep the same numerotation
             !
             iface=lpsur(ipoin)
             ifirst=iface
             nfloc=0_ip

             ipold=ipoin  
             iturn=1_ip

             do 

                if(lfold(1,iface)==ipold)then
                   iview=1_ip 
                else if(lfold(2,iface)==ipold)then
                   iview=2_ip 
                else if(lfold(3,iface)==ipold)then
                   iview=3_ip 
                else
                   write(*,*)'error corgeo 3'
                endif

                if(lsurfold(iface)==isurf)then
                   !
                   !     Did we consider this face already ?
                   !
                   do ifloc=1,nfloc
                      if(lfloc(ifloc)==iface)then
                         exit
                      endif
                   enddo
                   !
                   !     Do we have to add this face ?
                   !  
                   if(ifloc>nfloc)then
                      nfloc=nfloc+1
                      lfloc(nfloc)=iface
                   endif

                endif

                if(iturn==1)then
                   ifnew=eltoelold(ltab(1,iview),iface)
                   if(ifnew==ifirst)exit
                   if(ifnew==0)then
                      ifnew=eltoelold(ltab(2,iview),iface)
                      if(ifnew==0)exit
                      iturn=0_ip
                   endif
                else
                   ifnew=eltoelold(ltab(2,iview),iface)
                   if(ifnew==0)exit
                endif
                iface=ifnew

             enddo
             !
             !     Compute normal
             !
             rnx=c00   
             rny=c00   
             rnz=c00   
             do ifloc=1,nfloc
                iface=lfloc(ifloc)
                rnx=rnx+rnofaold(1,iface)
                rny=rny+rnofaold(2,iface)
                rnz=rnz+rnofaold(3,iface)
             enddo
             rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
             rnl=c10/rnl
             rnopo(1,ipoin)=rnx*rnl
             rnopo(2,ipoin)=rny*rnl
             rnopo(3,ipoin)=rnz*rnl
             lpsur(ipoin)=lfloc(1)   

          endif

       enddo
    enddo
    !
    !     DBG
    !
    !do iface=1,nfold
    !   if(lelemold(iface)/=0)then
    !      write(*,*)'Error in  corgeo, lelemold not clean'
    !      stop
    !   endif
    !enddo

  end subroutine corgeo

  subroutine resizep(npnew,coor,rnopo,lcart,lpsur,rsize,lptri,lmark,lpofa,lptype,lpsid)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: npnew
    integer(ip),pointer    :: lcart(:),lpsur(:),lptri(:),lmark(:)
    integer(ip),pointer    :: lptype(:,:),lpsid(:),lpofa(:)
    real(rp),pointer       :: rnopo(:,:),coor(:,:),rsize(:)

    call memrea(npnew,memor_msh,'COOR','resizep',coor)
    call memrea(npnew,memor_msh,'RNOPO','resizep',rnopo)
    call memrea(npnew,memor_msh,'LCART','resizep',lcart)
    call memrea(npnew,memor_msh,'LPSUR','resizep',lpsur)
    call memrea(npnew,memor_msh,'RSIZE','resizep',rsize)
    call memrea(npnew,memor_msh,'LPTRI','resizep',lptri)
    call memrea(npnew,memor_msh,'LMARK','resizep',lmark)
    call memrea(npnew,memor_msh,'LPOFA','resizep',lpofa)
    call memrea(npnew,memor_msh,'LPTYPE','resizep',lptype)
    call memrea(npnew,memor_msh,'LPSID','resizep',lpsid)

  end subroutine resizep

  subroutine rescueside(ie,rnofa,nnofa,nface,d1,d2,d3,pnew,ndim,coor,npoin,lface,ierr)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)     :: ie,nnofa,nface,ndim,npoin
    integer(ip),intent(in)     :: lface(nnofa,nface)
    integer(ip),intent(inout)  :: ierr
    real(rp),intent(in)        :: rnofa(nnofa,nface),pnew(ndim),coor(ndim,npoin)
    real(rp),intent(inout)     :: d1,d2,d3
    real(rp)                   :: rscal,pproj(ndim),c00,p1(ndim),p2(ndim),p3(ndim)
    real(rp)                   :: rnl,rx,ry,rz,rface(ndim),c10,dtot
    integer(ip)                :: ip1,ip2,ip3,icont

    c00=0.0d+00
    c10=1.0d+00
    !
    !     Transcribe element containing the new point
    !
    ip1=lface(1,ie)
    ip2=lface(2,ie)
    ip3=lface(3,ie)

    rface(1)=rnofa(1,ie)
    rface(2)=rnofa(2,ie)
    rface(3)=rnofa(3,ie)

    rx=pnew(1)-coor(1,ip1)
    ry=pnew(2)-coor(2,ip1)
    rz=pnew(3)-coor(3,ip1)

    rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
    pproj(1)= pnew(1)-rscal*rface(1)
    pproj(2)= pnew(2)-rscal*rface(2)
    pproj(3)= pnew(3)-rscal*rface(3)
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
    p1(1)=coor(1,ip2)-pproj(1)
    p1(2)=coor(2,ip2)-pproj(2)
    p1(3)=coor(3,ip2)-pproj(3)
    p2(1)=coor(1,ip3)-pproj(1)
    p2(2)=coor(2,ip3)-pproj(2)
    p2(3)=coor(3,ip3)-pproj(3)

    call orient3D(p1,p2,rface,d1,ndim)
    d1=d1/dtot
    !
    !     Side (ip3,ip1)
    !
    p1(1)=coor(1,ip3)-pproj(1)
    p1(2)=coor(2,ip3)-pproj(2)
    p1(3)=coor(3,ip3)-pproj(3)
    p2(1)=coor(1,ip1)-pproj(1)
    p2(2)=coor(2,ip1)-pproj(2)
    p2(3)=coor(3,ip1)-pproj(3)

    call orient3D(p1,p2,rface,d2,ndim)
    d2=d2/dtot
    !
    !     Side (ip1,ip2)
    !
    p1(1)=coor(1,ip1)-pproj(1)
    p1(2)=coor(2,ip1)-pproj(2)
    p1(3)=coor(3,ip1)-pproj(3)
    p2(1)=coor(1,ip2)-pproj(1)
    p2(2)=coor(2,ip2)-pproj(2)
    p2(3)=coor(3,ip2)-pproj(3)

    call orient3D(p1,p2,rface,d3,ndim)
    d3=d3/dtot
    !
    !     DBG 
    !
    icont=0_ip
    if(d1>c00)icont=icont+1
    if(d2>c00)icont=icont+1
    if(d3>c00)icont=icont+1

    if(icont/=2)then
       write(*,*)'Strange configuration in rescueside'
       ierr=1 
       return
    endif

    !
    !     Project on the correct side
    ! 
    if(d1<c00)then
       p3(1)=coor(1,ip2)
       p3(2)=coor(2,ip2)
       p3(3)=coor(3,ip2)
       p1(1)=coor(1,ip3)-p3(1)
       p1(2)=coor(2,ip3)-p3(2)
       p1(3)=coor(3,ip3)-p3(3)
    else if(d2<c00)then 
       p3(1)=coor(1,ip3)
       p3(2)=coor(2,ip3)
       p3(3)=coor(3,ip3)
       p1(1)=coor(1,ip1)-p3(1)
       p1(2)=coor(2,ip1)-p3(2)
       p1(3)=coor(3,ip1)-p3(3)
    else  
       p3(1)=coor(1,ip1)
       p3(2)=coor(2,ip1)
       p3(3)=coor(3,ip1)
       p1(1)=coor(1,ip2)-p3(1)
       p1(2)=coor(2,ip2)-p3(2)
       p1(3)=coor(3,ip2)-p3(3)
    endif

    p2(1)=pnew(1)-p3(1)
    p2(2)=pnew(2)-p3(2)
    p2(3)=pnew(3)-p3(3)

    rnl=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
    rnl=c10/rnl
    p1(1)=p1(1)*rnl
    p1(2)=p1(2)*rnl
    p1(3)=p1(3)*rnl

    rscal=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3) 

    p3(1)=p3(1)+rscal*p1(1)
    p3(2)=p3(2)+rscal*p1(2)
    p3(3)=p3(3)+rscal*p1(3)
    !
    !     This will be the shape functions
    !
    !
    !     Side (ip2,ip3)
    !
    p1(1)=coor(1,ip2)-p3(1)
    p1(2)=coor(2,ip2)-p3(2)
    p1(3)=coor(3,ip2)-p3(3)
    p2(1)=coor(1,ip3)-p3(1)
    p2(2)=coor(2,ip3)-p3(2)
    p2(3)=coor(3,ip3)-p3(3)

    call orient3D(p1,p2,rface,d1,ndim)
    d1=d1/dtot
    !
    !     Side (ip3,ip1)
    !
    p1(1)=coor(1,ip3)-p3(1)
    p1(2)=coor(2,ip3)-p3(2)
    p1(3)=coor(3,ip3)-p3(3)
    p2(1)=coor(1,ip1)-p3(1)
    p2(2)=coor(2,ip1)-p3(2)
    p2(3)=coor(3,ip1)-p3(3)

    call orient3D(p1,p2,rface,d2,ndim)
    d2=d2/dtot
    !
    !     Side (ip1,ip2)
    !
    p1(1)=coor(1,ip1)-p3(1)
    p1(2)=coor(2,ip1)-p3(2)
    p1(3)=coor(3,ip1)-p3(3)
    p2(1)=coor(1,ip2)-p3(1)
    p2(2)=coor(2,ip2)-p3(2)
    p2(3)=coor(3,ip2)-p3(3)

    call orient3D(p1,p2,rface,d3,ndim)
    d3=d3/dtot

  end subroutine rescueside

  subroutine updatepipe(ibegin,iend,lpipe,npipe,ipbegin,ipnew,rnofa,nnofa,nface,coor,ndim,&
       npoin,lfmark,pnew,eltoel,lface,ierr,rnx,rny,rnz,ipa,mpipe)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)     :: ibegin,iend,ipbegin,nface,nnofa,ipnew
    integer(ip),intent(in)     :: ndim,npoin,ipa,mpipe 
    integer(ip),intent(in)     :: lfmark(nface),eltoel(nnofa,nface),lface(nnofa,nface)
    integer(ip),intent(inout)  :: npipe,lpipe(mpipe),ierr
    real(rp),intent(in)        :: rnofa(nnofa,nface),coor(ndim,npoin),pnew(ndim)
    real(rp),intent(in)        :: rnx,rny,rnz
    integer(ip)                :: ieold,ie,lploc(100),nploc,iview1,ip1,ip2,ip3,iploc 
    integer(ip)                :: ielemend,ipos,npipold 
    real(rp)                   :: rface(ndim),pproj(ndim),rscal,p1(ndim),p2(ndim)
    real(rp)                   :: rx,ry,rz,dtot,d1,d2,epsil,csca
    integer(ip)                :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This subroutine updates the pipe from the position ibegin-1 to the position iend
    !     deleting all elements between ibegin and iend-1 included
    !
    epsil=1.0d-09
    !
    !     Copy the pipe from lpipe to lploc
    !
    do ipos=1,ibegin-1
       lploc(ipos)=lpipe(ipos)
    enddo
    nploc=ibegin-1
   
    npipold=npipe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Go through the pipe
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    !     Does iend make sense?
    !
    if(iend==npipe+1)then
       ielemend=-1
    else
       ielemend=lpipe(iend)
    endif      
    !
    !     Get element to begin the update
    !
    ieold=lpipe(ibegin-1)
    if(lface(1,ieold)==ipbegin)then
       iview1=1_ip 
    else if(lface(2,ieold)==ipbegin)then
       iview1=2_ip 
    else
       iview1=3_ip 
    endif

    ie=eltoel(iview1,ieold)

    do

       if(eltoel(1,ie)==ieold)then
          iview1=1_ip
       else if(eltoel(2,ie)==ieold)then
          iview1=2_ip
       else
          iview1=3_ip
       endif

       ip1=lface(iview1,ie)
       ip2=lface(ltab(1,iview1),ie)
       ip3=lface(ltab(2,iview1),ie)
       !
       !     Did we cross the front
       !
       if(lfmark(ie)==1)then
          write(*,*)'Element already marked in update'
          ierr=1_ip
          return
       endif
       !
       !     Remember for later
       !
       nploc=nploc+1
       lploc(nploc)=ie
       !
       !     Does ipnew belong to this element ?
       !
       if(ipnew==ip1)exit
       !
       !     Did we reach the last element ?
       ! 
       if(ie==ielemend)exit
       !
       !     Find which side is crossed
       ! 
       rface(1)=rnofa(1,ie)
       rface(2)=rnofa(2,ie)
       rface(3)=rnofa(3,ie)
       !
       !     Which side between (ip1,ip2) and (ip1,ip3) is being crossed ?
       !
       rx=coor(1,ip1)-coor(1,ipa)
       ry=coor(2,ip1)-coor(2,ipa)
       rz=coor(3,ip1)-coor(3,ipa)

       csca=rx*rnx+ry*rny+rz*rnz

       if(csca>epsil)then
          ieold=ie
          ie=eltoel(ltab(1,iview1),ie)
       else 
          ieold=ie
          ie=eltoel(ltab(2,iview1),ie)
       endif

    enddo
    !
    !     We must have only one element
    !     NOOOOO  we may swap and still have both elements in the pipe
    !
    !if(nploc/=1 .and. nploc/=2)then
    !   write(*,*)'Error in updatepipe, nploc=',nploc
    !   ierr=1_ip
    !   return
    !   !stop
    !endif
    
    !
    !     Now copy the last element
    !
    do ipos=iend+1,npipe
       nploc=nploc+1
       lploc(nploc)=lpipe(ipos)
    enddo
    !
    !     Finally transfer to lpipe
    !
    npipe=nploc
    do ipos=1,npipe
       lpipe(ipos)=lploc(ipos)
    enddo   
  
    if(npipe>npipold)then
       write(*,*)'Error in recover, npipe>npipold:',npipe,npipold    
       ierr=1
    endif

  end subroutine updatepipe

  subroutine horder(nnofa,npoin,lfold,nfold,coor,coorold,npold,rnofaold,ndim,&
       lpsur,eltoelold,lstotr2old,lsurfold,nsurf)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnofa,npold,ndim,nsurf
    integer(ip),intent(in)     :: npoin,nfold
    integer(ip),intent(in)     :: lfold(nnofa,nfold),lpsur(npoin),lstotr2old(nsurf+1)
    integer(ip),intent(in)     :: lsurfold(nfold)
    integer(ip),pointer        :: eltoelold(:,:)
    integer(ip),pointer        :: ptoel1(:),ptoel2(:),lrenu(:)
    integer(ip),pointer        :: lftoed(:,:),ledge(:,:),lptosu2(:),lptosu1(:),lfloc(:,:)
    real(rp),intent(inout)     :: coor(ndim,npoin)
    real(rp),intent(in)        :: coorold(ndim,npold)
    real(rp),intent(out)       :: rnofaold(ndim,nfold)
    real(rp),pointer           :: redg(:,:),rnopoold(:,:)
    real(rp)                   :: s(ndim),svec1(ndim),svec2(ndim) 
    real(rp)                   :: r1(ndim),r2(ndim),st1(ndim),st2(ndim) 
    real(rp)                   :: n1(ndim),n2(ndim),dtot,d1,d2,d3 
    real(rp)                   :: snl,snli,stl1,stl2,c05,c125,p3(ndim),c10,c20,c40
    real(rp)                   :: p1(ndim),p2(ndim),xnew(ndim),rface(ndim),pproj(ndim)
    real(rp)                   :: rsvec1,rsvec2,rsvec1i,rsvec2i
    integer(ip)                :: iedg1,iedg2,iedg3,ip1,ip2,ip3,ifold,iedge,ipoin,nedge
    integer(ip)                :: nfloc,isurf,nsto,iplace,iface,isto,ifloc
    integer(4)                 :: istat

    c05=0.5d+00
    c125=0.125d+00
    c10=1.0d+00
    c20=2.0d+00
    c40=4.0d+00
    !
    !     Allocate rnopold
    !
    allocate(rnopoold(ndim,npold),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOPOOLD','horder',rnopoold) 
    allocate(lfloc(nnofa,nfold),stat=istat)
    call memchk(zero,istat,memor_msh,'LFLOC','horder',lfloc) 
    allocate(lptosu2(nsurf+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTOSU2','horder',lptosu2) 
    allocate(lrenu(nfold),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','horder',lrenu) 
    allocate(redg(ndim,nfold),stat=istat)
    call memchk(zero,istat,memor_msh,'REDG','ptoelm',redg)
    !
    !     Get the points associated to each surface
    !
    do ipoin=1,npoin
       ifold=lpsur(ipoin)
       isurf=lsurfold(ifold)+1_ip
       lptosu2(isurf)=lptosu2(isurf)+1_ip
    enddo
    !
    !    Accumulate
    ! 
    lptosu2(1)=1_ip
    do isurf=2,nsurf+1
       lptosu2(isurf)=lptosu2(isurf)+lptosu2(isurf-1)
    enddo
    nsto=lptosu2(nsurf+1)
    allocate(lptosu1(nsto),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTOSU2','horder',lptosu1) 
    !
    !     Get points
    !
    do ipoin=1,npoin
       ifold=lpsur(ipoin)
       isurf=lsurfold(ifold)
       iplace=lptosu2(isurf)
       lptosu1(iplace)=ipoin
       lptosu2(isurf)=iplace+1_ip
    enddo
    !
    !     Reshuffle
    !
    do isurf=nsurf+1,2,-1
       lptosu2(isurf)=lptosu2(isurf-1)
    enddo
    lptosu2(1)=1_ip
    !
    !     Loop on surfaces
    !
    do isurf=1,nsurf
       !
       !     Copy global to local faces
       !
       nfloc=0_ip
       do iface=lstotr2old(isurf),lstotr2old(isurf+1)-1
          nfloc=nfloc+1_ip
          lfloc(1,nfloc)=lfold(1,iface) 
          lfloc(2,nfloc)=lfold(2,iface) 
          lfloc(3,nfloc)=lfold(3,iface)
          lrenu(iface)=nfloc
       enddo
       !
       !     Get the faces surrounding the points
       !
       call ptoelm(lfloc,nfloc,npold,nnofa,ptoel1,ptoel2)
       !
       !     Get the face normals
       !
       call gtfnrl(lfloc,nfloc,nnofa,ndim,coorold,npold,rnofaold)
       !
       !     Get the point normals
       !
       call gtpnrl(nfloc,rnopoold,npold,ndim,ptoel1,ptoel2,rnofaold)
       !
       !     Get the faces surrounding faces
       !
       call trtotr(lfloc,nnofa,nfloc,ptoel1,ptoel2,npold,eltoelold)
       !
       !     For high order recovery, get the face to edge pointer
       !
       call fatoed(nedge,nfloc,lfloc,nnofa,eltoelold,ledge,lftoed)
       !
       !     Allocate normal for quadratic recovery
       !
       call memrea(nedge,memor_msh,'REDG','horder',redg) 
       !
       !     Compute Hermite interpolation 
       !
       do iedge=1,nedge
          !
          !     Edge end points
          !
          ip1=ledge(1,iedge) 
          ip2=ledge(2,iedge) 

          s(1)=coorold(1,ip2)-coorold(1,ip1)
          s(2)=coorold(2,ip2)-coorold(2,ip1)
          s(3)=coorold(3,ip2)-coorold(3,ip1)
          snl=sqrt(s(1)*s(1)+s(2)*s(2)+s(3)*s(3))
          snli=c10/snl
          s(1)=s(1)*snli
          s(2)=s(2)*snli
          s(3)=s(3)*snli
          !
          !     Point normals
          !
          n1(1)=rnopoold(1,ip1)
          n1(2)=rnopoold(2,ip1)
          n1(3)=rnopoold(3,ip1)
          n2(1)=rnopoold(1,ip2)
          n2(2)=rnopoold(2,ip2)
          n2(3)=rnopoold(3,ip2)
          !
          !     Bitangent vectors
          !
          svec1(1)= s(2)*n1(3)-s(3)*n1(2)
          svec1(2)=-s(1)*n1(3)+s(3)*n1(1)
          svec1(3)= s(1)*n1(2)-s(2)*n1(1)
          rsvec1=sqrt(svec1(1)*svec1(1)+svec1(2)*svec1(2)+svec1(3)*svec1(3))
          rsvec1i=c10/rsvec1
          svec1(1)=rsvec1i*svec1(1)
          svec1(2)=rsvec1i*svec1(2)
          svec1(3)=rsvec1i*svec1(3)

          svec2(1)= s(2)*n2(3)-s(3)*n2(2)
          svec2(2)=-s(1)*n2(3)+s(3)*n2(1)
          svec2(3)= s(1)*n2(2)-s(2)*n2(1)
          rsvec2=sqrt(svec2(1)*svec2(1)+svec2(2)*svec2(2)+svec2(3)*svec2(3))
          rsvec2i=c10/rsvec2
          svec2(1)=rsvec2i*svec2(1)
          svec2(2)=rsvec2i*svec2(2)
          svec2(3)=rsvec2i*svec2(3)
          !
          !     Tangent vectors
          !
          st1(1)= n1(2)*svec1(3)-n1(3)*svec1(2)
          st1(2)=-n1(1)*svec1(3)+n1(3)*svec1(1)
          st1(3)= n1(1)*svec1(2)-n1(2)*svec1(1)
          stl1=sqrt(st1(1)*st1(1)+st1(2)*st1(2)+st1(3)*st1(3))
          stl1=c10/stl1
          st1(1)=st1(1)*stl1
          st1(2)=st1(2)*stl1
          st1(3)=st1(3)*stl1

          st2(1)= n2(2)*svec2(3)-n2(3)*svec2(2)
          st2(2)=-n2(1)*svec2(3)+n2(3)*svec2(1)
          st2(3)= n2(1)*svec2(2)-n2(2)*svec2(1)
          stl2=sqrt(st2(1)*st2(1)+st2(2)*st2(2)+st2(3)*st2(3))
          stl2=c10/stl2
          st2(1)=st2(1)*stl2
          st2(2)=st2(2)*stl2
          st2(3)=st2(3)*stl2

          r1(1)=snl*st1(1)
          r1(2)=snl*st1(2)
          r1(3)=snl*st1(3)

          r2(1)=snl*st2(1)
          r2(2)=snl*st2(2)
          r2(3)=snl*st2(3)

          redg(1,iedge)=c05*coorold(1,ip1)+c125*r1(1)+c05*coorold(1,ip2)-c125*r2(1)
          redg(2,iedge)=c05*coorold(2,ip1)+c125*r1(2)+c05*coorold(2,ip2)-c125*r2(2)
          redg(3,iedge)=c05*coorold(3,ip1)+c125*r1(3)+c05*coorold(3,ip2)-c125*r2(3)

          !write(*,*)iedge,redg(1,iedge),redg(2,iedge),redg(3,iedge)

       enddo
       !
       !     Now interpolate the new points
       !
       do isto=lptosu2(isurf),lptosu2(isurf+1)-1
          !
          !     Get point
          !
          ipoin=lptosu1(isto)
          !
          !     Compute shape functions
          ! 
          ifold=lpsur(ipoin)

          ip1=lfold(1,ifold)
          ip2=lfold(2,ifold)
          ip3=lfold(3,ifold)
          !
          !     The point is on the surface  -->  pproj=coor
          !
          pproj(1)=coor(1,ipoin)
          pproj(2)=coor(2,ipoin)
          pproj(3)=coor(3,ipoin)

          rface(1)=rnofaold(1,ifold)
          rface(2)=rnofaold(2,ifold)
          rface(3)=rnofaold(3,ifold)
          !
          !     Compute area to normalize
          !
          p1(1)=coorold(1,ip2)-coorold(1,ip1)
          p1(2)=coorold(2,ip2)-coorold(2,ip1)
          p1(3)=coorold(3,ip2)-coorold(3,ip1)
          p2(1)=coorold(1,ip3)-coorold(1,ip1)
          p2(2)=coorold(2,ip3)-coorold(2,ip1)
          p2(3)=coorold(3,ip3)-coorold(3,ip1)
          call orient3D(p1,p2,rface,dtot,ndim) 

          p1(1)=coorold(1,ip1)-pproj(1)
          p1(2)=coorold(2,ip1)-pproj(2)
          p1(3)=coorold(3,ip1)-pproj(3)
          p2(1)=coorold(1,ip2)-pproj(1)
          p2(2)=coorold(2,ip2)-pproj(2)
          p2(3)=coorold(3,ip2)-pproj(3)
          p3(1)=coorold(1,ip3)-pproj(1)
          p3(2)=coorold(2,ip3)-pproj(2)
          p3(3)=coorold(3,ip3)-pproj(3)

          call orient3D(p2,p3,rface,d1,ndim)
          call orient3D(p3,p1,rface,d2,ndim)
          call orient3D(p1,p2,rface,d3,ndim)

          d1=d1/dtot
          d2=d2/dtot
          d3=d3/dtot

          ifloc=lrenu(ifold)
          iedg1=lftoed(1,ifloc)
          iedg2=lftoed(2,ifloc)
          iedg3=lftoed(3,ifloc)
          !write(*,*)ledge(1,iedg1),ledge(2,iedg1)
          !write(*,*)ledge(1,iedg2),ledge(2,iedg2)
          !write(*,*)ledge(1,iedg3),ledge(2,iedg3)
          !
          !     Quadratic recovery
          !
          xnew(1)=d1*(c20*d1-c10)*coorold(1,ip1)+d2*(c20*d2-c10)*coorold(1,ip2)&
               +d3*(c20*d3-c10)*coorold(1,ip3)+c40*d1*d2*redg(1,iedg3)+c40*d2*d3*redg(1,iedg1)+c40*d3*d1*redg(1,iedg2)
          xnew(2)=d1*(c20*d1-c10)*coorold(2,ip1)+d2*(c20*d2-c10)*coorold(2,ip2)&
               +d3*(c20*d3-c10)*coorold(2,ip3)+c40*d1*d2*redg(2,iedg3)+c40*d2*d3*redg(2,iedg1)+c40*d3*d1*redg(2,iedg2)
          xnew(3)=d1*(c20*d1-c10)*coorold(3,ip1)+d2*(c20*d2-c10)*coorold(3,ip2)&
               +d3*(c20*d3-c10)*coorold(3,ip3)+c40*d1*d2*redg(3,iedg3)+c40*d2*d3*redg(3,iedg1)+c40*d3*d1*redg(3,iedg2)
          !
          !     Limit the recovery
          !

          !
          !     Apply recover
          !
          coor(1,ipoin)=xnew(1)
          coor(2,ipoin)=xnew(2)
          coor(3,ipoin)=xnew(3)
       enddo

    enddo

    call memchk(2_ip,istat,memor_msh,'LRENU','horder',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','horder',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTOSU1','horder',lptosu1)
    deallocate(lptosu1,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTOSU1','horder',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTOSU2','horder',lptosu2)
    deallocate(lptosu2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTOSU2','horder',0_ip)
    call memchk(2_ip,istat,memor_msh,'REDG','horder',redg)
    deallocate(redg,stat=istat)
    if(istat/=0) call memerr(2_ip,'REDG','horder',0_ip)
    call memchk(2_ip,istat,memor_msh,'RNOPOOLD','horder',rnopoold)
    deallocate(rnopoold,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOPOOLD','horder',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','horder',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','horder',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','horder',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','horder',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGE','horder',ledge)
    deallocate(ledge,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGE','horder',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFTOED','horder',lftoed)
    deallocate(lftoed,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFTOED','horder',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFLOC','horder',lfloc)
    deallocate(lfloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFLOC','horder',0_ip)

  end subroutine horder

  subroutine optglo(lface,nnofa,nface,lptype,npoin,coor,ndim,rnopo,eltoelold,npold,&
       nfold,coorold,rnofaold,lpsur,rsize,lfold,nsurf,lsurf,lsurfold,ptoel1,ptoel2,&
       eltoel,lblmsh,lstotr2,ptosi1,ptosi2,lside,nnosi,nside)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH,memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)          :: ndim,nnofa,npoin,nnosi,nside
    integer(ip),intent(in)          :: nfold,npold,nsurf
    integer(ip),intent(inout)       :: nface
    integer(ip),intent(in)          :: lptype(2,npoin),lside(nnosi,nside)
    real(rp),intent(inout)          :: coor(ndim,npoin),rnopo(ndim,npoin)    
    integer(ip),intent(inout)       :: lblmsh(nsurf),lstotr2(nsurf+1)
    integer(ip),intent(in)          :: eltoelold(nnofa,nfold),lfold(nnofa,nfold)
    integer(ip),intent(in)          :: lsurfold(nfold)
    integer(ip),intent(inout)       :: lpsur(npoin),lsurf(nface)
    real(rp),intent(in)             :: coorold(ndim,npold)
    real(rp),intent(in)             :: rnofaold(ndim,nfold),rsize(npoin)
    integer(ip)                     :: iter,nitermax,ipoin,nfloc,ipnt,ipo,isurf,iface
    integer(ip)                     :: nfnew
    integer(ip),pointer             :: lface(:,:),ptosi1(:),ptosi2(:)
    integer(ip),pointer             :: lfloc(:,:),lelemold(:),ptoel1(:),ptoel2(:)
    integer(ip),pointer             :: lptri(:),eltoel(:,:),lfnew(:,:),lstotr2new(:) 
    integer(4)                      :: istat
    !
    !      This subroutine optimizes the surface triangulation
    !
    allocate(lfnew(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFLNEW','optglo',lfnew) 
    allocate(lfloc(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFLOC','optglo',lfloc) 
    allocate(lelemold(nfold),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEMOLD','optglo',lelemold) 
    allocate(lptri(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTRI','optglo',lptri)
    allocate(lstotr2new(nsurf+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTOTR2NEW','optglo',lstotr2new)
    !
    !     Number of optimization iterations
    !
    nitermax=5_ip 
    !
    !     Initialize new face number for this round
    ! 
    nfnew=0_ip
    !
    !     Get the sides surrounding the points
    !
    call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)
    !
    !     Loop on surface patch
    !
    do isurf=1,nsurf

       write(*,*)'Optimizing surface:',isurf
       !
       !     Clean up lpoin
       !      
       do ipoin=1,npoin
          lptri(ipoin)=0_ip
       enddo
       !
       !     Copy global to local faces
       !
       nfloc=0_ip
       do iface=lstotr2(isurf)+lblmsh(isurf),lstotr2(isurf+1)-1
          nfloc=nfloc+1_ip
          !write(*,*)lface(1,iface),lface(2,iface),lface(3,iface)
          lfloc(1,nfloc)=lface(1,iface) 
          lfloc(2,nfloc)=lface(2,iface) 
          lfloc(3,nfloc)=lface(3,iface)
          lptri(lfloc(1,nfloc))=nfloc
          lptri(lfloc(2,nfloc))=nfloc
          lptri(lfloc(3,nfloc))=nfloc
       enddo
       !
       !     Get the faces surrounding the points
       !
       call ptoelm(lfloc,nfloc,npoin,nnofa,ptoel1,ptoel2)
       ! 
       !     Get the faces surrounding faces
       !
       call trtotr(lfloc,nnofa,nfloc,ptoel1,ptoel2,npoin,eltoel)
       !
       !     Swap the whole patch
       !
       call optsmo(lfloc,nnofa,nfloc,coor,npoin,ndim,rnopo,eltoel,nitermax,coorold,nfold,&
            npold,rnofaold,lfold,eltoelold,lptri,lpsur,lelemold,rsize,lptype,lsurfold,isurf,&
            ptosi1,ptosi2,lside,nnosi,nside)
       !
       !     Copy bl mesh
       !
       do iface=lstotr2(isurf),lstotr2(isurf)+lblmsh(isurf)-1
          nfnew=nfnew+1_ip
          lfnew(1,nfnew)=lface(1,iface) 
          lfnew(2,nfnew)=lface(2,iface) 
          lfnew(3,nfnew)=lface(3,iface)
          lsurf(nfnew)=isurf
       enddo
       ! 
       !     Copy local to global
       !
       do iface=1,nfloc
          nfnew=nfnew+1
          lfnew(1,nfnew)=lfloc(1,iface) 
          lfnew(2,nfnew)=lfloc(2,iface) 
          lfnew(3,nfnew)=lfloc(3,iface)
          lsurf(nfnew)=isurf
       enddo

    enddo
    !
    !     Copy lfnew in lface
    !
    nface=nfnew
    call memrea(nface,memor_msh,'LFACE','refglo',lface) 
    do iface=1,nfnew
       lface(1,iface)=lfnew(1,iface)
       lface(2,iface)=lfnew(2,iface)
       lface(3,iface)=lfnew(3,iface)
    enddo
    !
    !     Copy lstotr2new
    !
    do isurf=2,nsurf+1        
       lstotr2(isurf)=lstotr2new(isurf)
    enddo

    call memchk(2_ip,istat,memor_msh,'LSTOTR2NEW','optglo',lstotr2new)
    deallocate(lstotr2new,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTOTR2NEW','optglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTRI','optglo',lptri)
    deallocate(lptri,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTRI','optglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LELEMOLD','optglo',lelemold)
    deallocate(lelemold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEMOLD','optglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFLOC','optglo',lfloc)
    deallocate(lfloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFLOC','optglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFNEW','optglo',lfnew)
    deallocate(lfnew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFNEW','optglo',0_ip)

  end subroutine optglo

  subroutine swapglo(lface,nnofa,nface,npoin,coor,ndim,nsurf,lsurf,niter,&
       ptoel1,ptoel2,eltoel,rnofa,lside,nnosi,nside,ptosi1,ptosi2)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH,memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)          :: ndim,nnofa,npoin,niter,nnosi,nside
    integer(ip),intent(in)          :: nsurf,nface,lside(nnosi,nside)
    real(rp),intent(in)             :: coor(ndim,npoin)   
    real(rp),intent(inout)          :: rnofa(ndim,nface)   
    integer(ip),intent(inout)       :: lface(nnofa,nface)
    integer(ip),intent(inout)       :: lsurf(nface)
    integer(ip),pointer             :: ptoel1(:),ptoel2(:),eltoel(:,:)
    integer(ip),pointer             :: ptosi1(:),ptosi2(:)
    integer(ip)                     :: iter,ipoin,nfloc,nfloc2,ipnt,ipo,isurf,iface
    integer(ip),pointer             :: lfloc(:,:),lmark(:)
    real(rp),pointer                :: rnopo(:,:)
    integer(4)                      :: istat
    !
    !      This subroutine swaps the whole surface triangulation
    !
    !
    !      The point normal used is the one of the actual triangulation
    !
    allocate(lfloc(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFLOC','swapglo',lfloc) 
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','swapglo',lmark) 
    allocate(rnopo(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOPO','swapglo',rnopo) 
    !
    !     Get the sides surrounding the points
    !
    call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)
    !
    !     Loop on surface patch
    !
    do isurf=1,nsurf
       !
       !     Copy global to local faces
       !
       nfloc=0_ip
       ipnt=0_ip
       do iface=1,nface
          if(lsurf(iface)==isurf)then 
             nfloc=nfloc+1
             lfloc(1,nfloc)=lface(1,iface) 
             lfloc(2,nfloc)=lface(2,iface) 
             lfloc(3,nfloc)=lface(3,iface)
             lmark(lface(1,iface))=1_ip     
             lmark(lface(2,iface))=1_ip     
             lmark(lface(3,iface))=1_ip     
          else
             ipnt=ipnt+1
             lface(1,ipnt)=lface(1,iface)
             lface(2,ipnt)=lface(2,iface)
             lface(3,ipnt)=lface(3,iface)
             lsurf(ipnt)=lsurf(iface)
          endif
       enddo
       nfloc2=ipnt
       !
       !     Get the points of this surface
       !
       ipnt=0_ip
       do ipoin=1,npoin
          if(lmark(ipoin)==1)then
             ipnt=ipnt+1
             lmark(ipnt)=ipoin
          endif
       enddo
       !
       !     Get the faces surrounding the points
       !
       call ptoelm(lfloc,nfloc,npoin,nnofa,ptoel1,ptoel2)
       !
       !     Get the face normals
       !
       call gtfnrl(lfloc,nfloc,nnofa,ndim,coor,npoin,rnofa)
       !
       !     Get the point normals
       !
       call gtpnrl(nfloc,rnopo,npoin,ndim,ptoel1,ptoel2,rnofa)
       ! 
       !     Get the faces surrounding faces
       !
       call trtotr(lfloc,nnofa,nfloc,ptoel1,ptoel2,npoin,eltoel)
       !
       !     Swap the whole patch only by swapping
       !
       call swapsurf(lfloc,nnofa,nfloc,coor,npoin,ndim,rnopo,eltoel,niter,ptosi1,&
            ptosi2,lside,nnosi,nside)
       ! 
       !     Copy local to global   (The mesh size is the same)
       !
       do iface=1,nfloc
          nfloc2=nfloc2+1
          lface(1,nfloc2)=lfloc(1,iface) 
          lface(2,nfloc2)=lfloc(2,iface) 
          lface(3,nfloc2)=lfloc(3,iface)
          lsurf(nfloc2)=isurf
       enddo

    enddo

    call memchk(2_ip,istat,memor_msh,'RNOPO','swapglo',rnopo)
    deallocate(rnopo,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOPO','swapglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','swapglo',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','swapglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFLOC','swapglo',lfloc)
    deallocate(lfloc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFLOC','swapglo',0_ip)

  end subroutine swapglo

  subroutine moveopt(ipnew,coor,ndim,npoin,lface,nnofa,nface,rnopo,coorold,nfold,npold,&
       rnofaold,lfold,eltoelold,eltoel,lptri,lpsur,lelemold,rsize,lptype,lsurfold,isurf)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)          :: nface,ndim,npoin,nnofa,ipnew,isurf
    integer(ip),intent(in)          :: npold,nfold
    integer(ip),intent(in)          :: lface(nnofa,nface),lptri(npoin)
    integer(ip),intent(in)          :: lfold(nnofa,nfold),eltoel(nnofa,nface)
    integer(ip),intent(in)          :: eltoelold(nnofa,nfold),lptype(2,npoin),lsurfold(nfold)
    integer(ip),intent(inout)       :: lpsur(npoin),lelemold(nfold)
    real(rp),intent(inout)          :: coor(ndim,npoin),rnopo(ndim,npoin)
    real(rp),intent(in)             :: coorold(ndim,npold)
    real(rp),intent(in)             :: rnofaold(ndim,nfold),rsize(npoin)
    integer(ip)                     :: ip1,ip2,ip3,iefirst,ienext,ipsur,ierr
    integer(ip)                     :: ie,iter,maxiter,lball(100),nball,iface,lfacn(nnofa,100),ielem,ihost 
    real(rp)                        :: pold(ndim),rnopoold(ndim),c10,wei1,wei2,Qgeo,Qgeotol,rface(ndim),pnew(ndim) 
    real(rp)                        :: d1,d2,d3,rpoin(ndim),rfloc(ndim,100),modul(100),c00
    real(rp)                        :: nrmal1(ndim),nrmal2(ndim),optdir(3,4),lenloc,epsil,rl 
    real(rp)                        :: Qini,Qiniloc,pinit(ndim),QTOL,coef1,coef2,TOLGAIN,c05
    real(rp)                        :: coefmin,q,lenabs,Qini0,ploc(ndim),q1,Qinigeo,q1geo
    real(rp)                        :: rpoint(ndim),rinit(ndim),rloc(ndim),rsiz
    integer(ip)                     :: chkimp,idir,idiropt,ldir(4),check,ihostt,ipsinit,isloc,iimpr,tab(4) 
    !
    !     This subroutine optimizes the point position by improving the quality of the
    !     patch around the point
    !
    c10=1.0d+00
    c00=0.0d+00
    c05=0.5d+00
    Qgeotol=0.0d+00
    epsil=1.0d-12
    QTOL=1.2d+00
    coef1=0.1d+00
    coef2=2.0d+00
    coefmin=coef1*coef1*coef1
    TOLGAIN=0.1d+00
    tab(1)=2_ip
    tab(2)=1_ip
    tab(3)=4_ip
    tab(4)=3_ip
    !
    !     The point must be a smooth point
    !
    if(lptype(1,ipnew)/=ID_SMOOTH)return
    !
    !     First loop, get the ball of the point. Check if it is worth the trouble 
    !
    ie=lptri(ipnew)
    !
    !     Fix local size
    !
    rsiz=rsize(ipnew)
    !
    !     DBG
    !
    if(lface(1,ie)==ipnew)then
    else if(lface(2,ie)==ipnew)then
    else if(lface(3,ie)==ipnew)then
    else
       write(*,*)'Error in moveopt, ipnew not found'
       stop
    endif

    lfacn(1,1)=lface(1,ie)
    lfacn(2,1)=lface(2,ie)
    lfacn(3,1)=lface(3,ie)
    iefirst=ie
    nball=1_ip
    lball(1)=ie

    do 
       if(lface(1,ie)==ipnew)then
          ienext=eltoel(2,ie)
       else if(lface(2,ie)==ipnew)then
          ienext=eltoel(3,ie)
       else
          ienext=eltoel(1,ie)
       endif

       if(ienext==iefirst)exit

       nball=nball+1
       lball(nball)=ienext
       lfacn(1,nball)=lface(1,ienext)
       lfacn(2,nball)=lface(2,ienext)
       lfacn(3,nball)=lface(3,ienext)
       ie=ienext

    enddo
    !
    !     Compute initial quality
    !
    Qini=c00

    do ie=1,nball
       call gtfnr2(lfacn,nball,nnofa,ndim,coor,npoin,rfloc(1,ie),ie,modul(ie))
       call chkgeo(lfacn,nball,nnofa,ie,Qgeo,rfloc(1,ie),rnopo,npoin,ndim)
       modul(ie)=modul(ie)*c05
       call quality(lfacn,nball,ie,nnofa,coor,ndim,npoin,q,modul(ie))
       if(q>Qini)Qini=q
    enddo
    !
    !     Is it worth the trouble?
    !
    if(Qini<QTOL)return
    !
    !     Remember the position, normal and lpsur of ipnew
    !
    pold(1)=coor(1,ipnew)
    pold(2)=coor(2,ipnew)
    pold(3)=coor(3,ipnew)
    rnopoold(1)=rnopo(1,ipnew)
    rnopoold(2)=rnopo(2,ipnew)
    rnopoold(3)=rnopo(3,ipnew)
    ipsur=lpsur(ipnew)
    !
    !     Get tangent plane
    !
    nrmal2(1)=c00
    nrmal2(2)=c00
    nrmal2(3)=c00

    if(abs(rnopoold(1))<epsil)then
       nrmal2(1)=c10
    else if(abs(rnopoold(2))<epsil)then
       nrmal2(2)=c10
    else 
       nrmal2(3)=c10
    endif

    nrmal1(1)= rnopoold(2)*nrmal2(3)-rnopoold(3)*nrmal2(2)
    nrmal1(2)=-rnopoold(1)*nrmal2(3)+rnopoold(3)*nrmal2(1)
    nrmal1(3)= rnopoold(1)*nrmal2(2)-rnopoold(2)*nrmal2(1)

    rl=sqrt(nrmal1(1)*nrmal1(1)+nrmal1(2)*nrmal1(2)+nrmal1(3)*nrmal1(3))
    rl=c10/rl
    nrmal1(1)=rl*nrmal1(1)
    nrmal1(2)=rl*nrmal1(2)
    nrmal1(3)=rl*nrmal1(3)

    nrmal2(1)= rnopoold(2)*nrmal1(3)-rnopoold(3)*nrmal1(2)
    nrmal2(2)=-rnopoold(1)*nrmal1(3)+rnopoold(3)*nrmal1(1)
    nrmal2(3)= rnopoold(1)*nrmal1(2)-rnopoold(2)*nrmal1(1)

    rl=sqrt(nrmal2(1)*nrmal2(1)+nrmal2(2)*nrmal2(2)+nrmal2(3)*nrmal2(3))
    rl=c10/rl
    nrmal2(1)=rl*nrmal2(1)
    nrmal2(2)=rl*nrmal2(2)
    nrmal2(3)=rl*nrmal2(3)
    !
    !     Fill optimal direction
    !
    optdir(1,1)= nrmal1(1)
    optdir(2,1)= nrmal1(2)
    optdir(3,1)= nrmal1(3)
    optdir(1,2)=-nrmal1(1)
    optdir(2,2)=-nrmal1(2)
    optdir(3,2)=-nrmal1(3)
    optdir(1,3)= nrmal2(1)
    optdir(2,3)= nrmal2(2)
    optdir(3,3)= nrmal2(3)
    optdir(1,4)=-nrmal2(1)
    optdir(2,4)=-nrmal2(2)
    optdir(3,4)=-nrmal2(3)

    lenloc=rsize(ipnew)
    !
    !     Initialize optimization
    !
    lenloc=lenloc*coef1
    lenabs=coef1
    Qini0=Qini
    pinit(1)=pold(1)
    pinit(2)=pold(2)
    pinit(3)=pold(3)
    rinit(1)=rnopoold(1)
    rinit(2)=rnopoold(2)
    rinit(3)=rnopoold(3)
    ipsinit=ipsur

    check=0_ip
    !
    !     Loop iteratively
    !   
    maxiter=10_ip

    do iter=1,maxiter

       !
       !     Reset optimization flag
       !
       chkimp=0_ip
       Qiniloc=Qini

       !
       !     Loop on direction
       !
       do idir=1,4
          !
          !     Compute new position
          !

          pnew(1)=pinit(1)+lenloc*optdir(1,idir)
          pnew(2)=pinit(2)+lenloc*optdir(2,idir)
          pnew(3)=pinit(3)+lenloc*optdir(3,idir)
          !
          !     Find the host face in the old mesh
          !
          call gthost(pnew,lpsur(ipnew),ihost,rnofaold,ndim,nfold,npold,lfold,nnofa,coorold,&
               d1,d2,d3,eltoelold,lelemold,lsurfold,isurf,ierr,rsiz)
          lpsur(ipnew)=ihost
          !
          !     Get the projected point
          !
          ip1=lfold(1,ihost)
          ip2=lfold(2,ihost)
          ip3=lfold(3,ihost)

          pnew(1)=d1*coorold(1,ip1)+d2*coorold(1,ip2)+d3*coorold(1,ip3)
          pnew(2)=d1*coorold(2,ip1)+d2*coorold(2,ip2)+d3*coorold(2,ip3)
          pnew(3)=d1*coorold(3,ip1)+d2*coorold(3,ip2)+d3*coorold(3,ip3)
          rpoin(1)=rnofaold(1,ihost)
          rpoin(2)=rnofaold(2,ihost)
          rpoin(3)=rnofaold(3,ihost)
          !
          !     Transfer data to ipnew
          !
          coor(1,ipnew)=pnew(1)
          coor(2,ipnew)=pnew(2)
          coor(3,ipnew)=pnew(3)
          rnopo(1,ipnew)=rpoin(1)
          rnopo(2,ipnew)=rpoin(2)
          rnopo(3,ipnew)=rpoin(3)
          !
          !     Verify correctness and quality 
          !
          q1=c00

          do ie=1,nball
             call gtfnr2(lfacn,nball,nnofa,ndim,coor,npoin,rfloc(1,ie),ie,modul(ie))
             if(modul(ie)==c00)exit
             call chkgeo(lfacn,nball,nnofa,ie,Qgeo,rfloc(1,ie),rnopo,npoin,ndim)
             if(Qgeo<Qgeotol)then 
                exit
             endif
             modul(ie)=modul(ie)*c05
             call quality(lfacn,nball,ie,nnofa,coor,ndim,npoin,q,modul(ie))
             if(q>Qini)exit
             if(q>q1)q1=q
          enddo

          if(ie>nball)then
             !
             !     Successful, remember the point
             !
             ploc(1)=pnew(1) 
             ploc(2)=pnew(2) 
             ploc(3)=pnew(3)
             rloc(1)=rpoin(1)
             rloc(2)=rpoin(2)
             rloc(3)=rpoin(3)
             isloc=ihost
             Qini=q1 
             chkimp=1_ip
             idiropt=idir 
          endif

       enddo
       !
       !     Did we optimize something?
       !
       if(chkimp==0)then
          !
          !     Optimization unsuccessfull for this round 
          !

          !
          !     Initialize ldir
          !
          ldir(1)=0_ip
          ldir(2)=0_ip
          ldir(3)=0_ip
          ldir(4)=0_ip
          !
          !      Reduce the steps
          !
          lenloc=lenloc*coef1
          lenabs=lenabs*coef1
          !
          !     Did we reduce too much ?
          !
          if(lenabs<coefmin)exit

       else
          !
          !     Optimization successfull
          !
          !
          !     Is it worth the trouble?
          !
          !if((Qiniloc-Qini)<TOLGAIN)exit
          !
          !     Mark the oposite direction
          !
          ldir(tab(idiropt))=0_ip
          !
          !     Assign the valid new point
          !
          pinit(1)=ploc(1)
          pinit(2)=ploc(2)
          pinit(3)=ploc(3)
          rinit(1)=rloc(1)
          rinit(2)=rloc(2)
          rinit(3)=rloc(3)
          ipsinit=isloc 
          check=1_ip
       endif
    enddo

    coor(1,ipnew)=pinit(1)
    coor(2,ipnew)=pinit(2)
    coor(3,ipnew)=pinit(3)
    rnopo(1,ipnew)=rinit(1)
    rnopo(2,ipnew)=rinit(2)
    rnopo(3,ipnew)=rinit(3)
    lpsur(ipnew)=ipsinit

  end subroutine moveopt

  subroutine renucorner(nface,ndim,npoin,nnofa,nnosi,nside,lface,coor,lside,&
       nboup,lptype,lline,lsurf,tolsca)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH,memor_msh
    implicit none
    integer(ip),intent(in)          :: nface,ndim,npoin,nnofa,nside,nnosi
    integer(ip),intent(inout)       :: nboup
    integer(ip),intent(inout)       :: lface(nnofa,nface),lside(nnosi,nside),lptype(2,npoin)
    integer(ip),intent(in)          :: lline(nside),lsurf(nface) 
    real(rp),intent(in)             :: tolsca 
    real(rp),intent(inout)          :: coor(ndim,npoin)
    integer(ip),pointer             :: lrenu(:),lptypet(:),ptosi1(:),ptosi2(:) 
    real(rp),pointer                :: coort(:,:)
    integer(ip)                     :: ipoin,isurf,iline,ipnew,ncont,inosi,iside
    integer(ip)                     :: inofa,iface,iline1,iline2,iside1,iside2,ip1,ip2
    real(rp)                        :: rtx1,rty1,rtz1,rtx2,rty2,rtz2,rtl,csca,c10
    integer(4)                      :: istat
    !
    !     This subroutine renumbers the corner points at the beginning of the array 
    !
    c10=1.0d+00 
    allocate(lrenu(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','renucorner',lrenu) 
    !
    !     Get the sides surrounding points
    !
    call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)
    !
    !     Count the boundary points
    !
    nboup=0_ip
    ncont=0_ip

    do ipoin=1,npoin
       iside=ptosi2(ipoin+1)-ptosi2(ipoin)

       if(iside==0)then
          !
          !     No side --> smooth point
          !  
          lptype(1,ipoin)=ID_SMOOTH

       else if(iside==1)then    
          !
          !     One side --> cusp point
          !  
          lptype(1,ipoin)=ID_CUSP
          nboup=nboup+1_ip
          ncont=ncont+1_ip
          lrenu(ipoin)=ncont

       else if(iside==2)then
          !
          !     Two sides --> Is it a corner point?
          !  
          iline1=lline(ptosi1(ptosi2(ipoin)))
          iline2=lline(ptosi1(ptosi2(ipoin)+1))
          if(iline1==iline2)then 
             !
             !     Check scalar product
             !
             iside1=ptosi1(ptosi2(ipoin)) 
             iside2=ptosi1(ptosi2(ipoin)+1)
             !
             !     Find ipoin
             ! 
             if(lside(1,iside1)==ipoin)then 
                ip1=ipoin
                ip2=lside(2,iside1)
             else
                ip1=ipoin
                ip2=lside(1,iside1)
             endif
             rtx1=coor(1,ip2)-coor(1,ip1) 
             rty1=coor(2,ip2)-coor(2,ip1) 
             rtz1=coor(3,ip2)-coor(3,ip1) 
             rtl=sqrt(rtx1*rtx1+rty1*rty1+rtz1*rtz1)
             rtl=c10/rtl
             rtx1=rtx1*rtl             
             rty1=rty1*rtl             
             rtz1=rtz1*rtl             


             if(lside(1,iside2)==ipoin)then 
                ip1=ipoin
                ip2=lside(2,iside2)
             else
                ip1=ipoin
                ip2=lside(1,iside2)
             endif
             rtx2=coor(1,ip1)-coor(1,ip2) 
             rty2=coor(2,ip1)-coor(2,ip2) 
             rtz2=coor(3,ip1)-coor(3,ip2) 
             rtl=sqrt(rtx2*rtx2+rty2*rty2+rtz2*rtz2)
             rtl=c10/rtl
             rtx2=rtx2*rtl             
             rty2=rty2*rtl             
             rtz2=rtz2*rtl             

             csca=rtx1*rtx2+rty1*rty2+rtz1*rtz1
             !
             !     And the test
             !
             if(csca>tolsca) then

                lptype(1,ipoin)=ID_RIDGE
                nboup=nboup+1_ip

             else 

                lptype(1,ipoin)=ID_CORNER
                nboup=nboup+1_ip
                ncont=ncont+1_ip
                lrenu(ipoin)=ncont

             endif

          else      
             lptype(1,ipoin)=ID_CORNER
             nboup=nboup+1_ip
             ncont=ncont+1_ip
             lrenu(ipoin)=ncont

          endif

       else

          lptype(1,ipoin)=ID_CORNER  
          nboup=nboup+1_ip
          ncont=ncont+1_ip
          lrenu(ipoin)=ncont

       endif
    enddo
    !
    !     Do we have something to do?
    !
    if(nboup==0)then
       call memchk(2_ip,istat,memor_msh,'LRENU','renucorner',lrenu)
       deallocate(lrenu,stat=istat)
       if(istat/=0) call memerr(2_ip,'LRENU','renucorner',0_ip)
       return 
    endif
    !
    !     Renumber the other points at the end of the array
    !
    do ipoin=1,npoin
       if(lrenu(ipoin)==0)then
          ncont=ncont+1_ip 
          lrenu(ipoin)=ncont
       endif
    enddo
    !
    !     Allocate temporal arrays
    !
    allocate(coort(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'COORT','renucorner',coort) 
    allocate(lptypet(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTYPET','renucorner',lptypet) 
    !
    !     Reorder points and types
    !     lptype(2:ipoin) will be updated outside
    !
    do ipoin=1,npoin
       ipnew=lrenu(ipoin)
       coort(1,ipnew)=coor(1,ipoin)
       coort(2,ipnew)=coor(2,ipoin)
       coort(3,ipnew)=coor(3,ipoin)
       lptypet(ipnew)=lptype(1,ipoin)
    enddo

    do ipoin=1,npoin
       coor(1,ipoin)=coort(1,ipoin)
       coor(2,ipoin)=coort(2,ipoin)
       coor(3,ipoin)=coort(3,ipoin)
       lptype(1,ipoin)=lptypet(ipoin)
    enddo

    call memchk(2_ip,istat,memor_msh,'LPTYPET','renucorner',lptypet)
    deallocate(lptypet,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTYPET','renucorner',0_ip)
    call memchk(2_ip,istat,memor_msh,'COORT','renucorner',coort)
    deallocate(coort,stat=istat)
    if(istat/=0) call memerr(2_ip,'COORT','renucorner',0_ip)
    !
    !     Renumber faces
    ! 
    do iface=1,nface
       do inofa=1,nnofa
          lface(inofa,iface)=lrenu(lface(inofa,iface))
       enddo
    enddo
    !
    !     Renumber sides
    ! 
    do iside=1,nside
       do inosi=1,nnosi
          lside(inosi,iside)=lrenu(lside(inosi,iside))
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LRENU','renucorner',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','renucorner',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI1','renucorner',ptosi1)
    deallocate(ptosi1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI1','renucorner',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI2','renucorner',ptosi2)
    deallocate(ptosi2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI2','renucorner',0_ip)

    write(*,*)'Dumping old mesh in outface.msh'        
    call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)

  end subroutine renucorner

  subroutine renuface(nface,nnofa,nsurf,lface,lsurf,lstotr2)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH,memor_msh
    implicit none
    integer(ip),intent(in)          :: nface,nnofa,nsurf
    integer(ip),intent(inout)       :: lface(nnofa,nface),lstotr2(nsurf+1)
    integer(ip),intent(inout)       :: lsurf(nface) 
    integer(ip),pointer             :: lrenu(:),lfacet(:,:) 
    integer(ip)                     :: isurf,iface,iplace,ifnew
    integer(4)                      :: istat
    !
    !     This subroutine renumbers the faces with respect to surfaces 
    !
    allocate(lrenu(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','renuface',lrenu) 
    allocate(lfacet(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFACET','renuface',lfacet) 
    !
    !     Clean up lstotr2
    !
    do isurf=1,nsurf+1 
       lstotr2(isurf)=0_ip
    enddo
    !
    !     Get pointer
    !
    do iface=1,nface  
       isurf=lsurf(iface)+1
       lstotr2(isurf)=lstotr2(isurf)+1
    enddo
    !
    !     Sum up
    !
    lstotr2(1)=1_ip
    do isurf=2,nsurf+1
       lstotr2(isurf)=lstotr2(isurf)+lstotr2(isurf-1)
    enddo
    !
    !     Get renumbering
    !
    do iface=1,nface
       isurf=lsurf(iface) 
       iplace=lstotr2(isurf)
       lrenu(iface)=iplace
       lstotr2(isurf)=iplace+1
    enddo
    !
    !    Clean up 
    !
    do isurf=nsurf+1,2,-1 
       lstotr2(isurf)=lstotr2(isurf-1)
    enddo

    lstotr2(1)=1_ip  
    !
    !     Now copy
    !  
    do iface=1,nface
       ifnew=lrenu(iface)
       lfacet(1,ifnew)=lface(1,iface)
       lfacet(2,ifnew)=lface(2,iface)
       lfacet(3,ifnew)=lface(3,iface)
    enddo

    do iface=1,nface
       lface(1,iface)=lfacet(1,iface)
       lface(2,iface)=lfacet(2,iface)
       lface(3,iface)=lfacet(3,iface)
    enddo

    do isurf=1,nsurf
       do iface=lstotr2(isurf),lstotr2(isurf+1)-1
          lsurf(iface)=isurf  
       enddo
    enddo


    call memchk(2_ip,istat,memor_msh,'LFACET','renuface',lfacet)
    deallocate(lfacet,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFACET','renuface',0_ip)
    call memchk(2_ip,istat,memor_msh,'LRENU','renuface',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','renuface',0_ip)

  end subroutine renuface

  subroutine renufacecol(nface,nnofa,nsurf,lface,lsurf,lstotr2,lelem)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH,memor_msh
    implicit none
    integer(ip),intent(in)          :: nface,nnofa,nsurf
    integer(ip),intent(inout)       :: lface(nnofa,nface),lstotr2(nsurf+1)
    integer(ip),intent(inout)       :: lsurf(nface),lelem(nface) 
    integer(ip),pointer             :: lrenu(:),lfacet(:,:),lelemt(:) 
    integer(ip)                     :: isurf,iface,iplace,ifnew
    integer(4)                      :: istat
    !
    !     This subroutine renumbers the faces with respect to surfaces 
    !
    allocate(lrenu(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','renufacecol',lrenu) 
    allocate(lfacet(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFACET','renufacecol',lfacet) 
    allocate(lelemt(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEMT','renufacecol',lelemt) 
    !
    !     Clean up lstotr2
    !
    do isurf=1,nsurf+1 
       lstotr2(isurf)=0_ip
    enddo
    !
    !     Get pointer
    !
    do iface=1,nface  
       isurf=lsurf(iface)+1
       lstotr2(isurf)=lstotr2(isurf)+1
    enddo
    !
    !     Sum up
    !
    lstotr2(1)=1_ip
    do isurf=2,nsurf+1
       lstotr2(isurf)=lstotr2(isurf)+lstotr2(isurf-1)
    enddo
    !
    !     Get renumbering
    !
    do iface=1,nface
       isurf=lsurf(iface) 
       iplace=lstotr2(isurf)
       lrenu(iface)=iplace
       lstotr2(isurf)=iplace+1
    enddo
    !
    !    Clean up 
    !
    do isurf=nsurf+1,2,-1 
       lstotr2(isurf)=lstotr2(isurf-1)
    enddo

    lstotr2(1)=1_ip  
    !
    !     Now copy
    !  
    do iface=1,nface
       ifnew=lrenu(iface)
       lfacet(1,ifnew)=lface(1,iface)
       lfacet(2,ifnew)=lface(2,iface)
       lfacet(3,ifnew)=lface(3,iface)
       lelemt(ifnew)=lelem(iface)
    enddo

    do iface=1,nface
       lface(1,iface)=lfacet(1,iface)
       lface(2,iface)=lfacet(2,iface)
       lface(3,iface)=lfacet(3,iface)
       lelem(iface)=lelemt(iface)
    enddo

    do isurf=1,nsurf
       do iface=lstotr2(isurf),lstotr2(isurf+1)-1
          lsurf(iface)=isurf  
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LELEMT','renufacecol',lelemt)
    deallocate(lelemt,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEMT','renufacecol',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFACET','renufacecol',lfacet)
    deallocate(lfacet,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFACET','renufacecol',0_ip)
    call memchk(2_ip,istat,memor_msh,'LRENU','renufacecol',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','renufacecol',0_ip)

  end subroutine renufacecol

  subroutine chkorient(nface,nnofa,npoin,lface)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none 
    integer(ip),intent(in)     :: npoin,nface,nnofa
    integer(ip),intent(in)     :: lface(nnofa,nface)
    integer(ip),pointer        :: lmark(:),lstack(:)
    integer(ip)                :: ip1,ip2,ipa,ipb,iface,istack,nstack,ineigh,ipnt,j,k
    integer(4)                 :: istat
    integer(ip)                :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    integer(ip),pointer        :: ptoel1(:),ptoel2(:),eltoel(:,:) 
    !
    !     This subroutine verifies that the faces have a consistent orientation
    !
    nullify(ptoel1,ptoel2,eltoel)
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface, nface , npoin, nnofa ,ptoel1,ptoel2 )
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Allocate lmark && lstack
    !
    allocate(lstack(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','chkorient',lstack) 
    allocate(lmark(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','chkorient',lmark) 

    ipnt=1_ip

    do 
       !
       !     Find a seed
       !
       do iface=ipnt,nface
          if(lmark(iface)==0)exit
       enddo

       if(iface==nface+1)exit

       ipnt=iface+1_ip
       !
       !     Loop on neighbors
       ! 
       lstack(1)=iface
       lmark(iface)=1_ip 
       nstack=1_ip
       istack=0_ip

       do

          if(istack==nstack)exit 
          istack=istack+1_ip

          iface=lstack(istack)
          do j=1,3
             ineigh=eltoel(j,iface)
             if(ineigh==0)cycle
             if(lmark(ineigh)==0)then
                lmark(ineigh)=1_ip   
                nstack=nstack+1_ip
                lstack(nstack)=ineigh
                ip1=lface(ltab(1,j),iface)
                ip2=lface(ltab(2,j),iface)

                do k=1,3
                   if(eltoel(k,ineigh)==iface)exit
                enddo
                if(k==4)then
                   write(*,*)'Error strange in chkorient'
                   stop
                endif

                ipa=lface(ltab(1,k),ineigh)
                ipb=lface(ltab(2,k),ineigh)

                if(ip1==ipb .and. ip2==ipa)then
                   continue
                else
                   write(*,*)'Error orientation in chkorient'
                   stop
                endif

             endif

          enddo

       enddo

    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','chkorient',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','chkorient',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','chkorient',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','chkorient',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','chkorient',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','chkorient',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','chkorient',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','chkorient',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','chkorient',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','chkorient',0_ip)


  end subroutine chkorient

  subroutine chkline(npoin,nside,nnosi,lptype,ptosi1,ptosi2,lline,lside)
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    use def_kintyp, only       :  ip,rp,lg
    implicit none 
    integer(ip),intent(in) ::  npoin,nside,nnosi
    integer(ip),intent(in) ::  lptype(2,npoin),ptosi1(*),ptosi2(*),lline(nside),lside(nnosi,nside)
    integer(ip)            :: ipoin,is,iside,iline,lstack(500),istack,nstack

    !
    !     This sub checks that the corners have different lines
    !
    do ipoin=1,npoin
       if(lptype(1,ipoin)==ID_CORNER)then
          nstack=0_ip
          do is=ptosi2(ipoin),ptosi2(ipoin+1)-1
             iside=ptosi1(is)
             iline=lline(iside)
             do istack=1,nstack
                if(lstack(istack)==iline)then
                   write(*,*)'Error in data at corner:',ipoin
                   stop    
                endif
             enddo
             nstack=nstack+1_ip
             lstack(nstack)=iline 
          enddo
       endif
    enddo

  end subroutine chkline

  subroutine distEdg(lside,nside,nnosi,iside,pnew,coor,ndim,npoin,rdist,d1,d2,ipoin)
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    use def_kintyp, only       :  ip,rp,lg
    implicit none 
    integer(ip),intent(in)    :: npoin,nside,nnosi,ndim,iside
    integer(ip),intent(in)    :: lside(nnosi,nside)
    real(rp),intent(in)       :: coor(ndim,npoin),pnew(ndim)
    real(rp),intent(inout)    :: rdist,d1,d2
    integer(ip),intent(inout) :: ipoin
    integer(ip)               :: ip1,ip2
    real(rp)                  :: rtx,rty,rtz,rnl,rx,ry,rz,csca,cscal  
    real(rp)                  :: rpx,rpy,rpz,rtol0,rtol1,c10   

    c10=1.0d+00
    rtol0=-1.0d-05
    rtol1=c10+1.0d-05
    !
    !     Initialize ipoin
    ! 
    ipoin=0_ip
    !
    !     Host points 
    !
    ip1=lside(1,iside) 
    ip2=lside(2,iside) 
    !
    !     Host tangent
    !
    rtx=coor(1,ip2)-coor(1,ip1)
    rty=coor(2,ip2)-coor(2,ip1)
    rtz=coor(3,ip2)-coor(3,ip1)
    rnl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
    rnl=c10/rnl 
    rtx=rnl*rtx  
    rty=rnl*rty  
    rtz=rnl*rtz
    !
    !     Project
    !
    rx=pnew(1)-coor(1,ip1) 
    ry=pnew(2)-coor(2,ip1) 
    rz=pnew(3)-coor(3,ip1) 
    csca=rtx*rx+rty*ry+rtz*rz
    cscal=csca*rnl
    d2=cscal
    d1=c10-d2

    if(cscal>rtol0)then

       if(cscal<rtol1)then

          rpx=coor(1,ip1)+csca*rtx-pnew(1)
          rpy=coor(2,ip1)+csca*rty-pnew(2)
          rpz=coor(3,ip1)+csca*rtz-pnew(3)

          rdist=sqrt(rpx*rpx+rpy*rpy+rpz*rpz)

       else 

          rx=pnew(1)-coor(1,ip2) 
          ry=pnew(2)-coor(2,ip2) 
          rz=pnew(3)-coor(3,ip2) 

          rdist=sqrt(rx*rx+ry*ry+rz*rz) 
          ipoin=ip2 

       endif

    else

       rdist=sqrt(rx*rx+ry*ry+rz*rz) 
       ipoin=ip1

    endif

  end subroutine distEdg

  subroutine optridg(npoin,lptype,coor,ndim,nside,lside,nnosi,coorold,&
       nfold,npold,rnofaold,lcell,ncell,lelemold,&
       sitosiold,lsold,nsold,lstof,ptosi2old,ptosi1old,rsuni,&
       rsize,lpsid,rnopo,lcart,lpsur,rtol,ptosi1,ptosi2,ierr,lline,llinold)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH,memor_msh
    use mod_memchk
    implicit none 
    integer(ip),intent(in)    :: npoin,ndim,nside,nnosi,ncell
    integer(ip),intent(in)    :: npold,nsold,nfold
    integer(ip),intent(in)    :: lptype(2,npoin),lside(2,nside),lline(nside)
    real(rp),intent(inout)    :: coor(ndim,npoin),rsize(npoin),rnopo(ndim,npoin)
    integer(ip),intent(in)    :: lstof(nsold),llinold(nsold)
    integer(ip),intent(in)    :: sitosiold(nnosi,nsold),lsold(nnosi,nsold),ptosi2old(npold+1)
    integer(ip),intent(in)    :: ptosi1old(*)
    integer(ip),intent(inout) :: lelemold(nfold),lcart(npoin),lpsur(npoin),lpsid(npoin),ierr
    type(cell),intent(in)     :: lcell(ncell) 
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold),rsuni,rtol 
    integer(ip),pointer       :: lpsmoo(:),lneigh(:,:)
    integer(ip),pointer       :: ptosi1(:),ptosi2(:)
    integer(ip)               :: ipoin,iter,nitermax,ipsmoo,npsmoo,nsid,isid1,isid2
    integer(ip)               :: ipa,ipb,ipc,ipn1,ipn2,is,iguess,iside,ihost,ihostf 
    integer(ip)               :: ip1o,ip2o,icart 
    real(rp)                  :: rsiz1,rsiz2,rstot,c10,pnew(ndim),d1,d2
    integer(4)                :: istat
    !
    !     This sub smoothes all the ridges
    !
    c10=1.0d+00
    !
    !     Do we have something to do?
    !
    if(nside==0)return
    !
    !     Allocate lpsmoo
    ! 
    allocate(lpsmoo(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPSMOO','chkorient',lpsmoo) 
    !
    !     Get the ridge points to keep them always 
    !
    npsmoo=0_ip
    do ipoin=1,npoin
       if(lptype(1,ipoin)==ID_RIDGE)then
          npsmoo=npsmoo+1
          lpsmoo(npsmoo)=ipoin
       endif
    enddo
    !
    !     Allocate lneigh
    ! 
    allocate(lneigh(2,npsmoo),stat=istat)
    call memchk(zero,istat,memor_msh,'LPSMOO','chkorient',lneigh) 
    !
    !     Fix the maximum iteration number
    !
    nitermax=5_ip
    !
    !    Get the points surrounding points
    ! 
    do ipsmoo=1,npsmoo     

       ipoin=lpsmoo(ipsmoo)

       nsid=ptosi2(ipoin+1)-ptosi2(ipoin)

       if(nsid/=2)then
          write(*,*)'Error in optridg'
          stop
       endif
       !
       !     Get the surrounding edges
       !
       isid1=ptosi1(ptosi2(ipoin))
       isid2=ptosi1(ptosi2(ipoin)+1)

       ipa=lside(1,isid1)
       ipb=lside(2,isid1)

       if(ipa==ipoin)then
          ipn1=ipb
       else
          ipn1=ipa
       endif

       ipa=lside(1,isid2)
       ipb=lside(2,isid2)

       if(ipa==ipoin)then
          ipn2=ipb
       else
          ipn2=ipa
       endif
       !
       !     Store the surrouding points
       !
       lneigh(1,ipsmoo)=ipn1         
       lneigh(2,ipsmoo)=ipn2         

    enddo
    !
    !     Loop on smoothing iterations
    !
    do iter=1,nitermax

       do ipsmoo=1,npsmoo     

          ipoin=lpsmoo(ipsmoo)
          ipn1=lneigh(1,ipsmoo)
          ipn2=lneigh(2,ipsmoo)
          !
          !     Compute optimal point coordinate
          !
          rsiz1=rsize(ipn1) 
          rsiz2=rsize(ipn2) 
          rstot=rsiz1+rsiz2
          rstot=c10/rstot
          rsiz1=rsiz1*rstot
          rsiz2=rsiz2*rstot

          pnew(1)=rsiz2*coor(1,ipn1)+rsiz1*coor(1,ipn2)
          pnew(2)=rsiz2*coor(2,ipn1)+rsiz1*coor(2,ipn2)
          pnew(3)=rsiz2*coor(3,ipn1)+rsiz1*coor(3,ipn2)

          if(lpsid(ipoin)/=0)then
             iguess=lpsid(ipoin)
          else if(lpsid(ipn1)/=0)then
             iguess=lpsid(ipn1) 
          else if(lpsid(ipn2)/=0)then
             iguess=lpsid(ipn2) 
          else
             !
             !     Strange case
             !  
             write(*,*)'Strange case in optridg'
             stop
             !
             !     Check scalar product
             ! 
             do is=ptosi2old(ipn1),ptosi2old(ipn1+1)-1
                iside=ptosi1old(is)
                ipa=lsold(1,iside) 
                ipb=lsold(2,iside)
                if(ipa==ipn1)then
                   ipc=ipb
                else
                   ipc=ipa
                endif
                if(ipc==ipn2)goto 10
             enddo

             write(*,*)'Error in optridg, ip3 not found'
             stop     

10           continue

             iguess=iside
          endif
          call gthostsid(pnew,nsold,nnosi,npold,iguess,ihost,ndim,lsold,sitosiold,&
               coorold,d1,d2,ptosi2old,ptosi1old,npoin,lelemold,ierr,llinold,llinold(iguess))
          if(ierr==1)then
             write(*,*)'Error in optridg'
             return
          endif
          !
          !     Get a face attached to this side
          !
          ihostf=lstof(ihost)
          !
          !     Get the projected point
          !
          ip1o=lsold(1,ihost)
          ip2o=lsold(2,ihost)
          !
          !     Should check geometry...
          !
          coor(1,ipoin)=d1*coorold(1,ip1o)+d2*coorold(1,ip2o)
          coor(2,ipoin)=d1*coorold(2,ip1o)+d2*coorold(2,ip2o)
          coor(3,ipoin)=d1*coorold(3,ip1o)+d2*coorold(3,ip2o)
          rnopo(1,ipoin)=rnofaold(1,ihostf)
          rnopo(2,ipoin)=rnofaold(2,ihostf)
          rnopo(3,ipoin)=rnofaold(3,ihostf)
          icart=lcart(ipn1)

          call gtelem(ipoin,coor,npoin,ndim,lcell,ncell,icart,rtol)
          lcart(ipoin)=icart
          call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipoin,ndim,coor,rsuni)
          lpsur(ipoin)=ihostf
          lpsid(ipoin)=ihost 

       enddo

    enddo

    call memchk(2_ip,istat,memor_msh,'LNEIGH','optridg',lneigh)
    deallocate(lneigh,stat=istat)
    if(istat/=0) call memerr(2_ip,'LNEIGH','optridg',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPSMOO','optridg',lpsmoo)
    deallocate(lpsmoo,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPSMOO','optridg',0_ip)

  end subroutine optridg

  subroutine presmoo(eltoel,lface,nface,nnofa,ndim,coor,npoin,rnofa,ptoel1,ptoel2,&
       nboup,lsurf,lptype,rnopo,nsurf,lside,nside,lline,nline)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH,memor_msh
    use mod_memchk
    implicit none 
    integer(ip),intent(in)    :: nface,nnofa,ndim,npoin,nboup,nsurf,nside,nline 
    integer(ip),intent(inout) :: lface(nnofa,nface)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),eltoel(:,:)
    integer(ip),intent(inout) :: lsurf(nface)
    integer(ip),intent(in)    :: lptype(2,npoin),lside(2,nside),lline(nside)
    real(rp),intent(out)      :: rnofa(nnofa,nface),rnopo(ndim,npoin)
    real(rp),intent(inout)    :: coor(ndim,npoin)
    integer(ip)               :: iface,inofa,ineigh,iridg,nridg,nridg0,jnofa
    integer(ip)               :: ipf1,ipf2,ipoin1,ipoin2,leloc(2),neloc,iface1,iface2 
    integer(ip)               :: ip1,ip2,ip3,ie,icont,ipoin,np1,np2,irid,ipa,ipb,ipc
    integer(ip)               :: iploc,jsurf,isurf,nstack,istack,isurf1,isurf2
    integer(ip)               :: isid,isi,maxlevel,ilevel,nitermax,iter,jridg,ncont 
    integer(ip),parameter     :: mpoin=100,mstack=1000 
    integer(ip)               :: niter,jface,jp1,jp2,jp3,kface,knofa,iview,ipt,jpoin
    integer(ip)               :: lp1(mpoin),lp2(mpoin),lstack(mstack),level(mstack)
    integer(ip)               :: itergeo,nitergeo,kneigh 
    integer(ip),pointer       :: lridg(:,:),lmark(:),lmarkf(:)
    real(rp)                  :: tolscal,cscal,v(ndim),h 
    real(rp)                  :: nrmal1(ndim),nrmal2(ndim),nrmal(ndim),d1(ndim),d2(ndim) 
    real(rp)                  :: x1new(ndim),x2new(ndim),x1s(ndim),x2s(ndim),vect(ndim)
    real(rp)                  :: dxs1(ndim),dxs2(ndim),dx(ndim),vdxs1(ndim),relax
    real(rp)                  :: coef1,coef2,coef3,rnl,c10,c00,csca1,csca2,csca3,tolnorm 
    real(rp)                  :: rnl1,rnl2,rweight,rdist1t,rdist2t,rdist1(mpoin),rdist2(mpoin) 
    real(rp)                  :: rx,ry,rz,rdist12,rx1,ry1,rz1,rx2,ry2,rz2,c05,rtl,lambda 
    real(rp)                  :: rnx1,rny1,rnz1,rnx2,rny2,rnz2,rnx,rny,rnz,rtx,rty,rtz 
    integer(4)                :: istat
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))
    !
    !     This subroutine smoothes out the initial surface on ridges not 
    !     recognized as ridges
    !
    !     The surface is smoothed in a volume conserving fashion as described in: 
    ! 
    !     'Volume conserving smoothing for piecewise linear curves, surfaces and triple lines',
    !      Kuprat et al, JCP 172, 99-118 (2001) 
    !
    !
    c00=0.0d+00
    c05=0.5d+00
    c10=1.0d+00
    tolnorm=1.0d-10
    tolscal=0.95d+00
    relax=0.8d+00
    nitermax=1000_ip
    nitergeo=1_ip
    !
    !     First swap the old face to improve geometry
    !
    !niter=5_ip
    !call swapglo(lface,nnofa,nface,npoin,coor,ndim,rnopo,nsurf,lsurf,niter,&
    !   ptoel1,ptoel2,eltoel)
    !write(*,*)'Dumping swapped old mesh in outface.msh'        
    !call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)
    !
    !     Get the faces surrounding the points
    !
    !call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get the faces surrounding faces
    !
    !call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)

    nridg=0_ip
    goto 100
    !
    !   Allocate lmark 
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','presmoo',lmark) 
    !
    !     First identify the riges not detected as ridges
    !
    nridg=0_ip
    do iface=1,nface
       isurf=lsurf(iface)
       do inofa=1,nnofa
          ineigh=eltoel(inofa,iface)
          if(ineigh<=0)cycle
          jsurf=lsurf(ineigh)
          if(isurf/=jsurf)cycle 
          cscal=rnofa(1,iface)*rnofa(1,ineigh)+rnofa(2,iface)*rnofa(2,ineigh)+&
               rnofa(3,iface)*rnofa(3,ineigh)
          if(cscal<tolscal)then
             nridg=nridg+1_ip
             !write(*,*)nridg,iface,inofa
          endif
          eltoel(inofa,iface)=-ineigh 

          if(eltoel(1,ineigh)==iface)then
             eltoel(1,ineigh)=-iface
          else if(eltoel(2,ineigh)==iface)then
             eltoel(2,ineigh)=-iface
          else
             eltoel(3,ineigh)=-iface
          endif

       enddo
    enddo
    !
    !   Allocate lridg 
    !
    allocate(lridg(4,nridg),stat=istat)
    call memchk(zero,istat,memor_msh,'LRIDG','presmoo',lridg) 
    !
    !     Store the ridges
    !
    nridg=0_ip
    do iface=1,nface
       isurf=lsurf(iface)
       do inofa=1,nnofa
          ineigh=-eltoel(inofa,iface)
          if(ineigh<=0)cycle
          jsurf=lsurf(ineigh)
          if(isurf/=jsurf)cycle 
          cscal=rnofa(1,iface)*rnofa(1,ineigh)+rnofa(2,iface)*rnofa(2,ineigh)+&
               rnofa(3,iface)*rnofa(3,ineigh)
          if(cscal<tolscal)then
             nridg=nridg+1_ip
             !write(*,*)nridg,iface,inofa
             !write(*,*)'nridg=',nridg
             lridg(1,nridg)=lface(ltab(1,inofa),iface)
             lridg(2,nridg)=lface(ltab(2,inofa),iface)
             lridg(3,nridg)=iface
             lridg(4,nridg)=ineigh
          endif
          eltoel(inofa,iface)=ineigh 

          if(eltoel(1,ineigh)==-iface)then
             eltoel(1,ineigh)=iface
          else if(eltoel(2,ineigh)==-iface)then
             eltoel(2,ineigh)=iface
          else
             eltoel(3,ineigh)=iface
          endif
       enddo
    enddo
    !
    !     DBG
    !
    do iface=1,nface
       do inofa=1,nnofa
          ineigh=eltoel(inofa,iface)
          if(ineigh<0)then
             write(*,*)'Error in presmoo: eltoel not clean'
             stop
          endif
       enddo
    enddo
    call outpresmoo(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,lline,&
         lridg,nridg,nline)
    !
    !   Allocate lmarkf 
    !
    allocate(lmarkf(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKF','presmoo',lmarkf) 
    !
    !     Fix the number of layers around each ridge
    !
    maxlevel=10_ip
    !
    !     Mark the current ridges
    !     
    do iridg=1,nridg

       iface1=lridg(3,iridg)
       iface2=lridg(4,iridg)

       if(eltoel(1,iface1)==iface2)then
          eltoel(1,iface1)=-iface2 
       else if(eltoel(2,iface1)==iface2)then
          eltoel(2,iface1)=-iface2 
       else   
          eltoel(3,iface1)=-iface2 
       endif


       if(eltoel(1,iface2)==iface1)then
          eltoel(1,iface2)=-iface1 
       else if(eltoel(2,iface2)==iface1)then
          eltoel(2,iface2)=-iface1 
       else   
          eltoel(3,iface2)=-iface1 
       endif

    enddo
    !
    !     Add the edges surrounding the ridges not already detected
    !     to also smooth the surrounding region arround the ridges
    !
    nridg0=nridg

    do iridg=1,nridg0

       iface1=lridg(3,iridg)
       iface2=lridg(4,iridg)

       lstack(1)=iface1
       level(1)=1_ip
       lstack(2)=iface2
       level(2)=1_ip
       lmarkf(iface1)=iridg
       lmarkf(iface2)=iridg
       nstack=2_ip
       istack=0_ip
       !
       !     Loop on stack
       !
       do  
          if(istack==nstack)exit
          istack=istack+1_ip 
          iface=lstack(istack)
          isurf=lsurf(iface)  
          ilevel=level(istack)        

          if(ilevel==maxlevel)exit
          !
          !     Loop on neighbors
          !
          do inofa=1,nnofa
             ineigh=eltoel(inofa,iface)
             !
             !     Do we have a neighbor or has this edge already been marked?
             !
             if(ineigh<=0)cycle
             !
             !     Get surface number
             ! 
             jsurf=lsurf(ineigh)         
             !
             !     Are we still on the same surface?
             ! 
             if(isurf/=jsurf)cycle 
             !
             !     Is ineigh already in the stack?
             !
             if(lmarkf(ineigh)/=iridg)then
                nstack=nstack+1_ip
                lstack(nstack)=ineigh
                level(nstack)=ilevel+1_ip         
                lmarkf(ineigh)=iridg
             endif
             !
             !     Now add this edge
             ! 
             nridg=nridg+1_ip 
             call memrea(nridg,memor_msh,'LRIDG','presmoo',lridg)
             lridg(1,nridg)=lface(ltab(1,inofa),iface)            
             lridg(2,nridg)=lface(ltab(2,inofa),iface)            
             lridg(3,nridg)=iface
             lridg(4,nridg)=ineigh

             eltoel(inofa,iface)=-ineigh

             if(eltoel(1,ineigh)==iface)then
                eltoel(1,ineigh)=-iface
             else if(eltoel(2,ineigh)==iface)then
                eltoel(2,ineigh)=-iface
             else
                eltoel(3,ineigh)=-iface
             endif
          enddo

       enddo

    enddo
    !
    !     Clean up eltoel
    !
    do iridg=1,nridg
       isurf1=lridg(3,iridg)
       eltoel(1,isurf1)=abs(eltoel(1,isurf1)) 
       eltoel(2,isurf1)=abs(eltoel(2,isurf1)) 
       eltoel(3,isurf1)=abs(eltoel(3,isurf1)) 
       isurf2=lridg(4,iridg)
       eltoel(1,isurf2)=abs(eltoel(1,isurf2)) 
       eltoel(2,isurf2)=abs(eltoel(2,isurf2)) 
       eltoel(3,isurf2)=abs(eltoel(3,isurf2)) 
    enddo
    !
    !     DBG
    !
    do iface=1,nface
       do inofa=1,nnofa
          ineigh=eltoel(inofa,iface)
          if(ineigh<0)then
             write(*,*)'Error in presmoo 2,eltoel not clean'
             stop
          endif
       enddo
    enddo
    !
    !     Filter the ridges touching corner, cusp or ridge points
    !
    nridg0=nridg
    nridg=0_ip

    do iridg=1,nridg0 

       ip1=lridg(1,iridg)
       ip2=lridg(2,iridg)

       if(lptype(1,ip1)==ID_SMOOTH .and. lptype(1,ip2)==ID_SMOOTH)then
          nridg=nridg+1
          lridg(1,nridg)=lridg(1,iridg)
          lridg(2,nridg)=lridg(2,iridg)   
          lridg(3,nridg)=lridg(3,iridg)
          lridg(4,nridg)=lridg(4,iridg)   
       endif
    enddo
    !
    !     DBG
    ! 
    do iridg=1,nridg 
       ip1=lridg(1,iridg)
       ip2=lridg(2,iridg)
       do jridg=iridg+1,nridg 
          ipa=lridg(1,jridg)
          ipb=lridg(2,jridg)

          if((ip1==ipa .and.ip2==ipb) .or. (ip1==ipb .and. ip2==ipa))then 
             write(*,*)'Error in presmoo, edge appears twice'
          endif
       enddo
    enddo
    call outpresmoo(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,lline,&
         lridg,nridg,nline)
    !
    !     Loop on smoothing iterations
    !

    do iter=1,nitermax
       !
       !     Smooth these ridges 
       ! 
       do iridg=1,nridg


          ip1=lridg(1,iridg)
          ip2=lridg(2,iridg)
          rx=coor(1,ip1)-coor(1,ip2)
          ry=coor(2,ip1)-coor(2,ip2)
          rz=coor(3,ip1)-coor(3,ip2)
          rdist12=sqrt(rx*rx+ry*ry+rz*rz)
          !
          !     Compute v by finding both elements containing ip1 and ip2
          !
          neloc=0_ip 
          do ie=ptoel2(ip1),ptoel2(ip1+1)-1
             iface=ptoel1(ie)
             icont=0_ip
             do inofa=1,nnofa
                ipoin=lface(inofa,iface)
                if(ipoin==ip1)then
                   icont=icont+1_ip
                else if(ipoin==ip2)then
                   icont=icont+1_ip
                endif
             enddo
             if(icont==2)then
                neloc=neloc+1  
                leloc(neloc)=iface
             endif
          enddo
          iface1=leloc(1)
          iface2=leloc(2)
          !
          !     Orient iface1 and iface2
          !
          if(lface(1,iface1)==ip1)then
             ipf1=1_ip
          else if(lface(2,iface1)==ip1)then
             ipf1=2_ip
          else
             ipf1=3_ip
          endif
          !
          !     Do we have to switch iface1 and iface2?
          !
          if(lface(ltab(1,ipf1),iface1)==ip2)then

             iface1=leloc(2)
             iface2=leloc(1)

             ipf2=ltab(2,ipf1)

             if(lface(1,iface1)==ip1)then
                ipf1=2_ip
             else if(lface(2,iface1)==ip1)then
                ipf1=3_ip
             else
                ipf1=1_ip
             endif

          else

             ipf1=ltab(1,ipf1)

             if(lface(1,iface2)==ip1)then
                ipf2=3_ip
             else if(lface(2,iface2)==ip1)then
                ipf2=1_ip
             else
                ipf2=2_ip
             endif

          endif
          !
          !     Find the points that are not ip1 or ip2
          !
          ipoin1=lface(ipf1,iface1)
          ipoin2=lface(ipf2,iface2)

          v(1)=coor(1,ipoin2)-coor(1,ipoin1) 
          v(2)=coor(2,ipoin2)-coor(2,ipoin1) 
          v(3)=coor(3,ipoin2)-coor(3,ipoin1) 
          !
          !     Get relevant neighbors and non normalized point normals
          !
          np1=0_ip
          np2=0_ip
          lmark(ip1)=iridg
          lmark(ip2)=iridg
          lmark(ipoin1)=iridg
          lmark(ipoin2)=iridg

          nrmal1(1)=c00
          nrmal1(2)=c00
          nrmal1(3)=c00

          do ie=ptoel2(ip1),ptoel2(ip1+1)-1
             iface=ptoel1(ie)
             do inofa=1,nnofa
                ipoin=lface(inofa,iface)
                if(lmark(ipoin)/=iridg)then
                   np1=np1+1_ip
                   lp1(np1)=ipoin
                   lmark(ipoin)=iridg          
                endif
             enddo

             ipa=lface(1,iface)
             ipb=lface(2,iface)
             ipc=lface(3,iface)

             d1(1)=coor(1,ipb)-coor(1,ipa)
             d1(2)=coor(2,ipb)-coor(2,ipa)
             d1(3)=coor(3,ipb)-coor(3,ipa)
             d2(1)=coor(1,ipc)-coor(1,ipa)
             d2(2)=coor(2,ipc)-coor(2,ipa)
             d2(3)=coor(3,ipc)-coor(3,ipa)

             nrmal1(1)=nrmal1(1)+d1(2)*d2(3)-d1(3)*d2(2)
             nrmal1(2)=nrmal1(2)-d1(1)*d2(3)+d1(3)*d2(1)
             nrmal1(3)=nrmal1(3)+d1(1)*d2(2)-d1(2)*d2(1)

          enddo

          nrmal2(1)=c00
          nrmal2(2)=c00
          nrmal2(3)=c00

          do ie=ptoel2(ip2),ptoel2(ip2+1)-1
             iface=ptoel1(ie)
             do inofa=1,nnofa
                ipoin=lface(inofa,iface)
                if(lmark(ipoin)/=iridg)then
                   np2=np2+1_ip
                   lp2(np2)=ipoin          
                   lmark(ipoin)=iridg          
                endif

                ipa=lface(1,iface)
                ipb=lface(2,iface)
                ipc=lface(3,iface)

                d1(1)=coor(1,ipb)-coor(1,ipa)
                d1(2)=coor(2,ipb)-coor(2,ipa)
                d1(3)=coor(3,ipb)-coor(3,ipa)
                d2(1)=coor(1,ipc)-coor(1,ipa)
                d2(2)=coor(2,ipc)-coor(2,ipa)
                d2(3)=coor(3,ipc)-coor(3,ipa)

                nrmal2(1)=nrmal2(1)+d1(2)*d2(3)-d1(3)*d2(2)
                nrmal2(2)=nrmal2(2)-d1(1)*d2(3)+d1(3)*d2(1)
                nrmal2(3)=nrmal2(3)+d1(1)*d2(2)-d1(2)*d2(1)

             enddo
          enddo
          !
          !     DBG
          !
          !rnl1=sqrt(nrmal1(1)*nrmal1(1)+nrmal1(2)*nrmal1(2)+nrmal1(3)*nrmal1(3))
          !rnl1=c10/rnl1
          !rnl2=sqrt(nrmal2(1)*nrmal2(1)+nrmal2(2)*nrmal2(2)+nrmal2(3)*nrmal2(3))
          !rnl2=c10/rnl2
          !write(*,*)rnl1*nrmal1(1),rnl1*nrmal1(2),rnl1*nrmal1(3)
          !write(*,*)rnl2*nrmal2(1),rnl2*nrmal2(2),rnl2*nrmal2(3)
          !
          !     Add the two common points
          !
          np1=np1+1_ip
          lp1(np1)=ipoin1
          np1=np1+1_ip
          lp1(np1)=ipoin2

          np2=np2+1_ip
          lp2(np2)=ipoin1
          np2=np2+1_ip
          lp2(np2)=ipoin2
          !
          !     Compute distance
          !
          rdist1t=c00
          do iploc=1,np1
             ipoin=lp1(iploc)
             rx=coor(1,ip1)-coor(1,ipoin) 
             ry=coor(2,ip1)-coor(2,ipoin) 
             rz=coor(3,ip1)-coor(3,ipoin) 
             rnl=sqrt(rx*rx+ry*ry+rz*rz)  
             rdist1(iploc)=rnl
             rdist1t=rdist1t+rnl
          enddo
          rdist1t=rdist1t+rdist12
          rdist1t=c10/rdist1t

          rdist2t=c00
          do iploc=1,np2
             ipoin=lp2(iploc)
             rx=coor(1,ip2)-coor(1,ipoin) 
             ry=coor(2,ip2)-coor(2,ipoin) 
             rz=coor(3,ip2)-coor(3,ipoin) 
             rnl=sqrt(rx*rx+ry*ry+rz*rz)  
             rdist2(iploc)=rnl
             rdist2t=rdist2t+rnl
          enddo
          rdist2t=rdist2t+rdist12
          rdist2t=c10/rdist2t
          !
          !     Compute Laplacian smoothing
          !
          x1new(1)=c00
          x1new(2)=c00
          x1new(3)=c00

          do iploc=1,np1
             ipoin=lp1(iploc)
             !rweight=rdist1(iploc)*rdist1t
             rweight=1_ip
             x1new(1)=x1new(1)+rweight*coor(1,ipoin) 
             x1new(2)=x1new(2)+rweight*coor(2,ipoin) 
             x1new(3)=x1new(3)+rweight*coor(3,ipoin) 
          enddo

          x2new(1)=c00
          x2new(2)=c00
          x2new(3)=c00

          do iploc=1,np2
             ipoin=lp2(iploc)
             !rweight=rdist2(iploc)*rdist2t
             rweight=1_ip
             x2new(1)=x2new(1)+rweight*coor(1,ipoin) 
             x2new(2)=x2new(2)+rweight*coor(2,ipoin) 
             x2new(3)=x2new(3)+rweight*coor(3,ipoin) 
          enddo

          !write(*,*)(x1new(1)+coor(1,ip2))/real(np1+1),(x1new(2)+coor(2,ip2))/real(np1+1),(x1new(3)+coor(3,ip2))/real(np1+1)
          !write(*,*)(x2new(1)+coor(1,ip1))/real(np2+1),(x2new(2)+coor(2,ip1))/real(np2+1),(x2new(3)+coor(3,ip1))/real(np2+1)
          !
          !     Compute smoothed coordinates
          !
          !coef1=c10/real(np1*np2-1_ip)
          !coef2=coef1*real(np2)
          !x1s(1)=coef1*x2new(1)+coef2*x1new(1)
          !x1s(2)=coef1*x2new(2)+coef2*x1new(2)
          !x1s(3)=coef1*x2new(3)+coef2*x1new(3)
          !rweight=rdist12*rdist1t
          rweight=1_ip
          x1s(1)=(x1new(1)+rweight*coor(1,ip2))/real(np1+1_ip) 
          x1s(2)=(x1new(2)+rweight*coor(2,ip2))/real(np1+1_ip) 
          x1s(3)=(x1new(3)+rweight*coor(3,ip2))/real(np1+1_ip) 

          !coef3=c10/real(np2)
          !x2s(1)=coef3*(x1s(1)+x2new(1))
          !x2s(2)=coef3*(x1s(2)+x2new(2))
          !x2s(3)=coef3*(x1s(3)+x2new(3))
          !rweight=rdist12*rdist2t
          rweight=1_ip
          x2s(1)=(x2new(1)+rweight*coor(1,ip1))/real(np2+1_ip) 
          x2s(2)=(x2new(2)+rweight*coor(2,ip1))/real(np2+1_ip) 
          x2s(3)=(x2new(3)+rweight*coor(3,ip1))/real(np2+1_ip) 
          !
          !     Compute increments
          !
          dxs1(1)=relax*(x1s(1)-coor(1,ip1)) 
          dxs1(2)=relax*(x1s(2)-coor(2,ip1)) 
          dxs1(3)=relax*(x1s(3)-coor(3,ip1)) 

          dxs2(1)=relax*(x2s(1)-coor(1,ip2)) 
          dxs2(2)=relax*(x2s(2)-coor(2,ip2)) 
          dxs2(3)=relax*(x2s(3)-coor(3,ip2)) 
          !
          !     Compute area vector
          !
          dx(1)=dxs1(1)-dxs2(1)
          dx(2)=dxs1(2)-dxs2(2)
          dx(3)=dxs1(3)-dxs2(3)

          vect(1)= v(2)*dx(3)-v(3)*dx(2) 
          vect(2)=-v(1)*dx(3)+v(3)*dx(1) 
          vect(3)= v(1)*dx(2)-v(2)*dx(1) 

          nrmal(1)=nrmal1(1)+nrmal2(1)+vect(1)
          nrmal(2)=nrmal1(2)+nrmal2(2)+vect(2)
          nrmal(3)=nrmal1(3)+nrmal2(3)+vect(3)

          rnl=sqrt(nrmal(1)*nrmal(1)+nrmal(2)*nrmal(2)+nrmal(3)*nrmal(3))
          !
          !     Test
          !
          if(rnl>tolnorm)then

             rnl=c10/rnl
             nrmal(1)=nrmal(1)*rnl
             nrmal(2)=nrmal(2)*rnl
             nrmal(3)=nrmal(3)*rnl

             csca1=dxs1(1)*nrmal1(1)+dxs1(2)*nrmal1(2)+dxs1(3)*nrmal1(3)  
             csca2=dxs2(1)*nrmal2(1)+dxs2(2)*nrmal2(2)+dxs2(3)*nrmal2(3)  

             vdxs1(1)= v(2)*dxs1(3)-v(3)*dxs1(2) 
             vdxs1(2)=-v(1)*dxs1(3)+v(3)*dxs1(1) 
             vdxs1(3)= v(1)*dxs1(2)-v(2)*dxs1(1) 

             csca3=vdxs1(1)*dxs2(1)+vdxs1(2)*dxs2(2)+vdxs1(3)*dxs2(3)
             h=-(csca1+csca2+csca3)*rnl

             coor(1,ip1)=coor(1,ip1)+dxs1(1)!+h*nrmal(1)
             coor(2,ip1)=coor(2,ip1)+dxs1(2)!+h*nrmal(2)
             coor(3,ip1)=coor(3,ip1)+dxs1(3)!+h*nrmal(3)

             coor(1,ip2)=coor(1,ip2)+dxs2(1)!+h*nrmal(1)
             coor(2,ip2)=coor(2,ip2)+dxs2(2)!+h*nrmal(2)
             coor(3,ip2)=coor(3,ip2)+dxs2(3)!+h*nrmal(3)

          endif

       enddo

    enddo

100 continue
    !
    !     Loop on iteration
    !

    !do iter=1,nitermax
    iter=0_ip
    do
       !
       !     Loop on geometry optimization
       !
       do itergeo=1,nitergeo 

          iter=iter+1_ip
          ncont=0_ip
          !
          !     Loop on faces
          ! 
          do iface=1,nface

             isurf=lsurf(iface)
             ip1=lface(1,iface)          
             ip2=lface(2,iface)          
             ip3=lface(3,iface)          

             rx1=coor(1,ip2)-coor(1,ip1) 
             ry1=coor(2,ip2)-coor(2,ip1) 
             rz1=coor(3,ip2)-coor(3,ip1) 
             rx2=coor(1,ip3)-coor(1,ip1) 
             ry2=coor(2,ip3)-coor(2,ip1) 
             rz2=coor(3,ip3)-coor(3,ip1) 

             rnx1= ry1*rz2-rz1*ry2
             rny1=-rx1*rz2+rz1*rx2
             rnz1= rx1*ry2-ry1*rx2
             rnl=sqrt(rnx1*rnx1+rny1*rny1+rnz1*rnz1)
             rnl=c10/rnl
             rnx1=rnx1*rnl
             rny1=rny1*rnl
             rnz1=rnz1*rnl
             !
             !     Loop on neighbors
             ! 
             do inofa=1,nnofa

                jface=eltoel(inofa,iface)
                !
                !     Do we have a neighbor?
                ! 
                if(jface<=0)cycle
                !
                !     Find the viewer
                !
                if(eltoel(1,jface)==iface)then
                   iview=1 
                else if(eltoel(2,jface)==iface)then
                   iview=2 
                else 
                   iview=3 
                endif
                !
                !     Are we on the same surface?
                !
                jsurf=lsurf(jface)
                if(jsurf/=isurf)then
                   eltoel(inofa,iface)=-eltoel(inofa,iface)
                   eltoel(iview,jface)=-eltoel(iview,jface) 
                   cycle
                endif
                !
                !     It should be a smooth ridge
                !
                ip1=lface(ltab(1,inofa),iface)
                ip2=lface(ltab(2,inofa),iface)
                !if(lptype(1,ip1)/=ID_SMOOTH .or. lptype(1,ip2)/=ID_SMOOTH)cycle
                !
                !     Is it a ridge?
                !
                jp1=lface(1,jface)          
                jp2=lface(2,jface)          
                jp3=lface(3,jface)          

                rx1=coor(1,jp2)-coor(1,jp1) 
                ry1=coor(2,jp2)-coor(2,jp1) 
                rz1=coor(3,jp2)-coor(3,jp1) 
                rx2=coor(1,jp3)-coor(1,jp1) 
                ry2=coor(2,jp3)-coor(2,jp1) 
                rz2=coor(3,jp3)-coor(3,jp1) 

                rnx2= ry1*rz2-rz1*ry2
                rny2=-rx1*rz2+rz1*rx2
                rnz2= rx1*ry2-ry1*rx2
                rnl=sqrt(rnx2*rnx2+rny2*rny2+rnz2*rnz2)
                rnl=c10/rnl
                rnx2=rnx2*rnl
                rny2=rny2*rnl
                rnz2=rnz2*rnl

                cscal=rnx1*rnx2+rny1*rny2+rnz1*rnz2
                !
                !     And the test
                ! 
                if(cscal>tolscal)then
                   !
                   !     Mark the side 
                   !
                   eltoel(inofa,iface)=-eltoel(inofa,iface) 
                   eltoel(iview,jface)=-eltoel(iview,jface) 
                   cycle
                endif
                !
                !     Update counter
                !
                ncont=ncont+1
                !
                !     Get ridge
                !
                rx=coor(1,ip2)-coor(1,ip1)
                ry=coor(2,ip2)-coor(2,ip1)
                rz=coor(3,ip2)-coor(3,ip1)
                rnl=sqrt(rx*rx+ry*ry+rz*rz)
                rnl=c10/rnl
                rx=rx*rnl 
                ry=ry*rnl 
                rz=rz*rnl

                ipoin1=lface(inofa,iface)
                ipoin2=lface(iview,jface)
                !
                !     Get bisectrice of the ridge
                ! 
                rnx=rnx1+rnx2          
                rny=rny1+rny2          
                rnz=rnz1+rnz2          
                rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz) 
                rnl=c10/rnl
                rnx=rnx*rnl
                rny=rny*rnl
                rnz=rnz*rnl
                !
                !     Get mid points
                !
                rx1=c05*(coor(1,ipoin1)+coor(1,ip1))
                ry1=c05*(coor(2,ipoin1)+coor(2,ip1))
                rz1=c05*(coor(3,ipoin1)+coor(3,ip1))
                rx2=c05*(coor(1,ipoin2)+coor(1,ip1))
                ry2=c05*(coor(2,ipoin2)+coor(2,ip1))
                rz2=c05*(coor(3,ipoin2)+coor(3,ip1))
                !
                !     Get tangent
                !
                rtx=rx2-rx1
                rty=ry2-ry1
                rtz=rz2-rz1
                rtl=sqrt(rtx*rtx+rty*rty+rtz*rtz) 
                rtl=c10/rtl
                rtx=rtx*rtl
                rty=rty*rtl
                rtz=rtz*rtl
                !
                !     Get normal to optimal plane passing through both mid points
                !     and ridge
                !
                rnx1= rty*rz-rtz*ry
                rny1=-rtx*rz+rtz*rx
                rnz1= rtx*ry-rty*rx
                rnl=sqrt(rnx1*rnx1+rny1*rny1+rnz1*rnz1)
                rnl=c10/rnl
                rnx1=rnx1*rnl
                rny1=rny1*rnl
                rnz1=rnz1*rnl
                !
                !     Get intersection between this plane and the bisectrice direction
                !  
                rx=rx1-coor(1,ip1)
                ry=ry1-coor(2,ip1)
                rz=rz1-coor(3,ip1)

                lambda=(rx*rnx1+ry*rny1+rz*rnz1)/(rnx*rnx1+rny*rny1+rnz*rnz1)
                !
                !     Now move the four points if possible
                !
                if(lptype(1,ip1)==ID_SMOOTH)then

                   coor(1,ip1)=coor(1,ip1)+lambda*rnx
                   coor(2,ip1)=coor(2,ip1)+lambda*rny
                   coor(3,ip1)=coor(3,ip1)+lambda*rnz
                   !
                   !     Mark all the neighbors
                   !
                   do ie=ptoel2(ip1),ptoel2(ip1+1)-1
                      kface=ptoel1(ie)
                      eltoel(1,kface)=abs(eltoel(1,kface))
                      eltoel(2,kface)=abs(eltoel(2,kface))
                      eltoel(3,kface)=abs(eltoel(3,kface))
                      do knofa=1,nnofa
                         kneigh=abs(eltoel(knofa,kface))
                         if(kneigh==0)cycle 
                         eltoel(1,kneigh)=abs(eltoel(1,kneigh))
                         eltoel(2,kneigh)=abs(eltoel(2,kneigh))
                         eltoel(3,kneigh)=abs(eltoel(3,kneigh))
                      enddo
                   enddo
                endif

                if(lptype(1,ip2)==ID_SMOOTH)then

                   coor(1,ip2)=coor(1,ip2)+lambda*rnx
                   coor(2,ip2)=coor(2,ip2)+lambda*rny
                   coor(3,ip2)=coor(3,ip2)+lambda*rnz
                   !
                   !     Mark all the neighbors
                   !
                   do ie=ptoel2(ip1),ptoel2(ip1+1)-1
                      kface=ptoel1(ie)
                      eltoel(1,kface)=abs(eltoel(1,kface))
                      eltoel(2,kface)=abs(eltoel(2,kface))
                      eltoel(3,kface)=abs(eltoel(3,kface))
                      do knofa=1,nnofa
                         kneigh=abs(eltoel(knofa,kface))
                         if(kneigh==0)cycle 
                         eltoel(1,kneigh)=abs(eltoel(1,kneigh))
                         eltoel(2,kneigh)=abs(eltoel(2,kneigh))
                         eltoel(3,kneigh)=abs(eltoel(3,kneigh))
                      enddo
                   enddo
                endif

                if(lptype(1,ipoin1)==ID_SMOOTH)then

                   coor(1,ipoin1)=coor(1,ipoin1)-lambda*rnx
                   coor(2,ipoin1)=coor(2,ipoin1)-lambda*rny
                   coor(3,ipoin1)=coor(3,ipoin1)-lambda*rnz
                   !
                   !     Mark all the neighbors
                   !
                   do ie=ptoel2(ip1),ptoel2(ip1+1)-1
                      kface=ptoel1(ie)
                      eltoel(1,kface)=abs(eltoel(1,kface))
                      eltoel(2,kface)=abs(eltoel(2,kface))
                      eltoel(3,kface)=abs(eltoel(3,kface))
                      do knofa=1,nnofa
                         kneigh=abs(eltoel(knofa,kface))
                         if(kneigh==0)cycle 
                         eltoel(1,kneigh)=abs(eltoel(1,kneigh))
                         eltoel(2,kneigh)=abs(eltoel(2,kneigh))
                         eltoel(3,kneigh)=abs(eltoel(3,kneigh))
                      enddo
                   enddo
                endif

                if(lptype(1,ipoin2)==ID_SMOOTH)then

                   coor(1,ipoin2)=coor(1,ipoin2)-lambda*rnx
                   coor(2,ipoin2)=coor(2,ipoin2)-lambda*rny
                   coor(3,ipoin2)=coor(3,ipoin2)-lambda*rnz
                   !
                   !     Mark all the neighbors
                   !
                   do ie=ptoel2(ip1),ptoel2(ip1+1)-1
                      kface=ptoel1(ie)
                      eltoel(1,kface)=abs(eltoel(1,kface))
                      eltoel(2,kface)=abs(eltoel(2,kface))
                      eltoel(3,kface)=abs(eltoel(3,kface))
                      do knofa=1,nnofa
                         kneigh=abs(eltoel(knofa,kface))
                         if(kneigh==0)cycle 
                         eltoel(1,kneigh)=abs(eltoel(1,kneigh))
                         eltoel(2,kneigh)=abs(eltoel(2,kneigh))
                         eltoel(3,kneigh)=abs(eltoel(3,kneigh))
                      enddo
                   enddo
                endif

             enddo

          enddo
          !
          !    Output
          !
          write(*,*)'For iter:',iter,' ridges found:',ncont 
       enddo
       if(ncont==0)exit 
       if(iter==nitermax)exit 
       cycle
       !
       !     Now perform a Laplacian smoothing
       !
       do ipoin=1,npoin
          lmark(ipoin)=0_ip
       enddo

       do ipoin=1,npoin

          if(lptype(1,ipoin)/=ID_SMOOTH)cycle
          lmark(ipoin)=ipoin
          ncont=0_ip
          rx=c00
          ry=c00
          rz=c00
          do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
             iface=ptoel1(ie)  
             do inofa=1,nnofa
                jpoin=lface(inofa,iface)
                if(lmark(jpoin)/=ipoin)then
                   lmark(jpoin)=ipoin  
                   rx=rx+coor(1,jpoin)
                   ry=ry+coor(2,jpoin)
                   rz=rz+coor(3,jpoin)
                   ncont=ncont+1 
                endif
             enddo
          enddo
          ncont=max(1_ip,ncont)
          rx=rx/real(ncont)
          ry=ry/real(ncont)
          rz=rz/real(ncont)

          coor(1,ipoin)=rx 
          coor(2,ipoin)=ry 
          coor(3,ipoin)=rz 

       enddo

    enddo
    !
    !     Clean up eltoel
    !
    do iface=1,nface 
       do inofa=1,nnofa
          eltoel(inofa,iface)=abs(eltoel(inofa,iface))
       enddo
    enddo
    !
    !     Now check for anisotropy
    !






    !call outpresmoo(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,lline,&
    !     lridg,nridg,nline)

    !call memchk(2_ip,istat,memor_msh,'LMARKF','presmoo',lmarkf)
    !deallocate(lmarkf,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LMARKF','presmoo',0_ip)
    !call memchk(2_ip,istat,memor_msh,'LRIDG','presmoo',lridg)
    !deallocate(lridg,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LRIDG','presmoo',0_ip)
    !call memchk(2_ip,istat,memor_msh,'LMARK','presmoo',lmark)
    !deallocate(lmark,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LMARK','presmoo',0_ip)

    write(*,*)'Dumping smoothed old mesh in outface.msh'        
    call outface(nnofa,nface,npoin,ndim,lface,coor,lsurf)
    !
    !     Get the face normals
    !
    call gtfnrl(lface,nface,nnofa,ndim,coor,npoin,rnofa)
    !
    !     Get the point normals
    !
    call gtpnrl(nface,rnopo,npoin,ndim,ptoel1,ptoel2,rnofa)

  end subroutine presmoo

  subroutine outpresmoo(lface,nface,nnofa,ndim,coor,npoin,lsurf,nsurf,lside,nside,&
       lline,lridg,nridg,nline)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only       :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH,memor_msh
    use mod_memchk
    implicit none 
    integer(ip),intent(in)    :: nface,nnofa,ndim,npoin,nsurf,nside,nridg,nline 
    integer(ip),intent(in)    :: lface(nnofa,nface),lside(2,nside),lline(nside)
    integer(ip),intent(in)    :: lsurf(nface),lridg(4,nridg)
    real(rp),intent(in)       :: coor(ndim,npoin)
    integer(ip)               :: i,icont

    open(unit=50,file='outpresmoo.msh',status='unknown')
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

    write(50,2)
    write(50,3)
    write(50,4)
    write(50,11)
    write(50,5)

    icont=0_ip
    do i=1, nside
       icont=icont+1
       write(50,400)icont,lside(1,i),lside(2,i),lline(i)
    enddo

    do i=1, nridg
       icont=icont+1
       write(50,400)icont,lridg(1,i),lridg(2,i),nline+1
    enddo



    write(50,6)

    close(50)

1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
11  format('MESH dimension 3 ElemType  Linear Nnode 2')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
300 format(5i10)
400 format(4i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')


  end subroutine outpresmoo

  subroutine outinter(nnofa,nface,npoin,ndim,lface,coor,lsurf,lside,lline,nside,nnosi,nline)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim,nside,nnosi,nline
    integer(ip), intent(in)      :: lface(nnofa,nface),lsurf(nface)
    integer(ip), intent(in)      :: lside(nnosi,nside),lline(nside)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: i,icont
    !
    !     This sub displays a triangular surface mesh with its surface patch number 
    !

    open(unit=50,file='interface.msh',status='unknown')
    rewind 50

    write(50,*)'npoin ',npoin

    do i=1,npoin
       write(50,*)i,coor(1,i),coor(2,i),coor(3,i)
    enddo

    write(50,*)'nline ',nline
    write(50,*)'nside ',nside

    do i=1, nside
       write(50,*)i,lside(1,i),lside(2,i),lline(i)
    enddo

    write(50,*)'nface ',nface

    do i=1, nface
       write(50,*)i,lface(1,i),lface(2,i),lface(3,i),lsurf(i)
    enddo

    close(50)

  end subroutine outinter

  subroutine brutefside(lside,nnosi,nside,pnew,ihost,coor,ndim,npoin,lline,iline)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnosi,nside,ndim,npoin,iline
    integer(ip), intent(in)      :: lside(nnosi,nside),lline(nside)
    integer(ip),intent(inout)    :: ihost
    real(rp), intent(in)         :: coor(ndim,npoin),pnew(ndim)
    integer(ip)                  :: ip1,ip2,iside
    real(rp)                     :: rdist,rdistt,p1(3),p2(3),c10,csca,cscal,epsil0,epsil1,rl
    !
    !     Brute force on sides
    !
    c10=1.0d+00
    rdist=1.0d+12
    epsil1=1.0d+00+1.0d-04
    epsil0=1.0d-04

    do iside=1,nside

       if(lline(iside)/=iline)cycle 

       ip1=lside(1,iside)
       ip2=lside(2,iside)
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

          rdistt=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))

       else if(cscal>epsil1)then

          p2(1)=pnew(1)-coor(1,ip2)
          p2(2)=pnew(2)-coor(2,ip2)
          p2(3)=pnew(3)-coor(3,ip2)
          rdistt=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))

       else

          p2(1)=coor(1,ip1)+csca*p1(1)-pnew(1)
          p2(2)=coor(2,ip1)+csca*p1(2)-pnew(2)
          p2(3)=coor(3,ip1)+csca*p1(3)-pnew(3)
          rdistt=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))

       endif

       if(rdistt<rdist)then
          rdist=rdistt
          ihost=iside
       endif
    enddo

  end subroutine brutefside

end module mod_surf

