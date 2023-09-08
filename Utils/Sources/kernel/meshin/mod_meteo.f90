module mod_meteo

!  use mod_date

contains

  subroutine meteopre(npoin,nelem,ndim,nnofa,nnode,nface,lface,elem,coor,lsurf)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use mod_extrpr
    use mod_interp
    use def_meshin, only       : memor_msh
  use mod_ecoute, only :  ecoute
    use def_meteo
    implicit none
    integer(ip),intent(in)      :: ndim,nnofa,nnode
    integer(ip),intent(inout)   :: npoin,nelem,nface
    real(rp)                    :: rlat1,rlat2,rlon1,rlon2,d1,d2,d3,d4
    real(rp)                    :: thick,reason,rsize1,phi0,theta0,dphi0,dtheta0 
    real(rp)                    :: phi0new,theta0new,dphinew,c00,c10 
    integer(ip)                 :: ipoin,iy,ix,npbdr,invert,iboup
    integer(4)                  :: istat,idiff4
    integer(ip)                 :: ipold,nvar,ip1,ip2,ip3,ip4,ivar,ihost,npold0
    integer(ip)                 :: timi,idiff,niter,interval_seconds,nmod,iter
    type(meteo)                 :: ldata 
    character (len=256)         :: level,datestr,filestr,geog_data_path
    character (len=256)         :: start_date,end_date
    character (len=19)          :: valid_date,temp_date
    !
    !     Arrays for the new mesh
    !
    integer(ip)                 :: nblay,nboup,nx,ny,npsur 
    integer(ip),pointer         :: lface(:,:),elem(:,:),lboup(:),lelem(:)
    integer(ip),pointer         :: lsurf(:),lsurfold(:)
    real(rp),pointer            :: coor(:,:),rblay(:),rblold(:)
    real(rp),pointer            :: var(:,:),pres(:),temp(:),height(:)
    real(rp),pointer            :: uvel(:),vvel(:),wvel(:),rint(:,:)
    !
    !     Arrays for the old mesh
    !
    integer(ip)                 :: npold,neold,nxold,nyold,npsold,nfold,nblold,npsold0
    integer(ip)                 :: nfold0,nface0 
    integer(ip),pointer         :: elold(:,:),lfold(:,:),lmark(:)
    real(rp),pointer            :: coold(:,:)
    real(rp),pointer            :: varold(:,:)
    !
    !     This routine is the main routine for treating the pre-processing 
    !     of the meteo files in Alya
    ! 
    !
    c00=0.0d+00
    c10=1.0d+00
    !
    !     Call openfile Meteo
    !
    call openfi(9_ip)
    !
    !     Read the bl distribution
    !
    call readbl(nblay,rblay)
    !
    !     First read user file to define the area under study
    ! 
    call readmto(rlat1,rlat2,rlon1,rlon2,level,npbdr,filestr,geog_data_path,start_date,end_date,interval_seconds)
    !
    !     Then read topology file and fill coor
    ! 
    call readgeog(nx,ny,npsur,ndim,rlat1,rlat2,rlon1,rlon2,level,npbdr,coor,geog_data_path,phi0new,theta0new,dphinew)
    !
    !     Generate new surface mesh
    !
    call quatri(nx,ny,ndim,npsur,nface,nnofa,coor,lface)
    npoin=npsur
    !
    !     Allocate lsurf
    !
    allocate(lsurf(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSURF','meteo',lsurf)
    !
    !     The face number will change in blmesh, remember it
    !
    nface0=nface
    !
    !     Generate new BL mesh
    !
    call blmesh( nface,nnofa,nnode,ndim,nelem,npoin,nblay,npsur,&
         thick,reason,rsize1,lface,coor,elem,rblay,lsurf)
    call outfac(nnofa,nface,npoin,ndim,lface,coor)
    !
    !     Renumber the new mesh in a finite difference like way
    !
    call renustr(elem,nelem,nnode,npsur,nface0,npoin,nblay,coor,ndim,lface,nnofa,nface)
    !call outgidvolume(npoin,ndim,coor,elem,nelem,nnode)  
    !
    !     Get the boundary points
    !     
    call gtbpoin(elem,nelem,nnode,npoin,nboup,lboup)
    !
    !     Then read the variables and fill coold 
    !
    !datestr='2008-09-15_00'
    !filestr='/home/bsc-user/Soft/WPS/data/FILE'
    valid_date=trim(start_date)
    if (mod(interval_seconds,3600_ip) == 0) then
       write(temp_date,'(a13)') valid_date(1:10)//'_'//valid_date(12:13)
    else if (mod(interval_seconds,60_ip) == 0) then
       write(temp_date,'(a16)') valid_date(1:10)//'_'//valid_date(12:16)
    else
       write(temp_date,'(a19)') valid_date(1:10)//'_'//valid_date(12:19)
    end if
    write(*,*)'valid_date=',temp_date 

    datestr=temp_date 
    write(*,*)'datestr=',datestr 


    call gtmtodata(uvel,vvel,wvel,temp,height,pres,coold,nxold,nyold,&
         nblold,phi0,theta0,ndim,npold,rblold,invert,datestr,filestr,dphi0,dtheta0)
    npsold=nxold*nyold
    !
    !     Account for the changes at the surface in temperature and velocity
    ! 
    call corsurf(nxold,nyold,npsold,npold,ndim,coold,uvel,vvel,temp,phi0,theta0,dphi0,dtheta0)
    !
    !     Allocate old variables
    !
    nvar=5_ip
    allocate(varold(nvar,npold),stat=istat)
    call memchk(zero,istat,memor_msh,'VAROLD','meteo',varold)
    !
    !     Copy from single arrays to varold
    !
    do ipold=1,npold
       varold(1,ipold)=uvel(ipold)
       varold(2,ipold)=vvel(ipold)
       varold(3,ipold)=wvel(ipold)
       varold(4,ipold)=temp(ipold)
       varold(5,ipold)=pres(ipold)
    enddo
    !
    !     Build the old surface mesh 
    !
    call quatri(nxold,nyold,ndim,npsold,nfold,nnofa,coold,lfold)
    !
    !     Do we have to invert the orientation?
    !
    if(invert==1)then
       call invfac(lfold,nfold,nnofa)
    endif
    !
    !     Allocate lmark
    ! 
    allocate(lmark(npsold),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','filmet',lmark)
    !
    !     Remember npsold for filmet2
    !
    npsold0=npsold
    !
    !     Filter the old surface mesh
    !
    call filmet(lfold,nfold,nnofa,coold,ndim,npold,npsold,varold,nvar,height,coor,npoin,lmark,&
         nxold,nyold,phi0,theta0,dphi0,dtheta0)
    !call outfac2(nnofa,nface,npsur,ndim,lface,coor,lfold,nfold,coold,npold)
    !call outfac(nnofa,nfold,npold,ndim,lfold,coold)
    !
    !     Allocate lsurfold
    !
    allocate(lsurfold(nfold),stat=istat)
    call memchk(zero,istat,memor_msh,'LSURFOLD','meteo',lsurfold)
    !
    !     The face number will change in blmesh, remember it
    !
    nfold0=nfold
    !
    !     Generate the old BL mesh
    !
    call blmesh(nfold,nnofa,nnode,ndim,neold,npold,nblold,npsold,&
         thick,reason,rsize1,lfold,coold,elold,rblay,lsurfold)
    !
    !     Renumber the old mesh in a finite difference like way
    !
    call renustr(elold,neold,nnode,npsold,nfold0,npold,nblold,coold,ndim,lfold,nnofa,nfold)
    !call outgidvolume(npold,ndim,coold,elold,neold,nnode)  
    !
    !     Modify the old mesh as the coordinates are imposed
    !     for the points outside the surface 
    ! 
    call modmsh(coold,ndim,npold,height,nblold,npsold)
    !
    !     Output old mesh with results
    !
    !call outgidmeteo(npold,ndim,coold,elold,neold,nnode,varold,nvar)
    !
    !     Output both meshes
    !
    !call outgidvolume2(npoin,ndim,coor,elem,nelem,nnode,coold,npold,elold,neold)
    !
    !     Before interpolating put meshes on a plane
    !
    call sph2pla(coold,npold,phi0,theta0,dphi0,dtheta0,nxold,nyold,npsold,ndim,nblold)
    call sph2pla(coor,npoin,phi0new,theta0new,dphinew,dphinew,nx,ny,npsur,ndim,nblay)
    !
    !     Output old mesh with results
    !
    !call outgidmeteo(npold,ndim,coold,elold,neold,nnode,varold,nvar)
    !
    !     Output both meshes
    !
    !call outgidvolume2(npoin,ndim,coor,elem,nelem,nnode,coold,npold,elold,neold)
    !
    !     Allocate interpolation arrays
    ! 
    allocate(lelem(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEM','meteo',lelem)
    allocate(rint(4,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'RINT','meteo',rint)
    !
    !     Interpolate from the old mesh to the new mesh
    ! 
    call inter(elem,elold,nelem,nnode,neold,npoin,npold,coor,coold,ndim,lboup,&
         nboup,lelem,rint)
    !
    !     Put the new mesh back on earth
    ! 
    call pla2sph(coor,npoin,phi0new,theta0new,dphinew,dphinew,nx,ny,npsur,ndim,nblay)
    ! 
    !     Allocate new variables
    !
    allocate(var(nvar,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'VAR','meteo',var)
    !
    !     Now interpolate the values
    !
    do ipoin=1,npoin

       ihost=lelem(ipoin)
       !
       !     Do we have a valid element?
       !
       if(ihost>0)then
          !
          !     We interpolate from the old element
          !    
          ip1=elold(1,ihost)
          ip2=elold(2,ihost)
          ip3=elold(3,ihost)
          ip4=elold(4,ihost)

          d1=min(c10,max(c00,rint(1,ipoin)))
          d2=min(c10,max(c00,rint(2,ipoin)))
          d3=min(c10,max(c00,rint(3,ipoin)))
          d4=min(c10,max(c00,rint(4,ipoin)))

          do ivar=1,nvar 
             var(ivar,ipoin)=d1*varold(ivar,ip1)+d2*varold(ivar,ip2) &
                  &                 +d3*varold(ivar,ip3)+d4*varold(ivar,ip4)
          enddo

       else
          !
          !     We interpolate from the nearest neighbor
          !
          ipold=-lelem(ipoin)

          do ivar=1,nvar 
             var(ivar,ipoin)=varold(ivar,ipold)
          enddo

       endif

    enddo
    !
    !     Output new mesh with results
    !
    call outgidmeteo(npoin,ndim,coor,elem,nelem,nnode,var,nvar)
    !
    !     Output meteo
    !
    call outmeteo(npoin,ndim,coor,elem,nelem,nnode,var,nvar,nx,ny,npsur,nblay,lface,nface,nnofa,lsurf)
    !
    !     Check input
    !  
    idiff4 = idiff
   !!!!!!!!!!!!!!!!!!!!!!!!!!call geth_idts(end_date, start_date, idiff4)
    if(idiff<0)then
       write(*,*)'Error in dates, Ending date is earlier than starting date'
       stop
    endif
    niter = idiff / interval_seconds
    !
    !     Should be exact
    !
    nmod=niter*interval_seconds-idiff
    if(nmod/=0)then
       write(*,*)'Computation time is not a multiple of the interval seconds '
       stop
    endif
    !
    !     Now create the boundary conditions by looping in time
    !
    timi=interval_seconds

    do iter=1,niter 
       !
       !     Remember old npold
       !
       npold0=npold  
       !!call geth_newdate(valid_date,trim(start_date), time)
       write(*,*)'valid_date=',valid_date 
       if (mod(interval_seconds,3600_ip) == 0) then
          write(temp_date,'(a13)') valid_date(1:10)//'_'//valid_date(12:13)
       else if (mod(interval_seconds,60_ip) == 0) then
          write(temp_date,'(a16)') valid_date(1:10)//'_'//valid_date(12:16)
       else
          write(temp_date,'(a19)') valid_date(1:10)//'_'//valid_date(12:19)
       end if
       write(*,*)'valid_date=',temp_date 

       datestr=temp_date 
       write(*,*)'datestr=',datestr 
       call gtmtodata(uvel,vvel,wvel,temp,height,pres,coold,nxold,nyold,&
            nblold,phi0,theta0,ndim,npold,rblold,invert,datestr,filestr,dphi0,dtheta0)
       npsold=nxold*nyold
       if(npsold/=npsold0)then
          write(*,*)'Error in gtmtodata, files have different sizes, npsold/=npsold0'
          stop 
       endif
       !
       !     Account for the changes at the surface in temperature and velocity
       ! 
       call corsurf(nxold,nyold,npsold,npold,ndim,coold,uvel,vvel,temp,phi0,theta0,dphi0,dtheta0)
       !
       !     Copy from single arrays to varold
       !
       do ipold=1,npold
          varold(1,ipold)=uvel(ipold)
          varold(2,ipold)=vvel(ipold)
          varold(3,ipold)=wvel(ipold)
          varold(4,ipold)=temp(ipold)
          varold(5,ipold)=pres(ipold)
       enddo
       !
       !     Filter the points with the previous pattern
       !
       call filmet2(coold,ndim,npold,varold,nvar,height,lmark,npsold0)

       if(npold/=npold0)then
          write(*,*)'Error in gtmtodata, files have different sizes after filmet2'
          stop 
       endif
       !
       !     Now interpolate the boundary values
       !
       do iboup=1,nboup

          ipoin=lboup(iboup)
          ihost=lelem(ipoin)
          !
          !     Do we have a valid element?
          !
          if(ihost>0)then
             !
             !     We interpolate from the old element
             !    
             ip1=elold(1,ihost)
             ip2=elold(2,ihost)
             ip3=elold(3,ihost)
             ip4=elold(4,ihost)

             d1=rint(1,ipoin)
             d2=rint(2,ipoin)
             d3=rint(3,ipoin)
             d4=rint(4,ipoin)

             do ivar=1,nvar 
                var(ivar,ipoin)=d1*varold(ivar,ip1)+d2*varold(ivar,ip2) &
                     &                 +d3*varold(ivar,ip3)+d4*varold(ivar,ip4)
             enddo

          else
             !
             !     We interpolate from the nearest neighbor
             !
             ipold=-lelem(ipoin)

             do ivar=1,nvar 
                var(ivar,ipoin)=varold(ivar,ipold)
             enddo

          endif

       enddo
       !
       !     Output for gid
       !
       call outgidmeteobc(npoin,ndim,coor,var,nvar,lboup,nboup)
       !
       !     Output the boundary conditions
       !
       call outmetbc(npoin,var,nvar,nboup,lboup,iter)

       timi=timi+interval_seconds

    enddo
    !
    !     Deallocate
    !
    call memchk(2_ip,istat,memor_msh,'LSURFOLD','meteopre',lsurfold)
    deallocate(lsurfold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSURFOLD','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','meteopre',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'HEIGHT','meteopre',height)
    deallocate(height,stat=istat)
    if(istat/=0) call memerr(2_ip,'HEIGHT','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'LELEM','meteopre',lelem)
    deallocate(lelem,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEM','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'RINT','meteopre',rint)
    deallocate(rint,stat=istat)
    if(istat/=0) call memerr(2_ip,'RINT','meteopre',0_ip)
    !call memchk(2_ip,istat,memor_msh,'ELEM','meteopre',elem)
    !deallocate(elem,stat=istat)
    !if(istat/=0) call memerr(2_ip,'ELEM','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELOLD','meteopre',elold)
    deallocate(elold,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELOLD','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFOLD','meteopre',lfold)
    deallocate(lfold,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFOLD','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'COOLD','meteopre',coold)
    deallocate(coold,stat=istat)
    if(istat/=0) call memerr(2_ip,'COOLD','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'VAR','meteopre',var)
    deallocate(var,stat=istat)
    if(istat/=0) call memerr(2_ip,'VAR','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'VAROLD','meteopre',varold)
    deallocate(varold,stat=istat)
    if(istat/=0) call memerr(2_ip,'VAROLD','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'LBOUP','meteopre',lboup)
    deallocate(lboup,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBOUP','meteopre',0_ip)
    !call memchk(2_ip,istat,memor_msh,'LFACE','meteopre',lface)
    !deallocate(lface,stat=istat)
    !if(istat/=0) call memerr(2_ip,'LFACE','meteopre',0_ip)
    !call memchk(2_ip,istat,memor_msh,'COOR','meteopre',coor)
    !deallocate(coor,stat=istat)
    !if(istat/=0) call memerr(2_ip,'COOR','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'RBLAY','meteopre',rblay)
    deallocate(rblay,stat=istat)
    if(istat/=0) call memerr(2_ip,'RBLAY','meteopre',0_ip)
    call memchk(2_ip,istat,memor_msh,'RBLOLD','meteopre',rblold)
    deallocate(rblold,stat=istat)
    if(istat/=0) call memerr(2_ip,'RBLOLD','meteopre',0_ip)

  end subroutine meteopre

  subroutine readgeog(nx,ny,npoin,ndim,rlat1,rlat2,rlon1,rlon2,level,npbdr,coor,&
       geog_data_path,phi0new,theta0new,dphinew)
    use def_parame, only       :  pi
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_master
    use def_meteo, only           : EARTH_RADIUS_M
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(inout)         :: nx,ny,npoin
    integer(ip),intent(in)            :: ndim,npbdr
    real(rp),intent(inout)            :: rlat1,rlat2,rlon1,rlon2,phi0new,dphinew,theta0new
    real(rp),pointer                  :: coor(:,:)     
    character (len=256),intent(in)    :: level,geog_data_path

    integer(ip)               :: istatus,local_endian,local_wordsize,strlen,iplace,jplace
    integer(ip)               :: xdim,ydim,zdim,sign_convention,ix,iy,iz,ipoin
    integer(ip)               :: ilat1,ilat2,ilon1,ilon2,idx,ndx,ndy,i,j,ndxf,ndyf
    integer(ip)               :: ilat1f,ilat2f,ilon1f,ilon2f,ibegin,iend,jbegin,jend
    integer(ip)               :: iname1,iname2,jname1,jname2,indx,jndx,nxl,nyl
    real(rp)                  :: radio,radioloc,phi,theta,dtheta,dphi,c00,theta0,phi0,c180,c360,c90
    real(rp),pointer          :: height(:,:,:)
    real                      :: scalefac,dx,dxf,resol,rminute
    integer(4)                :: istat
    integer(ip),pointer       :: test(:)

    character (len=256)    :: file_name
    character (len=256)    :: file_root

    !
    !   This sub initialize the arrays of points given the user specifications 
    !
    c00=0.0d+00
    c360=360.0d+00
    c180=180.0d+00
    c90=90.0d+00
    !
    !     Create root filename
    !   
    file_root=geog_data_path
    !file_root='/home/bsc-user/Soft/WPS/geog/topo'
    !
    !     Set level: four possibilities  10m, 5m, 1m 30s 
    !     Initialize constants for each level
    !
    if(trim(level)=='10m')then
       idx=120_ip
       dxf=20.0d+00
       dx=360.0d+00/2160.0d+00
       ndxf=18_ip
       ndyf=9_ip
       ndx=2160_ip
       ndy=1080_ip
       file_root=adjustl(trim(file_root))//'_10m/'
       resol=10.0d+00
    else if(trim(level)=='5m')then
       idx=240_ip
       dxf=20.0d+00
       dx=360.0d+00/4320.0d+00
       ndxf=18_ip
       ndyf=9_ip
       ndx=4320_ip
       ndy=2160_ip
       file_root=adjustl(trim(file_root))//'_5m/'
       resol=5.0d+00
    else if(trim(level)=='2m')then
       idx=600_ip
       dxf=20.0d+00
       dx=360.0d+00/10800.0d+00
       ndxf=18_ip
       ndyf=9_ip
       ndx=10800_ip
       ndy=5400_ip
       file_root=adjustl(trim(file_root))//'_2m/'
       resol=1.0d+00
    else if(trim(level)=='30s')then
       idx=1200_ip
       dxf=10.0d+00
       dx=360.0d+00/43200.0d+00
       ndxf=36_ip
       ndyf=18_ip
       ndx=43200_ip
       ndy=21600_ip
       file_root=adjustl(trim(file_root))//'_30s/'
       resol=1.0d+00/60.0d+00*30.0d+00
    else
       write(*,*)'Error in readgeog, topography not defined'
       stop  
    endif
    !
    !     Change Latitude from [-90:90] to [0:180] 
    !
    rlat1=rlat1+c90
    rlat2=rlat2+c90
    !
    !     Change longitude as the 0 of the data corresponds to -180 ??? 
    !
    rlon1=rlon1+c180
    rlon2=rlon2+c180
    !
    !     Check input
    !
    if(rlat1<c00 .or. rlat2<c00 .or. rlon1<c00 .or. rlon2<c00)then
       write(*,*)'Error meteo input 1'
       stop
    endif

    if(rlat1>c180 .or. rlat2>c180 .or. rlon1>c360 .or. rlon2>c360)then
       write(*,*)'Error meteo input 2'
       stop
    endif
    !
    !     Compute position in bin of points
    !
    ilat1=int(rlat1/dx)+1_ip
    ilat2=int(rlat2/dx)+1_ip
    ilon1=int(rlon1/dx)+1_ip
    ilon2=int(rlon2/dx)+1_ip
    !
    !    Compute position in bin of files
    !
    ilat1f=(ilat1-1_ip)/idx+1_ip
    ilat2f=(ilat2-1_ip)/idx+1_ip
    ilon1f=(ilon1-1_ip)/idx+1_ip
    ilon2f=(ilon2-1_ip)/idx+1_ip
    !
    !     Arguments for readgeogrid
    !
    xdim=idx+2_ip*npbdr
    ydim=xdim
    zdim=1_ip
    local_endian=0_ip
    local_wordsize=2_ip
    sign_convention=1_ip
    istatus=1_ip 
    scalefac = 1.0d+00
    !
    !     Allocate height of the points
    !
    allocate(height(xdim,ydim,zdim),stat=istat)
    !
    !     Cartesian to spherical
    !
    rminute=21600.0d+00
    dphi=resol*2.0d+00*pi/rminute     ! 21600 mn ---> 2 Pi ,  resol --->  dphi
    dtheta=dphi
    !
    !     Remember dphinew, theta0new, phi0new
    !
    dphinew=dphi
    theta0new=(ilon1-1)*dphi  
    phi0new=(ilat1-1)*dphi-pi     
    !
    !     Earth radius 
    !
    radio=EARTH_RADIUS_M
    iz=1_ip
    !
    !     Area dimension in bin of points needed to allocate coor
    !
    nx=ilon2-ilon1+1_ip
    ny=ilat2-ilat1+1_ip
    npoin=nx*ny
    !
    !     Reallocate coor
    !
    if(.not.associated(coor)) then
       allocate(coor(ndim,npoin),stat=istat)
       call memchk(zero,istat,memor_msh,'COOR','readgeo',coor)
    else
       call memrea(npoin,memor_msh,'COOR','readgeo',coor)
    endif
    allocate(test(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'TEST','readgeo',test)
    !
    !     Loop on files to be read
    !
    do i=ilon1f,ilon2f
       do j=ilat1f,ilat2f 
          !
          !     Create the file name
          ! 
          iname1=(i-1)*idx+1
          iname2=i*idx 
          jname1=(j-1)*idx+1
          jname2=j*idx
          write(file_name,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)')trim(file_root),iname1,'-',iname2,'.',jname1,'-',jname2
          write(*,*)'Reading file:',trim(file_name)
          !
          !     Get the data slice
          !
          strlen=len(trim(file_name))

!!!!!     WARNING: THIS CALL IS SKIPPED BECAUSE THE ENDIAN AFFAIRE.
!!          call read_geogrid(file_name, strlen, height, xdim, ydim, zdim, &
!!               sign_convention, local_endian, scalefac, local_wordsize, istatus)
!!  HEMOS VOLADO LA SUBRUTINE READ_GEOGRID: ESTA PUESTA AL FINAL DE ESTE MODULO POR SI LAS MOSCAS

!!!!
          call runend("MOD_METEO: SOLVE THE READ_GEOGRID ENDIAN PROBLEM FIRST!!")
!!!!
          !
          !     Could we read the file?
          !
          if(istatus/=0)then
             write(*,*)'File not found in readgeog'
             write(*,*)trim(file_name)
             stop
          endif
          !
          !     Compute box intersection between file and area chosen
          !
          if(iname1<ilon1)then
             ibegin=ilon1-iname1+1_ip
          else
             ibegin=1_ip
          endif

          if(iname2>ilon2)then
             iend=ilon2-iname1+1_ip
          else
             iend=idx
          endif

          if(jname1<ilat1)then
             jbegin=ilat1-jname1+1_ip
          else
             jbegin=1_ip
          endif

          if(jname2>ilat2)then
             jend=ilat2-jname1+1_ip
          else
             jend=idx
          endif
          !
          !     Absolute index of the beginning of this slice
          ! 
          indx=(i-1)*idx
          jndx=(j-1)*idx 
          !
          !     Initial angles of the current file
          !
          theta0=(indx+ibegin-1)*dtheta
          phi0=(jndx+jbegin-1)*dphi-pi
          !
          !     Switch from spherical to cartesian coordinate
          ! 
          phi=phi0
          !
          !     Reassign the points in their respective place
          !  
          do iy=jbegin,jend

             theta=theta0

             do ix=ibegin,iend 

                radioloc=radio+height(ix+npbdr,iy+npbdr,iz)
                !
                !     Position in the bin of points
                !
                iplace=indx+ix
                jplace=jndx+iy
                !
                !     Position in the array
                !
                ipoin=(iplace-ilon1+1)+(jplace-ilat1)*nx
                test(ipoin)=1_ip
                !write(*,*)ipoin,iplace-ilon1+1,jplace-ilat1,i,j
                coor(1,ipoin)=cos(theta)*sin(phi)*radioloc
                coor(2,ipoin)=sin(theta)*sin(phi)*radioloc
                coor(3,ipoin)=cos(phi)*radioloc
                theta=theta+dtheta

             enddo

             phi=phi+dphi

          enddo
       enddo
    enddo

    !
    !     DBG 
    !
    do ipoin=1,npoin
       if(test(ipoin)==0)then
          write(*,*)'Error test resdgeog'
          stop
       endif
    enddo

    deallocate(test,stat=istat)


!!$  file_name="/home/bsc-user/Soft/WPS/geog/topo_10m/00481-00600.00601-00720"
!!$  c00=0.0d+00
!!$  xdim=126
!!$  ydim=126
!!$  zdim=1
!!$  local_endian=0
!!$  local_wordsize=2
!!$  sign_convention=1
!!$  istatus=1 
!!$  scalefac = 1.0
!!$
!!$  !
!!$  !     Allocate height of the points
!!$  !
!!$
!!$  allocate(height(xdim,ydim,zdim),stat=istat)
!!$  !
!!$  !     Read in binary file
!!$  !
!!$  strlen=len(trim(file_name))
!!$  call read_geogrid(file_name, strlen, height, xdim, ydim, zdim, &
!!$       sign_convention, local_endian, scalefac, local_wordsize, istatus)
!!$
!!$  !
!!$  !     Earth radius 
!!$  !
!!$
!!$  radio=6356.00d+03
!!$  iz=1_ip
!!$  nx=120_ip
!!$  ny=120_ip
!!$  ipoin=0_ip
!!$  npoin=nx*ny
!!$  !
!!$  !     Reallocate coor
!!$  !
!!$
!!$  if(associated(coor)==.false.)then
!!$     allocate(coor(ndim,npoin),stat=istat)
!!$     call memchk(zero,istat,memor_msh,'COOR','readgeo',coor)
!!$  else
!!$     call memrea(npoin,memor_msh,'COOR','readgeo',coor)
!!$  endif
!!$
!!$
!!$  !
!!$  !     Switch from spherical to cartesian coorinate
!!$  ! 
!!$  !
!!$  !     For the moment hard coded
!!$  !  
!!$
!!$  dphi=10*2*pi/21600.0d+00     ! 21600 m ---> 2 Pi ,  10mn --->
!!$  dtheta=dphi
!!$
!!$  theta0=c00
!!$  phi0=c00
!!$
!!$  phi=-pi+601*dphi
!!$  do iy=1,ny 
!!$
!!$     theta=theta0
!!$
!!$     do ix=1,nx 
!!$
!!$        ipoin=ipoin+1 
!!$        radioloc=radio+height(ix+3,iy+3,iz)
!!$        !write(*,*)ix,iy,height(ix+3,iy+3,iz),ipoin
!!$        !radioloc=radio
!!$        coor(1,ipoin)=cos(theta)*sin(phi)*radioloc
!!$        coor(2,ipoin)=sin(theta)*sin(phi)*radioloc
!!$        coor(3,ipoin)=cos(phi)*radioloc
!!$        theta=theta+dtheta
!!$     enddo
!!$     phi=phi+dphi
!!$  enddo
!!$



    deallocate(height,stat=istat)


  end subroutine readgeog

  subroutine readmto(rlat1,rlat2,rlon1,rlon2,level,npbdr,filestr,geog_data_path,start_date,end_date,interval_seconds)
    use def_kintyp, only       :  ip,rp,lg
    use      def_master
    use      def_inpout
    implicit none
    real(rp),intent(inout)    :: rlat1,rlat2,rlon1,rlon2
    integer(ip),intent(inout) :: npbdr,interval_seconds
    character (len=256),intent(inout)    :: level,filestr,geog_data_path,start_date,end_date
    character (len=256)       :: text
    !
    lispa         = 0
    !lisda         = lun_meteo                              ! Reading file
    lisre         = lun_outpu                              ! Writing file


    !
    ! Reach the section
    !
    !call ecoute('RMETEO')
    !
    ! Read data
    !
    !do while(words(1)/='ENDME')
    !   call ecoute('RMETEO')

    !   if(words(1)=='LATIT') then                          ! Latitude

    !      rlat1=param(1)
    !      rlat2=param(2)

    !   else if(words(1)=='LONGI') then                     ! Longitude

    !      rlon1=param(1)
    !      rlon2=param(2)

    !   else if(words(1)=='TOPOL') then                     ! Topology level

    !      level=words(2)

    !   else if(words(1)=='BOUND') then                     ! Boundary ghost points in the file

    !      npbdr=param(1)

    !   else if(words(1)=='FILEN')then

    !      filestr=words(2) 

    !   else if(words(1)=='GEOGD')then

    !     geog_data_path=words(2) 

    !   else if(words(1)=='START')then

    !      start_date=words(2)

    !   else if(words(1)=='ENDDA')then

    !      end_date=words(2) 

    !   endif
    !enddo

    !read(lun_meteo,*)text
    !read(lun_meteo,*)text,rlat1,rlat2
    !read(lun_meteo,*)text,rlon1,rlon2
    !read(lun_meteo,*)text,level 
    !read(lun_meteo,*)text,npbdr
    !read(lun_meteo,*)text,filestr
    !read(lun_meteo,*)text,geog_data_path
    !read(lun_meteo,*)text,start_date
    !read(lun_meteo,*)text,end_date
    !read(lun_meteo,*)text,interval_seconds
    !read(lun_meteo,*)text

  end subroutine readmto

  subroutine gtmtodata(uvel,vvel,wvel,temp,height,pres,coor,npx,npy,nblay,phi0,theta0,ndim,npoin,&
       rblay,invert,datestr,filestr,dphi0,dtheta0)
    use def_kintyp, only       :  ip,rp,lg
    use def_parame, only       :  pi
    use def_meteo, only          :  meteo,memor_meteo,EARTH_RADIUS_M  
    use mod_memchk
    implicit none

    ! Arguments
    real(rp),pointer            :: uvel(:),vvel(:),wvel(:),rblay(:)
    real(rp),pointer            :: temp(:),pres(:),coor(:,:),height(:)
    integer(ip),intent(inout)   :: nblay,npoin,invert
    integer(ip),intent(in)      :: ndim
    real(rp),intent(inout)      :: phi0,theta0,dphi0,dtheta0 
    character (len=256),intent(in):: datestr,filestr
    integer(ip)                 :: iter,istatus,npx,npxt,npy,npyt,npsur,npsurt
    integer(ip),parameter       :: npdim=6_ip
    integer(ip)                 :: np(npdim),ipdim,npoint,npnew,ix,iy,ilay,ipoin
    integer(ip)                 :: npcoor,iblay,iheight,input_unit,iostatus
    real(rp)                    :: radio,radioloc,theta,phi,dtheta,dphi,uv,vv,press 
    real(rp),pointer            :: pressurf(:)
    real(rp)                    :: pres0,temp0,gamar,z0,c10,R,coef1,coef2,grav,c00,deg2rad 
    type (meteo)                :: ldata
    character (len=256)         :: filename
    integer(4)                :: istat
    !
    !     Initialize earth radius
    !
    radio=EARTH_RADIUS_M
    !
    !     Initialize nblay 
    !
    nblay=0_ip 
    !
    !     Initialize invert 
    !
    invert=0_ip 
    !
    !     Initialize iheight 
    !
    iheight=0_ip 
    !
    !     Initialize npcoor 
    !
    npcoor=0_ip 
    !
    !     Initialize physical constants
    !
    c00=0.0d+00
    c10=1.0d+00
    pres0=1.0d+06
    temp0=300
    gamar=1.4d+00
    R=287.0d+00
    grav=9.81d+00
    z0=R*temp0/grav
    coef1=z0*gamar/(gamar-c10)
    coef2=(gamar-c10)/gamar
    deg2rad=pi/real(180)
    !
    !     Open file
    !    
    write(filename, '(a)') trim(filestr)//':'//trim(datestr)
    write(*,*)filename
    input_unit=53
    !open(unit=input_unit,convert='big_endian',file=trim(filename),status='old',form='unformatted',iostat=iostatus)
    if(iostatus>0)then
       write(*,*)'File not found in read_next_met_field'
       write(*,*)trim(filestr) 
       stop
    endif
    !
    !     Loop on the data until EOF found 
    !
    istatus = 0_ip 
    iter = 0_ip
    do ipdim=1,npdim
       np(ipdim)=0_ip
    enddo

    do while(istatus==0)

       iter=iter+1_ip

       call read_next_met_field(ldata,istatus,filename,input_unit)

       write(*,*)ldata%field


       if(istatus==0)then
          !
          !     What kind of data do we have?
          !
          !
          !     For the moment, only consider:
          !
          !     For 3D data
          ! 
          !        pressure:   
          !        height: HGT
          !        U velocity: UU 
          !        V velocity: VV 
          !        Temperature: TT 
          !
          !
          !    For 2D data
          !  
          !    Surface Pressure: PSFC
          !    Terrain height: SOILHGT
          !    Temperature at 2m
          !    U velocity at 10m 
          !    V velocity at 10m 
          !
          !     The first U and V velocities are measured  at 10 m
          !     The first temperatures are measured at 2 m
          !
          npxt=ldata%nx
          npyt=ldata%ny
          npsurt=npxt*npyt
          press=ldata%xlvl

          if(iter==1)then

             npx=npxt
             npy=npyt
             npsur=npsurt

          else

             if(npxt/=npx)then
                write(*,*)'Error in reading meteo data, npx/=npxt '
                stop
             endif
             if(npyt/=npy)then
                write(*,*)'Error in reading meteo data, npx/=npxt '
                stop
             endif

          endif
          !
          !     Store velocity in X
          !
          if(ldata%field == 'UU')then
             !
             !     Increment nblay
             ! 
             nblay=nblay+1

             npoin=np(1)
             npnew=npoin+npsur    
             if(.not.associated(uvel))then
                allocate(uvel(npnew),stat=istat)
                call memchk(zero,istat,memor_meteo,'UVEL','gtmtodata',uvel)
             else
                call memrea(npnew,memor_meteo,'UVEL','gtmtodata',uvel)
             endif
             !
             !     Store the slab
             !
             do iy=1,npy
                do ix=1,npx
                   npoin=npoin+1
                   uvel(npoin)=ldata%slab(ix,iy)
                enddo
             enddo
             np(1)=npnew       

             !
             !     Store pressure
             !
             npoin=np(5)
             npnew=npoin+npsur
             if(.not.associated(pres))then
                allocate(pres(npnew),stat=istat)
                call memchk(zero,istat,memor_meteo,'PRES','gtmtodata',pres)
             else
                call memrea(npnew,memor_meteo,'PRES','gtmtodata',pres)
             endif

             do iy=1,npy
                do ix=1,npx
                   npoin=npoin+1
                   pres(npoin)=press
                enddo
             enddo
             np(5)=npnew


          else if(ldata%field == 'VV')then
             !
             !     Store velocity in Y
             !
             npoin=np(2)
             npnew=npoin+npsur
             if(.not.associated(vvel))then
                allocate(vvel(npnew),stat=istat)
                call memchk(zero,istat,memor_meteo,'VVEL','gtmtodata',vvel)
             else
                call memrea(npnew,memor_meteo,'VVEL','gtmtodata',vvel)
             endif
             !
             !     Store the slab
             !

             do iy=1,npy
                do ix=1,npx
                   npoin=npoin+1
                   vvel(npoin)=ldata%slab(ix,iy)
                enddo
             enddo
             np(2)=npnew

          else if(ldata%field == 'HGT' .or. ldata%field == 'GHT')then
             !
             !     metgrid change HGT with GHT  ????????
             !


             !
             !     Store height
             !
             npoin=np(3)
             npnew=npoin+npsur
             if(.not.associated(height))then
                allocate(height(npnew),stat=istat)
                call memchk(zero,istat,memor_meteo,'HEIGHT','gtmtodata',height)
             else
                call memrea(npnew,memor_meteo,'HEIGHT','gtmtodata',height)
             endif
             !
             !     Store the slab
             !
             do iy=1,npy
                do ix=1,npx
                   npoin=npoin+1
                   height(npoin)=ldata%slab(ix,iy)
                enddo
             enddo
             np(3)=npnew

          else if(ldata%field == 'TT')then
             !
             !     Store temperature
             !
             npoin=np(4)
             npnew=npoin+npsur
             if(.not.associated(temp))then
                allocate(temp(npnew),stat=istat)
                call memchk(zero,istat,memor_meteo,'TEMP','gtmtodata',temp)
             else
                call memrea(npnew,memor_meteo,'TEMP','gtmtodata',temp)
             endif
             !
             !     Store the slab
             !
             do iy=1,npy
                do ix=1,npx
                   npoin=npoin+1
                   temp(npoin)=ldata%slab(ix,iy)
                enddo
             enddo
             np(4)=npnew

          else if(ldata%field == 'SOILHGT')then
             !
             !     Store orography
             !
             if(ldata%iproj/=0)then
                write(*,*)'Alya not yet ready for this projection'
                stop
             endif
             !
             !     Store only first time
             !
             if(npcoor==0)then

                phi=deg2rad * (ldata % startlat+90.0)-pi
                theta=deg2rad * ldata % startlon
                dphi=deg2rad * ldata % deltalat
                dtheta=deg2rad * ldata % deltalon 
                !write(*,*)ldata%startlat,ldata%startlon
                !write(*,*)ldata%deltalat,ldata%deltalon
                !write(*,*)ldata%startloc
                !
                !     Has the map been inverted?
                !
                if(ldata%deltalat<c00)invert=1_ip

                npcoor=npsur
                npoin=0_ip

                phi0=phi
                theta0=theta
                dphi0=dphi
                dtheta0=dtheta

                if(.not.associated(coor))then
                   allocate(coor(ndim,npnew),stat=istat)
                   call memchk(zero,istat,memor_meteo,'COOR','gtmtodata',coor)
                else
                   call memrea(npnew,memor_meteo,'COOR','gtmtodata',coor)
                endif


                do iy=1,npy

                   theta=theta0
                   !
                   !     Store the slab
                   !
                   do ix=1,npx
                      npoin=npoin+1
                      !radioloc=radio
                      radioloc=radio+ldata%slab(ix,iy)
                      coor(1,npoin)=cos(theta)*sin(phi)*radioloc
                      coor(2,npoin)=sin(theta)*sin(phi)*radioloc
                      coor(3,npoin)=cos(phi)*radioloc
                      theta=theta+dtheta
                   enddo

                   phi=phi+dphi

                enddo


             else

                write(*,*)'Error reading SOILHGT TWICE'
                stop

             endif

          else if(ldata%field == 'PSFC')then
             !
             !     Store surface pressure
             !
             npoin=np(6)
             if(npoin/=0)then
                write(*,*)'Error reading PSFC'
                stop
             endif

             npnew=npoin+npsur
             if(.not.associated(pressurf))then
                allocate(pressurf(npnew),stat=istat)
                call memchk(zero,istat,memor_meteo,'PRESSURF','gtmtodata',pressurf)
             else
                call memrea(npnew,memor_meteo,'PRESSURF','gtmtodata',pressurf)
             endif
             !
             !     Store the slab
             !
             do iy=1,npy
                do ix=1,npx
                   npoin=npoin+1
                   pressurf(npoin)=ldata%slab(ix,iy)
                enddo
             enddo
             np(6)=npnew

          endif

          deallocate(ldata%slab)

       endif

    enddo

    write(*,*)np(1),np(2),np(3),np(4)
    !
    !     Get total point form x velocity
    !
    npoin=np(1)
    !
    !     Check y velocity
    !
    if(np(2)/=npoin)then
       write(*,*)'Error reading meteo data np(2)/=npoin'
       stop
    endif
    !
    !     Check height
    ! 
    if(np(3)/=(npoin-npsur))then
       write(*,*)'Error reading meteo data np(3)/=(npoin-npsur)'
       stop
    endif
    !
    !     Check temperature
    !
    if(np(4)/=npoin)then
       write(*,*)'Error reading meteo data np(4)/=npoin'
       stop
    endif
    !
    !     Check surface pressure
    !
    if(np(6)==0)then
       write(*,*)'Error reading meteo data, PSFC not read'
       stop
    endif
    !
    !     Reallocate coor and copy height
    ! 
    call memrea(npoin,memor_meteo,'COOR','gtmtodata',coor)

    npoin=npsur

    do iblay=2,nblay

       phi=phi0

       do iy=1,npy

          theta=theta0
          !
          !     Store the slab
          !
          do ix=1,npx
             npoin=npoin+1
             radioloc=radio+height(npoin-npsur)
             coor(1,npoin)=cos(theta)*sin(phi)*radioloc
             coor(2,npoin)=sin(theta)*sin(phi)*radioloc
             coor(3,npoin)=cos(phi)*radioloc
             theta=theta+dtheta
          enddo

          phi=phi+dphi

       enddo
    enddo
    !
    !     Correct pressure for the first layer
    !
    do ipoin=1,npsur
       pres(ipoin)=pressurf(ipoin)
    enddo
    !
    !     Allocate vertival velocity
    !             
    allocate(wvel(npoin),stat=istat)
    call memchk(zero,istat,memor_meteo,'WVEL','gtmtodata',wvel)
    !
    !     Rotate velocity
    !

    npoint=0_ip
    do ilay=1,nblay

       phi=phi0
       do iy=1,npy

          theta=theta0         

          do ix=1,npx
             npoint=npoint+1
             uv=uvel(npoint)
             vv=vvel(npoint)
             !
             !     Rotate the velocity
             ! 
             uvel(npoint)=-sin(theta)*uv-cos(theta)*cos(phi)*vv
             vvel(npoint)= cos(theta)*uv-sin(theta)*cos(phi)*vv
             wvel(npoint)=                          sin(phi)*vv

             theta=theta+dtheta

          enddo

          phi=phi+dphi

       enddo

    enddo
    !
    !     Allocate a fictitious rblay, it will be corrected later
    !
    allocate(rblay(nblay),stat=istat)
    call memchk(zero,istat,memor_meteo,'RBLAY','gtmtodata',rblay)

    rblay(1)=c00
    do iblay=2,nblay
       rblay(iblay)=height(npsur*(iblay-1_ip)+1_ip)
    enddo
    nblay=nblay-1_ip
    !
    !     Compute the delta
    !
    do iblay=1,nblay
       rblay(iblay)=rblay(iblay+1)-rblay(iblay)
    enddo

    close(unit=input_unit) 
    !
    !     Deallocate pressurf
    !  
    call memchk(2_ip,istat,memor_meteo,'PRESSURF','gtmtodata',pressurf)
    deallocate(pressurf,stat=istat)
    if(istat/=0) call memerr(2_ip,'PRESSURF','gtmtodata',0_ip)

  end subroutine gtmtodata

  subroutine corsurf(npx,npy,npsurf,npoin,ndim,coor,uvel,vvel,temp,phi0,theta0,dphi0,dtheta0)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo, only        :  EARTH_RADIUS_M  
    implicit none
    integer(ip),intent(in) :: npx,npy,npsurf,npoin,ndim
    real(rp),intent(inout) :: coor(ndim,npoin),uvel(npoin),vvel(npoin),temp(npoin)
    real(rp),intent(in)    :: phi0,theta0,dphi0,dtheta0
    integer(ip)            :: npoint,iy,ix,ipsurf,ipvol
    real(rp)               :: uvsurf,vvsurf,tempsurf,c10,c20,phi,theta,dphi,dtheta,cten 
    real(rp)               :: uvnew,vvnew,tempnew,uvvol,vvvol,tempvol,heightsurf,heightvol
    real(rp)               :: dheight1,dheight2,ratio1,ratio2,radio 
    real(rp)               :: newtemp,newuv,newvv,rx,ry,rz,rnl 
    !
    !     Earth radius 
    !
    radio=EARTH_RADIUS_M
    !
    !     This subroutine corrects the velocity and temperature values at the surface
    !     as the temperature is given at 2m and the velocity at 10m 
    !
    c20=2.0d+00
    c10=1.0d+00 
    cten=1.0d+01

    npoint=0_ip

    phi=phi0
    dphi=dphi0
    dtheta=dtheta0
    !
    !     Loop on surface
    !
    do iy=1,npy

       theta=theta0

       do ix=1,npx
          npoint=npoint+1
          ipsurf=npoint
          !
          !     Get velocities and temperature at the surface 
          !
          uvsurf=uvel(ipsurf)
          vvsurf=vvel(ipsurf)
          tempsurf=temp(ipsurf) 
          !
          !     Get height of the point at the surface     
          !
          rx=coor(1,ipsurf)
          ry=coor(2,ipsurf)
          rz=coor(3,ipsurf)
          rnl=sqrt(rx*rx+ry*ry+rz*rz)
          heightsurf=rnl-radio 
          !
          !     Get the point of the next layer
          !
          ipvol=ipsurf+npsurf
          !
          !     Get velocities and temperature at the next layer 
          !
          uvvol=uvel(ipvol)
          vvvol=vvel(ipvol)
          tempvol=temp(ipvol) 
          !
          !     Get height of the point at the next layer     
          !
          rx=coor(1,ipvol)
          ry=coor(2,ipvol)
          rz=coor(3,ipvol)
          rnl=sqrt(rx*rx+ry*ry+rz*rz)
          heightvol=rnl-radio 
          !
          !     Compute difference in height for temperature
          !
          dheight1=heightsurf-c20
          dheight2=heightvol-c20 
          ratio1=dheight1/dheight2
          ratio2=c10-ratio1
          !
          !     Extrapolate (and limit) temperature
          !       
          newtemp=min(tempsurf*ratio2+tempvol*ratio1,tempsurf)
          !
          !     Compute difference in height for temperature
          !
          dheight1=heightsurf-cten
          dheight2=heightvol-cten
          ratio1=dheight1/dheight2
          ratio2=c10-ratio1
          !
          !     Extrapolate (and limit) velocities
          !       
          newuv=uvsurf*ratio2+uvvol*ratio1 
          newvv=vvsurf*ratio2+vvvol*ratio1 
          !
          !     Store velocities and temperature at the surface 
          !
          uvel(ipsurf)=newuv
          vvel(ipsurf)=newvv
          temp(ipsurf)=newtemp

          theta=theta+dtheta

       enddo

       phi=phi+dphi

    enddo

  end subroutine corsurf

  subroutine read_next_met_field(ldata,istatus,filename,input_unit)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo
    implicit none

    ! Arguments
    type (meteo), intent(inout) :: ldata
    integer(ip), intent(out) :: istatus
    integer(ip),intent(in)   :: input_unit
    character (len=256),intent(in) :: filename


    istatus = 0


    istatus = 1

    !  1) READ FORMAT VERSION
    read(unit=input_unit,err=1001,end=1001) ldata % version

    ! PREGRID
    if (ldata % version == 3) then

       read(unit=input_unit) ldata % hdate, &
            ldata % xfcst, &
            ldata % field, &
            ldata % units, &
            ldata % desc,  &
            ldata % xlvl,  &
            ldata % nx,    &
            ldata % ny,    &
            ldata % iproj

       ldata % map_source = ' '

       if (ldata % field == 'HGT      ') ldata % field = 'GHT      '

       ldata % starti = 1.0
       ldata % startj = 1.0

       ! Cylindrical equidistant
       if (ldata % iproj == 0) then
          ldata % iproj = PROJ_LATLON
          read(unit=input_unit,err=1001,end=1001) ldata % startlat, &
               ldata % startlon, &
               ldata % deltalat, &
               ldata % deltalon

          ! Mercator
       else if (ldata % iproj == 1) then
          ldata % iproj = PROJ_MERC
          read(unit=input_unit,err=1001,end=1001) ldata % startlat, &
               ldata % startlon, &
               ldata % dx,       &
               ldata % dy,       &
               ldata % truelat1

          ! Lambert conformal
       else if (ldata % iproj == 3) then
          ldata % iproj = PROJ_LC
          read(unit=input_unit,err=1001,end=1001) ldata % startlat, &
               ldata % startlon, &
               ldata % dx,       &
               ldata % dy,       &
               ldata % xlonc,    &
               ldata % truelat1, &
               ldata % truelat2

          ! Polar stereographic
       else if (ldata % iproj == 5) then
          ldata % iproj = PROJ_PS
          read(unit=input_unit,err=1001,end=1001) ldata % startlat, &
               ldata % startlon, &
               ldata % dx,       &
               ldata % dy,       &
               ldata % xlonc,    &
               ldata % truelat1

          ! ?????????
       else
          write(*,*)'Unrecognized projection code '
          stop     
       end if

       !ldata % earth_radius = EARTH_RADIUS_M / 1000.

       !#if (defined _GEOGRID) || (defined _METGRID)
       !         ldata % dx = ldata % dx * 1000.
       !         ldata % dy = ldata % dy * 1000.

       !         if (ldata % xlonc    > 180.) ldata % xlonc    = ldata%xlonc    - 360.

       !         if (ldata % startlon > 180.) ldata % startlon = ldata%startlon - 360.

       !         if (ldata % startlat < -90.) ldata % startlat = -90.
       !         if (ldata % startlat >  90.) ldata % startlat = 90.
       !#endif

       ldata % is_wind_grid_rel = .true.

       allocate(ldata % slab(ldata % nx, ldata % ny))
       read(unit=input_unit,err=1001,end=1001) ldata % slab

       istatus = 0 

       ! GRIB_PREP
    else if (ldata % version == 4) then

       read(unit=input_unit) ldata % hdate,      &
            ldata % xfcst,      &
            ldata % map_source, &
            ldata % field,      &
            ldata % units,      &
            ldata % desc,       &
            ldata % xlvl,       &
            ldata % nx,         &
            ldata % ny,         &
            ldata % iproj

       if (ldata % field == 'HGT      ') ldata % field = 'GHT      '

       ! Cylindrical equidistant
       if (ldata % iproj == 0) then
          ldata % iproj = PROJ_LATLON
          read(unit=input_unit,err=1001,end=1001) ldata%startloc, &
               ldata % startlat, &
               ldata % startlon, &
               ldata % deltalat, &
               ldata % deltalon

          ! Mercator
       else if (ldata % iproj == 1) then
          ldata % iproj = PROJ_MERC
          read(unit=input_unit,err=1001,end=1001) ldata%startloc, &
               ldata % startlat, &
               ldata % startlon, &
               ldata % dx,       &
               ldata % dy,       &
               ldata % truelat1

          ! Lambert conformal
       else if (ldata % iproj == 3) then
          ldata % iproj = PROJ_LC
          read(unit=input_unit,err=1001,end=1001) ldata%startloc, &
               ldata % startlat, &
               ldata % startlon, &
               ldata % dx,       &
               ldata % dy,       &
               ldata % xlonc,    &
               ldata % truelat1, &
               ldata % truelat2

          ! Polar stereographic
       else if (ldata % iproj == 5) then
          ldata % iproj = PROJ_PS
          read(unit=input_unit,err=1001,end=1001) ldata%startloc, &
               ldata % startlat, &
               ldata % startlon, &
               ldata % dx,       &
               ldata % dy,       &
               ldata % xlonc,    &
               ldata % truelat1

          ! ?????????
       else
          write(*,*)'Unrecognized projection code'
          stop

       end if

       if (ldata%startloc == 'CENTER  ') then
          ldata % starti = real(ldata % nx)/2.
          ldata % startj = real(ldata % ny)/2.
       else if (ldata%startloc == 'SWCORNER') then
          ldata % starti = 1.0
          ldata % startj = 1.0
       end if

       !ldata % earth_radius = EARTH_RADIUS_M / 1000.

       !#if (defined _GEOGRID) || (defined _METGRID)
       !         ldata % dx = ldata % dx * 1000.
       !         ldata % dy = ldata % dy * 1000.

       !         if (ldata % xlonc    > 180.) ldata % xlonc    = ldata % xlonc    - 360.

       !         if (ldata % startlon > 180.) ldata % startlon = ldata % startlon - 360.

       !         if (ldata % startlat < -90.) ldata % startlat = -90.
       !         if (ldata % startlat >  90.) ldata % startlat = 90.
       !#endif

       ldata % is_wind_grid_rel = .true.

       allocate(ldata % slab(ldata % nx, ldata % ny))
       read(unit=input_unit,err=1001,end=1001) ldata % slab

       istatus = 0

       ! WPS
    else if (ldata % version == 5) then

       read(unit=input_unit) ldata % hdate,      &
            ldata % xfcst,      &
            ldata % map_source, &
            ldata % field,      &
            ldata % units,      &
            ldata % desc,       &
            ldata % xlvl,       &
            ldata % nx,         &
            ldata % ny,         &
            ldata % iproj

       write(*,*)ldata % hdate      
       write(*,*)ldata % xfcst      
       write(*,*)ldata % map_source 
       write(*,*)ldata % field      
       write(*,*)ldata % units      
       write(*,*)ldata % desc       
       write(*,*)ldata % xlvl       
       write(*,*)ldata % nx         
       write(*,*)ldata % ny         
       write(*,*)ldata % iproj


       if (ldata % field == 'HGT      ') ldata % field = 'GHT      '

       ! Cylindrical equidistant
       if (ldata % iproj == 0) then
          ldata % iproj = PROJ_LATLON
          read(unit=input_unit,err=1001,end=1001) ldata%startloc, &
               ldata % startlat, &
               ldata % startlon, &
               ldata % deltalat, &
               ldata % deltalon, &
               ldata % earth_radius

          ! Mercator
       else if (ldata % iproj == 1) then
          ldata % iproj = PROJ_MERC
          read(unit=input_unit,err=1001,end=1001) ldata%startloc, &
               ldata % startlat, &
               ldata % startlon, &
               ldata % dx,       &
               ldata % dy,       &
               ldata % truelat1, &
               ldata % earth_radius

          ! Lambert conformal
       else if (ldata % iproj == 3) then
          ldata % iproj = PROJ_LC
          read(unit=input_unit,err=1001,end=1001) ldata%startloc, &
               ldata % startlat, &
               ldata % startlon, &
               ldata % dx,       &
               ldata % dy,       &
               ldata % xlonc,    &
               ldata % truelat1, &
               ldata % truelat2, &
               ldata % earth_radius

          ! Gaussian
       else if (ldata % iproj == 4) then
          ldata % iproj = PROJ_GAUSS
          read(unit=input_unit,err=1001,end=1001) ldata%startloc, &
               ldata % startlat, &
               ldata % startlon, &
               ldata % deltalat, &
               ldata % deltalon, &
               ldata % earth_radius

          ! Polar stereographic
       else if (ldata % iproj == 5) then
          ldata % iproj = PROJ_PS
          read(unit=input_unit,err=1001,end=1001) ldata%startloc, &
               ldata % startlat, &
               ldata % startlon, &
               ldata % dx,       &
               ldata % dy,       &
               ldata % xlonc,    &
               ldata % truelat1, &
               ldata % earth_radius

          ! ?????????
       else
          write(*,*)'Unrecognized projection code'
          stop

       end if

       if (ldata%startloc == 'CENTER  ') then
          ldata % starti = real(ldata % nx)/2.
          ldata % startj = real(ldata % ny)/2.
       else if (ldata%startloc == 'SWCORNER') then
          ldata % starti = 1.0
          ldata % startj = 1.0
       end if

       !#if (defined _GEOGRID) || (defined _METGRID)
       !         ldata % dx = ldata % dx * 1000.
       !         ldata % dy = ldata % dy * 1000.

       !         if (ldata % xlonc    > 180.) ldata % xlonc    = ldata % xlonc    - 360.

       !         if (ldata % startlon > 180.) ldata % startlon = ldata % startlon - 360.

       !         if (ldata % startlat < -90.) ldata % startlat = -90.
       !         if (ldata % startlat >  90.) ldata % startlat =  90.
       !#endif

       read(unit=input_unit,err=1001,end=1001) ldata % is_wind_grid_rel

       allocate(ldata % slab(ldata % nx, ldata % ny))
       read(unit=input_unit,err=1001,end=1001) ldata % slab

       istatus = 0

    else
       write(*,*)'Didn''t recognize format version of data'
       stop
    end if

    return

1001 return


  end subroutine read_next_met_field

  subroutine modmsh(coor,ndim,npoin,height,nblay,npsur)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo
    implicit none
    integer(ip),intent(in)    :: npoin,ndim,nblay,npsur 
    real(rp),intent(in)       :: height(npoin)
    real(rp),intent(inout)    :: coor(ndim,npoin)
    integer(ip)               :: ipoin,iblay,ipnew,ipheight
    real(rp)                  :: rnl,c10,rnx,rny,rnz,radio,radioloc

    c10=1.0d+00
    !
    !     Correct coordinates of the points
    ! 

    radio=EARTH_RADIUS_M

    do ipoin=1,npsur
       rnx=coor(1,ipoin)
       rny=coor(2,ipoin)
       rnz=coor(3,ipoin)
       rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
       rnl=c10/rnl
       rnx=rnx*rnl
       rny=rny*rnl
       rnz=rnz*rnl
       !
       !     Get the height of the points
       !     The height array has npsur points less than npoin 
       ! 
       do iblay=1,nblay
          ipnew=ipoin+iblay*npsur
          ipheight=ipnew-npsur 
          radioloc=radio+height(ipheight)
          coor(1,ipnew)=rnx*radioloc
          coor(2,ipnew)=rny*radioloc
          coor(3,ipnew)=rnz*radioloc
       enddo

    enddo

  end subroutine modmsh

  subroutine filmet(lface,nface,nnofa,coold,ndim,npold,npsur,var,nvar,height,coor,npoin,lmark,npx,npy,phi0,theta0,dphi,dtheta)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo
    use mod_cart
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,nvar,npoin
    integer(ip),intent(inout) :: npold,nface,npsur,npx,npy
    integer(ip),intent(inout) :: lface(nnofa,nface),lmark(npold) 
    real(rp),intent(inout)    :: coold(ndim,npold),var(nvar,npold),height(npold),phi0,theta0
    real(rp),intent(in)       :: coor(ndim,npoin),dphi,dtheta
    real(rp)                  :: bbox(ndim,2),xmin,xmax,ymin,ymax,zmin,zmax  
    integer(ip)               :: p1,p2,p3,iface,nface0,ivar,ipoin,npold0,iplace,npsur0     
    integer(ip)               :: isx,isy,isxmin,isymin,isymax,isxmax,jface     
    integer(ip)               :: ipa,ipb,ipc,ip1,ip2,ip3,ichk,ncont,jchk     
    integer(4)                :: istat
    !
    !     Compute bbox of the points of the surface triangulation
    !
    call boxbin(coor,ndim,npoin,bbox)
    !
    !     Filter the triangles 
    !
    nface0=nface
    nface=0_ip

    do iface=1,nface0

       p1=lface(1,iface)
       p2=lface(2,iface)
       p3=lface(3,iface)
       !
       !    Initialize min,max
       !
       xmin=coold(1,p1)
       ymin=coold(2,p1)
       zmin=coold(3,p1)  
       xmax=xmin
       ymax=ymin
       zmax=zmin
       !
       !     Compare with p2
       !
       if(coold(1,p2)<xmin)xmin=coold(1,p2)
       if(coold(1,p2)>xmax)xmax=coold(1,p2)

       if(coold(2,p2)<ymin)ymin=coold(2,p2)
       if(coold(2,p2)>ymax)ymax=coold(2,p2)

       if(coold(3,p2)<zmin)zmin=coold(3,p2)
       if(coold(3,p2)>zmax)zmax=coold(3,p2)
       !
       !     Compare with p3
       !
       if(coold(1,p3)<xmin)xmin=coold(1,p3)
       if(coold(1,p3)>xmax)xmax=coold(1,p3)

       if(coold(2,p3)<ymin)ymin=coold(2,p3)
       if(coold(2,p3)>ymax)ymax=coold(2,p3)

       if(coold(3,p3)<zmin)zmin=coold(3,p3)
       if(coold(3,p3)>zmax)zmax=coold(3,p3)

       if(xmin>bbox(1,2))cycle
       if(ymin>bbox(2,2))cycle
       if(zmin>bbox(3,2))cycle
       if(xmax<bbox(1,1))cycle
       if(ymax<bbox(2,1))cycle
       if(zmax<bbox(3,1))cycle
       !
       !     Mark the points to be kept
       ! 
       lmark(p1)=1_ip
       lmark(p2)=1_ip
       lmark(p3)=1_ip

    enddo
    !
    !     Take a bin pattern, find the max and the min in both direction 
    !  
    isxmin=npx
    isymin=npy
    isxmax=0_ip
    isymax=0_ip

    do isy=1,npy-1
       do isx=1,npx-1
          iface=2*(isx-1_ip)+2*(isy-1)*(npx-1)+1_ip
          jface=iface+1_ip          

          ip1=lface(1,iface)
          ip2=lface(2,iface)
          ip3=lface(3,iface)
          ipa=lface(1,jface)
          ipb=lface(2,jface)
          ipc=lface(3,jface)
          ichk=0_ip 


          if(lmark(ip1)+lmark(ip2)+lmark(ip3)==3)then
             lmark(ipa)=1_ip
             lmark(ipb)=1_ip
             lmark(ipc)=1_ip
             ichk=1_ip
          endif

          if(lmark(ipa)+lmark(ipb)+lmark(ipc)==3)then
             lmark(ip1)=1_ip
             lmark(ip2)=1_ip
             lmark(ip3)=1_ip
             ichk=1_ip
          endif

          if(ichk==1)then
             if(isx<isxmin)isxmin=isx
             if(isx>isxmax)isxmax=isx
             if(isy<isymin)isymin=isy
             if(isy>isymax)isymax=isy
          endif
       enddo
    enddo
    !
    !     Remark all the faces in the rank
    ! 
    do isy=isymin,isymax
       do isx=isxmin,isxmax
          iface=2*(isx-1_ip)+2*(isy-1)*(npx-1)+1_ip
          jface=iface+1_ip          

          ip1=lface(1,iface)
          ip2=lface(2,iface)
          ip3=lface(3,iface)
          ipa=lface(1,jface)
          ipb=lface(2,jface)
          ipc=lface(3,jface)
          lmark(ipa)=1_ip
          lmark(ipb)=1_ip
          lmark(ipc)=1_ip
          lmark(ip1)=1_ip
          lmark(ip2)=1_ip
          lmark(ip3)=1_ip
       enddo
    enddo
    !
    !     DBG
    !
    ncont=0_ip  
    do isy=1,npy-1
       do isx=1,npx-1
          iface=2*(isx-1_ip)+2*(isy-1)*(npx-1)+1_ip
          jface=iface+1_ip          

          ip1=lface(1,iface)
          ip2=lface(2,iface)
          ip3=lface(3,iface)
          ipa=lface(1,jface)
          ipb=lface(2,jface)
          ipc=lface(3,jface)
          ichk=0_ip
          jchk=0_ip


          if(lmark(ip1)+lmark(ip2)+lmark(ip3)==3)then
             ichk=1_ip
          endif
          if(lmark(ipa)+lmark(ipb)+lmark(ipc)==3)then
             jchk=1_ip
          endif

          if(ichk*jchk==0 .and. ichk+jchk>0)then
             write(*,*)'Error in filmet'
             stop
          endif

          if(ichk==1)then
             ncont=ncont+1_ip
          endif
       enddo
    enddo
    !
    !     Modify phi0 and theta0
    !
    phi0=phi0+(isymin-1)*dphi
    theta0=theta0+(isxmin-1)*dtheta
    !
    !     Remember npx and npy
    !
    npx=isxmax-isxmin+2_ip
    npy=isymax-isymin+2_ip

    do iface=1,nface0
       ip1=lface(1,iface)
       ip2=lface(2,iface)
       ip3=lface(3,iface)

       if(lmark(ip1)+lmark(ip2)+lmark(ip3)==3)then
          nface=nface+1_ip
          lface(1,nface)=ip1
          lface(2,nface)=ip2
          lface(3,nface)=ip3
       endif

    enddo
    !
    !     Get the renumbering for the surface points
    !
    npsur0=npsur
    npsur=0_ip
    do ipoin=1,npsur0
       if(lmark(ipoin)==1)then
          npsur=npsur+1_ip
          lmark(ipoin)=npsur 
          coold(1,npsur)=coold(1,ipoin)
          coold(2,npsur)=coold(2,ipoin)
          coold(3,npsur)=coold(3,ipoin)
       endif
    enddo
    !
    !     Renumber the faces
    ! 
    do iface=1,nface
       p1=lface(1,iface) 
       p2=lface(2,iface) 
       p3=lface(3,iface) 
       p1=lmark(p1)
       p2=lmark(p2)
       p3=lmark(p3)
       lface(1,iface)=p1 
       lface(2,iface)=p2 
       lface(3,iface)=p3 
    enddo
    !
    !     Compact the points
    !
    npold0=npold
    npold=0_ip

    do ipoin=1,npold0 
       iplace=mod(ipoin-1_ip,npsur0)+1_ip

       if(lmark(iplace)>0)then
          npold=npold+1
          height(npold)=height(ipoin)
          do ivar=1,nvar          
             var(ivar,npold)=var(ivar,ipoin) 
          enddo
       endif
    enddo

  end subroutine filmet

  subroutine filmet2(coold,ndim,npold,var,nvar,height,lmark,npsur)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo
    use mod_cart
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)    :: ndim,nvar,npsur
    integer(ip),intent(inout) :: npold
    integer(ip),intent(inout) :: lmark(npold) 
    real(rp),intent(inout)    :: coold(ndim,npold),var(nvar,npold),height(npold)
    real(rp)                  :: bbox(ndim,2),xmin,xmax,ymin,ymax,zmin,zmax  
    integer(ip)               :: p1,p2,p3,iface,nface0,ivar,ipoin,npold0,iplace     
    integer(4)                :: istat
    !
    !     Compact the points
    !
    npold0=npold
    npold=0_ip

    do ipoin=1,npold0 
       iplace=mod(ipoin-1_ip,npsur)+1_ip

       if(lmark(iplace)>0)then
          npold=npold+1
          height(npold)=height(ipoin)
          do ivar=1,nvar          
             var(ivar,npold)=var(ivar,ipoin) 
          enddo
          coold(1,npold)=coold(1,ipoin)
          coold(2,npold)=coold(2,ipoin)
          coold(3,npold)=coold(3,ipoin)
       endif
    enddo

  end subroutine filmet2

  subroutine invfac(lface,nface,nnofa)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo
    use mod_cart
    use mod_memchk
    implicit none
    integer(ip),intent(in)    :: nnofa,nface
    integer(ip),intent(inout) :: lface(nnofa,nface)  
    integer(ip)               :: iface,ip2,ip3

    do iface=1,nface
       ip2=lface(2,iface) 
       ip3=lface(3,iface)
       lface(2,iface)=ip3
       lface(3,iface)=ip2
    enddo

  end subroutine invfac

  subroutine renustr(elem,nelem,nnode,npsur,nface0,npoin,nblay,coor,ndim,lface,nnofa,nface)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo
    use mod_cart
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnode,nelem,npsur,ndim  
    integer(ip),intent(in)     :: nface0,npoin,nblay,nnofa,nface  
    integer(ip),intent(inout)  :: elem(nnode,nelem),lface(nnofa,nface)
    real(rp),intent(inout)     :: coor(ndim,npoin)
    integer(ip)                :: ielem,iface,ip1,ip2,ip3
    integer(ip)                :: ita,itb,itc,iblay,ipoin,ipnew
    integer(ip),pointer        :: lrenu(:)   
    real(rp),pointer           :: coort(:,:)   
    integer(4)                  :: istat
    !
    !     Allocate lmark
    ! 
    allocate(lrenu(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','renustr',lrenu)

    do iface=1,nface0

       ielem=3*nblay*(iface-1_ip)+1_ip
       !
       !     Get the point of the upper part
       !
       ip1=elem(1,ielem)
       ip2=elem(2,ielem)
       ip3=elem(3,ielem)

       lrenu(ip1)=ip1 
       lrenu(ip2)=ip2 
       lrenu(ip3)=ip3 


       ita=elem(1,ielem+2)
       itb=elem(4,ielem+2)
       itc=elem(3,ielem+2)
       ip1=ip1+npsur
       ip2=ip2+npsur
       ip3=ip3+npsur

       do iblay=1,nblay
          lrenu(ita)=ip1
          lrenu(itb)=ip2
          lrenu(itc)=ip3
          ita=ita+1_ip 
          itb=itb+1_ip 
          itc=itc+1_ip 
          ip1=ip1+npsur
          ip2=ip2+npsur
          ip3=ip3+npsur
       enddo
    enddo

    !
    !     DBG
    !
    do ipoin=1,npoin
       if(lrenu(ipoin)==0)then
          write(*,*)'Error in renustr',ipoin
          stop
       endif
    enddo
    !
    !     Renumber elements
    !
    do ielem=1,nelem
       elem(1,ielem)=lrenu(elem(1,ielem))
       elem(2,ielem)=lrenu(elem(2,ielem))
       elem(3,ielem)=lrenu(elem(3,ielem))
       elem(4,ielem)=lrenu(elem(4,ielem))
    enddo
    !
    !     Renumber faces
    !
    do iface=1,nface
       lface(1,iface)=lrenu(lface(1,iface))
       lface(2,iface)=lrenu(lface(2,iface))
       lface(3,iface)=lrenu(lface(3,iface))
    enddo
    !
    !     Allocate coort
    ! 
    allocate(coort(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'COORT','renustr',coort)

    do ipoin=1,npoin
       ipnew=lrenu(ipoin)
       coort(1,ipnew)=coor(1,ipoin)
       coort(2,ipnew)=coor(2,ipoin)
       coort(3,ipnew)=coor(3,ipoin)
    enddo

    do ipoin=1,npoin
       coor(1,ipoin)=coort(1,ipoin)
       coor(2,ipoin)=coort(2,ipoin)
       coor(3,ipoin)=coort(3,ipoin)
    enddo

    call memchk(2_ip,istat,memor_msh,'COORT','renustr',coort)
    deallocate(coort,stat=istat)
    if(istat/=0) call memerr(2_ip,'COORT','renustr',0_ip)
    call memchk(2_ip,istat,memor_msh,'LRENU','renustr',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','renustr',0_ip)


  end subroutine renustr

  subroutine gtbpoin(elem,nelem,nnode,npoin,nboup,lboup)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo
    use mod_cart
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnode,nelem,npoin  
    integer(ip),intent(inout)  :: nboup  
    integer(ip),intent(in)     :: elem(nnode,nelem)
    integer(ip),pointer        :: lboup(:)    
    integer(ip),pointer        :: ptoel1(:),ptoel2(:),eltoel(:,:),lmark(:)    
    integer(4)                 :: istat
    integer(ip)                :: ielem,inode,ineigh,ip1,ip2,ip3,ipoin
    integer(ip)                :: ltab(3,4) 

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
    !     Get the elements surrounding the points for the new mesh
    !
    call ptoelm(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !     Get the elements surrounding elements for the new mesh
    !
    call tetote(elem,nnode,nelem,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Allocate lmark
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','meteo',lmark)
    !
    !     Get the boundary points
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
          endif
       enddo
    enddo
    nboup=0_ip
    do ipoin=1,npoin
       if(lmark(ipoin)==1)then
          nboup=nboup+1_ip
       endif
    enddo
    allocate(lboup(nboup),stat=istat)
    call memchk(zero,istat,memor_msh,'LBOUP','meteo',lboup)
    nboup=0_ip
    do ipoin=1,npoin
       if(lmark(ipoin)==1)then
          nboup=nboup+1_ip
          lboup(nboup)=ipoin
       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','gtbpoin',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','gtbpoin',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','gtbpoin',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','gtbpoin',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','gtbpoin',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','gtbpoin',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','gtbpoin',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','gtbpoin',0_ip)

  end subroutine gtbpoin

  subroutine outgidmeteo(npoin,ndimn,coor,intmat,nelem,nnode,var,nvar)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)  :: nnode,nelem,ndimn,npoin,nvar    
    integer(ip),intent(in)  :: intmat(nnode,nelem)
    real(rp),intent(in)     ::  coor(ndimn,npoin),var(nvar,npoin)
    integer(ip)             :: ipoin,ip1,ip2,ip3,ip4,icont,iele
    real(rp)                :: rx,ry,rz,timer

    timer=0.0d+00


    open(unit=50,file='outgidmet.msh',status='unknown')
    rewind 50
1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(50,1)
    write(50,2)
    write(50,3)
    do ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(50,100)ipoin,rx,ry,rz
    end do
    write(50,4)
    write(50,5)
    icont=0
    do iele=1,nelem
       ip1=intmat(1,iele)
       ip2=intmat(2,iele)
       ip3=intmat(3,iele)
       ip4=intmat(4,iele)
       icont=icont+1
       write(50,200)icont,ip1,ip2,ip3,ip4
    end do
    write(50,6)
    close(50)

    open(unit=60,file='outgidmet.res',status='unknown')
    rewind 60

10  format('GID Post Results File 1.0')
11  format('Result "Velocity" "Analysis/time" ',1e20.10,' Vector OnNodes')
12  format('ComponentNames "Vx", "Vy" , "Vz" ')
13  format('Values')
14  format('End Values')
15  format('   ')
18  format('Result "Pressure" "Analysis/time" ',1e20.10,' Scalar OnNodes')
19  format('Result "Temperature" "Analysis/time" ',1e20.10,' Scalar OnNodes')
300 format(i10,1e20.10)


    write(60,10)

    write(60,11)timer
    write(60,12)
    write(60,13)
    do  ipoin=1,npoin
       write(60,100)ipoin,var(1,ipoin),var(2,ipoin),var(3,ipoin)
    enddo
    write(60,14)


    write(60,19)timer
    write(60,13)
    do  ipoin=1,npoin
       write(60,300)ipoin,var(4,ipoin)
    enddo
    write(60,14)


    write(60,18)timer
    write(60,13)
    do  ipoin=1,npoin
       write(60,300)ipoin,var(5,ipoin)
    enddo
    write(60,14)

    write(60,15)

    close(60)


  end subroutine outgidmeteo

  subroutine outgidmeteobc(npoin,ndimn,coor,var,nvar,lboup,nboup)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)  :: ndimn,npoin,nvar,nboup    
    integer(ip),intent(in)  :: lboup(nboup)
    real(rp),intent(in)     ::  coor(ndimn,npoin),var(nvar,npoin)
    integer(ip)             :: ipoin,ip1,ip2,ip3,ip4,icont,iele,iboup
    real(rp)                :: rx,ry,rz,timer

    timer=0.0d+00


    open(unit=50,file='outgidmetbc.msh',status='unknown')
    rewind 50
1   format('MESH dimension 3 ElemType Point Nnode 1')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(2i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(50,1)
    write(50,2)
    write(50,3)
    do iboup=1,nboup
       ipoin=lboup(iboup)
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(50,100)iboup,rx,ry,rz
    end do
    write(50,4)
    write(50,5)

    do iboup=1,nboup
       write(50,200)iboup,iboup
    end do


    write(50,6)
    close(50)

    open(unit=60,file='outgidmetbc.res',status='unknown')
    rewind 60

10  format('GID Post Results File 1.0')
11  format('Result "Velocity" "Analysis/time" ',1e20.10,' Vector OnNodes')
12  format('ComponentNames "Vx", "Vy" , "Vz" ')
13  format('Values')
14  format('End Values')
15  format('   ')
18  format('Result "Pressure" "Analysis/time" ',1e20.10,' Scalar OnNodes')
19  format('Result "Temperature" "Analysis/time" ',1e20.10,' Scalar OnNodes')
300 format(i10,1e20.10)


    write(60,10)

    write(60,11)timer
    write(60,12)
    write(60,13)
    do  iboup=1,nboup
       ipoin=lboup(iboup)
       write(60,100)iboup,var(1,ipoin),var(2,ipoin),var(3,ipoin)
    enddo
    write(60,14)


    write(60,19)timer
    write(60,13)
    do  iboup=1,nboup
       ipoin=lboup(iboup)
       write(60,300)iboup,var(4,ipoin)
    enddo
    write(60,14)


    write(60,18)timer
    write(60,13)
    do  iboup=1,nboup
       ipoin=lboup(iboup)
       write(60,300)iboup,var(5,ipoin)
    enddo
    write(60,14)

    write(60,15)

    close(60)


  end subroutine outgidmeteobc

  subroutine outmeteo(npoin,ndimn,coor,intmat,nelem,nnode,var,nvar,nx,ny,npsur&
       ,nblay,lface,nface,nnofa,lsurf)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)  :: nnode,nelem,ndimn,npoin,nvar,nx,ny,npsur,nblay,nnofa,nface    
    integer(ip),intent(in)  :: intmat(nnode,nelem),lface(nnofa,nface),lsurf(nface)
    real(rp),intent(in)     ::  coor(ndimn,npoin),var(nvar,npoin)
    integer(ip)             :: ipoin,ip1,ip2,ip3,ip4,icont,iele,ilay,cont,iface
    real(rp)                :: rx,ry,rz,timer
    real(rp)                :: u,v,w,rho,rhoe,rhou,rhov,rhow,pres,temp,c10,gam,R,c05,gam1

    timer=0.0d+00
    c10=1.0d+0    
    c05=0.5d+000

    open(unit=50,file='outmet.dat',status='unknown')
    rewind 50
100 format(i10,3e20.10)
200 format(6i10)
300 format(2i10)
    !
    !     Output point coordinates
    !
    write(50,*)'npoin',npoin
    do ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(50,100)ipoin,rx,ry,rz
    end do
    !
    !     Output elements
    !
    write(50,*)'nelem',nelem
    icont=0
    do iele=1,nelem
       ip1=intmat(1,iele)
       ip2=intmat(2,iele)
       ip3=intmat(3,iele)
       ip4=intmat(4,iele)
       icont=icont+1
       write(50,200)icont,ip1,ip2,ip3,ip4,0
    end do
    !
    !     Output boundary faces 
    !
    write(50,*)'nface',nface
    do iface=1,nface
       write(50,500)iface,lface(1,iface),lface(2,iface),lface(3,iface),lsurf(iface) 
    enddo
    !
    !     Output boundary points 
    !     First the points on the ground 
    !    
    write(50,*)'nboup',(2*npsur+2*nx+2*(ny-1))
    do ipoin=1,npsur 
       write(50,300)ipoin,1  
    enddo
    !
    !     Then sides
    !
    do ilay=2,nblay-1

       icont=npsur*(ilay-1) 

       do ipoin=1,nx
          write(50,300)icont+ipoin,2  
          write(50,300)icont+nx*(ny-1)+ipoin,2 
       enddo

       do ipoin=2,ny-1
          write(50,300)icont+nx*(ipoin-1)+1,2
          write(50,300)icont+nx*ipoin,2
       enddo
    enddo
    !
    !     Then ceiling
    !
    icont=(nblay-1)*npsur

    do ipoin=1,npsur
       write(50,300)icont+ipoin,3  
    enddo

    close(50)
    !
    !     Output results
    !
    R=287.0d+00
    gam=1.40027894
    gam1=gam-c10

    open(unit=60,file='outmet.init',status='unknown')
    rewind 60

    do ipoin=1,npoin
       u=var(1,ipoin) 
       v=var(2,ipoin) 
       w=var(3,ipoin)
       temp=var(4,ipoin)
       pres=var(5,ipoin)

       rho=pres/(R*temp)
       rhou=rho*u
       rhov=rho*v
       rhow=rho*w
       rhoe=pres/(gam1)+c05*(rho*(u*u+v*v+w*w))
       write(60,400)rho,rhou,rhov,rhow,rhoe 
    enddo

400 format(5e20.10)
500 format(5i10)

    close(60)


  end subroutine outmeteo

  subroutine outmetbc(npoin,var,nvar,nboup,lboup,iter)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)  :: npoin,nvar,nboup,iter    
    integer(ip),intent(in)  :: lboup(nboup) 
    real(rp),intent(in)     :: var(nvar,npoin)
    integer(ip)             :: ipoin,iboup
    character (len=256)     :: filename
    real(rp)                :: u,v,w,rho,rhoe,rhou,rhov,rhow,pres,temp,c10,gam,R,c05,gam1

    c10=1.0d+0    
    c05=0.5d+000


    write(filename,'(a,i5.5,a)')'outmet',iter,'.bc'

    open(unit=50,file=filename,status='unknown')
    rewind 50
    
    R=287.0d+00
    gam=1.40027894
    gam1=gam-c10

    do iboup=1,nboup
       ipoin=lboup(iboup)
       u=var(1,ipoin) 
       v=var(2,ipoin) 
       w=var(3,ipoin)
       temp=var(4,ipoin)
       pres=var(5,ipoin)

       rho=pres/(R*temp)
       rhou=rho*u
       rhov=rho*v
       rhow=rho*w
       rhoe=pres/(gam1)+c05*(rho*(u*u+v*v+w*w))
       write(50,100)ipoin,rho,rhou,rhov,rhow,rhoe
    enddo

100 format(i10,5e20.10)

    close(50)

  end subroutine outmetbc

  subroutine sph2pla(coor,npoin,phi0,theta0,dphi0,dtheta0,npx,npy,npsurf,ndim,nlay)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo, only        :  EARTH_RADIUS_M  
    implicit none
    integer(ip),intent(in) :: npx,npy,npsurf,npoin,ndim,nlay
    real(rp),intent(inout) :: coor(ndim,npoin)
    real(rp),intent(in)    :: phi0,theta0,dphi0,dtheta0
    integer(ip)            :: npoint,iy,ix,ipsurf,ipvol,ilay
    real(rp)               :: uvsurf,vvsurf,tempsurf,c10,c20,cten 
    real(rp)               :: uvnew,vvnew,tempnew,uvvol,vvvol,tempvol,heightsurf,heightvol
    real(rp)               :: dheight1,dheight2,ratio1,ratio2,radio,x0,y0,height,dx,dy 
    real(rp)               :: newtemp,newuv,newvv,rx,ry,rz,rnl 
    !
    !     Earth radius 
    !
    radio=EARTH_RADIUS_M
    !
    !     This subroutine puts a spherical mesh on a plane 
    !
    c20=2.0d+00
    c10=1.0d+00 
    cten=1.0d+01
    npoint=0_ip
    !
    !     Initial positions computed as l=r*theta as the earth is considered as spherical
    !
    x0=theta0*EARTH_RADIUS_M
    y0=phi0*EARTH_RADIUS_M
    dx=dtheta0*EARTH_RADIUS_M
    dy=dphi0*EARTH_RADIUS_M
    !
    !     Loop on surface
    !
    do iy=1,npy

       do ix=1,npx

          npoint=npoint+1
          ipsurf=npoint
          !
          !     Get height of the point at the surface     
          !
          rx=coor(1,ipsurf)
          ry=coor(2,ipsurf)
          rz=coor(3,ipsurf)
          rnl=sqrt(rx*rx+ry*ry+rz*rz)
          heightsurf=rnl-radio 
          !
          !     Change coordinates
          !
          coor(1,ipsurf)=dx*(ix-1_ip)+x0
          coor(2,ipsurf)=dy*(iy-1_ip)+y0
          coor(3,ipsurf)=heightsurf
          !
          !     Loop on height
          !
          do ilay=1,nlay
             !
             !     Get the point of the next layer
             !
             ipvol=ipsurf+(npsurf*ilay)
             !
             !     Get height of the point at the next layer     
             !
             rx=coor(1,ipvol)
             ry=coor(2,ipvol)
             rz=coor(3,ipvol)
             rnl=sqrt(rx*rx+ry*ry+rz*rz)
             heightvol=rnl-radio 
             !
             !     Change coordinates
             !
             coor(1,ipvol)=coor(1,ipsurf)
             coor(2,ipvol)=coor(2,ipsurf)
             coor(3,ipvol)=heightvol

          enddo

       enddo

    enddo

  end subroutine sph2pla

  subroutine pla2sph(coor,npoin,phi0,theta0,dphi0,dtheta0,npx,npy,npsurf,ndim,nlay)
    use def_kintyp, only       :  ip,rp,lg
    use def_meteo, only        :  EARTH_RADIUS_M  
    implicit none
    integer(ip),intent(in) :: npx,npy,npsurf,npoin,ndim,nlay
    real(rp),intent(inout) :: coor(ndim,npoin)
    real(rp),intent(in)    :: phi0,theta0,dphi0,dtheta0
    integer(ip)            :: npoint,iy,ix,ipsurf,ipvol,ilay
    real(rp)               :: uvsurf,vvsurf,tempsurf,c10,c20,phi,theta,dphi,dtheta,cten 
    real(rp)               :: uvnew,vvnew,tempnew,uvvol,vvvol,tempvol,heightsurf,heightvol
    real(rp)               :: dheight1,dheight2,ratio1,ratio2,radio,x0,y0,height,dx,dy 
    real(rp)               :: newtemp,newuv,newvv,radioloc 
    !
    !     Earth radius 
    !
    radio=EARTH_RADIUS_M
    !
    !     This subroutine put a spherical mesh on a plane 
    !
    c20=2.0d+00
    c10=1.0d+00 
    cten=1.0d+01

    npoint=0_ip

    phi=phi0
    dphi=dphi0
    dtheta=dtheta0

    !
    !     Loop on surface
    !
    do iy=1,npy

       theta=theta0

       do ix=1,npx
          npoint=npoint+1
          ipsurf=npoint
          !
          !     Get height of the point at the surface     
          !
          radioloc=radio+coor(3,ipsurf)
          !
          !     Change coordinates
          !
          coor(1,ipsurf)=cos(theta)*sin(phi)*radioloc
          coor(2,ipsurf)=sin(theta)*sin(phi)*radioloc
          coor(3,ipsurf)=cos(phi)*radioloc
          !
          !     Loop on height
          !
          do ilay=1,nlay
             !
             !     Get the point of the next layer
             !
             ipvol=ipsurf+(npsurf*ilay)
             radioloc=radio+coor(3,ipvol)
             coor(1,ipvol)=cos(theta)*sin(phi)*radioloc
             coor(2,ipvol)=sin(theta)*sin(phi)*radioloc
             coor(3,ipvol)=cos(phi)*radioloc

          enddo

          theta=theta+dtheta

       enddo

       phi=phi+dphi

    enddo

  end subroutine pla2sph

end module mod_meteo



!!$/* File: read_geogrid.c
!!$
!!$   Sample subroutine to read an array from the geogrid binary format.
!!$
!!$   Notes: Depending on the compiler and compiler flags, the name of 
!!$   the read_geogrid() routine may need to be adjusted with respect
!!$   to the number of trailing underscores when calling from Fortran.
!!$
!!$   Michael G. Duda, NCAR/MMM
!!$*/
!!$
!!$#include <stdlib.h>
!!$#include <stdio.h>
!!$#include <string.h>
!!$
!!$/* 
!!$#ifdef _UNDERSCORE
!!$#define read_geogrid read_geogrid_
!!$#endif
!!$#ifdef _DOUBLEUNDERSCORE
!!$#define read_geogrid read_geogrid__
!!$#endif
!!$
!!$#define BIG_ENDIAN    0
!!$#define LITTLE_ENDIAN 1
!!$*/
!!$
!!$int read_geogrid_(
!!$      char * fname,            /* The name of the file to read from */
!!$      int * len,               /* The length of the filename */
!!$      float * rarray,          /* The array to be filled */
!!$      int * nx,                /* x-dimension of the array */
!!$      int * ny,                /* y-dimension of the array */
!!$      int * nz,                /* z-dimension of the array */
!!$      int * isigned,           /* 0=unsigned data, 1=signed data */
!!$      int * endian,            /* 0=big endian, 1=little endian */
!!$      float * scalefactor,     /* value to multiply array elements by before truncation to integers */
!!$      int * wordsize,          /* number of bytes to use for each array element */
!!$      int * status)
!!${
!!$   int i, ival, cnt, narray;
!!$   int A2, B2;
!!$   int A3, B3, C3;
!!$   int A4, B4, C4, D4;
!!$   unsigned char * c;
!!$   char local_fname[1024];
!!$   FILE * bfile;
!!$
!!$   return 0;
!!$}
