!****************************************************************************
!
!  PROGRAM: Read_Chann
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program Read_Cav3D

  implicit none

  character(40):: fil_dummi, file, label, fil_meanv, fil_avden, fil_fvvel, fil_avtan, fil_avpre, fil_avpr2
  character(40):: fil_avvel, fil_avve2, fil_avvxy, fil_avtem, fil_avte2, fil_avtev, fil_coord
  integer       :: index, ifile, ipoin, inifile
  ! discretization 
  integer    :: nex, ney, nez, npx, npy, npz, npoin, ndime, ninte
  integer    :: i, j, iz, jz,kz,jj, jpoin, ip, jp, kp, iinte
  integer    :: ifile_g, ifile_p, ifile_b 
  integer    :: je, idime, NFILE
  integer    :: nfile_avvel, nfile_avve2, nfile_avvxy, nfile_avtem, nfile_avte2, nfile_avtev
  integer    :: nfile_avtan, nfile_avpre, nfile_avpr2
  integer    :: lun_coord
  integer    :: lun_y1plu, lun_homog, lun_outmn, lun_outwc, lun_outwh
  logical    :: iex, iex_avvel,iex_avve2,iex_avvxy,iex_avtem,iex_avte2,iex_avtev, kfl_reaho
  logical    :: iex_avden, iex_fvvel, kfl_incom, walla, iex_avtan, iex_avpre, iex_avpr2
  real       :: facto, zinte(10), Rayle, Prand, Rebul, nusse(2), cfmax(2), bqwal(2), ustar(2)
  real       :: rhbul, vebul, tebul, rhois, rhvis, temis,tstar(2), retau(2), qwall(2), tauwa(2)
  real       :: grate(2), grave(2), densi(2), cp, visco(2), condu(2), R, S, tehot, tecol, tauwa_av, ustar_av
  real       :: teref, delta, muref, dy, velis, maxve, rhref, Cfbul(2), wmult, walld(2)
  real, allocatable ::  &
       coord(:,:),      &   !Coordinates !ndime, npoin
       coorx(:),        &   !X -Coor
       coory(:),        &   !Y -Coor
       coorz(:),        &   !Z -Coor  
       !
       !  HOMOGENEIZED VARIABLES
       !
       hom_avvel(:,:),    &   !homogeneous V
       hom_avve2(:,:),    &   !homogeneous V**2
       hom_avverms(:,:),  &   !homogeneous rms V
       hom_avvxy(:,:),    &   !homogeneous Vx*Vy Vx*Vz Vy*Vz
       hom_avtem(:),      &   !homogeneous T
       hom_avte2(:),      &   !homogeneous T**2
       hom_avterms(:),    &   !homogeneous rms T
       hom_avtev(:,:),    &   !homogeneous V*T
       hom_avden(:),      &   !homogeneous densi      
       hom_fvvel(:,:),    &   !homogeneous favre velocity
       hom_avtan(:,:),    &   !homogeneous av tangent stress
       hom_avpre(:),      &   !homogeneous av pressure
       hom_avpr2(:),      &   !homogeneous av pressure**2
       hom_avprerms(:),   &   !homogeneous rms pressure
       yplus(:),          &   ! all nodes y plus 
       vplus(:),          &   ! all nodes v plus
       tplus(:)               ! all node t plus

  
  ! DATA , REFERENCE VISCOSITY
  muref = 1.0d0/180.0d0
  muref = 1.14d0*1.788d-5
  ! Reference density
  rhref = 1.14d0 
  ! WALL LAW
  walla = .true. ! exist
  wmult = 2.0d0  ! factor (multi)
  kfl_incom = .false.
  !
  !   NUMBER OF FILES (CONTAINING STEPS) THAT WILL BE MEANED
  !
  NFILE=15

  write(*,'(a,$)') ' Number of ensi files of results each to be read '
  read (*,*) NFILE
  index= 0
  do while (index/=1.and.index/=2)
     write(*,'(a,$)') 'Read all ensi files (1) or a previous homogeineized file (2) '
     read (*,*) index
  end do
  if (index==1) then
     kfl_reaho = .false.
  else if (index==2) then
     kfl_reaho = .true.
  else
     print *, 'ERROR IN INDEX'
     stop
  end if
  !
  !   INITIALIZATIONS
  !    
  ndime = 3
  iex=.false.
  index=0
  lun_outwc = 19
  lun_outwc = 21
  lun_coord=  23
  lun_outmn = 25
  lun_homog = 27


  ! existing files 

  iex_avvel = .false.
  iex_avve2 = .false.
  iex_avvxy = .false.
  iex_avtem = .false.
  iex_avte2 = .false.
  iex_avtev = .false.
  iex_avtan = .false.
  iex_avpre = .false.
  iex_avpr2 = .false.

  !
  !   LOOKS THE FILE NAME

  CALL GETARG(1, fil_meanv)  
  
  !
  !   outputs ON THE SCREEN the found input data name 
  !    
  write (*,*) 'READ DATA FOR PROBLEM : ', fil_meanv(1:8)
  !
  !   Loads problem discretization 
  !    
  read(fil_meanv(5:6),*)  nex
  if (nex.lt.20) read(fil_meanv(5:7),*)  nex
  ney = nex
  nez = nex
  npx = nex +1 ! number of nodes in x direction
  npy = ney +1
  npz = nez +1    
  print *, 'NEX,NEY,NEZ,= ', NEX,NEY,NEZ,NPX,NPY,NPZ
  
  

  file =trim(fil_meanv)//'.ensi'
  ! 

  fil_coord=  trim(file)//'.geo'
  fil_avvel=  trim(file)//'.AVVEL-000'
  fil_avve2=  trim(file)//'.AVVE2-000'
  fil_avvxy=  trim(file)//'.AVVXY-000'
  fil_avtem=  trim(file)//'.AVTEM-000'
  fil_avte2=  trim(file)//'.AVTE2-000'
  fil_avtev=  trim(file)//'.AVTEV-000'
  fil_avden=  trim(file)//'.AVDEN-000'
  fil_fvvel=  trim(file)//'.FVVEL-000'
  fil_avtan=  trim(file)//'.AVTAN-000'
  fil_avpre=  trim(file)//'.AVPRE-000'
  fil_avpr2=  trim(file)//'.AVPR2-000'

  !
  !    looks which files exists
  !

  !    Exists coord?

  inquire(file=fil_coord, EXIST=iex) 
  if (iex==.false.) then  
     print *, 'NOT FOUND coordinates file DATA',fil_coord
     stop
  end if
  
  !    Exists AVVEL ?
  fil_dummi = trim(fil_avvel)//'001'  
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avvel = .true.
  if (.not.iex) print *, fil_dummi, 'not found'

  !    Exists AVVE2 ?
  fil_dummi = trim(fil_avve2)//'001' 
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avve2 = .true.
  if (.not.iex) print *, fil_avve2, 'not found'

  !    Exists AVVXY ?
  fil_dummi = trim(fil_avvxy)//'001'  
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avvxy = .true.
  if (.not.iex) print *, fil_avvxy, 'not found'

  !    Exists AVTEM ?
  fil_dummi = trim(fil_avtem)//'001'  
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avtem = .true.
  if (.not.iex) print *, fil_avtem, 'not found'

  !    Exists AVTE2 ?
  fil_dummi = trim(fil_avte2)//'001'  
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avte2 = .true.
  if (.not.iex) print *, fil_avte2, 'not found'

  !    Exists AVTEV ?
  fil_dummi = trim(fil_avtev)//'001'  
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avtev = .true.
  if (.not.iex) print *, fil_avtev, 'not found'

  !    Exists AVDEN ?
  fil_dummi = trim(fil_avden)//'001'  
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avden = .true.
  if (.not.iex) print *, fil_avden, 'not found'

  !    Exists FVVEL ?y
  fil_dummi = trim(fil_fvvel)//'001'
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_fvvel = .true.
  if (.not.iex) print *, fil_fvvel, 'not found'
 
   !    Exists AVTAN ?
  fil_dummi = trim(fil_avtan)//'001'  
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avtan = .true.
  if (.not.iex) print *, fil_dummi, 'not found'

   !    Exists AVPRE ?
  fil_dummi = trim(fil_avpre)//'001'  
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avpre = .true.
  if (.not.iex) print *, fil_dummi, 'not found'
  
  !    Exists AVPR2 ?
  fil_dummi = trim(fil_avpr2)//'001'  
  inquire(file=fil_dummi, EXIST=iex) 
  if (iex)  iex_avpr2 = .true.
  if (.not.iex) print *, fil_dummi, 'not found'

  npoin = npx*npy*npz
 
  !******************************************************************************
  !  OPEN FILES
  !******************************************************************************  
  !
  !  READ COORDINATES FILE  name.ensi.geo
  !
  IF  (.not.kfl_reaho) THEN ! if reads from ensi files

     OPEN (UNIT= lun_coord, FILE= fil_coord, IOSTAT= JJ)  ! opens coordinates    
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*) npoin
     print *, 'npoin',npoin
     
     
     if (npoin /= npx*npy*npz) then
        print *, 'error calculating npoin'
        stop
     end if
     allocate(coord(ndime, npoin))
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,*)
     READ(lun_coord,'(i8)') jpoin
     do while (jpoin.ne.npoin)      
        READ(lun_coord,*) jpoin
     end do
     do idime =1, 3
        do ipoin =1, npoin
           read(lun_coord,*) coord(idime, ipoin)                                  
        end do
     end do
     close (lun_coord)      
     !
     ! save homogeneous sizes
     !
     allocate (coorx(npx))
     allocate (coory(npy))
     allocate (coorz(npz))
     !    ipoin = (i-1)*ny*nz + (j-1)*nz + k
     do ip =1, npx
        ipoin =  (ip-1)*npy*npz + 1
        coorx(ip) = coord(1,ipoin)
     end do
     do jp =1, npy
        ipoin =   jp*npz 
        coory(jp) = coord(2,ipoin)
     end do
     do kp =1, npz
        ipoin =  kp
        coorz(kp) = coord(3,ipoin)
     end do

     !
     ! Allocates structures
     !
     allocate (hom_avvel(ndime, npy))
     allocate (hom_avve2(ndime, npy))
     allocate (hom_avvxy(ndime, npy))
     allocate (hom_avtem(npy))
     allocate (hom_avte2(npy))
     allocate (hom_avtev(ndime, npy))
     allocate (hom_avden(npy))
     allocate (hom_fvvel(ndime, npy))
     allocate (hom_avtan(ndime, npy))
     allocate (hom_avpre(npy))
     allocate (hom_avpr2(npy))
     !
     !  READS FILES CONTAINING AVERGAED NODAL VARIABLES DURING some STEPS

     if (.not.iex_avtem) kfl_incom = .true.
     if (iex_avvel)   call readavfile(npoin, ndime,npx,npy,npz,fil_avvel,NFILE,nfile_avvel,hom_avvel)
     if (iex_avve2)   call readavfile(npoin, ndime,npx,npy,npz,fil_avve2,NFILE,nfile_avvel,hom_avve2)
     if (iex_avvxy)   call readavfile(npoin, ndime,npx,npy,npz,fil_avvxy,NFILE,nfile_avvel,hom_avvxy)
     if (iex_avtem)   call readavfile(npoin,     1,npx,npy,npz,fil_avtem,NFILE,nfile_avvel,hom_avtem)
     if (iex_avte2)   call readavfile(npoin,     1,npx,npy,npz,fil_avte2,NFILE,nfile_avte2,hom_avte2)
     if (iex_avtev)   call readavfile(npoin, ndime,npx,npy,npz,fil_avtev,NFILE,nfile_avtev,hom_avtev)
     if (iex_avden)   call readavfile(npoin,     1,npx,npy,npz,fil_avtev,NFILE,nfile_avtev,hom_avden)
     if (iex_fvvel)   call readavfile(npoin, ndime,npx,npy,npz,fil_avtev,NFILE,nfile_avtev,hom_fvvel)
     if (iex_avtan)   call readavfile(npoin, ndime,npx,npy,npz,fil_avtan,NFILE,nfile_avtan,hom_avtan)
     if (iex_avpre)   call readavfile(npoin,     1,npx,npy,npz,fil_avpre,NFILE,nfile_avpre,hom_avpre)
     if (iex_avpr2)   call readavfile(npoin,     1,npx,npy,npz,fil_avpr2,NFILE,nfile_avpr2,hom_avpr2)
!!$
!!$   Write a file with the resulting homogeneous results
!!$
     fil_dummi = trim(file)//'.HOMOG'
     OPEN (UNIT= lun_homog, FILE= fil_dummi, STATUS='replace')  
     write(lun_homog,'(a,i,a,i)') '# calculated homogeneous results from inifile ', inifile, ' to file ', nfile_avtev
     write(lun_homog,'(3i7)')   nex, ney, nez
     do jp =1, npy
        do kp =1, npz         
           write(lun_homog,'(2(f15.8))')   coory(jp)
           if (iex_avvel)  write(lun_homog,'(3(f15.8))')   hom_avvel(1:3, jp)
           if (iex_avve2)  write(lun_homog,'(3(f15.8))')   hom_avve2(1:3, jp)
           if (iex_avvxy)  write(lun_homog,'(3(f15.8))')   hom_avvxy(1:3, jp)
           if (iex_avtev)  write(lun_homog,'(3(f15.8))')   hom_avtev(1:3, jp)
           if (iex_avtem)  write(lun_homog,'(f15.8)')   hom_avtem(jp)
           if (iex_avte2)  write(lun_homog,'(f15.8)')   hom_avte2(jp)
           if (iex_avden)  write(lun_homog,'(f15.8)')   hom_avden(jp)
           if (iex_fvvel)  write(lun_homog,'(f15.8)')   hom_fvvel(1:3,jp)
           if (iex_avtan)  write(lun_homog,'(f15.8)')   hom_avtan(1:3,jp)
           if (iex_avpre)  write(lun_homog,'(f15.8)')   hom_avpre(jp)
           if (iex_avpr2)  write(lun_homog,'(f15.8)')   hom_avpr2(jp)
        end do
     end do
     close (lun_homog)

   !
   ! READS RESULTS FROM PREVIOUSLY HOMOGENEIZED FILE
   !
     
  ELSE IF(KFL_REAHO) then ! reads from previous homogeneized file
     
     fil_dummi = trim(file)//'.HOMOG'
     OPEN (UNIT= lun_homog, FILE= fil_dummi)  
     read(lun_homog,*)      
     read(lun_homog,*) nex, ney, nez
     npx = nex +1 ! number of nodes in x direction
     npy = ney +1
     npz = nez +1
     npoin = npx*npy*npz

     allocate (coorx(npx))
     allocate (coory(npy))
     allocate (coorz(npz))
     allocate (hom_avvel(ndime, npy))
     allocate (hom_avve2(ndime, npy))
     allocate (hom_avvxy(ndime, npy))
     if (.not.kfl_incom) then
        allocate (hom_avtev(ndime, npy))
        allocate (hom_avtem(npy))
        allocate (hom_avte2(npy))
        allocate (hom_avden(npy))
        allocate (hom_fvvel(ndime,npy))
     end if
     allocate (hom_avtan(ndime, npy))
     allocate (hom_avpre(npy))
     allocate (hom_avpr2(npy))

     do jp =1, npy
        do kp =1, npz
           read(lun_homog,*)   coory(jp)
           if (iex_avvel)   read(lun_homog,*)   hom_avvel(1:3, jp)
           if (iex_avve2)   read(lun_homog,*)   hom_avve2(1:3, jp)
           if (iex_avden)   read(lun_homog,*)   hom_avvxy(1:3, jp)
           if (iex_avvxy)   read(lun_homog,*)   hom_avtev(1:3, jp)
           if (iex_avtem)   read(lun_homog,*)   hom_avtem(jp)
           if (iex_avte2)   read(lun_homog,*)   hom_avte2(jp)
           if (iex_avden)   read(lun_homog,*)   hom_avden(jp)
           if (iex_fvvel)   read(lun_homog,*)   hom_fvvel(1:3, jp)
           if (iex_avtan)   read(lun_homog,*)   hom_avtan(1:3, jp)
           if (iex_avpre)   read(lun_homog,*)   hom_avpre(jp)
           if (iex_avpr2)   read(lun_homog,*)   hom_avpr2(jp)
        end do
     end do
     close (lun_homog)
  END IF
  
  !*******************************************************************************
  ! OBTAINS PHYSICAL PROPERTIES
  !*******************************************************************************
  delta = (coory(npy) - coory(1))*0.5d0
  if(.not.kfl_incom) then
     tehot = hom_avtem(npy)
     tecol = hom_avtem(1)

     teref = tecol

  
     Prand = 0.71d0
     cp = 1004.5d0
     R  = 287.0d0
     S  = 0.368d0 ! Sutherland's temper
     ! Density, viscosity and conductivity over the walls
     ! visco cold (Sutherland

     visco(1)= (tecol/teref)**1.5d0*(teref+S)/(tecol+S)*muref
     ! visco_hot
     visco(2)= (tehot/teref)**1.5d0*(teref+S)/(tehot+S)*muref

     condu(1) = visco(1)*cp/Prand
     condu(2) = visco(2)*cp/Prand
  else
     tehot =1.0d0
     tecol =1.0d0
     teref =1.0d0
     visco(1) = muref
     visco(2) = muref
  end if
  if (walla) then
     walld(1) = wmult*0.5d0*(coory(2)- coory(1))
     walld(2) = wmult*0.5d0*(coory(npy)- coory(npy-1))
  else
     walld(1) =0.0d0
     walld(2) =0.0d0
  end if

  ! approximated densities over walls 
  
!  densi(1) = hom_avden(1)
!  densi(2) = hom_avden(npy)
  if (tehot/tecol.gt.1.1) then
     prthe = 424.5d0
     densi(1) = prthe/(R*tecol)*rhref
     densi(2) = prthe/(R*tehot)*rhref
  else
     densi(1) = 1.0*rhref
     densi(2) = tecol/tehot *rhref
  end if
  
  print *, 'densi(1), densi(2)', densi(1), densi(2), kfl_incom

  allocate(yplus(npy))
  allocate(vplus(npy))
  allocate(tplus(npy))

!
!  WALL VALUES (1): Cold, (2): hot
!
  !  wall gradients
  grave(1) = (hom_avvel(1,2)-hom_avvel(1,1))/(coory(2)-coory(1))
  grave(2) = (hom_avvel(1,npy-1)-hom_avvel(1,npy))/(coory(npy)-coory(npy-1))

  ! temper gradient
  if (.not.kfl_incom) then
     grate(1) =  (hom_avtem(2)-tecol)/(coory(2)+1.0d0)
     grate(2) =  (hom_avtem(npy-1)-tehot)/(coory(npy-1)-1.0d0)
  end if

  ! tangent stress (simplified)
  if (walla) then
     if (.not.iex_avtan) then
        print *, 'AVTAN FILES NOT FOUND, NECESSARY WITH WALL LAW'
        stop
     end if
     tauwa(1) = abs(hom_avtan(1, 1))
     tauwa(2) = abs(hom_avtan(1,npy))
     print *, 'avtan', hom_avtan(1,1), hom_avtan(1,npy)
  else
     tauwa(1) = visco(1)*grave(1)
     tauwa(2) = visco(2)*grave(2)
  end if
     tauwa_av=0.5*(tauwa(1)+tauwa(2)) 

  if (.not.kfl_incom) then
     qwall(1) = abs( condu(1)*grate(1))
     qwall(2) = abs( condu(2)*grate(2))
  end if
  
  ustar(1) = sqrt(tauwa(1)/densi(1))
  ustar(2) = sqrt(tauwa(2)/densi(2))
  ustar_av=0.5*(ustar(1)+ustar(2))
  print *, 'walld(1), walld(2)', walld(1), walld(2)  
  print *, 'ustar(1), ustar(2)', ustar(1), ustar(2)  

! Reytau cold and hot ()

  Retau(1) = densi(1)*ustar(1)*delta/visco(1)
  Retau(2) = densi(2)*ustar(2)*delta/visco(2)
  if (.not.kfl_incom) then
     tstar(1) = qwall(1)/densi(1)/cp/ustar(1)
     tstar(2) = qwall(2)/densi(2)/cp/ustar(2)
  end if
  
  !
  ! Bulk values (integration) 
  !
  rhbul =0.0d0
  vebul =0.0d0
  tebul =0.0d0
  if (kfl_incom) then
     do je = 1, ney ! element loop 
        dy = coory(je+1) - coory(je)
        rhois =  rhref
        velis =  0.5d0*(hom_avvel(1,je)+hom_avvel(1,je+1))
        rhbul = rhbul + rhois*dy
        vebul = vebul + velis*dy
        print *, 'velis, vebul',velis, vebul
     end do
  rhbul= 0.5d0*rhbul/delta
  vebul =0.5d0*vebul/(delta*rhbul)
  print *, 'bulks (delta, rho_b, v_b):', delta, rhbul, vebul
!!$  else
!!$     do je = 1, ney ! element loop 
!!$        dy = coory(je+1) - coory(je)
!!$     rhois =  0.5d0*(hom_avden(je)+hom_avden(je+1))
!!$     velis =  0.5d0*(hom_avvel(1,je)+hom_avvel(1,je+1))
!!$     temis =  0.5d0*(hom_avtem(je)+hom_avtem(je+1))
!!$     rhvis =  0.5d0*(hom_fvvel(1,je)+hom_fvvel(1,je+1)) !rho vel
!!$     rhbul = rhbul + rhois
!!$     vebul = vebul + rhvis
!!$     tebul = tebul + rhvis*temis !This is incorrect!  should be used favre average!
!!$     end do
!!$  rhbul= 0.5d0*rhbul/delta
!!$  vebul =0.5d0*vebul/delta/rhbul  
!!$     tebul =0.5d0*tebul/delta/rhbul/vebul
  end if



  
  ! Reynolds bulk (flow rate)
  Rebul = rhbul*vebul*delta/muref 
 
  if (.not.kfl_incom) then
     ! Nussel numbers  
     nusse(1) = 4.0d0*grate(1)/(tebul-tecol)
     nusse(2) = 4.0d0*grate(2)/(tebul-tehot)
     !  Heat flux parameter Bq
     Bqwal(1) = qwall(1)/(rhbul*cp*ustar(1))
     Bqwal(2) = qwall(2)/(rhbul*cp*ustar(2))
  end if

  ! maximum velocity
  maxve = hom_avvel(1,2)
  do jp =3, npy -1
     maxve = max(maxve, hom_avvel(1,jp))
  end do

  ! Friction factors
  Cfmax(1) = 2.0d0*tauwa(1)/rhbul/maxve/maxve
  Cfmax(2) = 2.0d0*tauwa(2)/rhbul/maxve/maxve
  Cfbul(1) = 2.0d0*tauwa(1)/rhbul/vebul/vebul  ! shoud be v*v !
  Cfbul(2) = 2.0d0*tauwa(2)/rhbul/vebul/vebul

  !
  ! Obtains fluctuation <a'b'>=<ab> -<a><b>
  ! 
 
  allocate (hom_avprerms(npy))
  allocate (hom_avverms(ndime,npy))
  allocate (hom_avterms(npy))
  do jp =1, npy
     do idime =1, ndime
!        hom_avve2(idime,jp)= sqrt( hom_avve2(idime,jp)-  hom_avvel(idime,jp)*hom_avvel(idime,jp) )
        hom_avverms(idime,jp)= sqrt( hom_avve2(idime,jp)-  hom_avvel(idime,jp)*hom_avvel(idime,jp) )
     end do
     hom_avvxy(1,jp)  = hom_avvxy(1,jp)  -  hom_avvel(1,jp)*hom_avvel(2,jp) 
     hom_avvxy(2,jp)  = hom_avvxy(2,jp)  -  hom_avvel(2,jp)*hom_avvel(3,jp) 
     hom_avvxy(3,jp)  = hom_avvxy(3,jp)  -  hom_avvel(1,jp)*hom_avvel(3,jp) !!! changed component 2 for 3 rus
     if (.not.kfl_incom) then
!        hom_avte2(jp)    = sqrt(hom_avte2(jp) -  hom_avtem(jp)* hom_avtem(jp))
        hom_avterms(jp)    = sqrt(hom_avte2(jp) -  hom_avtem(jp)* hom_avtem(jp))
        do idime =1, ndime
           hom_avtev(idime,jp)= hom_avtev(idime,jp)- hom_avtem(jp) *hom_avvel(idime,jp)
        end do
     end if
     ! pressure rms
        hom_avprerms(jp)    = sqrt(hom_avpr2(jp) -  hom_avpre(jp)* hom_avpre(jp))
  end do
  
!  wall profiles 
  !  Cold wall
  do jp =1, (npy+1)/2-1
     yplus(jp) = (coory(jp)-coory(1)+walld(1))*densi(1)*ustar(1) /muref
     vplus(jp) = hom_avvel(1,jp)/ustar(1)
     if (.not.kfl_incom)  tplus(jp) = (hom_avtem(jp)-tecol)/tstar(1)
  end do
  !  Hot wall
  do jp =(npy+1)/2, npy        
     yplus(jp) = (coory(npy)-coory(jp)+walld(2))*densi(2)*ustar(2) /muref
     vplus(jp) = hom_avvel(1,jp)/ustar(2)
     if (.not.kfl_incom) tplus(jp) = (tehot-hom_avtem(jp))/tstar(2)
  end do

  
  !
  !  WRITES RESULTS 
  !

  fil_dummi = trim(file)//'_MEANPROFILE'
  OPEN (UNIT= lun_outmn, FILE= fil_dummi)  
  fil_dummi = trim(file)//'_MEAN_DW'
  OPEN (UNIT= lun_outwc, FILE =fil_dummi)
  fil_dummi = trim(file)//'_MEAN_UP'
  OPEN (UNIT= lun_outwh, FILE =fil_dummi)
  write(lun_outmn,'(a)')'# ycoor 2:avvex 3:avvey 4:avvez  5:avvex2  6:avvey2  7:avvez2  8:avxvy  9:avyvz  10:avxvz  11:avtem  12:avte2   13:avtevelx  14: avtevely  15: avtevelz  16: avpre  17:  avpr2  18: avprerms'      

  do jp = 1 , npy
     write(lun_outmn,110)  coory(jp), hom_avvel(1:3, jp), hom_avve2(1:3, jp), &
             hom_avvxy(1:3, jp), hom_avtem(jp), hom_avte2(jp), hom_avtev(1:3, jp),&
             hom_avpre(jp), hom_avpr2(jp), hom_avprerms(jp)
  end do
 
 
  ! COLD WALL
  write(lun_outwc,'(a)')'# yplus 2:vxplus  3: tplus 4: vxprms  5:  vyprms  6: vzprms 7: tprms'
  do jp = 1 , (npy+1)/2-1
     write(lun_outwc,110)  yplus(jp), vplus(jp), tplus(jp), hom_avverms(1:3,jp)/ustar(1), hom_avterms(jp)/tstar(1)
  end do
  ! HOT WALL
  write(lun_outwh,'(a)')'# yplus 2:vxplus  3: tplus 4: vxprms  5:  vyprms  6: vzprms 7: tprms'
  do jp = (npy+1)/2, npy
     write(lun_outwh,110)  yplus(jp), vplus(jp), tplus(jp), hom_avverms(1:3,jp)/ustar(2), hom_avterms(jp)/tstar(2)
  end do
  
  close(lun_outmn)
  close(lun_outwc)
  close(lun_outwh)

     deallocate (hom_avprerms)
     deallocate (hom_avverms)
     deallocate (hom_avterms)

  ! writes y1plus and nusselt distribution along wall
  fil_dummi = trim(file)//'_MEANGLOBAL'
  jp = 1
  if (walla) jp = 0
  OPEN (UNIT= lun_y1plu, FILE= fil_dummi)  
  write(lun_y1plu,'(a)')  '# 	1:Reycol 	2:Reyhot 	3:Rebulk 	 4:Nucol  	 5:Nuhot 	 6:Cfcol 	 7:Cfhot 	8:Bqcold 	 9:Bqhot 	10: Cf_bul_dw 	11: Cd_bul_up 	 12: y+_dw  	 13: y+_up 	14: maxvel  	  15:tauwal_col	    16:tauwal_hot'
  write(lun_y1plu,110)  Retau(1), Retau(2), Rebul, Nusse(1), Nusse(2),  Cfmax(1), Cfmax(2), Bqwal(1), Bqwal(2), Cfbul(1), Cfbul(2), yplus(1+jp), yplus(npy-jp), maxve, tauwa(1), tauwa(2) 
  close(lun_y1plu)


110 format(18(f17.8))
111 format(20(f17.8))

end program Read_Cav3D


subroutine readavfile(npoin,ndime,npx,npy,npz,fil_avvar,nfile, nfile_avvar,hom_avvar)
  !-----------------------------------------------------------------------
  !
  !   reads avverage file named fil_avvar
  !   the homogeneized variable will be  hom_avvar
  !   INPUT :
  !     npoin:  number of computational points
  !     ndime:  number of dimensions
  !     npx, npy, npz:   number of nodes in x, y and z
  !     fil_avvar:  name of the file containing data
  !     nfile  : max number of files
  !   OUTPUT :
  !     nfile_avvar : number of averaged files
  !     hom_avvar   : vector containing aveaged data
  !-----------------------------------------------------------------------
  implicit none
  integer, intent(in)        :: npoin, ndime, npx, npy, npz,nfile
  real,   intent(out)        :: hom_avvar(ndime, npy)
  integer, intent(out)       :: nfile_avvar
  character(40), intent(in)  :: fil_avvar
  ! local variables
  logical :: iex
  integer :: ifile, idime, ipoin, ip, jp, kp, lun_avvar, inifile
  character(40) :: label, fil_dummi
  real, allocatable :: avvar(:,:,:)

  lun_avvar = 299 
  print *, 'file, npoin,ndime, nfile', fil_avvar, npoin, ndime, nfile
  allocate (avvar(ndime, npoin, 300))
  iex=.true.
  ifile = 1  
  write (label,'(i3)') ifile
  label = '00'//trim(adjustl(label))

  DO WHILE (IEX)  ! looks if exists new file
     fil_dummi = trim(fil_avvar)//trim(adjustl(label))
     OPEN (UNIT= lun_avvar, FILE= fil_dummi) 
!     print *, 'ifile=',ifile
     READ(lun_avvar,*)
     READ(lun_avvar,*)
     do while (label(1:5).ne.'coord')
        READ(lun_avvar,*) label
     end do
     do idime =1, ndime
        DO ipoin =1, npoin
           read(lun_avvar,*)  avvar(idime,ipoin,ifile) 
        end do
     END DO
     CLOSE(lun_avvar)
     ifile = ifile + 1
     write (label,'(i3)') ifile
     
     if (ifile.lt.10) then
        label = '00'//trim(adjustl(label))
     elseif(ifile.lt.100) then
        label = '0'//trim(adjustl(label))
    
     end if
     fil_dummi = trim(fil_avvar)//trim(adjustl(label))
     inquire(file=fil_dummi, EXIST=iex) 
  END DO
  nfile_avvar = ifile -1
  print *,fil_avvar, '# AVVAL FILES READED=', nfile_avvar


  ! homogeneization along periodic direction X ,Z and interpolation

  inifile = max(nfile_avvar - nfile +1,1)
  print *, 'reading from file:', inifile
  do jp=1, npy
     do idime =1, ndime
        hom_avvar(idime,jp) = 0.0d0
        do kp=1, npz
           do ip =1, npx
              ipoin =  (ip-1)*npy*npz + (jp-1)*npz + kp
              do ifile =inifile, nfile_avvar
                 hom_avvar(idime,jp)= hom_avvar(idime,jp) + avvar(idime, ipoin, ifile)
              end do
           end do
        end do
        hom_avvar(idime,jp) =  hom_avvar(idime,jp)/dfloat((nfile_avvar-inifile+1)*npx*npy)
     end do
  end do
  print *, 'inifile, nfile_avvar', inifile, nfile_avvar 
  deallocate(avvar)

end subroutine readavfile
