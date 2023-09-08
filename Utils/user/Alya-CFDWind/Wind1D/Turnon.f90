  subroutine turnon
  ! This routines reads initial and boundary condition

  use def_master
  implicit none

  real(8)      :: z,  ustar1, dz, distz, rageo, tnume, tdeno,  sqkey, gpcod
  real(8)      :: lm, kz
  real(8)      :: heights(2000), lam, ztoh, LADEN, n
  integer      :: flag, kpoin
  character*20 :: label, file
  logical      :: wind = .true. ! realizable for wind

  !
  !*** Sets model constants
  !
  if  (kfl_model==0) then !CFDW2
     !********** Original Jones-Launder , 1972 **********
     ! write(*,*) 'Note: Using Jones-Launder (1972) constants'
     ! kar =0.41d0
     ! cmu=0.09d0
     ! c1= 1.44d0
     ! c2=1.92d0
     ! sigka=1.0d0
     ! sigep=1.11d0
     !
     !********** Crespo et al  ! From Crespo et al (Used in Alya) **********
      write(*,*) 'Note: Using Crespo et al. constants (Used in Alya)'
      kar =0.41d0
      cmu=0.0333d0
      c1= 1.176d0
      c2= 1.92d0
      sigka= 1.0d0
      sigep= kar*kar/(c2-c1)/sqrt(cmu) 
     !
     !********** Detering **********
     ! kar = 0.4
     ! cmu=0.0256d0
     ! c1= 1.13d0
     ! c2=1.90d0
     ! sigka=0.7407d0
     ! sigep=1.2987d0
     !
     !********** Sogachev **********
    ! write(*,*) 'Note: Using Sogachev constants'
     ! kar =0.4d0
     ! cmu=0.09d0
     ! c1= 1.52d0
     ! c2= 1.833d0
     ! sigka= 1.0d0
     ! sigep= kar*kar/(c2-c1)/sqrt(cmu)
     !
     !********** Koblitz / GABLS3 and Meso-Micro benchmark **********
     ! write(*,*) 'Note: Using Koblitz constants'
     ! kar =0.4d0
     ! cmu= 0.03d0
     ! c1=  1.52d0
     ! c2=  1.833d0
     ! sigep= kar*kar/(c2-c1)/sqrt(cmu)
     ! sigka= sigep
     !
     if (kfl_canmo == 2.and..false.) then ! SANZ model
        ! Atmospheric (Lopes da Costa)
        kar =0.40d0
        cmu= 0.0333d0
        c1=  1.44d0
        c2=  1.92d0
        sigka = 1.0d0
        sigep = kar*kar/(c2-c1)/sqrt(cmu)
     end if
     !     
  else if (kfl_model==1) then ! RNG
     kar =0.41d0
     cmu=0.085d0
     c1= 1.42d0
     c2=1.68d0
     sigka=0.7179d0
     sigep=0.7179d0  
  else if (kfl_model==2) then ! REALIZABLE
     kar =0.41d0
     cmu=0.09d0
     c1= 0.0d0
     c2=1.9d0
     sigka=1.0d0
     sigep=1.2d0  
     if (wind) then
        cmu=0.0333d0      
        sigep= kar*kar/(c2*sqrt(cmu)-0.43*sqrt(cmu/0.09d0))    
     end if
  else
     write(*,*) 'ERROR, NO CFWIND1 NO CFWIND2'
     stop
  end if
 
  ! Monin Obukhov wall law (Kansas 1)
  mo_prand = 0.74d0
  mo_gamm1 = 15.0d0 
  mo_gamm2 = 9.0d0
  mo_beta1 = 4.7d0
  mo_beta2 = 4.7d0  

  ! Monin Obukhov wall law (Kansas 2) (Alinot Masson)
  ! mo_prand = 1.0d0
  ! mo_gamm1 = 16.0d0 
  ! mo_gamm2 = 16.0d0
  ! mo_beta1 = 5.0d0
  ! mo_beta2 = 5.0d0  
  sigte = 0.74d0
  cmu0=cmu 

  !
  !*** Mesh Coordinates
  !
  allocate (coord(npoin))
  dz = dz1
  rageo = 1.01d0
  ielem = int(nelem*1.0)
  do while ((rageo**(ielem)-1.0d0)/(rageo-1.0d0).lt.length/dz1) 
     rageo = rageo + 0.002d0
  end do
  write (*,*) 'geometrical ratio (r+2)=', rageo+0.02
  flag=0
  ! creates mesh
  ! Loop along each node 
  do ipoin=1, npoin
     ! Defines each node height
     if (ipoin.eq.1) z=dwall
     if (ipoin.eq.2) z=dwall +  dz 
     if (ipoin.ne.1.and.ipoin.ne.2) then        
        if (z.lt.6.0*dz1) then 
           dz = dz*(rageo+0.02d0) !max(dz*1.5d0, 2.0d0*c/dfloat(nez))
        else if (z.lt.40.0*dz1) then 
           dz = dz*(rageo+0.02d0) !max(dz*1.5d0, 2.0d0*c/dfloat(nez))
        else
           dz=  dz*(rageo+0.02d0)
        end if
        distz= (length-z)/dfloat(npoin-ipoin+1)
        if (dz.gt.distz.or.flag==1) then           
           if (flag==0) then
              write (*,*)  'uniform length factor=', (z+distz)/length, &
                   'at z= ', z+distz, 'ipoin=', ipoin,'unilength=',distz, dz
           end if
           dz = distz      
           flag=1
        end if
        z = z + dz
     end if
     ! Sets node vertical coordinate
     coord(ipoin)= z 
  end do !ipoin
  !Check
  if (flag==0) then
     write(*,*) 'ERROR constructing mesh, correct input parameters \\ bad length, factor=', z/length, 'rageo', rageo
     stop
  end if


  !
  !*** GAUSS POINTS AND WEIGHTS
  !
  if (kfl_order==2) then  
     ngaus= 3
     nnode= 3
     npoin = 2*nelem +1
  end if
  allocate(weigp(ngaus))
  allocate(posgp(ngaus))
  if (ngaus==2.and.kfl_close==1) then
     write(*,*) 'closed integration rule'
     posgp(1)=-1.0d0
     posgp(2)= 1.0d0
     weigp(1)= 1.0d0
     weigp(2)= 1.0d0
  else if (ngaus==2) then
     posgp(1)=-0.577350269189626d0
     posgp(2)= 0.577350269189626d0
     weigp(1)= 1.0d0
     weigp(2)= 1.0d0
  else if(ngaus==3) then
     posgp(1)=-0.774596669241483d0
     posgp(2)= 0.0d0
     posgp(3)= 0.774596669241483d0
     weigp(1)= 0.555555555555556d0
     weigp(2)= 0.888888888888889d0
     weigp(3)= 0.555555555555556d0
  else if(ngaus==4)  then
     posgp(1)=-0.861136311594053d0
     posgp(2)=-0.339981043584856d0
     posgp(3)= 0.339981043584856d0
     posgp(4)= 0.861136311594053d0
     weigp(1)= 0.347854845137454d0
     weigp(2)= 0.652145154862546d0
     weigp(3)= 0.652145154862546d0
     weigp(4)= 0.347854845137454d0
  end if
  ! Definition of shape functions
  allocate(shape(nnode, ngaus))
  do igaus=1, ngaus
     z= posgp(igaus)
     shape(1,igaus)= -0.5d0*(z-1.0d0)
     shape(2,igaus)=  0.5d0*(z+1.0d0)
  end do
  if (kfl_order==2) then
     allocate(deriv(nnode, ngaus))
     allocate(deri2(nnode, ngaus))
     do igaus=1, ngaus
        z= posgp(igaus)
        shape(1,igaus)=  0.5d0*z*(z-1.0d0)            !  1    3    2
        shape(2,igaus)=  0.5d0*z*(z+1.0d0)
        shape(3,igaus)= -(z+1.0d0)*(z-1.0d0)
        deriv(1,igaus)= z - 0.5d0
        deriv(2,igaus)= z + 0.5d0
        deriv(3,igaus)=-2.0d0*z
        deri2(1,igaus)= 1.0d0
        deri2(2,igaus)= 1.0d0
        deri2(3,igaus)= -2.0d0
     end do
  end if
  

  !
  !*** Allocates structures
  !
  allocate (lnods(nnode,nelem))
  allocate (unkno(npoin))
  allocate (rhsid(npoin))
  allocate (avect(npoin))
  allocate (diago(npoin))
  allocate (cvect(npoin))
  allocate (veloc(npoin, 2, 2))
  allocate (keyva(npoin, 2))
  allocate (epsil(npoin, 2))
  allocate (veloc_ini(npoin, 2))
  allocate (keyva_ini(npoin))
  allocate (epsil_ini(npoin))
  allocate (ia(npoin+1))
  allocate (tempe(npoin,2))
  allocate (tempe_ini(npoin))

  !
  !*** Matrix connectivities
  !
  if (kfl_order==1) then 
     nzdom = 3*npoin -2  
     allocate (ja(nzdom))
     allocate(amatr(nzdom))
     ! conectivity matrix
     do ielem = 1, nelem 
        do inode = 1, nnode
           lnods(inode, ielem) = ielem +inode -1      
        end do
     end do
     ! matrix index (csr)
     ia(1) = 1  
     ja(1) = 1
     ja(2) = 2
     izdom = 3

     do ipoin = 2, npoin -1
        ia(ipoin)= izdom
        ja(izdom)= ipoin -1
        ja(izdom+1)= ipoin
        ja(izdom+2)= ipoin +1 
        izdom =izdom + 3
     end do
     ia(npoin)= izdom
     ja(izdom)= npoin -1
     ja(izdom+1)= npoin 
     nzdom =izdom + 1
     ia(npoin+1) = nzdom +1
     
  else if (kfl_order==2) then ! second order conectivity matrix
     nzdom =5*(nelem-1) + 3*(nelem+2)
     allocate (ja(nzdom))
     allocate(amatr(nzdom))
     ! conectivity matrix
     ipoin =1
     do ielem =1, nelem
        lnods(1,ielem)= ipoin
        lnods(2,ielem)= ipoin +2
        lnods(3,ielem)= ipoin +1
        ipoin = ipoin +2
     end do
     ! matrix index (csr)
     ia(1) = 1  
     ja(1) = 1
     ja(2) = 2
     ja(3) = 3
     izdom = 4
     do ipoin = 2, npoin -1    
        ia(ipoin)= izdom
        if (mod(ipoin,2)==0) then !even
           ja(izdom)= ipoin -1
           ja(izdom+1)= ipoin
           ja(izdom+2)= ipoin +1 
           izdom =izdom + 3
        else  !odd
           ja(izdom)  = ipoin -2
           ja(izdom+1)= ipoin -1
           ja(izdom+2)= ipoin  
           ja(izdom+3)= ipoin +1
           ja(izdom+4)= ipoin +2
           izdom =izdom + 5
        end if
     end do
     ! last row
     ia(npoin)= izdom
     ja(izdom)= npoin -2
     ja(izdom+1)= npoin  -1
     ja(izdom+2)= npoin 
     nzdom =izdom + 2
     ia(npoin+1) = nzdom +1
  else
     write(*,*) 'ERROR: KFL_ORDER SHOULD BE ONE OR TWO'
  end if

  !
  !*** Initial conditions
  !
  ustar1 = vegeo*kar/log(1.0d0+length/rough)
  if (kfl_case.ne.0) ustar = ustar1

  ! Veloc, key, eps, and tempe initial profiles
  if (.not.restart_in) then
     do ipoin=1, npoin
        z= coord(ipoin)

        ! Default
        veloc_ini(ipoin,1) = ustar1/kar*log(1.0d0+coord(ipoin)/rough)        
        veloc_ini(ipoin,2) = 0.0d0    
        keyva_ini(ipoin) = ustar*ustar/sqrt(cmu)
        !     lm =kar*(z+rough)*l_max/(kar*(z+rough)+l_max)
        !     epsil_ini(ipoin) =  ((cmu*keyva_ini(ipoin)*keyva_ini(ipoin))**(0.75d0))/lm
        epsil_ini(ipoin) = ustar*ustar*ustar/(kar*(z+rough))
        ! temperature initizalization
        if (z.lt.ztemin) then
           tempe_ini(ipoin) =  tewal +  gradbo*z
        else
           tempe_ini(ipoin) =  tewal +  gradbo*ztemin + gradto*(z-ztemin)
        end if
        
        ! GABLS2 and Cycle problem type
        if (kfl_case==1.or.kfl_case==2) then
           call gabls2_ini
        end if

        ! GABLS3 problem type and TKE and EPS creation
        if (kfl_case.eq.3) then
           call gabls3_ini
        end if

        ! GABLS1
        if (kfl_case.eq.4) then  !initial conditions
           veloc_ini(ipoin,1) = vegeo
           keyva_ini(ipoin)   = max(keyam, 0.4 * (1.0d0 - z/250.0 )**3.0) 
           kz =  kar*(z+rough)
           lm = kz/(1.0+ kz/l_max)
           epsil_ini(ipoin)   = max(epsam, ((cmu*keyva_ini(ipoin)*keyva_ini(ipoin))**0.75)/ lm)
        end if
       
!
        ! Logarithmic
        if (kfl_logva) then
           keyva_ini(ipoin) = log(keyva_ini(ipoin))
           epsil_ini(ipoin) = log(epsil_ini(ipoin) )
        end if
     end do  ! ipoin
     
     ! Calculates initial Mellor-Yosida length
     tnume = 0.0d0
     tdeno = 0.0d0
     do ielem =1, nelem
        ipoin = ielem
        jpoin = ielem + 1
        dz = coord(jpoin)-coord(ipoin)
        sqkey = 0.50d0*(sqrt(keyva_ini(ipoin))+sqrt(keyva_ini(jpoin)))
        gpcod = coord(ipoin) +0.50*dz
        tnume = tnume + gpcod*sqkey*dz
        tdeno = tdeno + sqkey*dz
     end do
     lenmy = 0.075*tnume/tdeno
     ustar2=ustar

  else ! Restart Initial conditions
     call read_restart
  end if

  ! For fixed temperature at top
  tetop = tempe_ini(npoin)

  ! Velocity components on top for GABLS3 case type
  if (kfl_case.eq.3) then
     ugeos(1) =  veloc_ini(npoin,1)
     ugeos(2) =  veloc_ini(npoin,2)
  else
     ugeos(1) =  vegeo
     ugeos(2) =  0.0d0 ! geostrophic paralel to x direction
  end if
  
  !
  !*** Initialize unknown
  !
  do ipoin =1, npoin
     veloc(ipoin, 1, 1) = veloc_ini(ipoin,1)
     veloc(ipoin, 2, 1) = veloc_ini(ipoin,2)
     keyva(ipoin, 1)    = keyva_ini(ipoin)
     epsil(ipoin, 1)    = epsil_ini(ipoin)
     tempe(ipoin, 1)    = tempe_ini(ipoin)
     veloc(ipoin, 1 ,2) = veloc(ipoin, 1,1)
     veloc(ipoin, 2 ,2) = veloc(ipoin, 2,1)
     keyva(ipoin, 2)    = keyva(ipoin,1)
     epsil(ipoin, 2)    = epsil(ipoin,1)     
     tempe(ipoin ,2)    = tempe(ipoin,1)
  end do

  !
  !*** open files
  !
  if (.not.restart_in) then
     open(lun_conve(1), FILE='Wind2.ns1.cvg',status='unknown')
     open(lun_conve(2), FILE='Wind2.ns2.cvg',status='unknown')
     open(lun_conve(3), FILE='Wind2.tur1.cvg',status='unknown')
     open(lun_conve(4), FILE='Wind2.tur2.cvg',status='unknown')
     if (kfl_trtem) then
        open(lun_conve(5), FILE='Wind2.tem.cvg',status='unknown')
        open(13, FILE='Surface_temperatures.txt',status='unknown')
     end if
     write(lun_conve(1),101)
     write(lun_conve(2),102)
     write(lun_conve(3),201)
     write(lun_conve(4),202)
     if (kfl_trtem) then
        write(lun_conve(5),301)
        write(13,'(a)') '#ctime,tewal,tempe(1,1),T_2,T_sfc,qw,L,ustar,logft,Stab'
     end if
     if (kfl_trtem) then
        open(lun_globa, FILE='plotglobal', status='unknown' )
        write(lun_globa,'(a)')  '#1: time          2:ustar          3:heatfl       4:tstar    5: tewal    6: maxle  7: MOlength  8: htpbl  9: ugeos(1)  10: ugeos(2)  11: ugeos_module  12:veloc(npoin,1)   13: veloc(npoin,1)  14: veloc_module(npoin)'
     end if
  else
     open(lun_conve(1), FILE='Wind2.ns1.cvg',status='unknown',access="append")
     open(lun_conve(2), FILE='Wind2.ns2.cvg',status='unknown',access="append")
     open(lun_conve(3), FILE='Wind2.tur1.cvg',status='unknown',access="append")
     open(lun_conve(4), FILE='Wind2.tur2.cvg',status='unknown',access="append")
     if (kfl_trtem) then
        open(lun_conve(5), FILE='Wind2.tem.cvg',status='unknown',access="append")
        open(13, FILE='Surface_temperatures.txt',status='unknown',access="append")
     end if
     if (kfl_trtem) then
        open(lun_globa, FILE='plotglobal', status='unknown') !,access="append")
     end if
  end if
  
     
  !
  !*** CUTS AND TRANSIENT FILES
  !
  do icuts = 1, ncuts        
     write(label,'(f7.2)') cutpr(icuts)
     file = 'plotcut_'//trim(adjustl(label))  ! plotcut_zz.zz
   !  if (.not.restart_in) then
        open(lun_cutre(icuts), FILE=file, status='unknown')
        write(lun_cutre(icuts),'(a)')  '#1: time          2:velocx          3:velocy      &
             4:key     5:eps     6:temp   7:htflx '
  !   else
  !      open(lun_cutre(icuts), FILE=file, status='unknown') !,access="append")
  !   end if
     
     ! INTERPOLATION VALUES
     ipoin =1
     jpoin =npoin
     kpoin =npoin/2
     do while((jpoin-ipoin).gt.1)
        if(cutpr(icuts).lt.coord(kpoin)) then
           jpoin = kpoin
        else
           ipoin = kpoin
        end if
        kpoin =(ipoin + jpoin)/2   
     end do
     ielec(icuts) = ipoin ! element to which icuts belongs
  end do

#ifdef _NETCDF_      ! This allows to compile without NETCDF library
  if (kfl_netCDF) then
     call netCDF_define
     if (.not.restart_in.and..not.append_netCDF) &
          call netCDF_Postpr
  end if
#endif
  if (kfl_trtem)  call postpr
  !
  !*** Formats
  !
  if (kfl_canop.and.kfl_candi==1) then ! write LAD distribution
     open(lun_lad, FILE='LAD_prof',status='unknown')
     write (lun_lad, *) '# z/h   LAD'
     lam =2.0d0*LAD
     do ipoin =1, 1000
        ztoh = dfloat(ipoin)/1000.0d0
        
        if (ztoh.lt.0.70d0) then
           n =6.0
        else                      
           n =0.5
        end if
        LADEN = lam *( ( 0.30 /max(1.0d0 -ztoh,0.0001d0))**n) *exp(n*(1.0d0- 0.30d0 /max(1.0d0 - ztoh, 0.0001d0) ))
        write (lun_lad, '(2(e15.7))') ztoh, LADEN
     end do
  end if

101 format('$ ','       Time','     Global','      Inner',&
       &      '       Current','   x- Velocity','  x - Velocity'/,&          
       & '$ ','       step','  iteration','  iteration',&
       &      '          time','      residual', '   norm')

201 format('$ ','       Time','     Global','      Inner',&
       &      '       Current','      kinetic_tur',' kinet'/,&          
       & '$ ','       step','  iteration','  iteration',&
       &      '          time','      residual','      residual', '   norm')

102 format('$ ','       Time','     Global','      Inner',&
       &      '       Current','   y- Velocity','   y - Veloc'/,&          
       & '$ ','       step','  iteration','  iteration',&
       &      '          time','      residual', '   norm')

202 format('$ ','       Time','     Global','      Inner',&
       &      '       Current','      epsil_tur','   epsil'/,&          
       & '$ ','       step','  iteration','  iteration',&
       &      '          time','      residual', '   norm')

301 format('$ ','       Time','     Global','      Inner',&
       &      '       Current','     Temper','   Temper'/,&          
       & '$ ','       step','  iteration','  iteration',&
       &      '          time','      residual', '   norm')
end subroutine turnon
