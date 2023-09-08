subroutine ibm_iniunk()
  !-----------------------------------------------------------------------
  !****f* ibm_iniunk/ibm_iniunk
  ! NAME
  !    ibm_iniunk
  ! DESCRIPTION
  !    This routines computes for each particle:
  !    1. The volume 
  !    2. The center of gravity
  !    3. The inertia tensor
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_elmtyp
  use def_immbou
  use mod_kdtree
  use mod_memchk
  use mod_cutele
  use mod_messages, only : livinf

  implicit none
  integer(ip)          :: iboib,igaub,pblty,pgaub,iimbo,ipoib,inoib,isize
  integer(ip)          :: idime,iiner,nsize,pnodb,dummi,ielem
  integer(ip)          :: posit(3),posix,posiy,posiz,ipoin
  integer(ip)          :: nfaib,nboxe(3)
  integer(4)           :: istat
  real(rp)             :: baloc(ndime,ndime),bocod(ndime,max(mnoib,mnodb)),xfact
  real(rp)             :: gbsur,eucta,bouno(3),gpcib(3),rx,ry,rz,nx,ny,nz
  real(rp)             :: rx3,ry3,rz3,r3,onov3,onovn,maxdi
  real(rp),    pointer :: inerc(:,:)
  integer(ip), pointer :: kfl_calgr(:)
  real(rp)             :: delta(3)
  character(20)        :: messa
  !
  ! Allocate memory: solver, etc.
  !
  call ibm_memall(0_ip)

  if( kfl_rstar /= 0 ) then

     !----------------------------------------------------------------------
     !
     ! Read restart file
     !
     !----------------------------------------------------------------------

     call ibm_restar(1_ip)   

  else

     !----------------------------------------------------------------------
     !
     ! Define body fitted 
     !
     !---------------------------------------------------------------------- 

     call ibm_bodfit()

     !----------------------------------------------------------------------
     !
     ! Initialization
     !
     !---------------------------------------------------------------------- 

     call livinf(69_ip,' ',0_ip)
     dtime_ibm = dtime
     do iimbo = 1,nimbo
        do idime = 1,ndime
           imbou(iimbo) % posil(idime,2) = imbou(iimbo) % posil(idime,1)
           imbou(iimbo) % velol(idime,2) = imbou(iimbo) % velol(idime,1)
           imbou(iimbo) % accel(idime,2) = imbou(iimbo) % accel(idime,1)

           imbou(iimbo) % posia(idime,2) = imbou(iimbo) % posia(idime,1)
           imbou(iimbo) % veloa(idime,2) = imbou(iimbo) % veloa(idime,1)
           imbou(iimbo) % accea(idime,2) = imbou(iimbo) % accea(idime,1)
        end do
     end  do
     nsize = 5 * ndime - 9                                           ! =1 in 2D, 6 in 3D
     onov3 = 1.0_rp / 3.0_rp
     onovn = 1.0_rp / real(ndime)
     allocate(inerc(nsize,nimbo),stat=istat)
     call memchk(zero,istat,memor_dom,'INERC','ibm_iniunk',inerc)
     allocate(kfl_calgr(nimbo),stat=istat)
     call memchk(zero,istat,memor_dom,'CALGR','ibm_iniunk',kfl_calgr)
     do iimbo = 1,nimbo   ! flag to decide if posgr is calculated or read, it is set based on the initial posgr
        kfl_calgr(iimbo) = 1_ip
     end do

     do iimbo = 1,nimbo
        if(imbou(iimbo) % posgr(1) < -0.5e12_rp) then
           do idime = 1,ndime
              imbou(iimbo) % posgr(idime) = 0.0_rp
           end do
        else
           kfl_calgr(iimbo) = 0_ip
        end if
        imbou(iimbo) % volum = 0.0_rp
        imbou(iimbo) % maxdi = 0.0_rp
     end do

     !----------------------------------------------------------------------
     !
     ! Center of gravity and volume
     !
     !----------------------------------------------------------------------

     do iimbo = 1,nimbo
        do iboib = 1,imbou(iimbo) % nboib
           pblty = imbou(iimbo) % ltyib(iboib) 
           pnodb = nnode(pblty)
           pgaub = ngaib(pblty)
           !
           ! BOCOD: Gather
           !
           do inoib = 1,pnodb
              ipoib = imbou(iimbo) % lnoib(inoib,iboib)
              do idime = 1,ndime
                 bocod(idime,inoib) = imbou(iimbo) % cooib(idime,ipoib)                 
              end do
           end do
           !
           ! Loop over Gauss points
           !
           do igaub = 1,pgaub
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty) % derib(1,1,igaub),&  
                   bocod,baloc,eucta)                 
              gbsur = elmar(pblty) % weiib(igaub)*eucta 
              !
              ! GPCIB: Coordinates of Gauss point
              !
              do idime = 1,ndime
                 gpcib(idime) = 0.0_rp
                 do inoib = 1,pnodb
                    gpcib(idime) = gpcib(idime) + &
                         bocod(idime,inoib) * elmar(pblty) % shaib(inoib,igaub)
                 end do
              end do
              !
              ! BOUNO: Exterior normal
              !
              call exteib(iimbo,pnodb,iboib,bouno,0_ip)
              !
              ! Center of gravity: Integrate xi^2.ni
              ! Volume: Integrate xi.ni
              !
              do idime = 1,ndime
                 xfact                       = gbsur * gpcib(idime) * bouno(idime)
                 if (kfl_calgr(iimbo) == 1_ip)   imbou(iimbo) % posgr(idime) = imbou(iimbo) % posgr(idime) + xfact * gpcib(idime) 
                 imbou(iimbo) % volum        = imbou(iimbo) % volum        + xfact
              end do
           end do

        end do

     end do
     !
     ! For body fitted type bodies, sum in parallel
     ! Take absolute value of volume as for body fitted boundaries, the normal points outwards
     !
     do iimbo = 1,nimbo
        if( imbou(iimbo) % kfl_typeb >= 1 ) then
           if (kfl_calgr(iimbo) == 1_ip)  call pararr('SUM',0_ip,ndime,imbou(iimbo) % posgr)
           call pararr('SUM',0_ip, 1_ip,imbou(iimbo) % volum)
           imbou(iimbo) % volum = -imbou(iimbo) % volum
           if (kfl_calgr(iimbo) == 1_ip)  then
              do idime = 1,ndime
                 imbou(iimbo) % posgr(idime) = -imbou(iimbo) % posgr(idime)
              end do
           end if
        end if
     end do

     !
     ! Correct cog and volume
     !
     do iimbo = 1,nimbo
        if (kfl_calgr(iimbo) == 1_ip)  then
           do idime = 1,ndime
              imbou(iimbo) % posgr(idime) = imbou(iimbo) % posgr(idime) * 0.5_rp
           end do
        end if
        imbou(iimbo) % volum = imbou(iimbo) % volum * onovn
     end do
     !
     ! Set posil = posgr if posil has not been read in ibm.dat. This is the logical choice for body fitted cases
     !
     do iimbo = 1,nimbo
        if(imbou(iimbo) % posil(1,1) < -0.5e12_rp) then
           do idime = 1,ndime
              imbou(iimbo) % posil(idime,1) = imbou(iimbo) % posgr(idime)
              imbou(iimbo) % posil(idime,2) = imbou(iimbo) % posgr(idime)
           end do
        end if
     end do

     !----------------------------------------------------------------------
     !
     ! Tensor of inertia
     !
     !----------------------------------------------------------------------

     do iimbo = 1,nimbo

        if(imbou(iimbo) % momin(1) < -0.5_rp) then  ! momin has been initialized to -1.0 to indicate that it must be calculated if 
           ! it is not read 

           do iboib = 1,imbou(iimbo) % nboib

              pblty = imbou(iimbo) % ltyib(iboib) 
              pnodb = nnode(pblty)
              pgaub = ngaib(pblty)
              !
              ! BOCOD: Gather
              !
              do inoib = 1,pnodb
                 ipoib = imbou(iimbo) % lnoib(inoib,iboib)
                 do idime = 1,ndime
                    bocod(idime,inoib) = imbou(iimbo) % cooib(idime,ipoib)
                 end do
              end do
              !
              ! Loop over Gauss points
              !
              do igaub = 1,pgaub

                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty) % derib(1,1,igaub),&  
                      bocod,baloc,eucta)                 
                 gbsur = elmar(pblty) % weiib(igaub)*eucta 
                 !
                 ! GPCIB: Coordinates of Gauss point
                 !
                 do idime = 1,ndime
                    gpcib(idime) = 0.0_rp
                    do inoib = 1,pnodb
                       gpcib(idime) = gpcib(idime) + &
                            bocod(idime,inoib) * elmar(pblty) % shaib(inoib,igaub)
                    end do
                 end do
                 !
                 ! BOUNO: Exterior normal
                 !
                 call exteib(iimbo,pnodb,iboib,bouno,0_ip)
                 !
                 ! Inertia tensor
                 !
                 if ( ndime == 2 ) then

                    rx = gpcib(1) - imbou(iimbo) % posgr(1)
                    ry = gpcib(2) - imbou(iimbo) % posgr(2)
                    nx = bouno(1) 
                    ny = bouno(2) 

                    inerc(1,iimbo) = inerc(1,iimbo) + gbsur * ( rx**3 * nx + ry**3 * ny ) * onov3

                 else if ( ndime == 3 ) then

                    rx  = gpcib(1) - imbou(iimbo) % posgr(1)
                    ry  = gpcib(2) - imbou(iimbo) % posgr(2)
                    rz  = gpcib(3) - imbou(iimbo) % posgr(3)
                    rx3 = rx * rx * rx * onov3
                    ry3 = ry * ry * ry * onov3
                    rz3 = rz * rz * rz * onov3
                    nx  = bouno(1) * gbsur
                    ny  = bouno(2) * gbsur
                    nz  = bouno(3) * gbsur
                    r3  = rx3 * nx + ry3 * ny + rz3 * nz 

                    inerc(1,iimbo) = inerc(1,iimbo) + ry3 * ny + rz3 * nz     ! I11
                    inerc(2,iimbo) = inerc(2,iimbo) + rx3 * nx + rz3 * nz     ! I22
                    inerc(3,iimbo) = inerc(3,iimbo) + rx3 * nx + ry3 * ny     ! I33
                    inerc(4,iimbo) = inerc(4,iimbo) - rx*rx*ry * 0.5_rp * nx  ! I12
                    inerc(5,iimbo) = inerc(5,iimbo) - rz*rz*rx * 0.5_rp * nz  ! I13
                    inerc(6,iimbo) = inerc(6,iimbo) - ry*ry*rz * 0.5_rp * ny  ! I23

                 end if

              end do

           end do

        end if

     end do
     !
     ! For body fitted type bodies, sum in parallel
     !
     do iimbo = 1,nimbo
        if(imbou(iimbo) % momin(1) < -0.5_rp) then  ! momin has been initialized to -1.0 to indicate that it must be calculated if 
           ! it is not read
           if( imbou(iimbo) % kfl_typeb >= 1 ) then
              call pararr('SUM',0_ip,nsize,inerc(1,iimbo))
              do isize = 1,nsize
                 inerc(isize,iimbo) = -inerc(isize,iimbo)
              end do
           end if
        end if
     end do

     !----------------------------------------------------------------------
     !
     ! Compute density form mass or mass from density
     !
     !----------------------------------------------------------------------

     do iimbo = 1,nimbo
        if( imbou(iimbo) % densi < 0.0_rp .and. imbou(iimbo) % massa < 0.0_rp ) then
           if( imbou(iimbo) % kfl_typeb == 0 .or. imbou(iimbo) % kfl_typeb == -2 ) then
              call runend('IBM_INIUNK: PARTICLE MASS OR DENSITY MUST BE PRESCRIBED')
           else
              imbou(iimbo) % densi = 1.0_rp
              imbou(iimbo) % massa = imbou(iimbo) % densi * imbou(iimbo) % volum
           end if
        else if( imbou(iimbo) % densi < 0.0_rp ) then
           imbou(iimbo) % densi = imbou(iimbo) % massa / imbou(iimbo) % volum
        else if( imbou(iimbo) % massa < 0.0_rp ) then
           imbou(iimbo) % massa = imbou(iimbo) % densi * imbou(iimbo) % volum
        else
           !           now we allow to prescribe both mas and density even if they are not consistent with the volume. A higher mass is a 
           !           way of relaxing the movement to obtain the same stationary solution.
        end if
        if( imbou(iimbo) % volum <= 0.0_rp ) then
           messa = intost(iimbo)                     
           call runend('IBM_INIUNK: NEGATIVE VOLUME FOR PARTICLE '//trim(messa)//'. CHECK BOUNDARY ORIENTATION')
        end if
     end do

     !----------------------------------------------------------------------
     !
     ! Multiply inertia tensor by density
     !
     !----------------------------------------------------------------------

     do iimbo = 1,nimbo
        if(imbou(iimbo) % momin(1) < -0.5_rp) then  ! momin has been initialized to -1.0 to indicate that it must be calculated if 
           ! it is not read
           do iiner = 1,nsize
              imbou(iimbo) % momin(iiner) = inerc(iiner,iimbo) * imbou(iimbo) % densi
           end do
        end if
     end do

     !----------------------------------------------------------------------
     !
     ! Modify COOIB so that center of gravity is at prescribed position
     ! 
     !----------------------------------------------------------------------

     do iimbo = 1,nimbo
        !
        ! Allocate the number of bouding boxes faces
        !
        allocate( imbou(iimbo) % fabox(2,ndime,imbou(iimbo) % nboib ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'FABOX','ibm_iniunk',imbou(iimbo) % fabox)

        imbou(iimbo) % maxdi = 0.0_rp          
        do ipoib = 1,imbou(iimbo) % npoib
           maxdi = 0.0_rp
           do idime = 1,ndime
              imbou(iimbo) % cooin(idime,ipoib) = imbou(iimbo) % cooib(idime,ipoib) - imbou(iimbo) % posgr(idime)
              imbou(iimbo) % cooib(idime,ipoib) = imbou(iimbo) % cooin(idime,ipoib) + imbou(iimbo) % posil(idime,1)
              imbou(iimbo) % cooi2(idime,ipoib) = imbou(iimbo) % cooib(idime,ipoib)              
              maxdi = maxdi + imbou(iimbo) % cooin(idime,ipoib)*imbou(iimbo) % cooin(idime,ipoib)              
           end do
           imbou(iimbo) % maxdi = max(imbou(iimbo) % maxdi,sqrt(maxdi))
        end do
     end do

     call memchk(two,istat,memor_dom,'INERC','ibm_iniunk',inerc)
     deallocate(inerc,stat=istat)
     if(istat/=0) call memerr(two,'INERC','ibm_iniunk',0_ip)
     call memchk(two,istat,memor_dom,'CALGR','ibm_iniunk',kfl_calgr)
     deallocate(kfl_calgr,stat=istat)
     if(istat/=0) call memerr(two,'CALGR','ibm_iniunk',0_ip)

     !---------------------------------------------------------------
     !
     ! Allocate and build the variables for tree structure
     !
     !---------------------------------------------------------------     

     do iimbo = 1,nimbo   
        nfaib = imbou(iimbo) % nboib
        !allocate(imbou(iimbo) % sabox(2,ndime,2_ip*nfaib-1_ip) ,stat=istat)
        !call memchk(zero,istat,mem_modul(1:2,modul),'SABOX','ibm_iniunk',imbou(iimbo) % sabox)
        !allocate(imbou(iimbo) % blink(2_ip*nfaib-1_ip)          ,stat=istat)
        !call memchk(zero,istat,mem_modul(1:2,modul),'BLINK','ibm_iniunk',imbou(iimbo) % blink)    
        !
        ! Create the bounding box for each particle and for each face in a particle
        ! Also, create tree structures used by find points inside the particle
        !      
        if( IMASTER .or. imbou(iimbo) % kfl_typeb >= 1 ) then
           continue
        else
           call kdtree(&
                1_ip,mnoib,imbou(iimbo) % npoib,imbou(iimbo) % nboib,&
                imbou(iimbo) % cooib,imbou(iimbo) % lnoib,imbou(iimbo) % ltyib,&
                imbou(iimbo) % fabox,imbou(iimbo) % bobox,imbou(iimbo) % sabox,&
                imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                imbou(iimbo) % lnele)     
        end if
     end do

     !---------------------------------------------------------------
     !
     ! Allocate and build the variables for cut elements
     !
     !---------------------------------------------------------------          

     if( INOTMASTER ) call buicut(1_ip)

     !---------------------------------------------------------------
     !
     ! Put the nodes of background mesh (fluid mesh) in differet boxes
     ! for a bucket sort
     !
     !---------------------------------------------------------------


     if( INOTMASTER ) then

        allocate(xmima_ibm(2,ndime),stat=istat)

        !
        ! Determine the number of boxes in an axis
        !
        nubox_ibm = int( real(npoin_2,rp)**(1.0_rp/real(ndime,rp)) / 10.0_rp )
        nubox_ibm = max(2_ip,nubox_ibm)

        do idime = 1,ndime
           xmima_ibm(1,idime) =  1.0e12_rp
           xmima_ibm(2,idime) = -1.0e12_rp
        end do

        do ipoin = 1,npoin_2
           do idime = 1,ndime
              xmima_ibm(1,idime) = min( xmima_ibm(1,idime) , coord(idime,ipoin) )
              xmima_ibm(2,idime) = max( xmima_ibm(2,idime) , coord(idime,ipoin) )
           end do
        end do

        nboxe(3) = 1
        do idime = 1,ndime
           nboxe(idime) = nubox_ibm
           delta(idime) = (xmima_ibm(2,idime)-xmima_ibm(1,idime)) / real(nboxe(idime),rp) 
           delta(idime) = delta(idime) + (delta(idime)/1000.0_rp)
        end do


        allocate(boxes_ibm(nboxe(1),nboxe(2),nboxe(3)),stat=istat)

        !
        ! Initialize the number of nodes for each box
        !
        do posix = 1,nboxe(1)
           do posiy = 1,nboxe(2)
              do posiz = 1,nboxe(3)
                 boxes_ibm(posix,posiy,posiz) % nnode = 0_ip
              end do
           end do
        end do

        !
        ! Determine the number of nodes in each box 
        !
        posit(3) = 1
        do ipoin = 1,npoin_2
           do idime = 1,ndime
              posit(idime) = int ( (coord(idime,ipoin) - xmima_ibm(1,idime))/delta(idime) ) + 1_ip
           end do
           boxes_ibm(posit(1),posit(2),posit(3)) % nnode = boxes_ibm(posit(1),posit(2),posit(3)) % nnode + 1
        end do

        !
        ! Allocate the space for each box
        !

        do posix = 1,nboxe(1)
           do posiy = 1,nboxe(2)
              do posiz = 1,nboxe(3)
                 dummi = boxes_ibm(posix,posiy,posiz) % nnode
                 allocate(boxes_ibm(posix,posiy,posiz) % nodes(dummi),stat=istat)          
              end do
           end do
        end do

        !
        ! Again, reinitialize the number of nodes for each box
        !
        do posix = 1,nboxe(1)
           do posiy = 1,nboxe(2)
              do posiz = 1,nboxe(3)
                 boxes_ibm(posix,posiy,posiz) % nnode = 0_ip
              end do
           end do
        end do
        !
        ! Store each node number in the fluid mesh in the corresponding box 
        !

        posit(3) = 1
        do ipoin = 1,npoin_2
           do idime = 1,ndime
              posit(idime) = int ( (coord(idime,ipoin) - xmima_ibm(1,idime))/delta(idime) ) + 1_ip
           end do
           boxes_ibm(posit(1),posit(2),posit(3)) % nnode = boxes_ibm(posit(1),posit(2),posit(3)) % nnode + 1
           dummi = boxes_ibm(posit(1),posit(2),posit(3)) % nnode
           boxes_ibm(posit(1),posit(2),posit(3)) % nodes(dummi) = ipoin     
        end do

     end if

     !---------------------------------------------------------------------
     !
     ! LNNIB, DISPM, LNTIB, LNDIB, LETIB, LNFIB, LNTRA, MASSC
     !
     !---------------------------------------------------------------------

     if( INOTMASTER ) then

        !allocate(lndib_ibm(npoin),stat=istat)
        !call memchk(zero,istat,mem_modul(1:2,modul),'LNDIB_IBM','ibm_iniunk',lndib_ibm)
        !allocate(lnfib_ibm(npoin),stat=istat)
        !call memchk(zero,istat,mem_modul(1:2,modul),'LNFIB_IBM','ibm_iniunk',lnfib_ibm)

        allocate(lntib(npoin_2),stat=istat)
        call memchk(zero,istat,memor_dom,'LNTIB','ibm_iniunk',lntib)
        allocate(lnti2(npoin_2),stat=istat)
        call memchk(zero,istat,memor_dom,'LNTI2','ibm_iniunk',lnti2)

        allocate(letib(nelem_2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LETIB','ibm_iniunk',letib)

        allocate(lelch_ibm(nelem_2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LELCH_IBM','ibm_iniunk',lelch_ibm)

        allocate(lnoch_ibm(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LNOCH_IBM','ibm_iniunk',lnoch_ibm)

        if( kfl_diric_ibm == 1 ) then
           allocate(lntra(npoin),stat=istat)
           call memchk(zero,istat,memor_dom,'LNTRA','ibm_iniunk',lntra)     

           allocate(massc(ndime+1,npoin),stat=istat)
           call memchk(zero,istat,memor_dom,'MASSC','ibm_iniunk',massc)    

           allocate(lnint(npoin),stat=istat)
           allocate(lnin2(npoin),stat=istat)
           !
           ! LNINT: Number of nodes to interpolate
           !
           do ipoin = 1,npoin
              lnint(ipoin) % limit = 0_ip
              lnin2(ipoin) % limit = 0_ip
           end do
        end if
        !
        ! LELCH: save initial value 
        !
        do ielem = 1,nelem_2
           lelch_ibm(ielem) = lelch(ielem)
        end do
        !
        ! LNOCH: save initial value
        !
        do ipoin = 1,npoin
           lnoch_ibm(ipoin) = lnoch(ipoin)
        end do        
     end if
  end if

  npoib = 0
  nboib = 0
  do iimbo = 1,nimbo
     npoib = npoib + imbou(iimbo) % npoib
     nboib = nboib + imbou(iimbo) % nboib
  end do

  !---------------------------------------------------------------
  !
  ! Initial condition
  !
  !---------------------------------------------------------------

  bacdt_ibm = dtime
  !
  ! call ibm_wallib()  ! tiene bug en L245  : redifine mnodb!!!!!!!!!!!!!
  !
  ! Build the parallel graph
  !
  call ibm_domgra()

  call ibm_solite(0_ip)
  call ibm_updunk(1_ip)

  !
  ! Redefine MGAUS
  !
  if( ndime == 2 ) then
     mgaus = max(mgaus,6*ngaus(TRI03))
  else
     if( lexis(HEX08) /= 0 ) then
        mgaus = max(mgaus,24*ngaus(TET04))
     else
        mgaus = max(mgaus,4*ngaus(TET04))
     end if
  end if

end subroutine ibm_iniunk
