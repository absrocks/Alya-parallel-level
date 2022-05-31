subroutine nsa_apply_physics

  !-----------------------------------------------------------------------
  !****f* Nastin/nsa_apply_physics
  ! NAME 
  !    nsa_apply_physics
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    nsa_updunk
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_elmtyp
  use def_domain
  use def_nastal
  use mod_memory
  implicit none

  integer(ip) :: iz, ipoin, nz, ncomp, irestart, iloop, ifnp, j, istep
  !  real(rp)    :: z_col(nz_nsa)
  !  real(rp)    :: rainnc(ncol_nsa), rainncv(ncol_nsa)
  real(rp)    :: z

  if (nzone > 1) call runend("NSA_APPLY_PHYSICS: THIS SUB IS NOT PREPARED TO RUN WITH ZONES.")

  irestart = vtkrestart_nsa

  !Local variables:
  nz = nz_nsa
  ncomp=min(3_ip,ncomp_nsa)
  !
  ! Store local working arrays (q, q_ref):
  !
  if( INOTMASTER ) then
     if(ndime < 3) then
        !
        ! 2D
        ! 
        do ipoin=1,npoin

           !Reference/hydrostatic variables:
           qref_nsa(1,ipoin) = rekee_nsa(ndime+1,ipoin) !densi
           qref_nsa(2,ipoin) = bvess_nsa(1,ipoin,1)     !u-velo
           qref_nsa(3,ipoin) = bvess_nsa(ndime,ipoin,1) !v-velo
           qref_nsa(4,ipoin) = rekee_nsa(ndime+2,ipoin) !theta

           qref_nsa(5,ipoin) = bvess_nsa(ndofn_nsa+1,ipoin,1) !q_vapor
           qref_nsa(6,ipoin) = bvess_nsa(ndofn_nsa+2,ipoin,1) !q_cloud
           qref_nsa(7,ipoin) = bvess_nsa(ndofn_nsa+3,ipoin,1) !q_rain

           !Solution variables (for us they are total so that we need to subtract the ref):
           !so that, to use them in kessler_col, we need to subtract
           !the reference values:
           q_nsa(1,ipoin) = densi(ipoin,ncomp)       - qref_nsa(1,ipoin)
           q_nsa(2,ipoin) = veloc(1,ipoin,ncomp)     - qref_nsa(2,ipoin)
           q_nsa(3,ipoin) = veloc(ndime,ipoin,ncomp) - qref_nsa(3,ipoin)
           q_nsa(4,ipoin) = tempe(ipoin,ncomp)       - qref_nsa(4,ipoin)

           q_nsa(5,ipoin) = conce(ipoin,1,ncomp)
           q_nsa(6,ipoin) = conce(ipoin,2,ncomp)
           q_nsa(7,ipoin) = conce(ipoin,3,ncomp)

           q_nsa(8,ipoin) = press(ipoin,ncomp) ! + press_ref(ipoin) 
        end do

     else
        !
        ! 3D
        !
        do ipoin=1,npoin

           !Reference/hydrostatic variables:
           qref_nsa(1,ipoin)  = rekee_nsa(ndime+1,ipoin)        !densi  ok
           qref_nsa(2,ipoin)  = bvess_nsa(1,ipoin,1)            !u-velo ok
           qref_nsa(3,ipoin)  = bvess_nsa(ndime-1,ipoin,1)      !v-velo ok
           qref_nsa(4,ipoin)  = bvess_nsa(ndime,ipoin,1)        !w-velo ok
           qref_nsa(5,ipoin)  = rekee_nsa(ndime+2,ipoin)        !theta  ok

           qref_nsa(6, ipoin) = bvess_nsa(ndofn_nsa+1,ipoin,1)  !qv ok
           qref_nsa(7, ipoin) = bvess_nsa(ndofn_nsa+2,ipoin,1)  !qc ok
           qref_nsa(8, ipoin) = bvess_nsa(ndofn_nsa+3,ipoin,1)  !qr ok

           !Solution variables (for us they are total so that we need to subtract the ref):
           !so that, to use them in kessler_col, we need to subtract
           !the reference values:
           q_nsa(1,ipoin) = densi(ipoin,ncomp)         - qref_nsa(1,ipoin)
           q_nsa(2,ipoin) = veloc(1,ipoin,ncomp)       - qref_nsa(2,ipoin)
           q_nsa(3,ipoin) = veloc(ndime-1,ipoin,ncomp) - qref_nsa(3,ipoin)
           q_nsa(4,ipoin) = veloc(ndime,ipoin,ncomp)   - qref_nsa(4,ipoin)
           q_nsa(5,ipoin) = tempe(ipoin,ncomp)         - qref_nsa(5,ipoin)

           q_nsa(6,ipoin) = conce(ipoin,1,ncomp)
           q_nsa(7,ipoin) = conce(ipoin,2,ncomp)
           q_nsa(8,ipoin) = conce(ipoin,3,ncomp)

           q_nsa(9,ipoin) = press(ipoin,ncomp) ! + press_ref(ipoin)

        end do
     end if
  end if
  !
  ! Do Kessler Physics on each column of data
  !
  call nsa_elmcol

  !
  ! Put data
  !
  if( INOTMASTER ) then
     if(ndime < 3) then
        !
        ! 2D
        !
        do ipoin = 1,npoin

           !Restore the total solution variables to continue in alya: 
           densi(ipoin,ncomp)   = q_nsa(1,ipoin) + qref_nsa(1,ipoin)
           tempe(ipoin,ncomp)   = q_nsa(4,ipoin) + qref_nsa(4,ipoin)
           conce(ipoin,1,ncomp) = q_nsa(5,ipoin)
           conce(ipoin,2,ncomp) = q_nsa(6,ipoin)
           conce(ipoin,3,ncomp) = q_nsa(7,ipoin)

        end do
     else
        !
        ! 3D
        !
        do ipoin = 1,npoin

           !Restore the total solution variables to continue in alya: 
           densi(ipoin,ncomp)   = q_nsa(1,ipoin) + qref_nsa(1,ipoin)
           tempe(ipoin,ncomp)   = q_nsa(5,ipoin) + qref_nsa(5,ipoin)
           conce(ipoin,1,ncomp) = q_nsa(6,ipoin)
           conce(ipoin,2,ncomp) = q_nsa(7,ipoin)
           conce(ipoin,3,ncomp) = q_nsa(8,ipoin)

        end do
     end if
  end if

  !-------------------------------------------------------------------------
  !Print the max of qv,qc,qr on screen
  !-------------------------------------------------------------------------
!!$  print*,''
!!$  print *, "nsa_kessler: Vapor qv: ", maxval(q_nsa(ndime+3,:)), minval(q_nsa(ndime+3,:))
!!$  !if (maxval(q_nsa(ndime+4,:)) > 0.0_rp .or. minval(q_nsa(ndime+4,:)) > 0.0_rp) then
!!$  print *, "nsa_kessler: Cloud qc: ", maxval(q_nsa(ndime+4,:)), minval(q_nsa(ndime+4,:))
!!$  !end if
!!$  !if (maxval(q_nsa(ndime+5,:)) > 0.0_rp .or. minval(q_nsa(ndime+5,:)) > 0.0_rp)then
!!$  print *, "nsa_kessler: Rain qr: ", maxval(q_nsa(ndime+5,:)), minval(q_nsa(ndime+5,:))
!!$  !end if
!!$  !-------------------------------------------------------------------------
end subroutine nsa_apply_physics

subroutine nsa_elmcol
!subroutine nsa_elmcol(ipoin,z)
!subroutine nsa_elmcol(ipoin,rainnc,rainncv,z)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsa_elmcol
  ! NAME 
  !    nsa_elmcol
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    nsa_apply_physics
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_elmtyp
  use def_domain
  use def_nastal
  use mod_memory
  implicit none
  real(rp)             :: elmat(mnode,mnode)
  real(rp)             :: elrhs(mnode)
  real(rp)             :: eladv(mnode)
  real(rp)             :: elpro(mnode)
  real(rp)             :: elcod(ndime,mnode)
  real(rp)             :: elden(mnode)
  real(rp)             :: gpcar(ndime,mnode,mgaus)      
  real(rp)             :: gphes(ntens,mnode,mgaus)      
  real(rp)             :: gpvol(mgaus)                  
  real(rp)             :: tragl(9),chale(2),hleng(3)  
  integer(ip)          :: ielem,pnode,pgaus,idime,inode,ipoin
  integer(ip)          :: pelty,pmate,plapl,porde,ptopo

  integer(ip)          :: ifiel,i,j,k
  integer(ip)          :: nfall
  integer(ip)          :: n, nz
  integer(ip)          :: nfall_new
  
  real(rp),    pointer :: vt(:)
  real(rp),    pointer :: vtden(:)
  real(rp),    pointer :: prod(:)
  real(rp),    pointer :: prodk(:)
  real(rp),    pointer :: rho(:)
  real(rp),    pointer :: rhok(:)
  real(rp),    pointer :: t(:)
  real(rp),    pointer :: qv(:)
  real(rp),    pointer :: qc(:)
  real(rp),    pointer :: qr(:)

  real(rp)             :: dtfall
  real(rp)             :: ppt
  real(rp)             :: crmax
  real(rp)             :: qrr
  real(rp)             :: time_sediment
  real(rp)             :: rhok1
  real(rp)             :: xlv
  real(rp)             :: ep2
  real(rp)             :: svp1
  real(rp)             :: svp2
  real(rp)             :: svp3
  real(rp)             :: svpt0
  real(rp)             :: rhowater
  real(rp)             :: pii
  real(rp)             :: z

  real(rp)             :: f5, rdz, product
  real(rp)             :: max_heating, max_condense, max_rain, maxqrp
  real(rp)             :: vtmax, ernmax, factorn
  real(rp)             :: qcr
  real(rp)             :: factorr
  real(rp)             :: qrprod
  real(rp)             :: ern
  real(rp)             :: gam
  real(rp)             :: rcgs
  real(rp)             :: rcgsi
  real(rp)             :: pressure
  real(rp)             :: temp
  real(rp)             :: es
  real(rp)             :: qvs
  real(rp)             :: dz
  
  integer(ip)          :: icolo,icolo_glo,icolo_max,iposi_max,ncolo
  integer(ip)          :: iposi
  integer(ip), pointer :: ifirs(:,:)
  integer(ip), pointer :: kfl_fixpr_nsa(:)
  real(rp),    pointer :: xfirs(:)

  !----------------------------------------------------------------
  ! parameters from the original module_mp_kessler.F
  !----------------------------------------------------------------
  real(rp), parameter  :: c1                   = 0.001_rp
  real(rp), parameter  :: c2                   = 0.001_rp
  real(rp), parameter  :: c3                   = 2.2_rp
  real(rp), parameter  :: c4                   = 0.875_rp
  real(rp), parameter  :: fudge                = 1.0_rp
  real(rp), parameter  :: mxfall               = 10.0_rp
  real(rp), parameter  :: max_cr_sedimentation = 0.75_rp
  !----------------------------------------------------------------
   
  !
  !  input values for squall-line test
  !
  xlv      = 2500000.0_rp        !Latent heat of vaporization
  ep2      = 0.6217504_rp
  svp1     = 0.6112000_rp   
  svp2     = 17.67000_rp
  svp3     = 29.65000_rp
  svpt0    = 273.1500_rp
  rhowater = 1000.000_rp  !density of rainwater
  
  !
  !  input arrays
  !
  !  t - potential temperature
  !  qv, qc, qr  - mixing ratio (g/g dry air) of water vapor, cloud water
  !                and rain water
  !  pii         - exner function
  !  dt_in - timestep
  !  z  - height of (t,qv,p,rho) points in meters
  !  dz8w - delta z between t points.
  !
  !  See Klemp and Wilhelmson (1978) Journal of the Atmospheric Sciences
  !  Vol 35, pp 1070-1096 for more details
  !
  f5     =  svp2*(svpt0-svp3)*xlv/cpcoe_nsa
  ernmax =    0.0_rp
  maxqrp = -100.0_rp
  
  !--------------------------------------------------------------------------
  ! parameters for the time split terminal advection
  !--------------------------------------------------------------------------
  max_heating  = 0.0_rp
  max_condense = 0.0_rp
  max_rain     = 0.0_rp
  
  if( INOTMASTER ) then
     call memory_alloca(mem_modul(1:2,modul),'VT','nsa_elmcol',vt,npoin)   
     call memory_alloca(mem_modul(1:2,modul),'VTDEN','nsa_elmcol',vtden,npoin)   
     call memory_alloca(mem_modul(1:2,modul),'PROD','nsa_elmcol',prod,npoin)    
     !call memory_alloca(mem_modul(1:2,modul),'PRODK','nsa_elmcol',prodk,npoin)
     call memory_alloca(mem_modul(1:2,modul),'RHO','nsa_elmcol',rho,npoin)    
     call memory_alloca(mem_modul(1:2,modul),'RHOK','nsa_elmcol',rhok,npoin)
     call memory_alloca(mem_modul(1:2,modul),'T','nsa_elmcol', t,npoin)     
     call memory_alloca(mem_modul(1:2,modul),'QV','nsa_elmcol',qv,npoin)      
     call memory_alloca(mem_modul(1:2,modul),'QC','nsa_elmcol',qc,npoin)     
     call memory_alloca(mem_modul(1:2,modul),'QR','nsa_elmcol',qr,npoin)  
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXPR_NSA','nsa_elmcol',kfl_fixpr_nsa,npoin)  
  end if

  ! are all the variables ready for the microphysics?
  ! start the microphysics
  ! do the sedimentation first
  crmax = 0.0_rp
  
  !----------------------------------------------------------------------------
  ! Terminal velocity calculation and advection, set up coefficients and
  ! compute stable timestep
  !----------------------------------------------------------------------------
    
  ifiel = 5 ! Change this with the number of the field set in .dom.dat. 
            ! E.g., if the COLUMN field is set as FIELD=5, set ifiel=5 here.
  !
  ! unkno in the solver is the PRODK of the original nsa_kessler.f90
  ! UNKNO <= PRODK
  !
  ! copy into 2D arrays (convert from element-based notation to one 2D array)
  if( INOTMASTER ) then
     if(ndime < 3) then
        !
        ! 2D, OK
        !
        do ipoin=1,npoin
           rho(ipoin) = q_nsa(1,ipoin) + qref_nsa(1,ipoin) !total density ok
           t(ipoin)   = q_nsa(4,ipoin) + qref_nsa(4,ipoin) !total theta   ok
           qv(ipoin)  = q_nsa(5,ipoin)                 !total vapor   ok
           qc(ipoin)  = q_nsa(6,ipoin)                 !total cloud   ok
           qr(ipoin)  = q_nsa(7,ipoin)                 !total rain    ok
        end do

     else
        !
        ! 3D
        !
        do ipoin=1,npoin
           rho(ipoin) = q_nsa(1,ipoin) + qref_nsa(1,ipoin) !total density
           t(ipoin)   = q_nsa(5,ipoin) + qref_nsa(5,ipoin) !total theta
           !qv(ipoin) = q_nsa(6,ipoin) + qref_nsa(6,ipoin)
           qv(ipoin)  = q_nsa(6,ipoin)                 !total vapor
           qc(ipoin)  = q_nsa(7,ipoin)                 !total cloud
           qr(ipoin)  = q_nsa(8,ipoin)                 !total rain

        end do !ipoin
     end if
  end if

  !
  ! Prepare data structures for column and global numbering.
  ! It is read from the FIELDS file '_COLUMN.dat' in *.dom.dat
  !
  if( INOTMASTER ) then
     call memgen(1_ip,npoin,0_ip)
     icolo_max = 0
     iposi_max = 0
     do ipoin = 1,npoin
        icolo_glo = int(xfiel(ifiel) % a(1,ipoin,1),ip)
        iposi     = int(xfiel(ifiel) % a(2,ipoin,1),ip)
        icolo_max = max(icolo_max,icolo_glo)
        iposi_max = max(iposi_max,iposi)
     end do
  end if
  call parari('MAX',0_ip,1_ip,icolo_max)
  call parari('MAX',0_ip,1_ip,iposi_max)

  allocate( xfirs(icolo_max) )
  if( INOTMASTER ) then
     allocate( ifirs(2,icolo_max) )
     icolo = 0
     do ipoin = 1,npoin
        icolo_glo = int(xfiel(ifiel) % a(1,ipoin,1),ip)
        iposi     = int(xfiel(ifiel) % a(2,ipoin,1),ip)
        if( iposi == 1 .and. ( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) ) then
           icolo          = icolo + 1
           ifirs(1,icolo) = icolo_glo
           ifirs(2,icolo) = ipoin
        end if
     end do
     ncolo = icolo !global number of columns         (SM)
     nz    = iposi !global number of vertical levels (SM)         
  end if

  do icolo = 1,icolo_max
     xfirs(icolo) = 0.0_rp
  end do
  
  if( INOTMASTER ) then
     do icolo = 1,ncolo
        icolo_glo        = ifirs(1,icolo)
        ipoin            = ifirs(2,icolo)
        xfirs(icolo_glo) = rho(ipoin)             !SACA rho in every subdomain. ESTA BIEN ESTO?
     end do
  end if
  call pararr('SUM',0_ip,icolo_max,xfirs)
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        icolo_glo    = int(xfiel(ifiel) % a(1,ipoin,1),ip)
        rhok1        = xfirs(icolo_glo)                           !rk      ok

        prodk_nsa(ipoin) = qr(ipoin)
        rhok(ipoin)  = rho(ipoin)

        vtden(ipoin) = sqrt( rhok1 / rhok(ipoin) )                !vtden   ok
        qrr = max(0.0_rp,qr(ipoin)*0.001_rp*rhok(ipoin))          !qrr     ok
        vt(ipoin) = 36.34*(qrr**0.1364_rp)*vtden(ipoin)        !vt      ok
     end do
  end if
  !
  ! Fix boundary condition
  !
  if( INOTMASTER ) then
    do ipoin = 1,npoin
        iposi     = int(xfiel(ifiel) % a(2,ipoin,1),ip)
        if( iposi_max == iposi ) kfl_fixpr_nsa(ipoin) = 1
     end do
  end if
  !
  ! Time step: CFL and dtfall:
  !
  if( INOTMASTER ) then
     do ielem = 1,nelem
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        porde = lorde(pelty)
        ptopo = ltopo(pelty)
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
           end do
        end do
        call elmlen(&
             ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
             hnatu(pelty),hleng)
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           crmax = max( vt(ipoin)*dtime/hleng(ndime),crmax )
        end do
     end do
  end if
  
  call pararr('MAX',0_ip,1_ip,crmax)
  ! OJO
  !nfall  = max(1_ip,nint(0.5_rp + crmax / max_cr_sedimentation,ip )) !nfall
  nfall  = max(1_ip,int(0.5_rp + crmax / max_cr_sedimentation,ip )) !nfall
  dtfall = dtime / real(nfall,rp)                                 !dtfall
  time_sediment = dtime /10                                       !time_sedimentation
 
  !----------------------------------------------------------------------
  !
  ! Assemble matrix and RHS
  !
  ! Terminal velocity calculation and advection
  ! Do a time split loop on this for stability.
  !----------------------------------------------------------------------
  column_sedimentation: do while( nfall > 0 )
     
     time_sediment = time_sediment - dtfall
     
     ! SM in the original kessler factor is computed, but we don't and I can't see how
     ! the equation is written (rhsid) SM
     !     do k = 1, kte-1
     !       factor(k) = dtfall*rdzk(k)/rhok(k)
     !    enddo
     !    factor(kte) = dtfall*rdzk(kte)
     if(INOTMASTER) then
        do ipoin = 1,npoin
           rhsid(ipoin) = 0.0_rp
        end do
     end if
     
     !------------------------------------------------------------------------------
     ! Time split loop, Fallout done with flux upstream
     !------------------------------------------------------------------------------
     if(INOTMASTER) then
        elements: do ielem=1,nelem
           !
           ! Element properties and dimensions
           !
           pelty = ltype(ielem)
           
           if( pelty > 0 ) then
              
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              plapl = 0
              porde = lorde(pelty)
              ptopo = ltopo(pelty)
              pmate = 1
              !
              ! Gather operations: ELCOD
              !
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 elden(inode) =  rhok(ipoin)
                 elpro(inode) =  prodk_nsa(ipoin)
                 eladv(inode) = -vt(ipoin)
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do
              !
              ! Cartesian derivatives, Hessian matrix and volume: GPCAR, PGVOL
              !
              call elmcar(&
                   pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                   elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
                   gphes,ielem)
              !
              ! HLENG and TRAGL at center of gravity and CHALE
              !
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                   hnatu(pelty),hleng)
              !
              ! Assemble elemental matrix
              !
              call nsa_elmco1(&
                   pgaus,pnode,lnods(1,ielem),gpvol,gpcar,&
                   elmar(pelty) % shape,chale,hleng,elden,&
                   eladv,elpro,elmat,elrhs,ielem)
              !
              ! Assembly: AMATR and RHSID
              !
              call assrhs(&
                   1_ip,pnode,lnods(1,ielem),elrhs,rhsid)
           end if
           
        end do elements
     end if
     
     if(INOTMASTER) then
        !----------------------------------------------------------------------
        !
        ! Solve system
        !
        !----------------------------------------------------------------------

        call rhsmod(1_ip,rhsid) 
        !
        ! Boundary condition
        !
        do ipoin = 1,npoin
           if( kfl_fixpr_nsa(ipoin) == 1 ) rhsid(ipoin) = 0.0_rp
        end do
        !SM : no creo que lo esta resolviendo bien. Donde esta el equivalente de factor(k) auqi?
        !unkno here is prodk_nsa() in the original column-wise code:
        do ipoin = 1,npoin
           prodk_nsa(ipoin) = prodk_nsa(ipoin) + dtfall / ( rhok(ipoin) * vmass(ipoin) ) * rhsid(ipoin)
        end do

!        prodk_nsa = unkno
     end if
     !SM donde esta el quivalente de 
     !k = kte
     !prodk_nsa(k) = prodk_nsa(k) - factor(k)*prodk_nsa(k)*vt(k)
     
     !------------------------------------------------------------------------------
     !
     ! compute new sedimentation velocity, and check/recompute new
     ! sedimentation timestep if this isn't the last split step.
     !
     !------------------------------------------------------------------------------

     if( nfall > 1_ip ) then! this wasn't the last split sedimentation timestep

        nfall = nfall - 1_ip
        crmax = 0.0_rp
        if( INOTMASTER ) then
           do ipoin = 1,npoin
              qrr       = max(0.0_rp,prodk_nsa(ipoin)*0.001_rp*rhok(ipoin)) !qrr
              vt(ipoin) = 36.34_rp*(qrr**0.1364_rp)*vtden(ipoin)        !vt
           end do
           !
           ! Time step: CFL and dtfall:
           !
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              porde = lorde(pelty)
              ptopo = ltopo(pelty)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                   hnatu(pelty),hleng)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 crmax = max( vt(ipoin)*time_sediment/hleng(ndime),crmax )
              end do
           end do
        end if !INOTMASTER
        
        call pararr('MAX',0_ip,1_ip,crmax)
        nfall_new = max(1_ip,int(0.5_rp+crmax/max_cr_sedimentation,ip))
        !nfall_new = max(1_ip,nint(0.5_rp+crmax/max_cr_sedimentation))
        if (nfall_new /= nfall ) then
           nfall  = nfall_new
           dtfall = time_sediment/real(nfall,rp)
        end if

     else

        if( INOTMASTER ) then
           do ipoin = 1,npoin
              !prod(ipoin) = unkno(ipoin)
              prod(ipoin) = prodk_nsa(ipoin)
           end do
        end if
        nfall = 0_ip  ! exit condition for sedimentation loop
     end if

  end do column_sedimentation
  
  ! now the conversion processes
  !----------------------------------------------------------------------------
  ! Production of rain and deletion of qc
  ! Production of qc from supersaturation
  ! Evaporation of QR
  !----------------------------------------------------------------------------

  if(INOTMASTER) then
     do ipoin = 1,npoin
        factorn = 1.0_rp / (1.0_rp+c3*dtime*max(0.0_rp,qr(ipoin))**c4)

        qrprod = qc(ipoin) * (1.0_rp - factorn) &
             + factorn*c1*dtime*max(qc(ipoin)-c2,0.0_rp)

        rcgs = 0.001_rp*rho(ipoin)

        qc(ipoin) = max(qc(ipoin) - qrprod, 0.0_rp)
        qr(ipoin) = (qr(ipoin) + prod(ipoin) - qr(ipoin))
        qr(ipoin) = max(qr(ipoin) + qrprod, 0.0_rp)

        pii      = (q_nsa(ndime+6,ipoin)/1.e5)**(287.0_rp/1004.0_rp) !ok
        temp     = t(ipoin)*pii                                      !ok
        pressure = q_nsa(ndime+6,ipoin)                              !ok

        gam = 2.5e+06/(1004.0_rp*pii)
        !      qvs       = 380.*exp(17.27*(temp-273.)/(temp- 36.))/pressure
        es        = 1000.0_rp*svp1*exp(svp2*(temp-svpt0)/(temp-svp3))
        qvs       = ep2*es/(pressure-es)
        !      prod(i,k) = (qv(i,k)-qvs) / (1.+qvs*f5/(temp-36.)**2)
        prod(ipoin) = (qv(ipoin)-qvs) / (1.0_rp+pressure/(pressure-es)*qvs*f5/ &
             (temp-svp3)**2)

        ern  = min(dtime * (((1.6_rp + 124.9_rp*(rcgs*qr(ipoin))**0.2046_rp) &
             * (rcgs*qr(ipoin))**0.525_rp)/(2.55e8/(pressure*qvs) &      
             +5.4e5))*(dim(qvs,qv(ipoin))/(rcgs*qvs)),  &           
             max(-prod(ipoin)-qc(ipoin), 0.0_rp),qr(ipoin))
       
        !
        ! Update all variables
        !
        product = max(prod(ipoin),-qc(ipoin))
        ! if(product > 1e-8 ) then
        !    PRINT*,'PRODUCTION', product, t(ipoin)
        ! end if
        t (ipoin) = t(ipoin) + gam*(product - ern)
        qv(ipoin) = max(qv(ipoin) - product + ern,0.0_rp)
        qc(ipoin) = qc(ipoin) + product
        qr(ipoin) = qr(ipoin) - ern

     end do

     !
     ! transform modified variables back to the element-based notation
     ! make sure only perturbations get modified, not the base state
     !
     if(ndime < 3) then
        !
        ! 2D
        !
        do ipoin=1,npoin   
           q_nsa(1,ipoin) = rho(ipoin) - qref_nsa(1,ipoin)
           q_nsa(4,ipoin) = t(ipoin)   - qref_nsa(4,ipoin)
           !q_nsa(5,ipoin) = qv(ipoin)  - qref_nsa(5,ipoin)
           q_nsa(5,ipoin) = qv(ipoin)  
           q_nsa(6,ipoin) = qc(ipoin)
           q_nsa(7,ipoin) = qr(ipoin)

           !qr(ipoin) = unkno(ipoin)
        end do !ipoin
     else
        !
        ! 3D
        !
        do ipoin=1,npoin
           q_nsa(1,ipoin) = rho(ipoin) - qref_nsa(1,ipoin)
           q_nsa(5,ipoin) = t(ipoin)   - qref_nsa(5,ipoin)
           !q_nsa(6,ipoin) = qv(ipoin)  - qref_nsa(6,ipoin)
           q_nsa(6,ipoin) = qv(ipoin)  
           q_nsa(7,ipoin) = qc(ipoin)
           q_nsa(8,ipoin) = qr(ipoin)

           !qr(ipoin) = unkno(ipoin)
        end do !ipoin
     end if
     
  end if !INOTMASTER

  !----------------------------------------------------------------------------
  ! Deallocate before exiting.
  !----------------------------------------------------------------------------
  
  if( INOTMASTER ) then   
     call memory_deallo(mem_modul(1:2,modul),'VT','nsa_elmcol',vt ) 
     call memory_deallo(mem_modul(1:2,modul),'VTDEN','nsa_elmcol',vtden )   
     call memory_deallo(mem_modul(1:2,modul),'PROD','nsa_elmcol', prod  )   
     !call memory_deallo(mem_modul(1:2,modul),'PRODK','nsa_elmcol', prodk)    
     call memory_deallo(mem_modul(1:2,modul),'RHO','nsa_elmcol',  rho   )   
     call memory_deallo(mem_modul(1:2,modul),'RHOK','nsa_elmcol', rhok  )
     call memory_deallo(mem_modul(1:2,modul),'T','nsa_elmcol',    t     )     
     call memory_deallo(mem_modul(1:2,modul),'QV','nsa_elmcol',   qv    )     
     call memory_deallo(mem_modul(1:2,modul),'QC','nsa_elmcol',   qc    )    
     call memory_deallo(mem_modul(1:2,modul),'QR','nsa_elmcol',   qr    ) 
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXPR_NSA','nsa_elmcol',kfl_fixpr_nsa) 
  end if
  deallocate( xfirs )
  if( INOTMASTER ) then
     deallocate( ifirs )
  end if

end subroutine nsa_elmcol

subroutine nsa_elmco1(&
     pgaus,pnode,lnods,gpvol,gpcar,gpsha,chale,&
     hleng,elden,eladv,elpro,elmat,elrhs,ielem)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsa_elmhy1
  ! NAME 
  !    nsa_elmhy1
  ! DESCRIPTION
  !    Assemble the elemental matrix from Gauss point contributions
  !   
  ! USES
  ! USED BY
  !    nsa_elmhyd
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_elmtyp, only     :  ELEXT
  use def_domain, only     :  mnode,ndime,coord,lpoty,lelch
  implicit none
  integer(ip), intent(in)  :: pgaus,pnode,ielem
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpvol(pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: chale(2)
  real(rp),    intent(in)  :: hleng(ndime)
  real(rp),    intent(in)  :: elden(pnode)
  real(rp),    intent(in)  :: eladv(pnode)
  real(rp),    intent(in)  :: elpro(pnode)
  real(rp),    intent(out) :: elmat(pnode,pnode)
  real(rp),    intent(out) :: elrhs(pnode)
  integer(ip)              :: igaus,inode,jnode,ipoin,ibopo,idime,dummi
  real(rp)                 :: fact1,fact2,tau,ggrap(mnode),vtest,dummr(2)
  real(rp)                 :: gpadv(3),gpden,gpunk
  real(rp)                 :: gradv(3,3),grden(3),grunk(3),bvess
  
  !----------------------------------------------------------------------
  !
  ! Compute LHS and RHS
  !
  !----------------------------------------------------------------------


  do inode = 1,pnode  
     elrhs(inode) = 0.0_rp
     do jnode = 1,pnode
        elmat(jnode,inode) = 0.0_rp
     end do
  end do

  do igaus = 1,pgaus
     
     gpadv = 0.0_rp
     gpden = 0.0_rp
     gpunk = 0.0_rp
     gradv = 0.0_rp
     grden = 0.0_rp
     grunk = 0.0_rp
     ggrap = 0.0_rp

     do inode = 1,pnode
        gpadv(ndime) = gpadv(ndime) + eladv(inode) * gpsha(inode,igaus) 
        gpden        = gpden        + elden(inode) * gpsha(inode,igaus) 
        gpunk        = gpunk        + elpro(inode) * gpsha(inode,igaus) 
        do idime = 1,ndime
           ggrap(inode)       = ggrap(inode)       + gpadv(idime) * gpcar(idime,inode,igaus)
           gradv(idime,ndime) = gradv(idime,ndime) + eladv(inode) * gpcar(idime,inode,igaus)
           grden(idime)       = grden(idime)       + elden(inode) * gpcar(idime,inode,igaus)
           grunk(idime)       = grunk(idime)       + elpro(inode) * gpcar(idime,inode,igaus)
        end do
     end do

     tau = 0.0_rp
     if( abs(gpadv(ndime)) > 1.0e-12_rp ) tau = hleng(ndime) / abs(gpadv(ndime))

     do inode = 1,pnode 
        vtest = ( gpsha(inode,igaus) + tau * ggrap(inode) ) * gpvol(igaus)
        do idime = 1,ndime
           elrhs(inode) = elrhs(inode) - vtest              &
                & * (   gpden * gpadv(idime) * grunk(idime) &
                &     + gpden * gpunk * gradv(idime,idime)  &
                &     + gpunk * gpadv(idime) * grden(idime) )
        end do
     end do

  end do

  !----------------------------------------------------------------------
  !
  ! Extension elements
  !
  !----------------------------------------------------------------------

  if( lelch(ielem) == ELEXT ) then
     call elmext(&
          4_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
          dummr,elrhs,dummr)
  end if

end subroutine nsa_elmco1
