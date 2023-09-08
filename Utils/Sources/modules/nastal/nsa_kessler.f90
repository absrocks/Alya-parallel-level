!!$!!NEW VERSION FOR PARALLELIZATION:
!!$!
!!$! 3D:
!!$!
!!$!-------------------------------------------------------------------------------------
!!$!From Jim original, Simone, Dec 22, 2011
!!$! Performs Kessler Physics on a col. of data
!!$! Iq(nputs:
!!$! q_col = col of prognostic variables (rho', u, w, theta', qv', qc, qr)
!!$! qref_col = col. of reference variables (rhoref, 0, 0, theteref, qvref, 0, 0)
!!$! p = col of total pressure ( p = p' + pref)
!!$! nvar = number of variables = 7
!!$! Outputs
!!$! rainnc = accumulation of rain at the surface
!!$! z = col. of z-values
!!$! nz = number of z-values in the col.
!!$! dt = time-step
!!$! James F. Kelly
!!$! 29 January 2010
!!$!-------------------------------------------------------------------------------------
!!$
!!$!------------------------------------------------------------
!!$! subroutine apply_physics
!!$!------------------------------------------------------------
!!$subroutine apply_physics(dt,rainnc,rainncv)
!!$
!!$  use      def_master
!!$  use      def_domain
!!$  use      mod_postpr
!!$  use      def_nastal
!!$
!!$  implicit none
!!$
!!$!!  real(rp)   ::  q_nsa(nvar_nsa,npoin),  q_ref_nsa(nvar_nsa,npoin)
!!$!!  real(rp)   ::  q_col_nsa(nvar_nsa,nz_nsa), qref_col_nsa(nvar_nsa,nz_nsa)
!!$!!  real(rp)   ::  press_col_nsa(nz_nsa)
!!$
!!$  integer(ip)::  icol, iz, ipoin, nz, ncomp, irestart, iloop, ifnp, j, istep
!!$!!  real(rp)   ::  z_col_nsa(nz_nsa)
!!$  !real(rp)   ::  press(npoin)
!!$  real(rp)   ::  dt
!!$!!  real(rp)   ::  rainnc_nsa(ncol_nsa), rainncv_nsa(ncol_nsa)
!!$  real(rp)   ::  rainnc(ncol_nsa), rainncv(ncol_nsa)
!!$
!!$  character               :: fnp1*24, fnp*72
!!$
!!$  irestart = vtkrestart_nsa
!!$
!!$  !Local variables:
!!$  nz = nz_nsa
!!$  ncomp=min(3,ncomp_nsa)
!!$  !
!!$  ! Store local working arrays (q, q_ref):
!!$  !
!!$  if(ndime < 3) then
!!$     !
!!$     ! 2D
!!$     ! 
!!$     do ipoin=1,npoin
!!$        
!!$        !Reference/hydrostatic variables:
!!$        q_ref_nsa(1,ipoin) = rekee_nsa(ndime+1,ipoin) !densi
!!$        q_ref_nsa(2,ipoin) = bvess_nsa(1,ipoin,1)     !u-velo
!!$        q_ref_nsa(3,ipoin) = bvess_nsa(ndime,ipoin,1) !v-velo
!!$        q_ref_nsa(4,ipoin) = rekee_nsa(ndime+2,ipoin) !theta
!!$    
!!$        q_ref_nsa(5,ipoin) = bvess_nsa(ndofn_nsa+1,ipoin,1)       !q_vapor
!!$        q_ref_nsa(6,ipoin) = bvess_nsa(ndofn_nsa+2,ipoin,1)       !q_cloud
!!$        q_ref_nsa(7,ipoin) = bvess_nsa(ndofn_nsa+3,ipoin,1)       !q_rain
!!$        
!!$        !Solution variables (for us they are total so that we need to subtract the ref):
!!$        !so that, to use them in kessler_col, we need to subtract
!!$        !the reference values:
!!$        q_nsa(1,ipoin) = densi(ipoin,ncomp)       - q_ref_nsa(1,ipoin)
!!$        q_nsa(2,ipoin) = veloc(1,ipoin,ncomp)     - q_ref_nsa(2,ipoin)
!!$        q_nsa(3,ipoin) = veloc(ndime,ipoin,ncomp) - q_ref_nsa(3,ipoin)
!!$        q_nsa(4,ipoin) = tempe(ipoin,ncomp)       - q_ref_nsa(4,ipoin)
!!$   
!!$        q_nsa(5,ipoin) = conce(ipoin,1,ncomp)
!!$        q_nsa(6,ipoin) = conce(ipoin,2,ncomp)
!!$        q_nsa(7,ipoin) = conce(ipoin,3,ncomp)
!!$        
!!$     end do
!!$  else
!!$     !
!!$     ! 3D
!!$     !
!!$     do ipoin=1,npoin
!!$
!!$        !Reference/hydrostatic variables:
!!$        q_ref_nsa(1,ipoin) = rekee_nsa(ndime+1,ipoin)        !densi  ok
!!$        q_ref_nsa(2,ipoin) = bvess_nsa(1,ipoin,1)            !u-velo ok
!!$        q_ref_nsa(3,ipoin) = bvess_nsa(ndime-1,ipoin,1)      !v-velo ok
!!$        q_ref_nsa(4,ipoin) = bvess_nsa(ndime,ipoin,1)        !w-velo ok
!!$        q_ref_nsa(5,ipoin) = rekee_nsa(ndime+2,ipoin)        !theta  ok
!!$
!!$        q_ref_nsa(6, ipoin) = bvess_nsa(ndofn_nsa+1,ipoin,1)  !qv ok
!!$        q_ref_nsa(7, ipoin) = bvess_nsa(ndofn_nsa+2,ipoin,1)  !qc ok
!!$        q_ref_nsa(8, ipoin) = bvess_nsa(ndofn_nsa+3,ipoin,1)  !qr ok
!!$
!!$!        q_ref_nsa(6, ipoin) = rekee_nsa(ndime+3 + 1, ipoin)  !qv ok testing
!!$!        q_ref_nsa(7, ipoin) = rekee_nsa(ndime+3 + 2, ipoin)  !qc ok testing
!!$!        q_ref_nsa(8, ipoin) = rekee_nsa(ndime+3 + 3, ipoin)  !qr ok testing
!!$
!!$        !Solution variables (for us they are total so that we need to subtract the ref):
!!$        !so that, to use them in kessler_col, we need to subtract
!!$        !the reference values:
!!$        q_nsa(1,ipoin) = densi(ipoin,ncomp)         - q_ref_nsa(1,ipoin)
!!$        q_nsa(2,ipoin) = veloc(1,ipoin,ncomp)       - q_ref_nsa(2,ipoin)
!!$        q_nsa(3,ipoin) = veloc(ndime-1,ipoin,ncomp) - q_ref_nsa(3,ipoin)
!!$        q_nsa(4,ipoin) = veloc(ndime,ipoin,ncomp)   - q_ref_nsa(4,ipoin)
!!$        q_nsa(5,ipoin) = tempe(ipoin,ncomp)         - q_ref_nsa(5,ipoin)
!!$       
!!$        q_nsa(6,ipoin) = conce(ipoin,1,ncomp)
!!$        q_nsa(7,ipoin) = conce(ipoin,2,ncomp)
!!$        q_nsa(8,ipoin) = conce(ipoin,3,ncomp)
!!$
!!$     end do
!!$  end if
!!$  
!!$  do icol = 1,ncol_nsa
!!$     ! Extract a column of data
!!$     if(ndime < 3) then
!!$        !
!!$        ! 2D
!!$        !
!!$        do iz = 1,nz_nsa
!!$           
!!$           ipoin = node_column_nsa(icol,iz)        !Mapping from column to global ipoin
!!$           
!!$           q_col_nsa(1,iz)    = q_nsa(1,ipoin)     !perturbed densi
!!$           q_col_nsa(2,iz)    = q_nsa(2,ipoin)     !u-velo
!!$           q_col_nsa(3,iz)    = q_nsa(3,ipoin)     !v-velo
!!$           q_col_nsa(4,iz)    = q_nsa(4,ipoin)     !perturbed theta
!!$           q_col_nsa(5,iz)    = q_nsa(5,ipoin)     !total vapor
!!$           q_col_nsa(6,iz)    = q_nsa(6,ipoin)     !total cloud
!!$           q_col_nsa(7,iz)    = q_nsa(7,ipoin)     !total rain
!!$
!!$           qref_col_nsa(1,iz) = q_ref_nsa(1,ipoin) !reference densi
!!$           qref_col_nsa(2,iz) = q_ref_nsa(2,ipoin) !reference u-velo
!!$           qref_col_nsa(3,iz) = q_ref_nsa(3,ipoin) !reference v-velo
!!$           qref_col_nsa(4,iz) = q_ref_nsa(4,ipoin) !reference theta
!!$           qref_col_nsa(5,iz) = q_ref_nsa(5,ipoin) !reference vapor
!!$           qref_col_nsa(6,iz) = q_ref_nsa(6,ipoin) !reference cloud
!!$           qref_col_nsa(7,iz) = q_ref_nsa(7,ipoin) !reference rain
!!$
!!$           press_col_nsa(iz)  = press(ipoin,ncomp) ! + press_ref(ipoin) 
!!$           !This muyst be TOTAL pressure.
!!$           !Add the hydrostatic value 
!!$           !if your press() is the perturbed pressure.
!!$           z_col_nsa(iz) = coord(ndime,ipoin)
!!$
!!$           !print*,'d', q_col_nsa(1,iz),q_col_nsa(4,iz), q_col_nsa(5,iz), press_col_nsa(iz)
!!$        end do
!!$     else
!!$        !
!!$        ! 3D
!!$        !
!!$        do iz = 1,nz_nsa
!!$           ipoin = node_column_nsa(icol,iz)
!!$           q_col_nsa(1,iz)    = q_nsa(1,ipoin)     !perturbed densi
!!$           q_col_nsa(2,iz)    = q_nsa(2,ipoin)     !u-velo
!!$           q_col_nsa(3,iz)    = q_nsa(3,ipoin)     !v-velo
!!$           q_col_nsa(4,iz)    = q_nsa(4,ipoin)     !w-velo
!!$           q_col_nsa(5,iz)    = q_nsa(5,ipoin)     !perturbed theta
!!$
!!$           q_col_nsa(6,iz)    = q_nsa(6,ipoin)     !total qv
!!$           q_col_nsa(7,iz)    = q_nsa(7,ipoin)     !total qc
!!$           q_col_nsa(8,iz)    = q_nsa(8,ipoin)     !total qr
!!$
!!$           qref_col_nsa(1,iz) = q_ref_nsa(1,ipoin) !reference densi
!!$           qref_col_nsa(2,iz) = q_ref_nsa(2,ipoin) !reference u-velo
!!$           qref_col_nsa(3,iz) = q_ref_nsa(3,ipoin) !reference v-velo
!!$           qref_col_nsa(4,iz) = q_ref_nsa(4,ipoin) !reference w-velo
!!$           qref_col_nsa(5,iz) = q_ref_nsa(5,ipoin) !reference theta
!!$
!!$           qref_col_nsa(6,iz) = q_ref_nsa(6,ipoin) !reference vapor
!!$           qref_col_nsa(7,iz) = q_ref_nsa(7,ipoin) !reference cloud
!!$           qref_col_nsa(8,iz) = q_ref_nsa(8,ipoin) !reference rain
!!$
!!$           press_col_nsa(iz)  = press(ipoin,ncomp) ! + press_ref(ipoin) 
!!$           !This must be TOTAL pressure.
!!$           !Add the hydrostatic value 
!!$           !if your press() is the perturbed pressure.
!!$           z_col_nsa(iz) = coord(ndime,ipoin)
!!$
!!$        end do
!!$     end if
!!$
!!$     !
!!$     ! Do Kessler Physics on each column of data
!!$     !
!!$     !print *, minval(press_col), maxval(press_col)
!!$     call kessler_col(icol,rainnc_nsa(icol),rainncv_nsa(icol),dt)
!!$     
!!$     !
!!$     ! Put data
!!$     !
!!$     if(ndime < 3) then
!!$        !
!!$        ! 2D
!!$        !
!!$        do iz = 1,nz
!!$
!!$           ipoin = node_column_nsa(icol,iz)
!!$           q_nsa(1,ipoin) = q_col_nsa(1,iz) !perturbed densi
!!$           q_nsa(4,ipoin) = q_col_nsa(4,iz) !perturbed theta
!!$           q_nsa(5,ipoin) = q_col_nsa(5,iz) !total qv
!!$           q_nsa(6,ipoin) = q_col_nsa(6,iz) !total qc
!!$           q_nsa(7,ipoin) = q_col_nsa(7,iz) !total qr
!!$
!!$           !Restore the total solution variables to continue in alya: 
!!$           densi(ipoin,ncomp)   = q_nsa(1,ipoin) + q_ref_nsa(1,ipoin)
!!$           tempe(ipoin,ncomp)   = q_nsa(4,ipoin) + q_ref_nsa(4,ipoin)
!!$           conce(ipoin,1,ncomp) = q_nsa(5,ipoin)
!!$           conce(ipoin,2,ncomp) = q_nsa(6,ipoin)
!!$           conce(ipoin,3,ncomp) = q_nsa(7,ipoin)
!!$
!!$        end do
!!$     else
!!$        !
!!$        ! 3D
!!$        !
!!$        do iz = 1,nz
!!$
!!$           ipoin = node_column_nsa(icol,iz)
!!$           q_nsa(1,ipoin) = q_col_nsa(1,iz) !perturbed densi
!!$           q_nsa(5,ipoin) = q_col_nsa(5,iz) !perturbed theta
!!$           q_nsa(6,ipoin) = q_col_nsa(6,iz) !total qv
!!$           q_nsa(7,ipoin) = q_col_nsa(7,iz) !total qc
!!$           q_nsa(8,ipoin) = q_col_nsa(8,iz) !total qr
!!$
!!$           !Restore the total solution variables to continue in alya: 
!!$           densi(ipoin,ncomp)   = q_nsa(1,ipoin) + q_ref_nsa(1,ipoin)
!!$           tempe(ipoin,ncomp)   = q_nsa(5,ipoin) + q_ref_nsa(5,ipoin)
!!$           conce(ipoin,1,ncomp) = q_nsa(6,ipoin)
!!$           conce(ipoin,2,ncomp) = q_nsa(7,ipoin)
!!$           conce(ipoin,3,ncomp) = q_nsa(8,ipoin)
!!$
!!$        end do
!!$     end if
!!$     
!!$  end do
!!$
!!$!  print*,''
!!$!  print *, "nsa_kessler: Vapor qv: ", maxval(q_nsa(ndime+3,:)), minval(q_nsa(ndime+3,:))
!!$!  if (maxval(q_nsa(ndime+4,:)) > 0.0_rp .or. minval(q_nsa(ndime+4,:)) > 0.0_rp) then
!!$!     print *, "nsa_kessler: Cloud qc: ", maxval(q_nsa(ndime+4,:)), minval(q_nsa(ndime+4,:))
!!$!  end if
!!$!  if (maxval(q_nsa(ndime+5,:)) > 0.0_rp .or. minval(q_nsa(ndime+5,:)) > 0.0_rp)then
!!$!     print *, "nsa_kessler: Rain qr: ", maxval(q_nsa(ndime+5,:)), minval(q_nsa(ndime+45,:))
!!$!  end if
!!$
!!$end subroutine apply_physics
!!$
!!$
!!$!!subroutine kessler_col(icol, q_col_nsa,qref_col_nsa, &
!!$!!     pcol_nsa,rainnc_nsa,rainncv_nsa,z_nsa,dt)
!!$subroutine kessler_col(icol,rainnc,rainncv,dt)
!!$
!!$  use      def_master
!!$  use      def_domain
!!$  use      mod_postpr
!!$  use      def_nastal
!!$  
!!$  implicit none
!!$
!!$  !global arrays
!!$  integer(ip)::  nz, nvar, icol
!!$!  real(rp)   ::  q_col_nsa(nvar_nsa,nz_nsa), qref_col_nsa(nvar_nsa,nz_nsa), 
!!$!                 z_nsa(nz_nsa)
!!$!  real(rp)   ::  pcol_nsa(nz_nsa)
!!$
!!$  !rain arrays
!!$  real(rp)   ::  rainnc, rainncv
!!$
!!$  ! local variables
!!$
!!$  integer(ip)  ::  i,j,k,ie
!!$  integer(ip)  ::  icount,jcount
!!$  integer(ip)  ::  its,ite,kts,kte
!!$
!!$  ! variables
!!$
!!$ !! real(rp), pointer :: rhocol_nsa(:),tcol_nsa(:),qvcol_nsa(:),qccol_nsa(:),qrcol_nsa(:)
!!$
!!$  ! local variables
!!$
!!$  real(rp)   ::  xlv,ep2,svp1,svp2,svp3,svpt0,rhowater, pii
!!$
!!$  ! local variables from the original module_mp_kessler.F
!!$
!!$  real(rp)   :: qrprod, ern, gam, rcgs, rcgsi
!!$!!  real(rp), pointer :: prodcol_nsa(:)
!!$!!  real(rp), pointer :: vtcol_nsa(:), prodkcol_nsa(:), vtdencol_nsa(:), rdzkcol_nsa(:), rhokcol_nsa(:), factorcol_nsa(:), rdzwcol_nsa(:)
!!$  integer(ip)  :: nfall, n, nfall_new
!!$  real(rp)     :: qrr, pressure, temp, es, qvs, dz, dt
!!$  real(rp)     :: f5, dtfall, rdz, product
!!$  real(rp)     :: max_heating, max_condense, max_rain, maxqrp
!!$  real(rp)     :: vtmax, ernmax, crmax, factorn, time_sediment
!!$  real(rp)     :: qcr, factorr, ppt
!!$
!!$  real(rp), parameter :: max_cr_sedimentation = 0.75_rp
!!$  
!!$  ! parameters from the original module_mp_kessler.F
!!$  !----------------------------------------------------------------
!!$  real(rp), parameter ::  c1     = 0.001_rp
!!$  real(rp), parameter ::  c2     = 0.001_rp
!!$  real(rp), parameter ::  c3     = 2.2_rp
!!$  real(rp), parameter ::  c4     = 0.875_rp
!!$  real(rp), parameter ::  fudge  = 1.0_rp
!!$  real(rp), parameter ::  mxfall = 10.0_rp
!!$  !----------------------------------------------------------------
!!$  
!!$  !Local constants:
!!$  nz   = nz_nsa
!!$  nvar = nvar_nsa
!!$
!!$  !
!!$  !  input values for squall-line test
!!$  !
!!$  xlv      = 2500000.0_rp        !Latent heat of vaporization
!!$  ep2      = 0.6217504_rp
!!$  svp1     = 0.6112000_rp   
!!$  svp2     = 17.67000_rp
!!$  svp3     = 29.65000_rp
!!$  svpt0    = 273.1500_rp
!!$  rhowater = 1000.000_rp  !density of rainwater
!!$  
!!$  !     retain kts,kte for the future work (splitting among parallel CPUs)
!!$  kts=1
!!$  kte=nz_nsa
!!$
!!$  !     allocate arrays in a normal fashion for a 2D array;
!!$  !     will have to copy element-based vars into 2D array first
!!$  !     and back to the element based syntax just before exiting
!!$
!!$!!  allocate ( rhocol_nsa(kts:kte),tcol_nsa(kts:kte),qvcol_nsa(kts:kte), &
!!$!!      qccol_nsa(kts:kte),qrcol_nsa(kts:kte))
!!$!!
!!$!!  allocate ( prodcol_nsa(kts:kte))
!!$!!  allocate ( vtcol_nsa(kts:kte), prodkcol_nsa(kts:kte), vtdencol_nsa(kts:kte),rdzkcol_nsa(kts:kte), &
!!$!!       rhokcol_nsa(kts:kte), factorcol_nsa(kts:kte), rdzwcol_nsa(kts:kte))
!!$
!!$  !  input arrays
!!$  !
!!$  !  t - potential temperature
!!$  !  qv, qc, qr  - mixing ratio (g/g dry air) of water vapor, cloud water
!!$  !                and rain water
!!$  !  pii         - exner function
!!$  !  dt_in - timestep
!!$  !  z  - height of (t,qv,p,rho) points in meters
!!$  !  dz8w - delta z between t points.
!!$  !
!!$  !  See Klemp and Wilhelmson (1978) Journal of the Atmospheric Sciences
!!$  !  Vol 35, pp 1070-1096 for more details
!!$  !
!!$
!!$  !   f5 = 237.3 * 17.27 * 2.5e6 / cpcoe_nsa
!!$  f5 = svp2*(svpt0-svp3)*xlv/cpcoe_nsa
!!$  ernmax = 0.0_rp
!!$  maxqrp = -100.0_rp
!!$
!!$  !------------------------------------------------------------------------------
!!$  ! parameters for the time split terminal advection
!!$  !------------------------------------------------------------------------------
!!$  max_heating = 0.
!!$  max_condense = 0.
!!$  max_rain = 0.
!!$
!!$  ! copy into 2D arrays (convert from element-based notation to one 2D array)
!!$  if(ndime < 3) then
!!$     !
!!$     ! 2D
!!$     !
!!$     do j=1,nz_nsa
!!$        rhocol_nsa(j) = q_col_nsa(1,j) + qref_col_nsa(1,j) !total density
!!$        tcol_nsa(j)   = q_col_nsa(4,j) + qref_col_nsa(4,j) !total theta
!!$        !qvcol_nsa(j) = q_col_nsa(5,j) + qref_col_nsa(5,j)
!!$        qvcol_nsa(j)  = q_col_nsa(5,j)                 !total vapor
!!$        qccol_nsa(j)  = q_col_nsa(6,j)                 !total cloud
!!$        qrcol_nsa(j)  = q_col_nsa(7,j)                 !total rain
!!$
!!$     end do !j
!!$  else
!!$     !
!!$     ! 3D
!!$     !
!!$     do j=1,nz_nsa
!!$        rhocol_nsa(j) = q_col_nsa(1,j) + qref_col_nsa(1,j) !total density
!!$        tcol_nsa(j)   = q_col_nsa(5,j) + qref_col_nsa(5,j) !total theta
!!$        !qvcol_nsa(j) = q_col_nsa(6,j) + qref_col_nsa(6,j)
!!$        qvcol_nsa(j)  = q_col_nsa(6,j)                 !total vapor
!!$        qccol_nsa(j)  = q_col_nsa(7,j)                 !total cloud
!!$        qrcol_nsa(j)  = q_col_nsa(8,j)                 !total rain
!!$
!!$     end do !j
!!$
!!$  end if
!!$
!!$  !
!!$  ! are all the variables ready for the microphysics?
!!$  !
!!$  ! start the microphysics
!!$  !
!!$  ! do the sedimentation first
!!$  !
!!$  crmax = 0.
!!$
!!$  !------------------------------------------------------------------------------
!!$  ! Terminal velocity calculation and advection, set up coefficients and
!!$  ! compute stable timestep
!!$  !------------------------------------------------------------------------------
!!$  do k = 1, kte-1
!!$     rdzkcol_nsa(k) = 1.0_rp/(z_nsa(k+1) - z_nsa(k) + 0.1_rp)
!!$  enddo
!!$  rdzkcol_nsa(kte) = 1.0_rp/(z_nsa(kte) - z_nsa(kte-1) + 0.1_rp)
!!$
!!$  do k = 1, kte
!!$     prodkcol_nsa(k) = qrcol_nsa(k)
!!$     rhokcol_nsa(k)  = rhocol_nsa(k)
!!$
!!$     qrr = max(0.0_rp,qrcol_nsa(k)*0.001_rp*rhokcol_nsa(k))  !ok
!!$     vtdencol_nsa(k) = sqrt(rhokcol_nsa(1)/rhokcol_nsa(k))   !ok
!!$     vtcol_nsa(k) = 36.34*(qrr**0.1364) * vtdencol_nsa(k)    !ok
!!$
!!$     !       vtmax = amax1(vtcol_nsa(k), vtmax)
!!$     crmax = max(vtcol_nsa(k)*dt*rdzkcol_nsa(k),crmax)
!!$  enddo
!!$
!!$  nfall = max(1_ip,nint(0.5_rp + crmax/max_cr_sedimentation)) ! courant number for big timestep.
!!$  dtfall = dt / float(nfall) ! splitting so courant number for sedimentation
!!$  time_sediment = dt      ! is stable
!!$
!!$  !------------------------------------------------------------------------------
!!$  ! Terminal velocity calculation and advection
!!$  ! Do a time split loop on this for stability.
!!$  !------------------------------------------------------------------------------
!!$
!!$  column_sedimentation: do while ( nfall > 0 )
!!$
!!$     time_sediment = time_sediment - dtfall
!!$     do k = 1, kte-1
!!$        factorcol_nsa(k) = dtfall*rdzkcol_nsa(k)/rhokcol_nsa(k)
!!$     enddo
!!$     factorcol_nsa(kte) = dtfall*rdzkcol_nsa(kte)
!!$
!!$     ppt=0.0_rp
!!$
!!$     k = 1
!!$     ppt=rhokcol_nsa(k)*prodkcol_nsa(k)*vtcol_nsa(k)*dtfall/rhowater
!!$     rainncv =ppt*1000.0_rp
!!$     rainnc =rainnc +  ppt*1000.0_rp ! unit = mm
!!$
!!$     !if(icol == ncol_nsa) &
!!$     if(rainnc > 0 ) &
!!$          print*,'kessler: rainnc mm', icol, rainnc
!!$     
!!$     !------------------------------------------------------------------------------
!!$     ! Time split loop, Fallout done with flux upstream
!!$     !------------------------------------------------------------------------------
!!$
!!$     do k = kts, kte-1
!!$        prodkcol_nsa(k) = prodkcol_nsa(k) - factorcol_nsa(k)  &   
!!$             * (rhokcol_nsa(k)*prodkcol_nsa(k)*vtcol_nsa(k) &  
!!$             -rhokcol_nsa(k+1)*prodkcol_nsa(k+1)*vtcol_nsa(k+1))
!!$     enddo
!!$
!!$     k = kte
!!$
!!$     prodkcol_nsa(k) = prodkcol_nsa(k) - factorcol_nsa(k)*prodkcol_nsa(k)*vtcol_nsa(k)
!!$
!!$     !------------------------------------------------------------------------------
!!$     ! compute new sedimentation velocity, and check/recompute new
!!$     ! sedimentation timestep if this isn't the last split step.
!!$     !------------------------------------------------------------------------------
!!$
!!$     if ( nfall > 1_ip ) then ! this wasn't the last split sedimentation timestep
!!$
!!$        nfall = nfall - 1_ip
!!$        crmax = 0.0_rp
!!$        do k = kts,kte
!!$           qrr = max(0.0_rp,prodkcol_nsa(k)*0.001_rp*rhokcol_nsa(k))
!!$           vtcol_nsa(k) = 36.34_rp*(qrr**0.1364_rp) * vtdencol_nsa(k)
!!$           !          vtmax = amax1(vtcol_nsa(k), vtmax)
!!$           crmax = max(vtcol_nsa(k)*time_sediment*rdzwcol_nsa(k),crmax)
!!$        enddo
!!$
!!$        nfall_new = max(1_ip,nint(0.5_rp+crmax/max_cr_sedimentation))
!!$        if (nfall_new /= nfall ) then
!!$           nfall = nfall_new
!!$           dtfall = time_sediment/nfall
!!$        end if
!!$
!!$     else  ! this was the last timestep
!!$
!!$        do k=kts,kte
!!$           prodcol_nsa(k) = prodkcol_nsa(k)
!!$        enddo
!!$        nfall = 0_ip  ! exit condition for sedimentation loop
!!$
!!$     endif
!!$
!!$  enddo column_sedimentation
!!$
!!$
!!$  ! now the conversion processes
!!$
!!$  !------------------------------------------------------------------------------
!!$  ! Production of rain and deletion of qc
!!$  ! Production of qc from supersaturation
!!$  ! Evaporation of QR
!!$  !------------------------------------------------------------------------------
!!$
!!$  do k = kts, kte
!!$     
!!$     factorn = 1.0_rp / (1.0_rp+c3*dt*max(0.0_rp,qrcol_nsa(k))**c4)
!!$
!!$     qrprod = qccol_nsa(k) * (1.0_rp - factorn) &
!!$          + factorn*c1*dt*max(qccol_nsa(k)-c2,0.0_rp)
!!$     
!!$     rcgs = 0.001_rp*rhocol_nsa(k)
!!$
!!$     qccol_nsa(k) = max(qccol_nsa(k) - qrprod,0.0_rp)
!!$     qrcol_nsa(k) = (qrcol_nsa(k) + prodcol_nsa(k)-qrcol_nsa(k))
!!$     qrcol_nsa(k) = max(qrcol_nsa(k) + qrprod,0.0_rp)
!!$
!!$     pii=(pcol_nsa(k)/1.e5)**(287.0_rp/1004.0_rp)
!!$     temp=tcol_nsa(k)*pii
!!$     pressure=pcol_nsa(k)
!!$     
!!$     gam = 2.5e+06/(1004.0_rp*pii)
!!$     !      qvs       = 380.*exp(17.27*(temp-273.)/(temp- 36.))/pressure
!!$     es        = 1000.0_rp*svp1*exp(svp2*(temp-svpt0)/(temp-svp3))
!!$     qvs       = ep2*es/(pressure-es)
!!$     !      prodcol_nsa(i,k) = (qvcol_nsa(i,k)-qvs) / (1.+qvs*f5/(temp-36.)**2)
!!$     prodcol_nsa(k) = (qvcol_nsa(k)-qvs) / (1.0_rp+pressure/(pressure-es)*qvs*f5/ &
!!$          (temp-svp3)**2)
!!$     
!!$     ern  = min(dt*(((1.6_rp+124.9_rp*(rcgs*qrcol_nsa(k))**0.2046_rp) &
!!$          *(rcgs*qrcol_nsa(k))**0.525_rp)/(2.55e8/(pressure*qvs) &      
!!$          +5.4e5))*(dim(qvs,qvcol_nsa(k))/(rcgs*qvs)),  &           
!!$          max(-prodcol_nsa(k)-qccol_nsa(k),0.0_rp),qrcol_nsa(k))
!!$
!!$     !
!!$     ! Update all variables
!!$     !
!!$     product = max(prodcol_nsa(k),-qccol_nsa(k))
!!$    ! if(product > 1e-8 ) then
!!$    !    PRINT*,'PRODUCTION', product, tcol_nsa(k)
!!$    ! end if
!!$     tcol_nsa(k)  = tcol_nsa(k) + gam*(product - ern)
!!$     qvcol_nsa(k) = max(qvcol_nsa(k) - product + ern,0.0_rp)
!!$     qccol_nsa(k) = qccol_nsa(k) + product
!!$     qrcol_nsa(k) = qrcol_nsa(k) - ern
!!$
!!$  enddo
!!$
!!$  !
!!$  ! transform modified variables back to the element-based notation
!!$  ! make sure only perturbations get modified, not the base state
!!$  !
!!$  if(ndime < 3) then
!!$     !
!!$     ! 2D
!!$     !
!!$     do j=1,nz_nsa       
!!$        q_col_nsa(1,j) = rhocol_nsa(j) - qref_col_nsa(1,j)
!!$        q_col_nsa(4,j) = tcol_nsa(j)   - qref_col_nsa(4,j)
!!$        !q_col_nsa(5,j) = qvcol_nsa(j)  - qref_col_nsa(5,j)
!!$        q_col_nsa(5,j) = qvcol_nsa(j)  
!!$        q_col_nsa(6,j) = qccol_nsa(j)
!!$        q_col_nsa(7,j) = qrcol_nsa(j)
!!$     end do !j
!!$  else
!!$     !
!!$     ! 3D
!!$     !
!!$     do j=1,nz_nsa     
!!$        q_col_nsa(1,j) = rhocol_nsa(j) - qref_col_nsa(1,j)
!!$        q_col_nsa(5,j) = tcol_nsa(j)   - qref_col_nsa(5,j)
!!$        !q_col_nsa(6,j) = qvcol_nsa(j)  - qref_col_nsa(6,j)
!!$        q_col_nsa(6,j) = qvcol_nsa(j)  
!!$        q_col_nsa(7,j) = qccol_nsa(j)
!!$        q_col_nsa(8,j) = qrcol_nsa(j)
!!$     end do !j
!!$  end if
!!$
!!$end subroutine kessler_col


!
! As version 1338 (working)
!
!-------------------------------------------------------------------------------------
!From Jim original, Simone, Dec 22, 2011
! Performs Kessler Physics on a col. of data
! Inputs:
! q_col = col of prognostic variables (rho', u, w, theta', qv', qc, qr)
! qref_col = col. of reference variables (rhoref, 0, 0, theteref, qvref, 0, 0)
! p = col of total pressure ( p = p' + pref)
! nvar = number of variables = 7
! Outputs
! rainnc = accumulation of rain at the surface
! z = col. of z-values
! nz = number of z-values in the col.
! dt = time-step
! James F. Kelly
! 29 January 2010
!-------------------------------------------------------------------------------------

!------------------------------------------------------------
! subroutine apply_physics
!------------------------------------------------------------
subroutine apply_physics(dt,rainnc,rainncv)

  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal

  implicit none

  real(rp)    ::  q(nvar_nsa,npoin),      q_ref(nvar_nsa,npoin)
  real(rp)    ::  q_col(nvar_nsa,nz_nsa), qref_col(nvar_nsa,nz_nsa)
  real(rp)    ::  press_col(nz_nsa)

  integer(ip) ::  icol, iz, ipoin, nz, ncomp, irestart, iloop, ifnp, j, istep
  real(rp)    ::  z_col(nz_nsa)
  !real(rp)   ::  press(npoin)
  real(rp)    ::  dt
  real(rp)    ::  rainnc(ncol_nsa), rainncv(ncol_nsa)

  character   :: fnp1*24, fnp*72

  if (nzone > 1) call runend("NSA_INITIAL_CONDITIONS: THIS SUB IS NOT PREPARED TO RUN WITH ZONES.")
  

  irestart = vtkrestart_nsa

  !Local variables:
  nz = nz_nsa
  ncomp=min(3_ip,ncomp_nsa)
  !
  ! Store local working arrays (q, q_ref):
  !
  if(ndime < 3) then
     !
     ! 2D
     ! 
     do ipoin=1,npoin
        
        !Reference/hydrostatic variables:
        q_ref(1,ipoin) = rekee_nsa(ndime+1,ipoin) !densi
        q_ref(2,ipoin) = bvess_nsa(1,ipoin,1)     !u-velo
        q_ref(3,ipoin) = bvess_nsa(ndime,ipoin,1) !v-velo
        q_ref(4,ipoin) = rekee_nsa(ndime+2,ipoin) !theta
    
        q_ref(5,ipoin) = bvess_nsa(ndofn_nsa+1,ipoin,1)       !q_vapor
        q_ref(6,ipoin) = bvess_nsa(ndofn_nsa+2,ipoin,1)       !q_cloud
        q_ref(7,ipoin) = bvess_nsa(ndofn_nsa+3,ipoin,1)       !q_rain
        
        !Solution variables (for us they are total so that we need to subtract the ref):
        !so that, to use them in kessler_col, we need to subtract
        !the reference values:
        q(1,ipoin) = densi(ipoin,ncomp)       - q_ref(1,ipoin)
        q(2,ipoin) = veloc(1,ipoin,ncomp)     - q_ref(2,ipoin)
        q(3,ipoin) = veloc(ndime,ipoin,ncomp) - q_ref(3,ipoin)
        q(4,ipoin) = tempe(ipoin,ncomp)       - q_ref(4,ipoin)
   
        q(5,ipoin) = conce(ipoin,1,ncomp)
        q(6,ipoin) = conce(ipoin,2,ncomp)
        q(7,ipoin) = conce(ipoin,3,ncomp)
        
     end do
  else
     !
     ! 3D
     !
     do ipoin=1,npoin

        !Reference/hydrostatic variables:
        q_ref(1,ipoin) = rekee_nsa(ndime+1,ipoin)        !densi  ok
        q_ref(2,ipoin) = bvess_nsa(1,ipoin,1)            !u-velo ok
        q_ref(3,ipoin) = bvess_nsa(ndime-1,ipoin,1)      !v-velo ok
        q_ref(4,ipoin) = bvess_nsa(ndime,ipoin,1)        !w-velo ok
        q_ref(5,ipoin) = rekee_nsa(ndime+2,ipoin)        !theta  ok

        q_ref(6, ipoin) = bvess_nsa(ndofn_nsa+1,ipoin,1)  !qv ok
        q_ref(7, ipoin) = bvess_nsa(ndofn_nsa+2,ipoin,1)  !qc ok
        q_ref(8, ipoin) = bvess_nsa(ndofn_nsa+3,ipoin,1)  !qr ok

!        q_ref(6, ipoin) = rekee_nsa(ndime+3 + 1, ipoin)  !qv ok testing
!        q_ref(7, ipoin) = rekee_nsa(ndime+3 + 2, ipoin)  !qc ok testing
!        q_ref(8, ipoin) = rekee_nsa(ndime+3 + 3, ipoin)  !qr ok testing

        !Solution variables (for us they are total so that we need to subtract the ref):
        !so that, to use them in kessler_col, we need to subtract
        !the reference values:
        q(1,ipoin) = densi(ipoin,ncomp)         - q_ref(1,ipoin)
        q(2,ipoin) = veloc(1,ipoin,ncomp)       - q_ref(2,ipoin)
        q(3,ipoin) = veloc(ndime-1,ipoin,ncomp) - q_ref(3,ipoin)
        q(4,ipoin) = veloc(ndime,ipoin,ncomp)   - q_ref(4,ipoin)
        q(5,ipoin) = tempe(ipoin,ncomp)         - q_ref(5,ipoin)
       
        q(6,ipoin) = conce(ipoin,1,ncomp)
        q(7,ipoin) = conce(ipoin,2,ncomp)
        q(8,ipoin) = conce(ipoin,3,ncomp)

     end do
  end if
  
  do icol = 1,ncol_nsa
     ! Extract a column of data
     if(ndime < 3) then
        !
        ! 2D
        !
        do iz = 1,nz_nsa
           ipoin = node_column_nsa(icol,iz)
           q_col(1,iz)    = q(1,ipoin)     !perturbed densi
           q_col(2,iz)    = q(2,ipoin)     !u-velo
           q_col(3,iz)    = q(3,ipoin)     !v-velo
           q_col(4,iz)    = q(4,ipoin)     !perturbed theta
           q_col(5,iz)    = q(5,ipoin)     !total vapor
           q_col(6,iz)    = q(6,ipoin)     !total cloud
           q_col(7,iz)    = q(7,ipoin)     !total rain

           qref_col(1,iz) = q_ref(1,ipoin) !reference densi
           qref_col(2,iz) = q_ref(2,ipoin) !reference u-velo
           qref_col(3,iz) = q_ref(3,ipoin) !reference v-velo
           qref_col(4,iz) = q_ref(4,ipoin) !reference theta
           qref_col(5,iz) = q_ref(5,ipoin) !reference vapor
           qref_col(6,iz) = q_ref(6,ipoin) !reference cloud
           qref_col(7,iz) = q_ref(7,ipoin) !reference rain

           press_col(iz)  = press(ipoin,ncomp) ! + press_ref(ipoin) 
           !This muyst be TOTAL pressure.
           !Add the hydrostatic value 
           !if your press() is the perturbed pressure.
           z_col(iz) = coord(ndime,ipoin)

           !print*,'d', q_col(1,iz),q_col(4,iz), q_col(5,iz), press_col(iz)
        end do
     else
        !
        ! 3D
        !
        do iz = 1,nz_nsa
           ipoin = node_column_nsa(icol,iz)
           q_col(1,iz)    = q(1,ipoin)     !perturbed densi
           q_col(2,iz)    = q(2,ipoin)     !u-velo
           q_col(3,iz)    = q(3,ipoin)     !v-velo
           q_col(4,iz)    = q(4,ipoin)     !w-velo
           q_col(5,iz)    = q(5,ipoin)     !perturbed theta

           q_col(6,iz)    = q(6,ipoin)     !total qv
           q_col(7,iz)    = q(7,ipoin)     !total qc
           q_col(8,iz)    = q(8,ipoin)     !total qr

           qref_col(1,iz) = q_ref(1,ipoin) !reference densi
           qref_col(2,iz) = q_ref(2,ipoin) !reference u-velo
           qref_col(3,iz) = q_ref(3,ipoin) !reference v-velo
           qref_col(4,iz) = q_ref(4,ipoin) !reference w-velo
           qref_col(5,iz) = q_ref(5,ipoin) !reference theta

           qref_col(6,iz) = q_ref(6,ipoin) !reference vapor
           qref_col(7,iz) = q_ref(7,ipoin) !reference cloud
           qref_col(8,iz) = q_ref(8,ipoin) !reference rain

           press_col(iz)  = press(ipoin,ncomp) ! + press_ref(ipoin) 
           !This muyst be TOTAL pressure.
           !Add the hydrostatic value 
           !if your press() is the perturbed pressure.
           z_col(iz) = coord(ndime,ipoin)

        end do
     end if

     !
     ! Do Kessler Physics on each column of data
     !
     !print *, minval(press_col), maxval(press_col)
     call kessler_col(icol, q_col,qref_col,press_col,&
          rainnc(icol),rainncv(icol),z_col,dt)

     !
     ! Put data
     !
     if(ndime < 3) then
        !
        ! 2D
        !
        do iz = 1,nz

           ipoin = node_column_nsa(icol,iz)
           q(1,ipoin) = q_col(1,iz) !perturbed densi
           q(4,ipoin) = q_col(4,iz) !perturbed theta
           q(5,ipoin) = q_col(5,iz) !total qv
           q(6,ipoin) = q_col(6,iz) !total qc
           q(7,ipoin) = q_col(7,iz) !total qr

           !Restore the total solution variables to continue in alya: 
           densi(ipoin,ncomp)   = q(1,ipoin) + q_ref(1,ipoin)
           tempe(ipoin,ncomp)   = q(4,ipoin) + q_ref(4,ipoin)
           conce(ipoin,1,ncomp) = q(5,ipoin)
           conce(ipoin,2,ncomp) = q(6,ipoin)
           conce(ipoin,3,ncomp) = q(7,ipoin)

        end do
     else
        !
        ! 3D
        !
        do iz = 1,nz

           ipoin = node_column_nsa(icol,iz)
           q(1,ipoin) = q_col(1,iz) !perturbed densi
           q(5,ipoin) = q_col(5,iz) !perturbed theta
           q(6,ipoin) = q_col(6,iz) !total qv
           q(7,ipoin) = q_col(7,iz) !total qc
           q(8,ipoin) = q_col(8,iz) !total qr

           !Restore the total solution variables to continue in alya: 
           densi(ipoin,ncomp)   = q(1,ipoin) + q_ref(1,ipoin)
           tempe(ipoin,ncomp)   = q(5,ipoin) + q_ref(5,ipoin)
           conce(ipoin,1,ncomp) = q(6,ipoin)
           conce(ipoin,2,ncomp) = q(7,ipoin)
           conce(ipoin,3,ncomp) = q(8,ipoin)

        end do
     end if
     
  end do

  print*,''
  print *, "nsa_kessler: Vapor qv: ", maxval(q(ndime+3,:)), minval(q(ndime+3,:))
  if (maxval(q(ndime+4,:)) > 0.0_rp .or. minval(q(ndime+4,:)) > 0.0_rp) then
     print *, "nsa_kessler: Cloud qc: ", maxval(q(ndime+4,:)), minval(q(ndime+4,:))
  end if
  if (maxval(q(ndime+5,:)) > 0.0_rp .or. minval(q(ndime+5,:)) > 0.0_rp)then
     print *, "nsa_kessler: Rain qr: ", maxval(q(ndime+5,:)), minval(q(ndime+5,:))
  end if

  !
  !Write vtk every n-steps
  !CAREFUL HERE: CHECK IT. Uncommenting this in 3D gives me error on DIVIDE BY ZERO!
!!$  if(ittim == 1 .or. mod(ittim,irestart) == 0 .or. ittim == mitim) then
!!$     istep_nsa = istep_nsa + 1
!!$     write(fnp1,'(i6)')istep_nsa
!!$     iloop=2 - int(log10(real(istep_nsa)))
!!$     do j=1,iloop
!!$        fnp1(j:j)='0000'
!!$     end do
!!$
!!$     !VTK:
!!$     fnp=trim('OUTVTK_alya') // '_' // trim(fnp1) // '.vtk'
!!$     call nsa_outvtk_kessler(fnp,q)
!!$
!!$     !Matlab:
!!$     fnp=trim('Qv_OUTMATLAB_alya') // '_' // trim(fnp1) // '.txt'
!!$     call nsa_outmatlab(fnp,q(5,:))
!!$
!!$     fnp=trim('Qc_OUTMATLAB_alya') // '_' // trim(fnp1) // '.txt'
!!$     call nsa_outmatlab(fnp,q(6,:))
!!$
!!$     fnp=trim('Qr_OUTMATLAB_alya') // '_' // trim(fnp1) // '.txt'
!!$     call nsa_outmatlab(fnp,q(7,:))
!!$
!!$     fnp=trim('Uvelo_OUTMATLAB_alya') // '_' // trim(fnp1) // '.txt'
!!$     call nsa_outmatlab(fnp,q(2,:))
!!$
!!$     fnp=trim('Wvelo_OUTMATLAB_alya') // '_' // trim(fnp1) // '.txt'
!!$     call nsa_outmatlab(fnp,q(ndime,:))
!!$
!!$     !  !Open sedimentation file to store the amount of 
!!$     !  !precipitated rain per column:
!!$     open(1, file="rain_accumulation.txt", position="append")
!!$!     write(1,*) "TIME"
!!$!     write(1,'(i7,E16.8)') ittim, cutim
!!$!     do icol =1,ncol_nsa
!!$!        write(1,'(i7,(E16.8,1x))') icol, rainncv(icol)
!!$!     end do
!!$     write(1,'(i7,9999999(E16.8,1x))') istep_nsa, rainncv(1:ncol_nsa)
!!$     close(1)
!!$
!!$  end if

end subroutine apply_physics

subroutine kessler_col(icol, q_col,qref_col,p_col,rainnc,rainncv,z,dt)
  !subroutine kessler_col(q_col,qref_col,p_col,rainnc,rainncv,z,dt)
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal
  
  implicit none

  !global arrays
  integer(ip)::  nz, nvar, icol
  real(rp)   ::  q_col(nvar_nsa,nz_nsa), qref_col(nvar_nsa,nz_nsa), z(nz_nsa)
  real(rp)   ::  p_col(nz_nsa)

  !rain arrays
  real(rp)   ::  rainnc, rainncv

  ! local variables

  integer(ip)  ::  i,j,k,ie
  integer(ip)  ::  icount,jcount
  integer(ip)  ::  its,ite,kts,kte

  ! variables

  real(rp), pointer :: rho(:),t(:),qv(:),qc(:),qr(:)

  ! local variables

  real(rp)   ::  xlv,ep2,svp1,svp2,svp3,svpt0,rhowater, pii

  ! local variables from the original module_mp_kessler.F

  real(rp)   :: qrprod, ern, gam, rcgs, rcgsi
  real(rp), pointer :: prod(:)
  real(rp), pointer :: vt(:), prodk(:), vtden(:), rdzk(:), &
                       rhok(:), factor(:), rdzw(:)
  integer(ip)  :: nfall, n, nfall_new
  real(rp)     :: qrr, pressure, temp, es, qvs, dz, dt
  real(rp)     :: f5, dtfall, rdz, product
  real(rp)     :: max_heating, max_condense, max_rain, maxqrp
  real(rp)     :: vtmax, ernmax, crmax, factorn, time_sediment
  real(rp)     :: qcr, factorr, ppt

  !----------------------------------------------------------------
  ! parameters from the original module_mp_kessler.F
  !----------------------------------------------------------------
  real(rp), parameter ::  c1     = 0.001_rp
  real(rp), parameter ::  c2     = 0.001_rp
  real(rp), parameter ::  c3     = 2.2_rp
  real(rp), parameter ::  c4     = 0.875_rp
  real(rp), parameter ::  fudge  = 1.0_rp
  real(rp), parameter ::  mxfall = 10.0_rp
  real(rp), parameter :: max_cr_sedimentation = 0.75_rp
  !----------------------------------------------------------------
  
  !Local constants:
  nz   = nz_nsa
  nvar = nvar_nsa

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
  
  !     retain kts,kte for the future work (splitting among parallel CPUs)
  kts=1
  kte=nz_nsa

  !     allocate arrays in a normal fashion for a 2D array;
  !     will have to copy element-based vars into 2D array first
  !     and back to the element based syntax just before exiting

  allocate ( rho(kts:kte),t(kts:kte),qv(kts:kte), &
       qc(kts:kte),qr(kts:kte))

  allocate ( prod(kts:kte))
  allocate ( vt(kts:kte), prodk(kts:kte), vtden(kts:kte),rdzk(kts:kte), &
       rhok(kts:kte), factor(kts:kte), rdzw(kts:kte))

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

  !   f5 = 237.3 * 17.27 * 2.5e6 / cpcoe_nsa
  f5 = svp2*(svpt0-svp3)*xlv/cpcoe_nsa
  ernmax = 0.0_rp
  maxqrp = -100.0_rp

  !------------------------------------------------------------------------------
  ! parameters for the time split terminal advection
  !------------------------------------------------------------------------------
  max_heating = 0.0_rp
  max_condense = 0.0_rp
  max_rain = 0.0_rp

  ! copy into 2D arrays (convert from element-based notation to one 2D array)
  if(ndime < 3) then
     !
     ! 2D
     !
     do j=1,nz_nsa
        rho(j) = q_col(1,j) + qref_col(1,j) !total density
        t(j)   = q_col(4,j) + qref_col(4,j) !total theta
        !qv(j) = q_col(5,j) + qref_col(5,j)
        qv(j)  = q_col(5,j)                 !total vapor
        qc(j)  = q_col(6,j)                 !total cloud
        qr(j)  = q_col(7,j)                 !total rain

     end do !j
  else
     !
     ! 3D
     !
     do j=1,nz_nsa
        rho(j) = q_col(1,j) + qref_col(1,j) !total density
        t(j)   = q_col(5,j) + qref_col(5,j) !total theta
        !qv(j) = q_col(6,j) + qref_col(6,j)
        qv(j)  = q_col(6,j)                 !total vapor
        qc(j)  = q_col(7,j)                 !total cloud
        qr(j)  = q_col(8,j)                 !total rain

     end do !j

  end if

  ! are all the variables ready for the microphysics?

  ! start the microphysics

  ! do the sedimentation first
  crmax = 0.0_rp

  !------------------------------------------------------------------------------
  ! Terminal velocity calculation and advection, set up coefficients and
  ! compute stable timestep
  !------------------------------------------------------------------------------

  do k = 1, kte-1
     rdzk(k) = 1.0_rp/(z(k+1) - z(k))
    enddo

  rdzk(kte) = 1.0_rp/(z(kte) - z(kte-1))

  do k = 1, kte
     prodk(k) = qr(k)
     rhok(k)  = rho(k)

     vtden(k) = sqrt(rhok(1)/rhok(k))          !ok
     qrr = max(0.0_rp,qr(k)*0.001_rp*rhok(k))  !ok
     vt(k) = 36.34*(qrr**0.1364) * vtden(k)    !ok

     !       vtmax = amax1(vt(k), vtmax)
     crmax = max(vt(k)*dt*rdzk(k),crmax)
  enddo

  !nfall = max(1_ip,nint(0.5_rp + crmax/max_cr_sedimentation)) ! courant number for big timestep.
  nfall = max(1_ip,int(0.5_rp + crmax/max_cr_sedimentation,ip)) ! courant number for big timestep.
  dtfall = dt / real(nfall,rp) ! splitting so courant number for sedimentation
  time_sediment = dt      ! is stable

  !---------------------------------------------------------------------------
  ! Terminal velocity calculation and advection
  ! Do a time split loop on this for stability.
  !---------------------------------------------------------------------------
  column_sedimentation: do while ( nfall > 0 )

     time_sediment = time_sediment - dtfall
     do k = 1, kte-1
        factor(k) = dtfall*rdzk(k)/rhok(k)
     enddo
     factor(kte) = dtfall*rdzk(kte)

     ppt=0.0_rp

     k = 1! !NOTE!!! In the original k = 1 because it uses the ground node.
     ppt=rhok(k)*prodk(k)*vt(k)*dtfall/rhowater
     rainncv =ppt*1000.0_rp
     rainnc =rainnc +  ppt*1000.0_rp ! unit = mm

     !if(icol == ncol_nsa) &
     if(rainnc > 0 ) &
          print*,'kessler: rainnc mm', icol, rainnc
     
     !------------------------------------------------------------------------------
     ! Time split loop, Fallout done with flux upstream
     !------------------------------------------------------------------------------

     do k = kts, kte-1
        prodk(k) = prodk(k) - factor(k)  &   
             * (rhok(k)*prodk(k)*vt(k) &  
             -rhok(k+1)*prodk(k+1)*vt(k+1))
     enddo

     k = kte
     prodk(k) = prodk(k) - factor(k)*prodk(k)*vt(k)

     !------------------------------------------------------------------------------
     ! compute new sedimentation velocity, and check/recompute new
     ! sedimentation timestep if this isn't the last split step.
     !------------------------------------------------------------------------------

     if ( nfall > 1_ip ) then ! this wasn't the last split sedimentation timestep

        nfall = nfall - 1_ip
        crmax = 0.0_rp
        do k = kts,kte
           qrr = max(0.0_rp,prodk(k)*0.001_rp*rhok(k))
           vt(k) = 36.34_rp*(qrr**0.1364_rp) * vtden(k)
           !          vtmax = amax1(vt(k), vtmax)
           crmax = max(vt(k)*time_sediment*rdzw(k),crmax)
        enddo

        !nfall_new = max(1_ip,nint(0.5_rp+crmax/max_cr_sedimentation))
        nfall_new = max(1_ip,int(0.5_rp+crmax/max_cr_sedimentation,ip))
        if (nfall_new /= nfall ) then
           nfall = nfall_new
           dtfall = time_sediment/nfall
        end if

     else  ! this was the last timestep

        do k=kts,kte
           prod(k) = prodk(k)
        enddo
        nfall = 0_ip  ! exit condition for sedimentation loop

     endif

  enddo column_sedimentation


  ! now the conversion processes
  !------------------------------------------------------------------------------
  ! Production of rain and deletion of qc
  ! Production of qc from supersaturation
  ! Evaporation of QR
  !------------------------------------------------------------------------------

  do k = kts, kte
     
     factorn = 1.0_rp / (1.0_rp+c3*dt*max(0.0_rp,qr(k))**c4)

     qrprod = qc(k) * (1.0_rp - factorn) &
          + factorn*c1*dt*max(qc(k)-c2,0.0_rp)
     
     rcgs = 0.001_rp*rho(k)

     qc(k) = max(qc(k) - qrprod,0.0_rp)
     qr(k) = (qr(k) + prod(k)-qr(k))
     qr(k) = max(qr(k) + qrprod,0.0_rp)

     pii=(p_col(k)/1.e5)**(287.0_rp/1004.0_rp)
     temp=t(k)*pii
     pressure=p_col(k)
     
     gam = 2.5e+06/(1004.0_rp*pii)
     !      qvs       = 380.*exp(17.27*(temp-273.)/(temp- 36.))/pressure
     es        = 1000.0_rp*svp1*exp(svp2*(temp-svpt0)/(temp-svp3))
     qvs       = ep2*es/(pressure-es)
     !      prod(i,k) = (qv(i,k)-qvs) / (1.+qvs*f5/(temp-36.)**2)
     prod(k) = (qv(k)-qvs) / (1.0_rp+pressure/(pressure-es)*qvs*f5/ &
          (temp-svp3)**2)
     
     ern  = min(dt*(((1.6_rp+124.9_rp*(rcgs*qr(k))**0.2046_rp) &
          *(rcgs*qr(k))**0.525_rp)/(2.55e8/(pressure*qvs) &      
          +5.4e5))*(dim(qvs,qv(k))/(rcgs*qvs)),  &           
          max(-prod(k)-qc(k),0.0_rp),qr(k))

     !
     ! Update all variables
     !
     product = max(prod(k),-qc(k))
    ! if(product > 1e-8 ) then
    !    PRINT*,'PRODUCTION', product, t(k)
    ! end if
     t (k) = t(k) + gam*(product - ern)
     qv(k) = max(qv(k) - product + ern,0.0_rp)
     qc(k) = qc(k) + product
     qr(k) = qr(k) - ern

  enddo

  !
  ! transform modified variables back to the element-based notation
  ! make sure only perturbations get modified, not the base state
  !
  if(ndime < 3) then
     !
     ! 2D
     !
     do j=1,nz_nsa       
        q_col(1,j) = rho(j) - qref_col(1,j)
        q_col(4,j) = t(j)   - qref_col(4,j)
        !q_col(5,j) = qv(j)  - qref_col(5,j)
        q_col(5,j) = qv(j)  
        q_col(6,j) = qc(j)
        q_col(7,j) = qr(j)
     end do !j
  else
     !
     ! 3D
     !
     do j=1,nz_nsa     
        q_col(1,j) = rho(j) - qref_col(1,j)
        q_col(5,j) = t(j)   - qref_col(5,j)
        !q_col(6,j) = qv(j)  - qref_col(6,j)
        q_col(6,j) = qv(j)  
        q_col(7,j) = qc(j)
        q_col(8,j) = qr(j)
     end do !j
  end if

  ! just before exiting deallocate
  deallocate (rho,t,qv,qc,qr)
  deallocate (prod)
  deallocate (vt, prodk, vtden, rdzk, rhok, factor, rdzw)
end subroutine kessler_col
