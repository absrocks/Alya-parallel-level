!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_checkp.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966   
!> @brief   Inject a particle
!> @details Define initial conditions for particle ILAGR
!>          ILAGR: accumulated particle number
!>          KLAGR: local particle number
!>          LAGRTYP(KLAGR) % ILAGR ....... Absolute particle number
!>          LAGRTYP(KLAGR) % ITYPE ....... Particle type
!>          LAGRTYP(KLAGR) % IELEM ....... Host element of particle
!>          LAGRTYP(KLAGR) % KFL_EXIST ... -5 (particle just injected, will be put to 1 right after)
!>          LAGRTYP(KLAGR) % COORD(:) .... Coordinates
!>          LAGRTYP(KLAGR) % DT .......... Current time step dt^n+1
!>          LAGRTYP(KLAGR) % DTO ......... Last time step dt^n
!>          LAGRTYP(KLAGR) % V_FLUID_K ....... Fluid velocity
!>
!>          LAGRTYP(KLAGR) % VELOC ....... Particle velocity at t^n+1
!>
!>          Example with two types:
!>
!>          1. First injection:
!> 
!>          NLAGR_POS =  9
!>          NLAGR_TOT = 18
!>          NLAGR_NEW = 14
!>          o o o   o o o
!>          o o x   o o x
!>          o o x   o o x        
!>          type1   type2
!>
!>          => NLACC_PTS=14
!>
!>          1. Second injection:
!> 
!>          NLAGR_POS =  9
!>          NLAGR_TOT = 18
!>          NLAGR_NEW = 16
!>          x o o   x o o
!>          o o o   o o o
!>          o o o   o o o        
!>          type1   type2
!>
!>          => NLACC_PTS=30
!>
!>
!>
!> @} 
!-----------------------------------------------------------------------

subroutine pts_checkp(nlagr_pos,nlagr_new,particle_position,particle_injector,particle_place,particle_type)
  
  use def_parame,    only : xmaxint4, pi
  use def_kintyp
  use def_master
  use def_kermod
  use def_domain
  use def_partis 
  use mod_ker_proper 
  use mod_memory,    only : memory_alloca
  use mod_memory,    only : memory_deallo
  use mod_elsest,    only : elsest_host_element
  use mod_maths,     only : maths_vector_to_new_basis
  use mod_maths,     only : maths_local_orthonormal_basis
  
  implicit none
  integer(ip),  intent(in)    :: nlagr_pos                   !< Number of injected particles
  integer(ip),  intent(out)   :: nlagr_new                   !< Number of new particles owned by myself
  integer(ip),  intent(in)    :: particle_injector(*)        !< (nlagr_pos)
  real(rp),     intent(in)    :: particle_position(ndime,*)  !< (ndime,nlagr_pos)
  integer(ip),  intent(inout) :: particle_place(*)           !< (nlagr_pos)
  integer(ip),  intent(inout) :: particle_type(*)            !< (nlagr_pos)

  integer(ip)                 :: new_size
  integer(ip)                 :: ielem,inode,ilagr,jlagr,klagr_last
  integer(ip)                 :: ipoin,iinj,itype
  integer(ip)                 :: iprop,pnode,dumm0,klagr,nlagr_free
  real(rp)                    :: coloc(3),deriv(ndime,mnode),xx(3),nn(3),basis(ndime,ndime)
  real(rp)                    :: xa,ya,za,norma
  real(rp)                    :: dummr(3),shapf(mnode),relse_sav,rad,v_tmp(3),nside
  real(rp)                    :: grafo,buofo,denfl,denpa,dista,theta,phi,tfact
  real(rp)                    :: tau_p,diame,visfl            
  real(rp),     pointer       :: host_shapf(:,:)
  integer(ip),  pointer       :: host_element(:)

  nullify(host_shapf)
  nullify(host_element)
  call memory_alloca(mem_modul(1:2,modul),'host_element','pts_checkp',host_element,nlagr_pos)
  call memory_alloca(mem_modul(1:2,modul),'host_shapf'  ,'pts_checkp',host_shapf,mnode,nlagr_pos)

  !----------------------------------------------------------------------
  ! 
  ! Look for host elements
  !
  !----------------------------------------------------------------------
  relse_sav   = relse(1)
  relse(1)    = 0.0_rp
  nlagr_new   = 0
  xx          = 0.0_rp
  xa          = 0.0_rp
  ya          = 0.0_rp
  za          = 0.0_rp

  !$OMP PARALLEL DO SCHEDULE (DYNAMIC,500)                          & 
  !$OMP DEFAULT   (NONE)                                            &    
  !$OMP PRIVATE   (ilagr,xx,dummr,ielem,deriv,coloc,dista)          &
  !$OMP SHARED    (nlagr_pos,particle_position,ielse,relse,         &
  !$OMP            kfl_exacs_pts,meshe,ndivi,                       &
#ifndef NDIMEPAR
  !$OMP            ndime,                                           &
#endif
  !$OMP            host_shapf,host_element)                         &
  !$OMP REDUCTION (+:nlagr_new)
  do ilagr = 1,nlagr_pos
     xx(1:ndime) = particle_position(1:ndime,ilagr)
     if( kfl_exacs_pts /= 0 ) then
        call pts_exacso(1_ip,0.0_rp,dummr,dummr,xx)
     end if
     call elsest_host_element(&
          ielse,relse,1_ip,meshe(ndivi),xx,ielem,&
          host_shapf(:,ilagr),deriv,coloc,dista)
     if( ielem > 0 ) then
        host_element(ilagr) = ielem
        nlagr_new           = nlagr_new + 1 ! Number of new particles injected successfully
     end if
  end do
  !$OMP END PARALLEL DO

  relse(1)  = relse_sav

  !----------------------------------------------------------------------
  !
  ! Reallocate LAGRTYP if necessary
  !
  !----------------------------------------------------------------------
  !
  ! Counter number of free places in LAGRTYP
  !
  nlagr_free = 0
  do klagr = 1,mlagr
     if( lagrtyp(klagr) % kfl_exist == 0 ) then
        nlagr_free = nlagr_free + 1 
     end if
  end do
  !
  ! Reallocate LAGRTYP if necessary
  ! I currently have MLAGR places
  ! NLAGR_FREE are available
  ! I need NLAGR_NEW new positions
  !
  if( nlagr_new > nlagr_free ) then
     dista = 1.2_rp*real(mlagr+nlagr_new-nlagr_free,rp)
     if( ip == 4 .and. dista > xmaxint4 ) then
        call runend('PTS_CHECKP: GOT TO 8 BYTES INTEGERS')
     else
        new_size = int(dista,ip)
        call pts_reallocate(new_size)
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Fill in LAGRTYP with new particles with host element
  !
  !----------------------------------------------------------------------

  klagr_last = 0

  !-$OMP PARALLEL DO SCHEDULE (DYNAMIC,500)                                          & 
  !-$OMP DEFAULT       (NONE)                                                        &    
  !-$OMP FIRSTPRIVATE  (klagr_last)                                                  &
  !-$OMP PRIVATE       (ilagr,ielem,xx,shapf,dumm0,klagr,itype,inode,ipoin,vv,iinj,  &
  !-$OMP                nx,ny,nz,sigma,xc,yc,zc,pnode,dummr,denfl,denpa,grafo,buofo) &       
  !-$OMP SHARED        (nlagr_pos,host_element,number_types_pts,particle_position,   &
  !-$OMP                host_shapf,mlagr,lagrtyp,particle_place,parttyp,gravi,       &
  !-$OMP                lnnod,lnods,kfl_injve_pts,particle_injector,ndime,grnor,     &
  !-$OMP                kfl_exacs_pts,mnode,nlacc_pts,kfl_adapt_pts,dtime,dtmin_pts, &
  !-$OMP                advec,parla_pts,parla2_pts)                                  &
  do ilagr = 1,nlagr_pos

     ielem = host_element(ilagr)

     if( ielem > 0 ) then

        xx(1:ndime)    = particle_position(1:ndime,ilagr)
        itype          = particle_type(ilagr)
        shapf(1:mnode) = host_shapf(1:mnode,ilagr)
        !
        ! Look for new available particle position: 
        !
        ! OPTIMIZE: klagr could be 0 at the begining of injection
        ! so that we do not start from 0 each time!!!!!
        !
        klagr = klagr_last
        loop_klagr: do while( klagr < mlagr )
           klagr = klagr + 1
           if( lagrtyp(klagr) % kfl_exist == 0 ) then
              exit loop_klagr
           end if
        end do loop_klagr
        klagr_last            = klagr
        particle_place(ilagr) = klagr                    ! To compare with my neighbor
        !
        ! Initialize particle
        !
        lagrtyp(klagr)             =  lagrtyp_init 
        lagrtyp(klagr) % ilagr     =  ilagr + nlacc_pts  ! Particle absolute ID

        lagrtyp(klagr) % itype     =  itype
        lagrtyp(klagr) % ielem     =  ielem
        lagrtyp(klagr) % ittim     =  0
        lagrtyp(klagr) % kfl_exist = -5
        lagrtyp(klagr) % coord     =  xx
        lagrtyp(klagr) % t         =  cutim - dtime

        ! Injection time
        lagrtyp(klagr) % t_inject  =  cutim - dtime

        ! Initial local & effective Stokes number
        diame                      = parttyp(itype) % diame      ! Particle diameter
        denpa                      = parttyp(itype) % denpa      ! Particle density
        pnode                      = lnnod(ielem)
        call ker_proper('VISCO','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)   ! Fluid viscosity
        visfl                      = max(zeror,dummr(1))
        tau_p                      = (denpa * diame * diame)/(18.0_rp*visfl)
        call pts_stk_local(ielem,tau_p, lagrtyp(klagr) % t_inject, lagrtyp(klagr) % Stk(1), lagrtyp(klagr) % Stk(2))    
        !
        ! Initial time step
        !
        if( parttyp(itype) % kfl_modla == 2 ) then
           if(      parttyp(itype) % kfl_tstep == 0 ) then
              lagrtyp(klagr) % dt_k = dtime
           else if( parttyp(itype) % kfl_tstep <  0 ) then
              lagrtyp(klagr) % dt_k = parttyp(itype) % dtime
           else
              lagrtyp(klagr) % dt_k = dtmin_pts
           end if
        else
           lagrtyp(klagr) % dt_k =  dtime
        end if
        lagrtyp(klagr) % dt_km1 = lagrtyp(klagr) % dt_k
        lagrtyp(klagr) % dt_km2 = lagrtyp(klagr) % dt_k
        !
        ! Fluid velocity
        !
        lagrtyp(klagr) % v_fluid_k =  0.0_rp
        do inode = 1,lnnod(ielem)
           ipoin = lnods(inode,ielem)  
           lagrtyp(klagr) % v_fluid_k(1:ndime) = lagrtyp(klagr) % v_fluid_k(1:ndime) &
                + shapf(inode) * advec(1:ndime,ipoin,3)
        end do
        lagrtyp(klagr) % v_fluid_km1(1:ndime) = lagrtyp(klagr) % v_fluid_k(1:ndime)
        lagrtyp(klagr) % v_fluid_km2(1:ndime) = lagrtyp(klagr) % v_fluid_k(1:ndime)
        !
        ! Initial particle velocity
        !
        ! If there is a drag force, initial particle velocity is fluid velocity
        !   
        if (     kfl_injve_pts == -1 ) then
           !
           ! Zero velocity
           !
           lagrtyp(klagr) % veloc(1:ndime) = 0.0_rp 

        else if ( kfl_injve_pts == 0 ) then
           !
           ! Fluid velocity
           !
           lagrtyp(klagr) % veloc(1:ndime) = lagrtyp(klagr) % v_fluid_k(1:ndime)

        else if ( kfl_injve_pts == 1 ) then
           !
           ! With respect to normal
           !
           vv   = parla2_pts(1)
           iinj = particle_injector(ilagr)
           nx   = parla_pts(iinj,5)
           ny   = parla_pts(iinj,6)
           nz   = parla_pts(iinj,7)
           lagrtyp(klagr) % veloc(1) =  vv*nx
           lagrtyp(klagr) % veloc(2) = -vv*ny
           if( ndime == 3 ) then
              lagrtyp(klagr) % veloc(3) = vv*nz
           end if

        else if ( kfl_injve_pts == 2 ) then
           !
           ! Gaussian velocity injector
           ! f(x) = 1/sqrt(2*pi*sgima) * exp[(x-mu)^2/(2*sigma)]
           !
           vv    = parla2_pts(1)
           iinj  = particle_injector(ilagr)
           if(  kfl_injla_pts(iinj) /= 4 ) call runend('PTS_CHECKP: GAUSSIAN VELOCITY ONLY VALID IF CIRCLE INJECTOR')
           xc    = parla_pts(iinj,1)
           yc    = parla_pts(iinj,2)
           zc    = parla_pts(iinj,3)
           nx    = parla_pts(iinj,5)
           ny    = parla_pts(iinj,6)
           nz    = parla_pts(iinj,7)                 
           sigma = parla2_pts(2)
           lagrtyp(klagr) % veloc(1) = vv * nx * exp(-((xx(1)-xc)**2 + (xx(2)-yc)**2 + (xx(3)-zc)**2) / (2*sigma)) 
           lagrtyp(klagr) % veloc(2) = vv * ny * exp(-((xx(1)-xc)**2 + (xx(2)-yc)**2 + (xx(3)-zc)**2) / (2*sigma))
           if( ndime == 3 ) then
              lagrtyp(klagr) % veloc(3) = vv * nz * exp(-((xx(1)-xc)**2 + (xx(2)-yc)**2 + (xx(3)-zc)**2) / (2*sigma))
           end if

        else if ( kfl_injve_pts == 3 ) then
           !
           ! Conic velocity normal
           !
           iinj  = particle_injector(ilagr)
           if(  kfl_injla_pts(iinj) /= 4 ) then
              call runend('PTS_CHECKP: CONIC VELOCITY ONLY VALID IF CIRCLE INJECTOR')
           end if
           vv    = parla2_pts(1)
           sigma = parla2_pts(2)*(pi/180.0_rp)
           xc    = parla_pts(iinj,1)
           yc    = parla_pts(iinj,2)
           zc    = parla_pts(iinj,3)
           rad   = parla_pts(iinj,4)
           nn(1)    = parla_pts(iinj,5)
           nn(2)    = parla_pts(iinj,6)
           nn(3)    = parla_pts(iinj,7)
           nside    = parla_pts(iinj,8)
           basis(1:ndime,1) = nn(1:3)

           if( ndime == 3 ) then    

              xa = xc + rad*cos(sigma)/sin(sigma)*nn(1)
              ya = yc + rad*cos(sigma)/sin(sigma)*nn(2)
              za = zc + rad*cos(sigma)/sin(sigma)*nn(3)

              v_tmp(1) = xa - xx(1)
              v_tmp(2) = ya - xx(2)
              v_tmp(3) = za - xx(3)

              norma    = v_tmp(1)**2 +  v_tmp(2)**2 + v_tmp(3)**2
              norma    = 1.0_rp / sqrt(norma)
              v_tmp = norma * v_tmp

              lagrtyp(klagr) % veloc(1:3) = - vv * v_tmp(1:3)

              ! tfact = sigma*(sqrt((xx(1)-xc)**2+(xx(2)-yc)**2+(xx(3)-zc)**2))/rad
              ! v_tmp(1) = vv * cos(tfact)
              ! v_tmp(2) = vv * sin(tfact) * sin(phi)
              ! v_tmp(3) = vv * sin(tfact) * cos(phi)

              ! v_tmp(2) = vv * sin(tfact) * cos(phi)
              ! v_tmp(3) = vv * sin(tfact) * sin(phi)

              ! call maths_local_orthonormal_basis(ndime,basis)
              ! call maths_vector_to_new_basis(ndime,basis,v_tmp)
              ! lagrtyp(klagr) % veloc(1:3) = v_tmp(1:3)
              !  print*,v_tmp
           end if
        else if ( kfl_injve_pts == 4 ) then
           !
           ! Spray velocity normal
           !
           iinj  = particle_injector(ilagr)
           if(  kfl_injla_pts(iinj) /= 4 ) then
              call runend('PTS_CHECKP: CONIC VELOCITY ONLY VALID IF CIRCLE INJECTOR')
           end if
           vv    = parla2_pts(1)
           sigma = parla2_pts(2)*(pi/180.0_rp)
           xc    = parla_pts(iinj,1)
           yc    = parla_pts(iinj,2)
           zc    = parla_pts(iinj,3)
           rad   = parla_pts(iinj,4)
           nn(1)    = -parla_pts(iinj,5)
           nn(2)    = -parla_pts(iinj,6)
           nn(3)    = -parla_pts(iinj,7)
           nside    = parla_pts(iinj,8)
           basis(1:ndime,1) = nn(1:3)

           if( ndime == 3 ) then    

              xa = xc + rad*cos(sigma)/sin(sigma)*nn(1)
              ya = yc + rad*cos(sigma)/sin(sigma)*nn(2)
              za = zc + rad*cos(sigma)/sin(sigma)*nn(3)
              !print*,'x',xa,xx(1)
              !print*,'y',ya,yy(1)
              !print*,'z',za,zz(1)
              v_tmp(1) = xa - xx(1)
              v_tmp(2) = ya - xx(2)
              v_tmp(3) = za - xx(3)

              norma    = v_tmp(1)**2 +  v_tmp(2)**2 + v_tmp(3)**2
              norma    = 1.0_rp / sqrt(norma)
              v_tmp = norma * v_tmp
              lagrtyp(klagr) % veloc(1:3) = - vv * v_tmp(1:3)
              !print*,'hola-spray',lagrtyp(klagr) % veloc(1:3)
           end if

        else if ( kfl_injve_pts == 5 ) then
           !
           ! Constant velocity 
           !
           iinj                      = particle_injector(ilagr)
           lagrtyp(klagr) % veloc(1) = parla_pts(iinj,4) 
           lagrtyp(klagr) % veloc(2) = parla_pts(iinj,5) 
           if( ndime == 3 ) then
              lagrtyp(klagr) % veloc(3) = parla_pts(iinj,6) 
           end if
           
        end if
        !
        ! Initial acceleration
        !
        if( parttyp(itype) % kfl_modla == 2 ) then
           !pnode    = lnnod(ielem)
           call ker_proper('DENSI','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)   ! Fluid density
           denfl    = dummr(1)          
           !denpa    = parttyp(itype) % denpa                                     ! Particle density
           grafo    = real( parttyp(itype) % kfl_grafo, rp )                     ! Gravity  force = 1.0
           buofo    = real( parttyp(itype) % kfl_buofo, rp )                     ! Buoyancy force = 1.0  
           lagrtyp(klagr) % accel(1:ndime) = grnor * gravi(1:ndime) * ( grafo - buofo * denfl / denpa )
        end if
        !
        ! Injected particles receive default value of properties from particle type
        !
        lagrtyp(klagr) % prope(1:mlapr) = parttyp(itype) % prope(1:mlapr)
        !
        ! Exact solution
        !
        if( kfl_exacs_pts /= 0 ) then
           call pts_exacso(2_ip,0.0_rp,lagrtyp(klagr) % accel,lagrtyp(klagr) % veloc,lagrtyp(klagr) % coord)
        end if
     end if
  end do
  !-$OMP END PARALLEL DO

  !
  ! Deallocate
  !
  call memory_deallo(mem_modul(1:2,modul),'host_element','pts_checkp',host_element)
  call memory_deallo(mem_modul(1:2,modul),'host_shapf'  ,'pts_checkp',host_shapf)

end subroutine pts_checkp
