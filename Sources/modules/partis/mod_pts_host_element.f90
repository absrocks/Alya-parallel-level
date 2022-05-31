!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    mod_pts_host_element.f90
!> @author  houzeaux
!> @date    2020-07-06
!> @brief   Check host element
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
!-----------------------------------------------------------------------

module mod_pts_host_element

  use def_kintyp_basic,            only : ip,rp
  use def_parame,                  only : xmaxint4,pi
  use def_master,                  only : zeror
  use def_master,                  only : kfl_paral
  use def_master,                  only : dtime
  use def_master,                  only : cutim
  use def_master,                  only : therm
  use def_master,                  only : advec
  use def_master,                  only : modul
  use def_master,                  only : mem_modul
  use def_master,                  only : IPARALL
  use def_master,                  only : INOTMASTER
  use def_master,                  only : INOTEMPTY
  use def_domain,                  only : ndime
  use def_domain,                  only : mnode
  use def_domain,                  only : lnods
  use def_domain,                  only : lnnod
  use def_domain,                  only : meshe
  use def_kermod,                  only : ielse,relse
  use def_kermod,                  only : ndivi
  use def_kermod,                  only : grnor,gravi
  use mod_ker_proper,              only : ker_proper
  use mod_memory,                  only : memory_alloca
  use mod_memory,                  only : memory_deallo
  use mod_elsest,                  only : elsest_host_element
  use mod_maths,                   only : maths_vector_to_new_basis
  use mod_maths,                   only : maths_local_orthonormal_basis
  use mod_physics,                 only : physics_sphere_mass
  use mod_pts_thermodynamic,       only : pts_thermodynamic_properties
  use mod_ker_space_time_function, only : ker_space_time_function
  use mod_communications,          only : PAR_POINT_TO_POINT_ARRAY_OPERATION
  use mod_communications,          only : PAR_MAX
  use mod_communications,          only : PAR_SUM
  use def_partis 
  implicit none
  private

  character(20), parameter :: vacal = 'mod_pts_host_element'
  
  public :: pts_host_element_and_initialize
  public :: pts_host_element
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-06
  !> @brief   Find host elements
  !> @details Find the host elements for a list of particles
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_host_element_parall(&
       nlagr_new,nlagr_inj,particle_place)

    integer(ip),            intent(out)   :: nlagr_new         !< Number of new particles owned by myself
    integer(ip),            intent(inout) :: nlagr_inj         !< New particles
    integer(ip),  pointer,  intent(inout) :: particle_place(:) !< Place of new particles
    integer(ip)                           :: ii,jj
        
    if( IPARALL ) then
       call PAR_POINT_TO_POINT_ARRAY_OPERATION(particle_place, 'IN MY ZONE', 'MIN RANK OR NEGATIVE')
       if (INOTMASTER) then
          nlagr_new = 0
          do ii = 1,nlagr_inj
             
             if(      particle_place(ii) < 0 ) then
                !
                ! One of my neighbors is in charge of this particle
                !
                jj = -particle_place(ii)
                lagrtyp(jj) % kfl_exist = 0
                
             else if( particle_place(ii) > 0 ) then
                !
                ! I take this particle
                !
                nlagr_new = nlagr_new + 1
                
             else if( particle_place(ii) == 0 ) then
                !
                ! Particle not found
                !
                continue
                
             end if
          end do
       end if
       call PAR_MAX(nlagr_inj,'IN MY CODE')
       call PAR_SUM(nlagr_new,'IN MY CODE')
    end if

  end subroutine pts_host_element_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-06
  !> @brief   Find host elements
  !> @details Find the host elements for a list of particles
  !> 
  !-----------------------------------------------------------------------
  
  subroutine pts_host_element_and_initialize(&
       nlagr_pos,nlagr_new,nlagr_inj,particle_position,particle_injector,&
       particle_type,particle_diameter)

    integer(ip),                     intent(in)    :: nlagr_pos                   !< Number of injected particles
    integer(ip),                     intent(inout) :: nlagr_new                   !< Number of new particles owned by myself
    integer(ip),                     intent(inout) :: nlagr_inj                   !< Number of injected particles
    real(rp),               pointer, intent(in)    :: particle_position(:,:)      !< (ndime,nlagr_pos)
    integer(ip),            pointer, intent(in)    :: particle_injector(:)        !< (nlagr_pos)
    integer(ip),            pointer, intent(in)    :: particle_type(:)            !< (nlagr_pos)
    real(rp),     optional, pointer, intent(in)    :: particle_diameter(:)        !< (nlagr_pos)
    real(rp),               pointer                :: host_shapf(:,:)          
    integer(ip),            pointer                :: host_element(:)
    integer(ip),            pointer                :: particle_place(:) 

    nlagr_new = 0
    nullify(host_shapf)
    nullify(host_element)
    nullify(particle_place)
    !
    ! Alllocate
    !
    call memory_alloca(mem_modul(1:2,modul),'HOST_ELEMENT'  ,vacal,host_element  ,nlagr_pos)
    call memory_alloca(mem_modul(1:2,modul),'HOST_SHAPF'    ,vacal,host_shapf    ,mnode,nlagr_pos)
    call memory_alloca(mem_modul(1:2,modul),'PARTICLE_PLACE',vacal,particle_place,max(1_ip,nlagr_inj))

    if( INOTEMPTY .and. nlagr_pos > 0 ) then
       !
       ! Find host elements
       !
       call pts_host_element(&
            nlagr_pos,nlagr_new,particle_position,&
            host_shapf,host_element)
       !
       ! Initialize particles
       !
       call pts_initial_condition(&
            nlagr_pos,nlagr_new,particle_position,particle_injector,&
            particle_place,particle_type,host_shapf,host_element,&
            particle_diameter)
    end if
    !
    ! Parallelization, treat particles with multiple hosts
    !
    call pts_host_element_parall(&
         nlagr_new,nlagr_inj,particle_place)
    !
    ! Deallocate
    !
    call memory_deallo(mem_modul(1:2,modul),'HOST_ELEMENT'  ,vacal,host_element)
    call memory_deallo(mem_modul(1:2,modul),'HOST_SHAPF'    ,vacal,host_shapf)
    call memory_deallo(mem_modul(1:2,modul),'PARTICLE_PLACE',vacal,particle_place)

  end subroutine pts_host_element_and_initialize
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-06
  !> @brief   Find host elements
  !> @details Find the host elements for a list of particles
  !> 
  !-----------------------------------------------------------------------
  
  subroutine pts_initial_condition(&
       nlagr_pos,nlagr_new,particle_position,particle_injector,&
       particle_place,particle_type,host_shapf,host_element,&
       particle_diameter)

    integer(ip),                      intent(in)    :: nlagr_pos                   !< Number of injected particles
    integer(ip),                      intent(in)    :: nlagr_new                   !< Number of new particles owned by myself
    real(rp),               pointer,  intent(in)    :: particle_position(:,:)      !< (ndime,nlagr_pos)
    integer(ip),            pointer,  intent(in)    :: particle_injector(:)        !< (nlagr_pos)
    integer(ip),            pointer,  intent(inout) :: particle_place(:)           !< (nlagr_pos)
    integer(ip),            pointer,  intent(in)    :: particle_type(:)            !< (nlagr_pos)
    real(rp),               pointer,  intent(in)    :: host_shapf(:,:)             !< Shape functions
    integer(ip),            pointer,  intent(in)    :: host_element(:)             !< Host elements
    real(rp),     optional,           intent(in)    :: particle_diameter(:)        !< (nlagr_pos)

    integer(ip)                                     :: new_size
    integer(ip)                                     :: ielem,inode,ilagr,klagr_last
    integer(ip)                                     :: ipoin,iinj,itype
    integer(ip)                                     :: pnode,dumm0,klagr,nlagr_free
    real(rp)                                        :: xx(3),nn(3),basis(ndime,ndime)
    real(rp)                                        :: xa,ya,za,norma,v_tmp(3),nside
    real(rp)                                        :: dummr(3),shapf(mnode),rad
    real(rp)                                        :: grafo,buofo,denpa,dista
    real(rp)                                        :: confl,sphfl,denfl,visfl,Dvg_m,Yv_surf,Yv_fluid_k
    real(rp)                                        :: Therm_fluid_k,T_fluid_k,conce_seen(nclas_pts)
    real(rp)                                        :: Pr_m, Sc_m, LK_m, xvap, w_nonf_seen
    real(rp)                                        :: tau_p,diame

    xa  = 0.0_rp
    ya  = 0.0_rp
    za  = 0.0_rp
    xx  = 0.0_rp
    
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
          call runend('PTS_HOST_ELEMENT: GOT TO 8 BYTES INTEGERS')
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
    !-$OMP                advec,parla_pts,kfl_injVelocFunction_pts)                    &
    do ilagr = 1,nlagr_pos

       ielem = host_element(ilagr)

       if( ielem > 0 ) then

          xx(1:ndime)    = particle_position(1:ndime,ilagr)
          itype          = particle_type(ilagr)
          iinj           = particle_injector(ilagr)
          shapf(1:mnode) = host_shapf(1:mnode,ilagr)
          pnode          = lnnod(ielem)
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
          call lagrtyp(klagr) % init()                     ! Initial value
          lagrtyp(klagr) % ilagr     =  ilagr + nlacc_pts  ! Particle absolute ID
          lagrtyp(klagr) % itype     =  itype              ! Type
          lagrtyp(klagr) % ielem     =  ielem              ! Host element
          lagrtyp(klagr) % ittim     =  0                  ! Time step
          lagrtyp(klagr) % kfl_exist = -5                  ! Just been injected
          lagrtyp(klagr) % coord     =  xx                 ! Initial coordinates
          lagrtyp(klagr) % t         =  cutim - dtime      ! Current time
          lagrtyp(klagr) % t_inject  =  cutim - dtime      ! Injection time
          lagrtyp(klagr) % mpi_rank  =  kfl_paral          ! My rank                

          !----------------------------------------------------------------
          !
          ! Size distribution model will go here 
          !
          !----------------------------------------------------------------
          if (present(particle_diameter)) then
             diame                   = particle_diameter(ilagr)                 ! Particle diameter
          else
             diame                   = parttyp(itype) % diame                   ! Particle diameter
          endif



          !----------------------------------------------------------------
          !
          ! Thermodynamic model
          !
          !---------------------------------------------------------------- 
          if( parttyp(itype) % kfl_therm /= 0 ) then
             !
             ! Temperature
             !
             if(      kfl_injte_pts(iinj) == 5 ) then
                lagrtyp(klagr) % tempe_k = param_tempe_pts(1,iinj)
             else if( kfl_injte_pts(iinj) == 0 ) then
                if( .not. associated(therm) ) call runend('PTS_HOST_ELEMENT: TEMPERATURE FIELD IS REQUIRED')
                lagrtyp(klagr) % tempe_k = 0.0_rp        
                do inode = 1,lnnod(ielem)
                   ipoin = lnods(inode,ielem)  
                   lagrtyp(klagr) % tempe_k = lagrtyp(klagr) % tempe_k &
                        + shapf(inode) * therm(ipoin,3)
                end do
             end if
             !
             ! Mass
             !
             call pts_thermodynamic_properties(itype, lagrtyp(klagr) % tempe_k, ielem, pnode, lnods(1:pnode,ielem), shapf, &
                  confl, sphfl, denfl, visfl, Dvg_m, Pr_m, Sc_m, LK_m, Yv_surf, Yv_fluid_k, &
                  Therm_fluid_k, T_fluid_k, xvap, w_nonf_seen, conce_seen)
             denpa     = parttyp(itype) % liq % rho               ! Particle density 
             lagrtyp(klagr) % diam_k = diame
             lagrtyp(klagr) % mass_k = physics_sphere_mass(diame,denpa)
             !
             ! Previous solutions
             !
             lagrtyp(klagr) % tempe_km1 = lagrtyp(klagr) % tempe_k
             lagrtyp(klagr) % tempe_km2 = lagrtyp(klagr) % tempe_k
             lagrtyp(klagr) % mass_km1  = lagrtyp(klagr) % mass_k
             lagrtyp(klagr) % mass_km2  = lagrtyp(klagr) % mass_k
          else
             !
             ! Constant temperature model
             !
             call ker_proper('VISCO','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)   ! Fluid viscosity
             visfl = max(zeror,dummr(1))                                           ! Fluid viscosity
             call ker_proper('DENSI','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)   ! Fluid density
             denfl = dummr(1)           
             denpa = parttyp(itype) % denpa                                        ! Particle density
          endif
          !
          ! Initial local & effective Stokes number
          !
          tau_p  = (denpa * diame * diame)/(18.0_rp*visfl)
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
          if( associated(advec) ) then
             do inode = 1,lnnod(ielem)
                ipoin = lnods(inode,ielem)  
                lagrtyp(klagr) % v_fluid_k(1:ndime) = lagrtyp(klagr) % v_fluid_k(1:ndime) &
                     + shapf(inode) * advec(1:ndime,ipoin,3)
             end do
          end if
          lagrtyp(klagr) % v_fluid_km1 = lagrtyp(klagr) % v_fluid_k
          lagrtyp(klagr) % v_fluid_km2 = lagrtyp(klagr) % v_fluid_k

          !----------------------------------------------------------------
          !
          ! Initial particle velocity
          !
          !----------------------------------------------------------------        

          if ( kfl_injVelocFunction_pts(iinj) /= 0 ) then
             call ker_space_time_function(&
                  kfl_injVelocFunction_pts(iinj),xx(1),xx(2),xx(ndime),cutim,vv)
             vv = vv * param_veloc_pts(1,iinj)
          else
             vv = param_veloc_pts(1,iinj)
          endif

          if (     kfl_injve_pts(iinj) == -1 ) then
             !
             ! Zero velocity
             !
             lagrtyp(klagr) % veloc(1:ndime) = 0.0_rp 

          else if ( kfl_injve_pts(iinj) == 0 ) then
             !
             ! Fluid velocity
             !
             lagrtyp(klagr) % veloc(1:ndime) = lagrtyp(klagr) % v_fluid_k(1:ndime)

          else if ( kfl_injve_pts(iinj) == 1 ) then
             !
             ! With respect to normal
             !
             iinj = particle_injector(ilagr)
             nx   = parla_pts(iinj,5)
             ny   = parla_pts(iinj,6)
             nz   = parla_pts(iinj,7)
             lagrtyp(klagr) % veloc(1) =  vv*nx
             lagrtyp(klagr) % veloc(2) = -vv*ny
             if( ndime == 3 ) then
                lagrtyp(klagr) % veloc(3) = vv*nz
             end if

          else if ( kfl_injve_pts(iinj) == 2 ) then
             !
             ! Gaussian velocity injector
             ! f(x) = 1/sqrt(2*pi*sgima) * exp[(x-mu)^2/(2*sigma)]
             !
             if( kfl_injla_pts(iinj) /= 4 ) call runend('PTS_HOST_ELEMENT: GAUSSIAN VELOCITY ONLY VALID IF CIRCLE INJECTOR')
             xc    = parla_pts(iinj,1)
             yc    = parla_pts(iinj,2)
             zc    = parla_pts(iinj,3)
             nx    = parla_pts(iinj,5)
             ny    = parla_pts(iinj,6)
             nz    = parla_pts(iinj,7)                 
             sigma = param_veloc_pts(2,iinj)
             lagrtyp(klagr) % veloc(1) = vv * nx * exp(-((xx(1)-xc)**2 + (xx(2)-yc)**2 + (xx(3)-zc)**2) / (2*sigma)) 
             lagrtyp(klagr) % veloc(2) = vv * ny * exp(-((xx(1)-xc)**2 + (xx(2)-yc)**2 + (xx(3)-zc)**2) / (2*sigma))
             if( ndime == 3 ) then
                lagrtyp(klagr) % veloc(3) = vv * nz * exp(-((xx(1)-xc)**2 + (xx(2)-yc)**2 + (xx(3)-zc)**2) / (2*sigma))
             end if

          else if ( kfl_injve_pts(iinj) == 3 ) then
             !
             ! Conic velocity normal
             !
             if(  (kfl_injla_pts(iinj) /= 4) .and. (kfl_injla_pts(iinj) /= 11)  ) then
                call runend('PTS_HOST_ELEMENT: CONIC VELOCITY ONLY VALID IF CIRCLE OR ANNULAR INJECTOR')
             end if
             sigma    = param_veloc_pts(2,iinj)*(pi/180.0_rp)
             xc       = parla_pts(iinj,1)
             yc       = parla_pts(iinj,2)
             zc       = parla_pts(iinj,3)
             rad      = parla_pts(iinj,4)
             if (kfl_injla_pts(iinj) == 4) then
                !
                ! Circle injector: on whole surface of circle
                !
                nside = parla_pts(iinj, 8)
                nn(1) = parla_pts(iinj, 5)
                nn(2) = parla_pts(iinj, 6)
                nn(3) = parla_pts(iinj, 7)
             elseif (kfl_injla_pts(iinj) == 11) then
                !
                ! Annular injector: on surface of an annulus
                !
                nside = parla_pts(iinj, 9)
                nn(1) = parla_pts(iinj, 6)
                nn(2) = parla_pts(iinj, 7)
                nn(3) = parla_pts(iinj, 8)
             endif
             basis(1:ndime, 1) = nn(1:3)

             if( ndime == 3 ) then    

                xa       = xc + rad*cos(sigma)/sin(sigma)*nn(1)
                ya       = yc + rad*cos(sigma)/sin(sigma)*nn(2)
                za       = zc + rad*cos(sigma)/sin(sigma)*nn(3)
                v_tmp(1) = xa - xx(1)
                v_tmp(2) = ya - xx(2)
                v_tmp(3) = za - xx(3)
                norma    = v_tmp(1)**2 +  v_tmp(2)**2 + v_tmp(3)**2
                norma    = 1.0_rp / sqrt(norma)
                v_tmp    = norma * v_tmp

                lagrtyp(klagr) % veloc(1:3) = - vv * v_tmp(1:3)
             end if

          else if ( kfl_injve_pts(iinj) == 4 ) then
             !
             ! Spray velocity normal
             !
             if(  (kfl_injla_pts(iinj) /= 4) .and. (kfl_injla_pts(iinj) /= 11)  ) then
                call runend('PTS_HOST_ELEMENT: CONIC VELOCITY ONLY VALID IF CIRCLE OR ANNULAR INJECTOR')
             end if
             sigma            =  param_veloc_pts(2,iinj)*(pi/180.0_rp)
             xc               =  parla_pts(iinj,1)
             yc               =  parla_pts(iinj,2)
             zc               =  parla_pts(iinj,3)
             rad              =  parla_pts(iinj,4)
             if (kfl_injla_pts(iinj) == 4) then
                !
                ! Circle injector: on whole surface of circle
                !
                nside = parla_pts(iinj, 8)
                nn(1) = parla_pts(iinj, 5)
                nn(2) = parla_pts(iinj, 6)
                nn(3) = parla_pts(iinj, 7)
             elseif (kfl_injla_pts(iinj) == 11) then
                !
                ! Annular injector: on surface of an annulus
                !
                nside = parla_pts(iinj, 9)
                nn(1) = parla_pts(iinj, 6)
                nn(2) = parla_pts(iinj, 7)
                nn(3) = parla_pts(iinj, 8)
             endif
             basis(1:ndime, 1) = nn(1:3)

             if( ndime == 3 ) then    

                xa       = xc + rad*cos(sigma)/sin(sigma)*nn(1)
                ya       = yc + rad*cos(sigma)/sin(sigma)*nn(2)
                za       = zc + rad*cos(sigma)/sin(sigma)*nn(3)
                v_tmp(1) = xa - xx(1)
                v_tmp(2) = ya - xx(2)
                v_tmp(3) = za - xx(3)

                norma    = v_tmp(1)**2 +  v_tmp(2)**2 + v_tmp(3)**2
                norma    = 1.0_rp / sqrt(norma)
                v_tmp    = norma * v_tmp
                lagrtyp(klagr) % veloc(1:3) = - vv * v_tmp(1:3)
             end if

          else if ( kfl_injve_pts(iinj) == 5 ) then
             !
             ! Constant velocity 
             !
             lagrtyp(klagr) % veloc(1) = vv           
             lagrtyp(klagr) % veloc(2) = param_veloc_pts(2,iinj) 
             if( ndime == 3 ) then
                lagrtyp(klagr) % veloc(3) = param_veloc_pts(3,iinj) 
             end if

          end if

          !----------------------------------------------------------------
          !
          ! Initial acceleration
          !
          !----------------------------------------------------------------

          if( parttyp(itype) % kfl_modla == 2 ) then
             ! AB calculated before call ker_proper('DENSI','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)   ! Fluid density
             ! AB calculated before denfl = dummr(1)           
             grafo = real( parttyp(itype) % kfl_grafo, rp )                        ! Gravity  force = 1.0
             buofo = real( parttyp(itype) % kfl_buofo, rp )                        ! Buoyancy force = 1.0  
             lagrtyp(klagr) % accel(1:ndime) = grnor * gravi(1:ndime) * ( grafo - buofo * denfl / denpa )
          end if

          !
          ! Exact solution
          !
          if( kfl_exacs_pts /= 0 ) then
             call pts_exacso(2_ip,0.0_rp,lagrtyp(klagr) % accel,lagrtyp(klagr) % veloc,lagrtyp(klagr) % coord)
          end if
          !
          ! Injected particles receive default value of properties from particle type
          !
          lagrtyp(klagr) % prope(1:mlapr) = parttyp(itype) % prope(1:mlapr)
       end if
    end do
    !-$OMP END PARALLEL DO

  end subroutine pts_initial_condition

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-06
  !> @brief   Find host elements
  !> @details Find the host elements for a list of particles
  !> 
  !-----------------------------------------------------------------------
  
  subroutine pts_host_element(nlagr_pos,nlagr_new,particle_position,host_shapf,host_element)

    integer(ip),           intent(in)    :: nlagr_pos               !< Number of injected particles
    integer(ip),           intent(out)   :: nlagr_new               !< Number of new particles owned by myself
    real(rp),     pointer, intent(in)    :: particle_position(:,:)  !< (ndime,nlagr_pos)
    real(rp),     pointer, intent(inout) :: host_shapf(:,:)
    integer(ip),  pointer, intent(inout) :: host_element(:)
 
    integer(ip)                          :: ielem,ilagr
    real(rp)                             :: coloc(3),deriv(ndime,mnode),xx(3)
    real(rp)                             :: dummr(3),relse_sav,dista

    !----------------------------------------------------------------------
    ! 
    ! Look for host elements
    !
    !----------------------------------------------------------------------

    relse_sav   = relse(1)
    relse(1)    = 0.0_rp
    nlagr_new   = 0
    xx          = 0.0_rp

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

  end subroutine pts_host_element

end module mod_pts_host_element
!> @}
