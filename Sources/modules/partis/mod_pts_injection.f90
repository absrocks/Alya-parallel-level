!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    mod_pts_injection.f90
!> @author  bsc21240
!> @date    2018-03-08
!> @brief   Module for injection
!> @details Inject particles
!> For the circle injector you can add two more options:
!> DISTRIBUTION: UNICA|UNIPO  (uniform spatially distrib. particles in cartesian coordinates 
!> or polar coordinates[this one is default previous implementation])
!> NPARTICLES: ASIS -- by default the number of particles provided in the injector parameters is squared, cubed or 
!> transformed in some other way. This parameter disables these transformations (only for CIRCLE)
!-----------------------------------------------------------------------

module mod_pts_injection
    use def_parame
    use def_master
    use def_kermod
    use def_partis
    use def_domain
    use mod_memory,           only : memory_alloca
    use mod_memory,           only : memory_deallo
    use mod_memory,           only : memory_copy
    use mod_communications,   only : PAR_BROADCAST
    use mod_parall,           only : PAR_COMM_MY_CODE_WM4
    use mod_messages,         only : livinf
    use mod_messages,         only : messages_live
    use mod_pts_host_element, only : pts_host_element_and_initialize
    use mod_strings,          only : integer_to_string
    implicit none

    private

    public :: pts_injection_initialization
    public :: pts_injection_injectors
    public :: pts_injection_parallelization

contains

    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21240
    !> @date    2018-03-08
    !> @brief   Initialize injectors
    !> @details Initialize injectors
    !> 
    !-----------------------------------------------------------------------

  subroutine pts_injection_parallelization()

    integer(ip) :: iinj

    do iinj = 1, size(kfl_injla_pts)
       !
       ! Only if injection is read from file
       !
       if (kfl_injla_pts(iinj) == 10) then
          call PAR_BROADCAST(injection_pts(iinj) % number_particles)
          if (injection_pts(iinj) % number_particles > 0) then
             call pts_injection_memory(iinj)
             call PAR_BROADCAST(injection_pts(iinj) % coord_particles)
             if (IMASTER) &
                  call memory_deallo(mem_modul(1:2, modul), 'COORD', 'pts_injection_memory', injection_pts(iinj) % coord_particles)
          end if
       end if
    end do

  end subroutine pts_injection_parallelization

    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21240
    !> @date    2018-03-08
    !> @brief   Allocate memory
    !> @details Allocate memory
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_memory(iinj, particle_place, particle_injector)

        integer(ip), intent(in) :: iinj !< Injector number
        integer(ip), intent(in), pointer, optional :: particle_place(:) !< Particle place
        integer(ip), intent(in), pointer, optional :: particle_injector(:) !< Particle place
        integer(ip) :: nn, ii, jj
        real(rp), pointer :: coord_sav(:)

        if (INOTMASTER) then
            nullify(coord_sav)
            !
            ! Only if injection is read from file
            !
            if (kfl_injla_pts(iinj) == 10) then
                nn = injection_pts(iinj) % number_particles
                if (nn > 0) then
                    call memory_alloca(mem_modul(1:2, modul), 'COORD', 'pts_injection_memory', injection_pts(iinj) % coord_particles, ndime * nn)
                    if (present(particle_place) .and. present(particle_injector)) then
                        injection_pts(iinj) % number_particles = count(particle_place > 0 .and. particle_injector == iinj)
                        call memory_copy(mem_modul(1:2, modul), 'COORD', 'pts_injection_memory', injection_pts(iinj) % coord_particles, coord_sav)
                        call memory_alloca(mem_modul(1:2, modul), 'COORD', 'pts_injection_memory', injection_pts(iinj) % coord_particles, ndime * nn)
                        jj = 0
                        do ii = 1, injection_pts(iinj) % number_particles
                            if (particle_place(ii) > 0) then
                                jj = jj + 1
                                injection_pts(iinj) % coord_particles(1:ndime * jj) = coord_sav(1:ndime * ii)
                            end if
                        end do
                        call memory_deallo(mem_modul(1:2, modul), 'COORD_SAV', 'pts_injection_memory', coord_sav)
                    end if
                end if
            end if
        end if

    end subroutine pts_injection_memory

    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21240
    !> @date    2018-03-08
    !> @brief   Initialize injectors
    !> @details Initialize injectors
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_initialization()
        implicit none
        integer(ip) :: ii

        injection_pts(:) % number_particles = 0
        do ii = 1, mpala
            nullify(injection_pts(ii) % coord_particles)
            nullify(injection_pts(ii) % parameters)
        end do

    end subroutine pts_injection_initialization


    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21304
    !> @date    2019-11-08
    !> @brief   Get particle size distribution
    !> @details Based on the mass to be injected, get a particle size
    !>          dictribution, that fits a given cummulative density function 
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_distribution(iinj, mass, denpa, len_list, d_list, mass_rem)
       use mod_random,          only : random_generate_number
       use mod_physics,         only : physics_sphere_mass
       use mod_physics,         only : physics_sphere_diameter
       use def_master,          only : parii, npari, nparr, &
                                       nparc, parin, parre, &
                                       parch, nparc, mem_modul, & 
                                       modul, ISLAVE, IMASTER

       use mod_memory,          only : memory_alloca
       use mod_memory,          only : memory_deallo
       use mod_communications,  only : PAR_BROADCAST
       use mod_communications,  only : PAR_EXCHANGE
       implicit none

       integer(ip),       intent(in)  :: iinj
       real(rp),          intent(in)  :: mass
       real(rp),          intent(in)  :: denpa
       integer(ip),       intent(out) :: len_list
       real(rp), pointer, intent(out) :: d_list(:)
       real(rp),          intent(out) :: mass_rem
       real(rp)                       :: U1, diame, mass_i, mass_all
       real(rp), pointer              :: d_list_swp(:)
       integer(ip)                    :: ii, len_old, jj, itype, isafety
       logical                        :: list_full

       real(rp)                       :: d_avg, d_std
       real(rp)                       :: rr_dbar, rr_K, rr_n
       real(rp)                       :: dmax, dmin
       real(rp), save                 :: d_next(pts_minj) = 0.0_rp
                  

       itype = max(1_ip,kfl_injty_pts(iinj))
       if ( INOTSLAVE ) then

          !
          ! Initialize
          !
          len_list = 100_ip
          nullify(d_list)
          call memory_alloca(mem_modul(1:2, modul), 'D_LIST', 'pts_injection_distribution', d_list, len_list)

          ii             = 0_ip
          isafety        = 0_ip
          mass_all       = 0.0_rp
          mass_rem       = mass
          list_full= .false.
          if (mass == 0.0_rp) list_full = .true.

          do while( .not. list_full)
             !
             ! Break if for some reason too many particles would be injected
             !
             isafety = isafety + 1_ip
             if ( isafety > 99999 ) then
                list_full = .true.
                print '(A)', '--| ALYA WARNING: Number of injected particles reached limit.'
             endif
             !
             ! Select diameter 
             !
             if (d_next(iinj) > zeror) then
                !
                ! Last selection was stored in d_next
                !
                diame        = d_next(iinj)
                d_next(iinj) = 0.0_rp
             else
                select case(kfl_injsd_pts(iinj)) 
                  case(0)
                     !
                     ! CONST
                     !
                     diame =  parttyp(itype) % diame
                  case(1)
                     !
                     ! UNIFORM
                     !
                     dmin    = param_sized_pts(1,iinj)
                     dmax    = param_sized_pts(2,iinj)
                     U1      = random_generate_number(UNIQUE_SEED=.true.)
                     diame   = dmin + (dmax - dmin) * U1
                  case(2)
                     !
                     ! ROSIN_RAMMLER
                     !
                     dmin    = param_sized_pts(1,iinj)
                     dmax    = param_sized_pts(2,iinj)
                     rr_dbar = param_sized_pts(3,iinj) - dmin
                     rr_n    = param_sized_pts(4,iinj)
                     rr_K    = 1.0_rp - exp( -1.0_rp*((dmax - dmin)/rr_dbar)**rr_n )
                     U1      = random_generate_number(UNIQUE_SEED=.true.)
                     diame   = dmin + rr_dbar * (-1.0_rp*log(1.0_rp - U1*rr_K))**(1.0_rp/rr_n)
                     diame   = min(dmax,max(dmin,diame))
                end select
             endif
            
             mass_i  = parttyp(itype) % n_drop * physics_sphere_mass(diame, denpa)

             if (mass <= (mass_all + mass_i + 1.0e-20_rp) ) then
                !
                ! Exceeded mass to be injected, limit particle to maximum allowable size
                !
                mass_rem      = mass - mass_all
                d_next(iinj)  = diame
                diame         = 0.0_rp 
                list_full     = .true.
             endif

             !
             ! Increment counters
             !
             if (diame > 0.0_rp) then
                mass_all = mass_all + mass_i
                ii = ii + 1_ip

                if (ii > size(d_list,KIND=ip)) then
                   !
                   ! Need to reallocate list, because it's too short
                   !
                   len_old = size(d_list,KIND=ip)  

                   !
                   ! Save list in temporary array:
                   !
                   nullify(d_list_swp)
                   call memory_alloca(mem_modul(1:2, modul), 'D_LIST_SWP', 'pts_injection_distribution', d_list_swp, len_old)
                   
                   do jj = 1,len_old
                      d_list_swp(jj) = d_list(jj)
                   enddo 

                   !
                   ! Allocate longer array
                   ! 
                   call memory_deallo(mem_modul(1:2, modul), 'D_LIST', 'pts_injection_distribution', d_list)
                   len_list = len_old + 100_ip
                   nullify(d_list)
                   call memory_alloca(mem_modul(1:2, modul), 'D_LIST', 'pts_injection_distribution', d_list, len_list)

                   do jj = 1,len_old
                      d_list(jj) = d_list_swp(jj)
                   enddo   
                endif

                d_list(ii) = diame 
             endif
          end do

          !
          ! Update real length:
          !
          len_list = ii

          
          d_avg = 0.0_rp
          do jj = 1,len_list
            d_avg = d_avg + d_list(jj)
          enddo
          d_avg = d_avg / max(real(len_list,rp),1.0_rp) 

          d_std = 0.0_rp
          do jj = 1,len_list
            d_std = d_std + (d_list(jj)-d_avg)**2
          enddo
          d_std = sqrt( d_std / max(real(len_list-1,rp),1.0_rp) )
          

          if (mass > 0.0_rp) print*, 'm, n, d_avg, d_std, d_next', mass, len_list, d_avg, d_std, d_next(iinj)
       endif

       !
       ! Communicate particle size distribution
       !
       if (IPARALL) then  
          

          nullify(parin)
          !
          ! Exchange dimensions
          !  
          do parii=1,2 
             npari=0
 
             call PAR_EXCHANGE(len_list, parin,npari,parii)           
             
             if( parii == 1 ) then
                call memory_alloca(mem_modul(1:2,modul),'PARIN','tab_par_exchange',parin,npari)
                if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
             else
                if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
             end if

          enddo
          call memory_deallo(mem_modul(1:2,modul),'PARIN','tab_par_exchange',parin)


          !
          ! Allocate for slaves
          !
          if ( ISLAVE ) then 
             nullify(d_list)
             call memory_alloca(mem_modul(1:2, modul), 'D_LIST', 'pts_injection_distribution', d_list, max(1_ip,len_list))
          endif


          !
          ! Exchange diameters
          !
          nullify(parre)
          do parii=1,2 
             nparr=0

             
             do jj = 1_ip, len_list
                call PAR_EXCHANGE(d_list(jj), parre, nparr, parii) 
             enddo
             
             if( parii == 1 ) then
                call memory_alloca(mem_modul(1:2,modul),'PARRE','tab_par_exchange',parre,nparr)
                if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
             else
                if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
             end if

          enddo
          call memory_deallo(mem_modul(1:2,modul),'PARRE','tab_par_exchange',parre)
       endif

    end subroutine pts_injection_distribution


    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21240
    !> @date    2018-03-08
    !> @brief   Inject particles
    !> @details Inject particles according to the different injectors
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_injectors()

      use mod_communications, only : PAR_SUM
      use mod_communications, only : PAR_MAX
      use mod_communications, only : PAR_POINT_TO_POINT_ARRAY_OPERATION
      use mod_maths,          only : maths_vector_from_new_basis
      use mod_maths,          only : maths_local_orthonormal_basis
      use mod_ker_timeline,   only : ker_timeline
      use mod_random,         only : random_generate_number
      use mod_physics,        only : physics_set_liquid_temperature
      use mod_physics,        only : physics_sphere_mass
      use mod_ker_space_time_function

      implicit none
      integer(ip)          :: ii, jj, kk, iinj
      integer(ip)          :: i1, i2, i3, i4, i5, i6, j1, j2, j3
      real(rp)             :: zfact, tfact, rfact
      real(rp)             :: xx(3), xmini, xmaxi, ymini, ymaxi, zmini, zmaxi, nn(3), basis(ndime, ndime)
      real(rp)             :: x1, y1, z1, x2, y2, z2, x3, y3, z3, p1(3), p2(3), p3(3)
      real(rp)             :: ux, uy, uz, wx, wy, wz, d1, d2, d3, U1, U2, U3, U4
      real(rp)             :: t_aux, t, z, r, radius, rad, rad_min, xe, ye, ze, rr, hh, alpha, rr_aux
      real(rp)             :: mass, cutla_loc_pts, denpa
      integer(ip)          :: nside, nlagr_new, nlagr_pos, nlagr_pos2, itype,dummi
      integer(ip)          :: number_injected_pts, num_types_injected, num_particles_per_injection
      
      real(rp),    pointer :: particle_position(:,:)
      integer(ip), pointer :: particle_injector(:)
      integer(ip), pointer :: particle_type(:)
      integer(ip), pointer :: particle_place(:)

      type(r1p),   pointer :: particle_diameter_inj(:)
      real(rp),    pointer :: particle_diameter(:)
      integer(ip), pointer :: n_inj(:)

      real(rp),    save    :: mass_remaining = 0.0_rp

      !
      ! By default NO particle is injected
      !
      kfl_injec = 0
      nlagr_pos = 0
      nlagr_pos2 = 0
      number_injected_pts = 0
      cutla_loc_pts = 0.0_rp
      !
      ! Decide if particles should be injected:
      !
      if (cutim >= tinla_pts) then
         cutla_pts = cutla_pts + dtime ! accumulate time since last injection
         if (cutla_pts >= tpela_pts - zeror) then
            if (cutla_pts > 1.0e11_rp) then
               cutla_loc_pts = dtime
            else
               cutla_loc_pts = cutla_pts
            endif
            kfl_injec = 1
            cutla_pts = 0.0_rp
         end if
      end if
      if (cutim >= tfila_pts) then
         kfl_injec = 0
      end if

      !
      ! Break execution if the injection period is used up.
      !
      if (kfl_injec == 0) return


      call ker_timeline('INI_ASSEMBLY')
      !
      ! Nullify pointers
      !
      nullify(particle_position)
      nullify(particle_injector)
      nullify(particle_type)
      nullify(particle_place)
      nullify(particle_diameter)
      nullify(particle_diameter_inj)
      nullify(n_inj)
      call memory_alloca(mem_modul(1:2,modul),'PARTICLE_DIAMETER_INJ' ,'pts_injection_injectors',particle_diameter_inj ,size(kfl_injla_pts,KIND=ip))
      call memory_alloca(mem_modul(1:2,modul),'N_INJ'                 ,'pts_injection_injectors',n_inj ,size(kfl_injla_pts,KIND=ip))
      !
      ! Number of injected particles when injection
      !
      do iinj = 1, size(kfl_injla_pts)
         n_inj(iinj) = 0
         num_types_injected = 0
         if (kfl_injla_pts(iinj) /= 0) then
            if (kfl_injty_pts(iinj) == 0) then
               num_types_injected = number_types_pts
            else
               num_types_injected = 1
            end if
         end if

         !
         ! Particle density
         !
         itype = max(1_ip,kfl_injty_pts(iinj))
         if( parttyp(itype) % kfl_therm /= 0 ) then
            call physics_set_liquid_temperature( parttyp(itype) % liq , param_tempe_pts(1,iinj))
            denpa = parttyp(itype) % liq % rho
         else
            denpa = parttyp(itype) % denpa  
         endif
         !
         ! Number of injected particles is determined by mass flow rate and particle size distribution
         !
         if (injector_npts_asis(iinj) >= 2) then
            !
            ! Flow rate imposed by space & time function 
            !
            if ( kfl_injMassFunction_pts(iinj) > 0 ) then
               xx = 0.0_rp    !! Time dependent function only
               call ker_space_time_function(&
                          kfl_injMassFunction_pts(iinj),xx(1),xx(2),xx(ndime),cutim,param_mass_pts(1, iinj))
            end if

            mass           = cutla_loc_pts * param_mass_pts(1, iinj) + mass_remaining
            mass_remaining = 0.0_rp

            call pts_injection_distribution(iinj, mass, denpa, nside, particle_diameter_inj(iinj)%a, mass_remaining )
            n_inj(iinj) = nside
         !
         ! Injection based on the specification of number of particles 
         !
         else

            if (kfl_injla_pts(iinj) == 1) then
               !
               ! Square
               !
               nside = int(parla_pts(iinj, 7), ip)
               n_inj(iinj) = nside * nside

            else if (kfl_injla_pts(iinj) == 2) then
               !
               ! Sphere injector
               !
               nside = int(parla_pts(iinj, 5), ip)
               n_inj(iinj) =  nside ** ndime

            else if (kfl_injla_pts(iinj) == 3) then
               !
               ! Semi-Sphere injector
               !
               nside = int(parla_pts(iinj, 8), ip)
               n_inj(iinj) =  2 * nside ** ndime
            else if (kfl_injla_pts(iinj) == 4) then
               !
               ! Circle injector
               !
               nside = int(parla_pts(iinj, 8), ip)

               if (injector_npts_asis(iinj) == 1) then !inject the requested number of particles
                  n_inj(iinj) =  nside
               else
                  n_inj(iinj) =  nside * nside
               end if

            else if (kfl_injla_pts(iinj) == 5) then
               !
               ! Rectangle injector
               !
               if (ndime == 3) then
                  nside = int(parla_pts(iinj, 10), ip)
               n_inj(iinj) =  nside * nside
               end if

            else if (kfl_injla_pts(iinj) == 6) then
               !
               ! Pointwise injector
               !
               nside = max(1_ip, int(parla_pts(iinj, 4), ip))
               n_inj(iinj) =  nside

            else if (kfl_injla_pts(iinj) == 7) then
               !
               ! Cone injector
               !
               nside = int(parla_pts(iinj, 9), ip)
               n_inj(iinj) =  nside ** 3

            else if (kfl_injla_pts(iinj) == 8) then
               !
               ! Random injector in a cube/square box
               !
               nside = int(parla_pts(iinj, 7), ip)
               n_inj(iinj) =  nside

            else if (kfl_injla_pts(iinj) == 9) then
               !
               ! Segment injector
               !
               nside = int(parla_pts(iinj, 7), ip)
               n_inj(iinj) =   nside

            else if (kfl_injla_pts(iinj) == 10) then
               !
               ! Read from file
               !
               n_inj(iinj) =  injection_pts(iinj) % number_particles
            else if (kfl_injla_pts(iinj) == 11) then
               !
               ! Annulus injector
               !
               nside = int(parla_pts(iinj, 9), ip)

               if (injector_npts_asis(iinj) == 1) then !inject the requested number of particles
                  n_inj(iinj) =  nside
               else
                  n_inj(iinj) = nside * nside 
               end if

            end if
            !!
            !! Allocate size distribution
            !!
            !nullify(particle_diameter_inj(iinj)%a)
            !call memory_alloca(mem_modul(1:2,modul),'PARTICLE_DIAMETER_INJ' ,'pts_injection_injectors',particle_diameter_inj(iinj)%a, max(1_ip,n_inj(iinj)))
            
            !! if (n_inj(iinj) > 0) then
            !!    mass = physics_sphere_mass(parttyp(itype) % diame,denpa) * n_inj(iinj)   
            !! else
            !!    mass = 0.0_rp
            !! endif


            nullify(particle_diameter_inj(iinj)%a)
            call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_DIAMETER_INJ % A', 'pts_injection_injectors', particle_diameter_inj(iinj)%a, max(1_ip,n_inj(iinj))* num_types_injected)

            if ( num_types_injected == 1 ) then
               itype = max(1_ip,kfl_injty_pts(iinj))
               do jj = 1, n_inj(iinj) * num_types_injected
                  particle_diameter_inj(iinj) % a(jj) = parttyp(itype) % diame
               enddo
            else
               do jj = 1, n_inj(iinj) * num_types_injected
                  itype = mod(jj,num_types_injected) + 1_ip
                  particle_diameter_inj(iinj) % a(jj) = parttyp(itype) % diame
               enddo
            endif

         endif
         !
         ! Total number of injected particles
         !
         number_injected_pts = number_injected_pts + n_inj(iinj) * number_types_pts

      end do

      if (INOTMASTER) then
         !
         ! Every subdomain does the same
         ! Later it is decided to which subdomain do they belong.
         !
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_POSITION', 'pts_injection_injectors', particle_position, ndime, max(1_ip,number_injected_pts))
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_INJECTOR', 'pts_injection_injectors', particle_injector, max(1_ip,number_injected_pts))
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_TYPE',     'pts_injection_injectors', particle_type,     max(1_ip,number_injected_pts))
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_PLACE'   , 'pts_injection_injectors', particle_place,    max(1_ip,number_injected_pts))
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_DIAMETER', 'pts_injection_injectors', particle_diameter, max(1_ip,number_injected_pts))

         do iinj = 1, size(kfl_injla_pts)

            if (kfl_injty_pts(iinj) == 0) then
               num_types_injected = number_types_pts
            else
               num_types_injected = 1
            end if

            !
            ! Particle diameter from individual injectors
            !
            do jj = 1, n_inj(iinj) * num_types_injected
               nlagr_pos2 = nlagr_pos2 + 1
               particle_diameter(nlagr_pos2) = particle_diameter_inj(iinj)%a(jj)
            enddo


            if (kfl_injla_pts(iinj) == 1) then

               !----------------------------------------------------------
               !
               ! Square injector given the corners
               !
               !----------------------------------------------------------

               nside = int(parla_pts(iinj, 7), ip)
               xmini = parla_pts(iinj, 1)
               ymini = parla_pts(iinj, 2)
               zmini = parla_pts(iinj, 3)
               xmaxi = parla_pts(iinj, 4)
               ymaxi = parla_pts(iinj, 5)
               zmaxi = parla_pts(iinj, 6)

               if (     zmini == zmaxi) then
                  i1 = 1
                  i2 = 2
                  i3 = 3
               else if (xmini == xmaxi) then
                  i1 = 3
                  i2 = 2
                  i3 = 1
               else if (ymini == ymaxi) then
                  i1 = 1
                  i2 = 3
                  i3 = 2
               else
                  call runend('LAGRAN: WRONG BOX')
               end if
               i4 = i1 + 3_ip
               i5 = i2 + 3_ip
               i6 = i3 + 3_ip
               j1 = i1
               j2 = i2
               j3 = i3

               if (xmini == xmaxi .or. ymini == ymaxi .or. zmini == zmaxi) then

                  xx(j3) = parla_pts(iinj, i3)

                  do jj = 1, nside

                     if (nside == 1) then
                        xx(j2) = parla_pts(iinj, i2)
                     else
                        xx(j2) = real(jj - 1, rp)/real(nside - 1, rp)*(parla_pts(iinj, i5) - parla_pts(iinj, i2)) + parla_pts(iinj, i2)
                     end if

                     do ii = 1, nside

                        if (nside == 1) then
                           xx(j1) = parla_pts(iinj, i1)
                        else
                           xx(j1) = real(ii - 1, rp)/real(nside - 1, rp)*(parla_pts(iinj, i4) - parla_pts(iinj, i1)) + parla_pts(iinj, i1)
                        end if

                        if (kfl_injty_pts(iinj) > 0) then
                           nlagr_pos = nlagr_pos + 1
                           particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                           particle_injector(nlagr_pos) = iinj
                           particle_type(nlagr_pos) = kfl_injty_pts(iinj)

                        else
                           do itype = 1, ntyla_pts
                              if (parttyp(itype) % kfl_exist /= 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = itype
                              end if
                           end do
                        end if

                     end do
                  end do

               end if

            else if (kfl_injla_pts(iinj) == 8) then

               !----------------------------------------------------------
               !
               ! Random injector in a square/cubic box
               !
               !----------------------------------------------------------

               nside = int(parla_pts(iinj, 7), ip)
               xmini = parla_pts(iinj, 1)
               ymini = parla_pts(iinj, 2)
               zmini = parla_pts(iinj, 3)
               xmaxi = parla_pts(iinj, 4)
               ymaxi = parla_pts(iinj, 5)
               zmaxi = parla_pts(iinj, 6)

               do ii = 1, n_inj(iinj)
                  U1 = random_generate_number(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                  U2 = random_generate_number(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                  U3 = random_generate_number(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                  xx(1) = xmini + (xmaxi - xmini) * U1
                  xx(2) = ymini + (ymaxi - ymini) * U2
                  if (ndime == 3) xx(3) = zmini + (zmaxi - zmini) * U3

                  if (kfl_injty_pts(iinj) > 0) then
                     nlagr_pos = nlagr_pos + 1
                     particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                     particle_injector(nlagr_pos) = iinj
                     particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                  else
                     do itype = 1, ntyla_pts
                        if (parttyp(itype) % kfl_exist /= 0) then
                           nlagr_pos = nlagr_pos + 1
                           particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                           particle_injector(nlagr_pos) = iinj
                           particle_type(nlagr_pos) = itype
                        end if
                     end do
                  end if

               end do

            else if (kfl_injla_pts(iinj) == 5) then

               !----------------------------------------------------------
               !
               ! Rectangle injector
               !
               !----------------------------------------------------------

               nside = int(parla_pts(iinj, 10), ip)
               if (ndime == 3) then

                  x1 = parla_pts(iinj, 1)
                  y1 = parla_pts(iinj, 2)
                  z1 = parla_pts(iinj, 3)
                  x2 = parla_pts(iinj, 4)
                  y2 = parla_pts(iinj, 5)
                  z2 = parla_pts(iinj, 6)
                  x3 = parla_pts(iinj, 7)
                  y3 = parla_pts(iinj, 8)
                  z3 = parla_pts(iinj, 9)

                  d1 = sqrt((x1 - x2)**2_ip + (y1 - y2)**2_ip + (z1 - z2)**2_ip)
                  d2 = sqrt((x3 - x2)**2_ip + (y3 - y2)**2_ip + (z3 - z2)**2_ip)
                  d3 = sqrt((x3 - x1)**2_ip + (y3 - y1)**2_ip + (z3 - z1)**2_ip)
                  if (d1 - d2 >= 0 .and. d1 - d3 >= 0)then
                     p1(1) = x1
                     p1(2) = y1
                     p1(3) = z1
                     p2(1) = x3
                     p2(2) = y3
                     p2(3) = z3
                     p3(1) = x2
                     p3(2) = y2
                     p3(3) = z2
                  else if (d2 - d1 >= 0 .and. d2 - d3 >= 0)then
                     p1(1) = x2
                     p1(2) = y2
                     p1(3) = z2
                     p2(1) = x1
                     p2(2) = y1
                     p2(3) = z1
                     p3(1) = x3
                     p3(2) = y3
                     p3(3) = z3
                  else if (d3 - d1 >= 0 .and. d3 - d2 >= 0)then
                     p1(1) = x3
                     p1(2) = y3
                     p1(3) = z3
                     p2(1) = x2
                     p2(2) = y2
                     p2(3) = z2
                     p3(1) = x1
                     p3(2) = y1
                     p3(3) = z1
                  end if
                  ux = x3 - x1
                  uy = y3 - y1
                  uz = z3 - z1
                  wx = x2 - x1
                  wy = y2 - y1
                  wz = z2 - z1
                  nn(1) = uy * wz - uz * wy
                  nn(2) = wx * uz - wz * ux
                  nn(3) = ux * wy - wx * uy
                  if (nn(1) == 0.0_rp)then
                     if (nn(3) == 0.0_rp)then
                        i1 = 1
                        i2 = 3
                        i3 = 2
                     else
                        i1 = 1
                        i2 = 2
                        i3 = 3
                        if (p1(i1) == p2(i1) .or. p2(i2) == p3(i2))then
                           i1 = 2
                           i2 = 1
                           i3 = 3
                        end if
                     end if
                  else
                     i1 = 3
                     i2 = 2
                     i3 = 1
                     if (p1(i1) == p2(i1) .or. p2(i2) == p3(i2))then
                        i1 = 2
                        i2 = 3
                        i3 = 1
                     end if

                  end if

                  do ii = 1, nside
                     if (nside == 1) then
                        xx(i1) = p2(i1)
                     else
                        if ((p2(i1) - p1(i1)) >= 0)then
                           xx(i1) = real(ii - 1, rp)/real(nside - 1, rp)*(p2(i1) - p1(i1)) + p1(i1)
                        else
                           xx(i1) = real(ii - 1, rp)/real(nside - 1, rp)*(p1(i1) - p2(i1)) + p2(i1)
                        end if
                     end if
                     do jj = 1, nside
                        if (nside == 1) then
                           xx(i2) = p2(i2)
                        else
                           if ((p3(i2) - p2(i2)) >= 0)then
                              xx(i2) = real(jj - 1, rp)/real(nside - 1, rp)*(p3(i2) - p2(i2)) + p2(i2)
                           else
                              xx(i2) = real(jj - 1, rp)/real(nside - 1, rp)*(p2(i2) - p3(i2)) + p3(i2)
                           end if
                        end if
                        xx(i3) = (-(xx(i1) - p1(i1)) * nn(i1) - (xx(i2) - p1(i2)) * nn(i2))/nn(i3) + p1(i3)

                        if (kfl_injty_pts(iinj) > 0) then
                           nlagr_pos = nlagr_pos + 1
                           particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                           particle_injector(nlagr_pos) = iinj
                           particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                        else
                           do itype = 1, ntyla_pts
                              if (parttyp(itype) % kfl_exist /= 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = itype
                              end if
                           end do
                        end if

                     end do
                  end do
               end if

            else if ((kfl_injla_pts(iinj) == 4) .or. (kfl_injla_pts(iinj) == 11)) then
               !----------------------------------------------------------
               !
               ! Circle or Annular injector
               !
               !----------------------------------------------------------

               !
               ! INPUT
               !
               xc    = parla_pts(iinj, 1) ! center x
               yc    = parla_pts(iinj, 2) ! center y
               zc    = parla_pts(iinj, 3) ! center z
               rad   = parla_pts(iinj, 4) ! outter radius
               if (kfl_injla_pts(iinj) == 4) then
                   !
                   ! Circle injector: on whole surface of circle
                   !
                   nside = int(parla_pts(iinj, 8), ip)
                   rad_min = 0.0_rp 

                   if (ndime == 3) then
                       nn(1) = parla_pts(iinj, 5)
                       nn(2) = parla_pts(iinj, 6)
                       nn(3) = parla_pts(iinj, 7)
                       basis(1:ndime, 1) = nn(1:3)
                   else
                       nx = parla_pts(iinj, 5)
                       ny = parla_pts(iinj, 6)
                       nz = parla_pts(iinj, 7)
                   endif
               elseif (kfl_injla_pts(iinj) == 11) then
                   !
                   ! Annular injector: on surface of an annulus
                   !
                   nside   = int(parla_pts(iinj, 9), ip)
                   rad_min = parla_pts(iinj, 5)

                   if (ndime == 3) then
                       nn(1) = parla_pts(iinj, 6)
                       nn(2) = parla_pts(iinj, 7)
                       nn(3) = parla_pts(iinj, 8)
                       basis(1:ndime, 1) = nn(1:3)
                   else
                       nx = parla_pts(iinj, 6)
                       ny = parla_pts(iinj, 7)
                       nz = parla_pts(iinj, 8)
                   endif
               endif
               



               if (ndime == 3) then !3d mesh
                  if (kfl_random_pts(iinj) == 1) then

                     if (injector_particle_distribution(iinj) == 0) then !uniform points in polar coordinates
                        call livinf(-7_ip, "Injecting particles unformly distributed in polar coordinates for injector ", iinj)
                        do jj = 1, nside
                           do kk = 1, nside
                              U1 = random_generate_number(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                              U2 = random_generate_number(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                              t = 2.0_rp * pi * U2
                              radius = rad_min + (rad-rad_min) * U1   ! Don't worry, rad_min is 0.0 for circle

                              xx(1) = 0.0_rp
                              xx(2) = radius * cos(t)
                              xx(3) = radius * sin(t)

                              call maths_local_orthonormal_basis(ndime, basis)
                              !call maths_vector_to_new_basis(ndime, basis, xx)
                              call maths_vector_from_new_basis(ndime,basis,xx)

                              xx(1) = xx(1) + xc
                              xx(2) = xx(2) + yc
                              xx(3) = xx(3) + zc

                              if (kfl_injty_pts(iinj) > 0) then
                                 nlagr_pos                             = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos)          = iinj
                                 particle_type(nlagr_pos)              = kfl_injty_pts(iinj)
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos                             = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos)          = iinj
                                       particle_type(nlagr_pos)              = itype
                                    end if
                                 end do
                              end if

                           end do
                        end do

                     else if (injector_particle_distribution(iinj) == 1) then !uniform points in cartesian coordinates

                        call livinf(-7_ip, "Injecting particles unformly distributed in cartesian coordinates for injector ", iinj)

                        do jj = 1, n_inj(iinj) 

                           U1     = random_generate_number(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4) 
                           U2     = random_generate_number(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4) 
                           
                           !
                           ! Take a random angle
                           !
                           t      = 2.0_rp * pi * U2

                           !
                           ! The idea is to give the outter radii higher
                           ! probability, becuse they have a bigger area.
                           ! The area in a drdt segment is exactly r*dr*dt.
                           ! U1 is the probability that the particle is found
                           ! within an area: 
                           ! A(radius) / Atot = (radius**2-rad_min**2) / (rad**2-rad_min**2) = U1
                           !                     -
                           !             -       |                       
                           !      -      |        |                      
                           ! o    |       |       |                      
                           !      -      |        |                      
                           !             -       |                       
                           !                     -                       
                           ! For just the circle, see:  
                           ! !https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
                           
                           radius = sqrt( rad_min**2 + (rad**2-rad_min**2) * U1 )

                           xx(1) = 0.0_rp
                           xx(2) = radius * cos(t)
                           xx(3) = radius * sin(t)

                           call maths_local_orthonormal_basis(ndime,basis)
                           !!!call maths_vector_to_new_basis(ndime,basis,xx)
                           call maths_vector_from_new_basis(ndime,basis,xx)

                           xx(1) = xx(1) + xc
                           xx(2) = xx(2) + yc
                           xx(3) = xx(3) + zc

                           if( kfl_injty_pts(iinj) > 0 ) then
                              nlagr_pos                             = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos)          = iinj
                              particle_type(nlagr_pos)              = kfl_injty_pts(iinj)
                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos                             = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos)          = iinj
                                    particle_type(nlagr_pos)              = itype
                                 end if
                              end do
                           end if

                        end do


                     end if

                  else  !not kfl_random_pts(iinj) == 1

                     tfact = 2.0_rp * pi / real(nside, rp)
                     rfact = (rad-rad_min) / real(nside, rp)

                     if (kfl_injla_pts(iinj) == 11) then
                        call runend('mod_pts_injection: deterministic particle palcement is not implemented for annulus.')
                     endif


                     if(kfl_random_pts(iinj)==2 )then  !PIPE=Injected following mass flow rate
                        tfact = 2.0_rp * pi / real(10, rp)
                        rfact = rad / real(10, rp)
                        num_particles_per_injection =  nside / 10  !!!nside siempre es divisible entre 10!!!!
                        do kk = 1,10
                           do ii = 1, num_particles_per_injection
                              rr = (real(kk, rp) - 1.0_rp) * rfact
                              t = (real(kk, rp) - 1.0_rp) * tfact
                              radius = (rr + rfact * (real(kk, rp) - 1.0_rp)/num_particles_per_injection) * rad

                              xx(1) = 0.0_rp
                              xx(2) = radius**0.5_rp * cos(t) 
                              xx(3) = radius**0.5_rp * sin(t) 
                              call maths_local_orthonormal_basis(ndime, basis)
                             !!! call maths_vector_to_new_basis(ndime, basis, xx)
                              call maths_vector_from_new_basis(ndime,basis,xx)

                              xx(1) = xx(1) + xc
                              xx(2) = xx(2) + yc
                              xx(3) = xx(3) + zc
                              if (kfl_injty_pts(iinj) > 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos) = iinj
                                       particle_type(nlagr_pos) = itype
                                    end if
                                 end do
                              end if
                           end do
                        end do
                     else
                        do jj = 1, nside
                           rr = (real(jj, rp) - 1.0_rp) * rfact
                           do kk = 1, nside
                              t = (real(kk, rp) - 1.0_rp) * tfact
                              radius = (rr + rfact * (real(kk, rp) - 1.0_rp)/nside) * rad

                              xx(1) = 0.0_rp
                              xx(2) = radius**0.5_rp * cos(t)
                              xx(3) = radius**0.5_rp * sin(t)

                              call maths_local_orthonormal_basis(ndime, basis)
                              !!!!call maths_vector_to_new_basis(ndime, basis, xx)
                              call maths_vector_from_new_basis(ndime,basis,xx)

                              xx(1) = xx(1) + xc
                              xx(2) = xx(2) + yc
                              xx(3) = xx(3) + zc

                              if (kfl_injty_pts(iinj) > 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos) = iinj
                                       particle_type(nlagr_pos) = itype
                                    end if
                                 end do
                              end if

                           end do
                        end do
                     end if
                  end if
               else !not 3d mesh
                  zfact = 2.0_rp / real(nside + 1, rp)
                  tfact = 2.0_rp * pi / real(nside + 1, rp)
                  rfact = rad / real(nside, rp)

                  if (kfl_injla_pts(iinj) == 11) then
                     call runend('mod_pts_injection: 2D model is not implemented for annulus.')
                  endif

                  if (kfl_random_pts(iinj) == 1) then
                     do jj = 1, nside
                        t_aux = real(jj, rp) * tfact
                        do kk = 1, nside
                           U1 = random_generate_number(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                           U2 = random_generate_number(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                           t = t_aux * U1
                           radius = real(kk, rp) * rfact * U2
                           xx(1) = radius * cos(t) + xc
                           xx(2) = radius * sin(t) + yc

                           if (kfl_injty_pts(iinj) > 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = itype
                                 end if
                              end do
                           end if

                        end do
                     end do
                  else if(kfl_random_pts(iinj)==2 )then  !PIPE=Injected following mass flow rate
                     tfact = 2.0_rp * pi / real(10, rp)
                     rfact = rad / real(10, rp)
                     num_particles_per_injection =  nside / 10  !!!nside siempre es divisible entre 10!!!!
                     do kk = 1,10
                        do ii = 1, num_particles_per_injection
                           rr = real(kk, rp) * rfact
                           t =  real(kk, rp) * tfact
                           radius = rr
                           xx(1) = radius * cos(t) + xc
                           xx(2) = radius * sin(t) + yc
                          

                           if (kfl_injty_pts(iinj) > 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = itype
                                 end if
                              end do
                           end if
                        end do
                     end do
                  else
                     do jj = 1, nside
                        t = real(jj, rp) * tfact
                        do kk = 1, nside
                           radius = real(kk, rp) * rfact
                           xx(1) = radius * cos(t) + xc
                           xx(2) = radius * sin(t) + yc

                           if (kfl_injty_pts(iinj) > 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = itype
                                 end if
                              end do
                           end if

                        end do
                     end do

                  end if
               end if

            else if (kfl_injla_pts(iinj) == 2) then

               !----------------------------------------------------------
               !
               ! Sphere injector
               !
               !----------------------------------------------------------

               nside = int(parla_pts(iinj, 5), ip)

               if (ndime == 3) then

                  xc = parla_pts(iinj, 1)
                  yc = parla_pts(iinj, 2)
                  zc = parla_pts(iinj, 3)
                  rad = parla_pts(iinj, 4)

                  zfact = 2.0_rp / real(nside + 1, rp)
                  tfact = 2.0_rp * pi / real(nside + 1, rp)
                  rfact = rad / real(nside, rp)

                  do ii = 1, nside

                     z = -1.0_rp + real(ii, rp) * zfact
                     r = sqrt(1.0_rp - z * z)

                     do jj = 1, nside

                        t = real(jj, rp) * tfact

                        do kk = 1, nside

                           radius = real(kk, rp) * rfact

                           xx(1) = radius * r * cos(t) + xc
                           xx(2) = radius * r * sin(t) + yc
                           xx(3) = radius * z + zc

                           if (kfl_injty_pts(iinj) > 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = itype
                                 end if
                              end do
                           end if

                        end do

                     end do
                  end do

               else

                  xc = parla_pts(iinj, 1)
                  yc = parla_pts(iinj, 2)
                  zc = parla_pts(iinj, 3)
                  rad = parla_pts(iinj, 4)

                  zfact = 2.0_rp / real(nside + 1, rp)
                  tfact = 2.0_rp * pi / real(nside + 1, rp)
                  rfact = rad / real(nside, rp)

                  do jj = 1, nside

                     t = real(jj, rp) * tfact

                     do kk = 1, nside

                        radius = real(kk, rp) * rfact
                        xx(1) = radius * cos(t) + xc
                        xx(2) = radius * sin(t) + yc

                        if (kfl_injty_pts(iinj) > 0) then
                           nlagr_pos = nlagr_pos + 1
                           particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                           particle_injector(nlagr_pos) = iinj
                           particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                        else
                           do itype = 1, ntyla_pts
                              if (parttyp(itype) % kfl_exist /= 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = itype
                              end if
                           end do
                        end if

                     end do

                  end do

               end if

            else if (kfl_injla_pts(iinj) == 3) then

               !----------------------------------------------------------
               !
               ! Semi-Sphere injector
               !
               !----------------------------------------------------------

               nside = int(parla_pts(iinj, 8), ip)

               if (ndime == 3) then

                  xc = parla_pts(iinj, 1)
                  yc = parla_pts(iinj, 2)
                  zc = parla_pts(iinj, 3)
                  rad = parla_pts(iinj, 4)
                  nx = parla_pts(iinj, 5)
                  ny = parla_pts(iinj, 6)
                  nz = parla_pts(iinj, 7)

                  zfact = 2.0_rp / real(nside + 1, rp)
                  tfact = 2.0_rp * pi / real(nside + 1, rp)
                  rfact = rad / real(nside, rp) / 2.0_rp

                  do ii = 1, nside

                     z = -1.0_rp + real(ii, rp) * zfact
                     r = sqrt(1.0_rp - z * z)

                     do jj = 1, nside

                        t = real(jj, rp) * tfact

                        do kk = 1, nside * 2

                           radius = real(kk, rp) * rfact

                           xx(1) = radius * r * cos(t) + xc
                           xx(2) = radius * r * sin(t) + yc
                           xx(3) = radius * z + zc

                           if ((xx(1) - xc) * nx + (xx(2) - yc) * ny + (xx(3) - zc) * nz >= 0.0_rp) then
                              if (kfl_injty_pts(iinj) > 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos) = iinj
                                       particle_type(nlagr_pos) = itype
                                    end if
                                 end do
                              end if

                           end if

                        end do

                     end do
                  end do

               else

                  xc = parla_pts(iinj, 1)
                  yc = parla_pts(iinj, 2)
                  zc = parla_pts(iinj, 3)
                  rad = parla_pts(iinj, 4)
                  nx = parla_pts(iinj, 5)
                  ny = parla_pts(iinj, 6)
                  nz = parla_pts(iinj, 7)

                  zfact = 2.0_rp / real(nside + 1, rp)
                  tfact = 2.0_rp * pi / real(nside + 1, rp)
                  rfact = rad / real(nside, rp) / 2.0_rp

                  do jj = 1, nside

                     t = real(jj, rp) * tfact

                     do kk = 1, nside * 2

                        radius = real(kk, rp) * rfact

                        xx(1) = radius * cos(t) + xc
                        xx(2) = radius * sin(t) + yc

                        if ((xx(1) - xc) * nx + (xx(2) - yc) * ny >= 0.0_rp) then
                           if (kfl_injty_pts(iinj) > 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = itype
                                 end if
                              end do
                           end if

                        end if

                     end do

                  end do

               end if


            else if (kfl_injla_pts(iinj) == 7) then

               !----------------------------------------------------------
               !
               ! Cone injector
               !
               !           _ _     _ _
               !    /\      |       |  z
               !   /r \     |  hh  _|_     zfact
               !  /<-  \    |       |  z
               ! /______\  _|_     _|_
               ! rad       radius = [0,rad]
               ! <--       t      = [0,2pi]
               !----------------------------------------------------------

               nside = int(parla_pts(iinj, 9), ip)

               if (ndime == 3) then

                  xc = parla_pts(iinj, 1)
                  yc = parla_pts(iinj, 2)
                  zc = parla_pts(iinj, 3)
                  rad = parla_pts(iinj, 4)
                  hh = parla_pts(iinj, 5)
                  nn(1) = parla_pts(iinj, 6)
                  nn(2) = parla_pts(iinj, 7)
                  nn(3) = parla_pts(iinj, 8)

                  basis(1:ndime, 1) = nn(1:3)
                  zfact = hh/real(nside, rp)

                  do ii = 1, nside

                     z = (-1.0_rp + real(ii, rp)) * zfact
                     r = (hh - z) / (hh/rad)
                     tfact = 2 * pi /real(nside, rp) + 1
                     rfact = r / real(nside, rp)

                     do jj = 1, nside

                        t = real(jj, rp) * tfact

                        do kk = 1, nside
                           radius = real(kk, rp) * rfact

                           xx(1) = radius * cos(t)
                           xx(2) = radius * sin(t)
                           xx(3) = z
                           ! Rotate axis in function of given normal (input)
                           call maths_local_orthonormal_basis(ndime, basis)
                           !!!call maths_vector_to_new_basis(ndime, basis, xx)
                           call maths_vector_from_new_basis(ndime,basis,xx)

                           xx(1) = xx(1) + xc
                           xx(2) = xx(2) + yc
                           xx(3) = xx(3) + zc

                           if (kfl_injty_pts(iinj) > 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = itype
                                 end if
                              end do
                           end if

                        end do

                     end do
                  end do
               else
                  call runend('PTS_INJECT: CONE INJECTOR ONLY WORKS IF 3D')
               end if

            else if (kfl_injla_pts(iinj) == 6) then

               !----------------------------------------------------------
               !
               ! Pointwise injector
               !
               !----------------------------------------------------------

               xx(1:3) = parla_pts(iinj, 1:3)

               do ii = 1, n_inj(iinj)
                  if (kfl_injty_pts(iinj) > 0) then
                     nlagr_pos = nlagr_pos + 1
                     particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                     particle_injector(nlagr_pos) = iinj
                     particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                  else
                     do itype = 1, ntyla_pts
                        if (parttyp(itype) % kfl_exist /= 0) then
                           nlagr_pos = nlagr_pos + 1
                           particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                           particle_injector(nlagr_pos) = iinj
                           particle_type(nlagr_pos) = itype
                        end if
                     end do
                  end if
               end do

            else if (kfl_injla_pts(iinj) == 9) then

               !----------------------------------------------------------
               !
               ! Segment injector
               !
               !----------------------------------------------------------

               nside = int(parla_pts(iinj, 7), ip)

               xmini = parla_pts(iinj, 1)
               ymini = parla_pts(iinj, 2)
               zmini = parla_pts(iinj, 3)
               xmaxi = parla_pts(iinj, 4)
               ymaxi = parla_pts(iinj, 5)
               zmaxi = parla_pts(iinj, 6)

               if( kfl_injty_pts(iinj) > 0 ) then
                  do ii = 1,nside
                     nlagr_pos = nlagr_pos + 1
                     particle_position(1, nlagr_pos) = xmini + real(ii-1,rp)/real(nside-1,rp)*(xmaxi - xmini)
                     particle_position(2, nlagr_pos) = ymini + real(ii-1,rp)/real(nside-1,rp)*(ymaxi - ymini)
                     if (ndime == 3) particle_position(3, nlagr_pos) = zmini + real(ii-1,rp)/real(nside-1,rp)*(zmaxi - zmini)
                     particle_injector(nlagr_pos) = iinj
                     particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                  end do
                  !print*,'hasta donde voy',iinj,size(kfl_injla_pts)
               else
                  do itype = 1, ntyla_pts
                     if (parttyp(itype) % kfl_exist /= 0) then
                        do ii = 1, nside
                           nlagr_pos = nlagr_pos + 1
                           particle_position(1, nlagr_pos) = xmini + real(ii-1,rp)/real(nside-1,rp) * (xmaxi - xmini)
                           particle_position(2, nlagr_pos) = ymini + real(ii-1,rp)/real(nside-1,rp) * (ymaxi - ymini)
                           if (ndime == 3) particle_position(3, nlagr_pos) = zmini + real(ii-1,rp)/real(nside-1,rp) * (zmaxi - zmini)
                           particle_injector(nlagr_pos) = iinj
                           particle_type(nlagr_pos) = itype
                        end do
                     end if
                  end do
               end if

            else if (kfl_injla_pts(iinj) == 10) then

               !----------------------------------------------------------
               !
               ! Read from file
               !
               !----------------------------------------------------------
                do ii = 1, n_inj(iinj)
                  if (kfl_injty_pts(iinj) > 0) then
                     nlagr_pos = nlagr_pos + 1
                     particle_position(1:ndime, nlagr_pos) = injection_pts(iinj) % coord_particles(ii * ndime - (ndime-1):ii * ndime)
                     particle_injector(nlagr_pos) = iinj
                     particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                  else
                    do itype = 1, ntyla_pts
                     if (parttyp(itype) % kfl_exist /= 0) then
                           nlagr_pos = nlagr_pos + 1
                           particle_position(1:ndime, nlagr_pos )= injection_pts(iinj) % coord_particles(ii * ndime - 2:ii * ndime)
                           particle_injector(nlagr_pos) = iinj
                           particle_type(nlagr_pos) = itype
                        end if
                     end do
                  end if
                 end do
               end if
         end do

      end if
      !
      ! Loop over particles to find host elements
      ! NLAGR_POS = Total number of injected particles
      ! NUMBER_INJECTED_PTS = Total number of particles (NLAGR_POS)
      !
      call pts_host_element_and_initialize(&
           nlagr_pos,nlagr_new,number_injected_pts,particle_position,particle_injector,&
           particle_type,particle_diameter)
      !
      ! Output injection info
      !
      call messages_live(&
           'INJECT '//integer_to_string(nlagr_new)//&
           ' LAGRANGIAN PARTICLES OUT OF '//integer_to_string(number_injected_pts),&
           'MODULE')
      if( nlagr_new < number_injected_pts ) then
         dummi = number_injected_pts - nlagr_new
         call messages_live(&
              'PARTICLES INJECTED OUT OF THE DOMAIN: '//integer_to_string(dummi),&
              'MODULE')
      end if
      !
      ! Actual total number of particles (over all partitions)
      !
      nlacc_pts = nlacc_pts + nlagr_new
      !
      ! Deallocate memory
      !
      call memory_deallo(mem_modul(1:2,modul), 'PARTICLE_POSITION'     ,'pts_injection_injectors', particle_position)
      call memory_deallo(mem_modul(1:2,modul), 'PARTICLE_INJECTOR'     ,'pts_injection_injectors', particle_injector)
      call memory_deallo(mem_modul(1:2,modul), 'PARTICLE_TYPE'         ,'pts_injection_injectors', particle_type)
      call memory_deallo(mem_modul(1:2,modul), 'PARTICLE_PLACE'        ,'pts_injection_injectors', particle_place)

      call memory_deallo(mem_modul(1:2,modul), 'PARTICLE_DIAMETER'     ,'pts_injection_injectors', particle_diameter)
      call memory_deallo(mem_modul(1:2,modul), 'PARTICLE_DIAMETER_INJ' ,'pts_injection_injectors', particle_diameter_inj)
      call memory_deallo(mem_modul(1:2,modul), 'N_INJ'                 ,'pts_injection_injectors', n_inj)
      !
      ! Output injected particles and the recover normal flag kfl_exist
      !
      ittyp = ITASK_BEGSTE  ! AB: Manually set this, because Filter(itask) is not called before Begste
#ifdef DBPARTICLES
      call pts_output_parall_db()
      call pts_output()
#else
      call pts_output()
#endif          
      if (INOTMASTER) where( lagrtyp(:) % kfl_exist == -5) lagrtyp(:) % kfl_exist = -1
      kfl_injec = 0

      call ker_timeline('END_ASSEMBLY', nlagr_new)
    end subroutine pts_injection_injectors

end module mod_pts_injection
!> @}
