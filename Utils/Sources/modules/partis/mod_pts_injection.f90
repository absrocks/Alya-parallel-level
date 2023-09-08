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
    use mod_memory,         only : memory_alloca
    use mod_memory,         only : memory_deallo
    use mod_memory,         only : memory_copy
    use mod_communications, only : PAR_BROADCAST
    use mod_parall,         only : PAR_COMM_MY_CODE_WM4
    use mod_messages,       only : livinf
    implicit none

    private

    public :: pts_injection_initialization
    public :: pts_injection_injectors
    public :: pts_injection_parallelization
    public :: pts_injection_reallocate

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
    !> @brief   Reallocate memory
    !> @details Reallocate memory
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_reallocate(particle_place, particle_injector)
        integer(ip), intent(in), pointer :: particle_place(:) !< Particle place
        integer(ip), intent(in), pointer :: particle_injector(:) !< Particle place
        integer(ip) :: iinj
        if (INOTMASTER) then
            do iinj = 1, size(kfl_injla_pts)
                if (kfl_injla_pts(iinj) == 10) then
                    call pts_injection_memory(iinj, particle_place, particle_injector)
                end if
            end do
        end if

    end subroutine pts_injection_reallocate

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
      use mod_random,         only : random_grnd
      implicit none
      integer(ip)          :: ii, jj, kk, iinj
      integer(ip)          :: i1, i2, i3, i4, i5, i6, j1, j2, j3
      real(rp)             :: zfact, tfact, rfact
      real(rp)             :: xx(3), xmini, xmaxi, ymini, ymaxi, zmini, zmaxi, nn(3), basis(ndime, ndime)
      real(rp)             :: x1, y1, z1, x2, y2, z2, x3, y3, z3, p1(3), p2(3), p3(3)
      real(rp)             :: ux, uy, uz, wx, wy, wz, d1, d2, d3, U1, U2, U3, U4
      real(rp)             :: t_aux, t, z, r, radius, rad, xe, ye, ze, rr, hh, alpha, rr_aux
      integer(ip)          :: nside, nlagr_new, nlagr_pos, itype
      integer(ip)          :: number_injected_pts, num_types_injected, num_particles_per_injection
      real(rp),    pointer :: particle_position(:,:)
      integer(ip), pointer :: particle_injector(:)
      integer(ip), pointer :: particle_place(:)
      integer(ip), pointer :: particle_type(:)
      !
      ! Initialization of the random number generator
      !

      kfl_injec = 0
      nlagr_pos = 0
      number_injected_pts = 0
      if (cutim >= tinla_pts) then
         cutla_pts = cutla_pts + dtime
         if (cutla_pts >= tpela_pts - zeror) then
            kfl_injec = 1
            cutla_pts = 0.0_rp
         end if
      end if
      if (cutim >= tfila_pts) then
         kfl_injec = 0
      end if

      if (kfl_injec == 0) return


      call ker_timeline('INI_ASSEMBLY')
      !
      ! Nullify pointers
      !
      nullify(particle_position)
      nullify(particle_injector)
      nullify(particle_place)
      nullify(particle_type)
      !
      ! Number of injected particles when injection
      !
      do iinj = 1, size(kfl_injla_pts)
         num_types_injected = 0
         if (kfl_injla_pts(iinj) /= 0) then
            if (kfl_injty_pts(iinj) == 0) then
               num_types_injected = number_types_pts
            else
               num_types_injected = 1
            end if
         end if
         if (kfl_injla_pts(iinj) == 1) then
            !
            ! Square
            !
            nside = int(parla_pts(iinj, 7), ip)
            number_injected_pts = number_injected_pts + nside * nside * num_types_injected

         else if (kfl_injla_pts(iinj) == 2) then
            !
            ! Sphere injector
            !
            nside = int(parla_pts(iinj, 5), ip)
            number_injected_pts = number_injected_pts + nside ** ndime * num_types_injected

         else if (kfl_injla_pts(iinj) == 3) then
            !
            ! Semi-Sphere injector
            !
            nside = int(parla_pts(iinj, 8), ip)
            !! print*,'holaaaaaaa-semi-esfere-num_inject',  nside ,ndime ,num_types_injected,2 * nside ** ndime * num_types_injected
            number_injected_pts = number_injected_pts + 2 * nside ** ndime * num_types_injected
         else if (kfl_injla_pts(iinj) == 4) then
            !
            ! Circle injector
            !
            nside = int(parla_pts(iinj, 8), ip)

            if (injector_npts_asis(iinj) == 1) then !inject the requested number of particles
               number_injected_pts = number_injected_pts + nside * num_types_injected
            else
               number_injected_pts = number_injected_pts + nside * nside * num_types_injected
            end if

         else if (kfl_injla_pts(iinj) == 5) then
            !
            ! Rectangle injector
            !
            if (ndime == 3) then
               nside = int(parla_pts(iinj, 10), ip)
               number_injected_pts = number_injected_pts + nside * nside * num_types_injected
            end if

         else if (kfl_injla_pts(iinj) == 6) then
            !
            ! Pointwise injector
            !
            nside = max(1_ip, int(parla_pts(iinj, 4), ip))
            number_injected_pts = number_injected_pts + nside * num_types_injected

         else if (kfl_injla_pts(iinj) == 7) then
            !
            ! Cone injector
            !
            nside = int(parla_pts(iinj, 9), ip)
            number_injected_pts = number_injected_pts + nside ** 3 * num_types_injected

         else if (kfl_injla_pts(iinj) == 8) then
            !
            ! Random injector in a cube/square box
            !
            nside = int(parla_pts(iinj, 7), ip)
            number_injected_pts = number_injected_pts + nside * num_types_injected

         else if (kfl_injla_pts(iinj) == 9) then
            !
            ! Segment injector
            !
            nside = int(parla_pts(iinj, 7), ip)
            number_injected_pts = number_injected_pts +  nside * num_types_injected

         else if (kfl_injla_pts(iinj) == 10) then
            !
            ! Read from file
            !
            number_injected_pts = number_injected_pts + injection_pts(iinj) % number_particles * num_types_injected
         end if

      end do

      if (INOTMASTER) then
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_POSITION', 'pts_inject', particle_position, ndime, number_injected_pts)
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_INJECTOR', 'pts_inject', particle_injector, number_injected_pts)
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_INJECTOR', 'pts_inject', particle_place,    number_injected_pts)
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_INJECTOR', 'pts_inject', particle_type,     number_injected_pts)

         do iinj = 1, size(kfl_injla_pts)

            if (kfl_injty_pts(iinj) == 0) then
               num_types_injected = number_types_pts
            else
               num_types_injected = 1
            end if

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

               do ii = 1, nside
                  U1 = random_grnd(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                  U2 = random_grnd(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                  U3 = random_grnd(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
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

            else if (kfl_injla_pts(iinj) == 4) then
               !----------------------------------------------------------
               !
               ! Circle injector
               !
               !----------------------------------------------------------
               nside = int(parla_pts(iinj, 8), ip)

               if (ndime == 3) then !3d mesh
                  xc    = parla_pts(iinj, 1)
                  yc    = parla_pts(iinj, 2)
                  zc    = parla_pts(iinj, 3)
                  rad   = parla_pts(iinj, 4)
                  nn(1) = parla_pts(iinj, 5)
                  nn(2) = parla_pts(iinj, 6)
                  nn(3) = parla_pts(iinj, 7)

                  basis(1:ndime, 1) = nn(1:3)
                  tfact = 2.0_rp * pi / real(nside, rp)
                  rfact = rad / real(nside, rp)

                  if (kfl_random_pts(iinj) == 1) then

                     if (injector_particle_distribution(iinj) == 0) then !uniform points in polar coordinates
                        call livinf(-7_ip, "Injecting particles unformly distributed in polar coordinates for injector ", iinj)
                        do jj = 1, nside
                           do kk = 1, nside
                              U1 = random_grnd(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                              U2 = random_grnd(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                              t = 2.0_rp * pi * U2
                              radius = rad * U1

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

                        num_particles_per_injection = 0
                        if( injector_npts_asis(iinj) == 1 ) then !inject the requested number of particles
                           num_particles_per_injection = nside
                        else
                           num_particles_per_injection = nside * nside
                        end if

                        !open (unit = 1977, file = "particles.csv")

                        do jj = 1, num_particles_per_injection

                           !U1     = random_grnd()
                           !U2     = random_grnd()
                           U1     = random_grnd(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4) 
                           U2     = random_grnd(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4) 
                           t      = 2.0_rp * pi * U2
                           radius = rad * rad * U1 !https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly

                           xx(1) = 0.0_rp
                           xx(2) = sqrt(radius) * cos(t)
                           xx(3) = sqrt(radius) * sin(t)

                           call maths_local_orthonormal_basis(ndime,basis)
                           !!!call maths_vector_to_new_basis(ndime,basis,xx)
                           call maths_vector_from_new_basis(ndime,basis,xx)

                           xx(1) = xx(1) + xc
                           xx(2) = xx(2) + yc
                           xx(3) = xx(3) + zc

                           !write (1977, *) xx(1), ',', xx(2), ',', xx(3)

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

                        !close(1977)

                        !End of generation of uniform cartesian particles

                     end if

                  else  !not kfl_random_pts(iinj) == 1
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

                  xc = parla_pts(iinj, 1)
                  yc = parla_pts(iinj, 2)
                  zc = parla_pts(iinj, 3)
                  rad = parla_pts(iinj, 4)
                  nx = parla_pts(iinj, 5)
                  ny = parla_pts(iinj, 6)
                  nz = parla_pts(iinj, 7)
                  zfact = 2.0_rp / real(nside + 1, rp)
                  tfact = 2.0_rp * pi / real(nside + 1, rp)
                  rfact = rad / real(nside, rp)

                  if (kfl_random_pts(iinj) == 1) then
                     do jj = 1, nside
                        t_aux = real(jj, rp) * tfact
                        do kk = 1, nside
                           U1 = random_grnd(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
                           U2 = random_grnd(UNIQUE_SEED=.true.,PAR_COMM_IN4=PAR_COMM_MY_CODE_WM4)
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

               nside = max(1_ip, int(parla_pts(iinj, 4), ip))
               xx(1:3) = parla_pts(iinj, 1:3)

               do ii = 1, nside
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
               if (kfl_injty_pts(iinj) > 0) then
                  do ii = 1, injection_pts(iinj) % number_particles
                     nlagr_pos = nlagr_pos + 1
                     particle_position(1:ndime, nlagr_pos) = injection_pts(iinj) % coord_particles(ii * ndime - (ndime-1):ii * ndime)
                     particle_injector(ii) = iinj
                     particle_type(nlagr_pos) = kfl_injty_pts(iinj)
                  end do
               else
                  do itype = 1, ntyla_pts
                     if (parttyp(itype) % kfl_exist /= 0) then
                        do ii = 1, injection_pts(iinj) % number_particles
                           nlagr_pos = nlagr_pos + 1
                           particle_position(1:ndime, nlagr_pos )= injection_pts(iinj) % coord_particles(ii * ndime - 2:ii * ndime)
                           particle_injector(ii) = iinj
                           particle_type(nlagr_pos) = itype
                        end do
                     end if
                  end do
               end if
            end if
         end do
         !
         ! Loop over particles to find host elements
         ! NLAGR_POS = Total number of injected particles
         ! NUMBER_INJECTED_PTS = Total number of particles (NLAGR_POS)
         !      
         call pts_checkp(nlagr_pos, nlagr_new, particle_position, particle_injector, particle_place, particle_type)
         !!print*,'salgo de ptj_checkp:nlagr_pos',nlagr_pos,nlagr_new
      end if
      !print*,'NUMBER_INJECTED_PTS /= NLAGR_POS',NUMBER_INJECTED_PTS,NLAGR_POS
      !! if( NUMBER_INJECTED_PTS /= NLAGR_POS .and. INOTMASTER ) call runend('HOUSTON, WE HAVE A PROBLEM')
      !
      ! Check injection. In Parallel, decide who owns the particle
      !
      if( IPARALL ) then
         call PAR_POINT_TO_POINT_ARRAY_OPERATION(particle_place, 'IN MY ZONE', 'MIN RANK OR NEGATIVE')
         if (INOTMASTER) then
            nlagr_new = 0
            do ii = 1,number_injected_pts

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
         call PAR_MAX(number_injected_pts, 'IN MY CODE')
         call PAR_SUM(nlagr_new, 'IN MY CODE')
      end if

      coutp(1) = 'INJECT'
      ioutp(1) = nlagr_new
      coutp(2) = 'LAGRANGIAN PARTICLES OUT OF'
      ioutp(2) = number_injected_pts
      call livinf(-10_ip, ' ', nlagr_pos)

      if( nlagr_new < number_injected_pts ) then
         coutp(1) = 'PARTICLES INJECTED OUT OF THE DOMAIN:'
         ioutp(1) = number_injected_pts - nlagr_new
         call livinf(-15_ip, ' ', 0_ip)
      end if
      !
      ! Actual total number of particles (over all partitions)
      !
      nlacc_pts = nlacc_pts + nlagr_new
      !
      ! Deallocate memory
      !
      call memory_deallo(mem_modul(1:2, modul), 'PARTICLE_POSITION', 'pts_inject', particle_position)
      call memory_deallo(mem_modul(1:2, modul), 'PARTICLE_INJECTOR', 'pts_inject', particle_injector)
      call memory_deallo(mem_modul(1:2, modul), 'PARTICLE_INJECTOR', 'pts_inject', particle_place)
      call memory_deallo(mem_modul(1:2, modul), 'PARTICLE_TYPE', 'pts_inject', particle_type)
      !
      ! Output injected particles and the recover normal flag kfl_exist
      !
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
