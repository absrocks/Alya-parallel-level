!-------------------------------------------------------------------------------
!> @addtogroup variablepressure boundary conditions
!> @{
!> @authors Adria Quintanas-Corominas : adria.quintanas@bsc.es
!> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
!> @date    May, 2020
!> @}
!------------------------------------------------------------------------------
module mod_sld_cardiac_cycle
    ! ----------------------------------------
    use def_kintyp_basic,               only :  ip, rp, lg
    ! ----------------------------------------
    implicit none
    ! ----------------------------------------
    integer(ip),   parameter                :: max_cavities   = 4_ip
    integer(ip),   parameter                :: max_cycles = 12_ip
    integer(ip)                             :: ncavi
    integer(ip)                             :: ncycles
    ! ----------------------------------------
    type heart_cavity
        logical(lg)                         :: initialized = .False.
        integer(ip)                         :: bouset
        integer(ip), dimension(2)           :: plane_nodset
        real(rp),    dimension(3)           :: plane_origin
        real(rp),    dimension(3)           :: plane_vector
        real(rp),    dimension(3,3)         :: plane_points
        real(rp),    dimension(5)           :: volume
        real(rp),    dimension(3)           :: dvol
        real(rp),    dimension(3)           :: ddvol
        real(rp)                            :: ini_vol
        real(rp)                            :: end_dia_vol
        real(rp)                            :: end_sys_vol
        real(rp),    dimension(3)           :: prest = 0.0_rp
        real(rp),    dimension(3)           :: pres = 0.0_rp
        integer(ip)                         :: phase = 0_ip
    end type
    ! ----------------------------------------
    type windkessel
        real(rp)                            :: r
        real(rp)                            :: c
        real(rp)                            :: p_ini
        real(rp), dimension(3)              :: pres=0.0_rp
    end type
    ! ----------------------------------------
    type cycle_control
        logical(lg)                         :: defined = .False.
        logical(lg)                         :: active = .False.
        integer(ip)                         :: cycle_id
        integer(ip)                         :: cavity_controled
        integer(ip)                         :: n_beats = 0_ip
        integer(ip)                         :: max_beats = 0_ip
        integer(ip)                         :: prevcycle = 0_ip
        integer(ip)                         :: nextcycle = 0_ip
        real(rp)                            :: last_phase_change=0.0_rp
        real(rp)                            :: tzero
        real(rp)                            :: tpstr
        real(rp)                            :: pzero
        real(rp)                            :: pstr0
        real(rp), dimension(2)              :: gains_contraction=0.0_rp
        type(windkessel)                    :: assoc_wdk
        real(rp), dimension(2)              :: gains_relaxation=0.0_rp
        real(rp)                            :: ppost
    end type
    ! ----------------------------------------
    type(heart_cavity), dimension(max_cavities)    :: cavities
    type(cycle_control), dimension(max_cycles) :: cycles
    ! ----------------------------------------
    
    contains 

        ! ---------------------------------------------------------------------
        !> @author   Adrià Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date    
        !> @brief    Initialise variables
        !> @details  
        !> @todo     
        ! ---------------------------------------------------------------------
        subroutine sld_cardiac_cycle_initialise_variables ( &
            )
            use def_solidz,           only :  kfl_cycle_sld, kfl_volca_sld
            use def_kintyp_basic,     only :  ip, rp, lg
            ! -------------------------------
            implicit none
            ! -------------------------------
            kfl_cycle_sld = 0_ip     ! Flag for cardiac cycle management (by default is OFF)
            kfl_volca_sld = 0_ip     ! Flag for compute cavity volumes
            ncavi         = 0_ip     ! number of cavities (ventricles)
            ncycles       = 0_ip     ! number of cavities (ventricles)
            ! -------------------------------
            
        end subroutine sld_cardiac_cycle_initialise_variables

        ! ---------------------------------------------------------------------
        !> @author   Adrià Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date    
        !> @brief    Write information
        !> @details  
        !> @todo     
        ! ---------------------------------------------------------------------
        subroutine sld_cardiac_cycle_write_res( &
            itask )
            use def_solidz, only            :  lun_carcy_res_sld
            use def_master, only            :  ittim, cutim
            use def_master, only            :  TIME_N
            use def_kintyp_basic, only      :  ip, rp, lg
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip), intent(in)        :: itask
            real(rp)                       :: auxr1, auxr2
            integer(ip)                    :: icycle, cycid1, cycid2
            ! -------------------------------
            select case(itask)

            case( 1_ip )

                write(lun_carcy_res_sld,*) '# --| ALYA Cardiac Cycle Computations'
                write(lun_carcy_res_sld,*) '# --| '
                write(lun_carcy_res_sld,*) '# --|  1. Step increment'
                write(lun_carcy_res_sld,*) '# --|  2. Step time'
                write(lun_carcy_res_sld,*) '# --|  3. Cavity #1'
                write(lun_carcy_res_sld,*) '# --|  4. Cavity #1: cycle'
                write(lun_carcy_res_sld,*) '# --|  5. Cavity #1: phase'
                write(lun_carcy_res_sld,*) '# --|  6. Cavity #1: volume'
                write(lun_carcy_res_sld,*) '# --|  7. Cavity #1: cavity pressure'
                write(lun_carcy_res_sld,*) '# --|  8. Cavity #1: outflow windkessel pressure'
                write(lun_carcy_res_sld,*) '# --|  9. Cavity #2'
                write(lun_carcy_res_sld,*) '# --|  10. Cavity #2: cycle'
                write(lun_carcy_res_sld,*) '# --|  11. Cavity #2: phase'
                write(lun_carcy_res_sld,*) '# --|  12. Cavity #2: volume'
                write(lun_carcy_res_sld,*) '# --|  13. Cavity #2: cavity pressure'
                write(lun_carcy_res_sld,*) '# --|  14. Cavity #2: outflow windkessel pressure'
                write(lun_carcy_res_sld,*) '# --| '
                
            case( 2_ip )
                auxr1=0.0_rp
                auxr2=0.0_rp
                cycid1=0_ip
                cycid2=0_ip
                do icycle=1, max_cycles
                    if(  cycles(icycle) % active .and. cycles(icycle) % cavity_controled .eq. 1) then
                     auxr1 = cycles(icycle) % assoc_wdk % pres(TIME_N)
                     cycid1= cycles(icycle) % cycle_id
                    elseif(  cycles(icycle) % active .and. cycles(icycle) % cavity_controled .eq. 2) then
                     auxr2 = cycles(icycle) % assoc_wdk % pres(TIME_N)
                     cycid2= cycles(icycle) % cycle_id
                    endif
                enddo

                write(lun_carcy_res_sld, 101)                              &
                    &   ittim,                                 &  ! Step increment
                    &   cutim,                                 &  ! Step time
                    &   1_ip,                                  &  ! Cavity #1
                    &   cycid1,                                &  ! Cavity #1: cycle
                    &   cavities(1) % phase,                   &  ! Cavity #1: Phase.
                    &   cavities(1) % volume(TIME_N),          &  ! Cavity #1: volume
                    &   cavities(1) % pres(TIME_N),            &  ! Cavity #1: cavity pressure
                    &   auxr1,                                 &  ! Cavity #1: outflow pressure
                    &   2_ip,                                  &  ! Cavity #2
                    &   cycid2,                                &  ! Cavity #2: cycle
                    &   cavities(2) % phase,                   &  ! Cavity #2: Phase.
                    &   cavities(2) % volume(TIME_N),          &  ! Cavity #2: volume
                    &   cavities(2) % pres(TIME_N),            &  ! Cavity #1: cavity pressure
                    &   auxr2                                     ! Cavity #1: outflow pressure

                flush( lun_carcy_res_sld )
                
            end select
            ! -------------------------------
            101 format(2x, i9, 2x, e16.8e3, 4(3(2x,i2) ,3(2x,e16.8e3)))
            ! -------------------------------
        end subroutine sld_cardiac_cycle_write_res

        ! ---------------------------------------------------------------------
        !> @author   Adrià Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date    
        !> @brief    Write information
        !> @details  
        !> @todo     
        ! ---------------------------------------------------------------------
        subroutine sld_cardiac_cycle_write_cvg( &
            itask )
            use def_solidz, only            :  lun_carcy_cvg_sld
            use def_master, only            :  itinn, ittim, modul
            use def_master, only            :  ITER_K
            use def_kintyp_basic, only      :  ip, rp, lg
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip), intent(in)        :: itask
            real(rp)                       :: auxr1, auxr2
            integer(ip)                    :: icycle, cycid1, cycid2
            ! -------------------------------
            select case(itask)

            case( 1_ip )

                write(lun_carcy_cvg_sld,*) '# --| ALYA Cardiac Cycle Computations'
                write(lun_carcy_cvg_sld,*) '# --| '
                write(lun_carcy_cvg_sld,*) '# --|  1. Step increment'
                write(lun_carcy_cvg_sld,*) '# --|  2. Iteration number'
                write(lun_carcy_cvg_sld,*) '# --|  3. Cavity #1'
                write(lun_carcy_cvg_sld,*) '# --|  4. Cavity #1: cycle'
                write(lun_carcy_cvg_sld,*) '# --|  5. Cavity #1: phase'
                write(lun_carcy_cvg_sld,*) '# --|  6. Cavity #1: volume'
                write(lun_carcy_cvg_sld,*) '# --|  7. Cavity #1: cavity pressure'
                write(lun_carcy_cvg_sld,*) '# --|  8. Cavity #1: outflow windkessel pressure'
                write(lun_carcy_cvg_sld,*) '# --|  9. Cavity #2'
                write(lun_carcy_cvg_sld,*) '# --|  10. Cavity #2: cycle'
                write(lun_carcy_cvg_sld,*) '# --|  11. Cavity #2: phase'
                write(lun_carcy_cvg_sld,*) '# --|  12. Cavity #2: volume'
                write(lun_carcy_cvg_sld,*) '# --|  13. Cavity #2: cavity pressure'
                write(lun_carcy_cvg_sld,*) '# --|  14. Cavity #2: outflow windkessel pressure'
                write(lun_carcy_cvg_sld,*) '# --| '
                
            case( 2_ip )
                auxr1=0.0_rp
                auxr2=0.0_rp
                cycid1=0_ip
                cycid2=0_ip
                do icycle=1, max_cycles
                    if(  cycles(icycle) % active .and. cycles(icycle) % cavity_controled .eq.1) then
                     auxr1 = cycles(icycle) % assoc_wdk % pres(ITER_K)
                     cycid1= cycles(icycle) % cycle_id
                    elseif(  cycles(icycle) % active .and. cycles(icycle) % cavity_controled .eq.2) then
                     auxr2 = cycles(icycle) % assoc_wdk % pres(ITER_K)
                     cycid2= cycles(icycle) % cycle_id
                    endif
                enddo
                write(lun_carcy_cvg_sld, 101)                              &
                    &   ittim,                                 &  ! Step increment
                    &   itinn(modul),                          &  ! Iteration number
                    &   1_ip,                                  &  ! Cavity #1
                    &   cycid1,                                &  ! Cavity #1: cycle
                    &   cavities(1) % phase,                 &  ! Cavity #1: Phase.
                    &   cavities(1) % volume(ITER_K),        &  ! Cavity #1: volume
                    &   cavities(1) % pres(ITER_K),          &  ! Cavity #1: cavity pressure
                    &   auxr1,                                 &  ! Cavity #1: outflow pressure
                    &   2_ip,                                  &  ! Cavity #2
                    &   cycid2,                                &  ! Cavity #1: cycle
                    &   cavities(2) % phase,                 &  ! Cavity #2: Phase.
                    &   cavities(2) % volume(ITER_K),        &  ! Cavity #2: volume
                    &   cavities(2) % pres(ITER_K),          &  ! Cavity #1: cavity pressure
                    &   auxr2                                     ! Cavity #1: outflow pressure

                flush( lun_carcy_cvg_sld )
                
            end select
            ! -------------------------------
            101 format(2x, i9, 2x, i9, 2(3(2x,i2) ,3(2x,e16.8e3)))
            ! -------------------------------
        end subroutine sld_cardiac_cycle_write_cvg

        ! ---------------------------------------------------------------------
        !> @author   Adrià Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date    
        !> @brief    Read data
        !> @details  
        !> @todo     
        ! ---------------------------------------------------------------------
        subroutine sld_cardiac_cycle_read_data( &
            )
            use def_kintyp_basic,     only :  ip, rp, lg
            use def_domain,           only :  ndime
            use def_inpout,           only :  words, param, exists, getint, getrea, getcha
            use mod_ecoute,           only :  ecoute
            use mod_messages,         only :  livinf
            use def_solidz,           only :  kfl_cycle_sld, kfl_volca_sld
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip)                    :: icavi
            integer(ip)                    :: id_aux
            real(rp)                       :: p1(3), p2(3), v1(3), nrm
            ! -------------------------------

            call ecoute('sld_cardiac_cycle_read_data')

            !CARDIAC_CYCLE section
            do while(words(1)/='ENDCA')

                if( words(1) == 'CAVIT' )then
                    ncavi=ncavi+1
                    kfl_volca_sld = 1_ip

                    !
                    ! Get cavity identification and check it was initialized
                    !
                    icavi = int(param(1))
                    if( cavities(icavi) % initialized )then
                        call runend('MOD_SLD_CARDIAC_CYCLE: Two cavities with the same ID')
                    else
                        cavities(icavi) % initialized  = .True.
                        cavities(icavi) % bouset       = 0_ip
                        cavities(icavi) % plane_nodset = 0_ip
                        cavities(icavi) % plane_origin = 0.0_rp
                        cavities(icavi) % plane_vector = 0.0_rp
                        cavities(icavi) % volume       = 0.0_rp
                    end if

                    call ecoute('sld_cardiac_cycle_read_data')
                    do while(words(1)/='ENDCA')

                        !
                        ! Boundary defining the cavity
                        !
                        if (words(1)=='BOUND') then 
                            cavities(icavi) % bouset = int(param(1))

                        endif
                        
                        !
                        ! Plane definition for open cavity
                        !
                        if(  words(1) == 'PLANE' )then

                            if( exists('POINT') )then

                                ! Plane defined by 2 points
                                p1(1:ndime) = param(2:4)
                                p2(1:ndime) = param(5:7)
                                v1(:)= p1(:) - p2(:)
                                nrm = sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
                                if( nrm < 0.0_rp ) call runend('MOD_SLD_CARDIAC_CYCLE: ||n|| of the plane is <= 0')
                                cavities(icavi) % plane_vector(:) = v1(:)/nrm
                                
                            elseif( exists('NODSE') .or. exists('SETNO') )then

                                ! Plane defined by a node set
                                cavities(icavi) % plane_nodset = (/int(param(2)), int(param(3))/)

                            elseif( exists('ORIGI') )then

                                ! Plane origin
                                cavities(icavi) % plane_origin = param(2:4)
                                
                            else
                                ! Plane defined by the normal vector
                                call runend('MOD_SLD_CARDIAC_CYCLE: PLANE POINT / OIRIGIN / NODSET')
                                
                            endif
                        endif

                        call ecoute('sld_cardiac_cycle_read_data')
                    end do

                else if( words(1) == 'VARIA') then
                    ncycles=ncycles+1
                    kfl_cycle_sld = 1_ip
                    do while( words(1) /= 'ENDVA' )
                      call ecoute('sld_cardiac_cycle_read_data')
                      if( words(1) == 'CYCLE') then
                        id_aux=getint('CYCLE',0_ip,'#CYCLE_ID')
                        cycles(id_aux) % cycle_id = id_aux
                        cycles(id_aux) % defined  = .True.
                        do while( words(1) /= 'ENDCY' )
                            if ( words(1) == 'FORCA') cycles(id_aux) %  cavity_controled = param(1)
                            if ( words(1) == 'MAXBE') cycles(id_aux) %  max_beats = param(1)
                            if ( words(1) == 'PREVC') cycles(id_aux) %  prevcycle = param(1)
                            if ( words(1) == 'NEXTC') cycles(id_aux) %  nextcycle = param(1)
                            if ( words(1) == 'TZERO') cycles(id_aux) %  tzero = param(1)
                            if ( words(1) == 'TPSTR') cycles(id_aux) %  tpstr = param(1)
                            if ( words(1) == 'GAINC') then
                                cycles(id_aux) %  gains_contraction(1) = param(1)
                                cycles(id_aux) %  gains_contraction(2) = param(2)
                            endif
                            if ( words(1) == 'PZERO') cycles(id_aux) %  pzero = param(1)
                            if ( words(1) == 'PSTR0') cycles(id_aux) %  pstr0 = param(1)
                            if ( words(1) == 'PART0') cycles(id_aux) %  assoc_wdk % p_ini = param(1)
                            if ( words(1) == 'CPRES') cycles(id_aux) %  assoc_wdk % c     = param(1)
                            if ( words(1) == 'RPRES') cycles(id_aux) %  assoc_wdk % r     = param(1)
                            if ( words(1) == 'GAINR') then
                                cycles(id_aux) %  gains_relaxation(1) = param(1)
                                cycles(id_aux) %  gains_relaxation(2) = param(2)
                            endif
                            if ( words(1) == 'PPOST') cycles(id_aux) %  ppost = param(1)

                           call ecoute('sld_cardiac_cycle_read_data')
                        end do
                      endif
                    enddo

                end if

                call ecoute('sld_cardiac_cycle_read_data')
            end do
            ! -------------------------------
        end subroutine sld_cardiac_cycle_read_data

        ! ---------------------------------------------------------------------
        !> @author   Adrià Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date    
        !> @brief    Interchange input data
        !> @details  
        !> @todo     
        ! ---------------------------------------------------------------------
        subroutine sld_cardiac_cycle_send_data( &
            itask )
            use def_kintyp_basic,     only :  ip, rp, lg
            use def_solidz,           only :  kfl_cycle_sld, kfl_volca_sld
            ! -------------------------------
            implicit none
            integer(ip),  intent(in)       :: itask
            integer(ip)                    :: icavi, icycle, ii,jj
            ! ------------------------------- 
            select case(itask)

            case( 1_ip )
                !
                ! Activation flags
                !
                call iexcha(kfl_volca_sld)
                call iexcha(kfl_cycle_sld)
                call iexcha(ncavi)
                call iexcha(ncycles)
            
            case( 2_ip )
                !
                ! Read information
                !
                if( kfl_volca_sld > 0_ip )then
                
                do icavi = 1_ip, max_cavities
                        call lexcha(cavities(icavi) % initialized )
                        call iexcha(cavities(icavi) % phase )
                        call iexcha(cavities(icavi) % bouset )
                        do ii = 1, 2
                            call iexcha(cavities(icavi) % plane_nodset(ii) )
                        enddo
                        do ii = 1, 3
                            call rexcha(cavities(icavi) % plane_vector(ii) )
                            call rexcha(cavities(icavi) % plane_origin(ii) )
                            do jj=1,3
                                call rexcha(cavities(icavi) % plane_points(ii,jj) )
                            enddo
                        end do
                    end do
                endif
                ! ------------------------------- 
                if( kfl_cycle_sld > 0_ip )then
                    do icycle = 1, max_cycles
                            call lexcha(cycles(icycle) % defined )
                            call lexcha(cycles(icycle) % active )
                            call iexcha(cycles(icycle) % cycle_id )
                            call iexcha(cycles(icycle) % cavity_controled )
                            call rexcha(cycles(icycle) % cavity_controled )
                            call iexcha(cycles(icycle) % n_beats )
                            call iexcha(cycles(icycle) % max_beats )
                            call iexcha(cycles(icycle) % prevcycle )
                            call iexcha(cycles(icycle) % nextcycle )
                            call rexcha(cycles(icycle) % tzero)
                            call rexcha(cycles(icycle) % tpstr)
                            call rexcha(cycles(icycle) % pzero)
                            call rexcha(cycles(icycle) % pstr0)
                            call rexcha(cycles(icycle) % gains_contraction(1))
                            call rexcha(cycles(icycle) % gains_contraction(2))
                            call rexcha(cycles(icycle) % assoc_wdk % p_ini)
                            call rexcha(cycles(icycle) % assoc_wdk % c)
                            call rexcha(cycles(icycle) % assoc_wdk % r)
                            call rexcha(cycles(icycle) % gains_relaxation(1))
                            call rexcha(cycles(icycle) % gains_relaxation(2))
                            call rexcha(cycles(icycle) % ppost)
                    end do
                endif

            end select
            ! -------------------------------
            contains 
                subroutine lexcha( lg_vari )
                    use def_kintyp_basic,     only :  ip, lg
                    logical(lg), intent(inout)    ::  lg_vari
                    integer(ip)                   ::  ip_vari
                    if( lg_vari )then
                        ip_vari = 1_ip
                    else
                        ip_vari = 0_ip
                    endif
                    call iexcha( ip_vari )
                    if( ip_vari == 1_ip )then
                        lg_vari = .True.
                    else
                        lg_vari = .False.
                    endif
                end subroutine lexcha
            ! -------------------------------
        end subroutine sld_cardiac_cycle_send_data


        ! ---------------------------------------------------------------------
        !> @author   Adrià Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date    
        !> @brief    Write information
        !> @details  
        !> @todo     
        ! ---------------------------------------------------------------------
        subroutine sld_cardiac_cycle_manage_restart( &
            itask )
            use def_master, only            :  INOTSLAVE, IPARALL, momod, modul
            use def_master, only            :  ITASK_READ_RESTART, ITASK_WRITE_RESTART
            use def_kintyp_basic, only      :  ip, rp, lg
            use mod_communications, only    : PAR_BROADCAST
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip), intent(in)        :: itask
            ! -------------------------------
            integer(ip)                    :: icavi, icycle
            ! -------------------------------

            if( itask == ITASK_READ_RESTART ) then
                !
                ! Read restart
                !
                if( INOTSLAVE )then
                    
                    do icavi = 1, max_cavities
                        read(momod(modul) % lun_rstar) cavities(icavi) % volume(:)
                        read(momod(modul) % lun_rstar) cavities(icavi) % dvol(:)
                        read(momod(modul) % lun_rstar) cavities(icavi) % ddvol(:)
                        read(momod(modul) % lun_rstar) cavities(icavi) % ini_vol
                        read(momod(modul) % lun_rstar) cavities(icavi) % end_dia_vol
                        read(momod(modul) % lun_rstar) cavities(icavi) % end_sys_vol
                        read(momod(modul) % lun_rstar) cavities(icavi) % prest(:)
                        read(momod(modul) % lun_rstar) cavities(icavi) % pres(:)
                        read(momod(modul) % lun_rstar) cavities(icavi) % phase
                    end do 
 
                    do icycle = 1, max_cycles
                        read(momod(modul) % lun_rstar) cycles(icycle) % last_phase_change
                        read(momod(modul) % lun_rstar) cycles(icycle) % assoc_wdk % pres 
                    end do

                end if

                if( IPARALL )then

                    do icavi = 1, max_cavities
                        call ARR_BROADCAST(cavities(icavi) % volume)
                        call ARR_BROADCAST(cavities(icavi) % dvol)
                        call ARR_BROADCAST(cavities(icavi) % ddvol)
                        call PAR_BROADCAST(cavities(icavi) % ini_vol)
                        call PAR_BROADCAST(cavities(icavi) % end_dia_vol)
                        call PAR_BROADCAST(cavities(icavi) % end_sys_vol)
                        call ARR_BROADCAST(cavities(icavi) % prest)
                        call ARR_BROADCAST(cavities(icavi) % pres)
                        call PAR_BROADCAST(cavities(icavi) % phase)
                    end do 

                    do icycle = 1, max_cycles
                        call PAR_BROADCAST(cycles(icycle) % last_phase_change)
                        call ARR_BROADCAST(cycles(icycle) % assoc_wdk % pres) 
                    end do

                end if

           
             else if( itask == ITASK_WRITE_RESTART ) then 
                !
                ! Write restart file
                !
                if( INOTSLAVE )then
                    
                    do icavi = 1, max_cavities
                        write(momod(modul) % lun_rstar) cavities(icavi) % volume(:)
                        write(momod(modul) % lun_rstar) cavities(icavi) % dvol(:)
                        write(momod(modul) % lun_rstar) cavities(icavi) % ddvol(:)
                        write(momod(modul) % lun_rstar) cavities(icavi) % ini_vol
                        write(momod(modul) % lun_rstar) cavities(icavi) % end_dia_vol
                        write(momod(modul) % lun_rstar) cavities(icavi) % end_sys_vol
                        write(momod(modul) % lun_rstar) cavities(icavi) % prest(:)
                        write(momod(modul) % lun_rstar) cavities(icavi) % pres(:)
                        write(momod(modul) % lun_rstar) cavities(icavi) % phase
                    end do 

                    do icycle = 1, max_cycles
                        write(momod(modul) % lun_rstar) cycles(icycle) % last_phase_change
                        write(momod(modul) % lun_rstar) cycles(icycle) % assoc_wdk % pres 
                    end do

                end if

            endif

            contains

                subroutine ARR_BROADCAST( array )
                    implicit none
                    integer(ip) :: ii
                    real(rp)    :: array(:)
                    do ii = 1, size(array,dim=1)
                        call PAR_BROADCAST(array(ii))
                    end do
                end subroutine ARR_BROADCAST

            ! -------------------------------
        end subroutine sld_cardiac_cycle_manage_restart

        ! ---------------------------------------------------------------------
        !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @author   Mariano Vázquez (mariano.vazquez@bsc.es)
        ! ---------------------------------------------------------------------
        subroutine sld_manage_cycles( &
            itask )
            use def_kintyp_basic,       only :  ip, rp, lg
            !--------------------------------
            implicit none
            ! -------------------------------
            integer(ip),  intent(in)       :: itask
            integer(ip)                    :: icycle
            ! -------------------------------
            do icycle=1, max_cycles
                if(cycles(icycle) % defined) then
                    cycles(icycle) % active=.True.
                    call sld_cardiac_phases(itask, cycles(icycle) % cycle_id)
                endif
            enddo
            ! -------------------------------
        endsubroutine sld_manage_cycles

        ! ---------------------------------------------------------------------
        !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @author   Mariano Vázquez (mariano.vazquez@bsc.es)
        ! ---------------------------------------------------------------------
        subroutine sld_cardiac_phases( &
            itask, id )
            use def_kintyp_basic,     only :  ip, rp, lg
            use def_master,           only :  ittim, itinn, cutim, modul, dtime
            use def_master,           only :  ITER_K, TIME_N, TIME_N_MINUS_1, ITASK_ENDINN, ITASK_ENDITE
            !--------------------------------
            implicit none
            ! -------------------------------
            integer(ip),  intent(in)       :: itask
            integer(ip),  intent(in)       :: id
            ! -------------------------------
            integer(ip)                    :: icavi
            real(rp)                       :: t_init, t_pres, preload_press, prestr_press, c_wdk, r_wdk, ini_wdk, p_fill, fill_rate
            real(rp)                       :: dvol_aux
            real(rp)                       :: delta_phase_change
            real(rp)                       :: Cvol, Cvel          
            real(rp)                       :: gain_err_c, gain_derr_c, gain_err_r, gain_derr_r
            ! -------------------------------

            ! ------------------------< VARIABLE DEFINITION >-----------------------------------
            icavi         = cycles(id) %  cavity_controled
            t_init        = cycles(id) %  tzero
            t_pres        = cycles(id) %  tpstr
            preload_press = cycles(id) %  pzero
            prestr_press  = cycles(id) %  pstr0
            gain_err_c    = cycles(id) %  gains_contraction(1)
            gain_derr_c   = cycles(id) %  gains_contraction(2)
            ini_wdk       = cycles(id) %  assoc_wdk % p_ini
            c_wdk         = cycles(id) %  assoc_wdk % c
            r_wdk         = cycles(id) %  assoc_wdk % r
            gain_err_r    = cycles(id) %  gains_relaxation(1)
            gain_derr_r   = cycles(id) %  gains_relaxation(2)
            p_fill        = cycles(id) %  ppost
            


            ! -------< INITIALISATION. FIRST TIME ITERATION AND FIRST INNER ITERATION>----------
            ! ----------------------------------------------------------------------------------
            if ((ittim == 1) .and. (itinn(modul) == 1)) then
                cycles(id) % assoc_wdk % pres(ITER_K) = ini_wdk
                cycles(id) % assoc_wdk % pres(TIME_N) = ini_wdk
                cavities(icavi) % ini_vol  = cavities(icavi) % volume(ITER_K)
            endif



            select case(itask)

                ! -------------< WHAT HAPPENS AT THE END OF EACH INNER ITERATION>-----------------
                ! -------------------------------------------------------------------------------
                case (ITASK_ENDINN)

                    ! Calculate volume increment and rate change
                    cavities(icavi) % dvol(ITER_K) = cavities(icavi) % volume(ITER_K)   - cavities(icavi) % volume(TIME_N)
                    cavities(icavi) % ddvol(ITER_K) = cavities(icavi) % dvol(ITER_K) / dtime

                    ! -------------------< EVOLVE WINDKESSEL (aortic pressure)  MODEL >---------------------
                    if(cavities(icavi) % phase .ne. 0) then
                        dvol_aux=0.0_rp
                        if(cavities(icavi) % phase .eq. 2) dvol_aux = cavities(icavi) % ddvol(ITER_K)
                        cycles(id) % assoc_wdk % pres(ITER_K) = cycles(id) % assoc_wdk % pres(TIME_N) - (dtime / c_wdk) *  &
                        (dvol_aux + (cycles(id) % assoc_wdk % pres(ITER_K) / r_wdk))
                    endif

                   select case(cavities(icavi) % phase)
                   ! -----------------------< PHASES OF THE CYCLE  >------------------------------
                   !
                   ! Phases:
                   !   
                   !    0) Preload and prestress (phase = 0)
                   !        * PARAMETERS:
                   !            + pzero_sld
                   !            + tzero_sld
                   !            + pstr0_sld
                   !            + tpstr_sld

                   !    1) Isovolumetric contraction:  Keeps volume constant
                   !        * PARAMETERS:
                   !            + gains_contraction
                   !    2) Ejection: Ventricular pressure equal to aorta
                   !
                   !    3) Isovolumetric relaxation: Keeps volume constant
                   !        * PARAMETERS:
                   !            + gains_relaxation
                   !
                   !    4) Filling: fills ventricle
                   !        * PARAMETERS:
                   !            + none
                   !
                   ! -------------------------------------------------------------------------------
                   ! --------------------< FLAG=0. PRELOAD AND PRE-STRESS  >------------------------
                   case(0_ip)
                        ! ------------------< ADDING PRELOAD CONTRIBUTION >-------------------------
                        cavities(icavi) % pres(ITER_K) = preload_press * (cutim / t_init)

                
                        ! ------------------< ADDING PRESTRES CONTRIBUTION >-------------------------
                        if(prestr_press .gt. 0.1_rp) then
                          if (cutim < t_pres) then
                              cavities(icavi) % prest(ITER_K) = prestr_press * (cutim / t_pres)
                          else
                              Cvol = cavities(icavi) % volume(ITER_K) / cavities(icavi) % prest(ITER_K)
                              Cvel = 1.0_rp / 1.0_rp ! Acceptable ranges are from 1.0 to 1.5

                              dvol_aux = cavities(icavi) % volume (ITER_K) - cavities(icavi) % ini_vol

                              cavities(icavi) % prest(ITER_K) = cavities(icavi) % prest(TIME_N) + dvol_aux / Cvol + cavities(icavi) % ddvol(ITER_K) / Cvel

                              if (cavities(icavi) % prest(ITER_K) < 1.0_rp) then
                                call runend('sld_windk: negative pre-stress!')
                                cavities(icavi) % prest(ITER_K) = 1.0_rp
                              endif
                          end if
                        endif

                        ! ------------------< SAVE LAST COMPUTED VALUE >-------------------------
                        cavities(icavi) % end_dia_vol = cavities(icavi) % volume(ITER_K)
                   ! ----------------------------------------------------------------------------
            
                   ! -------------< FLAG=1. ISOVOLUMETRIC CONTRACTION  >------------------------
                   case(1_ip)
                        !gain_err = cavities(icavi) % volume(ITER_K) / cavities(icavi) % pres(ITER_K)

                        dvol_aux = cavities(icavi) % volume (ITER_K) - cavities(icavi) % ini_vol
                        gain_err_c = cavities(icavi) % volume (ITER_K) / cavities(icavi) % pres(ITER_K)

                        cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N) - gain_err_c * dvol_aux - gain_derr_c * cavities(icavi) % ddvol(ITER_K)
                        if(cavities(icavi) % pres(ITER_K) .lt. cavities(icavi) % pres(TIME_N))then
                            cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N)
                        endif

                        ! Save new value in end diastolic in case oscilations or inertia make it larger
                        if(cavities(icavi) % volume(TIME_N) .gt. cavities(icavi) % end_dia_vol)then
                           cavities(icavi) % end_dia_vol = cavities(icavi) % volume(TIME_N)
                        endif

                   ! ----------------------------------------------------------------------------

                   ! -------------< FLAG=2. EJECTION  >------------------------------------------
                   case(2_ip)
                       ! -------------------< VENTRICLE PRESSURE EQUAL TO WDK >---------------------
                       cavities(icavi) % pres(ITER_K) = cycles(id) % assoc_wdk % pres(ITER_K)
                       cavities(icavi) % end_sys_vol = cavities(icavi) % volume(ITER_K)     
                   ! ----------------------------------------------------------------------------

                   ! -------------< FLAG=3. ISOVOLUMETRIC RELAXATION  >---------------------------
                   case(3_ip)
                       ! Fixed pressure is overrided and replaced with this automatic value
                       gain_err_r =  cavities(icavi) % pres(ITER_K)/cavities(icavi) % volume(ITER_K)
                       
                       dvol_aux = cavities(icavi) % volume(ITER_K) - cavities(icavi) % end_sys_vol

                       cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N) - gain_err_r*dvol_aux  - gain_derr_r * cavities(icavi) % ddvol(ITER_K)

                        ! clipping works well with automatic gain
                        if(cavities(icavi) % pres(ITER_K) .gt. cavities(icavi) % pres(TIME_N))then
                            cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N)
                        endif

                   ! ----------------------------------------------------------------------------

                   ! -------------< FLAG=4. DIASTOLIC FILLING > ---------------------------------
                   case(4_ip)

                       fill_rate=(preload_press-p_fill)/(cavities(icavi) % ini_vol - cavities(icavi) % end_sys_vol)
                       if(fill_rate.lt.-500.0_rp) fill_rate=-500.0_rp

                       cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N) + fill_rate* cavities(icavi) % dvol(ITER_K)
                        
                      if(cavities(icavi) % pres(ITER_K).lt.preload_press)then
                            cavities(icavi) % pres(ITER_K)=preload_press
                      endif
              
                   end select
            
            
                ! -------------< WHAT HAPPENS AT THE END OF EACH TIME ITERATION>-----------------
                ! -------------------------------------------------------------------------------
                case (ITASK_ENDITE)

                    ! -------------< VARIABLE UPDATE >-------------------------------------------   
                    cavities(icavi) % dvol(TIME_N)  = cavities(icavi) % volume(TIME_N) - cavities(icavi) % volume(TIME_N_MINUS_1)
                    cavities(icavi) % ddvol(TIME_N) = cavities(icavi) % dvol(TIME_N) / dtime
                    cavities(icavi) % pres(TIME_N) = cavities(icavi) % pres(ITER_K)
                    cavities(icavi) % prest(TIME_N) = cavities(icavi) % prest(ITER_K)
                    cycles(id) % assoc_wdk % pres(TIME_N) = cycles(id) % assoc_wdk % pres(ITER_K)
                    delta_phase_change= cutim-cycles(id) % last_phase_change

                    ! ----------------------------------------------------------------------------
                
                    ! -------------< PHASE TRANSITION >-------------------------------------------
                    ! Phase description:
                    !   -> 0:  initialisation
                    !   -> 1:  isovol contraction
                    !   -> 2:  ejection
                    !   -> 3:  isovol relax
                    !   -> 4:  filling
                    !   
                    ! -------------< preload(0)->isolov contraction(1) >--------------------------
                    if ((cavities(icavi) % phase == 0) .and. &
                         cutim .ge. t_init) then
                       cavities(icavi) % phase = 1
                       cycles(id) % last_phase_change=cutim
                    ! --------------------------------------------------------------------------
                    ! -------------< isovol contract(1)->ejection(2) >--------------------------
                    ! CONDITIONS:
                    !   * in phase 1
                    !   * wdk pressure greater than windkessel pressure 
                    !
                    elseif ((cavities(icavi) % phase == 1) .and. &
                        (cavities(icavi) % pres(TIME_N) > cycles(id) % assoc_wdk % pres(TIME_N))) then

                       cavities(icavi) % phase = 2
                       cycles(id) % last_phase_change=cutim
                    ! --------------------------------------------------------------------------
                    ! -------------< ejection(2)->isovol relax(2) >-----------------------------
                    ! CONDITIONS:
                    !   * in phase 2
                    !   * Volume is increasing (geometry is relaxing)
                    !   * It's still contracted (required for oscilations)
                    !
                    elseif ((cavities(icavi) % phase == 2)                 .and. &
                            (cavities(icavi) % dvol(TIME_N) > 0.00001_rp) .and. &
                            delta_phase_change.gt.0.05_rp .and.    &
                            (cavities(icavi) % volume(TIME_N) < 0.99_rp * cavities(icavi) % end_dia_vol)) then

                       cavities(icavi) % phase = 3
                       cycles(id) % last_phase_change=cutim
                    ! --------------------------------------------------------------------------
                    ! -------------< isovol relax(3)->filling(4) >------------------------------
                    ! CONDITIONS:
                    !   * in phase 3
                    !   * Pressure is under threashold
                    !
                    elseif ((cavities(icavi) % phase == 3) .and. &
                            (cavities(icavi) % pres(TIME_N) < p_fill))then

                       cavities(icavi) % phase = 4
                       cycles(id) % last_phase_change=cutim
                    ! --------------------------------------------------------------------------
                    ! -------------< filling(4)->isovol contract(1) >---------------------------
                    ! CONDITIONS:
                    !   * in phase 4
                    !   * Geometry is contracting
                    !
                    elseif ((cavities(icavi) % phase == 4) .and. &
                            delta_phase_change.gt.0.01_rp .and.    &
                            (cavities(icavi) % dvol(TIME_N) < -0.000001_rp) )then

                       cavities(icavi) % phase = 1
                       cycles(id) % last_phase_change=cutim
                       cycles(id) % n_beats = cycles(id) % n_beats +1
                    end if

            end select
            ! -------------------------------
        end subroutine sld_cardiac_phases

        ! ---------------------------------------------------------------------
        !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @author   Mariano Vázquez (mariano.vazquez@bsc.es)
        !> @details  Calculation of the volume of a cavity cutted by the vasal plane:
        !>           V = Int_{dOmega_0} (x \cdot e_1) (e_1 \cdot n) d dOmega_0    
        !>           where, x : current coordinates
        !>                  e : vector in the basal plane
        !>                  n : vector normal to the cavity surface
        !>           Refs: Levrero-Florencio et al (2020) DOI: 10.1016/j.cma.2019.112762
        !>                 Quarteroni et al. (2017)       DOI: 10.1016/j.cma.2016.05.031 
        ! ---------------------------------------------------------------------
        subroutine sld_cardiac_cycle_compute_cavity_volume( &
            itask )
            use def_kintyp_basic,     only :  ip, rp, lg
            use def_master,           only :  ittim, displ
            use def_master,           only :  INOTMASTER, ITER_K, TIME_N, TIME_N_MINUS_1, TIME_N_MINUS_2, ITASK_ENDITE
            use def_domain,           only :  mnode, mnodb, ltypb, nnode, ngaus, lelbo, ltype, nboun, npoin_own
            use def_domain,           only :  elmar, ndimb, lbset, lnodb, coord, lnods
            use def_domain,           only :  htable_lninv_loc
            use mod_htable,           only :  htalid
            use mod_communications,   only :  PAR_SUM  
            !-----------------------------
            implicit none
            !-----------------------------
            ! -------------------------------
            integer(ip),  intent(in)       :: itask
            ! -------------------------------
            integer(ip)                    :: icavi, ielem, ipoin, idime, iboun, igaub, inodb, ii
            integer(ip)                    :: pelty, pblty, pnodb, pgaub, pnode, inode
            real(rp)                       :: gpel(3), v1(3), nrm
            real(rp)                       :: baloc(3,3), bocod(3,mnodb), elcod(3,mnode), eucta
            ! -------------------------------
            !
            ! Initialise  
            !
            do icavi= 1, max_cavities
               cavities(icavi) % volume(ITER_K) = 0.0_rp
               cavities(icavi) % plane_points(1:3,1:3) = 0.0_rp
            enddo

            !
            ! Locate the nodes of the plane in the domain
            !
            if( INOTMASTER ) then
                
                do icavi= 1, max_cavities
                    if( cavities(icavi) % initialized .and. cavities(icavi) % plane_nodset(1) > 0_ip )then
                        do ii = 1, 2
                            ipoin = htalid(htable_lninv_loc,cavities(icavi) % plane_nodset(ii))
                            if( ipoin > 0_ip .and. ipoin <= npoin_own ) then
                                cavities(icavi) % plane_points(:,ii) = coord(:,ipoin) + displ(:,ipoin,1)
                            end if
                        enddo
                    endif
                end do

            end if

            ! mod_parall global_to_local

            do icavi = 1, max_cavities
                do ii = 1, 3
                    call PAR_SUM( 3_ip, cavities(icavi) % plane_points(:,ii) )
                enddo
           enddo

            !
            ! Compute the new origin and normal vector of the plane
            !
            do icavi= 1, max_cavities
                if( cavities(icavi) % initialized .and. cavities(icavi) % plane_nodset(1) > 0_ip  )then
                    v1(:) = cavities(icavi) % plane_points(:,1) - cavities(icavi) % plane_points(:,2) 
                    nrm = sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
                    cavities(icavi) % plane_vector(:) = v1(:)/nrm
                endif
            enddo

            !
            ! Compute the volume of the cavity
            !
            if( INOTMASTER ) then

                do iboun = 1, nboun
                    do icavi= 1, max_cavities
                        if( cavities(icavi) % initialized )then
                            if( lbset(iboun) == cavities(icavi) % bouset )then

                                pblty=ltypb(iboun) 
                                pnodb=nnode(pblty)
                                pgaub=ngaus(pblty)
                                ielem=lelbo(iboun)
                                pelty=ltype(ielem)
                                pnode=nnode(pelty)

                                ! Gather bondary and element coordiantes
                                do inodb=1,pnodb
                                    ipoin=lnodb(inodb,iboun)
                                    do idime=1, 3
                                        bocod(idime,inodb) = coord(idime,ipoin) + displ(idime,ipoin,1)
                                    end do
                                end do

                                do inode=1,pnode
                                    ipoin=lnods(inode,ielem)
                                    do idime=1, 3
                                        elcod(idime,inode) = coord(idime,ipoin) + displ(idime,ipoin,1)
                                    end do
                                end do

                                ! Divergence theorem
                                do igaub = 1, pgaub

                                    call bouder(pnodb,3_ip,ndimb,elmar(pblty)%deriv(1_ip,1_ip,igaub),bocod,baloc,eucta)    
                                    call chenor(pnode,baloc,bocod,elcod)

                                    gpel = 0.0_rp
                                    do inodb = 1_ip,pnodb
                                        do idime= 1_ip, 3_ip
                                            gpel(idime) = gpel(idime) + elmar(pblty)%shape(inodb,igaub) * &
                                            bocod(idime,inodb)
                                        end do
                                    end do

                                    cavities(icavi) % volume(ITER_K) =  cavities(icavi) % volume(ITER_K) &
                                        &    - dot_product(gpel(:), cavities(icavi) % plane_vector(:)) &
                                        &    * dot_product(cavities(icavi) % plane_vector(:), baloc(:,3)) &
                                        &    * elmar(pblty) % weigp(igaub) * eucta

                                enddo
                            end if
                        endif
                    end do
                end do 

            endif

            do icavi = 1, max_cavities
                call PAR_SUM( cavities(icavi) % volume(ITER_K) )
            enddo
            
            !
            ! From solidz to master arrays
            !
            do icavi = 1, max_cavities 

                !update volume, only at the end of the iterations loop
                if (itask == ITASK_ENDITE) then
                    cavities(icavi) % volume(TIME_N_MINUS_2) = cavities(icavi) % volume(TIME_N_MINUS_1)  !V^n-2
                    cavities(icavi) % volume(TIME_N_MINUS_1) = cavities(icavi) % volume(TIME_N)          !V^n-1
                    cavities(icavi) % volume(TIME_N)         = cavities(icavi) % volume(ITER_K)          !V^n

                    if (ittim == 1) then
                        cavities(icavi) % volume(TIME_N_MINUS_2) = cavities(icavi) % volume (TIME_N)
                        cavities(icavi) % volume(TIME_N_MINUS_1) = cavities(icavi) % volume (TIME_N)      !initially V^n-1 = V^n = V^i 
                    end if
                end if

            end do
            ! -------------------------------
        end subroutine sld_cardiac_cycle_compute_cavity_volume

end module mod_sld_cardiac_cycle
