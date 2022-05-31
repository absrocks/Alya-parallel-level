!-------------------------------------------------------------------------------
!> @addtogroup Phase Field
!> @{
!> @authors Adria Quintanas-Corominas : adria.quintanas@bsc.es
!> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
!> @date    March, 2020
!> @}
!------------------------------------------------------------------------------
module mod_exm_sld_eccoupling
    ! =========================================================================
    use def_kintyp_basic, only                :  ip, rp, lg
    use def_domain, only                      :  ndime, nmate
    use def_master, only                      :  cutim, ittim, dtime
    ! ----------------------------------------
    implicit none
    ! ----------------------------------------
    ! -< DIMENSIONS >-------------------------
    integer(ip), parameter                    :: kfl_exm_max_nmaterials = 10_ip
    integer(ip), parameter                    :: nmate_ecc = 10_ip
    integer(ip), parameter                    :: nprop_ecc = 17_ip
    integer(ip), parameter                    :: nsdvs_ecc =  7_ip
    integer(ip), parameter                    :: ncomp_ecc =  2_ip

    ! -< CELL MODELS >------------------------
    integer(ip), parameter                    :: EXMSLD_CELL_NO_CELL  = 0_ip
    integer(ip), parameter                    :: EXMSLD_CELL_FENTON   = 2_ip
    integer(ip), parameter                    :: EXMSLD_CELL_TT2006   = 4_ip
    integer(ip), parameter                    :: EXMSLD_CELL_OHARA    = 5_ip
    integer(ip), parameter                    :: EXMSLD_CELL_SCATRIA  = 6_ip
    integer(ip), parameter                    :: EXMSLD_CELL_SCVENTRI = 7_ip

    ! -< ELECTRO-MECHANICAL MODELS >----------
    integer(ip), parameter                    :: EXMSLD_EMEC_NO_EMEC    = 0_ip
    integer(ip), parameter                    :: EXMSLD_EMEC_HUNTER     = 1_ip
    integer(ip), parameter                    :: EXMSLD_EMEC_LAND       = 3_ip
    integer(ip), parameter                    :: EXMSLD_EMEC_LAND_BIDIR = 4_ip
    
    ! -< ELECTRO-MECHANICAL OPTIONS >---------
    integer(ip), parameter                    :: EMEC_COPT_ISOTROPIC = 1_ip
    integer(ip), parameter                    :: EMEC_COPT_TRANSISOT = 2_ip

    ! ----< GENERAL VARIABLES >---------------
    logical(lg)                               :: kfl_exmsld_3Dcou_ecc
    logical(lg)                               :: has_land
    logical(lg)                               :: kfl_print

    ! -----< WHICH MODULE IS CALLING >--------
    integer(ip), parameter                    :: SOLIDZ = 1_ip
    integer(ip), parameter                    :: EXMEDI = 2_ip  

    ! ----< POINTERS TO VARIABLES >-----------
    integer(ip), dimension(:),       pointer  :: law_cell_ecc
    integer(ip), dimension(:),       pointer  :: law_emec_ecc
    integer(ip), dimension(:),       pointer  :: kfl_copt_ecc
    real(rp),    dimension(:),       pointer  :: cocof_ecc 
    real(rp),    dimension(:),       pointer  :: timec_ecc 
    real(rp),    dimension(:),       pointer  :: hillc_ecc
    real(rp),    dimension(:),       pointer  :: cal50_ecc
    real(rp),    dimension(:),       pointer  :: ortk1_ecc
    real(rp),    dimension(:),       pointer  :: ortk2_ecc
    real(rp),    dimension(:),       pointer  :: trans_ecc
    real(rp),    dimension(:,:),     pointer  :: props_ecc
    real(rp),    dimension(:,:),     pointer  :: ortho_ecc
    real(rp),    dimension(:,:,:),   pointer  :: state_ecc
    real(rp),    dimension(:,:,:),   pointer  :: strch_ecc    
    real(rp),    dimension(:,:),     pointer  :: troponin_ecc 
    real(rp),    dimension(:,:),     pointer  :: calcium_ecc 
    real(rp),    dimension(:,:),     pointer  :: fiber_ecc 
    real(rp),    dimension(:,:),     pointer  :: cell_ca0_ecc
    ! ----------------------------------------
    private                                   :: get_active_cauchy_stress_tensor, &
    &                                            get_active_second_elasticity_tensor, &
    &                                            get_active_traction_HUNTER, &
    &                                            ROUNDOFF, &
    &                                            MULT_MxV, &
    &                                            MULT_MxT, &
    &                                            MULT_VxMxV, &
    &                                            OUT_VxV, &
    &                                            OUT_VxVxVxV, &
    &                                            INV3x3
    public                                    :: set_has_land, &
                                                 has_exmsld_coupling, &
                                                 get_active_traction_LANDS

    contains

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief   Manage arrays 
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine exm_sld_ecc_manage_arrays( &
        itask, narray, wtask )
        ! -------------------------------
        use def_domain,           only :  npoin, nelem, mgaus, nmate
        use mod_arrays,           only :  arrays, arrays_number, arrays_register 
        use def_master,           only :  kfl_eccty
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip),      intent(in)           :: itask
        integer(ip),      intent(in), optional :: narray
        character(len=*), intent(in), optional :: wtask
        ! -------------------------------
        integer(ip),      parameter            :: SOLIDZ_REGISTER = 0_ip
        integer(ip),      parameter            :: SOLIDZ_ARRAYS   = 1_ip
        integer(ip),      parameter            :: EXMEDI_REGISTER = 2_ip
        integer(ip),      parameter            :: EXMEDI_ARRAYS   = 3_ip
        ! -------------------------------
        ! -------------------------------
        integer(ip)                      :: imate
        select case(itask)

            case( SOLIDZ_REGISTER )
                call arrays_register( narray + 1_ip, (/ 'CALCI', 'SCALA', 'NPOIN', 'SECON' /))
                call arrays_register( narray + 2_ip, (/ 'EMSVD', 'SCALA', 'NELEM', 'PRIMA' /), state_ecc,    ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
                call arrays_register( narray + 3_ip, (/ 'EMSTR', 'SCALA', 'NELEM', 'PRIMA' /), strch_ecc,    ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
                call arrays_register( narray + 4_ip, (/ 'EMTRO', 'SCALA', 'NPOIN', 'PRIMA' /), troponin_ecc, ENTITY_POSITION= 1_ip, TIME_POSITION= 2_ip )
                call arrays_register( narray + 5_ip, (/ 'EMFIB', 'SCALA', 'NELEM', 'PRIMA' /), fiber_ecc,    ENTITY_POSITION= 2_ip)

           case( SOLIDZ_ARRAYS )
               call arrays( arrays_number('EMSVD'), wtask, state_ecc,    nsdvs_ecc*mgaus, nelem, ncomp_ecc)
               call arrays( arrays_number('EMSTR'), wtask, strch_ecc,    mgaus,           nelem, ncomp_ecc)
               call arrays( arrays_number('EMTRO'), wtask, troponin_ecc, npoin,           ncomp_ecc)
               
               ! TODO: Unify gpfib and fiber_ecc by creating a new way to register (nsdv, ngaus, nelem) arrays.
               if( wtask == 'WRITE RESTART') call set_fiber_ecc() 
               call arrays( arrays_number('EMFIB'), wtask, fiber_ecc,    ndime*mgaus,     nelem)
               if( wtask == 'READ RESTART')  call get_fiber_ecc()

            case( EXMEDI_REGISTER )
                call set_has_land()
                
                if(kfl_exmsld_3Dcou_ecc .and. has_land) then
                   call arrays_register( narray + 1_ip, (/ 'EMSVD', 'SCALA', 'NELEM', 'PRIMA' /), state_ecc,    ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
                   call arrays_register( narray + 3_ip, (/ 'EMTRO', 'SCALA', 'NPOIN', 'PRIMA' /), troponin_ecc, ENTITY_POSITION= 1_ip, TIME_POSITION= 2_ip )
                endif
            case( EXMEDI_ARRAYS )

        end select
        ! -------------------------------
        contains 

            subroutine set_fiber_ecc()
                use def_master, only : gpfib
                use def_domain, only : ndime, nelem, mgaus
                implicit none
                integer(ip)         :: ielem, ii, jj, kk 
                do ielem = 1, nelem
                    do ii = 1, mgaus
                        do kk = 1, ndime 
                           jj = (ii - 1_ip)*nsdvs_ecc + kk
                           fiber_ecc(jj,ielem) = gpfib(kk,ii,ielem)
                        enddo
                    enddo
                enddo    
            end subroutine set_fiber_ecc

            subroutine get_fiber_ecc()
                use def_master, only : gpfib 
                use def_domain, only : ndime, nelem, mgaus
                implicit none
                integer(ip)         :: ielem, ii, jj, kk
                do ielem = 1, nelem
                    do ii = 1, mgaus
                        do kk = 1, ndime
                           jj = (ii - 1_ip)*nsdvs_ecc + kk
                           gpfib(kk,ii,ielem) = fiber_ecc(jj,ielem)
                        enddo
                    enddo
                enddo
            end subroutine get_fiber_ecc

    end subroutine exm_sld_ecc_manage_arrays


    ! ---------------------------------------------------------------------
    !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief    Write information
    !> @details  
    !> @todo     
    ! ---------------------------------------------------------------------
    subroutine exm_sld_ecc_manage_restart( &
        itask )
        use def_master, only            :  INOTSLAVE, IPARALL, momod, modul, kfl_paral
        use def_master, only            :  ITASK_READ_RESTART, ITASK_WRITE_RESTART
        use mod_communications, only    :  PAR_BROADCAST
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: itask
        ! -------------------------------
        integer(ip)                    :: ii, jj
        ! -------------------------------

        if( itask == ITASK_READ_RESTART ) then
            !
            ! Read restart
            !
            if( INOTSLAVE )then

                do ii = 1, nmate_ecc
                    read(momod(modul) % lun_rstar) cell_ca0_ecc(1:3,ii) 
                enddo 
 
            endif

            if( IPARALL )then
             
                do ii = 1, nmate_ecc
                    do jj = 1, 3
                        call PAR_BROADCAST(cell_ca0_ecc(jj,ii))
                    enddo
                enddo
                
            endif


        else if( itask == ITASK_WRITE_RESTART ) then
            !
            ! Write restart file
            !
            if( INOTSLAVE )then

                do ii = 1, nmate_ecc
                    write(momod(modul) % lun_rstar) cell_ca0_ecc(1:3,ii) 
                enddo

            endif

        endif

    end subroutine exm_sld_ecc_manage_restart


    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief    Allocate memmory
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine exm_sld_ecc_allocate_memmory( &
        itask )
        ! -------------------------------
        use def_master,           only :  modul, mem_modul, IMASTER, INOTMASTER
        use def_master,           only :  kfl_eccty
        use def_domain,           only :  npoin, nelem, mgaus, nmate
        use def_master,           only :  gpfib 
        use mod_memory,           only :  memory_alloca
        use mod_memchk,           only :  zero, memchk
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: itask
        integer(ip)                    :: istat
        integer(ip), parameter         :: SOLIDZ_MEMPHY = 100_ip
        integer(ip), parameter         :: SOLIDZ_SOLMEM = 101_ip
        integer(ip), parameter         :: EXMEDI_MEMPHY = 102_ip
        integer(ip), parameter         :: EXMEDI_SOLMEM = 200_ip
        integer(ip), parameter         :: EXMEDI = 2_ip
        ! -------------------------------
        integer(ip)                      :: imate

        select case(itask)
            
            case( EXMEDI_SOLMEM )
       
                call memory_alloca(mem_modul(1:2,modul), 'CELL_CA0_ECC', 'exm_sld_ecc', cell_ca0_ecc, 3_ip, nmate_ecc) 

            case( SOLIDZ_MEMPHY )

                call memory_alloca(mem_modul(1:2,modul), 'KFL_COPT_ECC', 'exm_sld_ecc', kfl_copt_ecc, nmate_ecc)
                call memory_alloca(mem_modul(1:2,modul), 'PROPS_ECC'   , 'exm_sld_ecc', props_ecc   , nprop_ecc, nmate_ecc)
                call memory_alloca(mem_modul(1:2,modul), 'ORTHO_ECC'   , 'exm_sld_ecc', ortho_ecc   , 3_ip     , nmate_ecc)

            case( SOLIDZ_SOLMEM )

                if( IMASTER ) then
                    call memory_alloca(mem_modul(1:2,modul), 'CALCIUM_ECC' , 'exm_sld_ecc', calcium_ecc , 1_ip     , 1_ip)
                    call memory_alloca(mem_modul(1:2,modul), 'GPFIB', 'exm_sld_ecc', gpfib, ndime, mgaus, 1_ip)
                    if(kfl_exmsld_3Dcou_ecc) then
                        call memory_alloca(mem_modul(1:2,modul), 'CELL_CA0_ECC', 'exm_sld_ecc', cell_ca0_ecc, 3_ip, nmate_ecc) 
                    endif

                else
                    call memory_alloca(mem_modul(1:2,modul), 'CALCIUM_ECC' , 'exm_sld_ecc', calcium_ecc , 1_ip     , npoin)
                    call memory_alloca(mem_modul(1:2,modul), 'GPFIB', 'exm_sld_ecc', gpfib, ndime, mgaus, nelem)
                    if(kfl_exmsld_3Dcou_ecc) then
                        call memory_alloca(mem_modul(1:2,modul), 'CELL_CA0_ECC', 'exm_sld_ecc', cell_ca0_ecc, 3_ip, nmate_ecc) 
                    endif
                end if

                if( INOTMASTER )then
                    strch_ecc(:,:,:)   = 1.0_rp
                endif

            case( EXMEDI_MEMPHY )

        end select
        ! -------------------------------
    end subroutine exm_sld_ecc_allocate_memmory

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief    Read input variables concerning the coupling
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine exm_sld_ecc_read_data( &
        itask, imate )
        use def_master,           only :  kfl_eccty
        use def_inpout,           only :  words, param, exists, getint, getrea, getcha
        use mod_ecoute,           only :  ecoute
        use mod_messages,         only :  livinf
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: itask
        integer(ip), intent(in)        :: imate
        real(rp)                       :: cocof_ecc
        real(rp)                       :: timec_ecc
        real(rp)                       :: hillc_ecc
        real(rp)                       :: cal50_ecc
        real(rp)                       :: ortk1_ecc
        real(rp)                       :: ortk2_ecc
        real(rp)                       :: trans_ecc
        ! -------------------------------
        select case(itask)

            case( SOLIDZ )

                ! Initialise some variables
                kfl_copt_ecc(imate) = EMEC_COPT_ISOTROPIC
                cocof_ecc = 1.00_rp
                timec_ecc = 0.06_rp
                hillc_ecc = 3.00_rp
                cal50_ecc = 0.50_rp
                ortk1_ecc = 0.00_rp
                ortk2_ecc = 0.00_rp
                trans_ecc = 0.00_rp

                call ecoute('exm_sld_ecc')

                coupling_sld_exm: &
                do while( words(1)/='ENDCO' )

                    ! Electro-mechanical model
                    if( words(1)=='MODEL' )then

                        if(    words(2)=='HUNTE' )then
                            kfl_eccty(imate) = EXMSLD_EMEC_HUNTER

                        elseif( words(2)=='LAND')then
                            kfl_eccty(imate) = EXMSLD_EMEC_LAND
                            if( words(3)=='BIDIR' )then
                                kfl_eccty(imate)= EXMSLD_EMEC_LAND_BIDIR
                            endif

                       end if

                    end if

                    ! Control coupling factor
                    if( words(1)=='CONTR' )then
                        cocof_ecc = param(1)
                    endif

                    ! Time constant for the [Ca] fct
                    if( words(1)=='TIMEC' )then
                        timec_ecc = param(1)
                    endif

                    ! Hill coefficient
                    if( words(1)=='HILLC' )then
                        hillc_ecc = param(1)
                    endif

                    ! Ca 50 coefficient
                    if( words(1)=='CAL50' )then
                        cal50_ecc = param(1)
                    endif

                    ! First coefficient for the orthotropic activation
                    if( words(1)=='ORTK1' )then
                        ortk1_ecc= param(1)
                    endif

                    ! First coefficient for the orthotropic activation
                    if( words(1)=='ORTK2' )then
                        ortk2_ecc = param(1)
                    endif

                    ! Transversal coupling force ratio
                    if(words(1)=='TRISO') then
                        kfl_copt_ecc(imate) = EMEC_COPT_TRANSISOT
                        trans_ecc = param(1) 
                    end if

                    call ecoute('exm_sld_ecc')

                end do &
                coupling_sld_exm

                ! Properties of the coupling models
                props_ecc(:,imate) = 0.0_rp
                if(     kfl_eccty(imate) == EXMSLD_EMEC_HUNTER )then
                    props_ecc( 1,imate) = 1.00_rp            ! Cam 
                    props_ecc( 2,imate) = cal50_ecc
                    props_ecc( 3,imate) = timec_ecc
                    props_ecc( 4,imate) = 1.45_rp            ! eta
                    props_ecc( 5,imate) = hillc_ecc
                    props_ecc( 6,imate) = 0.30_rp            ! k
                    props_ecc( 7,imate) = cocof_ecc*1.0e6_rp

                elseif( kfl_eccty(imate) == EXMSLD_EMEC_LAND        .or. &
                &       kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR )then
                    props_ecc( 1,imate) = 100.0_rp           ! k_TRPN    
                    props_ecc( 2,imate) = 2.0_rp             ! n_TRPN   
                    props_ecc( 3,imate) = 0.805_rp           ! Ca50_ref  
                    props_ecc( 4,imate) = 1000.0_rp          ! k_u      
                    props_ecc( 5,imate) = 5.0_rp             ! n_tm     
                    props_ecc( 6,imate) = 0.35_rp            ! TRPN_50     
                    props_ecc( 7,imate) = 182.0_rp           ! k_uw         
                    props_ecc( 8,imate) = 12.0_rp            ! k_ws       
                    props_ecc( 9,imate) = 0.5_rp             ! r_w        
                    props_ecc(10,imate) = 0.25_rp            ! r_s         
                    props_ecc(11,imate) = 0.0085_rp          ! y_s          
                    props_ecc(12,imate) = 0.615_rp           ! y_w          
                    props_ecc(13,imate) = 2.23_rp            ! phi         
                    props_ecc(14,imate) = 25.0_rp            ! A_eff       
                    props_ecc(15,imate) = 2.3_rp             ! beta_0     
                    props_ecc(16,imate) = -2.4_rp            ! beta_1      
                    props_ecc(17,imate) = cocof_ecc*120.0_rp ! T_ref    
                endif

                ! Properties of the orthotropic coupling
                ortho_ecc( 1,imate) = trans_ecc
                ortho_ecc( 2,imate) = ortk1_ecc
                ortho_ecc( 3,imate) = ortk2_ecc

        end select 
        ! -------------------------------
    end subroutine exm_sld_ecc_read_data

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief    Send data
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine exm_sld_ecc_send_data( &
        itask )
        use def_master,            only : kfl_eccty
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: itask
        integer(ip)                    :: ii, jj
        ! -------------------------------

        call lexcha(kfl_exmsld_3Dcou_ecc)

        select case(itask)

            case(SOLIDZ)
            
                do ii = 1, nmate_ecc
                    call iexcha(kfl_eccty(ii)) 
                    call iexcha(kfl_copt_ecc(ii))
                    do jj = 1, nprop_ecc
                        call rexcha(props_ecc(jj,ii))
                    enddo
                    do jj = 1, 3
                        call rexcha(ortho_ecc(jj,ii))
                    enddo
                enddo
                call set_has_land()
                call lexcha(has_land)

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
    end subroutine exm_sld_ecc_send_data

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief   
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine exm_sld_ecc_manage_state_variables( &
        itask )
        use def_master,            only : ITASK_ENDSTE
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: itask
        ! -------------------------------
        select case(itask)

            case( ITASK_ENDSTE )

                ! ITER_AUX = 2; ITER_K = 1
                state_ecc(:,:,1) = state_ecc(:,:,2) 
                strch_ecc(:,:,1) = strch_ecc(:,:,2)

        end select
        ! -------------------------------
    end subroutine exm_sld_ecc_manage_state_variables

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief   
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine exm_sld_ecc_initialize_troponin( &
        )
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip)                    :: ii
        ! -------------------------------
        if(associated(troponin_ecc))then
            do ii = 1, size(troponin_ecc,DIM=1_ip,KIND=ip)
                troponin_ecc(:,1) = 0.0_rp
            enddo
        endif
        ! -------------------------------
    end subroutine exm_sld_ecc_initialize_troponin

    subroutine exm_sld_ecc_interchange_troponin( &
        )
        use def_domain,            only : vmass
        ! -------------------------------
        implicit none
        ! -------------------------------
        call rhsmod(1_ip,troponin_ecc(:,1))
        troponin_ecc(:,1) = troponin_ecc(:,1)/vmass(:)
        ! -------------------------------
    end subroutine exm_sld_ecc_interchange_troponin

    subroutine exm_sld_ecc_get_troponin_at_gp( &
        ielem, pgaus, gp_trop )
        use mod_matrix,            only : matrix_assemble_element_RHS
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: ielem
        integer(ip), intent(in)        :: pgaus
        real(rp),    intent(out)       :: gp_trop(pgaus)
        ! -------------------------------
        integer(ip)                    :: ii, jj, igaus
        ! -------------------------------
        do ii = 1, pgaus
            jj = (ii - 1_ip)*nsdvs_ecc + 3_ip 
            gp_trop(ii) = state_ecc(jj,ielem,2)        
        enddo
        ! -------------------------------
    end subroutine exm_sld_ecc_get_troponin_at_gp

    subroutine exm_sld_ecc_assemble_troponin( &
        ielem, pnode, pgaus, lnods, gp_N, gp_vol, gp_Trop )
        use mod_matrix,            only : matrix_assemble_element_RHS
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: ielem
        integer(ip), intent(in)        :: pnode
        integer(ip), intent(in)        :: pgaus
        integer(ip), intent(in)        :: lnods(pnode)
        real(rp),    intent(in)        :: gp_vol(pgaus)
        real(rp),    intent(in)        :: gp_N(pnode,pgaus)
        real(rp),    intent(in)        :: gp_trop(pgaus)
        ! -------------------------------
        integer(ip)                    :: inode, igaus
        real(rp)                       :: el_trop(pnode)
        ! -------------------------------
        !
        ! Integrate troponin
        !
        el_trop(:) = 0.0_rp
        do inode = 1, pnode
            do igaus = 1, pgaus
                el_trop(inode) = el_trop(inode) + gp_N(inode,igaus)*gp_vol(igaus)*gp_trop(igaus)
            end do
        end do
        !
        ! Assemble troponin
        !
        call matrix_assemble_element_RHS(1_ip,1_ip,pnode,lnods(:),el_trop(:),troponin_ecc(:,1))
        ! -------------------------------
    end subroutine exm_sld_ecc_assemble_troponin

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief 
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine exm_sld_ecc_set_fibers_at_gp( &
        ielem, gp_fib )
        use def_master,            only : gpfib
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: ielem 
        real(rp),    intent(in)        :: gp_fib(:,:)
        ! -------------------------------
        gpfib(:,:,ielem) = gp_fib(:,:)
        ! -------------------------------        
    end subroutine exm_sld_ecc_set_fibers_at_gp

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief 
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine exm_sld_ecc_get_calcium_at_gp( &
        ielem, pnode, pgaus, lnods, gp_N, gp_Ca )
        use def_master,            only : vconc
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: ielem
        integer(ip), intent(in)        :: pnode
        integer(ip), intent(in)        :: pgaus
        integer(ip), intent(in)        :: lnods(pnode)
        real(rp),    intent(in)        :: gp_N(pnode,pgaus)
        real(rp),    intent(out)       :: gp_Ca(pgaus)
        ! -------------------------------
        integer(ip)                    :: ipoin, igaus, inode
        real(rp)                       :: el_Ca(pnode)
        ! -------------------------------
        !
        ! Gather calciuml
        !
        do inode = 1, pnode
            ipoin = lnods(inode)
            if(kfl_exmsld_3Dcou_ecc) then
                el_Ca(inode) = calcium_ecc(1,ipoin)*1000_rp
            else
                el_Ca(inode) = vconc(1,ipoin,1)*1000_rp
            endif
        enddo

        !
        ! Interpolate calcium
        !
        gp_Ca(:) = 0.0_rp
        do igaus = 1, pgaus
            do inode = 1, pnode
                gp_Ca(igaus) = gp_Ca(igaus) + gp_N(inode,igaus)*el_Ca(inode)
            end do
        end do
        ! -------------------------------
    end subroutine exm_sld_ecc_get_calcium_at_gp

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief    
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine exm_sld_ecc_add_active_stress_and_moduli( &
        itask, ielem, imate, pgaus, &
        id_cell, id_ecct, id_copt, &
        gp_det_J, gp_F, gp_C, gp_Ca, &
        gp_fbr_0, gp_sht_0, gp_nrm_0, &
        gp_str, gp_pk1, gp_dds, gp_tmo )
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: itask                 ! 
        integer(ip), intent(in)        :: ielem                 ! Number of element
        integer(ip), intent(in)        :: imate                 ! Number of material
        integer(ip), intent(in)        :: pgaus                 ! Number of gauss points
        integer(ip), intent(in)        :: id_cell               ! Id for the cell model
        integer(ip), intent(in)        :: id_ecct               ! Id for the electro-mechanical model
        integer(ip), intent(in)        :: id_copt               ! Id for the coupling model
        real(rp),    intent(in)        :: gp_det_J(pgaus)       ! Determinant of the Jacobian at gauss points
        real(rp),    intent(in)        :: gp_F(3,3,pgaus)       ! Deformation Gradient at gauss points
        real(rp),    intent(in)        :: gp_C(3,3,pgaus)       ! Right-Cauchy tensor at gauss points
        real(rp),    intent(in)        :: gp_Ca(pgaus)          ! Calcium concentration at gauss points
        real(rp),    intent(in)        :: gp_fbr_0(3,pgaus)     ! Fiber direction in the reference config. at gauss points
        real(rp),    intent(in)        :: gp_sht_0(3,pgaus)     ! Sheat direction in the reference config. at gauss points
        real(rp),    intent(in)        :: gp_nrm_0(3,pgaus)     ! Norma direction in the reference config. at gauss points
        real(rp),    intent(inout)     :: gp_str(3,3,pgaus)     ! Second Piola-Kirchhoff at gauss points
        real(rp),    intent(inout)     :: gp_pk1(3,3,pgaus)     ! First Piola-Kirchhoff at gauss points
        real(rp),    intent(inout)     :: gp_dds(3,3,3,3,pgaus) ! Second elasticity tangent moduli at gauss points
        real(rp),    intent(inout)     :: gp_tmo(3,3,3,3,pgaus) ! First elasticity tangent moduli at gauss points
        ! -------------------------------
        integer(ip)                    :: igaus                 ! Counters for the loops
        integer(ip)                    :: i, j, k, l, m, n, ii, jj 
        real(rp)                       :: i4f, i4s, i4n         ! Kinematics invariants
        real(rp)                       :: gp_lmd(3,pgaus)       ! Stretches at gauss points
        real(rp)                       :: gp_fbr_c(3)           ! Fiber direction in the current config at the gauss points
        real(rp)                       :: gp_sht_c(3)
        real(rp)                       :: gp_nrm_c(3)
        real(rp)                       :: gp_fbr_t(3,pgaus)     ! Unit fiber direction in the current config. at the gauss points 
        real(rp)                       :: gp_sht_t(3,pgaus)
        real(rp)                       :: gp_nrm_t(3,pgaus)
        real(rp)                       :: statvar(10,2)         ! Local state variables
        real(rp)                       :: gp_F_inv(3,3,pgaus)   ! Inverse of the deformation gradient tesnor at the gauss points
        real(rp)                       :: gp_Tac_f(pgaus)       ! Active tractions at the guass points
        real(rp)                       :: gp_dTac_f(pgaus)      ! Active tractions at the guass points
        real(rp)                       :: gp_sig(3,3)           ! Cauchy stress tensor at the gauss points
        real(rp)                       :: gp_dsig(3,3,3,3)
        real(rp)                       :: gp_str_act(3,3,pgaus)    
        real(rp)                       :: gp_pk1_act(3,3,pgaus)    
        real(rp)                       :: gp_dds_act(3,3,3,3,pgaus)
        real(rp)                       :: gp_tmo_act(3,3,3,3,pgaus)
        ! -------------------------------
        integer(ip),  parameter        :: ASSEMBLE_MODULI = 2_ip ! It comes from solid
        real(rp), parameter, dimension(3,3) :: I_3x3 = reshape((/ &
                                        &   1.0, 0.0, 0.0,   &
                                        &   0.0, 1.0, 0.0,   &
                                        &   0.0, 0.0, 1.0    &
                                        &   /), (/3 , 3/))
        ! -------------------------------
        ! Compute stretches
        loop_gp_stretches: &
        do igaus = 1, pgaus
            
            ! Compute the invariants
            i4f = MULT_VxMxV(gp_fbr_0(:,igaus),3_ip,gp_C(:,:,igaus),gp_fbr_0(:,igaus),3_ip)
            i4s = MULT_VxMxV(gp_sht_0(:,igaus),3_ip,gp_C(:,:,igaus),gp_sht_0(:,igaus),3_ip)
            i4n = MULT_VxMxV(gp_nrm_0(:,igaus),3_ip,gp_C(:,:,igaus),gp_nrm_0(:,igaus),3_ip)

            ! Avoid round-off errors
            i4f = ROUNDOFF(i4f,1.0_rp)
            i4s = ROUNDOFF(i4s,1.0_rp)
            i4n = ROUNDOFF(i4n,1.0_rp)

            ! Calculate the stretch
            gp_lmd(1,igaus) = sqrt(i4f)
            gp_lmd(2,igaus) = sqrt(i4s)
            gp_lmd(3,igaus) = sqrt(i4n)

        end do &
        loop_gp_stretches

        ! Compute unitary current fiber direction
        loop_gp_update_fibers: &
        do igaus = 1, pgaus

            ! Current directions (f=F'*f0 ; s=F'*s0 ; n=F'*n0)
            gp_fbr_c(:) = MULT_MxV(gp_F(:,:,igaus),3_ip,3_ip,gp_fbr_0(:,igaus),3_ip)
            gp_sht_c(:) = MULT_MxV(gp_F(:,:,igaus),3_ip,3_ip,gp_sht_0(:,igaus),3_ip)
            gp_nrm_c(:) = MULT_MxV(gp_F(:,:,igaus),3_ip,3_ip,gp_nrm_0(:,igaus),3_ip)

            ! Updated fiber direction: lambda_f * f_true = F * f_0 = f
            ! => here f is the UNIT fiber direction in the current configuration
            if( gp_lmd(1,igaus) /= 0.0_rp ) then
                gp_fbr_t(:,igaus) = gp_fbr_c(:)/gp_lmd(1,igaus)
            else
                gp_fbr_t(:,igaus) = 0.0_rp
            end if

            if( gp_lmd(2,igaus) /= 0.0_rp ) then
                gp_sht_t(:,igaus) = gp_sht_c(:)/gp_lmd(2,igaus)
            else
                gp_sht_t(:,igaus) = 0.0_rp
            end if

            if( gp_lmd(3,igaus) /= 0.0_rp ) then
                gp_nrm_t(:,igaus) = gp_nrm_c(:)/gp_lmd(3,igaus)
            else
                gp_nrm_t(:,igaus) = 0.0_rp
            end if

        end do &
        loop_gp_update_fibers

        ! Active stress tensor 
        gp_str_act = 0.0_rp 
        gp_pk1_act = 0.0_rp

        loop_gp_active_stresses: &
        do igaus = 1, pgaus
            
            ! Auxiliary indices
            ii = (igaus - 1_ip)*nsdvs_ecc + 1_ip
            jj = ii + (nsdvs_ecc - 1_ip)
            
            ! Get state variables
            statvar(1:7,1) = state_ecc(ii:jj,ielem,1) 
            statvar(  8,1) = strch_ecc(igaus,ielem,1)
            statvar(  9,1) = sum(cell_ca0_ecc(:,imate))*1000.0_rp/3.0_rp
            statvar(:,2)=0.0_rp

            ! Compute active traction according to the electro-mechanical
            ! model
            if(     id_ecct == EXMSLD_EMEC_HUNTER )then
                call get_active_traction_HUNTER( id_cell, props_ecc(:,imate), gp_lmd(1,igaus), &
                &   gp_Ca(igaus), gp_Tac_f(igaus), gp_dTac_f(igaus), statvar(:,:) )
    
            elseif( id_ecct == EXMSLD_EMEC_LAND  .or. &
                    id_ecct == EXMSLD_EMEC_LAND_BIDIR )then
                call get_active_traction_LANDS( props_ecc(:,imate), gp_lmd(1,igaus), gp_Ca(igaus), &
                &   gp_Tac_f(igaus), gp_dTac_f(igaus), statvar(:,:) )
            endif
                    
            ! Set state variables
            state_ecc(ii:jj,ielem,2) = statvar(1:7,2)
            strch_ecc(igaus,ielem,2) = gp_lmd(1,igaus)

            ! Get Cauchy stress tensor from the active tractions
            gp_sig(:,:) = 0.0_rp
            call get_active_cauchy_stress_tensor( &
            &   gp_det_J(igaus), ortho_ecc(:,imate), &
            &   gp_fbr_t(:,igaus), gp_sht_t(:,igaus), gp_nrm_t(:,igaus), &
            &   gp_Tac_f(igaus), gp_dTac_f(igaus), &
            &   gp_sig(:,:), gp_dsig(:,:,:,:) ) 

            ! From Cauchy to Second Piola-Kirchhoff (PK2) stress tensor
            gp_F_inv(:,:,igaus) = INV3x3(gp_F(:,:,igaus))
            do i = 1, 3
                do j = 1, 3
                    do k = 1, 3
                        do l = 1, 3
                            gp_str_act(i,j,igaus) = gp_str_act(i,j,igaus) + &
                            &   gp_det_J(igaus)*gp_F_inv(i,k,igaus)*gp_sig(k,l)*gp_F_inv(j,l,igaus)
                        enddo
                    enddo
                enddo
            enddo

            ! From Second to First Piola-Kirchhoff (PK1) stress tensor
            !   P_jK = F_jI * S_IK
            !   {Ref : Belytschko book}
            do k = 1, 3
                do j = 1, 3
                    do i = 1, 3
                        gp_pk1_act(j,k,igaus) = gp_pk1_act(j,k,igaus) + gp_F(j,i,igaus) * gp_str_act(i,k,igaus)
                    enddo
                enddo
            enddo

        enddo &
        loop_gp_active_stresses

        ! Add active contribution to the passive one
        gp_str = gp_str + gp_str_act
        gp_pk1 = gp_pk1 + gp_pk1_act

        ! Active material moduli
        gp_dds_act = 0.0_rp
        gp_tmo_act = 0.0_rp

        if( itask == ASSEMBLE_MODULI )then

            do igaus = 1, pgaus

                ! Get second elasticity tensor (dS/dE)
                call get_active_second_elasticity_tensor( &
                &   ortho_ecc(:,imate), &
                &   gp_fbr_0(:,igaus), gp_sht_0(:,igaus), gp_nrm_0(:,igaus), &
                &   gp_lmd(:,igaus), gp_Tac_f(igaus), gp_dTac_f(igaus), &
                &   gp_dds_act(:,:,:,:,igaus) )
                
                ! From second to first elasticity tensor (dS/dE -> dP/dF)
                !   dPdF_iJkL = delta_ik * S_JL + F_iM * FkN * dSdE_MJNL
                !   {Ref : Belytschko book}
                do i = 1, 3
                    do j = 1, 3
                        do k = 1, 3 
                            do l = 1, 3
                                do m = 1, 3
                                    do n = 1, 3
                                        gp_tmo_act(i,j,k,l,igaus)= gp_tmo_act(i,j,k,l,igaus) + gp_dds_act(m,j,n,l,igaus) * gp_F(i,m,igaus) * gp_F(k,n,igaus)
                                    enddo
                                enddo
                                gp_tmo_act(i,j,k,l,igaus) = gp_tmo_act(i,j,k,l,igaus) + I_3x3(i,k) * gp_str_act(j,l,igaus)
                            enddo
                        enddo
                    enddo
                enddo

            end do

            ! Add active contribution to the passive one
            gp_dds = gp_dds + gp_dds_act
            gp_tmo = gp_tmo + gp_tmo_act

        endif
        ! -------------------------------
    end subroutine exm_sld_ecc_add_active_stress_and_moduli

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief   Decompose traction according to a basis
    !> @details 
    ! -------------------------------------------------------------------------
    subroutine get_active_cauchy_stress_tensor( &
        det_J, ratios, fib, sht, nrm, tract, dtract, sig, dsig )
        ! -------------------------------
        implicit none
        ! -------------------------------
        real(rp),    intent(in)        :: det_J           ! Det(J)
        real(rp),    intent(in)        :: ratios(3)       ! Transversal ratio
        real(rp),    intent(in)        :: fib(3)          ! Basis {f, s, n}
        real(rp),    intent(in)        :: sht(3)
        real(rp),    intent(in)        :: nrm(3)     
        real(rp),    intent(in)        :: tract           ! Active traction along the fiber direction
        real(rp),    intent(in)        :: dtract          ! Derivative of the active traction with respect to the stretch along the fiber direciton
        real(rp),    intent(out)       :: sig(3,3)        ! Active Cauchy stress tensor
        real(rp),    intent(out)       :: dsig(3,3,3,3)   ! Active Cauchy stress tensor
        ! -------------------------------
        integer(ip)                    :: ii
        real(rp),    dimension(3,3)    :: fxf
        real(rp),    dimension(3,3)    :: sxs
        real(rp),    dimension(3,3)    :: nxn
        real(rp)                       :: tract_f        ! Active traction : fiber direction
        real(rp)                       :: tract_s        ! Active traction : shear direction
        real(rp)                       :: tract_n        ! Active traction : normal direction
        real(rp)                       :: tract_i        ! Active traction : isotropic
        ! -------------------------------
        ! Coupling of the traction according to the direction
        ! - isotropic coupling
        tract_f = tract

        ! - transversally isotropic coupling
        tract_f = tract*(1.0_rp - ratios(1))
        tract_i = tract*ratios(1)

        ! -  orthotropic activation
        ! REF: Rossi et al. 2014
        tract_s = ratios(2)*tract_f
        tract_n = ratios(3)*tract_f

        ! Outer products 
        fxf = OUT_VXV(fib(:),3_ip,fib(:),3_ip)
        sxs = OUT_VXV(sht(:),3_ip,sht(:),3_ip)
        nxn = OUT_VXV(nrm(:),3_ip,nrm(:),3_ip)

        ! Active Cauchy stress tensor
        sig(:,:) = (1.0_rp/det_J)*(tract_f*fxf(:,:) + tract_s*sxs(:,:) + tract_n*nxn(:,:))

        ! Add isotropic contribution to the active Cauchy stress tensor 
        ! (only transversally isotropic case)
        sig(1,1) = sig(1,1) + tract_i
        sig(2,2) = sig(2,2) + tract_i
        sig(3,3) = sig(3,3) + tract_i

        ! Active CAUCHY 
        dsig(:,:,:,:) = 0.0_rp
        ! -------------------------------
    end subroutine get_active_cauchy_stress_tensor

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief   Decompose traction according to a basis
    !> @details 
    ! -------------------------------------------------------------------------
    subroutine get_active_second_elasticity_tensor( &
        ratios, fib, sht, nrm, lamda, tract, dtract, dtang )
        ! -------------------------------
        implicit none
        ! -------------------------------
        real(rp),    intent(in)        :: ratios(3)       ! Transversal ratio
        real(rp),    intent(in)        :: fib(3)          ! Basis {f, s, n}
        real(rp),    intent(in)        :: sht(3)
        real(rp),    intent(in)        :: nrm(3)     
        real(rp),    intent(in)        :: lamda(3)        ! Stretches
        real(rp),    intent(in)        :: tract           ! Active traction along the fiber direction
        real(rp),    intent(in)        :: dtract          ! Derivative of the active traction with respect to the stretch along the fiber direciton
        real(rp),    intent(inout)     :: dtang(:,:,:,:)  ! Active Cauchy stress tensor
        ! -------------------------------
        real(rp),   dimension(3,3,3,3) :: fxfxfxf, sxsxsxs, nxnxnxn, sxsxfxf, nxnxfxf
        ! -------------------------------
        ! Outer products 
        fxfxfxf = OUT_VxVxVxV(fib(:),3_ip,fib(:),3_ip,fib(:),3_ip,fib(:),3_ip)
        sxsxsxs = OUT_VxVxVxV(sht(:),3_ip,sht(:),3_ip,sht(:),3_ip,sht(:),3_ip)
        nxnxnxn = OUT_VxVxVxV(nrm(:),3_ip,nrm(:),3_ip,nrm(:),3_ip,nrm(:),3_ip)
        sxsxfxf = OUT_VxVxVxV(sht(:),3_ip,sht(:),3_ip,fib(:),3_ip,fib(:),3_ip)
        nxnxfxf = OUT_VxVxVxV(nrm(:),3_ip,nrm(:),3_ip,fib(:),3_ip,fib(:),3_ip)
            
        ! Fiber direction contribution
        dtang = (1.0_rp/(lamda(1)**3))*(dtract - 2.0_rp*(tract/lamda(1)))*fxfxfxf

        ! Sheet direction activation
        if( lamda(2) > 0.0_rp )then
            dtang = dtang + ratios(2)*((dtract/((lamda(2)**2)*lamda(1)))*sxsxfxf - &
                (2.0_rp/(lamda(2)**4))*tract*sxsxsxs) 
        endif
 
        ! Norma direction activation
        if( lamda(3) > 0.0_rp )then
            dtang = dtang + ratios(3)*((dtract/((lamda(3)**2)*lamda(1)))*nxnxfxf - &
                (2.0_rp/(lamda(3)**4))*tract*nxnxnxn) 
        endif
        ! -------------------------------
    end subroutine get_active_second_elasticity_tensor

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @author  Pierre Lafortune 
    !> @date    
    !> @brief   Calculate the active stress tensor
    !> @details Ref_A : P.J. Hunter  et al. 1998, DOI: 10.1016/S0079-6107(98)00013-3
    !           Ref_B : P. Lafortune et al. 2012, DOI: 10.1002/cnm.1494
    ! -------------------------------------------------------------------------
    subroutine get_active_traction_HUNTER( &
        id_cell, props, lamda, Ca, tract, dtract, sdv )
        use def_master,           only :  dtime
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: id_cell   ! ID cellular model
        real(rp),    intent(in)        :: props(:)  ! Properties
        real(rp),    intent(in)        :: lamda     ! Lambda (stretches)
        real(rp),    intent(in)        :: Ca        ! Calcium
        real(rp),    intent(inout)     :: tract     ! Active traction
        real(rp),    intent(inout)     :: dtract    ! Derivative active traction w.r.t lamda
        real(rp),    intent(inout)     :: sdv(:,:)  ! State variables
        ! -------------------------------
        real(rp)                       :: Ca0, Cam, Ca50, Ca2p
        real(rp)                       :: thau, eta, n, k
        real(rp)                       :: tmax
        ! -------------------------------
        ! Assign properties
        Cam  = props(1)
        Ca50 = props(2)
        thau = props(3)
        eta  = props(4)
        n    = props(5)
        k    = props(6)
        tmax = props(7)

        ! Compute Ca2p according to the cellular model
        if(     id_cell == EXMSLD_CELL_TT2006   .or.  &
                id_cell == EXMSLD_CELL_OHARA    .or.  &
                id_cell == EXMSLD_CELL_NO_CELL  )then 
                
            Ca2p = (Ca - SDV(9,1))

        else 

            call runend ("SLD-EXM-ECC : CELL MODEL NOT INCLUDED IN THE COUPLING WITH SOLIDZ MODULE")

        end if

        ! Compute active traction
        ! { Ref_B : Eq. 21}
        tract = (Ca2p**n)*tmax*(1.0_rp + eta*(lamda - 1.0_rp))/((Ca2p**n) + (Ca50**n))

        ! Compute derivative active traction w.r.t. lambda along the fiber direction
        dtract = (Ca2p**n)*tmax*eta/((Ca2p**n) + (Ca50)**n)
        ! -------------------------------
    end subroutine get_active_traction_HUNTER

    ! -------------------------------------------------------------------------
    !> @author  Francesc Levrero-Florencio
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief   Calculate the active stress tensor
    !> @details 
    ! -------------------------------------------------------------------------
    subroutine get_active_traction_LANDS( &
        props, lamda, Ca, tract, dtract, sdv )
        use def_master,           only :  dtime
        ! -------------------------------
        implicit none
        ! -------------------------------
        real(rp),    intent(in)        :: props(:)  ! Properties
        real(rp),    intent(in)        :: lamda     ! Lambda (stretches)
        real(rp),    intent(in)        :: Ca        ! Calcium
        real(rp),    intent(inout)     :: tract     ! Active traction
        real(rp),    intent(inout)     :: dtract    ! Derivative active traction w.r.t lamda
        real(rp),    intent(inout)     :: sdv(:,:)  ! State variables
        ! -------------------------------
        integer(ip)                    :: i, j
        real(rp)                       :: k_TRPN, n_TRPN, TRPN_50, Ca50_ref
        real(rp)                       :: n_tm
        real(rp)                       :: k_u, k_u_2, k_uw, k_ws, k_wu, k_su
        real(rp)                       :: r_s, r_w
        real(rp)                       :: y_s, y_w, y_su, y_wu
        real(rp)                       :: zeta_s, zeta_w
        real(rp)                       :: c_s, c_w
        real(rp)                       :: phi
        real(rp)                       :: A_eff, A_s
        real(rp)                       :: beta_0, beta_1
        real(rp)                       :: T_ref
        real(rp)                       :: lamda_curr, lamda_auxi, lamda_prev, lamda_rate
        real(rp)                       :: B, S, W, h, Ca50, CaTRPN, term
        integer(ip)                    :: ii
        integer(ip), parameter         :: old = 1_ip, new = 2_ip
        real(rp),    parameter         :: lamda_thr = 1.2_rp
        integer(ip), parameter         :: NR_ite_max = 100_ip
        real(rp),    parameter         :: NR_res_tol = 1.0e-10_rp
        integer(ip)                    :: NR_ite
        real(rp)                       :: NR_res_nrm
        real(rp),    dimension(6,2)    :: NR_sdv
        real(rp),    dimension(6)      :: NR_sol, NR_res, ODE_ydot
        real(rp),    dimension(6,6)    :: NR_jac
        real(rp),    dimension(6)      :: dRdlmd, dSdlmd
        real(rp)                       :: dhdlmd
        ! -------------------------------
        ! Assign properties
        k_TRPN   = props( 1)
        n_TRPN   = props( 2)
        Ca50_ref = props( 3)
        k_u      = props( 4)
        n_tm     = props( 5)
        TRPN_50  = props( 6)
        k_uw     = props( 7)
        k_ws     = props( 8)
        r_w      = props( 9)
        r_s      = props(10)
        y_s      = props(11)
        y_w      = props(12)
        phi      = props(13)
        A_eff    = props(14)
        beta_0   = props(15)
        beta_1   = props(16)
        T_ref    = props(17)

        ! Compute lambda values
        lamda_curr = lamda
        lamda_auxi = min(lamda_thr,lamda_curr)
        lamda_prev = sdv(8,old)                       ! Previous lambda from state variable
        lamda_rate = (lamda_curr - lamda_prev)/dtime  ! Lambda rate throught finite difference

        ! Compute Ca_50
        Ca50 = Ca50_ref + beta_1*(lamda_auxi - 1.0_rp)

        ! Compute h
        h = max(0.0_rp, 1.0_rp + beta_0*(lamda_auxi + min(0.87_rp, lamda_auxi) - 1.87_rp))

        ! Compute dependent parameters 
        k_wu = k_uw*((1.0_rp/r_w)-1.0_rp)-k_ws
        k_su = k_ws*((1.0_rp/r_s)-1.0_rp)*r_w
        A_s = (A_eff*r_s)/((1.0_rp-r_s)*r_w+r_s)
        c_s = phi*k_ws*(1.0_rp-r_s)*r_w/r_s
        c_w = phi*k_uw*((1.0_rp-r_s)*(1.0_rp-r_w))/((1.0_rp-r_s)*r_w)
        k_u_2 = (k_u*(TRPN_50**n_tm))/(1.0_rp-r_s-(1.0_rp-r_s)*r_w)

        ! Loop for the Newton-Raphson of the Backward Euler integration scheme
        ! - initialization
        NR_sdv(1:6,old) = sdv(1:6,old)
        NR_sdv(1:6,new) = sdv(1:6,old)
        S      = NR_sdv(1,new)
        W      = NR_sdv(2,new)
        CaTRPN = NR_sdv(3,new)
        B      = NR_sdv(4,new)
        zeta_s = NR_sdv(5,new)
        zeta_w = NR_sdv(6,new)
        ! - loop
        NR_loop: do NR_ite = 1, NR_ite_max

            ! Regularisation of some state variable 
            ! - compute y_su
            y_su = 0.0_rp
            if(     zeta_s + 1.0_rp < 0.0_rp )then
                y_su = -y_s*(zeta_s + 1.0_rp )
            elseif( zeta_s + 1.0_rp > 1.0_rp )then
                y_su = y_s*zeta_s
            endif
            ! - compute y_wu
            y_wu = y_w*abs(zeta_w)

            ! Compute the RHS of the ODE system
            ODE_ydot(1) = k_ws*W - k_su*S - y_su*S
            ODE_ydot(2) = k_uw*(1.0_rp - B - S - W) - k_wu*W - k_ws*W - y_wu*W
            ODE_ydot(3) = k_TRPN*(((Ca/Ca50)**n_TRPN)*(1.0_rp - CaTRPN) -CaTRPN)
            ODE_ydot(4) = k_u_2*min((CaTRPN**(-n_tm/2.0_rp)),100.0_rp)*(1.0_rp - B - S - W) - &
            &   k_u*(CaTRPN**(n_tm/2.0_rp))*B
            ODE_ydot(5) = A_s*lamda_rate - c_s*zeta_s
            ODE_ydot(6) = A_s*lamda_rate - c_w*zeta_w

            ! Compute Jacobian of the linearised system of equations
            NR_jac(:,:) = 0.0_rp
            ! jac(1,:)
            NR_jac(1,1) = -k_su - y_su
            NR_jac(1,2) =  k_ws
            if(     zeta_s + 1.0_rp < 0.0_rp )then
                NR_jac(1,5) = -y_s
            elseif( zeta_s + 1.0_rp > 1.0_rp )then
                NR_jac(1,5) = y_s
            endif
            ! jac(2,:)
            NR_jac(2,1) = -k_uw
            NR_jac(2,2) = -k_uw - k_wu -k_ws - y_wu
            NR_jac(2,4) = -k_uw
            if( zeta_w < 0.0_rp )then
                NR_jac(2,6) = -y_w
            else
                NR_jac(2,6) = y_w
            endif
            ! jac(3,:)
            NR_jac(3,3) = -k_TRPN*(((Ca/Ca50)**n_TRPN) + 1.0_rp)
            ! jac(4,:)
            if (CaTRPN**(-n_tm/2.0_rp) <= 100.0_rp) then
                term = CaTRPN**(-n_tm/2.0_rp)
            else
                term = 100.0_rp
            end if
            NR_jac(4,1) = -k_u_2*term
            NR_jac(4,2) = -k_u_2*term
            if( CaTRPN**(-n_tm/2.0_rp) <= 100.0_rp )then
                NR_jac(4,3) = (-n_tm/2.0_rp)*(k_u_2*(1.0_rp - B - S - W)*(CaTRPN**(-(n_tm/2.0_rp) - 1.0_rp)) + &
                &   k_u*B*(CaTRPN**((n_tm/2.0_rp)-1.0_rp)))
            else
                NR_jac(4,3) = (-n_tm/2.0_rp)*k_u*B*(CaTRPN**((n_tm/2.0_rp)-1.0_rp))
            end if
            NR_jac(4,4) = -k_u_2*term - k_u*(CaTRPN**(n_tm/2.0_rp))
            ! jac(5,:)
            NR_jac(5,5) = -c_s
            ! jac(6,:)
            NR_jac(6,6) = -c_w
            ! jac = I_6x6 - jac()
            NR_jac(:,:) = -NR_jac(:,:)*dtime
            do ii = 1, 6
                NR_jac(ii,ii) = 1.0_rp + NR_jac(ii,ii)
            enddo

            ! Compute residual of the linearised system of equations 
            ! (previous solution are the state variables)
            NR_res(:) = NR_sdv(:,new) - NR_sdv(:,old) - dtime*ODE_ydot(:)

            ! Solve the system of the equations
            call invert(NR_jac,6_ip,6_ip)
            NR_sdv(:,new) = NR_sdv(:,new) - matmul(NR_jac(:,:),NR_res(:))
            S      = NR_sdv(1,new)
            W      = NR_sdv(2,new)
            CaTRPN = NR_sdv(3,new)
            B      = NR_sdv(4,new)
            zeta_s = NR_sdv(5,new)
            zeta_w = NR_sdv(6,new)

            ! Check the convergence
            ! - compute norm of the residual of the linearised system pf equations
            NR_res_nrm = 0.0_rp
            do ii = 1, 6
                NR_res_nrm = NR_res_nrm + NR_res(ii)**2
            enddo
            NR_res_nrm = sqrt(NR_res_nrm)
            ! - apply convergence criterion
            if( NR_res_nrm < NR_res_tol )then
                exit
            endif

        enddo NR_loop

        ! Compute the active tractions
        ! - S      = NR_sdv(1,2)
        ! - W      = NR_sdv(2,2)
        ! - CaTRN  = NR_sdv(3,2)
        ! - B      = NR_sdv(4,2)
        ! - zeta_s = NR_sdv(5,2)
        ! - zeta_w = NR_sdv(6,2)
        tract = h*(T_ref/r_s)*(NR_sdv(1,2)*(NR_sdv(5,2)+1.0_rp)+NR_sdv(2,2)*NR_sdv(6,2))
        ! - convert from kPa to GS units 
        tract = tract*10000.0_rp

        ! Compute the derivative of the active tractions w.r.t.  lambda
        ! - dR/dlamda
        dRdlmd(:) = 0.0_rp
        if( lamda_curr < lamda_thr )then
            dRdlmd(3) = beta_1*n_TRPN*k_TRPN*dtime*(1.0_rp - CaTRPN)* &
            &   ((Ca**n_TRPN)*(Ca50**(-n_TRPN - 1.0_rp)))
        endif
        dRdlmd(5) = - A_s
        dRdlmd(6) = - A_s
        ! - dSdlamda
        ! NR_jac = NR_jac_inv (see above)
        dSdlmd(:) = -matmul(NR_jac,dRdlmd)
        ! - dhdlambda
        dhdlmd = 0.0_rp
        if(     lamda_curr < ((1.87_rp*beta_0 - 1.0_rp)/(2.0_rp*beta_0)) )then
            dhdlmd = 0.0_rp
        elseif( lamda_curr < 0.87_rp) then
            dhdlmd = 2.0_rp*beta_0
        elseif( lamda_curr < 1.2_rp )then
            dhdlmd = beta_0
        end if
        ! - DtractDlmd
        dtract = (T_ref*10000_rp/r_s)*(dhdlmd*((zeta_s + 1.0_rp)*S + zeta_w*W) + &
        &   h*((dSdlmd(5)*S + zeta_s*dSdlmd(1)) + (dSdlmd(6)*W + zeta_w*dSdlmd(2))))

        ! Assign new state variables
        sdv(1:6,new) = NR_sdv(1:6,new) ! Computed through the Backeuler
        sdv(  8,new) = sdv(8,old)      !  We didnt modify lamda previus
        ! -------------------------------
    end subroutine get_active_traction_LANDS

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief   
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    pure function ROUNDOFF(val,tgt) result(res)
        ! -------------------------------
        implicit none
        ! -------------------------------
        real(rp), intent(in)           :: val
        real(rp), intent(in)           :: tgt
        real(rp)                       :: res
        real(rp), parameter            :: thr = 1.0e-12_rp
        ! -------------------------------
        if( abs((val - tgt)) < thr )then
            res = tgt
        else
            res = val
        endif
        ! -------------------------------
    end function ROUNDOFF

    pure function MULT_MxV(M,mi,mj,v,vj) result(r)
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)                :: mi, mj, vj
        real(rp), dimension(mi,mj), intent(in) :: M
        real(rp), dimension(vj),    intent(in) :: v
        real(rp), dimension(mi)                :: r
        integer(ip)                            :: i, j
        ! -------------------------------
        r(:) = 0.0_rp
        do j = 1, mi
            do i = 1, mj
                r(i) =  r(i) + M(i,j)*v(j)
            enddo
        enddo
        ! -------------------------------
    end function MULT_MxV

    pure function MULT_MxT(M,mi,mj,T,ti,tj) result(R)
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)                :: mi, mj, ti, tj
        real(rp), dimension(mi,mj), intent(in) :: M
        real(rp), dimension(ti,tj), intent(in) :: T
        real(rp), dimension(mi,ti)             :: R
        integer(ip)                            :: i, j, k
        ! -------------------------------
        R(:,:) = 0.0_rp
        do i = 1, mi
            do j = 1, ti
                do k = 1, mj
                    R(i,j) = R(i,j) + M(i,k)*T(j,k)
                end do
            end do
        end do
        ! -------------------------------
    end function MULT_MxT

    pure function MULT_VxMxV(a,ai,M,b,bj) result(r)
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)                :: ai, bj
        real(rp), dimension(ai),    intent(in) :: a
        real(rp), dimension(ai,bj), intent(in) :: M
        real(rp), dimension(bj),    intent(in) :: b
        real(rp)                               :: r
        integer(ip)                            :: i, j
        ! -------------------------------
        r = 0.0_rp
        do j = 1, ai
            do i = 1, bj
                r = r + a(i)*M(i,j)*b(j)
            enddo
        enddo
        ! -------------------------------
    end function MULT_VxMxV

    pure function OUT_VxV(a,ai,b,bi) result(R)
        ! EQUIV to dot_product(a,b)
        !   where A and B are matrices
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)             :: ai, bi
        real(rp), dimension(ai), intent(in) :: a
        real(rp), dimension(bi), intent(in) :: b
        real(rp), dimension(ai,bi)          :: R
        integer(ip)                         :: i, j
        ! -------------------------------
        do i = 1, ai
            do j = 1, bi
                R(i,j) = a(i)*b(j)
            enddo
        enddo
        ! -------------------------------
    end function OUT_VxV

    pure function OUT_VxVxVxV(a,ai,b,bi,c,ci,d,di) result(R)
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)             :: ai, bi, ci, di
        real(rp), dimension(ai), intent(in) :: a
        real(rp), dimension(bi), intent(in) :: b
        real(rp), dimension(bi), intent(in) :: c
        real(rp), dimension(bi), intent(in) :: d
        real(rp), dimension(ai,bi,ci,di)    :: R
        integer(ip)                         :: i, j, k, l
        ! -------------------------------
        do i = 1, ai
            do j = 1, bi
                do k = 1, ci
                    do l = 1, di
                        R(i,j,k,l) = a(i)*b(j)*c(k)*d(l)
                    enddo
                enddo
            enddo
        enddo
        ! -------------------------------
    end function OUT_VxVxVxV

    pure function INV3x3(M) result(Mi)
        ! -------------------------------
        implicit none
        ! -------------------------------
        real(rp), dimension(3,3), intent(in) :: M
        real(rp), dimension(3,3)             :: Mi
        real(rp)                             :: det
        ! -------------------------------
        ! Determinant
        det = M(1,1)*M(2,2)*M(3,3) + M(1,3)*M(2,1)*M(3,2) + &
            M(3,1)*M(1,2)*M(2,3) - M(3,1)*M(2,2)*M(1,3) - &
            M(3,3)*M(1,2)*M(2,1) - M(1,1)*M(2,3)*M(3,2)
        ! Invert matrix
        Mi(1,1) =  (M(2,2)*M(3,3) - M(3,2)*M(2,3))/det
        Mi(1,2) = -(M(1,2)*M(3,3) - M(1,3)*M(3,2))/det
        Mi(1,3) =  (M(1,2)*M(2,3) - M(2,2)*M(1,3))/det
        Mi(2,1) = -(M(2,1)*M(3,3) - M(3,1)*M(2,3))/det
        Mi(2,2) =  (M(1,1)*M(3,3) - M(1,3)*M(3,1))/det
        Mi(2,3) = -(M(1,1)*M(2,3) - M(2,1)*M(1,3))/det
        Mi(3,1) =  (M(2,1)*M(3,2) - M(3,1)*M(2,2))/det
        Mi(3,2) = -(M(1,1)*M(3,2) - M(3,1)*M(1,2))/det
        Mi(3,3) =  (M(1,1)*M(2,2) - M(1,2)*M(2,1))/det
        ! -------------------------------
    end function INV3x3
    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief   
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    subroutine set_has_land()
      use def_master,           only :  kfl_eccty
      implicit none
        integer(ip)   :: imate
        has_land=.False.
        do imate = 1,nmate
         if( kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR ) then
             has_land=.True.
         endif
        end do
    endsubroutine set_has_land
    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @date    
    !> @brief   
    !> @details 
    !> @todo
    ! -------------------------------------------------------------------------
    function has_exmsld_coupling() result(it_has)
        use def_master,  only : coupling
        use def_master,  only :  kfl_eccty
        implicit none
        integer(ip)   :: imate
        logical(lg)   :: it_has
        it_has=.False.
        if((coupling('SOLIDZ','EXMEDI') >= 1_ip  .or. coupling('EXMEDI','SOLIDZ') >= 1_ip) ) it_has=.True.

    end function has_exmsld_coupling
 
    ! =========================================================================
end module mod_exm_sld_eccoupling
