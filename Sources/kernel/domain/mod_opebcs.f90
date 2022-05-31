module mod_opebcs
  !-----------------------------------------------------------------------
  !****f* outrut/mod_opebcs
  ! NAME
  !   mod_opebcs
  ! DESCRIPTION
  !   This routine manages the opebcsocess
  ! USES
  ! USED BY
  !   output_***
  !   outvar_***
  !***
  !-----------------------------------------------------------------------

  use def_master,         only : rp
  use def_master,         only : ip,lg
  use def_master,         only : parii
  use def_master,         only : npari
  use def_master,         only : nparc
  use def_master,         only : nparr
  use def_master,         only : nparx
  use def_master,         only : nparh
  use def_master,         only : parin
  use def_master,         only : parre
  use def_master,         only : parhh
  use def_master,         only : IMASTER
  use def_master,         only : ISLAVE
  use def_master,         only : IPARALL
  use def_master,         only : bc_nodes
  use def_master,         only : bc_bound
  use def_master,         only : coutp
  use def_master,         only : mem_modul
  use def_master,         only : modul
  use def_domain,         only : mcodb
  use def_domain,         only : mcono
  use def_domain,         only : kfl_icodb
  use def_domain,         only : kfl_icodn
  use def_domain,         only : kfl_geome
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_outfor,         only : outfor
  use mod_communications, only : PAR_EXCHANGE
  use mod_communications, only : PAR_BROADCAST
  implicit none

  private

  integer(ip) :: mcodc

  interface opebcs_initialization_structure
     module procedure     opebcs_initialization_structure_boundary,&
          &               opebcs_initialization_structure_node
  end interface opebcs_initialization_structure

  interface opebcs_initialization_variable
     module procedure     opebcs_initialization_variable_boundary,&
          &               opebcs_initialization_variable_boundary_s,&
          &               opebcs_initialization_variable_node_s,&
          &               opebcs_initialization_variable_node_1
  end interface opebcs_initialization_variable

  public :: boundary_conditions_read_node_codes  
  public :: boundary_conditions_read_boundary_codes   
  public :: boundary_conditions_impose_node_codes   ! Impose conditions on nodes 
  public :: opebcs_initialization_structure         ! Allocate bc structure
  public :: opebcs_initialization_variable          ! Initialize bc structure
  public :: opbbcs
  public :: opnbcs
  public :: spbbcs
  public :: spnbcs
  public :: spgbcs
  public :: cpybcs
  public :: cpybcs_boundaries

contains

  subroutine opebcs_initialization()

    mcodc = 3_ip * mcodb

  end subroutine

  subroutine opnbcs(itask,nvari,ndofn,kfl_ibopo,tncod_xxx)

    !--------------------------------------------------------------------
    !
    ! NODE
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),    intent(in)             :: itask,nvari,ndofn,kfl_ibopo
    integer(ip)                            :: icode,icono,idofn,ivari
    type(bc_nodes), intent(inout), pointer :: tncod_xxx(:)

    call opebcs_initialization()

    select case ( itask )

    case ( 0_ip )
       !
       ! Allocate type for NVARI variables
       !
       allocate( tncod_xxx(nvari) )
       do ivari = 1,nvari
          allocate( tncod_xxx(ivari) % l(mcodc) )
          tncod_xxx(ivari) % kfl_ibopo = kfl_ibopo
          tncod_xxx(ivari) % ndofn = ndofn  
          tncod_xxx(ivari) % ncode = 0
          do icode = 1,mcodc 
             allocate( tncod_xxx(ivari) % l(icode) % lcode(mcono) )
             allocate( tncod_xxx(ivari) % l(icode) % bvess(ndofn) )
             tncod_xxx(ivari) % l(icode) % kfl_fixno  = 0
             tncod_xxx(ivari) % l(icode) % kfl_value  = 0
             tncod_xxx(ivari) % l(icode) % kfl_funno  = 0
             tncod_xxx(ivari) % l(icode) % kfl_funtyp = 0
             tncod_xxx(ivari) % l(icode) % kfl_fixrs  = 0
             tncod_xxx(ivari) % l(icode) % tag        = ''
             tncod_xxx(ivari) % l(icode) % fname      = ''
             do icono = 1,mcono
                tncod_xxx(ivari) % l(icode) % lcode(icono) =  mcodb+1
             end do
             do idofn = 1,ndofn
                tncod_xxx(ivari) % l(icode) % bvess(idofn) =  0.0_rp
             end do
          end do
       end do

    case ( 1_ip )
       !
       ! Allocate type for NVARI variables
       !
       allocate( tncod_xxx(nvari) )

    case ( 2_ip )
       !
       ! Allocate memory for each code ICODE
       !
       allocate( tncod_xxx(nvari) % l(mcodc) )
       tncod_xxx(nvari) % kfl_ibopo = kfl_ibopo
       tncod_xxx(nvari) % ndofn = ndofn     
       tncod_xxx(nvari) % ncode = 0     
       do icode = 1,mcodc
          allocate( tncod_xxx(nvari) % l(icode) % lcode(mcono) )
          allocate( tncod_xxx(nvari) % l(icode) % bvess(ndofn) )
          tncod_xxx(nvari) % l(icode) % kfl_fixno  = 0
          tncod_xxx(nvari) % l(icode) % kfl_value  = 0
          tncod_xxx(nvari) % l(icode) % kfl_funno  = 0
          tncod_xxx(nvari) % l(icode) % kfl_funtyp = 0
          tncod_xxx(nvari) % l(icode) % kfl_fixrs  = 0
          tncod_xxx(nvari) % l(icode) % tag        = ''
          tncod_xxx(nvari) % l(icode) % fname      = ''
          do icono = 1,mcono
             tncod_xxx(nvari) % l(icode) % lcode(icono) =  mcodb+1
          end do
          do idofn = 1,ndofn
             tncod_xxx(nvari) % l(icode) % bvess(idofn) =  0.0_rp
          end do
       end do

    case ( 3_ip )
       !
       ! Deallocate whole type
       !
       do ivari = 1,size(tncod_xxx,KIND=ip)
          do icode = 1,mcodc
             deallocate( tncod_xxx(ivari) % l(icode) % lcode )
             deallocate( tncod_xxx(ivari) % l(icode) % bvess )
          end do
          deallocate( tncod_xxx(ivari) % l )
       end do
       deallocate( tncod_xxx )

    end select

  end subroutine opnbcs

  !--------------------------------------------------------------------
  !
  ! Initialization
  !
  !--------------------------------------------------------------------

  subroutine opebcs_initialization_structure_boundary(nvari,tbcod_xxx)
    implicit none
    integer(ip),             intent(in)    :: nvari
    integer(4)                             :: istat
    type(bc_bound), pointer, intent(inout) :: tbcod_xxx(:)

    if( associated(tbcod_xxx) ) call runend('STRUCTURE TBCOD ALREADY ASSOCIATED')
    allocate( tbcod_xxx(nvari) , stat = istat )
    if( istat /= 0 ) call runend('COULD NOT ALLOCATE STRUCTURE')

  end subroutine opebcs_initialization_structure_boundary

  subroutine opebcs_initialization_structure_node(nvari,tncod_xxx)
    implicit none
    integer(ip),             intent(in)    :: nvari
    integer(4)                             :: istat
    type(bc_nodes), pointer, intent(inout) :: tncod_xxx(:)

    if( associated(tncod_xxx) ) call runend('STRUCTURE TNCOD ALREADY ASSOCIATED')
    allocate( tncod_xxx(nvari) , stat = istat )
    if( istat /= 0 ) call runend('COULD NOT ALLOCATE STRUCTURE')

  end subroutine opebcs_initialization_structure_node

  subroutine opebcs_initialization_variable_boundary(ndofn,tbcod_xxx)
    implicit none
    integer(ip),             intent(in)    :: ndofn
    type(bc_bound), pointer, intent(inout) :: tbcod_xxx(:)
    integer(ip)                            :: icode,idofn,ivari


    call opebcs_initialization()

    if( .not. associated(tbcod_xxx) ) call runend('opebcs_initialization_variable_boundary: STRUCTURE TBCOD NOT ASSOCIATED')

    do ivari = 1,size(tbcod_xxx,KIND=ip)
       allocate( tbcod_xxx(ivari) % l(mcodc) )
       tbcod_xxx(ivari) % ndofn = ndofn     
       tbcod_xxx(ivari) % ncode = 0
       do icode = 1,mcodc
          allocate( tbcod_xxx(ivari) % l(icode) % bvnat(ndofn) )
          tbcod_xxx(ivari) % l(icode) % kfl_fixbo  =  0
          tbcod_xxx(ivari) % l(icode) % kfl_value  =  0
          tbcod_xxx(ivari) % l(icode) % lcode      =  mcodb+1
          tbcod_xxx(ivari) % l(icode) % kfl_funtyp =  0
          tbcod_xxx(ivari) % l(icode) % kfl_funbo  =  0
          do idofn = 1,ndofn
             tbcod_xxx(ivari) % l(icode) % bvnat(idofn) =  0.0_rp
          end do
       end do
    end do

  end subroutine opebcs_initialization_variable_boundary

  subroutine opebcs_initialization_variable_boundary_s(ndofn,tbcod_xxx)
    implicit none
    integer(ip),    intent(in)    :: ndofn
    type(bc_bound), intent(inout) :: tbcod_xxx
    integer(ip)                   :: icode,idofn


    call opebcs_initialization()

    allocate( tbcod_xxx % l(mcodc) )
    tbcod_xxx % ndofn = ndofn     
    tbcod_xxx % ncode = 0
    do icode = 1,mcodc
       allocate( tbcod_xxx % l(icode) % bvnat(ndofn) )
       tbcod_xxx % l(icode) % kfl_fixbo =  0
       tbcod_xxx % l(icode) % kfl_value =  0
       tbcod_xxx % l(icode) % lcode     =  mcodb+1
       tbcod_xxx % l(icode) % kfl_funtyp =  0
       tbcod_xxx % l(icode) % kfl_funbo =  0
       do idofn = 1,ndofn
          tbcod_xxx % l(icode) % bvnat(idofn) =  0.0_rp
       end do
    end do

  end subroutine opebcs_initialization_variable_boundary_s

  subroutine opebcs_initialization_variable_node_1(ndofn,tncod_xxx,kfl_ibopo)
    implicit none
    integer(ip),              intent(in)    :: ndofn
    type(bc_nodes), pointer,  intent(inout) :: tncod_xxx(:)
    integer(ip),    optional, intent(in)    :: kfl_ibopo
    integer(ip)                             :: kbopo
    integer(ip)                             :: icode,icono,idofn,ivari

    call opebcs_initialization()

    if( present(kfl_ibopo) ) then
       kbopo = kfl_ibopo
    else
       kbopo = 0
    end if

    if( .not. associated(tncod_xxx) ) call runend('opebcs_initialization_variable_node: STRUCTURE TNCOD NOT ASSOCIATED')
    do ivari = 1,size(tncod_xxx,KIND=ip)
       allocate( tncod_xxx(ivari) % l(mcodc) )
       tncod_xxx(ivari) % kfl_ibopo = kbopo
       tncod_xxx(ivari) % ndofn = ndofn     
       tncod_xxx(ivari) % ncode = 0     
       do icode = 1,mcodc
          allocate( tncod_xxx(ivari) % l(icode) % lcode(mcono) )
          allocate( tncod_xxx(ivari) % l(icode) % bvess(ndofn) )
          tncod_xxx(ivari) % l(icode) % kfl_fixno  = 0
          tncod_xxx(ivari) % l(icode) % kfl_value  = 0
          tncod_xxx(ivari) % l(icode) % kfl_funno  = 0
          tncod_xxx(ivari) % l(icode) % kfl_funtyp = 0
          tncod_xxx(ivari) % l(icode) % kfl_fixrs  = 0
          tncod_xxx(ivari) % l(icode) % tag        = ''
          tncod_xxx(ivari) % l(icode) % fname      = ''
          do icono = 1,mcono
             tncod_xxx(ivari) % l(icode) % lcode(icono) =  mcodb+1
          end do
          do idofn = 1,ndofn
             tncod_xxx(ivari) % l(icode) % bvess(idofn) =  0.0_rp
          end do
       end do
    end do

  end subroutine opebcs_initialization_variable_node_1

  subroutine opebcs_initialization_variable_node_s(ndofn,tncod_xxx,kfl_ibopo)
    implicit none
    integer(ip),              intent(in)    :: ndofn
    type(bc_nodes),           intent(inout) :: tncod_xxx
    integer(ip),    optional, intent(in)    :: kfl_ibopo
    integer(ip)                             :: kbopo,icode,idofn,icono



    call opebcs_initialization()

    if( present(kfl_ibopo) ) then
       kbopo = kfl_ibopo
    else
       kbopo = 0
    end if

    allocate( tncod_xxx % l(mcodc) )
    tncod_xxx % kfl_ibopo = kbopo
    tncod_xxx % ndofn = ndofn     
    tncod_xxx % ncode = 0     
    do icode = 1,mcodc
       allocate( tncod_xxx % l(icode) % lcode(mcono) )
       allocate( tncod_xxx % l(icode) % bvess(ndofn) )
       tncod_xxx % l(icode) % kfl_fixno  = 0
       tncod_xxx % l(icode) % kfl_value  = 0
       tncod_xxx % l(icode) % kfl_funno  = 0
       tncod_xxx % l(icode) % kfl_funtyp = 0
       tncod_xxx % l(icode) % kfl_fixrs  = 0
       tncod_xxx % l(icode) % tag        = ''
       tncod_xxx % l(icode) % fname      = ''
       do icono = 1,mcono
          tncod_xxx % l(icode) % lcode(icono) =  mcodb+1
       end do
       do idofn = 1,ndofn
          tncod_xxx % l(icode) % bvess(idofn) =  0.0_rp
       end do
    end do

  end subroutine opebcs_initialization_variable_node_s

  subroutine opbbcs(itask,nvari,ndofn,tbcod_xxx)

    !--------------------------------------------------------------------
    !
    ! BOUNDARY
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),             intent(in)    :: itask,nvari,ndofn
    type(bc_bound), pointer, intent(inout) :: tbcod_xxx(:)
    integer(ip)                            :: icode,ivari,idofn



    call opebcs_initialization()

    select case ( itask )

    case ( 0_ip )
       !
       ! Allocate type for NVARI variables
       !
       allocate( tbcod_xxx(nvari) )
       do ivari = 1,nvari
          allocate( tbcod_xxx(ivari) % l(mcodc) )
          tbcod_xxx(ivari) % ndofn = ndofn  
          tbcod_xxx(ivari) % ncode = 0
          do icode = 1,mcodc
             allocate( tbcod_xxx(ivari) % l(icode) % bvnat(ndofn) )
             tbcod_xxx(ivari) % l(icode) % kfl_fixbo  = 0
             tbcod_xxx(ivari) % l(icode) % kfl_value  = 0
             tbcod_xxx(ivari) % l(icode) % kfl_funtyp = 0
             tbcod_xxx(ivari) % l(icode) % kfl_funbo  = 0
             do idofn = 1,ndofn
                tbcod_xxx(ivari) % l(icode) % bvnat(idofn) =  0.0_rp
             end do
          end do
       end do

    case ( 1_ip )
       !
       ! Allocate type for NVARI variables
       !
       allocate( tbcod_xxx(nvari) )

    case ( 2_ip )
       !
       ! Allocate memory for each code ICODE
       !
       allocate( tbcod_xxx(nvari) % l(mcodc) )
       tbcod_xxx(nvari) % ndofn = ndofn     
       tbcod_xxx(nvari) % ncode = 0
       do icode = 1,mcodc
          allocate( tbcod_xxx(nvari) % l(icode) % bvnat(ndofn) )
          tbcod_xxx(nvari) % l(icode) % kfl_fixbo  =  0
          tbcod_xxx(nvari) % l(icode) % kfl_value  =  0
          tbcod_xxx(nvari) % l(icode) % lcode      =  mcodb+1
          tbcod_xxx(nvari) % l(icode) % kfl_funtyp =  0
          tbcod_xxx(nvari) % l(icode) % kfl_funbo  =  0
          do idofn = 1,ndofn
             tbcod_xxx(nvari) % l(icode) % bvnat(idofn) =  0.0_rp
          end do
       end do

    case ( 3_ip )
       !
       ! Deallocate whole type
       !
       do ivari = 1,size(tbcod_xxx,KIND=ip)
          do icode = 1,mcodc
             deallocate( tbcod_xxx(ivari) % l(icode) % bvnat )
          end do
          deallocate( tbcod_xxx(ivari) % l )
       end do
       deallocate( tbcod_xxx )

    end select

  end subroutine opbbcs

  subroutine spnbcs(tncod_xxx,INCLUDE_CHARACTER)

    !--------------------------------------------------------------------
    !
    ! PARALL for TNCOD_XXX
    !
    !--------------------------------------------------------------------

    type(bc_nodes), pointer, intent(inout) :: tncod_xxx(:)
    logical(lg),    optional, intent(in)   :: INCLUDE_CHARACTER
    integer(ip)                            :: dummi,idofn,nvari
    integer(ip)                            :: ivari,icono,ndofn,icode
    logical(lg)                            :: if_include_character

    if_include_character = .true.
    if( present(INCLUDE_CHARACTER) ) if_include_character = INCLUDE_CHARACTER

    call opebcs_initialization()

    if( IPARALL .and. associated(tncod_xxx) ) then

       nvari = size(tncod_xxx,KIND=ip)

       do parii = 1,2 
          npari = 0
          nparr = 0
          nparh = 0
          
          do ivari = 1,nvari
             ndofn = tncod_xxx(ivari) % ndofn
             call PAR_EXCHANGE(tncod_xxx(ivari) % kfl_ibopo,parin,npari,parii)
             call PAR_EXCHANGE(tncod_xxx(ivari) % ndofn    ,parin,npari,parii)
             call PAR_EXCHANGE(tncod_xxx(ivari) % ncode    ,parin,npari,parii)
             do icode = 1,mcodc
                call PAR_EXCHANGE(tncod_xxx(ivari) % l(icode) % kfl_fixno, parin,npari,parii)
                call PAR_EXCHANGE(tncod_xxx(ivari) % l(icode) % kfl_value, parin,npari,parii)
                call PAR_EXCHANGE(tncod_xxx(ivari) % l(icode) % kfl_funno, parin,npari,parii)
                call PAR_EXCHANGE(tncod_xxx(ivari) % l(icode) % kfl_funtyp,parin,npari,parii)
                call PAR_EXCHANGE(tncod_xxx(ivari) % l(icode) % kfl_fixrs, parin,npari,parii)   
                if( if_include_character ) then
                   call PAR_EXCHANGE(5_ip,tncod_xxx(ivari) % l(icode) % tag  ,parhh,nparh,parii)
                   call PAR_EXCHANGE(5_ip,tncod_xxx(ivari) % l(icode) % fname,parhh,nparh,parii)
                end if
                do icono = 1,mcono
                   call PAR_EXCHANGE(tncod_xxx(ivari) % l(icode) % lcode(icono),parin,npari,parii)
                end do
                do idofn = 1,ndofn
                   call PAR_EXCHANGE(tncod_xxx(ivari) % l(icode) % bvess(idofn),parre,nparr,parii) 
                end do
             end do
          end do
          !
          ! Allocate memory for the first pass
          !
          if( parii == 1 ) then
             call memory_alloca(mem_modul(1:2,modul),'PARIN','ker_parall',parin,npari,'DO_NOT_INITIALIZE')
             call memory_alloca(mem_modul(1:2,modul),'PARRE','ker_parall',parre,nparr,'DO_NOT_INITIALIZE')
             call memory_alloca(mem_modul(1:2,modul),'PARHH','ker_parall',parhh,nparh,'DO_NOT_INITIALIZE')
          end if
          if( ( parii == 1 .and. ISLAVE ) .or. ( parii == 2 .and. IMASTER ) ) then
             call PAR_BROADCAST(parin,'IN MY CODE')
             call PAR_BROADCAST(parre,'IN MY CODE')
             call PAR_BROADCAST(parhh,'IN MY CODE')
          end if
          
       end do

       call memory_deallo(mem_modul(1:2,modul),'PARIN','ker_parall',parin)
       call memory_deallo(mem_modul(1:2,modul),'PARRE','ker_parall',parre)
       call memory_deallo(mem_modul(1:2,modul),'PARHH','ker_parall',parhh)
       !
       ! Master deallocates memory
       !
       if( IMASTER ) then
          call opnbcs(3_ip,dummi,dummi,dummi,tncod_xxx)
       end if

    end if

  end subroutine spnbcs

  subroutine spgbcs(tncod_xxx)

    !--------------------------------------------------------------------
    !
    ! PARALL for TGCOD_XXX
    !
    !--------------------------------------------------------------------

    type(bc_nodes), pointer, intent(inout) :: tncod_xxx(:)
    integer(ip)                            :: dummi,idofn,nvari
    integer(ip)                            :: ivari,icono,ndofn,icode


    call opebcs_initialization()

    if( kfl_geome /= 0 .and. IPARALL .and. associated(tncod_xxx) ) then

       nvari = size(tncod_xxx,KIND=ip)

       do parii = 1,2 
          npari = 0
          nparr = 0
          nparh = 0
          do ivari = 1,nvari
             ndofn = tncod_xxx(ivari) % ndofn
             call PAR_EXCHANGE(tncod_xxx(ivari) % kfl_ibopo,parin,npari,parii)
             call PAR_EXCHANGE(tncod_xxx(ivari) % ndofn    ,parin,npari,parii)
             call PAR_EXCHANGE(tncod_xxx(ivari) % ncode    ,parin,npari,parii)           
             do icode = 1,mcodc                  
                call PAR_EXCHANGE(     tncod_xxx(ivari) % l(icode) % kfl_fixno, parin,npari,parii)           
                call PAR_EXCHANGE(     tncod_xxx(ivari) % l(icode) % kfl_value, parin,npari,parii)           
                call PAR_EXCHANGE(     tncod_xxx(ivari) % l(icode) % kfl_funno, parin,npari,parii)           
                call PAR_EXCHANGE(     tncod_xxx(ivari) % l(icode) % kfl_funtyp,parin,npari,parii)           
                call PAR_EXCHANGE(     tncod_xxx(ivari) % l(icode) % kfl_fixrs, parin,npari,parii)           
                call PAR_EXCHANGE(5_ip,tncod_xxx(ivari) % l(icode) % tag      , parhh,nparh,parii)
                call PAR_EXCHANGE(5_ip,tncod_xxx(ivari) % l(icode) % fname    , parhh,nparh,parii)

                do icono = 1,mcono
                   call PAR_EXCHANGE(tncod_xxx(ivari) % l(icode) % lcode(icono),parin,npari,parii)    
                end do
                do idofn = 1,ndofn
                   call PAR_EXCHANGE(tncod_xxx(ivari) % l(icode) % bvess(idofn),parre,nparr,parii) 
                end do
             end do
          end do
          !
          ! Allocate memory for the first pass
          !
          if( parii == 1 ) then
             call memory_alloca(mem_modul(1:2,modul),'PARIN','ker_parall',parin,npari,'DO_NOT_INITIALIZE')
             call memory_alloca(mem_modul(1:2,modul),'PARRE','ker_parall',parre,nparr,'DO_NOT_INITIALIZE')
             call memory_alloca(mem_modul(1:2,modul),'PARHH','ker_parall',parhh,nparh,'DO_NOT_INITIALIZE')
          end if
          if( ( parii == 1 .and. ISLAVE ) .or. ( parii == 2 .and. IMASTER ) ) then
             call PAR_BROADCAST(parin,'IN MY CODE')
             call PAR_BROADCAST(parre,'IN MY CODE')
             call PAR_BROADCAST(parhh,'IN MY CODE')
          end if

       end do
       call memory_deallo(mem_modul(1:2,modul),'PARIN','ker_parall',parin)
       call memory_deallo(mem_modul(1:2,modul),'PARRE','ker_parall',parre)
       call memory_deallo(mem_modul(1:2,modul),'PARHH','ker_parall',parhh)
       !
       ! Master deallocates memory
       !
       if( IMASTER ) then
          call opnbcs(3_ip,dummi,dummi,dummi,tncod_xxx)
       end if

    end if

  end subroutine spgbcs

  subroutine spbbcs(tbcod_xxx)

    !--------------------------------------------------------------------
    !
    ! PARALL for TBCOD_XXX
    !
    !--------------------------------------------------------------------

    type(bc_bound), pointer, intent(inout) :: tbcod_xxx(:)
    integer(ip)                            :: dummi,idofn,nvari
    integer(ip)                            :: ivari,ndofn,icode

    call opebcs_initialization()

    if( kfl_icodb /= 0 .and. IPARALL .and. associated(tbcod_xxx) ) then

       nvari = size(tbcod_xxx,KIND=ip)

       do parii = 1,2 
          npari = 0
          nparr = 0
          nparh = 0

          do ivari = 1,nvari
             ndofn = tbcod_xxx(ivari) % ndofn
             call PAR_EXCHANGE(tbcod_xxx(ivari) % ndofn,parin,npari,parii)
             call PAR_EXCHANGE(tbcod_xxx(ivari) % ncode,parin,npari,parii)
             do icode = 1,mcodc
                call PAR_EXCHANGE(tbcod_xxx(ivari) % l(icode) % kfl_fixbo,parin,npari,parii)
                call PAR_EXCHANGE(tbcod_xxx(ivari) % l(icode) % kfl_value,parin,npari,parii)
                call PAR_EXCHANGE(tbcod_xxx(ivari) % l(icode) % kfl_funtyp,parin,npari,parii)
                call PAR_EXCHANGE(tbcod_xxx(ivari) % l(icode) % kfl_funbo,parin,npari,parii)
                call PAR_EXCHANGE(tbcod_xxx(ivari) % l(icode) % lcode    ,parin,npari,parii)
                do idofn = 1,ndofn
                   call PAR_EXCHANGE(tbcod_xxx(ivari) % l(icode) % bvnat(idofn),parre,nparr,parii) 
                end do
             end do
          end do
          !
          ! Allocate memory for the first pass
          !
          if( parii == 1 ) then
             call memory_alloca(mem_modul(1:2,modul),'PARIN','ker_parall',parin,npari,'DO_NOT_INITIALIZE')
             call memory_alloca(mem_modul(1:2,modul),'PARRE','ker_parall',parre,nparr,'DO_NOT_INITIALIZE')
          end if
          if( ( parii == 1 .and. ISLAVE ) .or. ( parii == 2 .and. IMASTER ) ) then
             call PAR_BROADCAST(parin,      'IN MY CODE')
             call PAR_BROADCAST(parre,      'IN MY CODE')
          end if
       end do
       call memory_deallo(mem_modul(1:2,modul),'PARIN','ker_parall',parin)
       call memory_deallo(mem_modul(1:2,modul),'PARRE','ker_parall',parre)
       call memory_deallo(mem_modul(1:2,modul),'PARHH','ker_parall',parhh)
       !
       ! Master deallocates memory
       !
       if( IMASTER ) then
          call opbbcs(3_ip,dummi,dummi,tbcod_xxx)
       end if

    end if

  end subroutine spbbcs

  subroutine cpybcs(ivari,ivaro,tncod_xxx)

    !--------------------------------------------------------------------
    !
    ! Copy TNCOD_XXX to TNCOD_YYY
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),    intent(in)    :: ivari,ivaro
    type(bc_nodes), intent(inout) :: tncod_xxx(:)
    integer(ip)                   :: icono,ndofn,icode,idofn


    call opebcs_initialization()

    if( kfl_icodn /= 0) then

       ndofn = tncod_xxx(ivari) % ndofn
       tncod_xxx(ivaro) % kfl_ibopo = tncod_xxx(ivari) % kfl_ibopo 
       tncod_xxx(ivaro) % ndofn     = tncod_xxx(ivari) % ndofn      
       tncod_xxx(ivaro) % ncode     = tncod_xxx(ivari) % ncode      
       do icode = 1,mcodc                  
          tncod_xxx(ivaro) % l(icode) % kfl_fixno  = tncod_xxx(ivari) % l(icode) % kfl_fixno
          tncod_xxx(ivaro) % l(icode) % kfl_value  = tncod_xxx(ivari) % l(icode) % kfl_value
          tncod_xxx(ivaro) % l(icode) % kfl_funno  = tncod_xxx(ivari) % l(icode) % kfl_funno
          tncod_xxx(ivaro) % l(icode) % kfl_funtyp = tncod_xxx(ivari) % l(icode) % kfl_funtyp
          tncod_xxx(ivaro) % l(icode) % kfl_fixrs  = tncod_xxx(ivari) % l(icode) % kfl_fixrs
          tncod_xxx(ivaro) % l(icode) % tag        = tncod_xxx(ivari) % l(icode) % tag
          tncod_xxx(ivaro) % l(icode) % fname      = tncod_xxx(ivari) % l(icode) % fname
          do icono = 1,mcono
             tncod_xxx(ivaro) % l(icode) % lcode(icono) = tncod_xxx(ivari) % l(icode) % lcode(icono)
          end do
          do idofn = 1,ndofn
             tncod_xxx(ivaro) % l(icode) % bvess(idofn) = tncod_xxx(ivari) % l(icode) % bvess(idofn)
          end do
       end do

    end if

  end subroutine cpybcs

  subroutine cpybcs_boundaries(ivari,ivaro,tbcod_xxx)

    !--------------------------------------------------------------------
    !
    ! Copy TNCOD_XXX to TNCOD_YYY
    !
    !--------------------------------------------------------------------

    implicit none
    integer(ip),    intent(in)    :: ivari,ivaro
    type(bc_bound), intent(inout) :: tbcod_xxx(:)
    integer(ip)                   :: idofn,ndofn,icode


    call opebcs_initialization()

    if( kfl_icodb /= 0) then
       ndofn =  tbcod_xxx(ivari) % ndofn 
       tbcod_xxx(ivaro) % ndofn = tbcod_xxx(ivari) % ndofn
       tbcod_xxx(ivaro) % ncode = tbcod_xxx(ivari) % ncode  
       do icode = 1,mcodc
          tbcod_xxx(ivaro) % l(icode) % kfl_fixbo  = tbcod_xxx(ivari) % l(icode) % kfl_fixbo
          tbcod_xxx(ivaro) % l(icode) % kfl_value  = tbcod_xxx(ivari) % l(icode) % kfl_value
          tbcod_xxx(ivaro) % l(icode) % kfl_funtyp = tbcod_xxx(ivari) % l(icode) % kfl_funtyp
          tbcod_xxx(ivaro) % l(icode) % kfl_funbo  = tbcod_xxx(ivari) % l(icode) % kfl_funbo
          tbcod_xxx(ivaro) % l(icode) % lcode      = tbcod_xxx(ivari) % l(icode) % lcode
          do idofn = 1,ndofn
             tbcod_xxx(ivaro) % l(icode) % bvnat(idofn) = tbcod_xxx(ivari) % l(icode) % bvnat(idofn) 
          end do
       end do

    end if

  end subroutine cpybcs_boundaries

  subroutine boundary_conditions_read_node_codes(message)

    character(*), intent(in), optional :: message

    if( present(message) ) then
       coutp(1) = trim(message)
    else
       coutp(1) ='UNDEFINED VARIABLE'
    end if
    call reacod(100_ip)

  end subroutine boundary_conditions_read_node_codes

  subroutine boundary_conditions_read_boundary_codes(message)

    character(*), intent(in), optional :: message

    if( present(message) ) then
       coutp(1) = trim(message)
    else
       coutp(1) ='UNDEFINED VARIABLE'
    end if
    call reacod(200_ip)

  end subroutine boundary_conditions_read_boundary_codes

  !-------------------------------------------------------------------
  ! 
  ! Impose node or edge codes (done only for untagged nodes)
  !
  !-------------------------------------------------------------------

  ! alf alf
  ! debug debug
  ! kfl_funtyp should also be included here?
  subroutine boundary_conditions_impose_node_codes(tncod,kfl_fixno,bvess,kfl_funno,kfl_fixrs)

    use def_master, only : lun_outpu
    use def_master, only : intost
    use def_master, only : IMPOSE_NODE_CODES
    use def_master, only : IMPOSE_BOUNDARY_CODES  
    use def_master, only : IMPOSE_EDGE_CODES     
    use def_domain, only : kfl_codno
    use def_domain, only : kfl_coded
    use def_domain, only : mcono,npoin
    use def_domain, only : lpoty
    use def_domain, only : meshe,xfiel
    use def_kermod, only : ndivi
    use mod_chktyp, only : check_type
  
    type(bc_nodes), intent(in)                       :: tncod
    integer(ip),    intent(inout), pointer           :: kfl_fixno(:,:)
    real(rp),       intent(inout), pointer, optional :: bvess(:,:)
    integer(ip),    intent(inout), pointer, optional :: kfl_funno(:)
    integer(ip),    intent(inout), pointer, optional :: kfl_fixrs(:)

    integer(ip)                                      :: ipoin,iboun,icode,idofn,iparb,ibopo,dummi,iauxi
    integer(ip)                                      :: ndofn,ncode,ibves,kpoin,ierro,nnand,iword,itask
    integer(ip)                                      :: jcode,icono,nmcod,kcode(mcono),ivcod,ivcob
    integer(ip)                                      :: kfl_funno_tmp,kfl_funbo_tmp, kfl_funtyp_tmp, ifunc,ktest_size
    integer(ip)                                      :: ntotn,mcono_tmp
    character(20)                                    :: messa
    integer(ip)                                      :: kfl_fixrs_tmp
    character(5)                                     :: wfixrs,wfname,wtag
    integer(ip), pointer                             :: kfl_codes_tmp(:,:)

    itask = IMPOSE_NODE_CODES
    
    if( itask == IMPOSE_NODE_CODES ) then
       ntotn         =  npoin
       kfl_codes_tmp => kfl_codno
       mcono_tmp     =  mcono
    else if( itask == IMPOSE_EDGE_CODES ) then
       ntotn         =  meshe(ndivi) % nedge
       kfl_codes_tmp => kfl_coded
       mcono_tmp     =  2
    end if
    
    ndofn         = memory_size(kfl_fixno,1_ip)
    if( present(bvess) ) ibves = memory_size(bvess,1_ip)
    ibves = ntotn + 1

    do ncode = 1,tncod % ncode  
       !
       ! Do it only when the code is untagged
       !   
       if( itask == IMPOSE_NODE_CODES ) then          
          if(      tncod % l(ncode) % lcode(1) == mcodb+1 ) then
             nmcod = 0
          else if( tncod % l(ncode) % lcode(2) == mcodb+1 ) then
             nmcod = 1
          else if( tncod % l(ncode) % lcode(3) == mcodb+1 ) then
             nmcod = 2
          else 
             nmcod = 3
          end if
       else
          if(      tncod % l(ncode) % lcode(1) == mcodb+1 ) then
             nmcod = 0
          else if( tncod % l(ncode) % lcode(2) == mcodb+1 ) then
             nmcod = 1
          else 
             nmcod = 2
          end if
       end if
       !
       ! Check if value function exist
       !
       if( itask == IMPOSE_NODE_CODES ) then
          if( tncod % l(ncode) % kfl_value > 0 ) then
             ivcod = tncod % l(ncode) % kfl_value
             call check_type(xfiel,ivcod,ndofn,npoin)
          end if
       end if
       !
       ! Order codes
       !
       do jcode = 1,mcono_tmp
          kcode(jcode) = tncod % l(ncode) % lcode(jcode)
       end do
       call heapsorti1(2_ip,mcono_tmp,kcode)
       icode = tncod % l(ncode) % lcode(1)

       do ipoin = 1,ntotn
          icono = 0
          do jcode = 1,mcono_tmp
             if( kfl_codes_tmp(jcode,ipoin) == abs(kcode(jcode)) ) icono = icono + 1
          end do

          if( icono == mcono_tmp ) then

             kpoin = ipoin

             if( kpoin == 0 ) then

                messa = intost(ipoin)
                ierro = ierro + 1
                call outfor(2_ip,lun_outpu,&
                     'BOUNDARY CONDITION CANNOT BE IMPOSED ON INTERIOR NODE: '//trim(messa))
             else

                kfl_fixno(1,kpoin) = tncod % l(ncode) % kfl_fixno

                call codfix(ndofn,kfl_fixno(1,kpoin))

                if( present(bvess) ) then
                   if( ibves == ntotn ) then
                      if( tncod % l(ncode) % kfl_value == 0 ) then
                         bvess(kpoin,1) = tncod % l(ncode) % bvess(1)
                      else                          
                         ivcod = tncod % l(ncode) % kfl_value
                         bvess(kpoin,1) = xfiel(ivcod) % a(1,ipoin,1)
                      end if
                   else
                      if( tncod % l(ncode) % kfl_value == 0 ) then 
                         do idofn = 1,ndofn
                            bvess(idofn,kpoin) = tncod % l(ncode) % bvess(idofn)
                         end do
                      else
                         ivcod = tncod % l(ncode) % kfl_value
                         do idofn = 1,ndofn
                            bvess(idofn,kpoin) = xfiel(ivcod) % a(idofn,ipoin,1)
                         end do
                      end if
                   end if
                end if

                if( present(kfl_funno) ) then
                   kfl_funno(kpoin) = tncod % l(ncode) % kfl_funno                       
                   if( present(kfl_fixrs) ) then 
                      ibopo = lpoty(kpoin)
                      if( ibopo /= 0 ) then                         
                         kfl_fixrs(ibopo) = tncod % l(ncode) % kfl_fixrs
                         !if( kfl_fixrs(ibopo) == -2 .and. nskew > 0 ) call geofix(kpoin,ibopo)
                      end if
                   end if
                else if( present(kfl_fixrs) ) then
                   ibopo = lpoty(kpoin)
                   if( ibopo /= 0 ) then
                      kfl_fixrs(ibopo) = tncod % l(ncode) % kfl_fixrs
                      !if(kfl_fixrs(ibopo)==-2.and.nskew>0) call geofix(kpoin,ibopo)
                   end if
                end if
             end if

          end if

       end do

       !end if

    end do

  end subroutine boundary_conditions_impose_node_codes

end module mod_opebcs

