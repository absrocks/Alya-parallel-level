
!-----------------------------------------------------------------------
!
!> @defgroup Properties_Toolbox
!> Toolbox for creating, updating and computing properties
!> @{
!> @name    ToolBox for property creations
!> @file    mod_ker_proper.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for properties
!> @details To create a new property:
!>          1. ker_readat.f90: Add the reading of the property
!>          2. ker_allaws.f90: Decide in where to computer the property
!>          3. ker_updpro.f90: Compute the property and its gradients if
!>             needed
!
!------------------------------------------------------------------------

module mod_ker_proper

  use def_kintyp,         only : ip,rp 
#ifndef VECTOR_SIZE
  use def_master,         only : VECTOR_SIZE
#endif
  use def_master,         only : modul,mem_modul
  use def_domain,         only : nmate,nelem,npoin,ltypb
  use def_domain,         only : ltype,lmate,ngaus,nboun
  use def_domain,         only : lgaus,ndime,mnode
  use mod_parall,         only : par_omp_nelem_chunk
  use mod_ker_polynomial, only : ker_polynomial_proper
  use mod_ker_polynomial, only : ker_polynomial_name 
  use def_kermod
  implicit none

  interface ker_proper
     module procedure ker_proper_scalar_000, &
          &           ker_proper_scalar_00,  &
          &           ker_proper_scalar_0,   &
          &           ker_proper_scalar_1,   &
          &           ker_proper_scalar_2,   &
          &           ker_proper_vector_0,   &
          &           ker_proper_vector_1,   &
          &           ker_proper_vector_2,   &
          &           ker_proper_vector_3
  end interface ker_proper

  public :: ker_proper_initialization
  public :: ker_proper
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-08
  !> @brief   Initialization
  !> @details Give properties default values
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_initialization() 

    densi_ker % wlaws         = 'CONST'
    visco_ker % wlaws         = 'CONST'
    poros_ker % wlaws         = 'CONST'
    condu_ker % wlaws         = 'CONST'
    sphea_ker % wlaws         = 'CONST'
    dummy_ker % wlaws         = 'CONST'
    turmu_ker % wlaws         = 'CONST'
    absor_ker % wlaws         = 'CONST'
    scatt_ker % wlaws         = 'CONST'
    mixin_ker % wlaws         = 'CONST'

    densi_ker % rlaws         = 0.0_rp
    visco_ker % rlaws         = 0.0_rp
    poros_ker % rlaws         = 0.0_rp
    condu_ker % rlaws         = 0.0_rp
    sphea_ker % rlaws         = 0.0_rp
    dummy_ker % rlaws         = 0.0_rp
    turmu_ker % rlaws         = 0.0_rp
    absor_ker % rlaws         = 0.0_rp
    scatt_ker % rlaws         = 0.0_rp
    mixin_ker % rlaws         = 0.0_rp
    
    densi_ker % value_default = 0.0_rp
    visco_ker % value_default = 0.0_rp
    poros_ker % value_default = 0.0_rp
    condu_ker % value_default = 0.0_rp
    sphea_ker % value_default = 0.0_rp
    dummy_ker % value_default = 0.0_rp
    turmu_ker % value_default = 0.0_rp
    absor_ker % value_default = 0.0_rp
    scatt_ker % value_default = 0.0_rp
    mixin_ker % value_default = 0.0_rp

  end subroutine ker_proper_initialization
  
  subroutine ker_allpro(prope_ker)
    !------------------------------------------------------------------------
    !****f* Master/ker_allpro
    ! NAME
    !    ker_allpro
    ! DESCRIPTION
    !    Allocate memory
    ! USES
    !    postpr
    !    memgen
    ! USED BY
    !    ker_output
    !***
    !------------------------------------------------------------------------
    use mod_memchk
    use mod_memory
    use def_domain, only :  ndime,nnode
    use def_domain, only :  mgaus
    use def_master, only :  nturb
    implicit none
    type(typ_valpr_ker)  :: prope_ker
    integer(ip)          :: kfl_ipoin,kfl_ielem,kfl_const
    integer(ip)          :: ielem,imate,ilaws,pgaus,pgaub
    integer(ip)          :: iboun,kfl_grele,kfl_grpoi,kfl_drele,kfl_gdele,pelty,pnode
    integer(ip)          :: kfl_drele_tur,kfl_drele_vel
    integer(4)           :: istat
    !
    ! Nullify pointers
    !
    nullify( prope_ker % value_ipoin   )
    nullify( prope_ker % grval_ipoin   )
    nullify( prope_ker % value_ielem   )
    nullify( prope_ker % grval_ielem   )
    nullify( prope_ker % value_iboun   )
    nullify( prope_ker % value_const   )
    nullify( prope_ker % drval_ielem   )
    nullify( prope_ker % gdval_ielem   )
    nullify( prope_ker % drval_tur_ielem   )
    nullify( prope_ker % drval_vel_ielem   )
    !
    ! Define some flags
    !
    kfl_ipoin = 0
    kfl_ielem = 0
    kfl_const = 0

    kfl_grele = 0
    kfl_grpoi = 0
    kfl_drele = 0
    kfl_drele_tur = 0
    kfl_drele_vel = 0
    kfl_gdele = 0

    do imate = 1,nmate
       ilaws = prope_ker % ilaws(imate)
       if(      prope_ker % llaws(ilaws) % where == 'IPOIN' ) then
          kfl_ipoin = 1
          if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
             kfl_grpoi = 1
          end if
       else if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
          kfl_ielem = 1
          if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
             kfl_grele = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_deriv == 1 ) then
             kfl_drele = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_deriv_tur == 1 ) then
             kfl_drele_tur = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_deriv_vel == 1 ) then
             kfl_drele_vel = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_grder == 1 ) then
             kfl_gdele = 1
          end if
       else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then
          kfl_const = 1
       end if
    end do
    !
    ! Allocate memory
    !
    if ( (kfl_ipoin == 1 ) .or. ( kfl_ielem == 1 ) ) then
       !
       ! On nodes
       !
       ! even if the law is stored at ielem it might be smoothed if at some moment it is need at ipoin on npoin, so we allocate mem
       ! in order to save memory one could find a way to check if it will need to be smoothed -- needs to be thought
       !
       call memory_alloca(mem_modul(1:2,modul),'value_ipoin','ker_allpro',prope_ker % value_ipoin,npoin)
    end if
    if( ( kfl_ipoin == 1 ) .and. ( kfl_grpoi == 1 ) ) then
       call memory_alloca(mem_modul(1:2,modul),'grval_ipoin','ker_allpro',prope_ker % grval_ipoin,ndime,npoin)
    end if

    if( kfl_ielem == 1 ) then
       !
       ! On elements
       !
       call memory_alloca(mem_modul(1:2,modul),'value_ielem','ker_allpro',prope_ker % value_ielem,nelem)
       do ielem = 1,nelem
          pgaus = ngaus(abs(ltype(ielem)))
          if( kfl_cutel /= 0 ) pgaus = mgaus
          call memory_alloca(mem_modul(1:2,modul),'value_ielem(ielem)%a','ker_allpro',&
               prope_ker % value_ielem(ielem)%a,pgaus)
       end do
       if( kfl_grele == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'grval_ielem','ker_allpro',prope_ker % grval_ielem,nelem)
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             if( kfl_cutel /= 0 ) pgaus = mgaus
             call memory_alloca(mem_modul(1:2,modul),'grval_ielem(ielem)%a','ker_allpro',&
                  prope_ker % grval_ielem(ielem)%a,ndime,pgaus)
          end do
       end if


       if( kfl_drele == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'drval_ielem','ker_allpro',prope_ker % drval_ielem,nelem)
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             call memory_alloca(mem_modul(1:2,modul),'drval_ielem(ielem)%a','ker_allpro',&
                  prope_ker % drval_ielem(ielem)%a,pnode,pgaus)
          end do
       end if

       if( kfl_drele_tur == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'drval_tur_ielem','ker_allpro',prope_ker % drval_tur_ielem,nelem)
          nturb = 1
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             call memory_alloca(mem_modul(1:2,modul),'drval_tur_ielem(ielem)%a','ker_allpro',&
                  prope_ker % drval_tur_ielem(ielem)%a,nturb,pnode,pgaus)
          end do
       end if

       if( kfl_drele_vel == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'drval_vel_ielem','ker_allpro',prope_ker % drval_vel_ielem,nelem)
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             call memory_alloca(mem_modul(1:2,modul),'drval_vel_ielem(ielem)%a','ker_allpro',&
                  prope_ker % drval_vel_ielem(ielem)%a,ndime,pnode,pgaus)
          end do
       end if

       if( kfl_gdele == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'gdval_ielem','ker_allpro',prope_ker % gdval_ielem,nelem)
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             call memory_alloca(mem_modul(1:2,modul),'gdval_ielem(ielem)%a','ker_allpro',&
                  prope_ker % gdval_ielem(ielem)%a,ndime,pnode,pgaus)
          end do
       end if

       !
       ! On boundaries
       !
       call memory_alloca(mem_modul(1:2,modul),'value_iboun','ker_allpro',prope_ker % value_iboun,nboun)
       do iboun = 1,nboun
          pgaub = ngaus(abs(ltypb(iboun)))
          call memory_alloca(mem_modul(1:2,modul),'value_iboun(iboun)%a','ker_allpro',&
               prope_ker % value_iboun(iboun)%a,pgaub)
       end do

    end if
    !
    ! Constant properties
    !
    if( kfl_const == 1 ) then
       call memory_alloca(mem_modul(1:2,modul),'value_const','ker_allpro',prope_ker % value_const,nmate)
    end if

  end subroutine ker_allpro

  subroutine ker_wiprop(prope_ker)
    !------------------------------------------------------------------------
    !****f* Master/ker_wippro
    ! NAME
    !    ker_wippro
    ! DESCRIPTION
    !    Correspondance number vs name of properties
    ! USES
    !    postpr
    !    memgen
    ! USED BY
    !    ker_output
    !***
    !------------------------------------------------------------------------
    implicit none
    type(typ_valpr_ker)  :: prope_ker
    integer(ip)          :: imate,ilaws

    do imate = 1,nmate
       ilaws = 1
       do while( ilaws <=  mlaws_ker )
          if( prope_ker % llaws(ilaws) % wname == prope_ker % wlaws(imate) ) then
             prope_ker % ilaws(imate) = ilaws
             ilaws = mlaws_ker + 1
          end if
          ilaws = ilaws + 1
       end do
       if( ilaws /= mlaws_ker + 2 ) then
          write(*,*) 'KER_WIPROP: PROPERTY LAW DOES NOT EXIST:imate,trim(prope_ker % wlaws(imate))',&
                                                              imate,trim(prope_ker % wlaws(imate))
          ! this can happen if in your geo file you have more materias than those you define in ker.dat
          call runend('KER_WIPROP: PROPERTY LAW DOES NOT EXIST')
       end if
    end do

  end subroutine ker_wiprop

  subroutine ker_smopro(prope_ker,xvalu)
    !------------------------------------------------------------------------
    !****f* Master/ker_smopro
    ! NAME
    !    ker_smopro
    ! DESCRIPTION
    !    Smooth a property
    ! USES
    !    postpr
    !    memgen
    ! USED BY
    !    ker_output
    !***
    !------------------------------------------------------------------------
    use def_kintyp, only              :  ip,rp,r1p
    use def_kermod, only              :  typ_valpr_ker
    use def_master, only              :  INOTMASTER
    use def_domain, only              :  ndime,npoin,nelem,nnode,mnode,ntens
    use def_domain, only              :  lnods,ltype,coord,vmass,elmar,vmasc
    use def_domain, only              :  lexis,ngaus,kfl_naxis,mgaus,lelch
    use def_domain, only              :  lmate,lmatn,nmate
    use def_elmtyp
    use mod_memory
    implicit none
    type(typ_valpr_ker), pointer      :: prope_ker
    real(rp),     intent(out), target :: xvalu(*)
    integer(ip)                       :: ipoin,idime,inode,ielem,igaus
    integer(ip)                       :: pnode,pelty,pgaus,ilaws,imate
    integer(ip)                       :: jnode,knode
    integer(4)                        :: istat
    real(rp)                          :: detjm,gpvol,gpcar(ndime,mnode)
    real(rp)                          :: elcod(ndime,mnode)
    real(rp)                          :: xjaci(9),xjacm(9),xfact
    integer(ip),  pointer             :: lmele(:)
    real(rp),     pointer             :: vmass_ker(:)

    if( INOTMASTER ) then
       !
       ! Smooth property
       !
       nullify( lmele )
       nullify( vmass_ker )
       call memory_alloca(mem_modul(1:2,modul),'LMELE'    ,'ker_smopro',lmele,nelem)
       call memory_alloca(mem_modul(1:2,modul),'VMASS_KER','ker_smopro',vmass_ker,npoin)
       do imate = 1,nmate
          ilaws = prope_ker % ilaws(imate)
          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
             do ielem = 1,nelem
                if( lmate(ielem) == imate ) lmele(ielem) = 1
             end do
          end if
       end do
       !
       ! Initialization
       !
       do ipoin = 1,npoin
          vmass_ker(ipoin) = 0.0_rp
          xvalu(ipoin)     = 0.0_rp
       end do
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          pelty = ltype(ielem)

          if( pelty > 0 .and. lmele(ielem) == 1 ) then
             pnode = nnode(pelty)
             pgaus = ngaus(pelty)
             !
             ! Gather vectors
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                do idime = 1,ndime
                   elcod(idime,inode) = coord(idime,ipoin)
                end do
             end do
             !
             ! Loop over Gauss points
             !
             gauss_points: do igaus = 1,pgaus
                call elmder(&
                     pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                     elcod,gpcar,detjm,xjacm,xjaci)
                gpvol = elmar(pelty) % weigp(igaus) * detjm
                if( kfl_naxis == 1 ) then
                   call runend('MOD_GRADIE: NOT CODED')
                end if
                !
                ! Extension
                !
                if( lelch(ielem) == ELEXT ) then
                   knode = 1
                else
                   knode = pnode
                end if
                !
                ! Assemble
                !
                do inode = 1,knode
                   ipoin        = lnods(inode,ielem)
                   xfact        = gpvol * elmar(pelty) % shape(inode,igaus)
                   xvalu(ipoin) = xvalu(ipoin) + xfact * prope_ker % value_ielem(ielem) % a(igaus)
                   do jnode = 1,pnode
                      vmass_ker(ipoin) = vmass_ker(ipoin) + xfact * elmar(pelty) % shape(jnode,igaus)
                   end do
                end do

             end do gauss_points
          end if
       end do elements
       !
       ! Parallelization
       !
       call rhsmod(1_ip,xvalu)
       call rhsmod(1_ip,vmass_ker)
       !memory_alloca
       ! Solve diagonal system
       !
       do ipoin = 1,npoin
          if( vmass_ker(ipoin) /= 0.0_rp ) &
               xvalu(ipoin) = xvalu(ipoin) / vmass_ker(ipoin)
       end do
       !
       ! Deallocate memory
       !
       call memory_deallo(mem_modul(1:2,modul),'VMASS_KER','ker_smopro',vmass_ker)
       call memory_deallo(mem_modul(1:2,modul),'LMELE'    ,'ker_smopro',lmele)

    end if

  end subroutine ker_smopro

  subroutine ker_grapro(prope_ker,xvalu)
    !------------------------------------------------------------------------
    !****f* Master/ker_grapro
    ! NAME
    !    ker_grapro
    ! DESCRIPTION
    !    Smooth the gradient of a property
    ! USES
    !    postpr
    !    memgen
    ! USED BY
    !    ker_output
    !***
    !------------------------------------------------------------------------
    use def_kintyp, only              :  ip,rp,r1p
    use def_kermod, only              :  typ_valpr_ker
    use def_master, only              :  INOTMASTER
    use def_domain, only              :  ndime,npoin,nelem,nnode,mnode,ntens
    use def_domain, only              :  lnods,ltype,coord,vmass,elmar,vmasc
    use def_domain, only              :  lexis,ngaus,kfl_naxis,mgaus,lelch
    use def_domain, only              :  lmate,lmatn,nmate
    use def_elmtyp
    use mod_memory
    implicit none
    type(typ_valpr_ker), pointer      :: prope_ker
    real(rp),     intent(out), target :: xvalu(ndime,npoin)
    integer(ip)                       :: ipoin,idime,inode,ielem,igaus
    integer(ip)                       :: pnode,pelty,pgaus,ilaws,imate
    integer(ip)                       :: jnode,knode
    integer(4)                        :: istat
    real(rp)                          :: detjm,gpvol,gpcar(ndime,mnode)
    real(rp)                          :: elcod(ndime,mnode)
    real(rp)                          :: grunk(ndime)
    real(rp)                          :: xjaci(9),xjacm(9),xfact
    integer(ip),  pointer             :: lmele(:)
    real(rp),     pointer             :: vmass_ker(:)

    if( INOTMASTER ) then
       !
       ! Smooth property
       !
       nullify( lmele )
       nullify( vmass_ker )
       call memory_alloca(mem_modul(1:2,modul),'LMELE'    ,'ker_grapro',lmele,nelem)
       call memory_alloca(mem_modul(1:2,modul),'VMASS_KER','ker_grapro',vmass_ker,npoin)
       do imate = 1,nmate
          ilaws = prope_ker % ilaws(imate)
          if(    ( prope_ker % llaws(ilaws) % where == 'IELEM' .and.&
               &   prope_ker % llaws(ilaws) % kfl_gradi == 1 ) .or. &
               & ( prope_ker % llaws(ilaws) % where == 'IPOIN' .and.&
               &   prope_ker % llaws(ilaws) % kfl_gradi == 0 ) ) then
             do ielem = 1,nelem
                if( lmate(ielem) == imate ) lmele(ielem) = ilaws
             end do
          end if
       end do
       !
       ! Initialization
       !
       do ipoin = 1,npoin
          vmass_ker(ipoin) = 0.0_rp
       end do
       do ipoin = 1,npoin
          do idime = 1,ndime
             xvalu(idime,ipoin) = 0.0_rp
          end do
       end do
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          pelty = ltype(ielem)

          if( pelty > 0 .and. lmele(ielem) /= 0 ) then
             pnode = nnode(pelty)
             pgaus = ngaus(pelty)
             ilaws = lmele(ielem)
             !
             ! Gather vectors
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                do idime = 1,ndime
                   elcod(idime,inode) = coord(idime,ipoin)
                end do
             end do
             !
             ! Loop over Gauss points
             !
             gauss_points: do igaus = 1,pgaus
                call elmder(&
                     pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                     elcod,gpcar,detjm,xjacm,xjaci)
                gpvol = elmar(pelty) % weigp(igaus) * detjm

                if( kfl_naxis == 1 ) then
                   call runend('MOD_GRADIE: NOT CODED')
                end if
                !
                ! Extension
                !
                if( lelch(ielem) == ELEXT ) then
                   knode = 1
                else
                   knode = pnode
                end if
                !
                ! Assemble
                !
                if(        prope_ker % llaws(ilaws) % where == 'IELEM' &
                     .and. prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                      do idime = 1,ndime
                         xvalu(idime,ipoin) = xvalu(idime,ipoin) &
                              + xfact * prope_ker % grval_ielem(ielem) % a(idime,igaus)
                      end do
                      do jnode = 1,knode
                         vmass_ker(ipoin) = vmass_ker(ipoin) + xfact * elmar(pelty) % shape(jnode,igaus)
                      end do
                   end do

                else if(   prope_ker % llaws(ilaws) % where == 'IPOIN' &
                     .and. prope_ker % llaws(ilaws) % kfl_gradi == 0 ) then

                   grunk = 0.0_rp
                   do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      do idime = 1,ndime
                         grunk(idime) = grunk(idime) + prope_ker % value_ipoin(ipoin) * gpcar(idime,inode)
                      end do
                   end do
                   do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                      do idime = 1,ndime
                         xvalu(idime,ipoin) = xvalu(idime,ipoin) + xfact * grunk(idime)
                      end do
                      do jnode = 1,pnode
                         vmass_ker(ipoin) = vmass_ker(ipoin) + xfact * elmar(pelty) % shape(jnode,igaus)
                      end do
                   end do

                end if

             end do gauss_points
          end if
       end do elements
       !
       ! Parallelization
       !
       call rhsmod(ndime,xvalu)
       call rhsmod( 1_ip,vmass_ker)
       !
       ! Solve diagonal system
       !
       do ipoin = 1,npoin
          if( vmass_ker(ipoin) /= 0.0_rp ) then
             do idime = 1,ndime
                xvalu(idime,ipoin) = xvalu(idime,ipoin) / vmass_ker(ipoin)
             end do
          end if
       end do
       !
       ! Deallocate memory
       !
       call memory_deallo(mem_modul(1:2,modul),'VMASS_KER','ker_grapro',vmass_ker)
       call memory_deallo(mem_modul(1:2,modul),'LMELE'    ,'ker_grapro',lmele)

    end if

  end subroutine ker_grapro

  subroutine ker_proper_scalar(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    !-----------------------------------------------------------------------
    !****f* Kermod/ker_proper
    ! NAME
    !    ker_proper
    ! DESCRIPTION
    !    Update properties
    !
    !    Property can be defined with 3 formats
    !    It can be required at 8 different places
    !
    !       WHEREIN   IPOSI  IELBO
    !       --------------------
    !    1. IPOIN   IPOIN  NONSENSE
    !    2. NPOIN   DUMMI  Do smoothing => Loop over IELEM
    !    3. IGAUS   IGAUS  IELEM
    !    4. PGAUS   DUMMI  IELEM
    !    5. IGAUB   IGAUB  IBOUN        => IELEM = LBOEL(,IBOUN)
    !    6. PGAUB   DUMMI  IBOUN        => IELEM = LBOEL(,IBOUN)
    !    7. COG     DUMMI  IELEM
    !    8. PNODE   DUMMI  IELEM
    !
    !    GPSHA is dimensioned as a vector in order to use it as optional
    !    argument
    !
    ! USED BY
    !    many subroutines
    !***
    !-----------------------------------------------------------------------

    use def_kintyp,      only                 :  ip,rp
    use def_elmtyp,      only                 :  ELCUT
    use def_master,      only                 :  mem_modul,modul,nturb
    use def_domain,      only                 :  nmate,npoin,nelem,mnode,lelch
    use def_domain,      only                 :  lnods,lnnod,lmate,lmatn,nmatn
    use def_domain,      only                 :  ngaus,ltype,elmar,ndime,mgaus
    use def_domain,      only                 :  lelbo,ltypb,lnodb,nnode
    use def_kermod,      only                 :  typ_valpr_ker
    use def_kermod,      only                 :  densi_ker,visco_ker
    use def_kermod,      only                 :  poros_ker,condu_ker
    use def_kermod,      only                 :  sphea_ker,dummy_ker
    use def_kermod,      only                 :  absor_ker,scatt_ker
    use def_kermod,      only                 :  turmu_ker,mixin_ker
    use mod_memory
    implicit none
    character(5),        intent(in)           :: wname,wherein
    integer(ip),         intent(in)           :: iposi,ielbo
    real(rp),            intent(out)          :: xvalu(*)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(*)
    real(rp),            intent(in), optional :: gpcar(ndime,mnode,*)
    real(rp),            intent(in), optional :: gpcor(ndime,*)
    type(typ_valpr_ker), pointer              :: prope_ker
    integer(ip)                               :: imate,ielem,inode,ipoin,igaus,kauxi,iturb
    integer(ip)                               :: kpoin,ilaws,pelty,pgaus,pnode,jpoin
    integer(ip)                               :: pnodb,igaub,inodb,pblty,iboun
    integer(ip)                               :: pgaub,kposi,kfl_gradi,idime,kfl_deriv,kfl_grder
    integer(ip)                               :: kfl_deriv_tur,kfl_deriv_vel
    integer(ip)                               :: lposi
    real(rp)                                  :: elpro(mnode)
    real(rp)                                  :: elgrp(ndime,max(mgaus,mnode))
    real(rp),   pointer                       :: auxva(:)

    nullify(auxva)

    !----------------------------------------------------------------------
    !
    ! Some definitions
    !
    !----------------------------------------------------------------------

    kfl_gradi     = 0
    kfl_deriv     = 0
    kfl_grder     = 0
    kfl_deriv_tur = 0
    kfl_deriv_vel = 0

    if(      wname(1:2) == 'GR'  ) then
       kfl_gradi = 1
    else if( wname(1:2) == 'DR'  ) then
       kfl_deriv = 1
    else if( wname(1:2) == 'GD'  ) then
       kfl_grder = 1
    else if( wname(1:3) == 'TDR' ) then
       kfl_deriv_tur = 1
    else if( wname(1:3) == 'VDR' ) then
       kfl_deriv_vel = 1
    end if

    if(      wname == 'DENSI' .or. wname == 'GRDEN' .or. wname == 'DRDEN' .or. wname == 'GDDEN' ) then
       prope_ker => densi_ker
    else if( wname == 'VISCO' .or. wname == 'GRVIS' .or. wname == 'DRVIS' .or. wname == 'GDVIS' ) then
       prope_ker => visco_ker
    else if( wname == 'MIXIN' ) then
       prope_ker => mixin_ker
    else if( wname == 'POROS' .or. wname == 'GRPOS' ) then
       prope_ker => poros_ker
    else if( wname == 'CONDU' .or. wname == 'GRCON' ) then
       prope_ker => condu_ker
    else if( wname == 'SPHEA' .or. wname == 'GRSPE' ) then
       prope_ker => sphea_ker
    else if( wname == 'DUMMY' .or. wname == 'GRDUM' ) then
       prope_ker => dummy_ker
    else if( wname == 'TURBU' .or. wname == 'TURVI' .or. wname == 'GRTUR' .or. wname == 'TDRTU' .or. wname == 'VDRTU' ) then
       prope_ker => turmu_ker
    else if( wname == 'ABSOR' ) then
       prope_ker => absor_ker
    else if( wname == 'SCATT' ) then
       prope_ker => scatt_ker
    else
       print *,'--- I DONT UNDERSTAND ',wname
       call runend('KER_PROPER CALLED WITH WRONG PROPERTY') ! This avoids the bug of calling wrongly the ker_proper routine
    end if

    if( prope_ker % kfl_exist == 0 ) then
       !
       ! Put value to zero
       !
       if( wherein == 'IPOIN' .or. wherein == 'IGAUS' .or. wherein == 'IGAUB' .or. wherein == 'COG' ) then
          kposi = 1
       else if( wherein == 'NPOIN' ) then
          kposi = npoin
       else if( wherein == 'PNODE' ) then
          kposi = lnnod(ielbo)
       else if( wherein == 'PGAUS' ) then
          if( present(gpsha) ) then
             if( present(qgaus) ) then
                kposi = qgaus
             else
                call runend('MOD_KER_PROPER: ARGUMENT IS MISSING')
             end if
          else
             kposi = ngaus(abs(ltype(ielbo)))
          end if
       else if( wherein == 'PGAUB' ) then
          if( present(gpsha) ) then
             if( present(qgaus) ) then
                kposi = qgaus
             else
                call runend('MOD_KER_PROPER: ARGUMENT IS MISSING')
             end if
          else
             kposi = ngaus(abs(ltypb(ielbo)))
          end if
       end if
       if( kfl_gradi == 1 ) kposi = kposi * ndime
       do lposi = 1,kposi
          xvalu(lposi) = 0.0_rp
       end do
       
    else

       !----------------------------------------------------------------------
       !
       ! 1. Property required on IPOIN
       !
       !----------------------------------------------------------------------

       if( wherein == 'IPOIN' ) then

          if( kfl_gradi == 1 ) then
             call runend('KER_PROPER: DONT KNOW WHAT TO DO WITH THIS OPTION')
          end if

          ipoin = iposi
          !
          ! Smooth value - If in some material it is stored by ielem it is smoothed
          ! In the materials where it is not stored by ielem it will step over with erroneous values
          ! Therefore if is then recalculated in materials where it is stored by IPOIN or CONST
          !
          kauxi = 0
          if ( prope_ker % kfl_nedsm == 1 ) then
             imat2_loop: do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)
                if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
                   call memory_alloca(mem_modul(1:2,modul),'auxva','ker_proper',auxva,npoin)
                   call ker_smopro(prope_ker,auxva)
                   prope_ker % kfl_nedsm = 0
                   kauxi = 1
                   exit imat2_loop
                end if
             end do imat2_loop
          end if

          if (kauxi == 1) then  ! correct the value (auxva) in points where it is calculated by ipoin or const
             ! and has been steped over when smoothing due to materials whith law at ielem
             do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)

                if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

                   do kpoin = 1,nmatn(imate)
                      jpoin = lmatn(imate) % l(kpoin)
                      auxva(jpoin) = prope_ker % value_ipoin(jpoin)
                   end do

                else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                   do kpoin = 1,nmatn(imate)
                      jpoin = lmatn(imate) % l(kpoin)
                      auxva(jpoin) = prope_ker % value_const(imate)
                   end do

                end if
             end do
             !
             ! Save values obtained in auxva to prope_ker % value_ipoin(jpoin)
             !
             do kpoin =1,npoin
                prope_ker % value_ipoin(kpoin) = auxva(kpoin)
             end do
             call memory_deallo(mem_modul(1:2,modul),'auxva','ker_proper',auxva)
          end if


          do imate = 1,nmate
             ilaws = prope_ker % ilaws(imate)

             if( ( prope_ker % llaws(ilaws) % where == 'IPOIN' ) .or. ( prope_ker % llaws(ilaws) % where == 'IELEM' ) )then

                xvalu(1) = prope_ker % value_ipoin(ipoin)

             else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                xvalu(1) = prope_ker % value_const(imate)

             end if
          end do
       end if

       !----------------------------------------------------------------------
       !
       ! 2. Property required on NPOIN
       !
       !----------------------------------------------------------------------

       if( wherein == 'NPOIN' ) then

          if( kfl_gradi == 1 ) then
             !
             ! Gradient
             !
             imat0_loop: do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)
                if(    ( prope_ker % llaws(ilaws) % where == 'IELEM' .and.&
                     &   prope_ker % llaws(ilaws) % kfl_gradi == 1 ) .or. &
                     & ( prope_ker % llaws(ilaws) % where == 'IPOIN' .and.&
                     &   prope_ker % llaws(ilaws) % kfl_gradi == 0 ) ) then
                   call ker_grapro(prope_ker,xvalu)
                   exit imat0_loop
                end if
             end do imat0_loop

             do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)

                if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

                   if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                      do kpoin = 1,nmatn(imate)
                         ipoin = lmatn(imate) % l(kpoin)
                         kposi = (ipoin-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % grval_ipoin(idime,ipoin)
                         end do
                      end do

                   end if

                else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      kposi = (ipoin-1) * ndime
                      do idime = 1,ndime
                         kposi = kposi + 1
                         xvalu(kposi) = 0.0_rp
                      end do
                   end do

                end if
             end do

          else
             !
             ! Value - If in some material it is stored by ielem it is smoothed
             ! In the materials where it is not stored by ielem it will step over with erroneous values
             ! Therefore if is then recalculated in materials where it is stored by IPOIN or CONST
             !
             kauxi = 0
             if ( prope_ker % kfl_nedsm == 1 ) then
                imat1_loop: do imate = 1,nmate
                   ilaws = prope_ker % ilaws(imate)
                   if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
                      call memory_alloca(mem_modul(1:2,modul),'auxva','ker_proper',auxva,npoin)
                      call ker_smopro(prope_ker,auxva)
                      prope_ker % kfl_nedsm = 0
                      kauxi = 1
                      exit imat1_loop
                   end if
                end do imat1_loop
             end if

             if (kauxi == 1) then  ! correct the value (auxva) in points where it is calculated by ipoin or const
                ! and has been steped over when smoothing due to materials whith law at ielem
                do imate = 1,nmate
                   ilaws = prope_ker % ilaws(imate)

                   if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

                      do kpoin = 1,nmatn(imate)
                         jpoin = lmatn(imate) % l(kpoin)
                         auxva(jpoin) = prope_ker % value_ipoin(jpoin)
                      end do

                   else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                      do kpoin = 1,nmatn(imate)
                         jpoin = lmatn(imate) % l(kpoin)
                         auxva(jpoin) = prope_ker % value_const(imate)
                      end do

                   end if
                end do
                !
                ! Save values obtained in auxva to prope_ker % value_ipoin(jpoin)
                !
                do kpoin =1,npoin
                   prope_ker % value_ipoin(kpoin) = auxva(kpoin)
                end do
                call memory_deallo(mem_modul(1:2,modul),'auxva','ker_proper',auxva)
             end if

             do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)

                if( ( prope_ker % llaws(ilaws) % where == 'IPOIN' ) .or. ( prope_ker % llaws(ilaws) % where == 'IELEM' ) )then

                   do kpoin = 1,nmatn(imate)
                      jpoin = lmatn(imate) % l(kpoin)
                      xvalu(jpoin) = prope_ker % value_ipoin(jpoin)
                   end do

                else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                   do kpoin = 1,nmatn(imate)
                      jpoin = lmatn(imate) % l(kpoin)
                      xvalu(jpoin) = prope_ker % value_const(imate)
                   end do

                end if
             end do
          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 3. Property required on IGAUS
       !
       !----------------------------------------------------------------------

       if( wherein == 'IGAUS' ) then

          igaus = iposi
          ielem = ielbo
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   do idime = 1,ndime
                      xvalu(idime) = prope_ker % grval_ielem(ielem) % a(idime,igaus)
                   end do
                else
                   xvalu(1:ndime) = 0.0_rp
                end if
             else
                !
                ! Value
                !
                if( present(gpsha) ) then
                   pelty = abs(ltype(ielem))
                   pnode = lnnod(ielem)
                   pgaus = ngaus(pelty)
                   xvalu(1) = 0.0_rp
                   do inode = 1,pnode
                      do igaus = 1,pgaus
                         xvalu(1) = xvalu(1) + gpsha(inode) * elmar(pelty) % shaga(igaus,inode) * prope_ker % value_ielem(ielem) % a(igaus)
                      end do
                   end do
                else
                   xvalu(1) = prope_ker % value_ielem(ielem) % a(igaus)
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             pelty = abs(ltype(ielem))
             pnode = lnnod(ielem)
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elpro(inode) = prope_ker % value_ipoin(ipoin)
             end do

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                xvalu(1:ndime) = 0.0_rp

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      do idime = 1,ndime
                         elgrp(idime,inode) = prope_ker % grval_ipoin(idime,ipoin)
                      end do
                   end do
                   if( present(gpsha) ) then
                      kposi = (igaus-1) * qnode
                      do inode = 1,pnode
                         kposi = kposi + 1
                         do idime = 1,ndime
                            xvalu(idime) = xvalu(idime) + elgrp(idime,inode) * gpsha(kposi)
                         end do
                      end do
                   else
                      do inode = 1,pnode
                         do idime = 1,ndime
                            xvalu(idime) = xvalu(idime) + elgrp(idime,inode) * elmar(pelty) % shape(idime,igaus)
                         end do
                      end do
                   end if

                else

                   if( present(gpcar) ) then
                      do inode = 1,pnode
                         do idime = 1,ndime
                            xvalu(idime) = xvalu(idime) + elpro(inode) * gpcar(idime,inode,1)
                         end do
                      end do
                   else
                      call runend('MOD_KER_PROPER: GPCAR NEEDED')
                   end if
                end if

             else
                !
                ! Value
                !
                xvalu(1) = 0.0_rp
                if( present(gpsha) ) then
                   kposi    = (igaus-1) * qnode
                   do inode = 1,pnode
                      kposi    = kposi + 1
                      xvalu(1) = xvalu(1) + elpro(inode) * gpsha(kposi)
                   end do
                else
                   do inode = 1,pnode
                      xvalu(1) = xvalu(1) + elpro(inode) * elmar(pelty) % shape(inode,igaus)
                   end do
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 ) then
                xvalu(1:ndime) = 0.0_rp
             else
                xvalu(1) = prope_ker % value_const(imate)
             end if
          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 4. Property required on IGAUS=1,PGAUS
       !
       !----------------------------------------------------------------------

       if( wherein == 'PGAUS' ) then

          ielem = ielbo
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)
          pelty = abs(ltype(ielem))
          pgaus = ngaus(pelty)
          if( lelch(ielem) == ELCUT ) pgaus = lgaus(ielem)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus .and. kfl_cutel == 0 ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % grval_ielem(ielem) % a(idime,igaus)
                         end do
                      end do
                   else
                      if( qgaus /= pgaus .and. kfl_cutel == 0 ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % grval_ielem(ielem) % a(idime,igaus)
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*ndime
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_deriv == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1) * qnode
                         do inode = 1,qnode
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % drval_ielem(ielem) % a(inode,igaus)
                         end do
                      end do
                   else
                      do igaus = 1,pgaus
                         kposi = (igaus - 1) * qnode
                         do inode = 1,qnode
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % drval_ielem(ielem) % a(inode,igaus)
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_grder == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_grder == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1)* qnode*ndime
                         do inode = 1,qnode
                            do idime = 1, ndime
                                kposi = kposi + 1
                                xvalu(kposi) = prope_ker % gdval_ielem(ielem) % a(idime,inode,igaus)
                            enddo
                         end do
                      end do
                   else
                      do igaus = 1,pgaus
                         do inode = 1,qnode
                            kposi = (igaus - 1)*(inode-1) * qnode*ndime
                            do idime = 1, ndime
                                kposi = kposi + 1
                                xvalu(kposi) = prope_ker % gdval_ielem(ielem) % a(idime,inode,igaus)
                            enddo
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*ndime
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv_tur == 1 ) then

                nturb = 1
                if( prope_ker % llaws(ilaws) % kfl_deriv_tur == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1)* qnode*nturb
                         do inode = 1,qnode
                            do iturb = 1, nturb
                                kposi = kposi + 1
                                xvalu(kposi) = prope_ker % drval_tur_ielem(ielem) % a(iturb,inode,igaus)
                            enddo
                         end do
                      end do
                   else
                      do igaus = 1,pgaus
                         do inode = 1,qnode
                            kposi = (igaus - 1)*(inode-1) * qnode*nturb
                            do iturb = 1, nturb
                                kposi = kposi + 1
                                xvalu(kposi) = prope_ker % drval_tur_ielem(ielem) % a(iturb,inode,igaus)
                            enddo
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*nturb
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv_vel == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_deriv_vel == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1)* qnode*ndime
                         do inode = 1,qnode
                            do idime = 1, ndime
                                kposi = kposi + 1
                                xvalu(kposi) = prope_ker % drval_vel_ielem(ielem) % a(idime,inode,igaus)
                            enddo
                         end do
                      end do
                   else
                      do igaus = 1,pgaus
                         do inode = 1,qnode
                            kposi = (igaus - 1)*(inode-1) * qnode*ndime
                            do idime = 1, ndime
                                kposi = kposi + 1
                                xvalu(kposi) = prope_ker % drval_vel_ielem(ielem) % a(idime,inode,igaus)
                            enddo
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*ndime
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else

                if( present(gpsha) .and. kfl_cutel == 0 ) then
                   if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                   do igaus = 1,qgaus
                      xvalu(igaus) = prope_ker % value_ielem(ielem) % a(igaus)
                   end do
                else
                   do igaus = 1,pgaus
                      xvalu(igaus) = prope_ker % value_ielem(ielem) % a(igaus)
                   end do
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             pelty = abs(ltype(ielem))
             pnode = lnnod(ielem)
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elpro(inode) = prope_ker % value_ipoin(ipoin)
             end do

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do kposi = 1,pgaus*ndime
                   xvalu(kposi) = 0.0_rp
                end do

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      do idime = 1,ndime
                         elgrp(idime,inode) = prope_ker % grval_ipoin(idime,ipoin)
                      end do
                   end do
                   if( present(gpsha) ) then
                      do igaus = 1,qgaus
                         kposi = (igaus-1) * qnode
                         do inode = 1,pnode
                            kposi = kposi + 1
                            lposi = (igaus-1) * ndime
                            do idime = 1,ndime
                               lposi = lposi + 1
                               xvalu(lposi) = xvalu(lposi) + elgrp(idime,inode) * gpsha(kposi)
                            end do
                         end do
                      end do
                   else
                      do igaus = 1,qgaus
                         do inode = 1,pnode
                            lposi = (igaus-1) * ndime
                            do idime = 1,ndime
                               lposi = lposi + 1
                               xvalu(lposi) = xvalu(lposi) + elgrp(idime,inode) * elmar(pelty) % shape(inode,igaus)
                            end do
                         end do
                      end do
                   end if

                else

                   if( present(gpcar) ) then
                      do igaus = 1,qgaus
                         kposi = (igaus-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            do inode = 1,pnode
                               xvalu(kposi) = xvalu(kposi) + elpro(inode) * gpcar(idime,inode,igaus)
                            end do
                         end do
                      end do
                   end if

                end if

             else
                !
                ! Value
                !
                if( present(gpsha) ) then
                   do igaus = 1,qgaus
                      xvalu(igaus) = 0.0_rp
                      kposi        = (igaus-1) * qnode
                      do inode = 1,pnode
                         kposi = kposi + 1
                         xvalu(igaus) = xvalu(igaus) + elpro(inode) * gpsha(kposi)
                      end do
                   end do
                else
                   do igaus = 1,ngaus(pelty)
                      xvalu(igaus) = 0.0_rp
                      do inode = 1,pnode
                         xvalu(igaus) = xvalu(igaus) + elpro(inode) * elmar(pelty) % shape(inode,igaus)
                      end do
                   end do
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do kposi = 1,pgaus*ndime
                   xvalu(kposi) = 0.0_rp
                end do

             else
                !
                ! Value
                !
                if( present(gpsha) ) then
                   do igaus = 1,qgaus
                      xvalu(igaus) = prope_ker % value_const(imate)
                   end do
                else
                   do igaus = 1,ngaus(abs(ltype(ielem)))
                      xvalu(igaus) = prope_ker % value_const(imate)
                   end do
                end if

             end if

          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 5. Property required on IGAUB
       !
       !----------------------------------------------------------------------

       if( wherein == 'IGAUB' ) then

          if( kfl_gradi == 1 )  call runend('MOD_KER_PROPER: GRADIENT NOT AVILABLE ON BOUNDARY')
          igaub = iposi
          iboun = ielbo
          pblty = abs(ltypb(iboun))
          pnodb = nnode(pblty)
          ielem = lelbo(iboun)
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( present(gpsha) ) then
                if( qgaus /= ngaus(pblty) ) call runend('MOD_KER_PROPER: IGAUB NOT CODED')
             end if
             xvalu(1) = prope_ker % value_iboun(iboun) % a(igaub)

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                elpro(inodb) = prope_ker % value_ipoin(ipoin)
             end do
             xvalu(1) = 0.0_rp
             if( present(gpsha) ) then
                kposi = (igaub-1) * qnode
                do inodb = 1,pnodb
                   kposi = kposi + 1
                   xvalu(1) = xvalu(1) + elpro(inodb) * gpsha(kposi)
                end do
             else
                do inodb = 1,pnodb
                   xvalu(1) = xvalu(1) + elpro(inodb) * elmar(pblty) % shape(inodb,igaub)
                end do
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             xvalu(1) = prope_ker % value_const(imate)

          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 6. Property required on IGAUB=1,PGAUB
       !
       !----------------------------------------------------------------------

       if( wherein == 'PGAUB' ) then

          if( kfl_gradi == 1 ) call runend('MOD_KER_PROPER: GRADIENT NOT AVILABLE ON BOUNDARY')

          iboun = ielbo
          pblty = abs(ltypb(iboun))
          pnodb = nnode(pblty)
          ielem = lelbo(iboun)
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( present(gpsha) ) then
                if( qgaus /= ngaus(pblty) ) call runend('MOD_KER_PROPER: IGAUB NOT CODED')
             end if
             do igaub = 1,ngaus(pblty)
                xvalu(igaub) = prope_ker % value_iboun(iboun) % a(igaub)
             end do

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                elpro(inodb) = prope_ker % value_ipoin(ipoin)
             end do

             if( present(gpsha) ) then
                do igaub = 1,qgaus
                   xvalu(igaub) = 0.0_rp
                   kposi = (igaub-1) * qnode
                   do inodb = 1,pnodb
                      kposi        = kposi + 1
                      xvalu(igaub) = xvalu(igaub) + elpro(inodb) * gpsha(kposi)
                   end do
                end do
             else
                pgaub = ngaus(pblty)
                do igaub = 1,pgaub
                   xvalu(igaub) = 0.0_rp
                   do inodb = 1,pnodb
                      xvalu(igaub) = xvalu(igaub) + elpro(inodb) * elmar(pblty) % shape(inodb,igaub)
                   end do
                end do
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             do igaub = 1,ngaus(abs(ltypb(iboun)))
                xvalu(igaub) = prope_ker % value_const(imate)
             end do

          end if

       else if( wherein(1:3) == 'COG' ) then

          !----------------------------------------------------------------------
          !
          ! 7. Property required on C.O.G.
          !
          !----------------------------------------------------------------------

          ielem = ielbo
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)
          pelty = abs(ltype(ielem))
          pgaus = ngaus(pelty)
          pnode = nnode(pelty) 
          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 )  then
                !
                ! Gradient
                !
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   xvalu(1:ndime) = 0.0_rp
                   do igaus = 1,pgaus
                      do idime = 1,ndime
                         xvalu(idime) = xvalu(idime) + prope_ker % grval_ielem(ielem) % a(idime,igaus)
                      end do
                   end do
                   do idime = 1,ndime
                      xvalu(idime) = xvalu(idime) / real(pgaus,rp)
                   end do
                else
                   call runend('MOD_KER_PROPER: GRDAIENT ON COG NOT CODED')
                end if

             else
                !
                ! Value
                !
                xvalu(1) = 0.0_rp
                if( present(gpsha) ) then
                   do igaus = 1,qgaus
                      xvalu(1) = xvalu(1) + prope_ker % value_ielem(ielem) % a(igaus)
                   end do
                   xvalu(1) = xvalu(1) / real(qgaus,rp)
                else
                   do igaus = 1,pgaus
                      xvalu(1) = xvalu(1) + prope_ker % value_ielem(ielem) % a(igaus)
                   end do
                   xvalu(1) = xvalu(1) / real(pgaus,rp)
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             if( kfl_gradi == 1 )  then
                !
                ! Gradient
                !
                xvalu(1:ndime) = 0.0_rp

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      do idime = 1,ndime
                         xvalu(idime) = xvalu(idime) + prope_ker % grval_ipoin(idime,ipoin)
                      end do
                   end do
                   do idime = 1,ndime
                      xvalu(idime) = xvalu(idime) / real(pnode,rp)
                   end do

                else

                   if( present(gpcar) ) then
                      do igaus = 1,qgaus
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            do idime = 1,ndime
                               xvalu(idime) = xvalu(idime) &
                                    + prope_ker % value_ipoin(ipoin) * gpcar(idime,inode,1)
                            end do
                         end do
                      end do
                      do idime = 1,ndime
                         xvalu(idime) = xvalu(idime) / real(qgaus,rp)
                      end do
                   else
                      call runend('MOD_KER_PROPER: GRADIENT NOT AVAILABLE AT COG')
                   end if

                end if

             else
                !
                ! Value
                !
                xvalu(1) = 0.0_rp
                do inode = 1,pnode
                   ipoin    = lnods(inode,ielem)
                   xvalu(1) = xvalu(1) + prope_ker % value_ipoin(ipoin)
                end do
                xvalu(1) = xvalu(1) / real(pnode,rp)
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 )  then
                xvalu(1:ndime) = 0.0_rp
             else
                !if( .not. associated(prope_ker) ) call runend('PROB1')
                !if( .not. associated(prope_ker % value_const) ) call runend('PROB2')
                xvalu(1) = prope_ker % value_const(imate)
             end if

          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 8. Property required on PNODE
       !
       !----------------------------------------------------------------------

       if( wherein == 'PNODE' ) then

          ielem = ielbo
          pelty = abs(ltype(ielem))
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)
          pgaus = ngaus(pelty)
          pnode = nnode(pelty)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( present(gpsha) ) then
                if( pgaus /= qgaus ) call runend('MOD_KER_PROPER: COG NOT CODED')
             end if

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do inode = 1,ndime*pnode
                   xvalu(inode) = 0.0_rp
                end do

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         kposi = (inode-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = xvalu(kposi) + prope_ker % grval_ielem(ielem) % a(idime,igaus) &
                                 * elmar(pelty) % shaga(igaus,inode)
                         end do
                      end do
                   end do

                else

                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      elpro(inode) = 0.0_rp
                      do igaus = 1,pgaus
                         elpro(inode) = elpro(inode) + prope_ker % value_ielem(ielem) % a(igaus) &
                              * elmar(pelty) % shaga(igaus,inode)
                      end do
                   end do
                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         do idime = 1,ndime
                            elgrp(idime,igaus) = elgrp(idime,igaus) + elpro(inode) * elmar(pelty) % shape(inode,igaus)
                         end do
                      end do
                   end do
                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         kposi = (inode-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = xvalu(kposi) + elgrp(idime,igaus) * elmar(pelty) % shaga(igaus,inode)
                         end do
                      end do
                   end do

                end if

             else
                !
                ! Value
                !
                do inode = 1,pnode
                   xvalu(inode) = 0.0_rp
                   do igaus = 1,pgaus
                      xvalu(inode) = xvalu(inode) + prope_ker % value_ielem(ielem) % a(igaus) &
                           * elmar(pelty) % shaga(igaus,inode)
                   end do
                end do

             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do inode = 1,ndime*pnode
                   xvalu(inode) = 0.0_rp
                end do
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   do inode = 1,pnode
                      kposi = (inode-1) * ndime
                      do idime = 1,ndime
                         kposi = kposi + 1
                         xvalu(kposi) = prope_ker % grval_ipoin(idime,ipoin)
                      end do
                   end do
                else
                   if( .not. present(gpcar) ) call runend('MOD_KER_PROPER: GPCAR MISSING')
                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         ipoin = lnods(inode,ielem)
                         do idime = 1,ndime
                            elgrp(idime,igaus) = elgrp(idime,igaus) + prope_ker % value_ipoin(ipoin) * gpcar(idime,inode,igaus)
                         end do
                      end do
                   end do
                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         kposi = (inode-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = xvalu(kposi) + elgrp(idime,igaus) * elmar(pelty) % shaga(igaus,inode)
                         end do
                      end do
                   end do
                end if
             else
                !
                ! Value
                !
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   xvalu(inode) = prope_ker % value_ipoin(ipoin)
                end do
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 ) then
                do inode = 1,ndime*pnode
                   xvalu(inode) = 0.0_rp
                end do
             else
                do inode = 1,pnode
                   xvalu(inode) = prope_ker % value_const(imate)
                end do
             end if

          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 9. Property required Anywhere
       !
       !----------------------------------------------------------------------

       if( wherein == 'ANYWH' ) then

          igaus = iposi
          ielem = ielbo
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             call runend('MOD_KER_PROPER: TEST THIS WITH LAGRANGE')

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             pelty = abs(ltype(ielem))
             pnode = lnnod(ielem)
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elpro(inode) = prope_ker % value_ipoin(ipoin)
             end do
             xvalu(1) = 0.0_rp
             if( present(gpsha) ) then
                kposi    = (igaus-1) * qnode
                do inode = 1,pnode
                   kposi    = kposi + 1
                   xvalu(1) = xvalu(1) + elpro(inode) * gpsha(kposi)
                end do
             else
                call runend('MOD_KER_PROPER: NOT CODED')
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             xvalu(1) = prope_ker % value_const(imate)

          end if

       end if

    end if
    !nullify(prope_ker) ! This avoids the bug of calling wrongly the ker_proper routine

  end subroutine ker_proper_scalar

  subroutine ker_proper_vector_0(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_DIM)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE,1,*)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE,1,1,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE,1,*)
    integer(ip),                     optional :: VECTOR_DIM
    integer(ip)                               :: VECTOR_SIZE_LOC

    if( present(VECTOR_DIM) ) then
       VECTOR_SIZE_LOC = VECTOR_DIM
    else
       VECTOR_SIZE_LOC = VECTOR_SIZE       
    end if

    call ker_proper_vector(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)
    
  end subroutine ker_proper_vector_0

  subroutine ker_proper_vector_1(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_DIM)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE,*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE,1,*)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE,1,1,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE,1,*)
    integer(ip),                     optional :: VECTOR_DIM
    integer(ip)                               :: VECTOR_SIZE_LOC

    if( present(VECTOR_DIM) ) then
       VECTOR_SIZE_LOC = VECTOR_DIM
    else
       VECTOR_SIZE_LOC = VECTOR_SIZE       
    end if
    
    call ker_proper_vector(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)
    
  end subroutine ker_proper_vector_1

  subroutine ker_proper_vector_2(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_DIM)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE,1,*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE,1,*)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE,1,1,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE,1,*)
    integer(ip),                     optional :: VECTOR_DIM
    integer(ip)                               :: VECTOR_SIZE_LOC
    
    if( present(VECTOR_DIM) ) then
       VECTOR_SIZE_LOC = VECTOR_DIM
    else
       VECTOR_SIZE_LOC = VECTOR_SIZE       
    end if
    
    call ker_proper_vector(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)
    
  end subroutine ker_proper_vector_2

  subroutine ker_proper_vector_3(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_DIM)

    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE,1,1,*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE,1,*)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE,1,1,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE,1,*)
    integer(ip),                     optional :: VECTOR_DIM
    integer(ip)                               :: VECTOR_SIZE_LOC
    
    if( present(VECTOR_DIM) ) then
       VECTOR_SIZE_LOC = VECTOR_DIM
    else
       VECTOR_SIZE_LOC = VECTOR_SIZE       
    end if

    call ker_proper_vector(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)
    
  end subroutine ker_proper_vector_3

  subroutine ker_proper_scalar_00(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    ! For example: only on one Gauss point (IGAUS)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu(*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in)           :: gpsha(*)
    real(rp),            intent(in), optional :: gpcar(1,*)
    real(rp),            intent(in), optional :: gpcor(*)
    call ker_proper_scalar(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
  end subroutine ker_proper_scalar_00

  subroutine ker_proper_scalar_000(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    ! For example: only on one Gauss point (IGAUS)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in)           :: gpsha(*)
    real(rp),            intent(in), optional :: gpcar(1,*)
    real(rp),            intent(in), optional :: gpcor(*)
    real(rp)                                  :: xval2(2)
    call ker_proper_scalar(&
       wname,wherein,iposi,ielbo,xval2,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    xvalu = xval2(1)
  end subroutine ker_proper_scalar_000
 
  subroutine ker_proper_scalar_0(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu(*)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(1,*)
    real(rp),            intent(in), optional :: gpcar(1,1,*)
    real(rp),            intent(in), optional :: gpcor(1,*)

    if(  .not. present(qnode) .and. &
         .not. present(qgaus) .and. &
         .not. present(gpsha) .and. &
         .not. present(gpcar) .and. &
         .not. present(gpcor) ) then
      call ker_proper_scalar(&
            wname,wherein,iposi,ielbo,xvalu)
    else
       call ker_proper_scalar(&
            wname,wherein,iposi,ielbo,xvalu,&
            qnode,qgaus,gpsha,gpcar,gpcor)
    end if
  end subroutine ker_proper_scalar_0

  subroutine ker_proper_scalar_1(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu(1,*)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(1,*)
    real(rp),            intent(in), optional :: gpcar(1,1,*)
    real(rp),            intent(in), optional :: gpcor(1,*)
    call ker_proper_scalar(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
  end subroutine ker_proper_scalar_1

  subroutine ker_proper_scalar_2(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu(1,1,*)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(1,*)
    real(rp),            intent(in), optional :: gpcar(1,1,*)
    real(rp),            intent(in), optional :: gpcor(1,*)
    call ker_proper_scalar(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
  end subroutine ker_proper_scalar_2

  subroutine ker_proper_vector(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)
    !-----------------------------------------------------------------------
    !****f* Kermod/ker_proper
    ! NAME
    !    ker_proper
    ! DESCRIPTION
    !    Update properties
    !
    !    Property can be defined with 3 formats
    !    It can be required at 8 different places
    !
    !       WHEREIN   IPOSI  IELBO
    !       --------------------
    !    1. IPOIN   IPOIN  NONSENSE
    !    2. NPOIN   DUMMI  Do smoothing => Loop over IELEM
    !    3. IGAUS   IGAUS  IELEM
    !    4. PGAUS   DUMMI  IELEM
    !    5. IGAUB   IGAUB  IBOUN        => IELEM = LBOEL(,IBOUN)
    !    6. PGAUB   DUMMI  IBOUN        => IELEM = LBOEL(,IBOUN)
    !    7. COG     DUMMI  IELEM
    !    8. PNODE   DUMMI  IELEM
    !
    !    GPSHA is dimensioned as a vector in order to use it as optional
    !    argument
    !
    ! USED BY
    !    many subroutines
    !***
    !-----------------------------------------------------------------------

    use def_kintyp,      only                 :  ip,rp
    use def_elmtyp,      only                 :  ELCUT
    use def_master,      only                 :  mem_modul,modul,nturb
    use def_domain,      only                 :  nmate,npoin,nelem,mnode,lelch
    use def_domain,      only                 :  lnods,lnnod,lmate,lmatn,lnnod
    use def_domain,      only                 :  ngaus,ltype,elmar,ndime,mgaus
    use def_domain,      only                 :  ltypb,nnode
    use def_kermod,      only                 :  typ_valpr_ker
    use def_kermod,      only                 :  densi_ker,visco_ker
    use def_kermod,      only                 :  poros_ker,condu_ker
    use def_kermod,      only                 :  sphea_ker,dummy_ker
    use def_kermod,      only                 :  turmu_ker,mixin_ker
    use mod_memory
    implicit none
    integer(ip),         intent(in)           :: VECTOR_SIZE_LOC
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE_LOC)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE_LOC,*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE_LOC,qnode)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE_LOC,ndime,mnode,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE_LOC,ndime,*)
    type(typ_valpr_ker), pointer              :: prope_ker
    integer(ip)                               :: imate
    integer(ip)                               :: ielem(VECTOR_SIZE_LOC)
    integer(ip)                               :: pgaus
    integer(ip)                               :: inode,ipoin,igaus,iturb,kelem
    integer(ip)                               :: ilaws,pelty,pnode
    integer(ip)                               :: iboun,ivect
    integer(ip)                               :: pgaub,kposi,kfl_gradi,idime,kfl_deriv,kfl_grder
    integer(ip)                               :: kfl_deriv_tur,kfl_deriv_vel,jposi
    integer(ip)                               :: lposi
    real(rp)                                  :: elpro(VECTOR_SIZE_LOC,mnode)
    real(rp)                                  :: elgrp(VECTOR_SIZE_LOC,ndime,max(mgaus,mnode))
    real(rp),   pointer                       :: auxva(:)

    nullify(auxva)

    !----------------------------------------------------------------------
    !
    ! What to compute
    !
    !----------------------------------------------------------------------

    if( wname(1:2) == 'GR' ) then
       kfl_gradi = 1
    else
       kfl_gradi = 0
    end if

    if( wname(1:2) == 'DR' ) then
       kfl_deriv = 1
    else
       kfl_deriv = 0
    end if

    if( wname(1:2) == 'GD' ) then
       kfl_grder = 1
    else
       kfl_grder = 0
    end if

    if( wname(1:3) == 'TDR' ) then
       kfl_deriv_tur = 1
    else
       kfl_deriv_tur = 0
    end if

    if( wname(1:3) == 'VDR' ) then
       kfl_deriv_vel = 1
    else
       kfl_deriv_vel = 0
    end if

    !----------------------------------------------------------------------
    !
    ! Which property to compute
    !
    !----------------------------------------------------------------------

    if(      wname == 'DENSI' .or. wname == 'GRDEN' .or. wname == 'DRDEN' .or. wname == 'GDDEN') then
       prope_ker => densi_ker
    else if( wname == 'VISCO' .or. wname == 'GRVIS' .or. wname == 'DRVIS' .or. wname == 'GDVIS') then
       prope_ker => visco_ker
    else if( wname == 'MIXIN' ) then
       prope_ker => mixin_ker
    else if( wname == 'POROS' .or. wname == 'GRPOS' ) then
       prope_ker => poros_ker
    else if( wname == 'CONDU' .or. wname == 'GRCON' ) then
       prope_ker => condu_ker
    else if( wname == 'SPHEA' .or. wname == 'GRSPE' ) then
       prope_ker => sphea_ker
    else if( wname == 'DUMMY' .or. wname == 'GRDUM' ) then
       prope_ker => dummy_ker
    else if( wname == 'TURBU' .or. wname == 'TURVI' .or. wname == 'GRTUR' .or. wname == 'TDRTU' .or. wname == 'VDRTU' ) then
       prope_ker => turmu_ker
    else if( wname == 'ABSOR' ) then
       prope_ker => absor_ker
    else if( wname == 'SCATT' ) then
       prope_ker => scatt_ker
    else
       print *,'--- I DONT UNDERSTAND ',wname
       call runend('KER_PROPER CALLED WITH WRONG PROPERTY') ! This avoids the bug of calling wrongly the ker_proper routine
    end if

    if( prope_ker % kfl_exist == 0 ) then
       !
       ! Put value to zero
       !
       if( wherein == 'IPOIN' .or. wherein == 'IGAUS' .or. wherein == 'IGAUB' .or. wherein == 'COG' ) then
          kposi = 1
       else if( wherein == 'NPOIN' ) then
          kposi = npoin
       else if( wherein == 'PNODE' ) then
          kposi = lnnod(ielbo(1))
       else if( wherein == 'PGAUS' ) then
          if( present(gpsha) ) then
             kposi = qgaus
          else
             kposi = ngaus(abs(ltype(ielbo(1))))
          end if
       else
          call runend('KER_PROPER_VECTOR: NOT CODED 2')
       end if
       if( kfl_gradi == 1 ) kposi = kposi * ndime
       do kelem = 1,kposi
          xvalu(1:VECTOR_SIZE_LOC,kelem) = 0.0_rp
       end do

    else

       if( wherein == 'IGAUS' ) then

          !----------------------------------------------------------------------
          !
          ! 3. Property required on IGAUS
          !    Note: it has not been validated 03/03/2017
          !----------------------------------------------------------------------

          igaus = iposi
          ielem = ielbo(1)
          imate = lmate(ielem(1))
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   do ivect = 1,VECTOR_SIZE_LOC
                      if( ielem(ivect) > 0 ) then
                         do idime = 1,ndime
                            xvalu(ivect,idime) = prope_ker % grval_ielem(ielem(ivect)) % a(idime,igaus)
                         end do
                      end if
                   end do
                else
                   xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp
                end if

             else
                !
                ! Value
                !
                if( present(gpsha) ) then
                   do ivect = 1,VECTOR_SIZE_LOC
                      kelem = ielem(ivect)
                      if( kelem > 0 ) then
                         pelty = abs(kelem)
                         pnode = lnnod(kelem)
                         pgaus = ngaus(pelty)
                         xvalu(ivect,1) = 0.0_rp
                         do inode = 1,pnode
                            do igaus = 1,pgaus
                               xvalu(ivect,1) = xvalu(ivect,1) + gpsha(ivect,inode) * elmar(pelty) % shaga(igaus,inode) &
                                    * prope_ker % value_ielem(kelem) % a(igaus)
                            end do
                         end do
                      end if
                   end do
                else
                   do ivect = 1,VECTOR_SIZE_LOC
                      if( ielem(ivect) > 0 ) then
                         xvalu(ivect,1) = prope_ker % value_ielem(ielem(ivect)) % a(igaus)
                      end if
                   end do
                end if

             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             call runend('KER_PROPER: NOT AVAILABLE IGAUS 2')

!!$             pelty = abs(ltype(ielem))
!!$             pnode = lnnod(ielem)
!!$             do inode = 1,pnode
!!$                ipoin = lnods(inode,ielem)
!!$                elpro(inode) = prope_ker % value_ipoin(ipoin)
!!$             end do
!!$
!!$             if( kfl_gradi == 1 ) then
!!$                !
!!$                ! Gradient
!!$                !
!!$                xvalu(1:ndime) = 0.0_rp
!!$
!!$                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
!!$
!!$                   do inode = 1,pnode
!!$                      ipoin = lnods(inode,ielem)
!!$                      do idime = 1,ndime
!!$                         elgrp(idime,inode) = prope_ker % grval_ipoin(idime,ipoin)
!!$                      end do
!!$                   end do
!!$                   if( present(gpsha) ) then
!!$                      kposi = (igaus-1) * qnode
!!$                      do inode = 1,pnode
!!$                         kposi = kposi + 1
!!$                         do idime = 1,ndime
!!$                            xvalu(idime) = xvalu(idime) + elgrp(idime,inode) * gpsha(kposi)
!!$                         end do
!!$                      end do
!!$                   else
!!$                      do inode = 1,pnode
!!$                         do idime = 1,ndime
!!$                            xvalu(idime) = xvalu(idime) + elgrp(idime,inode) * elmar(pelty) % shape(idime,igaus)
!!$                         end do
!!$                      end do
!!$                   end if
!!$
!!$                else
!!$
!!$                   if( present(gpcar) ) then
!!$                      do inode = 1,pnode
!!$                         do idime = 1,ndime
!!$                            xvalu(idime) = xvalu(idime) + elpro(inode) * gpcar(idime,inode,1)
!!$                         end do
!!$                      end do
!!$                   else
!!$                      call runend('MOD_KER_PROPER: GPCAR NEEDED')
!!$                   end if
!!$                end if
!!$
!!$          else
!!$             !
!!$             ! Value
!!$             !
!!$             call runend('KER_PROPER: NOT AVAILABLE IGAUS 3')
!!$
!!$                xvalu(1) = 0.0_rp
!!$                if( present(gpsha) ) then
!!$                   kposi    = (igaus-1) * qnode
!!$                   do inode = 1,pnode
!!$                      kposi    = kposi + 1
!!$                      xvalu(1) = xvalu(1) + elpro(inode) * gpsha(kposi)
!!$                   end do
!!$                else
!!$                   do inode = 1,pnode
!!$                      xvalu(1) = xvalu(1) + elpro(inode) * elmar(pelty) % shape(inode,igaus)
!!$                   end do
!!$                end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 ) then
                xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp
             else
                xvalu(1:VECTOR_SIZE_LOC,1) = prope_ker % value_const(imate)
             end if

          end if

       else if( wherein == 'PGAUS' ) then

          !----------------------------------------------------------------------
          !
          ! 4. Property required on IGAUS=1,PGAUS
          !    Note: it has been validated 03/03/2017
          !----------------------------------------------------------------------

          ielem(1:VECTOR_SIZE_LOC) = ielbo(1:VECTOR_SIZE_LOC)
          imate                = lmate(ielem(1))
          ilaws                = prope_ker % ilaws(imate)
          pelty                = abs(ltype(ielem(1)))
          pgaus                = ngaus(pelty)
          if( lelch(ielem(1)) == ELCUT ) pgaus = lgaus(ielem(1))

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
             !
             ! Property is computed element-wise
             !
             if( kfl_gradi == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus .and. kfl_cutel == 0 ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            jposi = (igaus - 1) * ndime
                            do idime = 1,ndime
                               kposi = jposi + idime
                               xvalu(ivect,kposi) = prope_ker % grval_ielem(ielem(ivect)) % a(idime,igaus)
                            end do
                         end do
                      end do
                   else
                      if( qgaus /= pgaus .and. kfl_cutel == 0 ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            jposi = (igaus - 1) * ndime
                            do idime = 1,ndime
                               kposi = jposi + idime
                               xvalu(ivect,kposi) = prope_ker % grval_ielem(ielem(ivect)) % a(idime,igaus)
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*ndime
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv == 1 ) then

                call runend('KER_PROPER_VECTOR: NOT CODED 1')

                if( prope_ker % llaws(ilaws) % kfl_deriv == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            kposi = (igaus - 1) * qnode
                            do inode = 1,qnode
                               kposi = kposi + 1
                               xvalu(ivect,kposi) = prope_ker % drval_ielem(ielem(ivect)) % a(inode,igaus)
                            end do
                         end do
                      end do
                   else
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,pgaus
                            kposi = (igaus - 1) * qnode
                            do inode = 1,qnode
                               kposi = kposi + 1
                               xvalu(ivect,kposi) = prope_ker % drval_ielem(ielem(ivect)) % a(inode,igaus)
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_grder == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_grder == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            kposi = (igaus - 1)* qnode*ndime
                            do inode = 1,qnode
                               do idime = 1, ndime
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % gdval_ielem(ielem(ivect)) % a(idime,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   else
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,pgaus
                            do inode = 1,qnode
                               kposi = (igaus - 1)*(inode-1) * qnode*ndime
                               do idime = 1, ndime
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % gdval_ielem(ielem(ivect)) % a(idime,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*ndime
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv_tur == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_deriv_tur == 1 ) then

                   if( present(gpsha) ) then
                      nturb = 2
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            kposi = (igaus - 1)* qnode*ndime
                            do inode = 1,qnode
                               do iturb = 1, nturb
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % drval_tur_ielem(ielem(ivect)) % a(iturb,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   else
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,pgaus
                            do inode = 1,qnode
                               kposi = (igaus - 1)*(inode-1) * qnode*ndime
                               do iturb = 1, nturb
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % drval_tur_ielem(ielem(ivect)) % a(iturb,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*nturb
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv_vel == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_deriv_vel == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            kposi = (igaus - 1)* qnode*ndime
                            do inode = 1,qnode
                               do idime = 1, ndime
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % drval_vel_ielem(ielem(ivect)) % a(idime,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   else
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,pgaus
                            do inode = 1,qnode
                               kposi = (igaus - 1)*(inode-1) * qnode*ndime
                               do idime = 1, ndime
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % drval_vel_ielem(ielem(ivect)) % a(idime,inode,igaus)
                               enddo
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*ndime
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else


                if( present(gpsha) .and. kfl_cutel == 0 ) then
                   if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                   do ivect = 1,VECTOR_SIZE_LOC
                      do igaus = 1,qgaus
                         xvalu(ivect,igaus) = prope_ker % value_ielem(ielem(ivect)) % a(igaus)
                      end do
                   end do
                else
                   do ivect = 1,VECTOR_SIZE_LOC
                      do igaus = 1,pgaus
                         xvalu(ivect,igaus) = prope_ker % value_ielem(ielem(ivect)) % a(igaus)
                      end do
                   end do
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then
             !
             ! Property is computed node-wise
             ! Note: it has not been validated 03/03/2017
             do ivect = 1,VECTOR_SIZE_LOC
                pelty = abs(ltype(ielem(ivect)))
                pnode = nnode(pelty)
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem(ivect))
                   elpro(ivect,inode) = prope_ker % value_ipoin(ipoin)
                end do
             end do

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do kposi = 1,pgaus*ndime
                   xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                end do

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do ivect = 1,VECTOR_SIZE_LOC
                      do inode = 1,pnode
                         ipoin = lnods(inode,ielem(ivect))
                         do idime = 1,ndime
                            elgrp(ivect,idime,inode) = prope_ker % grval_ipoin(idime,ipoin)
                         end do
                      end do
                   end do

                   if( present(gpsha) ) then
                      do igaus = 1,qgaus
                         kposi = (igaus-1) * qnode
                         do inode = 1,pnode
                            kposi = kposi + 1
                            lposi = (igaus-1) * ndime
                            do idime = 1,ndime
                               lposi = lposi + 1
                               xvalu(1:VECTOR_SIZE_LOC,lposi) = xvalu(1:VECTOR_SIZE_LOC,lposi) + elgrp(1:VECTOR_SIZE_LOC,idime,inode) * gpsha(1:VECTOR_SIZE_LOC,kposi)
                            end do
                         end do
                      end do
                   else
                      do igaus = 1,qgaus
                         do inode = 1,pnode
                            lposi = (igaus-1) * ndime
                            do idime = 1,ndime
                               lposi = lposi + 1
                               xvalu(1:VECTOR_SIZE_LOC,lposi) = xvalu(1:VECTOR_SIZE_LOC,lposi) + elgrp(1:VECTOR_SIZE_LOC,idime,inode) * elmar(pelty) % shape(inode,igaus)
                            end do
                         end do
                      end do
                   end if

                else

                   if( present(gpcar) ) then
                      do igaus = 1,qgaus
                         kposi = (igaus-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            do inode = 1,pnode
                               xvalu(1:VECTOR_SIZE_LOC,kposi) = xvalu(1:VECTOR_SIZE_LOC,kposi) + elpro(1:VECTOR_SIZE_LOC,inode) * gpcar(1:VECTOR_SIZE_LOC,idime,inode,igaus)
                            end do
                         end do
                      end do
                   end if

                end if

             else
                !
                ! Value
                !
                xvalu(1:VECTOR_SIZE_LOC,1:qgaus) = 0.0_rp

                if( present(gpsha) ) then
                   do igaus = 1,qgaus
                      kposi        = (igaus-1) * qnode
                      do inode = 1,pnode
                         kposi = kposi + 1
                         xvalu(1:VECTOR_SIZE_LOC,igaus) = xvalu(1:VECTOR_SIZE_LOC,igaus) + elpro(1:VECTOR_SIZE_LOC,inode) * gpsha(1:VECTOR_SIZE_LOC,kposi)
                      end do
                   end do
                else
                   do igaus = 1,ngaus(pelty)
                      do inode = 1,pnode
                         xvalu(1:VECTOR_SIZE_LOC,igaus) = xvalu(1:VECTOR_SIZE_LOC,igaus) + elpro(1:VECTOR_SIZE_LOC,inode) * elmar(pelty) % shape(inode,igaus)
                      end do
                   end do
                end if

             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then
             !
             ! Property is constant
             !
             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do kposi = 1,pgaus*ndime
                   xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                end do

             else
                !
                ! Value
                !
                if( present(gpsha) ) then
                   do igaus = 1,qgaus
                      xvalu(1:VECTOR_SIZE_LOC,igaus) = prope_ker % value_const(imate)
                   end do
                else
                   do igaus = 1,ngaus(abs(ltype(ielem(1))))
                      xvalu(1:VECTOR_SIZE_LOC,igaus) = prope_ker % value_const(imate)
                   end do
                end if

             end if

          end if

          !----------------------------------------------------------------------
          !
          ! 7. Property required on C.O.G.
          !    Note: it has not been validated 03/03/2017
          !----------------------------------------------------------------------

       else if( wherein(1:3) == 'COG' ) then

          ielem(1:VECTOR_SIZE_LOC) = ielbo(1:VECTOR_SIZE_LOC)
          imate = lmate(ielem(1))
          ilaws = prope_ker % ilaws(imate)
          pelty = abs(ltype(ielem(1)))
          pgaus = ngaus(pelty)
          pnode = nnode(pelty)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 )  then
                !
                ! Gradient
                !
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp
                   do ivect = 1,VECTOR_SIZE_LOC
                      if( ielem(ivect) > 0 ) then
                         do igaus = 1,pgaus
                            do idime = 1,ndime
                               xvalu(ivect,idime) = xvalu(ivect,idime) + prope_ker % grval_ielem(ielem(ivect)) % a(idime,igaus)
                            end do
                         end do
                      end if
                   end do
                   do idime = 1,ndime
                      do ivect = 1,VECTOR_SIZE_LOC
                         xvalu(ivect,idime) = xvalu(ivect,idime) / real(pgaus,rp)
                      end do
                   end do
                else
                   call runend('MOD_KER_PROPER: GRDAIENT ON COG NOT CODED')
                end if

             else
                !
                ! Value
                !
                xvalu(1:VECTOR_SIZE_LOC,1) = 0.0_rp
                if( present(gpsha) ) then
                   do ivect = 1,VECTOR_SIZE_LOC
                      do igaus = 1,qgaus
                         xvalu(ivect,1) = xvalu(ivect,1) + prope_ker % value_ielem(abs(ielem(1))) % a(igaus)
                      end do
                   end do
                   xvalu(1:VECTOR_SIZE_LOC,1) = xvalu(1:VECTOR_SIZE_LOC,1) / real(qgaus,rp)
                else
                   do ivect = 1,VECTOR_SIZE_LOC
                      do igaus = 1,pgaus
                         xvalu(ivect,1) = xvalu(ivect,1) + prope_ker % value_ielem(abs(ielem(1))) % a(igaus)
                      end do
                   end do
                   xvalu(1:VECTOR_SIZE_LOC,1) = xvalu(1:VECTOR_SIZE_LOC,1) / real(pgaus,rp)
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             if( kfl_gradi == 1 )  then
                !
                ! Gradient
                !
                xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do ivect = 1,VECTOR_SIZE_LOC
                      if( ielem(ivect) > 0 ) then
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem(ivect))
                            do idime = 1,ndime
                               xvalu(ivect,idime) = xvalu(ivect,idime) + prope_ker % grval_ipoin(idime,ipoin)
                            end do
                         end do
                      end if
                   end do
                   do idime = 1,ndime
                      xvalu(1:VECTOR_SIZE_LOC,idime) = xvalu(1:VECTOR_SIZE_LOC,idime) / real(pnode,rp)
                   end do

                else

                   if( present(gpcar) ) then
                      do ivect = 1,VECTOR_SIZE_LOC
                         if( ielem(ivect) > 0 ) then
                            do igaus = 1,qgaus
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem(ivect))
                                  do idime = 1,ndime
                                     xvalu(ivect,idime) = xvalu(ivect,idime) &
                                          + prope_ker % value_ipoin(ipoin) * gpcar(ivect,idime,inode,1)
                                  end do
                               end do
                            end do
                         end if
                      end do
                      do idime = 1,ndime
                         xvalu(1:VECTOR_SIZE_LOC,idime) = xvalu(1:VECTOR_SIZE_LOC,idime) / real(qgaus,rp)
                      end do
                   else
                      call runend('MOD_KER_PROPER: GRADIENT NOT AVAILABLE AT COG')
                   end if

                end if

             else
                !
                ! Value
                !
                xvalu(1:VECTOR_SIZE_LOC,1) = 0.0_rp
                do ivect = 1,VECTOR_SIZE_LOC
                   if( ielem(ivect) > 0 ) then
                      do inode = 1,pnode
                         ipoin          = lnods(inode,ielem(ivect))
                         xvalu(ivect,1) = xvalu(ivect,1) + prope_ker % value_ipoin(ipoin)
                      end do
                   end if
                end do
                xvalu(1:VECTOR_SIZE_LOC,1) = xvalu(1:VECTOR_SIZE_LOC,1) / real(pnode,rp)
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 )  then
                xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp
             else
                xvalu(1:VECTOR_SIZE_LOC,1) = prope_ker % value_const(imate)
             end if

          end if

       else

          call runend('KER_PROPER_VECTOR: NOT CODED 3')

       end if

    end if
    !
    ! This fails with OpenMP
    !nullify(prope_ker) ! This avoids the bug of calling wrongly the ker_proper routine

  end subroutine ker_proper_vector

  subroutine ker_updpro(itask_xxx)

    !-----------------------------------------------------------------------
    !****f* Kermod/ker_updpro
    ! NAME
    !    ker_updpro
    ! DESCRIPTION
    !    Update properties
    !-----------------------------------------------------------------------

    use def_domain
    use def_master
    use def_kermod    ! wallcoupling_extr_boun, ptb_to_use, is_interior and others
    use def_parame, only :  pi,zero,one
    use def_elmtyp, only :  ELCUT
    use def_elmtyp, only :  TET04
    use def_elmtyp, only :  TRI03
    use mod_ker_regularization, only : regul_k, regul_e, kfl_regularization
    use mod_ker_ILSA_sgsmodel, only : ker_ILSA_sgs_viscosity
    use mod_memory,         only :  memory_alloca, memory_deallo
    use mod_interpolation,  only : COU_GET_INTERPOLATE_POINTS_VALUES

    implicit none
    integer(ip), intent(in), optional  :: itask_xxx
    integer(ip), save            ::  ipass(100) = 0
    integer(ip)                  :: itask, kfl_force
    integer(ip)                  :: imate,pnode,ielem,ipoin,kpoin
    integer(ip)                  :: pgaus,inode,pelty,inodb,pblty
    integer(ip)                  :: iprop,ilaws,pnodb,iboun,pgaub,igaub
    integer(ip)                  :: iresp,kfl_updpr(10),idime,igaus
    integer(ip)                  :: ispec,jspec,jdime
    real(rp)                     :: gpsha(mnode,mgaus)
    real(rp)                     :: gpder(ndime,mnode,mgaus)
    real(rp)                     :: gpcar(ndime,mnode,mgaus)
    real(rp)                     :: gphes(ntens,mnode,mgaus)
    real(rp)                     :: gpvol(mgaus)
    real(rp)                     :: dummr(ndime,mgaus)
    real(rp)                     :: elcod(ndime,mnode)
    real(rp)                     :: elvel(ndime,mnode), eltur(2, mnode), elwal(mnode)
    real(rp)                     :: gpcor(ndime,mgaus)
    real(rp)                     :: gp_temperature(mgaus)
    real(rp)                     :: gp_grad_tem(ndime,mgaus), gp_grad_wmean(ndime,mgaus)
    real(rp)                     :: gpvel(ndime,mgaus)
    real(rp)                     :: gptan(ndime,mgaus),gbtan(ndime)
    real(rp)                     :: gpgve(ndime,ndime,mgaus)
    real(rp)                     :: gptur

    real(rp)                     :: ellev(mnode),eltem(mnode),elwmean(mnode)
    real(rp)                     :: prope_wat,prope_air
    real(rp)                     :: hleng(3),tragl(ndime,ndime)

    real(rp)                     :: eps,ooeps,oovpi,f,phi,kap
    real(rp)                     :: densi_phase(2),z0,rho,ustar

    real(rp)                     :: T,T0,mu0,mu,S0,const, heigh, cd, cdlad, auxva , lad
    real(rp)                     :: lscale,tscale

    real(rp)                     :: phiij,sumphi(nspec),gptem,gbtem,gbwme,gphog(mgaus),zmaxi,zmaxf,gpcah(mgaus)
    real(rp)                     :: lamax, zz, n, shaaa, hoveg, gpwmean, gpden(mgaus), gpvis(1), gpnut, gpmut,gpnu,X,cv1,fv1
    real(rp)                     :: F2, kinen, gpwal, omega, gpvor, a1, as, epsil, relat, gbwal, gbden(1), gbvis(1),velmo,tanmo
    real(rp)                     :: dgpkin(mnode),dgpome(mnode), arg2,darg2_dkin(mnode),darg2_dome(mnode), dgpnut_du(mnode,ndime)
    real(rp)                     :: gpvor_mat(ndime,ndime), gpnuILSA(mgaus)
    real(rp)                     :: dF2_dkin(mnode), dF2_dome(mnode), dgpnut_dkin(mnode), dgpnut_dome(mnode), dgpvor_du(mnode,ndime)
    real(rp)                     :: dX_dnut,dfv1_dnut,dgpmut_dnut(mnode)
    real (rp),save               :: integ

    real (rp),parameter          :: h_can_critical=7.0_rp ! beware if this value is changed it should also be changed in Windmesh

    integer(ip)                  :: which_time, icoun, ncoun, kdime
    ! ke models
    real(rp)                     :: A0, As_r , W, phi_r, cmu, divve, simgr(3,3),  Wsqr6 
    real(rp)                     :: uaste,seci4, Cr,f0, f02
    real(rp)                     :: regue, reguk, sigmr

    type(typ_valpr_ker), pointer :: prope_ker
    type(lis_prope_ker), pointer :: listp_ker(:)

    real(rp), pointer :: conce_save(:,:,:)  !!! value to save conce for adjoint
    real(rp), pointer :: tempe_save(:,:)    !!! value to save tempe for adjoint
    real(rp), pointer :: untur_save(:,:,:)  !!! value to save untur for adjoint
    real(rp), pointer :: veloc_save(:,:,:)  !!! value to save veloc for adjoint

    !
    ! for shaw canopy model
    !
    integer(ip),parameter        :: ncoef_shaw=11
    real(rp)                     :: value,xcoef(2,ncoef_shaw),aux_dvari
    !
    ! for mixing length model 
    !
    real(rp)                     :: gbvel(ndime),gbgve(ndime,ndime),eltan(ndime,mnode)
    real(rp)    , pointer        :: tau_wall_boun(:,:)

    real(rp)                     :: gpmve,darcy,forch


    if( IMASTER .or. kfl_prope == 0 ) return

    nullify(prope_ker)
    nullify(listp_ker)

    if (present(itask_xxx)) then
       itask=itask_xxx
       kfl_force=0
    else
       itask=-1000
       kfl_force=1
    endif

    if (kfl_reset == 1) then
       which_time=3
    else
       which_time=1
    endif

    !----------------------------------------------------------------------
    !
    ! Assign forward values for adjoint case to temper and conce
    !
    !----------------------------------------------------------------------

    if( kfl_adj_prob == 1_ip) then
       conce_save => conce
       tempe_save => tempe
       untur_save => untur
       veloc_save => veloc
       conce      => conce_forw
       tempe      => tempe_forw
       untur      => untur_forw
       veloc      => veloc_forw
    end if

    allocate(listp_ker(mprop_ker))

    listp_ker( 1) % prop => densi_ker
    listp_ker( 2) % prop => visco_ker
    listp_ker( 3) % prop => poros_ker
    listp_ker( 4) % prop => condu_ker
    listp_ker( 5) % prop => sphea_ker
    listp_ker( 6) % prop => dummy_ker
    listp_ker( 7) % prop => turmu_ker
    listp_ker( 8) % prop => absor_ker
    listp_ker( 9) % prop => scatt_ker
    listp_ker(10) % prop => mixin_ker

    !----------------------------------------------------------------------
    !
    ! Constant properties: compute only at beginning
    !
    !----------------------------------------------------------------------

    if( itask == ITASK_INIUNK ) then
       do iprop = 1,mprop_ker
          prope_ker =>  listp_ker(iprop) % prop
          if( prope_ker % kfl_exist == 1 ) then
             do imate = 1,nmate
                if( prope_ker % wlaws(imate) == 'CONST' ) &
                     prope_ker % value_const(imate) = prope_ker % rlaws(1,imate)
             end do
          end if
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! Check if a property depends on the last module solved: LASTM_KER
    !
    !----------------------------------------------------------------------

    do iprop = 1,mprop_ker
       prope_ker => listp_ker(iprop)%prop
       kfl_updpr(iprop) = 0
       if( prope_ker % kfl_exist == 1 ) then
          do imate = 1,nmate
             if( prope_ker % wlaws(imate) /= 'CONST' ) then
                ilaws = prope_ker % ilaws(imate)
                do iresp = 1,mresp_ker
                   if( prope_ker % llaws(ilaws) % lresp(iresp) == lastm_ker&
                        .and.(prope_ker % update(1, imate)==itask&
                        .or.prope_ker % update(2, imate)==itask)) then
                      kfl_updpr(iprop) = 1
                   else if ( prope_ker % llaws(ilaws) % lresp(iresp) == -1 ) then
                      kfl_updpr(iprop) = 1
                   end if
                end do
             end if
          end do
          if ( ( kfl_force == 1 ) .or. ( itask == ITASK_INIUNK ) )  then
             kfl_updpr(iprop) = 1
          endif
       end if
    end do

    !----------------------------------------------------------------------
    !
    ! Variable properties
    !
    !----------------------------------------------------------------------

    do iprop = 1,mprop_ker

       if( kfl_updpr(iprop) == 1 ) then

          prope_ker => listp_ker(iprop) % prop

          if( prope_ker % kfl_exist == 1 ) then
                          
             do imate = 1,nmate

                if( itask == ITASK_INIUNK .and. prope_ker % value_default(imate) > 0.0_rp ) then

                   !----------------------------------------------------------------------
                   !
                   ! Run is starting (INIUNK): impose a default value if it exists (> 0)
                   !
                   !----------------------------------------------------------------------

                   ilaws = prope_ker % ilaws(imate)

                   if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

                      do ipoin = 1,npoin
                         prope_ker % value_ipoin(ipoin) = prope_ker % value_default(imate)
                      end do
                      if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                         do ipoin = 1,npoin
                            prope_ker % grval_ipoin(:,ipoin) = 0.0_rp
                         end do
                      end if

                   else if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

                      do ielem = 1,nelem
                         prope_ker % value_ielem(ielem) % a = prope_ker % value_default(imate)
                      end do
                      do iboun = 1,nboun
                         prope_ker % value_iboun(iboun) % a = prope_ker % value_default(imate)
                      end do

                      if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                         do ielem = 1,nelem
                            prope_ker % grval_ielem(ielem) % a = 0.0_rp
                         end do
                      end if

                   else if( prope_ker % wlaws(imate) == 'CONST' ) then

                      continue ! Nothing to do, it was already imposed
                      !prope_ker % value_const(imate) = prope_ker % value_default(imate)

                   end if

                else if( prope_ker % wlaws(imate) == 'IMMER' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Immersed boundary mixing
                   !
                   !----------------------------------------------------------------------
                   
                   if( iprop == 10 ) then

                      prope_wat = prope_ker % rlaws(1,imate)
                      prope_air = prope_ker % rlaws(2,imate)

                      do ielem = 1,nelem
                         pelty = ltype(ielem)
                         if( lmate(ielem) == imate .and. pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               ellev(inode) = fleve(ipoin,which_time)
                               elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                            end do
                            call elmca2(&
                                 pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                                 elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                                 gpder,gpcar,gphes,ielem)
                            call ker_biflui(&
                                 3_ip,pgaus,pnode,gpsha,gpcar,ellev,prope_wat,prope_air, &
                                 thicl,hleng,prope_ker % value_ielem(ielem) % a,         &
                                 dummr)            
                         end if
                      end do
                      
                      prope_ker % kfl_nedsm = 1

                   else

                      call runend('THIS PROPERTY DOES NOT MAKE SENSE FOR IMMERSED BOUNDARY')
                      
                   end if

                else if( prope_ker % wlaws(imate) == 'BIFLU' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Bifluid model
                   !
                   !----------------------------------------------------------------------

                   prope_wat = prope_ker % rlaws(1,imate)
                   prope_air = prope_ker % rlaws(2,imate)

                   do ielem = 1,nelem
                      pelty = ltype(ielem)
                      if( lmate(ielem) == imate .and. pelty > 0 ) then
                         pgaus = ngaus(pelty)
                         pnode = lnnod(ielem)
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            ellev(inode) = fleve(ipoin,which_time)
                            elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                         end do
                         call elmca2(&
                              pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                              elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                              gpder,gpcar,gphes,ielem)
                         call ker_biflui(&
                              1_ip,pgaus,pnode,gpsha,gpcar,ellev,prope_wat,prope_air, &
                              thicl,hleng,prope_ker % value_ielem(ielem) % a,         &
                              dummr)
                      end if
                   end do

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               ipoin = lnodb(inodb,iboun)
                               ellev(inodb) = fleve(ipoin,which_time)
                            end do
                            call ker_biflui(&
                                 2_ip,pgaub,pnodb,elmar(pblty) % shape,gpcar,ellev,prope_wat,prope_air,&
                                 thicl,hleng,prope_ker % value_iboun(iboun) % a,dummr)
                         end if
                      end if
                   end do
                   prope_ker % kfl_nedsm = 1

                else if( prope_ker % wlaws(imate) == 'BIFL2' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Bifluid model - including gradients
                   !
                   !----------------------------------------------------------------------

                   prope_wat = prope_ker % rlaws(1,imate)
                   prope_air = prope_ker % rlaws(2,imate)
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               ellev(inode) = fleve(ipoin,which_time)
                               elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                            end do
                            call elmca2(&
                                 pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                                 elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                                 gpder,gpcar,gphes,ielem)
                            call ker_biflu2(&
                                 1_ip,pgaus,pnode,gpsha,gpcar,ellev,prope_wat,prope_air, &
                                 thicl,hleng,prope_ker % value_ielem(ielem) % a,         &
                                 prope_ker % grval_ielem(ielem) % a)
                         end if
                      end if
                   end do

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               ipoin = lnodb(inodb,iboun)
                               ellev(inodb) = fleve(ipoin,which_time)
                            end do
                            call ker_biflu2(&
                                 2_ip,pgaub,pnodb,elmar(pblty) % shape,gpcar,ellev,prope_wat,prope_air,&
                                 thicl,hleng,prope_ker % value_iboun(iboun) % a,dummr)
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                else if( prope_ker % wlaws(imate) == 'SUTHE' ) then

                   !----------------------------------------------------------------------
                   !
                   !                                        T0 + S
                   ! Sutherland's law: mu/mu_0 = (T/T0)^1.5 --------
                   !                                        T  + S
                   !
                   !----------------------------------------------------------------------

                   mu0   = prope_ker % rlaws(1,imate)
                   T0    = prope_ker % rlaws(2,imate)
                   S0    = prope_ker % rlaws(3,imate)
                   if (S0.lt. epsilon(1.0_rp)) S0  = 110.0_rp
                   if (T0.lt. epsilon(1.0_rp)) T0  = 298.0_rp
                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time)
                      prope_ker % value_ipoin(ipoin) = mu0 * (T/T0)**1.5_rp * (T0+S0) / (T+S0)
                   end do
                else if( prope_ker % wlaws(imate) == 'CANOP' ) then
                   !----------------------------------------------------------------------
                   !
                   ! Canopy model
                   !
                   !----------------------------------------------------------------------
                   cd    = prope_ker % rlaws(2,imate) ! cd (generally =0.2)
                   lad   = prope_ker % rlaws(3,imate) ! leaf area density
                   zmaxf = prope_ker % rlaws(4,imate) ! zmaxf = zmax/heigh

                   if ( ipass(imate)==0 .and. kfl_canhe == 0) then   ! When a constant canopy height is used
                      do ielem = 1,nelem
                         if( lmate(ielem) == imate ) then
                            pelty = ltype(ielem)
                            if( pelty > 0 ) then
                               pnode = lnnod(ielem)
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  canhe(ipoin) = prope_ker % rlaws(1,imate)
                               end do
                            end if
                         end if
                      end do
                   end if

                   if ( ipass(imate)==0.and.zmaxf.gt.0.01_rp) then ! variable canopy density with height with a prescribed vertical density function
                      !
                      ! Obtain the integral for a height = 1   then we will just multiply by the real height
                      !
                      ncoun =1000
                      integ = 0.0_rp
                      heigh = 1.0_rp
                      zmaxi = zmaxf*heigh

                      do icoun = 1, ncoun -1
                         zz = real(icoun,rp)/real(ncoun,rp)*heigh
                         if (zz.lt.zmaxi) then
                            n = 6.0_rp
                         else
                            n = 0.5_rp
                         end if
                         integ = integ + (((heigh -zmaxi) /(heigh - zz))**n )*exp(n*(1.0_rp-(heigh -zmaxi) /(heigh - zz)))
                      end do

                      integ = integ  + 0.5_rp*(((heigh -zmaxi) /(heigh ))**6.0_rp )*exp(6.0_rp*(1.0_rp-(heigh -zmaxi) /(heigh )))

                      integ = integ  + 0.5_rp*(((heigh -zmaxi) /(1.0e-4_rp))**0.5_rp )*exp(0.5_rp*(1.0_rp-(heigh -zmaxi) /(1.0e-4_rp)))


                      integ = integ * heigh / real(ncoun,rp)
                   end if
                   ipass(imate) = 1

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            gpvel = 0.0_rp
                            gphog = 0.0_rp
                            gpcah = 0.0_rp
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)

                               do idime = 1,ndime
                                  elvel(idime, inode) = veloc(idime, ipoin,1)
                                  do igaus =1, pgaus
                                     gpvel(idime, igaus) = gpvel(idime,igaus) + &
                                          elvel(idime, inode) * elmar(pelty) % shape(inode,igaus)
                                  end do
                               end do
                               do igaus =1, pgaus
                                  shaaa = elmar(pelty) % shape(inode,igaus)
                                  gpcah(igaus) = gpcah(igaus) + canhe(ipoin)* shaaa
                                  gphog(igaus) = gphog(igaus) + heiov(ipoin)* shaaa
                               end do
                            end do


                            do igaus =1, pgaus
                               heigh = gpcah(igaus)
                               zmaxi = zmaxf*heigh
                               hoveg = gphog(igaus)
                               if ((hoveg.lt.heigh).and.(heigh>h_can_critical)) then ! below forest height and over the critical heig
                                  auxva =0.0_rp
                                  do idime =1, ndime
                                     auxva = auxva + gpvel(idime, igaus)*gpvel(idime, igaus)
                                  end do
                                  if (kfl_canla>0) then !leaf area index from field
                                     lad =0.0_rp
                                     do inode = 1,pnode
                                        lad = lad + elmar(pelty) % shape(inode,igaus)*canla(lnods(inode,ielem))
                                     end do
                                     cdlad = cd*lad
                                  else if (zmaxf.le.0.01_rp) then ! uniform canopy
                                     cdlad = cd*lad
                                  else
                                     if(1==1) then
                                        !non uniform canopy distribution from LALIC and MIHAILOVIC,
                                        !An Empirical Relation Describing Leaf-Area Density inside the Forest
                                        !http://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%282004%29043%3C0641%3AAERDLD%3E2.0.CO%3B2
                                        lamax = lad / integ
                                        if (hoveg.gt.0.999_rp*heigh ) then
                                           cdlad = 0.0_rp
                                        else if (hoveg.lt.zmaxi) then
                                           cdlad =   cd*lamax * (((heigh -zmaxi) /(heigh -hoveg))**6.0_rp )*&
                                                exp(6.0_rp*(1.0_rp-(heigh -zmaxi) /(heigh -hoveg)))
                                        else
                                           cdlad =   cd*lamax * (((heigh -zmaxi) /(heigh -hoveg))**0.5_rp )*&
                                                exp(0.5_rp*(1.0_rp-(heigh -zmaxi) /(heigh -hoveg)))
                                        end if
                                     else
                                        ! shaw and Schumann 92 - Large-eddy simulation of turbulent flow above and within a forest
                                        ! values for LAI=5
                                        if (hoveg.gt.0.999_rp*heigh ) then
                                           cdlad = 0.0_rp
                                        else
                                           ! digitized from fig 2 shaw
                                           xcoef(1,:) = [ 0.0_rp,    0.09899_rp, 0.20018_rp, 0.29857_rp, 0.3991_rp, 0.59839_rp, 0.69855_rp, 0.79769_rp, 0.7985_rp, 0.89736_rp, 1.0_rp ]
                                           xcoef(2,:) = [ 2.1446_rp, 2.3359_rp,  2.8216_rp,  3.7382_rp,  5.5695_rp,  7.5605_rp, 7.4996_rp,  6.8604_rp,  6.8646_rp, 5.149_rp,   0.0_rp ]
                                           call linint(ncoef_shaw, hoveg/heigh, xcoef, cdlad, aux_dvari)

                                           cdlad = cdlad * 0.15_rp
                                        end if
                                     end if
                                     
                                  end if
                                  prope_ker % value_ielem(ielem) % a(igaus) = densi_ker % value_const(imate)&
                                       *cdlad* sqrt(auxva)
                               else  ! above forest height
                                  prope_ker % value_ielem(ielem) % a(igaus) = 0.0_rp
                               end if
                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'DAFOR' ) then
                   !----------------------------------------------------------------------
                   !
                   ! Darcy-Forchheimer porosity
                   ! S_{i}=-(\mu D+ 1/2 \rho|u_{jj}|F) u_{i}
                   ! beware it only works with contant density and viscosity (visco_ker % value_const(imate))
                   !
                   !----------------------------------------------------------------------
                   darcy = prope_ker % rlaws(1,imate)  ! linear term
                   forch = prope_ker % rlaws(2,imate)  ! quadratic term
                   
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            gpvel = 0.0_rp
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)

                               do idime = 1,ndime
                                  elvel(idime, inode) = veloc(idime, ipoin,1)
                                  do igaus =1, pgaus
                                     gpvel(idime, igaus) = gpvel(idime,igaus) + &
                                          elvel(idime, inode) * elmar(pelty) % shape(inode,igaus)
                                  end do
                               end do
                            end do

                            do igaus =1, pgaus
                               gpmve = sqrt( dot_product(gpvel(1:ndime, igaus), gpvel(1:ndime, igaus)))
                               prope_ker % value_ielem(ielem) % a(igaus) =  darcy * visco_ker % value_const(imate)  &
                                    + 0.5_rp * forch * densi_ker % value_const(imate) * gpmve 
                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'GAUSS' ) then
                   !----------------------------------------------------------------------
                   !
                   ! Gaussian filter to generate turbulent fluctuations by diffusion
                   ! A. Kempf, M. Klein, J. Janicka, Flow Turbul. Combust. (2005)
                   !
                   !       D = C_D * L**2 / T
                   !----------------------------------------------------------------------
                   lscale = prope_ker % rlaws(1,imate)                   ! Target length scale
                   tscale = 1.0_rp / prope_ker % rlaws(2,imate)          ! Temporal integration
                   const  = 0.15915_rp                                   ! C_D = 1 / ( 2*pi)

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)

                            do igaus = 1,pgaus
                               gptur = const * lscale * lscale * tscale              ! D = C_D * L**2 / T
                               prope_ker % value_ielem(ielem) % a(igaus) = gptur

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = 0.0_rp
                               enddo
                            end do
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)

                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do igaub = 1,pgaub
                               gptur = const * lscale * lscale * tscale              ! D = C_D * L**2 / T
                               prope_ker % value_iboun(iboun) % a(igaub) =  gptur
                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'GPSUT' ) then
                   !----------------------------------------------------------------------
                   !
                   !                                                         T0 + S
                   ! Sutherland's law in Gauss points : mu/mu_0 = (T/T0)^1.5 --------
                   !                                                         T  + S
                   !
                   !----------------------------------------------------------------------

                   mu0   = prope_ker % rlaws(1,imate) ! reference value
                   T0    = prope_ker % rlaws(2,imate) ! referentemperature
                   ! Sutherland Temperature S0
                   S0    = prope_ker % rlaws(3,imate)
                   ! if not given, load default value (air)
                   if (S0.lt. epsilon(1.0_rp)) S0    = 110.0_rp

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltem(inode) = tempe(ipoin,which_time)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)
                            gp_temperature=0.0_rp
                            gp_grad_tem=0.0_rp
                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  gp_temperature(igaus) = gp_temperature(igaus) + eltem(inode) * elmar(pelty) % shape(inode,igaus)
                                  do idime = 1,ndime
                                     gp_grad_tem(idime,igaus) = gp_grad_tem(idime,igaus) + eltem(inode) * gpcar(idime,inode,igaus)
                                  end do
                               end do
                               if (gp_temperature(igaus) <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference temperature
                                  gp_temperature(igaus) = T0
                                  gp_grad_tem(1:ndime,igaus) = 0.0_rp
                               endif
                               mu =  mu0 * (gp_temperature(igaus)/T0)**1.5_rp * (T0+S0) / (gp_temperature(igaus)+S0)
                               prope_ker % value_ielem(ielem) % a(igaus) = mu
                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = gp_grad_tem(idime,igaus) * &
                                       (3.0_rp*S0 + gp_temperature(igaus)) / (2.0_rp * gp_temperature(igaus)* (gp_temperature(igaus)+S0))*mu
                               enddo

                               ! derivative of mu (at gauss points) respect to nodal tempe
                               do inode = 1,pnode
                                  prope_ker % drval_ielem(ielem) % a(inode,igaus) = elmar(pelty) % shape(inode,igaus) * &
                                       (3.0_rp*S0 + gp_temperature(igaus)) / (2.0_rp * gp_temperature(igaus)* (gp_temperature(igaus)+S0))*mu
                               enddo

                               ! derivative of grad(mu) (at gauss points) respect to nodal tempe
                               do inode = 1,pnode
                                  do idime = 1, ndime
                                     prope_ker % gdval_ielem(ielem) % a(idime, inode,igaus) = gpcar(idime,inode,igaus) * &
                                          (3.0_rp*S0 + gp_temperature(igaus)) / (2.0_rp * gp_temperature(igaus)* (gp_temperature(igaus)+S0))*mu + &
                                          gp_grad_tem(idime,igaus) * elmar(pelty) % shape(inode,igaus) * &
                                          ( (3.0_rp*S0 + gp_temperature(igaus))**2 / (2.0_rp * gp_temperature(igaus)* (gp_temperature(igaus)+S0))**2*mu + &
                                          ( -2*gp_temperature(igaus)**2 -12*gp_temperature(igaus)*S0 -6*S0**2 ) / &
                                          (4*gp_temperature(igaus)**2*(gp_temperature(igaus)+S0)**2 ) * mu )
                                  enddo
                               enddo

                            end do
                         end if
                      end if
                   end do
                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               eltem(inodb) = tempe(lnodb(inodb,iboun),1)
                            end do
                            do igaub = 1,pgaub
                               gbtem = 0.0_rp
                               do inodb = 1,pnodb
                                  gbtem = gbtem + eltem(inodb) * elmar(pblty) % shape(inodb,igaub)
                               end do
                               if (gbtem <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference temperature
                                  gbtem = 300.0_rp
                               endif
                               prope_ker % value_iboun(iboun) % a(igaub) =   mu0 * (gbtem/T0)**1.5_rp * (T0+S0) / (gbtem+S0)
                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'SMAGO' ) then

                   !----------------------------------------------------------------------
                   !
                   !
                   !  SMAGORINSKY MODEL nu_t = (C_s * Vol)^2 * sqrt(2*Sij*Sij)
                   !
                   !
                   !----------------------------------------------------------------------


                   const = prope_ker % rlaws(1,imate)

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                                  elvel(idime,inode) = veloc(idime,ipoin,which_time)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)

                            gpvel = 0.0_rp
                            gpgve = 0.0_rp
                            !
                            ! HLENG and TRAGL at center of gravity
                            !
                            call elmlen(&
                                 ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                                 hnatu(pelty),hleng)

                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  do idime = 1,ndime
                                     gpvel(idime,igaus) = gpvel(idime,igaus) + elvel(idime,inode) * elmar(pelty) % shape(inode,igaus)
                                     do jdime = 1,ndime
                                        gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) + elvel(idime,inode) * gpcar(jdime,inode,igaus)
                                     end do
                                  end do
                               end do
                               !
                               ! Computing turbulent viscosity mu_t at gauss points
                               ! according to Smagorinsky: gptur
                               !
                               call ker_turbul(1_ip,dummr,dummr,const,dummr,gpgve(:,:,igaus),hleng,gptur,dummr,dummr)

                               prope_ker % value_ielem(ielem) % a(igaus) = gptur

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = 0.0_rp
                               enddo
                            end do
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         pnode = nnode(pelty)
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            do idime = 1,ndime
                               elcod(idime,inode) = coord(idime,ipoin)
                            end do
                         end do
                         !
                         ! HLENG and TRAGL at center of gravity
                         !
                         call elmlen(&
                              ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                              hnatu(pelty),hleng)

                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               do idime = 1,ndime
                                  elvel(idime,inodb) = veloc(idime,lnodb(inodb,iboun),which_time)
                               end do
                            end do
                            do igaub = 1,pgaub
                               gbvel = 0.0_rp
                               gbgve = 0.0_rp
                               do inodb = 1,pnodb
                                  do idime = 1,ndime
                                     gbvel(idime) = gbvel(idime) + elvel(idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                                     do jdime = 1,ndime
                                        gbgve(jdime,idime) = gbgve(jdime,idime) + elvel(idime,inodb) * gpcar(jdime,inodb,igaub)
                                     end do
                                  end do
                               end do
                               !
                               ! Computing turbulent viscosity mu_t at gauss points
                               ! according to Smagorinsky: gptur
                               !
                               call ker_turbul(1_ip,dummr,dummr,const,dummr,gbgve(:,:),hleng,gptur,dummr,dummr)

                               prope_ker % value_iboun(iboun) % a(igaub) =  gptur

                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'SIGMA' ) then

                   !----------------------------------------------------------------------
                   !
                   !
                   !  Sigma - Baya et al., 2010, Proceeding of teh Summer
                   !  Program
                   !
                   !----------------------------------------------------------------------


                   const = prope_ker % rlaws(1,imate)

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                                  elvel(idime,inode) = veloc(idime,ipoin,which_time)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)

                            gpvel = 0.0_rp
                            gpgve = 0.0_rp
                            !
                            ! HLENG and TRAGL at center of gravity
                            !
                            call elmlen(&
                                 ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                                 hnatu(pelty),hleng)

                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  do idime = 1,ndime
                                     gpvel(idime,igaus) = gpvel(idime,igaus) + elvel(idime,inode) * elmar(pelty) % shape(inode,igaus)
                                     do jdime = 1,ndime
                                        gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) + elvel(idime,inode) * gpcar(jdime,inode,igaus)
                                     end do
                                  end do
                               end do
                               !
                               ! Computing turbulent viscosity mu_t at gauss points
                               ! according to Sigma: gptur
                               !
                               call ker_turbul(4_ip,dummr,dummr,const,dummr,gpgve(:,:,igaus),hleng,gptur,dummr,dummr)

                               prope_ker % value_ielem(ielem) % a(igaus) = gptur

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = 0.0_rp
                               enddo
                            end do
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         pnode = nnode(pelty)
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            do idime = 1,ndime
                               elcod(idime,inode) = coord(idime,ipoin)
                            end do
                         end do
                         !
                         ! HLENG and TRAGL at center of gravity
                         !
                         call elmlen(&
                              ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                              hnatu(pelty),hleng)

                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               do idime = 1,ndime
                                  elvel(idime,inodb) = veloc(idime,lnodb(inodb,iboun),which_time)
                               end do
                            end do
                            do igaub = 1,pgaub
                               gbvel = 0.0_rp
                               gbgve = 0.0_rp
                               do inodb = 1,pnodb
                                  do idime = 1,ndime
                                     gbvel(idime) = gbvel(idime) + elvel(idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                                     do jdime = 1,ndime
                                        gbgve(jdime,idime) = gbgve(jdime,idime) + elvel(idime,inodb) * gpcar(jdime,inodb,igaub)
                                     end do
                                  end do
                               end do
                               !
                               ! Computing turbulent viscosity mu_t at gauss points
                               ! according to Sigma: gptur
                               !
                               call ker_turbul(4_ip,dummr,dummr,const,dummr,gbgve(:,:),hleng,gptur,dummr,dummr)

                               prope_ker % value_iboun(iboun) % a(igaub) =  gptur

                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'WALE ' ) then

                   !----------------------------------------------------------------------
                   !
                   !  WALE  - Wall-adapting local eddy-viscosity - Nicoud-Ducros 99
                   !
                   !----------------------------------------------------------------------


                   const = prope_ker % rlaws(1,imate)

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                                  elvel(idime,inode) = veloc(idime,ipoin,which_time)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)

                            gpvel = 0.0_rp
                            gpgve = 0.0_rp

                            !                            call elmgeo_element_node_length(ndime,pnode,elcod,hleng)  ! for ker_turbul
                            !
                            ! HLENG and TRAGL at center of gravity
                            !
                            call elmlen(&
                                 ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                                 hnatu(pelty),hleng)
                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  do idime = 1,ndime
                                     gpvel(idime,igaus) = gpvel(idime,igaus) + elvel(idime,inode) * elmar(pelty) % shape(inode,igaus)
                                     do jdime = 1,ndime
                                        gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) + elvel(idime,inode) * gpcar(jdime,inode,igaus)
                                     end do
                                  end do
                               end do
                               !
                               ! Computing turbulent viscosity mu_t at gauss points
                               ! according to WALE: gptur
                               !
                               call ker_turbul(3_ip,dummr,dummr,const,dummr,gpgve(:,:,igaus),hleng,gptur,dummr,dummr)

                               prope_ker % value_ielem(ielem) % a(igaus) = gptur

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = 0.0_rp
                               end do
                            end do
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         pnode = nnode(pelty)
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            do idime = 1,ndime
                               elcod(idime,inode) = coord(idime,ipoin)
                            end do
                         end do
                         !
                         ! HLENG and TRAGL at center of gravity
                         !
                         call elmlen(&
                              ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                              hnatu(pelty),hleng)

                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               do idime = 1,ndime
                                  elvel(idime,inodb) = veloc(idime,lnodb(inodb,iboun),which_time)
                               end do
                            end do

                            do igaub = 1,pgaub
                               gbvel = 0.0_rp
                               gbgve = 0.0_rp
                               do inodb = 1,pnodb
                                  do idime = 1,ndime
                                     gbvel(idime) = gbvel(idime) + elvel(idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                                     do jdime = 1,ndime
                                        gbgve(jdime,idime) = gbgve(jdime,idime) + elvel(idime,inodb) * gpcar(jdime,inodb,igaub)
                                     end do
                                  end do
                               end do
                               !
                               ! Computing turbulent viscosity mu_t at gauss points
                               ! according to WALE: gptur
                               !
                               call ker_turbul(3_ip,dummr,dummr,const,dummr,gbgve(:,:),hleng,gptur,dummr,dummr)

                               prope_ker % value_iboun(iboun) % a(igaub) =  gptur

                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'VRMAN' ) then

                   !----------------------------------------------------------------------
                   !
                   !  VRMAN - Vreman SGS model static version
                   !
                   !----------------------------------------------------------------------

                   const = prope_ker % rlaws(1,imate)

                   !--------------------------------------------------------------------------     
                   !$OMP PARALLEL DO                                                         &
                   !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                      & 
                   !$OMP SHARED       ( par_omp_nelem_chunk )                                &
                   !--------------------------------------------------------------------------     
                   !--------------------------------------------------------------------------     
                   !$OMP DEFAULT      ( NONE )                                               &
                   !$OMP PRIVATE      ( pnode,pgaus,ielem,inode,ipoin,elcod,elvel,idime,     &
                   !$OMP                gpvol,gpcar,dummr,gpvel,gpgve,tragl,hleng,gptur,     &
                   !$OMP                pelty )                                              &
                   !$OMP SHARED       ( ngaus,lnods,lnnod,which_time,elmar,nelem,imate,      &
#ifndef NDIMEPAR
                   !$OMP                ndime,                                               &
#endif
                   !$OMP                prope_ker,lmate,hnatu,veloc,coord,ltype,const        )   
                   !--------------------------------------------------------------------------
                   
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                                  elvel(idime,inode) = veloc(idime,ipoin,which_time)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)
                            ! call elmgeo_element_node_length(ndime,pnode,elcod,hleng)  ! for ker_turbul
                            !
                            ! HLENG and TRAGL at center of gravity
                            !
                            call elmlen(&
                                 ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                                 hnatu(pelty),hleng)
                            gpvel = 0.0_rp
                            gpgve = 0.0_rp                          
                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  do idime = 1,ndime
                                     gpvel(idime,igaus) = gpvel(idime,igaus) + elvel(idime,inode) * elmar(pelty) % shape(inode,igaus)
                                     do jdime = 1,ndime
                                        gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) + elvel(idime,inode) * gpcar(jdime,inode,igaus)
                                     end do
                                  end do
                               end do
                               !
                               ! Computing turbulent viscosity mu_t at gauss points
                               ! according to VRMAN: gptur
                               !
                               call ker_turbul(5_ip,dummr,dummr,const,dummr,gpgve(:,:,igaus),hleng,gptur,dummr,dummr)
                               prope_ker % value_ielem(ielem) % a(igaus) = gptur

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = 0.0_rp
                               enddo
                            end do
                         end if
                      end if
                   end do
                   !$OMP END PARALLEL DO

                   prope_ker % kfl_nedsm = 1
                   
                   do iboun = 1,nboun
                      prope_ker % value_iboun(iboun) % a(:) =  0.0_rp
                   end do
 
                else if( prope_ker % wlaws(imate) == 'ILSA ' ) then

                   !---------------------------------------------------------------------------------
                   !
                   !  ILSA
                   !  Dynamic subfilter-scale stress model for large-eddy simulations,
                   !  A. Rouhi, U. Piomelli and B.J. Geurts, PHYSICAL REVIEW FLUIDS 1, 044401 (2016)
                   ! 
                   !---------------------------------------------------------------------------------

                   const = prope_ker % rlaws(1,imate)

                   !--------------------------------------------------------------------------     
                   !$OMP PARALLEL DO                                                         &
                   !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                      & 
                   !$OMP SHARED       ( par_omp_nelem_chunk )                                &
                   !--------------------------------------------------------------------------     
                   !--------------------------------------------------------------------------     
                   !$OMP DEFAULT      ( NONE )                                               &
                   !$OMP PRIVATE      ( pnode,pgaus,ielem,inode,ipoin,elcod,elvel,idime,     &
                   !$OMP                gpvol,gpcar,dummr,gpvel,gpgve,tragl,hleng,           &
                   !$OMP                pelty,gpnuILSA,gpden,gpvis )                         &
                   !$OMP SHARED       ( ngaus,lnods,lnnod,which_time,elmar,nelem,imate,      &
#ifndef NDIMEPAR
                   !$OMP                ndime,                                               &
#endif
                   !$OMP                prope_ker,lmate,hnatu,veloc,coord,ltype,const        )   
                   !--------------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)

                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                                  elvel(idime,inode) = veloc(idime,ipoin,which_time)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)

                            call elmlen(&
                                 ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                                 hnatu(pelty),hleng)
                            do igaus = 1,pgaus
                               do idime = 1,ndime
                                  gpvel(idime,igaus) = 0.0_rp
                                  do jdime = 1,ndime
                                     gpgve(jdime,idime,igaus) = 0.0_rp
                                  end do
                               end do
                               do inode = 1,pnode
                                  do idime = 1,ndime
                                     gpvel(idime,igaus) = gpvel(idime,igaus) + elvel(idime,inode) * elmar(pelty) % shape(inode,igaus)
                                     do jdime = 1,ndime
                                        gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) &
                                           + gpcar(jdime,inode,igaus) * elvel(idime,inode)
                                     end do
                                  end do
                               end do
                            end do

                            do igaus = 1,pgaus
                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               gpnuILSA(igaus) = (gpvis(1))/(gpden(1)+zeror)
                            end do

                            call ker_ILSA_sgs_viscosity(pgaus,1_ip,pgaus,gpnuILSA,& 
                               prope_ker % value_ielem(ielem) % a,gpvel,gpgve,ielem, prope_ker % rlaws(1,imate),& 
                               prope_ker % rlaws(2,imate),hleng)
                            
                            do idime = 1,ndime
                               prope_ker % grval_ielem(ielem) % a(idime,1:pgaus) = 0.0_rp
                            end do
                            
                         end if
                      end if
                   end do
                   !$OMP END PARALLEL DO

                   prope_ker % kfl_nedsm = 1
               
                   do iboun = 1,nboun
                      prope_ker % value_iboun(iboun) % a(:) =  0.0_rp
                   end do

                else if( prope_ker % wlaws(imate) == 'MIXIN' ) then

                   !----------------------------------------------------------------------
                   !
                   !
                   !  ALGEBRAIC MIXING LENGTH MODEL -- nu_t = (kappa * yplus)^2 * sqrt(2*Sij*Sij) * (1-exp(yplus/A))^2
                   !
                   !
                   !----------------------------------------------------------------------

                   const = prope_ker % rlaws(1,imate)  ! Here comes A+ = 26 
                   !
                   ! 1) Obtain tange - the tangential componenet of the shear stress at the wall - idem postpr
                   !
                   ! TANGE: Tangential stress
                   !
                   if( INOTMASTER ) then
                      call memgen(zero,ndime,npoin)
                      call runend (' can not call nsi_outtan here -- need to find a solution') 
!                      call nsi_outtan()              ! generates non zero values in gevec(idime,ipoin)
                   end if
                   !
                   ! 2) Extrapolate it to the interior
                   ! Allocate and obtain tau_wall_boun  -  the shear stress at wall nodes
                   !
                   kpoin = wallcoupling_extr_boun % wet % npoin_wet  ! temporary use of kpoin
                   call memory_alloca(mem_modul(1:2,modul),'tau_wall_boun','mod_ker_proper',tau_wall_boun, ndime, kpoin)

                   call COU_GET_INTERPOLATE_POINTS_VALUES(gevec,tau_wall_boun,wallcoupling_extr_boun)    ! tau_wall_boun(kount)  , tau_wall is the nodal shear stress
                   !
                   ! The extrapolated shear stress at an interior node is obtained - for boundary nodes I already have the correct value in gevec
                   !
                   kpoin =0
                   do ipoin= 1,npoin
                      if(is_interior(ipoin) /= 0_ip) then ! interior
                         kpoin = kpoin + 1
                         gevec(1:ndime,ipoin) = tau_wall_boun( 1:ndime, ptb_to_use(kpoin) )
                      end if
                   end do
                   call memory_deallo(mem_modul(1:2,modul),'tau_wall_boun','mod_ker_proper',tau_wall_boun)   ! ojo tal ves no hacerlo cada vez
                   !
                   ! 3) Prepare everything that is needed for ker_turbul(6_ip
                   !
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                                  elvel(idime,inode) = veloc(idime,ipoin,which_time)
                                  eltan(idime,inode) = gevec(idime,ipoin)
                               end do
                               elwal(inode) = walld(ipoin)
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)

                            gpvel = 0.0_rp
                            gptan = 0.0_rp
                            gpgve = 0.0_rp
                            gpwal = 0.0_rp
                            !
                            ! HLENG and TRAGL at center of gravity
                            !
                            call elmlen(&
                                 ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                                 hnatu(pelty),hleng)

                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  do idime = 1,ndime
                                     gpwal = gpwal + elmar(pelty) % shape(inode,igaus) * elwal(inode)
                                     gpvel(idime,igaus) = gpvel(idime,igaus) + elvel(idime,inode) * elmar(pelty) % shape(inode,igaus)
                                     gptan(idime,igaus) = gptan(idime,igaus) + eltan(idime,inode) * elmar(pelty) % shape(inode,igaus)
                                     do jdime = 1,ndime
                                        gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) + elvel(idime,inode) * gpcar(jdime,inode,igaus)
                                     end do
                                  end do
                               end do

                               velmo = 0.0_rp
                               tanmo = 0.0_rp
                               do idime= 1,ndime
                                  velmo = velmo + gpvel(idime,igaus) * gpvel(idime,igaus)
                                  tanmo = tanmo + gptan(idime,igaus) * gptan(idime,igaus)
                               end do

                               velmo = sqrt(velmo)
                               tanmo = sqrt(tanmo)

                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               !
                               ! Computing turbulent viscosity mu_t at gauss points
                               ! according to mixing length model: gptur
                               !
                               call ker_turbul(6_ip,gpwal,velmo,const,gpvis(1),gpgve(:,:,igaus),hleng,gptur,gpden(1),tanmo)   

                               prope_ker % value_ielem(ielem) % a(igaus) = gptur

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = 0.0_rp
                               enddo
                            end do
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         pnode = nnode(pelty)
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            do idime = 1,ndime
                               elcod(idime,inode) = coord(idime,ipoin)
                            end do
                         end do
                         !
                         ! HLENG and TRAGL at center of gravity
                         !
                         call elmlen(&
                              ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                              hnatu(pelty),hleng)

                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               do idime = 1,ndime
                                  elvel(idime,inodb) = veloc(idime,lnodb(inodb,iboun),which_time)
                                  eltan(idime,inodb) = gevec(idime,lnodb(inodb,iboun))
                               end do
                               elwal(inodb) = walld(lnodb(inodb,iboun))
                            end do
                            do igaub = 1,pgaub
                               gbvel = 0.0_rp
                               gbtan = 0.0_rp
                               gbgve = 0.0_rp
                               gbwal = 0.0_rp
                               do inodb = 1,pnodb
                                  do idime = 1,ndime
                                     gbvel(idime) = gbvel(idime) + elvel(idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                                     gbtan(idime) = gbvel(idime) + eltan(idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                                     do jdime = 1,ndime
                                        gbgve(jdime,idime) = gbgve(jdime,idime) + elvel(idime,inodb) * gpcar(jdime,inodb,igaub)
                                     end do
                                     gbwal = gbwal + elmar(pblty) % shape(inodb,igaub) * elwal(inodb)
                                  end do
                               end do

                               velmo = 0.0_rp
                               tanmo = 0.0_rp
                               do idime= 1,ndime
                                  velmo = velmo + gbvel(idime) * gbvel(idime)
                                  tanmo = tanmo + gbtan(idime) * gbtan(idime)
                               end do

                               velmo = sqrt(velmo)
                               tanmo = sqrt(tanmo)

                               call ker_proper('DENSI','IGAUB',igaub,iboun,gbden)
                               call ker_proper('VISCO','IGAUB',igaub,iboun,gbvis)
                               !
                               ! Computing turbulent viscosity mu_t at gauss points
                               ! according to mixing length model: gptur
                               !
                               call ker_turbul(6_ip,gbwal,velmo,const,gbvis(1),gbgve(:,:),hleng,gptur,gbden(1),tanmo)

                               prope_ker % value_iboun(iboun) % a(igaub) =  gptur

                            end do
                         end if
                      end if
                   end do

                   if( iavec == 1 ) call memgen(2_ip,one,one)
                   nullify(gevec)
                   
                else if( prope_ker % wlaws(imate) == 'SSTKO' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Turbulent viscosity calculated at gauss points for SST K-Omega model
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltur(1:2,inode) = untur(1:2,ipoin,which_time)
                               elwal(inode) = walld(ipoin)
                               elvel(1:ndime, inode) = veloc(1:ndime, ipoin,1)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)
                            do igaus = 1,pgaus
                               kinen = 0.0_rp
                               omega = 0.0_rp
                               gpwal = 0.0_rp

                               do idime = 1,ndime
                                  do jdime = 1,ndime
                                     gpgve(jdime,idime, igaus) = 0.0_rp
                                  end do
                               end do
                               do inode = 1,pnode
                                  !gpsha(inode, igaus) =  elmar(pelty) % shape(inode,igaus)
                                  kinen = kinen + elmar(pelty) % shape(inode,igaus) * eltur(1,inode)
                                  omega = omega + elmar(pelty) % shape(inode,igaus) * eltur(2,inode)
                                  gpwal = gpwal + elmar(pelty) % shape(inode,igaus) * elwal(inode)
                                  do idime = 1,ndime
                                     do jdime = 1,ndime
                                        gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) &
                                             + gpcar(jdime,inode,igaus) * elvel(idime,inode)
                                     end do
                                  end do
                               end do
                               if (kinen <= epsilon(1.0_rp) ) kinen = 1.0e10_rp
                               if (omega <= epsilon(1.0_rp) ) omega = 1.0e10_rp
                               !
                               ! ERROR in calculation of gpvor, gpvor = 2 W_ij W_ij where W_ij = sqrt(0.5 (dui_dxj - duj_dxi))
                               !
                               !                                gpvor = 0.0_rp   ! 2 S_ij : S_ij
                               !                                do idime = 1,ndime
                               !                                   do jdime = 1,ndime
                               !                                      gpvor = gpvor + gpgve(idime,jdime,igaus) &    ! 2 S_ij : S_ij
                               !                                           *(gpgve(idime,jdime,igaus)        &
                               !                                           +gpgve(jdime,idime,igaus))
                               !                                   end do
                               !                                end do

                               gpvor_mat = 0.0_rp
                               do idime = 1,ndime
                                  do jdime = 1,ndime
                                     gpvor_mat(idime,jdime) = 0.5_rp*(gpgve(jdime, idime,igaus) - gpgve(idime, jdime,igaus))
                                  end do
                               end do
                               gpvor = 0.0_rp
                               do idime = 1,ndime
                                  do jdime = 1,ndime
                                     gpvor = gpvor + 2.0_rp * gpvor_mat(idime,jdime) * gpvor_mat(idime,jdime)
                                  end do
                               end do

                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               gpvor = sqrt(gpvor)
                               a1 =0.31_rp
                               gpnut = 0.0_rp
                               if (gpwal.gt.1.0e-8_rp) then
!!! faltaba una gpden(1)
                                  F2 = max( 2.0_rp *sqrt(max(0.0_rp,kinen)) /(0.09_rp*omega*gpwal), 500.0_rp * gpvis(1) / (gpden(1)*gpwal*gpwal*omega))
                                  !                                   F2 = max( 2.0_rp *sqrt(max(0.0_rp,kinen)) /(0.09_rp*omega*gpwal), 500.0_rp * gpvis(1) / (gpwal*gpwal*omega)) !!!!!
                                  F2 = tanh( F2 * F2 )
                                  as = max( a1 * omega, gpvor * F2 )
                                  gpnut = a1*kinen/as
                               end if
                               ! check the limits
                               gpmut = max( gpden(1)*gpnut, 0.1_rp*gpvis(1))
                               gpmut = min(gpmut,1.0_rp * 1.0e20_rp)
                               prope_ker % value_ielem(ielem) % a(igaus) = gpmut

                               !
                               ! derivative of gpmut (at gauss points) respect to nodal kin and ome
                               ! arg2 = max (c1 , c2)                 F2 = TANH(arg2 * arg2)
                               !
                               ! dF2/dk = 2*arg2*d(arg2)/dk * (1-F2*F2)
                               ! dF2/dw = 2*arg2*d(arg2)/dw * (1-F2*F2)
                               !
                               do inode=1,pnode
                                  dgpkin(inode) = elmar(pelty) % shape(inode,igaus)
                                  dgpome(inode) = elmar(pelty) % shape(inode,igaus)
                               enddo

                               arg2 = max( 2.0_rp *sqrt(max(0.0_rp,kinen)) /(0.09_rp*omega*gpwal), 500.0_rp * gpvis(1) / (gpden(1)*gpwal*gpwal*omega))

                               ! calculation of d(arg2)/dk and d(arg2)/dw
                               darg2_dkin = 0.0_rp
                               darg2_dome = 0.0_rp
                               if (kinen > 0.0_rp) then
                                  if (2.0_rp *sqrt(max(0.0_rp,kinen)) /(0.09_rp*omega*gpwal) > 500.0_rp * gpvis(1) / (gpden(1)*gpwal*gpwal*omega) ) then
                                     do inode=1,pnode
                                        darg2_dkin(inode) = dgpkin(inode)/(0.09_rp*omega*gpwal*sqrt(max(0.0_rp,kinen)))
                                        darg2_dome(inode) = -2.0_rp*sqrt(max(0.0_rp,kinen))*dgpome(inode)/(0.09_rp*omega*omega*gpwal)
                                     enddo
                                  else
                                     do inode=1,pnode
                                        darg2_dkin(inode) = 0.0_rp
                                        darg2_dome(inode) = -500.0_rp * gpvis(1) * dgpome(inode)/ (gpden(1)*gpwal*gpwal*omega*omega)
                                     enddo
                                  endif
                               else
                                  if (0.0_rp > 500.0_rp * gpvis(1) / (gpden(1)*gpwal*gpwal*omega) ) then
                                     do inode=1,pnode
                                        darg2_dkin(inode) = 0.0_rp
                                        darg2_dome(inode) = 0.0_rp
                                     enddo
                                  else
                                     do inode=1,pnode
                                        darg2_dkin(inode) = 0.0_rp
                                        darg2_dome(inode) = -500.0_rp * gpvis(1) * dgpome(inode)/ (gpden(1)*gpwal*gpwal*omega*omega)
                                     enddo
                                  endif
                               endif

                               ! calculation of dF2/dk and dF2/dw
                               do inode = 1,pnode
                                  dF2_dkin(inode) = 2.0_rp*arg2*darg2_dkin(inode)*(1.0_rp-F2*F2)
                                  dF2_dome(inode) = 2.0_rp*arg2*darg2_dome(inode)*(1.0_rp-F2*F2)
                               enddo

                               ! calculation of d(nu_t)/dk and d(nu_t)/dw
                               dgpnut_dkin = 0.0_rp
                               dgpnut_dome = 0.0_rp
                               if (a1 * omega > gpvor * F2) then
                                  do inode = 1,pnode
                                     dgpnut_dkin(inode) = dgpkin(inode)/omega
                                     dgpnut_dome(inode) = -kinen*dgpome(inode)/(omega*omega)
                                  enddo
                               else
                                  do inode = 1,pnode
                                     dgpnut_dkin(inode) = a1*dgpkin(inode)/(gpvor * F2) - a1*kinen*dF2_dkin(inode)/(gpvor*F2*F2)
                                     dgpnut_dome(inode) = -a1*kinen*dF2_dome(inode)/(gpvor*F2*F2)
                                  enddo
                               endif

                               ! check limits
                               if ( gpden(1)*gpnut > 1.0e20_rp .or. gpden(1)*gpnut < 0.1_rp*gpvis(1) ) then
                                  dgpnut_dkin = 0.0_rp
                                  dgpnut_dome = 0.0_rp
                               endif
                               ! pass d(nu_t)/dk and d(nu_t)/dw to the kernel variables
                               do inode = 1,pnode
                                  prope_ker % drval_tur_ielem(ielem) % a(1,inode,igaus) = gpden(1) * dgpnut_dkin(inode)
                                  prope_ker % drval_tur_ielem(ielem) % a(2,inode,igaus) = gpden(1) * dgpnut_dome(inode)
                               end do

                               !
                               ! derivative of gpmut (at gauss points) respect to nodal veloc
                               ! arg2 = max (c1 , c2)                 F2 = TANH(arg2 * arg2)
                               ! dF2/dU = 2*arg2*d(arg2)/du * (1-F2*F2)
                               ! d(arg2)/du = d(gpvor)/du
                               ! gpvor = sqrt[ (du1/dx2 - du2_dx1)**2 + (du1/dx3 - du3/dx1)**2 + (du3/dx2 - du2/dx3)**2 ]
                               !
                               if (ndime == 2) then
                                  do inode = 1,pnode
                                     dgpvor_du(inode,1) = gpcar(2,inode,igaus)*(gpgve(2,1,igaus)-gpgve(1,2,igaus))/gpvor
                                     dgpvor_du(inode,2) = -gpcar(1,inode,igaus)*(gpgve(2,1,igaus)-gpgve(1,2,igaus))/gpvor
                                  enddo
                               elseif (ndime == 3) then
                                  do inode = 1,pnode
                                     dgpvor_du(inode,1) = ( gpcar(2,inode,igaus)*(gpgve(2,1,igaus)-gpgve(1,2,igaus)) &
                                          + gpcar(3,inode,igaus)*(gpgve(3,1,igaus)-gpgve(1,3,igaus)) ) /gpvor
                                     dgpvor_du(inode,2) = (-gpcar(1,inode,igaus)*(gpgve(2,1,igaus)-gpgve(1,2,igaus)) &
                                          + gpcar(3,inode,igaus)*(gpgve(3,2,igaus)-gpgve(2,3,igaus)) ) /gpvor
                                     dgpvor_du(inode,3) = ( gpcar(1,inode,igaus)*(gpgve(3,1,igaus)-gpgve(1,3,igaus)) &
                                          - gpcar(2,inode,igaus)*(gpgve(3,2,igaus)-gpgve(2,3,igaus)) ) /gpvor
                                  enddo
                               endif

                               ! calculation of d(nu_t)/du
                               dgpnut_du = 0.0_rp
                               if (a1 * omega < gpvor * F2) then
                                  do inode = 1,pnode
                                     do idime = 1,ndime
                                        dgpnut_du(inode,idime) = -a1*kinen*dgpvor_du(inode,idime)/(F2*gpvor*gpvor)
                                     enddo
                                  enddo
                               endif

                               ! check limits
                               if ( gpden(1)*gpnut > 1.0e20_rp .or. gpden(1)*gpnut < 0.1_rp*gpvis(1) ) then
                                  dgpnut_du = 0.0_rp
                               endif

                               ! pass d(nu_t)/du to the kernel variables
                               do inode = 1,pnode
                                  do idime = 1,ndime
                                     prope_ker % drval_vel_ielem(ielem) % a(idime,inode,igaus) = gpden(1) * dgpnut_du(inode,idime)
                                  enddo
                               end do

                            end do ! igaus

                            ! gradient
                            prope_ker % grval_ielem(ielem) % a = 0.0_rp

                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'STDKE' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Turbulent viscosity calculated at gauss points for standard K-Epsilon model
                   !
                   !----------------------------------------------------------------------

                   if (kfl_logva==1) then
                      relat =1.0_rp
                   else
                      relat =1.0_rp
                   end if
                   ! Realizable constant depending on cmu value
                   if (kfl_kemod_ker.ne.0) &
                        A0 = (1.0_rp - 3.0_rp*sqrt(0.5_rp*cmu_st))/cmu_st
                               
                   element_loop:  do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltur(1:2,inode) = untur(1:2,ipoin,which_time)
                               elvel(1:ndime, inode) = veloc(1:ndime, ipoin,1)
                               elcod(1:ndime,inode) = coord(1:ndime,ipoin)    
                            end do
                            gauss_loop:  do igaus = 1,pgaus
                               kinen = 0.0_rp
                               epsil = 0.0_rp
                               do inode = 1,pnode
                                  kinen = kinen + eltur(1,inode) * elmar(pelty)%shape(inode,igaus)
                                  epsil = epsil + eltur(2,inode) * elmar(pelty)%shape(inode,igaus)
                               end do
                               if (kfl_regularization) then
                                  reguk  = regul_k(kinen)
                                  regue  = regul_e(epsil)
                               else
                                  reguk  = kinen
                                  regue  = epsil
                               end if
                               if (kfl_kemod_ker==0.or.ittim==0) then ! k eps standard
                                  cmu = cmu_st
                                  
                               else if (kfl_kemod_ker==2.or.kfl_kemod_ker==3) then  ! Realizable or kefp models, modify cmu
                                  call elmca2(&
                                       pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                                       elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                                       gpder,gpcar,gphes,ielem)
         
                                  gpgve(1:ndime,1:ndime,igaus) = 0.0_rp
                                  do inode = 1,pnode
                                     do idime = 1,ndime
                                        do jdime = 1,ndime
                                           gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) &
                                                + gpcar(jdime,inode,igaus) * elvel(idime,inode)
                                        end do
                                     end do
                                  end do
                                  divve =0.0_rp       ! divergence of u
                                  do idime = 1,ndime  ! symmetric gradient (strain tensor)
                                     divve = divve + gpgve(idime, idime, igaus)
                                     do jdime = 1,ndime
                                        simgr(idime,jdime) = 0.5_rp*(gpgve(idime, jdime,igaus) &
                                             +  gpgve(jdime, idime,igaus))
                                     end do
                                  end do
                                  seci4 = 0.0_rp
                                  W=0.0_rp
                                  uaste =0.0_rp
                                  do idime = 1,ndime !deviator, necessary to compute W
                                     simgr(idime, idime ) = simgr(idime, idime) - divve/3.0_rp 
                                  end do
                                  do idime = 1,ndime
                                     do jdime = 1,ndime
                                        seci4 = seci4 + simgr(idime, jdime) &        ! S_ij : S_ij
                                             *simgr(idime, jdime)
                                        uaste = uaste +gpgve(idime,jdime,igaus) &    ! D_ij : D_ij
                                             *(gpgve(idime,jdime,igaus))   

                                        do kdime =1, ndime
                                           W = W +  simgr(idime,jdime)* &
                                                simgr(jdime,kdime)* &
                                                simgr(kdime,idime)
                                        end do
                                     end do
                                  end do
                                  uaste = sqrt(uaste -divve*divve /3.0_rp)
                                  if (seci4.gt.epsilon(1.0_rp)) W = W/sqrt(seci4*seci4*seci4)

                                  Wsqr6 = sqrt(6.0_rp)*W     

                                  Wsqr6 = min(Wsqr6,  1.0_rp)  ! corrects limits
                                  Wsqr6 = max(Wsqr6, -1.0_rp)  ! it should not happen

                                  phi_r = 1.0_rp/3.0_rp*acos(Wsqr6)
                                  As_r= sqrt(6.0_rp)*cos(phi_r)

                                  ! calculates turb viscosity
                                  if (kfl_kemod_ker==2) then ! REALIZABLE MODEL
                                     ! effective cmu
                                     cmu =1.0_rp/(A0+As_r*reguk/regue*uaste)
                                  else if (kfl_kemod_ker==3) then ! KEFP
                                     Cr  = 4.5_rp
                                     f0  = Cr/(Cr-1.0_rp)
                                     f02 = 4.0_rp*f0*(f0-1.0_rp)
                                     sigmr = cmu_st *(uaste*reguk/regue)*(uaste*reguk/regue)
                                     ! effective cmu
                                     cmu   = cmu_st *2.0_rp*f0 /(1.0_rp+ sqrt(1.0_rp + f02*sigmr))
                                  end if

 !                                 Gpmut(igaus)= max(gpden(igaus)*cmu*kinen*kinen/epsil, gpvis(igaus))

!                                  if (kfl_regularization) &
!                                       gpmut(igaus)  = cmu*gpden(igaus)*reguk*reguk/regue

                               end if ! end ke model
                               
                               ! density and viscosity
                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,elmar(pelty)%shape)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,elmar(pelty)%shape)
                               if (kfl_logva==0) then
                                  if (regue.gt.0) then
                                     gpmut = gpden(1)* cmu *  reguk*reguk/regue
                                  else
                                     gpmut =0.0_rp
                                  end if
                                  if (.not.kfl_regularization) &
                                       gpmut = max(gpmut,gpvis(1))
                               else
                                  gpmut = gpden(1)* cmu_st * exp(2.0_rp*reguk- regue)
                               end if
                               prope_ker % value_ielem(ielem) % a(igaus) = gpmut
                               !  USE RELAXATION
                               !                                    max(relat*gpden(1)*gpnut + (1.0_rp-relat)*  &
                               !                                    prope_ker % value_ielem(ielem) % a(igaus) , gpvis(1))
                               
                            end do gauss_loop
                            ! gradients! convergence problems
                            prope_ker % grval_ielem(ielem) % a = 0.0_rp
                         end if ! pelty > 0
                      end if
                   end do element_loop

                else if( prope_ker % wlaws(imate) == 'SPALA' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Turbulent viscosity calculated at gauss points for Spalart-Allmaras model
                   !
                   !----------------------------------------------------------------------
                   cv1 = 7.1_rp
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltur(1,inode) = untur(1,ipoin,which_time)
                               !                                eltur(1,inode) = untur_fix(ipoin)
                            end do
                            do igaus = 1,pgaus

                               gpnut = 0.0_rp
                               do inode = 1,pnode
                                  gpsha(inode,igaus) = elmar(pelty)%shape(inode,igaus)
                                  gpnut = gpnut + gpsha(inode, igaus)*eltur(1,inode)
                               end do
                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,gpsha)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,gpsha)

                               gpnu  = gpvis(1)/gpden(1)
                               X     = gpnut/gpnu
                               fv1 = X*X*X/(X*X*X + cv1*cv1*cv1)
                               gpmut = gpden(1)*fv1*gpnut

                               ! pass mu_t to the kernel variables
                               prope_ker % value_ielem(ielem) % a(igaus) = gpmut

                               ! calculation of d(mu_t)/d(nu_t)
                               do inode = 1,pnode
                                  dX_dnut = gpsha(inode,igaus)/gpnu
                                  dfv1_dnut = ( 3.0_rp*X*X*dX_dnut*(X*X*X + cv1*cv1*cv1) - 3.0_rp*X*X*dX_dnut*X*X*X )/(X*X*X + cv1*cv1*cv1)**2.0_rp
                                  dgpmut_dnut(inode) = gpden(1)*dfv1_dnut*gpnut + gpden(1)*fv1*gpsha(inode, igaus)
                               enddo

                               ! pass d(mu_t)/d(nu_t) to the kernel variables
                               do inode = 1,pnode
                                  prope_ker % drval_tur_ielem(ielem) % a(1,inode,igaus) = dgpmut_dnut(inode)
                               end do

                            end do

                         end if
                      end if
                   end do


                else if( prope_ker % wlaws(imate) == 'LOWMA' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Low Mach approximation, density = p(th) / R T
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)

                      if (tempe(ipoin,which_time) /= 0.0_rp ) then
                         prope_ker % value_ipoin(ipoin) = prthe(which_time) /(gasco * tempe(ipoin,which_time))
                      else
                         prope_ker % value_ipoin(ipoin) = 0.0_rp
                      endif
                   end do

                else if( prope_ker % wlaws(imate) == 'LOWMG' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Low Mach approximation, density = p(th) / R T IN GAUSS POINTS
                   !
                   !----------------------------------------------------------------------
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            if (associated(tesgs)) then
                               do igaus = 1,pgaus
                                  gptem = 0.0_rp
                                  do inode = 1,pnode
                                     ipoin = lnods(inode,ielem)
                                     gptem = gptem + tempe(ipoin,which_time) * elmar(pelty) % shape(inode,igaus)
                                  end do
                                  gptem = gptem+ tesgs(ielem)%a(1,igaus,1)
                                  if (gptem <= epsilon(1.0_rp) ) then
                                     ! For initialization we default to reference temperature
                                     gptem = 300.0_rp
                                  endif
                                  prope_ker % value_ielem(ielem) % a(igaus) = prthe(which_time)/gasco/gptem
                               end do
                            else
                               do igaus = 1,pgaus
                                  gptem = 0.0_rp
                                  do inode = 1,pnode
                                     ipoin = lnods(inode,ielem)
                                     gptem = gptem + tempe(ipoin,which_time) * elmar(pelty) % shape(inode,igaus)
                                  end do
                                  if (gptem <= epsilon(1.0_rp) ) then
                                     ! For initialization we default to reference temperature
                                     gptem = 300.0_rp
                                  endif
                                  prope_ker % value_ielem(ielem) % a(igaus) = prthe(which_time)/gasco/gptem
                               end do
                            end if

                         endif
                      end if
                   end do

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               eltem(inodb) = tempe(lnodb(inodb,iboun),1)
                            end do
                            do igaub = 1,pgaub
                               gbtem = 0.0_rp
                               do inodb = 1,pnodb
                                  gbtem = gbtem + eltem(inodb) * elmar(pblty) % shape(inodb,igaub)
                               end do
                               if (gbtem <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference temperature
                                  gbtem = 300.0_rp
                               endif
                               prope_ker % value_iboun(iboun) % a(igaub) =  prthe(1)/gasco/gbtem
                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'KLOWM' .or. prope_ker % wlaws(imate) == 'TLOWM' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Low Mach approximation, density = p(th) W / R T  where W is mean molecular density
                   !
                   ! Low Mach approximation for the syncronized CFI combustion model,
                   ! density is only updated in temper
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      if (tempe(ipoin,which_time) > 0.0_rp ) then
                         prope_ker % value_ipoin(ipoin) = prthe(which_time) * wmean(ipoin,which_time) /(gasco * tempe(ipoin,which_time))
                      else
                         prope_ker % value_ipoin(ipoin) = 0.0_rp
                      endif
                   end do

                else if( prope_ker % wlaws(imate) == 'GKLOW' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Low Mach approximation, density = p(th) W / R T at gauss_points where W is mean molecular density
                   !                         d(density)/dx = p(th) dW/dx / (R T) - p(th) W R/ (R T)^2
                   !
                   !----------------------------------------------------------------------
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)

                            gp_grad_tem = 0.0_rp
                            gp_grad_wmean = 0.0_rp

                            do igaus = 1,pgaus

                               gptem = 0.0_rp
                               gpwmean = 0.0_rp
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  gptem = gptem + tempe(ipoin,which_time) * elmar(pelty) % shape(inode,igaus)
                                  if (kfl_coupl(ID_TEMPER,ID_CHEMIC) /= 0 ) then
                                     gpwmean = gpwmean + wmean(ipoin,which_time) * elmar(pelty) % shape(inode,igaus)
                                  else
                                     gpwmean = gpwmean + 1.0_rp * elmar(pelty) % shape(inode,igaus)
                                  end if
                                  do idime = 1,ndime
                                     gp_grad_tem(idime,igaus)   = gp_grad_tem(idime,igaus) + tempe(ipoin,which_time) * gpcar(idime,inode,igaus)
                                     if (kfl_coupl(ID_TEMPER,ID_CHEMIC) /= 0 ) then
                                        gp_grad_wmean(idime,igaus) = gp_grad_wmean(idime,igaus) + wmean(ipoin,which_time) * gpcar(idime,inode,igaus)
                                     else
                                        gp_grad_wmean(idime,igaus) = 0.0_rp
                                     end if
                                  end do
                               end do

                               if (associated(tesgs)) then
                                  gptem = gptem + tesgs(ielem)%a(1,igaus,1)
                               end if
                               if (gptem <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference temperature
                                  gptem = 300.0_rp
                               end if
                               rho = prthe(which_time)*gpwmean/(gasco*gptem)
                               prope_ker % value_ielem(ielem) % a(igaus) = rho

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = prthe(which_time) * gp_grad_wmean(idime,igaus)/(gasco*gptem) &
                                       - prthe(which_time)*gpwmean*gp_grad_tem(idime,igaus)/(gasco*gptem**2)
                               end do

                               ! derivative of rho (at gauss points) respect to nodal tempe
                               do inode = 1,pnode
                                  prope_ker % drval_ielem(ielem) % a(inode,igaus) = -elmar(pelty) % shape(inode,igaus) * rho/gptem
                               end do

                               ! derivative of d(rho)/dx (at gauss points) respect to nodal tempe
                               do inode = 1,pnode
                                  do idime = 1, ndime
                                     prope_ker % gdval_ielem(ielem) % a(idime, inode,igaus) = &
                                          - prthe(which_time) * gp_grad_wmean(idime,igaus)/(gasco*gptem**2)*elmar(pelty) % shape(inode,igaus) &
                                          - rho/gptem * gpcar(idime,inode,igaus) &
                                          + 2*rho/(gptem**2)*gp_grad_tem(idime,igaus) * elmar(pelty) % shape(inode,igaus)
                                  end do
                               end do

                            end do !igaus
                         end if
                      end if
                   end do
                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               eltem(inodb) = tempe(lnodb(inodb,iboun),which_time)
                               if (kfl_coupl(ID_TEMPER,ID_CHEMIC) /= 0 ) then
                                  elwmean(inodb) = wmean(lnodb(inodb,iboun),which_time)
                               else
                                  elwmean(inodb) = 1.0_rp
                               end if
                            end do
                            do igaub = 1,pgaub
                               gbtem = 0.0_rp
                               gbwme = 0.0_rp
                               do inodb = 1,pnodb
                                  gbtem = gbtem + eltem(inodb) * elmar(pblty) % shape(inodb,igaub)
                                  gbwme = gbwme + elwmean(inodb) * elmar(pblty) % shape(inodb,igaub)
                               end do
                               if (gbtem <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference
                                  ! temperature
                                  gbtem = 300.0_rp
                               endif
                               prope_ker % value_iboun(iboun) % a(igaub) = prthe(1)*gbwme/(gasco*gbtem)
                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'MIXTU' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Mixture of constant density fluids, 1/density = Sum_k  Y_k / rho_k where Y_k is the mass fraction
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = 0.0_rp
                      do ispec = 1,nspec
                         prope_ker % value_ipoin(ipoin) = prope_ker % value_ipoin(ipoin) + conce(ipoin,ispec,which_time)/speci(ispec)%densi(1)
                      enddo
                      if ( prope_ker % value_ipoin(ipoin) .gt. 0.0_rp) then
                         prope_ker % value_ipoin(ipoin) = 1.0_rp /  prope_ker % value_ipoin(ipoin)
                      else
                         call runend("KERMOD: Density evaluated to <= 0, check species properties")
                      endif
                   end do

                else if( prope_ker % wlaws(imate) == 'BIPHA' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Mixture of constant density fluids in two phases
                   ! 1/density = Sum_k  Y_k / rho_k where Y_k is the mass fraction
                   ! If level set is negative, density corresponds to first phase, positive to second phase
                   !
                   !----------------------------------------------------------------------

                   if( thicl > 0.0_rp ) then
                      eps = thicl
                   else
                      call runend('Really need hleng in Ker_updpro.')
                      !eps = -thicl * hleng(1)
                   end if
                   ooeps = 1.0_rp/eps
                   oovpi = 1.0_rp/pi

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      ! Compute density of phases
                      !
                      densi_phase(1) = 0.0_rp
                      densi_phase(2) = 0.0_rp
                      do ispec = 1,nspec
                         densi_phase(1) = densi_phase(1) + conce(ipoin,ispec,1)/speci(ispec)%densi(1)
                         densi_phase(2) = densi_phase(2) + conce(ipoin,ispec,1)/speci(ispec)%densi(2)
                      enddo
                      if ( densi_phase(1) .gt. 0.0_rp) then
                         densi_phase(1) = 1.0_rp /  densi_phase(1)
                      else
                         call runend("KERMOD: Density in phase 1 evaluated to <= 0, check species properties")
                      endif
                      if ( densi_phase(2) .gt. 0.0_rp) then
                         densi_phase(2) = 1.0_rp /  densi_phase(2)
                      else
                         call runend("KERMOD: Density in phase 2 evaluated to <= 0, check species properties")
                      endif
                      ! Compute density at point
                      phi=fleve(ipoin,which_time)
                      if( phi < -eps ) then
                         ! First phase
                         prope_ker % value_ipoin(ipoin)=densi_phase(1)
                      else if( phi > eps ) then
                         ! Second phase
                         prope_ker % value_ipoin(ipoin)=densi_phase(2)
                      else
                         f = 0.5_rp*(1.0_rp+phi*ooeps+sin(pi*phi*ooeps)*oovpi)
                         prope_ker % value_ipoin(ipoin)=densi_phase(1) + (densi_phase(2)-densi_phase(1))*f
                      endif

                   end do

                else if( prope_ker % wlaws(imate) == 'TBIPH' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Mixture of temperature dependent density fluids in two phases
                   ! 1/density = Sum_k  Y_k / rho_k where Y_k is the mass fraction
                   ! If level set is negative, density corresponds to first phase, positive to second phase
                   ! Temperature dependence is assumed linear, rho = rho_a * T + rho_0
                   !
                   !----------------------------------------------------------------------

                   if( thicl > 0.0_rp ) then
                      eps = thicl
                   else
                      call runend('Really need hleng in Ker_updpro.')
                      !eps = -thicl * hleng(1)
                   end if
                   ooeps = 1.0_rp/eps
                   oovpi = 1.0_rp/pi

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      ! Compute density of phases
                      !
                      densi_phase(1) = 0.0_rp
                      densi_phase(2) = 0.0_rp
                      do ispec = 1,nspec
                         densi_phase(1) = densi_phase(1) + conce(ipoin,ispec,which_time)/ &
                              ( speci(ispec)%densi(1)*tempe(ipoin,which_time)+speci(ispec)%densi(2) )
                         densi_phase(2) = densi_phase(2) + conce(ipoin,ispec,which_time)/ &
                              ( speci(ispec)%densi(3)*tempe(ipoin,which_time)+speci(ispec)%densi(4) )
                      enddo
                      if ( densi_phase(1) .gt. 0.0_rp) then
                         densi_phase(1) = 1.0_rp /  densi_phase(1)
                      else
                         call runend("KERMOD: Density in phase 1 evaluated to <= 0, check species properties")
                      endif
                      if ( densi_phase(2) .gt. 0.0_rp) then
                         densi_phase(2) = 1.0_rp /  densi_phase(2)
                      else
                         call runend("KERMOD: Density in phase 2 evaluated to <= 0, check species properties")
                      endif
                      ! Compute density at point
                      phi=fleve(ipoin,which_time)
                      if( phi < -eps ) then
                         ! First phase
                         prope_ker % value_ipoin(ipoin)=densi_phase(1)
                      else if( phi > eps ) then
                         ! Second phase
                         prope_ker % value_ipoin(ipoin)=densi_phase(2)
                      else
                         f = 0.5_rp*(1.0_rp+phi*ooeps+sin(pi*phi*ooeps)*oovpi)
                         prope_ker % value_ipoin(ipoin)=densi_phase(1) + (densi_phase(2)-densi_phase(1))*f
                      endif

                   end do

                else if( prope_ker % wlaws(imate) == 'DNAST' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Density is controlled by NASTAL module
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = densi(ipoin,which_time)
                   end do

                else if( prope_ker % wlaws(imate) == 'MUTAB' .or. prope_ker % wlaws(imate) == 'TMUTA') then

                   !----------------------------------------------------------------------
                   !
                   ! Viscosity is read from table in CHEMIC module (CFI combustion model)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = visck(ipoin,1)
                   end do

                else if( prope_ker % wlaws(imate) == 'KTABL' .or. prope_ker % wlaws(imate) == 'TKTAB') then

                   !----------------------------------------------------------------------
                   !
                   ! Heat conductivity is read from table in CHEMIC module (CFI combustion model)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = condk(ipoin,1)
                   end do

                else if( prope_ker % wlaws(imate) == 'KQUAR') then

                   !----------------------------------------------------------------------
                   !
                   ! Heat conductivity for quartz glass is computed from polynomial function
                   ! based on temperature
                   ! Source: O. Sergeev, A. Shashkov, A. Umanskii, Thermophysical properties of
                   !         quartz glass, Journal of engineering physics 43 (1982) 1375–1383
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time)
                      prope_ker % value_ipoin(ipoin) = 2.14749_rp - 298.76_rp * T**(-1) + 20720.0_rp * T**(-2) - 540000.0_rp * T**(-3)
                   end do

                else if( prope_ker % wlaws(imate) == 'CPQUA') then

                   !----------------------------------------------------------------------
                   !
                   ! Specific heat for quartz glass is computed from polynomial function
                   ! based on temperature
                   ! Source: O. Sergeev, A. Shashkov, A. Umanskii, Thermophysical properties of
                   !         quartz glass, Journal of engineering physics 43 (1982) 1375–1383
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time)
                      prope_ker % value_ipoin(ipoin) = 931.3_rp + 0.256_rp * T - 24000000.0_rp * T**(-2)
                   end do

                else if( prope_ker % wlaws(imate) == 'KSTAI') then

                   !----------------------------------------------------------------------
                   !
                   ! Heat conductivity for stainless steel is computed from polynomial
                   ! function based on temperature
                   ! Source: EN1993 1.2 Annex C
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time) - 273.15_rp
                      prope_ker % value_ipoin(ipoin) = 14.6_rp + 0.0127_rp * T
                   end do

                else if( prope_ker % wlaws(imate) == 'CPSTA') then

                   !----------------------------------------------------------------------
                   !
                   ! Specific heat for stainless steel is computed from polynomial
                   ! function based on temperature
                   ! Source: EN1993 1.2 Annex C
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time) - 273.15_rp
                      prope_ker % value_ipoin(ipoin) = 450.0_rp + 0.28_rp * T - 0.000291_rp * T**2_rp + 0.000000134_rp * T**3_rp
                   end do

                else if( prope_ker % wlaws(imate) == 'CPTAB' .or. prope_ker % wlaws(imate) == 'TCPTA') then

                   !----------------------------------------------------------------------
                   !
                   ! Specific heat is read from table in CHEMIC module (CFI combustion model)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = sphek(ipoin,1)
                   end do

                else if( prope_ker % wlaws(imate) == 'MUMIX' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Viscosity given by a mixture of species (individual viscosities are computed
                   ! inside CHEMIC, here only the mixture is done)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      ! Mixture coefficients see Combustion Physics, Chung K. Law, page 155
                      !
                      ! phi(i,j) = 1/sqrt(8) 1/sqrt(1 + Wi/Wj) [1+ sqtr(mui/muj) (Wj/Wi)^1/4 ]^2
                      ! phi(i,i) = 1.0
                      !
                      ! Took out wmean(ipoin,1) from all terms because it simplifies in the final expression
                      !
                      do ispec = 1,nspec
                         sumphi(ispec) = 0.0_rp
                         do jspec = 1,nspec
                            if (jspec /= ispec) then
                               if (visck(ipoin,jspec) /= 0.0_rp) then
                                  phiij = sqrt(0.125_rp/(1.0_rp+speci(ispec)%weigh/speci(jspec)%weigh)) * &
                                       (1.0_rp+sqrt( visck(ipoin,ispec)&
                                       / visck(ipoin,jspec))*(speci(jspec)%weigh/speci(ispec)%weigh)**0.25_rp)**2_rp
                               else
                                  phiij =  sqrt(0.125_rp/(1.0_rp+speci(ispec)%weigh/speci(jspec)%weigh))
                               endif
                               sumphi(ispec) = sumphi(ispec) + phiij * conce(ipoin,jspec,which_time) / speci(jspec)%weigh
                            endif
                         enddo
                      enddo
                      !
                      prope_ker % value_ipoin(ipoin) = 0.0_rp
                      do ispec = 1,nspec
                         if (conce(ipoin,ispec,which_time) /= 0.0_rp) then
                            prope_ker % value_ipoin(ipoin) = prope_ker % value_ipoin(ipoin) + visck(ipoin,ispec) * &
                                 (conce(ipoin,ispec,which_time) / speci(ispec)%weigh) &
                                 / (conce(ipoin,ispec,which_time) / speci(ispec)%weigh + sumphi(ispec))
                         endif

                      enddo
                   end do

                else if( prope_ker % wlaws(imate) == 'KMIXT' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Heat conductivity given by a mixture of species (individual conductivities are computed
                   ! inside CHEMIC, here only the mixture is done)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      ! Mixture coefficients see Combustion Physics, Chung K. Law, page 155
                      ! Notice that it's ok that viscosities are used...
                      !
                      ! phi(i,j) = 1/sqrt(8) 1/sqrt(1 + Wi/Wj) [1+ sqtr(mui/muj) (Wj/Wi)^1/4 ]^2
                      ! phi(i,i) = 1.0
                      !
                      ! Took out wmean(ipoin,1) from all terms because it simplifies in the final expression
                      !
                      do ispec = 1,nspec
                         sumphi(ispec) = 0.0_rp
                         do jspec = 1,nspec
                            if (jspec /= ispec) then
                               if (visck(ipoin,jspec) /= 0.0_rp) then
                                  phiij = sqrt(0.125_rp/(1.0_rp+speci(ispec)%weigh/speci(jspec)%weigh)) * &
                                       (1.0_rp+sqrt( visck(ipoin,ispec)&
                                       / visck(ipoin,jspec))*(speci(jspec)%weigh/speci(ispec)%weigh)**0.25_rp)**2
                               else
                                  phiij =  sqrt(0.125_rp/(1.0_rp+speci(ispec)%weigh/speci(jspec)%weigh))
                               endif
                               sumphi(ispec) = sumphi(ispec) + phiij * conce(ipoin,jspec,which_time) / speci(jspec)%weigh
                            endif
                         enddo
                      enddo
                      !
                      prope_ker % value_ipoin(ipoin) = 0.0_rp
                      do ispec = 1,nspec
                         if (conce(ipoin,ispec,which_time) /= 0.0_rp) then
                            prope_ker % value_ipoin(ipoin) = prope_ker % value_ipoin(ipoin) + condk(ipoin,ispec) * &
                                 (conce(ipoin,ispec,which_time) / speci(ispec)%weigh) &
                                 / (conce(ipoin,ispec,which_time) / speci(ispec)%weigh + 1.065_rp * sumphi(ispec))
                         endif
                      enddo
                   end do

                else if( prope_ker % wlaws(imate) == 'CPMIX' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Specific heat of a mixture of species (individual speheas are computed
                   ! inside CHEMIC, here only the mixture is done)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      prope_ker % value_ipoin(ipoin) = 0.0_rp
                      do ispec = 1,nspec
                         prope_ker % value_ipoin(ipoin) = prope_ker % value_ipoin(ipoin) + sphek(ipoin,ispec) * conce(ipoin,ispec,which_time)
                      enddo
                   end do

                else if( prope_ker % wlaws(imate) == 'TEST1' ) then

                   !----------------------------------------------------------------------
                   !
                   ! TEST1
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            gpcor = 0.0_rp
                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  do idime = 1,ndime
                                     gpcor(idime,igaus) = gpcor(idime,igaus) + coord(idime,ipoin) * elmar(pelty) % shape(inode,igaus)
                                  end do
                               end do
                               prope_ker % value_ielem(ielem) % a(igaus)   =  2.0_rp * gpcor(1,igaus) - 3.0_rp * gpcor(2,igaus)
                               prope_ker % grval_ielem(ielem) % a(1,igaus) =  2.0_rp
                               prope_ker % grval_ielem(ielem) % a(2,igaus) = -3.0_rp
                            end do
                         end if
                      end if
                   end do
                   prope_ker % kfl_nedsm = 1

                else if( prope_ker % wlaws(imate) == 'TEST2' ) then

                   !----------------------------------------------------------------------
                   !
                   ! TEST2
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = 2.0_rp * coord(1,ipoin) - 3.0_rp * coord(2,ipoin)
                   end do

                else if( prope_ker % wlaws(imate) == 'TEST3' ) then

                   !----------------------------------------------------------------------
                   !
                   ! TEST3
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin)   =  2.0_rp * coord(1,ipoin) - 3.0_rp * coord(2,ipoin)
                      prope_ker % grval_ipoin(1,ipoin) =  2.0_rp
                      prope_ker % grval_ipoin(2,ipoin) = -3.0_rp
                   end do

                else if( prope_ker % wlaws(imate) == 'TEST4' ) then

                   !----------------------------------------------------------------------
                   !
                   ! TEST4
                   !
                   !----------------------------------------------------------------------

                   if( ndime == 2 ) then
                      do kpoin = 1,nmatn(imate)
                         ipoin = lmatn(imate) % l(kpoin)
                         prope_ker % value_ipoin(ipoin)   = 1.0_rp + 2.0_rp * coord(1,ipoin) + 3.0_rp * coord(2,ipoin)
                         prope_ker % grval_ipoin(1,ipoin) = 2.0_rp
                         prope_ker % grval_ipoin(2,ipoin) = 3.0_rp
                      end do
                   else
                      do kpoin = 1,nmatn(imate)
                         ipoin = lmatn(imate) % l(kpoin)
                         prope_ker % value_ipoin(ipoin)   = &
                              1.0_rp + 2.0_rp * coord(1,ipoin) + 3.0_rp * coord(2,ipoin) + 4.0_rp * coord(3,ipoin)
                         prope_ker % grval_ipoin(1,ipoin) = 2.0_rp
                         prope_ker % grval_ipoin(2,ipoin) = 3.0_rp
                         prope_ker % grval_ipoin(3,ipoin) = 4.0_rp
                      end do
                   end if

                else if( prope_ker % wlaws(imate) == 'ABL  ' ) then

                   !----------------------------------------------------------------------
                   !
                   ! ABL: mu = rho/Cmu*kap*(z+z0)*U_*
                   !
                   !----------------------------------------------------------------------

                   prope_ker % kfl_nedsm = 1
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = abs(ltype(ielem))
                         pgaus = ngaus(pelty)
                         pnode = lnnod(ielem)
                         gpcor = 0.0_rp
                         do igaus = 1,pgaus
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  gpcor(idime,igaus) = gpcor(idime,igaus) + coord(idime,ipoin) * elmar(pelty) % shape(inode,igaus)
                               end do
                            end do
                            if( kfl_rough > 0 ) then
                               z0 = 0.0_rp
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  z0    = z0 + rough(ipoin) * elmar(pelty) % shape(inode,igaus)
                               end do
                            else
                               z0 = rough_dom
                            end if
                            rho   = 1.0_rp
                            ustar = 0.23_rp
                            kap   = 0.41_rp
                            prope_ker % value_ielem(ielem) % a(      igaus) =  rho*ustar*kap*(gpcor(ndime,igaus)+z0+delta_dom)
                            prope_ker % grval_ielem(ielem) % a(    1,igaus) =  0.0_rp
                            prope_ker % grval_ielem(ielem) % a(    2,igaus) =  0.0_rp
                            prope_ker % grval_ielem(ielem) % a(ndime,igaus) =  rho*ustar*kap
                         end do
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == ker_polynomial_name ) then
                   call ker_polynomial_proper(prope_ker, imate, which_time)

                end if

             end do

          end if

       end if

    end do

    deallocate(listp_ker)

    !----------------------------------------------------------------------
    !
    ! Assign previous values for adjoint case to temper and conce
    !
    !----------------------------------------------------------------------

    if( kfl_adj_prob == 1_ip ) then
       conce => conce_save
       tempe => tempe_save
       untur => untur_save
       veloc => veloc_save
       nullify(conce_save)
       nullify(tempe_save)
       nullify(untur_save)
       nullify(veloc_save)
    end if

  end subroutine ker_updpro

end module mod_ker_proper
!> @}
