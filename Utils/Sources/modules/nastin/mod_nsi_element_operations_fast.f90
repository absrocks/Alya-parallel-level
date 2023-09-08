!------------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    mod_nsi_element_operations_fast.f90  - started from Alya 6954.
!> @author  Guillaume Houzeaux & Herbert Owen
!> @brief   Navier-Stokes system element assembly only of Auu part and send it to RHS
!>          The matrices are created at the elemental level but then they are multiplid by u_n and sent to RHS (explicit only).
!>          Moreover, boundary conditions will be applied after calling to the solver as is usually done in the explicit case. MISSING hhh
!>          For the moment I leave it as it is in mod_nsi_element_operations - when guillaume uploads dirichlet algorithmic I will change it
!> @details Elemental operations. Now there are no tasks.
!>
!>          It is only valid for LES (no RANS), no thermal coupling, emac , divergence form , force just gravity , no stablization, 3d
!>          We only have calls to:
!>                 1) element_shape_function_derivatives_jacobian - could be eliminated if I leave only tetras
!>                 2) cputim
!>                 3) ker_proper - could also be eliminated
!>
!>          OTHER SIMPLICATIONS:
!>
!>         1) When I inlined nsi_rhodt_rhotau_nu_vector   I only included the part for porde==1  (linear elements)   !! added in nsi_outerr
!>                    I am still missing add an error in nsi_outerr  if there are more than linear elements
!>                    I did this to avoid an if and to avid analysing whta is done for high order elements.
!>
!>         2) local time step does not work --  I set dtinv_loc(1:VECTOR_SIZE) = dtinv_nsi added
!>
!>         3) This only works for kfl_matdi_nsi ==2 (NSI_DIRICHLET_ALGORITHM) -added
!>
!>         4) added runened if nor diverg and emac
!>
!>
!>
!>          \verbatim
!>
!>          1 ........ Element calculations and assembly of global system:
!>                     b <= b^(e) - A^(e)*u^(e): RHS ................ RHSID
!>
!>          \endverbatim
!>
!>          CORRESPONDANCE OLD TO NEW SUBROUTINES: hhh delete
!>          --------------------------------------
!>
!>          nsi_elmma4         <= nsi_element_assembly_split_oss
!>          nsi_elmga3         <= nsi_element_operations_gather
!>          nsi_elmlen         <= elmgeo_element_characteristic_length
!>          elmchl             <= elmgeo_element_length
!>          elmca2             <= element_shape_function_derivatives_jacobian 
!>          nsi_elmtss         <= nsi_element_time_step
!>          nsi_elmres         <= nsi_element_residual
!>          nsi_updsgs         <= nsi_element_subgrid_scale
!>          nsi_elmsgs         <= nsi_element_stabilization
!>          nsi_elmort         <= nsi_assembly_projections
!>          nsi_elmope_omp     <= nsi_element_operations
!>          nsi_elmmat         <= nsi_element_assembly_asgs_oss
!>                                nsi_element_assembly_asgs_oss_old
!>          nsi_elmdi3         <= nsi_element_dirichlet
!>                                nsi_element_schur
!>          nsi_assemble_schur <= nsi_assembly_schur_method
!>          nsi_elmext         <= nsi_element_extension
!>          nsi_elmexa         <= nsi_element_manufactured_solution
!>          nsi_elmexf         <= nsi_element_external_force
!>
!>
!> @}
!------------------------------------------------------------------------

module mod_nsi_element_operations_fast

  use def_kintyp,                     only : ip,rp   !in_const
#ifndef VECTOR_SIZE
  use def_master,                     only : VECTOR_SIZE   !in_const
#endif
  use def_domain,                     only : ndime  !in_const_scalar  ! actually paarmeter if ndimepar
#ifndef SUPER_FAST
  use mod_element_integration,        only : element_shape_function_derivatives_jacobian
#endif

#ifdef _OPENACC
  use openacc
#endif

  use def_master,                     only : kfl_paral
  
  implicit none
  private
  public :: nsi_element_operations_fast
  public :: nsi_element_operations_fast5
  public :: nsi_element_operations_fast8



contains

  subroutine nsi_element_operations_fast(&
       pnode,pgaus,list_elements,time1)

    use def_kintyp,            only : ip,rp  ! in_const_scalar
    use def_master,            only : rhsid  ! out       ! real(rp), pointer     :: rhsid(:)
    use def_master,            only : veloc  ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
    use def_domain,            only : coord  ! in_const  ! real(rp), pointer     :: coord(:,:)
    use def_domain,            only : ltype  ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
    use def_domain,            only : lnods  ! in_const  ! integer(ip), pointer  :: lnods(:,:)
    use def_domain,            only : mnode  ! in_const scalar
#ifndef SUPER_FAST
    use def_domain,            only : elmar
    use mod_ker_proper,        only : ker_proper
#endif
    use def_domain,            only : ntens      ! in_const scalar
    use def_nastin,            only : dtinv_nsi  ! in_var   scalar
    use def_nastin,            only : dt_rho_nsi ! out  ! real(rp), pointer     :: dt_rho_nsi(:)
    use def_nastin,            only : mu_rho_nsi ! out  ! real(rp), pointer     :: mu_rho_nsi(:)
    use def_nastin,            only : grnor_nsi  ! in_const_scalar
    use def_nastin,            only : gravi_nsi  ! in_const  ! real(rp)         ::gravi_nsi(3)
    !
    ! For skews this part will be eliminated when guillaume uploads dirichlet algorithmic
#ifdef SUPER_FAST
    use def_kermod, only       :  densi_ker,visco_ker  ! in_const type(typ_valpr_ker), target   :: visco_ker
    ! to avoid problems we usa a real densi_aux = densi_ker % rlaws(1,1)
    ! once it works we could revert this to see if derived data types work fine
#endif

    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in)          :: pnode                        !< Number of nodes
    integer(ip), intent(in)          :: pgaus                        !< Number of Gauss points
    integer(ip), intent(in)          :: list_elements(VECTOR_SIZE)   !< List of elements

    real(rp),    intent(inout)       :: time1(10)                    ! Timings
    !
    ! Element matrices and vectors (stiffness and preconditioner)
    !
    real(rp)    :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)    :: elrbu(VECTOR_SIZE,ndime,pnode)                    ! bu
    real(rp)    :: eldtrho(VECTOR_SIZE,pnode)                        ! Projection of rho/dt
    real(rp)    :: elmurho(VECTOR_SIZE,pnode)                        ! Projection of mu/rho
    !
    ! Gather
    !
    real(rp)    :: elvel(VECTOR_SIZE,ndime,pnode)                       ! u
    real(rp)    :: elcod(VECTOR_SIZE,ndime,pnode)                         ! x
    !
    ! Indices and dimensions
    !
    integer(ip) :: ielem,inode,ivect
    integer(ip) :: pevat
    integer(ip) :: pelty,plapl
    integer(ip) :: ipoin,igaus
    integer(ip) :: lnods_loc(VECTOR_SIZE,pnode)
    integer(ip) :: list_elements_p(VECTOR_SIZE)                      ! List of elements (always positive)
    !
    ! Gauss point values
    !
    real(rp)    :: gpsha(VECTOR_SIZE,pnode,pgaus)                    ! N
    real(rp)    :: gpder(VECTOR_SIZE,ndime,pnode,pgaus)              ! dN/dsi
    real(rp)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)              ! dN/dxi
    real(rp)    :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)              ! d2N/dxidxj
    real(rp)    :: gpvol(VECTOR_SIZE,pgaus)                          ! w*|J|, |J|
    real(rp)    :: gpvis(VECTOR_SIZE,pgaus)                          ! Viscosity
    real(rp)    :: gpmut(VECTOR_SIZE,pgaus)                          ! mut
    real(rp)    :: gpden(VECTOR_SIZE,pgaus)                          ! Density
    real(rp)    :: gpadv(VECTOR_SIZE,ndime,pgaus)                    ! u+u'
    real(rp)    :: gprhs(VECTOR_SIZE,ndime,pgaus)                    ! RHS
    real(rp)    :: gpvel(VECTOR_SIZE,ndime,pgaus)                    ! u
    real(rp)    :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)              ! grad(u)
    !
    ! Internal
    !
#ifdef OPENACCHHH

#define FACT0X     fact0
#define FACT1X     fact1
#define FACT2X     fact2
#define FACT4X     fact4
#define FACT5X     fact5
#define FACT6X     fact6
#define DTINV_LOCX dtinv_loc
#define T1X        t1
#define T2X        t2
#define T3X        t3
#define DENOMX     denom

#else

#define FACT0X     fact0(1:VECTOR_SIZE)
#define FACT1X     fact1(1:VECTOR_SIZE)
#define FACT2X     fact2(1:VECTOR_SIZE)
#define FACT4X     fact4(1:VECTOR_SIZE)
#define FACT5X     fact5(1:VECTOR_SIZE)
#define FACT6X     fact6(1:VECTOR_SIZE)
#define DTINV_LOCX dtinv_loc(1:VECTOR_SIZE)
#define T1X        t1(1:VECTOR_SIZE)
#define T2X        t2(1:VECTOR_SIZE)
#define T3X        t3(1:VECTOR_SIZE)
#define DENOMX     denom(1:VECTOR_SIZE)

#endif

    real(rp)    :: FACT0X
    real(rp)    :: FACT1X
    real(rp)    :: FACT2X
    real(rp)    :: FACT4X
    real(rp)    :: FACT5X
    real(rp)    :: FACT6X
    real(rp)    :: DTINV_LOCX
#ifdef SUPER_FAST
    real(rp)    :: T1X
    real(rp)    :: T2X
    real(rp)    :: T3X
    real(rp)    :: DENOMX
#endif

    real(rp)    :: timea,timeb
    integer(ip) :: idime,jdime,jnode,idofv,jdofv,ievat,jevat,iauxi

    ! Local arrays - previously in element_assembly
    real(rp)    :: wgrgr(VECTOR_SIZE,pnode,pnode,pgaus)
    real(rp)    :: agrau(VECTOR_SIZE,pnode,pgaus)
    !
    ! For shape when SUPER_FAST
    !
#ifdef SUPER_FAST
    real(rp)                             :: weigp(4)
    real(rp)                             :: gpdet(VECTOR_SIZE,pgaus)
    real(rp)                             :: densi_aux
    real(rp)                             :: visco_aux
    real(rp)                             :: xjaci(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),parameter                   :: epsilgeo_div = epsilon(1.0_rp) !epsilgeo usad to avoid divisions by zero
    real(rp)                             :: sha_aux(4,4) = reshape((/ 0.585410196624969D+00 ,  0.138196601125011D+00 , &
         0.138196601125011D+00 ,  0.138196601125011D+00 ,  0.138196601125010D+00 ,  &
         0.585410196624969D+00 ,  0.138196601125011D+00 ,  0.138196601125011D+00 ,  &
         0.138196601125010D+00 ,  0.138196601125011D+00 ,  0.585410196624969D+00 , &
         0.138196601125011D+00 ,  0.138196601125010D+00 ,&
         0.138196601125011D+00 ,  0.138196601125011D+00 ,  0.585410196624969D+00/), (/4,4/))
#endif


#ifndef SUPER_FAST
    integer(ip) :: dummi
#endif


#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------
#ifdef SUPER_FAST
    densi_aux = densi_ker % rlaws(1,1)
    visco_aux = visco_ker % rlaws(1,1)
#endif

    call cputim(timea)
    DTINV_LOCX = dtinv_nsi
    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    pevat = ndime * pnode
    plapl = 0
#ifdef SUPER_FAST
    weigp(1:4)= 1.0_rp/24.0_rp
#endif

    list_elements_p = list_elements

    !$acc enter data create(agrau , lnods_loc , gpvol , elauu ,      &
    !$acc    eldtrho , elvel , gpmut     , gpadv , gpvel ,      &
    !$acc    gpgve   , elcod , gpsha     , gpcar ,              &
    !$acc    gpden   , xjaci , wgrgr     ,                      &
    !$acc    gpdet   , gprhs , elmurho   , elrbu , gpvis     )  &
    !$acc    copyin( weigp  , list_elements , sha_aux )


#ifndef OPENACCHHH
    call cputim(timea)
#endif

    !
    ! This do starts here both in the openacc version and in the nonopenacc
    ! for the non opencc it ends 30  lines later,
    ! in the openacc case it covers all the subroutine.
    ! Similarly the scatter in the nonopenacc case needs a do ivect
    !
    !$acc parallel loop gang vector default(present) async
    do ivect = 1,VECTOR_SIZE
       ielem = abs(list_elements(ivect))
       if ( ielem /= 0 ) then
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
       else
          list_elements_p(ivect)   = list_elements(1)
          lnods_loc(ivect,1:pnode) = 0
       end if

       ielem = list_elements(ivect)
       if ( ielem > 0 ) then
          !
          ! Transient
          !
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
             do idime = 1,ndime
                elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
                elcod(ivect,idime,inode) = coord(idime,ipoin)
             end do
          end do
       else
          !
          ! Element number is null
          !
          elcod(ivect,:,:)   = 0.0_rp
          elvel(ivect,:,:) = 0.0_rp
       end if

#ifndef OPENACCHHH
    end do
    call cputim(timeb)
    time1(1) = time1(1) + timeb - timea
    call cputim(timea)
#endif


    !--------------------------------------------------------------------
    !
    ! Element shape functions and derivatives and properties
    ! Here I coded the superfast version taht is with tetras inlined and without using ker_proper
    ! The non superfast version is not ready but we should try to recover it at some point at least teh sahpe functions part
    !
    !--------------------------------------------------------------------

#ifdef SUPER_FAST

    do igaus=1,pgaus
       do inode=1,pnode
          gpsha(DEF_VECT,inode,igaus) = sha_aux(inode,igaus)   !elmar(pelty) % shape(inode,igaus)
       end do
    end do
    !
    ! GPCAR (from elmgeo_cartesian_derivatives_jacobian_vector), and GPVOL
    !
    !
    ! 3D P1 element
    !
    gpcar(DEF_VECT,1,1,1) =  elcod(DEF_VECT,1,2)   - elcod(DEF_VECT,1,1)
    gpcar(DEF_VECT,1,2,1) =  elcod(DEF_VECT,1,3)   - elcod(DEF_VECT,1,1)
    gpcar(DEF_VECT,1,3,1) =  elcod(DEF_VECT,1,4)   - elcod(DEF_VECT,1,1)
    gpcar(DEF_VECT,2,1,1) =  elcod(DEF_VECT,2,2)   - elcod(DEF_VECT,2,1)
    gpcar(DEF_VECT,2,2,1) =  elcod(DEF_VECT,2,3)   - elcod(DEF_VECT,2,1)
    gpcar(DEF_VECT,2,3,1) =  elcod(DEF_VECT,2,4)   - elcod(DEF_VECT,2,1)
    gpcar(DEF_VECT,3,1,1) =  elcod(DEF_VECT,3,2)   - elcod(DEF_VECT,3,1)
    gpcar(DEF_VECT,3,2,1) =  elcod(DEF_VECT,3,3)   - elcod(DEF_VECT,3,1)
    gpcar(DEF_VECT,3,3,1) =  elcod(DEF_VECT,3,4)   - elcod(DEF_VECT,3,1)
    T1X          =  gpcar(DEF_VECT,2,2,1) * gpcar(DEF_VECT,3,3,1) - gpcar(DEF_VECT,3,2,1) * gpcar(DEF_VECT,2,3,1)
    T2X          = -gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,3,3,1) + gpcar(DEF_VECT,3,1,1) * gpcar(DEF_VECT,2,3,1)
    T3X          =  gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,3,2,1) - gpcar(DEF_VECT,3,1,1) * gpcar(DEF_VECT,2,2,1)
    gpdet(DEF_VECT,1)     =  gpcar(DEF_VECT,1,1,1) * T1X + gpcar(DEF_VECT,1,2,1) * T2X + gpcar(DEF_VECT,1,3,1) * T3X

    DENOMX       =  1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT,1))*max(abs(gpdet(DEF_VECT,1)),epsilgeo_div))

    xjaci(DEF_VECT,1,1,1) =  T1X * DENOMX
    xjaci(DEF_VECT,2,1,1) =  T2X * DENOMX
    xjaci(DEF_VECT,3,1,1) =  T3X * DENOMX
    xjaci(DEF_VECT,2,2,1) = ( gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,3,3,1) - gpcar(DEF_VECT,3,1,1) * gpcar(DEF_VECT,1,3,1)) * DENOMX
    xjaci(DEF_VECT,3,2,1) = (-gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,3,2,1) + gpcar(DEF_VECT,1,2,1) * gpcar(DEF_VECT,3,1,1)) * DENOMX
    xjaci(DEF_VECT,3,3,1) = ( gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,2,2,1) - gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,1,2,1)) * DENOMX
    xjaci(DEF_VECT,1,2,1) = (-gpcar(DEF_VECT,1,2,1) * gpcar(DEF_VECT,3,3,1) + gpcar(DEF_VECT,3,2,1) * gpcar(DEF_VECT,1,3,1)) * DENOMX
    xjaci(DEF_VECT,1,3,1) = ( gpcar(DEF_VECT,1,2,1) * gpcar(DEF_VECT,2,3,1) - gpcar(DEF_VECT,2,2,1) * gpcar(DEF_VECT,1,3,1)) * DENOMX
    xjaci(DEF_VECT,2,3,1) = (-gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,2,3,1) + gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,1,3,1)) * DENOMX

    gpcar(DEF_VECT,1,1,1) = -xjaci(DEF_VECT,1,1,1) - xjaci(DEF_VECT,2,1,1) - xjaci(DEF_VECT,3,1,1)
    gpcar(DEF_VECT,1,2,1) =  xjaci(DEF_VECT,1,1,1)
    gpcar(DEF_VECT,1,3,1) =  xjaci(DEF_VECT,2,1,1)
    gpcar(DEF_VECT,1,4,1) =  xjaci(DEF_VECT,3,1,1)
    gpcar(DEF_VECT,2,1,1) = -xjaci(DEF_VECT,1,2,1) - xjaci(DEF_VECT,2,2,1) - xjaci(DEF_VECT,3,2,1)
    gpcar(DEF_VECT,2,2,1) =  xjaci(DEF_VECT,1,2,1)
    gpcar(DEF_VECT,2,3,1) =  xjaci(DEF_VECT,2,2,1)
    gpcar(DEF_VECT,2,4,1) =  xjaci(DEF_VECT,3,2,1)
    gpcar(DEF_VECT,3,1,1) = -xjaci(DEF_VECT,1,3,1) - xjaci(DEF_VECT,2,3,1) - xjaci(DEF_VECT,3,3,1)
    gpcar(DEF_VECT,3,2,1) =  xjaci(DEF_VECT,1,3,1)
    gpcar(DEF_VECT,3,3,1) =  xjaci(DEF_VECT,2,3,1)
    gpcar(DEF_VECT,3,4,1) =  xjaci(DEF_VECT,3,3,1)

    do igaus = 2,pgaus
       gpdet(DEF_VECT,igaus) = gpdet(DEF_VECT,1)
       do idime = 1,3
          do jdime = 1,3
             xjaci(DEF_VECT,idime,jdime,igaus) = xjaci(DEF_VECT,idime,jdime,1)
          end do
          do inode = 1,4
             gpcar(DEF_VECT,idime,inode,igaus) = gpcar(DEF_VECT,idime,inode,1)
          end do
       end do
    end do

    do igaus = 1,pgaus
       gpvol(DEF_VECT,igaus) = weigp(igaus) * gpdet(DEF_VECT,igaus)
    end do
    !
    ! No Hessian
    !
    gphes(DEF_VECT,1:ntens,1:pnode,1:pgaus) = 0.0_rp


    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    gpden(DEF_VECT,:) = densi_aux
    gpvis(DEF_VECT,:) = visco_aux
    gpmut(DEF_VECT,:) = 0.0_rp

    gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)
    gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpmut(DEF_VECT,:)  ! Effective viscosity <= mu+mut

#else   ! NON SUPERFAST part missing

#ifdef OPENACCHHH
    call runend('nsi_element_operations_fast: non superfast version not ready for openacc case')
#endif

    call element_shape_function_derivatives_jacobian(&
         pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
         elmar(pelty) % deriv,elmar(pelty) % heslo,&
         elcod(DEF_VECT,:,:),gpvol(DEF_VECT,:),gpsha(DEF_VECT,:,:),gpder(DEF_VECT,:,:,:),gpcar(DEF_VECT,:,:,:),gphes(DEF_VECT,:,:,:),&
         list_elements(DEF_VECT))

    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    call ker_proper('DENSI','PGAUS',dummi,list_elements_p(DEF_VECT),gpden(DEF_VECT,:),  &
         pnode,pgaus,gpsha(DEF_VECT,:,:),gpcar(DEF_VECT,:,:,:))     ! rho
    call ker_proper('VISCO','PGAUS',dummi,list_elements_p(DEF_VECT),gpvis(DEF_VECT,:),  &
         pnode,pgaus,gpsha(DEF_VECT,:,:),gpcar(DEF_VECT,:,:,:))     ! mu
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p(DEF_VECT),gpmut(DEF_VECT,:),  &
         pnode,pgaus,gpsha(DEF_VECT,:,:),gpcar(DEF_VECT,:,:,:))     ! mut

    gpmut = gpden * gpmut
    gpvis = gpvis + gpmut  ! Effective viscosity <= mu+mut

#endif  ! SUPERFAST

    !----------------------------------------------------------------------
    !
    ! Gauss point values & eldtrho, elmurho
    !
    !----------------------------------------------------------------------
    gpgve(DEF_VECT,:,:,:)   = 0.0_rp
    gprhs(DEF_VECT,:,:)     = 0.0_rp
    gpvel(DEF_VECT,:,:) = 0.0_rp
    eldtrho(DEF_VECT,:) = 0.0_rp
    elmurho(DEF_VECT,:) = 0.0_rp

    do igaus = 1,pgaus
       FACT2X =  gpden(DEF_VECT,igaus)  * grnor_nsi
       do idime = 1,ndime
          do inode = 1,pnode
             gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) + elvel(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)

             do jdime = 1,ndime
                gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                     + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode)
             end do

          end do
          gpadv(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus)
          gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime)
       end do
       FACT2X = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
       FACT4X = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus)  )
       do inode = 1,pnode
          eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT2X
          elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT4X
       end do
    end do

    !----------------------------------------------------------------------
    !
    ! Element matrices
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    agrau(DEF_VECT,:,:)   = 0.0_rp
    wgrgr(DEF_VECT,:,:,:) = 0.0_rp

    !
    ! AGRAU = rho * (a.grad) Ni
    ! WGRGR = grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             agrau(DEF_VECT,inode,igaus) =  agrau(DEF_VECT,inode,igaus) + &
                  &                         gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
          end do
          agrau(DEF_VECT,inode,igaus) =  gpden(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus)
          do jnode = 1,pnode
             do idime = 1,ndime
                wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + &
                     &                              gpcar(DEF_VECT,idime,inode,igaus)*gpcar(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do
    end do

    !
    !  rho*(a.grad)Nj ) Ni + mu * grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus

       FACT6X = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode
          do idime = 1,ndime

             idofv           = (inode-1)*ndime+idime

             do jnode = 1,pnode
                jdofv           = (jnode-1)*ndime+idime
                FACT4X = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                FACT5X = FACT4X * ( agrau(DEF_VECT,jnode,igaus) ) &
                     &           + FACT6X *   wgrgr(DEF_VECT,inode,jnode,igaus)
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT5X
             end do
          end do
       end do
    end do


    do igaus = 1,pgaus

       FACT0X = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       FACT2X = 0.0_rp
       do idime = 1,ndime
          FACT2X = FACT2X + gpgve(DEF_VECT,idime,idime,igaus)
       end do
       FACT2X = FACT2X * FACT0X

       do inode = 1,pnode
          do idime = 1,ndime
             idofv = (inode-1)*ndime+idime
             do jnode = 1,pnode
                jdofv = (jnode-1)*ndime+idime
                !
                ! rho * (div u) u
                !
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                     + FACT2X * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                !
                ! rho * u.grad(u)^t
                !
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                        + FACT0X * gpsha(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) * gpvel(DEF_VECT,jdime,igaus)
                end do

             end do
             !
             ! 0.5*grad(u**2)
             !
             do jdime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     +  FACT0X * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus)*gpgve(DEF_VECT,idime,jdime,igaus)
             end do
          end do
          !
          ! ( mu*duj/dxi , dv/dxj ) (only div form)
          !
          do idime = 1,ndime
             idofv = (inode-1)*ndime + idime
             do jnode = 1,pnode
                FACT1X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                do jdime = 1,ndime
                   jdofv                       = (jnode-1)*ndime + jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT1X * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
          !
          ! bu = ( f , v )
          !
          FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus)
          end do
       end do
    end do
    !
    ! Send matrix to RHS
    !
    do jnode = 1,pnode
       do jdime = 1,ndime
          jevat = (jnode-1)*ndime+jdime
          do inode = 1,pnode
             do idime = 1,ndime
                ievat = (inode-1)*ndime+idime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     - elauu(DEF_VECT,ievat,jevat) * elvel(DEF_VECT,jdime,jnode)
             end do
          end do
       end do
    end do

    !
    ! Scatter element matrix to global one
    !
#ifndef OPENACCHHH
    call cputim(timeb)
    time1(7) = time1(7) + timeb - timea
    call cputim(timea)
    do ivect = 1,VECTOR_SIZE
#endif
       ielem = list_elements(ivect)
       if ( ielem > 0 ) then
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
             do idime = 1,ndime
                iauxi = idime + (ipoin-1) * ndime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
               !$acc atomic update
                rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect,idime,inode)

             end do
             !$acc atomic update
#ifdef NO_COLORING
             !$OMP ATOMIC
# endif
             dt_rho_nsi(ipoin) = dt_rho_nsi(ipoin) + eldtrho(ivect,inode)
             !$acc atomic update
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             mu_rho_nsi(ipoin) = mu_rho_nsi(ipoin) + elmurho(ivect,inode)

          end do
       end if
    end do
    !$acc end parallel loop
    !$acc wait
    !!$acc end data


#ifndef OPENACCHHH
    call cputim(timeb)
    time1(8) = time1(8) + timeb - timea
#endif

  end subroutine nsi_element_operations_fast






  subroutine nsi_element_operations_fast5(&
       VECTOR_DIM,pnode,pgaus,list_elements,time1)
!$acc routine(vecnor,frivel) seq
    use def_kintyp,            only : ip,rp                                              ! in_const_scalar
    use def_master,            only : rhsid                                              ! out       ! real(rp), pointer     :: rhsid(:)
    use def_master,            only : veloc                                              ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
    use def_domain,            only : coord                                              ! in_const  ! real(rp), pointer     :: coord(:,:)
    use def_domain,            only : ltype                                              ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
    use def_domain,            only : lnods                                              ! in_const  ! integer(ip), pointer  :: lnods(:,:)
    use def_domain,            only : lmate                                              ! in_const  ! integer(ip), pointer  :: lmate(:)
    use def_domain,            only : mnode                                              ! in_const scalar
    use def_domain,            only : elmar
    use def_domain,            only : elmda
    use def_domain,            only : ntens                                              ! in_const scalar
    use def_domain,            only : kfl_savda
    use def_domain,            only : elmda_gpvol
    use def_domain,            only : elmda_gpcar
    use def_domain,            only : lorde  
    use def_nastin,            only : dtinv_nsi                                          ! in_var   scalar
    use def_nastin,            only : dt_rho_nsi                                         ! out  ! real(rp), pointer     :: dt_rho_nsi(:)
    use def_nastin,            only : mu_rho_nsi                                         ! out  ! real(rp), pointer     :: mu_rho_nsi(:)
    use def_nastin,            only : grnor_nsi                                          ! in_const_scalar
    use def_nastin,            only : gravi_nsi                                          ! in_const  ! real(rp)         ::gravi_nsi(3)
    use def_nastin,            only : kfl_force_nsi 
    use def_nastin,            only : lforc_material_nsi
    use def_nastin,            only : xforc_material_nsi
    use mod_ker_proper,        only : ker_proper


    use def_kermod,            only : kfl_noslw_ker,avupo_ker,kfl_delta
    use def_kermod,            only : kfl_nswel_ker
    use def_kermod,            only : normal_nsw_ker
    use def_domain,            only : lpoty
    use def_domain,            only : ywale
    use mod_ker_proper,        only : delta_dom
    !
    ! include trucho
    !
    use def_kermod,                   only : lnsw_exch,kfl_rough,rough_dom,avta1_nsw_ker,fact_nsw_ker
    use def_domain,                   only : rough
    !
    !
    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in)          :: VECTOR_DIM                                       !< Number of nodes
    integer(ip), intent(in)          :: pnode                                            !< Number of nodes
    integer(ip), intent(in)          :: pgaus                                            !< Number of Gauss points
    integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                        !< List of elements

    real(rp),    intent(inout)       :: time1(10)                                        ! Timings
    !
    ! Element matrices and vectors (stiffness and preconditioner)
    !
    real(rp)                         :: elauu(VECTOR_DIM,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)                         :: elrbu(VECTOR_DIM,ndime,pnode)                    ! bu
    real(rp)                         :: eldtrho(VECTOR_DIM,pnode)                        ! Projection of rho/dt
    real(rp)                         :: elmurho(VECTOR_DIM,pnode)                        ! Projection of mu/rho
    !
    ! Gather
    !
    real(rp)                         :: elvel(VECTOR_DIM,ndime,pnode)                    ! u
    real(rp)                         :: elcod(VECTOR_DIM,ndime,pnode)                    ! x
    !
    ! Indices and dimensions
    !
    integer(ip)                      :: ielem,inode,ivect
    integer(ip)                      :: pelty,j,k,porde
    integer(ip)                      :: ipoin,igaus
    integer(ip)                      :: lnods_loc(VECTOR_DIM,pnode)
    integer(ip)                      :: list_elements_p(VECTOR_DIM)                      ! List of elements (always positive)
    !
    ! Gauss point values
    !
    real(rp)                         :: gpsha(VECTOR_DIM,pnode,pgaus)                    ! N
    real(rp)                         :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)              ! dN/dxi
    real(rp)                         :: gpvol(VECTOR_DIM,pgaus)                          ! w*|J|, |J|
    real(rp)                         :: gpvis(VECTOR_DIM,pgaus)                          ! Viscosity
    real(rp)                         :: gpvis_nsw(VECTOR_SIZE,pgaus)                     ! Viscosity for no slip wall
    real(rp)                         :: gpmut(VECTOR_DIM,pgaus)                          ! mut
    real(rp)                         :: gpden(VECTOR_DIM,pgaus)                          ! Density
    real(rp)                         :: gpadv(VECTOR_DIM,ndime,pgaus)                    ! u+u'
    real(rp)                         :: gprhs(VECTOR_DIM,ndime,pgaus)                    ! RHS
    real(rp)                         :: gpvel(VECTOR_DIM,ndime,pgaus)                    ! u
    real(rp)                         :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)              ! grad(u)
    real(rp)                         :: xjacm(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: xjaci(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: gpdet(VECTOR_DIM,pgaus)
    real(rp)                         :: wgrgr(VECTOR_DIM,pnode,pnode,pgaus)
    real(rp)                         :: agrau(VECTOR_DIM,pnode,pgaus)
    real(rp)                         :: T_dtrho(VECTOR_DIM) 
    real(rp)                         :: d_dtrho(VECTOR_DIM) 
    real(rp)                         :: T_murho(VECTOR_DIM) 
    real(rp)                         :: d_murho(VECTOR_DIM) 
    real(rp)                         :: timea,timeb
    integer(ip)                      :: idime,pmate
    integer(ip)                      :: jdime,jnode,idofv,jdofv
    integer(ip)                      :: ievat,jevat,iauxi,dummi
    !
    ! Gather No slip wall law
    !
    real(rp)    :: elibopo(VECTOR_SIZE,pnode)
    real(rp)    :: elnnsw(VECTOR_SIZE,ndime)
    real(rp)    :: elavv(VECTOR_SIZE,ndime,pnode)
    real(rp)    :: elywal(VECTOR_SIZE)
    !
    ! Internal
    !
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#ifdef VECTOR_SIZE_VARIABLE
#define DEF_VECT 1:VECTOR_DIM
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
#endif

#ifdef OPENACCHHH

#define FACT0X     fact0
#define FACT1X     fact1
#define FACT2X     fact2
#define FACT4X     fact4
#define FACT5X     fact5
#define FACT6X     fact6
#define DTINV_LOCX dtinv_loc
 #define T1X        t1
#define T2X        t2
#define T3X        t3
#define DENOMX     denom

#else

#ifdef VECTOR_SIZE_VARIABLE
#define FACT0X     fact0(1:VECTOR_DIM)
#define FACT1X     fact1(1:VECTOR_DIM)
#define FACT2X     fact2(1:VECTOR_DIM)
#define FACT4X     fact4(1:VECTOR_DIM)
#define FACT5X     fact5(1:VECTOR_DIM)
#define FACT6X     fact6(1:VECTOR_DIM)
#define DTINV_LOCX dtinv_loc(1:VECTOR_DIM)
#define T1X        t1(1:VECTOR_DIM)
#define T2X        t2(1:VECTOR_DIM)
#define T3X        t3(1:VECTOR_DIM)
#define DENOMX     denom(1:VECTOR_DIM)
#else
#define FACT0X     fact0(1:VECTOR_SIZE)
#define FACT1X     fact1(1:VECTOR_SIZE)
#define FACT2X     fact2(1:VECTOR_SIZE)
#define FACT4X     fact4(1:VECTOR_SIZE)
#define FACT5X     fact5(1:VECTOR_SIZE)
#define FACT6X     fact6(1:VECTOR_SIZE)
#define DTINV_LOCX dtinv_loc(1:VECTOR_SIZE)
#define T1X        t1(1:VECTOR_SIZE)
#define T2X        t2(1:VECTOR_SIZE)
#define T3X        t3(1:VECTOR_SIZE)
#define DENOMX     denom(1:VECTOR_SIZE)
#endif
    
#endif

    real(rp)    :: FACT0X
    real(rp)    :: FACT1X
    real(rp)    :: FACT2X
    real(rp)    :: FACT4X
    real(rp)    :: FACT5X
    real(rp)    :: FACT6X
    real(rp)    :: DTINV_LOCX
    real(rp)    :: T1X
    real(rp)    :: T2X
    real(rp)    :: T3X
    real(rp)    :: DENOMX
    !
    ! include
    !
    real(rp)         :: gpgnavv(VECTOR_SIZE,ndime,pgaus)  ! Gauss Point Gradient in Normal dir of AVerage Velocity
    real(rp)         :: elgnavv(VECTOR_SIZE,ndime)        ! ELement Gradient in Normal dir of AVerage Velocity
    real(rp)         :: elgnavvt(VECTOR_SIZE)             ! ELement Gradient in Normal dir of AVerage Tangent Velocity

    real(rp)         :: auxvi(VECTOR_SIZE)
    real(rp)         :: tveno_aux(VECTOR_SIZE)
    real(rp)         :: auxde(VECTOR_SIZE)
    real(rp)         :: auxmut(VECTOR_SIZE)
    real(rp)         :: avelavv(VECTOR_SIZE,ndime)        ! AVerage ELement AVerage Velocity
    real(rp)         :: avta1_aux(VECTOR_SIZE,ndime)
    real(rp)         :: auxi(VECTOR_SIZE)
    real(rp)         :: velfr(VECTOR_SIZE)
    real(rp)         :: fact_aux(VECTOR_SIZE)

    real(rp)         :: avtan_fric_grad_based(VECTOR_SIZE)    ! just the magnitude
    real(rp)         :: av_mu_mut(VECTOR_SIZE)



    real(rp)         :: rough_aux,vikin,tveno,auxi2,kount1
    integer(ip)      :: kk,ibopo
    integer(ip),parameter      :: imethod = 4   !now the default will be method 4 that is teh one taht is working best

    !end include
    
    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------

    call cputim(timea)
    
    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    porde = lorde(pelty)
    pmate = lmate(ielem)
    
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          list_elements_p(ivect)   = list_elements(ivect)
       else
          list_elements_p(ivect)   = list_elements(1)
       end if
    end do
    do ivect = 1,VECTOR_DIM      
       ielem = list_elements_p(ivect)
       gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
    end do

    call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! rho
    call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mut



!hh    ---------------
!hh        call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,gpcar)     ! mut  !OJO en el de fast no le apsa gpcar se arreglaar de alguna manera

!hh    if (kfl_cotur_nsi < 0) gpmut = gpden * gpmut                                             ! smago, wale, vreman etc   ! esto esta hecho mas abajo como gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)
!yo creo que se puede subir aca es mas no veo paar que le ponene explitas las dimensiones  
!  bueno asi DEF_VECT es ivect paar open acc y eso lo vectoriza la GPU no creoq ue haga mayor diferencia comparado con el coste de llamar a ker_proper


!hh    if ( kfl_noslw_ker /= 0_ip) call ker_nsw_visc(ndime,pnode,pgaus,list_elements,elavv,gpcar,elnnsw,gpvis,&    ! ojo tiene que ir mas abaj porque necesita elnnsw elavv elibop elywal
    !hh         gpden,gpmut,elywal,gpvis_nsw)  ! obtains gpvis_nsw
    ! con lo cual se va a la parte vectorizada  o sea que deberia meterlo con un include
    ! como hizo guillaume paar frivel en bound operatios fast  -- supongo que habra que usar lo de include 



    
    
!hh    gpvis = gpvis + gpmut                                                                    ! Effective viscosity <= mu+mut
!    ------------------



    
    
    DTINV_LOCX = dtinv_nsi

    call cputim(timeb)
    time1(3) = time1(3) + timeb - timea
    call cputim(timea)
    
    !--------------------------------------------------------------------
    !
    ! Shape function and derivatives   !ver si gpvis_nsw tine que entrar aca
    ! added avelavv  becaise it was giving me error at run time not sure if it is correct
    ! tambien avupo_ker,elavv
    ! Aunque estoy corriendo un caso sin no slip walllaw y esas solo aparcen dentro de if las quiere igual
    !
    !--------------------------------------------------------------------
    
    !$acc data create( agrau   , lnods_loc , gpvol     , elauu ,          &   
    !$acc              eldtrho , elvel     , gpmut     , gpadv , gpvel ,  &
    !$acc              gpgve   , elcod     , gpsha     , gpcar ,          &
    !$acc              gpden   , wgrgr     , gpdet     ,                  &
    !$acc              xjacm   , xjaci     ,                              &
    !$acc              elmar,                                             &
    !$acc              avelavv,avupo_ker,elavv,elnnsw,av_mu_mut,gpvis_nsw,&
    !$acc              avtan_fric_grad_based,elgnavvt,auxi,auxvi,auxde,   &
    !$acc              auxmut,elgnavv,gpgnavv,velfr,elywal,elibopo,       &
    !$acc              lpoty,avta1_aux,                                   &
    
    !$acc              elmar(pelty)%shape,                                &
    !$acc              elmar(pelty)%deriv,                                &
    !$acc              elmar(pelty)%weigp,                                &
    !$acc              T_dtrho, d_dtrho, T_murho, d_murho,                &
    !$acc              gprhs   , elmurho   , elrbu     , gpvis )          &
    !$acc copyin(      list_elements        ,                             &
    !$acc              xforc_material_nsi,lforc_material_nsi,             &
    !$acc              gpsha ,                                            &
    !$acc              gpden ,                                            &
    !$acc              gpvis ,                                            &
    !$acc              gpmut ,                                            &
    !$acc              elmar,                                             &
    !$acc              elmar(pelty)%shape ,                               &
    !$acc              elmar(pelty)%deriv ,                               &
    !$acc              elmar(pelty)%weigp                                 )
    !
    ! This do starts here both in the openacc version and in the nonopenacc
    ! for the non opencc it ends 30  lines later,
    ! in the openacc case it covers all the subroutine.
    ! Similarly the scatter in the nonopenacc case needs a do ivect
    !    
    !$acc parallel loop gang vector default(present)
    !
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
          ielem                    = list_elements(ivect)
       else
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,list_elements(1))
          ielem                    = list_elements(1)
       end if
       !
       ! Transient
       !
       do inode = 1,pnode
          ipoin = lnods_loc(ivect,inode)
          do idime = 1,ndime
             elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
             elcod(ivect,idime,inode) = coord(idime,ipoin)
          end do
       end do
       !
       ! no slip wall law - I could add flag so that i does it only on those elements where it is needed kfl_nswel_ker(ielem)
       !
       if ( kfl_noslw_ker /= 0_ip ) then
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
             if (lpoty(ipoin) > 0 ) then
                elibopo(ivect,inode) = 1.0_rp
             else
                elibopo(ivect,inode) = 0.0_rp
             end if
             do idime = 1,ndime
                elavv(ivect,idime,inode) = avupo_ker(idime,ipoin)
             end do
          end do

          if (kfl_nswel_ker(ielem) > 0_ip ) then
             do idime = 1,ndime
                elnnsw(ivect,idime) = normal_nsw_ker(idime,kfl_nswel_ker(ielem)) 
             end do
          else
             elnnsw(ivect,:) = 0.0_rp
          end if
          if (kfl_delta == 0) then
             elywal(ivect) = delta_dom
          else
             elywal(ivect) = ywale(ielem)
          end if
          !
          ! This comes directly in lnsw_exch(ielem)%velav(1:ndime)  - obtained in nsi_wallav
          ! for the moment it is not vectorized - not sure how easy/practical it would be
          ! these line where previously 250 forward - I belive it is better to do it here I leave them theer commented out just for the case
          ! when I want to get that out into a subroutine
          if (ielem /= 0_ip) then
             avelavv(ivect,1:ndime) = lnsw_exch(ielem)%velav(1:ndime)
          else
             avelavv(ivect,1:ndime) = 0.0_rp
          end if
       end if
       
       
#ifndef OPENACCHHH
    end do
    call cputim(timeb)
    time1(1) = time1(1) + timeb - timea
    call cputim(timea)
#endif
    !
    ! Why not trying this which can be vectorized easily
    !
    !do inode = 1,pnode
    !   do ivect = 1,VECTOR_DIM                      
    !      ipoin = lnods_loc(ivect,inode)
    !      elvel(ivect,1:ndime,inode) = veloc(1:ndime,ipoin,1)
    !      elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)          
    !   end do
    !end do
    
    !--------------------------------------------------------------------
    !
    ! Element Cartesian derivatives and Jacobian: GPCAR, GPVOL
    !
    !--------------------------------------------------------------------

    if( kfl_savda == 2 ) then
#ifndef OPENACCHHH
       do ivect = 1,VECTOR_DIM      
#endif
          ielem = abs(list_elements(ivect))
          if( ielem > 0 ) then
             gpvol(ivect,1:pgaus)                 = elmda_gpvol(1:pgaus,ielem)
             gpcar(ivect,1:ndime,1:mnode,1:pgaus) = elmda_gpcar(1:ndime,1:mnode,1:pgaus,ielem)
          else
             gpvol(ivect,1:pgaus)                 = 0.0_rp
             gpcar(ivect,1:ndime,1:mnode,1:pgaus) = 0.0_rp
          end if
#ifndef OPENACCHHH
       end do
#endif
    else 

       if(      ndime == 2 ) then
          
          do igaus = 1,pgaus
             xjacm(DEF_VECT,1,1)    =  0.0_rp
             xjacm(DEF_VECT,1,2)    =  0.0_rp
             xjacm(DEF_VECT,2,1)    =  0.0_rp
             xjacm(DEF_VECT,2,2)    =  0.0_rp
             do k = 1,pnode
                xjacm(DEF_VECT,1,1) =  xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,1,2) =  xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,2,1) =  xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,2,2) =  xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(2,k,igaus)
             end do
             gpdet(DEF_VECT,igaus)  =  xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)
             denom                  =  1.0_rp / gpdet(DEF_VECT,igaus)
             xjaci(DEF_VECT,1,1)    =  xjacm(DEF_VECT,2,2) * denom
             xjaci(DEF_VECT,2,2)    =  xjacm(DEF_VECT,1,1) * denom
             xjaci(DEF_VECT,2,1)    = -xjacm(DEF_VECT,2,1) * denom
             xjaci(DEF_VECT,1,2)    = -xjacm(DEF_VECT,1,2) * denom
             do j = 1, pnode
                gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,1) * elmar(pelty) % deriv(2,j,igaus)
                
                gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,2) * elmar(pelty) % deriv(2,j,igaus)
             end do
             gpvol(DEF_VECT,igaus) = elmar(pelty) % weigp(igaus) * gpdet(DEF_VECT,igaus)
          end do
       
       else if( ndime == 3 ) then
          
          do igaus = 1,pgaus
             xjacm(DEF_VECT,1:3,1:3)   = 0.0_rp
             do k = 1,pnode 
                xjacm(DEF_VECT,1,1)    =  xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,1,2)    =  xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,1,3)    =  xjacm(DEF_VECT,1,3) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(3,k,igaus)
                xjacm(DEF_VECT,2,1)    =  xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,2,2)    =  xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,2,3)    =  xjacm(DEF_VECT,2,3) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(3,k,igaus)
                xjacm(DEF_VECT,3,1)    =  xjacm(DEF_VECT,3,1) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,3,2)    =  xjacm(DEF_VECT,3,2) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,3,3)    =  xjacm(DEF_VECT,3,3) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(3,k,igaus)
             end do
             t1                        =  xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,2,3)
             t2                        = -xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,3)
             t3                        =  xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,2) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,2)
             gpdet(DEF_VECT,igaus)     =  xjacm(DEF_VECT,1,1) * t1 + xjacm(DEF_VECT,1,2) * t2 + xjacm(DEF_VECT,1,3) * t3
             denom                     =  1.0_rp / gpdet(DEF_VECT,igaus)
             xjaci(DEF_VECT,1,1)       =  t1 * denom
             xjaci(DEF_VECT,2,1)       =  t2 * denom
             xjaci(DEF_VECT,3,1)       =  t3 * denom
             xjaci(DEF_VECT,2,2)       =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,1,3)) * denom
             xjaci(DEF_VECT,3,2)       =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,2) + xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,1)) * denom
             xjaci(DEF_VECT,3,3)       =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)) * denom
             xjaci(DEF_VECT,1,2)       =  (-xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,1,3)) * denom
             xjaci(DEF_VECT,1,3)       =  ( xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,2,3) - xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,1,3)) * denom
             xjaci(DEF_VECT,2,3)       =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,3) + xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,3)) * denom
             
             do j = 1, pnode
                gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,1) * elmar(pelty) % deriv(2,j,igaus) &
                     &                      + xjaci(DEF_VECT,3,1) * elmar(pelty) % deriv(3,j,igaus)
                
                gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,2) * elmar(pelty) % deriv(2,j,igaus) &
                     &                      + xjaci(DEF_VECT,3,2) * elmar(pelty) % deriv(3,j,igaus)
                
                gpcar(DEF_VECT,3,j,igaus) =   xjaci(DEF_VECT,1,3) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,3) * elmar(pelty) % deriv(2,j,igaus) &
                     &                      + xjaci(DEF_VECT,3,3) * elmar(pelty) % deriv(3,j,igaus)
             end do
             gpvol(DEF_VECT,igaus) = elmar(pelty) % weigp(igaus) * gpdet(DEF_VECT,igaus)
          end do
       end if


    end if

    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)
    !
    if ( kfl_noslw_ker /= 0_ip ) then   ! I put it her becausa in the non fast version it goes before gpvis=gpvis+gpmut
       !call ker_nsw_visc(ndime,pnode,pgaus,list_elements,elavv,gpcar,elnnsw,gpvis,gpden,gpmut,elywal,gpvis_nsw)
       ! no de deja con include porque tiene ifdefs   asi que lo meto a mano aca  -- muy feo
          !----------------------------------------------------------------------
          !
          ! Gauss point values
          !
          !----------------------------------------------------------------------
          !
          ! GPGNAVV = dj uav_i nj =  dj N_I_i nj  Uav_I  ! gauss point gradient in normal direction of the average velocity
          !
          gpgnavv(DEF_VECT,:,:)   = 0.0_rp

          do igaus = 1,pgaus
             do inode = 1,pnode
                do idime = 1,ndime
                   do jdime = 1,ndime
                      gpgnavv(DEF_VECT,idime,igaus) = gpgnavv(DEF_VECT,idime,igaus) + elavv(DEF_VECT,idime,inode) * &
                           gpcar(DEF_VECT,jdime,inode,igaus) * elnnsw(DEF_VECT,jdime)
                   end do
                end do
             end do
          end do
          !
          ! Obtain the average over gauss points of the normal gradient of the average velocity for all gauss points
          ! I could merge this loop with the upper one and obtain elgnavv directly
          ! Also average density and viscosity
          !
          elgnavv  = 0.0_rp
          auxvi = 0.0_rp
          auxde = 0.0_rp
          auxmut = 0.0_rp
          do igaus = 1,pgaus
             do idime = 1,ndime
                elgnavv(DEF_VECT,idime)  = elgnavv(DEF_VECT,idime)  + gpgnavv(DEF_VECT,idime,igaus)
             end do
             auxvi(DEF_VECT)  = auxvi(DEF_VECT)  + gpvis(DEF_VECT,igaus)
             auxde(DEF_VECT)  = auxde(DEF_VECT)  + gpden(DEF_VECT,igaus)
             auxmut(DEF_VECT) = auxmut(DEF_VECT) + gpmut(DEF_VECT,igaus)
          end do
          do idime = 1,ndime
             elgnavv(DEF_VECT,idime) = elgnavv(DEF_VECT,idime) / real(pgaus,rp)
          end do
          auxvi(DEF_VECT)  = auxvi(DEF_VECT)  / real(pgaus,rp)
          auxde(DEF_VECT)  = auxde(DEF_VECT)  / real(pgaus,rp)
          auxmut(DEF_VECT) = auxmut(DEF_VECT) / real(pgaus,rp)
          !
          ! Substract normal component to keep only tangential one
          !
          !
          auxi(DEF_VECT) = 0.0_rp
          do idime=1,ndime
             auxi(DEF_VECT) = auxi(DEF_VECT) + elgnavv(DEF_VECT,idime) * elnnsw(DEF_VECT,idime)
          end do
          do idime=1,ndime
             elgnavv(DEF_VECT,idime) = elgnavv(DEF_VECT,idime) - auxi(DEF_VECT) * elnnsw(DEF_VECT,idime)
          end do
          !
          ! obtain modulus of tangential component
          !
          auxi(DEF_VECT) = 0.0_rp
          do idime=1,ndime
             auxi(DEF_VECT) = auxi(DEF_VECT) + elgnavv(DEF_VECT,idime) * elgnavv(DEF_VECT,idime)
          end do
          elgnavvt(DEF_VECT) = sqrt (auxi(DEF_VECT))


          !This comes directly in lnsw_exch(ielem)%velav(1:ndime)  - obtained in nsi_wallav
          ! for the moment it is not vectorized - not sure how easy/practical it would be
          !          do ivect = 1,VECTOR_SIZE
          !             ielem = list_elements(ivect)
          !             if (ielem/=0_ip) then
          !                avelavv(ivect,1:ndime) = lnsw_exch(ielem)%velav(1:ndime)
          !             else
          !                avelavv(ivect,1:ndime) = 0.0_rp
          !             end if
          !          end do
          !
          ! Substract normal component to keep only tangential one
          !
          !
          auxi(DEF_VECT) = 0.0_rp
          do idime=1,ndime
             auxi(DEF_VECT) = auxi(DEF_VECT) + avelavv(DEF_VECT,idime) * elnnsw(DEF_VECT,idime)
          end do
          do idime=1,ndime
             avelavv(DEF_VECT,idime) = avelavv(DEF_VECT,idime) - auxi(DEF_VECT) * elnnsw(DEF_VECT,idime)
          end do

          !
          ! All of the operations up to he are done for all elements when actually they could be done only for kfl_nswel_ker(ielem) >0
          ! RETHINK!!!
          !
          !
          ! Obtain tange from wall_law.
          ! Compute U*: VELFR   this part is not vectorized for the moment  I would need a vector frivel , not difficult
          ! also Time average of (mu+mut) d u_t / dn  - to be usad later
          !
#ifndef OPENACCHHH
          do ivect = 1,VECTOR_SIZE
#endif             
          ielem = list_elements(ivect)
          if( ielem > 0_ip ) then
             if(kfl_nswel_ker(ielem) >0) then
                vikin = auxvi(ivect) / auxde(ivect)                           ! nu
                call vecnor(avelavv(ivect,1:ndime),ndime,tveno,2_ip)          ! |u_tan-u_fix_tan|
                if(tveno > 1.0e-20_rp) then
                   tveno_aux(ivect) = tveno    ! I save it for later usage
                else
                   tveno_aux(ivect) = 1.0_rp
                end if

                if( kfl_rough == 0 ) then
                   !
                   ! Constant roughness
                   !
                   rough_aux = rough_dom

                else if( kfl_rough > 0 ) then
                   rough_aux = 0.0_rp
                   kount1 = 0_ip
                   do inode = 1,pnode
                      ipoin =  lnods(inode,ielem)
                      ibopo = lpoty(ipoin)
                      if (ibopo /= 0) then
                         kount1 = kount1 + 1_ip
                         rough_aux = rough_aux + rough(ipoin)
                      end if
                   end do
                   rough_aux = rough_aux / kount1
                end if

                call frivel(elywal(ivect),rough_aux,tveno,vikin,velfr(ivect))      ! U*
                !
                ! Time average of (mu+mut) d u_t / dn  - to be usad later
                !
                avta1_aux(ivect,:) = avta1_nsw_ker(:,kfl_nswel_ker(ielem))
                !
                ! fact_aux to be usad later
                !
                kk = 0_ip
                auxi2 = 0.0_rp
                do inode=1,pnode
                   ipoin = lnods(inode,ielem)
                   if (fact_nsw_ker(ipoin) > 0.0_rp) then    ! boundary node
                      kk = kk + 1_ip
                      auxi2 = auxi2 + fact_nsw_ker(ipoin)
                   end if
                   if (kk > 0_ip ) then
                      fact_aux(ivect) = auxi2 / real(kk,rp)
                   else
                      fact_aux(ivect) = 1.0_rp
                   end if
                end do

             else   ! Element is not boundary element
                velfr(ivect) = 0.0_rp
                tveno_aux(ivect) = 1.0_rp   ! just some value so that it does not give problems
                avta1_aux(ivect,:) = 0.0_rp
             end if
          else   ! Element number is null
             velfr(ivect) = 0.0_rp
             tveno_aux(ivect) = 1.0_rp      ! just some value so that it does not give problems
             avta1_aux(ivect,:) = 0.0_rp
          end if
#ifndef OPENACCHHH
       end do
#endif             


          !
          ! Tangential force that comes from the friction velocity - gradient based becausa it has already been transformed using fact
          !
          avtan_fric_grad_based(DEF_VECT) = velfr(DEF_VECT) * velfr(DEF_VECT) * auxde(DEF_VECT)
          gpvis_nsw(DEF_VECT,1) = ( avtan_fric_grad_based(DEF_VECT) / (abs(elgnavvt(DEF_VECT))+1.0e-30_rp) ) !!!- (auxmut(DEF_VECT) + auxvi(DEF_VECT))
          av_mu_mut(DEF_VECT) = avta1_aux(DEF_VECT,1) * avta1_aux(DEF_VECT,1) +  avta1_aux(DEF_VECT,2) * avta1_aux(DEF_VECT,2)
          if (ndime==3) av_mu_mut(DEF_VECT) = av_mu_mut(DEF_VECT) +  avta1_aux(DEF_VECT,3) * avta1_aux(DEF_VECT,3)
          av_mu_mut(DEF_VECT) = sqrt(av_mu_mut(DEF_VECT)) / (abs(elgnavvt(DEF_VECT))+1.0e-30_rp)
          gpvis_nsw(DEF_VECT,1) = gpvis_nsw(DEF_VECT,1) - av_mu_mut(DEF_VECT)
          do igaus=2,pgaus
             gpvis_nsw(DEF_VECT,igaus) = max ( - gpmut(DEF_VECT,igaus)  , gpvis_nsw(DEF_VECT,1) )
          end do
          gpvis_nsw(DEF_VECT,1) = max ( - gpmut(DEF_VECT,1)  , gpvis_nsw(DEF_VECT,1) )

    end if
    !
    !
    !

    
    gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpmut(DEF_VECT,:)  ! Effective viscosity <= mu+mut
    
    !----------------------------------------------------------------------
    !
    ! Gauss point values & eldtrho, elmurho
    !
    !----------------------------------------------------------------------

    gpgve(DEF_VECT,:,:,:) = 0.0_rp
    gprhs(DEF_VECT,:,:)   = 0.0_rp
    gpvel(DEF_VECT,:,:)   = 0.0_rp
    eldtrho(DEF_VECT,:)   = 0.0_rp
    elmurho(DEF_VECT,:)   = 0.0_rp

    do igaus = 1,pgaus
       FACT2X =  gpden(DEF_VECT,igaus) * grnor_nsi
       do idime = 1,ndime
          do inode = 1,pnode
             gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) + elvel(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)
             do jdime = 1,ndime
                gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                     + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode)
             end do
          end do
          gpadv(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus)
          gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime)
       end do
    end do
    !
    ! Material force
    !
    if( kfl_force_nsi == 1 ) then
       if( lforc_material_nsi(pmate) == 2 ) then
          do igaus = 1,pgaus
             do idime = 1,ndime
                gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) +  xforc_material_nsi(idime,pmate)
             end do
          end do          
       end if
    end if
 !
 ! Tau and Tim
 !
    if( porde == 1 ) then
       do igaus = 1,pgaus
          FACT2X = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
          FACT4X = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / gpden(DEF_VECT,igaus) 
          do inode = 1,pnode
             eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT2X
             elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT4X
          end do
       end do
    else

       T_dtrho(DEF_VECT) = 0.0_rp
       d_dtrho(DEF_VECT) = 0.0_rp
       T_murho(DEF_VECT) = 0.0_rp
       d_murho(DEF_VECT) = 0.0_rp
       do igaus = 1,pgaus
          FACT2X                 = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
          FACT4X                 = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / gpden(DEF_VECT,igaus)


          T_dtrho(DEF_VECT) = T_dtrho(DEF_VECT) + FACT2X
          T_murho(DEF_VECT) = T_murho(DEF_VECT) + FACT4X
          do inode = 1,pnode
             eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus)**2 * FACT2X
             elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus)**2 * FACT4X
          end do
       end do

       do inode = 1,pnode
          d_dtrho(DEF_VECT) = d_dtrho(DEF_VECT) + eldtrho(DEF_VECT,inode)
          d_murho(DEF_VECT) = d_murho(DEF_VECT) + elmurho(DEF_VECT,inode)
       end do

       do inode = 1,pnode
          eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) * T_dtrho(DEF_VECT) / d_dtrho(DEF_VECT)
          elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) * T_murho(DEF_VECT) / d_murho(DEF_VECT)
       end do
  

 end if

    !----------------------------------------------------------------------
    !
    ! Element matrices 
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    agrau(DEF_VECT,:,:)   = 0.0_rp
    wgrgr(DEF_VECT,:,:,:) = 0.0_rp
    !
    ! AGRAU = rho * (a.grad) Ni
    ! WGRGR = grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             agrau(DEF_VECT,inode,igaus) =  agrau(DEF_VECT,inode,igaus) + &
                  &                         gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
          end do
          agrau(DEF_VECT,inode,igaus) =  gpden(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus)
          do jnode = 1,pnode
             do idime = 1,ndime
                wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + &
                     &                              gpcar(DEF_VECT,idime,inode,igaus)*gpcar(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do
    end do
    !
    !  rho*(a.grad)Nj ) Ni + mu * grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus

       FACT6X = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode
          do idime = 1,ndime

             idofv           = (inode-1)*ndime+idime

             do jnode = 1,pnode
                jdofv           = (jnode-1)*ndime+idime
                FACT4X = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                FACT5X = FACT4X * ( agrau(DEF_VECT,jnode,igaus) ) &
                     &           + FACT6X *   wgrgr(DEF_VECT,inode,jnode,igaus)
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT5X
             end do
          end do
       end do
    end do

    do igaus = 1,pgaus

       FACT0X = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       FACT2X = 0.0_rp
       do idime = 1,ndime
          FACT2X = FACT2X + gpgve(DEF_VECT,idime,idime,igaus)
       end do
       FACT2X = FACT2X * FACT0X

       do inode = 1,pnode
          do idime = 1,ndime
             idofv = (inode-1)*ndime+idime
             do jnode = 1,pnode
                jdofv = (jnode-1)*ndime+idime
                !
                ! rho * (div u) u
                !
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                     + FACT2X * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                !
                ! rho * u.grad(u)^t
                ! 
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                        + FACT0X * gpsha(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) * gpvel(DEF_VECT,jdime,igaus)
                end do

             end do
             ! 
             ! 0.5 * grad(u**2)
             !
             do jdime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     +  FACT0X * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus)*gpgve(DEF_VECT,idime,jdime,igaus)
             end do
          end do
          !
          ! ( mu*duj/dxi , dv/dxj ) (only div form)
          !      
          do idime = 1,ndime
             idofv = (inode-1)*ndime + idime
             do jnode = 1,pnode
                FACT1X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                do jdime = 1,ndime
                   jdofv                       = (jnode-1)*ndime + jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT1X * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
          !
          ! bu = ( f , v ) 
          !
          FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus) 
          end do
       end do
    end do
    !
    ! Send matrix to RHS
    !
    do jnode = 1,pnode
       do jdime = 1,ndime
          jevat = (jnode-1)*ndime+jdime
          do inode = 1,pnode
             do idime = 1,ndime
                ievat = (inode-1)*ndime+idime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     - elauu(DEF_VECT,ievat,jevat) * elvel(DEF_VECT,jdime,jnode) 
             end do
          end do
       end do
    end do

#ifndef OPENACCHHH
    if ( kfl_noslw_ker /= 0_ip ) then 
       elauu(DEF_VECT,:,:)   = 0.0_rp   ! I can reusa elauu no need for elauu_nsw
       !
       ! WGRGR = grad(Ni) . grad(Nj)  Has already been calculated i do not need to recalculate it
       !
       !
       !   mu_nsw * grad(Ni) . grad(Nj)
       !
       do igaus = 1,pgaus

          FACT6X = gpvis_nsw(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

          do inode = 1,pnode
             do idime = 1,ndime
                idofv           = (inode-1)*ndime+idime
                do jnode = 1,pnode
                   jdofv           = (jnode-1)*ndime+idime
                   FACT5X =  FACT6X *   wgrgr(DEF_VECT,inode,jnode,igaus)
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT5X
                end do
             end do
          end do
       end do    ! no veo razon para tener 2 loops iagus ver como lo tenia yo en fast y si tiene alguna razon          
       do igaus = 1,pgaus
          !
          ! ( mu_nsw*duj/dxi , dv/dxj ) (only div form)
          !
          do inode = 1,pnode
             do idime = 1,ndime
                idofv = (inode-1)*ndime + idime
                do jnode = 1,pnode
                   FACT1X = gpvis_nsw(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                   do jdime = 1,ndime
                      jdofv                       = (jnode-1)*ndime + jdime
                      elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT1X * gpcar(DEF_VECT,jdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end do
       !
       ! Send matrix to RHS
       !
       do jnode = 1,pnode
          do jdime = 1,ndime
             jevat = (jnode-1)*ndime+jdime
             do inode = 1,pnode
                do idime = 1,ndime
                   ievat = (inode-1)*ndime+idime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                        - elauu(DEF_VECT,ievat,jevat) * elavv(DEF_VECT,jdime,jnode)
                end do
             end do
          end do
       end do
    end if
#endif
    !
    ! Scatter element matrix to global one 
    ! 
#ifndef OPENACCHHH
    call cputim(timeb)
    time1(7) = time1(7) + timeb - timea
    call cputim(timea)
    do ivect = 1,VECTOR_DIM                        
#endif 
       ielem = list_elements(ivect)
       if ( ielem > 0 ) then
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
             do idime = 1,ndime
                iauxi = idime + (ipoin-1) * ndime 
                !$acc atomic update
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect,idime,inode)

             end do
             !$acc atomic update
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             dt_rho_nsi(ipoin) = dt_rho_nsi(ipoin) + eldtrho(ivect,inode)
             !$acc atomic update
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             mu_rho_nsi(ipoin) = mu_rho_nsi(ipoin) + elmurho(ivect,inode)

          end do
       end if
    end do
    !$acc end parallel loop
    !$acc end data


#ifndef OPENACCHHH
    call cputim(timeb)
    time1(8) = time1(8) + timeb - timea                       
#endif

  end subroutine nsi_element_operations_fast5










    subroutine nsi_element_operations_fast6(&
       VECTOR_DIM,pnode,pgaus,list_elements,time1)

    use def_kintyp,            only : ip,rp                                              ! in_const_scalar
    use def_master,            only : rhsid                                              ! out       ! real(rp), pointer     :: rhsid(:)
    use def_master,            only : veloc                                              ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
    use def_domain,            only : coord                                              ! in_const  ! real(rp), pointer     :: coord(:,:)
    use def_domain,            only : ltype                                              ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
    use def_domain,            only : lnods                                              ! in_const  ! integer(ip), pointer  :: lnods(:,:)
    use def_domain,            only : mnode                                              ! in_const scalar
    use def_domain,            only : elmar
    use def_domain,            only : elmda
    use def_domain,            only : ntens                                              ! in_const scalar
    use def_domain,            only : kfl_savda
    use def_domain,            only : elmda_gpvol
    use def_domain,            only : elmda_gpcar
    use def_nastin,            only : dtinv_nsi                                          ! in_var   scalar
    use def_nastin,            only : dt_rho_nsi                                         ! out  ! real(rp), pointer     :: dt_rho_nsi(:)
    use def_nastin,            only : mu_rho_nsi                                         ! out  ! real(rp), pointer     :: mu_rho_nsi(:)
    use def_nastin,            only : grnor_nsi                                          ! in_const_scalar
    use def_nastin,            only : gravi_nsi                                          ! in_const  ! real(rp)         ::gravi_nsi(3)
    use mod_ker_proper,        only : ker_proper
    use def_kermod,            only : densi_ker 
    use def_kermod,            only : visco_ker 
    use def_kermod,            only : turmu_ker 
    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in)          :: VECTOR_DIM                                       !< Number of nodes
    integer(ip), intent(in)          :: pnode                                            !< Number of nodes
    integer(ip), intent(in)          :: pgaus                                            !< Number of Gauss points
    integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                        !< List of elements

    real(rp),    intent(inout)       :: time1(10)                                        ! Timings
    !
    ! Element matrices and vectors (stiffness and preconditioner)
    !
    real(rp)                         :: elauu(VECTOR_DIM,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)                         :: elrbu(VECTOR_DIM,ndime,pnode)                    ! bu
    real(rp)                         :: eldtrho(VECTOR_DIM,pnode)                        ! Projection of rho/dt
    real(rp)                         :: elmurho(VECTOR_DIM,pnode)                        ! Projection of mu/rho
    !
    ! Gather
    !
    real(rp)                         :: elvel(VECTOR_DIM,ndime,pnode)                    ! u
    real(rp)                         :: elcod(VECTOR_DIM,ndime,pnode)                    ! x
    !
    ! Indices and dimensions
    !
    integer                          :: ielem,inode,ivect
    integer(ip)                      :: pelty,j,k
    integer(ip)                      :: ipoin,igaus
    integer(ip)                      :: lnods_loc(VECTOR_DIM,pnode)
    integer(ip)                      :: list_elements_p(VECTOR_DIM)                      ! List of elements (always positive)
    !
    ! Gauss point values
    !
    real(rp)                         :: gpsha(VECTOR_DIM,pnode,pgaus)                    ! N
    real(rp)                         :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)              ! dN/dxi
    real(rp)                         :: gpvol(VECTOR_DIM,pgaus)                          ! w*|J|, |J|
    real(rp)                         :: gpvis(VECTOR_DIM,pgaus)                          ! Viscosity
    real(rp)                         :: gpmut(VECTOR_DIM,pgaus)                          ! mut
    real(rp)                         :: gpden(VECTOR_DIM,pgaus)                          ! Density
    real(rp)                         :: gpadv(VECTOR_DIM,ndime,pgaus)                    ! u+u'
    real(rp)                         :: gprhs(VECTOR_DIM,ndime,pgaus)                    ! RHS
    real(rp)                         :: gpvel(VECTOR_DIM,ndime,pgaus)                    ! u
    real(rp)                         :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)              ! grad(u)
    real(rp)                         :: xjacm(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: xjaci(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: gpdet(VECTOR_DIM,pgaus)
    real(rp)                         :: wgrgr(VECTOR_DIM,pnode,pnode,pgaus)
    real(rp)                         :: agrau(VECTOR_DIM,pnode,pgaus)
    real(rp)                         :: timea,timeb
    integer                          :: idime
    integer(ip)                      :: jdime,jnode,idofv,jdofv
    integer(ip)                      :: ievat,jevat,iauxi,dummi
    !
    ! Internal
    !
#ifdef OPENACCHHH

#define FACT0X     fact0
#define FACT1X     fact1
#define FACT2X     fact2
#define FACT4X     fact4
#define FACT5X     fact5
#define FACT6X     fact6
#define DTINV_LOCX dtinv_loc
#define T1X        t1
#define T2X        t2
#define T3X        t3
#define DENOMX     denom

#else

#define FACT0X     fact0(1:VECTOR_DIM)
#define FACT1X     fact1(1:VECTOR_DIM)
#define FACT2X     fact2(1:VECTOR_DIM)
#define FACT4X     fact4(1:VECTOR_DIM)
#define FACT5X     fact5(1:VECTOR_DIM)
#define FACT6X     fact6(1:VECTOR_DIM)
#define DTINV_LOCX dtinv_loc(1:VECTOR_DIM)
#define T1X        t1(1:VECTOR_DIM)
#define T2X        t2(1:VECTOR_DIM)
#define T3X        t3(1:VECTOR_DIM)
#define DENOMX     denom(1:VECTOR_DIM)

#endif

    real(rp)    :: FACT0X
    real(rp)    :: FACT1X
    real(rp)    :: FACT2X
    real(rp)    :: FACT4X
    real(rp)    :: FACT5X
    real(rp)    :: FACT6X
    real(rp)    :: DTINV_LOCX
    real(rp)    :: T1X
    real(rp)    :: T2X
    real(rp)    :: T3X
    real(rp)    :: DENOMX
    
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#ifdef VECTOR_SIZE_VARIABLE
#define DEF_VECT 1:VECTOR_DIM
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
#endif
    
    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------

    call cputim(timea)
    
    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          list_elements_p(ivect)   = list_elements(ivect)
       else
          list_elements_p(ivect)   = list_elements(1)
       end if
    end do
    do ivect = 1,VECTOR_DIM      
       ielem = list_elements_p(ivect)
       gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
    end do

    !call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! rho
    !call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu
    !call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mut
    
    DTINV_LOCX = dtinv_nsi

    !--------------------------------------------------------------------
    !
    ! Shape function and derivatives
    !
    !--------------------------------------------------------------------
    
    !$acc data create( agrau   , lnods_loc , gpvol     , elauu ,          &
    !$acc              eldtrho , elvel     , gpmut     , gpadv , gpvel ,  &
    !$acc              gpgve   , elcod     , gpsha     , gpcar ,          &
    !$acc              gpden   , wgrgr     , gpdet     ,                  &
    !$acc              xjacm   , xjaci     ,                              &
    !$acc              elmar,                                             &
    !$acc              elmar(pelty)%shape,                                &
    !$acc              elmar(pelty)%deriv,                                &
    !$acc              elmar(pelty)%weigp,                                &
    !$acc              gprhs   , elmurho   , elrbu     , gpvis )          &
    !$acc copyin(      list_elements        ,                             &
    !$acc              gpsha ,                                            &
    !$acc              gpden ,                                            &
    !$acc              gpvis ,                                            &
    !$acc              gpmut ,                                            &
    !$acc              elmar,                                             &
    !$acc              elmar(pelty)%shape ,                               &
    !$acc              elmar(pelty)%deriv ,                               &
    !$acc              elmar(pelty)%weigp                                 )

#ifndef OPENACCHHH
    call cputim(timea)
#endif
    !
    ! This do starts here both in the openacc version and in the nonopenacc
    ! for the non opencc it ends 30  lines later,
    ! in the openacc case it covers all the subroutine.
    ! Similarly the scatter in the nonopenacc case needs a do ivect
    !    
    !$acc parallel loop gang vector default(present)
    !
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
          ielem                    = list_elements(ivect)
       else
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,list_elements(1))
          ielem                    = list_elements(1)
       end if
       !
       ! Transient
       !
       do inode = 1,pnode
          ipoin = lnods_loc(ivect,inode)
          do idime = 1,ndime
             elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
             elcod(ivect,idime,inode) = coord(idime,ipoin)
          end do
       end do
       !
       ! Properties
       !
       gpden(ivect,1:pgaus) = densi_ker % value_const(1)
       gpvis(ivect,1:pgaus) = visco_ker % value_const(1)
       gpmut(ivect,1:pgaus) = turmu_ker % value_ielem(ielem) % a(1:pgaus)       

#ifndef OPENACCHHH
    end do

    call cputim(timeb)
    time1(1) = time1(1) + timeb - timea
    call cputim(timea)
#endif

    !--------------------------------------------------------------------
    !
    ! Element Cartesian derivatives and Jacobian: GPCAR, GPVOL
    !
    !--------------------------------------------------------------------

    if( kfl_savda == 2 ) then
#ifndef OPENACCHHH
       do ivect = 1,VECTOR_DIM      
#endif
          ielem = abs(list_elements(ivect))
          if( ielem > 0 ) then
             gpvol(ivect,1:pgaus)                 = elmda_gpvol(1:pgaus,ielem)
             gpcar(ivect,1:ndime,1:mnode,1:pgaus) = elmda_gpcar(1:ndime,1:mnode,1:pgaus,ielem)
          else
             gpvol(ivect,1:pgaus)                 = 0.0_rp
             gpcar(ivect,1:ndime,1:mnode,1:pgaus) = 0.0_rp
          end if
#ifndef OPENACCHHH
       end do
#endif
    else 

       if(      ndime == 2 ) then
          
          do igaus = 1,pgaus
             xjacm(DEF_VECT,1,1)    =  0.0_rp
             xjacm(DEF_VECT,1,2)    =  0.0_rp
             xjacm(DEF_VECT,2,1)    =  0.0_rp
             xjacm(DEF_VECT,2,2)    =  0.0_rp
             do k = 1,pnode
                xjacm(DEF_VECT,1,1) =  xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,1,2) =  xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,2,1) =  xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,2,2) =  xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(2,k,igaus)
             end do
             gpdet(DEF_VECT,igaus)  =  xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)
             denom                  =  1.0_rp / gpdet(DEF_VECT,igaus)
             xjaci(DEF_VECT,1,1)    =  xjacm(DEF_VECT,2,2) * denom
             xjaci(DEF_VECT,2,2)    =  xjacm(DEF_VECT,1,1) * denom
             xjaci(DEF_VECT,2,1)    = -xjacm(DEF_VECT,2,1) * denom
             xjaci(DEF_VECT,1,2)    = -xjacm(DEF_VECT,1,2) * denom
             do j = 1, pnode
                gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,1) * elmar(pelty) % deriv(2,j,igaus)
                
                gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,2) * elmar(pelty) % deriv(2,j,igaus)
             end do
             gpvol(DEF_VECT,igaus) = elmar(pelty) % weigp(igaus) * gpdet(DEF_VECT,igaus)
          end do
       
       else if( ndime == 3 ) then
          
          do igaus = 1,pgaus
             xjacm(DEF_VECT,1:3,1:3)   = 0.0_rp
             do k = 1,pnode 
                xjacm(DEF_VECT,1,1)    =  xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,1,2)    =  xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,1,3)    =  xjacm(DEF_VECT,1,3) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(3,k,igaus)
                xjacm(DEF_VECT,2,1)    =  xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,2,2)    =  xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,2,3)    =  xjacm(DEF_VECT,2,3) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(3,k,igaus)
                xjacm(DEF_VECT,3,1)    =  xjacm(DEF_VECT,3,1) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,3,2)    =  xjacm(DEF_VECT,3,2) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,3,3)    =  xjacm(DEF_VECT,3,3) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(3,k,igaus)
             end do
             t1                        =  xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,2,3)
             t2                        = -xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,3)
             t3                        =  xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,2) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,2)
             gpdet(DEF_VECT,igaus)     =  xjacm(DEF_VECT,1,1) * t1 + xjacm(DEF_VECT,1,2) * t2 + xjacm(DEF_VECT,1,3) * t3
             denom                     =  1.0_rp / gpdet(DEF_VECT,igaus)
             xjaci(DEF_VECT,1,1)       =  t1 * denom
             xjaci(DEF_VECT,2,1)       =  t2 * denom
             xjaci(DEF_VECT,3,1)       =  t3 * denom
             xjaci(DEF_VECT,2,2)       =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,1,3)) * denom
             xjaci(DEF_VECT,3,2)       =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,2) + xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,1)) * denom
             xjaci(DEF_VECT,3,3)       =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)) * denom
             xjaci(DEF_VECT,1,2)       =  (-xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,1,3)) * denom
             xjaci(DEF_VECT,1,3)       =  ( xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,2,3) - xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,1,3)) * denom
             xjaci(DEF_VECT,2,3)       =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,3) + xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,3)) * denom
             
             do j = 1, pnode
                gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,1) * elmar(pelty) % deriv(2,j,igaus) &
                     &                      + xjaci(DEF_VECT,3,1) * elmar(pelty) % deriv(3,j,igaus)
                
                gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,2) * elmar(pelty) % deriv(2,j,igaus) &
                     &                      + xjaci(DEF_VECT,3,2) * elmar(pelty) % deriv(3,j,igaus)
                
                gpcar(DEF_VECT,3,j,igaus) =   xjaci(DEF_VECT,1,3) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,3) * elmar(pelty) % deriv(2,j,igaus) &
                     &                      + xjaci(DEF_VECT,3,3) * elmar(pelty) % deriv(3,j,igaus)
             end do
             gpvol(DEF_VECT,igaus) = elmar(pelty) % weigp(igaus) * gpdet(DEF_VECT,igaus)
          end do
       end if


    end if

    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)
    gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpmut(DEF_VECT,:)  ! Effective viscosity <= mu+mut
    
    !----------------------------------------------------------------------
    !
    ! Gauss point values & eldtrho, elmurho
    !
    !----------------------------------------------------------------------

    gpgve(DEF_VECT,:,:,:) = 0.0_rp
    gprhs(DEF_VECT,:,:)   = 0.0_rp
    gpvel(DEF_VECT,:,:)   = 0.0_rp
    eldtrho(DEF_VECT,:)   = 0.0_rp
    elmurho(DEF_VECT,:)   = 0.0_rp

    do igaus = 1,pgaus
       FACT2X =  gpden(DEF_VECT,igaus) * grnor_nsi
       do idime = 1,ndime
          do inode = 1,pnode
             gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) + elvel(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)
             do jdime = 1,ndime
                gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                     + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode)
             end do
          end do
          gpadv(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus)
          gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime)
       end do
       FACT2X = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
       FACT4X = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus)  )
       do inode = 1,pnode
          eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT2X
          elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT4X
       end do
    end do

    !----------------------------------------------------------------------
    !
    ! Element matrices 
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    agrau(DEF_VECT,:,:)   = 0.0_rp
    wgrgr(DEF_VECT,:,:,:) = 0.0_rp
    !
    ! AGRAU = rho * (a.grad) Ni
    ! WGRGR = grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             agrau(DEF_VECT,inode,igaus) =  agrau(DEF_VECT,inode,igaus) + &
                  &                         gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
          end do
          agrau(DEF_VECT,inode,igaus) =  gpden(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus)
          do jnode = 1,pnode
             do idime = 1,ndime
                wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + &
                     &                              gpcar(DEF_VECT,idime,inode,igaus)*gpcar(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do
    end do
    !
    !  rho*(a.grad)Nj ) Ni + mu * grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus

       FACT6X = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode
          do idime = 1,ndime

             idofv           = (inode-1)*ndime+idime

             do jnode = 1,pnode
                jdofv           = (jnode-1)*ndime+idime
                FACT4X = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                FACT5X = FACT4X * ( agrau(DEF_VECT,jnode,igaus) ) &
                     &           + FACT6X *   wgrgr(DEF_VECT,inode,jnode,igaus)
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT5X
             end do
          end do
       end do
    end do

    do igaus = 1,pgaus

       FACT0X = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       FACT2X = 0.0_rp
       do idime = 1,ndime
          FACT2X = FACT2X + gpgve(DEF_VECT,idime,idime,igaus)
       end do
       FACT2X = FACT2X * FACT0X

       do inode = 1,pnode
          do idime = 1,ndime
             idofv = (inode-1)*ndime+idime
             do jnode = 1,pnode
                jdofv = (jnode-1)*ndime+idime
                !
                ! rho * (div u) u
                !
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                     + FACT2X * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                !
                ! rho * u.grad(u)^t
                ! 
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                        + FACT0X * gpsha(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) * gpvel(DEF_VECT,jdime,igaus)
                end do

             end do
             ! 
             ! 0.5 * grad(u**2)
             !
             do jdime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     +  FACT0X * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus)*gpgve(DEF_VECT,idime,jdime,igaus)
             end do
          end do
          !
          ! ( mu*duj/dxi , dv/dxj ) (only div form)
          !      
          do idime = 1,ndime
             idofv = (inode-1)*ndime + idime
             do jnode = 1,pnode
                FACT1X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                do jdime = 1,ndime
                   jdofv                       = (jnode-1)*ndime + jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT1X * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
          !
          ! bu = ( f , v ) 
          !
          FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus) 
          end do
       end do
    end do
    !
    ! Send matrix to RHS
    !
    do jnode = 1,pnode
       do jdime = 1,ndime
          jevat = (jnode-1)*ndime+jdime
          do inode = 1,pnode
             do idime = 1,ndime
                ievat = (inode-1)*ndime+idime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     - elauu(DEF_VECT,ievat,jevat) * elvel(DEF_VECT,jdime,jnode) 
             end do
          end do
       end do
    end do
    !
    ! Scatter element matrix to global one 
    ! 
#ifndef OPENACCHHH
    call cputim(timeb)
    time1(7) = time1(7) + timeb - timea
    call cputim(timea)
    do ivect = 1,VECTOR_DIM                        
#endif 
       ielem = list_elements(ivect)
       if ( ielem > 0 ) then
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
             do idime = 1,ndime
                iauxi = idime + (ipoin-1) * ndime 
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                !$acc atomic update
                rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect,idime,inode)

             end do
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             !$acc atomic update
             dt_rho_nsi(ipoin) = dt_rho_nsi(ipoin) + eldtrho(ivect,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             !$acc atomic update
             mu_rho_nsi(ipoin) = mu_rho_nsi(ipoin) + elmurho(ivect,inode)

          end do
       end if
    end do
    !$acc end parallel loop
    !$acc end data


#ifndef OPENACCHHH
    call cputim(timeb)
    time1(8) = time1(8) + timeb - timea                       
#endif

  end subroutine nsi_element_operations_fast6


    subroutine nsi_element_operations_fast8(&
       VECTOR_DIM,pnode,pgaus,list_elements,time1,streamid)

    use def_kintyp,            only : ip,rp                                              ! in_const_scalar
    use def_master,            only : rhsid                                              ! out       ! real(rp), pointer     :: rhsid(:)
    use def_master,            only : veloc                                              ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
    use def_domain,            only : coord                                              ! in_const  ! real(rp), pointer     :: coord(:,:)
    use def_domain,            only : ltype                                              ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
    use def_domain,            only : lnods                                              ! in_const  ! integer(ip), pointer  :: lnods(:,:)
    use def_domain,            only : mnode                                              ! in_const scalar
    use def_domain,            only : elmar
    use def_domain,            only : elmda
    use def_domain,            only : ntens                                              ! in_const scalar
    use def_domain,            only : kfl_savda
    use def_domain,            only : elmda_gpvol
    use def_domain,            only : elmda_gpcar
    use def_nastin,            only : dtinv_nsi                                          ! in_var   scalar
    use def_nastin,            only : dt_rho_nsi                                         ! out  ! real(rp), pointer     :: dt_rho_nsi(:)
    use def_nastin,            only : mu_rho_nsi                                         ! out  ! real(rp), pointer     :: mu_rho_nsi(:)
    use def_nastin,            only : grnor_nsi                                          ! in_const_scalar
    use def_nastin,            only : gravi_nsi                                          ! in_const  ! real(rp)         ::gravi_nsi(3)
    use mod_ker_proper,        only : ker_proper
    use def_kermod,            only : densi_ker 
    use def_kermod,            only : visco_ker 
    use def_kermod,            only : turmu_ker 
    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in)          :: VECTOR_DIM                                       !< Number of nodes
    integer(ip), intent(in)          :: pnode                                            !< Number of nodes
    integer(ip), intent(in)          :: pgaus                                            !< Number of Gauss points
    integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                        !< List of elements

    real(rp),    intent(inout)       :: time1(10)                                        ! Timings
    integer(ip)                      :: streamid
    !
    ! Element matrices and vectors (stiffness and preconditioner)
    !
    real(rp)                         :: elauu(VECTOR_DIM,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)                         :: elrbu(VECTOR_DIM,ndime,pnode)                    ! bu
    real(rp)                         :: eldtrho(VECTOR_DIM,pnode)                        ! Projection of rho/dt
    real(rp)                         :: elmurho(VECTOR_DIM,pnode)                        ! Projection of mu/rho
    !
    ! Gather
    !
    real(rp)                         :: elvel(VECTOR_DIM,ndime,pnode)                    ! u
    real(rp)                         :: elcod(VECTOR_DIM,ndime,pnode)                    ! x
    !
    ! Indices and dimensions
    !
    integer                          :: ielem,inode,ivect
    integer(ip)                      :: pelty,j,k
    integer(ip)                      :: ipoin,igaus
    integer(ip)                      :: lnods_loc(VECTOR_DIM,pnode)
    integer(ip)                      :: list_elements_p(VECTOR_DIM)                      ! List of elements (always positive)
    !
    ! Gauss point values
    !
    real(rp)                         :: gpsha(VECTOR_DIM,pnode,pgaus)                    ! N
    real(rp)                         :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)              ! dN/dxi
    real(rp)                         :: gpvol(VECTOR_DIM,pgaus)                          ! w*|J|, |J|
    real(rp)                         :: gpvis(VECTOR_DIM,pgaus)                          ! Viscosity
    real(rp)                         :: gpmut(VECTOR_DIM,pgaus)                          ! mut
    real(rp)                         :: gpden(VECTOR_DIM,pgaus)                          ! Density
    real(rp)                         :: gpadv(VECTOR_DIM,ndime,pgaus)                    ! u+u'
    real(rp)                         :: gprhs(VECTOR_DIM,ndime,pgaus)                    ! RHS
    real(rp)                         :: gpvel(VECTOR_DIM,ndime,pgaus)                    ! u
    real(rp)                         :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)              ! grad(u)
    real(rp)                         :: xjacm(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: xjaci(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: gpdet(VECTOR_DIM,pgaus)
    real(rp)                         :: wgrgr(VECTOR_DIM,pnode,pnode,pgaus)
    real(rp)                         :: agrau(VECTOR_DIM,pnode,pgaus)
    real(rp)                         :: timea,timeb
    integer                          :: idime
    integer(ip)                      :: jdime,jnode,idofv,jdofv
    integer(ip)                      :: ievat,jevat,iauxi,dummi
    real(rp)                         :: scaden,scavis   !
    ! Internal
    !
#ifdef OPENACCHHH

#define FACT0X     fact0
#define FACT1X     fact1
#define FACT2X     fact2
#define FACT4X     fact4
#define FACT5X     fact5
#define FACT6X     fact6
#define DTINV_LOCX dtinv_loc
#define T1X        t1
#define T2X        t2
#define T3X        t3
#define DENOMX     denom

#else

#define FACT0X     fact0(1:VECTOR_DIM)
#define FACT1X     fact1(1:VECTOR_DIM)
#define FACT2X     fact2(1:VECTOR_DIM)
#define FACT4X     fact4(1:VECTOR_DIM)
#define FACT5X     fact5(1:VECTOR_DIM)
#define FACT6X     fact6(1:VECTOR_DIM)
#define DTINV_LOCX dtinv_loc(1:VECTOR_DIM)
#define T1X        t1(1:VECTOR_DIM)
#define T2X        t2(1:VECTOR_DIM)
#define T3X        t3(1:VECTOR_DIM)
#define DENOMX     denom(1:VECTOR_DIM)

#endif

    real(rp)    :: FACT0X
    real(rp)    :: FACT1X
    real(rp)    :: FACT2X
    real(rp)    :: FACT4X
    real(rp)    :: FACT5X
    real(rp)    :: FACT6X
    real(rp)    :: DTINV_LOCX
    real(rp)    :: T1X
    real(rp)    :: T2X
    real(rp)    :: T3X
    real(rp)    :: DENOMX
    
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#ifdef VECTOR_SIZE_VARIABLE
#define DEF_VECT 1:VECTOR_DIM
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
#endif
    
    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------

    call cputim(timea)
    scaden=densi_ker % value_const(1)  
    scavis=visco_ker % value_const(1)
 
    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          list_elements_p(ivect)   = list_elements(ivect)
       else
          list_elements_p(ivect)   = list_elements(1)
       end if
    end do
    do ivect = 1,VECTOR_DIM      
       ielem = list_elements_p(ivect)
       gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
    end do

 !   call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! rho
!    call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mut
    
    DTINV_LOCX = dtinv_nsi
#ifdef OPENACCHHH

    !--------------------------------------------------------------------
    !
    ! Shape function and derivatives
    !
    !--------------------------------------------------------------------
    !$acc enter data create( agrau, lnods_loc , gpvol     , elauu ,          &
    !$acc              eldtrho , elvel     ,  gpadv , gpvel ,  &
    !$acc              gpgve   , elcod     ,  gpcar ,          &
    !$acc              wgrgr     , gpdet     ,gpden, gpvis,           &
    !$acc              xjacm   , xjaci     ,                              &
    !$acc              gprhs   , elmurho   , elrbu)          &
    !$acc copyin(      list_elements        ,coord,                             &
    !$acc              gpsha ,                                            &
    !$acc              gpmut ,                                            &
    !$acc              elmar,                                             &
    !$acc              elmar(pelty)%shape ,                               &
    !$acc              elmar(pelty)%deriv ,                               &
    !$acc              elmar(pelty)%weigp                 ) async(streamid)
#endif   

#ifndef OPENACCHHH
    call cputim(timea)
#endif
    !
    ! This do starts here both in the openacc version and in the nonopenacc
    ! for the non opencc it ends 30  lines later,
    ! in the openacc case it covers all the subroutine.
    ! Similarly the scatter in the nonopenacc case needs a do ivect
    !    
#ifdef OPENACCHHH
    !$acc parallel loop gang vector default(present) &
    !$acc        async(streamid)
    !
#endif
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
          ielem                    = list_elements(ivect)
       else
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,list_elements(1))
          ielem                    = list_elements(1)
       end if
       !
       ! Transient
       !
       do inode = 1,pnode
          ipoin = lnods_loc(ivect,inode)
          do idime = 1,ndime
             elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
             elcod(ivect,idime,inode) = coord(idime,ipoin)
          end do
       end do
       !
       gpden(ivect,1:pgaus) = scaden 
       gpvis(ivect,1:pgaus) = scavis 
!
#ifndef OPENACCHHH
    end do

    call cputim(timeb)
    time1(1) = time1(1) + timeb - timea
    call cputim(timea)
#endif

    !--------------------------------------------------------------------
    !
    ! Element Cartesian derivatives and Jacobian: GPCAR, GPVOL
    !
    !--------------------------------------------------------------------

    if( kfl_savda == 2 ) then
#ifndef OPENACCHHH
       do ivect = 1,VECTOR_DIM      
#endif
          ielem = abs(list_elements(ivect))
          if( ielem > 0 ) then
             gpvol(ivect,1:pgaus)                 = elmda_gpvol(1:pgaus,ielem)
             gpcar(ivect,1:ndime,1:mnode,1:pgaus) = elmda_gpcar(1:ndime,1:mnode,1:pgaus,ielem)
          else
             gpvol(ivect,1:pgaus)                 = 0.0_rp
             gpcar(ivect,1:ndime,1:mnode,1:pgaus) = 0.0_rp
          end if
#ifndef OPENACCHHH
       end do
#endif
    else 

       if(      ndime == 2 ) then
          
          do igaus = 1,pgaus
             xjacm(DEF_VECT,1,1)    =  0.0_rp
             xjacm(DEF_VECT,1,2)    =  0.0_rp
             xjacm(DEF_VECT,2,1)    =  0.0_rp
             xjacm(DEF_VECT,2,2)    =  0.0_rp
             do k = 1,pnode
                xjacm(DEF_VECT,1,1) =  xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,1,2) =  xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,2,1) =  xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,2,2) =  xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(2,k,igaus)
             end do
             gpdet(DEF_VECT,igaus)  =  xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)
             denom                  =  1.0_rp / gpdet(DEF_VECT,igaus)
             xjaci(DEF_VECT,1,1)    =  xjacm(DEF_VECT,2,2) * denom
             xjaci(DEF_VECT,2,2)    =  xjacm(DEF_VECT,1,1) * denom
             xjaci(DEF_VECT,2,1)    = -xjacm(DEF_VECT,2,1) * denom
             xjaci(DEF_VECT,1,2)    = -xjacm(DEF_VECT,1,2) * denom
             do j = 1, pnode
                gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,1) * elmar(pelty) % deriv(2,j,igaus)
                
                gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,2) * elmar(pelty) % deriv(2,j,igaus)
             end do
             gpvol(DEF_VECT,igaus) = elmar(pelty) % weigp(igaus) * gpdet(DEF_VECT,igaus)
          end do
       
       else if( ndime == 3 ) then
          
          do igaus = 1,pgaus
             xjacm(DEF_VECT,1:3,1:3)   = 0.0_rp
             do k = 1,pnode 
                xjacm(DEF_VECT,1,1)    =  xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,1,2)    =  xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,1,3)    =  xjacm(DEF_VECT,1,3) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(3,k,igaus)
                xjacm(DEF_VECT,2,1)    =  xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,2,2)    =  xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,2,3)    =  xjacm(DEF_VECT,2,3) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(3,k,igaus)
                xjacm(DEF_VECT,3,1)    =  xjacm(DEF_VECT,3,1) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(1,k,igaus)
                xjacm(DEF_VECT,3,2)    =  xjacm(DEF_VECT,3,2) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(2,k,igaus)
                xjacm(DEF_VECT,3,3)    =  xjacm(DEF_VECT,3,3) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(3,k,igaus)
             end do
             t1                        =  xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,2,3)
             t2                        = -xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,3)
             t3                        =  xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,2) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,2)
             gpdet(DEF_VECT,igaus)     =  xjacm(DEF_VECT,1,1) * t1 + xjacm(DEF_VECT,1,2) * t2 + xjacm(DEF_VECT,1,3) * t3
             denom                     =  1.0_rp / gpdet(DEF_VECT,igaus)
             xjaci(DEF_VECT,1,1)       =  t1 * denom
             xjaci(DEF_VECT,2,1)       =  t2 * denom
             xjaci(DEF_VECT,3,1)       =  t3 * denom
             xjaci(DEF_VECT,2,2)       =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,1,3)) * denom
             xjaci(DEF_VECT,3,2)       =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,2) + xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,1)) * denom
             xjaci(DEF_VECT,3,3)       =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)) * denom
             xjaci(DEF_VECT,1,2)       =  (-xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,1,3)) * denom
             xjaci(DEF_VECT,1,3)       =  ( xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,2,3) - xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,1,3)) * denom
             xjaci(DEF_VECT,2,3)       =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,3) + xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,3)) * denom
             
             do j = 1, pnode
                gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,1) * elmar(pelty) % deriv(2,j,igaus) &
                     &                      + xjaci(DEF_VECT,3,1) * elmar(pelty) % deriv(3,j,igaus)
                
                gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,2) * elmar(pelty) % deriv(2,j,igaus) &
                     &                      + xjaci(DEF_VECT,3,2) * elmar(pelty) % deriv(3,j,igaus)
                
                gpcar(DEF_VECT,3,j,igaus) =   xjaci(DEF_VECT,1,3) * elmar(pelty) % deriv(1,j,igaus) &
                     &                      + xjaci(DEF_VECT,2,3) * elmar(pelty) % deriv(2,j,igaus) &
                     &                      + xjaci(DEF_VECT,3,3) * elmar(pelty) % deriv(3,j,igaus)
             end do
             gpvol(DEF_VECT,igaus) = elmar(pelty) % weigp(igaus) * gpdet(DEF_VECT,igaus)
          end do
       end if


    end if

!!! ----
    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)
    gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpmut(DEF_VECT,:)  ! Effective viscosity <= mu+mut
    
    !----------------------------------------------------------------------
    !
    ! Gauss point values & eldtrho, elmurho
    !
    !----------------------------------------------------------------------

    gpgve(DEF_VECT,:,:,:) = 0.0_rp
    gprhs(DEF_VECT,:,:)   = 0.0_rp
    gpvel(DEF_VECT,:,:)   = 0.0_rp
    eldtrho(DEF_VECT,:)   = 0.0_rp
    elmurho(DEF_VECT,:)   = 0.0_rp
    do igaus = 1,pgaus
       FACT2X =  gpden(DEF_VECT,igaus) * grnor_nsi
       do idime = 1,ndime
          do inode = 1,pnode
             gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) + elvel(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)
             do jdime = 1,ndime
                gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                     + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode)
             end do
          end do
          gpadv(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus)
          gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime)
       end do
       FACT2X = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
       FACT4X = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus)  )
       do inode = 1,pnode
          eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT2X
          elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT4X
       end do
    end do



    !----------------------------------------------------------------------
    !
    ! Element matrices 
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    agrau(DEF_VECT,:,:)   = 0.0_rp
    wgrgr(DEF_VECT,:,:,:) = 0.0_rp
    !
    ! AGRAU = rho * (a.grad) Ni
    ! WGRGR = grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             agrau(DEF_VECT,inode,igaus) =  agrau(DEF_VECT,inode,igaus) + &
                  &                         gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
          end do
          agrau(DEF_VECT,inode,igaus) =  gpden(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus)
          do jnode = 1,pnode
             do idime = 1,ndime
                wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + &
                     &                              gpcar(DEF_VECT,idime,inode,igaus)*gpcar(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do
    end do

    !
    !  rho*(a.grad)Nj ) Ni + mu * grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus

       FACT6X = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode
          do idime = 1,ndime

             idofv           = (inode-1)*ndime+idime

             do jnode = 1,pnode
                jdofv           = (jnode-1)*ndime+idime
                FACT4X = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                FACT5X = FACT4X * ( agrau(DEF_VECT,jnode,igaus) ) &
                     &           + FACT6X *   wgrgr(DEF_VECT,inode,jnode,igaus)
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT5X
             end do
          end do
       end do
    end do

    do igaus = 1,pgaus

       FACT0X = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       FACT2X = 0.0_rp
       do idime = 1,ndime
          FACT2X = FACT2X + gpgve(DEF_VECT,idime,idime,igaus)
       end do
       FACT2X = FACT2X * FACT0X

       do inode = 1,pnode
          do idime = 1,ndime
             idofv = (inode-1)*ndime+idime
             do jnode = 1,pnode
                jdofv = (jnode-1)*ndime+idime
                !
                ! rho * (div u) u
                !
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                     + FACT2X * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                !
                ! rho * u.grad(u)^t
                ! 
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                        + FACT0X * gpsha(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) * gpvel(DEF_VECT,jdime,igaus)
                end do

             end do
             ! 
             ! 0.5 * grad(u**2)
             !
             do jdime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     +  FACT0X * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus)*gpgve(DEF_VECT,idime,jdime,igaus)
             end do
          end do
          !
          ! ( mu*duj/dxi , dv/dxj ) (only div form)
          !      
          do idime = 1,ndime
             idofv = (inode-1)*ndime + idime
             do jnode = 1,pnode
                FACT1X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                do jdime = 1,ndime
                   jdofv                       = (jnode-1)*ndime + jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT1X * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
          !
          ! bu = ( f , v ) 
          !
          FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus) 
          end do
       end do
    end do

!!!    !ultima version


    !
    ! Send matrix to RHS
    !
    do jnode = 1,pnode
       do jdime = 1,ndime
          jevat = (jnode-1)*ndime+jdime
          do inode = 1,pnode
             do idime = 1,ndime
                ievat = (inode-1)*ndime+idime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     - elauu(DEF_VECT,ievat,jevat) * elvel(DEF_VECT,jdime,jnode) 
             end do
          end do
       end do
    end do
    !
    ! Scatter element matrix to global one 
    ! 
#ifndef OPENACCHHH
    call cputim(timeb)
    time1(7) = time1(7) + timeb - timea
    call cputim(timea)
    do ivect = 1,VECTOR_DIM                        
#endif 
       ielem = list_elements(ivect)
       if ( ielem > 0 ) then
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
             do idime = 1,ndime
               iauxi = idime + (ipoin-1) * ndime 

#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
 
#ifdef OPENACCHHH
               !$acc atomic update
#endif
                rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect,idime,inode)

             end do
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
#ifdef OPENACCHHH
            !$acc atomic update
#endif
             dt_rho_nsi(ipoin) = dt_rho_nsi(ipoin) + eldtrho(ivect,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
#ifdef OPENACCHHH
           !$acc atomic update
#endif
             mu_rho_nsi(ipoin) = mu_rho_nsi(ipoin) + elmurho(ivect,inode)

          end do
       end if
    end do
#ifdef OPENACCHHH

    !$acc end parallel loop

    !$acc exit data delete (agrau, lnods_loc , gpvol, elauu, &
    !$acc      eldtrho , elvel     , gpmut     , gpadv , gpvel ,  &
    !$acc      gpgve   , elcod     , gpsha     , gpcar ,          &
    !$acc      gpden   , wgrgr     , gpdet     ,                  &
    !$acc      xjacm   , xjaci     ,                              &
    !$acc      elmar,                                             &
    !$acc      elmar(pelty)%shape,                                &
    !$acc      elmar(pelty)%deriv,                                &
    !$acc      elmar(pelty)%weigp,                                &
    !$acc      gprhs   , elmurho   , elrbu, gpvis,list_elements ) &
    !$acc    async(streamid)

#endif


#ifndef OPENACCHHH
    call cputim(timeb)
    time1(8) = time1(8) + timeb - timea                       
#endif

  end subroutine nsi_element_operations_fast8




end module mod_nsi_element_operations_fast
