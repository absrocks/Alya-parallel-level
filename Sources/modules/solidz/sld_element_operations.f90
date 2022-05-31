!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_element_operations.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Elemental operations
!> @details Elemental operations
!>          1. Compute elemental matrix and RHS
!>          2. Impose Dirichlet boundary conditions
!>          3. Assemble them
!> @}
!-----------------------------------------------------------------------

subroutine sld_element_operations(itask,ielem,ielem_min,gpdet_min)

  use def_master                                                               ! general global variables
  use def_elmtyp                                                               ! type of elements
  use def_domain                                                               ! geometry information
  use def_solidz                                                               ! general solidz module information
  use def_kermod,                 only : kfl_element_to_csr                    ! element to CSR database
  use mod_ker_detection,          only : ker_detection_min_max_value           ! Event detection output
  use mod_communications,         only : PAR_MAX
  use mod_matrix,                 only : matrix_assemble_element_matrix_to_CSR ! Matrix assembly
  use mod_matrix,                 only : matrix_assemble_element_RHS           ! RHS assembly
  use mod_sld_interface_element,  only : ELINT_elemental_operations            ! Cohesive element
  use mod_sld_contshell_elements, only : SHELL_elemental_operations
  use mod_sld_fibers
  use mod_exm_sld_eccoupling
  use mod_sld_cardiac_cycle

  implicit none

  integer(ip) :: itask                            !< What to do
  integer(ip) :: ielem                            !< Element to assemble
  integer(ip) :: ielem_min                        !< Element with minimum negative Jacobian
  real(rp)    :: gpdet_min                        !< Mminimum negative Jacobian
  integer(ip) :: pevat                            !< Indices and dimensions
  integer(ip) :: dummi
  integer(ip) :: pelty,pmate,pnode,pgaus,plapl,porde
  real(rp)    :: elrhs(ndofn_sld*mnode)
  real(rp)    :: elmat(ndofn_sld*mnode,ndofn_sld*mnode)
  real(rp)    :: elfin(ndofn_sld*mnode)
  real(rp)    :: elfex(ndofn_sld*mnode)
  real(rp)    :: elmac(ndofn_sld*mnode)
  real(rp)    :: elfda(ndofn_sld*mnode)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: eldis(ndime,mnode,ncomp_sld)
  real(rp)    :: eldip(ndime,mnode)
  real(rp)    :: elddis(ndime,mnode,ncomp_sld)
  real(rp)    :: eldix(ndime,mnode,ncomp_sld)
  real(rp)    :: elvel(ndime,mnode,ncomp_sld)
  real(rp)    :: elacc(ndime,mnode,3)
  real(rp)    :: elrst(    6,mnode)
  real(rp)    :: elvex(ndime,mnode,ncomp_sld)
  real(rp)    :: elepo(ndime,ndime,mnode)
  real(rp)    :: elepo_new_inverse(ndime,ndime,mnode)
  real(rp)    :: elepo_new_determinant(mnode)
  real(rp)    :: elfrx(ndofn_sld,mnode)
  real(rp)    :: elfco(ndofn_sld,mnode)
  real(rp)    :: gpvol(mgaus),gpdet(mgaus)
  real(rp)    :: gpdis(ndime,mgaus,ncomp_sld)           ! Displacement at time n
  real(rp)    :: gpvel(ndime,mgaus,ncomp_sld)           ! Velocity at time n
  real(rp)    :: gpacc(ndime,mgaus,3)                   ! Acceleration
  real(rp)    :: gprat(ndime,ndime,mgaus)               ! Rate of deformation Fdot at time n
  real(rp)    :: gpgdi(ndime,ndime,mgaus)               ! Displacement Gradient F
  real(rp)    :: gpigd(ndime,ndime,mgaus)               ! Displacement Gradient Inverse F^{-1}
  real(rp)    :: gpgdi_eps(ndime,ndime,mgaus)           ! Displacement Gradient F for the perturbed state
  real(rp)    :: gpigd_eps(ndime,ndime,mgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp)    :: gpdet_eps(mgaus)
  real(rp)    :: gppio_eps(ndime,ndime,mgaus)
  real(rp)    :: gpcau(ndime,ndime,mgaus)               ! Cauchy tensor C
  real(rp)    :: gpstr(ndime,ndime,mgaus,2)             ! Stress tensor S
  real(rp)    :: gprestre(ndime,ndime,mgaus)            ! Residual stress tensor
  real(rp)    :: gpfibe(ndime,mgaus)                    ! Fibers length vector (relative to its ini. length)
  real(rp)    :: gpene(mgaus)                           ! Energy
  real(rp)    :: gpsha(mnode,mgaus)                     ! Ni
  real(rp)    :: gpmof(mgaus)                           ! Modulator fields
  real(rp)    :: gphes(ntens,mnode,mgaus)               ! dNk/dxidxj
  real(rp)    :: gpder(ndime,mnode,mgaus)               ! dNk/dsj
  real(rp)    :: gpcar(ndime,mnode,mgaus)               ! dNk/dxj
  real(rp)    :: gpcod(ndime,mgaus)                     ! Gauss points' physical coordinates
  real(rp)    :: gphea(mnode,mgaus)                     ! Shifted Heaviside (XFEM)
  real(rp)    :: gptmo(ndime,ndime,ndime,ndime,mgaus)   ! Tangent moduli
  real(rp)    :: elmof(mnode)                           ! Modulating fields
  real(rp)    :: elmuu(mnode)
  real(rp)    :: elmaa(mnode)
  real(rp)    :: elmas(mnode,mnode)
  real(rp)    :: gpdds(ndime,ndime,ndime,ndime,mgaus)
  logical(lg) :: debugging
  real(rp)    :: hleng(3)
  real(rp)    :: tragl(9)
  real(rp)    :: gpcal(mgaus), gptrp(mgaus)                         ! Calcium   {SLD-EXM-ECC Coupling}
  real(rp)    :: gpfbr0(3,mgaus), gpsht0(3,mgaus), gpnrm0(3,mgaus)  ! Fiber directions  {SLD-EXM-ECC Coupling}

  !----------------------------------------------------------------------
  !
  ! RHS assembly (r) and Jacobian matrix (J)
  !
  !----------------------------------------------------------------------

  debugging = .false.
  !
  ! Element dimensions
  !
  pelty = ltype(ielem)
  !
  if ( pelty > 0 .and. pelty /= SHELL .and. pelty /= BAR3D .and. lawst_sld(lmate(ielem)) /= 200_ip .and. lelch(ielem) /= ELINT) then
     pnode = nnode(pelty)
     porde = lorde(pelty)
     pgaus = ngaus(pelty)
     plapl = llapl(pelty)
     pevat = ndofn_sld * pnode
     pmate = 1
     gpmof = 1.0_rp

     if( nmate > 1 ) then
        pmate = lmate_sld(ielem)
     end if

     !
     ! Gather operations
     !
     call sld_elmgat(&
          pnode,lnods(1,ielem),eldis,elddis,eldip,elvel,elacc,elcod,eldix,elvex,elmof,elrst)

     !
     ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, GPVOL
     !
     call elmca2(&
          pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
          elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
          gpder,gpcar,gphes,ielem)

     !
     ! Element characteristic lenght
     !
     call elmlen(&
          ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
          hnatu(pelty),hleng)
     celen_sld(ielem) = hleng(ndime)/real(porde,rp)

     !
     ! Variables at Gauss points
     !
     if (kfl_xfeme_sld == 1) then
        call gather(2_ip,pgaus,pnode,ndime,dummi,gpsha,elcod,gpcod)
        call sld_xfemic(2_ip,pnode,pgaus,gphea,gpcod,elcod,ielem,dummi)
     end if

     call sld_elmpre(&
          ielem,pnode,pgaus,pmate,eldis,eldip,elvel,elacc,eldix,elvex,elcod,&
          elmof,elrst,gpsha,gpcar,gphea,gpgdi_sld(ielem)%a,gpdet,gpvol,&
          elepo,elepo_new_inverse,elepo_new_determinant,&
          gpvel,gpacc,&
          gprat,gpdis,gpcau,gpigd,gpcod,gpgdi,&
          gpigd_eps,gpgdi_eps,gpdet_eps,&
          gpmof,gprestre,gpdet_min,ielem_min)

     !
     ! Constitutive law
     !
     call sld_elmcla(&
          itask,pgaus,pmate,gpvol,gpcau,gpgdi_sld(ielem)%a,gpgdi,gpene,gpstr,gppio_sld(ielem)%a,&
          gpdet,gprat,gpigd,gptmo,&
          gpfibe,ielem,elcod,pnode,lnods(1,ielem),gpsha,gpdds,&
          gpigd_eps,gpgdi_eps,gpdet_eps,gppio_eps,&
          gpmof)

     !
     ! Include pre-stress to active force
     !
     if(cycles(1)% pstr0 .gt.0.0_rp) then
       call runend("SLD_ELEMENT_OPERATIONS: SOLIDZ NOT PREPARED TO DEAL WITH PRE-STRESS")
       ! The pre stress should be added to the active contribution.
       ! Included pre-stresss should be:
       ! tactf=tactf+prest_sld(imate)
     endif

     !
     ! Adding Active Contribution (SLD-EXM-ECC Coupling)
     !
     if( has_exmsld_coupling() .or. kfl_exmsld_3Dcou_ecc) then
          call sld_fib_get_normalized_fiber_directions_at_gauss_points( & 
          &    ielem, pnode, pgaus, lnods, gpsha, gpfbr0, gpsht0, gpnrm0 ) 
          call exm_sld_ecc_set_fibers_at_gp( &
          &    ielem, gpfibe )
          call exm_sld_ecc_get_calcium_at_gp( &
          &    ielem, pnode, pgaus, lnods(1:pnode,ielem), gpsha, gpcal )
          call exm_sld_ecc_get_troponin_at_gp( &
          &    ielem, pgaus, gptrp )
          call exm_sld_ecc_add_active_stress_and_moduli( &
          &    itask, ielem, pmate, pgaus, &
          &    kfl_cellmod(pmate), kfl_eccty(pmate), kfl_copt_ecc(pmate), &
          &    gpdet, gpgdi, gpcau, gpcal, &
          &    gpfbr0, gpsht0, gpnrm0, &
          &    gpstr(:,:,:,1), gppio_sld(ielem)%a, gpdds, gptmo )
          if( itask <= 2 )then
               call exm_sld_ecc_assemble_troponin( &
               &    ielem, pnode, pgaus, lnods(1:pnode,ielem), gpsha, gpvol, gptrp )
          endif
     end if
     
     !
     ! This if statement is moved here, because elmpre, elmcla, etc are also
     ! required for post-processing of stress, strain and internal variables
     ! (the call to sld_bitcul with itask = 3)
     !
     if( itask <= 2 ) then
        !
        ! RHS and MATRIX
        !
        call sld_elmmat(itask,&
             ielem,pgaus,pmate,pnode,lnods(1,ielem),gpgdi,gppio_sld(ielem)%a,gppio_eps,gpstr,&
             gpvol,gpvel,gpacc,gptmo,gprat,gpsha,gpcar,gphea,gpdis,gpdds,gprestre,eldis,elcod,&
             elrhs,elmat,elmuu,elmaa,elfrx,elfin,elfex,elmac,eldip,elfco)
        !
        ! Add the stabilization terms
        !
        if (1 == 2) then

           call sld_stabil(&
                itask,pnode,pevat,pgaus,pmate,gpvol,gpsha,eldis,elvel,&
                elmuu,elmat,elfda,elrhs)

        end if
        !
        ! Assemble push forward
        !
        if( kfl_gdepo /= 0 ) then
           call matrix_assemble_element_RHS(&
                ndime*ndime,ndime*ndime,pnode,lnods(:,ielem),elepo,gdepo)
           call matrix_assemble_element_RHS(&
                ndime*ndime,ndime*ndime,pnode,lnods(:,ielem),elepo_new_inverse,gdeinv)
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elepo_new_determinant,gdedet)
        end if

        !
        ! ASSEMBLY
        !
        ! Explicit/Implicit
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elrhs,rhsid)
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elfrx,frxid_sld)
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elfin,finte_sld)
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elfex,fexte_sld)
        if ( kfl_conta_sld /= 0 ) then
           call matrix_assemble_element_RHS(&
                ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elfco,fcont_sld)
        end if

        ! Dynamic only
        if( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elmuu,vmass_sld)
        end if

        ! Implicit only
        if( itask == 2_ip ) then
           call matrix_assemble_element_RHS(&
                ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elmac,macce_sld)
           call matrix_assemble_element_matrix_to_CSR(&
                kfl_element_to_csr,solve(1) % ndofn,pnode,pevat,&
                ielem,lnods(:,ielem),elmat,r_dom,c_dom,amatr,lezdo)
        end if

        ! X-FEM
        if( itask == 1_ip .and. kfl_xfeme_sld == 1_ip ) then
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elmaa,vmasx_sld)
        end if

     else if( itask == 3 ) then
        !
        ! Compute caust, green and lepsi, called in sld_fotens
        !
        call sld_bitcul(&
             itask,ielem,pnode,lnods(1,ielem),pgaus,pmate,gpsha,gpvol,&
             gpgdi_sld(ielem)%a,gppio_sld(ielem)%a,gpfibe,dedef_sld(ielem)%a)

     else if( itask == 4 ) then
        !
        ! Fibers on element,
        ! needs gpfib, computed in sld_elmcla
        !
        call sld_bitcul(&
             itask,ielem,pnode,lnods(1,ielem),pgaus,pmate,gpsha,gpvol,&
             gpgdi_sld(ielem)%a,gppio_sld(ielem)%a,gpfibe,dedef_sld(ielem)%a)

     else if( itask == 5 ) then
        !
        ! Fibers on nodes,
        ! needs gpfib, computed in sld_elmcla
        !
        call sld_bitcul(&
             itask,ielem,pnode,lnods(1,ielem),pgaus,pmate,gpsha,gpvol,&
             gpgdi_sld(ielem)%a,gppio_sld(ielem)%a,gpfibe,dedef_sld(ielem)%a)

     else if( itask == 6 ) then
        !
        ! Update crack elements based on maximum tensile stress criteria
        !
        call sld_crapro(&
             ielem,pgaus,pnode,gpcod,elcod,gpgdi_sld(ielem)%a,&
             gppio_sld(ielem)%a,dedef_sld(ielem)%a)

     else if( itask == 7 ) then
        !
        ! Compute nopio_sld (first Piola--Kirchhoff stress) at nodes
        !
        call sld_bitcul(&
             itask,ielem,pnode,lnods(1,ielem),pgaus,pmate,gpsha,gpvol,&
             gpgdi_sld(ielem)%a,gppio_sld(ielem)%a,gpfibe,dedef_sld(ielem)%a)

     else if( itask == 8 ) then
        !
        ! Exact solution: compute errors
        !
        if( kfl_exacs_sld /= 0 ) then
           call sld_elmexa(&
                2_ip,pgaus,pnode,elcod,gpsha,gpcar,gpdds,gpvol,&
                eldis,gpdis,gpgdi,elrhs)
        end if

     else if( itask == 9 ) then
        !
        ! Exact solution: Assemble force term
        !
        if( kfl_exacs_sld /= 0 ) then
           elrhs = 0.0_rp
           call sld_elmexa(&
                1_ip,pgaus,pnode,elcod,gpsha,gpcar,gpdds,gpvol,&
                eldis,gpdis,gpgdi,elrhs)
           call matrix_assemble_element_RHS(&
                ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elrhs,rhsid)
           !call assrhs(ndofn_sld,pnode,lnods(1,ielem),elrhs,rhsid)
        end if

     else if( itask == 10 ) then
        !
        ! Compute element-wise cauchy Schwarz
        !
        call sld_bitcul(&
             itask,ielem,pnode,lnods(1,ielem),pgaus,pmate,gpsha,gpvol,&
             gpgdi_sld(ielem)%a,gppio_sld(ielem)%a,gpfibe,dedef_sld(ielem)%a)

     else if( itask == 11 ) then
        !
        ! Compute recovered stresses: element-wise cauchy stress
        !
        call sld_bitcul(&
             itask,ielem,pnode,lnods(1,ielem),pgaus,pmate,gpsha,gpvol,&
             gpgdi_sld(ielem)%a,gppio_sld(ielem)%a,gpfibe,dedef_sld(ielem)%a)

     else if( itask == 12 ) then
        !
        ! LOCAL CAUCHY STRESS TENSOR
        !
        call sld_bitcul(&
             itask,ielem,pnode,lnods(1,ielem),pgaus,pmate,gpsha,gpvol,&
             gpgdi_sld(ielem)%a,gppio_sld(ielem)%a,gpfibe,dedef_sld(ielem)%a)

     else if( itask == 13 ) then
        !
        ! INFINITESIMAL STRAIN TENSOR
        !
        call sld_bitcul(&
             itask,ielem,pnode,lnods(1,ielem),pgaus,pmate,gpsha,gpvol,&
             gpgdi_sld(ielem)%a,gppio_sld(ielem)%a,gpfibe,dedef_sld(ielem)%a)

     else if( itask == 14 ) then
        !
        ! INFINITESIMAL STRESS TENSOR (CAUCHY)
        !
        call sld_bitcul(&
             itask,ielem,pnode,lnods(1,ielem),pgaus,pmate,gpsha,gpvol,&
             gpgdi_sld(ielem)%a,gppio_sld(ielem)%a,gpfibe,dedef_sld(ielem)%a)

     end if

  else if ( lelch(ielem) == ELINT ) then
     !
     ! INTERFACE COHESIVE ELEMENTS
     !
     pnode = nnode(pelty)
     pevat = ndofn_sld*pnode
     pmate = lmate_sld(ielem)
     !
     ! Compatibility checks
     !
     if ( pelty /= QUA04 .and. pelty /= HEX08 ) then
        call runend('SLD_ELEMENT_OPERATIONS: WRONG ELEMENT TYPE FOR INTERFACE ELEMENTS')
     end if
     !
     ! Gather operations
     !
     call sld_elmgat(pnode,lnods(1,ielem),eldis,elddis,eldip,elvel,elacc,elcod,eldix,elvex,elmof,elrst)
     !
     ! Element characteristic lenght
     !
     celen_sld(ielem) = 1.0_rp
     !
     ! Call cohesive element elemental operations
     !
     call ELINT_elemental_operations( itask, dtime, ielem, ndime, nsvar_sld, pnode, pevat, elcod(:,:), &
          eldis(:,:,ITER_K), elrhs(:), elmat(:,:), elfrx(:,:) )
     !
     ! Assembly
     !
     if ( itask <=2 ) then
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elrhs,rhsid)
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elrhs,finte_sld)
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elfrx,frxid_sld)

        if( itask == 2 ) then
           call matrix_assemble_element_matrix_to_CSR(&
                kfl_element_to_csr,solve(1) % ndofn,pnode,pevat,&
                ielem,lnods(:,ielem),elmat,r_dom,c_dom,amatr,lezdo)
        end if

     end if

  elseif (lawst_sld(lmate(ielem)) == 200_ip) then
     !
     ! CONTINUUM SHELLS
     !
     pnode = nnode(abs(pelty))
     pevat = ndofn_sld * pnode
     pmate = lmate_sld(ielem)
     pgaus = ngaus(pelty)        ! we need it to calculate mass matrix
     !
     ! Gather operations
     !
     call sld_elmgat(pnode,lnods(1,ielem),eldis,elddis,eldip,elvel,elacc,elcod,eldix,elvex,elmof,elrst)
     !
     ! Cartesian derivatives (PROVISIONAL only for the calculation of the mass matrix)
     !
     call elmca2(&  ! we need it to calculate an artificial mass matrix for stabilize
          pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
          elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
          gpder,gpcar,gphes,ielem)
     !
     ! Element characteristic lenght
     !
     call elmlen(&
          ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
          hnatu(pelty),hleng)
     celen_sld(ielem) = hleng(ndime)

     if (itask <= 2) then
        !
        ! Call shell element elemental operations
        ! >> variables provisionals: pgaus, gpvol, gpsha
        call SHELL_elemental_operations(itask, ielem, elcod, eldis(:,:,ITER_K), elmas, elmuu, &
             elfex, elfin, elmac, elrhs, elfrx, elmat, pgaus, gpvol, gpsha)
        !
        ! Add stabilization terms
        !
        call sld_stabil(&
             itask,pnode,pevat,pgaus,pmate,gpvol,gpsha,eldis,elvel,&
             elmuu,elmat,elfda,elrhs)
        !
        ! Assembly
        !
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elrhs,rhsid)
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elfrx,frxid_sld)
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elfin,finte_sld)
        call matrix_assemble_element_RHS(&
             ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elfex,fexte_sld)

        if( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elmuu,vmass_sld)
        end if

        if( itask == SLD_IMPLICIT_SCHEME ) then
           call matrix_assemble_element_matrix_to_CSR(&
                kfl_element_to_csr,solve(1) % ndofn,pnode,pevat,&
                ielem,lnods(:,ielem),elmat,r_dom,c_dom,amatr,lezdo)
           if( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
              call matrix_assemble_element_RHS(&
                   ndofn_sld,ndofn_sld,pnode,lnods(:,ielem),elmac,macce_sld)
           end if
        end if

     end if

  end if

end subroutine sld_element_operations
