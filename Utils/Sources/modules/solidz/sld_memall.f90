!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_memall.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   General memory allocation
!> @details General memory allocation
!> @}
!-----------------------------------------------------------------------

subroutine sld_memall()

  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver, only : solve_sol
  use mod_memchk
  use mod_memory
  use def_solidz

  implicit none
  integer(ip) :: lodof,ielem,pgaus,pelty,ntaul,num_land
  integer(ip) :: ncouv,ncomo,idime,igaus,ivari,imate,npoin_roloc
  integer(ip) :: ipoin
  integer(4)  :: istat

  !----------------------------------------------------------------------
  !
  ! Solver
  !
  !----------------------------------------------------------------------
  !
  ! Memory
  !
  solve_sol => solve(1:2)
  call soldef(4_ip)
  !
  ! Boundary conditions
  !
  solve_sol => solve(1:)
  solve(1) % bvess     => bvess_sld(:,:,1)
  solve(1) % bvnat     => bvnat_sld(:,:,1)
  solve(1) % kfl_fixno => kfl_fixno_sld
  !
  ! Dirichlet value is 0 (we solve increments)
  !
  if( solve(1) % kfl_iffix /= 0 ) solve(1) % kfl_iffix = 2

  !----------------------------------------------------------------------
  !
  ! Arrays
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then
     !
     ! Displacement, velocity and acceleration
     !
     call memory_alloca(mem_modul(1:2,modul),'DISPL',    'sld_memall',    displ,ndime,npoin,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'VELOC_SLD','sld_memall',veloc_sld,ndime,npoin,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ACCEL_SLD','sld_memall',accel_sld,ndime,npoin,ncomp_sld)
     !
     ! Correction variables for Newton-Raphson
     !
     call memory_alloca(mem_modul(1:2,modul),'DUNKN_SLD',   'sld_memall',dunkn_sld,ndofn_sld*npoin)
     call memory_alloca(mem_modul(1:2,modul),'DDISP_SLD',   'sld_memall',ddisp_sld,ndime,npoin,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'VELOCTMP_SLD','sld_memall',veloctmp_sld,ndime,npoin)
     call memory_alloca(mem_modul(1:2,modul),'UNKNOTMP_SLD','sld_memall',unknotmp_sld,ndime*npoin)
     !
     ! Mass matrix
     !
     call memory_alloca(mem_modul(1:2,modul),'VMASS_SLD','sld_memall',vmass_sld,npoin)
     !
     ! Global/assembled forces
     !
     call memory_alloca(mem_modul(1:2,modul),'FRXID_SLD','sld_memall',frxid_sld,ndofn_sld*npoin)
     call memory_alloca(mem_modul(1:2,modul),'FINTE_SLD','sld_memall',finte_sld,ndofn_sld*npoin)
     call memory_alloca(mem_modul(1:2,modul),'FEXTE_SLD','sld_memall',fexte_sld,ndofn_sld*npoin)
     call memory_alloca(mem_modul(1:2,modul),'MACCE_SLD','sld_memall',macce_sld,ndofn_sld*npoin)
     if( kfl_conta_sld /= 0_ip ) then
        call memory_alloca(mem_modul(1:2,modul),'FCONT_SLD','sld_membcs',fcont_sld,ndime*npoin)
     end if
     call memory_alloca(mem_modul(1:2,modul),'FEXTT_SLD','sld_memall',fextt_sld,ndime,npoin,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'FINTT_SLD','sld_memall',fintt_sld,ndime,npoin,ncomp_sld)
     !
     ! Energies
     !
     call memory_alloca(mem_modul(1:2,modul),'ALLWK_SLD','sld_memall',allwk_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ALLIE_SLD','sld_memall',allie_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ALLKE_SLD','sld_memall',allke_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ETOTA_SLD','sld_memall',etota_sld,ncomp_sld)
     !
     ! State dependent variables (stored at Gauss Points for each Element)
     !
     if( kfl_sdvar_sld == 1_ip ) then
        call memory_alloca(mem_modul(1:2,modul),'SVEGM_SLD','sld_memall',svegm_sld,nelem)
        call memory_alloca(mem_modul(1:2,modul),'SVAUX_SLD','sld_memall',svaux_sld,nelem) ! Postprocess purposes for GID
        do ielem = 1, nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'SVEGM_SLD(IELEM)','sld_memall',svegm_sld(ielem)%a,nsvar_sld,pgaus,2_ip)
           call memory_alloca(mem_modul(1:2,modul),'SVAUX_SLD(IELEM)','sld_memall',svaux_sld(ielem)%a,     1_ip,pgaus,1_ip)
        end do
     end if
     !
     ! Element characteristic length
     !
     call memory_alloca(mem_modul(1:2,modul),'CELEN_SLD','sld_memall',celen_sld,nelem)
     if ( kfl_xfeme_sld == 1 .or. kfl_damag_sld == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'ELENG_SLD','sld_memall',eleng_sld,nelem)
     end if
     !
     ! Stress and strains vectors
     !
     call memory_alloca(mem_modul(1:2,modul),'EPSEL_SLD','sld_memall',epsel_sld,nelem)
     call memory_alloca(mem_modul(1:2,modul),'LEPSE_SLD','sld_memall',lepse_sld,nelem)
     do ielem = 1, nelem
        pgaus = ngaus(abs(ltype(ielem)))
        call memory_alloca(mem_modul(1:2,modul),'EPSEL_SLD(IELEM)','sld_memall',epsel_sld(ielem)%a,ndime,ndime,pgaus)
        call memory_alloca(mem_modul(1:2,modul),'LEPSE_SLD(IELEM)','sld_memall',lepse_sld(ielem)%a,ndime,ndime,pgaus)
     end do
     !
     ! Push forward operator
     !
     if( kfl_gdepo /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'GDEPO', 'sld_memall',gdepo, ndime,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GDEINV','sld_memall',gdeinv,ndime,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GDEDET','sld_memall',gdedet,npoin)
        !
        ! Initialize gdepo to the identity (useful for FSI coupling)
        !
        do ipoin = 1_ip, npoin
           do idime = 1_ip, ndime
              gdepo(idime,idime,ipoin)  = 1.0_rp
              gdeinv(idime,idime,ipoin) = 1.0_rp
           end do
           gdedet(ipoin) = 1.0_rp
        end do
     end if
     !
     ! Coupling with NASTIN
     !
     if( coupling('SOLIDZ','NASTIN') >= 1 ) then
        ! the mesh velocity is the solid deformation velocity
        !velom => veloc_sld(:,:,1)
        call memory_alloca(mem_modul(1:2,modul),'DDISM_SLD','sld_memall',ddism_sld,ndime,npoin)
     end if
     !
     ! Coupling with EXMEDI
     !
     ntaul= 1
     ncouv= 1
     ncomo= 1
     num_land = 1
     if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then
        if( kfl_gdepo == 0 ) call runend('SLD_MEMALL: ROTATION ON SHOULD BE INCLUDED IN SOLIDZ-EXMEDI SIMULATIONS ')
        ntaul= npoin
        do imate= 1,nmate_sld
           !only one variable coupling is enough to allocate covar_sld
           if (kfl_coupt_sld(imate) == 11) ncouv = nelem
           if ( kfl_eccty(imate) ==  3 .or. kfl_eccty(imate)==4 ) ncomo = nelem
        end do
     else
        if (kfl_modul(ID_EXMEDI) == 1) &
             call runend('SLD_MEMALL: COUPLING WITH EXMEDI NOW IN THE KER.DAT FILE!! CHECK DOCS.')
     end if
     !
     ! Stuff mostly related to exmedi coupling
     !
     call memory_alloca(mem_modul(1:2,modul),'KACTI_SLD','sld_memall',kacti_sld,ntaul)
     call memory_alloca(mem_modul(1:2,modul),'STRETLAM','sld_memall',stretlam,ntaul)
     call memory_alloca(mem_modul(1:2,modul),'COVAR_SLD','sld_memall',covar_sld,ncouv)
     call memory_alloca(mem_modul(1:2,modul),'STATELAND','sld_memall',stateland,6_ip,ncomo,mgaus,2_ip)
     call memory_alloca(mem_modul(1:2,modul),'EXM_LAMBDA_SLD','sld_memall',exm_lambda_sld,ncomo,mgaus,2_ip)
     covar_sld(1:ncouv) = 1.0_rp
     stateland(1:2,:,:,1) = 0.0_rp
     stateland(3,:,:,1)   = 0.00000001_rp
     stateland(4,:,:,1)   = 1.0_rp
     stateland(5:6,:,:,1) = 0.0_rp
     exm_lambda_sld(:,:,:) = 1.0_rp
     covar_sld(1:ncouv) = 1.0_rp
     num_land = 1
     do imate= 1,nmate
        if (kfl_eccty(imate) == 3 .or. kfl_eccty(imate)==4 ) num_land = npoin
     end do
     call memory_alloca(mem_modul(1:2,modul),'TROPONIN','sld_memall',troponin,num_land)
     call memory_alloca(mem_modul(1:2,modul),'TROPONIN_PREV','sld_memall',troponin_prev,num_land)
     !
     ! Others
     !
     call memory_alloca(mem_modul(1:2,modul),'PRESS','sld_memall',press,npoin,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'GPFIB','sld_memall',gpfib,ndime,mgaus,nelem)
     call memory_alloca(mem_modul(1:2,modul),'FORFI','sld_memall',forc_fsi,ndime,npoin)
     call memory_alloca(mem_modul(1:2,modul),'VDIAG_SLD','sld_memall',vdiag_sld,npoin)
     call memory_alloca(mem_modul(1:2,modul),'AUX','sld_memall',aux_sld,ndime,npoin,2_ip)
     aux_sld = 0.0_rp
     !
     ! disep_sld(ndime,npoin) is the perturbation displacement, only required for inexact newton
     !
     if (kfl_ninex_sld == 1) then
        call memory_alloca(mem_modul(1:2,modul),'DISEP_SLD','sld_memall',disep_sld,ndime,npoin)
     else
        call memory_alloca(mem_modul(1:2,modul),'DISEP_SLD','sld_memall',disep_sld,ndime, 1_ip)
     end if
     !
     ! Local pseudo time steps (tau)
     !
     call memory_alloca(mem_modul(1:2,modul),'DTTAU_SLD','sld_memall',dttau_sld,nelem)
     !
     ! Displacement enrichment function for XFEM / GFEM
     !
     if (kfl_xfeme_sld > 0) then
        call memory_alloca(mem_modul(1:2,modul),'DXFEM_SLD','sld_memall',dxfem_sld,ndime,npoin,ncomp_sld)
        call memory_alloca(mem_modul(1:2,modul),'VXFEM_SLD','sld_memall',vxfem_sld,ndime,npoin,ncomp_sld)
        call memory_alloca(mem_modul(1:2,modul),'AXFEM_SLD','sld_memall',axfem_sld,ndime,npoin,ncomp_sld)
        call memory_alloca(mem_modul(1:2,modul),'CRAPX_SLD','sld_memall',crapx_sld,ndime,nelem)
        call memory_alloca(mem_modul(1:2,modul),'CRANX_SLD','sld_memall',cranx_sld,ndime,nelem)
        call memory_alloca(mem_modul(1:2,modul),'SGMAX_SLD','sld_memall',sgmax_sld,nelem)
        call memory_alloca(mem_modul(1:2,modul),'LEENR_SLD','sld_memall',leenr_sld,nelem)
        call memory_alloca(mem_modul(1:2,modul),'LNENR_SLD','sld_memall',lnenr_sld,npoin)
        call memory_alloca(mem_modul(1:2,modul),'VMASX_SLD','sld_memall',vmasx_sld,npoin)
     end if
     !
     ! Postprocess of element wise Cauchy stress
     !
     ivari = 36
     call posdef(25_ip,ivari)
     if( ivari /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'CAUSE_SLD','sld_memall',cause_sld,nelem)
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'CAUSE_SLD(IELEM)','sld_memall',cause_sld(ielem)%a,ndime,pgaus,1_ip)
        end do
     end if
     !
     ! Cohesive law and friction/contact
     !
     if (kfl_cohes_sld > 0 .and. kfl_elcoh > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'DFRIC_SLD','sld_memall',dfric_sld,nelem,2_ip,mgaus,2_ip)
        call memory_alloca(mem_modul(1:2,modul),'DSLIP_SLD','sld_memall',dslip_sld,nelem,2_ip,mgaus,2_ip)
        call memory_alloca(mem_modul(1:2,modul),'DCEFF_SLD','sld_memall',dceff_sld,nelem,2_ip,mgaus,2_ip)
        call memory_alloca(mem_modul(1:2,modul),'DCMAX_SLD','sld_memall',dcmax_sld,nelem,2_ip,mgaus)
        call memory_alloca(mem_modul(1:2,modul),'NOPIO_SLD','sld_memall',nopio_sld,ndime*ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'LECOH_SLD','sld_memall',lecoh_sld,nelem)
        call memory_alloca(mem_modul(1:2,modul),'TREFF_SLD','sld_memall',treff_sld,npoin)
        call memory_alloca(mem_modul(1:2,modul),'NOCOH_SLD','sld_memall',nocoh_sld,npoin)
        call memory_alloca(mem_modul(1:2,modul),'COHNX_SLD','sld_memall',cohnx_sld,ndime,nelem)
     end if
     !
     ! Damage model
     !
     if (kfl_damag_sld > 0) then
        allocate(ledam_sld(nelem),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LEDAM_SLD','sld_memall',ledam_sld)
     else if (kfl_damag_sld == 0) then
        allocate(ledam_sld(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LEDAM_SLD','sld_memall',ledam_sld)
     end if

     !
     ! Allocate memory for stress (1st PK tensor)
     !
     allocate(gppio_sld(nelem),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'GPPIO_SLD','sld_memall',gppio_sld)
     if (kfl_xfeme_sld > 0) then
        do ielem = 1,nelem
           allocate(gppio_sld(ielem)%a(ndime,ndime,mgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'GPPIO_SLD','sld_memall',gppio_sld(ielem)%a)
        end do
     else
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           allocate(gppio_sld(ielem)%a(ndime,ndime,pgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'GPPIO_SLD','sld_memall',gppio_sld(ielem)%a)

        end do
     end if
     !
     ! Allocate memory for F (F1 stocked in gpgdi_sld)
     !
     allocate(gpgdi_sld(nelem),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'GPGDI_SLD','sld_memall',gpgdi_sld)
     if (kfl_xfeme_sld > 0) then
        do ielem = 1,nelem
           allocate(gpgdi_sld(ielem)%a(ndime,ndime,mgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'GPGDI_SLD','sld_memall',gpgdi_sld(ielem)%a)
        end do
     else
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           allocate(gpgdi_sld(ielem)%a(ndime,ndime,pgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'GPGDI_SLD','sld_memall',gpgdi_sld(ielem)%a)
        end do
     end if
     !
     ! Allocate memory for detF
     !
     allocate(dedef_sld(nelem),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DEDEF_SLD','sld_memall',dedef_sld)
     if (kfl_xfeme_sld > 0) then
        do ielem = 1,nelem
           allocate(dedef_sld(ielem)%a(mgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'DEDEF_SLD','sld_memall',dedef_sld(ielem)%a)
        end do
     else
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           allocate(dedef_sld(ielem)%a(pgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'DEDEF_SLD','sld_memall',dedef_sld(ielem)%a)
        end do
     end if
     !
     ! Dimensions required for the solver
     !
     lodof     = ndofn_sld
     nusol     = max(nusol,1_ip)
     mxdof     = max(mxdof,lodof)
     nusol     = 1
     nzmat_sld = 0
     nzrhs_sld = ndime*npoin
     if (kfl_xfeme_sld == 1) nzrhs_sld = nzrhs_sld + ndime*npoin

     !if(kfl_timet_sld==2) then
     !   if (kfl_xfeme_sld == 1) then
     !      solve_sol => solve(1:2) ! OJO CAMBIAR
     !      call soldef(4_ip)
     !   else
     !      solve_sol => solve(1:2) ! OJO CAMBIAR
     !      call soldef(4_ip)
     !   end if
     !end if
     !
     ! NEW
     !
     !solve_sol => solve(1:2) ! OJO CAMBIAR
     !call soldef(4_ip)
     !
     ! Initialized F to [1]
     !
     if (kfl_xfeme_sld > 0) then
        do ielem = 1,nelem
           do idime=1,ndime
              do igaus=1,mgaus
                 gpgdi_sld(ielem)%a(idime,idime,igaus)=1.0_rp
              end do
           end do
        end do
     else
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           do idime=1,ndime
              do igaus=1,pgaus
                 gpgdi_sld(ielem)%a(idime,idime,igaus)=1.0_rp
              end do
           end do
        end do
     end if

     !
     ! Cauchy stress and Green strain: CAUST_SLD, GREEN_SLD, LSFSN_SLD, SEQVM_SLD, INV1/2/3E_SLD' (for postprocessing stage)
     !
     if (kfl_foten_sld == 1) then
        allocate(caust_sld(nvoig_sld,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'CAUST_SLD','sld_memall',caust_sld)
        allocate(caunn_sld(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'CAUNN_SLD','sld_memall',caunn_sld)
        allocate(green_sld(nvoig_sld,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GREEN_SLD','sld_memall',green_sld)
        allocate(lepsi_sld(nvoig_sld,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LSFSN_SLD','sld_memall',lepsi_sld)
        allocate(seqvm_sld(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'SEQVM_SLD','sld_memall',seqvm_sld)
        allocate(grlst_sld(nelem,ndime*ndime),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GRLST_SLD','sld_memall',grlst_sld)
     else
        allocate(caust_sld(nvoig_sld,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'CAUST_SLD','sld_memall',caust_sld)
        allocate(caunn_sld(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'CAUNN_SLD','sld_memall',caunn_sld)
        allocate(green_sld(nvoig_sld,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GREEN_SLD','sld_memall',green_sld)
        allocate(lepsi_sld(nvoig_sld,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LSFSN_SLD','sld_memall',lepsi_sld)
        allocate(seqvm_sld(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'SEQVM_SLD','sld_memall',seqvm_sld)
        allocate(grlst_sld(1,ndime*ndime),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GRLST_SLD','sld_memall',grlst_sld)
     end if

     if (kfl_plast_sld == 1) then
        allocate(epsee_sld(nvoig_sld,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'EPSEE_SLD','sld_memall',epsee_sld)
     else
        allocate(epsee_sld(nvoig_sld,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'EPSEE_SLD','sld_memall',epsee_sld)
     end if

     if (kfl_restr_sld == 0) then
        ! No restr field read
        allocate(restr_sld(6,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'RESTR_SLD','sld_memall',restr_sld)
     end if

     if (kfl_fiber_sld < 4_ip) then
        allocate(fibde_sld(ndime,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'FIBDE_SLD','sld_memall',fibde_sld)
        if (kfl_fiber_sld == 3) then
           !
           ! sheet and normal fields to be computed, not read. when read, they are assigned
           ! in SLD_INIVAR
           !
           allocate(fibts_sld(ndime,npoin),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'FIBTS_SLD','sld_memall',fibts_sld)
           allocate(fibtn_sld(ndime,npoin),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'FIBTN_SLD','sld_memall',fibtn_sld)

        end if
     else
        allocate(fibde_sld(ndime,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'FIBDE_SLD','sld_memall',fibde_sld)
     end if

     !!     npoin_roloc = 1
     !!     if (kfl_rotei_sld == 1) npoin_roloc = npoin
     npoin_roloc = npoin
     allocate(sigei_sld(npoin_roloc),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'SIGEI_SLD','sld_memall',sigei_sld)
     allocate(roloc_sld(ndime,ndime,npoin_roloc),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ROLOC_SLD','sld_memall',roloc_sld)

     !
     ! Elemental fibers
     !
     ivari = 10
     call posdef(25_ip,ivari)
     if( ivari /= 0 ) then
        allocate(fibeg_sld(nelem),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'FIBEG_SLD','sld_memall',fibeg_sld)
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           allocate(fibeg_sld(ielem)%a(ndime,pgaus,1),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'FIBEG_SLD','sld_memall',fibeg_sld(ielem)%a)
        end do
     end if
     !
     ! Actualize maximum sizes of matrix and RHS
     !
     nzmat=max(nzmat,nzmat_sld)
     nzrhs=max(nzrhs,nzrhs_sld)

     !
     ! CAUCHY STRESS TENSOR
     !   The values are not average at nodes so they are stored for each
     !   element.
     ivari = 38
     call posdef(25_ip,ivari)
     if( ivari /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'CAUNP_SLD','sld_memall',caunp_sld,nelem)
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'CAUNP_SLD(IELEM)','sld_memall',caunp_sld(ielem)%a,nvoig_sld,pgaus,1_ip)
        end do
     end if

     !
     ! LOCAL CAUCHY STRESS TENSOR                                                  ( *AQU* )
     !   Local stress tensor according to the global coordinate system.
     !   The values are not average at nodes so they are stored for each
     !   element.
     ivari = 39
     call posdef(25_ip,ivari)
     if( ivari /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'CAULO_SLD','sld_memall',caulo_sld,nelem)
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'CAULO_SLD(IELEM)','sld_memall',caulo_sld(ielem)%a,nvoig_sld,pgaus,1_ip)
        end do
     end if
     !
     ! Isochrones
     !
     if(kfl_isoch_sld==1_ip) then
        allocate(isoch_sld(npoin,2), stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'isoch_sld','sld_memall',isoch_sld)
        isoch_sld=-1.0_rp

        allocate(iswav_sld(npoin,2), stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'iswav_sld','sld_memall',iswav_sld)
        iswav_sld=0_ip
     endif

  else
     !
     ! MASTER: Allocate minimum memory
     !
     !
     ! Displacement, velocity and acceleration
     !
     call memory_alloca(mem_modul(1:2,modul),'DISPL',    'sld_memall',displ,    1_ip,1_ip,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'VELOC_SLD','sld_memall',veloc_sld,1_ip,1_ip,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ACCEL_SLD','sld_memall',accel_sld,1_ip,1_ip,ncomp_sld)
     !
     ! Correction variables for Newton-Raphson
     !
     call memory_alloca(mem_modul(1:2,modul),'DDISP_SLD','sld_memall',ddisp_sld,1_ip,1_ip,ncomp_sld)
     call memory_alloca_min(dunkn_sld)
     call memory_alloca(mem_modul(1:2,modul),'VELOCTMP_SLD','sld_memall',veloctmp_sld,1_ip,1_ip)
     call memory_alloca_min(unknotmp_sld)
     !
     ! Force vectors
     !
     call memory_alloca_min(frxid_sld)
     call memory_alloca_min(finte_sld)
     call memory_alloca_min(fexte_sld)
     call memory_alloca_min(macce_sld)
     if( kfl_conta_sld /= 0 ) call memory_alloca_min(fcont_sld)
     call memory_alloca(mem_modul(1:2,modul),'FINTT_SLD','sld_memall',fintt_sld,1_ip,1_ip,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'FEXTT_SLD','sld_memall',fextt_sld,1_ip,1_ip,ncomp_sld)
     !
     ! Energies
     !
     call memory_alloca(mem_modul(1:2,modul),'ALLIE_SLD','sld_memall',allie_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ALLWK_SLD','sld_memall',allwk_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ALLKE_SLD','sld_memall',allke_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ETOTA_SLD','sld_memall',etota_sld,ncomp_sld)
     !
     ! Mass matrix
     !
     call memory_alloca_min(vmass_sld)
     !
     ! State dependent variables (stored at Gauss Points for each Element)
     !
     if( kfl_sdvar_sld == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'SVEGM_SLD','sld_memall',svegm_sld,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SVAUX_SLD','sld_memall',svaux_sld,1_ip)
     end if
     !
     ! Element characteristic length
     !
     call memory_alloca(mem_modul(1:2,modul),'CELEN_SLD','sld_memall',celen_sld,1_ip)
     !
     ! Push forward operator
     !
     if( kfl_gdepo /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'GDEPO','sld_memall',gdepo,ndime,ndime,1_ip)
        call memory_alloca_min(gdeinv)
        call memory_alloca_min(gdedet)
     end if
     !
     ! Other
     !
     ! Cohesive elements pandolfi
     !
     if( kfl_elcoh /= 0 )then
        allocate(dfric_sld(1,1,1,2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'DFRIC_SLD','sld_memall',dfric_sld)
        allocate(dslip_sld(1,1,1,2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'DSLIP_SLD','sld_memall',dslip_sld)
        allocate(dceff_sld(1,1,1,2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'DCEFF_SLD','sld_memall',dceff_sld)
        allocate(dcmax_sld(1,1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'DCMAX_SLD','sld_memall',dcmax_sld)
        allocate(nopio_sld(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'NOPIO_SLD','sld_memall',nopio_sld)
        allocate(lecoh_sld(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LECOH_SLD','sld_memall',lecoh_sld)
        allocate(treff_sld(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'TREFF_SLD','sld_memall',treff_sld)
        call memory_alloca(mem_modul(1:2,modul),'NOCOH_SLD','sld_memall',nocoh_sld,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'COHNX_SLD','sld_memall',cohnx_sld,1_ip,1_ip)
     end if
     !
     ! XFEM
     !
     if( kfl_xfeme_sld /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'DXFEM_SLD','sld_memall',dxfem_sld,1_ip,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'VXFEM_SLD','sld_memall',vxfem_sld,1_ip,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'AXFEM_SLD','sld_memall',axfem_sld,1_ip,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'CRAPX_SLD','sld_memall',crapx_sld,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'CRANX_SLD','sld_memall',cranx_sld,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SGMAX_SLD','sld_memall',sgmax_sld,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'ELENG_SLD','sld_memall',eleng_sld,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'LEENR_SLD','sld_memall',leenr_sld,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'LNENR_SLD','sld_memall',lnenr_sld,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'VDIAG_SLD','sld_memall',vdiag_sld,1_ip)
     end if

     allocate(gpfib(ndime,mgaus,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'GPFIB','sld_memall',gpfib)

     allocate(caust_sld(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CAUST_SLD','sld_memall',caust_sld)
     allocate(green_sld(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'GREEN_SLD','sld_memall',green_sld)
     allocate(lepsi_sld(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LSFSN_SLD','sld_memall',lepsi_sld)
     allocate(grlst_sld(1,ndime*ndime),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'grlst_SLD','sld_memall',grlst_sld)
     allocate(restr_sld(6,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RESTR_SLD','sld_memall',restr_sld)
     allocate(fibde_sld(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'FIBDE_SLD','sld_memall',fibde_sld)
     allocate(seqvm_sld(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'SEQVM_SLD','sld_memall',seqvm_sld)

     allocate(iswav_sld(1,2), stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'iswav_sld','sld_memall',iswav_sld)

     allocate(dttau_sld(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DTTAU_SLD','sld_memall',dttau_sld)

     !! Allocate minimum memory for master for the land model
     call memory_alloca(mem_modul(1:2,modul),'STRETLAM' ,'sld_memall',stretlam,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'STATELAND','sld_memall',stateland,6_ip,1_ip,1_ip,2_ip)
     call memory_alloca(mem_modul(1:2,modul),'TROPONIN','sld_memall',troponin,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'TROPONIN_PREV','sld_memall',troponin_prev,1_ip)

  end if

  !dxfem_sld=0.0_rp
  !vxfem_sld=0.0_rp
  !axfem_sld=0.0_rp
  !crapx_sld=0.0_rp
  !cranx_sld=0.0_rp
  sidev_sld=0.0_rp
  dltap_sld=0.0_rp

  ifase_sld=0_ip
  kfase_sld=9_ip

end subroutine sld_memall
