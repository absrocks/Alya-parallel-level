!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_solmem.f90
!> @date    14/06/2019
!> @author  hozeaux
!> @brief   Allocate memory
!> @details Allocate memory for solver and variables that should be
!>          reallocated and not redistributed in case of repartitioning
!> @}
!-----------------------------------------------------------------------

subroutine sld_solmem()

  use def_kintyp,               only : ip
  use def_master,               only : IMASTER
  use def_master,               only : ID_EXMEDI
  use def_master,               only : solve
  use def_master,               only : press
  use def_master,               only : modul, mem_modul
  use def_master,               only : gdepo,gdeinv,gdedet
  use def_master,               only : gpfib
  use def_master,               only : coupling
  use def_master,               only : kfl_modul
  use def_domain,               only : ndime,npoin,nelem,ngaus,nmate
  use def_domain,               only : ltype
  use def_domain,               only : mgaus
  use def_solver,               only : solve_sol
  use mod_memory,               only : memory_alloca
  use mod_output_postprocess,   only : output_postprocess_check_variable_postprocess
  use mod_exm_sld_eccoupling
  use def_solidz

  implicit none

  integer(ip) :: ielem,pgaus,idime,ipoin
  integer(ip) :: istat

  if( kfl_rigid_sld == 0 ) then

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
     ! Global variables
     !
     !----------------------------------------------------------------------
     !
     ! Force vectors
     !
     call memory_alloca(mem_modul(1:2,modul),'FINTE_SLD','sld_solmem',finte_sld,ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'FINTT_SLD','sld_solmem',fintt_sld,ndime,npoin,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'FEXTE_SLD','sld_solmem',fexte_sld,ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'FEXTT_SLD','sld_solmem',fextt_sld,ndime,npoin,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'MACCE_SLD','sld_solmem',macce_sld,ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'FRXID_SLD','sld_solmem',frxid_sld,ndime*npoin)
     if( kfl_conta_sld /= 0_ip ) then
        call memory_alloca(mem_modul(1:2,modul),'FCONT_SLD','sld_solmem',fcont_sld,ndime*npoin)
     end if
     !
     ! Mass matrix
     !
     call memory_alloca(mem_modul(1:2,modul),'VMASS_SLD',   'sld_solmem',vmass_sld,npoin)
     !
     ! Variables solution methods
     !
     call memory_alloca(mem_modul(1:2,modul),'DUNKN_SLD',   'sld_solmem',dunkn_sld,   ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'UNKNOTMP_SLD','sld_solmem',unknotmp_sld,ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'VELOCTMP_SLD','sld_solmem',veloctmp_sld,ndime,npoin)
     call memory_alloca(mem_modul(1:2,modul),'DDISP_SLD',   'sld_solmem',ddisp_sld,   ndime,npoin,ncomp_sld)
     !
     ! Element characteristic length
     !
     call memory_alloca(mem_modul(1:2,modul),'CELEN_SLD','sld_solmem',celen_sld,nelem)
     if ( kfl_damag_sld == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'ELENG_SLD','sld_solmem',eleng_sld,nelem)
     end if
     !
     ! Push forward operator
     !
     if( kfl_gdepo /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'GDEPO', 'sld_solmem',gdepo, ndime,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GDEINV', 'sld_solmem',gdeinv, ndime,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GDEDET','sld_solmem',gdedet,npoin)
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
     ! Energy
     !
     call memory_alloca(mem_modul(1:2,modul),'ALLWK_SLD','sld_solmem',allwk_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ALLIE_SLD','sld_solmem',allie_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ALLKE_SLD','sld_solmem',allke_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ETOTA_SLD','sld_solmem',etota_sld,ncomp_sld)

     !----------------------------------------------------------------------
     !
     ! Post-process
     !
     !----------------------------------------------------------------------
     !
     ! Postprocess: Stress and strains vectors
     !
     call memory_alloca(mem_modul(1:2,modul),'EPSEL_SLD','sld_solmem',epsel_sld,nelem)
     call memory_alloca(mem_modul(1:2,modul),'LEPSE_SLD','sld_solmem',lepse_sld,nelem)
     do ielem = 1, nelem
        pgaus = ngaus(abs(ltype(ielem)))
        call memory_alloca(mem_modul(1:2,modul),'EPSEL_SLD(IELEM)','sld_solmem',epsel_sld(ielem)%a,ndime,ndime,pgaus)
        call memory_alloca(mem_modul(1:2,modul),'LEPSE_SLD(IELEM)','sld_solmem',lepse_sld(ielem)%a,ndime,ndime,pgaus)
     end do
     !
     ! Postprocess: Stress and strains
     !
     if( kfl_foten_sld == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'CAUST_SLD','sld_solmem',caust_sld,nvoig_sld,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GREEN_SLD','sld_solmem',green_sld,nvoig_sld,npoin)
        call memory_alloca(mem_modul(1:2,modul),'LEPSI_SLD','sld_solmem',lepsi_sld,nvoig_sld,npoin)
        call memory_alloca(mem_modul(1:2,modul),'SEQVM_SLD','sld_solmem',seqvm_sld,npoin)
     else
        call memory_alloca(mem_modul(1:2,modul),'CAUST_SLD','sld_solmem',caust_sld,nvoig_sld,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'GREEN_SLD','sld_solmem',green_sld,nvoig_sld,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'LEPSI_SLD','sld_solmem',lepsi_sld,nvoig_sld,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SEQVM_SLD','sld_solmem',seqvm_sld,1_ip)
     end if
     !
     ! Postprocess: Fibers
     !
     call memory_alloca(mem_modul(1:2,modul),'FIBDE_SLD','sld_solmem',fibde_sld,ndime,npoin)

      !----------------------------------------------------------------------
      !
      ! Couplings
      !
      !----------------------------------------------------------------------
      !
      ! Coupling with NASTIN
      !
      if( coupling('SOLIDZ','NASTIN') >= 1 ) then
         ! the mesh velocity is the solid deformation velocity
         ! velom => veloc_sld(:,:,1)
         call memory_alloca(mem_modul(1:2,modul),'DDISM_SLD','sld_solmem',ddism_sld,ndime,npoin)
      end if

      call memory_alloca(mem_modul(1:2,modul),'PRESS',   'sld_solmem',press,npoin,ncomp_sld)
      call memory_alloca(mem_modul(1:2,modul),'FORC_FSI','sld_solmem',forc_fsi,ndime,npoin)
      !
      ! Coupling with EXMEDI
      !
      if( has_exmsld_coupling() .or. kfl_exmsld_3Dcou_ecc )then
          call exm_sld_ecc_allocate_memmory( 101_ip )
         if( kfl_gdepo == 0_ip )then
            call runend('SLD_SOLMEM: ROTATION ON SHOULD BE INCLUDED IN SOLIDZ-EXMEDI SIMULATIONS ')
         endif
      else
         if( kfl_modul(ID_EXMEDI) == 1_ip )then
            call runend('SLD_SOLMEM: COUPLING WITH EXMEDI NOW IN THE KER.DAT FILE!! CHECK DOCS.')
         endif
      end if


     !----------------------------------------------------------------------
     !
     ! Others: (Requires revision/deletion)
     !
     !----------------------------------------------------------------------

     if( IMASTER )then
         call memory_alloca(mem_modul(1:2,modul),'GRLST_SLD','sld_solmem',grlst_sld,1_ip,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'EPRIN_SLD','sld_solmem',eprin_sld,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'CAUNN_SLD','sld_solmem',caunn_sld,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'SIGEI_SLD','sld_solmem',sigei_sld,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'ROLOC_SLD','sld_solmem',roloc_sld,1_ip,1_ip,1_ip)
     else
         call memory_alloca(mem_modul(1:2,modul),'GPGDI_SLD','sld_solmem',gpgdi_sld,nelem)
         call memory_alloca(mem_modul(1:2,modul),'GPPIO_SLD','sld_solmem',gppio_sld,nelem)
         call memory_alloca(mem_modul(1:2,modul),'DEDEF_SLD','sld_solmem',dedef_sld,nelem)
         call memory_alloca(mem_modul(1:2,modul),'GRLST_SLD','sld_solmem',grlst_sld,nelem,ndime*ndime)
         call memory_alloca(mem_modul(1:2,modul),'EPRIN_SLD','sld_solmem',eprin_sld,npoin)
         call memory_alloca(mem_modul(1:2,modul),'CAUNN_SLD','sld_solmem',caunn_sld,npoin)
         call memory_alloca(mem_modul(1:2,modul),'SIGEI_SLD','sld_solmem',sigei_sld,npoin)
         call memory_alloca(mem_modul(1:2,modul),'ROLOC_SLD','sld_solmem',roloc_sld,ndime,ndime,npoin)
         do ielem = 1,nelem
            pgaus = ngaus(abs(ltype(ielem)))
            call memory_alloca(mem_modul(1:2,modul),'GPGDI_SLD(IELEM)','sld_solmem',gpgdi_sld(ielem)%a,ndime,ndime,pgaus)
            call memory_alloca(mem_modul(1:2,modul),'GPPIO_SLD(IELEM)','sld_solmem',gppio_sld(ielem)%a,ndime,ndime,pgaus)
            call memory_alloca(mem_modul(1:2,modul),'DEDEF_SLD','sld_solmem',dedef_sld(ielem)%a,pgaus)
         end do
     endif
     call memory_alloca(mem_modul(1:2,modul),'AUX','sld_solmem',aux_sld,ndime,npoin,2_ip)
     allocate(gpsl0_sld(ndime,ndime,30),stat=istat); gpsl0_sld=0.0_rp
     allocate(rstr0_sld(ndime,ndime,30),stat=istat); rstr0_sld=0.0_rp

  end if

end subroutine sld_solmem
