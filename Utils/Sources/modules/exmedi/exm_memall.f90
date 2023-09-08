subroutine exm_memall
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_memall
  ! NAME 
  !    exm_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    equations of EXMEDI
  !    Allocation is done for:
  !      - Potentials (intra- and extracellular, and transmembrane)
  !      - Subcell and no-subcell models of ionic currents
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !    exm_turnon
  !***
  !-----------------------------------------------------------------------
  !
  ! DESCRIPTION OF VARIABLES' MEANING 
  !
  !-----------------------------------------------------------------------
  !
  ! Variable: elmag
  ! 
  !     Potentials (i.e. the unknowns) are stored in elmag(npoin,ncomp_exm):
  !        - internal action potential (v) ---> elmag(:,:)
  !        - external action potential (u) ---> extac_exm(:,:) TO BE DONE
  !        - FHN recovery potential        ---> refhn_exm(:,:)
  !
  !-----------------------------------------------------------------------
  !
  ! Variable: vauxi_exm
  !
  !     vauxi_exm(nauxi_exm,npoin,ncomp_exm)
  !        It is used in exm_ceauxi.f90 subroutine for computation of activation/
  !        inactivation states corresponding to Hodgkin-Huxley formalism
  !
  !-----------------------------------------------------------------------
  !
  ! Variable: 

  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      def_solver
  use      def_exmedi
  use mod_memory

  implicit none
  integer(ip)    :: lodof,ncdof,nzrhs_exm,nauxi
  integer(4)     :: istat,imate

  ! FOR SUBCELLULAR IONIC CURRENTS MODELS (TT, LR, BR, ...)

  ncdof= 1
  nzrhs_exm     = 1
  ngate_exm = 1    ! Default for fhn

  !
  ! Derived dimensions according to the general algorithm type
  !
  if(kfl_genal_exm==1 .or. kfl_genal_exm==2)  then
     lodof     = 1                              ! explicit or decoupled
  else if(kfl_genal_exm==3)  then
     lodof     = ndofn_exm                      ! monolithic
  end if

  nauxi_exm = 1
  nconc_exm = 1
  nicel_exm = 1
  ncdof     = 1    

  do imate= 1,nmate

     if (kfl_cellmod(imate) == CELL_NOMOD_EXMEDI) then
        nauxi_exm = max(nauxi_exm,1_ip)    ! Default value of number of variables used for activation/inactivation 
        nconc_exm = max(nconc_exm,1_ip)    ! Default value of number of variables used for concentrations
        nicel_exm = max(nicel_exm,1_ip)
     else if (kfl_cellmod(imate) == CELL_FENTON_EXMEDI) then
        call runend("EXM_MEMALL: FENTON KARMA MODELS NOT UPDATED")
        nconc_exm = 2
        ngate_exm = 2
     else if (kfl_cellmod(imate) == CELL_FITZHUGH_EXMEDI) then
        nauxi_exm = max(nauxi_exm,12_ip)  
        nconc_exm = max(nconc_exm, 8_ip)  
        ncdof = npoin
     else if (kfl_cellmod(imate) == CELL_TT2006_EXMEDI) then
        nauxi_exm = max(nauxi_exm,12_ip)  
        nconc_exm = max(nconc_exm,10_ip)  
        nicel_exm = max(nconc_exm,18_ip)  
        ncdof = npoin
     else if (kfl_cellmod(imate) == CELL_OHARA_EXMEDI) then
        nauxi_exm = max(nauxi_exm,29_ip)  
        nconc_exm = max(nconc_exm,11_ip)  
        nicel_exm = max(nconc_exm,26_ip)  
        ncdof = npoin
     else if (kfl_cellmod(imate) == CELL_SCATRIA_EXMEDI) then
        nauxi_exm = max(nauxi_exm,14_ip)  
        nconc_exm = max(nconc_exm,3_ip)  
        nicel_exm = max(nconc_exm,15_ip)  
        ncdof = npoin
     else if (kfl_cellmod(imate) == CELL_SCVENTRI_EXMEDI) then
        nauxi_exm = max(nauxi_exm,14_ip)  
        nconc_exm = max(nconc_exm,3_ip)  
        nicel_exm = max(nconc_exm,15_ip)  
        ncdof = npoin
!     else if (kfl_cellmod(imate) == 6) then  !TTINA
!        nauxi_exm = max(nauxi_exm,14_ip)  
!        nconc_exm = max(nconc_exm,10_ip)  
!        nicel_exm = max(nconc_exm,18_ip)  
!        ncdof = npoin
     end if

  end do


  if(kfl_timei_exm==1) ncomp_exm = 2 + kfl_tiacc_exm

  ! kfl_gemod_exm = 0 --> No-Subcellular models for ionic currents: FHN, FENTON...
  ! kfl_gemod_exm = 1 --> Subcellular models for ionic currents: TT, LR, BR, ...

  if(INOTMASTER) then


     ! ALLOCATION of variables:
     !   - elmag
     !   - vconc
     !   - taulo
     !   - vauxi_exm
     !   - vicel_exm
     !   - ticel_exm
     !   - jicel_exm
     !   - rhsau_exm
     !   - nospm_exm

     !!     if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then
     call memory_alloca(mem_modul(1:2,modul),'TAULO','exm_memall',taulo,npoin)
     taulo= 0.0_rp

     call memory_alloca(mem_modul(1:2,modul),'FISOC','exm_memall',fisoc,npoin)
     fisoc = -1.0_rp
     if( thiso_exm(1)>0.0_rp ) then !if ISHOCH AUTOSAVE
        fisoc = 0.0_rp !here we will save both upstrokes with positive time and downstrokes with negative time. There will never be activation at 0 time, so it's safe
     else
        fisoc = -1.0_rp
     end if

     call memory_alloca(mem_modul(1:2,modul),'KWAVE_EXM','exm_memall',kwave_exm,npoin)
     kwave_exm = 0_ip

     allocate(isoch_modified(npoin), STAT=istat)
     if (istat /= 0) then
        call runend('EXMEDI: Insufficient memory for isoch_modified')
     end if
     isoch_modified = .FALSE.


     call memory_alloca(mem_modul(1:2,modul),'VDIAG_EXM','exm_memall',vdiag_exm,npoin)
     call memory_alloca(mem_modul(1:2,modul),'IDIMA_EXM','exm_memall',idima_exm,npoin)

     call memory_alloca(mem_modul(1:2,modul),'ELMAG','exm_memall',elmag,npoin,ncomp_exm)
     call memory_alloca(mem_modul(1:2,modul),'REFHN_EXM','exm_memall',refhn_exm,npoin,ncomp_exm)
     call memory_alloca(mem_modul(1:2,modul),'VAUXI_EXM','exm_memall',vauxi_exm,nauxi_exm,ncdof,ncomp_exm)
     call memory_alloca(mem_modul(1:2,modul),'VCONC'    ,'exm_memall',vconc    ,nconc_exm,ncdof,ncomp_exm)
     call memory_alloca(mem_modul(1:2,modul),'CELL_CA0_ECC','exm_memall',cell_ca0_ecc ,3_ip,nmate)
     call memory_alloca(mem_modul(1:2,modul),'VICEL_EXM','exm_memall',vicel_exm,nicel_exm,ncdof,ncomp_exm)

     call memory_alloca(mem_modul(1:2,modul),'TICEL_EXM','exm_memall',ticel_exm,npoin)
     call memory_alloca(mem_modul(1:2,modul),'JICEL_EXM','exm_memall',jicel_exm,npoin)
     call memory_alloca(mem_modul(1:2,modul),'QNETO','exm_memall',qneto_exm,npoin)
     qneto_exm = 0_rp

     if (kfl_appli_exm > 0) then
        !
        ! Source fields
        !
        !    Applied current field flag:
        call memory_alloca(mem_modul(1:2,modul),'LAPNO_EXM','exm_memall',lapno_exm,npoin)
        lapno_exm=0

        !    Applied current field:
        call memory_alloca(mem_modul(1:2,modul),'APPFI_EXM','exm_memall',appfi_exm,npoin)
        appfi_exm=0.0_rp
     end if


     if(kfl_algso_exm==0 .and. kfl_genal_exm > 1) then
        call runend('NOT CODED')
        !allocate(lpexm(npoin),stat=istat)
        !call memchk(zero,istat,memdi,'LPEXM','exm_memall',lpexm)
        !call mediso(lodof,npoin,lun_solve_exm,lpexm)
     end if

     !
     ! Residual size
     !
     nzrhs_exm=npoin
     !
     ! Conductivity, etc.
     !
     if( kfl_comdi_exm == 2 ) then

        call memory_alloca(mem_modul(1:2,modul),'CEDIF_EXM','exm_memall',cedif_exm,ndime,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAFI_EXM','exm_memall',grafi_exm,ndime,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'KGRFI_EXM','exm_memall',kgrfi_exm,nelem)

     else if( kfl_comdi_exm == 3 ) then

        call memory_alloca(mem_modul(1:2,modul),'KGRFI_EXM','exm_memall',kgrfi_exm,nelem)

     end if

  else

     !
     ! Master: allocate minimum memory
     !

     call memory_alloca(mem_modul(1:2,modul),'VDIAG_EXM','exm_memall',vdiag_exm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'ELMAG','exm_memall',elmag,1_ip,ncomp_exm)
     call memory_alloca(mem_modul(1:2,modul),'REFHN_EXM','exm_memall',refhn_exm,1_ip,ncomp_exm)
     call memory_alloca(mem_modul(1:2,modul),'VAUXI_EXM','exm_memall',vauxi_exm,nauxi_exm,1_ip,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'VCONC','exm_memall',vconc,nconc_exm,1_ip,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'CELL_CA0_ECC','exm_memall',cell_ca0_ecc ,3_ip,nmate)
     call memory_alloca(mem_modul(1:2,modul),'VICEL_EXM','exm_memall',vicel_exm,nicel_exm,1_ip,1_ip)

     call memory_alloca(mem_modul(1:2,modul),'TICEL_EXM','exm_memall',ticel_exm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'JICEL_EXM','exm_memall',jicel_exm,1_rp)
     call memory_alloca(mem_modul(1:2,modul),'QNETO','exm_memall',qneto_exm,1_rp)

     if (kfl_appli_exm > 0) then
        !
        ! Source fields
        !
        !    Applied current field flag:
        call memory_alloca(mem_modul(1:2,modul),'LAPNO_EXM','exm_memall',lapno_exm,1_ip)
        lapno_exm=0

        !    Applied current field:
        call memory_alloca(mem_modul(1:2,modul),'APPFI_EXM','exm_memall',appfi_exm,1_ip)
        appfi_exm=0.0_rp
     end if


     if(kfl_algso_exm==0) then
        call runend('NOT CODED')
        !allocate(lpexm(1),stat=istat)
        !call memchk(zero,istat,memdi,'LPEXM','exm_memall',lpexm)
        !call mediso(lodof,npoin,lun_solve_exm,lpexm)
     end if
     !
     ! Conductivity, etc.
     !
     call memory_alloca(mem_modul(1:2,modul),'CEDIF_EXM','exm_memall',cedif_exm,ndime,ndime,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'GRAFI_EXM','exm_memall',grafi_exm,ndime,ndime,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'KGRFI_EXM','exm_memall',kgrfi_exm,1_ip)

  end if

  !
  ! Actualize maximum sizes of RHSID
  !


  nzrhs=max(nzrhs,nzrhs_exm)
  !
  ! Solver memory for NVARI variables
  !
  solve_sol => solve(1:1)
  call soldef(4_ip)
  

end subroutine exm_memall
