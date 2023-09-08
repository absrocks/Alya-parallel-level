!-----------------------------------------------------------------------
!> @addtogroup NastalTurnon
!> @ingroup    Nastal
!> @{
!> @file    nsa_memall.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Allocate memory 
!> @details Allocate memory 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_memall
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      def_solver
  use      mod_memchk

  use      def_nastal


  implicit none
  integer(ip) :: knodb,iboun,lodof,ielem,inode,ipoin,nvara_nsa,nrhsa_nsa,ndire,iauxi,nauxi,i,pelty,pgaus
  integer(4)  :: istat

  !
  ! Boundary conditions
  !
  !solve(1) % bvess     => bvess_nsa(:,:,1)
  !solve(1) % kfl_fixno => kfl_fixno_nsa

  ! Element operations database, defined for all

  



  if(INOTMASTER) then
     nvara_nsa=1
     nrhsa_nsa=1
     if (ntomp_nsa>0) nvara_nsa= nunkn_nsa
     if (nromp_nsa>0) nrhsa_nsa= nunkn_nsa

     allocate(rhsax_nsa(nrhsa_nsa,nromp_nsa+1,nfrap_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RHSAX_NSA','nsa_memall',rhsax_nsa)

     allocate(rhsou_nsa(npoin*ndofn_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RHSOU_NSA','nsa_memall',rhsou_nsa)
     rhsou_nsa = 0.0_rp

     allocate(reafo_nsa(npoin*ndofn_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'REAFO_NSA','nsa_memall',reafo_nsa)
     reafo_nsa = 0.0_rp

     !RK:
     allocate(unkna_nsa(npoin*ndofn_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'UNKNA_NSA','nsa_memall',unkna_nsa)
     allocate(unkit_nsa(npoin*ndofn_nsa,nromp_nsa+1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'UNKIT_NSA','nsa_memall',unkit_nsa)

     !
     ! Residual smoothing
     !
     nauxi=1
     if (kfl_resmo_nsa == 1) nauxi= npoin*ndofn_nsa
     allocate(resmo_nsa(nauxi),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RESMO_NSA','nsa_memall',resmo_nsa)
     resmo_nsa= 0.0_rp

     !
     ! Shock capturing 
     !
     nauxi=1
     do i=1,ndofn_nsa
        if (kfl_shock_nsa(i) > 0) nauxi= nelem
     end do
     allocate(shocktau_nsa(nauxi),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'SHOCKTAU_NSA','nsa_memall',shocktau_nsa)
     do ielem = 1,nauxi
        pelty = abs(ltype(ielem))
        pgaus = ngaus(pelty)
        allocate(shocktau_nsa(ielem)%a(3,pgaus,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'SHOCKTAU_NSA','nsa_memall',shocktau_nsa(ielem)%a)
     end do
     
     !
     ! Delta unknown
     !
     allocate(dunkn_nsa(ndofn_nsa*npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DUNKN_NSA','nsa_memall',dunkn_nsa)
     dunkn_nsa = 0.0_rp
     
     ! Added for FSI
!      allocate(nsa_forfsi(ndime,npoin),stat=istat)
!      call memchk(zero,istat,mem_modul(1:2,modul),'NSA_FORFSI','nsa_forfsi',nsa_forfsi)
     allocate(veloc(ndime,npoin,ncomp_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VELOC','nsa_memall',veloc)
     allocate(press(npoin,ncomp_nsa)      ,stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PRESS','nsa_memall',press)         
     allocate(vorti(ndime+1,npoin)        ,stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VORTI','nsa_memall',vorti)

     allocate(umome(ndime,npoin,ncomp_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'UMOME','nsa_memall',umome)
     
     allocate(densi(npoin,ncomp_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DENSI','nsa_memall',densi)
     allocate(tempe(npoin,ncomp_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'TEMPE','nsa_memall',tempe)
     allocate(energ(npoin,ncomp_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ENERG','nsa_memall',energ)
     allocate(visco(npoin,ncomp_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VISCO','nsa_memall',visco)
     allocate(vmach_nsa(npoin) ,stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VMACH_NSA','nsa_memall',vmach_nsa)
  
     if (kfl_infun_nsa > 0 .and. kfl_benme_nsa >= 200) then
        !
        ! Meteo problems
        !
        ncol_nsa = npoin/(nelz_nsa+1)

        allocate(node_column_nsa(ncol_nsa,nelz_nsa+1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'node_column_nsa','nsa_memall',node_column_nsa)
        allocate(intma_column_nsa(ncol_nsa,nelz_nsa+1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'intma_column_nsa','nsa_memall',intma_column_nsa)

        !Kessler's variables:
        allocate(rhocol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rhocol_nsa','nsa_memall',rhocol_nsa)
        allocate(tcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'tcol_nsa','nsa_memall',tcol_nsa)
        allocate(qvcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qvcol_nsa','nsa_memall',qvcol_nsa)
        allocate(qccol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qccol_nsa','nsa_memall',qccol_nsa)
        allocate(qrcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qrcol_nsa','nsa_memall',qrcol_nsa)
        allocate(vtcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'vtcol_nsa','nsa_memall',vtcol_nsa)
        allocate(prodcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'prodcol_nsa','nsa_memall',prodcol_nsa)
        allocate(prodkcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'prodkcol_nsa','nsa_memall',prodkcol_nsa)
        allocate(vtdencol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'vtdencol_nsa','nsa_memall',vtdencol_nsa)
        allocate(rdzkcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rdzkcol_nsa','nsa_memall',rdzkcol_nsa)
        allocate(rdzwcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rdzwcol_nsa','nsa_memall',rdzwcol_nsa)
        allocate(rhokcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rhokcol_nsa','nsa_memall',rhokcol_nsa)
        allocate(factorcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'factorcol_nsa','nsa_memall',factorcol_nsa)

        allocate(q_nsa(ndime+6,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'q_nsa','nsa_memall',q_nsa)
        allocate(qref_nsa(ndime+6,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qref_nsa','nsa_memall',qref_nsa)
        allocate(q_col_nsa(ndime+6,nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'q_col_nsa','nsa_memall',q_col_nsa)
        allocate(qref_col_nsa(ndime+6,nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qref_col_nsa','nsa_memall',qref_col_nsa)

        allocate(pcol_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'pcol_nsa','nsa_memall',pcol_nsa)
        allocate(z_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'z_nsa','nsa_memall',z_nsa)
        allocate(z_col_nsa(nz_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'z_col_nsa','nsa_memall',z_col_nsa)

        allocate(rainnc_nsa(ncol_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rainnc_nsa','nsa_memall',rainnc_nsa)
        allocate(rainncv_nsa(ncol_nsa),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rainncv_nsa','nsa_memall',rainncv_nsa)

        allocate(prodk_nsa(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'prodk_nsa','nsa_memall',prodk_nsa)

     end if !METEO Variables

     ndire= 1
     if (kfl_relat_nsa == 1) ndire= npoin           ! Relativistic flow
        
     allocate(cdens(ndire,ncomp_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CDENS','nsa_memall',cdens)
     allocate(entha(ndire,ncomp_nsa),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ENTHA','nsa_memall',entha)        
     
     
     ! Allocate time increments vector
     
     allocate(dtieq_nsa(ndtdf_nsa,npoin,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DTIEQ_NSA','nsa_memall',dtieq_nsa)
     dtieq_nsa = 0.0_rp

     ! Allocate convergence fields (the fields that correspond to rensi)

     allocate(crens_nsa(ndtdf_nsa,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CRENS_NSA','nsa_memall',crens_nsa)
     crens_nsa = 0.0_rp

     if (kfl_pro2d_nsa > 0) then
        allocate(linod_nsa(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LINOD_NSA','nsa_memall',linod_nsa)        
     end if

     if (kfl_fasts_nsa /=0) then        
        allocate(vmacp_nsa(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'VMACP_NSA','nsa_memall',vmacp_nsa)
     end if

!     iauxi= kfl_dttyp_nsa(1)+kfl_dttyp_nsa(2)+kfl_dttyp_nsa(3)
!     if (kfl_diagi_nsa > 0 .or. kfl_lotim_nsa > 0 .or. iauxi > 0) then        
     allocate(vdiag_nsa(ndofn_nsa*npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VDIAG_NSA','nsa_memall',vdiag_nsa)
!     end if

     allocate(umoss_nsa(ndime,npoin,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'UMOSS_NSA','nsa_memall',umoss_nsa)
     allocate(denss_nsa(npoin,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DENSS_NSA','nsa_memall',denss_nsa)
     allocate(eness_nsa(npoin,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ENESS_NSA','nsa_memall',eness_nsa)
     
     allocate(umosg_nsa(ndime,nelem,mgaus,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'UMOSG_NSA','nsa_memall',umosg_nsa)
     allocate(densg_nsa(nelem,mgaus,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DENSG_NSA','nsa_memall',densg_nsa)
     allocate(enesg_nsa(nelem,mgaus,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ENESG_NSA','nsa_memall',enesg_nsa)
 
     allocate(ortpr_nsa(ndofn_nsa,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ORTPR_NSA','nsa_memall',ortpr_nsa)
    
     allocate(frequ_nsa(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'FREQU_NSA','nsa_memall',frequ_nsa)
     
     if (kfl_coupl(ID_NASTAL,ID_TURBUL) == 0 ) then
         allocate(turmu(npoin),stat=istat)
         call memchk(zero,istat,mem_modul(1:2,modul),'TURMU','nsa_memall',turmu)        
         turmu     = 0.0_rp
     endif

     if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
         allocate(wmean(npoin,1),stat=istat)
         call memchk(zero,istat,mem_modul(1:2,modul),'WMEAN','nsa_memall',wmean)        
         wmean     = 0.0_rp
     endif

     allocate(shecp_nsa(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'SHECP_NSA','nsa_memall',shecp_nsa)

     umosg_nsa = 0.0_rp
     densg_nsa = 0.0_rp
     enesg_nsa = 0.0_rp
     
     umoss_nsa = 0.0_rp
     denss_nsa = 0.0_rp
     eness_nsa = 0.0_rp
     
     ortpr_nsa = 0.0_rp
     frequ_nsa = 0.0_rp

     shecp_nsa = 0.0_rp

     allocate(avpre_nsa(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVPRE_NSA','nsa_memall',avpre_nsa)
     allocate(avtem_nsa(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVTEM_NSA','nsa_memall',avtem_nsa)

     allocate(avvel_nsa(ndime,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVVEL_NSA','nsa_memall',avvel_nsa)
     allocate(avve2_nsa(ndime,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVVE2_NSA','nsa_memall',avve2_nsa)
     allocate(avvxy_nsa(ndime,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVVXY_NSA','nsa_memall',avvxy_nsa)
     allocate(avmom_nsa(ndime,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVMOM_NSA','nsa_memall',avmom_nsa)

     avpre_nsa = 0.0_rp
     avtem_nsa = 0.0_rp
     avvel_nsa = 0.0_rp
     avve2_nsa = 0.0_rp
     avvxy_nsa = 0.0_rp
     avmom_nsa = 0.0_rp

     !------------------------------------------------------------------------
     ! Convection for convection-diffusion-reaction equation
     !------------------------------------------------------------------------
     !
     ! Dimensions required for the mesh adaptivity (when required)
     !
     if (nzerr > 1) nzerr= max(nzerr, nelem*ndofn_nsa)
     !
     ! Dimensions required for the solver (when required)
     !
     solve_sol => solve(1:1)
     call soldef(4_ip)
     !
     ! Boundary conditions
     !
     solve(1) % bvess     => bvess_nsa(:,:,1)
     solve(1) % kfl_fixno => kfl_fixno_nsa
     
  else
     !
     ! Master: allocate minimum memory
     !
     allocate(veloc(1,1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VELOC','nsa_memall',veloc)
     allocate(press(1,3)      ,stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PRESS','nsa_memall',press)         
     allocate(vorti(1,3)      ,stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VORTI','nsa_memall',vorti)

     allocate(umome(1,1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'UMOME','nsa_memall',umome)
     allocate(densi(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DENSI','nsa_memall',densi)
     allocate(tempe(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'TEMPE','nsa_memall',tempe)
     allocate(energ(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ENERG','nsa_memall',energ)
     allocate(visco(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VISCO','nsa_memall',visco)  
     allocate(vmach_nsa(1) ,stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VMACH_NSA','nsa_memall',vmach_nsa)
     allocate(resmo_nsa(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RESMO_NSA','nsa_memall',resmo_nsa)

     allocate(dunkn_nsa(1) ,stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DUNKN_NSA','nsa_memall',dunkn_nsa)
     dunkn_nsa = 0.0_rp

     allocate(reafo_nsa(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'REAFO_NSA','nsa_memall',reafo_nsa)
     reafo_nsa = 0.0_rp

     allocate(rhsou_nsa(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RHSOU_NSA','nsa_memall',rhsou_nsa)
     rhsou_nsa = 0.0_rp

     nauxi=1
     allocate(shocktau_nsa(nauxi),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'SHOCKTAU_NSA','nsa_memall',shocktau_nsa)
     allocate(shocktau_nsa(1)%a(3,1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'SHOCKTAU_NSA','nsa_memall',shocktau_nsa(1)%a)
     

     if (kfl_infun_nsa > 0 .and. kfl_benme_nsa >= 200) then
        !
        ! Meteo problems
        !
        
        allocate(node_column_nsa(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'node_column_nsa','nsa_memall',node_column_nsa)
        allocate(intma_column_nsa(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'intma_column_nsa','nsa_memall',intma_column_nsa)

        !Kessler's variables:
        allocate(rhocol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rhocol_nsa','nsa_memall',rhocol_nsa)
        allocate(tcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'tcol_nsa','nsa_memall',tcol_nsa)
        allocate(qvcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qvcol_nsa','nsa_memall',qvcol_nsa)
        allocate(qccol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qccol_nsa','nsa_memall',qccol_nsa)
        allocate(qrcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qrcol_nsa','nsa_memall',qrcol_nsa)
        allocate(vtcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'vtcol_nsa','nsa_memall',vtcol_nsa)
        allocate(prodcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'prodcol_nsa','nsa_memall',prodcol_nsa)
        allocate(prodkcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'prodkcol_nsa','nsa_memall',prodkcol_nsa)
        allocate(vtdencol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'vtdencol_nsa','nsa_memall',vtdencol_nsa)
        allocate(rdzkcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rdzkcol_nsa','nsa_memall',rdzkcol_nsa)
        allocate(rdzwcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rdzwcol_nsa','nsa_memall',rdzwcol_nsa)
        allocate(rhokcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rhokcol_nsa','nsa_memall',rhokcol_nsa)
        allocate(factorcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'factorcol_nsa','nsa_memall',factorcol_nsa)

        allocate(q_nsa(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'q_nsa','nsa_memall',q_nsa)
        allocate(qref_nsa(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qref_nsa','nsa_memall',qref_nsa)
        allocate(q_col_nsa(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'q_col_nsa','nsa_memall',q_col_nsa)
        allocate(qref_col_nsa(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'qref_col_nsa','nsa_memall',qref_col_nsa)

        allocate(pcol_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'pcol_nsa','nsa_memall',pcol_nsa)
        allocate(z_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'z_nsa','nsa_memall',z_nsa)
        allocate(z_col_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'z_col_nsa','nsa_memall',z_col_nsa)

        allocate(rainnc_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rainnc_nsa','nsa_memall',rainnc_nsa)
        allocate(rainncv_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rainncv_nsa','nsa_memall',rainncv_nsa)
        
        allocate(prodk_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'prodk_nsa','nsa_memall',prodk_nsa)

     end if !METEO Variables minimum allocation for master only
     

     allocate(cdens(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CDENS','nsa_memall',cdens)
     allocate(entha(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ENTHA','nsa_memall',entha) 
       
     allocate(dtieq_nsa(1,1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DTIEQ_NSA','nsa_memall',dtieq_nsa)
     allocate(crens_nsa(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CRENS_NSA','nsa_memall',crens_nsa)
     
     allocate(umoss_nsa(1,1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'UMOSS_NSA','nsa_memall',umoss_nsa)
     allocate(denss_nsa(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DENSS_NSA','nsa_memall',denss_nsa)
     allocate(eness_nsa(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ENESS_NSA','nsa_memall',eness_nsa)
     
     allocate(umosg_nsa(1,1,1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'UMOSG_NSA','nsa_memall',umosg_nsa)
     allocate(densg_nsa(1,1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DENSG_NSA','nsa_memall',densg_nsa)
     allocate(enesg_nsa(1,1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ENESG_NSA','nsa_memall',enesg_nsa)

     allocate(vdiag_nsa(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VDIAG_NSA','nsa_memall',vdiag_nsa)
     
     allocate(ortpr_nsa(3,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ORTPR_NSA','nsa_memall',ortpr_nsa)

     allocate(frequ_nsa(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'FREQU_NSA','nsa_memall',frequ_nsa)

     if (kfl_coupl(ID_NASTAL,ID_TURBUL) == 0 ) then
         allocate(turmu(1),stat=istat)
         call memchk(zero,istat,mem_modul(1:2,modul),'TURMU','nsa_memall',turmu)        
     endif

     if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
         allocate(wmean(npoin,1),stat=istat)
         call memchk(zero,istat,mem_modul(1:2,modul),'WMEAN','nsa_memall',wmean)        
     endif

     allocate(shecp_nsa(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'SHECP_NSA','nsa_memall',shecp_nsa)

     allocate(avpre_nsa(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVPRE_NSA','nsa_memall',avpre_nsa)
     allocate(avtem_nsa(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVTEM_NSA','nsa_memall',avtem_nsa)

     allocate(avvel_nsa(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVVEL_NSA','nsa_memall',avvel_nsa)
     allocate(avve2_nsa(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVVE2_NSA','nsa_memall',avve2_nsa)
     allocate(avvxy_nsa(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVVXY_NSA','nsa_memall',avvxy_nsa)
     allocate(avmom_nsa(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'AVMOM_NSA','nsa_memall',avmom_nsa)

  end if
  !
  ! Actualize maximum sizes of RHSID
  !
  !nzrhs=max(nzrhs,nzrhs_nsa)
 
  !
  ! Subscale coefficients (everyone)
  !
!!$  allocate(cosma_nsa(mgaus,mgaus,2*mfreq_nsa+1,2*mfreq_nsa+1,2*mfreq_nsa+1,nelty),stat=istat)
!!$!  call memchk(zero,istat,mem_modul(1:2,modul),'COSMA_NSA','nsa_memall',cosma_nsa)
!!$  allocate(sinma_nsa(mgaus,mgaus,2*mfreq_nsa+1,2*mfreq_nsa+1,2*mfreq_nsa+1,nelty),stat=istat)
!!$!  call memchk(zero,istat,mem_modul(1:2,modul),'SINMA_NSA','nsa_memall',sinma_nsa) 

  allocate(cosma_nsa(mgaus,mgaus,2*frmax_nsa+1,2*frmax_nsa+1,2*frmax_nsa+1,nelty),stat=istat)
!  call memchk(zero,istat,mem_modul(1:2,modul),'COSMA_NSA','nsa_memall',cosma_nsa)
  allocate(sinma_nsa(mgaus,mgaus,2*frmax_nsa+1,2*frmax_nsa+1,2*frmax_nsa+1,nelty),stat=istat)
!  call memchk(zero,istat,mem_modul(1:2,modul),'SINMA_NSA','nsa_memall',sinma_nsa) 


end subroutine nsa_memall
