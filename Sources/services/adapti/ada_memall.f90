subroutine ada_memall(itask)
!-----------------------------------------------------------------------
!****f* adapti/ada_memall
! NAME 
!    ada_memall
! DESCRIPTION
!    This routine allocates memory for the module arrays
! USES
!    memchk
!    mediso
! USED BY
!    ada_turnon
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      def_solver
  use      mod_memchk

  use      def_adapti


  implicit none
  integer(ip) :: istat,iimmo,nelau,npoau,itask,ielem,inode


  if (itask == 1) then
     
     nelto_ada=0
     npoto_ada=0
     do iimmo=1,nimmo_ada
        if (nelto_ada < nelei_ada(iimmo)) nelto_ada = nelei_ada(iimmo)
        if (npoto_ada < npoii_ada(iimmo)) npoto_ada = npoii_ada(iimmo)
     end do

     nenew_ada = 1
     nechi_ada = 1  ! max number of child elements per parent element

     if(kfl_paral/=0) then
        !
        ! Slaves and sequential run: allocate all the memory required
        !
        nzerr= 1
        npoau = npoin
        if (kfl_redgr_ada >= 1) then
           if (ndime == 2) then
              nechi_ada = 6
           else if (ndime == 3) then
              call runend('ADA_MEMALL: 3D NOT YET PROGRAMMED!')
           end if
           nzerr= nelem          ! <-- to be used in kernel's memunk
           nenew_ada = nelem * 4 * (ndime-1)
           npoau= nedge
           ! lsofa : to which father a son belongs
           allocate(lsofa_ada(nelem),stat=istat)       
           call memchk(zero,istat,mem_servi(1:2,modul),'LSOFA_ADA','ada_memall',lsofa_ada)
           ! lfaso : which are the nechi_ada sons that each father has
           allocate(lfaso_ada(nechi_ada,nelem),stat=istat)       
           call memchk(zero,istat,mem_servi(1:2,modul),'LFASO_ADA','ada_memall',lfaso_ada)
        end if
     
        allocate(conew_ada(ndime,npoau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'CONEW_ADA','ada_memall',conew_ada)
        allocate(lnnew_ada(mnode,nenew_ada),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNNEW_ADA','ada_memall',lnnew_ada)
        allocate(lnoed_ada(mnode,mnode+3,nelem),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNOED_ADA','ada_memall',lnoed_ada)

        allocate(lnseg_ada(nelem),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNSEG_ADA','ada_memall',lnseg_ada)
        allocate(lnodi_ada(nnodi_ada,nelto_ada,nimmo_ada),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNODI_ADA','ada_memall',lnodi_ada)
        allocate(coori_ada(ndime,npoto_ada,nimmo_ada),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'COORI_ADA','ada_memall',coori_ada)
        allocate(gpoic_ada(mnode,npoto_ada,nimmo_ada),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'GPOIC_ADA','ada_memall',gpoic_ada)
        allocate(gshac_ada(mnode,npoto_ada,nimmo_ada),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'GSHAC_ADA','ada_memall',gshac_ada)
        allocate(lelim_ada(npoto_ada,nimmo_ada),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LELIM_ADA','ada_memall',lelim_ada)
        
        allocate(lnori_ada(mnode,nelem),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNORI_ADA','ada_memall',lnori_ada)
        allocate(lsori_ada(nelem),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LSORI_ADA','ada_memall',lsori_ada)

        if (kfl_redgr_ada == 0) allocate(lenew_ada(nelem),stat=istat)
        
        !
        ! lnori initialization: it stores the original lnods
        !
        npori_ada = npoin
        neori_ada = nelem
        do ielem=1,nelem
           do inode=1,mnode
              lnori_ada(inode,ielem) = lnods(inode,ielem)
           end do
        end do
        
     else
        !
        ! Master: allocate minimum memory
        !
        nelau = 1
        npoau = 1
        allocate(conew_ada(ndime,npoau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'CONEW_ADA','ada_memall',conew_ada)
        allocate(lnnew_ada(mnode,nenew_ada),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNNEW_ADA','ada_memall',lnnew_ada)
        allocate(lnoed_ada(mnode,mnode+3,nelau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNOED_ADA','ada_memall',lnoed_ada)

        allocate(lsofa_ada(nelau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LSOFA_ADA','ada_memall',lsofa_ada)
        allocate(lfaso_ada(nechi_ada,nelau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LFASO_ADA','ada_memall',lfaso_ada)

        allocate(lnseg_ada(nelau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNSEG_ADA','ada_memall',lnseg_ada)
        allocate(lnodi_ada(nnodi_ada,nelau,nelau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNODI_ADA','ada_memall',lnodi_ada)
        allocate(coori_ada(ndimi_ada,npoau,npoau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'COORI_ADA','ada_memall',coori_ada)
        allocate(gpoic_ada(mnode,npoau,npoau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'GPOIC_ADA','ada_memall',gpoic_ada)
        allocate(gshac_ada(mnode,npoau,npoau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'GSHAC_ADA','ada_memall',gshac_ada)
        allocate(lelim_ada(npoau,npoau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LELIM_ADA','ada_memall',lelim_ada)
        
        allocate(lnori_ada(mnode,nelau),stat=istat)
        call memchk(zero,istat,mem_servi(1:2,modul),'LNORI_ADA','ada_memall',lnori_ada)

        if (kfl_redgr_ada == 0) allocate(lenew_ada(nelau),stat=istat)
        
     end if
     
     coori_ada = 0.0_rp
     lnodi_ada = 0
     gshac_ada = 0.0_rp
     gpoic_ada = 0
     lnseg_ada = 0

  else if (itask == 2) then
     !
     ! First set of auxiliary vectors
     !
     allocate(lnosu_ada(nnodb_ada+1,nelem),stat=istat)
     call memchk(zero,istat,mem_servi(1:2,modul),'LNOSU_ADA','ada_memall',lnosu_ada)
     allocate(lolne_ada(4,npoin),stat=istat)
     call memchk(zero,istat,mem_servi(1:2,modul),'LOLNE_ADA','ada_memall',lolne_ada)

  else if (itask == 3) then
     !
     ! Second set of auxiliary vectors
     !
     allocate(lnoad_ada(nnodb_ada+1,nelsu_ada),stat=istat)
     call memchk(zero,istat,mem_servi(1:2,modul),'LNOAD_ADA','ada_memall',lnoad_ada)

  else if (itask == 4) then
  
  else if (itask == 5) then
     !
     ! Interpolation
     !
     allocate(venew_ada(ndime,npnew_ada),stat=istat)
     allocate(scnew_ada(      npnew_ada),stat=istat)     

  end if

end subroutine ada_memall
