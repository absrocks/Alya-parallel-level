subroutine nsa_membcs(itask)
!-----------------------------------------------------------------------
!****f* Nastin/nsa_membcs
! NAME 
!    nsa_memcbs
! DESCRIPTION
!    This routine allocates memory for the boundary conditions arrays
! USES
!    ecoute
! USED BY
!    nsa_reabcs
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_nastal
  use      def_domain
  use      mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ndbve
  integer(4)              :: istat
  integer(ip)             :: ifunc
  
  select case (itask)
     
  case(1)
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------!

call allocate_lodi_bou() 
!!! MIGUELITO EL MAYA LO MIRARA
!!$    allocate(normal_nsa(npoin), stat=istat)
!!$    normal_nsa(1:npoin)%id    = -1
!!$    normal_nsa(1:npoin)%idime = -1
!!$    normal_nsa(1:npoin)%ichrc = -1
!!$    !
!!$    allocate( xschlrn_nsa(nelem,mgaus) ,stat=istat)
!!$    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', xschlrn_nsa)
!!$    !
!!$    allocate( schlrn_nsa(npoin) ,stat=istat)
!!$    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', schlrn_nsa)
!!$    !
!!$    allocate( xchrc_nsa(ndofn_nsa,nelem,mgaus,ndime), stat=istat)
!!$    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', xchrc_nsa)
!!$    ! 
!!$    allocate( chrc_nsa(ndofn_nsa,npoin,ndime), stat=istat)
!!$    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', chrc_nsa)
!!$    !
!!$    !
!!$    allocate( gamma_nsa(npoin) ,stat=istat)
!!$    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', gamma_nsa)
!!$    gamma_nsa(1:npoin) = adgam_nsa
!!$    !
!!$    allocate( sound_nsa(npoin) ,stat=istat)
!!$    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', sound_nsa)
!!$    sound_nsa(1:npoin) = 343.2_rp 
!!$    !
!!$    allocate( sxinv_nsa(ndofn_nsa,ndofn_nsa,ndime,mgaus,nelem), stat=istat)
!!$    allocate(    sx_nsa(ndofn_nsa,ndofn_nsa,ndime,mgaus,nelem), stat=istat)
!!$    !-----------------------------------------------------------------------!
     
     allocate(kfl_fixno_nsa(ndofn_nsa,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_NSA','nsa_membcs',kfl_fixno_nsa)
     kfl_fixno_nsa=-1     

     allocate(kfl_fixbo_nsa(nboun),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_NSA','nsa_membcs',kfl_fixbo_nsa)
     
     allocate(jacrot_du_dq_nsa(ndofn_nsa,ndofn_nsa,nbopo),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'JACROT_DU_DQ_NSA','nsa_membcs',jacrot_du_dq_nsa)
     jacrot_du_dq_nsa = 0.0_rp
     allocate(jacrot_dq_du_nsa(ndofn_nsa,ndofn_nsa,nbopo),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'JACROT_DQ_DU_NSA','nsa_membcs',jacrot_dq_du_nsa)
     jacrot_dq_du_nsa = 0.0_rp

     allocate(kfl_fixrs_nsa(nbopo),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXRS_NSA','nsa_membcs',kfl_fixrs_nsa)


  case(2)
     
     if (kfl_infun_nsa > 0) then
        !
        ! Hydrostatic reference field: keep speed, density, pressure and 
        ! temperature (ndofn_nsa + 1) and some other stuff...
        !
        if(kfl_benme_nsa >= 200) then

           !Allocate space for 3 additional water variables:
           allocate(rekee_nsa(nkeep_nsa+3,npoin),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'REKEE_NSA',    'nsa_membcs',rekee_nsa)        
        else
           allocate(rekee_nsa(nkeep_nsa,npoin),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'REKEE_NSA',    'nsa_membcs',rekee_nsa)
        end if
     else
        allocate(rekee_nsa(nkeep_nsa,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'REKEE_NSA',    'nsa_membcs',rekee_nsa)                
     end if

     if ( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or. kfl_benme_nsa == 5 .or. &
          kfl_benme_nsa == 15  .or. kfl_benme_nsa >= 200) then
        !
        ! Hydrostatic reference field: keep speed, density, pressure and 
        ! temperature (ndofn_nsa + 1)
        !
        allocate(bspon_nsa(ndofn_nsa+1,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BSPON_NSA',    'nsa_membcs',bspon_nsa)
        allocate(aa_nsa(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'AA_NSA','nsa_membcs',aa_nsa)
        allocate(bb_nsa(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BB_NSA','nsa_membcs',bb_nsa)
        
     else
        allocate(bspon_nsa(ndofn_nsa+1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BSPON_NSA',    'nsa_membcs',bspon_nsa)    

        allocate(aa_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'AA_NSA','nsa_membcs',aa_nsa)
        allocate(bb_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BB_NSA','nsa_membcs',bb_nsa)
            
     end if
     
     allocate(funpa_nsa(20),       stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'FUNPA_NSA','nsa_membcs',funpa_nsa)

     allocate(bvnat_nsa(ndofn_nsa,nboun),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'BVNAT_NSA',    'nsa_membcs',bvnat_nsa)
     if(kfl_conbc_nsa == 0) then                                   ! non-constant boundary conditions
        allocate(kfl_funno_nsa(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNNO_NSA','nsa_membcs',kfl_funno_nsa)
        allocate(kfl_funbo_nsa(nboun),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNBO_NSA','nsa_reabcs',kfl_funbo_nsa)
        !kfl_funno_nsa=-1     !!! COMMENTED OUT THIS LINE; WHY WAS IT HERE?
        allocate(bvess_nsa(ndofn_nsa,npoin,2),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_NSA',    'nsa_membcs',bvess_nsa)
     else                                                       ! constant boundary conditions
        allocate(kfl_funno_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNNO_NSA','nsa_membcs',kfl_funno_nsa)
        allocate(kfl_funbo_nsa(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FUNBO_NSA','nsa_reabcs',kfl_funbo_nsa)
        !kfl_funno_nsa=-1     !!! COMMENTED OUT THIS LINE; WHY WAS IT HERE?
        if(kfl_benme_nsa >= 200) then
           !
           ! Kessler variables (3 additional vars set in initial conditions):
           !
           allocate(bvess_nsa(ndofn_nsa+3,npoin,1),stat=istat)    
           call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_NSA',    'nsa_membcs',bvess_nsa)
                  
        else
           allocate(bvess_nsa(ndofn_nsa,npoin,1),stat=istat)    
           call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_NSA',    'nsa_membcs',bvess_nsa)
        end if

     end if
     !Initialize to zero
     bvess_nsa     = 0.0_rp

  case (3)

     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXNO_NSA','nsa_membcs',kfl_fixno_nsa)
     deallocate(kfl_fixno_nsa,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXNO_NSA','nsa_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXBO_NSA','nsa_membcs',kfl_fixbo_nsa)
     deallocate(kfl_fixbo_nsa,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXBO_NSA','nsa_membcs',0_ip)     
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXRS_NSA','nsa_membcs',kfl_fixrs_nsa)
     deallocate(kfl_fixrs_nsa,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXRS_NSA','nsa_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'BVNAT_NSA','nsa_membcs',bvnat_nsa)
     deallocate(bvnat_nsa,stat=istat)
     if(istat/=0) call memerr(two,'BVNAT_NSA','nsa_membcs',0_ip)     
     call memchk(two,istat,mem_modul(1:2,modul),'REKEE_NSA',    'nsa_membcs',rekee_nsa)                
     deallocate(rekee_nsa,stat=istat)
     if(istat/=0) call memerr(two,'REKEE_NSA','nsa_membcs',0_ip)     
     
     if(kfl_conbc_nsa == 0) then 
        call memchk(two,istat,mem_modul(1:2,modul),'bvess_nsa','nsa_membcs',bvess_nsa)
        deallocate(bvess_nsa,stat=istat)
        if(istat/=0) call memerr(two,'bvess_nsa','nsa_membcs',0_ip)
        call memchk(two,istat,mem_modul(1:2,modul),'kfl_funno_nsa','nsa_membcs',kfl_funno_nsa)
        deallocate(kfl_funno_nsa,stat=istat)
        if(istat/=0) call memerr(two,'kfl_funno_nsa','nsa_membcs',0_ip)
        call memchk(two,istat,mem_modul(1:2,modul),'funpa_nsa','nsa_membcs',funpa_nsa)
        deallocate(funpa_nsa,stat=istat)
        if(istat/=0) call memerr(two,'funpa_nsa','nsa_membcs',0_ip)
     else
        call memchk(two,istat,mem_modul(1:2,modul),'bvess_nsa','nsa_membcs',bvess_nsa)
        deallocate(bvess_nsa,stat=istat)
     end if

    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------!
    call deallocate_lodi_bou() 
!!$    call memchk(two,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', xschlrn_nsa)
!!$    deallocate( xschlrn_nsa, stat=istat)
!!$    !
!!$    call memchk(two,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', schlrn_nsa)
!!$    deallocate( schlrn_nsa ,stat=istat)
!!$    !
!!$    call memchk(two,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', xchrc_nsa)
!!$    deallocate( xchrc_nsa, stat=istat)
!!$    ! 
!!$    call memchk(two,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', chrc_nsa)
!!$    deallocate( chrc_nsa, stat=istat)
!!$    !
!!$    call memchk(two,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', gamma_nsa)
!!$    deallocate( gamma_nsa ,stat=istat)
!!$    !
!!$    deallocate( sxinv_nsa, stat=istat)
!!$    deallocate(    sx_nsa, stat=istat)
!!$    !-----------------------------------------------------------------------!

  case (4) 

     allocate(tload_nsa(10),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'TLOAD_NSA','nsa_reabcs',tload_nsa)

  case default
     !
     !Allocate memory for time dependent boundary conditions
     !
     if ((itask > 10_ip).and.(itask < 21_ip)) then 
        ifunc = itask - 10
        allocate(tload_nsa(ifunc)%a(ndofn_nsa+1,mtloa_nsa(ifunc)),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'TLOAD_NSA','nsa_reabcs',tload_nsa(ifunc)%a)
     endif

  end select

end subroutine nsa_membcs



