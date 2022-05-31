subroutine nsa_memphy
!-----------------------------------------------------------------------
!****f* Nastal/nsa_memphy
! NAME 
!    nsa_memphy
! DESCRIPTION
!    This routine allocates memory for physical arrays
! USES
!    ecoute
!    memchk
!    runend
! USED BY
!    nsa_reaphy
!***
!-----------------------------------------------------------------------
  use      def_parame 
  use      def_inpout
  use      def_master
  use      def_nastal
  use      def_domain
  use      mod_memchk
  implicit none
  integer(4) :: istat

  if(nmate_nsa>1) then
     allocate(lmate_nsa(nelem),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LMATE_NSA','nsa_memphy',lmate_nsa)
     allocate(lmatn_nsa(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LMATN_NSA','nsa_memphy',lmatn_nsa)
     lmate_nsa=1 ! Fluid by default
  else
     allocate(lmatn_nsa(1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LMATN_NSA','nsa_memphy',lmatn_nsa)
  end if

  if (kfl_infun_nsa == 3 .and. kfl_brunt_nsa == 2) then   ! brunt-vaisala frequency field
     allocate(brunt_nsa(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'BRUNT_NSA','nsa_memphy',brunt_nsa)     
  else
     allocate(brunt_nsa(    1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'BRUNT_NSA','nsa_memphy',brunt_nsa)
  end if

  
end subroutine nsa_memphy
