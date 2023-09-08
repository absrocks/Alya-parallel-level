subroutine exm_membcs(itask)
!-----------------------------------------------------------------------
!****f* Exmedi/exm_membcs
! NAME 
!    exm_membcs
! DESCRIPTION
!    This routine allocates memory for the boundary conditions
! USES
!    ecoute
! USED BY
!    exm_turnon
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk

  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: nfipo,ndobo
  integer(4)              :: istat
     
  if (itask <= 1) then
     ndobo= 1
     npoin= 1
     if (itask == 1) then
        ndobo = nboun
        nfipo = npoin
     end if
     !
     ! Allocate memory
     !
     allocate(kfl_fixno_exm(ndofn_exm,nfipo),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_EXM','exm_membcs',kfl_fixno_exm)
     allocate(kfl_fixbo_exm(ndobo),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXBO_EXM','exm_membcs',kfl_fixbo_exm)
     allocate(bvess_exm(ndofn_exm,nfipo),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_EXM',    'exm_membcs',bvess_exm)
     allocate(bvnat_exm(ndofn_exm,mnodb,ndobo),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'BVNAT_EXM',    'exm_membcs',bvnat_exm)

  else if (itask==2) then
     

     ! to be done...



  end if

end subroutine exm_membcs



