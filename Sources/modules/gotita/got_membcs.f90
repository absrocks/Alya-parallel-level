subroutine got_membcs(itask)
!-----------------------------------------------------------------------
!****f* Gotita/got_membcs
! NAME 
!    got_memcbs
! DESCRIPTION
!    This routine allocates memory for the boundary conditions arrays
! USES
!    ecoute
! USED BY
!    got_reabcs
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_gotita
  use      def_domain
  use      mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case (itask)

  case(1)
     !
     ! Allocate memory
     !
     call runend('GOT_MEMBCS: HAY QUE CAMBIAR FIXNO PARA QUE TENGA UNA PRIMERA DIMENSION')
     allocate(kfl_fixno_got(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_FIXNO_GOT','got_membcs',kfl_fixno_got)
     if(kfl_probl_got==1) then
        allocate(bvess_got(ndofn_got(3),npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_GOT','got_membcs',bvess_got)
     else if(kfl_probl_got==2) then
        allocate(bvess_got(ndime,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_GOT','got_membcs',bvess_got)
     else if(kfl_probl_got==3) then
        allocate(bvess_got(1,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_GOT','got_membcs',bvess_got)
     end if

  case(2)
     !
     ! Deallocate memory
     !
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXNO_GOT','got_membcs',kfl_fixno_got)
     deallocate(kfl_fixno_got,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXNO_GOT','nsi_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'BVESS_GOT','got_membcs',bvess_got)
     deallocate(bvess_got,stat=istat)
     if(istat/=0) call memerr(two,'BVESS_GOT','nsi_membcs',0_ip)

  end select

end subroutine got_membcs

