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
  use      mod_memory

  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: itask
     
  if (itask == 1) then
     !
     ! Allocate memory
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_EXM','exm_membcs',kfl_fixno_exm,ndofn_exm,npoin)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_EXM','exm_membcs',kfl_fixbo_exm,nboun)
     call memory_alloca(mem_modul(1:2,modul),'BVESS_EXM'    ,'exm_membcs',bvess_exm,ndofn_exm,npoin)
     call memory_alloca(mem_modul(1:2,modul),'BVNAT_EXM'    ,'exm_membcs',bvnat_exm,ndofn_exm,mnodb,nboun)
  
  else if (itask==2) then
     

     ! to be done...



  end if

end subroutine exm_membcs



