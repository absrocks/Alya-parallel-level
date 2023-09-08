subroutine exm_memphy(itask)
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_memphy
  ! NAME 
  !    exm_memcbs	
  ! DESCRIPTION	
  !    This routine allocates memory for the nodal/elementary properties
  ! USES	
  !    ecoute
  ! USED BY
  !    exm_reaphy
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  
  use      def_exmedi
	
  implicit none
  integer(4)              :: istat
  integer(ip)             :: nsize,itask,nsiel
  
  nsize = npoin	                            ! nodal-wise
  nsiel = nelem
  
  select case(itask)
     
  case( 1_ip)
     !
     ! Conductivity:  nodal-CEDIF 
     !
     allocate(cedif_exm(ndime,ndime,nsize),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CEDIF_EXM','exm_memphy',cedif_exm)
     allocate(grafi_exm(ndime,ndime,nsize),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'GRAFI_EXM','exm_memphy',grafi_exm)
     allocate(fiber_exm(ndime,ndime),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'FIBER_EXM','exm_memphy',fiber_exm)


  end select
  
end subroutine exm_memphy


