subroutine exm_turnof
!-----------------------------------------------------------------------
!****f* Exmedi/exm_turnof
! NAME 
!    exm_turnof
! DESCRIPTION
!    This routine closes the run for the current module
! USES
!    exm_outcpu
!    exm_output
! USED BY
!    Exmedi
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_exmedi
  use      def_postpr, only: kfl_ivari
  implicit none

  !if automatic isochrone saving is on, save the last set of isochrones
  if (thiso_exm(3)>0.0_rp) then
    kfl_ivari(1) = 26_ip
    call exm_outvar(26_ip,0_ip)
  end if

  !this array is now allocated/deallocated in mod_exm_arrays
  !if(INOTMASTER) then
  !  call memory_deallo(mem_modul(1:2,modul),'ISOCH_MODIFIED', 'exm_turnof', isoch_modified)
  !end if


end subroutine exm_turnof

