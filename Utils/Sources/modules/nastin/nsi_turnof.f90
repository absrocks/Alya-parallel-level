subroutine nsi_turnof
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_turnof
  ! NAME 
  !    nsi_turnof
  ! DESCRIPTION
  !    This routine closes NASTIN module
  ! USES
  !    nsi_output
  !    nsi_outlat
  !    nsi_openfi
  ! USED BY
  !    Nastin
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_nastin
  implicit none

  !-------------------------------------------------------------------
  !
  ! Interpolation from coarse to fine mesh
  !
  !-------------------------------------------------------------------
  !
  if(kfl_meshi_nsi /= 0_ip) call nsi_coarfine(2_ip)

  !
  ! Output CPU times
  !
  call nsi_outcpu()
  !
  ! Output latex file
  !
  call nsi_outlat(two)

!  if( INOTMASTER ) then
!  call memgen(0_ip,npoin,0_ip)
!  do ipoin = 1,npoin
!     gesca(ipoin)=real(lnlev(ipoin))
!  end do
!  call possla(2_ip,parin,gesca,parre)
!end if
!call runend('OESOS')

   call nsi_cadan(7_ip) !llamada a Alya ADAN para cierre de archivos

end subroutine nsi_turnof

