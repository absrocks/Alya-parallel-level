subroutine ada_modmsh(itask)
!-----------------------------------------------------------------------
!****f* adapti/ada_modmsh
! NAME 
!    ada_modmsh
! DESCRIPTION
!    This routine performs mesh modification tasks
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_adapti
  implicit none
  integer(ip), intent(in) :: itask


  if (nimmo_ada == 0) then
     !
     ! No immersed objects     
     !     
     if (kfl_remst_ada > 10) then
        call ada_coamsh(itask)
        call ada_divmsh(itask)
        call ada_newcon(itask)
     end if

  else
     !
     ! Immersed objects are present
     !
     ! Use elsest to identify the immersed object situation relative to
     ! the background mesh
     !
     call ada_elsest
     !
     ! Segmentate the background mesh
     ! 
     call ada_xtract
     call ada_sahole
     call ada_saptch

  end if
  
end subroutine ada_modmsh
