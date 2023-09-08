subroutine ada_livinf(message)
    
!-----------------------------------------------------------------------
!****f* adapti/ada_livinf
! NAME
!    ada_livinf
! DESCRIPTION
!    This routine write particular live information for this module/service
! USES
! USED BY
!    Reapro
!***
!-----------------------------------------------------------------------
    
  use      def_master
  use      def_adapti
  implicit none
  character(*)              :: message
  character(300)            :: messa,dumml

  if(kfl_paral<=0.and.kfl_livei_ada/=0.and.lun_livei/=0) then
     
     messa='--|       {  Adapti : '//adjustl(trim(message))//' } '
     write(6,'(a)') adjustl(trim(messa))

  end if

end subroutine ada_livinf
