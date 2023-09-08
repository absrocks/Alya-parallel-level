subroutine ada_reapro
!-----------------------------------------------------------------------
!****f* adapti/Adapti
! NAME 
!    ada_reapro
! DESCRIPTION
!    This routine reads adapti parameters
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_adapti
  use def_inpout
  use mod_ecoute, only :  ecoute
  implicit none
  
  !
  ! Reach the section
  !      
  rewind(lisda)
  do while(words(1)/='RUNDA')
     call ecoute('ADA_REAPRO')
  end do
  do while(words(1)/='ENDRU')
     call ecoute('ADA_REAPRO')
  end do
  do while(words(1)/='PROBL')
     call ecoute('ADA_REAPRO')
  end do
  !
  ! Read data
  !      
  do while(words(1)/='ENDPR')
     call ecoute('ADA_REAPRO')
     if(words(1)==naser(servi)(1:5)) then 
        if(exists('ON   ')) then
           kfl_livei_ada   = 0
           kfl_servi(servi)= 1
           do while(words(1)/='END'//naser(servi)(1:2))
              call ecoute('ADA_REAPRO')
              if (words(1) .eq.'VERBO') kfl_livei_ada = 1
           end do
        end if
     end if
  end do

end subroutine ada_reapro
