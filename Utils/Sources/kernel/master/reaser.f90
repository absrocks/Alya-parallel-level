subroutine reaser()
  !------------------------------------------------------------------------
  !****f* Dodeme/reaser
  ! NAME
  !    reaser
  ! DESCRIPTION
  !    This routine 
  ! OUTPUT
  ! USED BY
  !    Dodeme
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_inpout
  use mod_ecoute, only :  ecoute
  implicit none
  character(10) :: messa

  if( INOTSLAVE ) then
     !
     ! Reach the section
     !      
     lisda=lun_pdata
     messa=exser(servi)//'_REAPRO'
     rewind(lisda)
     do while(words(1)/='RUNDA')
        call ecoute(messa)
     end do
     do while(words(1)/='ENDRU')
        call ecoute(messa)
     end do
     do while(words(1)/='PROBL')
        call ecoute(messa)
     end do
     !
     ! Read data
     !      
     do while(words(1)/='ENDPR')
        call ecoute(messa)
        if(words(1)==naser(servi)(1:5)) then 
           if(exists('ON   ')) then
              kfl_servi(servi)=1
              do while(words(1)/='END'//naser(servi)(1:2))
                 call ecoute(messa)
              end do
           end if
        end if
     end do
     
  end if
  
end subroutine reaser

