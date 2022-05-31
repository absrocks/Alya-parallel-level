subroutine opealy(itask,wopos)

  !-----------------------------------------------------------------------
  !    
  !
  ! open case file and write some data on it
  !
  !
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use mod_iofile
  use def_postpr
  implicit none
  integer(ip),  intent(in) :: itask
  character(5), intent(in) :: wopos(2)
  character(150)           :: filva
  
  select case(itask)

  case(1)
     !
     ! Open file 
     !  
     if( ncoun_pos == 0 ) call wrialy(2_ip,'NULL')
     call wrialy(3_ip,wopos)
     ncoun_pos = ncoun_pos + 1
     
     filva = trim(fil_postp)//'.'//trim(wopos(1))//'-'//adjustl(trim(nunam_pos))//'.alya'
     call iofile(zero,lun_postp,filva,'ALYA ASCII'//trim(wopos(1))//' RESULTS FILE')

  case(2)
     !
     ! Close file
     !
     call iofile(two,lun_postp,'NULL','ALYS ASCII'//trim(wopos(1))//' RESULTS FILE')

  end select


end subroutine opealy
