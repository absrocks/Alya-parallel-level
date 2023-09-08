subroutine chm_outinf()
  !-----------------------------------------------------------------------
  !****f* partis/chm_outinf
  ! NAME 
  !    chm_outinf
  ! DESCRIPTION
  !    This routine writes informtation
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solver
  use def_chemic
  use mod_outfor, only : outfor
  implicit none
  character(60) :: equat
  integer(ip)   :: ierhs,imate,ispec
  character(2)  :: wcpcv

  if( INOTSLAVE ) then
     !
     ! Write information in Result file
     !
     if(kfl_rstar/=2) then
        if (kfl_model_chm == 4) then
           call outfor(25_ip,momod(modul) % lun_outpu,'MODEL: COMBUSTION (TRANSPORT OF SPECIES)')
           write(momod(modul) % lun_outpu,*) 'SPECIES NAMES:'
           write(momod(modul) % lun_outpu,*) '--------------'
           do ispec=1,nspec_chm
              write(momod(modul) % lun_outpu,*) ispec,':',speci(ispec)%name
           enddo
        endif
        flush(momod(modul) % lun_outpu)
     end if

  end if


end subroutine chm_outinf

