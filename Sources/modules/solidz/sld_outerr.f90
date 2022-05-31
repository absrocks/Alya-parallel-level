subroutine sld_outerr()
  !------------------------------------------------------------------------
  !****f* Solidz/sld_outerr
  ! NAME 
  !    sld_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    sld_turnon
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_kermod
  use def_domain
  use def_solidz
  use mod_outfor, only : outfor
  implicit none
  integer(ip)    :: ierro=0,iwarn=0,jerro
  character(200) :: wmess

  iwarn = 0
  ierro = 0

  if( nmate == 0 ) then
     ierro = ierro+1
     wmess = 'MATERIALS MUST BE DECLARED EXPLICITLY IN THE DOMAIN (*.DOM.DAT FILE)'
     call outfor(1_ip,momod(modul)%lun_outpu,trim(wmess))        
  end if

  if( kfl_cutel == 0 .and. kfl_xfeme_sld == 1 ) then
     ierro = ierro+1
     wmess = 'PUT CUT ELEMENTS ON IN THE KERNAL DATA FILE WHEN USING X-FEM'
     call outfor(1_ip,momod(modul)%lun_outpu,trim(wmess))        
  end if

  if( kfl_vofor_sld > 0 ) then
     jerro = 0
     if( INOTMASTER ) then
        if( size(xfiel) < kfl_vofor_sld ) then
           jerro = 1
        else if( .not. associated(xfiel(kfl_vofor_sld) % a)) then
           jerro = 1
        end if
     end if
     call parari('MAX',0_ip,1_ip,jerro)
     if( jerro == 1 ) then
        ierro = ierro + 1
        wmess = 'EXTERNAL VOLUME FORCE FIELD DOES NOT EXIST'
        call outfor(1_ip,momod(modul)%lun_outpu,trim(wmess))        
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------

  call errors(3_ip,ierro,iwarn,' ')

end subroutine sld_outerr
