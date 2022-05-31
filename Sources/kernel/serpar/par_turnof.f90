subroutine par_turnof
  !-----------------------------------------------------------------------
  !****f* parall/par_turnof
  ! NAME
  !    par_turnof
  ! DESCRIPTION
  !    This subroutine turns off service
  ! USES
  ! USED BY
  !    Turnof
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_parall
  use mod_parall
  use mod_outfor, only : outfor
  implicit none
  !
  ! Write CPU time heading and master's CPU time
  !
  call par_outcpu()
  !
  ! Write tail for formatted files
  !
  if(IMASTER) then
     call outfor( 26_ip,lun_outpu_par,' ')
  else if(ISLAVE.and.kfl_outpu_par==1.and.kfl_outpu==1) then
     call outfor(-26_ip,lun_outpu,' ')
  end if
  !
  ! Close files
  !
  call par_openfi(four)

end subroutine par_turnof
