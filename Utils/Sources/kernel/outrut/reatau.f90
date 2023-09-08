subroutine reatau(kfl_taust)
  !-----------------------------------------------------------------------
  !****f* outrut/reatau
  ! NAME 
  !    reatau
  ! DESCRIPTION
  !    This routine reads the tau strategy use to compute tau in tauadr
  ! INPUT
  ! OUTPUT
  ! USES
  ! USED BY
  !    ***_reanut
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip
  use def_inpout, only      :  words
  implicit none
  integer(ip),  intent(out) :: kfl_taust

  if(words(2)=='OFF  ') then
     kfl_taust=0 
  else if(words(2)=='CODIN') then
     kfl_taust=1
  else if(words(2)=='EXACT') then
     kfl_taust=2
  else if(words(2)=='SHAKI') then
     kfl_taust=3
  else if(words(2)=='INCLU') then
     kfl_taust=5
  else if(words(2)=='TIMES') then
     kfl_taust=6
  end if

end subroutine reatau
