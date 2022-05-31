subroutine realen(kfl_ellen)
  !-----------------------------------------------------------------------
  !****f* outrut/realen
  ! NAME 
  !    realen
  ! DESCRIPTION
  !    This routine reads thge natural length calculation
  ! INPUT
  ! OUTPUT
  ! USES
  !    ecoute
  ! USED BY
  !    ***_reanut
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip
  use def_inpout, only      :  words
  implicit none
  integer(ip),  intent(out) :: kfl_ellen

  if(words(2)=='MAXIM') then
     !
     ! Maximum length: hmax
     !
     kfl_ellen=1

  else if(words(2)=='MINIM') then
     !
     ! Minimum length: hmin
     !
     kfl_ellen=0

  else if(words(2)=='AVERA') then
     !
     ! Average length: 0.5*(hmin+hmax)
     !
     kfl_ellen=2

  else if(words(2)=='FLOWD') then
     !
     ! Length in flow direction
     !
     kfl_ellen=3

  else if(words(2)=='APPRO') then
     !
     ! Approximate parameter: sqrt(hmin*hmax)
     !
     kfl_ellen=4

  else if(words(2)=='REALF') then
     !
     ! Length in real flow direction
     !
     kfl_ellen=5
   else if(words(2)=='MIXLE') then
     !
     ! Mixed element length : hmin for tau1 and hmax for tau2
     !	     
     kfl_ellen=6

  end if

end subroutine realen
