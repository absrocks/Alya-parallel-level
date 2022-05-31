subroutine outvax(&
     ivari,jttim,dutim,wopos)
  !-----------------------------------------------------------------------
  !****f* outrut/outvax
  ! NAME
  !   outvax
  ! DESCRIPTION
  !    Output variables for:
  !    - Time step  jttim
  !    - Time value dutim
  ! USES
  !    postpr
  ! USED BY
  !    ***_outvar
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use mod_postpx
  
  implicit none
  
  integer(ip),  intent(in)    :: ivari
  integer(ip),  intent(in)    :: jttim
  real(rp),     intent(in)    :: dutim
  character(5), intent(in)    :: wopos(3)
  integer(ip)                 :: ioscx,iovex
  complex(rp)                 :: duscx(1),duvex(1,1)
  character(5)                :: wauxi(2,6)
  character(2)                :: waux2(6)

  ivapo = ivari
  if( mitim == 0 .and. ittyp == ITASK_ENDRUN ) goto 10

  ioscx = 0
  iovex = 0
  
  print *,kfl_paral,'outvax ',wopos(2)

  if( wopos(2) == 'SCALX' ) then
     !
     ! Postprocess complex scalar
     !
     ioscx = 1

  else if( wopos(2) == 'VECTX' ) then
     !
     ! Postprocess complex vector
     !
     iovex = 1
     
  end if

  if( ioscx == 1 ) then
     !
     ! Postprocess a complex scalar
     !    
     if( IMASTER ) then
      call postpx(duscx,wopos,jttim,dutim)
     else
      call postpx(gescx,wopos,jttim,dutim)
     end if
  end if

  if( iovex == 1 ) then
     !
     ! Postprocess a complex vector
     !    
     if( IMASTER ) then
      call postpx(duvex,wopos,jttim,dutim) 
     else
      call postpx(gevex,wopos,jttim,dutim)
     end if
  end if

10 continue

  ivapo = 0

end subroutine outvax
