subroutine qua_exacso(&
     itask,gpcod,gpden,gpdif,gprea,gpgrd,&
     gpvel,gpsou,extem,exteg,gprhs)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_exacso
  ! NAME 
  !    qua_exacso
  ! DESCRIPTION
  !    This routine computes the exact solution at the integration points
  !    ITASK = 1 ... Computes phion in and gradients 
  !    ITASK = 2 ... Computes force vector in to be applied as a source
  ! USES
  ! USED BY
  !    qua_elmope
  !    qua_exacso
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master, only       :  cutim
  use def_domain, only       :  ndime,kfl_naxis
  use def_quanty, only       :  kfl_exacs_qua,kfl_timei_qua !,expar_qua
  implicit none
  integer(ip), intent(in)    :: itask
  real(rp),    intent(in)    :: gpden,gpdif,gprea,gpgrd(ndime)
  real(rp),    intent(in)    :: gpsou,gpvel(ndime),gpcod(ndime)
  real(rp),    intent(out)   :: extem,exteg(ndime)
  real(rp),    intent(inout) :: gprhs
  real(rp)                   :: x,y,z,t,u,v,r,k,a,dtdx,dtdy,zequa
  real(rp)                   :: d2tdy2,d2tdx2,dtdt,Q,dkdx,dkdy
  real(rp)                   :: rmayo,rmeno
  !
  ! Initializations
  ! 
  if(itask==2) then

  end if

  if (kfl_exacs_qua==1)then

     !
     ! one exact solution! 
     !
     
  else if (kfl_exacs_qua==2)then

     !
     ! other exact solution !
     !
     ! and go on! 
  end if

  if(itask==1) then
     !
     ! Exact solution and gradients
     !

  else if(itask==2) then

  end if

end subroutine qua_exacso
