subroutine chm_updtsp(dtmin)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updtsp
  ! NAME 
  !    chm_updtsp
  ! DESCRIPTION
  !    This routine computes next timestep based on the accuracy reached
  !    with previous timestep
  ! USED BY
  !    chm_timste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  real(rp),   intent(inout) :: dtmin
  integer(ip)               :: ipoin,ispec,ii,ipoi1,ipoi2
  real(rp)                  :: numer,denom,kp,ki,kd,tol
  real(rp),    target       :: dtpar(1)
  real(rp),    save         :: e(3)

  if( kfl_paral >= 0 ) call runend('TO DO')

  kp  = 0.075_rp
  ki  = 0.175_rp
  kd  = 0.01_rp
  tol = 2.0_rp

  if( ittim == 0 ) e = 0.0_rp

  numer = 0.0_rp
  denom = 0.0_rp
  do ispec = 1,nspec_chm
     do ii = 1,2
        if( ii == 1 ) then
           ipoi1 = 1
           ipoi2 = npoi1
        else
           ipoi1 = npoi2
           ipoi2 = npoi3                 
        end if
        do ipoin = ipoi1,ipoi2
           numer = numer + ( conce(ipoin,ispec,3) - conce(ipoin,ispec,4) ) ** 2
           denom = denom + conce(ipoin,ispec,3)**2
        end do
     end do
  end do

  e(3) = e(2)
  e(2) = e(1)
  e(1) = sqrt(numer/denom)

  if( ittim <= 4 ) then

     dtmin = dtold(1)

  else

     dtmin = ( e(2) / e(1) )**kp * ( tol / e(1) )**ki * ( e(2) * e(2) / ( e(1) * e(3) ) )**kd
     dtmin = dtmin * dtold(1)

  end if
  
end subroutine chm_updtsp
