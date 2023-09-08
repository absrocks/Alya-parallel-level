 subroutine chm_updsou()
  !------------------------------------------------------------------------
  !****f* chemic/chm_updsou
  ! NAME 
  !    chm_updsou
  ! DESCRIPTION
  !    This routine performs the gather operations
  ! USES
  ! USED BY
  !    chm_begite
  !***
  !------------------------------------------------------------------------ 
  use def_kintyp
  use def_domain
  use def_master
  use def_chemic
  use def_kermod, only     : gasco
  implicit none 
  integer(ip) :: ipoin
  real(rp)    :: cpl,cpml,cpv,cvv,Lv0,T0,R,Rv,eps
  real(rp)    :: T,p,Lv,rrv,rrc,rrl,Rm,cvml,dt,es
  real(rp)    :: rvs,cp,cv,denom
 
  if( kfl_sourc_chm == -4 ) then

     !-------------------------------------------------------------------
     !
     ! Moisture model
     !
     !                ( rv - rvs )
     ! Q =  -----------------------------------
     !      dt * ( 1 + Lv^2 rvs / ( Cp Rv T^2 )
     !
     !-------------------------------------------------------------------

     cp   = cpcoe_chm
     cv   = cvcoe_chm
     cpl  = cpliq_chm  
     cpml = cpmli_chm  
     cpv  = cpvap_chm
     cvv  = cvvap_chm
     Lv0  = lhref_chm
     T0   = teref_chm
     R    = gasco
     Rv   = rgava_chm
     eps  = R / Rv
     dt   = dtime

     do ipoin = 1,npoin

        T     = tempe(ipoin,1)
        p     = press(ipoin,1)
        Lv    = Lv0 - ( cpl - cpv ) * ( T - T0 ) 
        rrv   = conce(ipoin,1,1)
        rrc   = conce(ipoin,2,1)
        rrl   = rrc
        es    = 611.2_rp * exp( 17.67 * ( T - 273.15_rp ) / ( T - 29.65_rp ) )
        rvs   = eps * es / ( p - es )
        Rm    = R  + Rv * rrv
        cvml  = cv + cvv * rrv + cpl * rrl           
        denom = 1.0_rp / ( dt * ( 1.0_rp + Lv * Lv * rvs / ( cp * Rv * T * T ) ) )
        !
        ! ICLAS = 1: rv
        !
        treac_chm(1,ipoin) =  1.0_rp * denom       
        tmrat_chm(1,ipoin) =  rvs    * denom
        !
        ! ICLAS = 2: rc
        !
        treac_chm(2,ipoin) = -1.0_rp * denom       
        tmrat_chm(2,ipoin) = -rvs    * denom

     end do

  end if
  
end subroutine chm_updsou
