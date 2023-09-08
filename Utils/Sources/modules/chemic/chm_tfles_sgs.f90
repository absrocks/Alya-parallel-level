subroutine chm_tfles_sgs(gpspe,gpthi,gptur,gpvol,sgsef,gpfac,omega) 
  !-----------------------------------------------------------------------
  ! ****f* partis/chm_tfles_sgs
  ! NAME 
  !    chm_tfles_sgs
  ! DESCRIPTION
  !    Interpolate the laminar flame speed and flame thickness for a given fuel
  !    from the equivalance ratio.
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_chemic
  use def_master 
  use mod_ker_proper 
  use def_kermod
  implicit none

  real(rp),     intent(in)  :: gpspe                 ! Flame speed
  real(rp),     intent(in)  :: gpthi                 ! Flame thickness
  real(rp),     intent(in)  :: gptur                 ! Turbulent viscosity
  real(rp),     intent(in)  :: gpvol                 ! Element volume at gauss points
  real(rp),     intent(in)  :: omega                 ! Dynamic sensor for DTFLES
  real(rp),    intent(out)  :: sgsef                 ! Efficiency function for the wrinkling of the flame front
  real(rp),    intent(out)  :: gpfac                 ! Dynamic thickening factor at gauss point

  real(rp)                  :: filter                 ! Filter size for TFLES (different from LES)
  real(rp)                  :: const                  ! Model constant Cs to extract the subgrid scale velocity from mut
  real(rp)                  :: alpha                  ! Model constant alpha
  real(rp)                  :: uprime                 ! Subgrid scale velocity
  real(rp)                  :: wratio_laminar_flame   ! wratio = <a>_s / (u'/ h) laminar flame
  real(rp)                  :: wratio_thicken_flame   ! wratio = <a>_s / (u'/ h) thickened flame
  real(rp)                  :: wrink_lam              ! Wrinkling of the sgs flame front (laminar)
  real(rp)                  :: wrink_thi              ! Wrinkling of the sgs flame front (thickened)
  
  !
  ! Computation TFLES filter size
  !
  filter = 5.0_rp * (gpvol) ** 0.333333_rp

  !
  ! Computation subgrid scale velocity
  !
  const  = 0.25_rp
  if (filter == 0.0_rp) then
     uprime = 0.0_rp
  else
     uprime = gptur / (const * filter ) 
  endif

  alpha = 0.19_rp
  sgsef = 1.0_rp

  if (gpspe > 0.0_rp .and. uprime > 0.0_rp) then
     !
     ! Ratio of strain rate over turbulent fluctuation wratio = <a>_s / (u'/ h)
     !
     wratio_laminar_flame = 0.75_rp * exp( -1.2_rp / (uprime / gpspe)**(0.3_rp) ) * (filter / gpthi )**0.6666667_rp
     wratio_thicken_flame = 0.75_rp * exp( -1.2_rp / (uprime / gpspe)**(0.3_rp) ) * (filter / gpthi / tfles_chm )**0.6666667_rp  
     !
     ! Wrinkling of the sgs flame front
     !
     wrink_lam = 1.0_rp + alpha * wratio_laminar_flame * (uprime / gpspe)
     wrink_thi = 1.0_rp + alpha * wratio_thicken_flame * (uprime / gpspe)     
     !
     ! Efficiency function
     !     
     sgsef = wrink_lam / wrink_thi

  endif
 
  sgsef = max(1.0_rp,sgsef)
  sgsef = min(tfles_chm**(0.666666_rp),sgsef)
  !
  ! Dynamic thickening factor F = 1 + (F_max - 1)*OMEGA
  !
  gpfac = 1.0_rp + (tfles_chm - 1.0_rp)*omega

end subroutine chm_tfles_sgs
