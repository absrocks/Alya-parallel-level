subroutine chm_flame_prop(gpfar,gpspe,gpthi) 
  !-----------------------------------------------------------------------
  ! ****f* partis/chm_flame_prop
  ! NAME 
  !    chm_flame_speed
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

  real(rp),     intent(in)  :: gpfar
  real(rp),     intent(out) :: gpspe
  real(rp),     intent(out) :: gpthi 

  real(rp),   dimension(10) :: flamsp  ! Flame speed database
  real(rp),   dimension(10) :: flamth  ! Flame thickness database
  real(rp),   dimension(10) :: flameq  ! Equivalence ratio database
  real(rp)                  :: mfspe
  real(rp)                  :: mfthi
  real(rp)                  :: gpeqr
  real(rp)                  :: far_st  ! Fuel / air stoichoimetric ratio
  integer(ip)               :: ival1,ival2,idata

  !
  ! Tabulated values of flame speed and flame thickness for methane
  ! C. K. Law B. H. Chao, F. N. Egolfopoulos, Structure and propagation of premixed flame in nozzlegenerated
  !                                           counterflow, Combustion and Flame 109 (1997), 620–638. 
  flameq = (/ 0.5_rp, 0.6_rp, 0.7_rp, 0.8_rp, 0.9_rp, 1.0_rp, 1.1_rp, 1.2_rp, 1.3_rp, 1.4_rp/)
  flamsp = (/ 0.08_rp, 0.12_rp, 0.2_rp, 0.25_rp, 0.33_rp, 0.365_rp, 0.37_rp, 0.33_rp, 0.21_rp, 0.08_rp /)
  flamth = (/ 0.0015_rp, 0.0013_rp, 0.001_rp, 0.0007_rp, 0.0005_rp, 0.0004_rp, 0.0006_rp, 0.00085_rp, 0.0009_rp, 0.0011_rp /)
  far_st = 17.16_rp                           ! Stoichoimetric fuel/air ratio for methane

  ! 
  ! Conversion to equivalence ratio ER = (F/A)real / (F/A)stoi
  ! 
  gpeqr = gpfar * far_st

  !
  ! Interpolation
  !
  ival1 = 0
  ival2 = 0

  if (gpeqr >= flameq(1) .and. gpeqr <= flameq(10)) then
     do idata =1,size(flameq)-1_ip
        if (gpeqr >= flameq(idata) .and. gpeqr <= flameq(idata+1_ip) ) then
           ival1 = idata
           ival2 = idata + 1_ip
        end if
     end do 

     mfspe = (flamsp(ival2) - flamsp(ival1) ) / (flameq(ival2) - flameq(ival1) )
     mfthi = (flamth(ival2) - flamth(ival1) ) / (flameq(ival2) - flameq(ival1) )
 
     gpspe = flamsp(ival1) + mfspe * (gpeqr - flameq(ival1)) 
     gpthi = flamth(ival1) + mfthi * (gpeqr - flameq(ival1))

  else
     gpspe = 0.0_rp
     gpthi = 0.0_rp
  endif
  
end subroutine chm_flame_prop
