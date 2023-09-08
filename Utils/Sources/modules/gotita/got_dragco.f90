subroutine got_dragco(ReD,CDReD)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_dragco
  ! NAME 
  !    got_dragco
  ! DESCRIPTION
  !    Compute the drag coefficient
  ! USES
  !
  ! USED BY
  !    got_elmrpo
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only  :  ip,rp
  use def_gotita, only  :  kfl_exacs_got,kfact_got
  implicit none
  real(rp), intent(in)  :: ReD
  real(rp), intent(out) :: CDReD
  real(rp)              :: w,CD
  integer(ip)           :: kfl_dragm_got

  if(kfl_exacs_got==0) then

     kfl_dragm_got=1

     select case(kfl_dragm_got)

     case(1)
        !
        ! Water - default
        !
        if(ReD<=1300.0_rp) then
           CDReD=24.0_rp*(1.0_rp+0.15_rp*ReD**0.687_rp)
        else
           CDReD=0.4_rp*ReD
        end if

     case(2)
        !
        ! Water  - extended Reynolds
        !
        w = log10(ReD)
        if(ReD<=0.01_rp) then
           CD = 3.16_rp+24.0_rp*ReD
        else if(ReD<=20.0_rp) then
           CD = 24.0_rp/ReD*(1.0_rp+0.1315_rp*ReD**(0.82_rp-0.05_rp*w))
        else if(ReD<=260.0_rp) then
           CD = 24.0_rp/ReD*(1.0_rp+0.1935_rp*ReD**0.6305_rp)
        else if(ReD<=1500.0_rp) then
           CD = 10.0_rp**(1.6435_rp-1.1242_rp*w+0.1558_rp*w**2.0_rp)
        else if(ReD<=1.2e4) then
           CD = 10.0_rp**(-2.4571_rp+2.5558_rp*w-0.9295_rp*w**2.0_rp+0.1049_rp*w**3.0_rp)
        else if(ReD<=4.4e4) then
           CD = 10.0_rp**(-1.9181_rp+0.6370_rp*w-0.0636_rp*w**2.0_rp)
        else if(ReD<=3.38e5) then
           CD = 10.0_rp**(-4.3390_rp+1.5809_rp*w-0.1546_rp*w**2.0_rp)
        else if(ReD<=4.0e5) then
           CD = 29.78_rp-5.3_rp*w
        else if(ReD<=1.0e6) then
           CD = -0.49_rp+0.1_rp*w
        else
           CD = 0.19_rp-80000.0_rp/ReD
        end if
        CDReD=CD*ReD

     end select

  else
     !
     ! Exact solution: make sig=1
     !
     CDReD=24.0_rp*kfact_got

  end if

end subroutine got_dragco
