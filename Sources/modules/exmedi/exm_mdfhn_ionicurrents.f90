  !------------------------------------------------------------------------
  !> @addtogroup ExmediIonicurrents
  !> @{
  !> @file    exm_mdfhn_ionicurrents.f90
  !> @date    17/10/2016
  !> @author  Jazmin A-S
  !> @brief   Computes the modified Fitzhugh-Nagumo EP model
  !> @details Implementation from Rogers & McCulloch 1994
  !> @}
  !
  ! Model solved in this subroutine:
  ! (see that the last term contains u, this is a modification suggested in R&M 1994)
  !                   Iion= C_1*u(u-A)(1-u)-C_2*v*u 
  !                   dv/dt=  B*(u-Dv)   
  ! 
  ! This subroutine computes Iion, which can be then added to the original equation
  ! 
  !                   du/dt= diffusion terms + (Istim+Iion)/Cm
  ! Where:
  !       u           => voltage (elmag(:,ITER_K))
  !       v           => recovery variable (refhn_exm)
  !       afhn=0.13_rp  => original values of paper: a=0.13;    old alpha
  !       bfhn=0.0035_rp => original values of paper: b=0.013;  old epsilon
  !       c1fhn=0.26_rp  => original values of paper: c1=0.26;  old cmemb
  !       c2fhn=0.1_rp   => original values of paper: c2=0.1;   
  !       dfhn=1.0_rp    => original values of paper: d=1.0;    old gamma
  !------------------------------------------------------------------------
subroutine exm_mdfhn_ionicurrents(ipoin,xioni,dioni)
  use      def_master
  use      def_domain
  use      def_elmtyp
  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: ipoin
  real(rp), intent(out)   :: xioni
  real(rp), intent(out)   :: dioni
  integer(ip)             :: imate 

  real(rp)   :: afhn, & ! A
       bfhn, & ! B
       c1fhn, & ! C1
       c2fhn, & ! C2
       dfhn, dv, k1, k2, k3, k4, dtfhn,  & !  D
       voltiter,refiter,reftime,refnew, &
       aauxi, bauxi

  imate= nodemat(ipoin)

  afhn  = xmopa_exm( 1,imate)           
  bfhn  = xmopa_exm( 2,imate)
  c1fhn = xmopa_exm( 3,imate)           
  c2fhn = xmopa_exm( 4,imate)           
  dfhn  = xmopa_exm( 5,imate)                           

  !
  !refhn_exm(ipoin,1) == new v
  !refhn_exm(ipoin,2) == used v
  !refhn_exn(ipoin,3) == old v

  ! voltage scale conversion:
  ! phi mV = (120_rp*x)-85.0_rp = (poref_fhn_exm(2)-poref_fhn_exm(1)) * x + poref_fhn_exm(1) = aauxi * x + bauxi
  ! x      = (phi - bauxi) / aauxi 
  ! time scale conversion:
  ! dtfhn  = dt * 1000

  dtfhn=  dtime * sms_conversion_currents_exm 

  bauxi = poref_fhn_exm(1)
  aauxi = poref_fhn_exm(2) - poref_fhn_exm(1)

  voltiter= (elmag(ipoin,ITER_K  ) - bauxi) / aauxi

  refiter = refhn_exm(ipoin,ITER_K)
  reftime = refhn_exm(ipoin,TIME_N)

  ! solve the refhn time equation ! dv/dt = b*(u - d*v)

  dv=bfhn*(voltiter - dfhn*reftime)   

  k1 = dv
  k2 = dv + 0.5_rp * dtfhn* k1
  k3 = dv + 0.5_rp * dtfhn * k2
  k4 = dv + dtfhn * k3

  ! update refhn_exm (NOT DONE IN UPDUNK)
  refnew = reftime + (dtfhn / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
  refhn_exm(ipoin,ITER_K) = refnew

  ! compute the new xioni: using REFITER instead of REFNEW, because it is more robust
  ! con la correccion de mcculloch: C_1*u(u-A)(1-u)-C_2*v*u 

  ! dioni is the ionic current derivative w.r.t. phi
  ! xioni is negative because it goes to the RHS

  !!    dioni = c1fhn*(voltiter-afhn)*(voltiter-1.0_rp)-c2fhn*refiter  ! antes
  !!    xioni = stimu + voltiter * dioni  ! antes

  ! o, como estaba antes, sin la correccion de mcculloch: C_1*u(u-A)(1-u)-C_2*v
  !    dioni = c1fhn*(voltiter-afhn)*(voltiter-1.0_rp)

  dioni = (c1fhn*(voltiter-afhn)*(1.0_rp - voltiter)-c2fhn*refnew)
  xioni = voltiter * dioni 

  xioni = -aauxi * xioni * sms_conversion_currents_exm
  dioni = aauxi * dioni

end subroutine exm_mdfhn_ionicurrents

