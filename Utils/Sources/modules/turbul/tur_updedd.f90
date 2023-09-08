subroutine tur_updedd()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_updedd
  ! NAME 
  !    tur_updedd
  ! DESCRIPTION
  !    This routine updates the eddy viscosity
  ! OUTPUT
  !    TURNU: Turbulent viscosity mut
  ! USED BY
  !    tur_updunk
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use def_kermod , only     :  kfl_prope, turmu_ker, kfl_logva,kfl_adj_prob
  use mod_ker_proper
  implicit none
  integer(ip) :: ipoin,dummi
  real(rp)    :: dumm1,dumm2,nut,nu,rho(1),mu(1)

  if (kfl_adj_prob == 1) return
  if( INOTMASTER.and.turmu_ker%kfl_exist==0 ) then ! and turmu not from proper (matias)

     if(kfl_delay(modul)/=0) then
        !
        ! TURBUL Module is delayed: put turbulent viscosity to zero
        !
        call tur_updunk(11_ip)

     else

        do ipoin=1,npoin
           if ( kfl_prope /= 0_ip ) then
              call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
              call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
           else
              call tur_nodpro(ipoin,rho(1),mu(1))
           end if
           nu = mu(1)/rho(1) 
           call tur_nut2nd(2_ip,ipoin,nu,nut,dumm1,dumm2)
           if( (postp(1) % npp_stepi (42)/=0.or.maxval(postp(1) % pos_times (1:nvart,42))>zetur)) then
             olded_tur(ipoin) = turmu(ipoin)
           end if
!            turmu(ipoin) = 1.0_rp*max(rho(1)*nut,mu(1)*1.0e-4_rp) + 0.0_rp*turmu(ipoin)
           turmu(ipoin) = rho(1)*nut
           if (kfl_logva==1) turmu(ipoin) = rho(1)*nut
        end do
     
     end if
  end if

end subroutine tur_updedd
