module Chemic_Utils
  !
  !  FIXED POINT ITERATION ALGORITHM
  !  For now maybe it works, but eventually we could upgrade this to a faster option like Newton-Raphson (unlikely, the derivative is not easy to compute)
  !  or some trapezoidal-like thing
  !
  use def_kintyp, only      :  ip,rp
  use def_chemic, only      :  nreac_chm,nspec_chm, &
       fixed_point_tolerance_chm, max_fixed_point_iterations_chm
  use def_master, only      :  speci, prthe
  use def_kermod, only      :  gasco

  implicit none

  type iter_vector
     real(rp),pointer  :: conce(:)
     real(rp)  :: temp
  end type iter_vector

  real(rp),pointer    :: difnu(:,:)  ! Stoichiometric coefficients: backward-forward
  real(rp),pointer    :: T_min_cp(:), T_max_cp(:)
  integer(ip),pointer :: i_max_cp(:,:)

real(rp) :: mierda

  interface iterate_fixed_point
     module procedure iterate_fixed_point_only_conce, iterate_fixed_point_with_tempe
  end interface iterate_fixed_point

contains
  !
  ! Utilities
  !
  function vector_distance(vec1,vec2) result(dist)
    type(iter_vector),intent(in) :: vec1,vec2
    real(rp) :: dist
    ! Weak convergence, all species have at most tolerance error
    dist = max(maxval(abs(vec1%conce - vec2%conce)), abs(vec1%temp - vec2%temp))

  end function vector_distance
  !
  !
  !
  subroutine iterate_fixed_point_only_conce(initial_value, dt, densi,conce_iplus1)
    type(iter_vector),intent(in) :: initial_value
    real(rp),intent(in)          :: dt, densi  ! Time step, density
    type(iter_vector),intent(inout):: conce_iplus1

    integer(ip)               :: ireac,ispec,iterations
    real(rp)                  :: Q(nreac_chm)  ! Progress rate of each reaction
    type(iter_vector)         :: conce_i       ! Mass ratios for fixed-point iterations
    real(rp)                  :: source_xxx(nspec_chm),convergence,dummy

    ! Initialize
    allocate(conce_i%conce(nspec_chm))
!!    allocate(conce_iplus1%conce(nspec_chm))
    conce_iplus1%conce=initial_value%conce
    conce_iplus1%Temp =initial_value%Temp
    conce_i%Temp      =initial_value%Temp
    iterations = 0
    convergence=1.0D6
    ! Start fixed point iterations
    do 
       conce_i%conce=conce_iplus1%conce
       ! Compute heat release for each reaction
       call chm_heat_of_reactions(initial_value%Temp,densi,conce_iplus1%conce,Q)
       do ispec= 1,nspec_chm
          source_xxx(ispec) = 0.0_rp
          do ireac = 1,nreac_chm
             source_xxx(ispec) = source_xxx(ispec) + speci(ispec)%weigh * difnu(ireac,ispec) * Q(ireac)
          enddo
       enddo
       ! Evolve       
       conce_iplus1%conce= initial_value%conce + dt * source_xxx / densi
       convergence = vector_distance(conce_iplus1, conce_i)
       iterations=iterations+1
       if (iterations >= max_fixed_point_iterations_chm .or. convergence < fixed_point_tolerance_chm) then
          exit
       endif
    end do
    deallocate(conce_i%conce)
  end subroutine iterate_fixed_point_only_conce



  subroutine iterate_fixed_point_with_tempe(initial_value, dt,conce_iplus1)
    type(iter_vector),intent(in) :: initial_value
    real(rp),intent(in)  :: dt  ! Time step
    type(iter_vector),intent(inout)    :: conce_iplus1

    integer(ip)          :: ireac,ispec,iterations,icoef
    real(rp)             :: Q(nreac_chm)  ! Progress rate of each reaction
    type(iter_vector)    :: conce_i       ! Mass ratios for fixed-point iterations
    real(rp)             :: source_xxx(nspec_chm),convergence
    real(rp)             :: specific_heat(nspec_chm), enthalpy(nspec_chm),heat_source,mol_weight_inv, total_cp,densi
    real(rp)             :: c(0:7), T, dummy

    ! Initialize
    allocate(conce_i%conce(nspec_chm))
    conce_iplus1%conce=initial_value%conce
    conce_iplus1%Temp =initial_value%Temp
    iterations = 0
    convergence=1.0D6

    ! Start fixed point iterations
    do 
       conce_i%conce=conce_iplus1%conce
       conce_i%temp =conce_iplus1%temp
       !Compute mean molecular weight
       mol_weight_inv = 0.0_rp
       do ispec = 1,nspec_chm
          mol_weight_inv = mol_weight_inv +  conce_iplus1%conce(ispec) / speci(ispec)%weigh 
       enddo
       !Compute density
       densi = prthe(1) /(gasco * conce_iplus1%temp * mol_weight_inv) 
       !Compute specific heat and enthalpies
       do ispec= 1_ip,nspec_chm
          T=conce_iplus1%temp/1000.0_rp  ! This 1000 is just from the Shomate eq.
          if ( T < T_min_cp(ispec)) then ! Temp below minimum, assume constant
             T=T_min_cp(ispec)
             icoef=1
          else if ( T > T_max_cp(ispec)) then ! Temp above max range, assume constant
             T=T_max_cp(ispec)
             icoef=i_max_cp(1,ispec)-1
          else
             ! Search the temperature range for the Cp coefficients
             icoef=1
             do while (  icoef < 6 .and. 1000.0_rp *T >= speci(ispec)%trang(icoef+1) ) !!We only have up to 4 ranges=5 temperatures
                icoef=icoef+1
             enddo
          endif
          c=speci(ispec)%cpcoe(1:8,icoef)
          specific_heat(ispec) = c(0) + c(1)*T + c(2)*T**2 + c(3)*T**3 + c(4)/T**2 !+ c(5)*sqrt(T) 
          enthalpy(ispec) = 1000.0_rp * &  ! For some reason Shomate eq. gives enthalpy in kJ so we transform to J
               (speci(ispec)%entha(1) + &
               c(0)*T + c(1)*T**2/2.0_rp + c(2)*T**3/3.0_rp + c(3)*T**4/4.0_rp - c(4)/T + c(6) - c(7) ) ! + (2.0_rp/3.0_rp)*c(5)*T**1.5
       enddo
       total_cp = 0.0_rp
       do ispec= 1_ip,nspec_chm
          total_cp = total_cp +specific_heat(ispec) * conce_iplus1%conce(ispec)
       enddo
       ! Compute heat release for each reaction
       call chm_heat_of_reactions(conce_iplus1%temp,densi,conce_iplus1%conce,Q)
       do ispec= 1,nspec_chm
          source_xxx(ispec) = 0.0_rp
          do ireac = 1,nreac_chm
             source_xxx(ispec) = source_xxx(ispec) + speci(ispec)%weigh * difnu(ireac,ispec) * Q(ireac)
          enddo
       enddo
       !compute heat source
       heat_source = 0.0_rp
       do ispec= 1_ip,nspec_chm
          heat_source = heat_source - enthalpy(ispec) * source_xxx(ispec)
       enddo
       ! Finally, EVOLVE
       conce_iplus1%conce= initial_value%conce + dt * source_xxx / densi
       conce_iplus1%temp = initial_value%temp  + dt * heat_source / (densi * total_cp)
       convergence = vector_distance(conce_iplus1, conce_i) 
       iterations=iterations+1
       if (iterations >= max_fixed_point_iterations_chm .or. convergence < fixed_point_tolerance_chm) then
         ! print *,'out by ',iterations >= max_fixed_point_iterations_chm, convergence < fixed_point_tolerance_chm,conce_iplus1%temp, dt ,heat_source , densi ,( total_cp)
          exit
       endif
    end do
    deallocate(conce_i%conce)
    mierda = mierda +  dt * heat_source / (densi * total_cp)
  end subroutine iterate_fixed_point_with_tempe


end module Chemic_Utils


subroutine chm_updrea(dt_total)
  !-----------------------------------------------------------------------
  !    This routine computes the time dependency of the concentration
  !    as a function of the reaction terms (source terms in ADR language)
  !    This is used in fractional combustion model, and assumes that the ADR part of the
  !    species equations has already been computed.  
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,npoin
  use def_master, only      :  tempe,conce,inotmaster,speci
  use def_chemic
  use mod_ker_proper 
  use def_kermod
  use Chemic_Utils

  implicit none
  real(rp),intent(in)       :: dt_total
  real(rp)                  :: dt_current,dt_next,dt_accumulated
  integer(ip)               :: ireac,ispec,ipoin,icoef
  real(rp)                  :: T  ! Temperature
  real(rp)                  :: densi_xxx(npoin) ! Density from kernel
  type(iter_vector)         :: conce_current,conce_next  ! Mass ratios at current and next time
  type(iter_vector)         :: conce_full_step, conce_half_step_1, conce_half_step_2  ! Mass ratios for adaptive time stepping
  integer(ip)               :: dummi
  real(rp)                  :: dummr(3),error
  ! real(rp),external         :: iterate_fixed_point
  ! Init

  if (INOTMASTER) then

     allocate(difnu(nreac_chm,nspec_chm))
     difnu = 0.0_rp
     do ireac = 1,nreac_chm
        do icoef=6, 5+lreac_chm(ireac)%l(1) !Forward
           ispec = lreac_chm(ireac)%l(  icoef )
           difnu(ireac,ispec) = difnu(ireac,ispec)-stoic_chm(ispec,ireac,1)
        enddo
        do icoef=6 +lreac_chm(ireac)%l(1),5 +lreac_chm(ireac)%l(1)+lreac_chm(ireac)%l(2) !Backward
           ispec = lreac_chm(ireac)%l(icoef )
           difnu(ireac,ispec) = difnu(ireac,ispec)+stoic_chm(ispec,ireac,2)
        enddo
     enddo
     allocate(conce_current%conce(nspec_chm))
     allocate(conce_next%conce(nspec_chm))
     allocate(conce_full_step%conce(nspec_chm))
     allocate(conce_half_step_1%conce(nspec_chm))
     allocate(conce_half_step_2%conce(nspec_chm))
     if (kfl_stagg_chm == 1) then
        call ker_proper('DENSI','NPOIN',dummi,dummi,densi_xxx)
        do ipoin = 1,npoin
           select case (lawte_chm)
           case(0)
              T=temma_chm(1)
           case(-2)
              T=tempe_chm(ipoin)
           end select
           !
           ! Time evolution: Backward Euler (implicit), because it is A-stable and L-stable
           ! and therefore good for stiff systems - -even though it is only second order in time.
           !
           dt_accumulated=0.0_rp
           dt_current=dt_total/initial_fraction_step_chm
           dt_next=dt_current
           conce_current%conce=conce(ipoin,:,1)
           conce_next%conce=   conce(ipoin,:,1)
           conce_current%temp = T
           conce_next%temp = T
           conce_full_step%temp = T
           conce_half_step_1%temp = T 
           conce_half_step_2%temp = T
           do while (dt_accumulated < dt_total)
              !
              ! We compute a full time step and then 2 half steps, to estimate the error
              !
              do ! Adaptive time step
                 ! Full step:
                 call iterate_fixed_point(conce_current,dt_current, densi_xxx(ipoin),conce_full_step)
                 ! Two half steps:
                 call iterate_fixed_point(conce_current,dt_current/2.0_rp, densi_xxx(ipoin), conce_half_step_1)
                 call iterate_fixed_point(conce_half_step_1, dt_current/2.0_rp, densi_xxx(ipoin), conce_half_step_2)
                 ! Now we estimate the error 
                 error=vector_distance(conce_half_step_2, conce_full_step) !Weak convergence
                 !if (error.ne.0.0_rp) dt_next = 0.9*dt_current*odeint_tolerance_chm/error ! Estimate next time step
                 if (dt_next < timestep_min_chm) then
                    dt_next=timestep_min_chm
                    conce_next%conce =2.0_rp* conce_half_step_2%conce - conce_full_step%conce  ! conce_half_step_2 !
                    exit
                 else if (error < odeint_tolerance_chm ) then ! Successful step
                    ! Richardson extrapolation, increases approximation to dt^3
                    dt_next=2.0*dt_current
                    conce_next%conce =2.0_rp* conce_half_step_2%conce - conce_full_step%conce  ! conce_half_step_2 !
                    exit
                 else
                    dt_current=0.5*dt_current
                 end if
              enddo
              conce_current%conce = conce_next%conce
              dt_accumulated=dt_accumulated+dt_current
              dt_current=min(dt_next, dt_total-dt_accumulated, dt_total/initial_fraction_step_chm)

           enddo ! Time stepping
           conce(ipoin,:,2)= conce_current%conce
        enddo ! Loop over nodes

     else 

        allocate(T_min_cp(nspec_chm))
        allocate(i_max_cp(1,nspec_chm))
        allocate(T_max_cp(nspec_chm))
        do ispec= 1_ip,nspec_chm
           T_min_cp(ispec) = speci(ispec)%trang(1)/1000.0_rp
           i_max_cp(:,ispec) = maxloc( speci(ispec)%trang ) 
           T_max_cp(ispec) = speci(ispec)%trang(i_max_cp(1,ispec)) /1000.0_rp
        enddo
        !
        mierda=0.0_rp
        do ipoin = 1,npoin
           T=tempe(ipoin,1)
           !
           ! Time evolution: Backward Euler (implicit), because it is A-stable and L-stable
           ! and therefore good for stiff systems - -even though it is only second order in time.
           !
           dt_accumulated=0.0_rp
           dt_current=dt_total/initial_fraction_step_chm
           dt_next=dt_current
           conce_current%conce=conce(ipoin,:,1)
           conce_next%conce=conce(ipoin,:,1)
           conce_current%temp = T
           conce_next%temp = T
           conce_full_step%temp = T
           conce_half_step_1%temp = T 
           conce_half_step_2%temp = T
           do while (dt_accumulated < dt_total)
              !
              ! We compute a full time step and then 2 half steps, to estimate the error
              !
              do ! Adaptive time step
                 ! Full step:
                 call iterate_fixed_point(conce_current,dt_current, conce_full_step)
!                 print *,'P:',ipoin,dt_current,conce_current%temp,conce_full_step%temp
                 ! Two half steps:
                 call iterate_fixed_point(conce_current,     dt_current/2.0_rp, conce_half_step_1)
                 call iterate_fixed_point(conce_half_step_1, dt_current/2.0_rp, conce_half_step_2)
                 ! Now we estimate the error 
                 error=vector_distance(conce_half_step_2, conce_full_step) !Weak convergence
                 !if (error.ne.0.0_rp) dt_next = 0.9*dt_current*odeint_tolerance_chm/error ! Estimate next time step
                 if (dt_next < timestep_min_chm) then
                    dt_next=timestep_min_chm
                    conce_next%conce =2.0_rp* conce_half_step_2%conce - conce_full_step%conce  ! conce_half_step_2 !
                    conce_next%temp  =2.0_rp* conce_half_step_2%temp  - conce_full_step%temp   ! conce_half_step_2 !
                    exit
                 else if (error < odeint_tolerance_chm ) then ! Successful step
                    ! Richardson extrapolation, increases approximation to dt^3
                    dt_next=2.0*dt_current
                    conce_next%conce =2.0_rp* conce_half_step_2%conce - conce_full_step%conce  ! conce_half_step_2 !
                    conce_next%temp  =2.0_rp* conce_half_step_2%temp  - conce_full_step%temp   ! conce_half_step_2 !
                    exit
                 else
                    dt_current=0.5*dt_current
                 end if
              enddo
              conce_current%conce = conce_next%conce
              conce_current%temp = conce_next%Temp
              dt_accumulated=dt_accumulated+dt_current
              dt_current=min(dt_next, dt_total-dt_accumulated, dt_total/initial_fraction_step_chm)

           enddo ! Time stepping
           conce(ipoin,:,2)= conce_current%conce
!!!           mierda=mierda+abs((tempe(ipoin,1)-conce_current%temp))
           tempe(ipoin,1)= conce_current%temp
        enddo ! Loop over nodes
        !
        deallocate(T_min_cp)
        deallocate(T_max_cp)
        deallocate(i_max_cp)
     endif

     deallocate(difnu)
     deallocate(conce_current%conce)
     deallocate(conce_next%conce)
     deallocate(conce_full_step%conce)
     deallocate(conce_half_step_1%conce)
     deallocate(conce_half_step_2%conce)

  endif


1 format(1(1x,i7),10(1x,e12.6))

end subroutine chm_updrea

