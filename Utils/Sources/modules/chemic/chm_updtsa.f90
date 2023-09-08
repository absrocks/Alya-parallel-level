subroutine chm_updtsa(dtmin)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updtsa
  ! NAME 
  !    chm_updtsa
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
  integer(ip)               :: ipoin,ispec,ipoi1,ipoi2,ii
  real(rp)                  :: dtmax,dtht,dttr,dcdt1,dcdt2,facc
  real(rp),    target       :: dtpar(1)

  if( ittim <= 4 ) then

     dtmin = dtold(1)

  else

     ! Adaptive timestep
     ! Prediction of next timestep dtn+1 based on accuracy
     ! reached with previous timestep dtold.

     dtmax = 0.0_rp

     if( INOTMASTER ) then

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

                 !  - Calculation of first derivative from values of concentrations at tn and tn-1.
                 !  dC/dt)n = (Cn - Cn-1)/dtn

                 dcdt1 = (conce(ipoin,ispec,3)-conce(ipoin,ispec,4))/dtold(1)

                 !  - Calculation of second derivative from values of concentrations
                 !  at previous times tn, tn-1 and tn-2.
                 !  d2C/dt2)n = ((Cn - Cn-1) - (Cn-1 - Cn-2)*dtn/dtn-1)/dtn**2

                 dcdt2=&
                      ((conce(ipoin,ispec,3) - conce(ipoin,ispec,4)) -                     &
                      (conce(ipoin,ispec,4)  - conce(ipoin,ispec,5)) * dtold(1)/dtold(2) ) &
                      / (dtold(1)*dtold(1))

                 !   - Criterion of influence of higher order terms on C:
                 !     1/2|dC2/dt2|*dt_ht**2=eps_trunc*|C| --> dt_ht=...

                 if( dcdt2 /= 0.0_rp ) then
                    dtht = sqrt(2.0_rp * epsht_chm * abs(conce(ipoin,ispec,3)) / abs(dcdt2))
                 else
                    dtht = dtold(1)
                 end if

                 !   - Criterion of truncation error due to time discretization
                 !     1/2|d2C/dt2|*dt_trunc=eps_HT*|dC/dt| --> dt_trunc=...

                 if( dcdt2 /= 0.0_rp ) then
                    dttr = 2.0_rp * epsht_chm * abs(dcdt1) / abs(dcdt2)
                 else
                    dttr = dtold(1)
                 end if

                 ! epstr_chm and epsht_chm to be defined by user in INPUT (for the moment)

                 ! Storing dtmax

                 dtmax=max(dtmax,dtht,dttr)

              end do
           end do
        end do
     end if
     !
     ! Look for maximum over whole mesh
     !
     if( IPARALL ) then
        dtpar =  dtmax
        nparr =  1
        parre => dtpar
        call Parall(10_ip)
        dtmax =  dtpar(1)
     end if

     !   - Calculation of accuracy reached with precious timestep dtold:
     !      facc=max(dt_trunc,dt_HT)/dtold

     facc = dtmax / dtold(1)

     !   - Prediction of next time step dtn+1 to be used in the solver:
     !   dtn+1=dtold * min(strech,ln(1+(damp-1)*facc)/ln(damp))
     !  
     !   strech: strech factor to limit timestep increase such thate dtnew <= P2.dtold.
     !   strech_m to be defined by the user in INPUT. strech >1. (eg: strech=1.05 means 5%).
     !   damp: damping factor to avoid oscillations of the timestep.
     !   Given by user in INPUT. damp > 1.0.

     dtmin = dtold(1) * min(strec_chm,log(1.0_rp + (dampi_chm-1.0_rp)*facc) / log(dampi_chm))
     dtmin = max(1.0e-8_rp,dtmin) 

  end if

end subroutine chm_updtsa
