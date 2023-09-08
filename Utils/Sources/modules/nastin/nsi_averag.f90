subroutine nsi_averag()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_averag
  ! NAME 
  !    nsi_averag
  ! DESCRIPTION
  !    Average velocity and pressure
  ! USES
  ! USED BY
  !    nsi_endste 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_gradie
  use mod_projec,         only : projec_elements_to_nodes
  use def_kermod,         only : kfl_wlaav_ker,nsteps_ensemble,kfl_twola_ker,kfl_noslw_ker,kfl_waexl_ker

  implicit none
  integer(ip) :: idime,ipoin,iauxi
  real(rp)    :: auxii
  real(rp),    pointer    :: velgr(:,:)
  !  logical(lg) :: kensem

  if ( (kfl_wlaav_ker == 1_ip) .and. ( (kfl_waexl_ker == 0_ip) .or. (kfl_noslw_ker == 1_ip) ) .and. (ittim > 0_ip) ) then
     call nsi_wallav(2_ip)
  end if
  
  if( cutim > avtim_nsi ) then
     
     
     if( INOTMASTER ) then

        !AVVEL
        if( postp(1) % npp_stepi(21) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,21)) > zensi ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 avvel_nsi(idime,ipoin)=avvel_nsi(idime,ipoin)&
                      +veloc(idime,ipoin,1)*dtime
              end do
           end do
        end if

        ! AVPRE
        if( postp(1) % npp_stepi(28) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,28)) > zensi ) then
           do ipoin = 1,npoin
              avpre_nsi(ipoin)=avpre_nsi(ipoin)&
                   +press(ipoin,1)*dtime
           end do
        end if

        ! AVTAN
        if( postp(1) % npp_stepi(46) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,46)) > zensi ) then
           call memgen(zero,ndime,npoin)
           call nsi_memphy(4_ip)
           call nsi_denvis()
           call nsi_outtan(1_ip)     ! for no slip wall law, variational forces are obtained - comes in gevec
           call nsi_memphy(-4_ip)
           do ipoin=1,npoin
              do idime=1,ndime
                 avtan_nsi(idime,ipoin) = avtan_nsi(idime,ipoin)&
                      + gevec(idime,ipoin)*dtime
              end do
           end do
           call memgen(2_ip,ndime,npoin)
        end if

        ! V**2
        if( postp(1) % npp_stepi(57) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,57)) > zensi ) then
           do ipoin=1,npoin
              do idime=1,ndime
                 auxii = veloc(idime,ipoin,1)                 
                 avve2_nsi(idime,ipoin)=avve2_nsi(idime,ipoin)&
                      +auxii*auxii*dtime
              end do
           end do
        end if

        !Vx*Vy
        if( postp(1) % npp_stepi(58) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,58)) > zensi ) then
           do ipoin=1,npoin
              avvxy_nsi(1,ipoin)=avvxy_nsi(1,ipoin)&
                   +veloc(1,ipoin,1)*&
                   veloc(2,ipoin,1)*dtime
              if (ndime==3) then
                 avvxy_nsi(2,ipoin)=avvxy_nsi(2,ipoin)&
                      +veloc(2,ipoin,1)*&
                      veloc(3,ipoin,1)*dtime
                 avvxy_nsi(3,ipoin)=avvxy_nsi(3,ipoin)&
                      +veloc(3,ipoin,1)*&
                      veloc(1,ipoin,1)*dtime
              end if
           end do
        end if

        !P**2
        if( postp(1) % npp_stepi(59) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,59)) > zensi ) then
           do ipoin = 1,npoin
              avpr2_nsi(ipoin)=avpr2_nsi(ipoin)&
                   +press(ipoin,1)*press(ipoin,1)*dtime
           end do
        end if

        ! AVMUT
        if( postp(1) % npp_stepi(75) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,75)) > zensi ) then
           call memgen(zero,npoin,zero)
           call projec_elements_to_nodes(turmu_nsi,gesca)
           do ipoin=1,npoin
              avmut_nsi(ipoin) = avmut_nsi(ipoin)&
                   + gesca(ipoin)*dtime
           end do
           call memgen(2_ip,npoin,zero)
        end if

        ! AVSTX: mu_t * grad(u)
        if( postp(1) % npp_stepi(76) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,76)) > zensi ) then
           allocate(velgr(ndime,npoin))
           call gradie(veloc(1,1:npoin,1),velgr)
           call memgen(zero,npoin,zero)
           call projec_elements_to_nodes(turmu_nsi,gesca)
           do ipoin=1,npoin
              do idime=1,ndime
                 avstx_nsi(idime,ipoin) = avstx_nsi(idime,ipoin)&
                      + gesca(ipoin) * velgr(idime,ipoin) * dtime
              end do
           end do
           call memgen(2_ip,npoin,zero)
           deallocate(velgr)
        end if

        ! AVSTY: mu_t * grad(v)
        if( postp(1) % npp_stepi(77) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,77)) > zensi ) then
           allocate(velgr(ndime,npoin))
           call gradie(veloc(2,1:npoin,1),velgr)
           call memgen(zero,npoin,zero)
           call projec_elements_to_nodes(turmu_nsi,gesca)
           do ipoin=1,npoin
              do idime=1,ndime
                 avsty_nsi(idime,ipoin) = avsty_nsi(idime,ipoin)&
                      + gesca(ipoin) * velgr(idime,ipoin) * dtime
              end do
           end do
           call memgen(2_ip,npoin,zero)
           deallocate(velgr)
        end if

        ! AVSTZ: mu_t * grad(w)
        if( postp(1) % npp_stepi(78) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,78)) > zensi ) then
           allocate(velgr(ndime,npoin))
           call gradie(veloc(3,1:npoin,1),velgr)
           call memgen(zero,npoin,zero)
           call projec_elements_to_nodes(turmu_nsi,gesca)
           do ipoin=1,npoin
              do idime=1,ndime
                 avstz_nsi(idime,ipoin) = avstz_nsi(idime,ipoin)&
                      + gesca(ipoin) * velgr(idime,ipoin) * dtime
              end do
           end do
           call memgen(2_ip,npoin,zero)
           deallocate(velgr)
        end if
        
        
        if (nsteps_ensemble > 0 ) then

           !ENVEL
           if (( postp(1) % npp_stepi(86) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,86)) > zensi ) ) then
              iauxi = mod(ittim, postp(1) % npp_stepi(86) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 do ipoin=1,npoin
                    do idime=1,ndime
                       envel_nsi(idime,ipoin)=envel_nsi(idime,ipoin)&
                            +veloc(idime,ipoin,1)*dtime
                    end do
                 end do
              end if
           end if

           ! ENPRE
           if(( postp(1) % npp_stepi(87) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,87)) > zensi ) )  then
              iauxi = mod(ittim, postp(1) % npp_stepi(87) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 do ipoin = 1,npoin
                    enpre_nsi(ipoin)=enpre_nsi(ipoin)&
                         +press(ipoin,1)*dtime
                 end do
              end if
           end if

           ! V**2
           if(( postp(1) % npp_stepi(88) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,88)) > zensi ) )  then
              iauxi = mod(ittim, postp(1) % npp_stepi(88) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 do ipoin=1,npoin
                    do idime=1,ndime
                       auxii = veloc(idime,ipoin,1)                 
                       enve2_nsi(idime,ipoin)=enve2_nsi(idime,ipoin)&
                            +auxii*auxii*dtime
                    end do
                 end do
              end if
           end if

           !Vx*Vy
           if(( postp(1) % npp_stepi(89) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,89)) > zensi ) )  then
              iauxi = mod(ittim, postp(1) % npp_stepi(89) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 do ipoin=1,npoin
                    envxy_nsi(1,ipoin)=envxy_nsi(1,ipoin)&
                         +veloc(1,ipoin,1)*&
                         veloc(2,ipoin,1)*dtime
                    if (ndime==3) then
                       envxy_nsi(2,ipoin)=envxy_nsi(2,ipoin)&
                            +veloc(2,ipoin,1)*&
                            veloc(3,ipoin,1)*dtime
                       envxy_nsi(3,ipoin)=envxy_nsi(3,ipoin)&
                            +veloc(3,ipoin,1)*&
                            veloc(1,ipoin,1)*dtime
                    end if
                 end do
              end if
           end if

           !P**2
           if(( postp(1) % npp_stepi(90) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,90)) > zensi ) )  then
              iauxi = mod(ittim, postp(1) % npp_stepi(90) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 do ipoin = 1,npoin
                    enpr2_nsi(ipoin)=enpr2_nsi(ipoin)&
                         +press(ipoin,1)*press(ipoin,1)*dtime
                 end do
              end if
           end if


           ! ENTAN

           if(( postp(1) % npp_stepi(91) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,91)) > zensi ) )  then
              iauxi = mod(ittim, postp(1) % npp_stepi(91) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 call memgen(zero,ndime,npoin)
                 call nsi_memphy(4_ip)
                 call nsi_denvis()
                 call nsi_outtan(1_ip)  ! for no slip wall law, variational forces are obtained - comes in gevec
                 call nsi_memphy(-4_ip)
                 do ipoin=1,npoin
                    do idime=1,ndime
                       entan_nsi(idime,ipoin) = entan_nsi(idime,ipoin)&
                            + gevec(idime,ipoin)*dtime
                    end do
                 end do
                 call memgen(2_ip,ndime,npoin)
              end if
           end if

           ! ENMUT

           if(( postp(1) % npp_stepi(92) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,92)) > zensi ) )  then
              iauxi = mod(ittim, postp(1) % npp_stepi(92) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 call memgen(zero,npoin,zero)
                 call projec_elements_to_nodes(turmu_nsi,gesca)
                 do ipoin=1,npoin
                    enmut_nsi(ipoin) = enmut_nsi(ipoin)&
                         + gesca(ipoin)*dtime
                 end do
                 call memgen(2_ip,npoin,zero)
              end if
           end if

           ! ENSTX: mu_t * grad(u)
           if(( postp(1) % npp_stepi(93) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,93)) > zensi ) )  then
              iauxi = mod(ittim, postp(1) % npp_stepi(93) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 allocate(velgr(ndime,npoin))
                 call gradie(veloc(1,1:npoin,1),velgr)
                 call memgen(zero,npoin,zero)
                 call projec_elements_to_nodes(turmu_nsi,gesca)
                 do ipoin=1,npoin
                    do idime=1,ndime
                       enstx_nsi(idime,ipoin) = enstx_nsi(idime,ipoin)&
                            + gesca(ipoin) * velgr(idime,ipoin) * dtime
                    end do
                 end do
                 call memgen(2_ip,npoin,zero)
                 deallocate(velgr)
              end if
           end if

           ! ENSTY: mu_t * grad(v)
           if(( postp(1) % npp_stepi(94) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,94)) > zensi ) )  then
              iauxi = mod(ittim, postp(1) % npp_stepi(94) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 allocate(velgr(ndime,npoin))
                 call gradie(veloc(2,1:npoin,1),velgr)
                 call memgen(zero,npoin,zero)
                 call projec_elements_to_nodes(turmu_nsi,gesca)
                 do ipoin=1,npoin
                    do idime=1,ndime
                       ensty_nsi(idime,ipoin) = ensty_nsi(idime,ipoin)&
                            + gesca(ipoin) * velgr(idime,ipoin) * dtime
                    end do
                 end do
                 call memgen(2_ip,npoin,zero)
                 deallocate(velgr)
              end if
           end if

           ! ENSTZ: mu_t * grad(w)
           if(( postp(1) % npp_stepi(95) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,95)) > zensi ) )  then
              iauxi = mod(ittim, postp(1) % npp_stepi(95) )
              if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
                 allocate(velgr(ndime,npoin))
                 call gradie(veloc(3,1:npoin,1),velgr)
                 call memgen(zero,npoin,zero)
                 call projec_elements_to_nodes(turmu_nsi,gesca)
                 do ipoin=1,npoin
                    do idime=1,ndime
                       enstz_nsi(idime,ipoin) = enstz_nsi(idime,ipoin)&
                            + gesca(ipoin) * velgr(idime,ipoin) * dtime
                    end do
                 end do
                 call memgen(2_ip,npoin,zero)
                 deallocate(velgr)
              end if
           end if

        end if ! ensemble

        if ( kfl_wlaav_ker == 1_ip ) then
           call nsi_wallav(1_ip)
        end if

        if ( kfl_twola_ker == 1_ip ) then
           call nsi_tluave(1_ip)
        end if

        if ( kfl_twola_ker == 2_ip ) then
           call nsi_coutan()
        end if

     end if

  end if

end subroutine nsi_averag
