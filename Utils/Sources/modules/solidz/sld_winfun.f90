  !-----------------------------------------------------------------------
  !> @addtogroup Solidz
  !> @{
  !> @file    sld_winfun.f90
  !> @author  Mariano Vazquez
  !> @date    08/01/2018
  !> @brief   Compute pressure Windkessel functions for boundary conditions 
  !> @details Compute pressure Windkessel functions for boundary conditions 
  !> @} 
  !-----------------------------------------------------------------------
  subroutine sld_winfun(itask)
    use def_parame
    use def_master
    use def_domain 
    use def_solidz

    implicit none

    integer(ip) :: itask !> itask 
    integer(ip) :: icavi 
    integer(ip), save :: flag_heart(4)
    integer(ip), save :: beat_num
    real(rp), save :: invol(4), esvol(4), edvol(4) !four cavities
    real(rp) :: deltavol, deltavol_prev, deltavol_sys(4), deltavol_dia(4),    &
     deltavol_pst(4), ddeltavol, Cvol, Cvel
    
    select case (itask)

    case (1_ip)

    ! Left Ventricle
    ! ======================================================================
    
       ! Calculate volume increments
       deltavol = volst(ITER_K,1) - volst(TIME_N,1)
       deltavol_prev = volst(TIME_N,1) - volst(TIME_N_MINUS_1,1)

       ! Calculate rate of volume change
       ddeltavol = (volst(ITER_K,1) - volst(TIME_N,1)) / dtime

       ! Three element Windkessel model
       if ((ittim == 1) .and. (itinn(modul) == 1)) then
          pwink(10:11,1) = part0_sld(1)
       else
       !elseif (flag_heart(1) == 2) then
          if (cutim > tzero_sld(1)) then

             ! Three element Windkessel
             !1 = 40.0_rp !40.0_rp
             !2 = 840.0_rp !840.0_rp
             !mod = 0.04_rp !0.04_rp
             !if (flag_heart(1) == 2) then
             !  pwink(10,1) = pwink(11,1) - ((pwink(11,1) * dtime) / (R2 *    &
             !   Cmod)) - (1.0_rp + (R1 / R2)) * (deltavol / Cmod) - R1 *     &
             !   ((deltavol - deltavol_prev) / dtime)
             !else
             !   pwink(10,1) = pwink(11,1) - ((pwink(11,1) * dtime) / (R2 *    &
             !    Cmod)) - (1.0_rp + (R1 / R2)) * (deltavol / Cmod) - R1 *     &
             !    (deltavol / dtime)
             !end if

             ! Two element Windkessel
             !if (pwink(TIME_N,1) > pwink_in(1)) then
             if (flag_heart(1) == 2) then
                !pwink(10,1) = pwink(11,1) + (dtime / C(1)) * (ddeltavol -              &
                ! (pwink(10,1) / R(1)))
                pwink(10,1) = pwink(11,1) - (dtime / cpres_sld(1)) *           &
                 (ddeltavol + (pwink(10,1) / rpres_sld(1)))
             else
                pwink(10,1) = part0_sld(1)
             end if
          end if
       end if

       ! First phase, initialisation
       if ((ittim == 1) .and. (itinn(modul) == 1)) then
         prest_sld(ITER_K) = 0.0_rp
         prest_sld(TIME_N) = 0.0_rp
         invol(1) = volst(ITER_K,1)
       end if

       if (cutim <= tzero_sld(1)) then

          ! Preload
          flag_heart(1) = 1   
          pwink(ITER_K,1) = pzero_sld(1) * (cutim / tzero_sld(1))
          edvol(1) = volst(ITER_K,1)

          ! Prestress
          deltavol_pst(1) = volst(ITER_K,1) - invol(1)
          if (cutim < tpstr_sld(1)) then
            prest_sld(ITER_K) = pstr0_sld(1) * (cutim / tpstr_sld(1))
          else
            Cvol = volst(TIME_N,1) / prest_sld(TIME_N)
            !Cvol = 1.0 / 10.0_rp
            !Cvol = volst(ITER_K,1) / prest_sld(ITER_K)
            Cvel = 1.0 / 2.5_rp !12.5_rp !2.0_rp ! Acceptable ranges are from 1.0 to 1.5
            !prest_sld(ITER_K) = prest_sld(TIME_N) + deltavol_pst(1) / Cvol +  &
            ! ddeltavol / Cvel
            prest_sld(ITER_K) = prest_sld(TIME_N) + deltavol_pst(1) / Cvol +  &
             ddeltavol / Cvel
            if (prest_sld(ITER_K) < 1.0_rp) prest_sld(ITER_K) = 1.0_rp

          end if

          ! First heart beat
          beat_num = 1

       ! Rest of the cardiac cycle
       else

          ! Isovolumetric contraction
          if (flag_heart(1) == 1) then
             !Cvol = volst(TIME_N,1) / pwink(TIME_N)
             Cvol = volst(ITER_K,1) / pwink(ITER_K,1) 
             Cvel = 1.0_rp / 10.0_rp ! Acceptable ranges are from 0.5 to 1.0
             deltavol_sys(1) = volst(ITER_K,1) - invol(1)
             !ddeltavol_sys(1) = (volst(ITER_K,1) - edvol(1)) / dtime
             pwink(ITER_K,1) = pwink(TIME_N,1) - deltavol_sys(1) / Cvol -     &
              ddeltavol / Cvel

          ! Ejection
          elseif (flag_heart(1) == 2) then

             ! Ventricular is equal to aortic pressure
             pwink(ITER_K,1) = pwink(10,1)
             esvol(1) = volst(ITER_K,1)     
 
          ! Isovolumetric relaxation
          elseif (flag_heart(1) == 3) then
             !Cvol = volst(TIME_N,1) / pwink(TIME_N)
             Cvol = volst(ITER_K,1) / pwink(ITER_K,1)
             Cvel = 1.0_rp / 10.0_rp ! Acceptable ranges are from 0.5 to 1.0
             deltavol_dia(1) = volst(ITER_K,1) - esvol(1)
             !ddeltavol_dia(1) = (volst(ITER_K,1) - esvol(1)) / dtime
             pwink(ITER_K,1) = pwink(TIME_N,1) - deltavol_dia(1) / Cvol -     &
              ddeltavol / Cvel

          ! Filling
          elseif (flag_heart(1) == 4) then
            !Cvol = volst(ITER_K,1) / pwink(ITER_K,1)
            !Cvol = 1.0_rp / 0.05_rp
            !Cvel = 1.0_rp / 10.0_rp  !0.5_rp
            !deltavol_sys(1) = volst(ITER_K,1) - invol(1)
            !pwink(ITER_K,1) = pwink(TIME_N,1) - deltavol_sys(1) / Cvol -      &
            !  ddeltavol / Cvel
            !K(1) = 500.0_rp !1.0_rp
            pwink(ITER_K,1) = pwink(TIME_N,1) - gfill_sld(1) * deltavol
            !pwink(ITER_K,1) = preload

          end if
       end if
       ! ======================================================================

       if (mcavi_sld > 1) then

          ! Right ventricle
          ! ======================================================================

          ! ======================================================================

       end if

    case (2_ip)
 
       ! update pwink for the next time step
       pwink(11,1) = pwink(10,1)
       pwink(TIME_N,1:4) = pwink(ITER_K,1:4)       
       prest_sld(TIME_N) = prest_sld(ITER_K)
       deltavol = volst(TIME_N,1) - volst(TIME_N_MINUS_1,1)

       ! Check for phase transition
       if (cutim > tzero_sld(1)) then
          if ((flag_heart(1) == 1) .and. (pwink(TIME_N,1) > pwink(11,1))) then
             flag_heart(1) = 2
          elseif ((flag_heart(1) == 2) .and. (deltavol > 0.000001_rp) .and.   &
           (volst(TIME_N,1) < 0.9_rp * edvol(1))) then
             flag_heart(1) = 3
          elseif ((flag_heart(1) == 3) .and. (pwink(TIME_N,1) < ppost_sld(1))) &
           then
             flag_heart(1) = 4
          !elseif ((flag_heart(1) == 4) .and. (cutim >= beat_num * 0.857_rp))  &
          ! then
          !   flag_heart(1) = 1
          !   beat_num = beat_num + 1
          end if
       end if

    end select
  end subroutine sld_winfun
