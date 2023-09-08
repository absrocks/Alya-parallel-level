!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_averag.f90
!> @author  Daniel Mira
!> @date    13/06/204
!> @brief   Average variables: temperature
!> @details Average variables: temperature
!> @} 
!-----------------------------------------------------------------------
subroutine chm_averag()
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip) :: ipoin
  real(rp)    :: xfact,zechm
  real(rp), pointer         :: prope_tmp(:)

  zechm = epsilon(0.0_rp)

  if( cutim > avtim_chm ) then
     xfact     = 0.5_rp*dtime

     if( INOTMASTER ) then

        !
        ! AVTEM: average temperature
        !
        if( postp(1) % npp_stepi(28) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,28)) > zechm ) then
          do ipoin=1,npoin
             avtem_chm(ipoin) = avtem_chm(ipoin)&
                                + (tempe(ipoin,1)+tempe(ipoin,3)) * xfact
          end do
        end if
        
        !
        ! AVCON: average concentration
        !
        if( postp(1) % npp_stepi(30) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,30)) > zechm ) then
          do ipoin=1,npoin
             avcon_chm(ipoin) = avcon_chm(ipoin)&
                                    + (conce(ipoin,1,1)+conce(ipoin,1,3)) * xfact
          end do
        end if
        
        !
        ! AVCO2: average squared of concentration c*c
        !
        if( postp(1) % npp_stepi(31) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,31)) > zechm ) then
          do ipoin=1,npoin
             avco2_chm(ipoin) = avco2_chm(ipoin)&
                                      + (conce(ipoin,1,1)*conce(ipoin,1,1) + conce(ipoin,1,3)*conce(ipoin,1,3)) * xfact
          end do
        end if

        !
        ! AVVAR: average variance of reaction progress variable
        !
        if( postp(1) % npp_stepi(32) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,32)) > zechm ) then
          do ipoin=1,npoin
             avvar_chm(ipoin) = avvar_chm(ipoin)&
                                      + (conce(ipoin,2,1) + conce(ipoin,2,3)) * xfact
          end do
        end if

        !
        ! AVIME: average enthalpy scalar (IMEAN)
        !
        if( postp(1) % npp_stepi(33) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,33)) > zechm ) then
          do ipoin=1,npoin
             avime_chm(ipoin) = avime_chm(ipoin)&
                                      + dtime*encfi(ipoin)
          end do
        end if

        !
        ! AVCHM: average chemical heat (CHEMICAL_HEAT)
        !
        if( postp(1) % npp_stepi(34) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,34)) > zechm ) then

          nullify ( prope_tmp )
          allocate( prope_tmp(npoin) )

          call smooth (chemical_heat, prope_tmp)

          do ipoin=1,npoin
             avchm_chm(ipoin) = avchm_chm(ipoin)&
                                      + dtime*prope_tmp(ipoin)
          end do

          deallocate( prope_tmp )

        end if

        !
        ! AVMIX: average mixture fraction
        !
        if( postp(1) % npp_stepi(35) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,35)) > zechm ) then
          do ipoin=1,npoin
             avmix_chm(ipoin) = avmix_chm(ipoin)&
                                    + (conce(ipoin,3,1)+conce(ipoin,3,3)) * xfact
          end do
        end if

        !
        ! AVMI2: average squared of mixture fraction f*f
        !
        if( postp(1) % npp_stepi(36) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,36)) > zechm ) then
          do ipoin=1,npoin
             avmi2_chm(ipoin) = avmi2_chm(ipoin)&
                                      + (conce(ipoin,3,1)*conce(ipoin,3,1) + conce(ipoin,3,3)*conce(ipoin,3,3)) * xfact
          end do
        end if

     end if

  end if

end subroutine chm_averag
