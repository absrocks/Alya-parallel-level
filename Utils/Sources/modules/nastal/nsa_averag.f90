!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_averag.f90
!> @author  Daniel Mira
!> @date    16/11/1966
!> @brief   Average velocity, normal and shear stress, momentum, pressure and temperature
!> @details Average velocity, normal and shear stress, momentum, pressure and temperature
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_averag()
  use def_parame
  use def_master
  use def_domain
  use def_nastal
  implicit none
  integer(ip) :: idime,kpoin,ipoin
  real(rp)    :: xfact,auxii

  if( cutim > avtim_nsa ) then
     xfact     = 0.5_rp*dtime

     if( INOTMASTER ) then

        !AVVEL: average velocity

        if( postp(1) % npp_stepi(55) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,55)) > zensa ) then
           do ipoin = 1,npoin
              do idime=1,ndime
                 avvel_nsa(idime,ipoin)=avvel_nsa(idime,ipoin)&
                      +(veloc(idime,ipoin,1)+veloc(idime,ipoin,3))*xfact
              end do
           end do
        endif

        ! AVVE2: average normal stress 

        if( postp(1) % npp_stepi(56) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,56)) > zensa ) then
           do ipoin = 1,npoin
              do idime=1,ndime
                 auxii                  = (veloc(idime,ipoin,1)+veloc(idime,ipoin,3))*0.5_rp                 
                 avve2_nsa(idime,ipoin) = avve2_nsa(idime,ipoin)&
                                          + auxii*auxii*dtime
              end do
           end do
        end if

        ! AVVXY: average shear stress

        if( postp(1) % npp_stepi(57) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,57)) > zensa ) then
           do ipoin = 1,npoin
              avvxy_nsa(1,ipoin)=avvxy_nsa(1,ipoin)&
                   +(veloc(1,ipoin,1)+veloc(1,ipoin,3))*0.5_rp*&
                   (veloc(2,ipoin,1)+veloc(2,ipoin,3))*xfact
              if (ndime==3) then
                 avvxy_nsa(2,ipoin)=avvxy_nsa(2,ipoin)&
                      +(veloc(2,ipoin,1)+veloc(2,ipoin,3))*0.5_rp*&
                      (veloc(3,ipoin,1)+veloc(3,ipoin,3))*xfact
                 avvxy_nsa(3,ipoin)=avvxy_nsa(3,ipoin)&
                      +(veloc(3,ipoin,1)+veloc(3,ipoin,3))*0.5_rp*&
                      (veloc(1,ipoin,1)+veloc(1,ipoin,3))*xfact
              end if
           end do
        end if

        !AVMOM: average momentum

        if( postp(1) % npp_stepi(58) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,58)) > zensa ) then
           do ipoin = 1,npoin
              do idime=1,ndime
                 avmom_nsa(idime,ipoin)=avmom_nsa(idime,ipoin)&
                      +(umome(idime,ipoin,1)+umome(idime,ipoin,3))*xfact
              end do
           end do
        endif

        ! AVPRE: average pressure

        if( postp(1) % npp_stepi(59) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,59)) > zensa ) then
           do ipoin = 1,npoin
              avpre_nsa(ipoin)=avpre_nsa(ipoin)&
                   +(press(ipoin,1)+press(ipoin,3))*xfact
           end do
        end if

        ! AVTEM: average temperature

        if( postp(1) % npp_stepi(60) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,60)) > zensa ) then
          do ipoin = 1,npoin
              avtem_nsa(ipoin)=avtem_nsa(ipoin)&
                   +(tempe(ipoin,1)+tempe(ipoin,3)) * xfact
           end do
        end if

     end if

  end if

end subroutine nsa_averag
