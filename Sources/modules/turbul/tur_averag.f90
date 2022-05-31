subroutine tur_averag()
  !------------------------------------------------------------------------
  !****f* Turbul/tur_averag
  ! NAME 
  !    tur_averag
  ! DESCRIPTION
  !    Average key, omega, turvi
  ! USES
  ! USED BY
  !    tur_endste 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  implicit none
  integer(ip) :: idime,ipoin
  real(rp)    :: xfact,auxii

  if( INOTMASTER ) then

     if( cutim > avtim_tur ) then
        xfact     = 0.5_rp*dtime 

        ! AVKEY: average key

        if( output_postprocess_check_variable_postprocess(40_ip) ) then
          do ipoin = 1,npoin
              avkey_tur(ipoin)=avkey_tur(ipoin)&
                   +(untur(1,ipoin,1)+untur(1,ipoin,3)) * xfact
           end do
        end if

        ! AVOME: average omega

        if( output_postprocess_check_variable_postprocess(41_ip) ) then
          do ipoin = 1,npoin
              avome_tur(ipoin)=avome_tur(ipoin)&
                   +(untur(2,ipoin,1)+untur(2,ipoin,3)) * xfact
           end do
        end if

        ! AVTVI: average turbulent viscosity

        if( output_postprocess_check_variable_postprocess(42_ip) ) then
          do ipoin = 1,npoin
              avtvi_tur(ipoin)=avtvi_tur(ipoin)&
                   +(turmu(ipoin)+olded_tur(ipoin)) * xfact

           end do
        end if

     end if

  end if

end subroutine tur_averag

