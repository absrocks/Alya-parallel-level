subroutine chm_scale_cfi()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_scale_cfi
  ! NAME 
  !    chm_scale_cfi
  ! DESCRIPTION
  !    Scale reaction progress variable with equilibrium mass fractions after initialization
  ! USES
  ! USED BY
  !    chm_iniunk
  !   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_master, only      : table_cfi,conce
  use def_domain, only      : npoin
  use def_chemic, only      : kfl_uncfi_chm
  implicit none
  integer(ip)               :: ipoin,loc
  real(rp)                  :: f,c,delta_f,delta_y_c,y_c_eq

  do ipoin = 1,npoin

     f = conce(ipoin,3,1)
     c = conce(ipoin,1,1)

     !
     ! Impose equilibrium values in the flamability range 
     !
     if ( f > table_cfi(1)%fmima(1) .and. f < table_cfi(1)%fmima(2) ) then

       if (kfl_uncfi_chm == 1_ip ) then
          if (abs (c - 1.0_rp) < 1.0e-4_rp ) then

             loc = minloc(abs(table_cfi(1) % ymass(:,1) - f), 1)

             if (f < table_cfi(1) % ymass(loc,1)) loc = loc - 1_ip
             !
             ! Linear interpolation: Y_c = Y_c_1 + (f - f_1) * Delta_Y_c / Delta_f 
             !
             delta_f   = table_cfi(1) % ymass(loc + 1,1) - table_cfi(1) % ymass(loc,1)
             delta_y_c = table_cfi(1) % ymass(loc + 1,2) - table_cfi(1) % ymass(loc,2)

             if (delta_f /= 0.0_rp ) then
                y_c_eq = table_cfi(1) % ymass(loc,2) + (f - table_cfi(1) % ymass(loc,1)) * delta_y_c / delta_f
             else
                y_c_eq = table_cfi(1) % ymass(loc,2)
             endif
             !
             ! Impose equilibrium value: y_c_eq 
             !
             conce(ipoin,1,1) = y_c_eq
             conce(ipoin,1,2) = y_c_eq
             conce(ipoin,1,3) = y_c_eq 
          end if

       endif 
 
     else
       conce(ipoin,1,1) = 0.0_rp
       conce(ipoin,1,2) = 0.0_rp
       conce(ipoin,1,3) = 0.0_rp
     end if

  end do
  
end subroutine chm_scale_cfi
