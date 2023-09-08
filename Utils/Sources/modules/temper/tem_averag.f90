subroutine tem_averag()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_averag
  ! NAME 
  !    tem_averag
  ! DESCRIPTION
  !    This routine averages the temperature along time
  ! USES
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod, only       :  gasco 
  implicit none
  integer(ip) :: ipoin, idime
  real(rp)    :: fact0
  real(rp), pointer         :: prope_tmp(:)

  if (cutim > avtim_tem) then

     
     if( INOTMASTER ) then
        
        ! AVTEM
        if( postp(1) % npp_stepi(5) /= 0 ) then        
           do ipoin=1,npoin
              avtem_tem(ipoin)=avtem_tem(ipoin)&
                   +dtime*tempe(ipoin,1)
           end do
        end if
        
        ! TEM**2
        if( postp(1) % npp_stepi(18) /= 0 ) then        
           do ipoin=1,npoin
              avte2_tem(ipoin)=avte2_tem(ipoin)&
                   +tempe(ipoin,1) &
                   *tempe(ipoin,1)*dtime 
           end do
        end if

        ! VEL*TEM
        if( postp(1) % npp_stepi(19) /= 0 ) then        
           do ipoin=1,npoin
              fact0 = tempe(ipoin,1)
              do idime =1, ndime
                 avtev_tem(idime, ipoin)=avtev_tem(idime, ipoin)&
                   +fact0*veloc(idime, ipoin,1)*dtime
              end do
           end do
        end if
        ! AVDEN
        if( postp(1) % npp_stepi(20) /= 0 ) then
           if (kfl_regim_tem>=3) then
              fact0 = prthe(1)/gasco
              do ipoin=1,npoin
                 avden_tem(ipoin)=avden_tem(ipoin)&
                      +dtime*fact0/tempe(ipoin,1)
              end do
           end if
        end if
         ! RHO*VEL
        if( postp(1) % npp_stepi(21) /= 0 ) then        
           if (kfl_regim_tem>=3) then
              do ipoin=1,npoin
                 fact0=prthe(1)/gasco/tempe(ipoin,1)
                 do idime =1, ndime
                    fvvel_tem(idime, ipoin)=fvvel_tem(idime, ipoin)&
                         +fact0*veloc(idime, ipoin,1)*dtime
                 end do
              end do
           end if
        end if

        ! AVRES
        if( postp(1) % npp_stepi(32) /= 0 ) then
           nullify ( prope_tmp )
           allocate( prope_tmp(npoin) )
           prope_tmp(1:npoin) =  solve(1)%reaction(1,1:npoin)
           do ipoin=1,npoin
              avres_tem(ipoin)=avres_tem(ipoin)&
                   +dtime*prope_tmp(ipoin)
           end do
           deallocate( prope_tmp )
        end if

     end if
  end if

end subroutine tem_averag
