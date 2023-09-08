subroutine wav_outset
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_outset
  ! NAME 
  !    wav_outset
  ! DESCRIPTION
  !    Calculation on sets
  ! USED BY
  !    wav_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  use mod_iofile
  implicit none
  integer(ip) :: inset,dummi

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if(maxval(postp(1) % npp_setsn)>0) then

     if( INOTMASTER ) then
        do inset=1,nnset
           if(lnsec(inset)/=0) then
              if(postp(1) % npp_setsn(1)/=0) postp(1) % vnset(1,inset)=wavam(lnsec(inset),1)
           end if
        end do
     end if
     call posdef(23_ip,dummi)

  end if

end subroutine wav_outset
