subroutine got_tistep
  !-----------------------------------------------------------------------
  !****f* Gotita/got_tistep
  ! NAME 
  !    got_tistep
  ! DESCRIPTION
  !    This routine sets the time step. 
  ! USES
  ! USED BY
  !    got_begite
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_gotita
  use mod_outfor, only : outfor
  implicit none

  if(kfl_timco/=2) then
     dtinv_got=dtinv
     if(kfl_timei_got==0) dtinv_got = 0.0_rp 
     if(kfl_stead_got==1) dtinv_got = 0.0_rp
     !
     ! Trapezoidal rule: Euler iterations
     !
     if(ittim<=neule_got) then
        kfl_tiacc_got=1
     else
        kfl_tiacc_got=kfl_tiaor_got
     end if
     if(kfl_tiacc_got==2) dtinv_got = 2.0_rp*dtinv_got
  end if

  routp(1)=dtcri_got
  ioutp(1)=kfl_timei_got
  ioutp(2)=kfl_stead_got
  routp(2)=timin_got
  routp(3)=timax_got
  call outfor(8_ip,lun_outpu,' ')

end subroutine got_tistep
