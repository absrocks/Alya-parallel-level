subroutine qua_tistep
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_tistep
  ! NAME 
  !    qua_tistep
  ! DESCRIPTION
  !    This routine sets the time step
  ! USES
  ! USED BY
  !    qua_begite
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_quanty
  use mod_outfor, only : outfor
  implicit none

  if(kfl_timco/=2) then
     
	 dtinv_qua=dtinv
     if(kfl_stead_qua==1) dtinv_qua = 0.0_rp
     if(kfl_timei_qua==0) dtinv_qua = 0.0_rp  

! no se que carajo hace aca!!!

!     if(kfl_tisch_qua==1) then
        !
        ! Trapezoidal rule: Euler iterations
        !
    !    if(ittim<=neule_qua) then
    !       kfl_tiacc_qua=1
    !    else
    !       kfl_tiacc_qua=kfl_tiaor_qua
    !    end if
    !    if(kfl_tiacc_qua==2) dtinv_qua = 2.0_rp*dtinv_qua
!     else
        !
        ! BDF scheme: increase integration order at each time step
        !
!        kfl_tiacc_qua=min(kfl_tiaor_qua,ittim)
!        call parbdf(kfl_tiacc_qua,pabdf_qua)
!     end if
  
  end if

  routp(1)=dtcri_qua
  routp(2)=0.0_rp
  routp(3)=0.0_rp
  ioutp(1)=kfl_timei_qua
  ioutp(2)=kfl_stead_qua

  call outfor(8_ip,lun_outpu,' ')

end subroutine qua_tistep
