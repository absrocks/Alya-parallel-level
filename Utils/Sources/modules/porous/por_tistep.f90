!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_tistep.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Sets the time step
!> @details Sets the time step
!> @} 
!------------------------------------------------------------------------
subroutine por_tistep
  use      def_parame
  use      def_master
  use      def_porous
  use mod_outfor, only : outfor
  implicit none

  if(kfl_timco/=2) then
     dtinv_por = dtinv
     if(kfl_stead_por==1) dtinv_por = 0.0_rp
!     if(kfl_timei_por==0) dtinv_por = 0.0_rp     ! esto deberái repensarlo ya que kfl_timei_por es un vector aqui

     if(kfl_tisch_por==1) then
        !
        ! Trapezoidal rule: Euler iterations
        !
!        if(ittim<=neule_por) then   we do not have euler iterations in porous

!           kfl_tiacc_por=1
!        else
!           kfl_tiacc_por=kfl_tiaor_por
!        end if
!        if(kfl_tiacc_por==2) dtinv_por = 2.0_rp*dtinv_por     ! esto deberái repensarlo ya que kfl_timei_por es un vector aqui
     else
        !
        ! BDF scheme: increase integration order at each time step
        !
!        kfl_tiacc_por=min(kfl_tiaor_por,ittim)
!        call parbdf(kfl_tiacc_por,pabdf_por)
     end if
  end if

  routp(1)=dtcri_por
  routp(2)=0.0_rp
  routp(3)=0.0_rp
  ioutp(1)=kfl_timei_por(1)
  ioutp(2)=kfl_stead_por
  call outfor(8_ip,lun_outpu,' ')

end subroutine por_tistep
