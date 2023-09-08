subroutine ale_begste
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_begste
  ! NAME 
  !    ale_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the ALE formulation
  !    equation      
  ! USES
  !    ale_updunk
  ! USED BY
  !    Alefor
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  implicit none
  integer(ip) :: idime

  if ( kfl_rigid_ale == 1 ) then
     !
     ! Put force to zero - not needed here it is done inside ale_forces
     !
!     do iimbo = 1,nrbod
!        do idime = 1,ndime
!           rbbou(iimbo) % force(idime,1) = 0.0_rp
!        end do
!     end do
     !
     ! Obtain xrota_ale and xline_ale from nstro_ale and nstli_ale
     !
     xrota_ale = 0.0_rp
     xline_ale = 0.0_rp
     do idime =1,3_ip   
        if ( ittim >= nstro_ale(idime) )  xrota_ale(idime) = 1.0_rp
        if ( ittim >= nstli_ale(idime) )  xline_ale(idime) = 1.0_rp
     end do
     !
     ! Initial guess for rigid body unknowns: a(n,0,*) <-- a(n-1,*,*).
     !
     call ale_updunk(21_ip)

  else
     !
     ! Initial guess for dispm: d(n,*,*) <-- d(n-1,*,*) (important in strong coupling schemes for FSI)
     !
     call ale_updunk(1_ip)
     !
     ! Use temporal predictor for zonal coupling (if any)
     !
     call ale_coupre()

  end if

end subroutine ale_begste

