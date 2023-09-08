subroutine tem_endste()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_endste
  ! NAME 
  !    tem_endste
  ! DESCRIPTION
  !    This routine ends a time step of the temperature equation.
  ! USES
  !    tem_cvgunk
  !    tem_updunk
  !    tem_output
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_temper
  use def_domain,    only : npoin
  implicit none
  integer(ip) :: itask,ipoin,ivalu
  real(rp)    :: dummr,tenew,cploc(6,2)
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  !
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if(kfl_stead_tem==0.and.kfl_timei_tem==1) then

     call tem_cvgunk(three)
     call tem_updunk(five)

     if (kfl_regim_tem == 4) then 
       ! 
       ! Re-compute the enthalpy with the new values of C_p to ensure the correct T at the boundary 
       !
       if (kfl_regim_tem==4 .and. kfl_plepp_tem /= 4) call tem_calcEnthalpyBC()
       ! 
       ! Compute current temperature from enthalpy 
       !
       call tem_comtem
       ! 
       ! Temperature storage (:,3) <= (:,1) and high-order temporal scheme update
       !
       call tem_store_tempe(two) 

    end if

    call tem_updthe(2_ip) ! PRTHE_NSI: Thermodynamic pressure
  end if
  !
  ! Compute averaged variables
  !  
  call tem_averag()
  !
  ! Write restart file
  !
  call tem_restar(2_ip)
  !
  ! Coupling with dynamic solver
  !
  call tem_dyncou(2_ip)
  !
  ! If not steady, go on
  !
  if(kfl_stead_tem==0.and.kfl_timei_tem==1.and.kfl_conve(modul)==1) kfl_gotim = 1

end subroutine tem_endste
