subroutine nsi_endite(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_endite
  ! NAME 
  !    nsi_endite
  ! DESCRIPTION
  !    This routine checks convergence and updates unknowns at:
  !    - itask=1 The end of an internal iteration
  !    - itask=2 The end of the internal loop iteration
  ! USES
  !    nsi_cvgunk
  !    nsi_updunk
  ! USED BY
  !    nsi_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_ker_proper
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: itask

  integer(ip) :: ielem
  
  select case ( itask )

  case ( ITASK_ENDINN )

     !-------------------------------------------------------------------
     !
     !  Compute convergence residual of the internal iteration 
     !
     !-------------------------------------------------------------------

     if( NSI_MONOLITHIC ) then
        call nsi_updrel()
        call nsi_updunk(1500_ip)                 ! Relax:    UNKNO
     end if
     call nsi_cvgunk(ITASK_ENDINN)               ! Residual: ||VELOC(:,1)-UNKNO||
     call nsi_updunk(ITASK_ENDINN)               ! Update:   VELOC(:,1)=UNKNO
     !
     ! Solve Subgrid scale equation
     !
     call nsi_solsgs(2_ip)
     !
     ! Convergence and timings
     !
     call nsi_cvgunk(0_ip)   
     !
     ! Output matrix
     !
     call nsi_outite()
     !
     ! Update properties if needed
     !
     call ker_updpro(ITASK_ENDITE)

  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     !  Compute convergence residual of the external iteration 
     !
     !-------------------------------------------------------------------

     if( NSI_MONOLITHIC ) call livinf(16_ip,' ',itinn(modul))

     call nsi_cvgunk(ITASK_ENDITE)               ! Residual: ||VELOC(:,1)-VELOC(:,2)||
     call nsi_updunk(ITASK_ENDITE)               ! Update:   VELOC(:,2)=VELOC(:,1)
     !
     ! Compute forces on IB (particles)
     !
     call nsi_coupli(ITASK_ENDITE)
     call nsi_updunk(1400_ip)                    ! Update:   xx(:,2)=xx(:,1) - for RB variables

     !
     ! Couple with ADAN if flag is present
     !
     call nsi_cadan(2_ip)
     call nsi_cadan(3_ip)

  case ( ITASK_INNITE )
     
     !-------------------------------------------------------------------
     !
     !  Compute convergence After Runge-Kutta substep
     !
     !-------------------------------------------------------------------

     call nsi_updunk(ITASK_INNITE)              ! Update velocity     
     
  end select

end subroutine nsi_endite
 
