subroutine tur_endite(itask)
!-----------------------------------------------------------------------
!****f* Turbul/tur_endite
! NAME 
!    tur_endite
! DESCRIPTION
!    This routine checks convergence and performs updates of the
!    turbulence variables at:
!    - itask=1 The end of an internal iteration
!    - itask=2 The end of the internal loop iteration
! USES
!    tur_cvgunk
!    tur_updunk
! USED BY
!    tur_doiter
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_turbul
  use def_kermod
  use mod_messages,   only : livinf
  use mod_ker_proper, only : ker_updpro
  implicit none
  integer(ip) :: itask

  select case(itask)

  case(one)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || f(n,i,j) - f(n,i,j-1)|| / ||f(n,i,j)||) and update unknowns:
     !  f(n,i,j-1) <-- f(n,i,j) 
     !
     call tur_updrel()     ! Compute relaxation factor
     call tur_updunk(8_ip) ! Relax UNKNO
     call tur_cvgunk(1_ip) ! Compute residual=UNTUR-UNKN
     call tur_updunk(3_ip) ! Actualize UNTUR=UNKNO and TURMU
     !
     ! Residual projection and subgrid scale equation
     !
     call tur_solsgs()
    
     call ker_updpro(ITASK_ENDITE)
  case(two)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || f(n,i,*) - f(n,i-1,*)|| / ||f(n,i,*)||) and update unknowns:
     !  f(n,i-1,*) <-- f(n,i,*) 
     !
     call livinf(16_ip,' ',itinn(modul))
     call tur_cvgunk(2_ip)
     call tur_updunk(4_ip)
     
  end select

end subroutine tur_endite
