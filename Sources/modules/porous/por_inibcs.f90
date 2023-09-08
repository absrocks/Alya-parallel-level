subroutine por_inibcs()
  !-----------------------------------------------------------------------
  !****f* Temper/por_inibcs
  ! NAME
  !    por_inibcs
  ! DESCRIPTION
  !    This routine applied boundary conditions
  ! OUTPUT 
  ! USES
  ! USED BY
  !    por_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_porous
  use mod_opebcs
  use mod_memory
  implicit none

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Allocate memory
     !
     !-------------------------------------------------------------

     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_POR','por_inibcs',kfl_fixno_por,1_ip,npoin,2_ip)
     call memory_alloca(mem_modul(1:2,modul),'BVESS_POR',    'por_inibcs',bvess_por,    1_ip,npoin,2_ip)
 
     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then

        iffun     =  0
        ifloc     =  0
        ifbop     =  0
        kfl_fixno => kfl_fixno_por(:,:,1)
        bvess     => bvess_por(:,:,1) 
        tncod     => tncod_por(1:)
        call reacod(10_ip) 
        kfl_fixno => kfl_fixno_por(:,:,2)
        bvess     => bvess_por(:,:,2)
        tncod     => tncod_por(2:)
        call reacod(10_ip)

     end if

  end if

end subroutine por_inibcs
