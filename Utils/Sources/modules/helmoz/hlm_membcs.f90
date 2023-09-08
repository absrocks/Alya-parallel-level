subroutine hlm_membcs(itask)

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_membcs.f90
  ! NAME
  !    hlm_membcs
  ! DESCRIPTION
  !    This routine allocates memory for boundary conditions.
  ! OUTPUT 
  ! USES
  ! USED BY
  !    hlm_inibcs
  !-----------------------------------------------------------------------

  use def_master
  use def_domain
  use def_kermod
  use def_helmoz
  use mod_memory

  implicit none

  integer(ip), intent(in) :: itask

  integer(ip)             :: ipoin,iedge

  select case (itask)

  case ( 1_ip )

     if( kfl_edges_hlm == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_HLM','hlm_membcs',kfl_fixno_hlm,1_ip,meshe(ndivi) % nedge)
        do iedge = 1,meshe(ndivi) % nedge
           kfl_fixno_hlm(1,iedge) = -1_ip
        end do
     else
        call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_HLM','hlm_membcs',kfl_fixno_hlm,1_ip,npoin)
        do ipoin = 1,npoin
           kfl_fixno_hlm(1,ipoin) = -1_ip
        end do
     end if

  case (-1_ip )

     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_HLM','hlm_membcs',kfl_fixno_hlm)

  end select

end subroutine hlm_membcs
