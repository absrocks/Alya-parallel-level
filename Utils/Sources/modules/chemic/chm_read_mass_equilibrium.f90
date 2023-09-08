subroutine chm_read_mass_equilibrium()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_reasca
  ! NAME 
  !    chm_reasca
  ! DESCRIPTION
  !    Read equilibrium mass factions for scaling RPV for non-premixed conditions 
  ! USES
  ! USED BY
  !    chm_reaphy
  !   
  !***
  !-----------------------------------------------------------------------
  use def_inpout, only      : param,words
  use def_kintyp, only      : ip,rp
  use def_master, only      : table_cfi
  use mod_ecoute, only      : ecoute
  implicit none
  integer(ip)               :: ii,icoef
  real(rp)                  :: aux(3)

  !
  ! Initialization equilibrium mass fractions
  !
  table_cfi(1) % ymass = 0.0_rp

  call ecoute('chm_reaphy')
  
  if (words(1) == 'EQUIL') then
     call ecoute('chm_reaphy')

     do ii=1,table_cfi(1)%nvcfi(3)
        do icoef = 1,3
           table_cfi(1)%ymass(ii,icoef) = param(icoef)
        end do
        call ecoute('chm_reaphy')
     end do  
  endif

end subroutine chm_read_mass_equilibrium
