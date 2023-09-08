subroutine nsa_vortic(iwopo,kvopo)
!-----------------------------------------------------------------------
!****f* nastal/nsa_vortex
! NAME 
!    nsa_vortex
! DESCRIPTION
!    Compute elemental vorticity and assemble it in the global vector
! USES
!    nsa_...
! USED BY
!    nsa_elcons
!    nsa_elchea
!***
!-----------------------------------------------------------------------
  use def_kintyp
  integer(ip) ::  kvopo,iwopo
  
  if (iwopo == 4 .or. iwopo == 29 .or. iwopo == 31 .or. iwopo == 32) then
     if (kvopo == 0) then
        kvopo = 1
        call vortic(1_ip)             ! compute vorticity
     end if
  end if
  
end subroutine nsa_vortic
