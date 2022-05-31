subroutine nsa_elfmom
!-----------------------------------------------------------------------
!****f* Nstinc/nsa_elfmom
! NAME 
!    nsa_elfmom
! DESCRIPTION
!    Fractional momentum per-element operations:
!    1. Compute elemental matrix and RHS 
!    2. Compute boundary contributions
!    3. Assemble 
! USES
!    nsa_...
! USED BY
!    nsa_gofmom
!***
!-----------------------------------------------------------------------
  implicit none

  call runend('NSA_ELFMOM: DEPRECATED ALGORITHM')

end subroutine nsa_elfmom
