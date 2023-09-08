module def_optsol
  !-----------------------------------------------------------------------
  !****f* Optsol/def_optsol
  ! NAME
  !    def_optsol
  ! DESCRIPTION
  !    Optsol module
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master

  real(rp) ::         &
       scalefactor,   &               ! factor used by the design update       d^{k+1} = d^k + scalefactor * stepj * descdir
       stepfactor                     ! factor used by the stepm length update alpha^{k+1} = alpha^k * stepfactor


end module def_optsol
