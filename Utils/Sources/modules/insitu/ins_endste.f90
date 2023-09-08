subroutine ins_endste()
  use def_insitu
  use def_parame
  use def_master
  use def_domain
  use mod_ker_proper

  implicit none

!  if(kfl_timei_exm==1)
  kfl_gotim = 1
end subroutine ins_endste
