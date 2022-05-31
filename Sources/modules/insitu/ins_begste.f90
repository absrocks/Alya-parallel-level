subroutine ins_begste()
  use def_insitu
  use def_master, only    : IMASTER,INOTMASTER
  implicit none

#ifdef INVIZ
  call updateframe()
#endif
  
end subroutine ins_begste
