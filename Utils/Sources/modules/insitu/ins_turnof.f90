subroutine ins_turnof()
  use def_insitu
  implicit none

  deallocate(advec)
  deallocate(elmag)
  deallocate(global_num)

#ifdef INVIZ
  call shutdownviz()
#endif
  
end subroutine ins_turnof
