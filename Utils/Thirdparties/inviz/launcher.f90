subroutine launchviz(grank,gsize,lrank,lsize,dataptr,NDglobal,V)
  use def_kintyp
  use mod_parall, only : PAR_COMM_CURRENT
  implicit none  
  integer(ip) :: grank,gsize,lrank,lsize
  integer(ip) :: NDglobal,V
  real(rp)    :: dataptr(*)
  
  call insitu_communicator_init(PAR_COMM_CURRENT)
  call index_viewer_main(grank,gsize,lrank,lsize)
  
end subroutine launchviz

subroutine shutdownviz()
  implicit none

!  call index_viewer_shutdown()

end subroutine shutdownviz

subroutine updatedata(globaln,data,ND,V)
  use def_kintyp
  implicit none
  real(rp) :: data(*)
  integer(ip) :: ND,V
  integer(ip) :: globaln(*)
  
  call update_framedata(globaln,data,ND,V)
end subroutine updatedata

subroutine updateframe()
  use def_kintyp
  implicit none
  
  call render_frame()

end subroutine updateframe

subroutine launchasync()
  use def_kintyp
  use mod_parall, only : PAR_COMM_CURRENT, PAR_COMM_MY_CODE, PAR_COMM_MY_CODE_WM4
  use def_master, only : IMASTER
  implicit none
  integer(ip) :: flag =0
  
  if ( IMASTER ) then
     
  else

  end if
  
end subroutine launchasync
