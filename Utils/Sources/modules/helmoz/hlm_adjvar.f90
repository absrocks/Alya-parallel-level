
subroutine hlm_adjvar
  !------------------------------------------------------------------------
  !****f* Helmoz/hlm_adjvar
  ! NAME
  !    hlm_adjvar
  ! DESCRIPTION
  ! OUTPUT
  ! USED BY
  !    Helmoz
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_inpout
  use      def_solver
  implicit none
  integer(ip) :: ii,nbnodes
  real(rp) :: normunknx,time1,time2
  
  !call moddef(10_ip)
  solad     => momod(modul) % solad
  solad_sol => solad
  solve_sol => solad_sol


  if (IMASTER) then
     nbnodes = 0_ip 
  else
     nbnodes = npoin
  endif

  do ii=1,size(aunknx)
     aunknx(ii)=cmplx(0.0_rp,0.0_rp,kind=rp)
  end do

  !normunknx=0.0_rp
  !call no2plx(nbnodes,4_ip,dcostx,normunknx)
  !write(*,*)'norm2(dcostx)=',normunknx

  call cputim(time1)
  call solvex(dcostx,aunknx,amatx,pmatx)
  call cputim(time2)


  call Parall(20_ip) !MPI_Barrier

  write(*,*)'kfl_paral=',kfl_paral,' iteration solver-adj time: ',time2-time1

  call moddef(9_ip)

end subroutine hlm_adjvar
