subroutine solshu(nzmat,ndofn,npoin,rhsid,unkno,amatr)
  !------------------------------------------------------------------------
  !****f* solshu/solshu
  ! NAME 
  !    solshu
  ! DESCRIPTION
  !    Bridge to Schur solver
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_master, only    :   INOTMASTER
  use def_master, only    :   aii,aib,abi,abb,bbi,bbb,xxi,xxb
  use def_domain, only    :   npoin_ii
  use def_solver 
  implicit none
  integer(ip), intent(in) :: ndofn
  integer(ip), intent(in) :: nzmat
  integer(ip), intent(in) :: npoin
  real(rp),    target     :: rhsid(npoin*ndofn)
  real(rp),    target     :: unkno(npoin*ndofn)
  real(rp),    target     :: amatr(nzmat*ndofn)
  integer(ip)             :: poaii,poaib,poabi,poabb
  !
  ! Bridge to Schur complement solver
  ! 
  if( INOTMASTER ) then
     poaii =  1
     poaib =  poaii + nzdom_aii * solve_sol(1) % ndof2
     poabi =  poaib + nzdom_aib * solve_sol(1) % ndof2
     poabb =  poabi + nzdom_abi * solve_sol(1) % ndof2
     aii   => amatr(poaii:)
     aib   => amatr(poaib:)
     abi   => amatr(poabi:)
     abb   => amatr(poabb:)
     bbi   => rhsid
     bbb   => rhsid(npoin_ii*ndofn+1:)
     xxi   => unkno
     xxb   => unkno(npoin_ii*ndofn+1:)
  end if
  call Parall(427_ip)

end subroutine solshu
