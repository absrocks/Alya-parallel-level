subroutine csrluf(ndof,nbrows,rhsid,unkno,amatr,ia,ja)

  !------------------------------------------------------------------------
  !****f* solite/solite
  ! NAME 
  !    solite
  ! DESCRIPTION
  !    Sparse direct solver
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------

  use def_kintyp, only : ip,rp
  use mod_csrdir
  implicit none 

  integer(ip), pointer      :: iL(:),jL(:)
  real(rp),    pointer      :: Ln(:)
  integer(ip), pointer      :: iU(:),jU(:)
  real(rp),    pointer      :: Un(:)

  real(rp),    intent(in)   :: rhsid(*)
  real(rp),    intent(out)  :: unkno(*)
  real(rp),    intent(in)   :: amatr(*)
  integer(ip), intent(in)   :: ia(*)
  integer(ip), intent(in)   :: ja(*)
  integer(ip)               :: nbrows
  integer(ip)               :: ndof
  integer(ip)               :: ii,sing
  integer(ip), pointer      :: invpR(:),invpC(:) 
  !
  ! Permutation arrays
  !
  nullify(Il,Jl,Ln,iU,jU,Un)
  nullify(invpR,invpC)
  !
  ! Initialization
  !
  do ii = 1, nbrows * ndof
     unkno(ii) = 0.0_rp 
  end do
  !
  ! CSR LU Factorization  
  !
  call CSR_LU_Factorization(&
       nbrows,ndof,ia,ja,amatr,IL,JL,&
       LN,IU,JU,UN,sing)
  !
  ! Errors
  !
  if( sing /= 0 ) write(*,*) 'Error!! Singularity detected.'
  !
  ! Solver for triangular linear systems of equations
  !
  call CSR_LUsol(&
       nbrows,ndof,invpR,invpC,IL,JL,&
       LN,IU,JU,UN,rhsid,unkno)
  !
  ! Deallocate memory
  !
  call CSR_LUfin(&
       IL,JL,LN,IU,JU,UN)

end subroutine csrluf
