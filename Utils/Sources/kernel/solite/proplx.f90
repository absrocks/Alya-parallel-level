subroutine proplx(nbnodes,nbvar,xx,yy,sumxx)

  !------------------------------------------------------------------------
  ! Sources/kernel/solite/proplx.f90
  ! NAME 
  !    proplx
  ! DESCRIPTION
  !    This routine computes the sclara product of two complex vectors:
  !                      SUMXX = <XX,YY>
  ! INPUT ARGUMENTS
  !    NBNODES .................. Number of nodes
  !    NBVAR .................... Number of variables per node
  !    XX(NBVAR*NBNODES) ........ Complex vector #1
  !    YY(NBVAR*NBNODES) ........ Complex vector #2
  ! OUTPUT ARGUMENTS
  !    SUMXX .................... Scalar (dot) product 
  ! USES
  ! USED BY 
  !------------------------------------------------------------------------


  use def_kintyp, only : ip,rp
  use def_master, only : kfl_paral,npoi1,npoi2,npoi3,parcx,npari,nparr,nparx,nparc
! Declaration statements
  implicit none

! Dummy arguments
  integer(ip), intent(in)  :: nbnodes,nbvar
  complex(rp), intent(in)  :: xx(*),yy(*)
  complex(rp), intent(out) :: sumxx

! Local variables
  integer(ip)              :: kk,temp
  complex(rp), target      :: dummx_par(1)

!--------------------------------------------------------------------------------------------------------
! Inner product: sumxx = <x, y> 
!--------------------------------------------------------------------------------------------------------

  sumxx = (0.0_rp,0.0_rp)
 
  if ( kfl_paral == -1_ip ) then
  !Sequential
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(kk)      &
! !$omp            schedule(static) &
! !$omp            reduction(+:sumxx)
	do kk = 1,nbvar*nbnodes
	  sumxx = sumxx + conjg(yy(kk)) * xx(kk) 
	end do
! !$omp end parallel do
  else if ( kfl_paral >= 1_ip ) then
  !Parallel: slaves
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(kk)      &
! !$omp            schedule(static) &
! !$omp            reduction(+:sumxx)
  	do kk = 1,nbvar*npoi1
      sumxx = sumxx + conjg(yy(kk)) * xx(kk) 
  	end do     
! !$omp end parallel do
    temp = nbvar*(npoi2 - 1_ip) + 1_ip
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(kk)      &
! !$omp            schedule(static) &
! !$omp            reduction(+:sumxx)
  	do kk = temp,nbvar*npoi3
      sumxx = sumxx + conjg(yy(kk)) * xx(kk) 
	  end do
! !$omp end parallel do
  end if

  !call pararr('SUM',0_ip,1_ip,sumxx)
  npari = 0_ip
  nparc = 0_ip
  nparr = 0_ip
  if ( kfl_paral >= 0_ip ) then
  !Parallel: reduce sum
    nparx        =  1_ip
    dummx_par(1) =  sumxx
    parcx        => dummx_par
    call Parall(9_ip) 
    sumxx        =  dummx_par(1)
  end if
 
end subroutine proplx


