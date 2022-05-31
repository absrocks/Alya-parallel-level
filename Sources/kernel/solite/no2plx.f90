subroutine no2plx(nbnodes,nbvar,xx,sumxx)

  !---------------------------------------------------------------------------------
  ! Sources/kernel/solite/no2plx.f90
  ! NAME 
  !    no2plx
  ! DESCRIPTION
  !    This routine computes the 2-norm (Eucledian norm) of a complex vector XX:
  !    SUMXX = ||XX||_2 = ( sum_i |XX_i|^2 )^{1/2}
  ! INPUT ARGUMENTS
  !    NBNODES .................. Number of nodes
  !    NBVAR .................... Number of variables per node
  !    XX(NBVAR*NBNODES) ........ Vector
  ! OUTPUT ARGUMENTS
  !    SUMXX .................... Eucledian norm
  ! USES
  ! USED BY 
  !***
  !---------------------------------------------------------------------------------


  use def_kintyp, only : ip,rp
  use def_master, only : kfl_paral,npoi1,npoi2,npoi3,parre,npari,nparr,nparx,nparc
  !Declaration statements
  implicit none

  !Dummy arguments
  integer(ip), intent(in)  :: nbnodes,nbvar
  complex(rp), intent(in)  :: xx(*)
  real(rp),    intent(out) :: sumxx

  !Local variables
  integer(ip)              :: ii,temp
  real(rp),    target      :: dummr_par(1)

  !--------------------------------------------------------------------------------------------------------
  ! Euclidean norm: sumxx = ||x||
  !--------------------------------------------------------------------------------------------------------
  sumxx = 0.0_rp

  if (kfl_paral == -1_ip) then
    !Sequential
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(ii)      &
! !$omp            schedule(static) &
! !$omp            reduction(+:sumxx)
    do ii = 1,nbvar*nbnodes
      sumxx = sumxx + conjg(xx(ii)) * xx(ii)       !Inner product: sumxx = <xx, xx> = ||x||^2
    enddo
! !$omp end parallel do
  else if (kfl_paral >= 1_ip) then
    !Parallel: slaves
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(ii)      &
! !$omp            schedule(static) &
! !$omp            reduction(+:sumxx)
    do ii = 1,nbvar*npoi1
      sumxx = sumxx + conjg(xx(ii)) * xx(ii)
    enddo
! !$omp end parallel do
    temp = nbvar*(npoi2 - 1_ip) + 1_ip
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(ii)      &
! !$omp            schedule(static) &
! !$omp            reduction(+:sumxx)
    do ii = temp,nbvar*npoi3
      sumxx = sumxx + conjg(xx(ii)) * xx(ii)
    enddo
! !$omp end parallel do
  endif
  
  !call pararr('SUM',0_ip,1_ip,sumxx)
  !
  ! Parallel: reduce sum
  !
  npari = 0_ip
  nparc = 0_ip
  nparx = 0_ip
  if( kfl_paral >= 0 ) then
    nparr        =  1
    dummr_par(1) =  sumxx
    parre        => dummr_par
    call par_operat(3_ip)           
    sumxx        =  dummr_par(1)
  end if

  sumxx = sqrt(sumxx)

end subroutine no2plx

