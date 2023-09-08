subroutine inimum()
!-----------------------------------------------------------------------
!
! This routine calls MUMPS for initialization and analysis
!
!-----------------------------------------------------------------------
!  use      def_master
!  use      def_solver
!  use      def_domain
!  implicit none
!  interface
!     subroutine dmumps(id)
!        include 'dmumps_struc.h'
!        type(dmumps_struc) id
!     end subroutine dmumps
!  end interface
!  include 'dmumps_struc.h'
!  include 'mpif.h'
!  integer(ip)       , intent(in)    :: ndofn,lusol
!  type(dmumps_struc), intent(inout) :: mumps
!  real(rp)                          :: cpu_refe1,cpu_refe2
!  integer(ip)                       :: ipoin,jpoin,idofn,jdofn,izdom,iztot
!!
!  call cputim(cpu_refe1)
!!
!! Call MUMPS for initialization
!!
!  mumps%JOB = -1
!  mumps%SYM = 0
!  mumps%PAR = 1
!  call dmumps(mumps)
!!
!! Build matrix structure
!!
!  mumps%N  = npoin*ndofn
!  mumps%NZ = nzsol*ndofn*ndofn
!  allocate(mumps%IRN(mumps%NZ))
!  allocate(mumps%JCN(mumps%NZ))
!  iztot=0
!  do ipoin=1,npoin
!     do izdom=r_sol(ipoin),r_sol(ipoin+1)-1
!        jpoin=c_sol(izdom)
!        do idofn=1,ndofn
!           do jdofn=1,ndofn
!              iztot=iztot+1
!              mumps%IRN(iztot)=(ipoin-1)*ndofn+idofn
!              mumps%JCN(iztot)=(jpoin-1)*ndofn+jdofn
!           end do
!        end do
!     end do
!  end do
!!
!! Call MUMPS for analysis
!!
!  mumps%JOB = 1
!  mumps%ICNTL(1)=lusol  ! output file 
!  mumps%ICNTL(3)=lusol  ! output file 
!  call dmumps(mumps)
!!
!! Compute CPU time 
!!
!  call cputim(cpu_refe2)
!  cpu_solve = cpu_solve + cpu_refe2 - cpu_refe1
!  
end subroutine inimum
