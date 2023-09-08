!------------------------------------------------------------------------
!> @addtogroup Mathematics
!> @{
!> @file    spcdec.f90
!> @author  Herbert Owen
!> @brief   Computes eigenval. and eigenvec. of real, symmetrix matrix A (3x3). Sorts them in accending order.
!> @details Output is a vector D
!!          containing the eigenvalues in ascending order, and a matrix V whose
!!          columns contain the corresponding eigenvectors.
!> @}
!------------------------------------------------------------------------
subroutine spcdec(A,D,V,NROT,kfl_wivec,callersub)
  use def_kintyp

  implicit none
  real(rp),    intent(in)  :: A(3,3)
  real(rp),    intent(out) :: D(3), V(3,3)
  integer(ip), intent(out) :: NROT
  integer(ip), intent(in)  :: kfl_wivec   ! also obtain eigenvectors

  real(rp)                 :: E(3,3)
  real(rp)                 :: daux(3,10)!,dauxi(3),error(3)
  integer(ip)              :: idime,jdime,imeth

  character(*)             :: callersub  ! who is calling spcdec

  if ( kfl_wivec == 0_ip ) then
     imeth = 2_ip
  else
     imeth = 1_ip  ! perhaps there are better option still needs to be checked
  end if

  select case(imeth)

  case(1_ip)
     do idime=1,3
        do jdime=1,3
           E(idime,jdime) = A(idime,jdime)
        end do
     end do
     !call eigenv(E,D,V,NROT)
     call jacobi_nre(E,D,V,NROT,callersub)
     call eigsrt(D,V)

  case(2_ip) ! cardano only eigenval

     if ( kfl_wivec /= 1_ip) then

!       to check if cardano is working ok
!        do idime=1,3
!           do jdime=1,3
!              E(idime,jdime) = A(idime,jdime)
!           end do
!        end do
!        call jacobi_nre(E,D,V,NROT,callersub)
!        call eigsrt(D,V)
!        dauxi = D

        call DSYEVC3(A, D)     ! actually I unsderstand it does not touch the matrix so A could be used directly
        call eigsrt(D,V)       ! Here I should order only eigenvalues

!        error = D - dauxi
!        do idime=1,3
!           if ( D(idime) /= 0.0_rp) then
! this I left if unfinished becuase I did not know exactly what to compare? May be Fernado can give me some help on what to use

     else

        do idime=1,3
           do jdime=1,3
              E(idime,jdime) = A(idime,jdime)
           end do
        end do
        call DSYEVV3(E, V, D)  ! needs E
        call eigsrt(D,V)

     end if

  case(3_ip) ! Jacobi  - needs E

     do idime=1,3
        do jdime=1,3
           E(idime,jdime) = A(idime,jdime)
        end do
     end do
     call DSYEVJ3(E, V, D)
     call eigsrt(D,V)

  case(4_ip) ! Cuppen's Divide & Conquer algorithm - A could be used directly

     call DSYEVD3(A, V, D)
     call eigsrt(D,V)

  case(5_ip) ! Hybrid - A could be used directly

     call DSYEVH3(A, V, D)
     call eigsrt(D,V)

  case(6_ip) ! QL - A could be used directly

     call DSYEVQ3(A, V, D)
     call eigsrt(D,V)

  case(7_ip) ! here I could put all to compare

     do idime=1,3
        do jdime=1,3
           E(idime,jdime) = A(idime,jdime)
        end do
     end do
     call jacobi_nre(E,D,V,NROT,callersub)
     call eigsrt(D,V)
     daux(:,1)=D

     call DSYEVC3(A, D)     ! actually I unsderstand it does not touch the matrix so A could be used directly
     call eigsrt(D,V)
     daux(:,2)=D

     call DSYEVD3(A, V, D)
     call eigsrt(D,V)
     daux(:,3)=D

     call DSYEVH3(A, V, D)
     call eigsrt(D,V)
     daux(:,4)=D

     call DSYEVQ3(A, V, D)
     call eigsrt(D,V)
     daux(:,5)=D

     do idime=1,3
        do jdime=1,3
           E(idime,jdime) = A(idime,jdime)
        end do
     end do
     call DSYEVJ3(E, V, D) ! at the end because it alters E
     call eigsrt(D,V)
     daux(:,6)=D

     print*,'daux(1,1:6)',daux(1,1:6)
     print*,'daux(2,1:6)',daux(2,1:6)
     print*,'daux(3,1:6)',daux(3,1:6)

     print*,'diffdaux(1,1:6)',daux(1,2)-daux(1,1), daux(1,3)-daux(1,1), daux(1,4)-daux(1,1), daux(1,5)-daux(1,1),daux(1,6)-daux(1,1)
     print*,'diffdaux(2,1:6)',daux(2,2)-daux(2,1), daux(2,3)-daux(2,1), daux(2,4)-daux(2,1), daux(2,5)-daux(2,1),daux(2,6)-daux(2,1)
     print*,'diffdaux(3,1:6)',daux(3,2)-daux(3,1), daux(3,3)-daux(3,1), daux(3,4)-daux(3,1), daux(3,5)-daux(3,1),daux(3,6)-daux(3,1)

     print*,'**********************************************************'

  end select

end subroutine spcdec

