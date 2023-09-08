!------------------------------------------------------------------------
!> @addtogroup Mathematics
!> @{
!> @file    jacobi_nre.f90
!> @author  Herbert Owen
!> @brief   computes all the eigenvalues and eigenvectors of real, symmetrix matrix A (3x3)
!> @details The eigenvalues go into matrix D, while the
!!          eigenvectors go into matrix V. NROT is the number of Jacobi rotations
!!          which were required.
!!          Beware A is altered.
!!          *This subroutine is taken from "Numerical Recipes F90", page 1225
!!          In Alya there was a similar subrotine called sld_jacobi thar I renamed eigenv
!> @}
!------------------------------------------------------------------------

SUBROUTINE jacobi_nre(a,d,v,nrot,callersub)
  USE def_kintyp
  USE mod_nrutil, ONLY : assert_eq,get_diag,nrerror,unit_matrix,upper_triangle
  IMPLICIT NONE

  INTEGER(ip), INTENT(OUT) :: nrot
  !  REAL(rp), DIMENSION(:), INTENT(OUT) :: d
  !  REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: a
  !  REAL(rp), DIMENSION(:,:), INTENT(OUT) :: v
  REAL(rp), DIMENSION(3), INTENT(OUT) :: d     ! for the moment I restrict it to dim 3 because asset_eq gives me trouble
  REAL(rp), DIMENSION(3,3), INTENT(INOUT) :: a
  REAL(rp), DIMENSION(3,3), INTENT(OUT) :: v
  !
  ! Computes all eigenvalues and eigenvectors of a real symmetric N × N matrix a. On output, elements of a
  ! above the diagonal are destroyed. d is a vector of length N that returns the eigenvalues of a. v is
  ! an N × N matrix whose columns contain, on output, the normalized eigenvectors of a. nrot
  ! returns the number of Jacobi rotations that were required.
  !
  INTEGER(ip) :: i,ipp,iq,n,j
  REAL(rp) :: c,g,h,s,sm,t,tau,theta,tresh,toler
  REAL(rp), DIMENSION(size(d)) :: b,z
  !  n=assert_eq((/size(a,1),size(a,2),size(d),size(v,1),size(v,2)/),'jacobi')

  character(*)             :: callersub  ! who is calling spcdec
  character(300)           :: messa

  n = 3_ip
  toler = 1.0e-10            ! Tolerance
  call unit_matrix(v(:,:))   ! Initialize v to the identity matrix.
  b(:)=get_diag(a(:,:))      ! Initialize b and d to the diagonal of a. This vector will accumulate terms of
  d(:)=b(:)
  z(:)=0.0_rp                ! This vector will accumulate terms of the form tapq as in eq. (11.1.14).
  nrot=0
  do i=1,50
     sm=sum(abs(a),mask=upper_triangle(n,n))    ! Sum off-diagonal elements.
     if (sm <= toler) RETURN      !The normal return, which relies on quadratic convergence to machine underflow.

     tresh=merge(0.2_rp*sm/n**2,0.0_rp, i < 4 ) ! On the first three sweeps, we will rotate only if tresh exceeded.
     do ipp=1,n-1
        do iq=ipp+1,n
           g=100.0_rp*abs(a(ipp,iq))
           !After four sweeps, skip the rotation if the off-diagonal element is small.
           if ((i > 4) .and. (abs(d(ipp))+g == abs(d(ipp)))  .and. (abs(d(iq))+g == abs(d(iq)))) then
              a(ipp,iq)=0.0
           else if (abs(a(ipp,iq)) > tresh) then
              h=d(iq)-d(ipp)
              if (abs(h)+g == abs(h)) then
                 t=a(ipp,iq)/h         !t = 1/(2tita)
              else
                 theta=0.5_rp*h/a(ipp,iq)
                 t=1.0_rp/(abs(theta)+sqrt(1.0_rp+theta**2))
                 if (theta < 0.0) t=-t
              end if
              c         = 1.0_rp/sqrt(1.0_rp+t**2)
              s         = t*c
              tau       = s/(1.0_rp+c)
              h         = t*a(ipp,iq)
              z(ipp)    = z(ipp)-h
              z(iq)     = z(iq)+h
              d(ipp)    = d(ipp)-h
              d(iq)     = d(iq)+h
              a(ipp,iq) = 0.0_rp

              !              call jrotate(a(1:ipp-1,ipp),a(1:ipp-1,iq),s,tau)
              !Case of rotations 1 <= j < p.
              !              call jrotate(a(ipp,ipp+1:iq-1),a(ipp+1:iq-1,iq),s,tau)
              !Case of rotations p < j < q.
              !              call jrotate(a(ipp,iq+1:n),a(iq,iq+1:n),s,tau)
              !Case of rotations q < j <= n.
              !              call jrotate(v(:,ipp),v(:,iq),s,tau)

              ! I am having problems with jrotate y prefer to do it as in eigenv

              ! case of rotations 1<= J < P

              do J=1, ipp-1
                 G = A(J,ipp)
                 H = A(J,iq)
                 A(J,ipp) = G - S*(H + G*TAU)
                 A(J,iq) = H + S*(G - H*TAU)
              end do

              ! case of rotations P < J < Q

              do J=ipp+1,iq-1
                 G = A(ipp,J)
                 H = A(J,iq)
                 A(ipp,J) = G - S*(H + G*TAU)
                 A(J,iq) = H + S*(G - H*TAU)
              end do

              ! case of rotations Q < J <= N

              do J=iq+1, 3
                 G = A(ipp,J)
                 H = A(iq,J)
                 A(ipp,J) = G - S*(H + G*TAU)
                 A(iq,J) = H + S*(G - H*TAU)
              end do

              do J=1,3
                 G = V(J,ipp)
                 H = V(J,iq)
                 V(J,ipp) = G - S*(H + G*TAU)
                 V(J,iq) = H + S*(G - H*TAU)
              end do

              nrot=nrot+1
           end if
        end do
     end do
     b(:)=b(:)+z(:)
     d(:)=b(:)
     z(:)=0.0_rp
  end do
  messa = 'JACOBI_NRE: TOO MANY ITERATIONS. CALLER SUBRUTINE IS  '//adjustl(trim(callersub))
  call nrerror(trim(messa))

END SUBROUTINE jacobi_nre

SUBROUTINE jrotate(a1,a2,s,tau)
  USE def_kintyp
  REAL(rp), DIMENSION(3), INTENT(INOUT) :: a1,a2
  REAL(rp), INTENT(IN)                  :: s,tau
  REAL(rp), DIMENSION(3) :: wk1

  wk1(:)=a1(:)
  a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
  a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
END SUBROUTINE jrotate
