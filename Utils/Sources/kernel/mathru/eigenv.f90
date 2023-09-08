!------------------------------------------------------------------------
!> @addtogroup Mathematics
!> @{
!> @file    eigenv.f90
!> @author  Herbert Owen
!> @brief   computes all the eigenvalues and eigenvectors of real, symmetrix matrix A (3x3)
!> @details The eigenvalues go into matrix D, while the
!!          eigenvectors go into matrix V. NROT is the number of Jacobi rotations
!!          which were required.
!!          Beware A is altered.
!!          *This subroutine is taken from "Numerical Recipes", page 346
!> @}
!------------------------------------------------------------------------

subroutine eigenv(A,D,V,NROT)

  use def_kintyp
  implicit none
  real(rp),intent(inout)   :: A(3,3)
  real(rp),intent(out)     :: D(3), V(3,3)
  integer(ip),intent(out)  :: NROT

  real(rp)                 :: B(3), Z(3), SM, THRESH
  real(rp)                 :: C, S, G, H, T, TAU, THETA, tolei, tolej, toletole, toler
  integer(ip)              :: J, idime, jdime, i

  toletole = 10000.0_rp
  toler = 1.0e-10           ! Tolerance

  do idime = 1,3
     do jdime = 1,3

        ! store identity matrix in V
        if (idime==jdime) then
           V(idime,jdime) = 1.0_rp
        else
           V(idime,jdime) = 0.0_rp
        end if

     end do

     ! store zero matrix in Z
     Z(idime) = 0.0_rp

     ! initialize B and D to the diagonal of A
     B(idime) = A(idime,idime)
     D(idime) = B(idime)

  end do

  NROT=0

  do i = 1,50

     ! sum off-diagonal elements of A
     SM = ABS(A(1,2)) + ABS(A(1,3)) + ABS(A(2,3))

     ! if sum is 0, then return
     if (SM <= toler) return

     ! in first three sweeps, carry out PQ rotation only if
     ! |A_PQ| > TRESH, a threshold value

     if (i .lt. 4) then
        THRESH = 0.20_rp*SM/3**2
     else
        THRESH = 0.0_rp
     end if

     do idime=1,2
        do jdime=idime+1,3
           G = 100.0_rp*ABS(A(idime,jdime))

           ! after four sweeps, skip rotation if the
           ! off-diagonal element is small
           tolei = G / abs(D(idime))
           tolej = G / abs(D(jdime))
           if ((i.gt.4) .and. (ABS(D(idime))+G .eq. ABS(D(idime))) .and. (ABS(D(jdime))+G .eq. ABS(D(jdime)))) then
!           if ((i.gt.4) .and. (tolei .lt. toletole) .and. (tolej .lt. toletole)) then
!           write(6,*) 'eeeeee', tolei, tolej
              A(idime,jdime)=0.D0
           else if (ABS(A(idime,jdime)) .gt. THRESH) then
              H = D(jdime) - D(idime)
!!!              tole= G / abs(H)
              if (ABS(H)+G .eq. ABS(H)) then
!              if (tole .lt. toletole) then
                 T = A(idime,jdime)/H
              else
                 THETA = 0.5_rp*H/A(idime,jdime)
                 T = 1.0_rp/(ABS(THETA)+SQRT(1.D0+THETA**2))
                 if (THETA .lt. 0.0_rp) then
                    T = -T
                 end if
              end if
              C = 1.0_rp/SQRT(1.0_rp + T**2)
              S = T*C
              TAU = S/(1.0_rp+C)
              H = T*A(idime,jdime)
              Z(idime) = Z(idime) - H
              Z(jdime) = Z(jdime) + H
              D(idime) = D(idime) - H
              D(jdime) = D(jdime) + H
              A(idime,jdime) = 0.0_rp

              ! case of rotations 1<= J < P

              do J=1, idime-1
                 G = A(J,idime)
                 H = A(J,jdime)
                 A(J,idime) = G - S*(H + G*TAU)
                 A(J,jdime) = H + S*(G - H*TAU)
              end do

              ! case of rotations P < J < Q

              do J=idime+1,jdime-1
                 G = A(idime,J)
                 H = A(J,jdime)
                 A(idime,J) = G - S*(H + G*TAU)
                 A(J,jdime) = H + S*(G - H*TAU)
              end do

              ! case of rotations Q < J <= N

              do J=jdime+1, 3
                 G = A(idime,J)
                 H = A(jdime,J)
                 A(idime,J) = G - S*(H + G*TAU)
                 A(jdime,J) = H + S*(G - H*TAU)
              end do

              do J=1,3
                 G = V(J,idime)
                 H = V(J,jdime)
                 V(J,idime) = G - S*(H + G*TAU)
                 V(J,jdime) = H + S*(G - H*TAU)
              end do

              NROT = NROT + 1

           end if
        end do
     end do

     ! update D with the sum of T*A(PQ), and reinitialize Z

     do idime=1,3
        B(idime) = B(idime) + Z(idime)
        D(idime) = B(idime)
        Z(idime) = 0.0_rp
     end do
  end do

end subroutine eigenv

