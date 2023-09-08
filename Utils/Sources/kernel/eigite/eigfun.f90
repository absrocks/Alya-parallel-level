real(rp) function eigdot(npoin,X,Y)
  use def_kintyp, only :  ip,rp
  implicit none
  integer(ip)         :: npoin
  integer(ip)         :: kk
  real(rp)            :: X(*),Y(*)

  eigdot = 0.0_rp

  DO KK=1,npoin
     EIGDOT=EIGDOT + X(KK)*Y(KK)
  ENDDO

end function eigDOT

SUBROUTINE  ORTOGONALIZA(NP,AUTOVEC,MM) 
  use def_kintyp, only  :  ip,rp
  use def_master, only  :  npoi1,npoi2,npoi3,parre,nparr,&
       &                   INOTMASTER,IMASTER
  implicit none
  integer(ip)            :: NP,MM,JJ,KK
  real(rp)               :: AUTOVEC(NP,*),NORMA1,NORMA2
  real(rp),  ALLOCATABLE :: X(:)
  real(rp)               :: dummr

  if( INOTMASTER ) then

     ALLOCATE(X(NP))
     DO JJ=1,NP
        X(JJ) = AUTOVEC(JJ,MM)
     ENDDO

     DO KK=1,MM-1

        call prodts(1_ip,np,AUTOVEC(1,KK),AUTOVEC(1,MM),norma1,norma2)

        !XMAX = 0.0_rp
        DO JJ=1,NP
           X(JJ)=X(JJ) - NORMA2* AUTOVEC(JJ,KK) / NORMA1
           !IF(XMAX.LT.ABS(X(JJ))) XMAX = ABS(X(JJ)) 
        ENDDO
     ENDDO

     DO JJ=1,NP
        AUTOVEC(JJ,MM) = X(JJ)
     ENDDO

     DEALLOCATE(X)

  else

     DO KK=1,MM-1        
        call prodts(1_ip,np,dummr,dummr,norma1,norma2)
     end DO

  end if

  RETURN
END SUBROUTINE ORTOGONALIZA


SUBROUTINE  ORTOGeig3(np,MM,napunte,AUTOVEC) 
  use def_kintyp, only  :  ip,rp
  use def_master, only  :  npoi1,npoi2,npoi3,parre,nparr,&
       &                   INOTMASTER,IMASTER
  implicit none
  integer(ip)            :: NP,MM,JJ,KK,II
  integer(ip)            :: napunte(MM)
  real(rp)               :: AUTOVEC(NP,*),NORMA1,NORMA2
  real(rp),  ALLOCATABLE :: X(:)
  real(rp)               :: dummr

  if( INOTMASTER ) then

     ALLOCATE(X(NP))

     do II = 2,MM

        IF(NAPUNTE(II).EQ.0) THEN

           DO JJ = 1,NP
              X(JJ) = AUTOVEC(JJ,II)
           ENDDO

           DO KK = 1,II-1

              call prodts(1_ip,np,AUTOVEC(1,KK),AUTOVEC(1,II),norma1,norma2)

              DO JJ = 1,NP
                 X(JJ) = X(JJ) - NORMA2 * AUTOVEC(JJ,KK) / NORMA1
              ENDDO
           ENDDO

           DO JJ = 1,NP
              AUTOVEC(JJ,II) = X(JJ)
           ENDDO

        ENDIF

     enddo
     DEALLOCATE(X)

  else

     do II = 2,MM

        IF(NAPUNTE(II).EQ.0) THEN


           DO KK = 1,II-1

              call prodts(1_ip,np,dummr,dummr,norma1,norma2)

           ENDDO


        ENDIF

     enddo

  end if

  RETURN
END SUBROUTINE ORTOGEIG3


!    SUBROUTINE  ortogGS(MM,napunte) 
!    implicit none

!    INTEGER  MM,ii,jj,kk
! INTEGER NAPUNTE(mm)
! DOUBLE PRECISION NORMA1,NORMA2,xmax
!    DOUBLE PRECISION, ALLOCATABLE :: X(:)

! ALLOCATE(X(nbnodes))

!      do ii=2,mm

!     IF(NAPUNTE(II).EQ.0) THEN
!          DO JJ=1,nbnodes
!         X(JJ) = eigen(JJ,ii)
!       ENDDO

!    DO KK=1,ii-1
!   NORMA1 = 0.0
!   NORMA2 = 0.0
!   DO JJ=1,nbnodes
!     NORMA1 = NORMA1 + eigen(JJ,KK)**2
!     NORMA2 = NORMA2 + eigen(JJ,KK)*eigen(JJ,ii)
!   ENDDO
!   XMAX = 0
!   DO JJ=1,nbnodes
!     X(JJ)=X(JJ) - NORMA2* eigen(JJ,KK) / NORMA1
!     IF(XMAX.LT.ABS(X(JJ))) XMAX = ABS(X(JJ)) 
!   ENDDO
!    ENDDO
!
!    DO JJ=1,nbnodes
!   eigen(JJ,ii) = X(JJ)
!    ENDDO
!       ENDIF

! enddo

!      DEALLOCATE(X)

!      RETURN
! END




subroutine NORMALIZA_AM(npoin,keig,kfl_massm,eigen,bmatr,c_sym,r_sym)
  use def_kintyp, only  :  ip,rp
  use def_master, only  :  kfl_paral,npoi1,npoi2,npoi3,parre,nparr,&
       &                   INOTMASTER,IMASTER
  use def_solver, only  :  solve_sol
  implicit none
  integer(ip)         :: npoin,keig,kfl_massm
  real(rp)            :: eigen(npoin,*)
  real(rp)            :: bmatr(*)
  integer(ip)         :: c_sym(*),r_sym(*)
  integer(ip)         :: kk
  real(rp)            :: norma,dummr
  real(rp),   pointer :: xt(:),x(:) 

  if( kfl_massm == 0 ) then
     !
     ! Diagonal matrix
     !
     if( IMASTER ) then
        call proxyz(1_ip,npoin,dummr,dummr,dummr,norma)
     else
        call proxyz(1_ip,npoin,eigen(1,keig),eigen(1,keig),bmatr,norma)
        norma = 1.0_rp / sqrt(abs(norma))
        do kk = 1,npoin
           eigen(kk,keig) = eigen(kk,keig) * norma
        end do
     end if

  else
     !
     ! Sparse matrix
     !
     if( INOTMASTER ) then
        allocate (xt(npoin),x(npoin))
        do kk = 1,npoin
           x(kk) = eigen(kk,keig)
        end do
        if( solve_sol(1)%kfl_symme == 1 ) then
           call bsymax(1_ip,npoin,1_ip,bmatr,c_sym,r_sym,x,xt)
        else
           call bcsrax(1_ip,npoin,1_ip,bmatr,c_sym,r_sym,x,xt)
        end if
        call prodxy(1_ip,npoin,eigen(1,keig),xt,norma)
        norma = 1.0_rp/sqrt(abs(norma))
        do kk=1,npoin
           eigen(kk,keig)= eigen(kk,keig) * norma
        end do
        deallocate(xt,x)
     else
        call prodxy(1_ip,npoin,dummr,dummr,norma)
     end if

  end if

end subroutine normaLIZA_AM

subroutine orto_AM(npoin,keig,kfl_massm,eigen,x,bmatr,c_sym,r_sym) 
  use def_kintyp, only  :  ip,rp
  use def_master, only  :  IMASTER,INOTMASTER
  use def_solver, only  :  solve_sol
  implicit none
  integer(ip)           :: npoin,keig,kfl_massm
  real(rp)              :: eigen(npoin,*),x(npoin)
  real(rp)              :: bmatr(*)
  integer(ip)           :: c_sym(*),r_sym(*)
  integer(ip)           :: jj,kk
  real(rp)              :: norma2,dummr
  real(rp),   pointer   :: xt(:)

  if( kfl_massm == 0 ) then
     !
     ! Diagonal matrix
     !
     do kk = 1,keig-1

        if( IMASTER ) then
           call proxyz(1_ip,npoin,dummr,dummr,dummr,norma2)
        else
           call proxyz(1_ip,npoin,eigen(1,kk),bmatr,x,norma2)
           do jj = 1,npoin
              x(jj) = x(jj) - norma2 * eigen(jj,kk)
           end do
        end if

     end do

  else
     !
     ! Sparse matrix
     !
     if( INOTMASTER ) then
        allocate(xt(npoin))
        if( solve_sol(1)%kfl_symme == 1 ) then
           call bsymax(1_ip,npoin,1_ip,bmatr,c_sym,r_sym,x,xt)
        else
           call bcsrax(1_ip,npoin,1_ip,bmatr,c_sym,r_sym,x,xt)
        end if
        do kk = 1,keig-1
           call prodxy(1_ip,npoin,eigen(1,kk),xt,norma2)
           do jj = 1,npoin
              x(jj) = x(jj) - norma2 * eigen(jj,kk) 
           end do
        end do
        deallocate(xt)
     else
        do kk = 1,keig-1
           call prodxy(1_ip,npoin,dummr,dummr,norma2)
        end do
     end if

  end if

end subroutine orto_AM

real(rp) function funmaxi(npoin,eigen,keig)
  use def_kintyp, only :  ip,rp
  use def_master, only :  INOTMASTER,IPARALL,parre,nparr
  implicit none
  integer(ip)          :: npoin,keig
  real(rp)             :: eigen(npoin,*)
  real(rp),    target  :: dummr(1) 
  integer(ip)          :: kk

  funmaxi = 0.0_rp

  if( INOTMASTER ) then
     do kk = 1,npoin
        if( funmaxi < ABS(eigen(kk,keig)) ) THEN
           funmaxi = ABS(eigen(kk,keig))
        end if
     end do
  end if
  if( IPARALL ) then
     nparr    =  1
     dummr(1) =  funmaxi
     parre    => dummr
     call Parall(10_ip)           
     funmaxi  =  dummr(1)
  end if

end function funmaxi

subroutine ortogGS_BIS(npoin,keig,kfl_massm,eigen,bmatr,c_sym,r_sym)
  use def_kintyp, only  :  ip,rp
  use def_solver, only  :  solve_sol
  use def_master, only  :  INOTMASTER
  implicit none 
  integer(ip)           :: npoin,keig,kfl_massm
  real(rp)              :: eigen(npoin,*)
  real(rp)              :: bmatr(*)
  integer(ip)           :: c_sym(*),r_sym(*)
  integer(ip)           :: jj,kk
  real(rp), ALLOCATABLE :: x(:),xt(:)
  real(rp)              :: norma2,dummr

  if( kfl_massm == 0 ) then

     if( INOTMASTER ) then
        allocate (xt(npoin))
        do kk = 1,npoin
           xt(kk) = bmatr(kk) * eigen(kk,keig)
        end do
        do kk = 1,keig-1
           call prodxy(1_ip,npoin,eigen(1,kk),xt,norma2)
           do jj = 1,npoin
              eigen(jj,keig) = eigen(jj,keig) - norma2 * eigen(jj,keig)
           end do
        end do
        deallocate(xt)
     else
        do kk = 1,keig-1
           call prodxy(1_ip,npoin,dummr,dummr,norma2)
        end do
     end if

  else

     if( INOTMASTER ) then
        allocate (xt(npoin),x(npoin))
        do kk = 1,npoin
           x(kk) = eigen(kk,keig)
        enddo
        if( solve_sol(1)%kfl_symme == 1 ) then
           call bsymax(1_ip,npoin,1_ip,bmatr,c_sym,r_sym,x,xt)
        else
           call bcsrax(1_ip,npoin,1_ip,bmatr,c_sym,r_sym,x,xt)
        end if
        do kk = 1,keig-1
           call prodxy(1_ip,npoin,eigen(1,kk),xt,norma2)
           do jj=1,npoin
              eigen(jj,keig)=eigen(jj,keig) - norma2 * eigen(jj,keig) !/norma2
           end do
        end do
        deallocate(xt,x)
     else
        do kk = 1,keig-1
           call prodxy(1_ip,npoin,dummr,dummr,norma2)
        end do
     end if

  end if

end subroutine ortogGS_BIS

real(rp) function funmaxmin(npoin,x)
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,IPARALL,parre,nparr,kfl_paral
  use mod_communications, only : PAR_MIN,PAR_MAX
  implicit none
  integer(ip)          :: npoin
  real(rp)             :: x(npoin)
  real(rp)             :: minva,maxva
  integer(ip)          :: kk

  if( INOTMASTER ) then

     minva =  1.0e6_rp
     maxva = -1.0e6_rp
     do kk = 1,npoin
        if( x(kk) < minva ) THEN
           minva = x(kk)
        end if
        if( x(kk) > maxva ) THEN
           maxva = x(kk)
        end if
     end do

  end if

  if( IPARALL ) then
     call PAR_MAX(maxva)
     call PAR_MIN(minva)
  end if

  if( abs(minva) > abs(maxva) ) then
     funmaxmin = minva
  else
     funmaxmin = maxva
  end if

end function funmaxmin

subroutine eigup1(nbnodes,keig,ALANDA,y,x,abstol,autovec)
  use def_kintyp, only      :  ip,rp
  use def_master, only      :  IPARALL,ISLAVE,ISEQUEN,&
       &                       parre,nparr,npoi1,npoi2,npoi3
  implicit none
  integer(ip)               :: nbnodes,keig
  real(rp),   intent(in)    :: ALANDA,y(*)
  real(rp),   intent(inout) :: x(*)
  real(rp),   intent(out)   :: abstol,autovec(nbnodes,*)
  integer(ip)               :: kk
  real(rp)                  :: ax,diff,denom,xn
  real(rp),   target        :: dummr(2)

  if( ISEQUEN ) then

     ABSTOL = 0.0_rp
     DENOM  = 0.0_rp     
     do KK = 1,nbnodes
        AX               = X(KK)
        X(KK)            = Y(KK)  * ALANDA 
        diff             = X(KK)  - AX
        ABSTOL           = ABSTOL + diff  * diff
        DENOM            = DENOM  + X(KK) * X(KK)
        autovec(KK,KEIG) = Y(KK)
     end do

  else if( ISLAVE ) then

     dummr(1) = 0.0_rp
     dummr(2) = 0.0_rp     
     do KK = 1,npoi1
        xn               = Y(KK)  * ALANDA
        diff             = xn - X(KK)
        dummr(1)         = dummr(1) + diff * diff
        dummr(2)         = dummr(2) + xn   * xn
     end do
     do KK = npoi2,npoi3
        xn               = Y(KK)  * ALANDA
        diff             = xn - X(KK)
        dummr(1)         = dummr(1) + diff * diff
        dummr(2)         = dummr(2) + xn   * xn
     end do

     do kk = 1,nbnodes
        X(KK)            = Y(KK)  * ALANDA 
        autovec(KK,KEIG) = Y(KK)
     end do

  end if

  if( IPARALL ) then
     nparr =  2
     parre => dummr
     call Parall(9_ip)           
     ABSTOL = dummr(1)
     DENOM  = dummr(2)
  end if

  if( DENOM == 0.0_rp ) DENOM = 1.0_rp
  ABSTOL = SQRT(ABSTOL/DENOM)

end subroutine eigup1
