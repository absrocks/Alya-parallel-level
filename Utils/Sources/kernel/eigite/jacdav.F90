      SUBROUTINE JD(ncomp_eig, nequa_sol, jd_mmax, jd_mmin, amatr, cmass, c_sol, r_sol, v0, eval, evec) 
	  implicit none
      integer          :: ncomp_eig, nequa_sol, jd_mmax, jd_mmin
      integer          ::  c_sol(*), r_sol(*)
      real*8           :: amatr(*), cmass(*), v0(*), eval(*), evec(nequa_sol,*)
      integer          :: i, k, m, m2, off, of2, iaux, l
      real             :: norm
      real,  allocatable :: t(:), bt(:),                                      &
                          MM(:), W(:), VR(:), uu(:), pp(:), ua(:),            &
                          rr(:), Z(:,:), VA(:,:), VB(:,:), VV(:,:),VAUX(:,:)
      real*8, external :: ddot
      real*8, parameter:: ONE=1.0, ZERO=0.0
      real*8           :: jac_eps = 1E-4
      integer          :: MGS_k   = 0.25
      real*8           :: gmr_eps = 1E-6
      integer          :: gmr_dim = 50
      integer          :: gmr_max = 200
      integer          :: jd_it    = 0
      integer          :: itsol_it = 0
	  logical          :: GOOD_EIG


      m2 = (jd_mmax*(jd_mmax+1))/2
!      nullify( t, bt, VA, VB, VV, MM, W, VR, uu, pp, ua, rr, Z, VAUX )
      allocate( t(nequa_sol), bt(nequa_sol) )
      allocate( VA(nequa_sol,jd_mmax), VB(nequa_sol,jd_mmax), VV(nequa_sol,jd_mmax) )
      allocate( MM(m2), W(jd_mmax), VR(jd_mmax*jd_mmax) )
      allocate( uu(nequa_sol), pp(nequa_sol), ua(nequa_sol), rr(nequa_sol) )
      allocate( Z(nequa_sol,ncomp_eig) )
      allocate( VAUX(nequa_sol,jd_mmax) )

      do k=1,nequa_sol
	    t(k) = v0(k)
	  enddo
      
	  k = 0
      m = 0
      l = 0
      do while ( k .lt. ncomp_eig )
        jd_it = jd_it + 1
        l     = l + 1
!        write(2,*) ' '
!        write(2,*) '  NEW JD_ITER:', jd_it, 'ITSOL_IT:', itsol_it,
!     &             'JD_dim:', m+1, '  EV found:', k
!       (3) - (5) Modified Gram Schmidt with refinement

        call MGS_REFINED( nequa_sol, m, t, VV, VB, MGS_k )

!       (6)  m=m+1; bt=B*t; norm = sqrt( t*bt )
        m = m + 1
!        call MATXVEC(ib,jb,b,t,bt,N)

        call bcsrax(0,nequa_sol,1,cmass,c_sol,r_sol,t,bt) 

        norm = sqrt(ddot( nequa_sol, t, 1, bt, 1 ))

!       (7)  VVm=t/norm; VAm=A*VVm;  VBm=bt/norm
        VV(1:nequa_sol,m) = t/norm
        
		call bcsrax(0,nequa_sol,1,amatr,c_sol,r_sol,VV(1,m),VA(1,m)) 


        VB(1:nequa_sol,m) = bt/norm

!       (8) - (10) Update one new row of MM = VV*A*VV
        off = (m-1)*m/2
        do i= 1, m
          MM(off+i) = ddot( nequa_sol, VV(1,i), 1, VA(1,m), 1 )
        enddo

!       (11) Compute eigenvalues and eigenvectors of matrix MM
        call COMPUTE_EIGENPAIRS( m, MM, W, VR )

!       (12) uu=VV*VR(:,1); pp=VB*VR(:,1); ua=VA*VR(:,1); rr=ua-W(1)*pp
        call DGEMV( 'N', nequa_sol, m, ONE, VV, nequa_sol, VR, 1, ZERO, uu, 1 ) ! reemplazar
        call DGEMV( 'N', nequa_sol, m, ONE, VB, nequa_sol, VR, 1, ZERO, pp, 1 ) ! reemplazar
        call DGEMV( 'N', nequa_sol, m, ONE, VA, nequa_sol, VR, 1, ZERO, ua, 1 ) ! reemplazar
        rr = ua - W(1)*pp

!       (13) while ||r||/||u|| < eps
        do while( GOOD_EIG( nequa_sol, rr, uu,jac_eps ) )
          write(2,*) '    GOOD EIGENVALUE(', k+1, ')=', W(1)
!         (14) eval(k+1)=W(1); EVEC=[EVEC,uu]; Z=[Z,pp]; k=k+1
          eval(k+1)     = W(1)
          EVEC(1:nequa_sol,k+1) = uu
          Z(1:nequa_sol,k+1)    = pp
          k             = k + 1
!         (15) Exit if we have enought eigenvalues/eigenvectors
          if (k.eq.ncomp_eig) goto 100

!         (16) m = m - 1; M = 0
          m         = m - 1
          off       = (m+1)*m/2
          if (m.gt.0) then
!           (18) VV_=VV*VR(:,2:m+1); VA_=VA*VR(:,2:m+1); VB_=VB*VR(:,2:m+1)
            MM(1:off) = 0.0
            call DGEMM( 'N', 'N', nequa_sol, m, m+1, ONE, VV, nequa_sol, VR(m+2), m+1,  &
                       ZERO, VAUX, nequa_sol )  ! reemplazar
            VV(1:nequa_sol,1:m) = VAUX(1:nequa_sol,1:m)
            call DGEMM( 'N', 'N', nequa_sol, m, m+1, ONE, VA, nequa_sol, VR(m+2), m+1,  &
                       ZERO, VAUX, nequa_sol )  ! reemplazar
            VA(1:nequa_sol,1:m) = VAUX(1:nequa_sol,1:m)
            call DGEMM( 'N', 'N', nequa_sol, m, m+1, ONE, VB, nequa_sol, VR(m+2), m+1,  &
                       ZERO, VAUX, nequa_sol )  ! reemplazar
            VB(1:nequa_sol,1:m) = VAUX(1:nequa_sol,1:m)
          endif

!         (19) M(i,i)=W(i+1); VR(i)=ei; W(i)=W(i+1)
          off = 0
          of2 = 0
          VR(1:m*m) = 0.0
          do i= 1, m
            off       = off + i
            MM(off)   = W(i+1)
            VR(of2+i) = 1.0
            W(i)      = W(i+1)
            of2       = of2 + m
          enddo

!         (21) uu=VV(1:n,1); pp=VB(1:n,1); rr=VA(1:n,1)-W(1)*pp
          uu = VV(1:nequa_sol,1)
          pp = VB(1:nequa_sol,1)
!!!          ua = VA(1:n,1)    <<<< ¿To Do or not to do?
          rr = VA(1:nequa_sol,1) - W(1)*pp

!         Reset l
          l = 1
        end do

        if (m.eq.jd_mmax) then
        write(2,*) '  Restart Krylov Space'
!         (24) M = 0
          off       = ((jd_mmin+1)*jd_mmin)/2
          MM(1:off) = ZERO

!         (25) - (28) VV = [uu,VV*VR(2:jd_mmin)]
          call DGEMM( 'N', 'N', nequa_sol, jd_mmin-1, m, ONE, VV, nequa_sol, VR(m+1),   &
                     m+1, ZERO, VAUX, nequa_sol )  !reemplazar

          VV(1:nequa_sol,1)           = uu
          VV(1:nequa_sol,2:jd_mmin) = VAUX(1:nequa_sol,1:jd_mmin-1)

!         VA = [ua,VA*VR(2:jd_mmin)]
          call DGEMM( 'N', 'N', nequa_sol, jd_mmin-1, m, ONE, VA, nequa_sol, VR(m+1),    &
                     m+1, ZERO, VAUX, nequa_sol )   !reemplazar

          VA(1:nequa_sol,1)           = ua
          VA(1:nequa_sol,2:jd_mmin) = VAUX(1:nequa_sol,1:jd_mmin-1)

!         VB = [pp,VB*VR(2:jd_mmin)]
          call DGEMM( 'N', 'N', nequa_sol, jd_mmin-1, m, ONE, VB, nequa_sol, VR(m+1),     &
                      m+1, ZERO, VAUX, nequa_sol )  !reemplazar

          VB(1:nequa_sol,1)           = pp
          VB(1:nequa_sol,2:jd_mmin) = VAUX(1:nequa_sol,1:jd_mmin-1)

!         M(i,i)=W(i+1)
          m   = jd_mmin
          off = 0
          do i= 1, m
            off     = off + i
            MM(off) = W(i)
          enddo
        endif

!       (30) Z=[Z,p]
        Z(1:nequa_sol,k+1) = pp

!       (31) JD Correction Equation
        call JD_CE( nequa_sol, k+1, Z, amatr, cmass, c_sol, r_sol, W(1), rr, t ,gmr_dim,gmr_max, itsol_it,gmr_eps )

      enddo

100   continue

      deallocate(t, bt, VA, VB, VV, MM, W, VR, uu, pp, ua, rr, Z, VAUX)

      END SUBROUTINE JD




      SUBROUTINE JD_CE( n, k, Z, amatr, cmass, c_sol, r_sol, teta, rr, t,gmr_dim,gmr_max, itsol_it,gmr_eps )
      implicit none
      integer,    intent(in) :: n, k, c_sol(*), r_sol(*),itsol_it
      real*8,     intent(in) :: Z(n,k), amatr(*), cmass(*), teta, rr(n)
      real*8,    intent(out) :: t(n)
      target                 :: Z, rr

      integer                :: info, nnz, nzmax, ierr, fil, i, j, jj, l,gmr_dim,gmr_max
      real*8                 :: tol, alfa,gmr_eps 
      integer,       pointer :: ipiv(:)
      real*8,        pointer :: Z_(:,:), M(:,:), r_(:), gama(:),alpha(:)
      real*8,      parameter :: ONE=1.0, ZERO=0.0, M_ONE=-1.0
      nullify( M, ipiv, gama, alpha )
      allocate( M(k,k), ipiv(k), gama(k), alpha(k) )

!     (1) K*Z_ => Z

      Z_ => Z(1:n,1:k)

!     (2) M = Z^(*)*Z_
      call DGEMM( 'T', 'N', k, k, n, ONE, Z, n, Z_, n, ZERO, M, k ) ! reemplzar

      ! (3)
!C     GETRF computes an LU factorization of a general M-by-N matrix A
!C     using partial pivoting with row interchanges.

      call DGETRF( k, k, M, k, ipiv, info )

      if (info.ne.0) STOP 'ERROR IN DGETRF'

!     (5) solve K r_ = r
      r_ => rr(1:n)

      ! (6) gama = Z r
      call DGEMV( 'T', n, k, ONE, Z, n, r_, 1, ZERO, gama, 1 )

      ! (7) (8)
!C     ZGETRS solves a system of linear equations
!C     A * X = B,  A**T * X = B,  or  A**H * X = B
!C     with a general N-by-N matrix A using the LU factorization 
!C     computed by DGETRF.

      alpha = gama
      call DGETRS( 'N', k, 1, M, k, ipiv, alpha, k, info )
      if (info.ne.0) STOP 'ERROR IN DGETRS'

      ! (9) r_ = Z_*alpha - r_
      call DGEMV( 'N', n, k, ONE, Z_, n, alpha, 1, M_ONE, r_, 1 )


      gmr_dim = n
      gmr_max = n
      call JD_GMRES( n, k, amatr, cmass, c_sol, r_sol, teta, r_, Z, Z_, M, ipiv, t ,gmr_max, itsol_it,gmr_dim,gmr_eps)

      deallocate( M, ipiv, gama, alpha )
      END SUBROUTINE JD_CE

      SUBROUTINE MGS_REFINED( n, m, t, VV, VB,MGS_k )
      implicit none
      integer,    intent(in) :: n, m,MGS_k
      real*8,     intent(in) :: VV(n,*), VB(n,*)
      real*8,  intent(inout) :: t(n)

      integer                :: i, j
      real*8                 :: taut_ini, taut_end, alpha
      real*8,       external :: ddot

      if (m.gt.0) then
        taut_ini = sqrt(ddot( n, t, 1, t, 1 ))
        do i= 1, m
          alpha = ddot( n, t, 1, VB(1,i), 1 )
          t     = t - alpha*VV(1:n,i)
        enddo
        taut_end = sqrt(ddot( n, t, 1, t, 1 ))
        if (taut_end/taut_ini .le. MGS_k) then
          do i= 1, m
            alpha = ddot( n, t, 1, VB(1,i), 1 )
            t     = t - alpha*VV(1:n,i)
          enddo
        endif
      endif

      END SUBROUTINE MGS_REFINED

      SUBROUTINE COMPUTE_EIGENPAIRS( m, VAV, W, VR )
      implicit none
!C     Input parameters
      integer,    intent(in) :: m
      real*8                 :: VAV((m*m+m)/2)
!C     Output parameters
      real*8                 :: W(m), VR(m*m)
!C     Local Variables
      integer                :: info
      real*8,        pointer :: A(:), work(:)

      ! LAPACK
      !
      ! DSPEV computes all the eigenvalues and, optionally, eigenvectors
      ! of a complex Hermitian matrix in packed storage
      !
      ! The right eigenvector v(j) of A satisfies
      !                  A * v(j) = lambda(j) * v(j)
      ! where lambda(j) is its eigenvalue.
      ! The left eigenvector u(j) of A satisfies
      !               u(j)**H * A = lambda(j) * u(j)**H
      ! where u(j)**H denotes the conjugate transpose of u(j).
      !
      ! The computed eigenvectors are normalized to have Euclidean norm
      ! equal to 1 and largest component real.
      
      nullify( A, work )
      allocate( A((m*m+m)/2), work(3*m) )
      A = VAV

      call DSPEV( 'V', 'U', m, A, W, VR, m, work, info )

      if (info.ne.0) STOP 'ERROR in compute_eigenpairs->DSPEV'
      write(2,*)'      CANDIDATE EIGENVALUES:', W
      deallocate( work, A )

      END SUBROUTINE COMPUTE_EIGENPAIRS

      logical FUNCTION GOOD_EIG( n, x, y, jac_eps )
      implicit none
      integer,     intent(in) :: n
      real*8,      intent(in) :: x(n), y(n),jac_eps
!      logical                 :: GOOD_EIG
      real*8                  :: rx, ry
      real*8,        external :: ddot

      rx = ddot( n, x, 1, x, 1 )
      ry = ddot( n, y, 1, y, 1 )
      GOOD_EIG = sqrt(rx/ry) .le. jac_eps
      
	  END FUNCTION GOOD_EIG



      subroutine JD_GMRES( n, k, amatr, cmass, c_sol, r_sol, teta, r, Z, Z_, M, ipiv, t ,gmr_max,itsol_it,gmr_dim,gmr_eps)
      implicit none

      integer,     intent(in) :: n, k, c_sol(*), r_sol(*), ipiv(k)
      real*8,      intent(in) :: amatr(*), cmass(*), teta, r(n), Z(n,k),Z_(n,k), M(k,k)

      real*8,     intent(out) :: t(n)

      logical                 :: conv, fin2
      integer                 :: iter, ii, jj, jj1, kk, idx, dim2,gmr_max,itsol_it,gmr_dim
      real*8                  :: raux, stopcri, gamma, tmp, gmr_eps
      real*8,         pointer :: kryl(:,:), rs(:), hh(:),cc(:), ss(:)
      real*8,        external :: ddot

      conv = .FALSE.
      t    = 0.0
      raux = SQRT(ddot( n, r, 1, r, 1 ))
      tmp  = raux
      if (raux.eq.0.0) return
      stopcri = gmr_eps*raux
      dim2    = gmr_dim*(gmr_dim+3)/2

      nullify( kryl, rs, hh, cc, ss )
      allocate( kryl(n,gmr_dim+1), rs(gmr_dim+1), hh(dim2), cc(gmr_dim), ss(gmr_dim) )


      iter = 0
      DO WHILE (.NOT. conv)
        
		call JD_MV( n, k,amatr, cmass, c_sol, r_sol, teta, t, Z, Z_, M, ipiv, kryl(:,1) )

        kryl(1:n,1) = r - kryl(1:n,1)

!C       raux = ||kryl(*,1)||
        raux = SQRT(ddot( n, kryl(1,1), 1, kryl(1,1), 1 ))
!c        write(2,*) '      GMRES raux:', iter, raux/tmp

        if (raux.le.stopcri) then
!C         The initial guess is the solution
          conv = .TRUE.
        else
!         Initialize 1-st term of the rhs of hessenberg system
          rs(1) = raux
!         Ortonormalize kryl(*,1)
          raux = 1.0/raux
          kryl(1:n,1) = kryl(1:n,1)*raux
        endif

        jj   = 0
        idx  = -1
        FIN2 = conv

!C       Inner loop. Restarted each kryldim iterations.
        DO WHILE (.NOT. FIN2)
          itsol_it = itsol_it + 1
          iter = iter + 1
          jj   = jj + 1
          jj1  = jj + 1
          idx  = idx + jj

          call JD_MV( n, k, amatr, cmass, c_sol, r_sol, teta, kryl(:,jj),     &
                     Z, Z_, M, ipiv, kryl(:,jj1) )
!C         Modified Gram-Schmidt
!C         For i= 1, j
!C           H(i,j) = <v_i, v_j1>
!C           v_j1   = v_j1 - H(i,j) * v_i
          DO ii= 1, jj
            raux = ddot( n, kryl(1:n,ii), 1, kryl(1:n,jj1), 1 )
            hh(idx+ii) = raux
            kryl(1:n,jj1) = kryl(1:n,jj1) - raux*kryl(1:n,ii)
          ENDDO

!C         H(jj1,jj) = ||kryl(*,jj1)||
          raux = SQRT(ddot( n, kryl(1:n,jj1), 1, kryl(1:n,jj1), 1 ))
          hh(idx+jj1) = raux

          if (raux.eq.0.0) then
            FIN2 = .TRUE.
            conv = .TRUE.
            idx  = idx - jj
            jj   = jj - 1
            iter = iter - 1
          else
!C           Ortonormalize kryl(*,jj1)
            raux = 1.0/raux
            kryl(1:n,jj1) = kryl(1:n,jj1)*raux

!C           Update factorization of H. Perform previous 
!C           transformations on jj-th column of H
            do ii= 1, jj-1
              kk   = ii + 1
              raux = hh(idx+ii)

              hh(idx+ii) =  cc(ii)*raux + ss(ii)*hh(idx+kk)
              hh(idx+kk) = -ss(ii)*raux + cc(ii)*hh(idx+kk)
            enddo

            gamma = hh(idx+jj)*hh(idx+jj) + hh(idx+jj1)*hh(idx+jj1)
            gamma = SQRT(gamma)
     
!C           if gamma is zero then take any small
!C           value will affect only residual estimate
            if (gamma.eq.0.0) gamma = gmr_eps

!C           Get next plane rotation
            cc(jj) = hh(idx+jj)/gamma
            ss(jj) = hh(idx+jj1)/gamma

            hh(idx+jj) = cc(jj)*hh(idx+jj) + ss(jj)*hh(idx+jj1)

!C           Update the rhs of the LS problem
            rs(jj1) = -ss(jj)*rs(jj)
            rs(jj)  =  cc(jj)*rs(jj)

!C           Convergence Test
            raux = ABS( rs(jj1) )
!c            write(2,*) '      GMRES raux:', iter, raux/tmp
            if (raux.le.stopcri) then
!              write(2,*) '      GMRES CONVERGED!!!!', iter
              conv = .TRUE.
              FIN2 = .TRUE.
            else
              if (iter.ge.gmr_max) then
!                write(2,*) '         GMRES ERROR: MAXITER REACHED'
                conv = .TRUE.
                FIN2 = .TRUE.
              else if (jj .ge. gmr_dim) then
                FIN2 = .TRUE.
              endif
            endif
          endif
        END DO   ! Krylov Loop End
        if (jj .gt. 0) then
!C         Compute y => Solve upper triangular system
          do ii= jj, 2, -1
            rs(ii) = rs(ii)/hh(idx+ii)
            raux   = rs(ii)

            do kk= 1, ii-1
              rs(kk) = rs(kk) - hh(idx+kk) * raux
            enddo

            idx = idx - ii
          enddo
          rs(1) = rs(1) / hh(1)

!C         Linear combination of kryl(*,jj)'s to get the solution.
          do ii= 1, jj
            raux = rs(ii)
            t = t + raux*kryl(1:n,ii)
          enddo
        endif

      ENDDO

      deallocate( kryl, rs, hh, cc, ss )
      end subroutine JD_GMRES

      subroutine JD_MV( n, k, amatr, cmass, c_sol, r_sol,teta, v, Z, Z_, M, ipiv, zz )
      implicit none

      integer,         intent(in) :: n, k, ipiv(k), c_sol(*), r_sol(*)
      real*8,          intent(in) :: amatr(*), cmass(*), teta, v(n), Z(n,k), Z_(n,k), M(k,k)

      real*8,         intent(out) :: zz(n)

      integer                     :: info
      real*8,             pointer :: y(:), y_(:), gamma(:), wa(:)
      real*8,           parameter :: ONE=1.0, ZERO=0.0, M_ONE=-1.0

      nullify( y, y_, gamma, wa )
      allocate( y(n), y_(n), gamma(k), wa(n) )

      ! (11)
      call bcsrax(0,n,1,amatr,c_sol,r_sol,v,y) 
      call bcsrax(0,n,1,cmass,c_sol,r_sol,v,y_) 

!      call MATXVEC(ia,ja,a,v,y,N)
!      call MATXVEC(ib,jb,b,v,y_,N)

      y = y - teta*y_

      ! (12) K*zz = y
      zz = y

      ! (13)
      call DGEMV( 'T', n, k, ONE, Z, n, y, 1, ZERO, gamma, 1 )

      ! (14) (15)
!C     DGETRS solves a system of linear equations
!C         A * X = B,  A**T * X = B,  or  A**H * X = B
!C     with a general N-by-N matrix A using the LU factorization
!C     computed by DGETRF.
      call DGETRS ( 'N', k, 1, M, k, ipiv, gamma, k, info )
      if (info.ne.0) STOP 'ERROR in jd_mv->DGETRS'

      ! (16)
      call DGEMV( 'N', n, k, M_ONE, Z_, n, gamma, 1, ONE, zz, 1 )

      deallocate( y, y_, gamma, wa )
      end subroutine JD_MV



