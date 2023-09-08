subroutine bcsrai(itask,nbnodes,nbvar,an,ja,ia,xx,yy)
  !----------------------------------------------------------------------
  !****f* mathru/bcsrai
  ! NAME
  !     bcsrai
  ! DESCRIPTION
  !     Multiply a non symmetric matrix stored in BCSR by a vector
  !     YY = A XX
  ! INPUT
  !    NBNODES .... Number of equations
  !    NBVAR ...... Number of variables
  !    AN ......... Matrix
  !    JA ......... List of elements
  !    IA ......... Pointer to list of elements
  !    XX ......... Vector
  ! OUTPUT
  !    YY ......... result vector
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp,         only     :  ip,rp
  use def_master,         only     :  INOTMASTER,npoi1,IPARALL
  use def_master,         only     :  NPOIN_TYPE,kfl_paral
  use def_solver,         only     :  solve_sol
  use mod_communications, only     :  PAR_INTERFACE_OWN_NODE_EXCHANGE
  implicit none
  integer(ip), intent(in)          :: itask,nbnodes,nbvar
  real(rp),    intent(in)          :: an(nbvar,nbvar,*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  real(rp),    intent(inout)       :: xx(nbvar,*)
  real(rp),    intent(out), target :: yy(nbvar,*)
  integer(ip)                      :: ii,jj,kk,ll,col
  real(rp)                         :: raux,raux1,raux2,raux3
  real(rp)                         :: time1,time2,time3,time4,time5,time_spmv

  if( INOTMASTER ) then

     call cputim(time1)

     if( solve_sol(1) % kfl_full_rows == 1 ) then

        !-------------------------------------------------------------------
        !
        ! Full row
        !
        !-------------------------------------------------------------------

        call cputim(time2)
        call PAR_INTERFACE_OWN_NODE_EXCHANGE(nbvar,xx,'SEND RECEIVE')
        call cputim(time3)
        !
        ! Interior nodes (square matrix)
        !
        if( nbvar == 1 ) then
           !
           ! NBVAR=1
           !
           !$OMP PARALLEL  DO                                     &
           !$OMP SCHEDULE ( STATIC )                              &
           !$OMP DEFAULT  ( NONE )                                &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, raux )
           !
           do ii = 1,nbnodes
              yy(1,ii) = 0.0_rp
              do jj   = ia(ii),ia(ii+1)-1
                 col  = ja(jj)
                 if( col <= npoi1 ) then
                    raux = xx(1,col)
                    yy(1,ii) = yy(1,ii) +an(1,1,jj) * raux
                 end if
              end do
           end do

        else if( nbvar == 2 ) then
           !
           ! NBVAR=2
           !
           !$OMP PARALLEL  DO                                     &
           !$OMP SCHEDULE ( STATIC )                              &
           !$OMP DEFAULT  ( NONE )                                &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
           !
           do ii = 1,nbnodes
              yy(1,ii) = 0.0_rp
              yy(2,ii) = 0.0_rp
              do jj       = ia(ii),ia(ii+1)-1
                 col      = ja(jj)
                 if( col <= npoi1 ) then
                    raux1    = xx(1,col)
                    raux2    = xx(2,col)
                    yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
                    yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
                    yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
                    yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
                 end if
              end do
           end do

        else if( nbvar == 3 ) then
           !
           ! NBVAR=3
           !
           !$OMP PARALLEL  DO                                       &
           !$OMP SCHEDULE ( STATIC )                                &
           !$OMP DEFAULT  ( NONE )                                  &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )    &
           !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
           !
           do ii = 1,nbnodes
              yy(1,ii) = 0.0_rp
              yy(2,ii) = 0.0_rp
              yy(3,ii) = 0.0_rp
              do jj       = ia(ii),ia(ii+1)-1
                 col      = ja(jj)
                 if( col <= npoi1 ) then
                    raux1    = xx(1,col)
                    raux2    = xx(2,col)
                    raux3    = xx(3,col)
                    yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
                    yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
                    yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux3
                    yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
                    yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
                    yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux3
                    yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux1
                    yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux2
                    yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux3
                 end if
              end do

           end do

        else if( nbvar == 4 ) then
           !
           ! NBVAR=4
           !
           !$OMP PARALLEL  DO                                     &
           !$OMP SCHEDULE ( STATIC )                              &
           !$OMP DEFAULT  ( NONE )                                &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, raux )
           !
           do ii = 1,nbnodes
              yy(1,ii) = 0.0_rp
              yy(2,ii) = 0.0_rp
              yy(3,ii) = 0.0_rp
              yy(4,ii) = 0.0_rp
              do jj       = ia(ii),ia(ii+1)-1
                 col      = ja(jj)
                 if( col <= npoi1 ) then
                    raux     = xx(1,col)
                    yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
                    yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux
                    yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux
                    yy(4,ii) = yy(4,ii) + an(1,4,jj) * raux
                    raux     = xx(2,col)
                    yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux
                    yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux
                    yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux
                    yy(4,ii) = yy(4,ii) + an(2,4,jj) * raux
                    raux     = xx(3,col)
                    yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux
                    yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux
                    yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux
                    yy(4,ii) = yy(4,ii) + an(3,4,jj) * raux
                    raux     = xx(4,col)
                    yy(1,ii) = yy(1,ii) + an(4,1,jj) * raux
                    yy(2,ii) = yy(2,ii) + an(4,2,jj) * raux
                    yy(3,ii) = yy(3,ii) + an(4,3,jj) * raux
                    yy(4,ii) = yy(4,ii) + an(4,4,jj) * raux
                 end if
              end do

           end do

        else
           !
           ! NBVAR = whatever
           !
           !$OMP PARALLEL  DO                                            &
           !$OMP SCHEDULE ( STATIC )                                     &
           !$OMP DEFAULT  ( NONE )                                       &
           !$OMP SHARED   ( an, ia, ja, nbnodes, nbvar, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
           !
           do ii = 1,nbnodes
              do kk = 1,nbvar
                 yy(kk,ii) = 0.0_rp
              end do
              do jj  = ia(ii),ia(ii+1)-1
                 col = ja(jj)
                 if( col <= npoi1 ) then
                    do ll = 1,nbvar
                       raux = xx(ll,col)
                       do kk = 1,nbvar
                          yy(kk,ii) = yy(kk,ii) + an(ll,kk,jj) * raux
                       end do
                    end do
                 end if
              end do
           end do

        end if
        !
        ! Wait and assemble
        !
        call cputim(time4)
        call PAR_INTERFACE_OWN_NODE_EXCHANGE(nbvar,xx,'WAIT AND ASSEMBLE')
        call cputim(time5)
        !
        ! Interior nodes (missing rectangle)
        !
        !$OMP PARALLEL  DO                                            &
        !$OMP SCHEDULE ( STATIC )                                     &
        !$OMP DEFAULT  ( NONE )                                       &
        !$OMP SHARED   ( an, ia, ja, nbnodes, nbvar, npoi1, xx, yy )  &
        !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
        !
        do ii = 1,nbnodes
           do jj  = ia(ii),ia(ii+1)-1
              col = ja(jj)
              if( col > npoi1 ) then
                 do ll = 1,nbvar
                    raux = xx(ll,col)
                    do kk = 1,nbvar
                       yy(kk,ii) = yy(kk,ii) + an(ll,kk,jj) * raux
                    end do
                 end do
              end if
           end do
        end do

     else

        !-------------------------------------------------------------------
        !
        ! Partial row
        !
        !-------------------------------------------------------------------
        !
        ! Boundary nodes
        !
        ! Compute boundary elements

        if( nbvar == 1 ) then
           !
           ! NBVAR=1
           !
           !$OMP PARALLEL  DO                                     &
           !$OMP SCHEDULE ( STATIC )                              &
           !$OMP DEFAULT  ( NONE )                                &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, raux )
           !
           do ii = npoi1+1,nbnodes
              yy(1,ii) = 0.0_rp
              do jj   = ia(ii),ia(ii+1)-1
                 col  = ja(jj)
                 raux = xx(1,col)
                 yy(1,ii) = yy(1,ii) +an(1,1,jj) * raux
              end do
           end do

        else if( nbvar == 2 ) then
           !
           ! NBVAR=2
           !
           !$OMP PARALLEL  DO                                     &
           !$OMP SCHEDULE ( STATIC )                              &
           !$OMP DEFAULT  ( NONE )                                &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
           !
           do ii = npoi1+1,nbnodes
              yy(1,ii) = 0.0_rp
              yy(2,ii) = 0.0_rp
              do jj       = ia(ii),ia(ii+1)-1
                 col      = ja(jj)
                 raux1    = xx(1,col)
                 raux2    = xx(2,col)
                 yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
                 yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
                 yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
                 yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
              end do
           end do

        else if( nbvar == 3 ) then
           !
           ! NBVAR=3
           !
           !$OMP PARALLEL  DO                                       &
           !$OMP SCHEDULE ( STATIC )                                &
           !$OMP DEFAULT  ( NONE )                                  &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )    &
           !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
           !
           do ii = npoi1+1,nbnodes
              yy(1,ii) = 0.0_rp
              yy(2,ii) = 0.0_rp
              yy(3,ii) = 0.0_rp
              do jj       = ia(ii),ia(ii+1)-1
                 col      = ja(jj)
                 raux1    = xx(1,col)
                 raux2    = xx(2,col)
                 raux3    = xx(3,col)
                 yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
                 yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
                 yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux3
                 yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
                 yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
                 yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux3
                 yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux1
                 yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux2
                 yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux3
              end do

           end do

        else if( nbvar == 4 ) then
           !
           ! NBVAR=4
           !
           !$OMP PARALLEL  DO                                     &
           !$OMP SCHEDULE ( STATIC )                              &
           !$OMP DEFAULT  ( NONE )                                &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, raux )
           !
           do ii = npoi1+1,nbnodes
              yy(1,ii) = 0.0_rp
              yy(2,ii) = 0.0_rp
              yy(3,ii) = 0.0_rp
              yy(4,ii) = 0.0_rp
              do jj       = ia(ii),ia(ii+1)-1
                 col      = ja(jj)
                 raux     = xx(1,col)
                 yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
                 yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux
                 yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux
                 yy(4,ii) = yy(4,ii) + an(1,4,jj) * raux
                 raux     = xx(2,col)
                 yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux
                 yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux
                 yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux
                 yy(4,ii) = yy(4,ii) + an(2,4,jj) * raux
                 raux     = xx(3,col)
                 yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux
                 yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux
                 yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux
                 yy(4,ii) = yy(4,ii) + an(3,4,jj) * raux
                 raux     = xx(4,col)
                 yy(1,ii) = yy(1,ii) + an(4,1,jj) * raux
                 yy(2,ii) = yy(2,ii) + an(4,2,jj) * raux
                 yy(3,ii) = yy(3,ii) + an(4,3,jj) * raux
                 yy(4,ii) = yy(4,ii) + an(4,4,jj) * raux
              end do

           end do

        else
           !
           ! NBVAR = whatever
           !
           !$OMP PARALLEL  DO                                            &
           !$OMP SCHEDULE ( STATIC )                                     &
           !$OMP DEFAULT  ( NONE )                                       &
           !$OMP SHARED   ( an, ia, ja, nbnodes, nbvar, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
           !
           do ii = npoi1+1,nbnodes
              do kk = 1,nbvar
                 yy(kk,ii) = 0.0_rp
              end do
              do jj  = ia(ii),ia(ii+1)-1
                 col = ja(jj)
                 do ll = 1,nbvar
                    raux = xx(ll,col)
                    do kk = 1,nbvar
                       yy(kk,ii) = yy(kk,ii) + an(ll,kk,jj) * raux
                    end do
                 end do
              end do
           end do

        end if

        ! Send/recv boundary elements contributions asynchronously

        !-------------------------------------------------------------------
        !
        ! Modify YY due do periodicity and Parall service
        !
        !-------------------------------------------------------------------

        call cputim(time2)
        if( itask == 1 .and. IPARALL ) then
           call pararr('SLA',NPOIN_TYPE,nbnodes*nbvar,yy)
        end if
        call cputim(time3)

        ! Compute interior elements async

        !
        ! Interior nodes
        !
        if( nbvar == 1 ) then
           !
           ! NBVAR=1
           !
           !$OMP PARALLEL  DO                                     &
           !$OMP SCHEDULE ( STATIC )                              &
           !$OMP DEFAULT  ( NONE )                                &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, raux )
           !
           do ii = 1,npoi1
              yy(1,ii) = 0.0_rp
              do jj   = ia(ii),ia(ii+1)-1
                 col  = ja(jj)
                 raux = xx(1,col)
                 yy(1,ii) = yy(1,ii) +an(1,1,jj) * raux
              end do
           end do

        else if( nbvar == 2 ) then
           !
           ! NBVAR=2
           !
           !$OMP PARALLEL  DO                                     &
           !$OMP SCHEDULE ( STATIC )                              &
           !$OMP DEFAULT  ( NONE )                                &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )
           !
           do ii = 1,npoi1
              yy(1,ii) = 0.0_rp
              yy(2,ii) = 0.0_rp
              do jj       = ia(ii),ia(ii+1)-1
                 col      = ja(jj)
                 raux1    = xx(1,col)
                 raux2    = xx(2,col)
                 yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
                 yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
                 yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
                 yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
              end do
           end do

        else if( nbvar == 3 ) then
           !
           ! NBVAR=3
           !
           !$OMP PARALLEL  DO                                       &
           !$OMP SCHEDULE ( STATIC )                                &
           !$OMP DEFAULT  ( NONE )                                  &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )    &
           !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )
           !
           do ii = 1,npoi1
              yy(1,ii) = 0.0_rp
              yy(2,ii) = 0.0_rp
              yy(3,ii) = 0.0_rp
              do jj       = ia(ii),ia(ii+1)-1
                 col      = ja(jj)
                 raux1    = xx(1,col)
                 raux2    = xx(2,col)
                 raux3    = xx(3,col)
                 yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
                 yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
                 yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux3
                 yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
                 yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
                 yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux3
                 yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux1
                 yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux2
                 yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux3
              end do

           end do

        else if( nbvar == 4 ) then
           !
           ! NBVAR=4
           !
           !$OMP PARALLEL  DO                                     &
           !$OMP SCHEDULE ( STATIC )                              &
           !$OMP DEFAULT  ( NONE )                                &
           !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, raux )
           !
           do ii = 1,npoi1
              yy(1,ii) = 0.0_rp
              yy(2,ii) = 0.0_rp
              yy(3,ii) = 0.0_rp
              yy(4,ii) = 0.0_rp
              do jj       = ia(ii),ia(ii+1)-1
                 col      = ja(jj)
                 raux     = xx(1,col)
                 yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
                 yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux
                 yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux
                 yy(4,ii) = yy(4,ii) + an(1,4,jj) * raux
                 raux     = xx(2,col)
                 yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux
                 yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux
                 yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux
                 yy(4,ii) = yy(4,ii) + an(2,4,jj) * raux
                 raux     = xx(3,col)
                 yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux
                 yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux
                 yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux
                 yy(4,ii) = yy(4,ii) + an(3,4,jj) * raux
                 raux     = xx(4,col)
                 yy(1,ii) = yy(1,ii) + an(4,1,jj) * raux
                 yy(2,ii) = yy(2,ii) + an(4,2,jj) * raux
                 yy(3,ii) = yy(3,ii) + an(4,3,jj) * raux
                 yy(4,ii) = yy(4,ii) + an(4,4,jj) * raux
              end do

           end do

        else
           !
           ! NBVAR = whatever
           !
           !$OMP PARALLEL  DO                                            &
           !$OMP SCHEDULE ( STATIC )                                     &
           !$OMP DEFAULT  ( NONE )                                       &
           !$OMP SHARED   ( an, ia, ja, nbnodes, nbvar, npoi1, xx, yy )  &
           !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )
           !
           do ii = 1,npoi1
              do kk = 1,nbvar
                 yy(kk,ii) = 0.0_rp
              end do
              do jj  = ia(ii),ia(ii+1)-1
                 col = ja(jj)
                 do ll = 1,nbvar
                    raux = xx(ll,col)
                    do kk = 1,nbvar
                       yy(kk,ii) = yy(kk,ii) + an(ll,kk,jj) * raux
                    end do
                 end do
              end do
           end do

        end if


        ! Wait boundary elements contributions async

        !-------------------------------------------------------------------
        !
        ! Wait all and sum up contributions of boundary nodes
        !
        !-------------------------------------------------------------------

        call cputim(time4)
        if( itask == 1 .and. IPARALL ) then
           call pararr('SLA',NPOIN_TYPE,nbnodes*nbvar,yy)
        end if
        call cputim(time5)

     end if

     time_spmv                  = time2 + time4 - ( time1 + time3 )
     solve_sol(1) % num_spmv    = solve_sol(1) % num_spmv    + 1
     solve_sol(1) % cpu_spmv(1) = solve_sol(1) % cpu_spmv(1) + time_spmv
     solve_sol(1) % cpu_spmv(2) = solve_sol(1) % cpu_spmv(2) + time5 - ( time1 + time_spmv )


     ! Finish


  end if

end subroutine bcsrai

