
subroutine bcsrax_2(itask,nbnodes,nbvar,an,ja,ia,xx,yy) 
  !----------------------------------------------------------------------
  !****f* mathru/bcsrax_2
  ! NAME 
  !     bcsrax_2
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
  use def_kintyp, only             :  ip,rp
  use def_master, only             :  INOTMASTER,kfl_async,IPARALL
  use def_master, only             :  NPOIN_TYPE
  implicit none
  integer(ip), intent(in)          :: itask,nbnodes,nbvar
  real(rp),    intent(in)          :: an(nbvar,nbvar,*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  real(rp),    intent(in)          :: xx(nbvar,*)
  real(rp),    intent(out), target :: yy(nbvar,*)
  integer(ip)                      :: ii,jj,kk,ll,col
  real(rp)                         :: raux,raux1,raux2,raux3

  !if( IPARALL .and. kfl_async == 1 ) then
  !
  !
  !   call bcsrai(itask,nbnodes,nbvar,an,ja,ia,xx,yy) 

  if( INOTMASTER ) then

     if( nbvar == 1 ) then
        !
        ! NBVAR=1
        !
        !*OMP   PARALLEL DO SCHEDULE (GUIDED)        & 
        !*OMP   DEFAULT (NONE)                       &
        !*OMP   PRIVATE ( ii, jj, col, raux)         &
        !*OMP   SHARED ( nbnodes, xx, yy, ia, ja, an)
        do ii = 1,nbnodes
           yy(1,ii) = 0.0_rp
           do jj   = ia(ii),ia(ii+1)-1
              col  = ja(jj)
              raux = xx(1,col)
              yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
           end do
        end do

     else if( nbvar == 2 ) then
        !
        ! NBVAR=2
        !
        !*OMP   PARALLEL DO SCHEDULE (GUIDED)         & 
        !*OMP   DEFAULT (NONE)                        &
        !*OMP   PRIVATE ( ii, jj, col, raux1, raux2)  &
        !*OMP   SHARED ( nbnodes, xx, yy, ia, ja, an)
        do ii = 1,nbnodes
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
        !*OMP   PARALLEL DO SCHEDULE (GUIDED)             & 
        !*OMP   DEFAULT (NONE)                            &
        !*OMP   PRIVATE ( ii, jj, col, raux1,raux2,raux3) &
        !*OMP   SHARED ( nbnodes, xx, yy, ia, ja, an)
        do ii = 1,nbnodes
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
        !*OMP   PARALLEL DO SCHEDULE (GUIDED)        & 
        !*OMP   DEFAULT (NONE)                       &
        !*OMP   PRIVATE ( ii, jj, col, raux)         &
        !*OMP   SHARED ( nbnodes, xx, yy, ia, ja, an)
        do ii = 1,nbnodes
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
        !*OMP   PARALLEL DO SCHEDULE (GUIDED)               & 
        !*OMP   DEFAULT (NONE)                              &
        !*OMP   PRIVATE ( ii, jj, kk, ll, col, raux)        &
        !*OMP   SHARED ( nbnodes, nbvar, xx, yy, ia, ja, an)
        do ii = 1,nbnodes
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
     !
     ! Modify YY due do periodicity and Parall service
     !
     call pararr('SSS',NPOIN_TYPE,nbvar*nbnodes,yy)

  end if

end subroutine bcsrax_2

