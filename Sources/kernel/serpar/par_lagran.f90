subroutine par_lagran(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_lagran
  ! NAME
  !    par_lagran
  ! DESCRIPTION
  !    This subroutine exchange arrays between slaves
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall, only : PAR_INTEGER,commd
  use mod_parall, only : PAR_COMM_WORLD,PAR_COMM_MY_CODE4
  use mod_parall, only : par_memor
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
  integer(4)               :: status(MPI_STATUS_SIZE)
#endif
  integer(ip), intent(in)  :: itask
  integer(ip)              :: ii,jj
  integer(4)               :: zero4,npari4,npasi4,npasr4,bsize4,dom_i
  integer(4)               :: istat
  integer(ip), pointer     :: loc_rpari1(:)

  if( IPARALL ) then

     select case ( itask )

     case ( 1_ip ) 
        !
        ! GATHER
        !
        zero4  = 0_4
        npari4 = int( npari , 4_4 )
        npasi4 = int( npasi , 4_4 )
#ifdef MPI_OFF
#else
        CALL MPI_Gather( paris(1:) , npasi4 , PAR_INTEGER, parig(1:) , npari4 , &
             PAR_INTEGER , zero4 , PAR_COMM_MY_CODE4 , status)
#endif
        npari = 0
        npasi = 0

     case ( 2_ip ) 
        !
        ! GATHERV
        !
        ! parrs: sendbuf
        ! npasr: sendcount
        ! parre: recvbuf
        ! parig: recvcount array 
        ! pari1: displacements
        ! zero:  root (master receives)
        !
        zero4  = 0_4
        npasr4 = int( npasr , 4_4 )
#ifdef MPI_OFF
#else
        call MPI_Gatherv( parrs, npasr4, MPI_DOUBLE_PRECISION, parre, parig, pari1, &
             MPI_DOUBLE_PRECISION, zero4, PAR_COMM_MY_CODE4, status )
#endif
        npasr = 0

     case ( 3_ip )
        !
        ! Check repeated particles
        !
        if( ISLAVE ) then
           allocate(loc_rpari1(npari),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)
           pard1 = 0

           do ii = 1, nneig

              dom_i = int(commd % neights(ii),4)            

#ifdef MPI_OFF
#else
              bsize4 = int(npari,4)
              call MPI_Sendrecv(                 &
                   pari1(1:),  bsize4,           &
                   PAR_INTEGER,  dom_i, 0_4,     &
                   loc_rpari1(1:), bsize4,       &
                   PAR_INTEGER, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif

              do jj = 1,npari
                 if( pari1(jj) > 0 .and. loc_rpari1(jj) > 0 ) then
                    if( kfl_paral > dom_i ) then
                       pard1 = 1
                       pari1(jj) = -abs(pari1(jj))
                    end if
                 end if
              end do

           end do

           call memchk(two,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)
           deallocate(loc_rpari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARI1','par_slexch',0_ip)

        end if

     case ( 4_ip ) 

        zero4  = 0_4
        npari4 = int( npari , 4_4 )
        npasi4 = int( npasi , 4_4 )

#ifdef MPI_OFF
#else
        !
        ! paris: sendbuf
        ! npasi: sendcount
        ! parin: recvbuf
        ! parig: recvcount array 
        ! pari1: displacements
        ! zero:  root (master receives)
        !
        call MPI_Gatherv( paris, npasi4, PAR_INTEGER, parin, parig, pari1, &
             PAR_INTEGER, zero4, PAR_COMM_MY_CODE4, status )
#endif

        npari = 0
        npasi = 0

      end select

  end if

end subroutine par_lagran
