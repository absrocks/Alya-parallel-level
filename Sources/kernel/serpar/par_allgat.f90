subroutine par_allgat(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_allgat
  ! NAME
  !    par_allgat
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
  use mod_parall, only : PAR_INTEGER,PAR_COMM_MY_CODE4
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
  integer(4)               :: status(MPI_STATUS_SIZE)
#endif
  integer(ip), intent(in)  :: itask
  integer(ip)              :: isize,ii
  integer(4)               :: zero4,npari4,npasi4,npasr4
  integer(4),  pointer     :: parig_4(:),pari1_4(:)

  if( IPARALL ) then

     select case ( itask )

     case ( 1_ip ) 
        !
        ! All gather an array PARIS(NPARI) in variable PARIN
        !
        npari4 = int( npari , 4_4 )
        
#ifdef MPI_OFF
#else
        CALL MPI_AllGather( &
             paris(1:) , npari4 , PAR_INTEGER , &
             parin(1:) , npari4 , PAR_INTEGER , &
             PAR_COMM_MY_CODE4 , status)
#endif

     case ( 2_ip ) 
        !
        ! All gather an array PARIS(NPARI) in variable PARIN
        !
        zero4  = 0_4
        npasr4 = int( npasr , 4_4 )
#ifdef MPI_OFF
#else
        !call MPI_AllGatherv( parrs, npasr4, MPI_DOUBLE_PRECISION, parre, parig, pari1, &
        !     MPI_DOUBLE_PRECISION, PAR_COMM_MYCODE, status )
        if( ip /= 4 ) then
           isize = size(parig)
           allocate( parig_4(isize) )
           do ii = 1,isize
              parig_4(ii) = int(parig(ii),4)
           end do
           isize = size(pari1)
           allocate( pari1_4(isize) )
           do ii = 1,isize
              pari1_4(ii) = int(pari1(ii),4)
           end do
           call MPI_AllGatherv( &
                parrs, npasr4, MPI_DOUBLE_PRECISION, &
                parre, parig_4, pari1_4, &
                MPI_DOUBLE_PRECISION, PAR_COMM_MY_CODE4, status )
           deallocate( parig_4 , pari1_4 )
        else
           call MPI_AllGatherv( &
                parrs, npasr4, MPI_DOUBLE_PRECISION, &
                parre, parig, pari1, &
                MPI_DOUBLE_PRECISION, PAR_COMM_MY_CODE4, status )
        end if
#endif

     case ( 3_ip ) 
        !
        ! All gather an array PARIS(NPARI) in variable PARIN
        !
        zero4  = 0_4
        npasi4 = int( npasi , 4_4 )
#ifdef MPI_OFF
#else
        if( ip /= 4 ) then
           isize = size(parig)
           allocate( parig_4(isize) )
           do ii = 1,isize
              parig_4(ii) = int(parig(ii),4)
           end do
           isize = size(pari1)
           allocate( pari1_4(isize) )
           do ii = 1,isize
              pari1_4(ii) = int(pari1(ii),4)
           end do
           call MPI_AllGatherv( &
                paris, npasi4, PAR_INTEGER, &
                parin, parig_4,  pari1_4,       &
                PAR_INTEGER, PAR_COMM_MY_CODE4, status )
           deallocate( parig_4 , pari1_4 )
        else
           call MPI_AllGatherv( &
                paris, npasi4, PAR_INTEGER, &
                parin, parig,  pari1,       &
                PAR_INTEGER, PAR_COMM_MY_CODE4, status )
        end if
#endif

     end select

  end if

end subroutine par_allgat
