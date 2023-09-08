subroutine par_sumset()
  !------------------------------------------------------------------------
  !****f* Parall/par_sumset
  ! NAME
  !    par_sumset
  ! DESCRIPTION
  !    This routine operates for sets
  !    ITASK=1 ... Minimum 
  !    ITASK=2 ... Maximum
  !    ITASK=3 ... Sum
  ! OUTPUT
  !    NPARI ..... Integer array
  !    NPARR ..... Real array
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_parall
  use def_master
  use mod_parall, only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall, only : commd,PAR_COMM_MY_CODE4
  use mod_parall, only : PAR_INTEGER
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
#endif
  integer(ip)           :: ii,ndim1,ndim2,idim1,idim2
  integer(4)            :: istat=0,nparr4
  real(rp)              :: time1,time2
  real(rp),    pointer  :: rwa(:)
  !
  ! Operations on real arrays
  !
  ndim1 = size(parr2,1)
  ndim2 = size(parr2,2)
  nparr = ndim1 * ndim2

  if( nparr > 0 ) then 

#ifdef MPI_OFF
#else

     call cputim(time1)
     allocate(rwa(nparr),   stat=istat)
     allocate(parre(nparr), stat=istat)

     if( ISLAVE ) then
        ii = 0
        do idim2 = 1,ndim2
           do idim1 = 1,ndim1
              ii = ii + 1
              rwa(ii) = parr2(idim1,idim2)
           end do
        end do
     else
        do ii = 1,nparr
           rwa(ii) = 0.0_rp
        end do        
     end if

     nparr4 = int(nparr,4)
     call MPI_AllReduce(rwa,parre,nparr4,MPI_DOUBLE_PRECISION,MPI_SUM,PAR_COMM_MY_CODE4,istat)        

     ii = 0
     do idim2 = 1,ndim2
        do idim1 = 1,ndim1
           ii = ii + 1
           parr2(idim1,idim2) = parre(ii)
        end do
     end do

     deallocate(rwa,  stat=istat)
     deallocate(parre,stat=istat)
     call cputim(time2)
     cpu_paral(22) = cpu_paral(22) + time2 - time1
     nparr = 0
#endif

  end if

  if( istat /= 0 ) call runend('PARALL: FUNCTION MPI_ALLREDUCE HAS FAILED')

end subroutine par_sumset
