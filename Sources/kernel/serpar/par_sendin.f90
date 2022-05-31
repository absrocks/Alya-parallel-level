subroutine par_sendin()
  !------------------------------------------------------------------------
  !****f* Parall/par_send
  ! NAME
  !    par_send
  ! DESCRIPTION
  !    This routine Send all buffers to process 'kfl_desti'
  ! OUTPUT
  !   
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_parall
  use def_master
  use mod_memchk
  use mod_iofile
  use mod_par_virfil
  use mod_parall,               only : PAR_INTEGER
  use mod_parall,               only : PAR_COMM_MY_CODE4
  use mod_communications_tools, only : PAR_MPI_RUNEND
  use mod_memory
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
#endif
  integer(4)             :: istat,kfl_desti_par4,npari4,nparr4,nparc4,iunit4
  integer(4)             :: nparr42
  integer(ip)            :: iunit,ivari,nparr2,iparx
  real(rp)               :: time1,time2  
  character(150)         :: cfile
  real(rp),      pointer :: parrx(:)

  call cputim(time1)
  kfl_desti_par4=int(kfl_desti_par,4)

  if( PART_AND_RUN() ) then

     !-------------------------------------------------------------------
     !
     ! MPI communication
     !
     !-------------------------------------------------------------------

#ifdef MPI_OFF
#else
     if( npari > 0 ) then
        npari4=int(npari,4)
        call MPI_Send( parin(1:npari), npari4, PAR_INTEGER, kfl_desti_par4, 0_4, PAR_COMM_MY_CODE4, istat )
        if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
        npari=0
     end if

     if( nparr > 0 ) then
        nparr4=int(nparr,4)
        call MPI_Send( parre(1:nparr), nparr4, MPI_DOUBLE_PRECISION, kfl_desti_par4, 0_4, PAR_COMM_MY_CODE4, istat )
        if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
        nparr=0
     end if

     if( nparc > 0 ) then
        nparc4=int(nparc,4)
        call MPI_Send( parch(1:nparc), nparc4, MPI_CHARACTER, kfl_desti_par4, 0_4, PAR_COMM_MY_CODE4, istat )
        if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
        nparc=0
     end if

     if( nparx > 0 ) then
        nparr2  = 2*nparx
        nparr42 = int(nparr2,4)
        allocate( parrx(nparr2) )
        ivari = 0
        do iparx = 1,nparx
           ivari = ivari + 1
           parrx(ivari) = real  ( parcx(iparx) )
           ivari = ivari + 1
           parrx(ivari) = aimag ( parcx(iparx) )
        end do
        call MPI_Send( parrx(1:nparr2), nparr42, MPI_DOUBLE_PRECISION, kfl_desti_par4, 0_4, PAR_COMM_MY_CODE4, istat )
        if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
        nparx=0
        deallocate( parrx )
     end if

#endif

  else if( ( IMASTER .or. nproc_par > 1 ) .and. PART_AND_WRITE() ) then

     !-------------------------------------------------------------------
     !
     ! File communication
     !
     !-------------------------------------------------------------------
     
     iunit  = lun_aonlp_par + kfl_desti_par
     iunit4 = int(iunit,4) 

     if( kfl_virfi_par == 1 ) then
        !
        ! Virtual file
        !
        call par_wribuf(kfl_desti_par)

     else 
        !
        ! Binary/ASCII format
        !
        if( kfl_filio_par == 1 ) then
           call par_filnam(1_ip,kfl_desti_par,fil_rstar_par,cfile)
           if( kfl_ascii_par == 0 ) then
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','unformatted','append')
           else
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','formatted',  'append')      
           end if
        end if

        call par_wrifil(iunit)

        if( kfl_filio_par == 1 )  close(iunit4)

     end if

     npari = 0
     nparr = 0
     nparc = 0
     nparx = 0

  end if

  call cputim(time2)
  cpu_paral(21)=cpu_paral(21)+time2-time1

end subroutine par_sendin
