subroutine par_receiv()
!------------------------------------------------------------------------
!****f* Parall/par_receiv
! NAME
!    par_receiv
! DESCRIPTION
!    This routine Receive all buffers from process 'kfl_desti'
! OUTPUT
!   
! USED BY
!    Parall
!***
!------------------------------------------------------------------------
  use def_parall
  use def_master
  use mod_iofile
  use mod_parall,               only : PAR_INTEGER
  use mod_parall,               only : PAR_COMM_MY_CODE4
  use mod_communications_tools, only : PAR_MPI_RUNEND
  use mod_memory
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
  integer(4)     :: status(MPI_STATUS_SIZE)
#endif
  integer(4)     :: istat,kfl_desti_par4,npari4,nparr4,nparc4,iunit4
  integer(ip)    :: ipari,iparr,iunit
  real(rp)       :: time1,time2
  character(300) :: messa
  character(20)  :: cnume

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
        npari4 = int(npari,4)
        call MPI_Recv( parin(1:npari), npari4, PAR_INTEGER, kfl_desti_par4, 0_4, PAR_COMM_MY_CODE4, status, istat )
        if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
        npari = 0
     end if
     
     if( nparr > 0 ) then
        nparr4 = int(nparr,4)
        call MPI_Recv( parre(1:nparr), nparr4, MPI_DOUBLE_PRECISION, kfl_desti_par4, 0_4, PAR_COMM_MY_CODE4, status, istat )
        if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
        nparr = 0
     end if

     if( nparc > 0 ) then
        nparc4 = int(nparc,4)
        call MPI_Recv( parch(1:nparc), nparc4, MPI_CHARACTER, kfl_desti_par4, 0_4, PAR_COMM_MY_CODE4, status, istat )
        if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
        nparc = 0
     end if
#endif

  else 

     !-------------------------------------------------------------------
     !
     ! Read from file
     !
     !-------------------------------------------------------------------

     iunit  = lun_aonlp_par !+ kfl_paral
     iunit4 = int(iunit,4)

     if(kfl_ascii_par==0) then
        read(iunit4,err=1) npari,nparr,nparc
        if( npari > 0 ) read(iunit4,err=1,end=1)  ( parin(ipari),   ipari=1,npari )
        if( nparr > 0 ) read(iunit4,err=1,end=1)  ( parre(iparr),   iparr=1,nparr )
        if( nparc > 0 ) read(iunit4,err=1,end=1)    parch(1:nparc)
     else
        read(iunit4,*,err=1) npari,nparr,nparc
        read(iunit4,*,err=1) strin,strre,strch
        if( npari > 0 ) read(iunit4,*,err=1,end=1) ( parin(ipari),  ipari=1,npari )
        if( nparr > 0 ) read(iunit4,*,err=1,end=1) ( parre(iparr),  iparr=1,nparr )
        if( nparc > 0 ) read(iunit4,*,err=1,end=1)   parch(1:nparc)      
     end if
     npari=0
     nparr=0
     nparc=0

  end if

  call cputim(time2)
  cpu_paral(22)=cpu_paral(22)+time2-time1
  return

1 cnume=intost(kfl_paral)
  messa='PARALL: ERROR WHILE SLAVE '//trim(cnume)&
       //' IS READING RESTART FILE. CHECK FILE FORMAT.'&
       //'TRYING TO READ: '//trim(strin)&
       //', '//trim(strre)//', '//trim(strch)
  call runend(messa)

end subroutine par_receiv
