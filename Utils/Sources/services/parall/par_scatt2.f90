subroutine par_scatt2(iloca)
  !-----------------------------------------------------------------------
  !****f* parall/par_scatt2
  ! NAME 
  !    par_scatt2
  ! DESCRIPTION
  ! USES
  !
  ! USED BY
  !
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_parall
  use mod_iofile
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : PAR_COMM_MY_CODE4
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
#endif
  integer(ip), intent(out) :: iloca
  integer(ip)              :: iunit,dummi
  integer(4)               :: istat,iloca4
  character(300)           :: messa
  character(150)           :: cfile
  character(20)            :: cnume

  iloca4=int(iloca,4)

  if(kfl_ptask==1) then
     !
     ! Normal run
     !
#ifdef MPI_OFF 
#else
     if(IMASTER) then
        call MPI_SCATTER(&
             parin, 1_4, PAR_INTEGER, iloca, 1_4, PAR_INTEGER, 0_4,&
             PAR_COMM_MY_CODE4, istat)
     else if(ISLAVE) then
        call MPI_SCATTER(&
             dummi, 1_4, PAR_INTEGER, iloca, 1_4, PAR_INTEGER, 0_4,&
             PAR_COMM_MY_CODE4, istat)        
     end if

#endif

  else if(kfl_ptask==0.and.IMASTER) then
     !
     ! Only preprocess run
     !
     do kfl_desti_par=1,npart_par
        iunit = lun_aonlp_par + kfl_desti_par
        if(kfl_ascii_par==0) then
           !
           ! ASCII Format
           !
           if(kfl_filio_par==1) then
              call par_filnam(1_ip,kfl_desti_par,fil_rstar_par,cfile)
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','unformatted','append')
           end if
           write(iunit) npari
           if(npari>0) write(iunit) parin(kfl_desti_par+1)
           if(kfl_filio_par==1) close(iunit)
        else
           !
           ! Binary format
           !
           if(kfl_filio_par==1) then
              call par_filnam(1_ip,kfl_desti_par,fil_rstar_par,cfile)
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','formatted','append')
           end if
           write(iunit,*) npari
           write(iunit,*) trim(strin)
           if(npari>0) write(iunit,*) parin(kfl_desti_par+1)
           if(kfl_filio_par==1) close(iunit)
        end if
     end do
     npari=0
     strin='NULL'

  else if(kfl_ptask==2.and.kfl_paral/=0) then
     !
     ! Only run
     !
     iunit= lun_aonlp_par + kfl_paral
     if(kfl_ascii_par==0) then
        read(iunit,err=1) npari
        if(npari>0) read(iunit,err=1,end=1) iloca
     else
        read(iunit,*,err=1) npari
        read(iunit,*,err=1) strin
        if(npari>0) read(iunit,*,err=1,end=1) iloca
     end if
     npari=0   
     strin='NULL'

  end if

  return

1 cnume=intost(kfl_paral)
  messa='PARALL: ERROR WHILE SLAVE '//trim(cnume)&
       //' IS READING RESTART FILE. CHECK FILE FORMAT.'&
       //'TRYING TO READ: '//trim(strin)&
       //', '//trim(strre)//', '//trim(strch)
  call runend(messa)

end subroutine par_scatt2
