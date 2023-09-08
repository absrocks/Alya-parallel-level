subroutine par_sen2ma()
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
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
#endif
  integer(4)             :: kfl_desti_par4
  integer(ip)            :: pleng
  real(rp)               :: time1,time2  

  call cputim(time1)

#ifdef MPI_OFF
#else
  
  kfl_desti_par = 0
  kfl_desti_par4 = int(kfl_desti_par,4)

  if( npari > 0 ) then
     pleng = npari
     if( pardi == 1 ) then
        call par_parari('SND',0_ip,pleng,pari1)
     else
        call par_parari('SND',0_ip,pleng,pari2)
     end if
  end if

  if( nparr > 0 ) then
     pleng = nparr
     if( pardi == 1 ) then
        call par_pararr('SND',0_ip,pleng,parr1)
     else
        call par_pararr('SND',0_ip,pleng,parr2)
     end if
  end if

  if( nparx > 0 ) then
     call runend('PAR_SEN2MA: NOT DONE')
  end if

#endif

  nparc = 0
  npari = 0
  nparr = 0
  nparx = 0

  call cputim(time2)
  cpu_paral(21)=cpu_paral(21)+time2-time1

end subroutine par_sen2ma
