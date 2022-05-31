subroutine par_sledif()
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall, only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall, only : commd,PAR_COMM_MY_CODE4
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : par_memor
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
  integer(4)               :: status(MPI_STATUS_SIZE)
#endif
  integer(ip)              :: ipoin,ii,jj,bsize,ji,poin,ini,dom_i,ibopo,ineig
  integer(ip)              :: kpoin,kk
  integer(4)               :: istat,bsize4
  real(rp)                 :: time1,time2
  integer(ip), allocatable :: loc_spari1(:),   loc_rpari1(:)
  real(rp),    allocatable :: loc_sparr1(:),   loc_rparr1(:)
  complex(rp), allocatable :: loc_sparx1(:),   loc_rparx1(:)

  call cputim(time1)

  if(kfl_paral>0) then

     !-------------------------------------------------------------
     !
     ! REAL(NPOIN)
     !
     !-------------------------------------------------------------

     allocate(loc_sparr1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)

     allocate(loc_rparr1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)

     do jj= 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        loc_sparr1(jj) = parr1(ipoin)
     enddo

     do ii= 1, nneig
        dom_i = commd%neights(ii)

        ini   = commd%bound_size(ii)
        bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
        bsize4=int(bsize,4)
        call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
             MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
             loc_rparr1(ini:), bsize4,              &
             MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
             PAR_COMM_MY_CODE4, status, istat )

#endif
     enddo

     do jj = 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        parr1(ipoin) = 0.0_rp
     end do

     do jj = 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        parr1(ipoin) = max( parr1(ipoin), abs(loc_sparr1(jj)-loc_rparr1(jj)))
     enddo

     call memchk(two,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)
     deallocate(loc_rparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

     call memchk(two,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)
     deallocate(loc_sparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

  endif

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_sledif
