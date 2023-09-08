subroutine par_slonly()
  !-----------------------------------------------------------------------
  !****f* Parall/par_slonly
  ! NAME
  !    par_slonly
  ! DESCRIPTION
  !    This subroutine exchange arrays between slaves.
  !    If the array is different from zero on different subdomains
  !    pout it to zero except on the one with lowest rank
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall, only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall, only : commd,PAR_COMM_MY_CODE4
  use mod_parall, only : PAR_INTEGER
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

  if( IPARALL ) then

     if( party == NPOIN_TYPE ) then
        !
        ! Node
        !
        if( pardi == 1 .and. parki == 1 ) then

           !-------------------------------------------------------------
           !
           ! INT(NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_spari1(commd % bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_slonly',loc_spari1)

           allocate(loc_rpari1(commd % bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_slonly',loc_rpari1)

           do jj = 1,commd % bound_dim
              ipoin = commd % bound_perm(jj)
              loc_spari1(jj) = pari1(ipoin)
           end do

           do ii = 1,nneig

              dom_i = commd % neights(ii)
              ini   = commd % bound_size(ii)
              bsize = commd % bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv(                 &
                   loc_spari1(ini:), bsize4,     &
                   PAR_INTEGER,  dom_i, 0_4,     &
                   loc_rpari1(ini:), bsize4,     &
                   PAR_INTEGER, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           end do

           do ii = 1,nneig
              dom_i = commd % neights(ii)
              do jj = commd % bound_size(ii),commd % bound_size(ii+1)-1
                 ipoin = commd % bound_perm(jj)
                 if( pari1(ipoin) /= 0 .and. loc_rpari1(jj) /= 0 ) then
                    if( kfl_paral > dom_i ) pari1(ipoin) = 0

                    

                 end if
              end do
           end do

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_slonly',loc_rpari1)
           deallocate(loc_rpari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARI1','par_slonly',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_slonly',loc_spari1)
           deallocate(loc_spari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARI1','par_slonly',0_ip)

        else 

           call runend('PAR_SLONLY: NOT CODED')

        end if

     else

        call runend('PAR_SLONLY: NOT CODED')

     end if
  end if

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slonly
