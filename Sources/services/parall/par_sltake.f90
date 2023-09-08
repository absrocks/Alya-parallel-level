subroutine par_sltake
  !-----------------------------------------------------------------------
  !****f* Parall/par_sltake
  ! NAME
  !    par_sltake
  ! DESCRIPTION
  !    This subroutine exchange arrays between master and slaves
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
  integer(ip)              :: ipoin,ii,jj,bsize,ini, dom_i
  integer(4)               :: istat,bsize4
  real(rp)                 :: time1,time2
  integer(ip), allocatable :: loc_spari1(:),   loc_rpari1(:)

  call cputim(time1)

  if(kfl_paral>0) then

     if(party==1) then
        !
        ! Element
        !
     else if(party==2) then
        !
        ! Boundary
        !
     else if(party==3) then
        !
        ! Node
        !
        if(pardi==1.and.parki==1) then

           !-------------------------------------------------------------
           !
           ! INT(NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_spari1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_sltake',loc_spari1)

           allocate(loc_rpari1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_sltake',loc_rpari1)

           do jj = 1,commd%bound_dim
              ipoin = commd%bound_perm(jj)
              loc_spari1(jj) = pari1(ipoin)
           end do
           
           do ii = 1,nneig

              dom_i = commd%neights(ii)
              ini   = commd%bound_size(ii)
              bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
              bsize4=int(bsize,4)
              call MPI_Sendrecv( loc_spari1(ini:), bsize4,&
                   PAR_INTEGER,  dom_i, 0_4,     &
                   loc_rpari1(ini:), bsize4,              &
                   PAR_INTEGER, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           end do

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              if( pari1(ipoin) == 0 .and. loc_rpari1(jj) > 0 ) pari1(ipoin) = 2
           end do

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_sltake',loc_rpari1)
           deallocate(loc_rpari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARI1','par_sltake',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_sltake',loc_spari1)
           deallocate(loc_spari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARI1','par_sltake',0_ip)

        else if(pardi==1.and.parki==2) then

        else if(pardi==1.and.parki==5) then

        else if(pardi==2.and.parki==1) then

        else if(pardi==2.and.parki==2) then

        end if
        !
        ! Boundary node: REAL ARRAY(NBOPO)
        !
     else if(party==4) then
     end if
  endif

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_sltake
