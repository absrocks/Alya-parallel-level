subroutine par_slequa()
  !-----------------------------------------------------------------------
  !****f* Parall/par_slequa
  ! NAME
  !    par_slequa
  ! DESCRIPTION
  !    This subroutine exchange arrays between slaves
  !    commd%bound_perm(jj):      permutation array
  !    loc_sparr1:                my local values
  !    loc_rparr1:                values given by neighbor ii
  !    commd%bound_dim:           size of communication array
  !    nneig:                     number of neighbors that share same group
  !    commd%neights(ii):         number of subdomain ii
  !    commd%bound_size(ii):      where my local arrays sart to exchange with ii
  !    commd%bound_size(ii+1)
  !        -commd%bound_size(ii): number of groups to exchange with ii
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
  integer(ip)              :: ipoin,ii,jj,bsize,ini,dom_i
  integer(ip)              :: kk
  integer(4)               :: istat,bsize4
  real(rp)                 :: time1,time2
  integer(ip), allocatable :: loc_spari1(:),   loc_rpari1(:)
  real(rp),    allocatable :: loc_sparr1(:),   loc_rparr1(:)

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
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_slequa',loc_spari1)

           allocate(loc_rpari1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_slequa',loc_rpari1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              loc_spari1(jj) = pari1(ipoin)
           enddo

           do ii= 1, nneig
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
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              pari1(ipoin) = pari1(ipoin) + loc_rpari1(jj)
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_slequa',loc_rpari1)
           deallocate(loc_rpari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARI1','par_slequa',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_slequa',loc_spari1)
           deallocate(loc_spari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARI1','par_slequa',0_ip)

        else if(pardi==1.and.parki==2) then

           !-------------------------------------------------------------
           !
           ! REAL(NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slequa',loc_sparr1)

           allocate(loc_rparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slequa',loc_rparr1)

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
 
           !funin(1) = 0.0_rp
           !funin(2) = 0.0_rp
           !do jj= 1, commd%bound_dim
           !   ipoin = commd%bound_perm(jj)
           !   funin(1) = funin(1) + ( parr1(ipoin) - loc_rparr1(jj) ) ** 2
           !   funin(2) = funin(2) + ( loc_rparr1(jj) ) ** 2
           !enddo

           do ii= 1, nneig
              dom_i = commd%neights(ii)
              if( kfl_paral < dom_i ) then
                 do jj = commd%bound_size(ii),commd%bound_size(ii+1)-1
                    ipoin = commd%bound_perm(jj)
                    parr1(ipoin) = loc_rparr1(jj)
                 end do
              end if
           end do

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slequa',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slequa',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slequa',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slequa',0_ip)

        else if(pardi>=1.and.parki==5) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN) => REAL(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slequa',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slequa',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 loc_sparr1(pard1*(jj-1)+ii) = parr1(pard1*(ipoin-1)+ii)
              enddo
           enddo

           do ii= 1, nneig
              dom_i = commd%neights(ii)

              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   loc_rparr1(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           end do

           do ii= 1, nneig
              dom_i = commd%neights(ii)
              if( kfl_paral < dom_i ) then
                 do jj = commd%bound_size(ii),commd%bound_size(ii+1)-1
                    ipoin = commd%bound_perm(jj)
                    do kk = 1,pard1
                       parr1(pard1*(ipoin-1)+kk) = loc_rparr1(pard1*(jj-1)+kk)
                    end do
                 end do
              end if
           end do

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slequa',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slequa',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slequa',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slequa',0_ip)

        else if(pardi==2.and.parki==1) then

           call runend('PAR_SLEQUA: NOT CODED')

        else if(pardi==2.and.parki==2) then

           call runend('PAR_SLEQUA: NOT CODED')

        else if(pardi==1.and.parki==6) then

           call runend('OBSOLETE')
 
       else if ( pardi==1 .and. parki==4 ) then

           call runend('PAR_SLEQUA: NOT CODED')

	else if ( pardi == 1 .and. parki == 7 ) then

           call runend('PAR_SLEQUA: NOT CODED')

        end if

     else if(party==4) then

        if( parki == 2 .and. pardi == 1 ) then

           call runend('PAR_SLEQUA: NOT CODED')

        else if( parki == 2 .and. pardi == 2 ) then

           call runend('PAR_SLEQUA: NOT CODED')

        else if(pardi>=1.and.parki==5) then

           call runend('PAR_SLEQUA: NOT CODED')

        endif

     end if
  endif

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slequa
