subroutine par_slegro
  !-----------------------------------------------------------------------
  !****f* Parall/par_slegro
  ! NAME
  !    par_slegro
  ! DESCRIPTION
  !    This subroutine exchange arrays between master and slaves
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
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
  integer(4)            :: status(MPI_STATUS_SIZE)
#endif
  integer(ip)           :: ipoin
  integer(ip)           :: ii,jj,bsize,ini, dom_i
  integer(4)            :: istat,bsize4,bsizepard14
  real(rp)              :: time1,time2
  real(rp), allocatable :: loc_sparr1(:),   loc_rparr1(:)
  real(rp), allocatable :: loc_sparr2(:,:), loc_rparr2(:,:)

  call cputim(time1)

  

  if(kfl_paral>0) then

     if(party==3) then
        !
        ! Node
        !
        if(pardi==1.and.parki==2) then
           !
           ! comle(icoml)%commd%bound_perm(jj):               permutation array
           ! loc_sparr1:                                      my local values
           ! loc_rparr1:                                      values given by neighbor ii
           ! comle(icoml)%commd%bound_dim:                    size of communication array
           ! nneig:                                           number of neighbors that share same group
           ! comle(icoml)%commd%neights(ii):                  number of subdomain ii
           ! comle(icoml)%commd%bound_size(ii):               where my local arrays sart to exchange with ii
           ! comle(icoml)%commd%bound_size(ii+1)-comle(icoml)%commd%bound_size(ii): number of groups to exchange with ii
           !
           allocate(loc_sparr1(comle(icoml)%commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slegro',loc_sparr1)

           allocate(loc_rparr1(comle(icoml)%commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slegro',loc_rparr1)

           do jj= 1, comle(icoml)%commd%bound_dim
              ipoin = comle(icoml)%commd%bound_perm(jj)
              loc_sparr1(jj) = parr1(ipoin)
           enddo

           do ii= 1, comle(icoml)%nneig
              dom_i = comle(icoml)%commd%neights(ii)

              ini   = comle(icoml)%commd%bound_size(ii)
              bsize = comle(icoml)%commd%bound_size(ii+1) - ini

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

           do jj= 1, comle(icoml)%commd%bound_dim
              ipoin = comle(icoml)%commd%bound_perm(jj)
              parr1(ipoin) = parr1(ipoin) + loc_rparr1(jj)
           enddo

           call memchk(two,istat,par_memor,'LOC_RPARR1','par_slegro',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slegro',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR1','par_slegro',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slegro',0_ip)

        else if(pardi==1.and.parki==5) then

           allocate(loc_sparr1(pard1*comle(icoml)%commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slegro',loc_sparr1)

           allocate(loc_rparr1(pard1*comle(icoml)%commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slegro',loc_rparr1)

           do jj= 1, comle(icoml)%commd%bound_dim
              ipoin = comle(icoml)%commd%bound_perm(jj)
              do ii= 1, pard1
                 loc_sparr1(pard1*(jj-1)+ii) = parr1(pard1*(ipoin-1)+ii)
              enddo
           enddo

           do ii= 1, comle(icoml)%nneig
              dom_i = comle(icoml)%commd%neights(ii)

              ini   = pard1*(comle(icoml)%commd%bound_size(ii)-1)   + 1
              bsize = pard1*(comle(icoml)%commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   loc_rparr1(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           enddo

           do jj= 1, comle(icoml)%commd%bound_dim
              ipoin = comle(icoml)%commd%bound_perm(jj)
              do ii= 1, pard1
                 parr1(pard1*(ipoin-1)+ii) = parr1(pard1*(ipoin-1)+ii) + loc_rparr1(pard1*(jj-1)+ii)
              enddo
           enddo
           call memchk(two,istat,par_memor,'LOC_RPARR1','par_slegro',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slegro',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR1','par_slegro',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slegro',0_ip)

        else if(pardi==2.and.parki==1) then

        else if(pardi==2.and.parki==2) then

           allocate(loc_sparr2(pard1,comle(icoml)%commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR2','par_slegro',loc_sparr2)

           allocate(loc_rparr2(pard1,comle(icoml)%commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR2','par_slegro',loc_rparr2)

           do jj= 1, comle(icoml)%commd%bound_dim
              ipoin = comle(icoml)%commd%bound_perm(jj)
              loc_sparr2(1:pard1,jj) = parr2(1:pard1,ipoin)
           enddo

           do ii= 1, comle(icoml)%nneig
              dom_i = comle(icoml)%commd%neights(ii)

              ini   = comle(icoml)%commd%bound_size(ii)
              bsize = comle(icoml)%commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
              bsize4      = int(bsize,4)
              bsizepard14 = bsize4*int(pard1,4)
              call MPI_Sendrecv( loc_sparr2(1:,ini), bsizepard14,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,            &
                   loc_rparr2(1:,ini), bsizepard14,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,             &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           enddo

           do jj= 1, comle(icoml)%commd%bound_dim
              ipoin = comle(icoml)%commd%bound_perm(jj)
              parr2(1:pard1,ipoin) = parr2(1:pard1,ipoin) + loc_rparr2(1:pard1,jj)
           enddo

           call memchk(two,istat,par_memor,'LOC_RPARR2','par_slegro',loc_rparr2)
           deallocate(loc_rparr2,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR2','par_slegro',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR2','par_slegro',loc_sparr2)
           deallocate(loc_sparr2,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR2','par_slegro',0_ip)
        end if
 
     end if
  endif

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slegro
