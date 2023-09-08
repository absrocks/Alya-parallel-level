subroutine par_slexca
  !-----------------------------------------------------------------------
  !****f* Parall/par_slexca
  ! NAME
  !    par_slexch
  ! DESCRIPTION
  !    This subroutine exchange arrays between master and slaves
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
  use mod_parall, only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall, only : commd,PAR_COMM_MY_CODE4
  use mod_parall, only : PAR_INTEGER
  use mod_memory
  use mod_memchk
  
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
  integer(4)                 :: status(MPI_STATUS_SIZE,2*commd % nneig)
#endif
  integer(ip)                :: ipoin,ii,jj,kk,bsize,ini,dom_i,ibopo
  integer(ip)                :: ndofi,ndofj
  integer(ip), save          :: count4
  integer(4)                 :: istat,bsize4
  real(rp)                   :: time1,time2
  integer(ip), allocatable   :: loc_spari1(:), loc_rpari1(:)
  integer(4)                 :: ierr4

  call cputim(time1)

  if( ISLAVE ) then

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
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_slexch',loc_spari1)

           allocate(loc_rpari1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_slexch',loc_rpari1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              loc_spari1(jj) = pari1(ipoin)
           enddo

           do ii= 1, commd%nneig
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

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_slexch',loc_rpari1)
           deallocate(loc_rpari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARI1','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_slexch',loc_spari1)
           deallocate(loc_spari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARI1','par_slexch',0_ip)

        else if(pardi==1.and.parki==2) then

           !-------------------------------------------------------------
           !
           ! REAL(NPOIN): AQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
           !
           !-------------------------------------------------------------
           
           if( ipass_par == 0 ) then

              ipass_par = 1
              allocate(ireq4(commd%nneig*2),stat=istat)
              call memchk(zero,istat,mem_servi(1:2,servi),'IREQ4','par_slexch',ireq4)

              allocate(parws(commd%bound_dim),stat=istat)
              call memchk(zero,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)

              allocate(parwr(commd%bound_dim),stat=istat)
              call memchk(zero,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)

              do jj= 1, commd%bound_dim
                 ipoin = commd%bound_perm(jj)
                 parws(jj) = parr1(ipoin)
              enddo

              kk = 0
              do ii = 1,commd % nneig
                 dom_i = commd % neights(ii)

                 ini   = commd % bound_size(ii)
                 bsize = commd % bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
                 
                 bsize4=int(bsize,4)
                 kk = kk + 1
                 call MPI_Isend( parws(ini:ini+bsize-1), bsize4, &
                      MPI_DOUBLE_PRECISION,  dom_i, 0_4,   &
                      PAR_COMM_MY_CODE4, ireq4(kk), istat )
                 kk = kk + 1
                 call MPI_Irecv( parwr(ini:ini+bsize-1), bsize4, &
                      MPI_DOUBLE_PRECISION,  dom_i, 0_4,   &
                      PAR_COMM_MY_CODE4, ireq4(kk), istat )
#endif
              enddo

           else


              ipass_par  = 0
              count4 = commd % nneig*2
#ifdef MPI_OFF
#else
              CALL MPI_WAITALL(count4,ireq4,status,ierr4)
#endif

              do jj = 1,commd % bound_dim
                 ipoin = commd % bound_perm(jj)
                 parr1(ipoin) = parr1(ipoin) + parwr(jj)
              enddo

              call memchk(two,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)
              deallocate(parwr,stat=istat)
              if(istat/=0) call memerr(two,'PARWR','par_slexch',0_ip)

              call memchk(two,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)
              deallocate(parws,stat=istat)
              if(istat/=0) call memerr(two,'PARWS','par_slexch',0_ip)

              call memchk(two,istat,mem_servi(1:2,servi),'IREQ4','par_slexch',ireq4)
              deallocate(ireq4,stat=istat)
              if(istat/=0) call memerr(two,'IREQ4','par_slexch',0_ip)

           end if


        else if(pardi>=1.and.parki==5) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN) => REAL(PARD1*NPOIN) AQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
           !
           !-------------------------------------------------------------

           if( ipass_par == 0 ) then

              ipass_par = 1

              allocate(ireq4(commd%nneig*2),stat=istat)
              call memchk(zero,istat,mem_servi(1:2,servi),'IREQ4','par_slexch',ireq4)
              
              allocate(parws(pard1*commd%bound_dim),stat=istat)
              call memchk(zero,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)
              
              allocate(parwr(pard1*commd%bound_dim),stat=istat)
              call memchk(zero,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)
              
              do jj= 1, commd%bound_dim
                 ipoin = commd%bound_perm(jj)
                 ndofi = pard1*(ipoin-1)
                 ndofj = pard1*(jj-1)
                 do ii= 1, pard1
                    ndofi = ndofi + 1
                    ndofj = ndofj + 1
                    parws(ndofj) = parr1(ndofi)
                 enddo
              enddo

              kk = 0
              do ii= 1, commd%nneig
                 dom_i = commd%neights(ii)
                 
                 ini   = pard1*(commd%bound_size(ii)-1) + 1
                 bsize = pard1*(commd%bound_size(ii+1)-commd%bound_size(ii))
                 
#ifdef MPI_OFF
#else
                 bsize4 = int(bsize,4)
                 kk = kk + 1
                 call MPI_Isend( parws(ini:ini+bsize-1), bsize4, &
                      MPI_DOUBLE_PRECISION,  dom_i, 0_4,   &
                      PAR_COMM_MY_CODE4, ireq4(kk), istat )
                 kk = kk + 1
                 call MPI_Irecv( parwr(ini:ini+bsize-1), bsize4, &
                      MPI_DOUBLE_PRECISION,  dom_i, 0_4,   &
                      PAR_COMM_MY_CODE4, ireq4(kk), istat )              
#endif
              end do

           else

              ipass_par  = 0
              count4 = commd%nneig*2

#ifdef MPI_OFF
#else
              CALL MPI_WAITALL(count4,ireq4,status,ierr4)
#endif
              do jj = 1,commd%bound_dim
                 ipoin = commd%bound_perm(jj)
                 ndofi = pard1*(ipoin-1)
                 ndofj = pard1*(jj-1)
                 do ii = 1,pard1
                    ndofi = ndofi + 1
                    ndofj = ndofj + 1
                    parr1(ndofi) = parr1(ndofi) + parwr(ndofj)
                 end do
              end do

              call memchk(two,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)
              deallocate(parwr,stat=istat)
              if(istat/=0) call memerr(two,'PARWR','par_slexch',0_ip)

              call memchk(two,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)
              deallocate(parws,stat=istat)
              if(istat/=0) call memerr(two,'PARWS','par_slexch',0_ip)

              call memchk(two,istat,mem_servi(1:2,servi),'IREQ4','par_slexch',ireq4)
              deallocate(ireq4,stat=istat)
              if(istat/=0) call memerr(two,'IREQ4','par_slexch',0_ip)

           end if

        else if(pardi==2.and.parki==1) then

        else if(pardi==2.and.parki==2) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN)
           !
           !-------------------------------------------------------------

           allocate(parws(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)

           allocate(parwr(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 parws(pard1*(jj-1)+ii) = parr2(ii,ipoin)
              enddo
           enddo

           do ii= 1, commd%nneig
              dom_i = commd%neights(ii)

              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( parws(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   parwr(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 parr2(ii,ipoin) = parr2(ii,ipoin) + parwr(pard1*(jj-1)+ii)
              enddo
           enddo
           call memchk(two,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)
           deallocate(parwr,stat=istat)
           if(istat/=0) call memerr(two,'PARWR','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)
           deallocate(parws,stat=istat)
           if(istat/=0) call memerr(two,'PARWS','par_slexch',0_ip)

        else if(pardi==1.and.parki==6) then

           call runend('OBSOLETE')

        end if

     else if(party==4) then

        if( parki == 2 .and. pardi == 1 ) then

           !-------------------------------------------------------------
           !
           ! REAL(NBOPO)
           !
           !-------------------------------------------------------------

           allocate(parws(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)

           allocate(parwr(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) &
                   parws(jj) = parr1(ibopo)
           enddo

           do ii= 1, commd%nneig
              dom_i = commd%neights(ii)

              ini   = commd%bound_size(ii)
              bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
              bsize4=int(bsize,4)
              call MPI_Sendrecv( parws(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   parwr(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) &
                   parr1(ibopo) = parr1(ibopo) + parwr(jj)
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)
           deallocate(parwr,stat=istat)
           if(istat/=0) call memerr(two,'PARWR','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)
           deallocate(parws,stat=istat)
           if(istat/=0) call memerr(two,'PARWS','par_slexch',0_ip)

        else if( parki == 2 .and. pardi == 2 ) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NBOPO)
           !
           !-------------------------------------------------------------

           allocate(parws(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)

           allocate(parwr(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    parws(pard1*(jj-1)+ii) = parr2(ii,ibopo)
                 enddo
              end if
           enddo

           do ii= 1, commd%nneig
              dom_i = commd%neights(ii)

              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( parws(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   parwr(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    parr2(ii,ibopo) = parr2(ii,ibopo) + parwr(pard1*(jj-1)+ii)
                 enddo
              end if
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)
           deallocate(parwr,stat=istat)
           if(istat/=0) call memerr(two,'PARWR','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)
           deallocate(parws,stat=istat)
           if(istat/=0) call memerr(two,'PARWS','par_slexch',0_ip)

        else if(pardi>=1.and.parki==5) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NBOPO) => REAL(PARD1*NBOPO)
           !
           !-------------------------------------------------------------

           allocate(parws(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)

           allocate(parwr(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    parws(pard1*(jj-1)+ii) = parr1(pard1*(ibopo-1)+ii)
                 enddo
              end if
           enddo

           do ii= 1, commd%nneig
              dom_i = commd%neights(ii)

              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( parws(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   parwr(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    parr1(pard1*(ibopo-1)+ii) = parr1(pard1*(ibopo-1)+ii) + parwr(pard1*(jj-1)+ii)
                 enddo
              end if
           enddo
           call memchk(two,istat,mem_servi(1:2,servi),'PARWR','par_slexch',parwr)
           deallocate(parwr,stat=istat)
           if(istat/=0) call memerr(two,'PARWR','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'PARWS','par_slexch',parws)
           deallocate(parws,stat=istat)
           if(istat/=0) call memerr(two,'PARWS','par_slexch',0_ip)

        endif

     end if
  endif

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slexca
