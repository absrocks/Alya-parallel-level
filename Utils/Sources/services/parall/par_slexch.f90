subroutine par_slexch()
  !-----------------------------------------------------------------------
  !****f* Parall/par_slexch
  ! NAME
  !    par_slexch
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
  use mod_maths,  only : maths_heap_sort
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_parall, only : commd,PAR_COMM_MY_CODE4
  use mod_parall, only : PAR_INTEGER
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
  integer(4)               :: status(MPI_STATUS_SIZE)
#endif
  integer(ip)              :: ipoin,ii,jj,bsize,ji,poin,ini,dom_i,ibopo,ineig
  integer(ip)              :: kpoin,kk,kfl_order
  integer(4)               :: istat,bsize4
  real(rp)                 :: time1,time2
  integer(ip), allocatable :: loc_spari1(:),   loc_rpari1(:)
  real(rp),    allocatable :: loc_sparr1(:),   loc_rparr1(:)
  complex(rp), allocatable :: loc_sparx1(:),   loc_rparx1(:)
  real(rp),    pointer     :: my_rparr1(:)
  integer(ip), pointer     :: neights_ordered(:)
  integer(ip), pointer     :: perm_ordered(:)
  
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

        else if( pardi >= 1 .and. parki == 6 ) then

           !-------------------------------------------------------------
           !
           ! INT(PARD1,NPOIN) => INT(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_spari1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_slexch',loc_spari1)

           allocate(loc_rpari1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_slexch',loc_rpari1)

           do jj = 1,commd%bound_dim
              ipoin = commd%bound_perm(jj)
              kpoin = pard1*(ipoin-1)
              kk    = pard1*(jj-1)
              do ii = 1,pard1
                 kk    = kk + 1
                 kpoin = kpoin + 1
                 loc_spari1(kk) = pari1(kpoin)
              enddo
           enddo

           do ii = 1,commd%nneig

              dom_i = commd%neights(ii)
              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( loc_spari1(ini:), bsize4,&
                   PAR_INTEGER,  dom_i, 0_4,     &
                   loc_rpari1(ini:), bsize4,              &
                   PAR_INTEGER, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           end do

           do jj = 1,commd%bound_dim
              ipoin = commd%bound_perm(jj)
              kpoin = pard1*(ipoin-1)
              kk    =  pard1*(jj-1)
              do ii = 1,pard1
                 kk    = kk + 1
                 kpoin = kpoin + 1
                 pari1(kpoin) = pari1(kpoin) + loc_rpari1(kk)
              enddo
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
           ! REAL(NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)

           allocate(loc_rparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)

           kfl_order = 0
           
           if( kfl_order == 1 ) then

              nullify(my_rparr1)
              nullify(neights_ordered)
              nullify(perm_ordered)
              call memory_alloca(mem_servi(1:2,servi),'MY_RPARR1'      ,'par_slexch',my_rparr1      ,npoin-npoi1  ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'NEIGHTS_ORDERED','par_slexch',neights_ordered,commd % nneig,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'PERM_ORDERED'   ,'par_slexch',perm_ordered   ,commd % nneig,'DO_NOT_INITIALIZE')
              do ipoin = npoi1+1,npoin
                 my_rparr1(ipoin-npoi1) = parr1(ipoin)
              end do
              neights_ordered(1:commd % nneig) = commd % neights(1:commd % nneig)
              do ii = 1,commd % nneig
                 perm_ordered(ii) = ii
              end do
              call maths_heap_sort(2_ip,commd % nneig,neights_ordered,'NORMAL',perm_ordered)
           end if
           
           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              loc_sparr1(jj) = parr1(ipoin)
           enddo

           do ii= 1, commd%nneig
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

           if( kfl_order == 1 ) then 
              do ipoin = npoi1+1,npoin
                 parr1(ipoin) = 0.0_rp 
              end do

              ii = 0
              dom_i = neights_ordered(1)
              
              do while( dom_i < kfl_paral .and. ii < commd % nneig )
                 ii    = ii + 1
                 dom_i = neights_ordered(ii)
                 kk    = perm_ordered(ii)

                 do jj = commd % bound_size(kk),commd % bound_size(kk+1)-1
                    ipoin        = commd % bound_perm(jj)
                    parr1(ipoin) = parr1(ipoin) + loc_rparr1(jj)
                 end do                 
              end do
           
              do ipoin = npoi1+1,npoin
                 parr1(ipoin) = parr1(ipoin) + my_rparr1(ipoin-npoi1) 
              end do              
            
              do while( ii < commd % nneig )
                 ii    = ii + 1
                 dom_i = neights_ordered(ii)
                 kk    = perm_ordered(ii)
                 do jj = commd % bound_size(kk),commd % bound_size(kk+1)-1
                    ipoin        = commd % bound_perm(jj)
                    parr1(ipoin) = parr1(ipoin) + loc_rparr1(jj)
                 end do                 
              end do

              call memory_deallo(mem_servi(1:2,servi),'PERM_ORDERED'   ,'par_slexch',perm_ordered)
              call memory_deallo(mem_servi(1:2,servi),'NEIGHTS_ORDERED','par_slexch',neights_ordered)
              call memory_deallo(mem_servi(1:2,servi),'MY_RPARR1'      ,'par_slexch',my_rparr1)

           else
              do jj= 1, commd%bound_dim
                 ipoin = commd%bound_perm(jj)
                 parr1(ipoin) = parr1(ipoin) + loc_rparr1(jj)
              end do
           end if
           
           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

        else if(pardi>=1.and.parki==5) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN) => REAL(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 loc_sparr1(pard1*(jj-1)+ii) = parr1(pard1*(ipoin-1)+ii)
              enddo
           enddo

           do ii= 1, commd%nneig
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

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 parr1(pard1*(ipoin-1)+ii) = parr1(pard1*(ipoin-1)+ii) + loc_rparr1(pard1*(jj-1)+ii)
              enddo
           enddo
           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

        else if(pardi==2.and.parki==1) then

        else if(pardi==2.and.parki==2) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 loc_sparr1(pard1*(jj-1)+ii) = parr2(ii,ipoin)
              enddo
           enddo

           do ii= 1, commd%nneig
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
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 parr2(ii,ipoin) = parr2(ii,ipoin) + loc_rparr1(pard1*(jj-1)+ii)
              enddo
           enddo
           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

        else if ( pardi==1 .and. parki==4 ) then

           !-------------------------------------------------------------
           !
           ! COMPLEX(NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparx1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARX1','par_slexch',loc_sparx1)

           allocate(loc_rparx1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARX1','par_slexch',loc_rparx1)

           do jj = 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              loc_sparx1(jj) = parx1(ipoin)
           enddo

           do ii = 1, commd%nneig
              dom_i = commd%neights(ii)

              ini   = commd%bound_size(ii)
              bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( loc_sparx1(ini:), bsize4,&
                   MPI_DOUBLE_COMPLEX, dom_i, 0_4,        &
                   loc_rparx1(ini:), bsize4,              &
                   MPI_DOUBLE_COMPLEX, dom_i, 0_4,        &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           enddo

           do jj = 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              parx1(ipoin) = parx1(ipoin) + loc_rparx1(jj)
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARX1','par_slexch',loc_rparx1)
           deallocate(loc_rparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARX1','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARX1','par_slexch',loc_sparx1)
           deallocate(loc_sparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARX1','par_slexch',0_ip)

        else if ( pardi == 1 .and. parki == 7 ) then

           !-------------------------------------------------------------
           !
           ! COMPLEX(PARD1,NPOIN) => COMPLEX(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparx1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARX1','par_slexch',loc_sparx1)

           allocate(loc_rparx1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARX1','par_slexch',loc_rparx1)

           do jj = 1,commd%bound_dim

              ji    = pard1 * (jj - 1)
              ipoin = commd%bound_perm(jj)
              poin  = pard1 * (ipoin - 1)
              do ii = 1,pard1

                 loc_sparx1(ji+ii) = parx1(poin+ii)

              enddo

           enddo

           do ii = 1,commd%nneig

              dom_i = commd%neights(ii)

              ini   = pard1 * (commd%bound_size(ii) - 1) + 1
              bsize = pard1 * (commd%bound_size(ii+1) - 1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv(loc_sparx1(ini:),bsize4,&
                   MPI_DOUBLE_COMPLEX,dom_i,0_4,        &
                   loc_rparx1(ini:),bsize4,             &
                   MPI_DOUBLE_COMPLEX,dom_i,0_4,        &
                   PAR_COMM_MY_CODE4,status,istat)
#endif
           enddo

           do jj = 1,commd%bound_dim

              ji = pard1 * (jj - 1)
              ipoin = commd%bound_perm(jj)
              poin = pard1 * (ipoin - 1)
              do ii = 1,pard1

                 parx1(poin+ii) = parx1(poin+ii) + loc_rparx1(ji+ii)

              enddo

           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARX1','par_slexch',loc_rparx1)
           deallocate(loc_rparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARX1','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARX1','par_slexch',loc_sparx1)
           deallocate(loc_sparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARX1','par_slexch',0_ip)


        end if

     else if(party==4) then

        if( parki == 2 .and. pardi == 1 ) then

           !-------------------------------------------------------------
           !
           ! REAL(NBOPO)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)

           allocate(loc_rparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) &
                   loc_sparr1(jj) = parr1(ibopo)
           enddo

           do ii= 1, commd%nneig
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

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) &
                   parr1(ibopo) = parr1(ibopo) + loc_rparr1(jj)
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

        else if( parki == 2 .and. pardi == 2 ) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NBOPO)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    loc_sparr1(pard1*(jj-1)+ii) = parr2(ii,ibopo)
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
              call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   loc_rparr1(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    parr2(ii,ibopo) = parr2(ii,ibopo) + loc_rparr1(pard1*(jj-1)+ii)
                 enddo
              end if
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

        else if(pardi>=1.and.parki==5) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NBOPO) => REAL(PARD1*NBOPO)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    loc_sparr1(pard1*(jj-1)+ii) = parr1(pard1*(ibopo-1)+ii)
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
              call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   loc_rparr1(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    parr1(pard1*(ibopo-1)+ii) = parr1(pard1*(ibopo-1)+ii) + loc_rparr1(pard1*(jj-1)+ii)
                 enddo
              end if
           enddo
           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

        endif

     end if
  endif

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slexch
