subroutine par_slexma()
  !-----------------------------------------------------------------------
  !****f* Parall/par_slexma
  ! NAME
  !    par_slexma
  ! DESCRIPTION
  !    This subroutine exchange arrays between slaves and takes the max
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
  integer(ip)              :: ipoin,ii,jj,bsize,ji,poin,ini,dom_i,ibopo,ineig
  integer(4)               :: istat,bsize4
  real(rp)                 :: time1,time2
  integer(ip), allocatable :: loc_spari1(:),   loc_rpari1(:)
  real(rp),    allocatable :: loc_sparr1(:),   loc_rparr1(:)
  complex(rp), allocatable :: loc_sparx1(:),   loc_rparx1(:)

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
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_slexma',loc_spari1)

           allocate(loc_rpari1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_slexma',loc_rpari1)

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
              pari1(ipoin) = max(pari1(ipoin),loc_rpari1(jj))
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARI1','par_slexma',loc_rpari1)
           deallocate(loc_rpari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARI1','par_slexma',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARI1','par_slexma',loc_spari1)
           deallocate(loc_spari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARI1','par_slexma',0_ip)

        else if(pardi==1.and.parki==2) then

           !-------------------------------------------------------------
           !
           ! REAL(NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)

           allocate(loc_rparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)

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

           goto 10
           call memgen(1_ip,nneig+1,0_ip)
           allocate(parre(npoin))
           do ipoin = 1,npoin
              parre(ipoin) = parr1(ipoin)
           end do
           do ii = 1,nneig
              gisca(ii) = commd%neights(ii)
           end do
           gisca(nneig+1)=kfl_paral
           call heapsorti1(2_ip,nneig+1,gisca)

           do jj = 1,commd%bound_dim
              ipoin = commd%bound_perm(jj)
              parr1(ipoin) = 0.0_rp
           end do

           do ineig = 1,nneig+1
              dom_i = gisca(ineig)
              if( dom_i == kfl_paral ) then
                 do ipoin= npoi1+1,npoin
                    parr1(ipoin) = parr1(ipoin) + parre(ipoin)
                 end do
              else
                 ii = 1
                 do while( commd%neights(ii) /= dom_i )
                    ii = ii + 1
                 end do
                 do jj = commd%bound_size(ii),commd%bound_size(ii+1)-1
                    ipoin = commd%bound_perm(jj)
                    parr1(ipoin) = parr1(ipoin) + loc_rparr1(jj)
                 end do
              end if
           end do
           call memgen(3_ip,nneig+1,0_ip)
           deallocate(parre)
10         continue 
           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              parr1(ipoin) = max(parr1(ipoin),loc_rparr1(jj))
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexma',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexma',0_ip)

        else if(pardi>=1.and.parki==5) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN) => REAL(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)

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

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 parr1(pard1*(ipoin-1)+ii) = max(parr1(pard1*(ipoin-1)+ii),loc_rparr1(pard1*(jj-1)+ii))
              enddo
           enddo
           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexma',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexma',0_ip)

        else if(pardi==2.and.parki==1) then
        
           call runend('PAR_SLEXMA: NOT CODED')

        else if(pardi==2.and.parki==2) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 loc_sparr1(pard1*(jj-1)+ii) = parr2(ii,ipoin)
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
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 parr2(ii,ipoin) = max(parr2(ii,ipoin), loc_rparr1(pard1*(jj-1)+ii))
              enddo
           enddo
           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexma',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexma',0_ip)

        else if(pardi==1.and.parki==6) then

           call runend('OBSOLETE')

        else if ( pardi==1 .and. parki==4 ) then

           !-------------------------------------------------------------
           !
           ! COMPLEX(NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparx1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARX1','par_slexma',loc_sparx1)

           allocate(loc_rparx1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARX1','par_slexma',loc_rparx1)

           do jj = 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              loc_sparx1(jj) = parx1(ipoin)
           enddo

           do ii = 1, nneig
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
              !parx1(ipoin) = max(parx1(ipoin), loc_rparx1(jj))
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARX1','par_slexma',loc_rparx1)
           deallocate(loc_rparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARX1','par_slexma',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARX1','par_slexma',loc_sparx1)
           deallocate(loc_sparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARX1','par_slexma',0_ip)

        else if ( pardi == 1 .and. parki == 7 ) then

           !-------------------------------------------------------------
           !
           ! COMPLEX(PARD1,NPOIN) => COMPLEX(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparx1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARX1','par_slexma',loc_sparx1)

           allocate(loc_rparx1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARX1','par_slexma',loc_rparx1)

           do jj = 1,commd%bound_dim

              ji    = pard1 * (jj - 1)
              ipoin = commd%bound_perm(jj)
              poin  = pard1 * (ipoin - 1)
              do ii = 1,pard1

                 loc_sparx1(ji+ii) = parx1(poin+ii)

              enddo

           enddo

           do ii = 1,nneig

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

                 !parx1(poin+ii) = max(parx1(poin+ii), loc_rparx1(ji+ii))

              enddo

           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARX1','par_slexma',loc_rparx1)
           deallocate(loc_rparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARX1','par_slexma',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARX1','par_slexma',loc_sparx1)
           deallocate(loc_sparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARX1','par_slexma',0_ip)


        else

           call runend('PAR_SLEXMA: NOT ESTA PROGRAMADA VAGO')

        end if

     else if(party==4) then

        if( parki == 2 .and. pardi == 1 ) then

           !-------------------------------------------------------------
           !
           ! REAL(NBOPO)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)

           allocate(loc_rparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) &
                   loc_sparr1(jj) = parr1(ibopo)
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

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) &
                   parr1(ibopo) = max(parr1(ibopo), loc_rparr1(jj))
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexma',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexma',0_ip)

        else if( parki == 2 .and. pardi == 2 ) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NBOPO)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    loc_sparr1(pard1*(jj-1)+ii) = parr2(ii,ibopo)
                 enddo
              end if
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
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    parr2(ii,ibopo) = max(parr2(ii,ibopo),loc_rparr1(pard1*(jj-1)+ii))
                 enddo
              end if
           enddo

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexma',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexma',0_ip)

        else if(pardi>=1.and.parki==5) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NBOPO) => REAL(PARD1*NBOPO)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    loc_sparr1(pard1*(jj-1)+ii) = parr1(pard1*(ibopo-1)+ii)
                 enddo
              end if
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
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do ii= 1, pard1
                    parr1(pard1*(ibopo-1)+ii) = max(parr1(pard1*(ibopo-1)+ii),loc_rparr1(pard1*(jj-1)+ii))
                 enddo
              end if
           enddo
           call memchk(two,istat,mem_servi(1:2,servi),'LOC_RPARR1','par_slexma',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexma',0_ip)

           call memchk(two,istat,mem_servi(1:2,servi),'LOC_SPARR1','par_slexma',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexma',0_ip)

        endif

     end if
  endif

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slexma
