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
  use mod_memory, only : memory_size
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_parall, only : commd,PAR_COMM_MY_CODE4
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : par_memor
  use mod_parall, only : sendbuff_rp
  use mod_parall, only : recvbuff_rp
  use mod_parall, only : sendbuff_ip
  use mod_parall, only : recvbuff_ip
  use mod_communications_tools, only : PAR_MPI_RUNEND
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
  integer(4)               :: status(MPI_STATUS_SIZE)
#endif
  integer(ip)              :: ipoin,ii,jj,bsize,ji,poin,ini,dom_i,ibopo
  integer(ip)              :: kpoin,kk
  integer(4)               :: istat,bsize4
  real(rp)                 :: time1,time2
  complex(rp), allocatable :: loc_sparx1(:),   loc_rparx1(:)
  real(rp),    pointer     :: my_rparr1(:)
  
  call cputim(time1)

  if( ISLAVE ) then

     if(party==3) then
        !
        ! Node
        !
        if(pardi==1.and.parki==1) then

           !-------------------------------------------------------------
           !
           ! INT(NPOIN)
           !
           !-------------------------------------------------------------

           if( memory_size(sendbuff_ip) <  commd%bound_dim ) then
              call memory_deallo(par_memor,'SENDBUFF_IP','par_slexca',sendbuff_ip)
              call memory_alloca(par_memor,'SENDBUFF_IP','par_slexca',sendbuff_ip,commd%bound_dim,'DO_NOT_INITIALIZE')
           end if
           if( memory_size(recvbuff_ip) <  commd%bound_dim ) then
              call memory_deallo(par_memor,'RECVBUFF_IP','par_slexca',recvbuff_ip)
              call memory_alloca(par_memor,'RECVBUFF_IP','par_slexca',recvbuff_ip,commd%bound_dim,'DO_NOT_INITIALIZE')
           end if

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              sendbuff_ip(jj) = pari1(ipoin)
           enddo

           do ii= 1, commd%nneig
              dom_i = commd%neights(ii)

              ini   = commd%bound_size(ii)
              bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
              bsize4=int(bsize,4)
              call MPI_Sendrecv( sendbuff_ip(ini:), bsize4,&
                   PAR_INTEGER,  dom_i, 0_4,     &
                   recvbuff_ip(ini:), bsize4,              &
                   PAR_INTEGER, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
              if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCH')
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_invp(jj)
              pari1(ipoin) = pari1(ipoin) + recvbuff_ip(jj)
           enddo


        else if( pardi >= 1 .and. parki == 6 ) then

           !-------------------------------------------------------------
           !
           ! INT(PARD1,NPOIN) => INT(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           if( memory_size(sendbuff_ip) <  commd%bound_dim ) then
              call memory_deallo(par_memor,'SENDBUFF_IP','par_slexca',sendbuff_ip)
              call memory_alloca(par_memor,'SENDBUFF_IP','par_slexca',sendbuff_ip,commd%bound_dim,'DO_NOT_INITIALIZE')
           end if
           if( memory_size(recvbuff_ip) <  commd%bound_dim ) then
              call memory_deallo(par_memor,'RECVBUFF_IP','par_slexca',recvbuff_ip)
              call memory_alloca(par_memor,'RECVBUFF_IP','par_slexca',recvbuff_ip,commd%bound_dim,'DO_NOT_INITIALIZE')
           end if

           do jj = 1,commd%bound_dim
              ipoin = commd%bound_perm(jj)
              kpoin = pard1*(ipoin-1)
              kk    = pard1*(jj-1)
              do ii = 1,pard1
                 kk    = kk + 1
                 kpoin = kpoin + 1
                 sendbuff_ip(kk) = pari1(kpoin)
              enddo
           enddo

           do ii = 1,commd%nneig

              dom_i = commd%neights(ii)
              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( sendbuff_ip(ini:), bsize4,&
                   PAR_INTEGER,  dom_i, 0_4,     &
                   recvbuff_ip(ini:), bsize4,              &
                   PAR_INTEGER, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
              if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCH')

#endif
           end do

           do jj = 1,commd%bound_dim
              ipoin = commd%bound_invp(jj)
              kpoin = pard1*(ipoin-1)
              kk    =  pard1*(jj-1)
              do ii = 1,pard1
                 kk    = kk + 1
                 kpoin = kpoin + 1
                 pari1(kpoin) = pari1(kpoin) + recvbuff_ip(kk)
              enddo
           enddo

        else if(pardi==1.and.parki==2) then

           !-------------------------------------------------------------
           !
           ! REAL(NPOIN)
           !
           !-------------------------------------------------------------

           if( memory_size(sendbuff_rp) <  commd%bound_dim ) then
              call memory_deallo(par_memor,'SENDBUFF_RP','par_slexca',sendbuff_rp)
              call memory_alloca(par_memor,'SENDBUFF_RP','par_slexca',sendbuff_rp,commd%bound_dim,'DO_NOT_INITIALIZE')
           end if
           if( memory_size(recvbuff_rp) <  commd%bound_dim ) then
              call memory_deallo(par_memor,'RECVBUFF_RP','par_slexca',recvbuff_rp)
              call memory_alloca(par_memor,'RECVBUFF_RP','par_slexca',recvbuff_rp,commd%bound_dim,'DO_NOT_INITIALIZE')
           end if
           
           if( kfl_order_exchange_par == 1 .and. commd % npoin-commd % npoi1 > 0 ) then

              nullify(my_rparr1)
              allocate(my_rparr1(commd % npoin-commd % npoi1))
              do ipoin = commd % npoi1+1,commd % npoin
                 my_rparr1(ipoin-commd % npoi1) = parr1(ipoin)
              end do  
           end if
           
           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              sendbuff_rp(jj) = parr1(ipoin)
           enddo          
           
           do ii= 1, commd%nneig
              dom_i = commd%neights(ii)

              ini   = commd%bound_size(ii)
              bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
              bsize4=int(bsize,4)
              call MPI_Sendrecv( sendbuff_rp(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   recvbuff_rp(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
              if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCH')

#endif
           enddo

           if( kfl_order_exchange_par == 1 .and. commd % npoin-commd % npoi1 > 0 ) then

              do ipoin = commd % npoi1+1,commd % npoin
                 parr1(ipoin) = 0.0_rp 
              end do

              do ii = 1,commd % nneig_1
                 dom_i = commd%neights_ordered(ii)
                 kk    = commd%perm_ordered(ii)                
                 do jj = commd % bound_size(kk),commd % bound_size(kk+1)-1
                    ipoin        = commd % bound_invp(jj)
                    parr1(ipoin) = parr1(ipoin) + recvbuff_rp(jj)
                 end do                 
              end do
           
              do ipoin = commd % npoi1+1,commd % npoin
                 parr1(ipoin) = parr1(ipoin) + my_rparr1(ipoin-commd % npoi1) 
              end do              
            
               do ii = commd % nneig_2,commd % nneig
                 dom_i = commd%neights_ordered(ii)
                 kk    = commd%perm_ordered(ii)
                 do jj = commd % bound_size(kk),commd % bound_size(kk+1)-1
                    ipoin        = commd % bound_invp(jj)
                    parr1(ipoin) = parr1(ipoin) + recvbuff_rp(jj)
                 end do                 
              end do

              deallocate(my_rparr1)

           else
              do jj= 1, commd%bound_dim
                 ipoin = commd%bound_invp(jj)
                 parr1(ipoin) = parr1(ipoin) + recvbuff_rp(jj)
              end do
           end if
           
        else if(pardi>=1.and.parki==5) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN) => REAL(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           if( memory_size(sendbuff_rp) <  pard1*commd%bound_dim ) then
              call memory_deallo(par_memor,'SENDBUFF_RP','par_slexca',sendbuff_rp)
              call memory_alloca(par_memor,'SENDBUFF_RP','par_slexca',sendbuff_rp,pard1*commd%bound_dim,'DO_NOT_INITIALIZE')
           end if
           if( memory_size(recvbuff_rp) <  pard1*commd%bound_dim ) then
              call memory_deallo(par_memor,'RECVBUFF_RP','par_slexca',recvbuff_rp)
              call memory_alloca(par_memor,'RECVBUFF_RP','par_slexca',recvbuff_rp,pard1*commd%bound_dim,'DO_NOT_INITIALIZE')
           end if
           
           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 sendbuff_rp(pard1*(jj-1)+ii) = parr1(pard1*(ipoin-1)+ii)
              enddo
           enddo
           
           if( kfl_order_exchange_par == 1 .and. commd % npoin-commd % npoi1 > 0 ) then
              nullify(my_rparr1)
              allocate(my_rparr1(pard1*(commd % npoin-commd % npoi1)))
              do ipoin = commd % npoi1+1,commd % npoin
                 do kk = 1,pard1
                    my_rparr1((ipoin-commd % npoi1-1)*pard1+kk) = parr1((ipoin-1)*pard1+kk)
                 end do
              end do  
           end if

           do ii= 1, commd%nneig
              dom_i = commd%neights(ii)

              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( sendbuff_rp(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   recvbuff_rp(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
              if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCH')
#endif
           end do

           if( kfl_order_exchange_par == 1 .and. commd % npoin-commd % npoi1 > 0 ) then

              do ipoin = (commd % npoi1)*pard1+1,commd % npoin*pard1
                 parr1(ipoin) = 0.0_rp 
              end do

              do ii = 1,commd % nneig_1
                 dom_i = commd%neights_ordered(ii)
                 kk    = commd%perm_ordered(ii)                
                 do jj = commd % bound_size(kk),commd % bound_size(kk+1)-1
                    ipoin        = commd % bound_invp(jj)
                    do kk = 1,pard1
                       parr1((ipoin-1)*pard1+kk) = parr1((ipoin-1)*pard1+kk) + recvbuff_rp((jj-1)*pard1+kk)
                    end do
                 end do                 
              end do
           
              do ipoin = commd % npoi1+1,commd % npoin
                 do kk = 1,pard1
                    parr1((ipoin-1)*pard1+kk) = parr1((ipoin-1)*pard1+kk) + my_rparr1((ipoin-commd % npoi1-1)*pard1+kk)
                 end do
              end do              
            
               do ii = commd % nneig_2,commd % nneig
                 dom_i = commd%neights_ordered(ii)
                 kk    = commd%perm_ordered(ii)
                 do jj = commd % bound_size(kk),commd % bound_size(kk+1)-1
                    ipoin        = commd % bound_invp(jj)
                    do kk = 1,pard1
                       parr1((ipoin-1)*pard1+kk) = parr1((ipoin-1)*pard1+kk) + recvbuff_rp((jj-1)*pard1+kk)
                    end do
                 end do                 
              end do

              deallocate(my_rparr1)

           else
              do jj= 1, commd%bound_dim
                 ipoin = commd%bound_invp(jj)
                 do ii= 1, pard1
                    parr1(pard1*(ipoin-1)+ii) = parr1(pard1*(ipoin-1)+ii) + recvbuff_rp(pard1*(jj-1)+ii)
                 enddo
              enddo
           end if
           
        else if(pardi==2.and.parki==2) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN)
           !
           !-------------------------------------------------------------

           if( memory_size(sendbuff_rp) <  pard1*commd%bound_dim ) then
              call memory_deallo(par_memor,'SENDBUFF_RP','par_slexca',sendbuff_rp)
              call memory_alloca(par_memor,'SENDBUFF_RP','par_slexca',sendbuff_rp,pard1*commd%bound_dim,'DO_NOT_INITIALIZE')
           end if
           if( memory_size(recvbuff_rp) <  pard1*commd%bound_dim ) then
              call memory_deallo(par_memor,'RECVBUFF_RP','par_slexca',recvbuff_rp)
              call memory_alloca(par_memor,'RECVBUFF_RP','par_slexca',recvbuff_rp,pard1*commd%bound_dim,'DO_NOT_INITIALIZE')
           end if           

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 sendbuff_rp(pard1*(jj-1)+ii) = parr2(ii,ipoin)
              enddo
           enddo

           do ii= 1, commd%nneig
              dom_i = commd%neights(ii)

              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( sendbuff_rp(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   recvbuff_rp(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
              if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCH')
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_invp(jj)
              do ii= 1, pard1
                 parr2(ii,ipoin) = parr2(ii,ipoin) + recvbuff_rp(pard1*(jj-1)+ii)
              enddo
           enddo

        else if ( pardi==1 .and. parki==4 ) then

           !-------------------------------------------------------------
           !
           ! COMPLEX(NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparx1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARX1','par_slexch',loc_sparx1)

           allocate(loc_rparx1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARX1','par_slexch',loc_rparx1)

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
              if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCH')
#endif
           enddo

           do jj = 1, commd%bound_dim
              ipoin = commd%bound_invp(jj)
              parx1(ipoin) = parx1(ipoin) + loc_rparx1(jj)
           enddo

           call memchk(two,istat,par_memor,'LOC_RPARX1','par_slexch',loc_rparx1)
           deallocate(loc_rparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARX1','par_slexch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARX1','par_slexch',loc_sparx1)
           deallocate(loc_sparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARX1','par_slexch',0_ip)

        else if ( pardi == 1 .and. parki == 7 ) then

           !-------------------------------------------------------------
           !
           ! COMPLEX(PARD1,NPOIN) => COMPLEX(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparx1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARX1','par_slexch',loc_sparx1)

           allocate(loc_rparx1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARX1','par_slexch',loc_rparx1)

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
              if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCH')
#endif
           enddo

           do jj = 1,commd%bound_dim

              ji = pard1 * (jj - 1)
              ipoin = commd%bound_invp(jj)
              poin = pard1 * (ipoin - 1)
              do ii = 1,pard1

                 parx1(poin+ii) = parx1(poin+ii) + loc_rparx1(ji+ii)

              enddo

           enddo

           call memchk(two,istat,par_memor,'LOC_RPARX1','par_slexch',loc_rparx1)
           deallocate(loc_rparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARX1','par_slexch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARX1','par_slexch',loc_sparx1)
           deallocate(loc_sparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARX1','par_slexch',0_ip)


        end if
     end if
  endif

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slexch
