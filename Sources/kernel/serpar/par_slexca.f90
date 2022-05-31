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
  use mod_parall,               only : commd,PAR_COMM_MY_CODE4
  use mod_parall,               only : PAR_INTEGER
  use mod_parall,               only : par_memor
  use mod_parall,               only : sendbuff_rp
  use mod_parall,               only : recvbuff_rp
  use mod_parall,               only : sendbuff_ip
  use mod_parall,               only : recvbuff_ip
  use mod_communications_tools, only : PAR_MPI_RUNEND
  use mod_memory
  use mod_memchk

  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
  integer(4)                 :: status(MPI_STATUS_SIZE,2*commd % nneig)
#endif
  integer(ip)                :: ipoin,ii,jj,kk,bsize,ini,dom_i
  integer(ip)                :: ndofi,ndofj
  integer(ip), save          :: count4
  integer(4)                 :: istat,bsize4
  real(rp)                   :: time1,time2
  integer(4)                 :: ierr4

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
              if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_invp(jj)
              pari1(ipoin) = pari1(ipoin) + recvbuff_ip(jj)
           enddo

        else if(pardi==1.and.parki==2) then

           !-------------------------------------------------------------
           !
           ! REAL(NPOIN)
           !
           !-------------------------------------------------------------

           if( ipass_par == 0 ) then

              ipass_par = 1
              allocate(ireq4(commd%nneig*2),stat=istat)
              call memchk(zero,istat,par_memor,'IREQ4','par_slexch',ireq4)

              if( memory_size(sendbuff_rp) <  commd%bound_dim ) then
                 call memory_deallo(par_memor,'SENDBUFF_RP','par_slexca',sendbuff_rp)
                 call memory_alloca(par_memor,'SENDBUFF_RP','par_slexca',sendbuff_rp,commd%bound_dim,'DO_NOT_INITIALIZE')
              end if
              if( memory_size(recvbuff_rp) <  commd%bound_dim ) then
                 call memory_deallo(par_memor,'RECVBUFF_RP','par_slexca',recvbuff_rp)
                 call memory_alloca(par_memor,'RECVBUFF_RP','par_slexca',recvbuff_rp,commd%bound_dim,'DO_NOT_INITIALIZE')
              end if

              do jj= 1, commd%bound_dim
                 ipoin = commd%bound_perm(jj)
                 sendbuff_rp(jj) = parr1(ipoin)
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
                 call MPI_Isend( sendbuff_rp(ini:ini+bsize-1), bsize4, &
                      MPI_DOUBLE_PRECISION,  dom_i, 0_4,   &
                      PAR_COMM_MY_CODE4, ireq4(kk), istat )
                 if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
                 kk = kk + 1
                 call MPI_Irecv( recvbuff_rp(ini:ini+bsize-1), bsize4, &
                      MPI_DOUBLE_PRECISION,  dom_i, 0_4,   &
                      PAR_COMM_MY_CODE4, ireq4(kk), istat )
                 if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
#endif
              enddo

           else


              ipass_par  = 0
              count4 = commd % nneig*2
#ifdef MPI_OFF
#else
              CALL MPI_WAITALL(count4,ireq4,status,ierr4)
              if( ierr4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(ierr4,'PAR_SLEXCA')
#endif

              do jj = 1,commd % bound_dim
                 ipoin = commd % bound_invp(jj)
                 parr1(ipoin) = parr1(ipoin) + recvbuff_rp(jj)
              enddo

              call memchk(two,istat,par_memor,'IREQ4','par_slexch',ireq4)
              deallocate(ireq4,stat=istat)
              if(istat/=0) call memerr(two,'IREQ4','par_slexch',0_ip)

           end if


        else if(pardi>=1.and.parki==5) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN) => REAL(PARD1*NPOIN) 
           !
           !-------------------------------------------------------------

           if( ipass_par == 0 ) then

              ipass_par = 1

              allocate(ireq4(commd%nneig*2),stat=istat)
              call memchk(zero,istat,par_memor,'IREQ4','par_slexch',ireq4)

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
                 ndofi = pard1*(ipoin-1)
                 ndofj = pard1*(jj-1)
                 do ii= 1, pard1
                    ndofi = ndofi + 1
                    ndofj = ndofj + 1
                    sendbuff_rp(ndofj) = parr1(ndofi)
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
                 call MPI_Isend( sendbuff_rp(ini:ini+bsize-1), bsize4, &
                      MPI_DOUBLE_PRECISION,  dom_i, 0_4,   &
                      PAR_COMM_MY_CODE4, ireq4(kk), istat )
                 if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
                 kk = kk + 1
                 call MPI_Irecv( recvbuff_rp(ini:ini+bsize-1), bsize4, &
                      MPI_DOUBLE_PRECISION,  dom_i, 0_4,   &
                      PAR_COMM_MY_CODE4, ireq4(kk), istat )              
                 if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
#endif
              end do

           else

              ipass_par  = 0
              count4 = commd%nneig*2

#ifdef MPI_OFF
#else
              CALL MPI_WAITALL(count4,ireq4,status,ierr4)
              if( ierr4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(ierr4,'PAR_SLEXCA')
#endif
              do jj = 1,commd%bound_dim
                 ipoin = commd%bound_invp(jj)
                 ndofi = pard1*(ipoin-1)
                 ndofj = pard1*(jj-1)
                 do ii = 1,pard1
                    ndofi = ndofi + 1
                    ndofj = ndofj + 1
                    parr1(ndofi) = parr1(ndofi) + recvbuff_rp(ndofj)
                 end do
              end do

              call memchk(two,istat,par_memor,'IREQ4','par_slexch',ireq4)
              deallocate(ireq4,stat=istat)
              if(istat/=0) call memerr(two,'IREQ4','par_slexch',0_ip)

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
              call MPI_Sendrecv( &
                   sendbuff_rp(ini:), bsize4,&
                   MPI_DOUBLE_PRECISION,  dom_i, 0_4,     &
                   recvbuff_rp(ini:), bsize4,              &
                   MPI_DOUBLE_PRECISION, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE4, status, istat )
              if( istat /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat,'PAR_SLEXCA')
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_invp(jj)
              do ii= 1, pard1
                 parr2(ii,ipoin) = parr2(ii,ipoin) + recvbuff_rp(pard1*(jj-1)+ii)
              enddo
           enddo

        else if(pardi==1.and.parki==6) then

           call runend('OBSOLETE')

        end if

     end if
  end if

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slexca
