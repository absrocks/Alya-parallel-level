subroutine par_sencom
  !------------------------------------------------------------------------
  !****f* Parall/par_sencom
  ! NAME
  !   par_sencom 
  ! DESCRIPTION
  !    This rouotine sends communication strategy to slaves
  ! OUTPUT
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_parall
  use def_master
  use def_domain
  use mod_memchk
  use mod_par_memchk
  use mod_memory
  use mod_parall, only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall, only : commd
  use mod_parall, only : PAR_COMM_MY_CODE
  use mod_parall, only : PAR_INITIALIZE_COMMUNICATION_ARRAY
  implicit none 
  type(comm_data_par), pointer :: cdata(:)
  integer(ip)                  :: domai,dom_j,dom_i,kk
  integer(ip)                  :: ipoin,ii,jj
  
  allocate( PAR_COMM_MY_CODE_ARRAY(1) )
  call PAR_INITIALIZE_COMMUNICATION_ARRAY(PAR_COMM_MY_CODE_ARRAY(1))
  PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD = PAR_COMM_MY_CODE
  PAR_COMM_MY_CODE_ARRAY(1) % RANK4          = int(kfl_paral,4)
 
  if( IMASTER ) then

     !-------------------------------------------------------------------
     !
     ! Master
     !
     !-------------------------------------------------------------------

     if( .not. READ_AND_RUN() ) then

        nullify(cdata)
        allocate(cdata(npart_par))
        do domai = 1,npart_par
           call PAR_INITIALIZE_COMMUNICATION_ARRAY(cdata(domai))
        end do

        do domai = 1,npart_par

           kfl_desti_par = domai
           nneig         = lneig_par(domai) 

           call memory_alloca(mem_servi(1:2,servi),'NEIGHTS'   ,'par_sencom',cdata(domai) % neights,nneig)
           call memory_alloca(mem_servi(1:2,servi),'BOUND_SIZE','par_sencom',cdata(domai) % bound_size,nneig+1)
           !
           ! List of neighbours
           !
           dom_j = 1
           do dom_i = 1, nbcol
              if (lcomm_par(dom_i,domai)/=0) then
                 cdata(domai) % neights(dom_j) = lcomm_par(dom_i,domai)
                 dom_j = dom_j + 1
              endif
           enddo
           !
           ! Boundary size with every neighbour
           !
           cdata(domai) % bound_size(1) = 1
           do dom_i = 1, nneig
              dom_j = cdata(domai) % neights(dom_i)
              if( domai > dom_j ) then
                 kk = (domai*(domai-1))/2 + dom_j
              else
                 kk = (dom_j*(dom_j-1))/2 + domai
              endif
              cdata(domai) % bound_size(dom_i+1) = cdata(domai) % bound_size(dom_i) + neighdom(kk)
           enddo
           cdata(domai) % bound_dim = cdata(domai) % bound_size(nneig+1) - 1

           call memory_alloca(mem_servi(1:2,servi),'BOUND_PERM','par_sencom',cdata(domai) % bound_perm,cdata(domai) % bound_dim)

        end do
        !
        ! boundary points
        !
        do ipoin= 1, gnb
           do ii= badj(ipoin), badj(ipoin+1)-1
              dom_i = bdom(ii)
              do jj= ii+1, badj(ipoin+1)-1
                 dom_j = bdom(jj)
                 !
                 ! Update dom_i
                 !
                 kk = 1
                 do while (cdata(dom_i)%neights(kk)/=dom_j)
                    kk = kk + 1
                 enddo
                 cdata(dom_i) % bound_perm(cdata(dom_i) % bound_size(kk)) = bpoin(ii)
                 cdata(dom_i) % bound_size(kk) = cdata(dom_i) % bound_size(kk) + 1
                 !
                 ! Update dom_j
                 !
                 kk = 1
                 do while( cdata(dom_j) % neights(kk) /= dom_i )
                    kk = kk + 1
                 end do
                 cdata(dom_j) % bound_perm(cdata(dom_j) % bound_size(kk)) = bpoin(jj)
                 cdata(dom_j) % bound_size(kk) = cdata(dom_j)%bound_size(kk) + 1
              end do
           end do
        end do

        do domai = 1,npart_par

           kfl_desti_par = domai
           nneig         = lneig_par(domai)
           !
           ! Recompute bound_size
           !
           do dom_i = nneig, 1, -1
              cdata(domai) % bound_size(dom_i+1) = cdata(domai) % bound_size(dom_i)
           enddo
           cdata(domai) % bound_size(1) = 1

           npari =  nneig
           parin => cdata(domai) % neights
           strin =  'CDATA(DOMAIN) % NEIGHTS'
           call par_sendin()

           npari =  nneig+1
           parin => cdata(domai) % bound_size
           strin =  'CDATA(DOMAIN) % BOUND_SIZE'
           call par_sendin()

           npari =  cdata(domai) % bound_dim 
           parin => cdata(domai) % bound_perm
           strin =  'CDATA(DOMAIN) % BOUND_PERM'
           call par_sendin()

           call memory_deallo(mem_servi(1:2,servi),'BOUND_PERM','par_sencom',cdata(domai) % bound_perm)
           call memory_deallo(mem_servi(1:2,servi),'BOUND_PERM','par_sencom',cdata(domai) % bound_size)
           call memory_deallo(mem_servi(1:2,servi),'BOUND_PERM','par_sencom',cdata(domai) % neights)

        end do

        deallocate(cdata)
        !
        ! Deallocate memory
        !
        call par_memory(five)
     end if
     !
     ! Create communication array for master
     !
     !allocate(PAR_COMM_MY_CODE_ARRAY(1))
     !PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD = PAR_COMM_MY_CODE
     commd     => PAR_COMM_MY_CODE_ARRAY(1)
     commd_glo => PAR_COMM_MY_CODE_ARRAY(1)

  else if ( ISLAVE ) then

     !-------------------------------------------------------------------
     !
     ! Slaves
     !
     !-------------------------------------------------------------------

     kfl_desti_par = 0

     call memory_alloca(mem_servi(1:2,servi),'NEIGHTS'   ,'par_sencom',PAR_COMM_MY_CODE_ARRAY(1) % neights,nneig)
     call memory_alloca(mem_servi(1:2,servi),'BOUND_SIZE','par_sencom',PAR_COMM_MY_CODE_ARRAY(1) % bound_size,nneig+1)

     npari =  nneig
     parin => PAR_COMM_MY_CODE_ARRAY(1) % neights
     call par_receiv()

     npari =  nneig+1
     parin => PAR_COMM_MY_CODE_ARRAY(1) % bound_size
     call par_receiv()

     PAR_COMM_MY_CODE_ARRAY(1) % bound_dim = PAR_COMM_MY_CODE_ARRAY(1) % bound_size(nneig+1) - 1

     call memory_alloca(mem_servi(1:2,servi),'BOUND_PERM','par_sencom',PAR_COMM_MY_CODE_ARRAY(1) % bound_perm,PAR_COMM_MY_CODE_ARRAY(1) % bound_dim)

     npari =  PAR_COMM_MY_CODE_ARRAY(1) % bound_dim
     parin => PAR_COMM_MY_CODE_ARRAY(1) % bound_perm
     call par_receiv()

     PAR_COMM_MY_CODE_ARRAY(1) % nneig = nneig 
     PAR_COMM_MY_CODE_ARRAY(1) % npoi1 = npoi1
     PAR_COMM_MY_CODE_ARRAY(1) % npoi2 = npoi2
     PAR_COMM_MY_CODE_ARRAY(1) % npoi3 = npoi3

     commd     => PAR_COMM_MY_CODE_ARRAY(1)
     commd_glo => PAR_COMM_MY_CODE_ARRAY(1)

  end if

  nullify(parin)
  nullify(parre)

end subroutine par_sencom
