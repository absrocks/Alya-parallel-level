subroutine pararr(wtask,ntype,ndimr,rvarr)
  !------------------------------------------------------------------------
  !****f* Parall/pararr
  ! NAME
  !    pararr
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only  :  ip,rp,lg
  use def_master,         only  :  npari,nparr,nparc,parin,parre,parr1,kfl_paral
  use def_master,         only  :  party,parki,pardi,IPARALL,pard1,pari1,npasi
  use def_master,         only  :  ISEQUEN,INOTMASTER
  use def_master,         only  :  NPOIN_TYPE,NFACE_TYPE,NBOPO_TYPE,lntra,nfacg,lninv_loc
  use def_domain,         only  :  npoin,nbopo,npoin_2
  use mod_couplings,      only  :  COU_INTERPOLATE_NODAL_VALUES
  use mod_communications, only  :  PAR_INTERFACE_NODE_EXCHANGE
  use def_coupli,         only  :  ncoup_implicit
  use def_coupli,         only  :  mcoup
  use def_coupli,         only  :  coupling_type
  use def_coupli,         only  :  ncoup_implicit_d
  use def_coupli,         only  :  lcoup_implicit_d
  use def_coupli,         only  :  ncoup_implicit_n
  use def_coupli,         only  :  lcoup_implicit_n
  use def_coupli,         only  :  RESIDUAL,UNKNOWN
  use def_coupli,         only  :  BETWEEN_SUBDOMAINS
  use def_master,         only  :  current_zone,kfl_paral,modul,ISLAVE,mmodu
  use def_solver,         only  :  solve_sol
  use mod_communications, only  :  PAR_MIN,PAR_MAX,PAR_SUM
  implicit none 
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimr
  real(rp),     target     :: rvarr(ndimr)

  integer(ip)              :: ii,kk,icoup,jcoup
  integer(ip)              :: ipoin,jj,ll,kcoup
  integer(ip)              :: kfl_mask,kpoin
  real(rp),     pointer    :: rvarr_tmp(:)
  real(rp),     pointer    :: rvarr_int(:)


  if( IPARALL .or. ( ISEQUEN .and. (wtask == 'SLX' .or. wtask == 'SLG' .or. wtask == 'SOL') ) ) then

     npari = 0
     nparc = 0

     if( wtask == 'BCT' ) then
        !
        ! Broadcast: master computes something that is broadcasted to the slaves
        !
        nparr =  ndimr
        parre => rvarr
        call Parall(2_ip)
         
     else if( wtask == 'SND' ) then

        nparr =  ndimr
        parre => rvarr
        call Parall(3_ip)

     else if( wtask == 'RCV' ) then

        nparr =  ndimr
        parre => rvarr
        call Parall(4_ip)

     else if( wtask == 'MIN' ) then

        call PAR_MIN(ndimr,rvarr)
        
     else if( wtask == 'SUM' ) then

        call PAR_SUM(ndimr,rvarr)

     else if( wtask == 'MAX' ) then

        call PAR_MAX(ndimr,rvarr)

     else if( wtask == 'S2M' ) then

        nparr =  ndimr
        parre => rvarr
        call Parall(24_ip)

     else if( wtask == 'GAT' ) then

        call runend('PARARR: GATHER NO LONGER HERE')

     else if( wtask == 'IBI' ) then 

        pari1 => lntra
        pard1 =  ndimr/npoin
        parr1 => rvarr(1:ndimr)
        call Parall(405_ip)

     else if( wtask == 'SLX' .or. wtask == 'SLG' .or. wtask == 'SOL' ) then
        !
        ! par_slexch
        !
        if( wtask == 'SOL' ) then
           kfl_mask = 1
        else
           kfl_mask = 0
        end if

        if( INOTMASTER ) then
           if( wtask == 'SLG' ) then
              call Parall(611_ip) 
           end if

           if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
              if( ntype == NPOIN_TYPE ) pard1 = ndimr/npoin
              if( ntype == NBOPO_TYPE ) pard1 = ndimr/nbopo
              party =  ntype
              if( pard1 == 1 ) then
                 parki =  2
                 pardi =  1
              else
                 parki =  5
                 pardi =  1
              end if

              if( ISLAVE ) then
                 parr1 => rvarr
                 call Parall(400_ip) 
              end if
              
           else if( ntype == NFACE_TYPE ) then

              pard1 = ndimr/nfacg
              party = ntype
              if( pard1 == 1 ) then
                 parki =  2
                 pardi =  1
              else
                 parki =  5
                 pardi =  1
              end if
              
              if( ISLAVE ) then
                 parr1 => rvarr
                 call Parall(607_ip)
              end if
              
           else
              call runend('PARARR: NOT CODED')
           end if

           if( wtask == 'SLG' ) then
              call Parall(612_ip)
           end if
           !
           ! Subdomain coupling 
           !
           if( ntype == NPOIN_TYPE .and. ncoup_implicit > 0 ) then
              
              allocate( rvarr_tmp(ndimr),rvarr_int(ndimr) )
              rvarr_tmp(1:ndimr) = rvarr(1:ndimr)
              rvarr_int(1:ndimr) = 0.0_rp
              kk                 = ndimr/npoin
              !
              ! Neumann transmission condition
              !
              do kcoup = 1,ncoup_implicit_n
                 icoup = lcoup_implicit_n(kcoup)
                 if( wtask == 'SOL' ) then
                    call COU_INTERPOLATE_NODAL_VALUES(icoup,kk,rvarr_int,rvarr_tmp,solve_sol(1) % kfl_fixno)
                 else
                    call COU_INTERPOLATE_NODAL_VALUES(icoup,kk,rvarr_int,rvarr_tmp)
                 end if                 
                 !!!call PAR_INTERFACE_NODE_EXCHANGE(kk,rvarr_int,'SUM')
                 rvarr(1:ndimr) = rvarr(1:ndimr) + rvarr_int(1:ndimr)    
              end do
              deallocate(rvarr_int)
              rvarr_tmp(1:ndimr) = rvarr(1:ndimr)
              !
              ! Dirichlet transmission condition
              !              
              do kcoup = 1,ncoup_implicit_d
                 icoup = lcoup_implicit_d(kcoup)
                 if( wtask == 'SOL' ) then
                    call COU_INTERPOLATE_NODAL_VALUES(icoup,kk,rvarr,rvarr_tmp,solve_sol(1) % kfl_fixno)
                 else
                    call COU_INTERPOLATE_NODAL_VALUES(icoup,kk,rvarr,rvarr_tmp)
                 end if
                 !!!call PAR_INTERFACE_NODE_EXCHANGE(kk,rvarr,'DISTRIBUTE') 
              end do
              deallocate(rvarr_tmp)

           end if
        end if

     else if( wtask == 'SLE' ) then
        !
        ! par_slequa: Equal values at interface
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimr/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimr/nbopo
           party =  ntype
           if( pard1 == 1 ) then
              parki =  2
              pardi =  1
           else
              parki =  5
              pardi =  1
           end if
           parr1 => rvarr
           call Parall(900_ip)
        else
           call runend('PARARR: NOT CODED')
        end if

     else if( wtask == 'SLA' ) then
        !
        ! par_slexca
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then

           if( ntype == NPOIN_TYPE ) pard1 = ndimr/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimr/nbopo

           party =  ntype
           if( pard1 == 1 ) then
              parki =  2
              pardi =  1
           else
              parki =  5
              pardi =  1
           end if
           parr1 => rvarr
           call Parall(407_ip)
        else
           call runend('PARARR: NOT CODED')
        end if

     else if( wtask == 'SNR' ) then 
        !
        ! par_slaves
        !
        npasi = 0
        npari = 0
        call Parall(411_ip)

     else if( wtask == 'SRA' ) then 
        !
        ! par_slaves: Send receive to any slave
        !
        npasi = 0
        npari = 0
        call Parall(408_ip)

     else if( wtask == 'SMI' ) then
        !
        ! par_slexch
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimr/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimr/nbopo
           party =  ntype
           if( pard1 == 1 ) then
              parki =  2
              pardi =  1
           else
              parki =  5
              pardi =  1
           end if
           parr1 => rvarr
           call Parall(424_ip)
        else
           call runend('PARARR: NOT CODED')
        end if

     else if( wtask == 'SMA' ) then
        !
        ! par_slexch
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimr/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimr/nbopo
           party =  ntype
           if( pard1 == 1 ) then
              parki =  2
              pardi =  1
           else
              parki =  5
              pardi =  1
           end if
           parr1 => rvarr
           call Parall(425_ip)
        else
           call runend('PARARR: NOT CODED')
        end if

     else if( wtask == 'FRI' ) then
        !
        ! par_slexfr: take fringe geometry nodes
        !
        party = ntype
        if( ntype == NPOIN_TYPE ) then
           pard1 = ndimr/npoin_2
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           else
              parki =  5
              pardi =  2
           end if
           parr1 => rvarr
           call Parall(429_ip)
        end if

     else if( wtask == 'SCH' ) then 

        parr1 => rvarr
        call Parall(426_ip)

     else if( wtask == 'SSS' )  then

        parr1 => rvarr
        call Parall(1005_ip)

     else

        call runend('PARARR: WRONG CASE')

     end if

     nparr = 0
     nullify(parre)
     nullify(parr1)

  end if

end subroutine pararr
