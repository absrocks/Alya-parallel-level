subroutine par_sengro
  !------------------------------------------------------------------------
  !****f* Parall/par_sengro
  ! NAME
  !   par_sengro 
  ! DESCRIPTION
  !    This rouotine sends communication strategy to slaves
  ! OUTPUT
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use      def_kintyp
  use      def_parame
  use      def_parall
  use      def_master
  use      def_domain
  use      mod_memchk
  use      mod_par_memchk
  use mod_parall, only : par_memor
  implicit none

  type(comm_data_par), pointer :: cdata(:)
  integer(ip)                  :: domai,dom_j,dom_i,kk
  integer(ip)                  :: ipoin,ii,jj
  integer(ip),         target  :: dummi(1)
  integer(4)                   :: istat

  if( IMASTER ) then

     !-------------------------------------------------------------------
     !
     ! Master
     !
     !-------------------------------------------------------------------

     allocate(cdata(npart_par),stat=istat)
     call par_memchk(zero,istat,par_memor,'CDATA','par_sengro',cdata)

     do domai= 1, npart_par
        kfl_desti_par = domai
        nneig         = comle(icoml)%lneig_par(domai) 

        allocate( cdata(domai)%neights(nneig),stat=istat)
        call memchk( zero, istat, par_memor, 'cdata(domai)%neights', &
             'par_sengro', cdata(domai)%neights )

        allocate( cdata(domai)%bound_size(nneig+1),stat=istat)
        call memchk( zero, istat, par_memor, 'cdata(domai)%bound_size', &
             'par_sengro', cdata(domai)%bound_size )
        !
        ! NEIGHTS:
        ! Compacted list of neighbours following communication ordering
        !
        dom_j = 1
        do dom_i = 1, comle(icoml)%nbcol
           if ( comle(icoml)%lcomm_par(dom_i,domai) /= 0 ) then
              cdata(domai)%neights(dom_j) = comle(icoml)%lcomm_par(dom_i,domai)
              dom_j = dom_j + 1
           endif
        enddo
        !
        ! cdata(domai)%bound_dim
        ! Boundary size with every neighbour communication of DOMAI with DOM_J
        !
        cdata(domai)%bound_size(1) = 1  
        do dom_i = 1, nneig
           dom_j = cdata(domai)%neights(dom_i)
           if (domai>dom_j) then
              kk = (domai*(domai-1))/2 + dom_j
           else
              kk = (dom_j*(dom_j-1))/2 + domai
           endif
           cdata(domai)%bound_size(dom_i+1) = cdata(domai)%bound_size(dom_i) + &
                comle(icoml)%neighdom(kk)
        enddo

        cdata(domai)%bound_dim = cdata(domai)%bound_size(nneig+1) - 1
        allocate( cdata(domai)%bound_perm(cdata(domai)%bound_dim),stat=istat)
        call memchk( zero, istat, par_memor, 'cdata(domai)%bound_perm', &
             'par_sengro', cdata(domai)%bound_perm )
     enddo
     !
     ! boundary points: boundary point IPOIN belongs to DOM_I and DOM_J
     !
     do ipoin = 1, comle(icoml)%gnb
        do ii = comle(icoml)%badj(ipoin), comle(icoml)%badj(ipoin+1)-1
           dom_i = comle(icoml)%bdom(ii)
           do jj = ii+1, comle(icoml)%badj(ipoin+1)-1
              dom_j = comle(icoml)%bdom(jj)
              !
              ! Update dom_i
              !
              kk = 1
              do while ( cdata(dom_i)%neights(kk) /= dom_j )
                 kk = kk + 1
              enddo
              cdata(dom_i)%bound_perm( cdata(dom_i)%bound_size(kk) ) = comle(icoml)%bpoin(ii) 
              cdata(dom_i)%bound_size( kk )                          = cdata(dom_i)%bound_size(kk) + 1
              !
              ! Update dom_j
              !
              kk = 1
              do while ( cdata(dom_j)%neights(kk) /= dom_i )
                 kk = kk + 1
              enddo
              cdata(dom_j)%bound_perm( cdata(dom_j)%bound_size(kk) ) = comle(icoml)%bpoin(jj)
              cdata(dom_j)%bound_size( kk )                          = cdata(dom_j)%bound_size(kk) + 1
           enddo
        enddo
     enddo

     do domai= 1, npart_par

        kfl_desti_par = domai
        nneig         = comle(icoml)%lneig_par(domai)
        !
        ! recompute bound_size
        !
        do dom_i = nneig, 1, -1
           cdata(domai)%bound_size(dom_i+1) = cdata(domai)%bound_size(dom_i)
        enddo
        cdata(domai)%bound_size(1) = 1

        npari    = 1
        dummi(1) = nneig
        parin => dummi
        strin =  'CDATA(DOMAIN)%NNEIG'
        call par_sendin()

        npari =  nneig
        parin => cdata(domai)%neights
        strin =  'CDATA(DOMAIN)%NEIGHTS'
        call par_sendin()

        npari =  nneig+1
        parin => cdata(domai)%bound_size
        strin =  'CDATA(DOMAIN)%BOUND_SIZE'
        call par_sendin()

        npari =  cdata(domai)%bound_dim 
        parin => cdata(domai)%bound_perm
        strin =  'CDATA(DOMAIN)%BOUND_PERM'
        call par_sendin()

        call memchk( two, istat, par_memor, 'cdata(domai)%bound_perm', &
             'par_sengro', cdata(domai)%bound_perm)
        deallocate( cdata(domai)%bound_perm, stat=istat )
        if(istat/=0) call memerr( two, 'cdata(domai)%bound_perm', 'par_sengro',0_ip)

        call memchk( two, istat, par_memor, 'cdata(domai)%bound_size', &
             'par_sengro', cdata(domai)%bound_size)
        deallocate( cdata(domai)%bound_size, stat=istat )
        if(istat/=0) call memerr( two, 'cdata(domai)%bound_size', 'par_sengro',0_ip)

        call memchk( two, istat, par_memor, 'cdata(domai)%neights', &
             'par_sengro', cdata(domai)%neights)
        deallocate( cdata(domai)%neights, stat=istat )
        if(istat/=0) call memerr( two, 'cdata(domai)%neights', 'par_sengro',0_ip)
     enddo

     call par_memchk(two,istat,par_memor,'cdata','par_sengro',cdata)
     deallocate( cdata, stat=istat )
     if(istat/=0) call memerr(two,'cdata','par_sengro',0_ip)

  else if ( ISLAVE ) then

     !-------------------------------------------------------------------
     !
     ! Slaves
     !
     !-------------------------------------------------------------------

     kfl_desti_par = 0

     npari = 1
     parin => dummi
     call par_receiv()
     comle(icoml)%nneig = dummi(1)

     allocate( cdata(1),stat=istat)

     allocate( cdata(1)%neights(comle(icoml)%nneig),stat=istat)
     call memchk( zero, istat, par_memor, 'cdata(1)%neights', &
          'par_sengro', cdata(1)%neights )

     allocate( cdata(1)%bound_size(comle(icoml)%nneig+1),stat=istat)
     call memchk( zero, istat, par_memor, 'cdata(1)%bound_size', &
          'par_sengro', cdata(1)%bound_size )

     npari =  comle(icoml)%nneig
     parin => cdata(1)%neights
     call par_receiv()
    
     npari =  comle(icoml)%nneig+1
     parin => cdata(1)%bound_size
     call par_receiv()
 
     cdata(1)%bound_dim = cdata(1)%bound_size(comle(icoml)%nneig+1) - 1
     allocate( cdata(1)%bound_perm(cdata(1)%bound_dim),stat=istat)
     call memchk( zero, istat, par_memor, 'cdata(1)%bound_perm', &
          'par_sengro', cdata(1)%bound_perm )

     npari =  cdata(1)%bound_dim
     parin => cdata(1)%bound_perm
     call par_receiv()

     comle(icoml)%commd => cdata(1)

  end if

end subroutine par_sengro
