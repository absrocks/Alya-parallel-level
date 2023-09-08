!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_comzon.f90
!> @author  Guillaume Houzeaux
!> @brief   Zone-wise communication array
!> @details Zone-wise communication array
!> @} 
!------------------------------------------------------------------------
subroutine par_comzon(itask)
  use def_parame
  use def_kintyp
  use def_domain
  use def_master
  use def_parall
  use mod_memory 
  use mod_postpr  
  use mod_parall, only    :  PAR_COMM_MY_CODE4,PAR_INTEGER
  use mod_parall, only    :  PAR_COMM_WORLD,commd
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
#endif
  integer(ip),      intent(in)       :: itask
  integer(ip)                        :: ii,jj,kk,ll,ineig,dom_i,jpart,imodu
  integer(ip)                        :: ipart,izone,kzone,myzon,ipoin
  integer(4),          pointer       :: icolor(:)
  integer(4),          pointer       :: ikey(:)
  integer(4)                         :: istat
  integer(ip),                  save :: PAR_COMM_WORLD_SAV
  type(comm_data_par), pointer, save :: commd_sav
  type(comm_data_par), pointer       :: my_commz(:) => null()
  logical(lg)                        :: enter_in_module
  
  nullify(icolor)
  nullify(ikey)

  select case ( itask ) 

  case ( 1_ip )

     !-------------------------------------------------------------
     !
     ! Construct zone-wise communication arrays
     !
     !-------------------------------------------------------------

     if( nzone_par > 1 ) then

        if( ISLAVE ) call memory_alloca(mem_servi(1:2,servi),'LSUBZ_PAR','par_comzon',lsubz_par,npart)

        call par_parari('BCT',0_ip,npart,lsubz_par)
        call par_parari('BCT',0_ip,10_ip,lzone_par)

        if( ISLAVE ) then

           allocate( my_commz(1) )
           nullify( my_commz(1) % neights    )
           nullify( my_commz(1) % bound_size )
           nullify( my_commz(1) % bound_perm )
           nullify( my_commz(1) % bound_scal )
           nullify( my_commz(1) % bface_size )
           nullify( my_commz(1) % bface_perm )
           !nullify( my_commz(1) % frins_size )
           !nullify( my_commz(1) % frins_perm )
           !nullify( my_commz(1) % frinr_size )
           !nullify( my_commz(1) % frinr_perm )

           my_commz(1) % nneig     = 0
           my_commz(1) % bound_dim = 0
           do ii = 1,commd % nneig
              dom_i = commd % neights(ii)        
              if( lsubz_par(dom_i) == lsubz_par(kfl_paral) ) then
                 my_commz(1) % nneig     = my_commz(1) % nneig + 1
                 my_commz(1) % bound_dim = my_commz(1) % bound_dim + commd % bound_size(ii+1) - commd % bound_size(ii)
              end if
           end do

           my_commz(1) % neights    => null()
           my_commz(1) % bound_size => null()
           my_commz(1) % bound_perm => null()

           call memory_alloca(mem_servi(1:2,servi),'MY_COMMZ(1) % NEIGHTS'   ,&
                'par_comzon',my_commz(1) % neights   , max(my_commz(1) % nneig,1_ip))
           call memory_alloca(mem_servi(1:2,servi),'MY_COMMZ(1) % BOUND_SIZE','par_comzon',&
                my_commz(1) % bound_size, my_commz(1) % nneig+1_ip)
           call memory_alloca(mem_servi(1:2,servi),'MY_COMMZ(1) % BOUND_PERM','par_comzon',&
                my_commz(1) % bound_perm, max(my_commz(1) % bound_dim,1_ip))

           kk = 0
           ineig = 0
           do ii = 1,commd % nneig
              dom_i = commd % neights(ii)        
              if( lsubz_par(dom_i) == lsubz_par(kfl_paral) ) then
                 ineig = ineig + 1
                 my_commz(1) % bound_size(ineig) = commd % bound_size(ii+1) - commd % bound_size(ii)
                 my_commz(1) % neights(ineig)    = dom_i
                 do jj =  commd % bound_size(ii),commd % bound_size(ii+1)-1
                    kk = kk + 1
                    my_commz(1) % bound_perm(kk) = commd % bound_perm(jj)
                 end do
              end if
           end do

           jj = my_commz(1) % bound_size(1)
           my_commz(1) % bound_size(1) = 1 
           do ii = 2,my_commz(1) % nneig+1
              ll = my_commz(1) % bound_size(ii)
              my_commz(1) % bound_size(ii) = my_commz(1) % bound_size(ii-1) + jj
              jj = ll
           end do
           commz => my_commz(1)

        else

           allocate( my_commz(1) )
           commz => my_commz(1)

        end if
        !
        ! Split MPI communicator
        !
        call memory_alloca(mem_servi(1:2,servi),'ICOLOR'        ,'par_comzon',icolor,npart+1_ip)
        call memory_alloca(mem_servi(1:2,servi),'IKEY'          ,'par_comzon',ikey,npart+1_ip)
        call memory_alloca(mem_servi(1:2,servi),'PAR_COMM_ZONES','par_comzon',PAR_COMM_ZONES,nzone_par)

        do kzone = 1,nzone_par
           icolor    = 0
           ikey      = 0
           jpart     = 0
           do ipart = 2,npart+1
              if( lsubz_par(ipart-1) == kzone ) then
                 jpart = jpart + 1
                 ikey(ipart) = int(jpart,4)
              end if
           end do
           if( IMASTER ) then
              icolor(1) = 1
              ikey(1)   = 0
           else if( lsubz_par(kfl_paral) == kzone ) then
              icolor(1) = 1        
              ikey(1)   = ikey(kfl_paral)
           else
              icolor(1) = 0
              ikey(1)   = 0
           end if
#ifdef MPI_OFF
#else
           CALL MPI_COMM_SPLIT(PAR_COMM_MY_CODE4, icolor, ikey, PAR_COMM_ZONES(kzone), istat )
#endif
        end do

        call memory_deallo(mem_servi(1:2,servi),'ICOLOR','par_comzon',icolor)
        call memory_deallo(mem_servi(1:2,servi),'IKEY',  'par_comzon',ikey)
        !
        ! My communicator
        !
        if( INOTMASTER ) myzon = lsubz_par(kfl_paral) 
        if( IMASTER ) then
           commz % PAR_COMM_WORLD = PAR_COMM_ZONES(1)
        else
           commz % PAR_COMM_WORLD = PAR_COMM_ZONES(myzon)
        end if
        !
        ! For scalar products, some own boundary nodes should be added if no one
        ! from my zone owns it
        !
        if( INOTMASTER ) then
           !
           ! Check if inside my zone we have lost own boundary nodes, that is
           ! nodes that are not own boundary nodes in any of the subdomains
           ! of my zone
           !
           call memgen(1_ip,npoin,0_ip)
           kzone = myzon
           do ipoin = npoi2,npoi3
              gisca(ipoin) = -1
           end do
           !
           ! Exchange in my zone only
           !
           commd_sav => commd
           commd     => commz 
           call par_parari('SLX',NPOIN_TYPE,npoin,gisca)
           !
           ! Decide who's got the lost boundary point: subdomain woth maximum rank
           !
           do ii = 1,commd_glo % nneig
              dom_i = commd_glo % neights(ii)
              do jj = commd_glo % bound_size(ii),commd_glo % bound_size(ii+1)-1
                 ipoin = commd_glo % bound_perm(jj)
                 if( gisca(ipoin) >= 0 ) then
                    if( lsubz_par(dom_i) == myzon ) then
                       gisca(ipoin) = max(kfl_paral,dom_i) 
                    else
                       gisca(ipoin) = max(kfl_paral,gisca(ipoin))
                    end if
                 end if
              end do
           end do
            
           commd => commd_sav

           call memgen(3_ip,npoin,0_ip)
        end if
        !
        ! From now on, LSUBZ_PAR(1) is my zone number
        !
        call memory_deallo(mem_servi(1:2,servi),'LSUBZ_PAR','par_comzon',lsubz_par)
        call memory_alloca(mem_servi(1:2,servi),'LSUBZ_PAR','par_comzon',lsubz_par,1_ip)
        if( IMASTER ) then
           lsubz_par(1) = -1
        else
           lsubz_par(1) =  myzon
        end if

     end if
     !
     ! Compute module communication arrays
     !m
     if( nzone_par <= 1 ) then
        do imodu = 1,mmodu
           if( kfl_modul(imodu) == 1 ) then
              momod(imodu) % commd => commd
           end if
        end do
     else
        do imodu = 1,mmodu
           if( kfl_modul(imodu) == 1 ) then
              kzone_loop: do kzone = 1,nzone_par
                 izone = lzone_par(kzone)
                 if( izone == lzone(imodu) ) then
                    momod(imodu) % commd => commz 
                    exit kzone_loop
                 else
                    momod(imodu) % commd => commd
                 end if
              end do kzone_loop
           end if
        end do
     end if

  case ( 2_ip )

     !-------------------------------------------------------------
     !
     ! Construct zone-wise communication arrays
     !
     !-------------------------------------------------------------

     if( igene == ITASK_DOITER .and. nzone_par > 1 .and. modul > 0 ) then
        enter_in_module = .false.
        izone = lzone(modul)                          ! Geometrical zone of current module
        myzon = lsubz_par(1)                          ! My process zone
        if( myzon == -1 ) then
           enter_in_module = .true.
        else if( lzone_par(myzon) == izone ) then
           enter_in_module = .true.
        else
           enter_in_module = .false.
        end if
        if( enter_in_module ) then
           commd => momod(modul) % commd
           if( IMASTER ) then
              kzone = lzone_par(1)                    ! Parallel zone of current process
              ii = 1
              do while( kzone /= izone )
                 ii = ii + 1
                 kzone = lzone_par(ii)
              end do
              PAR_COMM_WORLD = PAR_COMM_ZONES(kzone)
           else if( enter_in_module ) then
              PAR_COMM_WORLD = momod(modul) % commd % PAR_COMM_WORLD
           end if
        end if
        if( .not. enter_in_module ) modul = -2_ip
     end if

  case ( 3_ip ) 

     !-------------------------------------------------------------
     !
     ! Recover global communication arrays and PAR_COMM_MY_CODE4
     !
     !-------------------------------------------------------------

     if( igene == ITASK_DOITER .and. nzone_par > 1 ) then
        commd          => commd_glo
#ifdef MPI_OFF
#else
        PAR_COMM_WORLD =  PAR_COMM_MY_CODE4
#endif
     end if

  case ( 4_ip ) 

     !-------------------------------------------------------------
     !
     ! Force global communication arrays and PAR_COMM_MY_CODE4
     !
     !-------------------------------------------------------------

     commd_sav          => commd
     PAR_COMM_WORLD_SAV =  PAR_COMM_WORLD

     commd              => commd_glo
#ifdef MPI_OFF
#else
     PAR_COMM_WORLD     =  PAR_COMM_MY_CODE4
#endif

  case ( 5_ip ) 

     !-------------------------------------------------------------
     !
     ! Recover last communication arrays and MPI communicator
     !
     !-------------------------------------------------------------

     commd              => commd_sav
     PAR_COMM_WORLD     =  PAR_COMM_WORLD_SAV

  case ( 6_ip )

     !-------------------------------------------------------------
     !
     ! Check if my zone solves module modul
     !
     !-------------------------------------------------------------

     if( nzone_par > 1 ) then
        enter_in_module = .false.
        igene           = 1
        izone           = lzone(modul)                          ! Geometrical zone of current module
        myzon           = lsubz_par(1)                          ! My process zone
        if( myzon == -1 ) then
           enter_in_module = .true.
        else if( lzone_par(myzon) == izone ) then
           enter_in_module = .true.
        else
           enter_in_module = .false.
        end if
        if( .not. enter_in_module ) igene = 0
     end if

  end select


end subroutine par_comzon
