subroutine ker_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_inivar
  ! NAME
  !    ker_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables
  !                are not initialized before
  !    ITASK=3 ... When starting a time step (from ker_begste)
  ! USES
  ! USED BY
  !    ker_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_cutele
  use mod_ker_space_time_function
  use mod_ker_detection,          only : ker_events_directory_name
  use mod_ker_subdomain,          only : ker_subdomain_initialization
  use mod_ker_proper,             only : ker_proper_initialization
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ifunc

  select case ( itask )

  case ( 0_ip )
     !
     ! Postprocess
     !
     postp(1) % kfl_oonce     = 1
     postp(1) % npp_iniso     = 1

     postp(1) % wopos( 1, 1)  = 'EXNOR'
     postp(1) % wopos( 1, 2)  = '-----' ! Free
     postp(1) % wopos( 1, 3)  = 'LPOIN'
     postp(1) % wopos( 1, 4)  = 'SKCO1'
     postp(1) % wopos( 1, 5)  = 'SKCO2'
     postp(1) % wopos( 1, 6)  = 'SKCO3'
     postp(1) % wopos( 1, 7)  = 'HANGI'
     postp(1) % wopos( 1, 8)  = 'DISPM'
     postp(1) % wopos( 1, 9)  = 'LPOTY'
     postp(1) % wopos( 1,10)  = 'PONUM'
     postp(1) % wopos( 1,11)  = 'CODNO'
     postp(1) % wopos( 1,12)  = 'YWALP'
     postp(1) % wopos( 1,13)  = 'LNTIB'
     postp(1) % wopos( 1,14)  = 'HOLES'

     postp(1) % wopos( 1,15)  = 'DENSI'
     postp(1) % wopos( 1,16)  = 'VISCO'
     postp(1) % wopos( 1,17)  = 'POROS'
     postp(1) % wopos( 1,18)  = 'CONDU'
     postp(1) % wopos( 1,19)  = 'SPECI'
     postp(1) % wopos( 1,20)  = 'DUMMY'
     postp(1) % wopos( 1,22)  = 'LEVEL'
     postp(1) % wopos( 1,23)  = 'GROUP'
     postp(1) % wopos( 1,24)  = 'MASSM'

     postp(1) % wopos( 1,25)  = 'MASSC'
     postp(1) % wopos( 1,26)  = 'GEONO'
     postp(1) % wopos( 1,27)  = 'SUBDO'
     postp(1) % wopos( 1,28)  = 'WALLD'
     postp(1) % wopos( 1,29)  = 'ROUGH'
     postp(1) % wopos( 1,30)  = 'KEKET'
     postp(1) % wopos( 1,31)  = 'CODBO'
     postp(1) % wopos( 1,32)  = 'MATER'
     postp(1) % wopos( 1,33)  = 'LMATN'
     postp(1) % wopos( 1,34)  = 'LNOCH'! TO REMOVE
     postp(1) % wopos( 1,35)  = 'LELEV'
     postp(1) % wopos( 1,36)  = 'QUALI'
     postp(1) % wopos( 1,37)  = 'ELNUM'
     postp(1) % wopos( 1,38)  = 'CODBB'
     postp(1) % wopos( 1,39)  = 'LETIB'
     postp(1) % wopos( 1,40)  = 'LTYP2'
     postp(1) % wopos( 1,41)  = 'LELC2'
     postp(1) % wopos( 1,42)  = 'CONTA'
     postp(1) % wopos( 1,43)  = 'RDOM '
     postp(1) % wopos( 1,44)  = 'VORTX'
     postp(1) % wopos( 1,45)  = 'LESET'
     postp(1) % wopos( 1,46)  = 'DISMM'
     postp(1) % wopos( 1,47)  = 'VELOC'
     postp(1) % wopos( 1,48)  = 'LBSET'
     postp(1) % wopos( 1,49)  = 'LMESH'
     postp(1) % wopos( 1,50)  = 'TURBU'
     postp(1) % wopos( 1,51)  = 'LNSUB'
     postp(1) % wopos( 1,52)  = 'WETNO'
     postp(1) % wopos( 1,53)  = 'LESUB'
     postp(1) % wopos( 1,54)  = 'CANOP'
     postp(1) % wopos( 1,55)  = 'HEIGH'
     postp(1) % wopos( 1,56)  = 'BEATR'
     postp(1) % wopos( 1,57)  = 'COLOR'
     postp(1) % wopos( 1,58)  = 'WALLN'
     postp(1) % wopos( 1,59)  = 'LOCAL'
     postp(1) % wopos( 1,60)  = 'ADVEC'
     postp(1) % wopos( 1,62)  = 'RENEL'
     postp(1) % wopos( 1,63)  = 'RENPO'
     postp(1) % wopos( 1,64)  = 'ORIOL'
     postp(1) % wopos( 1,65)  = 'OMPSS'
     postp(1) % wopos( 1,66)  = 'OPENM'
     postp(1) % wopos( 1,67)  = 'GAUSS'
     postp(1) % wopos( 1,68)  = 'WALLO'
     postp(1) % wopos( 1,69)  = 'MIXIN'
     postp(1) % wopos( 1,70)  = 'FIELD'
     postp(1) % wopos( 1,71)  = 'VELOM'
     postp(1) % wopos( 1,72)  = 'OMPSB'
     postp(1) % wopos( 1,73)  = 'NSWVI'
     postp(1) % wopos( 1,74)  = 'NUMBE'

     postp(1) % wopos( 2, 1)  = 'VECTO'
     postp(1) % wopos( 2, 2)  = 'SCALA'
     postp(1) % wopos( 2, 3)  = 'SCALA'
     postp(1) % wopos( 2, 4)  = 'VECTO'
     postp(1) % wopos( 2, 5)  = 'VECTO'
     postp(1) % wopos( 2, 6)  = 'VECTO'
     postp(1) % wopos( 2, 7)  = 'SCALA'
     postp(1) % wopos( 2, 8)  = 'VECTO'
     postp(1) % wopos( 2, 9)  = 'SCALA'
     postp(1) % wopos( 2,10)  = 'SCALA'
     postp(1) % wopos( 2,11)  = 'VECTO'
     postp(1) % wopos( 2,12)  = 'SCALA'
     postp(1) % wopos( 2,13)  = 'SCALA'
     postp(1) % wopos( 2,14)  = 'SCALA'

     postp(1) % wopos( 2,15)  = 'SCALA'
     postp(1) % wopos( 2,16)  = 'SCALA'
     postp(1) % wopos( 2,17)  = 'SCALA'
     postp(1) % wopos( 2,18)  = 'SCALA'
     postp(1) % wopos( 2,19)  = 'SCALA'
     postp(1) % wopos( 2,20)  = 'SCALA'

     postp(1) % wopos( 2,22)  = 'SCALA'
     postp(1) % wopos( 2,23)  = 'SCALA'
     postp(1) % wopos( 2,24)  = 'SCALA'
     postp(1) % wopos( 2,25)  = 'SCALA'
     postp(1) % wopos( 2,26)  = 'SCALA'
     postp(1) % wopos( 2,27)  = 'SCALA'
     postp(1) % wopos( 2,28)  = 'SCALA'
     postp(1) % wopos( 2,29)  = 'SCALA'
     postp(1) % wopos( 2,30)  = 'SCALA'
     postp(1) % wopos( 2,31)  = 'SCALA'
     postp(1) % wopos( 2,32)  = 'SCALA'
     postp(1) % wopos( 2,33)  = 'SCALA'
     postp(1) % wopos( 2,34)  = 'SCALA'
     postp(1) % wopos( 2,35)  = 'SCALA'
     postp(1) % wopos( 2,36)  = 'SCALA'
     postp(1) % wopos( 2,37)  = 'SCALA'
     postp(1) % wopos( 2,38)  = 'SCALA'
     postp(1) % wopos( 2,39)  = 'SCALA'
     postp(1) % wopos( 2,40)  = 'SCALA'
     postp(1) % wopos( 2,41)  = 'SCALA'
     postp(1) % wopos( 2,42)  = 'SCALA'
     postp(1) % wopos( 2,43)  = 'SCALA'
     postp(1) % wopos( 2,44)  = 'VECTO'
     postp(1) % wopos( 2,45)  = 'SCALA'
     postp(1) % wopos( 2,46)  = 'VECTO'
     postp(1) % wopos( 2,47)  = 'VECTO'
     postp(1) % wopos( 2,48)  = 'SCALA'
     postp(1) % wopos( 2,49)  = 'SCALA'
     postp(1) % wopos( 2,50)  = 'SCALA'
     postp(1) % wopos( 2,51)  = 'SCALA'
     postp(1) % wopos( 2,52)  = 'SCALA'
     postp(1) % wopos( 2,53)  = 'SCALA'
     postp(1) % wopos( 2,54)  = 'SCALA'
     postp(1) % wopos( 2,55)  = 'SCALA'
     postp(1) % wopos( 2,56)  = 'SCALA'
     postp(1) % wopos( 2,57)  = 'SCALA'
     postp(1) % wopos( 2,58)  = 'VECTO'
     postp(1) % wopos( 2,59)  = 'SCALA'
     postp(1) % wopos( 2,60)  = 'VECTO'
     postp(1) % wopos( 2,64)  = 'VECTO'
     postp(1) % wopos( 2,67)  = 'SCALA'
     postp(1) % wopos( 2,68)  = 'SCALA'
     postp(1) % wopos( 2,69)  = 'SCALA'
     postp(1) % wopos( 2,70)  = 'SCALA'
     postp(1) % wopos( 2,71)  = 'VECTO'
     postp(1) % wopos( 2,72)  = 'SCALA'
     postp(1) % wopos( 2,73)  = 'SCALA'
     postp(1) % wopos( 2,74)  = 'SCALA'

     postp(1) % wopos( 3,22)  = 'NELEM'
     postp(1) % wopos( 3,31)  = 'NBOUN'
     postp(1) % wopos( 3,32)  = 'NELEM'
     postp(1) % wopos( 3,35)  = 'NELEM'
     postp(1) % wopos( 3,36)  = 'NELEM'
     postp(1) % wopos( 3,37)  = 'NELEM'
     postp(1) % wopos( 3,53)  = 'NELEM'
     postp(1) % wopos( 3,56)  = 'NELEM'
     postp(1) % wopos( 3,57)  = 'NELEM'
     postp(1) % wopos( 3,62)  = 'NELEM'
     postp(1) % wopos( 3,65)  = 'NELEM'
     postp(1) % wopos( 3,66)  = 'NELEM'
     postp(1) % wopos( 3,67)  = 'NELEM'
     postp(1) % wopos( 3,72)  = 'NBOUN'

     ! NOT ONLY ONCE
     postp(1) % kfl_oonce(4:6)   = 0
     postp(1) % kfl_oonce(8)     = 0
     postp(1) % kfl_oonce(13)    = 0
     postp(1) % kfl_oonce(14)    = 0
     postp(1) % kfl_oonce(15:20) = 0
     postp(1) % kfl_oonce(26)    = 0
     postp(1) % kfl_oonce(34)    = 0
     postp(1) % kfl_oonce(28)    = 0
     postp(1) % kfl_oonce(40)    = 0
     postp(1) % kfl_oonce(41)    = 0
     postp(1) % kfl_oonce(43)    = 0
     postp(1) % kfl_oonce(44)    = 0
     postp(1) % kfl_oonce(47)    = 0
     postp(1) % kfl_oonce(50)    = 0
     postp(1) % kfl_oonce(60)    = 0
     postp(1) % kfl_oonce(68)    = 0
     postp(1) % kfl_oonce(69)    = 0
     postp(1) % kfl_oonce(71)    = 0
     postp(1) % kfl_oonce(72)    = 0
     postp(1) % kfl_oonce(73)    = 0

     !postp(1) % kfl_oonce(24)    = 0
     !
     ! Nullify pointers
     !
     nullify(lnodb_mm)
     nullify(coord_mm)
     nullify(tncod_ker)
     nullify(tgcod_ker)
     nullify(tbcod_ker)
     nullify(cowit)
     nullify(uwall_ker)
     nullify(uwal2_ker)
     nullify(shwit)
     nullify(dewit)
     nullify(displ_ker)
     nullify(lewit)
     nullify(kfl_funno_walld_ker)
     nullify(kfl_funbo_walld_ker)
     nullify(kfl_fixno_walld_ker)
     nullify(kfl_fixbo_walld_ker)
     nullify(kfl_funty_walld_ker)
     nullify(kfl_fixno_walln_ker)
     nullify(kfl_fixbo_walln_ker)
     nullify(funpa_walld_ker)
     nullify(bvess_walld_ker)
     nullify(bvnat_walld_ker)
     nullify(bvess_walln_ker)
     nullify(kfl_funno_defor_ker)
     nullify(kfl_fixno_defor_ker)
     nullify(kfl_funty_defor_ker)
     nullify(bvess_defor_ker)
     nullify(kfl_funno_rough_ker)
     nullify(kfl_funbo_rough_ker)
     nullify(kfl_fixno_rough_ker)
     nullify(kfl_fixbo_rough_ker)
     nullify(kfl_fixbo_nsw_ker)
     nullify(kfl_funty_rough_ker)
     nullify(funpa_rough_ker)
     nullify(bvess_rough_ker)
     nullify(bvnat_rough_ker)
     nullify(kfl_fixno_suppo_ker)
     nullify(bvess_suppo_ker)
     nullify(fact_nsw_ker)
     nullify(avupo_ker)
     nullify(lnsw_exch)

     do ifunc = 1,max_space_time_function
        nullify( space_time_function(ifunc) % expression )
     end do
     do ifunc = 1,max_time_function
        nullify( time_function(ifunc) % parameters )
     end do
     nullify(subdomain)
     !
     ! Solvers - kernel problems
     !
     call moddef( 9_ip)
     call soldef(-6_ip)   ! Allocate memory for 6 solvers

     solve(1) % ndofn     = 1
     solve(1) % kfl_solve = 1
     solve(1) % wprob     = 'ROUGHNESS'

     solve(2) % ndofn     = 1
     solve(2) % kfl_solve = 1
     solve(2) % wprob     = 'WALL_DISTANCE'

     solve(3) % ndofn     = ndime
     solve(3) % kfl_solve = 1
     solve(3) % wprob     = 'SUPPORT_GEOMETRY'

     solve(4) % ndofn     = ndime
     solve(4) % kfl_solve = 1
     solve(4) % wprob     = 'MESH_DEFORMATION'

     solve(5) % ndofn     = 1
     solve(5) % kfl_solve = 1
     solve(5) % wprob     = 'WALL_NORMAL'

     solve(6) % ndofn     = 1
     solve(6) % kfl_solve = 1
     solve(6) % wprob     = 'GROUPS'
     solve(6) % kfl_iffix = 1
     !
     ! Laws for properties
     !
     call ker_allaws()
     call ker_mempro(1_ip)
     call ker_proper_initialization()
     !
     ! Others
     !
     number_event = 0
     !
     ! Get events directory name
     !
     call ker_events_directory_name()
     call ker_subdomain_initialization()

  case ( 1_ip )
     !
     ! Redefine MGAUS if there are cut elements
     !
     if( kfl_cutel == 1 ) then
        if( ndime == 2 ) then
           if( lexis(TRI03) == 0 ) call runend('CDERDA: WHEN USING CUT ELEMENTS, DECLARE TRI03 ELEMENTS')
           mgaus = max(mgaus,9*ngaus(TRI03))
        else
           if( lexis(HEX08) /= 0 ) then
              mgaus = max(mgaus,36*ngaus(TET04))
           else
              if( lexis(TET04) == 0 ) call runend('CDERDA: WHEN USING CUT ELEMENTS, DECLARE TET04 ELEMENTS')
              mgaus = max(mgaus,6*ngaus(TET04))
           end if
        end if
     end if
     !
     ! Initialize space/time functions
     !
     if( INOTMASTER ) call ker_init_space_time_function()
     !
     ! Vector size
     !
     if( kfl_vector_size /= 0 ) then
#ifdef VECTOR_SIZE_VARIABLE
        VECTOR_SIZE = kfl_vector_size
#endif
     end if

  case ( 2_ip )
     !
     ! Solver arrays
     !
     solve(1) % kfl_fixno => kfl_fixno_rough_ker
     solve(2) % kfl_fixno => kfl_fixno_walld_ker
     solve(3) % kfl_fixno => kfl_fixno_suppo_ker
     solve(5) % kfl_fixno => kfl_fixno_walln_ker
     !
     ! Allocate memory for cut elements
     !
     if( kfl_cutel == 1 ) then
        if( INOTMASTER ) call buicut(1_ip)
        mnoga = max(mnode,mgaus)
     end if

  end select

end subroutine ker_inivar
