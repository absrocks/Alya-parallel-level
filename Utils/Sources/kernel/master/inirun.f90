subroutine inirun()
  !-----------------------------------------------------------------------
  !****f* master/inirun
  ! NAME
  !    inirun
  ! DESCRIPTION
  !    This subroutine initializes and defines some run parameters
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_meshin
  use def_master
  use def_solver
  use def_postpr
  use def_inpout
  use def_coupli
  use mod_opebcs
  use def_coupli
  use mod_memory,             only : lun_memor
  use mod_memory,             only : lun_varcount
  use mod_memory,             only : memory_initialization
  use mod_postpr,             only : postpr_initialization
  use mod_elsest,             only : elsest_initialization
  use mod_random,             only : random_initialization
  use mod_couplings,          only : couplings_initialization
  use mod_htable,             only : htable_initialization
  use mod_alya_direct_solver, only : alya_direct_solver_initialization
  use mod_messages,           only : messages_initialization
  use mod_alya2metis,         only : alya2metis_initialization
  use mod_elmgeo,             only : elmgeo_element_type_initialization
  use mod_timings,            only : timings_initialization
  use mod_mpio_seq_postpr,    only : mod_mpio_seq_postpr_initialization

  implicit none
  integer(ip) :: dummi,imodu,jmodu,ii

  call cputim(cpu_initi)         ! Initial CPU time
  !
  ! Initializaiton of some modules
  !
  call memory_initialization()
  call postpr_initialization()
  call elsest_initialization()
  call random_initialization()
  call couplings_initialization()

  call alya_direct_solver_initialization()
  call htable_initialization(htable_lninv_loc)
  call messages_initialization()
  call alya2metis_initialization()
  call elmgeo_element_type_initialization()
  call timings_initialization()
  call mod_mpio_seq_postpr_initialization()
  !
  ! Units: General
  !
  lun_pdata     = 11             ! Data file unit
  lun_outpu     = 12             ! Output (log) file unit
  lun_memor     = 13             ! Memory file unit
  lun_conve     = 14             ! Convergence file unit
  lun_rstar     = 17             ! Restart file unit
  lun_latex     = 18             ! Latex file unit: text
  lun_gnupl     = 19             ! Latex file unit: gnuplot
  lun_commu     = 20             ! Communication with Alya
  lun_binar     = 24             ! Geometry binary file
  lun_pdata_dom = 21             ! Domain data file unit
  lun_elsta_dom = 25             ! Elsest statistics
  lun_elmsh_dom = 26             ! Elsest mesh
  lun_elres_dom = 27             ! Elsest results
  lun_syste     = 28             ! System info
  lun_tempo     = 29             ! Temporary unit
  lun_rstib     = 45             ! IB Restart file unit
  lun_detec     = 33             ! Automatic detection file
  lun_timeline  = 42             ! Timeline file
  lun_memory    = 34             ! Memory evolution
  lun_varcount  = 50             ! Memory variable counter file unit
  !  lun_time      = 50             ! Time values file unit
  !
  ! Units: Domain postprocess
  !
  lun_outpu_dom = 22             ! Output domain file unit
  lun_postp     = 15             ! Postprocess domain unit
  lun_posvx     = 16             ! Postprocess voxel unit
  lun_pos00     = 30             ! Additional output file (VU)
  lun_pos01     = 31             ! Additional output file (VU)
  lun_pos02     = 32             ! Additional output file (VU)
  !kfl_oumes     =  0             ! Do not output
  !kfl_oumes(1)  =  1             ! Output mesh
  !
  ! Units: Set postprocess
  !
  lun_quali     = 35             ! Mesh quality
  !
  ! Units: Set postprocess
  !
  lun_pos09     = 43             ! Additional output file (VU) (Filter)
  lun_pos10     = 44             ! Additional output file (VU) (Filter)
  !
  ! IB
  !
  lun_mshib     = 46             ! IB mesh
  lun_resib     = 47             ! IB results
  lun_mshi2     = 48             ! IB mesh (2)
  lun_resi2     = 49             ! IB results (2)
  !
  ! Lagrangian particles
  !
  lun_rstla     = 39
  lun_posla     = 40
  lun_cvgla     = 41
  !
  ! Coupling
  !
  lun_coupl_dat = 36
  lun_coupl_res = 37
  lun_coupl_cvg = 38
  !
  ! Read/write units
  !
  lun_pdata_dom = 21             ! Domain data file unit
  lun_outpu_dom = 22             ! Output domain file unit
  lispa         = 0              ! 0 passes through ecoute
  lisda         = lun_pdata      ! Temporary data file
  lisre         = lun_outpu      ! Results file
  endst         = 1              ! Stop Alya of end of file found
  kfl_split_plus= 0              ! By default the symbol + is not a separator
  !
  ! Memory
  !
  iar3p        = 0               ! ger3p not allocated
  iasca        = 0               ! gesca not allocated
  iavec        = 0               ! gevec not allocated
  iaten        = 0               ! geten not allocated
  !
  ! Master
  !
  cpu_outpu     = 0.0_rp         ! Output CPU
  cpu_other     = 0.0_rp         ! CPU's
  ittim         = 0              ! First time step
  itti2         = 0              ! First time step
  kfl_goopt     = 1              ! Go in optimization
  kfl_gotim     = 1              ! Go in time
  kfl_gocou     = 1              ! Go in coupling iterations
  kfl_reset     = -1             ! Reset step -1 is off, 0 is on but not required, 1 is do reset
  memke         = 0              ! Current and maximum memory
  dtold         = 0.0_rp         ! Old time step
  dtime         = 0.0_rp         ! Time step
  ittyp         = 0              ! We are in initial run
  isect         = 0              ! Live output sections
  inews         = 0              ! Live output New section
  file_opened   = .true.         ! File was openend successfully
  nspec         = 0              ! Species
  !
  ! Parall service
  !
  kfl_paral     = -1             ! Alya not initiated by MPI
  kfl_ptask     =  1             ! 0=Master ONLY partition
  kfl_outpu_par =  0             ! Do not output parall slave info
  kfl_postp_par =  1             ! Postprocess in Master postprocess file
  nzone_par =  1                 ! Only one zone
  icoml         =  1             ! Current level number
  npart         =  1             ! Sequential run
  npasi         =  0_ip
  npari         =  0_ip
  nparx         =  0_ip
  npasr         =  0_ip
  nparr         =  0_ip
  nparl         =  0_ip
  call vocabu(-1_ip,0_ip,0_ip)
  !
  ! Solvers
  !
  nusol         = 1              ! Number of solves
  mxdof         = 0              ! D.o.f. per node
  memit         = 0              ! Iterative solver memory
  smemo         = 0              ! Solver memory
  kfl_symgr     = 0              ! Symmetric graph not needed (define in ***_inivar)
  kfl_schur     = 0              ! No Schur complement solver exists
  kfl_aiipr     = 0              ! No Aii preconditioner exists

  nzmat         = 1              ! # Max matrix size over all modules
  nzmbt         = 1              ! # Max RHS eigen matrix size over all modules
  nzrhs         = 1              ! # Max RHS size over all modules
  nzpre         = 1              ! # Max Preconditioner size over all modules

  nzmax         = 1              ! = max(nzsol,nzsky,nzexp): Components of A
  nzrhx         = 1              ! = max(nzsol,nzsky,nzexp): Components of RHS
  nzprx         = 1              ! Components of Preconditioner

  neige         = 1              ! # Max eigenvalue vector size
  neiva         = 1              ! # Max eignvalues
  nzerr         = 1              ! # Max Error Estimator size over all modules
  !
  ! Domain
  !
  kfl_conma     = 0              ! Consistent mass not nedded
  nzsky         = 1              ! # of comp. in the skyline matrix of the graph
  nzsol         = 1              ! # of comp. in the CSR matrix = nzdom
  mpopo         = 0              ! # Max Node-element connectivity
  iffun         = 0              ! Read bc codes (no function)
  ifloc         = 0              ! Read bc codes (no local bc)
  ifbop         = 0              ! Read bc codes (not on boundary nodes)
  ifbes         = 1              ! Read boundary values
  !
  ! Variables read in readat
  !
  kfl_elses  = 0                 ! Elsest
  xscal(1)   = 1.0_rp            ! X scale factor
  xscal(2)   = 1.0_rp            ! Y scale factor
  xscal(3)   = 1.0_rp            ! Z scale factor
  ielse(1)   = 100               ! nx
  ielse(2)   = 100               ! ny
  ielse(3)   = 100               ! nz
  ielse(4)   = 0                 ! data format (0=type,1=list)
  ielse(5)   = 10                ! Maximum number of possible meshes
  ielse(6)   = 2                 ! Second try strategy: if box is not found
  ielse(7)   = 0                 ! Output unit
  ielse(8)   = 1                 ! Search strategy (0=bin,1=Quad)
  ielse(9)   = 100               ! Points per node for Quad/Octtree
  ielse(10)  = 0                 ! Result output frequency
  ielse(11)  = 0                 ! Neighboring-boxes-search radius, 0: search until all boxes finished
  ielse(12)  = 0                 ! Postprocess mesh
  ielse(13)  = 0                 ! Postprocess results
  relse      = 0.0_rp            ! Elsest
  relse(1)   = 0.01_rp           ! Tolerance for iteration
  !
  ! Global variables read by modules
  !
  kfl_coibm     = 0              ! Immbou coupling: Read by IMMBOU.
  kfl_advec     = 0              ! Mesh advection: Read by IMMBOU.
  kfl_async     = 1              ! Parall: asynchronous communications
  kfl_forca_res = 0              ! Forces in residual formulation for rigid body(alefor decides)
  thicl         = 0.0_rp         ! Interface thicknes. Read by LEVELS.
  mcono         = 3              ! Max # codes per nodes
  !
  ! Required arrays
  !
  kfl_lface     = 0              ! List of global faces not required: LFACG
  kfl_lelp2     = 0              ! List of extended node-element graph: PELPO_2, LELPO_2
  kfl_lelbf     = 0              ! List of element boundary faces: LELBF
  kfl_symgr     = 0              ! If symmetric graph is needed
  kfl_conma     = 0              ! If consistent mass is needed
  kfl_element_bin = 0            ! If element bin is required
  !
  ! Meshing
  !
  kfl_ifbox     = 0              ! If bounding box is prescribed
  !
  ! Modules
  !
  do modul = 0,mmodu
     momod(modul) % cpu_modul = 0.0_rp
     momod(modul) % mem_modul = 0
     momod(modul) % kfl_delay = 0
     momod(modul) % kfl_solve = AT_EACH_TIME_STEP
     nullify(momod(modul) % solve)
     nullify(momod(modul) % solad)
     nullify(momod(modul) % postp)
     nullify(momod(modul) % eigen)
  end do
  modul         = 0                       ! I am kernel
  kfl_timei     = 0                       ! Steady calculation
  cpu_modul     = 0.0_rp                  ! Module CPU
  cpu_modul(CPU_MINI_ASSEMBLY,:) = 1e6_rp ! Min CPU, to remove variablity effects
  mem_modul     = 0                       ! Module memory to zero
  namod         = '      '                ! Module names
  exmod         = '   '                   ! Module extension
  kfl_delay     = 0                       ! Delay module
  kfl_conve     = 1                       ! Module convergence required
  kfl_modul     = 0                       ! Module not used
  kfl_solve     = AT_EACH_TIME_STEP       ! When module is solved
  ndela         = 0                       ! Steps to delay module
  itinn         = 0                       ! First internal iteration
  kfl_modul(0)  = 1                       ! Kernel is always On!
  kfl_modul(mmodu) = 1                    ! Kermod is always On!
  do imodu = 0,mmodu
     lzone(imodu) = 1
     do jmodu = 0,mmodu
        kfl_coupl(jmodu,imodu) = 0
        kfl_cowhe(jmodu,imodu) = 0
     end do
     do ii = 1,20
        kfl_itask(ii,imodu) = 0
     end do
  end do
  namod(mmodu)  = 'KERMOD'
  exmod(mmodu)  = 'ker'
  namod(0)      = 'KERNEL'
  exmod(0)      = 'ke2'
  namod(1)      = 'NASTIN'
  exmod(1)      = 'nsi'
  namod(2)      = 'TEMPER'
  exmod(2)      = 'tem'
  namod(3)      = 'CODIRE'
  exmod(3)      = 'cdr'
  namod(4)      = 'TURBUL'
  exmod(4)      = 'tur'
  namod(5)      = 'EXMEDI'
  exmod(5)      = 'exm'
  namod(6)      = 'NASTAL'
  exmod(6)      = 'nsa'
  namod(7)      = 'ALEFOR'
  exmod(7)      = 'ale'
  namod(8)      = 'LATBOL'
  exmod(8)      = 'lat'
  namod(9)      = 'APELME'
  exmod(9)      = 'ape'
  namod(10)     = 'SOLIDZ'
  exmod(10)     = 'sld'
  namod(11)     = 'GOTITA'
  exmod(11)     = 'got'
  namod(12)     = 'WAVEQU'
  exmod(12)     = 'wav'
  namod(14)     = 'LEVELS'
  exmod(14)     = 'lev'
  namod(15)     = 'QUANTY'
  exmod(15)     = 'qua'
  namod(16)     = 'MAGNET'
  exmod(16)     = 'mag'
  namod(17)     = 'PARTIS'
  exmod(17)     = 'pts'
  namod(18)     = 'NASEDG'
  exmod(18)     = 'nsg'
  namod(19)     = 'CHEMIC'
  exmod(19)     = 'chm'
  namod(20)     = 'HELMOZ'
  exmod(20)     = 'hlm'
  namod(21)     = 'IMMBOU'
  exmod(21)     = 'ibm'
  namod(22)     = 'RADIAT'
  exmod(22)     = 'rad'
  namod(23)     = 'CASIMI'
  exmod(23)     = 'cas'
  namod(24)     = 'POROUS'
  exmod(24)     = 'por'
  namod(24)     = 'XXXXXX'
  exmod(24)     = 'xxx'
  namod(26)     = 'NEUTRO'
  exmod(26)     = 'neu'
  namod(27)     = 'INSITU'
  exmod(27)     = 'ins'
  namod(28)     = 'SOLFE2'
  exmod(28)     = 'fe2'
  !
  ! Services
  !
  cpu_servi     = 0.0_rp         ! Service CPU
  mem_servi     = 0              ! Service memory to zero
  naser         = '      '       ! Service names
  exser         = '   '          ! Service extension
  naser(1)      = 'SHAPAR'
  exser(1)      = 'shp'
  naser(2)      = 'SOLMUM'
  exser(2)      = 'slm'
  naser(3)      = 'GIDPOS'
  exser(3)      = 'gid'
  naser(4)      = 'HANDFP'
  exser(4)      = 'hdp'
  naser(5)      = 'PARALL'
  exser(5)      = 'par'
  naser(6)      = 'CGNS24'
  exser(6)      = 'cgn'
  naser(7)      = 'DODEME'
  exser(7)      = 'dod'
  naser(9)      = 'ADAPTI'
  exser(9)      = 'ada'
  naser(10)     = 'ARPACK'
  exser(10)     = 'apa'
  naser(11)     = 'LAPACK'
  exser(11)     = 'lpa'
  naser(12)     = 'HDFPOS'
  exser(12)     = 'hdf'
  naser(13)     = 'OPTSOL'
  exser(13)     = 'opt'
  !
  ! Postprocess
  !
  call posdef(0_ip,dummi)
  kfl_reawr             = 0       ! Postprocess mode
  ivapo                 = 0
  kfl_ivari             = 0
  modul                 = 0       ! Now we deal with kernel's postprocess
  call moddef(9_ip)
  postp(1) % kfl_oonce = 1        ! Kernel postprocess is only once by defalut
  postp(1) % npp_iniso = 1        ! Kernel postprocess always initial values
  !
  ! Geometry
  !
  postp(1) % wopos( 1,15)  = 'LNODS'
  postp(1) % wopos( 1,16)  = 'COORD'
  postp(1) % wopos( 1,17)  = 'LTYPE'
  postp(1) % wopos( 1,18)  = 'LNINV'
  postp(1) % wopos( 1,19)  = 'LELCH'
  postp(1) % wopos( 1,20)  = 'LNODB'
  postp(1) % wopos( 1,21)  = 'LTYPB'
  postp(1) % wopos( 1,22)  = 'LESUB'
  postp(1) % wopos( 1,23)  = 'LEINV'
  postp(1) % wopos( 1,24)  = 'LMATE'
  postp(1) % wopos( 1,25)  = 'LBINV'
  postp(1) % wopos( 1,26)  = 'LBOCH'
  postp(1) % wopos( 1,27)  = 'LNOCH'
  postp(1) % wopos( 1,33)  = 'LELBO'
  postp(1) % wopos( 1,75)  = 'LMAST'

  !
  ! Sets
  !
  postp(1) % wopos( 1,30)  = 'LESET'
  postp(1) % wopos( 1,31)  = 'LBSET'
  !
  ! Boundary conditions
  !
  postp(1) % wopos( 1,28)  = 'CODBO'
  postp(1) % wopos( 1,29)  = 'CODNO'
  !
  ! Fields
  !
  postp(1) % wopos( 1,32)  = 'XFIEL'


  postp(1) % wopos( 2,15)  = 'VECTO'
  postp(1) % wopos( 2,16)  = 'VECTO'
  postp(1) % wopos( 2,17)  = 'SCALA'
  postp(1) % wopos( 2,18)  = 'SCALA'
  postp(1) % wopos( 2,19)  = 'SCALA'
  postp(1) % wopos( 2,20)  = 'VECTO'
  postp(1) % wopos( 2,21)  = 'SCALA'
  postp(1) % wopos( 2,22)  = 'SCALA'
  postp(1) % wopos( 2,23)  = 'SCALA'
  postp(1) % wopos( 2,24)  = 'SCALA'
  postp(1) % wopos( 2,25)  = 'SCALA'
  postp(1) % wopos( 2,26)  = 'SCALA'
  postp(1) % wopos( 2,27)  = 'SCALA'

  postp(1) % wopos( 2,28)  = 'SCALA'
  postp(1) % wopos( 2,29)  = 'VECTO'

  postp(1) % wopos( 2,30)  = 'SCALA'
  postp(1) % wopos( 2,31)  = 'SCALA'

  postp(1) % wopos( 2,32)  = 'VECTO'
  postp(1) % wopos( 2,33)  = 'SCALA'
  postp(1) % wopos( 2,75)  = 'SCALA'

  postp(1) % wopos( 3,15)  = 'NELEM'
  postp(1) % wopos( 3,16)  = 'NPOIN'
  postp(1) % wopos( 3,17)  = 'NELEM'
  postp(1) % wopos( 3,18)  = 'NPOIN'
  postp(1) % wopos( 3,19)  = 'NELEM'
  postp(1) % wopos( 3,20)  = 'NBOUN'
  postp(1) % wopos( 3,21)  = 'NBOUN'
  postp(1) % wopos( 3,22)  = 'NELEM'
  postp(1) % wopos( 3,23)  = 'NELEM'
  postp(1) % wopos( 3,24)  = 'NELEM'
  postp(1) % wopos( 3,25)  = 'NBOUN'
  postp(1) % wopos( 3,26)  = 'NBOUN'
  postp(1) % wopos( 3,27)  = 'NPOIN'

  postp(1) % wopos( 3,28)  = 'NBOUN'
  postp(1) % wopos( 3,29)  = 'NPOIN'

  postp(1) % wopos( 3,30)  = 'NELEM'
  postp(1) % wopos( 3,31)  = 'NBOUN'

  postp(1) % wopos( 3,32)  = 'UNKNO'
  postp(1) % wopos( 3,33)  = 'NBOUN'
  postp(1) % wopos( 3,75)  = 'NPOIN'
  !
  ! nullify pointers
  !
  nullify(enthalpy_transport)
  nullify(chemical_heat)
  nullify(radiative_heat)
  nullify(table_cfi)
  nullify(div_enthalpy_transport)
  nullify(tfles_factor)
  nullify(tfles_sensor)
  nullify(tfles_sgseff)
  nullify(gefil)
  nullify(veloc)
  nullify(press)
  nullify(tempe)
  nullify(densi)
  nullify(energ)
  nullify(visco)
  nullify(umome)
  nullify(untur)
  nullify(uncdr)
  nullify(elmag)
  nullify(dispm)
  nullify(velom)
  nullify(displ)
  nullify(spins)
  nullify(disfu)
  nullify(vdrop)
  nullify(cdrop)
  nullify(wavam)
  nullify(fleve)
  nullify(erres)
  nullify(vorti)
  nullify(conce)
  nullify(cdens)
  nullify(entha)
  nullify(therm)
  nullify(phion)
  nullify(fiber)
  nullify(fisoc)
  nullify(gpfib)
  nullify(radso)
  nullify(vconc)
  nullify(taulo)
  nullify(kfl_fixno_ale)
  nullify(bvess_ale)
  nullify(gisca)
  nullify(givec)
  nullify(giscp)
  nullify(givep)
  nullify(geten)
  nullify(gevec)
  nullify(gesca)
  nullify(getep)
  nullify(gevep)
  nullify(gescp)
  nullify(getex)
  nullify(gevex)
  nullify(gescx)
  nullify(ger3p)
  nullify(turmu)
  nullify(rhoon)
  nullify(forcf)
  nullify(forca)
  nullify(sphea)
  nullify(kcond)
  nullify(wmean)
  nullify(visck)
  nullify(massk)
  nullify(lescl)
  nullify(encfi)
  nullify(momen)

  nullify(condk)
  nullify(sphek)
  nullify(advec)
  nullify(sphec)
  nullify(vesgs)
  nullify(tesgs)
  nullify(cosgs)
  nullify(veset)
  nullify(vbset)
  nullify(vnset)
  nullify(witne)
  nullify(unkno)
  nullify(eigen)
  nullify(eigva)
  nullify(amatr)
  nullify(bmatr)
  nullify(rhsid)
  nullify(pmatr)
  nullify(pschu)
  nullify(aii)
  nullify(aib)
  nullify(abi)
  nullify(abb)
  nullify(xxi)
  nullify(xxb)
  nullify(bbi)
  nullify(bbb)
  nullify(lumma)
  nullify(amatx)
  nullify(rhsix)
  nullify(unknx)
  nullify(pmatx)
  nullify(parin)
  nullify(pari1)
  nullify(pari2)
  nullify(pari3)
  nullify(paris)
  nullify(parig)
  nullify(lnbin)
  nullify(lgpar)
  nullify(lnwit)
  nullify(ledgg)
  nullify(lfacg)
  nullify(ledgb)
  nullify(lfacb)
  nullify(ledgp)
  nullify(pedgp)
  nullify(parlo)
  nullify(parhh)
  nullify(nelem_tot)
  nullify(npoin_tot)
  nullify(nboun_tot)
  nullify(npoia_tot)
  nullify(npoin_par)
  nullify(nelem_par)
  nullify(nboun_par)
  nullify(lninv_loc)
  nullify(leinv_loc)
  nullify(lpoi4)
  nullify(recvbuf_ip)
  nullify(sendbuf_ip)
  nullify(recvbuf_rp)
  nullify(sendbuf_rp)
  nullify(displs)
  nullify(recvcounts)

  nullify(iaren_par)
  nullify(jaren_par)
  nullify(permr_par)
  nullify(invpr_par)
  nullify(lelbf)
  nullify(lelfa)
  nullify(parre)
  nullify(parr1)
  nullify(parrs)
  nullify(parr2)
  nullify(parr3)
  nullify(parcx)
  nullify(parx1)
  nullify(parx2)
  nullify(parx3)
  nullify(par3p)
  nullify(pai1p)
  nullify(kfl_fixno_dod)
  nullify(lnsub_dod)
  nullify(lbsub_dod)
  nullify(gesca_hdf)
  nullify(gevec_hdf)
  nullify(ger3p_hdf)
  nullify(imbou)
  nullify(rbbou)
  nullify(lnint)
  nullify(lntib)
  nullify(lnti2)
  nullify(lndom)
  nullify(lntra)
  nullify(letib)
  nullify(massc)
  nullify(statt)

  nullify(gescar_mpio)
  nullify(gevecr_mpio)
  nullify(gescai_mpio)
  nullify(geveci_mpio)

  nullify(I_AM_IN_ZONE)
  nullify(I_AM_IN_SUBD)
  nullify(I_AM_IN_CODE)
  nullify(application_names)
  !
  ! Code/zone/subdomain (essentially used by service COUPLI)
  !
  !current_zone = lzone(ID_KERMOD)
  current_zone = 0
  current_code = 1
  current_subd = 1
  kfl_gozon    = 1
  !
  ! Others
  !
  IEMPTY    = .false.
  INOTEMPTY = .true.
  !
  ! Solvers - kernel problems
  !
  !  modul = 0
  !  call moddef( 9_ip)
  !  call soldef(-2_ip)                                ! Allocate memory
  !  solve(1) % kfl_solve = 1
  !  solve(1) % wprob     = 'ROUGHNESS'
  !  solve(2) % kfl_solve = 1
  !  solve(2) % wprob     = 'WALL_DISTANCE'
  !  !
  !  ! Boundary conditions
  !  !
  !  ! Node
  !  call opnbcs(1_ip,2_ip,dummi,dummi,tncod_ker)      ! Memory for structure
  !  call opnbcs(2_ip,1_ip, 1_ip, 0_ip,tncod_ker)      ! Memory for mesh deformation/smoothing
  !  call opnbcs(2_ip,2_ip, 1_ip, 1_ip,tncod_ker)      ! Memory for wall distance
  !  ! Boundary
  !  call opbbcs(0_ip,1_ip, 1_ip,tbcod_ker)            ! Memory for mesh deformation/smoothing
  !  call opbbcs(0_ip,2_ip, 1_ip,tbcod_ker)            ! Memory for wall distance
  !  ! Geometrical
  !  call opnbcs(0_ip,1_ip, 1_ip, 0_ip,tgcod_ker)      ! Memory for mesh deformation/smoothing
  !  call opnbcs(0_ip,2_ip, 1_ip, 0_ip,tgcod_ker)      ! Memory for wall distance

end subroutine inirun
