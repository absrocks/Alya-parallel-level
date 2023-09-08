subroutine nsa_reabcs
!-----------------------------------------------------------------------
!****f* Nastal/nsa_reabcs
! NAME 
!    nsa_reabcs
! DESCRIPTION
!    This routine reads the boundary conditions for the NS
!    equation.
!
!    For conditions on nodes, bvess_nsa(idofn,ipoin,1) contains the
!    imposed values for u, rho (or p) and T, 
!    which corresponds to the momentum, continuity and energy equations 
!    respectively. 
!    The codes for the momentum equations (u) are:
!    = 0 ... Free or initial
!    = 1 ... Fixing on global cartesian basis (GCB)
!    = 2 ... Fixing on local curvilinear basis (LCB)
!    and the special code for the first velocity component
!    = 5 ... SPECIAL
!    Iffix 5 determines the rest of the conditions of the node depending on the 
!    inflow velocity angle. Normal inflows and outflows are considered
!    as 5's when INLETS are set to REFERENCE in the dat file.
!    The codes for the continuity equation are (rho or p):
!    = 0 ... Free or initial
!    = 1 ... Fixing rho
!    = 2 ... Fixing p
!    The codes for the energy equation are (T):
!    = 0 ... Free or initial
!    = 1 ... Fixing T
!    Alternatively, the condition imposed on a boundary
!    element can modify a GCB (only a GCB!) condition, flipping it to LCB. 
!
!    It is also possible to assign bvess values (and then fixnos) to
!    two other variables: density (or pressure) and temperature. This is
!    is the reason why bvess_nsa is dimensioned as ndofn_nsa X npoin.
!    For instance, the stagnation nodes on Euler compressible flows:
!
!                  11101    0.0 0.0 0.0   0.0   Tstag
!
!    For conditions on boundaries, bvnat_nsa(iboun) is allocated ONLY if
!    a condition is imposed on it. Its length depends on the
!    condition type. 
!
!    Nastal module can accept nstinc-like fixities files, providing that
!    this fact is explicitly warned with the keyword NASTIN just before the
!    fixities themselves in the *.nsa.dat file. For instance:
!
!    BOUNDARY_CONDITIONS, NASTIN, CONSTANT
!      INCLUDE my-old-nstinc-file-of-fixities.fix
!    END_BOUNDARY_CONDITIONS
!
!    When prescribing boundary conditions on boundary faces, the following
!    codes are used:
!
!    1. INHERITED FROM NASTIN, USED IN FORCED INCOMPRESSIBLE FLOW (kfl_foreg = 1):
!    = 1 ... Dirichlet ............ u
!    = 2 ... Pressure imposed  .... sig.n=-p n
!    = 3 ... Wall law ............. sig.n.t=-rho*U_*^2*u/|u|
!    = 4 ... Symmetry/Slip wall ... sig.n.t=0, u.n=0
!    = 5 ... Dynamic pressure ..... sig.n=-1/2*rho*u^2
!    = 6 ... Open flow ............ Assemble 2*mu*Sym(grad(u).n
!    = 7 ... No slip wall ......... u=0
!
!    2. SPECIFIC FOR NASTAL, USED IN FORCED COMPRESSIBLE FLOW (kfl_foreg = 0):
!    = 20 ... Free, force to check . rho or p or nothing; forced to be set according 
!                                    to local Mach number. Sometimes for OUTFLOWS.
!    = 21 ... Dirich, check mach ... u,rho,T ; set according to local 
!                                    Mach number. Typical INFLOWS.
!    = 22 ... Pressure imposed  .... p
!    = 23 ... Wall law ............. sig.n.t=-rho*U_*^2*u/|u|
!    = 24 ... Symmetry/Slip wall ... sig.n.t=0, u.n=0
!    = 25 ... Dynamic pressure ..... sig.n=-1/2*rho*u^2
!    = 26 ... Open flow ............ Assemble 2*mu*Sym(grad(u).n
!    = 27 ... No slip wall ......... u= 0
!    = 31 ... Heat flux ............ q= k n.grad(T)
!    = ... etc ...
!
!    At the end of the subroutine conditions on boundaries are
!    transfered to conditions on nodes.
!
!    When a fixities file with nstinc-like faces-wise boundary conditions is 
!    used, the codes are converted to the proper ones when reading the file.
!
!
! USES
!    ecoute
!    memchk
!    runend
! USED BY
!    nsa_turnon
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  use      mod_opebcs
  use      def_nastal
  use      def_kermod
  use mod_ker_space_time_function
  use mod_ecoute, only      : ecoute
  use mod_messages, only : livinf

  implicit none
  integer(ip)              :: ipoin,idime,ifunc
  integer(ip)              :: ncodf,kfl_autlo,kposi
  integer(ip)              :: lbcod(50,2),licod(50,2)
  integer(ip)              :: nbcod,nicod,itloa

  integer(ip)              :: krang(50,3),nspbc(5)
       
  real(rp)                 :: rbcod(50,10),ricod(50,5),rauxi(2)
  character(300)           :: messa
  

  do ipoin=1,50
     krang(ipoin,1)=0
     krang(ipoin,2)=0
     krang(ipoin,3)=0
  end do
  kfl_dumbo_nsa=0
  kfl_funty_nsa=0
  fubcs_nsa=0.0_rp

  !
  ! Initialization
  !
  
  nfunc_nsa     = 0                                    ! No time-dependent boundary condition types
  nodpr_nsa     = 0                                    ! Node where to impose pressure
  kfl_confi_nsa = 0                                    ! Flow is not confined
  kfl_local_nsa = 0                                    ! No local system of reference
  kfl_conbc_nsa = 1                                    ! Constant boundary conditions
  kfl_bofty_nsa = 0                                    ! Nastal b.c. file type
  kfl_skewb_nsa = 0                                    ! No skew symetric in bvess
  kfl_spong_nsa = -1                                   ! No Rayleigh sponge
  kfl_inlet_nsa = 0                                    ! Activation of boundary condition to force the conservation of the flow rate (0=OFF,0/=ON)
  inlet_nsa     = -1.0_rp                              ! Target mass flow rate at the boundary
  zspon_nsa     = 0.0_rp
  

  kfl_bofie_nsa = 0                                    ! No boundary conditions from field values

  kfl_autlo = 0                                        ! Do not compute autloatic intersection

  kfl_nopen_nsa = 0                                    ! No penetration condition in strong form

  kfl_tredg_nsa = 0                                    ! Do not automatically correct trailing edges b.c.
  angle_tredg_nsa = 0.3_rp                             ! Scalar projection threshold (corresponding to 72 degrees)

  
  rtico_nsa     = 0.0_rp                               ! Initialize time-dependent b.c. time counter
  
  lbcod         = -1                                   ! Boundary codes
  nbcod         = -1
  rbcod         = 0.0_rp                               ! Boundary codes values
  
  licod         = -1                                   ! Initial cond. codes
  nicod         = -1
  ricod         = 0.0_rp                               ! Initial cond. codes values
  nspbc         = 0
  kposi         = 0                                    ! non-positional bc

  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opnbcs(0_ip,1_ip,ndofn_nsa,0_ip,tncod_nsa) 
  end if
  if( kfl_icodb > 0 ) then
     call opbbcs(0_ip,1_ip,ndofn_nsa,tbcod_nsa)      
  end if
  if( kfl_geome > 0 ) then
     call opnbcs(0_ip,1_ip,ndofn_nsa, 0_ip,tgcod_nsa)
  end if


  if( INOTSLAVE ) then
     messa = &
          '        READING BOUNDARY CONDITIONS...'
     call livinf(0_ip,messa,one)
     !
     ! Reach the nodal-wise section.
     !
     call ecoute('nsa_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('nsa_reabcs')
     end do
     !-----------------------------------------------------------------------
     ! ADOC[0]> $--------------------------------------
     ! ADOC[0]> $-- Module specific boundary conditions
     ! ADOC[0]> $--------------------------------------
     ! ADOC[0]> BOUNDARY_CONDITIONS
     !-----------------------------------------------------------------------
     !
     ! Read header.
     !
     if(exists('NASTI')) kfl_bofty_nsa = 1
     if(exists('FIELD')) then
        kfl_bofie_nsa = 1
       if (kfl_inifi_nsa(1) == 0 .or. kfl_inifi_nsa(2) == 0 .or. kfl_inifi_nsa(3) == 0 ) kfl_bofie_nsa= 0
     end if
     if(exists('FIXPR')) then
        kfl_confi_nsa =  1
        nodpr_nsa=getint('FIXPR',1_ip,'#Node where to impose pressure')
     else if(exists('DONTF')) then
        kfl_confi_nsa = -1
     end if
     if(exists('WALLD')) then
        delta_nsa=getrea('WALLD',0.0_rp,'#Distance to the wall')
     end if
     if(exists('NONCO').or.exists('TRANSI')) then               ! non-constant boundary conditions
        kfl_conbc_nsa = 0
     else                                                       ! constant boundary conditions
        kfl_conbc_nsa = 1
     end if

     if(exists('TRAIL')) then
        kfl_tredg_nsa = 1
!        angle_tredg_nsa = getrea('TRAIL',75.0_rp,'#Angle threshold for trailing edge detection')
     end if

     call ecoute('nsa_reabcs')

     if (exists('SKEWS')) then
        kfl_skewb_nsa    = 1          ! skcos links will be set in bvess
     end if
     if (exists('SPONG')) then
        kfl_spong_nsa = 1             ! rayleigh sponge
        zspon_nsa(1:3) = param(6:8)
     end if
          
     ! Adapt reading to the b.c. file origin (nastal)
     ncodf = ndofn_nsa
     if (kfl_bofty_nsa == 1) ncodf = ndime

     !
     ! Allocate the prescription time function vector of vectors
     !
     call nsa_membcs(4_ip)        

     do while(words(1)/='ENDBO')
        
        ! User-defined codes for BOUNDARY conditions
        if(words(1)=='CODES') then

           messa = &
                '        READING BOUNDARY CONDITION CODES...'
           call livinf(0_ip,messa,one)

           if(number_time_function > 0 .and. kfl_conbc_nsa /= 0 ) kfl_conbc_nsa = 0

           if (exists('NODES')) then
              !
              ! User-defined codes on nodes
              !
              messa = &
                   '        CODES DEFINED ON NODES (CODES, NODES)'
              call livinf(0_ip,messa,one)
              if (kfl_extra == 0) then
                 messa = &
                      '        **** WARNING: NO EXTRAPOLATE IS SET IN DOM.DAT FILE'
                 call livinf(0_ip,messa,one)
              end if

              if( exists('GEOME') ) then
                 !
                 ! Velocity: geometrical node code
                 !              
                 tgcod => tgcod_nsa(1:)
                 !tgcod => momod(ID_NASTAL) % tgcod(1:)
                 messa = &
                      '        CODES DEFINED ON NODES GEOMETRICALLY'
                 call livinf(0_ip,messa,one)
                 call reacod(4_ip)
              else                 
                 tncod => tncod_nsa(1:)
                 call reacod(1_ip)
              end if

           else if (exists('BOUND')) then
              !-------------------------------------------------------------
              !
              ! User-defined codes on boundaries
              !          
              !-------------------------------------------------------------
              
              messa = &
                   '        CODES DEFINED ON FACES (CODES, BOUND)'
              call livinf(0_ip,messa,one)

              tbcod => tbcod_nsa(1:)
              call reacod(2_ip)

           else 
              !-------------------------------------------------------------
              !
              ! Alya does not know what are these codes for...
              !          
              !-------------------------------------------------------------
              call runend('NSA_REABCS: CODES defined in nsa.dat file, are they for NODES or BOUNDARIES ?')

           end if

        else if(words(1)=='FUNCT') then
           call ecoute('nsa_reabcs')
           !
           ! OJOOOOOOO: A PARTIR DE AHORA, SI HAY FUNCION ENTONCES SON NON-CONSTANT, NO HACE FALTA PONERLO ARRIBA
           !
           kfl_conbc_nsa = 0                   ! nonconstant bc
           do while(words(1)/='ENDFU')
              if(kfl_conbc_nsa==0) then
                 if (words(1)=='FUNCT') then
                    ifunc=getint('FUNCT',1_ip,'#FUNCTION NUMBER')                       
                    if(ifunc<0.or.ifunc>10) then
                       call runend('nsa_reabcs: WRONG FUNCION NUMBER')
                    else
                       nfunc_nsa = nfunc_nsa + 1
                    end if
                 else if (words(1)=='TIMES') then      ! time shape
                    if(words(2)=='PARAB' .or. words(2)=='LINEA' .or. words(2)=='POLYN') then
                       kfl_funty_nsa(ifunc,1)=1
                       kfl_funty_nsa(ifunc,2)=20
                    else if(words(2)=='PERIO') then
                       kfl_funty_nsa(ifunc,1)=2
                       kfl_funty_nsa(ifunc,2)=20
                    else if(words(2)=='DISCR') then
                       kfl_funty_nsa(ifunc,1)=3
                       rauxi= 0.0_rp
                       call ecoute('nsa_reabcs')
                       if (words(1) == 'SHAPE') then   ! defining the shape by discrete points
                          rauxi(1)= getrea('REFER', 0.0_rp, 'Reference value, sometimes useful')
                          rauxi(2)= getrea('TSTAR', 0.0_rp, 'Startint time, sometimes useful')
                          call ecoute('nsa_reabcs')
                          kfl_funty_nsa(ifunc,2) = int(param(1)) - rauxi(2)
                          mtloa_nsa(ifunc) = kfl_funty_nsa(ifunc,2) - rauxi(1)
                          itloa = 0

                          call nsa_membcs(10_ip + ifunc)  ! allocate the prescription time function vector for ifunc
                          call ecoute('nsa_reabcs')

                          do while(words(1)/='ENDSH')
                             itloa= itloa + 1
                             tload_nsa(ifunc)%a(ndofn_nsa+1,itloa)= param(1)                 ! time
                             tload_nsa(ifunc)%a(1:ndofn_nsa,itloa)= param(2:ndofn_nsa+1)     ! prescribed values
                             call ecoute('nsa_reabcs') 
                          end do

                          if (itloa /= mtloa_nsa(ifunc) ) &
                               call runend('NASTAL: WRONG TOTAL OF TIME SHAPE NUMBER OF POINTS')
                          fubcs_nsa(ifunc,1)= tload_nsa(ifunc)%a(ndofn_nsa+1,1)
                          fubcs_nsa(ifunc,2)= tload_nsa(ifunc)%a(ndofn_nsa+1, mtloa_nsa(ifunc) )
                          fubcs_nsa(ifunc,3)= 1.0_rp * ifunc    ! needed to know which tload to use without passing extra parameters
                       end if
                    else
                       kfl_funty_nsa(ifunc,1)=0
                    end if

                    kfl_funty_nsa(ifunc,3:7)=1      ! default: apply to all equations 

                    if(kfl_funty_nsa(ifunc,1)>0 .and. kfl_funty_nsa(ifunc,1)/=3) then ! Read parameter functions
                       fubcs_nsa(ifunc,1:10)=param(2:11)
                    end if
                 else if (words(1)=='REFVA') then    ! reference values                       
                    fubcs_nsa(ifunc,15)=param(1)
                    fubcs_nsa(ifunc,16)=param(2)
                    fubcs_nsa(ifunc,17)=param(3)
                    fubcs_nsa(ifunc,18)=param(4)
                    fubcs_nsa(ifunc,19)=param(5)
                 else if (words(1)=='TIMEL') then    ! time lapse
                    fubcs_nsa(ifunc,11)   =param(1)  !   start
                    fubcs_nsa(ifunc,12)   =param(2)  !   end            
                    rtico_nsa(ifunc,1)    =param(1)  !   initial start time for all eqs.
                    rtico_nsa(ifunc,2)    =param(2)  !   initial final time for all eqs.
                 else if (words(1)=='TIMER') then    ! time repeat
                    fubcs_nsa(ifunc,13)=param(1)     !   when repeat
                 else if (words(1)=='EQUAT') then    ! which equation
                    kfl_funty_nsa(ifunc,3:7)= 0
                    if (exists('MOMEN')) then
                       do idime=1,ndime
                          kfl_funty_nsa(ifunc,2+idime) = 1
                       end do
                    end if
                    if (exists('XMOME')) kfl_funty_nsa(ifunc,2+1) = 1
                    if (exists('YMOME').and.ndime>1) kfl_funty_nsa(ifunc,2+2) = 1
                    if (exists('ZMOME').and.ndime>2) kfl_funty_nsa(ifunc,2+3) = 1
                    if (exists('ENERG')) kfl_funty_nsa(ifunc,2+ndime+2) = 1
                    if (exists('CONTI')) kfl_funty_nsa(ifunc,2+ndime+1) = 1
                    if (exists('ENERG')) kfl_funty_nsa(ifunc,2+ndime+2) = 1
                 else if (words(1)=='FIXIT') then    ! special behavior controlled by the fixity code
                    kfl_funty_nsa(ifunc,8)= 0
                    if (exists('PRESV')) then        ! pressure valve: alternates velocity and pressure
                       kfl_funty_nsa(ifunc,8) = 1                             ! non-slip
                       if (kfl_visco_nsa == 0) kfl_funty_nsa(ifunc,8) = 2     ! slip
                    else if (exists('VELOV')) then   ! velocity valve: alternates velocity and velocity
                       kfl_funty_nsa(ifunc,8) = 3                             ! non-slip
                       if (kfl_visco_nsa == 0) kfl_funty_nsa(ifunc,8) = 4     ! slip
                    end if
                 end if
              end if
              call ecoute('nsa_reabcs')
           end do

        else if(words(1)=='MASSF') then

           kfl_inlet_nsa = getint('BOUND',  -1_ip, 'Boundary patch')
           inlet_nsa     = getrea('FLOWR',  -0.25_rp, 'Flow rate target parameter') 

        else if(words(1)=='ONNOD') then  !DMM this is really needed?  
           !
           ! Allocate memory for the vectors needed to define the BC's 
           !
           call nsa_membcs(one)
           !
           ! iffix boundary conditions: reading the conditions
           !
           call ecoute('nsa_reabcs')
           do while(words(1)/='ENDON')
              ipoin                     = int(param(1))
              kfl_fixno_nsa(1,ipoin)    = int(param(2))
              call codfix(ncodf,kfl_fixno_nsa(1,ipoin))
              
              call ecoute('nsa_reabcs')
           end do
        end if
        !
        call ecoute('nsa_reabcs')
     end do
     !-----------------------------------------------------------------------
     ! ADOC[0]> END_BOUNDARY_CONDITIONS
     !-----------------------------------------------------------------------

     ! give a default size to non-defined tloads
     do ifunc= 1,10
        if (mtloa_nsa(ifunc) == 0) then
           mtloa_nsa(ifunc) = 1            ! done to allocate a default memory space
           call nsa_membcs(10_ip + ifunc)  ! allocate the prescription time function vector for ifunc           
        end if
     end do

     ! if at least a space-time function is given, deactivate the constant bc option
     if (number_space_time_function > 0) then
        kfl_conbc_nsa = 0
     end if


  end if


!
! Formats
!
100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)

  
end subroutine nsa_reabcs
