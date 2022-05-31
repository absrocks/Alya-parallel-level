!-----------------------------------------------------------------------
!> @addtogroup SolidzInput
!> @{
!> @file    sld_reabcs.f90
!> @author  Mariano Vazquez
!> @date    November, 2017-Adds smooth step
!> @brief   Read boundary conditions
!> @details Read boundary conditions
!> @}
!-----------------------------------------------------------------------

subroutine sld_reabcs()

  use def_kintyp, only : ip, rp
  use def_inpout, only : words, exists, param, getint, getrea, getcha
  use def_master, only : INOTSLAVE
  use def_domain, only : kfl_icodb, kfl_icodn
  use def_domain, only : ndime, tbcod, tncod
  use mod_opebcs, only : opnbcs, opbbcs
  use mod_opebcs, only : boundary_conditions_read_node_codes, boundary_conditions_read_boundary_codes
  use mod_ecoute, only : ecoute
  use def_solidz, only : tbcod_sld, tncod_sld
  use def_solidz, only : kfl_funty_sld, nfunc_sld, mtloa_sld, tload_sld
  use def_solidz, only : fubcs_sld
  use def_solidz, only : kfl_local_sld, kfl_csysl_sld
  use def_solidz, only : csysl_sld, SLD_CSYS_CARTESIAN, SLD_CSYS_CYLINDRICAL
  use def_solidz, only : SLD_CSYS_SPHERICAL, SLD_CSYS_EXNOR
  use def_solidz, only : fubcs_sld, rtico_sld
  use def_solidz, only : kfl_conta_sld
  use def_solidz, only : contactbou_sld, neumann_relax, coupling_contact_tol, coupling_contact_its
  use def_solidz, only : vect_proje_sld
  use def_solidz, only : kfl_bodyf_sld, kfl_conbc_sld
  use def_solidz, only : ncoef_sld
  use def_solidz, only : ncrak_sld, crkco_sld
  use def_solidz, only : kfl_conta_stent, conbo_sld, r_fin_stent, contact_friction_sld


  implicit none

  integer(ip)   :: itloa,ifunc,nauxi,icrak,idime,dummi,ipoin,icavi,imate,nposi
  real(rp)      :: rauxi(2)
  character(5)  :: wcsys

  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opnbcs(0_ip,1_ip,ndime,0_ip,tncod_sld)
  end if
  if( kfl_icodb > 0 ) then
     call opbbcs(0_ip,1_ip,ndime,tbcod_sld)
  end if

  if( INOTSLAVE ) then
     !
     ! Initializations global variables
     !
     kfl_conbc_sld  = 0                    ! Constant boundary conditions
     kfl_local_sld  = 0                    ! No local system of reference
     kfl_conta_sld  = 0                    ! PDN contact
     kfl_bodyf_sld  = 0                    ! Number of body forces (by default is OFF)
     kfl_conta_stent = 0                   ! Flag for stenting: 1 is crimpring, 2 i expansion, 3 is charge
     r_fin_stent = 0                       ! Initialization of final raidus of stent
     nfunc_sld      = 0                    ! Number time-dependent boundary condition types
     mtloa_sld(:)   = 0                    ! Number of discrete times per function
     contactbou_sld      = 0                    ! Identificator of contact boundary
     neumann_relax  = 0.9_rp               ! Neumann relaxation
     coupling_contact_tol = 1e-4_rp        ! Contact coupling tolerance
     coupling_contact_its = 20_ip          ! Contact coupling iterations
     vect_proje_sld  = (/ 0.0_rp,1.0_rp,0.0_rp /)  ! Projection direction
     csysl_sld      = 0.0_rp               ! Parameters local coordinate system for boundary conditions


     !
     ! Reach section BOUNDARY_CONDITIONS
     !
     call ecoute('sld_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('sld_reabcs')
     end do
     !
     ! Read header
     !
     if(exists('CONST')) then
        kfl_conbc_sld = 1
     else if(exists('TRANSI').or.(exists('NONCO'))) then
        kfl_conbc_sld = 0
     end if
     !
     ! Allocate the prescription time function vector of vectors
     !
     call sld_membcs(3_ip)
     !
     ! Read data
     !
     call ecoute('sld_reabcs')
     !
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> $ Boundary conditions
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> BOUNDARY_CONDITIONS, TRANSIENT | CONSTANT
     !
     do while (words(1) /= 'ENDBO')

        if (words(1) == 'PARAM' ) then

           !-------------------------------------------------------------
           !
           ! Parameters
           !
           !-------------------------------------------------------------
           !
           ! ADOC[1]> PARAMETERS
           ! ADOC[d]> <ul> </ul>
           ! ADOC[d]> Parameters:
           !
           call ecoute('sld_reabcs')
           do while(words(1)/='ENDPA')

              if ( words(1) == 'COORD' ) then
                 !
                 ! Coordinate sytem parameters
                 !
                 ! ADOC[2]> COORDINATE_SYSTEM: BASIS=CARTESIAN | CYLINDRICAL | EXNOR   $ Local coordinate system
                 ! ADOC[d]> COORDINATE_SYSTEM:
                 ! ADOC[d]> Local coordinate system for boundary conditions.
                 ! ADOC[d]> <ul>
                 ! ADOC[d]> <li> Set BASIS=CARTESIAN (default) <tt>rea1,rea2,rea3,rea4,rea5,rea6,rea7,rea8,rea9</tt>. It defines a cartesian
                 ! ADOC[d]> a coordinate system by three points. The first point <tt>c(rea1,rea2,rea3)</tt> corresponds to the center of
                 ! ADOC[d]> the coordinate system;
                 ! ADOC[d]> The second point  <tt>a(rea4,rea5,rea6)</tt> has to form the primary axis P; The last point
                 ! ADOC[d]> has to form the secondary axis S <tt>b(rea7,rea8,rea9)</tt>. The normal axis N is calculated internally doing the vectorial
                 ! ADOC[d]> product of N = P x S.</li>
                 ! ADOC[d]> <li> Set BASIS=CYLINDRICAL to define a cylindrical coordinate system <tt>rea1,rea2,rea3,rea4,rea5,rea6,rea7,rea8,rea9</tt>.
                 ! ADOC[d]> It defines a cylindrical coordinate sytem by three points. The first point <tt>c(rea1,rea2,rea3)</tt> corresponds to the
                 ! ADOC[d]> the center or the coordinate system; The second point  <tt>a(rea4,rea5,rea6)</tt> has to form the radial axis R;
                 ! ADOC[d]> The last point <tt>b(rea7,rea8,rea9)</tt> has to form the axial axis A; The tangencial axis T is calculated
                 ! ADOC[d]> internally doing the vectorial product of T = R x A.</li>
                 ! ADOC[d]> <li> Set BASIS=EXNOR to use Alya coordinate system based on exteriorn normal. No points has to be defined.</li>
                 ! ADOC[d]> </ul>
                 !
                 ! Type of coordinate system
                 !
                 if ( words(2) == 'BASIS' ) then
                    nposi = 2_ip
                    kfl_local_sld = 1_ip
                    wcsys = getcha('BASIS','     ','#Local coordinate system')
                    if (      wcsys == 'CARTE' ) then
                       nposi = nposi + 1_ip
                       kfl_csysl_sld = SLD_CSYS_CARTESIAN
                    else if ( wcsys == 'CYLIN' ) then
                       nposi = nposi + 1_ip
                       kfl_csysl_sld = SLD_CSYS_CYLINDRICAL
                    else if ( wcsys == 'SPHER' ) then
                       nposi = nposi + 1_ip
                       kfl_csysl_sld = SLD_CSYS_SPHERICAL
                    else if ( wcsys == 'EXNOR' ) then
                       kfl_csysl_sld = SLD_CSYS_EXNOR
                    end if
                    !
                    ! Read parameters
                    if ( wcsys /= 'EXNOR' ) csysl_sld(1:ndime*3) = param(nposi:ndime*3+nposi-1)
                 end if

              else if ( words(1) == 'CONTA' ) then

                 !----------------------------------------------------------
                 !
                 ! Contact parameters
                 !
                 !----------------------------------------------------------
                 !
                 ! ADOC[2]> CONTACT_PDN                      $ PDN Contact
                 ! ADOC[3]>   UNILATERAL | BILATERAL | RBO_DEFORMABLE
                 ! ADOC[3]>   BOUNDARY=   int1
                 ! ADOC[3]>   PROJECTION: X | Y | Z
                 ! ADOC[2]> END_CONTACT_PDN
                 ! ADOC[d]> CONTACT_PDN: Partial Dirichlet-Neumann contact algorithm (Rivero et. al 2018)
                 !
                 kfl_local_sld = 1_ip ! local axes
                 !
                 call ecoute('sld_reabcs')
                 do while( words(1) /= 'ENDCO' )

                    if (words(1) == 'CODES') then
                       conbo_sld(1) = int(param(1))  !EXTERIOR SURFACE
                       conbo_sld(2) = int(param(2))  !INTERIOR SURFACE

                    else if ( words(1) == 'STENT' ) then
                       if ( words(2) == 'CRIMP' ) then
                         kfl_conta_stent = 1_ip
                         if ( words(3) == 'RMIN' ) then
                           r_fin_stent = param(3)
                         else if ( words(3) == 'RELAX' ) then
                           kfl_conta_stent = -1_ip 
                         end if
                       else if ( words(2) == 'EXPAN' ) then
                         kfl_conta_stent = 2_ip
                         if ( words(3) == 'RMAX' ) then
                           r_fin_stent = param(3)
                         else if ( words(3) == 'RELAX' ) then
                           kfl_conta_stent = -2_ip 
                         end if
                       else if ( words(2) == 'CHARG' ) then
                           kfl_conta_stent = 3_ip
                       else if (words(2) == 'MOVEV' ) then
                           kfl_conta_stent = 4_ip
                       end if

                    else if ( words(1) == 'FRICT' ) then
                       contact_friction_sld = param(1)

                    else if ( words(1) == 'UNILA' ) then
                       kfl_conta_sld = 1

                    else if ( words(1) == 'BILAT' ) then
                       kfl_conta_sld = 2
                       if ( exists('NEUMA') ) neumann_relax        = getrea('NEUMA',0.9_rp,'#Neumann relaxation factor')
                       if ( exists('TOLER') ) coupling_contact_tol = getrea('TOLER',1.0E-4_rp,'#Coupling contact tolerance')
                       if ( exists('ITERA') ) coupling_contact_its = getint('ITERA',20_ip,'#Coupling contact iterations')

                    else if ( words(1) == 'RBODE' ) then
                       kfl_conta_sld = 3

                    else if ( words(1) == 'BOUND' ) then
                       contactbou_sld = int(param(1))

                    else if ( words(1) == 'PROJE' ) then
                       if (      words(2) == 'X    ' ) then
                          vect_proje_sld(1:3) = (/ 1.0_rp,0.0_rp,0.0_rp /)
                       else if ( words(2) == 'Y    ' ) then
                          vect_proje_sld(1:3) = (/ 0.0_rp,1.0_rp,0.0_rp /)
                       else if ( words(2) == 'Z    ' ) then
                          vect_proje_sld(1:3) = (/ 0.0_rp,0.0_rp,1.0_rp /)
                       else
                          vect_proje_sld(1:3) = (/ 0.0_rp,1.0_rp,0.0_rp /)
                       end if

                    end if
                    call ecoute('sld_reabcs')

                 end do

              end if
              call ecoute('sld_reabcs')

           end do
           ! ADOC[1]> END_PARAMETERS

        else if (words(1) == 'CODES' .and. exists('NODES') ) then

           !----------------------------------------------------------
           !
           ! User-defined codes on nodes
           !
           !----------------------------------------------------------
           !
           ! ADOC[d]> <ul> </ul>
           ! ADOC[d]> Codes:
           !
           ! ADOC[1]> CODES, NODES
           ! ADOC[2]> int1 int2 rea1 rea2 rea3 int3    $ int1=code, int2=fixity codes, rea1-3=DoFs int3=function num.
           ! ADOC[2]> int1 int2 rea1 rea2 rea3 int3 AXES=LOCAL
           ! ADOC[2]> int1 int2 rea1 rea2 rea3      SPACE_&_TIME_FUNCTIONS=func_name
           ! ADOC[2]> ...
           ! ADOC[1]> END_CODES
           ! ADOC[d]> <b>CODES, NODES</b>:
           ! ADOC[d]> Interpret the code on nodes to impose displacement for each degrees of freedom.
           ! ADOC[d]> <ul>
           ! ADOC[d]> <li>
           ! ADOC[d]>      <tt>int1</tt> is the node code, e.g. 1, 2, 3 etc.
           ! ADOC[d]>      If the node has multiple codes (e.g. if the code was extrapolated from boundary codes), the
           ! ADOC[d]>      Syntaxis is <tt>int1 = 1 & 2 & 4</tt>. It means that nodes with the three codes 1,2 and 4 are considered.
           ! ADOC[d]> </li>
           ! ADOC[d]> <li>
           ! ADOC[d]>      <tt>int2</tt> has the form f1f2f3 (11, 10, 01 in 2D and 101, 111, etc in 3D). 0 means free or initial condition, 1 means prescribed.
           ! ADOC[d]>      Values f1, f2 and f3 are fixity codes of the displacement for each degrees of freedom.
           ! ADOC[d]>      f1 refers to the first displacement component, f2 to the second and f3 to the third.
           ! ADOC[d]> </li>
           ! ADOC[d]> <li>
           ! ADOC[d]>      <tt>rea1</tt>, <tt>rea2</tt> and <tt>rea3</tt> are the corresponding values for each degree of freedom (first to third displacement component).
           ! ADOC[d]> </li>
           ! ADOC[d]> <li>
           ! ADOC[d]>      <tt>int3</tt> is the time function f(t) to applied to the nodes with code <tt>int1</tt>. The value of the nodal condition will be
           ! ADOC[d]>      <tt>f(t)*(real1,real2,real3)</tt>. The function has to be defined in sample.sld.dat file as FUNCTIONS section.
           ! ADOC[d]> </li>
           ! ADOC[d]> <li>
           ! ADOC[d]>      Set AXES=LOCAL to use local boundary conditions. In this case a COORDINATE_SYSTEM must be defined in the PARAMETERS section.
           ! ADOC[d]> </li>
           ! ADOC[d]> <li>
           ! ADOC[d]>      Set SPACE_&_TIME_FUNCTIONS=func_name to use a space and time function defined in sample.ker.dat file.
           ! ADOC[d]> </li>
           ! ADOC[d]> </ul>
           ! ADOC[d]> Examples:
           ! ADOC[d]> <ul>
           ! ADOC[d]> <li>
           ! ADOC[d]>  3 101 0.0 0.0 1.0 1 => Displacement is prescribed in z-direction, free in y-direction and constrained in x-direction.
           ! ADOC[d]>  Moreover, displacement is applied accoring to the space/time function coded 1.
           ! ADOC[d]> </li>
           ! ADOC[d]> <li>
           ! ADOC[d]>  4 & 6 111 0.0 0.0 0.0 => Prescribe value function number 2 as a Dirichlet boundary condition on all nodes.
           ! ADOC[d]> </li>
           ! ADOC[d]> <li>
           ! ADOC[d]>  4 & 6 111 -1.0 0.0 0.0 => Constant displacement in x-direction and constrained in y- and z-direction at nodes
           ! ADOC[d]>  from the intersection between boundaries coded 4 and 6. No space/time function (constant).
           ! ADOC[d]> </li>
           ! ADOC[d]> </ul>
           tncod => tncod_sld(:)
           call boundary_conditions_read_node_codes('DISPLACEMENT')

        else if( words(1) == 'CODES' .and. exists('BOUND') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on boundaries
           !
           !-------------------------------------------------------------
           !
           ! ADOC[1]> CODES, BOUNDARIES
           ! ADOC[2]> int1 int2 int3 int4 int5 int6    $ int1=code, int2=3, int3-5=DoFs, int6=function number
           ! ADOC[2]> ...
           ! ADOC[1]> END_CODES
           ! ADOC[d]> <b>CODES, BOUNDARIES</b>:
           ! ADOC[d]> Impose a natural boundary condition of type int2 on boundaries with code int1. The different boundary conditions available are
           ! ADOC[d]> <ul>
           ! ADOC[d]> <li>
           ! ADOC[d]>  int2 = 3 (Pressure has to be normal to the surface)
           ! ADOC[d]> </li>
           ! ADOC[d]> </ul>
           ! ADOC[d]> Examples:
           ! ADOC[d]> <ul>
           ! ADOC[d]> <li>
           ! ADOC[d]>  6 3 1 => Impose pressure using function number 1 on boundaries with code 6.
           ! ADOC[d]> </li>
           ! ADOC[d]> </ul>
           tbcod => tbcod_sld(1:)
           call boundary_conditions_read_boundary_codes('DISPLACEMENT')


        else if ( words(1) == 'CRACK' ) then

           !----------------------------------------------------------
           !
           ! Cracks
           !
           !----------------------------------------------------------

           if( exists('NUMBE') ) then
              ncrak_sld = getint('NUMBE',1_ip,'#Number of cracks')
           else
              call runend('SLD_RABCS: NUMBER OF CRACKS MISSING')
           end if
           call ecoute('sld_reabcs')
           call sld_membcs(31_ip)
           if( ndime == 2 ) then
              dummi = 2
           else
              dummi = 4
           end if
           do while( words(1) /= 'ENDCR' )
              if( words(1) == 'CRACK') then
                 icrak = getint('CRACK',1_ip,'#CRACK NUMBER')
                 if( icrak > ncrak_sld .or. icrak < 1 ) call runend('SLD_REABCS: WRONG NUMBER OF CRACKS')
                 call ecoute('sld_reabcs')
                 ipoin = 0
                 do while( words(1) /= 'ENDCR' )
                    ipoin = ipoin + 1
                    if( ipoin > dummi ) call runend('SLD_REABCS: WRONG CRACK DEFINITION')
                    do idime = 1,ndime
                       crkco_sld(idime,ipoin,icrak) = param(idime)
                    end do
                    call ecoute('sld_reabcs')
                 end do
              end if
              call ecoute('sld_reabcs')
           end do

        else if ( words(1) == 'BODYF' ) then

           !----------------------------------------------------------
           !
           ! Body forces
           !
           !----------------------------------------------------------

           kfl_bodyf_sld = int(param(1))

        else if ( words(1) == 'FUNCT' ) then

           !----------------------------------------------------------
           !
           ! Functions
           !
           !----------------------------------------------------------
           !
           ! ADOC[d]> <ul> </ul>
           ! ADOC[d]> Functions:
           !
           ! ADOC[1]> FUNCTIONS                            $ Functions (transient only)
           ! ADOC[2]> TOTAL_NUMBER: int1               $ Total number of functions
           ! ADOC[2]> CONDITION
           ! ADOC[3]> FUNCTION_NUMBER: int2        $ Function number
           ! ADOC[3]> TIME_SHAPE: DISCRETE, LINEAR | SMOOTH  $ Function type
           ! ADOC[3]> SHAPE_DEFINITION
           ! ADOC[4]> int3                     $ Number of lines
           ! ADOC[4]> rea1 rea2 rea3 rea4      $ rea1=time, rea2-4=DoFs
           ! ADOC[4]> ...
           ! ADOC[3]> END_SHAPE_DEFINITION
           ! ADOC[2]> END_CONDITION
           ! ADOC[1]> END_FUNCTIONS
           ! ADOC[d]> <b>FUNCTIONS</b>:
           ! ADOC[d]> Space/time functions applied to those nodes/boundaries.
           ! ADOC[d]> <ul>
           ! ADOC[d]> <li>
           ! ADOC[d]>  DISCRETE, LINEAR
           ! ADOC[d]>  Linear step (Ramp) is applied between discrete times.
           ! ADOC[d]> </li>
           ! ADOC[d]> <li>
           ! ADOC[d]>  DISCRETE, SMOOTH
           ! ADOC[d]>  This function is such that first and second derivatives are zero between two discrete times.
           ! ADOC[d]>  This definition is intended to ramp up or down smoothly from one amplitude value to another.
           ! ADOC[d]> </li>
           ! ADOC[d]> </ul>

           call ecoute('sld_reabcs')

           ! Initializations
           nauxi = 0
           ifunc = 0

           ! Total number of functions
           if (words(1) == 'TOTAL') nfunc_sld =  int(param(1))
           if (nfunc_sld == 0) then
              call runend('SLD_REABCS: PROVIDE A TOTAL_NUMBER OF FUNCTIONS')
           else if (nfunc_sld > 9_ip) then
              call runend('SLD_REABCS: THE TOTAL_NUMBER OF FUNCTIONS MUST BE LOWER THAN 10')
           end if

           do while(words(1) /= 'ENDFU')
              if(kfl_conbc_sld == 0) then                   ! non-constant (transient) BC

                 ! Condition
                 do while(words(1) /= 'ENDCO')
                    if (words(1) == 'FUNCT') then
                       ! Function number
                       ifunc = getint('FUNCT',1_ip,'#FUNCTION NUMBER')
                       kfl_funty_sld(8,ifunc) = 1_ip        ! default: only apply function to fixity codes 1
                       if (ifunc < 0 .or. ifunc > 10_ip) then
                          call runend('SLD_REABCS: WRONG FUNCION NUMBER, MUST BE GT.0 AND LT.10')
                       else
                          nauxi = nauxi + 1
                          if (nauxi > nfunc_sld) call runend('SLD_REABCS: MORE FUNCTIONS THAN THE TOTAL_NUMBER ARE GIVEN')
                       end if

                       ! Function type
                    else if ( words(1)=='TIMES' ) then      ! time shape

                       if ( words(2)=='PARAB' .or. words(2)=='LINEA' .or. words(2)=='POLYN' ) then
                          kfl_funty_sld(1,ifunc)=1
                          if (words(2) == 'LINEA') kfl_funty_sld(1,ifunc)=10
                          kfl_funty_sld(2,ifunc)=20
                          kfl_funty_sld(3:7,ifunc)=1        ! default: apply to all equations

                       else if ( words(2) == 'PERIO' ) then
                          kfl_funty_sld(1,ifunc)=2
                          kfl_funty_sld(2,ifunc)=20
                          kfl_funty_sld(3:7,ifunc)=1        ! default: apply to all equations

                       else if ( words(2) == 'DISCR' ) then
                          kfl_funty_sld(1,ifunc)=3_ip
                          if (words(3) == 'SMOOT') then
                             kfl_funty_sld(2,ifunc)= 1_ip
                          else
                             kfl_funty_sld(2,ifunc)= 0_ip
                          end if

                          ! Shape defintion
                          call ecoute('sld_reabcs')
                          rauxi(:) = 0.0_rp
                          if ( words(1) == 'SHAPE' ) then   ! defining the shape by discrete points
                             ! Reference time and value
                             rauxi(1)= getrea('REFER', 0.0_rp, 'Reference value, sometimes useful')
                             rauxi(2)= getrea('TSTAR', 0.0_rp, 'Startint time, sometimes useful')

                             ! Save total number of tabular lines
                             call ecoute('sld_reabcs')
                             mtloa_sld(ifunc) = int(param(1))

                             ! Allocate the prescription time function vector for ifunc
                             call sld_membcs(10_ip + ifunc)

                             ! Save discrete times and values
                             call ecoute('sld_reabcs')
                             itloa = 0
                             do while(words(1)/='ENDSH')
                                itloa= itloa + 1
                                tload_sld(ifunc)%a(ndime+1,itloa)= param(1) - rauxi(2)         ! time
                                tload_sld(ifunc)%a(1:ndime,itloa)= param(2:ndime+1) - rauxi(1) ! prescribed value
                                call ecoute('sld_reabcs')
                             end do
                          end if

                          kfl_funty_sld(3:7,ifunc)=1       ! default: apply to all equations

                       else
                          kfl_funty_sld(1,ifunc)=0

                       end if

                       !
                       ! OLD WAY: funcrion without an input file
                       !
                       if(kfl_funty_sld(1,ifunc)>0) then
                          if(kfl_funty_sld(1,ifunc) /= 3) then
                             fubcs_sld(1:10,ifunc)=param(2:11)
                          end if
                       end if

                    else if ( words(1)=='REFVA' ) then    ! reference values

                       fubcs_sld(15,ifunc)=param(1)
                       fubcs_sld(16,ifunc)=param(2)
                       fubcs_sld(17,ifunc)=param(3)
                       fubcs_sld(18,ifunc)=param(4)
                       fubcs_sld(19,ifunc)=param(5)

                    else if ( words(1)=='TIMEL' ) then    ! time lapse

                       fubcs_sld(11,ifunc)   =param(1)  !   start
                       fubcs_sld(12,ifunc)   =param(2)  !   end
                       fubcs_sld(13,ifunc)   =1.0e10_rp !   when repeat (default value, very high)
                       rtico_sld( 1,ifunc)    =param(1) !   initial start time for all eqs.
                       rtico_sld( 2,ifunc)    =param(2) !   initial final time for all eqs.
                       if (kfl_funty_sld(1,ifunc)==10) then   ! LINEAR: divide the load interval in equispaced steps
                          kfl_funty_sld(1,ifunc)=1
                          fubcs_sld(3:10,ifunc)= 0.0_rp
                          fubcs_sld(2,ifunc)   = 1.0_rp/(fubcs_sld(12,ifunc)-fubcs_sld(11,ifunc))
                       end if

                    else if ( words(1)=='TIMER' ) then    ! time repeat

                       fubcs_sld(13,ifunc)=param(1)     !   when repeat

                    else if ( words(1)=='FIXIT' ) then    ! special behavior controlled by the fixity code

                       if (exists('CONST')) then
                          kfl_funty_sld(8,ifunc) = 1    ! constrained in the direction with fixity code other than 1
                          call runend('SLD_REABCS: CONSTRAINED IS A DEPRECATED OPTION')
                       else if (exists('UNCON')) then   ! free in the direction with fixity code other than 1
                          kfl_funty_sld(8,ifunc) = 2
                          call runend('SLD_REABCS: UNCONSTRAINED IS A DEPRECATED OPTION')
                       else if (exists('ALLDI')) then   ! apply function to all dimensions
                          kfl_funty_sld(8,ifunc) = 3
                          call runend('SLD_REABCS: ALLDIMENSIONS IS A DEPRECATED OPTION')
                       end if

                    end if

                    call ecoute('sld_reabcs')
                 end do
                 ! End condition

              end if

              call ecoute('sld_reabcs')

           end do
           ! End functions

        end if

        call ecoute('sld_reabcs')

     end do
     ! ADOC[0]> END_BOUNDARY_CONDITIONS

     ! give a default size to non-defined tloads
     do ifunc= 1,10
        if (mtloa_sld(ifunc) == 0) then
           mtloa_sld(ifunc) = 1            ! done to allocate a default memory space
           call sld_membcs(10_ip + ifunc)  ! allocate the prescription time function vector for ifunc
        end if
     end do

  end if

end subroutine sld_reabcs
