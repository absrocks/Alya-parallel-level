!-----------------------------------------------------------------------
!> @addtogroup DomainInput
!> @{
!> @file    reastr.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Read mesh strategy
!> @details Read mesh strategy:       
!>          \verbatim
!>          - LQUAD(NELTY) ... Type of integration rule (1 = closed, 0 = open)
!>          - NGAUS(NELTY) ... Number of domain integration points
!>          - NGROU_DOM ...... Number of groups for deflation based solvers
!>          - SCALE(3) ....... Geometrical scale factor
!>          - TRANS(3) ....... Geometrical translation vector
!>          - KFL_GEOME ...... If geometrical conditions exist
!>          - KFL_CONVX ...... What to do with convex node (if geometrical conditions are applied)
!>          - KFL_ICODN ...... If there are codes on nodes
!>          - KFL_ICODB ...... If the are codes on boundaries
!>          - MCONO .......... Max # codes per nodes
!>          - KFL_EXTRA ...... If boundary codes are extrapolated to nodes
!>          - NPBCS .......... Max # of possibe geometrical conditions
!>          - GEOAN .......... Geometrical angle to decide if we have corners (if geometrical conditions are applied)
!>          - LSBCS(8,100) ... List of boundaries for a given geometrical code (first argument is condition)
!>          \endverbatim
!> @} 
!-----------------------------------------------------------------------
subroutine reastr()
  use def_parame
  use def_master 
  use def_domain
  use def_inpout
  use def_elmtyp
  use mod_memchk
  use mod_iofile
  use mod_ecoute, only : ecoute
  use mod_reabcs, only : reabcs_geometrical
  use mod_elmgeo, only : elmgeo_element_name_to_type
  implicit none
  integer(ip) :: ielty,jelty,imate,jmate,icode


  !if( ISEQUEN .or. ( IMASTER .and. kfl_ptask /= 2 ) ) then
  if( INOTSLAVE ) then

     lquad              =  0          ! Open integration rule
     ngaus              =  0          ! No Gauss points

     kfl_ngrou          =  0          ! Strategy to construct groups
     ngrou_dom          =  0          ! Number of groups (for deflated): -2 (slaves), -1 (automatic, master), >0 (master)
     ngrou_dom_target   =  0          ! Number of groups (for deflated): -2 (slaves), -1 (automatic, master), >0 (master)
     ngrou_boxes_coarse =  1          ! Number of coarse boxes for SFC
     ngrou_boxes_fine   =  128        ! Number of fines boxes for SFC

     xscal(1)           =  1.0_rp     ! X scale factor
     xscal(2)           =  1.0_rp     ! Y scale factor
     xscal(3)           =  1.0_rp     ! Z scale factor

     trans(1)           =  0.0_rp     ! X translation
     trans(2)           =  0.0_rp     ! Y translation 
     trans(3)           =  0.0_rp     ! Z translation

     kfl_extra          =  0          ! Extrapolate from boundary to node

     kfl_geome          =  0          ! Do not compute geometrical normals
     kfl_convx          =  1          ! Use exnor for convex angle nodes
     kfl_frees          =  0          ! Freestream criteria. 0: use wind angle to decide inflow/outflow
     npbcs              =  0          ! Max # of possibe geometrical conditions
     lsbcs              =  0          ! List of boundaries for a given geometrical code
     awind              =  0.0_rp     ! Wind angle for freestream
     tolan              =  5.0_rp*pi/180.0_rp ! Tolerance used to define inflow from freestream
     geoan              = 45.0_rp     ! Angle for geometrical normals

     kfl_chege          =  0          ! Don't check geometry
     kfl_naxis          =  0          ! Cartesian coordinate system
     kfl_spher          =  0          ! Cartesian coordinate system
     kfl_bouel          =  1          ! Element # connected to boundary is known
     kfl_divid          =  0          ! Divide elements into tetra
     curvatureDataField =  0          ! curvature data field id
     curvatureField     =  0          ! curvature data field id

     materials_nlaye    = 0           ! Automatic generation of materials
     materials_icode    = 0           ! Automatic generation of materials
     materials_imate    = 1           ! Automatic generation of materials
     !
     ! Reach the section 
     !
     imate = 1
     jmate = 1
     call ecoute('REASTR')
     do while( words(1) /= 'STRAT' )
        call ecoute('REASTR')
     end do
     !
     !
     !.md<module>kernel
     !.md<input>case.dom.dat
     !.md<pos>1
     !.md<sec>
     !.md<0># Strategy for the mesh related operations
     !.md<>
     !.md<code>
     !.md<0><b>STRATEGY</b>
     do while(words(1) /= 'ENDST')
        call ecoute('REASTR')

        if( words(1) == 'INTEG' ) then
           !
           !.md<1>INTEGRATION_RULE = Open/Close, Open/Close...              $ Integration rule (Open by default)
           !.md<field>INTEGRATION_RULE
           !.md<com>List of intgeration rules for each of the element type declared in DIMENSION field.
           !.md<com>The list should be written in the following ordering:
           !.md<com>    - 1D: BAR02,BAR03,BAR04
           !.md<com>    - 2D: TRI03,TRI06,QUA04,QUA08,QUA09,QUA16
           !.md<com>    - 3D: TET04,TET10,PYR05,PYR14,PEN06,PEN15,PEN18,HEX08,HEX20,HEX27,HEX64
           !.md<com>
           !.md<com>If only one option is given, it is applied to all element types.
           !
           if( words(3) == '' ) then
              if( words(2) == 'CLOSE' ) then
                 lquad = (/(1_ip,ielty=1,nelty)/)
              else
                 lquad = (/(0_ip,ielty=1,nelty)/)
              end if
           else              
              do ielty = 1,min(nelty,size(words,KIND=ip)-1)
                 if( words(ielty+1) == 'CLOSE' ) then
                    lquad(ielty) = 1_ip
                 else
                    lquad(ielty) = 0_ip
                 end if
              end do
           end if

        else if( words(1) == 'DOMAI' ) then 
           !
           !.md<1>DOMAIN_INTEGRATION_POINTS = int, int...                   $ Number of Integration points
           !.md<field>DOMAIN_INTEGRATION_POINTS
           !.md<com>Number of integration points (0 for automatic) for each of the element type declared in DIMENSION field.
           !.md<com>When automatic, the number of integration points COINCIDES with the number of element nodes.
           !.md<com>To explicitly give a value, the list should be written using the following element ordering:
           !.md<com>    - 1D: BAR02,BAR03,BAR04
           !.md<com>    - 2D: TRI03,TRI06,QUA04,QUA08,QUA09,QUA16
           !.md<com>    - 3D: TET04,TET10,PYR05,PYR14,PEN06,PEN15,PEN18,HEX08,HEX20,HEX27,HEX64
           !
           jelty = 0
           do ielty = 1,nelty
              if( lexis(ielty) > 0 ) then
                 jelty = jelty + 1
                 ngaus(ielty) = int(param(jelty))
                 if( ngaus(ielty) == 1 ) lquad(ielty) = 0      ! One gauss point forced to open rule 
              end if
           end do

        else if( words(1) == 'GAUSS' ) then 
           !
           ! Gauss points
           !
           call ecoute('REAST')
           do while( words(1) /= 'ENDGA' )
              call elmgeo_element_name_to_type(words(1),ielty)
              if( ielty < 1 .or. ielty > nelty ) call runend('REASTR: WRONG ELEMENT TYPE')
              ngaus(ielty) = int(param(1))
              call ecoute('REAST')
           end do

        else if( words(1) == 'GROUP' ) then
           !
           !.md<1>GROUPS = int, SEQUENTIAL_FRONTAL | PARALLEL_FRONTAL | PARTITION | SFC | FIELD $ Number of groups (used for deflation)
           !.md<field>GROUPS
           !.md<com>This option is the number of groups to be constructed for deflation based solvers.
           !.md<com>    - GROUPS = -1 ... The number of groups is computed in an automatic way.
           !.md<com>    - GROUPS = -2 ... In parallel, one group per subdomain.
           !.md<com>    - GROUPS >  0 ... Number of groups specified by the user (highly recommanded).
           !.md<com>    - Options:
           !.md<com>        - `SEQUENTIAL_FRONTAL` (default)
           !.md<com>        - `PARALLEL_FRONTAL` (whole parallel execution)
           !.md<com>        - `PARTITION`
           !.md<com>        - `SFC, COARSE=int | FINE=int`
           !.md<com>        - `FIELD`

           ngrou_dom = getint('GROUP',0_ip,'#NUMBER OF GROUPS')
           if(      exists('SEQUE') .or. exists('FRONT' ) ) then
              kfl_ngrou = -1
           else if( exists('PARAL') ) then
              kfl_ngrou = -2
           else if( exists('PARTI') ) then
              kfl_ngrou = -3
           else if( exists('SFC  ') ) then
              kfl_ngrou = -4
              if( exists('COARS') ) ngrou_boxes_coarse = getint('COARS',1_ip  ,'#Number of coarse bins')
              if( exists('FINE ') ) ngrou_boxes_fine   = getint('FINE ',128_ip,'#Number of fine bins')
           else if( exists('FIELD') ) then
              kfl_ngrou = getint('FIELD',1_ip,'#groups from field')
           else
              kfl_ngrou = -1
           end if
           ngrou_dom_target = ngrou_dom
           
        else if(words(1)=='SCALE') then 
           !
           !.md<1>SCALE : XSCALE = real, YSCALE = real, ZSCALE = real       $ Scaling factors for the geometry
           !.md<field>SCALE
           !.md<com>Multiply each coordinates by its corresponding scaling factor.
           !.md<com>In postprocess, the geometry is written using this scaling.
           !
           xscal(1) = getrea('XSCAL',1.0_rp,'#x-factor')
           xscal(2) = getrea('YSCAL',1.0_rp,'#y-factor')
           xscal(3) = getrea('ZSCAL',1.0_rp,'#z-factor')

        else if( words(1) == 'TRANS' ) then 
           !
           !.md<1>TRANSLATION = XTRANS = real, YTRANS = real, ZTRANS = real $ Translation of the geometry
           !.md<field>TRANSLATION
           !.md<com>Translate the geometry.
           !
           trans(1) = getrea('XTRAN',0.0_rp,'#x-translation')
           trans(2) = getrea('YTRAN',0.0_rp,'#y-translation')
           trans(3) = getrea('ZTRAN',0.0_rp,'#z-translation')

        else if( words(1) == 'GEOME' ) then
           !
           ! Read geometrical boundary conditions
           !
           call reabcs_geometrical()

        else if( words(1) == 'BOUND' ) then
           !
           ! Compute LELBO because it is unknown
           !
           if( words(2) == 'KNOWN' ) then
              kfl_bouel = 1
           else if( words(2) == 'UNKNO' ) then
              kfl_bouel = 0
           else if( option('BOUND') ) then
              kfl_bouel = 1
           else
              kfl_bouel = 0
           end if

        else if( words(1) == 'CHECK' ) then
           !
           ! Check geometry
           !
           if( option('CHECK') ) kfl_chege = 1

        else if( words(1) == 'CRVDA' ) then
           !
           ! Reading curvatureDataField
           !
           curvatureDataField = getint('CRVDA',1_ip,'#Curve dat field id')

        else if( words(1) == 'CRVGE' ) then
           !
           ! Reading curvatureField
           !
           curvatureField = getint('CRVGE',1_ip,'#Curved geometry field id')

        else if( words(1) == 'AXISY' ) then
           !
           ! Cylindrical coordinates
           !
           if( option('AXISY') ) kfl_naxis = 1

        else if( words(1) == 'SPHER' ) then
           !
           ! Spherical coordinates
           !
           if( option('SPHER') ) kfl_spher = 1

        else if( words(1) == 'DIVID' ) then
           !
           ! Divide mesh into tetras
           !
           if( words(2) == 'ON   ' .or. words(2) == 'YES  ' ) kfl_divid = 1

        else if( exists('EXTRA') ) then
           !
           ! Extrapolate boundary conditions
           !
           if( option('EXTRA') ) kfl_extra = 1

        else if( words(1) == 'MATER' ) then

           !-------------------------------------------------------------
           !
           ! Materials depending on boundary codes
           !
           !-------------------------------------------------------------
           !
           !.md<1>MATERIALS, BOUNDARIES, DEFAULT= int
           !.md<2>CODE =     int,int...                                                        $ From which codes material should be generated
           !.md<2>MATERIAL = int                                                               $ Material number
           !.md<2>LAYERS =   int                                                               $ Number of element layers
           !.md<1>END_MATERIALS
           !.md<field>MATERIALS, BOUNDARIES
           !.md<com>This option enables to generate automatically a certain number of element layers (LAYERS= int) of one material (MATERIAL= int)
           !.md<com>starting from boundaries with a given code (CODE= int). The option DEFAULT= int assign automatically this material
           !.md<com>to all the element before applying the layers. It can be useful for example in CFD to define a different material at outflows
           !.md<com>on which the viscosity will be higher.
           !
           if( words(2) == 'BOUND' ) then

              icode = 0
              call ecoute('reastr')
              do while(words(1)/='ENDMA')
                 if( words(1)== 'CODE ' ) then
                    jmate = imate + nnpar - 1
                    icode = nnpar
                    if( jmate > size(materials_icode,KIND=ip) ) call runend('REASTR: TOO MANY MATERIALS DEFINED FROM BOUNDARIES')                    
                    materials_icode(imate:jmate) = int(param(1:icode),ip)
                 else
                    if( icode == 0 ) then
                       call runend('REASTR: IN MATERIALS_FROM_BOUDNARIES, CODE SHOULD BE DEFINED FIRST')
                    else if(  words(1)== 'MATER' ) then
                       materials_imate(imate:jmate) = int(param(1),ip)
                    else if( words(1) == 'LAYER' ) then
                       materials_nlaye(imate:jmate) = int(param(1),ip)
                    end if
                 end if
                 call ecoute('reastr')
              end do
              imate = imate + icode

           end if

        end if

     end do
     !.md<0><b>END_STRATEGY</b>
     !.md</code>
     !.md<>

  end if

end subroutine reastr
