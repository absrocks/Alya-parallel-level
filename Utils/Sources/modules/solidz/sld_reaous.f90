!-----------------------------------------------------------------------
!> @addtogroup SolidzInput
!> @{
!> @file    sld_reaous.f90
!> @author  Solidz Team
!> @date    August, 2006
!>          - Subroutine written
!> @note    <GGU> Revisar estas lineas
!> @brief   Read post-process Solidz data
!>
!> @details Read post-process data
!>           - Array post-process (Field outputs)
!>           - Element, boundary and node sets
!>           - Witness points
!>           - Error (exact solution)
!>           - Matrix output
!>           - Post-process using a local coordinate system
!>           - Operations on reaction forces and displacement
!>           - Rotation
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_reaous()

  use def_kintyp,  only : ip, rp
  use def_inpout,  only : words, param, getint, getrea, getcha
  use def_master,  only : INOTSLAVE, mem_modul, modul
  use def_domain,  only : ndime
  use mod_memory,  only : memory_alloca
  use mod_ecoute,  only : ecoute
  use def_solidz,  only : kfl_exacs_sld, kfl_foten_sld, kfl_rotei_sld
  use def_solidz,  only : kfl_rsfor_sld, kfl_psmat_sld, psmat_sld
  use def_solidz,  only : kfl_csysp_sld
  use def_solidz,  only : csysp_sld, SLD_CSYS_CARTESIAN, SLD_CSYS_CYLINDRICAL
  use def_solidz,  only : idilimi_sld, resnode_sld, numnodfor_sld, reslimi_sld

  implicit none

  integer(ip)          :: dummi,ipoin,nposi
  integer(ip)          :: icsys
  character(5)         :: wcsys

  if( INOTSLAVE ) then
     !
     ! Initializations global variables
     !
     kfl_exacs_sld = 0                  ! Exact solution (FEM errors)
     kfl_foten_sld = 1                  ! Tensors of forces caust, green, lepsi (0 if none)
     kfl_rotei_sld = 0                  ! Rotate and correct sigma eigenvalues (ONLY FOR 3D PRISMATIC SHELLS)
     kfl_rsfor_sld = 0                  ! Do not compute reaction forces and displacement
     kfl_psmat_sld = 0                  ! Postscript file of the matrix
     kfl_csysp_sld = SLD_CSYS_CARTESIAN ! Type of local coordinate system for post-process
     csysp_sld(:,:) = 0.0_rp            ! Parameters local coordinate system for post-process

     !
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> $ Output and postprocess
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> OUTPUT_&_POST_PROCESS
     ! ADOC[1]> START_POSTPROCESS_AT: STEP=int | TIME=int   $ Start from step number or time
     ! ADOC[d]> START_POSTPROCESS_AT:
     ! ADOC[d]> This option is used to start the post-process at a certain step or time in the simulation.
     ! ADOC[d]> Set TIME=<tt>int</tt> to specify the initial time for post-process or STEP=<tt>int</tt> to specify the first time step for post-process.
     !
     ! ADOC[1]> POSTPROCESS XXXXX, STEPS=int                $ Variable for post-process
     ! ADOC[1]> ...
     ! ADOC[d]> POSTPROCESS:
     ! ADOC[d]> This option is to post-process a variable XXXXX for visualitzation. Set STEPS=<tt>int</tt> equal to the post-process frequency.
     ! ADOC[d]> Results are written in <tt>caseName-XXXXX-YYYYYYYY.post.alyabin</tt> file.
     !
     ! ADOC[1]> ELEMENT_SET
     ! ADOC[1]>  XXXXX                                      $ Variables computed in volume integrals
     ! ADOC[1]>  ...
     ! ADOC[1]> END_ELEMENT_SET
     ! ADOC[d]> ELEMENT_SET:
     ! ADOC[d]> Variables computed as volume integrals on element sets. The results are written in
     ! ADOC[d]> <tt>caseName-element.sld.set</tt> file.
     ! ADOC[d]> The following variables can be post-processed
     ! ADOC[d]> <ul>
     ! ADOC[d]> <li>LEPRT: Averaged logarithminc strain using polar-cylindrical coordiantes (eps_rr + eps_tt)</li>
     ! ADOC[d]> </ul>
     !
     ! ADOC[1]> BOUNDARY_SET
     ! ADOC[1]>  XXXXX                                      $ Variables computed in boundary integrals
     ! ADOC[1]>  ...
     ! ADOC[1]> END_BOUNDARY_SET
     ! ADOC[d]> BOUNDARY_SET:
     ! ADOC[d]> Variables computed as boundary integrals on boundary sets. The results are written in
     ! ADOC[d]> <tt>caseName-boundary.sld.set</tt> file.
     ! ADOC[d]> The following variables can be post-processed
     ! ADOC[d]> <ul>
     ! ADOC[d]> <li>DIBOX, DIBOY, DIBOZ, DIBOU (Magnitude): Averaged displacement</li>
     ! ADOC[d]> <li>FRBOX, FRBOY, FRBOZ, FRBOU (Magnitude): Sum of reaction forces</li>
     ! ADOC[d]> </ul>
     !
     ! ADOC[1]> NODE_SET
     ! ADOC[1]>  XXXXX                                      $ Variables computed on nodes
     ! ADOC[1]>  ...
     ! ADOC[1]> END_NODE_SET
     ! ADOC[d]> NODE_SET:
     ! ADOC[d]> Variables on node sets. The results are written
     ! ADOC[d]> in <tt>caseName-node.sld.set</tt> file.
     ! ADOC[d]> The following variables can be post-processed
     ! ADOC[d]> <ul>
     ! ADOC[d]> <li>DISPX, DISPY, DISPZ: Displacements</li>
     ! ADOC[d]> <li>VELOX, VELOY, VELOZ: Velocities</li>
     ! ADOC[d]> <li>ACCEX, ACCEY, ACCEZ: Accelerations</li>
     ! ADOC[d]> <li>FRXIX, FRXIY, FRXIZ: Reaction Forces</li>
     ! ADOC[d]> <li>SIGXX, SIGYY, SIGZZ, SIGYZ, SIGXZ, SIGXY: Global stresses</li>
     ! ADOC[d]> <li>EPSXX, EPSYY, EPSZZ, EPSYZ, EPSXZ, EPSXY: Global strains</li>
     ! ADOC[d]> <li>LEPXX, LEPYY, LEPZZ, LEPYZ, LEPXZ, LEPXY: Global logarithmic strains</li>
     ! ADOC[d]> <li>COORX,COORY,COORZ: Coordinates</li>
     ! ADOC[d]> <li>SEQVM: Von Mises Stress</li>
     ! ADOC[d]> </ul>
     ! ADOC[d]>
     ! ADOC[1]> WITNESS_POINTS
     ! ADOC[1]>  XXXXX                                      $ Variables computed witness points
     ! ADOC[1]>  ...
     ! ADOC[1]> END_WITNESS_POINTS
     ! ADOC[d]> WITNESS_POINTS:
     ! ADOC[d]> Variables on witness points. The results are stored
     ! ADOC[d]> <tt>caseName.sld.wit</tt>.
     ! ADOC[d]> The following variables can be postprocessed
     ! ADOC[d]> <ul>
     ! ADOC[d]> <li>DISPX, DISPY, DISPZ: Displacements</li>
     ! ADOC[d]> <li>VELOX, VELOY, VELOZ: Velocities</li>
     ! ADOC[d]> <li>ACCEX, ACCEY, ACCEZ: Accelerations</li>
     ! ADOC[d]> <li>FRXIX, FRXIY, FRXIZ: Reaction Forces</li>
     ! ADOC[d]> <li>SIGXX, SIGYY, SIGZZ, SIGYZ, SIGXZ, SIGXY: Global stresses</li>
     ! ADOC[d]> <li>EPSXX, EPSYY, EPSZZ, EPSYZ, EPSXZ, EPSXY: Global strains</li>
     ! ADOC[d]> <li>LEPXX, LEPYY, LEPZZ, LEPYZ, LEPXZ, LEPXY: Global Logarithmic strains</li>
     ! ADOC[d]> <li>COORX,COORY,COORZ: Coordinates</li>
     ! ADOC[d]> <li>SEQVM: Von Mises Stress</li>
     ! ADOC[d]> </ul>
     !
     !
     ! Reach the section
     !
     call ecoute('reaous')
     do while(words(1)/='OUTPU')
        call ecoute('reaous')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDOU')
        call ecoute('reaous')
        call posdef(2_ip,dummi)

        if(words(1)=='OUTPU') then

           !----------------------------------------------------------
           !
           ! Output (Exact solution)
           !
           !----------------------------------------------------------

           if(words(2)=='ERROR') then

              !----------------------------------------------------------
              !
              ! ERROR
              !
              !----------------------------------------------------------
              !
              ! ADOC[1]> OUTPUT, ERROR, SOLUTION = int               $ Manufactured solution
              ! ADOC[d]> OUTPUT, ERROR:
              ! ADOC[d]> Manufactured solutions used to carry out code typing checking and mesh convergence.
              ! ADOC[d]> The L2, L1 and Linf norm of displacement
              ! ADOC[d]> gradients can be found in file <tt>caseName.sld.log</tt>. Dirichlet boundary
              ! ADOC[d]> conditions are automatically imposed on all the boundaries.
              ! ADOC[d]> The manufactured solutions available are the following
              ! ADOC[d]> <ul>
              ! ADOC[d]> <li> SOLUTION = 2. Exact solution using Solidz Implicit integration scheme. </li>
              ! ADOC[d]> </ul>
              !
              kfl_exacs_sld = getint('SOLUT',1_ip,'#Exact solution')

           else if (words(2) == 'MATRI') then

              !----------------------------------------------------------
              !
              ! MATRIX
              !
              !----------------------------------------------------------
              !
              ! ADOC[1]> OUTPUT, MATRIX, STEP = int1, ITERA = int2   $ Matrix output
              ! ADOC[d]> OUTPUT, MATRIX:
              ! ADOC[d]> Postprocess global system matrix. Set STEP=<tt>int1</tt> to choose the time step and
              ! ADOC[d]> ITERA=<tt>int2</tt> to choose the N-R iteration. The available formats for matrix output are the following
              ! ADOC[d]> <ul>
              ! ADOC[d]> <li> Postcript file (.ps) </li>
              ! ADOC[d]> <li> Data file (.dat) </li>
              ! ADOC[d]> </ul>
              !
              kfl_psmat_sld = 1_ip
              psmat_sld(1) = getrea('STEP ',1.0_rp,'#Step to postprocess matrix')
              psmat_sld(2) = getrea('ITERA',1.0_rp,'#Iteration to postprocess matrix')

           end if

        else if (words(1) == 'ROTAT') then

           !----------------------------------------------------------
           !
           ! Rotation
           !
           !----------------------------------------------------------
           !
           ! ADOC[1]> ROTATION: ON | OFF                          $ Rotation
           ! ADOC[d]> ROTATION:
           ! ADOC[d]> To fill
           !
           kfl_rotei_sld = 1

        else if (words(1) == 'PARAM') then

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
                 ! ADOC[2]> COORDINATE_SYSTEM: BASIS=CARTESIAN | CYLINDRICAL         $ Local coordinate for post-process
                 ! ADOC[d]> COORDINATE_SYSTEM:
                 ! ADOC[d]> Local coordinate system for post-process purposes.
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
                 ! ADOC[d]> The last point <tt>b(rea7,rea8,rea9)</tt> has to form the axial axis; The tangencial axis T is calculated
                 ! ADOC[d]> internally doing the vectorial product of T = R x A.</li>
                 ! ADOC[d]> </ul>
                 !
                 ! Type of coordinate system
                 !
                 icsys = getint('COORD',1_ip,'#Code for coordinate system')
                 if ( words(2) == 'BASIS' ) then
                    nposi = 3_ip
                    wcsys = getcha('BASIS','     ','#Type of the coordinate system')
                    if (      wcsys == 'CARTE' ) then
                       kfl_csysp_sld = SLD_CSYS_CARTESIAN
                    else if ( wcsys == 'CYLIN' ) then
                       kfl_csysp_sld = SLD_CSYS_CYLINDRICAL
                    end if
                    !
                    ! Read parameters
                    csysp_sld(1:ndime*3,icsys) = param(nposi:ndime*3+nposi-1)
                 end if

              end if
              call ecoute('sld_reabcs')

           end do
           ! ADOC[1]> END_PARAMETERS

        else if (words(1) == 'RESID') then

           !----------------------------------------------------------
           !
           ! Operate on variables
           !
           !----------------------------------------------------------
           !
           ! ADOC[1]> RESIDUAL_FORCE:                             $ Operate on variables
           ! ADOC[d]> RESIDUAL_FORCE:
           ! ADOC[d]> To fill
           !
           kfl_rsfor_sld = 1
           idilimi_sld = 1
           if (words(2) == 'INTER') then
              kfl_rsfor_sld = 2
              reslimi_sld(1) = getrea('COORX',0.0_rp,'#Position')
              reslimi_sld(2) = getrea('COORY',0.0_rp,'#Position')
              reslimi_sld(3) = getrea('COORZ',0.0_rp,'#Position')
              idilimi_sld(1) = getint('DIMEN',1_ip,'#Dimension')

           else if (words(2) == 'SETOF') then
              ! <GGU> Revisar este caso
              kfl_rsfor_sld = 3

              call ecoute('sld_reaous')
              if (words(1)=='ONNOD') then
                 ipoin = 1
                 numnodfor_sld = getint('ONNOD', 0_ip,'#numberofnodes')
                 call ecoute('sld_reaous')
                 call memory_alloca(mem_modul(1:2,modul),'RESNODE_SLD','sld_reaous',resnode_sld,numnodfor_sld)
                 do while(words(1) /= 'ENDON')
                    resnode_sld(ipoin) = int(param(1))
                    ipoin = ipoin + 1_ip
                    call ecoute('sld_reaous')
                 end do
              end if

           else
              idilimi_sld(1) = getint('DIMEN',1_ip,'#Dimension')
              reslimi_sld(1) = getrea('POSIT',0.0_rp,'#Position')
           end if
        end if
     end do
     !
     ! ADOC[0]> END_OUTPUT_&_POST_PROCESS
     !
  end if

end subroutine sld_reaous
