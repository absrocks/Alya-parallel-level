subroutine rrudat()
  !------------------------------------------------------------------------
  !****f* master/rrudat
  ! NAME
  !    rrudat
  ! DESCRIPTION
  !    This routine reads run data
  ! OUTPUT
  !   KFL_CUSTO ...... Customer
  !   KFL_PRELI ...... If this is preliminary run
  !   NPRIT     ...... Preliminary frequence
  !   KFL_RSTAR ...... If this is a restart run
  !   KFL_OUTFO ...... Output format
  !   KFL_LATEX ...... If output info in latex file
  !   LUN_LIVEI ...... If live file is screen or file
  !   KFL_MEMOR ...... If output memory
  !   KFL_VARCOUNT ... If output variable memory counter
  !   KFL_LOTME ...... If we do not have a lot of memory
  !   KFL_LOGFI ...... If log file should be written
  !   CPU_LIMIT ...... CPU limit
  !   TITLE     ...... Problem title
  ! USES
  ! USED BY
  !    Reapro
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_inpout
  use def_kermod
  use mod_memory,         only : kfl_memor
  use mod_memory,         only : kfl_varcount
  use mod_coupling_timer
  use mod_ecoute, only :  ecoute

  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     current_code = 1_ip       ! Current code
     kfl_custo    = 0          ! There is no customer
     kfl_preli    = 0          ! Not a preliminary run
     nprit        = huge(1_ip) ! Preliminary frequence
     kfl_rstar    = 0          ! Not a restart run
     kfl_rsfil    = 0          ! Type of restart file name
     kfl_timeline = 0          ! Output of timeline
     kfl_commu    = 0          ! Communication with Alya
     kfl_outfo    = 1          ! Output format
     kfl_latex    = 0          ! No latex output file
     lun_livei    = 6          ! Log file is screen
     kfl_memor    = 0          ! Do not Output memory
     kfl_varcount = 0          ! Do not Output variable memory counter
     kfl_timin    = 0          ! Timing
     kfl_lotme    = 0          ! We do not have a lot of memory
     kfl_freme    = 1          ! Master deallocates mesh-related memory in memgeo
     kfl_outpu    = 1          ! Log file on
     cpu_limit    = 1.0e20_rp  ! Default CPU limit
     title        = ' '        ! Problem title
     kfl_vtk      = 0          ! VTK output format way
     kfl_rread    = 0          ! Not reread ini data in READ_AND_RUN mode
     !
     ! Begin
     !
     call ecoute('RRUDAT')

     if( words(1) /= 'RUNDA' ) call runend('RRUDAT: WRONG RUN_DATA CARD')

     do while( words(1) /= 'ENDRU' )

        if( words(1) == 'ALYA ' ) then
           !
           ! Read and write title
           !
           title = trim(wname)

        else if( words(1) == 'CODE ' ) then
           !
           ! Read and write title
           !
           current_code = getint('CODE ',1_ip,'#My code number')

        else if( words(1) == 'RUNTY' ) then
           !
           ! Read the type of run
           !
           if( exists('PRELI') ) then
              kfl_preli = 1                                           ! Preliminary run
              nprit = getint('FREQU',1_ip,'#Preliminary frequency')   ! Preliminary frequency
           end if
           if( exists('RESTA') ) kfl_rstar = 1                        ! Restart run: initial
           if( exists('INITI') ) kfl_rstar = 1                        ! Restart run: initial
           if( exists('CONTI') ) kfl_rstar = 2                        ! Restart run: continue
           if( exists('INTER') ) kfl_rstar = 3                        ! Restart run: interpolate and initial
           if( exists('APPEN') ) kfl_rsfil = 1                        ! Append time step to file names

        else if( words(1) == 'COMMU' ) then
           if( words(2) == 'YES  '.or.words(2) == 'ON   ' ) kfl_commu = 1 ! Communication with Alya

        else if( words(1) == 'CPULI' ) then
           !
           ! Read the CPU limit
           !
           if( param(1)/=0.0_rp) cpu_limit = param(1)

        else if( words(1) == 'CUSTO' ) then
           !
           ! Read the customer
           !
           if( words(2) == 'MAREK' ) then
              kfl_custo = 1                                         ! Marek
           else if( words(2) == 'CFDW1' ) then
              kfl_custo =  2                                        ! Iberdrola: CFDWind1
           else if( words(2) == 'CFDW2' ) then
              kfl_custo =  3                                        ! Iberdrola: CFDWind2
           else if( words(2) == 'CFDW0' ) then
              kfl_custo = -2                                        ! Iberdrola: CFDWind1 + base field
           end if

        else if( words(1) == 'LOGFI' ) then
           !
           ! Log file
           !
           if( words(2) == 'YES  '.or.words(2) == 'ON   ' ) then
              kfl_outpu = 1
           else
              kfl_outpu = 0
           end if

        else if( words(1) == 'OUTPU' ) then
           !
           ! Output format
           !
           if( words(2) == 'GID  ' ) then                                 ! GiD Ascii (our format)
              kfl_outfo = 1
           else if( words(2) == 'GIDAS' ) then                            ! GiD Ascii (use GiD lib)
              kfl_outfo = 3
           else if( words(2) == 'GIDBI' ) then                            ! GiD Binary (use GiD lib)
              kfl_outfo = 4
           else if( words(2) == 'FEMVI' ) then                            ! Femview
              kfl_outfo = 2
           else if( words(2) == 'CGNS ' ) then                            ! CGNS
              kfl_outfo = 5
           else if( words(2) == 'GNUPL' ) then                            ! Gnuplot
              kfl_outfo = 6
           else if( words(2) == 'ALYA ' ) then                            ! Alya (ASCII and BIN)
              kfl_outfo = 8
              if( exists('BIN  ') ) kfl_outfo = 7
           else if( words(2) == 'ENSIG' .or. words(2) == 'VISIT'  ) then  ! ENSIGHT / VISIT (from 10 to 19)
              kfl_outfo = 10                                              ! default ASCII
              if (exists('BINAR')) kfl_outfo = 15
           else if( words(2) == 'VU   ' ) then                            ! VU
              kfl_outfo = 30
           else if( words(2) == 'HDF  ' .or. words(2) == 'HDF5 ' ) then   ! HDF
              kfl_outfo = 50
           else if( words(2) == 'VTK  ' ) then                            ! VTK (default BINARY, writing /scratch and cp /gpfs)
              kfl_outfo = 40
              if      ( words(3) == 'ASCII' ) then                        ! VTK ASCII
                 kfl_outfo = 41
              else if ( words(3) == 'GPFS' ) then
                 kfl_vtk = 1                                              ! VTK (writing in the /gpfs)
              end if
!           else if( words(2) == 'MPIAL' ) then                            ! MPIO format
!              kfl_outfo = 324
           end if

!        else if( words(1) == 'INPUT' ) then
!           if( words(2) == 'MPIAL' ) then                            ! Alya (ASCII and BIN)
!              kfl_inpfo = 324
!           end if
!        else if( words(1) == 'RESTA' ) then
!           if( words(2) == 'MPIAL' ) then                            ! Alya (ASCII and BIN)
!              kfl_rstfo = 324
!           end if
        else if( words(1) == 'LIVEI' ) then
           !
           ! Live information
           !
           if( words(2) == 'SCREE' ) then
              lun_livei=6
           else if( words(2) == 'FILE ' ) THEN
              lun_livei=16
           end if

        else if( words(1) == 'TIMEL' ) then
           !
           ! Output timeline
           !
           if( option('TIMEL') ) kfl_timeline = 1

        else if( exists('STATI') ) then

           call coupling_readat()  !< 2017SEP22

        else if( words(1) == 'LATEX' ) then
           !
           ! Latex file
           !
           if( words(2) == 'YES  '.or.words(2) == 'ON   ') kfl_latex=1

        else if( words(1) == 'MEMOR' ) then
           !
           ! Memory count
           !
           if( words(2) == 'YES  '.or.words(2) == 'ON   ') then
              kfl_memor = -1
              if( exists('VARIA') ) kfl_varcount = -1
           end if
           
        else if( words(1) == 'TIMIN' ) then
           if( words(2) == 'YES  '.or.words(2) == 'ON   ') kfl_timin = 1

        else if( words(1) == 'LOTOF' ) then
           if( words(2) == 'YES  '.or.words(2) == 'ON   ') kfl_lotme = 1

        else if( words(1) == 'FREEM' ) then
           if( words(2) == 'NO   '.or.words(2) == 'OFF  ') kfl_freme = 0

        ! Reread ini data in READ_AND_RUN mode
        else if( words(1) == 'REREA' ) then
           if( words(2) == 'YES  '.or.words(2) == 'ON   ') kfl_rread = 1

        end if
        call ecoute('RRUDAT')
     end do

  else

     lun_livei    = 0
     kfl_memor    = 0
     kfl_varcount = 0

  end if

end subroutine rrudat
