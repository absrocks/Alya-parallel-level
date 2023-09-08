subroutine pts_reaous()
  !------------------------------------------------------------------------
  !****f* Partis/pts_reaous
  ! NAME 
  !    pts_reaous
  ! DESCRIPTION
  !    This routine reads data
  ! OUTPUT
  ! USES
  ! USED BY
  !    Reapro
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_kermod
  use def_inpout
  use def_domain
  use def_partis
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: dummi,ivari

  if( INOTSLAVE ) then

     kfl_posla_pts = 3                                    ! Lagrangian particles: postprocess
     kfl_oudep_pts = 0                                    ! Do not output deposition map     
     kfl_depos_surface_pts = 0                            ! Do not output deposition surface
     kfl_oufre_pts = 1                                    ! Output frequency 
     kfl_dbfre_pts = 0                                    ! Do not output DB frequency 
     kfl_exacs_pts = 0                                    ! Exact solution
     
     if( kfl_rstar == 0 ) then                            ! Independent Restart of particles
        kfl_rstar_pts = 0
     else
        kfl_rstar_pts = 1
     end if

     call ecoute('pts_reaous')
     do while( words(1) /= 'OUTPU' )
        call ecoute('pts_reaous')
     end do
     call ecoute('pts_reaous')

     do while(words(1) /= 'ENDOU' )

        call posdef(2_ip,dummi)

        !Add PTSRES BINARY in the postprocess section to save pts.res in binary rather than text
        if( words(1) == 'PTSRE' ) then
           if( exists('BINAR') ) then
              !save pts.res in binary. The integers will have IP precision + 255 char header
              kfl_ptsres_binary = 1
           end if
           if( exists('SPLIT') ) then
              !save pts.res separate files, one per timestep. Files will be pts.00000x.res. 
              kfl_ptsres_split = 1
           end if

        end if


        if( words(1) == 'LEVEL' ) then
           kfl_posla_pts = getint('LEVEL',0_ip,'#POSTPROCESS LEVEL FOR LAGRANGIAN PARTICLE')
        end if
        if( words(1) == 'FREQU' ) then
           kfl_oufre_pts = getint('FREQU',1_ip,'#OUTPUT FREQUENCY FOR PARTICLES')
           kfl_oufre_pts = max(1_ip,kfl_oufre_pts)
        end if

        if( words(1) == 'DBFRE' ) then
           kfl_dbfre_pts = getint('DBFRE',0_ip,'#OUTPUT FREQUENCY FOR DB PARTICLES OPTION')
        end if

        if( words(1) == 'OUTPU' ) then

           if( words(2) == 'DEPOS') then
              kfl_oudep_pts = 1
              if( words(3) == 'SURFA') then
                 kfl_depos_surface_pts = 1
              endif
           else if(words(2)=='ERROR') then
              kfl_exacs_pts = getint('SOLUT',1_ip,'#Exact solution')

           end if
        end if

        if( words(1) == 'DATAB' ) then
           !
           ! Data base
           !
           call  pts_readbp('pts_reaous')
           do while(words(1) /= 'ENDDA' )
              write(*,*) 'WORD1:', words(1)
              if( words(1) == 'CUBES' ) then
                 dbSettings%kfl_db_deep = getint('CUBES',1_ip,'#Number of leves on db')
              else if( words(1) == 'CONNE' ) then
                 dbSettings%kfl_db_url_conn = wordsdb(2)
              else if( words(1) == 'CPORT' ) then
                 dbSettings%port = getint('CPORT',1_ip,'#Number connection port')
              else if( words(1) == 'MAXIN' ) then
                 dbSettings%maxInserts = getint('MAXIN',1_ip,'#Number of inserts on db')
              else if( words(1) == 'BATCH' ) then
                 dbSettings%numParticlesBatch = getint('BATCH',1_ip,'#Number of batchs on db')
              end if
              call pts_readbp('pts_reaous')
           end do
           
        else if( words(1) == 'VARIA' ) then
           !
           ! List of variables to output
           !
              if( words(2) /= 'ADD  ' ) postprocess_var_pts = .false.
           postprocess_var_pts(pts_name_to_variable_number('ILAGR')) = .true.
           postprocess_var_pts(pts_name_to_variable_number('ITYPE')) = .true.
           postprocess_var_pts(pts_name_to_variable_number('EXIST')) = .true.
           postprocess_var_pts(pts_name_to_variable_number('T'))     = .true.           
           call ecoute('pts_reaous')
           do while( words(1) /= 'ENDVA' )
              ivari = pts_name_to_variable_number(wname)
              if( ivari /= 0 ) postprocess_var_pts(ivari) = .true.
              if( words(1) == 'COORD' ) then
                 postprocess_var_pts(pts_name_to_variable_number('COORX')) = .true.
                 postprocess_var_pts(pts_name_to_variable_number('COORY')) = .true.
                 postprocess_var_pts(pts_name_to_variable_number('COORZ')) = .true.
              end if
              if( words(1) == 'VELOC' ) then
                 postprocess_var_pts(pts_name_to_variable_number('VELOX')) = .true.
                 postprocess_var_pts(pts_name_to_variable_number('VELOY')) = .true.
                 postprocess_var_pts(pts_name_to_variable_number('VELOZ')) = .true.
              end if
              if( words(1) == 'ACCEL' ) then
                 postprocess_var_pts(pts_name_to_variable_number('ACCEX')) = .true.
                 postprocess_var_pts(pts_name_to_variable_number('ACCEY')) = .true.
                 postprocess_var_pts(pts_name_to_variable_number('ACCEZ')) = .true.
              end if
              if( words(1) == 'STK1' ) then
                 postprocess_var_pts(pts_name_to_variable_number('STK1')) = .true.
              end if
              if( words(1) == 'STK2' ) then
                 postprocess_var_pts(pts_name_to_variable_number('STK2')) = .true.
              end if
              call ecoute('pts_reaous')
           end do
           
        else if( words(1) == 'VARIA' ) then
           !
           ! List of variables to output
           !
           if( words(2) /= 'ADD  ' ) deposition_var_pts = .false.
           deposition_var_pts(pts_name_to_variable_number('ILAGR')) = .true.
           deposition_var_pts(pts_name_to_variable_number('ITYPE')) = .true.
           deposition_var_pts(pts_name_to_variable_number('EXIST')) = .true.
           deposition_var_pts(pts_name_to_variable_number('T'))     = .true.           
           call ecoute('pts_reaous')
           do while( words(1) /= 'ENDVA' )
              ivari = pts_name_to_variable_number(wname)
              if( ivari /= 0 ) deposition_var_pts(ivari) = .true.
              if( words(1) == 'COORD' ) then
                 deposition_var_pts(pts_name_to_variable_number('COORX')) = .true.
                 deposition_var_pts(pts_name_to_variable_number('COORY')) = .true.
                 deposition_var_pts(pts_name_to_variable_number('COORZ')) = .true.
              end if
              if( words(1) == 'VELOC' ) then
                 deposition_var_pts(pts_name_to_variable_number('VELOX')) = .true.
                 deposition_var_pts(pts_name_to_variable_number('VELOY')) = .true.
                 deposition_var_pts(pts_name_to_variable_number('VELOZ')) = .true.
              end if
              if( words(1) == 'ACCEL' ) then
                 deposition_var_pts(pts_name_to_variable_number('ACCEX')) = .true.
                 deposition_var_pts(pts_name_to_variable_number('ACCEY')) = .true.
                 deposition_var_pts(pts_name_to_variable_number('ACCEZ')) = .true.
              end if
              call ecoute('pts_reaous')
           end do
           
        else if( words(1) == 'RESTART' ) then
           !
           ! Restart
           !
           if(exists('ON   ')) kfl_rstar_pts = 1 
           if(exists('OFF  ')) kfl_rstar_pts = 0
        end if

        call ecoute('pts_reaous')
     end do

  end if
  
end subroutine pts_reaous

