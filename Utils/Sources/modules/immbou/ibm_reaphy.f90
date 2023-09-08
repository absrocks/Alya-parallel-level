subroutine ibm_reaphy()
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_reaphy
  ! NAME
  !    ibm_reaphy
  ! DESCRIPTION
  !    Read physical IB
  !    kfl_cofor = 0 ... No viscous no pressure forces
  !                1 ... Viscous forces 
  !                2 ... Pressure forces 
  !                3 ... Viscous and pressure forces 
  !    kfl_coibm = 0 ... No coupling with kernel
  !                1 ... Coupling with kernel
  ! OUTPUT
  ! USED BY
  !    ibm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use def_inpout
  use mod_ecoute, only :  ecoute
  implicit none
  real(rp) :: dummr

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_timei_ibm =  0             ! IB is moving
     kfl_rotib_ibm =  1             ! IB enable rotation
     kfl_linib_ibm =  1             ! IB enable linear motion
     nstro_ibm =  0                 ! IB start rotation from step 0
     nstli_ibm =  0                 ! IB start linear motion from step 0
     kfl_mvext_ibm =  0             ! IB do not move exterior
     kfl_ralei_ibm =  0             ! IB no Raleigh damping
     nstra_ibm =  1000000           ! IB start eliminating Raleigh damping at a very high step
     nenra_ibm =  2000000           ! IB end eliminating Raleigh damping at a very high step
     kfl_staib_ibm =  0             ! Starting IB
     kfl_colli_ibm =  0             ! Collisions
     kfl_coibm     =  0             ! Coupling with kernel  (ALE, interpolation) 
     kfl_cofor     = -1             ! Coupling with modules (forces)

     kfl_grafo_ibm =  0             ! Gravity force
     kfl_catfo_ibm =  0             ! Catamaran force
     kfl_buofo_ibm =  0             ! Buoyancy force
     kfl_drafo_ibm =  0             ! Drag force
     kfl_extfo_ibm =  0             ! External force 

     staib_ibm     =  0.0_rp        ! Starting time
     spher_ibm     =  1.0_rp        ! Sphericity 
     xline_ibm     =  0.0_rp        ! Linear motion - now no longer read but calculated fron nstli_ibm
     xrota_ibm     =  0.0_rp        ! Rotation motion - now no longer read but calculated fron nstro_ibm
     !
     ! Reach the section 
     !
     call ecoute('ibm_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('ibm_reaphy')
     end do
     call ecoute('ibm_reaphy')

     do while(words(1)/='ENDPH') 

        if(words(1)=='MOTIO' .or. words(1) == 'TIMEV' ) then
           if( words(2) == 'OFF  ' .or. words(2) == 'NO   ' ) then
              kfl_timei_ibm = 0
           else if( words(2) == 'ON   ' .or. words(2) == 'EXTRA' ) then
              kfl_timei_ibm = 1
           end if

        else if(words(1)=='ROTAT' ) then
           if( words(2) == 'OFF  ' .or. words(2) == 'NO   ' ) then
              kfl_rotib_ibm = 0
              xrota_ibm     = 0.0_rp
           end if
           if(words(2)=='START')   nstro_ibm(1:3) = int(param(2:4),ip)           


        else if(words(1)=='TRANS' ) then
           if( words(2) == 'OFF  ' .or. words(2) == 'NO   ' ) then
              kfl_linib_ibm = 0
              xline_ibm     = 0.0_rp
           end if
           if(words(2)=='START')   nstli_ibm(1:3) = int(param(2:4),ip)


        else if(words(1)=='MOVEE' ) then
              kfl_mvext_ibm = 1

        else if(words(1)=='RALEI' ) then
           kfl_ralei_ibm = 1
           ralei_ibm     = 1.0_rp
           if(words(2)=='VALUE')   ralei_ibm = param(2)
           if(words(3)=='START')   nstra_ibm = int(param(3),ip)
           if(words(4)=='END  ')   nenra_ibm = int(param(4),ip)

        else if( words(1) == 'FORCE' ) then
           !
           ! Forces
           !
           if( words(2) == 'ALL  ' ) then
              kfl_grafo_ibm = 1                       
              kfl_buofo_ibm = 1
              kfl_drafo_ibm = 2
              kfl_extfo_ibm = 1
              if( exists('FUNCT') ) &
                   kfl_extfo_ibm = getint('FUNCT',1_ip,'#EXTERNAL FORCE FUNCTION NUMBER')
              if( exists('SPHER') ) &
                   spher_ibm = getrea('SPHER',1.0_rp,'#SPHERICITY OF THE PARTICLE TYPE')
           else
              if( exists('CATAM') ) kfl_catfo_ibm = 1
              if( exists('GRAVI') ) kfl_grafo_ibm = 1
              if( exists('BUOYA') ) kfl_buofo_ibm = 1
              if( exists('DRAG ') ) then
                 spher_ibm  = getrea('SPHER',1.0_rp,'#SPHERICITY OF THE PARTICLE TYPE')
                 if( exists('GANSE') ) then
                    kfl_cofor     = -1 
                    kfl_drafo_ibm =  2
                 else if( exists('CHENG') ) then
                    kfl_cofor     = -1 
                    kfl_drafo_ibm =  1
                 else if( exists('ARAST') ) then
                    kfl_cofor     = -1 
                    kfl_drafo_ibm =  3
                 else if( exists('WILSO') ) then
                    kfl_cofor     = -1 
                    kfl_drafo_ibm =  4
                 else if( exists('MODUL') ) then
                    kfl_drafo_ibm =  0
                    if( exists('VISCO') .and. exists('PRESS') ) then
                       kfl_cofor =  3
                    else if( exists('VISCO') ) then
                       kfl_cofor =  1
                    else if( exists('PRESS') ) then
                       kfl_cofor =  2
                    end if
                 else
                    kfl_cofor     = -1 
                    kfl_drafo_ibm =  2
                 end if
              end if
              if( exists('EXTER') ) then
                 kfl_extfo_ibm = getint('EXTER',1_ip,'#EXTERNAL FORCE FUNCTION NUMBER')
              end if
           end if

        else if(words(1)=='KERNE' ) then
           if( words(2) == 'ON   ') then
              kfl_coibm = 1
           else if( words(2) == 'OFF  ') then
              kfl_coibm = 0
           end if

        else if(words(1)=='START' ) then
           if( words(2) == 'ATTIM' ) then
              staib_ibm = getrea('ATTIM',0.0_rp,'#Starting stim')
           end if
           if( words(2) == 'ATSTE' ) then
              kfl_staib_ibm = getint('ATSTE',0_ip,'#Starting stim')
           end if

        else if(words(1)=='COLLI' ) then
           if( words(2) == 'OFF  ' .or. words(2) == 'NO   ' ) then
              kfl_colli_ibm = 0
           else if( words(2) == 'ON   ' .or. words(2) == 'YES  ' ) then
              kfl_colli_ibm = 1
           end if

        end if
        call ecoute('ibm_reaphy')
     end do

  end if

end subroutine ibm_reaphy
