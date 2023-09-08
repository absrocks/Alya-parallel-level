!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_reaphy.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Read physical data
!> @details Read physical data
!> @} 
!------------------------------------------------------------------------
subroutine por_reaphy()
  use def_parame
  use def_inpout
  use def_master
  use def_porous
  use def_domain
  use def_kermod
  use mod_ker_space_time_function
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip)     :: irows,jmate,iwell,nposi,ntime
  real(rp)        :: dummr
  character(5)    :: wcond


  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     do kprsa_por=1,nprsa_por
        kfl_timei_por(kprsa_por) = 1                         ! Transient flow
     end do
     kfl_advec_por(1) = 0                                    ! Convection is off for pressure
     kfl_advec_por(2) = 1                                    ! Convection is on for saturation - for the moment they are not used
     denhy_por  = 0.0_rp                                     ! Do not substract hydrostatic component
     !
     ! Reach the section
     !
     call ecoute('por_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('por_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDPH')
        call ecoute('por_reaphy')
        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !
           call ecoute('por_reaphy')
           do while(words(1)/='ENDPR')

!              if(words(1)=='CONVE') then               ! Convective term
!                 if(exists('ON   ')) then
!                    kfl_advec_por = 1  ! This I can use for the pressure equation - the conv term taht I do not know whether to include or not 
!                 end if
!              end if
              call ecoute('por_reaphy')
           end do

        else if(words(1)=='PROPE') then
           !
           ! Allocate memory
           !
           call por_memphy(1_ip)

           call ecoute('por_reaphy')
           do while(words(1)/='ENDPR')

              if(words(1)=='ROCKC') then
                 !
                 ! ADOC[2]> ROCK_COMPRESSIBILITY:      real1                                           $ rock compressibility
                 ! ADOC[d]> ROCK_COMPRESSIBILITY:
                 ! ADOC[d]> Rock compressibility. Used to obtain the porosity from the initial porosity and the pressure.
                 !
                 comro_por = param(1)      

              else if(words(1)=='WATCO') then
                 !
                 ! ADOC[2]> WAT_COMPRESSIBILITY:      real1                                          $ water compressibility
                 ! ADOC[d]> WAT_COMPRESSIBILITY:
                 ! ADOC[d]> Water compressibility. Used to obtain Bw from Bw(Pref) and the pressure.
                 !
                 comwa_por = param(1)

              else if(words(1)=='OILCO') then
                 !
                 ! ADOC[2]> OIL_COMPRESSIBILITY:      real1                                          $ oil compressibility
                 ! ADOC[d]> OIL_COMPRESSIBILITY:
                 ! ADOC[d]> Oil compressibility. Used to obtain Bo from Bo(Pref) and the pressure.
                 !
                 comoi_por = param(1)

              else if(words(1)=='BWREF') then
                 !
                 ! ADOC[2]> BW_REFERENCE:      real1                                          $ Bw at reference Pressure
                 ! ADOC[d]> BW_REFERENCE:
                 ! ADOC[d]> Bw at reference Pressure. Used to obtain Bw from Bw(Pref) and the pressure.
                 !
                 bwref_por = param(1)

              else if(words(1)=='BOREF') then
                 !
                 ! ADOC[2]> BO_REFERENCE:      real1                                          $ Bo at reference Pressure
                 ! ADOC[d]> BO_REFERENCE:
                 ! ADOC[d]> Bo at reference Pressure. Used to obtain Bo from Bo(Pref) and the pressure.
                 !
                 boref_por = param(1)

              else if(words(1)=='PRREF') then
                 !
                 ! ADOC[2]> PR_REFERENCE:      real1                                          $ Reference Pressure
                 ! ADOC[d]> PR_REFERENCE:
                 ! ADOC[d]> Reference Pressure.
                 !
                 prref_por = param(1)

              else if(words(1)=='VISCW') then
                 !
                 ! ADOC[2]> VISC_WATER:      real1                                          $ Water Viscosity
                 ! ADOC[d]> VISC_WATER:
                 ! ADOC[d]> Water Viscosity.
                 !
                 muwat_por = param(1)

              else if(words(1)=='VISCO') then
                 !
                 ! ADOC[2]> VISC_OIL:      real1                                          $ Oil Viscosity
                 ! ADOC[d]> VISC_OIL:
                 ! ADOC[d]> Oil Viscosity.
                 !
                 muoil_por = param(1)

              else if(words(1)=='DENSW') then
                 !
                 ! ADOC[2]> DENS_WATER:      real1                                          $ Water Density
                 ! ADOC[d]> DENS_WATER:
                 ! ADOC[d]> Water Density.
                 !
                 denwa_por = param(1)

              else if(words(1)=='DENSO') then
                 !
                 ! ADOC[2]> DENS_OIL:      real1                                          $ Oil Density
                 ! ADOC[d]> DENS_OIL:
                 ! ADOC[d]> Oil Density.
                 !
                 denoi_por = param(1)

              else if(words(1)=='DENSH') then
                 !
                 ! ADOC[2]> DENS_HYDRO:      real1                                          $ Hydrostatic Density
                 ! ADOC[d]> DENS_HYDRO:
                 ! ADOC[d]> Density used to substract Hydrostatic component - th default is zero so that the Hydrostatic 
                 ! ADOC[d]> component is not substracted.
                 !
                 denhy_por = param(1)

              else if(words(1)=='PRINI') then
                 !
                 ! ADOC[2]> INITIAL_PRESSURE:      real1                                          $ Initial Pressure
                 ! ADOC[d]> INITIAL_PRESSURE:
                 ! ADOC[d]> Initial Pressure.
                 !
                 prini_por = param(1)

              else if(words(1)=='GRAVI') then
                 !
                 ! ADOC[2]> GRAVITY:              NORM= real, GX= real, GY= real, GZ= real                    $ Gravity accelaration
                 ! ADOC[d]> GRAVITY:
                 ! ADOC[d]> If this option is present, the gravity acceleration is added to the 
                 ! ADOC[d]> Navier-Stokes equations. Porous does not use Kermod's gravity, as some
                 ! ADOC[d]> modules could require gravity (Partis) but Porous could run without.
                 !
                 grnor_por    = getrea('NORM ',0.0_rp,'#Gravity norm')
                 gravi_por(1) = getrea('GX   ',0.0_rp,'#x-component of g')
                 gravi_por(2) = getrea('GY   ',0.0_rp,'#y-component of g')
                 gravi_por(3) = getrea('GZ   ',0.0_rp,'#z-component of g')
                 call vecuni(3_ip,gravi_por,dummr)

              else if(words(1)=='INFOR') then
                 !
                 ! ADOC[2]> INFORMATION_WELL ,    NUMBER_OF_WELLS= integer   , MAX_ROWS= integer,  MAX_HEIGHT= integer
                 ! ADOC[d]> WELL_NUMBER=iwell, WELL_RADIUS=wellr_por(iwell), NROWS=nroww_por(iwell), NHEIG=6, CONDITION:WATER_FLOW or 
                 !                                                  TOTAL_FLOW or PBH
                 ! ADOC[d]> wvalu_por(1,1,iwell)   wvalu_por(2,1,iwell)            *                         wvalu_por(nheig+1,1,iwell)
                 ! ADOC[d]> wvalu_por(1,2,iwell)   wvalu_por(2,2,iwell)            *                         wvalu_por(nheig+1,2,iwell)   
                 ! ADOC[d]> *                      * 
                 ! ADOC[d]> *                      *
                 ! ADOC[d]> wvalu_por(1,nroww_por(iwell),iwell)   wvalu_por(2,nroww_por(iwell),iwell)    *         *
                 ! ADOC[d]> END_WELL_NUMBER
                 ! ADOC[d]> WELL_NUMBER=iwell, WELL_RADIUS=wellr_por(iwell), NROWS=nroww_por(iwell) , NHEIG=6, CONDITION:PBH_MIDDLE
                 ! ADOC[d]>                                                                     (this was the only option in 2526)
                 ! ADOC[d]> wvalu_por(1,1,iwell)   wvalu_por(2,1,iwell) 
                 ! ADOC[d]> wvalu_por(1,2,iwell)   wvalu_por(2,2,iwell) 
                 ! ADOC[d]> *                      * 
                 ! ADOC[d]> *                      *
                 ! ADOC[d]> wvalu_por(1,nroww_por(iwell),iwell)   wvalu_por(2,nroww_por(iwell),iwell) 
                 ! ADOC[d]> END_WELL_NUMBER
                 ! ADOC[d]> END_INFORMATION_WELL
                 !
                 nwell_por = getint('NUMBE',1_ip,'#Number of wells')
                 mroww_por = getint('MAXRO',20_ip,'#Max. Number of rows for well')
                 mheiw_por = getint('MAXHE',20_ip,'#Max. Number of nodes per well') ! for
                 !
                 ! Allocate memory for wells
                 !
                 call por_memphy(4_ip)
                 call ecoute('por_reaphy')

                 do while(words(1)/='ENDIN')
                    !
                    ! Type 3, 4: pbh prescribed
                    ! Type 1, 2: q   prescribed
                    !
                    iwell                    = abs(getint('WELLN',1_ip,'#Number of well'))
                    tywel_por(iwell) % ntime = getint('NROWS',1_ip,'#Number of rows for the bottom hole pressure table')
                    tywel_por(iwell) % nposi = getint('NHEIG',1_ip,'#Number vertical points that form the well as read from data file')
                    tywel_por(iwell) % radiu = getrea('WELLR',1.0_rp,'#Well radius')
                    !
                    ! Allocate memory for this well
                    !
                    igene = iwell
                    call por_memphy(5_ip)
                    !
                    ! Read well
                    !
                    if( exists('CONDI') ) then    ! Borrowed from reacod AXES
                       wcond = getcha('CONDI','     ','#Condition')
                       if( wcond == 'WATER' ) then
                          tywel_por(iwell) % itype = 1
                       else if( wcond == 'TOTAL' ) then
                          tywel_por(iwell) % itype = 2
                       else if( wcond == 'OIL  ' ) then
                          tywel_por(iwell) % itype = 5
                       else if( wcond == 'PBH  ' ) then
                          tywel_por(iwell) % itype = 3
                       else if( wcond == 'PBHMI' ) then
                          tywel_por(iwell) % itype = 4
                       end if
                    else
                       kfl_wellc_por(iwell) = 4  ! This is the default - the only option prior to 2526
                       tywel_por(iwell) % itype = 4
                    end if
                    
                    nposi = tywel_por(iwell) % nposi
                    ntime = tywel_por(iwell) % ntime

                    if( tywel_por(iwell) % itype == 3 .or. tywel_por(iwell) % itype == 4 ) then
                       call ecoute('por_reaphy')
                       tywel_por(iwell) % pbh_coord(1:nposi) = param(2:1+nposi)
                    end if

                   do irows = 1,ntime
                       call ecoute('por_reaphy')
                       if( tywel_por(iwell) % itype == 3 .or. tywel_por(iwell) % itype == 4 ) then
                          tywel_por(iwell) % pbh_table(irows,1:nposi+1) = param(1:nposi+1) 
                       else
                          tywel_por(iwell) % q_table(irows,1:2) = param(1:2) 
                       end if
                     end do
                    call ecoute('por_reaphy')
                    if( (words(1)/='ENDWE') .and. (words(1)/='ENDIN') ) &
                         call runend('POR_REAPHY: SOMETHING WRONG IN WELL NUMBER '//intost(iwell))
                    if (words(1)/='ENDIN')  call ecoute('por_reaphy')

                 end do

!!$                 do while(words(1)/='ENDIN')
!!$                    iwell = abs(getint('WELLN',1_ip,'#Number of well'))
!!$                    rwell_por(iwell) = getrea('WELLR',1.0_rp,'#Well radius')
!!$                    nroww_por(iwell) = getint('NROWS',1_ip,'#Number of rows for the bottom hole pressure table')
!!$                    nheiw_por(iwell) = getint('NHEIG',1_ip,'#Number vertical points that form the well as read from data file')
!!$                    if( exists('CONDI') ) then    ! Borrowed from reacod AXES
!!$                       wcond = getcha('CONDI','     ','#Condition')
!!$                       if( wcond == 'WATER' ) then
!!$                          kfl_wellc_por(iwell) = 1
!!$                       else if( wcond == 'TOTAL' ) then
!!$                          kfl_wellc_por(iwell) = 2
!!$                       else if( wcond == 'PBH  ' ) then
!!$                          kfl_wellc_por(iwell) = 3
!!$                       else if( wcond == 'PBHMI' ) then
!!$                          kfl_wellc_por(iwell) = 4
!!$                       end if
!!$                    else
!!$                       kfl_wellc_por(iwell) = 4  ! This is the default - the only option prior to 2526
!!$                    end if
!!$                    !print*,'iwell,kfl_wellc_por(iwell),nheiw_por(iwell)',iwell,kfl_wellc_por(iwell),nheiw_por(iwell)
!!$                    do irows = 1,nroww_por(iwell)
!!$                       call ecoute('por_reaphy')
!!$                       if ( kfl_wellc_por(iwell) == 4) then
!!$                          wvalu_por(1:2,irows,iwell) = param(1:2)
!!$                       else
!!$                          iauxi = nheiw_por(iwell) + 1_ip
!!$                          wvalu_por(1:iauxi,irows,iwell) = param(1:iauxi)
!!$                       end if
!!$                    end do
!!$                    call ecoute('por_reaphy')
!!$                    if( (words(1)/='ENDWE') .and. (words(1)/='ENDIN') ) &
!!$                         call runend('POR_REAPHY: SOMETHING WRONG IN WELL_NUMBER')
!!$                    if (words(1)/='ENDIN')  call ecoute('por_reaphy')
!!$                 end do

              else if(words(1)=='KRELA') then
                 !
                 ! ADOC[2]> K_RELATIVE       NUM_MATER= integer, MAX_ROWS= integer 
                 ! ADOC[d]> K_RELATIVE:
                 ! ADOC[d]> krw & kwr as a table in function of Sw
                 !
                 ! ADOC[3]> MATERIAL=int1, N_ROWS=int2.                                  $  
                 ! ADOC[3]> Sw1(irow=1) krw(irow=1)  kro(irow=1) 
                 ! ADOC[3]> Sw1(irow=2) krw(irow=2)  kro(irow=2) 
                 ! ADOC[3]> ..... 
                 ! ADOC[3]> Sw1(irow=nrows_por) krw(irow=nrows_por)  kro(irow=nrows_por) 
                 ! ADOC[3]> END_MATERIAL
                 ! ADOC[d]> ... repeat for all the materials
                 ! ADOC[d]> END_K_RELATIVE
                 !
                 nmate_por = getint('NUMMA',7_ip,'#Number of material')
                 mrows_por = getint('MAXRO',20_ip,'#Max. Number of rows')
                 !
                 ! Allocate memory for tkrel_por
                 !
                 call por_memphy(3_ip)      !tkrel_por(2+1,mrows_por,nmate_por)  
                 call ecoute('por_reaphy')

                 do while(words(1)/='ENDKR')
                    jmate = getint('MATER',1_ip,'#Number of material')
                    nrows_por(jmate) = getint('NROWS',20_ip,'#Number of rows')
                    do irows=1,nrows_por(jmate)
                       call ecoute('por_reaphy')
                       tkrel_por(1:3,irows,jmate) = param(1:3)
                    end do
                    call ecoute('por_reaphy')
                    if( (words(1)/='ENDKR') .and. (words(1)/='ENDMA') ) &
                         call runend('POR_REAPHY: SOMETHING WRONG IN K_RELATIVE')
                    if (words(1)/='ENDKR')  call ecoute('por_reaphy')
                 end do

              end if
              call ecoute('por_reaphy')
           end do
        end if
     end do

  end if

end subroutine por_reaphy
