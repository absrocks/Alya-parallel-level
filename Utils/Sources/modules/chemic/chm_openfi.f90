subroutine chm_openfi(itask)
  !------------------------------------------------------------------------
  !****f* partis/chm_openfi
  ! NAME 
  !    chm_openfi
  ! DESCRIPTION
  !    This subroutine gets ALL the file names and open them to be used by 
  !    the module in two possible ways:
  ! 
  !    1. Recalling them from the environment, when Alya is launched
  !    encapsulated in a shell script, or
  ! 
  !    2. Composing the names out of the problem name which is given as argument
  !    when the binary file Alya is launched "naked".  
  ! USES
  ! USED BY
  !    chm_turnon
  !------------------------------------------------------------------------
  use def_chemic
  use def_parame
  use def_master
  use def_domain
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask 
  character(150)          :: fil_resu1_chm,fil_resu2_chm,fil_resu3_chm
  character(150)          :: fil_spcvg_chm
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if( INOTSLAVE ) then

     if( kfl_rstar == 2 ) then
        statu = 'old'
        forma = 'formatted'
        posit = 'append'
     else
        statu = 'unknown'
        forma = 'formatted'
        posit = 'asis'
     end if

     select case (itask)

     case (2_ip)
        !
        ! Open files needed occasionally
        !
        if( kfl_naked == 0 ) then
           call GET_ENVIRONMENT_VARIABLE('FOR1910',fil_sized_chm)     
           call GET_ENVIRONMENT_VARIABLE('FOR1911',fil_times_chm)     
           call GET_ENVIRONMENT_VARIABLE('FOR1912',fil_time2_chm)     
           call GET_ENVIRONMENT_VARIABLE('FOR1930',fil_remet_chm)     
           call GET_ENVIRONMENT_VARIABLE('FOR1931',fil_resou_chm)      
           call GET_ENVIRONMENT_VARIABLE('FOR1932',fil_spcvg_chm)      
        else
           fil_sized_chm = adjustl(trim(namda))//'-size-distribution.'//exmod(modul)//'.res'   
           fil_times_chm = adjustl(trim(namda))//'-time-step.'        //exmod(modul)//'.res'   
           fil_time2_chm = adjustl(trim(namda))//'-time-step-target.' //exmod(modul)//'.res'   
           fil_remet_chm = adjustl(trim(namda))//'.'                  //exmod(modul)//'.met'
           fil_resou_chm = adjustl(trim(namda))//'.'                  //exmod(modul)//'.src'
           fil_spcvg_chm = adjustl(trim(namda))//'-species.'          //exmod(modul)//'.cvg'
        end if
        !
        ! Size distribution
        !
        if( kfl_sized_chm == 1 ) then
           call iofile(zero,lun_sized_chm,fil_sized_chm,namod(modul)//' SIZE DISTRIBUTION') 
        end if
        !
        ! Time step strategy
        !
        call iofile(zero,lun_times_chm,fil_times_chm,namod(modul)//' TIME STEP INFORMATION') 
        if( kfl_dttar_chm >= 1 ) then
           call iofile(zero,lun_time2_chm,fil_time2_chm,namod(modul)//' TIME STEP TARGET') 
        end if
        !
        ! Specific problems
        !
        if( wprob_chm == 'COIN1' ) then

           if( kfl_naked == 0 ) then
              call GET_ENVIRONMENT_VARIABLE('FOR1921',fil_resu1_chm)     
              call GET_ENVIRONMENT_VARIABLE('FOR1922',fil_resu2_chm)     
           else
              fil_resu1_chm = adjustl(trim(namda))//'-integral-evolution.'   //exmod(modul)//'.res'   
              fil_resu2_chm = adjustl(trim(namda))//'-size-distribution-In.' //exmod(modul)//'.res'   
           end if
           call iofile(zero,lun_resu1_chm,fil_resu1_chm,namod(modul)//' RESULT 1') 
           call iofile(zero,lun_resu2_chm,fil_resu2_chm,namod(modul)//' RESULT 2') 

        else if( wprob_chm == 'INVN1' ) then

           if( kfl_naked == 0 ) then
              call GET_ENVIRONMENT_VARIABLE('FOR1921',fil_resu1_chm)     
              call GET_ENVIRONMENT_VARIABLE('FOR1922',fil_resu2_chm)     
              call GET_ENVIRONMENT_VARIABLE('FOR1923',fil_resu3_chm)     
           else
              fil_resu1_chm = adjustl(trim(namda))//'-integral-evolution.'   //exmod(modul)//'.res'   
              fil_resu2_chm = adjustl(trim(namda))//'-size-distribution-In.' //exmod(modul)//'.res'   
              fil_resu3_chm = adjustl(trim(namda))//'-size-distribution-Vn.' //exmod(modul)//'.res'   
           end if
           call iofile(zero,lun_resu1_chm,fil_resu1_chm,namod(modul)//' RESULT 1') 
           call iofile(zero,lun_resu2_chm,fil_resu2_chm,namod(modul)//' RESULT 2') 
           call iofile(zero,lun_resu3_chm,fil_resu3_chm,namod(modul)//' RESULT 3') 

        else if( wprob_chm == 'INTE1' ) then

           if( kfl_naked == 0 ) then 
              call GET_ENVIRONMENT_VARIABLE('FOR1921',fil_resu1_chm)    
           else
              fil_resu1_chm = adjustl(trim(namda))//'-defects-temperature.'//exmod(modul)//'.res'  
           end if
           call iofile(zero,lun_resu1_chm,fil_resu1_chm,namod(modul)//' RESULT 1') 

        end if
        !
        ! METEO model: Properties file
        !
        if( kfl_meteo_chm == 1 ) then         ! debug
           continue
        else if( kfl_meteo_chm == 2 ) then    ! ASCII
           fil_remet_chm = adjustl(trim(namda))//'.'//exmod(modul)//'.met'
           call iofile(zero,lun_remet_chm,fil_remet_chm,'CHEMIC ASCII METEO FILE',statu,forma,posit)  
        else if( kfl_meteo_chm == 3 ) then    ! netCDF
           fil_remet_chm = adjustl(trim(namda))//'.'//exmod(modul)//'.met.nc'
           call runend('ERROR: BIN CHEMIC METEO FILE TYPE NOT IMPLEMENTED YET')
        end if
        !
        ! METEO model: Source file
        !
        if( kfl_sourc_chm == -1 ) then        ! debug
           continue 
        else if( kfl_sourc_chm == -2 ) then   ! ASCII
           fil_resou_chm = adjustl(trim(namda))//'.'//exmod(modul)//'.src'  
           call iofile(zero,lun_resou_chm,fil_resou_chm,namod(modul)//' ASCII SOURCE FILE',statu,forma,posit)  
        else if( kfl_sourc_chm == -3 ) then
           call runend('Type of source not implemented yet')
        end if
        !
        ! Mechano-biological: convergence file for all species
        !
        if( kfl_model_chm == 3 ) then
           call iofile(zero,lun_spcvg_chm,fil_spcvg_chm,namod(modul)//' SEPCIES CONVERGENCE',statu,forma,posit)  
        end if

     case (4_ip)
        !
        ! Close files
        !
        if( kfl_sized_chm == 1 ) then
           call iofile(two,lun_sized_chm,fil_sized_chm,namod(modul)//' SIZE DISTRIBUTION') 
        end if
        call iofile(two,lun_times_chm,fil_times_chm,namod(modul)//' TIME STEP') 

     end select

  end if

end subroutine chm_openfi

