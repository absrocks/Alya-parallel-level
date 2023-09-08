subroutine qua_openfi(itask)
  !------------------------------------------------------------------------
  !****f* Quanty/qua_openfi
  ! NAME 
  !    qua_openfi
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
  !    GET_ENVIRONMENT_VARIABLE
  !    iofile
  ! USED BY
  !    qua_turnon
  !------------------------------------------------------------------------
  use def_quanty
  use def_master
  use def_domain
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask 
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if(kfl_paral<=0) then

     if(kfl_rstar==2) then
        statu='old'
        forma='formatted'
        posit='append'
     else
        statu='unknown'
        forma='formatted'
        posit='asis'
     end if

     select case (itask)

     case(1_ip)
        !
        ! Open files that are always needed
        !
        if (kfl_naked==0) then
           !
           ! Encapsulated, then get names from the environment   
           !  
           
		   call GET_ENVIRONMENT_VARIABLE('FOR1501',fil_pdata_qua)
           call GET_ENVIRONMENT_VARIABLE('FOR1502',fil_outpu_qua)
           call GET_ENVIRONMENT_VARIABLE('FOR1503',fil_conve_qua)
           call GET_ENVIRONMENT_VARIABLE('FOR1504',fil_solve_qua)
           call GET_ENVIRONMENT_VARIABLE('FOR1505',fil_rstar_qua)
           call GET_ENVIRONMENT_VARIABLE('FOR1506',fil_witne_qua)
           call GET_ENVIRONMENT_VARIABLE('FOR1507',fil_ppseu_qua)

        else if (kfl_naked==1) then
           !
           ! Naked, then compose the names     
           !
           fil_pdata_qua = adjustl(trim(namda))//'.'//exmod(modul)//'.dat'
           fil_outpu_qua = adjustl(trim(namda))//'.'//exmod(modul)//'.log'
           fil_conve_qua = adjustl(trim(namda))//'.'//exmod(modul)//'.cvg'
           fil_solve_qua = adjustl(trim(namda))//'.'//exmod(modul)//'.sol'
           fil_rstar_qua = adjustl(trim(namda))//'.'//exmod(modul)//'.rst'
           fil_witne_qua = adjustl(trim(namda))//'.'//exmod(modul)//'.wit'
           fil_ppseu_qua = adjustl(trim(namda))//'.'//exmod(modul)//'.seu'
        end if
        !
        ! Open files
        !
        call iofile(zero,lun_pdata_qua,fil_pdata_qua,'QUANTY DATA',         'old')
        call iofile(zero,lun_outpu_qua,fil_outpu_qua,'QUANTY OUTPUT',       statu,forma,posit)
        call iofile(zero,lun_conve_qua,fil_conve_qua,'QUANTY CONVERGENCE',  statu,forma,posit)
        call iofile(zero,lun_solve_qua,fil_solve_qua,'QUANTY SOLVER',       statu,forma,posit)
        call iofile(zero,lun_rstar_qua,fil_rstar_qua,'QUANTY RESTAR',       statu,forma,posit)
        call iofile(zero,lun_witne_qua,fil_witne_qua,'QUANTY WITNESS',      statu,forma,posit)
       ! call iofile(zero,lun_ppseu_qua,fil_ppseu_qua,'QUANTY PPSEUDO',      statu,forma,posit)

     case(2_ip)
        !
        ! Open files needed occasionally
        !
        !
        ! Solver convergence file
        !
!        if(solve_qua(1)%kfl_cvgso/=0) &
!             call iofile(zero,solve_qua(1)%lun_cvgso,fil_cvgso_qua,'QUANTY SOLVER CONVER',statu,forma,posit)  
        
     case(3_ip)
        !
        ! Close data file
        !
        call iofile(two,lun_pdata_qua,' ','QUANTY DATA')

     end select

  end if

end subroutine qua_openfi

