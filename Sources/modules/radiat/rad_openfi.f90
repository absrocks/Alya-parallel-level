subroutine rad_openfi(itask)
  !------------------------------------------------------------------------
  !****f* Radiat/rad_openfi
  ! NAME 
  !    rad_openfi
  ! DESCRIPTION
  !    This subroutine gets ALL the file names and open them to be used by 
  !    the module in 2_ip possible ways:
  ! 
  !    1. Recalling them from the environment, when Alya is launched
  !    encapsulated in a shell script, or
  ! 
  !    2. Composing the names out of the problem name which is given as argument
  !    when the binary file Alya is launched "naked".  
  ! USES
  ! USED BY
  !    rad_turnon
  !------------------------------------------------------------------------
  use def_radiat
  use def_master
  use def_domain
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask 
  character(150)          :: fil_pdata_rad,fil_outpu_rad
  character(150)          :: fil_solve_rad,fil_splot_rad,fil_psmat_rad
  character(150)          :: fil_setse_rad,fil_setsb_rad,fil_cvgso_rad
  character(150)          :: fil_setsn_rad,fil_bound_rad,fil_funck_rad
  character(150)          :: fil_funcc_rad,fil_witne_rad,fil_intbc_rad
  character(150)          :: fil_dynlo_rad,fil_dynre_rad
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit
  character(20)           :: wmate

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

     case (2_ip)
        !
        ! Open files needed occasionally
        !
        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR210',fil_bound_rad)   
           call GET_ENVIRONMENT_VARIABLE('FOR214',fil_psmat_rad)   
           call GET_ENVIRONMENT_VARIABLE('FOR221',fil_intbc_rad)      
           call GET_ENVIRONMENT_VARIABLE('FOR222',fil_dynin_rad)      
           call GET_ENVIRONMENT_VARIABLE('FOR223',fil_dynou_rad)      
           call GET_ENVIRONMENT_VARIABLE('FOR224',fil_dynlo_rad)      
           call GET_ENVIRONMENT_VARIABLE('FOR225',fil_dynre_rad)      
           call GET_ENVIRONMENT_VARIABLE('FOR232',fil_splot_rad)      
           call GET_ENVIRONMENT_VARIABLE('FOR233',fil_ramsh_rad)      
           call GET_ENVIRONMENT_VARIABLE('FOR234',fil_rares_rad)       
        else
           fil_bound_rad = adjustl(trim(namda))//'.'                 //exmod(modul)//'.bcs'
           fil_psmat_rad = adjustl(trim(namda))//'-matrix.'          //exmod(modul)//'.ps'    
           fil_intbc_rad = adjustl(trim(namda))//'-bcinterpolation.' //exmod(modul)//'.fix'   
           fil_dynin_rad = adjustl(trim(namda))//'-dynamic-coupling.'//exmod(modul)//'.in'     
           fil_dynou_rad = adjustl(trim(namda))//'-dynamic-coupling.'//exmod(modul)//'.out'    
           fil_dynlo_rad = adjustl(trim(namda))//'-dynamic-coupling.'//exmod(modul)//'.log'    
           fil_dynre_rad = adjustl(trim(namda))//'-dynamic-coupling.'//exmod(modul)//'.res'    
           fil_splot_rad = adjustl(trim(namda))//'-specificheat.'    //exmod(modul)//'.fun'     
           fil_ramsh_rad = adjustl(trim(namda))//'-radiation.'       //exmod(modul)//'.post.msh'     
           fil_rares_rad = adjustl(trim(namda))//'-radiation.'       //exmod(modul)//'.post.res'    
        end if
        !
        ! Surface plot file
        !
        if(kfl_splot_rad==1) &
             call iofile(0_ip,lun_splot_rad,fil_splot_rad,'RADIAT SURFACE PLOT ',statu,forma,posit)
        !
        ! Boundary conditions
        !
        if(npp_bound_rad>0) &
             call iofile(0_ip,lun_bound_rad,fil_bound_rad,'RADIAT BOUND. COND. ',statu,forma,posit)
        !
        ! Matrix profile
        !
        if(kfl_psmat_rad>0) &
             call iofile(0_ip,lun_psmat_rad,fil_psmat_rad,'RADIAT MATRIX',statu,forma,posit)
        !
        ! Bc interpolation file
        !
        if(kfl_intbc_rad/=0) &
             call iofile(0_ip,lun_intbc_rad,fil_intbc_rad,'RADIAT BC INTERPOLATION','old')  
!!$        !
!!$        ! Dynamic coupling
!!$        !
!!$        if(kfl_dynco_rad/=0) then
!!$           call iofile(0_ip,lun_dynlo_rad,fil_dynlo_rad,'RADIAT DYNAMIC COUPLING LOG')  
!!$           call iofile(0_ip,lun_dynre_rad,fil_dynre_rad,'RADIAT DYNAMIC COUPLING RESULT')  
!!$        end if

     case(4_ip)
        !
        ! Close matrix file
        !
        call iofile(2_ip,lun_psmat_rad,' ','RADIAT MATRIX PROFILE')

     case(5_ip)
        !
        ! Interpolation file for k
        !
!!$        if (kfl_naked==0) then 
!!$           call GET_ENVIRONMENT_VARIABLE('FOR215',fil_funck_rad)       
!!$        else
!!$           fil_funck_rad = adjustl(trim(namda))//'-conductivity.'//exmod(modul)//'.fun'   
!!$        end if
!!$        wmate=intost(imate_rad)
!!$        fil_funck_rad=trim(fil_funck_rad)//trim(wmate)
!!$        call iofile(0_ip,lun_funck_rad,trim(fil_funck_rad),'RADIAT K INTERPOLATION','old')

     case(6_ip)
        !
        ! Interpolation file for Cp
        !
!!$        if (kfl_naked==0) then 
!!$           call GET_ENVIRONMENT_VARIABLE('FOR216',fil_funcc_rad)       
!!$        else
!!$           fil_funcc_rad = adjustl(trim(namda))//'-specificheat.'//exmod(modul)//'.fun'   
!!$        end if
!!$        wmate=intost(imate_rad)
!!$        fil_funcc_rad=trim(fil_funcc_rad)//trim(wmate)
!!$        call iofile(0_ip,lun_funcc_rad,trim(fil_funcc_rad),'RADIAT CP INTERPOLATION','old')

     end select

  end if

end subroutine rad_openfi

