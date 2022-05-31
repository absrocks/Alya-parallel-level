subroutine nsa_openfi(itask)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_openfi
  ! NAME 
  !    nsa_openfi
  ! DESCRIPTION
  !    This subroutine gets ALL the file names and open them to be used by 
  !    the module in two possible ways:
  ! 
  !    1. Recalling them from the environment, when Alya is launched
  !       encapsulated in a shell script, or
  ! 
  !    2. Composing the names out of the problem name which is given as
  !       argument when the binary file Alya is launched "naked".
  ! USES
  !    iofile
  ! USED BY
  !    nsa_turnon
  !***
  !-----------------------------------------------------------------------
  use def_nastal
  use def_parame
  use def_master
  use def_domain
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ilcha 
  character(150)          :: fil_pdata_nsa,fil_conve_nsa,fil_force_nsa
  character(150)          :: fil_outpu_nsa,fil_solve_nsa,fil_maxmi_nsa
  character(150)          :: fil_rstar_nsa,fil_pro2d_nsa,fil_dumb1_nsa,fil_dumb2_nsa
  character(150)          :: fil_setse_nsa,fil_setsb_nsa,fil_setsn_nsa

  if(kfl_paral<=0) then

     select case(itask)

     case(1)
        !
        ! kfl_naked is set in the kernel subroutine getnam
        !
        if (kfl_naked==0) then
           !  encapsulated, then get names from the environment
           call GET_ENVIRONMENT_VARIABLE('FOR606',fil_force_nsa)
           call GET_ENVIRONMENT_VARIABLE('FOR607',fil_chkpo_nsa(1))     ! default in
           call GET_ENVIRONMENT_VARIABLE('FOR608',fil_chkpo_nsa(2))     ! default out
           call GET_ENVIRONMENT_VARIABLE('FOR609',fil_pro2d_nsa)
           call GET_ENVIRONMENT_VARIABLE('FOR610',fil_maxmi_nsa)
           call GET_ENVIRONMENT_VARIABLE('FOR620',fil_dumb1_nsa)
           call GET_ENVIRONMENT_VARIABLE('FOR621',fil_dumb2_nsa)
        else if (kfl_naked==1) then
           !  naked, then compose the names     
           fil_force_nsa = adjustl(trim(namda))//'.'//exmod(modul)//'.frc'
           fil_chkpo_nsa(1) = adjustl(trim(namda))//'.'//exmod(modul)//'.chk.in'        ! default in     
           fil_chkpo_nsa(2) = adjustl(trim(namda))//'.'//exmod(modul)//'.chk.out'       ! default out
           if (ndime==2) &
                fil_pro2d_nsa = adjustl(trim(namda))//'.'//exmod(modul)//'.p2d'     
           fil_maxmi_nsa = adjustl(trim(namda))//'.'//exmod(modul)//'.mxm'
           fil_dumb1_nsa = adjustl(trim(namda))//'.'//exmod(modul)//'-fixboundary.fix'
           fil_dumb2_nsa = adjustl(trim(namda))//'.'//exmod(modul)//'-fixboundary.geo'
        end if
        !
        ! Open files
        !
        call iofile(zero,lun_force_nsa,fil_force_nsa,'NASTAL FORCES AND COEFFICIENTS')
        if (ndime==2) &
             call iofile(zero,lun_pro2d_nsa,fil_pro2d_nsa,'NASTAL 2D-PROFILE DIST')
        call iofile(zero,lun_maxmi_nsa,fil_maxmi_nsa,'NASTAL MIN-MAX')

     case(2)
        !
        ! Dump derived boundary conditions
        !
        if (kfl_dumbo_nsa > 0) then
           call iofile(zero,lun_dumb1_nsa,fil_dumb1_nsa,'DERIVED BOUNDARY CONDITIONS')
           call iofile(zero,lun_dumb2_nsa,fil_dumb2_nsa,'DERIVED BOUNDARIES ')           
        end if

     case(3)

     case(4)
        !
        ! Close output files
        !
        call iofile(two,lun_maxmi_nsa,' ','NASTAL MAXMIN')

     end select

  end if

end subroutine nsa_openfi

