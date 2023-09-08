subroutine got_openfi(itask)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_openfi
  ! NAME 
  !    got_openfi
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
  !    got_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_gotita
  use      def_master
  use      def_postpr
  use      mod_iofile
  implicit none
  integer(ip), intent(in) :: itask  
  character(150)          :: fil_pdata_got,fil_outpu_got
  character(150)          :: fil_solve
  character(150)          :: fil_setse_got,fil_setsb_got
  character(150)          :: fil_setsn_got,fil_bound_got
  character(150)          :: fil_cvgso_got,fil_stasg_got
  character(150)          :: fil_cvgsg_got
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if(kfl_paral<=0) then
     !
     ! Define unit opening option if this is a restart run
     !
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

     case (2)
        !
        ! Open files needed occasionally
        !
        if(kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR1110',fil_bound_got)     
           call GET_ENVIRONMENT_VARIABLE('FOR1112',fil_stasg_got)
           call GET_ENVIRONMENT_VARIABLE('FOR1113',fil_cvgsg_got)
        else
           fil_bound_got = adjustl(trim(namda))//'.'         //exmod(modul)//'.bcs'
           fil_stasg_got = adjustl(trim(namda))//'.'         //exmod(modul)//'.sgs'
           fil_cvgsg_got = adjustl(trim(namda))//'.'         //exmod(modul)//'.csg'
        end if
        !
        ! Boundary conditions
        !
        if(npp_bound_got==1) &
             call iofile(0_ip,lun_bound_got,fil_bound_got,'GOTITA BOUND. COND',  statu,forma,posit)
        !
        ! Subgrid scales
        !
        if(kfl_sgsco_got==1) then
           call iofile(0_ip,lun_stasg_got,fil_stasg_got,'GOTITA SUBGRID SCALE STAT.',statu,forma,posit)
           call iofile(0_ip,lun_cvgsg_got,fil_cvgsg_got,'GOTITA SUBGRID SCALE CONV.',statu,forma,posit)
        end if

     end select

  end if

end subroutine got_openfi

