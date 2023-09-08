subroutine dod_openfi(itask)
  !-----------------------------------------------------------------------
  !****f* Dodeme/dod_openfi
  ! NAME 
  !    dod_openfi
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
  !    dod_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_dodeme
  use      def_master
  use      def_postpr
  use      mod_iofile
  implicit none
  integer(ip), intent(in) :: itask  
  character(150)          :: fil_pdata_dod,fil_outpu_dod
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

     case (1)
        !
        ! Open files that are always needed
        !
        if(kfl_naked==0) then
           !
           !  encapsulated, then get names from the environment
           !
           call GET_ENVIRONMENT_VARIABLE('FOR5701',fil_pdata_dod)    
           call GET_ENVIRONMENT_VARIABLE('FOR5702',fil_outpu_dod)    

        else if(kfl_naked==1) then
           !
           !  naked, then compose the names     
           !
           fil_pdata_dod = adjustl(trim(namda))//'.'//exser(servi)//'.dat'    
           fil_outpu_dod = adjustl(trim(namda))//'.'//exser(servi)//'.log'    
    
        end if

        call iofile(zero,lun_pdata_dod,fil_pdata_dod,'DODEME DATA       ','old')
        call iofile(zero,lun_outpu_dod,fil_outpu_dod,'DODEME OUTPUT     ',statu,forma,posit)

     case (2)
        !
        ! Close data file
        !       
        call iofile(two,lun_pdata_dod,' ','DODEME DATA')

     case (3)
        !
        ! Close result file
        !       
        call iofile(two,lun_outpu_dod,' ','DODEME OUTPUT')

     end select

  end if

end subroutine dod_openfi

