subroutine qua_turnon()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_turnon
  ! NAME 
  !    qua_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the temperature equation.
  !    - Write some info
  !    - Check if there are errrors
  !    - Allocate memory
  ! USES
  !    qua_inivar
  !    qua_openfi
  !    qua_reaphy
  !    qua_reabcs
  !    qua_reanut
  !    qua_reaous
  !    qua_outinf
  !    qua_memall
  !    qua_restar
  ! USED BY
  !    Quanty
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_quanty
  use      mod_iofile
  implicit none
  !
  ! Initial variables
  !
  call qua_inivar(zero) ! :o) 
  !
  ! Open files
  !
  call qua_openfi(one)  ! :o)
  !
  ! Initial variables
  !
  !call qua_inivar(zero) ! :o)
  !
  ! Read the physical problem
  !
  call qua_reaphy()     !  (?):|
  !
  ! Read the numerical treatment
  !
  call qua_reanut()     !  (?):|
  !
  ! Read the output strategy
  !
  call qua_reaous()     !  (?):|
  !
  ! Parall service
  !
  call qua_parall(one)
  !
  ! Read the boundary conditions
  !
  call qua_reabcs()       ! (?) :| 
  !
  ! If a deflated solver is used
  !
  call qua_inivar(5_ip)
  !
  ! Parall service
  !
  call qua_parall(two)
  !
  ! Initial variables
  !
  call qua_inivar(one)     ! :o)
  !
  ! Write info
  !
  call qua_outinf()        ! :o) 
  !
  ! Warnings and errors
  !
  call qua_outerr()        ! :o) 
  !
  ! Allocate memory
  !
  call qua_memall()        ! :| (ojo!!!?)  
  !
  ! Read restart file
  !
  call qua_restar(one)     ! :o) 
  !
  ! Open additional files
  !
  call qua_openfi(two)     ! :o) 
  !
  ! Close input data files
  !
  call qua_openfi(three)   ! :o) 

end subroutine qua_turnon
