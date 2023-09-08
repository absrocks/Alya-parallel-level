subroutine rad_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_inivar
  ! NAME 
  !    rad_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=0 ... When starting the run (from Turnon)
  !    ITASK=1 ... Solver initialization after some variables have been read
  ! USES
  ! USED BY
  !    rad_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  use def_solver
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(0_ip)  
     !
     ! Postprocess
     !
     postp(1) % wopos (1, 1) = 'RADAV'
     postp(1) % wopos (1, 2) = 'RADSO'
     postp(1) % wopos (1, 3) = 'ERROR'
     postp(1) % wopos (1, 4) = 'RESID'
     postp(1) % wopos (1, 5) = 'HEATF'
     postp(1) % wopos (1, 6) = 'TEMPR'

     postp(1) % wopos (2, 1) = 'SCALA'
     postp(1) % wopos (2, 2) = 'SCALA'
     postp(1) % wopos (2, 3) = 'SCALA'
     postp(1) % wopos (2, 4) = 'SCALA'
     postp(1) % wopos (2, 5) = 'VECTO'
     postp(1) % wopos (2, 6) = 'SCALA'
     !
     ! Set and witness variables 
     !
     postp(1) % woese(1)     = 'MEANG'
     postp(1) % wobse(1)     = 'MEANH'  
     postp(1) % wonse(1)     = 'RADAV'  
     postp(1) % wowit(1)     = 'RADAV'
     !
     ! Initialize general solver variables
     !     
     call soldef(-1_ip)
     solve(1) % wprob     = 'RADIATIONP1'          ! Equation name
     solve(1) % kfl_solve = 1                      ! Output flag

     !
     ! Check if we need our own # of species
     !
     if (kfl_modul(ID_CHEMIC)==1) then
        nspec_rad = nspec
     else
        nspec_rad=1
     end if

  case(1_ip)  

     ncomp_rad = 2     
     ittot_rad = 0                                  ! Total number of iterations
     !
     ! Solver 
     !
     solve_sol => solve(1:)
     if( solve_sol(1)%kfl_algso == 2 .or. solve_sol(1)%kfl_algso == 12 .or. solve_sol(1)%kfl_algso == 13 ) then
        if( INOTMASTER ) then
           solve_sol(1)%limpo => kfl_fixno_rad(1,1:)
        else
           allocate(solve_sol(1)%limpo(1))
           solve_sol(1)%limpo(1) = 0
        end if
        call cregro()
     end if

  case default

        call runend('RADIAT: Inivar called with wrong order')

  end select

end subroutine rad_inivar
