!------------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_reanut.f90
!> @author  Guillaume Houzeaux
!> @date    February, 2017
!> @brief   Read numerical treatment
!>
!> @details Read numerical treatment
!>
!> @}
!------------------------------------------------------------------------
subroutine ale_reanut

 !.md<module>alefor
 !.md<input>case.ale.dat
 !.md<pos>1
 !.md<sec>

  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_alefor
  use def_domain
  use mod_ecoute,             only : ecoute
  use mod_read_domain_arrays, only : read_domain_arrays_types
  implicit none
  integer(ip) :: iboun,ipoin,itype,ktype,lexis_ad(nelty)

  if( INOTSLAVE ) then
     ! 
     !  Initializations (defaults)
     !
     kfl_smoot_ale     =  0                                      ! No mesh smoothing
     kfl_timef_ale     =  0                                      ! ALEFOR running coupled, with no timefunction
     nsmoo_ale         =  1                                      ! Number of smoothing steps (loading steps)
     kfl_defor_ale     =  0                                      ! No mesh deformation
     ndefo_ale         =  1                                      ! Number of deformation steps (loading steps)
     kfl_smobo_ale     =  0                                      ! Boundary smoothing
     kfl_fixsm_ale     =  1                                      ! Fix boundary nodes by default
     nsmob_ale         =  0                                      ! Number of boundary smoothing iterations
     kfl_crist_ale     = 0_ip                                    ! Cristobal 1 , new 0  ! RIGID BODY
     kfl_foexo_ale     = 2_ip                                    ! Force ( & Torque) extrapolation order  ! RIGID BODY
     kfl_disor_ale     = 2_ip                                    ! Integration order for the RB displacements  ! RIGID BODY
     kfl_nforc_ale     = 2_ip                                    ! Use the average of the previous 2 forces to calculate accel ! RB
     kfl_suppo_ale     =  0                                      ! Support geometry
     ansmo_ale         = -1.0_rp                                 ! Sharp edge detection
     resmo_ale         =  1.0_rp                                 ! Relaxation factor
     mnodb_ad          =  2
     !
     ! Reach the section
     !
     !.md<0># Numerical Treatment
     !.md<code>
     !.mdNUMERICAL_TREATMENT
     solve_sol => solve   
     call ecoute('ale_reanut')
     do while( words(1) /= 'NUMER' )
        call ecoute('ale_reanut')
     end do
     !
     ! Begin to read data
     !
     do while( words(1) /= 'ENDNU' )
        call ecoute('ale_reanut')

        if( words(1) == 'SMOOT' ) then
           !
           ! Mesh smoothing
           !
           if( exists('GAUSS') ) then
              kfl_smoot_ale = 2
           elseif( exists('CONST') ) then
              kfl_smoot_ale = 3
           else
              kfl_smoot_ale = 1
           end if
           if( exists('STEPS') ) nsmoo_ale = getint('STEPS',1_ip,  '#Number of smoothing steps') 
           if( exists('RELAX') ) resmo_ale = getrea('RELAX',1.0_rp,'#Smoothing relaxation factor') 

        else if( words(1) == 'DEFOR' ) then
           !
           ! Mesh deformation
           !
           if( exists('UNIFO') ) then
              kfl_defor_ale = 1
           else if( exists('SMALL') ) then
              kfl_defor_ale = 2
           else if( exists('ONLYB') ) then
              kfl_defor_ale = 4
           else if( exists('ISOTR') ) then
              kfl_defor_ale = 5
           else if( exists('BOUND') ) then
              kfl_defor_ale = 6
              if( exists('ALLEL') ) kfl_defor_ale = 7
              if( exists('TIMEF') ) kfl_timef_ale = 1
           else
              kfl_defor_ale = 1
           end if
           if( exists('STEPS') ) ndefo_ale = getint('STEPS',1_ip,'#Number of deformation steps') 

        else if( words(1) == 'BOUND' ) then
           !
           ! Boundary smoothing
           !
           kfl_smobo_ale = 1
           kfl_fixsm_ale = 0
           if( exists('STEPS') ) nsmob_ale = getint('STEPS',1_ip,'#Number of boundary smoothing steps') 
           if( exists('SHARP') ) ansmo_ale = getrea('SHARP',45.0_rp,'#Sharp edge angle') 

        else if( words(1) == 'FIXBO' ) then
           !
           ! Fix boundaries
           !
           if( words(2) == 'YES  ' .or. words(2) == 'ON   ' ) then
              kfl_fixsm_ale = 1
           else if( words(2) == 'NO   ' .or. words(2) == 'OFF  ' ) then
              kfl_fixsm_ale = 0
           end if

        else if( words(1) == 'ALGEB' ) then
           call reasol(1_ip)

        else if( words(1) == 'PRECO' ) then 
           call reasol(2_ip)

        else if( words(1) == 'RBTIM' ) then
           !
           ! Rigid Body time integration
           !
           if( words(2) == 'CRIST') then
              kfl_crist_ale = 1_ip
           else if( words(2) == 'SECON') then
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 2_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 2_ip                  ! Integration order for the RB displacements
           else if( words(2) == 'FIRST') then
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 1_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 1_ip                  ! Integration order for the RB displacements
           else if( words(2) == 'VSDF ') then
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 2_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 1_ip                  ! Integration order for the RB displacements
           else if( words(2) == 'VFDS ') then
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 1_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 2_ip                  ! Integration order for the RB displacements             
           else   ! default 2 2
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 2_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 2_ip                  ! Integration order for the RB displacements
           end if

        else if(words(1)=='FORCE') then
           if( words(2) == 'ORDIN' ) then
              kfl_nforc_ale = 1
           elseif( words(2) == 'AVERA' ) then
              kfl_nforc_ale = 2
           end if

        else if( words(1) == 'SUPPO' ) then
           !
           !.md<1>SUPPORT_GEOMETRY, NODES=int1, BOUNDARIES=int2                              $ Number of mesh multiplication levels. 0 means do not multiply
           !.md<2>COORDINATES
           !.md<3>...
           !.md<2>END_COORDINATES
           !.md<2>BOUNDARIES
           !.md<3>...
           !.md<2>END_BOUNDARIES
           !.md<1>END_SUPPORT_GEOMETRY   
           !.md<field>SUPPORT_GEOMETRY
           !
           npoin_ad      = getint('NODES',0_ip,'#NUMBER OF NODES OF THE SUPPORT GEOMETRY')
           nboun_ad      = getint('BOUND',0_ip,'#NUMBER OF NODES OF THE BOUNDARIES GEOMETRY')
           ktype         = 0
           kfl_suppo_ale = 1
           if( npoin_ad + nboun_ad == 0 ) call runend('ALE_REANUT: WRONG DIMENSIONS OF THE SUPPORT GEOMETRY')  
           call ale_membcs(5_ip)

           call ecoute('ale_reanut')
           do while( words(1) /= 'ENDSU' )

              if( words(1) == 'COORD' ) then
                 !
                 ! Coordinate COORD_AD(:,:)
                 !
                 call ecoute('ale_reanut')
                 do while( words(1) /= 'ENDCO')
                    ipoin = int(param(1),ip)
                    if( ipoin < 0 .or. ipoin > npoin_ad ) then
                       call runend('ALE_REANUT: WRONG SUPPORT GEOMETRY COORDINATES')
                    end if
                    coord_ad(1:ndime,ipoin) = param(2:1+ndime)
                    call ecoute('ale_reanut')
                 end do

              else if( words(1) == 'TYPES' ) then
                 !
                 ! Connectivity LTYPB_AD(:,:)
                 !
                 call read_domain_arrays_types(0_ip,nboun_ad,ktype,ltypb_ad,lexis_ad)

              else if( words(1) == 'BOUND' ) then
                 !
                 ! Connectivity LNODB_AD(:,:)
                 !
                 call ecoute('ale_reanut')
                 do while( words(1) /= 'ENDBO')
                    iboun = int(param(1))
                    if( iboun < 0 .or. iboun > nboun_ad ) then
                       call runend('ALE_REANUT: WRONG SUPPORT BOUNDARY CONNECTIVITY')
                    end if
                    lnodb_ad(1:mnodb_ad,iboun) = int(param(2:1+mnodb_ad),ip)
                    call ecoute('ale_reanut')
                 end do

              else if( words(1) == 'ALGEB' ) then
                 !
                 ! Algebraic solver
                 !                       
                 solve_sol => solve(3:)
                 call reasol(1_ip)

              end if

              call ecoute('ale_reanut')
           end do

        end if

     !.mdEND_NUMERICAL_TREATMENT
     !.md</code>

     end do
  end if

end subroutine ale_reanut
