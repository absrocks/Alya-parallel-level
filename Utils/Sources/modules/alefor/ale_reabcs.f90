subroutine ale_reabcs()
  !------------------------------------------------------------------------
  !****f* Alefor/ale_reabcs
  ! NAME 
  !    ale_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions 
  ! USES
  !    ale_bcntoe
  !    ecoute
  !    memchk
  !    runend
  !    ale_bounod
  !    ale_autbcs
  ! USED BY
  !    ale_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_alefor 
  use def_kermod, only: number_space_time_function
  use mod_memchk
  use mod_opebcs
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: ifunc,dummi
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opnbcs(1_ip,1_ip,dummi,dummi,tncod_ale) ! Memory for structure
     call opnbcs(2_ip,1_ip,ndime, 0_ip,tncod_ale) ! Memory for variable
 end if
  if( kfl_icodb > 0 ) then
     call opbbcs(0_ip,1_ip,1_ip,tbcod_ale)      
  end if
  if( kfl_geome > 0 ) then
     call opnbcs(0_ip,1_ip,ndime, 0_ip,tgcod_ale)
  end if

  if( INOTSLAVE ) then
     !
     ! Initialization global variables
     !
     kfl_conbc_ale = 1                                       ! Constant boundary conditions
     if (number_space_time_function > 0) kfl_conbc_ale = 0   ! Non-constant boundary conditions coming from kernel functions
     !
     ! Reach the nodal-wise section
     !
     call ecoute('ale_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('ale_reabcs')
     end do
     !
     ! Loop over nodes and or boundaries
     !
     call ecoute('ale_reabcs')
     do while(words(1)/='ENDBO')

        if( words(1) == 'CODES' .and. exists('NODES') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on nodes
           !
           !-------------------------------------------------------------

           if(exists('GEOME')) then
              !
              ! Geometrical node code
              !              
              tgcod => tgcod_ale(1:)
              call reacod(4_ip)

           else
              !
              ! Velocity: node codes
              !
              tncod => tncod_ale(1:)
              call reacod(1_ip)

           end if

        else if( words(1) == 'CODES' .and. exists('BOUND') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on boundaries
           !          
           !-------------------------------------------------------------

           tbcod => tbcod_ale(1:)
           call reacod(2_ip)

        else if( words(1) == 'FUNCT' ) then

           !-------------------------------------------------------------
           !
           ! Function definitions
           !
           !-------------------------------------------------------------

           kfl_conbc_ale = 0
           call ale_membcs(47_ip)  

           call ecoute('ale_reabcs')
           do while( words(1) /= 'ENDFU' )
              
              if( words(1) == 'FUNCT' ) then

                 do while( words(1) /= 'ENDFU' )
                    ifunc = getint('FUNCT',1_ip,'#FUNCTION NUMBER')

                    if( ifunc < 0 .or. ifunc > 10 ) then
                       
                       call runend('readat: WRONG FUNCTION NUMBER')
                       
                    else
                       
                       if( words(2) == 'TRANS' ) then
                          kfl_funty_ale(ifunc,1) =  6
                          kfl_funty_ale(ifunc,2) =  4
                       else if( words(2) == 'ROTAT' ) then
                          kfl_funty_ale(ifunc,1) =  7
                          kfl_funty_ale(ifunc,2) =  4
                        else if( words(2) == 'SIN  ' ) then
                          kfl_funty_ale(ifunc,1) =  8
                          kfl_funty_ale(ifunc,2) =  4
                        else if( words(2) == 'LINEA' ) then
                          kfl_funty_ale(ifunc,1) =  9
                          kfl_funty_ale(ifunc,2) =  3
                       else
                          kfl_funty_ale(ifunc,1) =  0 
                       end if
                       
                       if( kfl_funty_ale(ifunc,1) > 0 ) then
                          igene = ifunc
                          call ale_membcs(46_ip)  
                          funpa_ale(ifunc) % a(1:kfl_funty_ale(ifunc,2)) = param(3:3+kfl_funty_ale(ifunc,2)-1)
                       end if
                    end if
                    call ecoute('READAT')
                          
                 end do
              end if
              call ecoute('READAT')
                          
           end do

        end if
        call ecoute('ale_reabcs')
     end do

  end if

end subroutine ale_reabcs
