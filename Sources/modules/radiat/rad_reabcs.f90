subroutine rad_reabcs()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_reabcs
  ! NAME
  !    rad_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions for the radiation module
  !
  !    For now the only thing that needs is the emmisivity or reflectivity of the wall
  !   
  !    for open space, use emmisivity = 0 !!F ???
  !
  !    * For conditions on boundaries, bvnat(iboun) contains the emmisivity
  !
  !         Emissivity ... fixbo(iboun)=1 -> bvnat(1,iboun)
  !         Reflectivity ..... fixbo(iboun)=2 -> bvnat(1,iboun)=1-reflectivity
  !
  !    * For conditions on nodes, values are stored
  !      temporarily in bvess_rad(ipoin) and then changed to conditions on 
  !      boundaries  in the routine rad_bcntoe (called at the end) 
  !
  !    * Conditions on nodes have priority over conditions on boundaries.
  !      This is done and explained in rad_bcntoe. 
  !
  ! OUTPUT 
  ! USES
  ! USED BY
  !    rad_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_radiat
  use mod_opebcs
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip)  :: ifunc,ipara,ncodf,nbcod,dummi
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opnbcs(1_ip,1_ip,dummi,dummi,tncod_rad) ! Memory for structure
     call opnbcs(2_ip,1_ip, 1_ip, 0_ip,tncod_rad) ! Memory for variable
  end if
  npnat_rad = 1                                   ! # only parameter needed is emmisivity/reflectivity
  if( kfl_icodb > 0 ) then
     call opbbcs(0_ip,npnat_rad,1_ip,tbcod_rad)      
  end if
  !if( kfl_geome > 0 ) then
  !   call opnbcs(0_ip,1_ip,1_ip,0_ip,tgcod_rad)
  !end if

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     kfl_conbc_rad = 1             ! Constant boundary conditions
     kfl_inico_rad = 1             ! Default initial condition function: Black body radiation
     kfl_intbc_rad = 0             ! Do not interpolate !!F What
     !
     ! Read flags
     !
     iknbo_rad     =  0            ! Boundary condition is known
     ncodf         =  1            ! Radiat-P1 has 1 degree of freedom 
     nbcod         = -1

     call ecoute('rad_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('rad_reabcs')
     end do
     if(exists('NONCO')) then
        kfl_conbc_rad=0
     else
        kfl_conbc_rad=1
     end if
     if(exists('UNKNO')) then
        iknbo_rad=1
     end if
     if(exists('TIMEI')) then
        kfl_intbc_rad=2
     end if
     !
     ! Read data
     !
     call ecoute('rad_reabcs')
     do while( words(1) /= 'ENDBO' )

        if(words(1)=='CODES'.and.exists('NODES')) then 
           !
           ! User-defined codes on nodes
           !
           tncod => tncod_rad
           call reacod(1_ip)

        else if(words(1)=='CODES'.and.exists('BOUND')) then
           !
           ! User-defined codes on boundaries
           !          
           kfl_fixbo => kfl_fixbo_rad
           bvnat     => bvnat_rad(:,:,1)
           tbcod     => tbcod_rad(1:)
           call reacod(2_ip)

        else if(words(1)=='CODES') then
           call runend('RAD_REABCS: SPECIFY IF CODES ARE APPLIED TO NODES OR BOUNDARIES')

        else if(words(1)=='INITI') then 
           !
           ! Initial condition
           !
           if(words(2)=='CONST') then        ! Constant
              kfl_inico_rad = 2
              do ipara=1,10
                 bvcoe_rad(ipara)=param(1+ipara) ! First coeff is for bulk, rest is for boundaries
              end do
           else if(words(2)=='BLACK') then   ! Black body radiation
              kfl_inico_rad = 1
           end if

        else if(words(1)=='FUNCT') then
           !
           ! Functions
           !
           call ecoute('rad_reabcs')
           do while(words(1)/='ENDFU')
              if(kfl_conbc_rad==0) then
                 ifunc=getint('FUNCT',1_ip,'#FUNCTION NUMBER')
                 if(ifunc<0.or.ifunc>10) then
                    call runend('rad_reabcs: WRONG FUNCION NUMBER')
                 else
                    if(words(2)=='LINEA') then
                       kfl_funty_rad(ifunc)=1
                    else if(words(2)=='PARAB') then
                       kfl_funty_rad(ifunc)=2
                    else if(words(2)=='PERIO') then
                       kfl_funty_rad(ifunc)=3
                    end if
                    if(kfl_funty_rad(ifunc)>0) funpa_rad(1:6,ifunc)=param(3:8)
                 end if
              end if
              call ecoute('rad_reabcs')
           end do

        end if

        call ecoute('rad_reabcs')

     end do

  end if

end subroutine rad_reabcs
