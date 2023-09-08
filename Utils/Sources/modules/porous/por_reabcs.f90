!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_outvar.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Reads the boundary conditions for porous equations.
!> @details Reads the boundary conditions for porous equations. For the moment we do not need any bcs for porous.
!> @} 
!------------------------------------------------------------------------
subroutine por_reabcs()
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_porous
  use mod_opebcs
  implicit none
  integer(ip) :: ii ,dummi
  !
  ! Allocate memory
  ! 
  if( kfl_icodn > 0 ) then
     call opnbcs(1_ip,2_ip,dummi,dummi,tncod_por)       ! Memory for structure
     do ii = 1,2
        call opnbcs(2_ip,ii,1_ip, 0_ip,tncod_por)       ! Memory for variable
     end do
  end if

  if( INOTSLAVE ) then

     call ecoute('tem_reabcs')
     do while( words(1) /= 'ENDBO' )

        if( words(1) == 'CODES' .and. exists('NODES') .and. kfl_icodn > 0 ) then 
           !
           ! User-defined codes on nodes
           !
           if( exists('VARIA') ) then
              ii     = getint('VARIA',1_ip,'*USER BOUNDARY CONDITIONS')
              tncod => tncod_por(ii:)
              call reacod(1_ip)
           else
              tncod => tncod_por(1:)
              call reacod(1_ip)
              do ii = 2,2
                 call cpybcs(1_ip,ii,tncod_por) 
              end do
           end if

        end if

        call ecoute('tem_reabcs')

     end do

  end if

end subroutine por_reabcs
