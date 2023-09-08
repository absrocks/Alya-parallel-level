subroutine hlm_reabcs()

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_reabcs.f90
  ! NAME
  !    hlm_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions. 
  ! USES
  ! USED BY
  !    hlm_turnon
  !-----------------------------------------------------------------------

  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_helmoz
  use mod_opebcs
  use mod_ecoute, only : ecoute
  implicit none

  integer(ip)  :: ipoin,pnodb,iboun,inodb,ifunc,ipara,ibsta,knodb(mnodb)
  integer(ip)  :: pblty,ncodf,nbcod,dummi
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opnbcs(1_ip,1_ip,dummi,dummi,tncod_hlm)       ! Memory for structure
     call opnbcs(2_ip,1_ip, 1_ip, 0_ip,tncod_hlm)       ! Memory for variable
  end if
  if( kfl_icodb > 0 ) then
     call opnbcs(1_ip,1_ip,dummi,dummi,tecod_hlm)       ! Memory for structure
     call opnbcs(2_ip,1_ip, 1_ip, 0_ip,tecod_hlm)       ! Memory for variable
  end if

  if( INOTSLAVE ) then
     !
     ! Reach section
     !
     call ecoute('hlm_reabcs')
     do while (words(1) /= 'BOUND')
        call ecoute('hlm_reabcs')
     enddo
     !
     ! Read data - codes on the nodes
     !
     call ecoute('hlm_reabcs')

     do while (words(1) /= 'ENDBO')

        if (     words(1) == 'CODES' .and. exists('NODES') ) then 
           !
           ! User-defined codes on nodes
           !
           tncod => tncod_hlm
           call reacod(READ_NODE_CODES)

        else if( words(1) == 'CODES' .and. exists('EDGES') ) then 
           !
           ! User-defined codes on edges
           !
           tncod => tecod_hlm
           call reacod(READ_EDGE_CODES)

        end if

        call ecoute('hlm_reabcs')
     end do

  end if

end subroutine hlm_reabcs
