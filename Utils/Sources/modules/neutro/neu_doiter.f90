!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_doiter.f90
!> @date    29/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Solve inner iterations
!> @details Solve inner iterations
!> @}
!------------------------------------------------------------------------

subroutine neu_doiter()

  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_neutro

  implicit none

  call neu_begite()

  do while( kfl_goite_neu == 1 )

     itinn(modul) = itinn(modul) + 1

     do current_energy_neu = 1,num_energies_neu
        do current_direction_neu = 1,num_directions_neu
           if( INOTSLAVE ) print*,'Solving=',current_energy_neu,current_direction_neu
           call neu_solite()
           call neu_endite(1_ip)
        end do
     end do
     call neu_endite(3_ip) ! Check global convergence

  end do

  call neu_endite(2_ip)

end subroutine neu_doiter

