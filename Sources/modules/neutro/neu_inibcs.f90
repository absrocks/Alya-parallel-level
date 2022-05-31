 !-----------------------------------------------------------------------
!> @addtogroup NeutroTurnon
!> @ingroup Neutro
!> @{
!> @file    neu_inibcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Impose boundary conditions
!> @details Impose boundary conditions
!> @} 
!-----------------------------------------------------------------------
subroutine neu_inibcs()
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_neutro 
  implicit none
  integer(ip) :: ienergy,idirection,iunkn

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Allocate memory for the vectors needed to define the BC's 
     !
     !-------------------------------------------------------------
     
     call neu_membcs(1_ip) 

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then !if there are codes

        do idirection = 1,num_directions_neu
           do ienergy = 1,num_energies_neu

              kfl_fixno => kfl_fixno_neu(ienergy,idirection) % l !we associate the nodal fixity
              bvess     => bvess_neu(ienergy,idirection)     % a(:,:,1) ! Essential velocity bc values
              iunkn     =  (idirection-1)*num_energies_neu + ienergy !unkown identifier
              tncod     => tncod_neu(iunkn:) !node code type
              call reacod(IMPOSE_NODE_CODES)

           end do
        end do

     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then

        do idirection = 1,num_directions_neu
           do ienergy = 1,num_energies_neu

              kfl_fixbo => kfl_fixbo_neu(ienergy,idirection) % l
              bvnat     => bvnat_neu(ienergy,idirection)     % a(:,:,1)
              iunkn     =  (idirection-1)*num_energies_neu + ienergy
              tbcod     => tbcod_neu(iunkn:)
              call reacod(IMPOSE_BOUNDARY_CODES)

           end do
        end do

     end if

  end if

end subroutine neu_inibcs
