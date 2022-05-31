!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_coefficients.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Coefficients
!> @details Here we calculate the source (rhs) term (including the scattered radiation from other directions), the current velocity (direction vector) and the total cross-section
!> @} 
!-----------------------------------------------------------------------

subroutine neu_coefficients(pnode,pgaus,gpsha,gpabs,gpsca,elrad,elunk,gpadv,gprea,gprhs,gpunk)

  use def_kintyp, only     :  ip,rp 
  use def_parame, only     :  pi,Stefan_Boltzmann,in4pi
  use def_domain, only     :  ndime
  use def_neutro, only     :  ADR_neu
  use def_neutro, only     :  num_energies_neu,num_directions_neu
  use def_neutro, only     :  current_energy_neu
  use def_neutro, only     :  current_direction_neu
  use def_neutro, only     :  direc_neu
  use def_neutro, only     :  weigd_neu
  use def_neutro, only     :  scattering_neu
  implicit none

  integer(ip), intent(in)  :: pnode 
  integer(ip), intent(in)  :: pgaus
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpabs(pgaus)
  real(rp),    intent(in)  :: gpsca(pgaus)
  real(rp),    intent(in)  :: elrad(num_energies_neu,num_directions_neu,pnode)
  real(rp),    intent(in)  :: elunk(pnode,ADR_neu % ntime)
  real(rp),    intent(out) :: gpadv(ndime,pnode)
  real(rp),    intent(out) :: gprea(pnode) 
  real(rp),    intent(out) :: gprhs(pnode)
  integer(ip)              :: inode,itime
  integer(ip)              :: igaus,idire,jdire
  real(rp)                 :: gprad(num_energies_neu,num_directions_neu) !> This variable will only be used to calculate the radiation scattered from other directions inthe RHS term
  real(rp)                 :: gpunk(pgaus,ADR_neu % ntime),gpsou(pgaus) !> Values for the unkown at the Gauss points for each time step
  
  !real(rp) :: tempe

  gpsou = 0.0_rp !> We keep the source in case we wanted it for the future
  idire = current_direction_neu !> We set the index for the direction into the current direction we consider
  gpunk = 0.0_rp 
  gprea = 0.0_rp

  do igaus = 1,pgaus !>Iterating over the Gauss points 
     !
     ! Source term
     !
     !tempe        = 0.0_rp
     !gpsou(igaus) = Stefan_Boltzmann*gpabs(igaus)*(tempe**4)/pi 

     gprad = 0.0_rp 
     do inode = 1,pnode !> Iterating over the number of nodes of the element considered
        gprad(1:num_energies_neu,1:num_directions_neu) = gprad(1:num_energies_neu,1:num_directions_neu) &
             + gpsha(inode,igaus) * elrad(1:num_energies_neu,1:num_directions_neu,inode) !> We add the contribution of each node of the element (radiation*shape at the Gauss point)
        do itime = 1,ADR_neu % ntime
           gpunk(igaus,itime) = gpunk(igaus,itime) + gpsha(inode,igaus) * elunk(inode,itime) !> We add the contribution of each element to the unkown for each time
        end do
     end do

     gpadv(1:ndime,igaus) = direc_neu(1:ndime,idire) !> The "velocity" is set equal to the direction we consider at the moment
     gprea(igaus)         = gpabs(igaus) + gpsca(igaus) * ( 1.0_rp - in4pi * weigd_neu(idire) * scattering_neu(idire,idire) ) !> We calculate the total cross section for the current direction MUST CHECK THIS
     gprhs(igaus)         = gpsou(igaus) !>The right-hand side of the equation is the source term

     do jdire = 1,num_directions_neu
        if( jdire /= idire ) then
           gprhs(igaus) = gprhs(igaus) + gpsca(igaus) * in4pi * weigd_neu(jdire) * scattering_neu(idire,jdire) * gprad(1,jdire) !> We add to RHS the contributions from scattered rad from other directions
        end if
     end do

  end do

end subroutine neu_coefficients
