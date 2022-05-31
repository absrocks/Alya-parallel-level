!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_outvar.f90
!> @date    01/04/2016
!> @author  Guillaume Houzeaux
!> @brief   Postprocess
!> @details Postprocess
!> @} 
!-----------------------------------------------------------------------

subroutine neu_outvar(ivari)

  use def_master
  use def_domain
  use def_neutro
  use mod_postpr
  use mod_outvar,         only : outvar
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: idire,idime,ipoin
  character(5)            :: wopos(3)
  character(20)           :: wdire

  select case ( ivari )  

  case ( 0_ip )
     !
     ! Do nothing
     !
     return

  case ( 1_ip )
     !
     ! CURRE: neutron current
     !
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)!> Allocate memory for a vector array of length ndime and element dimension npoin which will be allocated in gevec (for a scalar it would be gesca)
        do ipoin = 1,npoin  !> Iterating over the points of the mesh
           do idire = 1,num_directions_neu !> Iterating over the directions in each point
              do idime = 1,ndime !> Iterating over the spatial dimensions
                 gevec(idime,ipoin) = gevec(idime,ipoin) + neutr(1,idire,ipoin,1) * direc_neu(idime,idire) * weigd_neu(idire) !> The current is the vector sum of the flux*direction*weight
              end do
           end do
        end do
     end if

  case ( 2_ip )
     !
     ! FLUX: neutron flux
     !
     if( INOTMASTER ) then
        call memgen(0_ip,npoin,0_ip) !> Allocate memory for an npoin scalar array in gesca
        do ipoin = 1,npoin  
           do idire = 1,num_directions_neu
              gesca(ipoin) = gesca(ipoin) + neutr(1,idire,ipoin,1) * weigd_neu(idire) !> We add the contributions for each direction (scalarly) multiplying by the weight of that direction
           end do
        end do
     end if

  case ( 3_ip )
     !
     ! NEUTRONS: Neutrons (what is this??)
     !
     wopos(2:3) = postp(1) % wopos(2:3,3)

     if( INOTMASTER ) call memgen(0_ip,npoin,0_ip)
     do idire = 1,num_directions_neu
        wdire = adjustl(intost(idire))
        if( idire < 10 ) then
           wopos(1) = postp(1) % wopos(1,3)(1:3) // '0' // trim(wdire)
        else
           wopos(1) = postp(1) % wopos(1,3)(1:3) // trim(wdire)
        end if
        if( INOTMASTER ) gesca(1:npoin) = neutr(1,idire,1:npoin,1)
        call postpr(gesca,wopos,ittim,cutim)
     end do
     if( INOTMASTER ) call memgen(2_ip,npoin,0_ip)
     return

  end select
  !
  ! Output GESCA for a scalar or GEVEC for a vector
  !
  call outvar(ivari,ittim,cutim,postp(1) % wopos(:,ivari))

end subroutine neu_outvar
