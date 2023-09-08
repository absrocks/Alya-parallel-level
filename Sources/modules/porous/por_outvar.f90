!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_outvar.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Output a postprocess variable
!> @details Output a postprocess variable
!> @} 
!------------------------------------------------------------------------
subroutine por_outvar(ivari)
  use def_parame
  use def_master
  use def_domain
  use def_porous
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: ipoin,jpoin,iline,icont,ielem,inode

  select case (ivari)  

  case(0_ip)
     !
     ! Do nothing
     !
     return

  case(1_ip)
     !
     ! Velocity  - I guess I must be obtained after solving for pressure - for the moment a simple projection 
     !
     gevec => veloc(:,:,1)

  case(2_ip)
     !
     ! Pressure
     !
     gesca => press(:,1) 

  case(3_ip)
     !
     ! Water Saturation
     !
     gesca => wasat(:,1) 

  case(4_ip)
     !
     ! IJKNO
     !
     if( INOTMASTER ) then
        do ipoin=1,npoin
           ! rhsid(ipoin) = real(ijkno(ipoin))   ! fo the moment I still do not have a variable for this
           rhsid(ipoin) = 0.0_rp
        end do
        gesca => rhsid
     end if

  case(5_ip)
     !
     ! Porosity
     !
     gesca => nodpo_por  ! uses smoothed value

  case(6_ip)
     !
     ! Wells
     !
     if( INOTMASTER ) then
        call memgen(0_ip,npoin,0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = real(iwell_por(ipoin),rp)
        end do
     end if

  case(7_ip)
     !
     ! Linelets 
     !
     if( INOTMASTER ) then
        icont = 0
        do ipoin = 1,npoin
           rhsid(ipoin) = 0.0_rp
        end do 
        do iline = 1,solve(1) % nline
           icont = icont+1
           do ipoin = solve(1) % lline(iline),solve(1) % lline(iline+1)-1
              jpoin = solve(1) % lrenup(ipoin)
              rhsid(jpoin) = real(icont,rp)
           end do
        end do
        gesca => rhsid
     end if

  case(8_ip)
     !
     ! GROUP: GROUPS FOR DEFLATED pressure
     !
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           rhsid(ipoin) = real(solve(1) % lgrou(ipoin),rp)
        end do
        gesca => rhsid
     end if

  case(9_ip)
     !
     ! Permeability
     !
     gevec => nodpe_por  ! uses smoothed value

  end select

  call outvar(&
       ivari,&
       ittim,cutim,postp(1) % wopos(1,ivari))

end subroutine por_outvar
