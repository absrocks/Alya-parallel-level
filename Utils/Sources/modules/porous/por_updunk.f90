!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_updunk.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Performs several types of updates for the porous equation
!> @details Performs several types of updates for the porous equation
!> @} 
!------------------------------------------------------------------------
subroutine por_updunk(itask)
  use def_parame
  use def_master
  use def_domain
  use def_porous
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin


  if( INOTMASTER ) then

     select case (itask)

     case(1_ip)
        !
        ! Assign P(n,0,*) <-- P(n-1,*,*), initial guess for outer iterations - also S
        !
        do ipoin=1,npoin
           wasat(ipoin,2) = wasat(ipoin,ncomp_por)
           press(ipoin,2) = press(ipoin,ncomp_por)
        end do

     case(2_ip)
        !
        ! Assign P(n,i,0) <-- P(n,i-1,*), initial guess for inner iterations - also S
        !
        do ipoin=1,npoin
           wasat(ipoin,1) = wasat(ipoin,2)
           press(ipoin,1) = press(ipoin,2)
        end do

     case(3_ip)
        !
        ! Assign P(n,i,j-1) <-- P(n,i,j), update of the pressure
        !
        do ipoin=1,npoin
           press(ipoin,1) = unkno(ipoin)
        end do

     case(4_ip)
        !
        ! Assign P(n,i-1,*) <-- P(n,i,*) - also S
        !        
        do ipoin=1,npoin
           wasat(ipoin,2) = wasat(ipoin,1)
           press(ipoin,2) = press(ipoin,1)
        end do

     case(5_ip)
        !
        ! Obtain P(n,*,*) for the Crank-Nicolson method and assign
        ! P(n-1,*,*) <-- P(n,*,*)
        !        
        if(kfl_tisch_por==1.and.kfl_tiacc_por(1)==2) then
           !
           ! Crank-Nicolson method 
           !
           call runend('por_updunk: for the moment CN not ready')        
!           do ipoin=1,npoin
!              press(ipoin,1) = 2.0_rp*press(ipoin,1)-press(ipoin,3)
!           end do

        else if(kfl_tisch_por==2) then
           !
           ! BDF scheme
           !
           call runend('por_updunk: for the moment CN not ready')        
!           do ipoin=1,npoin
!              do itime=2+kfl_tiacc_por(1),4,-1
!                 press(ipoin,itime) = press(ipoin,itime-1)
!              end do
!           end do
        end if
!        call por_averag()         !not ready
        do ipoin=1,npoin
           wasat(ipoin,3) = wasat(ipoin,1)
           press(ipoin,3) = press(ipoin,1)
        end do

        if(kfl_sgsti_por(1)==1) then
           !
           ! Time tracking of the subscales
           !        
           call runend('por_updunk:tracking of the subscales not ready')

        end if

     case(6_ip) 
        !
        ! Assign P(n,1,*) <-- P(n-1,*,*), when readin from restart file
        ! 
        call runend('por_updunk: for moment restart not ready')        

!        icomp=min(3_ip,ncomp_por) 
!        do ipoin=1,npoin
!           press(ipoin,1) = press(ipoin,icomp)
!        end do

     case(7_ip)
        !
        ! Assign S(n,i,j-1) <-- S(n,i,j), update of the pressure
        !
        do ipoin=1,npoin
           wasat(ipoin,1) = unkno(ipoin)
        end do

     case(8_ip)
        !
        ! Initial guess for solver - water saturartion
        !
        do ipoin=1,npoin
           unkno(ipoin)   = wasat(ipoin,1)
        end do

     case(9_ip)
        !
        ! Initial guess for solver - pressure
        !
        do ipoin=1,npoin
           unkno(ipoin)   = press(ipoin,1)
        end do

     end select

  end if

end subroutine por_updunk

