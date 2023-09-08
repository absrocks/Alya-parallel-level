subroutine nsa_partim
!-----------------------------------------------------------------------
!****f* Nastal/nsa_partim
! NAME 
!    nsa_inivar
! DESCRIPTION
!    This routine sets the time scheme parameters
!    
!
!
! USES
! USED BY
!    nsa_inivar
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  implicit none

  parkb_nsa      = 0.0_rp
  parka_nsa      = 0.0_rp

!!  if (kfl_timul_nsa == 1) then

     if (kfl_tisch_nsa == 30) then
        !
        !   heun 3-stages, O(3)
        !
        !   parka is the A matrix
        !   parkb is the b vector
        !
        !   parka and parkb last indices label the fractional level
        !
        parka_nsa(1,1,1) = 1.0_rp/3.0_rp
        parka_nsa(1,2,1) = 0.0_rp
        parka_nsa(2,2,1) = 2.0_rp/3.0_rp

        parkb_nsa(1,1) = 1.0_rp / 4.0_rp   
        parkb_nsa(2,1) = 0.0_rp
        parkb_nsa(3,1) = 3.0_rp / 4.0_rp 

     else if (kfl_tisch_nsa == 31) then
        !
        ! explicit crank-nicholson with fixed point, 3-stages, O(2)
        !
        parka_nsa(1,1,1) = 1.0_rp
        parka_nsa(1,2,1) = 1.0_rp/2.0_rp        
        parka_nsa(2,2,1) = 1.0_rp/2.0_rp        

        parkb_nsa(1,1) = 1.0_rp / 2.0_rp   
        parkb_nsa(2,1) = 0.0_rp
        parkb_nsa(3,1) = 1.0_rp / 2.0_rp 

     else if (kfl_tisch_nsa == 41) then
        !
        ! explicit Runge-Kutta 3-stages
        !
        parka_nsa(1,1,1) = 1.0_rp
        parka_nsa(2,1,1) = 0.75_rp
        parka_nsa(2,2,1) = 0.25_rp
        parka_nsa(3,1,1) = 0.333333333333_rp
        parka_nsa(3,3,1) = 0.666666666666667_rp

        parkb_nsa(1,1) = 1.0_rp
        parkb_nsa(2,1) = 0.25_rp
        parkb_nsa(3,1) = 0.666666666666667_rp
     end if

!!!  end if

  
  safet_nsa = safet_nsa * safrk_nsa


end subroutine nsa_partim
