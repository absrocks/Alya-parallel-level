!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_onetor.f90
!> @date    12/04/2013
!> @author  Zhinuo Jenny Wang
!> @brief   Initial condition setup for ToR-ORd 2019 heterogeneous model\n
!!          The init_toggle controls whether the initialisation is run as single cell steady-state\n
!!          or during a 2D/3D simulation\n
!> @}
!------------------------------------------------------------------------

subroutine exm_initor(ipoin, mat)

   use def_master
   use def_domain
   use def_elmtyp
   use def_exmedi 
   
   implicit none 
   integer(ip), intent(in) :: ipoin, mat
   integer(ip) :: ituss_exm
   
   ! Initialise single cell states
   vicel_exm(1:21,ipoin,1) = 0.0_rp ! [muA/muF]

   ! Initialise each node 
   if(kfl_timei_exm==1_ip) then  
      ! Initialise to values after 30 minutes pacing at cycle length of 1000 ms  
      ituss_exm = int(celty_exm(1,ipoin))

      elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)
      
      if(ituss_exm == EXM_CELLTYPE_EPI) then   !epicardial
         elmag(ipoin,1:3) = vminimate_exm(3,mat)
         vconc(1:14, ipoin, 1) = vcoin_exm(1:14,3,mat)
         vconc(1:14, ipoin, 2) = vcoin_exm(1:14,3,mat)
         vconc(1:14, ipoin, 3) = vcoin_exm(1:14,3,mat)
         vauxi_exm(1:31, ipoin, 1) = vauin_exm(1:31,3,mat) 
         vauxi_exm(1:31, ipoin, 2) = vauin_exm(1:31,3,mat)
         vauxi_exm(1:31, ipoin, 3) = vauin_exm(1:31,3,mat)
         
      else if (ituss_exm == EXM_CELLTYPE_MID) then ! mid-myocardial
         elmag(ipoin,1:3) = vminimate_exm(2,mat)
         vconc(1:14, ipoin, 1) = vcoin_exm(1:14,2,mat)
         vconc(1:14, ipoin, 2) = vcoin_exm(1:14,2,mat)
         vconc(1:14, ipoin, 3) = vcoin_exm(1:14,2,mat)
         vauxi_exm(1:31, ipoin, 1) = vauin_exm(1:31,2,mat)
         vauxi_exm(1:31, ipoin, 2) = vauin_exm(1:31,2,mat)
         vauxi_exm(1:31, ipoin, 3) = vauin_exm(1:31,2,mat)
      
      else if (ituss_exm == EXM_CELLTYPE_ENDO) then ! endocardial
         elmag(ipoin,1:3) = vminimate_exm(1,mat)
         vconc(1:14, ipoin, 1) = vcoin_exm(1:14,1,mat)
         vconc(1:14, ipoin, 2) = vcoin_exm(1:14,1,mat)
         vconc(1:14, ipoin, 3) = vcoin_exm(1:14,1,mat)
         vauxi_exm(1:31, ipoin, 1) = vauin_exm(1:31,1,mat)
         vauxi_exm(1:31, ipoin, 2) = vauin_exm(1:31,1,mat)
         vauxi_exm(1:31, ipoin, 3) = vauin_exm(1:31,1,mat)
         
      end if
   end if

end subroutine exm_initor

