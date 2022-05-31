!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_torord_ionicurrents.f90
!> @author  Zhinuo Jenny Wang
!> @brief   Wrapper for ToR-ORd 2020 model in 3D
!> @date    14/FEB/2020
!> @details 
!> @}
!!-----------------------------------------------------------------------
subroutine exm_torord_ionicurrents(ipoin, xioni, dioni, cai)

   use      def_parame
   use      def_master
   use      def_elmtyp
   use      def_domain
   use      def_exmedi
   use      mod_exm_torord_model, only: exm_torord_model

   ! Definition of variables
   implicit none
   integer(ip), intent(in) :: ipoin !< node
   real(rp), intent(out) :: xioni   !< current
   real(rp), intent(out)   :: dioni !< current derivative
   real(rp), intent(out)  :: cai    !< Intracellular calcium concentration
   real(rp) :: a2bas                !< Apex-to-base scaling factor for Gks
   real(rp) :: dtimeEP              !< Time step for solving cell model
   real(rp) :: gkr_scaling          !< Possible GKr scaling for BZ
   real(rp) :: drugd(24)             !< Drug effects scaling
   integer(ip) :: i, n,imate,ituss_exm
   logical :: flag_land,flag_3D,flag_border,flag_drug
   real(rp) :: statvar(7,2)

   if (INOTMASTER) then

      ! Get apex to base gradient of gKs if it is defined
      if (kfl_atbhe_exm == 0_ip) then
         a2bas = 1.0_rp
      else
         a2bas = atbhe_exm(1,ipoin)
      end if

      if (kfl_cellmod(nodemat(ipoin)) == CELL_TORORD_EXMEDI) then

         ituss_exm = int(celty_exm(1,ipoin))

         n = nodemat(ipoin)

         if (n==0) then
            call runend('Elements with material 0 are present')
         end if

         dtimeEP = dtime * 1000.0_rp ! Convert to [ms]

         ! Provide stimulus
         vicel_exm(20,ipoin,1) = appfi_exm(ipoin)
         flag_land = .FALSE.
         do imate = 1,nmate_exm
            if (kfl_eccty(imate) == 4_ip) flag_land=.TRUE.
         end do
         flag_3D = .TRUE.

         !if (kfl_borde_exm(n)==1) then
         !   flag_border = .TRUE.
         !   gkr_scaling = border_gkrsc_exm(n)
         !else
            flag_border = .FALSE.
            gkr_scaling = 1.0_rp
         !end if

         if (kfl_drugsmate_exm(n) == 1) then
            flag_drug = .TRUE.
            drugd(:) = drugdmate_exm(:,n)
         else
            flag_drug = .FALSE.
            drugd = 1.0_rp
         end if

         statvar(:,:) = 0.0_rp

         call exm_torord_model(ituss_exm, elmag(ipoin,ITER_K),vconc(:,ipoin,:),vauxi_exm(:,ipoin,:),vicel_exm(:,ipoin,1),dtimeEP, a2bas, &
           & qneto_exm(ipoin),flag_land,flag_3D,flag_border,flag_drug,gkr_scaling,drugd,statvar,ipoin)

         dioni = 0.0_rp
         xioni = 0.0_rp
         cai = vconc(5,ipoin,2) ! Output intracellular calcium concentration
         do i=1,19
            xioni = xioni + vicel_exm(i,ipoin,1) ! Istim is added in mod_exm_ionicurrents
         end do
         xioni = xioni + vicel_exm(21,ipoin,1) ! Add stretch-activated current

         vconc(:,ipoin,1) = vconc(:,ipoin,2)
         vauxi_exm(:,ipoin,1) = vauxi_exm(:,ipoin,2)
      end if
   end if

end subroutine exm_torord_ionicurrents
