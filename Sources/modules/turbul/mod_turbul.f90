module mod_turbul

   use def_kintyp, only: ip, rp, lg
   use def_master
   use def_turbul
   implicit none
   private

   public :: turbul_main

contains

   subroutine turbul_main(order)
      implicit none
      integer(ip), intent(in) :: order
      call Turbul(order)
   end subroutine turbul_main

end module mod_turbul
