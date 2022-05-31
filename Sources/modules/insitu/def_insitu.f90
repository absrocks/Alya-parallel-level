!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    def_insitu.f90
!> @date    27/03/2017
!> @author  Vishal Mehta
!> @brief   Def file
!> @details Def file
!> @}
!------------------------------------------------------------------------
module def_insitu
  use def_kintyp
  integer(ip) :: PAR_CURRENT_RANK_INS,PAR_CURRENT_SIZE_INS
  !! Variables for coupling
  real(rp), pointer        :: advec(:,:)
  real(rp), pointer        :: elmag(:)
  real(rp), pointer        :: elmag_buf(:)
  integer(ip),pointer      :: global_num(:)
  integer(ip),pointer      :: npoin_all(:)
  integer(ip)              :: NDglobal
  
end module def_insitu
