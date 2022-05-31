!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_coupli.f90
!> @author  Simone Marras
!> @date    16/11/1966
!> @brief   Read initial condition from a restart file
!> @details Read initial condition from a restart file
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_rearst(outvar)
  use def_master
  use def_domain
  use def_nastal
  use def_parame

  implicit none
  integer(ip) :: outvar  !< Which is the variable

  integer(ip) :: ipoin,kpoin
  real(rp)    :: dummyr

  integer(ip) :: icase_npoin
  character   :: restart_file_name*24
  

  !Re-initiaize bvess_nsa to zero:
  bvess_nsa = 0.0_rp

  !
  ! Open and read file
  !
  open(1,file='restart-file.txt')
  !
  ! The case200 file is stored in npoin rows (obtained from the matlab routines):
  ! 
  if(outvar == 1) then
     !
     ! bvess_nsa(ndofn_nsa-1,ipoin,1) <-- DENSI
     !
     do ipoin = 1,npoin
        !
        !    press     densi      theta    u-velo w-velo 
        !
        read(1,*) &
             dummyr, & !THIS IS DUMMY BECAUSE bvess stores densi when outvar == 1
             bvess_nsa(ndofn_nsa-1,ipoin,1), & !DDHYD (rho')
             bvess_nsa(ndofn_nsa,  ipoin,1), & !THBAS (theta')
             bvess_nsa(1        ,  ipoin,1), & !XVELO
             bvess_nsa(ndime    ,  ipoin,1)    !YVELO
        
     end do

  else if(outvar == 2) then
     !
     ! bvess_nsa(ndofn_nsa-1,ipoin,1) <-- PRESS
     !
     do ipoin = 1,npoin
        !
        !    press     densi      theta    u-velo w-velo 
        !
        read(1,*) &
             bvess_nsa(ndofn_nsa-1,ipoin,1), & !DPHYD (p')
             dummyr, & !THIS IS DUMMY BECAUSE bvess stores pressure when outvar == 2
             bvess_nsa(ndofn_nsa,  ipoin,1), & !THBAS (theta')
             bvess_nsa(1        ,  ipoin,1), & !XVELO
             bvess_nsa(ndime    ,  ipoin,1)    !YVELO
        
     end do
  end if
  close(1)


  !
  ! The file contains the perturbation variables:
  ! add the background hydrsotatic state to store the total
  ! variables into bvess_nsa:
  !
  if(outvar == 1) then
     do ipoin = 1,npoin
        !Densi
        bvess_nsa(ndofn_nsa-1,ipoin,1) = bvess_nsa(ndofn_nsa-1,ipoin,1) + rekee_nsa(ndime+1,ipoin)
        !Tempe
        bvess_nsa(ndofn_nsa,ipoin,1)   = bvess_nsa(ndofn_nsa,ipoin,1)   + rekee_nsa(ndime+2,ipoin)
     end do

  else if(outvar == 2) then
     do ipoin = 1,npoin
        !Press
        bvess_nsa(ndofn_nsa-1,ipoin,1) = bvess_nsa(ndofn_nsa-1,ipoin,1) + rekee_nsa(ndime+3,ipoin)
        !Tempe
        bvess_nsa(ndofn_nsa,ipoin,1)   = bvess_nsa(ndofn_nsa,ipoin,1)   + rekee_nsa(ndime+2,ipoin)
     end do
  end if

end subroutine nsa_rearst
