module mod_sld_cou
  use def_kintyp,         only :  ip, rp
  use def_coupli,         only : mcoup
  use def_coupli,         only : coupling_type
  use mod_parall,         only :  PAR_MY_CODE_RANK
  use mod_communications, only :  PAR_SUM, PAR_BARRIER
  use def_master,         only :  kfl_eccty, kfl_cellmod
  use def_domain,         only: nmate
  use mod_exm_sld_eccoupling,  only :  cell_ca0_ecc
  use mod_exm_sld_eccoupling,  only : kfl_exmsld_3Dcou_ecc

  implicit none

  !private ::

  public :: mod_sld_cou_initvar
  public :: mod_sld_cou_initexchange
  public :: mod_sld_cou_physics_initialised
    


contains
    !---------------------------------------------------------------
    !>
    !> @author  Alfonso Santiago
    !> @date    24/09/2018
    !> @brief   
    !> @details
    !>
    !---------------------------------------------------------------
    subroutine mod_sld_cou_initvar()
      implicit none
      integer(ip)   :: i

      kfl_exmsld_3Dcou_ecc=.False.
      if(mcoup.ge.1_ip)then
        do i=1,mcoup
          if(coupling_type(i) % variable .eq. 'CALCI')then
             kfl_exmsld_3Dcou_ecc=.True.
          endif
         enddo
       endif
    endsubroutine mod_sld_cou_initvar
    !---------------------------------------------------------------
    !>
    !> @author  Alfonso Santiago
    !> @date    24/09/2018
    !> @brief   
    !> @details
    !>
    !---------------------------------------------------------------
    subroutine mod_sld_cou_initexchange()
      use def_master,         only :  kfl_exm_max_nmaterials
      implicit none

      integer(ip)                          :: i
      integer(ip), pointer, dimension(:)   :: pv

       allocate(pv(kfl_exm_max_nmaterials))

       ! Receive kfl_eccty
        kfl_eccty=0_ip
        pv=0_ip

       call PAR_SUM(pv, 'IN THE UNIVERSE')


       do i=1,kfl_exm_max_nmaterials
          kfl_eccty(i)=pv(i)
       enddo


        ! Send kfl_cellmod only one time
       
       do i=1,kfl_exm_max_nmaterials
          pv(i)=kfl_cellmod(i)
       enddo

       IF(PAR_MY_CODE_RANK.ne.1_ip) pv(:)=0_ip

        call PAR_SUM(pv, 'IN THE UNIVERSE')
 
        deallocate(pv)
    endsubroutine mod_sld_cou_initexchange

     !---------------------------------------------------------------
     !>
     !> @author  Alfonso Santiago
     !> @date    24/09/2018
     !> @brief
     !> @details
     !>
     !---------------------------------------------------------------
     subroutine mod_sld_cou_physics_initialised()
       use def_master,         only :  kfl_exm_max_nmaterials

       implicit none

       integer(ip)                          :: i,j
       real(rp), pointer, dimension(:,:)    :: pM

       allocate(pM(3_ip,kfl_exm_max_nmaterials))
        pM=0.0_rp


        ! Send intitial cell calcium cell_ca0_ecc

       do i=1,nmate
          do j=1_ip,3_ip
          pM(j,i)=cell_ca0_ecc(j,i)
          enddo
       enddo

       IF(PAR_MY_CODE_RANK.ne.1_ip) pM(:,:)=0.0_rp
    
        call PAR_SUM(pM, 'IN THE UNIVERSE')

       deallocate(pM)
     endsubroutine mod_sld_cou_physics_initialised

endmodule mod_sld_cou

