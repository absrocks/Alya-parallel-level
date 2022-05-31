module mod_exm_cou
  use def_kintyp,              only :  ip, rp
  use def_coupli,              only : mcoup
  use def_coupli,              only : coupling_type
  use mod_parall,              only :  PAR_MY_CODE_RANK
  use mod_communications,      only :  PAR_SUM, PAR_BARRIER
  use def_master,              only :  kfl_eccty, kfl_cellmod, kfl_exm_max_nmaterials
  use mod_exm_sld_eccoupling,  only : set_has_land
  use mod_exm_sld_eccoupling,  only : kfl_exmsld_3Dcou_ecc
  use mod_exm_sld_eccoupling,  only :  cell_ca0_ecc

  implicit none


contains
    !---------------------------------------------------------------
    !>
    !> @author  Alfonso Santiago
    !> @date    24/09/2018
    !> @brief   
    !> @details
    !>
    !---------------------------------------------------------------
    subroutine mod_exm_cou_initvar()
      implicit none
      integer(ip)       :: i
      kfl_exmsld_3Dcou_ecc=.False.
      if(mcoup.ge.1_ip)then
        do i=1,mcoup
          if(coupling_type(i) % variable .eq. 'CALCI')then
             kfl_exmsld_3Dcou_ecc=.True.
          endif
        enddo
      endif
    endsubroutine mod_exm_cou_initvar
    !---------------------------------------------------------------
    !>
    !> @author  Alfonso Santiago
    !> @date    24/09/2018
    !> @brief   
    !> @details
    !>
    !---------------------------------------------------------------
    subroutine mod_exm_cou_initexchange()
      use def_master,         only :  kfl_exm_max_nmaterials
      implicit none

      integer(ip)                          :: i
      integer(ip), pointer, dimension(:)   :: pv

      allocate(pv(kfl_exm_max_nmaterials))

      ! Send eccty

      do i=1,kfl_exm_max_nmaterials
         pv(i)=kfl_eccty(i)
      enddo

      IF(PAR_MY_CODE_RANK.ne.1_ip) pv=0_ip

      call PAR_SUM(pv,'IN THE UNIVERSE')


      ! Receive kfl_cellmod 
      pv=0_ip
      call PAR_SUM(pv,'IN THE UNIVERSE')

      do i=1,kfl_exm_max_nmaterials
         kfl_cellmod(i)=pv(i)
      enddo


      deallocate(pv)
    endsubroutine mod_exm_cou_initexchange

    !---------------------------------------------------------------
    !>
    !> @author  Alfonso Santiago
    !> @date    24/09/2018
    !> @brief   
    !> @details
    !>
    !---------------------------------------------------------------
    subroutine mod_exm_cou_physics_initialised()
      use def_master,         only :  kfl_exm_max_nmaterials

      implicit none

      integer(ip)                          :: i,j
      real(rp), pointer, dimension(:,:)    :: pM

      !
      ! Set if land is coupled
      !
      call set_has_land()


      allocate(pM(3_ip,kfl_exm_max_nmaterials))

      ! Receive initial cell calcium cell_ca0_ecc
      pM=0.0_rp
      call PAR_SUM(pM,'IN THE UNIVERSE')

      do i=1,kfl_exm_max_nmaterials
         do j=1,3_ip
            cell_ca0_ecc(j,i)=pM(j,i)
         enddo
      enddo

      deallocate(pM)
    endsubroutine mod_exm_cou_physics_initialised
endmodule mod_exm_cou

