!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_massma.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Obtains lumpped mass matrix from vmass and the nodal values of porosity and Bw - also sends Mass_mat * Sw_n to RHS
!> @details Obtains lumpped mass matrix from vmass and the nodal values of porosity and Bw - also sends Mass_mat * Sw_n to RHS
!!          Similar to nsi_assdia
!> @} 
!------------------------------------------------------------------------
subroutine por_massma()
  !
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  npoin
  use def_porous, only       :  nodpo_por,prref_por,comwa_por,comro_por,bwref_por,dtinv_por,denwa_por
  use def_master, only       :  pmatr,press,INOTMASTER

  implicit none
  real(rp)    :: deltp,xx_w,xx_r,Bw,poro
  integer(ip)                :: ipoin

  if( INOTMASTER ) then

     do ipoin = 1,npoin
        !
        ! pmatr
        !
        deltp = press(ipoin,1) - prref_por
        xx_w = comwa_por * deltp
        xx_r = comro_por * deltp
        Bw = bwref_por / (1.0_rp + xx_w + (0.5_rp*xx_w*xx_w) )
        poro = nodpo_por(ipoin) * (1.0_rp + xx_r + (0.5_rp*xx_r*xx_r) )
        pmatr(ipoin) =  Bw / (denwa_por * poro * dtinv_por)
     end do

  end if

end subroutine por_massma
