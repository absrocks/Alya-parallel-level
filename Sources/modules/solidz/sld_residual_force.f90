!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_residual_force.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute the residual force over a group of points
!> @details Compute the residual force over a group of points
!> @} 
!-----------------------------------------------------------------------
subroutine sld_residual_force
  use      def_parame
  use      def_master
  use      def_domain
  use      def_solidz
  use      mod_communications, only : PAR_SUM, PAR_MAX, PAR_MIN

  implicit none

  integer(ip)  :: ipoin, idime, idofn, inode
  real(rp)     :: dum1min, dum1max, dum2min, dum2max, dum3min, dum3max

  if (kfl_rsfor_sld == 0) return

  resdispl_sld = 0.0_rp
  resforce_sld = 0.0_rp
 
  !
  ! By position greater than a fixed value (one direction)
  !
  if (kfl_rsfor_sld == 1_ip) then
     do ipoin = 1,npoin
        if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
           idofn = (ipoin-1) * ndime
           if (coord(idilimi_sld(1),ipoin) > reslimi_sld(1)) then  
              do idime =1,ndime
                 if(kfl_fixno_sld(idime,ipoin) > 0) then
                    resforce_sld(idime,1) = resforce_sld(idime,1) + frxid_sld(idofn+idime) 
                    resdispl_sld(idime,1) = displ(idime,ipoin,1)
                 end if
              end do
           end if
        end if
     end do

  !
  ! By position greater than a fixed value (three directions)
  !
  else if (kfl_rsfor_sld == 2_ip) then 

     do ipoin=1,npoin
        if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
           idofn = (ipoin-1) * ndime
           if ((coord(1,ipoin) > reslimi_sld(1)) .and. (coord(2,ipoin) > reslimi_sld(2)) .and.  (coord(3,ipoin) > reslimi_sld(3))) then
              do idime=1,ndime
                 if(kfl_fixno_sld(idime,ipoin) > 0) then
                    resforce_sld(idime,1) = resforce_sld(idime,1) + frxid_sld(idofn+idime)
                    resdispl_sld(idime,1) = displ(idime,ipoin,1)
                 end if
              end do
           end if
        end if
     end do

  !
  ! By node set
  !
  else if (kfl_rsfor_sld == 3_ip) then

     do ipoin=1,size(resnode_sld,1)
        if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
           inode = resnode_sld(ipoin)
           idofn = (inode-1) * ndime
           do idime=1,ndime
              if(kfl_fixno_sld(idime,ipoin) > 0) then
                 resforce_sld(idime,1) = resforce_sld(idime,1) + frxid_sld(idofn+idime)
                 resdispl_sld(idime,1) = displ(idime,inode,1)
              end if
           end do
        end if
     end do

  end if
  
  !
  ! Parall exchange reaction force 
  !
  call PAR_SUM(resforce_sld(1,1),'IN MY CODE')
  call PAR_SUM(resforce_sld(2,1),'IN MY CODE')
  call PAR_SUM(resforce_sld(3,1),'IN MY CODE')

  !
  ! Parall exchange displacement
  !
  dum1max = resdispl_sld(1,1)
  dum1min = resdispl_sld(1,1)
  call PAR_MIN(dum1min, 'IN MY CODE')
  call PAR_MAX(dum1max, 'IN MY CODE')
  if( dum1min == 0.0_rp ) resdispl_sld(1,1) = dum1max
  if( dum1max == 0.0_rp ) resdispl_sld(1,1) = dum1min

  dum2max = resdispl_sld(2,1)
  dum2min = resdispl_sld(2,1)
  call PAR_MIN(dum2min, 'IN MY CODE')
  call PAR_MAX(dum2max, 'IN MY CODE')
  if( dum2min == 0.0_rp ) resdispl_sld(2,1) = dum2max
  if( dum2max == 0.0_rp ) resdispl_sld(2,1) = dum2min

  dum3max = resdispl_sld(3,1)
  dum3min = resdispl_sld(3,1)
  call PAR_MIN(dum3min, 'IN MY CODE')
  call PAR_MAX(dum3max, 'IN MY CODE')
  if( dum3min == 0.0_rp ) resdispl_sld(3,1) = dum3max 
  if( dum3max == 0.0_rp ) resdispl_sld(3,1) = dum3min     


end subroutine sld_residual_force
