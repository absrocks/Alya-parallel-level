subroutine chm_post_projection(ivari,ielem,pnode,pgaus,gpsha,gp_Chi,gpvol) 
   !------------------------------------------------------------------------
   ! NAME 
   !    chm_post_projection
   ! DESCRIPTION
   !    Projection from Gauss points to nodes
   ! USES
   ! USED BY
   !    chm_post_scalar_dissipation_rate
   !***
   !------------------------------------------------------------------------
   use def_kintyp, only      : ip,rp
   use def_domain, only      : lnods
   use def_chemic, only      : xYr_chm, &
                               xZr_chm, &
                               xYs_chm, &
                               xZs_chm

   implicit none

   integer(ip), intent(in) :: ivari,ielem,pnode,pgaus
   real(rp),    intent(in) :: gpsha(pnode,pgaus)
   real(rp),    intent(in) :: gpvol(pgaus)
   real(rp),    intent(in) :: gp_Chi(pgaus)

   integer(ip) :: ipoin,igaus,inode

   select case (ivari)

   case (46_ip) 

      do inode = 1,pnode
         ipoin = lnods(inode,ielem)
         do igaus = 1,pgaus
            xYr_chm(ipoin) = xYr_chm(ipoin) + gpsha(inode,igaus) * gp_Chi(igaus) * gpvol(igaus)            
         end do
      end do

   case (47_ip) 

      do inode = 1,pnode
         ipoin = lnods(inode,ielem)
         do igaus = 1,pgaus
            xZr_chm(ipoin) = xZr_chm(ipoin) + gpsha(inode,igaus) * gp_Chi(igaus) * gpvol(igaus)            
         end do
      end do

   case (48_ip) 

      do inode = 1,pnode
         ipoin = lnods(inode,ielem)
         do igaus = 1,pgaus
            xYs_chm(ipoin) = xYs_chm(ipoin) + gpsha(inode,igaus) * gp_Chi(igaus) * gpvol(igaus)            
         end do
      end do

   case (49_ip) 

      do inode = 1,pnode
         ipoin = lnods(inode,ielem)
         do igaus = 1,pgaus
            xZs_chm(ipoin) = xZs_chm(ipoin) + gpsha(inode,igaus) * gp_Chi(igaus) * gpvol(igaus)            
         end do
      end do

   end select

end subroutine chm_post_projection

