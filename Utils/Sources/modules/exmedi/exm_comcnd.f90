!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_comcnd.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute the diffusivity tensor
!> @details Compute the diffusivity tensor
!> @} 
!-----------------------------------------------------------------------
subroutine exm_comcnd(ielem,imate,difin,noion)
  use      def_master
  use      def_domain
  use      def_exmedi
  implicit none

  real(rp)    :: &
       difin(ndime,ndime,mgaus),xshap,xceco(mgaus),&
       auxic(3,3),difin_norm, aux
  integer(ip) :: ielem,jelem,idime,jdime,kdime,pdime,inode,igaus,ipoin,icond,&
       pelty,pnode,pgaus,imate, itest,noion

  itest= 0
  if (kfl_comdi_exm == 0 .or. ielem < 0) itest = 1

  if (itest == 1) then
     
     difin= 0.0_rp
     difin(1,1,1)= gdiff_exm(1,1,imate)
     difin(2,2,1)= gdiff_exm(1,2,imate)
     difin(ndime,ndime,1)= gdiff_exm(1,2,imate)
     

  else 

     ! with fibers

     jelem= ielem
     !!     if (ielem < 0) jelem= -ielem

     pelty = ltype(jelem)
     pnode = nnode(pelty)
     pgaus = ngaus(pelty)

     !!     if (ielem < 0) pgaus=1

     ! Gather operations     
     difin_norm=0.0_rp
     do igaus=1,pgaus
        difin(1,1,igaus) = 0.0_rp
        difin(2,2,igaus) = 0.0_rp
        difin(1,2,igaus) = 0.0_rp
        difin(2,1,igaus) = 0.0_rp
        do inode = 1,pnode
           ipoin = lnods(inode,jelem)
           xshap = elmar(pelty)%shape(inode,igaus)
           !!           if (ielem < 0) xshap= 1.0_rp / real(pnode)
           difin(1,1,igaus) = difin(1,1,igaus) + xshap * cedif_exm(1,1,ipoin)
           difin(2,1,igaus) = difin(2,1,igaus) + xshap * cedif_exm(2,1,ipoin)
           difin(1,2,igaus) = difin(1,2,igaus) + xshap * cedif_exm(1,2,ipoin)
           difin(2,2,igaus) = difin(2,2,igaus) + xshap * cedif_exm(2,2,ipoin)              
        end do

        difin_norm= difin_norm &
             + difin(1,1,igaus) * difin(1,1,igaus) &
             + difin(2,2,igaus) * difin(2,2,igaus) &
             + difin(1,2,igaus) * difin(1,2,igaus) &
             + difin(2,1,igaus) * difin(2,1,igaus)
        
     end do

     if (ndime == 3) then
        difin_norm = 0.0_rp
        do igaus=1,pgaus
           difin(3,1,igaus) = 0.0_rp
           difin(3,2,igaus) = 0.0_rp
           difin(1,3,igaus) = 0.0_rp
           difin(2,3,igaus) = 0.0_rp
           difin(3,3,igaus) = 0.0_rp
           do inode = 1,pnode
              ipoin = lnods(inode,jelem)
              xshap = elmar(pelty)%shape(inode,igaus)
              !!              if (ielem < 0) xshap= 1.0_rp / real(pnode)
              difin(3,1,igaus) = difin(3,1,igaus) + xshap * cedif_exm(3,1,ipoin)
              difin(3,2,igaus) = difin(3,2,igaus) + xshap * cedif_exm(3,2,ipoin)              
              difin(1,3,igaus) = difin(1,3,igaus) + xshap * cedif_exm(1,3,ipoin)
              difin(2,3,igaus) = difin(2,3,igaus) + xshap * cedif_exm(2,3,ipoin)
              difin(3,3,igaus) = difin(3,3,igaus) + xshap * cedif_exm(3,3,ipoin)              
           end do

           difin_norm= difin_norm &
                + difin(3,1,igaus) * difin(3,1,igaus) &
                + difin(3,2,igaus) * difin(3,2,igaus) &
                + difin(1,3,igaus) * difin(1,3,igaus) &
                + difin(2,3,igaus) * difin(2,3,igaus) &
                + difin(3,3,igaus) * difin(3,3,igaus) 

        end do
        
     end if

     noion= 0
     if (difin_norm < 1.0e-10) then
        !noion= 1
        aux= (gdiff_exm(1,1,imate)+gdiff_exm(1,2,imate))/2.0_rp
        difin(1,1,imate)= aux
        difin(2,2,imate)= aux
        difin(ndime,ndime,imate)= aux
     end if

     if (kfl_cemod_exm == 2) then
        call runend('EXM_COMCND: NO BIDOMAIN MODEL WITH FIBERS PROGRAMMED')
     end if

  end if


end subroutine exm_comcnd
