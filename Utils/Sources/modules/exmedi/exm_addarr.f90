!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_addarr.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine computes some additional arrays
!> @details This routine computes some additional arrays
!> @} 
!-----------------------------------------------------------------------
subroutine exm_addarr()
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_addarr
  ! NAME 
  !    exm_addarr
  ! DESCRIPTION
  !   subroutine that computes the conductivity tensor for 2D and 3D 
  ! USES
  !    exm_...
  ! USED BY
  !    exm_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_exmedi
  implicit none
  integer(ip) :: istat,npara,ipoin,kpoin,ielem,iipar,iauxi,imate,pelty,pnode,inode
  integer(ip) :: idime,icomo,igrou,ngrou,ivalu,icmod,nrema,irema,nrefi,irefi,igaus
  real(rp)    :: xauxi,xrefi(12),xbalo(3,3),xnorm(3),xlong,xtra1,xtra2,elcod(ndime,mnode),elfib(ndime,mnode)
  real(rp)    :: cartd(ndime,mnode)
  real(rp)    :: dvolu,detjm
  real(rp)    :: xjaci(ndime,ndime),xjacm(ndime,ndime),xshap


  if( INOTMASTER ) then
     !
     ! Fibers
     !
     if( modfi_exm < 0 ) then
        fiber_exm => xfiel(-modfi_exm) % a(:,:,1)
     else if ( modfi_exm > 0 ) then
        gfibe_exm= 0.0_rp
        gfibe_exm(modfi_exm) = 1.0_rp
     end if
     !
     ! Cell types
     !
     if( kfl_heter_exm == 1_ip ) then
        celty_exm => xfiel(-modce_exm) % a(:,:,1)
     end if
     if( kfl_atbhe_exm == 1_ip ) then
        atbhe_exm => xfiel(-modab_exm) % a(:,:,1)
     end if

     if( kfl_comdi_exm == 2_ip) then

        nrema = ndime
        nrefi = npoin

        cedif_exm= 0.0_rp

        elements: do ielem = 1,nelem
           
           ! Initialize
           pelty = ltype(ielem)
           
           if( pelty > 0 ) then
              pnode = nnode(pelty)
              
              imate = 1
              if( nmate > 1 ) then
                 imate = lmate_exm(ielem)
              end if

              do inode=1,pnode
                 ipoin= lnods(inode,ielem)
                 elcod(1:ndime,inode)=coord(1:ndime,ipoin)
              end do

              if (modfi_exm < 0) then
                 
                 do inode=1,pnode
                    ipoin= lnods(inode,ielem)
                    do idime=1,ndime
                       elfib(idime,inode)=fiber_exm(idime,ipoin)
                    end do
                 end do

              else if (modfi_exm > 0) then
                 do inode=1,pnode
                    elfib(1:ndime,inode)=gfibe_exm(1:ndime)                 
                 end do
              end if

              gauss_points: do igaus=1,ngaus(pelty)

                 ! Cartesian derivatives and Jacobian
                 call elmder(pnode,ndime,elmar(pelty)%deriv(1,1,igaus),elcod,cartd,detjm,xjacm,xjaci)
                 dvolu=elmar(pelty)%weigp(igaus)*detjm

                 xrefi= 0.0_rp
                 do inode=1,pnode
                    xshap= elmar(pelty)%shape(inode,igaus)
                    do idime= 1,ndime
                       xrefi(idime)= xrefi(idime)+elfib(idime,inode)*xshap
                    end do
                 end do
             
                 xrefi(ndime+1)       = gdiff_exm(1,1,imate)  ! longitudinal intra
                 xrefi(ndime+2)       = gdiff_exm(1,2,imate)  ! cross-wise (1) intra
                 xrefi(ndime+ndime)   = gdiff_exm(1,3,imate)  ! cross-wise (2) intra
                 
                 xauxi = xrefi(1) * xrefi(1) + xrefi(2) * xrefi(2)                 
                 if (ndime==3) xauxi= xauxi+xrefi(3)*xrefi(3)
                 xauxi= sqrt(xauxi)
                 
                 xbalo= 0.0_rp
                 if (xauxi > 1.0e-8) then              
                    xbalo(1,1) = xrefi(1)/xauxi
                    xbalo(2,1) = xrefi(2)/xauxi
                    if (ndime==3) xbalo(ndime,1) = xrefi(ndime)/xauxi
                    
                    if( ndime == 2 ) then
                       
                       xbalo(1,2)=   - xbalo(2,1)
                       xbalo(2,2)=     xbalo(1,1)
                       
                       do inode=1,pnode
                          ipoin= lnods(inode,ielem)
                          xshap= elmar(pelty)%shape(inode,igaus)
                          
                          cedif_exm(1,1,ipoin) = cedif_exm(1,1,ipoin) + &
                               xshap * (xrefi(ndime+1)*xbalo(1,1)*xbalo(1,1) &
                               + xrefi(ndime+2)*xbalo(2,1)*xbalo(2,1) ) * dvolu / vmass(ipoin)
                          
                          cedif_exm(2,2,ipoin) = cedif_exm(2,2,ipoin) + &
                               xshap * (xrefi(ndime+1)*xbalo(1,2)*xbalo(1,2) &
                               + xrefi(ndime+2)*xbalo(2,2)*xbalo(2,2) ) * dvolu / vmass(ipoin)
                          
                          cedif_exm(1,2,ipoin) = cedif_exm(1,2,ipoin) + &
                               xshap * (xrefi(ndime+1)*xbalo(1,2)*xbalo(1,1) &
                               + xrefi(ndime+2)*xbalo(2,1)*xbalo(2,2) ) * dvolu / vmass(ipoin)
                          
                          
                          cedif_exm(2,1,ipoin) = cedif_exm(1,2,ipoin)
                          
                       end do
                       
                    else if (ndime==3) then
                       
                       do inode= 1,pnode
                          
                          ipoin= lnods(inode,ielem)
                          xshap= elmar(pelty)%shape(inode,igaus)

                          xnorm(1)=0.0_rp
                          xnorm(2)=1.0_rp
                          xnorm(3)=0.0_rp
                          
                          xbalo(1,2) = xnorm(2)*xbalo(3,1)-xnorm(3)*xbalo(2,1)
                          xbalo(2,2) = xnorm(3)*xbalo(1,1)-xnorm(1)*xbalo(3,1)
                          xbalo(3,2) = xnorm(1)*xbalo(2,1)-xnorm(2)*xbalo(1,1)
                          
                          xauxi=xbalo(1,2)*xbalo(1,2)+xbalo(2,2)*xbalo(2,2)+xbalo(3,2)*xbalo(3,2)
                          xauxi=sqrt(xauxi)
                          
                          if (xauxi < 0.000001_rp) then
                             xnorm(1) = 1.0_rp
                             xnorm(2) = 0.0_rp
                             xnorm(3) = 0.0_rp
                             
                             xbalo(1,2) = xnorm(2)*xbalo(3,1)-xnorm(3)*xbalo(2,1)
                             xbalo(2,2) = xnorm(3)*xbalo(1,1)-xnorm(1)*xbalo(3,1)
                             xbalo(3,2) = xnorm(1)*xbalo(2,1)-xnorm(2)*xbalo(1,1)
                             
                             xauxi = xbalo(1,2)*xbalo(1,2)+xbalo(2,2)*xbalo(2,2)+xbalo(3,2)*xbalo(3,2)
                             xauxi = sqrt(xauxi)                                  
                          end if
                          
                          xbalo(1,2) = xbalo(1,2)/xauxi
                          xbalo(2,2) = xbalo(2,2)/xauxi
                          xbalo(3,2) = xbalo(3,2)/xauxi
                          
                          ! vector axial
                          
                          xbalo(1,3) = xbalo(2,1)*xbalo(3,2)-xbalo(3,1)*xbalo(2,2)
                          xbalo(2,3) = xbalo(3,1)*xbalo(1,2)-xbalo(1,1)*xbalo(3,2)
                          xbalo(3,3) = xbalo(1,1)*xbalo(2,2)-xbalo(2,1)*xbalo(1,2)
                          
                          xlong= xrefi(ndime+1)
                          xtra1= xrefi(ndime+2)
                          xtra2= xrefi(ndime+2)
                          
                          cedif_exm(1,1,ipoin) = cedif_exm(1,1,ipoin) + &
                               xshap*(xlong*xbalo(1,1)*xbalo(1,1) + &
                               xtra1*xbalo(1,2)*xbalo(1,2) + &
                               xtra2*xbalo(1,3)*xbalo(1,3)) * dvolu / vmass(ipoin)                
                          cedif_exm(2,2,ipoin) = cedif_exm(2,2,ipoin) +&
                               xshap*(xlong*xbalo(2,1)*xbalo(2,1) + &
                               xtra1*xbalo(2,2)*xbalo(2,2) + &
                               xtra2*xbalo(2,3)*xbalo(2,3)) * dvolu / vmass(ipoin) 
                          cedif_exm(3,3,ipoin) = cedif_exm(3,3,ipoin) + &
                               xshap*(xlong*xbalo(3,1)*xbalo(3,1) + &
                               xtra1*xbalo(3,2)*xbalo(3,2) + &
                               xtra2*xbalo(3,3)*xbalo(3,3)) * dvolu / vmass(ipoin)              
                          cedif_exm(1,2,ipoin) = cedif_exm(1,2,ipoin) + &
                               xshap*(xlong*xbalo(1,1)*xbalo(2,1) + &
                               xtra1*xbalo(1,2)*xbalo(2,2) + &
                               xtra2*xbalo(1,3)*xbalo(2,3)) * dvolu / vmass(ipoin)              
                          cedif_exm(1,3,ipoin) = cedif_exm(1,3,ipoin) + &
                               xshap*(xlong*xbalo(1,1)*xbalo(3,1) + &
                               xtra1*xbalo(1,2)*xbalo(3,2) + &
                               xtra2*xbalo(1,3)*xbalo(3,3)) * dvolu / vmass(ipoin)              
                          cedif_exm(2,3,ipoin) = cedif_exm(2,3,ipoin) + &
                               xshap*(xlong*xbalo(2,1)*xbalo(3,1) + &
                               xtra1*xbalo(2,2)*xbalo(3,2) + &
                               xtra2*xbalo(2,3)*xbalo(3,3)) * dvolu / vmass(ipoin)              

                          cedif_exm(3,2,ipoin) = cedif_exm(2,3,ipoin)
                          cedif_exm(3,1,ipoin) = cedif_exm(1,3,ipoin)
                          cedif_exm(2,1,ipoin) = cedif_exm(1,2,ipoin)
                          
                          
                       end do
                       
                    end if
                    
                 end if
              end do gauss_points
              
           end if
        end do elements
     end if

     call rhsmod(ndime*ndime,cedif_exm)



  end if
  
  

end subroutine exm_addarr
