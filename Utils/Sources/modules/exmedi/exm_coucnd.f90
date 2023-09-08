!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_comcnd.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute the diffusivity
!> @details Compute the diffusivity
!> @} 
!-----------------------------------------------------------------------
subroutine exm_coucnd(imate,difin,nfibe,noion)
  use      def_master
  use      def_domain
  use      def_exmedi
  implicit none

  integer(ip) :: idime,imate,noion
  real(rp)    :: xauxi,xrefi(12),xbalo(3,3),xnorm(3),difin(ndime,ndime),nfibe(ndime),&
       xlong,xtra1,xtra2


  noion= 0
  if (kfl_cellmod(imate) == 0) then
     noion= 1
     difin(1,1)= 0.01*gdiff_exm(1,1,1)
     difin(2,2)= 0.01*gdiff_exm(1,1,1)
     difin(ndime,ndime)= 0.01*gdiff_exm(1,1,1)     
     return
  end if
  
  xrefi= 0.0_rp
  do idime= 1,ndime
     xrefi(idime)= nfibe(idime)
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

        difin(1,1) = xrefi(ndime+1)*xbalo(1,1)*xbalo(1,1) &
             + xrefi(ndime+2)*xbalo(2,1)*xbalo(2,1)

        difin(2,2) = xrefi(ndime+1)*xbalo(1,2)*xbalo(1,2) &
             + xrefi(ndime+2)*xbalo(2,2)*xbalo(2,2)

        difin(1,2) = xrefi(ndime+1)*xbalo(1,2)*xbalo(1,1) &
             + xrefi(ndime+2)*xbalo(2,1)*xbalo(2,2)

        difin(2,1) = difin(1,2)
                          
     else if (ndime == 3) then

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
        
        difin(1,1) = xlong*xbalo(1,1)*xbalo(1,1) &
             + xtra1*xbalo(1,2)*xbalo(1,2) &
             + xtra2*xbalo(1,3)*xbalo(1,3)
        difin(2,2) = xlong*xbalo(2,1)*xbalo(2,1) &
             + xtra1*xbalo(2,2)*xbalo(2,2) &
             + xtra2*xbalo(2,3)*xbalo(2,3)
        difin(3,3) = xlong*xbalo(3,1)*xbalo(3,1) &
             + xtra1*xbalo(3,2)*xbalo(3,2) &
             + xtra2*xbalo(3,3)*xbalo(3,3)
        difin(1,2) = xlong*xbalo(1,1)*xbalo(2,1) &
             + xtra1*xbalo(1,2)*xbalo(2,2) &
             + xtra2*xbalo(1,3)*xbalo(2,3)
        difin(1,3) = xlong*xbalo(1,1)*xbalo(3,1) &
             + xtra1*xbalo(1,2)*xbalo(3,2) &
             + xtra2*xbalo(1,3)*xbalo(3,3)
        difin(2,3) = xlong*xbalo(2,1)*xbalo(3,1) &
             + xtra1*xbalo(2,2)*xbalo(3,2) &
             + xtra2*xbalo(2,3)*xbalo(3,3)
        difin(3,2) = difin(2,3)
        difin(3,1) = difin(1,3)
        difin(2,1) = difin(1,2)
        
     end if
     

  end if





end subroutine exm_coucnd
