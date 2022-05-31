subroutine element_grid_test(ndime,pnode,pelty,elcod,coglo,coloc)

  use def_kintyp, only : ip,rp
  use mod_elmgeo, only : elmgeo_shape2
  use mod_elmgeo, only : element_type
  implicit none

  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: pelty
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: coglo(*)
  real(rp),    intent(out) :: coloc(*)
  integer(ip)              :: nn(3),ii,jj,kk
  integer(ip)              :: inode,inod1,inod2,iedge
  integer(ip)              :: iimin,jjmin
  integer(ip)              :: iimax,jjmax,iiter
  real(rp)                 :: iso_min(3)
  real(rp)                 :: iso_max(3)
  real(rp)                 :: dn(3),xx
  real(rp)                 :: s(3),x(3)
  real(rp)                 :: dista,dista_min
  real(rp)                 :: gpsha(pnode)
  real(rp)                 :: smin(3),smax(3),sdif(3)
  real(rp)                 :: smin_sav(3),smax_sav(3),sdif_sav(3)
  real(rp)                 :: xnorm,xnorm_min

  iiter     =  0
  smin      = -1.0_rp
  smax      =  1.0_rp
  sdif      =  smax-smin
  nn        =  10
  dn        =  1.0_rp / real(nn-1_ip,rp)
  xnorm_min =  huge(1.0_rp)
  dista_min =  huge(1.0_rp)

  do iedge = 1,element_type(pelty) % number_edges
     inod1 = element_type(pelty) % list_edges(1,iedge)
     inod2 = element_type(pelty) % list_edges(2,iedge)
     xnorm = dot_product(elcod(1:ndime,inod1)-elcod(1:ndime,inod2),elcod(1:ndime,inod1)-elcod(1:ndime,inod2))
     if( xnorm < xnorm_min ) xnorm_min = xnorm
  end do

  do while( sqrt(dista_min/xnorm_min) > 1.0e-3_rp )

     iiter     = iiter + 1 
     smin_sav  = smin
     smax_sav  = smax
     sdif_sav  = sdif
     dista_min = huge(1.0_rp)
     
     do ii = 1,nn(1)
        s(1) = smin(1) + sdif(1) * real(ii-1,rp) * dn(1)
        do jj = 1,nn(2)
           s(2) = smin(2) + sdif(2) * real(jj-1,rp) * dn(2)
           x = 0.0_rp
           call elmgeo_shape2(s(1),s(2),pnode,gpsha)
           do inode = 1,pnode
              x(1:2) = x(1:2) + gpsha(inode) * elcod(1:2,inode)
           end do
           dista = dot_product(x(1:2)-coglo(1:2),x(1:2)-coglo(1:2))
           if( dista <= dista_min ) then
              dista_min  = dista
              coloc(1:2) = s(1:2)
              iimin      = max(1_ip ,ii-1)
              iimax      = min(nn(1),ii+1)
              jjmin      = max(1_ip ,jj-1)
              jjmax      = min(nn(2),jj+1)
           end if
        end do
     end do
     
     smin(1)    = smin_sav(1) + sdif_sav(1) * real(iimin-1,rp) * dn(1)
     smax(1)    = smin_sav(1) + sdif_sav(1) * real(iimax-1,rp) * dn(1)
     
     smin(2)    = smin_sav(2) + sdif_sav(2) * real(jjmin-1,rp) * dn(2)
     smax(2)    = smin_sav(2) + sdif_sav(2) * real(jjmax-1,rp) * dn(2)                
     sdif       = smax-smin
  end do

  print*,'iiter=',iiter
  print*,'a=',coloc(1:2)
  
end subroutine element_grid_test
