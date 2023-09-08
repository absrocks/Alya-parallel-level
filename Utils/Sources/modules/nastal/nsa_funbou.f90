subroutine nsa_funbou(itask,ifixi,xinpu,xretu)
  !------------------------------------------------------------------------
  !
  ! This subroutine computes transient boundary conditions for boundary
  ! elements coming from transient data files
  !
  ! counter ifixi can be a boundary element (iboun, itask=1) or a node (ipoin, itask=2)
  !
  !
  !-------------------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_nastal
  implicit none
  
  integer(ip) :: itask,idata,ifunc,ifixi,idofn
  real(rp)    :: t1,p1,t2,p2, b, m, cdeti, venew
  real(rp), INTENT(OUT) ::xretu(ndofn_nsa)
  real(rp), INTENT(IN)  ::xinpu(ndofn_nsa)
  real(rp), external    :: funcre



  if (itask == 1_ip) then
     !
     ! on boundary elements
     !
 
  else   if (itask == 2_ip) then
     !
     ! on boundary nodes
     !
     ifunc = kfl_funno_nsa(ifixi)

     if(kfl_funty_nsa(ifunc,1)/=3_ip) then
       !
       ! periodic time function
       !
       venew = funcre(fubcs_nsa(ifunc,1:8),8_ip,kfl_funty_nsa(ifunc,1),cutim)

       do idofn = 1,ndofn_nsa
         if (kfl_funty_nsa(ifunc,idofn+2) == 1) then
           xretu(idofn) = xinpu(idofn) * venew
         else
           xretu(idofn) = xinpu(idofn)
         end if
       end do


     else
       !
       ! discrete time function
       !
       do idata= 1,mtloa_nsa(ifunc)-1
  
          if (cutim >= tload_nsa(ifunc)%a(ndofn_nsa+1,idata)) then
             if (cutim < tload_nsa(ifunc)%a(ndofn_nsa+1,idata+1)) then
  
                t1= tload_nsa(ifunc)%a(ndofn_nsa+1,idata)
                t2= tload_nsa(ifunc)%a(ndofn_nsa+1,idata+1)
                do idofn=1,ndofn_nsa
                   if (kfl_fixno_nsa(idofn,ifixi) > 0) then
                      p1= tload_nsa(ifunc)%a(idofn,idata  )
                      p2= tload_nsa(ifunc)%a(idofn,idata+1)
                      m=(p2-p1)/(t2-t1)
                      b=p1-m*t1
                      xretu(idofn)=m*cutim+b                              ! return value
                   end if
                end do
                exit
             end if
          end if
       end do     

     end if

   else   if (itask == 3_ip) then
     !
     ! on every node (body forces)
     !

  end if


end subroutine nsa_funbou
