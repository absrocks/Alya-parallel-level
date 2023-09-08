subroutine ada_interp(wopad)
!-----------------------------------------------------------------------
!****f* adapti/ada_interp
! NAME 
!    ada_interp
! DESCRIPTION
!    This routine interpolates the values from the coarse to the 
!    refined mesh. For nested meshes, it does it linearly.
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_inpout
  use      mod_postpr
  use      def_adapti
  use      mod_iofile
  implicit none
  integer(ip)       :: &
       lunau,npoau,nelau,ielem,iface,ipoin,jpoin,inpoi,&
       idime,pelty,pnode,inose,inode
  real(rp)          :: rlifa,rlice
  character(5)      :: wopad(2)
  real(rp), pointer :: rpove(:,:),rposc(:)

  rlifa= 1.0_rp / 2.0_rp
  rlice= 1.0_rp / 3.0_rp
  if (ndime == 3) then
     rlifa= 1.0_rp / 3.0_rp
     rlice= 1.0_rp / 4.0_rp
  end if

  if (wopad(1)=='PRESS') then 
     wopad(2)='SCALA' 
     rposc    =>  press(:,1)
  else if (wopad(1)=='DENSI') then 
     wopad(2)='SCALA'      
     rposc    =>  densi(:,1)
  else if (wopad(1)=='TEMPE') then 
     wopad(2)='SCALA'      
     rposc    =>  tempe(:,1)
  else if (wopad(1)=='ENERG') then 
     wopad(2)='SCALA'      
     rposc    =>  energ(:,1)
  else if (wopad(1)=='VISCO') then 
     wopad(2)='SCALA'      
     rposc    =>  visco(:,1)
  else if (wopad(1)=='VELOC') then 
     wopad(2)='VECTO'      
     rpove    =>  veloc(:,:,1)
  else if (wopad(1)=='UMOME') then 
     wopad(2)='VECTO'      
     rpove    =>  umome(:,:,1)
  end if

  if (wopad(2)=='SCALA') then 
     do ipoin=1,npoin
        scnew_ada(ipoin)  =rposc(ipoin)
        venew_ada(1,ipoin)=1.0_rp
     end do
     do ipoin=npoin+1,npnew_ada
        scnew_ada(ipoin)  =0.0_rp
        venew_ada(1,ipoin)=0.0_rp
     end do
     
     do ielem=1,nelem
        pelty=ltype(ielem)
        pnode=nnode(pelty)           
        if (lnseg_ada(ielem) == 1) then
           do iface=  1,pnode
              if (lnoed_ada(iface,1,ielem) > 0 ) then
                 ipoin= lnoed_ada(iface,mnode+2,ielem)
                 do inose=1,ndime
                    inpoi= lnoed_ada(iface,inose,ielem) 
                    scnew_ada(ipoin)  = scnew_ada(ipoin) + rposc(inpoi)                       
                    venew_ada(1,ipoin)= venew_ada(1,ipoin)+1.0_rp
                 end do
              end if
           end do
        else if (lnseg_ada(ielem) < 0) then
           ipoin= - lnseg_ada(ielem)
           do inode= 1,pnode
              inpoi= lnods(inode,ielem)
              scnew_ada(ipoin)= scnew_ada(ipoin) + rposc(inpoi)
              venew_ada(1,ipoin)= venew_ada(1,ipoin)+1.0_rp
           end do
        end if
     end do

     do ipoin=1,npnew_ada
        scnew_ada(ipoin)  = scnew_ada(ipoin) /  venew_ada(1,ipoin)
     end do

     
  else if (wopad(2)=='VECTO') then 
     do ipoin=1,npoin
        scnew_ada(ipoin)  =1.0_rp
        do idime=1,ndime
           venew_ada(idime,ipoin)=rpove(idime,ipoin)
        end do
     end do
     do ipoin=npoin+1,npnew_ada
        scnew_ada(ipoin)  =0.0_rp
        do idime=1,ndime
           venew_ada(idime,ipoin)=0.0_rp
        end do
     end do

     do ielem=1,nelem
        pelty=ltype(ielem)
        pnode=nnode(pelty)           
        if (lnseg_ada(ielem) == 1) then
           do iface=  1,pnode
              if (lnoed_ada(iface,1,ielem) > 0 ) then
                 ipoin= lnoed_ada(iface,mnode+2,ielem)
                 do inose=1,ndime
                    inpoi= lnoed_ada(iface,inose,ielem) 
                    scnew_ada(ipoin)  =scnew_ada(ipoin)  + 1.0_rp
                    do idime=1,ndime
                       venew_ada(idime,ipoin)= venew_ada(idime,ipoin) + rpove(idime,inpoi)
                    end do
                 end do
              end if
           end do
        else if (lnseg_ada(ielem) < 0) then
           ipoin= - lnseg_ada(ielem)
           do inode= 1,pnode
              do idime=1,ndime
                 inpoi= lnods(inode,ielem)
                 venew_ada(idime,ipoin)= venew_ada(idime,ipoin) + rpove(idime,inpoi)
              end do
              scnew_ada(ipoin)  =scnew_ada(ipoin)  + 1.0_rp
           end do
        end if
     end do
     
     
     do ipoin=1,npnew_ada
        do idime=1,ndime
           venew_ada(idime,ipoin)= venew_ada(idime,ipoin) / scnew_ada(ipoin)
        end do        
     end do
     
  end if
  
end subroutine ada_interp
