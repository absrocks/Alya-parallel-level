subroutine ada_elsest
!-----------------------------------------------------------------------
!****f* adapti/ada_elsest
! NAME 
!    ada_modmsh
! DESCRIPTION
!    This routine use elsest to identify the immersed object situation 
!    relative to the background mesh
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain

  use      def_adapti
  implicit none

  integer(ip)   :: lopoi(7)
  integer(ip)   :: &
       ipoii,kpoin,jpoin,iimmo,ifoun,npoie,&
       lenin,lenre,inofo,iyefo,isobp,luels,iele1,ipoto
  real(rp)      :: toaux,doaux,topct,dtpct,xshac
  real(rp)      :: loarr(50),xbbox(3,3),xnter(3),rdote,rdone
  character(80) :: trafi
  
  lopoi(1)=1
  lopoi(2)=lopoi(1)+ndime
  lopoi(3)=lopoi(2)+mnode+1
  lopoi(4)=lopoi(3)+ndime*(mnode+1)
  lopoi(5)=lopoi(4)+ndime*ndime
  lopoi(6)=lopoi(5)+ndime*ndime
  lopoi(7)=lopoi(6)+ndime
  
  lenin = 8
  lenre = 8

  rdote=9.9999d0
  title=''
  topct= 1.0             !percentage off
  dtpct= 0.5

  luels= lun_outpu_ada

  do iimmo = 1,nimmo_ada     
     npoie= npoii_ada(iimmo)
     inofo=0
     iyefo=0
     isobp=2

     do ipoii= 1, npoie
        iele1 = 0
        ifoun = 0
        loarr(1:ndime) = coori_ada(1:ndime,ipoii,iimmo)        
        toaux= topct
        doaux= dtpct

        do while (ifoun.eq.0)            
           !call Elsest(                                            &
           !     1,1,0,1,iele1,ifoun,                               &
           !     npoin,mnode,nelem,2,mnode,ndime,                   &
           !     lenin,lenre,0,luels,0,0,                           &
           !     0,npoie,title,coord,lnods,coori_ada(1,ipoii,iimmo),&
           !     loarr(lopoi(2)),loarr(lopoi(3)),loarr(lopoi(4)),   &
           !     loarr(lopoi(5)),loarr(lopoi(6)),loarr(lopoi(7)),   &
           !     xnter,toaux,xbbox)
           
           if (ifoun.eq.0) then
              toaux= toaux+doaux
              npoie=npoie+1
              !                  if (toaux.gt.3.0d0) then
              if (toaux.gt.1.5d0) then
                 !  if (toaux.gt.5.5d0) then
                 ifoun = -1
              end if
           else if (ifoun.eq.1) then
              lelim_ada(ipoii,iimmo) = iele1
              iyefo=iyefo+1
           end if
           
        end do
        
        if (ifoun.eq.-1) then !out of the background mesh
           inofo=inofo+1
           lelim_ada(ipoii,iimmo) = ifoun
        else
           call ada_scagec(ipoii,iele1,iimmo,loarr(ndime+1))           
        end if
        
        rdone=float(ipoii*100/npoie)
        
        if (rdone.gt.rdote) then
           rdote=rdote+9.99999d0
!!           write (6,*) 'Done: ',rdone,'% of the patch mesh'
        end if
        
        
     end do
     
  end do

  


end subroutine ada_elsest
