subroutine nsa_aebody(itask,iboun,ielem,pnodb,ifibo)
!-----------------------------------------------------------------------
!****f* nastal/nsa_aebody
! NAME 
!    nsa_elcons
! DESCRIPTION
!    This subroutine compute aerodynamic body coefficients    
! USES
!    nsa_...
! USED BY
!    nsa_elcons
!    nsa_elfmom
!    nsa_elchea
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_nastal
  implicit none
  integer(ip)    :: &
       idime,itask,inode,iboun,pnodb,ielem,ipoin,ifibo,ibobo,inodb,kfoun
  real(rp)       :: &
       xcoef,tcoef(3,20)
  
  if (kfl_bodyf_nsa == 1) then
     
     call nsa_chkpar                  ! determine parallelization status 

     if (itask == 1) then
        
        ifibo=0
        ibobo= 0
        inodb= 1
        do while (ibobo == 0)
           kfoun= 0
           inode = lboel(inodb,iboun)
           ipoin = lnods(inode,ielem)                 
           if (coord(1,ipoin) > bodyr_nsa(1,1,1)) then
              if (coord(1,ipoin) < bodyr_nsa(1,2,1)) then
                 kfoun= kfoun+1
              end if
           end if
           if (kfoun == 1) then
              if (coord(2,ipoin) > bodyr_nsa(1,1,2)) then
                 if (coord(2,ipoin) < bodyr_nsa(1,2,2)) then
                    kfoun= kfoun+1
                 end if
              end if
           end if
           if ((kfoun == 2) .and. (ndime==3)) then
              if (coord(3,ipoin) > bodyr_nsa(1,1,3)) then
                 if (coord(3,ipoin) < bodyr_nsa(1,2,3)) then
                    kfoun= kfoun+1
                 end if
              end if
           end if
           if (kfoun == ndime) then
              ifibo= ifibo+1
              ibobo= 1
             end if
           inodb= inodb + 1
           if (inodb > pnodb) ibobo=1
        end do        

     else if (itask == 2) then
        
        sufro_nsa = -sufro_nsa     ! needed, because of the way we compute the frontal areas
        if (ndime==2) then
           ! in 2D problems, the wing area is the airfoil cord, which if X is the streamwise
           ! direction, then the cord is the frontal area in the Y direction
           surto_nsa= sufro_nsa(2) 
        end if

        ccoef_nsa(1)      = surto_nsa                                                     ! total area
        do idime=1,ndime
           ccoef_nsa(idime+1)      = fbody_nsa(1,idime)                                   ! press force
           ccoef_nsa(idime+4)      = fbody_nsa(2,idime)                                   ! visco force
           ccoef_nsa(idime+7)      = sufro_nsa(idime)                                     ! frontal area
        end do
        
     else if (itask == 3) then

        call nsa_parall(6)

        tcoef= 0.0_rp

        if (imaster .or. iloner) then
           xcoef=  1.0_rp / densi_nsa / speed_nsa / speed_nsa
           do idime=1,ndime
              tcoef(idime,1) = 2_rp * ccoef_nsa(idime+1) * xcoef / ccoef_nsa(1)       ! press cx total
              tcoef(idime,2) = 2_rp * ccoef_nsa(idime+4) * xcoef / ccoef_nsa(1)       ! visco cx total
              tcoef(idime,3) = 2_rp * ccoef_nsa(idime+1) * xcoef / ccoef_nsa(idime+7) ! press cx frontal
              tcoef(idime,4) = 2_rp * ccoef_nsa(idime+4) * xcoef / ccoef_nsa(idime+7) ! visco cx frontal
           end do

!#  --|  1. Time Step   2. Global Iteration   3. Inner Iteration   4. Current time
!#  --|  5-7.   Press CL/D (total area)  X,Y,Z
!#  --|  8-10.  Visco CL/D (total area)  X,Y,Z
!#  --|  11-13. Press CX (frontal area)  X,Y,Z
!#  --|  14-16. Visco CX (frontal area)  X,Y,Z
!#  --|  17-19. Press force  X,Y,Z
!#  --|  20-22. Visco force  X,Y,Z

           if ((ittim == 1) .and. (itcou ==1)) then
              write(lun_force_nsa,1010) '#  --|    Total area     =  ', ccoef_nsa( 1)  
              write(lun_force_nsa,1010) '#  --|    X-Frontal area =  ', ccoef_nsa( 8)
              write(lun_force_nsa,1010) '#  --|    Y-Frontal area =  ', ccoef_nsa( 9)
              write(lun_force_nsa,1010) '#  --|    Z-Frontal area =  ', ccoef_nsa(10)
           end if
           write(lun_force_nsa,1000) ittim,itcou,itinn(modul),cutim,&
                tcoef(1:3,1) , &
                tcoef(1:3,2) , &
                tcoef(1:3,3) , &
                tcoef(1:3,4) , &
                ccoef_nsa(2 ), ccoef_nsa(3 ) , ccoef_nsa(4 ), &
                ccoef_nsa(5 ), ccoef_nsa(6 ) , ccoef_nsa(7 )
           flush(lun_force_nsa)
        end if
        
     end if

  end if
  
1000 format(4x,i9,2x,i9,2x,i9,22(2x,e12.6))
1010 format(a,2x,e12.6)

end subroutine nsa_aebody
