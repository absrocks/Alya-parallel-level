subroutine hlm_elmopefirst()

  use def_parame
  use def_master
  use def_domain
  use def_helmoz

  implicit none

  complex(rp)         :: elmat(4*mnode,4*mnode),elrhs(4*mnode)   !Element matrix, element RHS
  complex(rp)         :: pvecpo(ndime,mnode)                     !Primary magnetic vector potential
  complex(rp)         :: pscapo(mnode)                           !Primary electric scalar potential
  integer(ip)         :: ielem,igaus,idime,inode,ii              !Indices and dimensions
  integer(ip)         :: pelty,pmate,pnode,ipoin,poin
  integer(ip)         :: pgaus,plapl,porde,ptopo
  real(rp)            :: elcod(ndime,mnode)
  real(rp)            :: gpvol(mgaus)                            !w * |J|(gaus point)
  complex(rp)         :: gprea(ncond_hlm,mgaus),dgprea(ncond_hlm,mgaus)  !Reaction terms
  complex(rp)         :: gprhs(4,mnode)                          !f (all terms)
  real(rp)            :: gpcar(ndime,mnode,mgaus)                !dNk/dxj ...
  real(rp)            :: gphes(ntens,mnode,mgaus)                !dNk/dxidx ...

  real(rp)            :: dmax,dmin,rcons

  real(rp)            :: dsigmatmp
  integer(ip)         :: counterInnerElems, isInnerElem
  real(rp)            :: centerX, centerY, centerZ
  real(rp)            :: minX, minY, minZ
  real(rp)            :: maxX, maxY, maxZ

  if(kfl_scheme_opt==4) then

     design_vars(:)=1.0E+030_rp
     design_vars_ref(:)=-1.0E+030_rp
  
     diffj_illum(:)=0.0_rp
     diffj_isInside(:)=0_ip

     counterInnerElems=0_ip
        
     !dsigmatmp = 96.7_rp  ! target  (sigma=100.0)  sigma_B=3.3
     !dsigmatmp = 46.7_rp  ! target  (sigma=50.0)  sigma_B=3.3
     !dsigmatmp = 21.7_rp  ! target  (sigma=25.0)  sigma_B=3.3
     !dsigmatmp = 6.7_rp   ! target  (sigma=10.0)  sigma_B=3.3
     !dsigmatmp = -2.1_rp  ! target  (sigma=1.2)  sigma_B=3.3
     !dsigmatmp = -2.2_rp  ! target  (sigma=1.1)  sigma_B=3.3
     !dsigmatmp = -2.8_rp  ! target  (sigma=0.5)  sigma_B=3.3
     !dsigmatmp = -3.2_rp  ! target  (sigma=0.1)  sigma_B=3.3
     !dsigmatmp = -3.225_rp  ! target  (sigma=0.075)  sigma_B=3.3
     !dsigmatmp = -3.229_rp  ! target  (sigma=0.071)  sigma_B=3.3
     dsigmatmp = -3.25_rp ! target  (sigma=0.05)  sigma_B=3.3
     !dsigmatmp = -3.26_rp ! target  (sigma=0.04)  sigma_B=3.3
     !dsigmatmp = -3.268377223_rp ! target  (sigma=0.0316227)  sigma_B=3.3
     !dsigmatmp = -3.27_rp ! target  (sigma=0.03)  sigma_B=3.3
     !dsigmatmp = -3.28_rp ! target  (sigma=0.02)  sigma_B=3.3
     !dsigmatmp = -3.29_rp ! target  (sigma=0.01)  sigma_B=3.3
     !dsigmatmp = -3.291_rp ! target  (sigma=0.009)  sigma_B=3.3
     !dsigmatmp = -3.292_rp ! target  (sigma=0.008)  sigma_B=3.3
     !dsigmatmp = -3.295_rp ! target  (sigma=0.005)  sigma_B=3.3
     !dsigmatmp = -3.299_rp ! target  (sigma=0.001)  sigma_B=3.3
     !dsigmatmp = -3.29982_rp ! target  (sigma=0.00018)  sigma_B=3.3
     !dsigmatmp = -3.2999_rp ! target  (sigma=0.0001)  sigma_B=3.3
     !dsigmatmp = -3.29999_rp ! target  (sigma=0.00001)  sigma_B=3.3
     !dsigmatmp = -3.299999_rp ! target  (sigma=0.000001)  sigma_B=3.3
     !dsigmatmp = -3.2999999_rp ! target  (sigma=0.0000001)  sigma_B=3.3

     elementsOpt: do ielem = 1,nelem

        centerX=0.0_rp
        centerY=0.0_rp
        centerZ=0.0_rp
        minX=1.0E+30_rp
        minY=1.0E+30_rp
        minZ=1.0E+30_rp
        maxX=-1.0E+30_rp
        maxY=-1.0E+30_rp
        maxZ=-1.0E+30_rp

        !print *,kfl_paral,': elementsOpt: ',ielem,' en elmope' 

        pelty = ltype(ielem)       !Element type	
        pnode = nnode(pelty)       !Number of element nodes 

        !Check material
        pmate = 1_ip
        if ( nmate > 1_ip ) then
           pmate = lmate(ielem)
        end if

        ! calculate center of tetra
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           centerX = centerX + coord(1, ipoin) 
           centerY = centerY + coord(2, ipoin) 
           centerZ = centerZ + coord(3, ipoin) 
           if(coord(1,ipoin)<minX)then
              minX=coord(1,ipoin)
           end if 
           if(coord(2,ipoin)<minY)then
              minY=coord(2,ipoin)
           end if 
           if(coord(3,ipoin)<minZ)then
              minZ=coord(3,ipoin)
           end if
           if(coord(1,ipoin)>maxX)then
              maxX=coord(1,ipoin)
           end if 
           if(coord(2,ipoin)>maxY)then
              maxY=coord(2,ipoin)
           end if 
           if(coord(3,ipoin)>maxZ)then
              maxZ=coord(3,ipoin)
           end if  
        end do
        centerX=centerX * 0.25_rp
        centerY=centerY * 0.25_rp
        centerZ=centerZ * 0.25_rp

        ! check if center of tetra is inside of anomaly
!        ! test A
!        if(     (centerX .lt. 3000.0_rp) .and. (centerX .gt. 1000.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 1000.0_rp) .and. (centerZ .gt. 500.0_rp) &
!        ) then
!        ! test B
!        if(     (centerX .lt. 2000.0_rp) .and. (centerX .gt. -2000.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 2000.0_rp) .and. (centerZ .gt. 1000.0_rp) &
!        ) then
!           do inode = 1,pnode
!              ipoin = lnods(inode,ielem)
!              design_vars(lninv_loc(ipoin))= log(dsigmatmp+bckco_hlm(1)) 
!           end do
!        ! test C/D
!        if(     (centerX .lt. 2000.0_rp) .and. (centerX .gt. -2000.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 2200.0_rp) .and. (centerZ .gt. 1500.0_rp) &
!        ) then
!        ! test E/F
!        if(     (centerX .lt. 1000.0_rp) .and. (centerX .gt. -1000.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 2500.0_rp) .and. (centerZ .gt. 2000.0_rp) &
!        ) then
!        ! test G
!        if(     (centerX .lt. 4000.0_rp) .and. (centerX .gt. 0.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 2500.0_rp) .and. (centerZ .gt. 2000.0_rp) &
!        ) then
!        ! test H
!        if(     (maxX .lt. 200.0_rp) .and. (minX .gt. -600.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1500.0_rp) .and. (minZ .gt. 1000.0_rp) &
!        ) then
!        ! test H2: starting model modified
!        if(     (maxX .lt. 600.0_rp) .and. (minX .gt. -600.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1500.0_rp) .and. (minZ .gt. 1000.0_rp) &
!        ) then
!        ! test H3: starting model modified
!        if(     (maxX .lt. 1200.0_rp) .and. (minX .gt. -1200.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1500.0_rp) .and. (minZ .gt. 1000.0_rp) &
!        ) then
!        ! test H4: starting model modified
!        if(     (maxX .lt. 2400.0_rp) .and. (minX .gt. -2400.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1500.0_rp) .and. (minZ .gt. 1000.0_rp) &
!        ) then
!        ! test H5: starting model modified
!        if(     (maxX .lt. 200.0_rp) .and. (minX .gt. -2400.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1500.0_rp) .and. (minZ .gt. 1000.0_rp) &
!        ) then
!        ! test H6: starting model modified
!        if(     (maxX .lt. 2400.0_rp) .and. (minX .gt. -200.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1500.0_rp) .and. (minZ .gt. 1000.0_rp) &
!        ) then
!        ! test I
!        if(     (maxX .lt. 200.0_rp) .and. (minX .gt. -600.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test I2: starting model modified
!        if(     (maxX .lt. 600.0_rp) .and. (minX .gt. -600.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test I3: starting model modified
!        if(     (maxX .lt. 1200.0_rp) .and. (minX .gt. -1200.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test I4: starting model modified
!        if(     (maxX .lt. 2400.0_rp) .and. (minX .gt. -2400.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test I5: starting model modified
!        if(     (maxX .lt. 200.0_rp) .and. (minX .gt. -2400.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test I6: starting model modified
!        if(     (maxX .lt. 2400.0_rp) .and. (minX .gt. -200.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test J
!        if(     (maxX .lt. 2200.0_rp) .and. (minX .gt. 1400.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test J2: starting model modified
!         if(     (maxX .lt. 600.0_rp) .and. (minX .gt. -600.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test J3: starting model modified
!        if(     (maxX .lt. 1200.0_rp) .and. (minX .gt. -1200.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test J4: starting model modified
!        if(     (maxX .lt. 2400.0_rp) .and. (minX .gt. -2400.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
!        ! test J5: starting model modified
!        if(     (maxX .lt. 200.0_rp) .and. (minX .gt. -2400.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 500.0_rp) &
!        ) then
!        ! test J6: starting model modified
!        if(     (maxX .lt. 2400.0_rp) .and. (minX .gt. -200.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
!        ) then
        ! test K: starting model modified
        ! 2K Inner (original)
        if(     (maxX .lt. 2400.0_rp) .and. (minX .gt. -2400.0_rp) .and. &
                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
                (maxZ .lt. 1200.0_rp) .and. (minZ .gt. 600.0_rp) &
        ! 100K Inner
        !if(     (maxX .lt. 5000.0_rp) .and. (minX .gt. -5000.0_rp) .and. &
        !        (maxY .lt. 1000.0_rp) .and. (minY .gt. -1000.0_rp) .and. & 
        !        (maxZ .lt. 3000.0_rp) .and. (minZ .gt. 0.0_rp) &
        !) then
        ! all
        !if(   1==1   &
        ) then


           do inode = 1,pnode
              ipoin = lnods(inode,ielem)

              if(diffj_isInside(lninv_loc(ipoin))==0_ip)then
                 design_vars(lninv_loc(ipoin))= log(dsigmatmp+bckco_hlm(1)) 
              end if

              diffj_illum(lninv_loc(ipoin)) = 1.0_rp * &
                 exp( -1.2_rp * ((3000_rp - abs(coord(3, ipoin)))/ (1452.8775_rp))**2.0_rp  )  ! skin depth
                 !exp( -1.8_rp * ((3000_rp - abs(coord(3, ipoin)))/ (1452.8775_rp))**2.0_rp  )  ! skin depth
                 !exp( -1.8_rp * ((2000.0_rp - abs(coord(3, ipoin)))/ (1452.8775_rp))**2.0_rp  )  ! skin depth


           end do


!           if( maxX < 0.0_rp ) then
!              do inode = 1,pnode
!                 ipoin = lnods(inode,ielem)
!                 design_vars(lninv_loc(ipoin))= log(dsigmatmp+bckco_hlm(1)) 
!                 diffj_illum(lninv_loc(ipoin)) = 1.0_rp !* &
!                 !   exp( -1.2_rp * ((3000.0_rp - abs(coord(3, ipoin)))/ (1452.8775_rp))**2.0_rp  )  ! skin depth
!              end do
!           else
!              do inode = 1,pnode
!                 ipoin = lnods(inode,ielem)
!                 design_vars(lninv_loc(ipoin))= log(0.1*(dsigmatmp+bckco_hlm(1))) 
!                 diffj_illum(lninv_loc(ipoin)) = 1.0_rp !* &
!                 !   exp( -1.2_rp * ((3000.0_rp - abs(coord(3, ipoin)))/ (1452.8775_rp))**2.0_rp  )  ! skin depth
!              end do
!           end if


           ! tag inside nodes
           counterInnerElems = counterInnerElems +1
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              diffj_isInside(lninv_loc(ipoin)) = 1_ip 
           end do


        !else
        !   do inode = 1,pnode
        !      ipoin = lnods(inode,ielem)
        !      !if(design_vars(lninv_loc(ipoin))==1.0E+030_rp)then
        !         design_vars(lninv_loc(ipoin))= log(dsigma_hlm(1,pmate)+bckco_hlm(1)) 
        !      !end if
        !      !diffj_illum(lninv_loc(ipoin)) = 0.0_rp
        !   end do
        end if







        ! reference models
        !! test J
        !if(     (maxX .lt. 2200.0_rp) .and. (minX .gt. 1400.0_rp) .and. &
        !        (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
        !        (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
        !) then
        ! test J perturbed
        !if(     (maxX .lt. 2200.0_rp) .and. (minX .gt. 1400.0_rp) .and. &
        !        (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
        !        (maxZ .lt. 1300.0_rp) .and. (minZ .gt. 800.0_rp) &
        !) then
        !   do inode = 1,pnode
        !      ipoin = lnods(inode,ielem)
        !      design_vars_ref(lninv_loc(ipoin))= log(dsigmatmp+bckco_hlm(1)) 
        !   end do

        !else
        !   do inode = 1,pnode
        !      ipoin = lnods(inode,ielem)
        !      if(design_vars_ref(lninv_loc(ipoin))==-1.0E+030_rp)then
        !         design_vars_ref(lninv_loc(ipoin))= log(dsigma_hlm(1,pmate)+bckco_hlm(1)) 
        !      end if
        !   end do
        !end if














        ! check if center of tetra is inside of interest window
!!        if(     (centerX .lt. 6000.0_rp) .and. (centerX .gt. -6000.0_rp) .and. &
!        if(     (centerX .lt. 2500.0_rp) .and. (centerX .gt. -2500.0_rp) .and. &
!        !if(     (centerX .lt. 400.0_rp) .and. (centerX .gt. -400.0_rp) .and. &
!        !if(     (centerY .lt. 4000.0_rp) .and. (centerY .gt. -4000.0_rp) .and. & 
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 2100.0_rp) .and. (centerZ .gt. 100.0_rp) &
!                !(centerZ .lt. 500.0_rp) .and. (centerZ .gt. 100.0_rp) &
!!                (centerZ .lt. 2000.0_rp) .and. (centerZ .gt. 1000.0_rp) &
!!                (centerZ .lt. 3000.0_rp) .and. (centerZ .gt. 1500.0_rp) &
!        ) then
!
!        ! test A
!        if(     (centerX .lt. 3000.0_rp) .and. (centerX .gt. 1000.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 1000.0_rp) .and. (centerZ .gt. 500.0_rp) &
!        ) then
        ! test B
!        if(     (centerX .lt. 2000.0_rp) .and. (centerX .gt. -2000.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 2000.0_rp) .and. (centerZ .gt. 1000.0_rp) &
!        ) then
        ! test C/D
!        if(     (centerX .lt. 2000.0_rp) .and. (centerX .gt. -2000.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 2200.0_rp) .and. (centerZ .gt. 1500.0_rp) &
!        ) then 
        ! test E/F
!        if(     (centerX .lt. 1000.0_rp) .and. (centerX .gt. -1000.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 2500.0_rp) .and. (centerZ .gt. 2000.0_rp) &
!        ) then
!        ! test G
!        if(     (centerX .lt. 4000.0_rp) .and. (centerX .gt. 0.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 2500.0_rp) .and. (centerZ .gt. 2000.0_rp) &
!        ) then
        ! test H
!        if(     (centerX .lt. 200.0_rp) .and. (centerX .gt. -600.0_rp) .and. &
!                (centerY .lt. 400.0_rp) .and. (centerY .gt. -400.0_rp) .and. & 
!                (centerZ .lt. 1500.0_rp) .and. (centerZ .gt. 1000.0_rp) &
!        ) then

!        if(     (maxX .lt. 2400.0_rp) .and. (minX .gt. -2400.0_rp) .and. &
!                (maxY .lt. 400.0_rp) .and. (minY .gt. -400.0_rp) .and. & 
!                (maxZ .lt. 1200.0_rp) .and. (minZ .gt. 600.0_rp) &
!        ) then
!                counterInnerElems = counterInnerElems +1
!                do inode = 1,pnode
!                        ipoin = lnods(inode,ielem)
!                        diffj_isInside(lninv_loc(ipoin)) = 1_ip 
!                end do
!        end if

     end do elementsOpt

     !write(*,*)kfl_paral,' counterInnerElems=',counterInnerElems  

  end if

end subroutine hlm_elmopefirst

