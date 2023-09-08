subroutine platri()
  !-----------------------------------------------------------------------
  !****f* Domain/platri
  ! NAME
  !    platri
  ! DESCRIPTION
  !    Inrestect a plane with an element
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_master
  use def_domain
  use mod_memchk
  use def_parame

  implicit none
  integer(ip)  :: ipoin,jpoin,izdom,ifoun,n,ii,ppoin
  integer(ip)  :: inode,pnode,ielem
  real(rp)     :: xp1(ndime),xp2(ndime)
  real(rp)     :: a,b,c,d,t


  if( INOTMASTER ) then
     if( ndime == 3 ) then
        !
        ! Equation of the plane: a x + b y + c z + d = 0
        !
        a = pafil(1)
        b = pafil(2)
        c = pafil(3)
        d = pafil(4)
        n = int(pafil(5),ip)  ! Depth of the intersection


        do ipoin = 1,npoin
           izdom = r_dom(ipoin)-1
           do while( izdom < r_dom(ipoin+1)-1 )

              izdom = izdom + 1
              jpoin = c_dom(izdom)
              ifoun = 0

              if( jpoin < ipoin ) then

                 xp1(1) = coord(1,ipoin)
                 xp1(2) = coord(2,ipoin)
                 xp1(3) = coord(3,ipoin)          
                 xp2(1) = coord(1,jpoin)
                 xp2(2) = coord(2,jpoin)
                 xp2(3) = coord(3,jpoin)   
                 !
                 ! Scalar product t = n.(P1,P2)
                 !
                 t = a*(xp2(1)-xp1(1)) + b*(xp2(2)-xp1(2)) + c*(xp2(3)-xp1(3))  

                 if( t /= 0.0_rp ) then
                    !
                    ! Compute parametric coordinate t on P1-P2
                    !
                    t = (- d - a*xp1(1) - b*xp1(2) - c*xp1(3) ) / t
                    if( t >= 0.0_rp .and. t <= 1.0_rp ) ifoun = 1

                 else if( abs( a * xp1(1) + b * xp1(2) + c * xp1(3) + d ) == 0.0_rp ) then
                    !
                    ! (P1,P2) is parallel to plane: check if P1 on plane
                    !
                    ifoun = 1
                 end if

                 if( ifoun == 1 ) then
                    !
                    ! Mark IPOIN and JPOIN
                    !
                    gefil(ipoin) = 1
                    gefil(jpoin) = 1
                    izdom        = r_dom(ipoin+1)-1                    
                 end if

              end if
           end do
        end do
        !
        ! Increase depth
        !
        do ii=2,n
           do ipoin = 1,npoin
              izdom = r_dom(ipoin)-1
              if( gefil(ipoin) <= ii-1 .and. gefil(ipoin) > 0 ) then
                 do while( izdom < r_dom(ipoin+1)-1 )              
                    izdom = izdom + 1
                    jpoin = c_dom(izdom)
                    if( gefil(jpoin) <= ii-1 ) then
                       gefil(jpoin) = ii
                    end if
                 end do
              end if
           end do
        end do

     end if
     !
     ! Exchange in Parallel
     !
     call parari('SLX',NPOIN_TYPE,npoin,gefil)
     do ipoin = 1,npoin
        gefil(ipoin) = min(1_ip,gefil(ipoin))
     end do
     !
     ! Mark nodes who's elements are marked...
     !
     do ielem = 1,nelem
        inode = 0
        pnode = nnode(ltype(ielem)) 
        do while( inode < pnode )
           inode = inode + 1
           ipoin = lnods(inode,ielem)
           if( gefil(ipoin) == 1 ) inode = pnode + 1
        end do
        if( inode == pnode + 1 ) then
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              if( gefil(ipoin) == 0 ) gefil(ipoin) = 2
           end do
        end if
     end do
     !
     ! Order GEFIL
     !
     ppoin = 0
     do ipoin = 1,npoin
        if( gefil(ipoin) > 0_ip ) then
           ppoin = ppoin + 1
           gefil(ipoin) = ppoin 
        end if
     end do
  end if

end subroutine platri

subroutine noseplan(a,b,c,d)
  !-----------------------------------------------------------------------
  !****f* Domain/noseplan
  ! NAME
  !    platri
  ! DESCRIPTION
  !    Inrestect a plane with an element
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_master
  use def_domain
  use mod_memchk
  use def_parame

  implicit none
  real(rp), intent(in) :: a,b,c,d
  integer(ip)  :: ipoin,jpoin,izdom,ifoun,n,ii,ppoin
  integer(ip)  :: inode,pnode,ielem
  real(rp)     :: xp1(ndime),xp2(ndime),t
  

  if( INOTMASTER ) then
     if( ndime == 3 ) then
        !
        ! Equation of the plane: a x + b y + c z + d = 0
        !
     
        n = 1  ! Depth of the intersection


        do ipoin = 1,npoin
           izdom = r_dom(ipoin)-1
           do while( izdom < r_dom(ipoin+1)-1 )

              izdom = izdom + 1
              jpoin = c_dom(izdom)
              ifoun = 0

              if( jpoin < ipoin ) then

                 xp1(1) = coord(1,ipoin)
                 xp1(2) = coord(2,ipoin)
                 xp1(3) = coord(3,ipoin)          
                 xp2(1) = coord(1,jpoin)
                 xp2(2) = coord(2,jpoin)
                 xp2(3) = coord(3,jpoin)   
                 !
                 ! Scalar product t = n.(P1,P2)
                 !
                 t = a*(xp2(1)-xp1(1)) + b*(xp2(2)-xp1(2)) + c*(xp2(3)-xp1(3))  

                 if( t /= 0.0_rp ) then
                    !
                    ! Compute parametric coordinate t on P1-P2
                    !
                    t = (- d - a*xp1(1) - b*xp1(2) - c*xp1(3) ) / t
                    if( t >= 0.0_rp .and. t <= 1.0_rp ) ifoun = 1

                 else if( abs( a * xp1(1) + b * xp1(2) + c * xp1(3) + d ) == 0.0_rp ) then
                    !
                    ! (P1,P2) is parallel to plane: check if P1 on plane
                    !
                    ifoun = 1
                 end if

                 if( ifoun == 1 .and.                                                 &
                      (((coord(1,ipoin)<0.116_rp .and. coord(1,ipoin)>0.0752_rp).and. &
                      (coord(2,ipoin)<0.105_rp .and. coord(2,ipoin)>0.0146_rp).and.  &
                      (coord(3,ipoin)<0.0539_rp .and. coord(3,ipoin)>0.00114_rp))    &
                      .or.                                                           &
                      ((coord(1,ipoin)<0.595_rp .and. coord(1,ipoin)>-0.396_rp).and. &
                      (coord(2,ipoin)<0.326_rp .and. coord(2,ipoin)>0.105_rp).and.    &
                      (coord(3,ipoin)<0.041_rp .and. coord(3,ipoin)>-0.488_rp)))     &
                      ) then
                    !
                    ! Mark IPOIN and JPOIN
                    !
                    gefil(ipoin) = 1
                    gefil(jpoin) = 1
                    izdom        = r_dom(ipoin+1)-1                    
                 end if

              end if
           end do
        end do
        !
        ! Increase depth
        !
        do ii=2,n
           do ipoin = 1,npoin
              izdom = r_dom(ipoin)-1
              if( gefil(ipoin) <= ii-1 .and. gefil(ipoin) > 0 ) then
                 do while( izdom < r_dom(ipoin+1)-1 )              
                    izdom = izdom + 1
                    jpoin = c_dom(izdom)
                    if( gefil(jpoin) <= ii-1 ) then
                       gefil(jpoin) = ii
                    end if
                 end do
              end if
           end do
        end do

     end if
  end if

end subroutine noseplan

subroutine noseplanclean()
  !-----------------------------------------------------------------------
  !****f* Domain/noseplanclean
  ! NAME
  !    platri
  ! DESCRIPTION
  !    clean the gefil for the seval plan for the nose
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_master
  use def_domain
  use mod_memchk
  use def_parame

  implicit none
  integer(ip)  :: ipoin,jpoin,izdom,ifoun,n,ii,ppoin
  integer(ip)  :: inode,pnode,ielem
  real(rp)     :: xp1(ndime),xp2(ndime),t


  if( INOTMASTER ) then
     !
     ! Exchange in Parallel
     !
     call parari('SLX',NPOIN_TYPE,npoin,gefil)
     do ipoin = 1,npoin
        gefil(ipoin) = min(1_ip,gefil(ipoin))
     end do
     !
     ! Mark nodes who's elements are marked...
     !
     do ielem = 1,nelem
        inode = 0
        pnode = nnode(ltype(ielem)) 
        do while( inode < pnode )
           inode = inode + 1
           ipoin = lnods(inode,ielem)
           if( gefil(ipoin) == 1 ) inode = pnode + 1
        end do
        if( inode == pnode + 1 ) then
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              if( gefil(ipoin) == 0 ) gefil(ipoin) = 2
           end do
        end if
     end do
     !
     ! Order GEFIL
     !
     ppoin = 0
     do ipoin = 1,npoin
        if( gefil(ipoin) > 0_ip ) then
           ppoin = ppoin + 1
           gefil(ipoin) = ppoin 
        end if
     end do
  end if

end subroutine noseplanclean
