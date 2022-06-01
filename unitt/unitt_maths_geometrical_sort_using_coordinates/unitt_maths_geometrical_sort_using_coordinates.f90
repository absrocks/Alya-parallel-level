program unitt_maths_geometrical_sort_using_coordinates 
  
  use def_kintyp_basic, only : ip, rp
  use mod_maths_sort,   only : maths_geometrical_sort_using_coordinates 
  
  implicit none
  real(rp)                 :: rr
  real(rp),    parameter   :: epsil=epsilon(1.0_8)
  integer(ip)              :: num_solutions,ii,iimax,itask,idime,ndime
  integer(ip)              :: numb
  integer(ip), allocatable :: ivin(:)
  real(rp),    allocatable :: xx(:,:)
  real(rp)                 :: rval(3)
  integer                  :: nums
  integer,     allocatable :: seed(:)
  integer                  :: values(8)
  
  call date_and_time(VALUES=values)
  call random_seed(size = nums)
  allocate(seed(nums))
  iimax = min(nums,size(values))
  seed(1:iimax) = values(1:iimax)
  call random_seed(put=seed)
  deallocate(seed)

  ndime=3
  numb=10000
  allocate(xx(ndime,numb))
  allocate(ivin(numb))

  do ii = 1,numb     
     call RANDOM_NUMBER(rval)     
     xx(:,ii) = rval(:)
  end do
  !
  ! Decreasing order
  !
  itask = 1
  call maths_geometrical_sort_using_coordinates(itask,3_ip,numb,xx,ivin,TOLERANCE=epsil)
  do ii = 1,numb-1
     if( xx(1,ii) > xx(1,ii+1) ) then
        continue
     else if( abs(xx(1,ii) - xx(1,ii+1)) <= epsil ) then
        if( xx(2,ii) > xx(2,ii+1) ) then
           continue
        else if( abs(xx(2,ii) - xx(2,ii+1)) <= epsil ) then
           if( xx(3,ii) > xx(3,ii+1) ) then
              continue
           else
              print*,'a=',xx(1:ndime,ii)
              print*,'b=',xx(1:ndime,ii+1)
              stop 1
           end if
        else
           print*,'a=',xx(1:ndime,ii)
           print*,'b=',xx(1:ndime,ii+1)
           stop 1
        end if
     else
        print*,'a=',xx(1:ndime,ii)
        print*,'b=',xx(1:ndime,ii+1)
        stop 1
     end if
  end do
  !
  ! Increasing order
  !
  itask = 2
  call maths_geometrical_sort_using_coordinates(itask,3_ip,numb,xx,ivin,TOLERANCE=epsil)
  do ii = 1,numb-1
     if( xx(1,ii) < xx(1,ii+1) ) then
        continue
     else if( abs(xx(1,ii) - xx(1,ii+1)) <= epsil ) then
        if( xx(2,ii) < xx(2,ii+1) ) then
           continue
        else if( abs(xx(2,ii) - xx(2,ii+1)) <= epsil ) then
           if( xx(3,ii) < xx(3,ii+1) ) then
              continue
           else
              print*,'a=',xx(1:ndime,ii)
              print*,'b=',xx(1:ndime,ii+1)
              stop 1
           end if
        else
           print*,'a=',xx(1:ndime,ii)
           print*,'b=',xx(1:ndime,ii+1)
           stop 1
        end if
     else
        print*,'a=',xx(1:ndime,ii)
        print*,'b=',xx(1:ndime,ii+1)
        stop 1
     end if
  end do

  !
  ! Simple test
  !
  !>          Unsorted list:
  !>          node      x     y     z
  !>             1   -3.0  -2.0   3.0
  !>             2   -3.0  -2.0   4.0
  !>             3    0.0  -2.0  -7.0
  !>             4    0.0  -2.0   2.0
  !>             5   -1.0  -2.0  -1.0
  !>             6   -1.0   1.0   2.0
  !>             7   -3.0   2.0  12.0
  !>             8    0.0  -2.0   8.0
  !>             9    0.0  -2.0   9.0
  !>
  !>          Sorted list:
  !>          node      x     y     z
  !>             1   -3.0  -2.0   3.0
  !>             2   -3.0  -2.0   4.0
  !>             7   -3.0   2.0  12.0
  !>             5   -1.0  -2.0  -1.0
  !>             6   -1.0   1.0   2.0
  !>             3    0.0  -2.0  -7.0
  !>             4    0.0  -2.0   2.0
  !>             8    0.0  -2.0   8.0
  !>             9    0.0  -2.0   9.0
  !
  1 continue
  numb = 9
  do ii = 1,9
     ivin(ii) = ii
  end do
  xx(1:3,1) = (/ -3.0_rp, -2.0_rp,  3.0_rp /)
  xx(1:3,2) = (/ -3.0_rp, -2.0_rp,  4.0_rp /)
  xx(1:3,3) = (/  0.0_rp, -2.0_rp, -7.0_rp /)
  xx(1:3,4) = (/  0.0_rp, -2.0_rp,  2.0_rp /)
  xx(1:3,5) = (/ -1.0_rp, -2.0_rp, -1.0_rp /)
  xx(1:3,6) = (/ -1.0_rp,  1.0_rp,  2.0_rp /)
  xx(1:3,7) = (/ -3.0_rp,  2.0_rp, 12.0_rp /)
  xx(1:3,8) = (/  0.0_rp, -2.0_rp,  8.0_rp /)
  xx(1:3,9) = (/  0.0_rp, -2.0_rp,  9.0_rp /)
  
  call maths_geometrical_sort_using_coordinates(2_ip,3_ip,numb,xx,ivin,TOLERANCE=0.0_rp)

  do ii = 1,9
     print*,ivin(ii),xx(:,ii)
  end do
  
  if(      ivin(1) == 1 .and. &
       &   ivin(2) == 2 .and. &   
       &   ivin(3) == 7 .and. &   
       &   ivin(4) == 5 .and. &   
       &   ivin(5) == 6 .and. &   
       &   ivin(6) == 3 .and. &   
       &   ivin(7) == 4 .and. &  
       &   ivin(8) == 8 .and. &  
       &   ivin(9) == 9 ) then
     continue
  else
     print*,'Simple test failed'
     stop 1
  end if
  
  deallocate(xx)

end program unitt_maths_geometrical_sort_using_coordinates

