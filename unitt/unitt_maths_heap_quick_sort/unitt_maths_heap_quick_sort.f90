program unitt_maths_heap_quick_sort
  !
  ! Sort reals and integers using quick and heap sort algorithms
  !
  use def_kintyp_basic, only : ip, rp
  use mod_maths_sort,   only : maths_heap_sort
  use mod_maths_sort,   only : maths_quick_sort

  implicit none

  integer(ip)           :: nn,ii,itask,imeth,ierro,jj,iimax
  real(rp),   pointer   :: aa(:)
  integer(ip),pointer   :: kk(:)
  real(rp)              :: time1,time2
  character(10)         :: method(2)
  character(10)         :: realint(2)
  character(10)         :: order(2)

  integer                  :: nums
  integer,     allocatable :: seed(:)
  integer                  :: values(8)

  call date_and_time(VALUES=values)
  call random_seed(size = nums)
  allocate(seed(nums))
  iimax = min(nums,size(values))
  seed(1:iimax) = values(1:iimax)
  call random_seed(put=seed)

  order(1)   = 'DECREASING'
  order(2)   = 'INCREASING'
  method(1)  = 'QUICK SORT'
  method(2)  = 'HEAP  SORT'
  realint(1) = 'REAL'
  realint(2) = 'INT.'

  nn = 600000
  allocate(aa(nn))
  allocate(kk(nn))
  aa = 0.0_rp
  kk = 0_ip

  do itask = 1,2
     !
     ! ITASK = 1 ... Decreasing order
     ! ITASK = 2 ... Increasing order
     !
     ierro = 0
     do imeth = 1,2
        do jj = 1,2
           call cpu_time(time1)
           do ii = 1,nn
              call RANDOM_NUMBER(aa(ii))     
              kk(ii) = -nn/2+FLOOR(nn*aa(ii))
           end do
           if( imeth == 1 ) then
              if( jj == 1 ) then
                 call maths_quick_sort(itask,nn,aa)
              else
                 call maths_quick_sort(itask,nn,kk)
              end if
           else
              if( jj == 1 ) then
                 call maths_heap_sort(itask,nn,aa)
              else
                 call maths_heap_sort(itask,nn,kk)
              end if
           end if
           call cpu_time(time2)
           print*,trim(method(imeth))//', '//trim(order(ITASK))//', '//trim(realint(jj))//', ','TIME= ',time2-time1
           if( jj == 1 ) then
              if( itask == 2 ) then
                 do ii = 1,nn-1
                    if( aa(ii) > aa(ii+1) ) then
                       ierro = 1
                       goto 1
                    end if
                 end do
              else
                 do ii = 1,nn-1
                    if( aa(ii) < aa(ii+1) ) then
                       ierro = 1
                       goto 1
                    end if
                 end do
              end if
           else
              if( itask == 2 ) then
                 do ii = 1,nn-1
                    if( kk(ii) > kk(ii+1) ) then
                       ierro = 1
                       goto 1
                    end if
                 end do
              else
                 do ii = 1,nn-1
                    if( kk(ii) < kk(ii+1) ) then
                       ierro = 1
                       goto 1
                    end if
                 end do
              end if
           end if
        end do
     end do

  end do

  deallocate(aa,kk)
  stop

1 print*,'Error: ',imeth,itask
  stop 1

end program unitt_maths_heap_quick_sort
