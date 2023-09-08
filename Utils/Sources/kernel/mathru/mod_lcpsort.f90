! Module for sorting and merging ASCII text strings
! Uses standard Mergesort, with two enhancements:
!
! A.(i) Along with each string the number of matching leading characters
!       compared to the preceding string is recorded and updated.
!   (ii) A full string comparison is needed only at these times:
!          (a) Initially
!          (b) When merging, and the two LCPL
!              (longest commmon prefix length)-s are equal
!          (c) When comparing the head of the lower part with the
!              tail of the upper part
!
!   (iii) At all other times, only the LCPL-s  of the
!         two strings are compared and the string associated with the
!         smaller LCPL is output before the other.
!
! B. If the head of the lower list is already in order with respect to the
!    tail of the upper list, merging is skipped. The lower list is copied
!    into scratch space and the output list overlays the composite list.
!
! Ref. See Ph.D. dissertation:
!    dspace.wul.waseda.ac.jp/dspace/bitstream/2065/28672/4/Honbun-4624_01.pdf
!
! Programmed 5/2009 by: N. Shamsundar, University of Houston, shamsundar@uh.edu
!
! Can be used for sorting text files with fixed length or variable length records
! For variable length records, make LINELEN equal to the longest possible
! record length, and change the file I/O in the example main program
! to perform variable length record I/O.
!
module mod_lcpsort

  use def_kintyp, only : ip
  implicit none

  private
  
  integer(ip),      parameter :: LINELEN=20          ! maximum length of strings
  character(len=1), parameter :: EOL=achar(10)

  type :: lptype
     integer(ip)                :: lcpl
     character(LINELEN),pointer :: str
  end type lptype

  integer(ip) :: nmerge,ncomp,nskip

  public :: lptype
  public :: lcpsort_sort
  public :: lcpsort_initialization
  public :: lcpsort_deallocate
  public :: lcpsort_concatenate
  
contains

  subroutine lcpsort_concatenate(X,n,nout)
    implicit none
    integer(ip),                intent(in)    :: n
    integer(ip),                intent(inout) :: nout
    type (lptype), allocatable, intent(inout) :: X(:)
    integer(ip)                               :: ii,jj

    nout = 1
    do ii = 2,n
      if( trim(X(ii)%str) /= trim(X(ii-1)%str) ) then
         nout = nout + 1
         X(nout)%str = X(ii)%str
       end if
    end do
    
  end subroutine lcpsort_concatenate
  
  subroutine lcpsort_initialization(X,n,Z)
    implicit none
    integer(ip),                intent(in)    :: n
    type (lptype), allocatable, intent(inout) :: X(:),Z(:)
    integer(ip)                               :: ii
    
    allocate(X(n))
    allocate(Z((n+1)/2))
    do ii = 1,n
       X(ii) % lcpl =  0
       nullify(X(ii) % str)
    end do
    
  end subroutine lcpsort_initialization
  
  subroutine lcpsort_deallocate(X,Z)
    implicit none
    type (lptype), allocatable, intent(inout) :: X(:),Z(:)
    
    deallocate(X)
    deallocate(Z)
    
  end subroutine lcpsort_deallocate
  !
  ! ***  SUBROUTINE MERGELCP  ***
  !
  ! Input: X(1:lX) and Y(1:lY) contain (i) arrays of derived type lptype,
  !        of sizes lX and lY, to be merged. X and Y may not overlap.
  ! Output: Z(1:lX+lY) will contained the merged array, with the
  !        lcpl components updated when string pointers are swapped.
  !        Z(lX+1::lx+lY) may overlap Y.
  !
  subroutine lcpsort_merge(X,lX,Y,lY,Z)
    implicit none
    integer(ip),  intent (in)     :: lX,lY
    type(lptype), intent (in out) :: X(lX),Y(lY),Z(lX+lY)
    integer(ip)                   :: i,j,k,ix

    logical :: debug=.FALSE.

    if(debug)then
       write(*,*)' List X'
       write(*,'(1x,I3,2x,A)')(X(i)%lcpl,X(i)%str,i=1,lX)
       write(*,*)
       write(*,*)' List Y'
       write(*,'(1x,I3,2x,A)')(Y(i)%lcpl,Y(i)%str,i=1,lY)
       write(*,*)
    end if

    i=1
    j=1
    k=1
    do while (i <= lX .and. j <= lY)
       if(X(i)%lcpl /= Y(j)%lcpl)then
          if(X(i)%lcpl > Y(j)%lcpl)then
             Z(k)=X(i)
             i=i+1
          else
             Z(k)=Y(j)
             j=j+1
          end if
          k=k+1
       else
          ix=X(i)%lcpl+1
          ncomp=ncomp+1
          do while(X(i)%str(ix:ix) == Y(j)%str(ix:ix) .and. &
               X(i)%str(ix:ix) /= EOL)
             ix=ix+1
          end do
          if(X(i)%str(ix:ix) > Y(j)%str(ix:ix))then
             X(i)%lcpl=ix-1
             Z(k)=Y(j)
             j=j+1
          else
             Y(j)%lcpl=ix-1
             Z(k)=X(i)
             i=i+1
          end if
          k=k+1
       end if
    end do
    do while(i <= lX)
       Z(k)=X(i)
       i=i+1
       k=k+1
    end do

    if(debug)then
       write(*,*)' Merged List'
       write(*,'(1x,I3,2x,A)')(Z(i)%lcpl,Z(i)%str,i=1,lX+lY)
       write(*,*)
    end if

    return
  end subroutine Lcpsort_merge

  ! ***   SUBROUTINE SortLCP   ***
  !
  ! Sorts an array of Ascii text strings using Longest Common Prefix
  ! Merge Sort. The parameters are:
  ! The array X(1:n) of derived type lptype (see top of file),
  ! with pointers to the strings to be sorted in the %str
  ! components, and with the %lcpl components initialized to zero.

  ! The array Z is used as scratch space, and must be at least half as long
  ! as X. Since the routine is recursive, the design decision was taken that
  ! allocation of Z was to be performed before calling the subroutine.
  !
  recursive subroutine lcpsort_sort(X,n,Z)
    implicit none
    integer(ip),   intent(in)    :: n
    type (lptype), intent(inout) :: X(n),Z((n+1)/2)
    type (lptype)                :: TMP
    integer(ip)                  :: i,m
    
    if( n <  2 ) return                  ! no need to sort array of size 1
    if( n == 2 ) then                    ! special handling for pairs
       ncomp=ncomp+1
       i=1
       do while(X(1)%str(i:i) == X(2)%str(i:i) .and. &
            X(1)%str(i:i) /= EOL)
          i=i+1
          if(i>LINELEN) then
             i = LINELEN
             exit
          end if
       end do
       if(X(1)%str(i:i) > X(2)%str(i:i))then
          TMP=X(2)
          X(2)=X(1)
          X(1)=TMP
       end if
       X(2)%lcpl=i-1
       return
    end if
    m=(n+1)/2                              ! divide and conquer for n > 2
    call lcpsort_sort(X,m,Z)               ! sort lower part
    call lcpsort_sort(X(m+1),n-m,Z)        ! sort upper part
    ncomp=ncomp+1
    i=1                                    ! compare head of lower with tail of upper
    do while(X(m)%str(i:i) == X(m+1)%str(i:i) .and. &
         X(m)%str(i:i) /= EOL)
       i=i+1
       if(i>LINELEN) then
          i = LINELEN
          exit
       end if
    end do
    if(X(m)%str(i:i) > X(m+1)%str(i:i))then
       DO I=1,m
          Z(I)=X(i)                         ! copy "left" array into scratch space
       END DO
       call lcpsort_merge(Z,m,X(m+1),n-m,X) ! merge the two parts
       nmerge=nmerge+1
    else
       nskip=nskip+1
       X(m+1)%lcpl=i-1                      ! skip merging, but update lcpl
    end if
    
  end subroutine lcpsort_sort
  
end module mod_lcpsort
