module mod_sort 

contains

  subroutine sort(criter,nleni)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)            :: nleni 
    integer(ip),intent(inout)         :: criter(nleni) 
    integer(ip)                       :: l,r,i,j,n,crit
    !
    !     This subroutine sorts with a heap sort the array criter in place
    !
    !     On input/output:   criter(nleni)
    !
    !
    n=nleni
    if(n<2)return 
    l= n/2 +1
    r=n

    do 

       if(l<=1) then
          crit=criter(r)
          criter(r)=criter(1)
          r=r-1
          if(r==1)then
             criter(1)=crit
             exit
          endif
       else
          l=l-1
          crit=criter(l)
       endif

       j=l

       do 
          i=j
          j=2*j
          if(j<r) then
             if(criter(j)<criter(j+1))then
                j=j+1
             endif
             if(crit>=criter(j)) then
                criter(i)=crit
                exit
             else   
                criter(i)=criter(j)
             endif
          else if(j==r) then
             if(crit>=criter(j)) then
                criter(i)=crit
                exit
             else  
                criter(i)=criter(j)
             endif
          else 
             criter(i)=crit
             exit
          endif
       enddo
    enddo


  end subroutine sort

  subroutine sort2(array,criter,nleni)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)            :: nleni 
    integer(ip),intent(inout)         :: criter(nleni),array(nleni) 
    integer(ip)                       :: l,r,i,j,n,crit,icrit
    !
    !     This subroutine sorts with a heap sort the array criter in place
    !     and modifies array accordingly
    !
    !     On input/output:   criter(nleni)
    !                        array(nleni)
    !
    n=nleni
    if(n<2)return 
    l= n/2 +1
    r=n

    do 

       if(l<=1) then
          crit=criter(r)
          icrit=array(r)
          criter(r)=criter(1)
          array(r)=array(1) 
          r=r-1
          if(r==1)then
             criter(1)=crit
             array(1)=icrit
             exit
          endif
       else
          l=l-1
          crit=criter(l)
          icrit=array(l)
       endif

       j=l

       do 
          i=j
          j=2*j
          if(j<r) then
             if(criter(j)<criter(j+1))then
                j=j+1
             endif
             if(crit>=criter(j)) then
                criter(i)=crit
                array(i)=icrit
                exit
             else   
                criter(i)=criter(j)
                array(i)=array(j) 
             endif
          else if(j==r) then
             if(crit>=criter(j)) then
                criter(i)=crit
                array(i)=icrit 
                exit
             else  
                criter(i)=criter(j)
                array(i)=array(j) 
             endif
          else 
             criter(i)=crit
             array(i)=icrit
             exit
          endif
       enddo
    enddo
  
  end subroutine sort2 

  subroutine sort3(array,criter,nleni)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)            :: nleni 
    integer(ip),intent(inout)         :: array(nleni) 
    real(rp),intent(inout)            :: criter(nleni) 
    real(rp)                          :: crit
    integer(ip)                       :: l,r,i,j,n,icrit
    !
    !     This subroutine sorts with a heap sort the array criter in place
    !     and modifies array accordingly
    !
    n=nleni
    if(n<2)return 
    l= n/2 +1
    r=n

    do 

       if(l<=1) then
          crit=criter(r)
          icrit=array(r)
          criter(r)=criter(1)
          array(r)=array(1) 
          r=r-1
          if(r==1)then
             criter(1)=crit
             array(1)=icrit
             exit
          endif
       else
          l=l-1
          crit=criter(l)
          icrit=array(l)
       endif

       j=l

       do 
          i=j
          j=2*j
          if(j<r) then
             if(criter(j)<criter(j+1))then
                j=j+1
             endif
             if(crit>=criter(j)) then
                criter(i)=crit
                array(i)=icrit
                exit
             else   
                criter(i)=criter(j)
                array(i)=array(j) 
             endif
          else if(j==r) then
             if(crit>=criter(j)) then
                criter(i)=crit
                array(i)=icrit 
                exit
             else  
                criter(i)=criter(j)
                array(i)=array(j) 
             endif
          else 
             criter(i)=crit
             array(i)=icrit
             exit
          endif
       enddo
    enddo
  
  end subroutine sort3 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !     This is the interface for the heap sort for the 
  !     advancing front like mesh generator
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gtsmall(lheap,nheap,rheap,ismall,lfront,nfront)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(inout)     :: nheap
    integer(ip),intent(in)        :: nfront
    integer(ip),intent(inout)     :: lheap(nheap),ismall,lfront(3,nfront)
    real(rp),intent(in)           :: rheap(nfront)
    integer(ip)                   :: ipson1,ipson2,ipfath,ieson1,ieson2,iefath
    real(rp)                      :: rtol
    ! 
    rtol=1.05d+00
    !
    !     Remove the smallest element
    !
    ismall=lheap(1)
    !
    !     Restore the tree structure
    !
    !     The last element will become the first one. For grids with the same
    !     size, we want the last element to stay the last one for grids to 
    !     be regular --> the test is done with a tolerance 
    !
    lheap(1)=lheap(nheap)
    lfront(3,lheap(1))=1_ip 
    nheap=nheap-1

    ipfath=1
    ipson1=2*ipfath
    ipson2=ipson1+1
    do while(nheap>=ipson1)
       if(nheap>=ipson2)then
          ieson1=lheap(ipson1)
          ieson2=lheap(ipson2)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(ieson2))then
             if(rheap(ieson1)<rheap(iefath)*rtol)then
                lheap(ipfath)=ieson1
                lheap(ipson1)=iefath
                lfront(3,iefath)=ipson1
                lfront(3,ieson1)=ipfath
                ipfath=ipson1
             else
                exit
             endif
          else
             if(rheap(ieson2)<rheap(iefath)*rtol)then
                lheap(ipfath)=ieson2
                lheap(ipson2)=iefath
                lfront(3,iefath)=ipson2
                lfront(3,ieson2)=ipfath
                ipfath=ipson2

             else
                exit
             endif
          endif

       else
          ieson1=lheap(ipson1)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(iefath)*rtol)then
             lheap(ipfath)=ieson1
             lheap(ipson1)=iefath
             lfront(3,iefath)=ipson1
             lfront(3,ieson1)=ipfath
             ipfath=ipson1
          else
             exit
          endif
       endif
       ipson1=2*ipfath
       ipson2=ipson1+1
    enddo

    !
    !     DBG
    !

!!$  do iheap=1,nheap
!!$     ipfath=iheap
!!$     ipson1=2*ipfath
!!$     ipson2=2*ipfath+1
!!$     iefath=lheap(ipfath)
!!$
!!$     if(ipson1<nheap)then
!!$        ieson1=lheap(ipson1)
!!$        if(rheap(iefath)>rheap(ieson1))then
!!$           write(*,*)'error 1 in gtsmall'
!!$           stop
!!$        endif
!!$     endif
!!$
!!$     if(ipson2<nheap)then
!!$        ieson2=lheap(ipson2)
!!$        if(rheap(iefath)>rheap(ieson2))then
!!$           write(*,*)'error 2 in gtsmall'
!!$           stop
!!$        endif
!!$     endif
!!$  enddo

  end subroutine gtsmall

  subroutine addheap(lheap,nheap,rheap,ie,lfront,nfront)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)        :: ie,nfront
    integer(ip),intent(inout)     :: nheap,lheap(nheap+1),lfront(3,nfront)
    real(rp),intent(inout)        :: rheap(nfront)
    integer(ip)                   :: ipson,ipfath,ieson,iefath
    real(rp)                      :: rtol
    !
    !     This subroutine adds an element at the end of the heap
    !     and updates the heap
    !     
    rtol=1.05d+00
    ! 
    !     For grids with uniform size, we want the last element to remain the 
    !     last element --> the test is done with a tolerance
    !
    !
    nheap=nheap+1
    lheap(nheap)=ie
    lfront(3,ie)=nheap

    ipson=nheap
    do while(ipson/=1)
       ipfath=ipson/2
       ieson=lheap(ipson)
       iefath=lheap(ipfath)
       if(rheap(ieson)*rtol<rheap(iefath))then
          lheap(ipson)=iefath
          lheap(ipfath)=ieson
          lfront(3,iefath)=ipson
          lfront(3,ieson)=ipfath
          ipson=ipfath
       else
          exit
       endif

    enddo

    !
    !     DBG
    !

!!$  do iheap=1,nheap
!!$     ipfath=iheap
!!$     ipson1=2*ipfath
!!$     ipson2=2*ipfath+1
!!$     iefath=lheap(ipfath)
!!$
!!$     if(ipson1<nheap)then
!!$        ieson1=lheap(ipson1)
!!$        if(rheap(iefath)>rheap(ieson1))then
!!$           write(*,*)'error 1 in addheap'
!!$           stop
!!$        endif
!!$     endif
!!$
!!$     if(ipson2<nheap)then
!!$        ieson2=lheap(ipson2)
!!$        if(rheap(iefath)>rheap(ieson2))then
!!$           write(*,*)'error 2 in addheap'
!!$           stop
!!$        endif
!!$     endif
!!$  enddo

  end subroutine addheap

  subroutine delheap(idel,lheap,nheap,lfront,nfront,rheap)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: idel,nfront
    integer(ip), intent(inout) :: nheap
    integer(ip), intent(inout) :: lfront(3,nfront),lheap(nheap)
    integer(ip)                :: ipfath,ipson1,ipson2,iefath,ieson1,ieson2
    real(rp),intent(inout)     :: rheap(nfront)
    real(rp)                   :: rtol
    !
    !     This subroutine deletes an element from the heap at place idel 
    !     and updates the heap   
    !
    rtol=1.05d+00
    ! 
    !     For grids with uniform size, we want the last element to remain the 
    !     last element --> the test is done with a tolerance
    !
    !
    !
    !     Fill the holes
    !
    ipson1=idel
    do while(ipson1/=1)
       ipfath=ipson1/2
       iefath=lheap(ipfath)
       lheap(ipson1)=iefath
       lfront(3,iefath)=ipson1
       ipson1=ipfath 
    enddo
    !
    !     Restore the tree
    !
    lheap(1)=lheap(nheap)
    lfront(3,lheap(1))=1_ip 
    nheap=nheap-1

    ipfath=1
    ipson1=2*ipfath
    ipson2=ipson1+1
    do while(nheap>=ipson1)
       if(nheap>=ipson2)then
          ieson1=lheap(ipson1)
          ieson2=lheap(ipson2)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(ieson2))then
             if(rheap(ieson1)<rheap(iefath)*rtol)then
                lheap(ipfath)=ieson1
                lheap(ipson1)=iefath
                lfront(3,iefath)=ipson1
                lfront(3,ieson1)=ipfath
                ipfath=ipson1
             else
                exit
             endif
          else
             if(rheap(ieson2)<rheap(iefath)*rtol)then
                lheap(ipfath)=ieson2
                lheap(ipson2)=iefath
                lfront(3,iefath)=ipson2
                lfront(3,ieson2)=ipfath
                ipfath=ipson2

             else
                exit
             endif
          endif

       else
          ieson1=lheap(ipson1)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(iefath)*rtol)then
             lheap(ipfath)=ieson1
             lheap(ipson1)=iefath
             lfront(3,iefath)=ipson1
             lfront(3,ieson1)=ipfath
             ipfath=ipson1
          else
             exit
          endif
       endif
       ipson1=2*ipfath
       ipson2=ipson1+1
    enddo

    !
    !     DBG
    !

!!$  do iheap=1,nheap
!!$     ipfath=iheap
!!$     ipson1=2*ipfath
!!$     ipson2=2*ipfath+1
!!$     iefath=lheap(ipfath)
!!$
!!$     if(ipson1<nheap)then
!!$        ieson1=lheap(ipson1)
!!$        if(rheap(iefath)>rheap(ieson1))then
!!$           write(*,*)'error 1 in delheap'
!!$           stop
!!$        endif
!!$     endif
!!$
!!$     if(ipson2<nheap)then
!!$        ieson2=lheap(ipson2)
!!$        if(rheap(iefath)>rheap(ieson2))then
!!$           write(*,*)'error 2 in delheap'
!!$           stop
!!$        endif
!!$     endif
!!$  enddo

  end subroutine delheap

  subroutine initheap(lheap,nheap,lfront,nfront)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)        :: nheap,nfront,lheap(nheap)
    integer(ip),intent(inout)     :: lfront(3,nfront)
    integer(ip)                   :: iheap,ifront

    do iheap=1,nheap
       ifront=lheap(iheap)
       lfront(3,ifront)=iheap
    enddo

  end subroutine initheap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !     This is the interface for the heap sort for the 
  !     3D advancing front like mesh generator
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gtsmall3d(lheap,nheap,rheap,ismall,lfront,nfront)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(inout)     :: nheap
    integer(ip),intent(in)        :: nfront
    integer(ip),intent(inout)     :: lheap(nheap),ismall,lfront(4,nfront)
    real(rp),intent(in)           :: rheap(nfront)
    integer(ip)                   :: ipson1,ipson2,ipfath,ieson1,ieson2,iefath
    real(rp)                      :: rtol
    ! 
    rtol=1.05d+00
    !
    !     Remove the smallest element
    !
    ismall=lheap(1)
    !
    !     Restore the tree structure
    !
    !     The last element will become the first one. For grids with the same
    !     size, we want the last element to stay the last one for grids to 
    !     be regular --> the test is done with a tolerance 
    !
    lheap(1)=lheap(nheap)
    lfront(4,lheap(1))=1_ip 
    nheap=nheap-1

    ipfath=1
    ipson1=2*ipfath
    ipson2=ipson1+1
    do while(nheap>=ipson1)
       if(nheap>=ipson2)then
          ieson1=lheap(ipson1)
          ieson2=lheap(ipson2)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(ieson2))then
             if(rheap(ieson1)<rheap(iefath)*rtol)then
                lheap(ipfath)=ieson1
                lheap(ipson1)=iefath
                lfront(4,iefath)=ipson1
                lfront(4,ieson1)=ipfath
                ipfath=ipson1
             else
                exit
             endif
          else
             if(rheap(ieson2)<rheap(iefath)*rtol)then
                lheap(ipfath)=ieson2
                lheap(ipson2)=iefath
                lfront(4,iefath)=ipson2
                lfront(4,ieson2)=ipfath
                ipfath=ipson2

             else
                exit
             endif
          endif

       else
          ieson1=lheap(ipson1)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(iefath)*rtol)then
             lheap(ipfath)=ieson1
             lheap(ipson1)=iefath
             lfront(4,iefath)=ipson1
             lfront(4,ieson1)=ipfath
             ipfath=ipson1
          else
             exit
          endif
       endif
       ipson1=2*ipfath
       ipson2=ipson1+1
    enddo

  end subroutine gtsmall3d

  subroutine addheap3d(lheap,nheap,rheap,ie,lfront,nfront)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)        :: ie,nfront
    integer(ip),intent(inout)     :: nheap,lheap(nheap+1),lfront(4,nfront)
    real(rp),intent(inout)        :: rheap(nfront)
    integer(ip)                   :: ipson,ipfath,ieson,iefath
    real(rp)                      :: rtol
    !
    !     This subroutine adds an element at the end of the heap
    !     and updates the heap
    !     
    rtol=1.05d+00
    ! 
    !     For grids with uniform size, we want the last element to remain the 
    !     last element --> the test is done with a tolerance
    !
    !
    nheap=nheap+1
    lheap(nheap)=ie
    lfront(4,ie)=nheap

    ipson=nheap
    do while(ipson/=1)
       ipfath=ipson/2
       ieson=lheap(ipson)
       iefath=lheap(ipfath)
       if(rheap(ieson)*rtol<rheap(iefath))then
          lheap(ipson)=iefath
          lheap(ipfath)=ieson
          lfront(4,iefath)=ipson
          lfront(4,ieson)=ipfath
          ipson=ipfath
       else
          exit
       endif

    enddo

    !
    !     DBG
    !

!!$  do iheap=1,nheap
!!$     ipfath=iheap
!!$     ipson1=2*ipfath
!!$     ipson2=2*ipfath+1
!!$     iefath=lheap(ipfath)
!!$
!!$     if(ipson1<nheap)then
!!$        ieson1=lheap(ipson1)
!!$        if(rheap(iefath)>rheap(ieson1))then
!!$           write(*,*)'error 1 in addheap'
!!$           stop
!!$        endif
!!$     endif
!!$
!!$     if(ipson2<nheap)then
!!$        ieson2=lheap(ipson2)
!!$        if(rheap(iefath)>rheap(ieson2))then
!!$           write(*,*)'error 2 in addheap'
!!$           stop
!!$        endif
!!$     endif
!!$  enddo

  end subroutine addheap3d

  subroutine delheap3d(idel,lheap,nheap,lfront,nfront,rheap)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: idel,nfront
    integer(ip), intent(inout) :: nheap
    integer(ip), intent(inout) :: lfront(4,nfront),lheap(nheap)
    integer(ip)                :: ipfath,ipson1,ipson2,iefath,ieson1,ieson2
    real(rp),intent(inout)     :: rheap(nfront)
    real(rp)                   :: rtol
    !
    !     This subroutine deletes an element from the heap at place idel 
    !     and updates the heap   
    !
    rtol=1.05d+00
    ! 
    !     For grids with uniform size, we want the last element to remain the 
    !     last element --> the test is done with a tolerance
    !
    !
    !
    !     Fill the holes
    !
    ipson1=idel
    do while(ipson1/=1)
       ipfath=ipson1/2
       iefath=lheap(ipfath)
       lheap(ipson1)=iefath
       lfront(4,iefath)=ipson1
       ipson1=ipfath 
    enddo
    !
    !     Restore the tree
    !
    lheap(1)=lheap(nheap)
    lfront(4,lheap(1))=1_ip 
    nheap=nheap-1

    ipfath=1
    ipson1=2*ipfath
    ipson2=ipson1+1
    do while(nheap>=ipson1)
       if(nheap>=ipson2)then
          ieson1=lheap(ipson1)
          ieson2=lheap(ipson2)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(ieson2))then
             if(rheap(ieson1)<rheap(iefath)*rtol)then
                lheap(ipfath)=ieson1
                lheap(ipson1)=iefath
                lfront(4,iefath)=ipson1
                lfront(4,ieson1)=ipfath
                ipfath=ipson1
             else
                exit
             endif
          else
             if(rheap(ieson2)<rheap(iefath)*rtol)then
                lheap(ipfath)=ieson2
                lheap(ipson2)=iefath
                lfront(4,iefath)=ipson2
                lfront(4,ieson2)=ipfath
                ipfath=ipson2

             else
                exit
             endif
          endif

       else
          ieson1=lheap(ipson1)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(iefath)*rtol)then
             lheap(ipfath)=ieson1
             lheap(ipson1)=iefath
             lfront(4,iefath)=ipson1
             lfront(4,ieson1)=ipfath
             ipfath=ipson1
          else
             exit
          endif
       endif
       ipson1=2*ipfath
       ipson2=ipson1+1
    enddo

    !
    !     DBG
    !

!!$  do iheap=1,nheap
!!$     ipfath=iheap
!!$     ipson1=2*ipfath
!!$     ipson2=2*ipfath+1
!!$     iefath=lheap(ipfath)
!!$
!!$     if(ipson1<nheap)then
!!$        ieson1=lheap(ipson1)
!!$        if(rheap(iefath)>rheap(ieson1))then
!!$           write(*,*)'error 1 in delheap'
!!$           stop
!!$        endif
!!$     endif
!!$
!!$     if(ipson2<nheap)then
!!$        ieson2=lheap(ipson2)
!!$        if(rheap(iefath)>rheap(ieson2))then
!!$           write(*,*)'error 2 in delheap'
!!$           stop
!!$        endif
!!$     endif
!!$  enddo

  end subroutine delheap3d

  subroutine initheap3d(lheap,nheap,lfront,nfront)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)        :: nheap,nfront,lheap(nheap)
    integer(ip),intent(inout)     :: lfront(4,nfront)
    integer(ip)                   :: iheap,ifront

    do iheap=1,nheap
       ifront=lheap(iheap)
       lfront(4,ifront)=iheap
    enddo

  end subroutine initheap3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !     This is the standard interface for the heap sort
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine buildheap(lheap,nheap,rheap)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)        :: nheap
    integer(ip),intent(inout)     :: lheap(nheap)
    real(rp),intent(in)           :: rheap(nheap)
    integer(ip)                   :: ipson,ipfath,ieson,iefath,ie
    !
    !     This subroutine builds the heap 
    !
    do ie=1,nheap

       ipson=ie
       lheap(ie)=ie
       do while(ipson/=1)
          ipfath=ipson/2
          ieson=lheap(ipson)
          iefath=lheap(ipfath)
          if(rheap(ieson)<rheap(iefath))then
             lheap(ipson)=iefath
             lheap(ipfath)=ieson
             ipson=ipfath
          else
             exit
          endif

       enddo

    enddo
    !
    !     DBG
    !

!!$  do iheap=1,nheap
!!$     ipfath=iheap
!!$     ipson1=2*ipfath
!!$     ipson2=2*ipfath+1
!!$     iefath=lheap(ipfath)
!!$
!!$     if(ipson1<nheap)then
!!$        ieson1=lheap(ipson1)
!!$        if(rheap(iefath)>rheap(ieson1))then
!!$           write(*,*)'error 1 in buildheap'
!!$           stop
!!$        endif
!!$     endif
!!$
!!$     if(ipson2<nheap)then
!!$        ieson2=lheap(ipson2)
!!$        if(rheap(iefath)>rheap(ieson2))then
!!$           write(*,*)'error 2 in buildheap'
!!$           stop
!!$        endif
!!$     endif
!!$  enddo

  end subroutine buildheap

  subroutine gtkey(lheap,nheap,rheap,ismall,nsize)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(inout)     :: nheap
    integer(ip),intent(in)        :: nsize
    integer(ip),intent(inout)     :: lheap(nheap),ismall
    real(rp),intent(in)           :: rheap(nsize)
    integer(ip)                   :: ipson1,ipson2,ipfath,ieson1,ieson2,iefath
    !
    !     Remove the smallest element
    !
    ismall=lheap(1)
    !
    !     Restore the tree structure
    !
    lheap(1)=lheap(nheap)
    nheap=nheap-1

    ipfath=1
    ipson1=2*ipfath
    ipson2=ipson1+1
    do while(nheap>=ipson1)
       if(nheap>=ipson2)then
          ieson1=lheap(ipson1)
          ieson2=lheap(ipson2)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(ieson2))then
             if(rheap(ieson1)<rheap(iefath))then
                lheap(ipfath)=ieson1
                lheap(ipson1)=iefath
                ipfath=ipson1
             else
                exit
             endif
          else
             if(rheap(ieson2)<rheap(iefath))then
                lheap(ipfath)=ieson2
                lheap(ipson2)=iefath
                ipfath=ipson2

             else
                exit
             endif
          endif

       else
          ieson1=lheap(ipson1)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(iefath))then
             lheap(ipfath)=ieson1
             lheap(ipson1)=iefath
             ipfath=ipson1
          else
             exit
          endif
       endif
       ipson1=2*ipfath
       ipson2=ipson1+1
    enddo

  end subroutine gtkey

  subroutine addkey(lheap,nheap,rheap,ikey,nsize)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        : memor_msh
    implicit none
    integer(ip),intent(inout)     :: nheap
    integer(ip),intent(in)        :: ikey,nsize
    integer(ip),pointer           :: lheap(:)
    real(rp),intent(in)           :: rheap(nsize)
    integer(ip)                   :: ipson,ipfath,ieson,iefath
    !
    !     Resize lheap and rheap 
    !
    nheap=nheap+1_ip
    call memrea(nheap,memor_msh,'LHEAP','addkey',lheap)
    lheap(nheap)=ikey
    !
    !     Initialize ipson
    !
    ipson=nheap

    do while(ipson/=1)
       ipfath=ipson/2
       ieson=lheap(ipson)
       iefath=lheap(ipfath)
       if(rheap(ieson)<rheap(iefath))then
          lheap(ipson)=iefath
          lheap(ipfath)=ieson
          ipson=ipfath
       else
          exit
       endif

    enddo

  end subroutine addkey

  subroutine buildheapi(lheap,nheap,rheap)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)        :: nheap
    integer(ip),intent(inout)     :: lheap(nheap)
    integer(ip),intent(in)        :: rheap(nheap)
    integer(ip)                   :: ipson,ipfath,ieson,iefath,ie
    !
    !     This subroutine builds the heap for the advancing front like grid
    !
    do ie=1,nheap

       ipson=ie
       lheap(ie)=ie
       do while(ipson/=1)
          ipfath=ipson/2
          ieson=lheap(ipson)
          iefath=lheap(ipfath)
          if(rheap(ieson)<rheap(iefath))then
             lheap(ipson)=iefath
             lheap(ipfath)=ieson
             ipson=ipfath
          else
             exit
          endif

       enddo

    enddo
   
  end subroutine buildheapi

  subroutine gtkeyi(lheap,nheap,rheap,ismall,nsize)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(inout)     :: nheap
    integer(ip),intent(in)        :: nsize
    integer(ip),intent(inout)     :: lheap(nheap),ismall
    integer(ip),intent(in)        :: rheap(nsize)
    integer(ip)                   :: ipson1,ipson2,ipfath,ieson1,ieson2,iefath
    !
    !     Remove the smallest element
    !
    ismall=lheap(1)
    !
    !     Restore the tree structure
    !
    lheap(1)=lheap(nheap)
    nheap=nheap-1

    ipfath=1
    ipson1=2*ipfath
    ipson2=ipson1+1
    do while(nheap>=ipson1)
       if(nheap>=ipson2)then
          ieson1=lheap(ipson1)
          ieson2=lheap(ipson2)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<rheap(ieson2))then
             if(rheap(ieson1)<=rheap(iefath))then
                lheap(ipfath)=ieson1
                lheap(ipson1)=iefath
                ipfath=ipson1
             else
                exit
             endif
          else
             if(rheap(ieson2)<=rheap(iefath))then
                lheap(ipfath)=ieson2
                lheap(ipson2)=iefath
                ipfath=ipson2

             else
                exit
             endif
          endif

       else
          ieson1=lheap(ipson1)
          iefath=lheap(ipfath)
          if(rheap(ieson1)<=rheap(iefath))then
             lheap(ipfath)=ieson1
             lheap(ipson1)=iefath
             ipfath=ipson1
          else
             exit
          endif
       endif
       ipson1=2*ipfath
       ipson2=ipson1+1
    enddo

  end subroutine gtkeyi

  subroutine addkeyi(lheap,nheap,rheap,ikey,nsize)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        : memor_msh
    implicit none
    integer(ip),intent(inout)     :: nheap
    integer(ip),intent(in)        :: ikey,nsize
    integer(ip),pointer           :: lheap(:)
    integer(ip),intent(in)        :: rheap(nsize)
    integer(ip)                   :: ipson,ipfath,ieson,iefath
    !
    !     Resize lheap and rheap 
    !
    nheap=nheap+1_ip
    call memrea(nheap,memor_msh,'LHEAP','addkey',lheap)
    lheap(nheap)=ikey
    !
    !     Initialize ipson
    !
    ipson=nheap

    do while(ipson/=1)
       ipfath=ipson/2
       ieson=lheap(ipson)
       iefath=lheap(ipfath)
       if(rheap(ieson)<rheap(iefath))then
          lheap(ipson)=iefath
          lheap(ipfath)=ieson
          ipson=ipfath
       else
          exit
       endif

    enddo

  end subroutine addkeyi

end module mod_sort
