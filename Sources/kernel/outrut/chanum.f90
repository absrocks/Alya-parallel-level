subroutine chanum(ptopo,pnode,nocha,noinv,westr)

!-----------------------------------------------------------------------
!****f* Mathru/chanum
! NAME 
!    chanum
! DESCRIPTION
!    This routine changes nodal numbering for stream-function algorithm
!    and sets weigths for central node
! USES
! USED BY
!    stream 
!***
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: ptopo,pnode
  integer(ip), intent(out) :: nocha(pnode), noinv(pnode)
  real(rp),    intent(out) :: westr(pnode-1)
  integer(ip)              :: istop,inode
!
! Initialization
!
  do inode=1,pnode
     nocha(inode)=inode
     noinv(inode)=inode
  end do
  istop=0
!
! Quads & bricks
!
  if(ptopo==0) then
     if(pnode==4) then
        continue
     else if(pnode==5) then
        nocha(pnode)=0
        noinv(pnode)=0
        do inode=1,pnode-1
           westr(inode)=0.25_rp
        end do
     else if((pnode==8).or.(pnode==9)) then
        nocha(2)=3
        nocha(3)=5
        nocha(4)=7
        nocha(5)=2
        nocha(6)=4
        nocha(7)=6
        noinv(2)=5
        noinv(3)=2
        noinv(4)=6
        noinv(5)=3
        noinv(6)=7
        noinv(7)=4
        if(pnode==9) then
           nocha(pnode)=0
           noinv(pnode)=0
           do inode=1,4
              westr(inode  )=-0.25_rp
              westr(inode+4)= 0.50_rp
           end do
        end if
     else if(pnode==16) then
        nocha(2 )=4
        nocha(3 )=7
        nocha(4 )=10
        nocha(5 )=2
        nocha(6 )=3
        nocha(7 )=5
        nocha(8 )=6
        nocha(9 )=8
        nocha(10)=9
        nocha(11)=11
        nocha(12)=12
        nocha(13)=0
        nocha(14)=0
        nocha(15)=0
        nocha(16)=0
        noinv(2 )=5
        noinv(3 )=6
        noinv(4 )=2
        noinv(5 )=7
        noinv(6 )=8
        noinv(7 )=3
        noinv(8 )=9
        noinv(9 )=10
        noinv(10)=4
        noinv(11)=11
        noinv(12)=12
        noinv(13)=0
        noinv(14)=0
        noinv(15)=0
        noinv(16)=0
        do inode=1,4
           westr(inode)=0.25_rp
        end do
        do inode=5,8
           westr(inode)=0.0_rp
        end do
     else
        istop=1
     end if
!
! Triangles & tetrahedra
!
  else if(ptopo==1) then
     if(pnode==3) then
        continue
     else if(pnode==4) then
        nocha(pnode)=0
        noinv(pnode)=0
        do inode=1,pnode-1
           westr(inode)=1.0_rp/3.0_rp
        end do
     else if((pnode==6).or.(pnode==7)) then
        nocha(2)=3
        nocha(3)=5
        nocha(4)=2
        nocha(5)=4
        noinv(2)=4
        noinv(3)=2
        noinv(4)=5
        noinv(5)=3
        if(pnode==7) then
           nocha(pnode)=0
           noinv(pnode)=0
           do inode=1,3
              westr(inode  )=-1.0_rp/9.0_rp
              westr(inode+3)= 4.0_rp/9.0_rp
           end do
        end if
     else if(pnode==10) then
        nocha(2 )=4
        nocha(3 )=7
        nocha(4 )=2
        nocha(5 )=3
        nocha(6 )=5
        nocha(7 )=6
        nocha(10)=0
        noinv(2 )=4
        noinv(3 )=5
        noinv(4 )=2
        noinv(5 )=6
        noinv(6 )=7
        noinv(7 )=3
        noinv(10)=0
        do inode=1,3
           westr(inode)=1.0_rp/3.0_rp
        end do
        do inode=4,8
           westr(inode)=0.0_rp
        end do
     else
        istop=1
     end if
  end if
!
! Error
!
  if(istop==1)&
       call runend('CHANUM: STREAMFUNCTION FACILITY NOT AVAILABLE')

end subroutine chanum
