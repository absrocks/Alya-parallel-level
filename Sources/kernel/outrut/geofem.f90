subroutine geofem

!------------------------------------------------------------------------
!
! This routine dumps the geometry in Femview format into the 
! postprocess file
!      
!------------------------------------------------------------------------
  use      def_kintyp
  use      def_parame
  use      def_domain
  use      def_master
  use      def_inpout
  use      mod_iofile
  implicit none
  integer(ip) :: lnodw(12,12)
  integer(ip) :: istop,nnodw,iflag,idime,ielty,inode
  integer(ip) :: ipoin,nelew,kelem,jelem,ielem
!
! Coordinates
!
  write(lun_postp,10) 1,'C',title(1:6)
  write(lun_postp,10) 2,'C','Coord.'
  do ipoin  = 1,npoin
     write(lun_postp,20) -1,ipoin,(coord(idime,ipoin),idime=1,ndime)
  enddo
  write(lun_postp,30) -3
  write(lun_postp,10) 3,'C','Conne.'
  
  kelem=0
  do ielty=1,nelty
     do ielem=1,nelem
        nnodw = 0
        istop = 1
        if(ltype(ielem)==ielty) then
           if(ndime==2) then
              if(nnode(ielty)==4.and.ltopo(ielty)==1) then             ! P1+bubble
                 call top004(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 7
                 istop = 0
              else if(nnode(ielty)==5) then                            ! Q1+bubble
                 call top05a(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 7
                 istop = 0
              else if(nnode(ielty)==7) then                            ! P2+bubble = P2^+
                 call top007(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 9
                 istop = 0
              else if(nnode(ielty)==10) then                           ! P3
                 call top10a(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 7
                 istop = 0 
              else if(nnode(ielty)==9) then                            ! Q2
                 call top09a(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 9
                 istop = 0
              else if(nnode(ielty)==16) then                           ! Q3
                 call top016(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 9
                 istop = 0
              endif
           else if(ndime==3) then
              if(nnode(ielty)==5) then                                 ! P1+bubble
                 call top05b(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 3
                 istop = 0
              else if(nnode(ielty)==9) then                            ! Q1+bubble
                 call top09b(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 3
                 istop = 0
              else if(nnode(ielty)==11) then                           ! P2+bubble
                 call top011(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 3
                 istop = 0
              else if(nnode(ielty)==15) then                           ! P2+bubbles = P2^+
                 call top015(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 1
                 istop = 0
              else if(nnode(ielty)==27) then                           ! Q2
                 call top027(lnods(1,ielem),lnodw,nelew,nnodw)
                 iflag = 1
                 istop = 0
              end if
           end if          
           if(nnodw==0) then
              if(ndime==2) then
                 if(nnode(ielty)==3) then
                    iflag = 7
                    istop = 0
                 else if(nnode(ielty)==4) then
                    iflag = 9
                    istop = 0
                 else if(nnode(ielty)==6) then
                    iflag = 8
                    istop = 0 
                 else if(nnode(ielty)==8) then
                    iflag = 10
                    istop = 0
                 end if
              else if(ndime==3) then
                 if(nnode(ielty)==4) then
                    iflag = 3
                    istop = 0
                 else if(nnode(ielty)==8) then
                    iflag = 1
                    istop = 0 
                 else if(nnode(ielty)==10) then
                    iflag = 6           
                    istop = 0
                 else if(nnode(ielty)==20) then
                    iflag = 17
                    istop = 0
                 end if
              end if
           end if           
           if(istop==1) call runend('GEOFEM: NOT AVAILABLE ELEMENT')     
           if(nnodw==0) then
              kelem=kelem+1
              write(lun_postp,40) -1,kelem,iflag,1,1
              write(lun_postp,40) -2,(lnods(inode,ielem),inode=1,nnode(ielty))
           else 
              do jelem=1,nelew
                 kelem=kelem+1
                 write(lun_postp,40) -1,kelem,iflag,1,1
                 write(lun_postp,40) -2,(lnodw(inode,jelem),inode=1,nnodw)
              end do
           end if
        end if
     end do
  end do

  write(lun_postp,30) -3
  flush(lun_postp)

10 format(1x,i4,a1,a6)
20 format(1x,i2,i5,3e12.5)
30 format(1x,i2)
40 format(1x,i2,20i5)

end subroutine geofem

subroutine top004(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 3
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(2)
  lnodw(3,melem) = lnods(4)
  melem = melem+1
  lnodw(1,melem) = lnods(2)
  lnodw(2,melem) = lnods(3)
  lnodw(3,melem) = lnods(4)
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(4)
  lnodw(3,melem) = lnods(3)

end subroutine top004

subroutine top007(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 4
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(4)
  lnodw(3,melem) = lnods(7)
  lnodw(4,melem) = lnods(6)
  melem = melem+1
  lnodw(1,melem) = lnods(4)
  lnodw(2,melem) = lnods(2)
  lnodw(3,melem) = lnods(5)
  lnodw(4,melem) = lnods(7)
  melem = melem+1
  lnodw(1,melem) = lnods(6)
  lnodw(2,melem) = lnods(7)
  lnodw(3,melem) = lnods(5)
  lnodw(4,melem) = lnods(3)

end subroutine top007

subroutine top10a(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 3
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1 )
  lnodw(2,melem) = lnods(4 )
  lnodw(3,melem) = lnods(9 )
  melem = melem+1
  lnodw(1,melem) = lnods(4 )
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(9 )
  melem = melem+1
  lnodw(1,melem) = lnods(4 )
  lnodw(2,melem) = lnods(5 )
  lnodw(3,melem) = lnods(10)
  melem = melem+1
  lnodw(1,melem) = lnods(5 )
  lnodw(2,melem) = lnods(6 )
  lnodw(3,melem) = lnods(10)
  melem = melem+1
  lnodw(1,melem) = lnods(5 )
  lnodw(2,melem) = lnods(2 )
  lnodw(3,melem) = lnods(6 )
  melem = melem+1
  lnodw(1,melem) = lnods(9 )
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(8 )
  melem = melem+1
  lnodw(1,melem) = lnods(10)
  lnodw(2,melem) = lnods(7 )
  lnodw(3,melem) = lnods(8 )
  melem = melem+1
  lnodw(1,melem) = lnods(10)
  lnodw(2,melem) = lnods(6 )
  lnodw(3,melem) = lnods(7 )
  melem = melem+1
  lnodw(1,melem) = lnods(8 )
  lnodw(2,melem) = lnods(7 )
  lnodw(3,melem) = lnods(3 )

end subroutine top10a

subroutine top011(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 4
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(7)
  lnodw(4,melem) = lnods(8)
  melem = melem+1
  lnodw(1,melem) = lnods(2)
  lnodw(2,melem) = lnods(6)
  lnodw(3,melem) = lnods(5)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(3)
  lnodw(2,melem) = lnods(7)
  lnodw(3,melem) = lnods(6)
  lnodw(4,melem) = lnods(10)
  melem = melem+1
  lnodw(1,melem) = lnods(4)
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(9)
  lnodw(4,melem) = lnods(8)
  melem = melem+1
  lnodw(1,melem) = lnods(8)
  lnodw(2,melem) = lnods(9)
  lnodw(3,melem) = lnods(5)
  lnodw(4,melem) = lnods(11)
  melem = melem+1
  lnodw(1,melem) = lnods(7)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(6)
  lnodw(4,melem) = lnods(11)
  melem = melem+1
  lnodw(1,melem) = lnods(7)
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(8)
  lnodw(4,melem) = lnods(11)
  melem = melem+1
  lnodw(1,melem) = lnods(6)
  lnodw(2,melem) = lnods(9)
  lnodw(3,melem) = lnods(10)
  lnodw(4,melem) = lnods(11)
  melem = melem+1
  lnodw(1,melem) = lnods(7)
  lnodw(2,melem) = lnods(8)
  lnodw(3,melem) = lnods(5)
  lnodw(4,melem) = lnods(11)
  melem = melem+1
  lnodw(1,melem) = lnods(5)
  lnodw(2,melem) = lnods(9)
  lnodw(3,melem) = lnods(6)
  lnodw(4,melem) = lnods(11)
  melem = melem+1
  lnodw(1,melem) = lnods(6)
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(7)
  lnodw(4,melem) = lnods(11)
  melem = melem+1
  lnodw(1,melem) = lnods(8)
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(9)
  lnodw(4,melem) = lnods(11)

end subroutine top011

subroutine top015(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 8
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(11)
  lnodw(4,melem) = lnods(7)
  lnodw(5,melem) = lnods(8)
  lnodw(6,melem) = lnods(12)
  lnodw(7,melem) = lnods(15)
  lnodw(8,melem) = lnods(14)
  melem = melem+1
  lnodw(1,melem) = lnods(5)
  lnodw(2,melem) = lnods(2)
  lnodw(3,melem) = lnods(6)
  lnodw(4,melem) = lnods(11)
  lnodw(5,melem) = lnods(12)
  lnodw(6,melem) = lnods(9)
  lnodw(7,melem) = lnods(13)
  lnodw(8,melem) = lnods(15)
  melem = melem+1
  lnodw(1,melem) = lnods(7)
  lnodw(2,melem) = lnods(11)
  lnodw(3,melem) = lnods(6)
  lnodw(4,melem) = lnods(3)
  lnodw(5,melem) = lnods(14)
  lnodw(6,melem) = lnods(15)
  lnodw(7,melem) = lnods(13)
  lnodw(8,melem) = lnods(10)
  melem = melem+1
  lnodw(1,melem) = lnods(8)
  lnodw(2,melem) = lnods(12)
  lnodw(3,melem) = lnods(15)
  lnodw(4,melem) = lnods(14)
  lnodw(5,melem) = lnods(4)
  lnodw(6,melem) = lnods(9)
  lnodw(7,melem) = lnods(13)
  lnodw(8,melem) = lnods(10)

end subroutine top015

subroutine top016(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 4
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1 )
  lnodw(2,melem) = lnods(5 )
  lnodw(3,melem) = lnods(13)
  lnodw(4,melem) = lnods(12)
  melem = melem+1
  lnodw(1,melem) = lnods(5 )
  lnodw(2,melem) = lnods(6 )
  lnodw(3,melem) = lnods(14)
  lnodw(4,melem) = lnods(13)
  melem = melem+1
  lnodw(1,melem) = lnods(6 )
  lnodw(2,melem) = lnods(2 )
  lnodw(3,melem) = lnods(7 )
  lnodw(4,melem) = lnods(14)
  melem = melem+1
  lnodw(1,melem) = lnods(12)
  lnodw(2,melem) = lnods(13)
  lnodw(3,melem) = lnods(16)
  lnodw(4,melem) = lnods(11)
  melem = melem+1
  lnodw(1,melem) = lnods(13)
  lnodw(2,melem) = lnods(14)
  lnodw(3,melem) = lnods(15)
  lnodw(4,melem) = lnods(16)
  melem = melem+1
  lnodw(1,melem) = lnods(14)
  lnodw(2,melem) = lnods(7 )
  lnodw(3,melem) = lnods(8 )
  lnodw(4,melem) = lnods(15)
  melem = melem+1
  lnodw(1,melem) = lnods(11)
  lnodw(2,melem) = lnods(16)
  lnodw(3,melem) = lnods(10)
  lnodw(4,melem) = lnods(4 )
  melem = melem+1
  lnodw(1,melem) = lnods(16)
  lnodw(2,melem) = lnods(15)
  lnodw(3,melem) = lnods(9 )
  lnodw(4,melem) = lnods(10)
  melem = melem+1
  lnodw(1,melem) = lnods(15)
  lnodw(2,melem) = lnods(8 )
  lnodw(3,melem) = lnods(3 )
  lnodw(4,melem) = lnods(9 )

end subroutine top016

subroutine top027(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 8
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(9)
  lnodw(3,melem) = lnods(21)
  lnodw(4,melem) = lnods(12)
  lnodw(5,melem) = lnods(13)
  lnodw(6,melem) = lnods(22)
  lnodw(7,melem) = lnods(27)
  lnodw(8,melem) = lnods(25)
  melem = melem+1
  lnodw(1,melem) = lnods(9)
  lnodw(2,melem) = lnods(2)
  lnodw(3,melem) = lnods(10)
  lnodw(4,melem) = lnods(21)
  lnodw(5,melem) = lnods(22)
  lnodw(6,melem) = lnods(14)
  lnodw(7,melem) = lnods(23)
  lnodw(8,melem) = lnods(27)
  melem = melem+1
  lnodw(1,melem) = lnods(21)
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(3)
  lnodw(4,melem) = lnods(11)
  lnodw(5,melem) = lnods(27)
  lnodw(6,melem) = lnods(23)
  lnodw(7,melem) = lnods(15)
  lnodw(8,melem) = lnods(24)
  melem = melem+1
  lnodw(1,melem) = lnods(12)
  lnodw(2,melem) = lnods(21)
  lnodw(3,melem) = lnods(11)
  lnodw(4,melem) = lnods(4)
  lnodw(5,melem) = lnods(25)
  lnodw(6,melem) = lnods(27)
  lnodw(7,melem) = lnods(24)
  lnodw(8,melem) = lnods(16)
  melem = melem+1
  lnodw(1,melem) = lnods(13)
  lnodw(2,melem) = lnods(22)
  lnodw(3,melem) = lnods(27)
  lnodw(4,melem) = lnods(25)
  lnodw(5,melem) = lnods(5)
  lnodw(6,melem) = lnods(17)
  lnodw(7,melem) = lnods(26)
  lnodw(8,melem) = lnods(20)
  melem = melem+1
  lnodw(1,melem) = lnods(22)
  lnodw(2,melem) = lnods(14)
  lnodw(3,melem) = lnods(23)
  lnodw(4,melem) = lnods(27)
  lnodw(5,melem) = lnods(17)
  lnodw(6,melem) = lnods(6)
  lnodw(7,melem) = lnods(18)
  lnodw(8,melem) = lnods(26)
  melem = melem+1
  lnodw(1,melem) = lnods(27)
  lnodw(2,melem) = lnods(23)
  lnodw(3,melem) = lnods(15)
  lnodw(4,melem) = lnods(24)
  lnodw(5,melem) = lnods(26)
  lnodw(6,melem) = lnods(18)
  lnodw(7,melem) = lnods(7)
  lnodw(8,melem) = lnods(19)
  melem = melem+1
  lnodw(1,melem) = lnods(25)
  lnodw(2,melem) = lnods(27)
  lnodw(3,melem) = lnods(24)
  lnodw(4,melem) = lnods(16)
  lnodw(5,melem) = lnods(20)
  lnodw(6,melem) = lnods(26)
  lnodw(7,melem) = lnods(19)
  lnodw(8,melem) = lnods(8)

end subroutine top027

subroutine top05a(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 3
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(2)
  lnodw(3,melem) = lnods(5)
  melem = melem+1
  lnodw(1,melem) = lnods(2)
  lnodw(2,melem) = lnods(3)
  lnodw(3,melem) = lnods(5)
  melem = melem+1
  lnodw(1,melem) = lnods(3)
  lnodw(2,melem) = lnods(4)
  lnodw(3,melem) = lnods(5)
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(4)

end subroutine top05a

subroutine top05b(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 4
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(4)
  lnodw(3,melem) = lnods(2)
  lnodw(4,melem) = lnods(5)
  melem = melem+1
  lnodw(1,melem) = lnods(2)
  lnodw(2,melem) = lnods(4)
  lnodw(3,melem) = lnods(3)
  lnodw(4,melem) = lnods(5)
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(3)
  lnodw(3,melem) = lnods(4)
  lnodw(4,melem) = lnods(5)
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(2)
  lnodw(3,melem) = lnods(3)
  lnodw(4,melem) = lnods(5)

end subroutine top05b

subroutine top09a(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 4
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(9)
  lnodw(4,melem) = lnods(8)
  melem = melem+1
  lnodw(1,melem) = lnods(5)
  lnodw(2,melem) = lnods(2)
  lnodw(3,melem) = lnods(6)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(9)
  lnodw(2,melem) = lnods(6)
  lnodw(3,melem) = lnods(3)
  lnodw(4,melem) = lnods(7)
  melem = melem+1
  lnodw(1,melem) = lnods(8)
  lnodw(2,melem) = lnods(9)
  lnodw(3,melem) = lnods(7)
  lnodw(4,melem) = lnods(4)

end subroutine top09a

subroutine top09b(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 4
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(2)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(5)
  lnodw(2,melem) = lnods(6)
  lnodw(3,melem) = lnods(2)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(2)
  lnodw(2,melem) = lnods(6)
  lnodw(3,melem) = lnods(7)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(7)
  lnodw(2,melem) = lnods(3)
  lnodw(3,melem) = lnods(2)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(3)
  lnodw(2,melem) = lnods(7)
  lnodw(3,melem) = lnods(4)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(7)
  lnodw(2,melem) = lnods(8)
  lnodw(3,melem) = lnods(4)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(4)
  lnodw(2,melem) = lnods(8)
  lnodw(3,melem) = lnods(5)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(4)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(1)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(6)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(7)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(7)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(8)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(2)
  lnodw(3,melem) = lnods(4)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(2)
  lnodw(2,melem) = lnods(3)
  lnodw(3,melem) = lnods(4)
  lnodw(4,melem) = lnods(9)

end subroutine top09b

subroutine top10b(lnods,lnodw,melem,nnodw)
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(out) :: lnodw(12,12),nnodw,melem

  nnodw = 4
  melem = 0
  melem = melem+1
  lnodw(1,melem) = lnods(1)
  lnodw(2,melem) = lnods(5)
  lnodw(3,melem) = lnods(7)
  lnodw(4,melem) = lnods(8)
  melem = melem+1
  lnodw(1,melem) = lnods(2)
  lnodw(2,melem) = lnods(6)
  lnodw(3,melem) = lnods(5)
  lnodw(4,melem) = lnods(9)
  melem = melem+1
  lnodw(1,melem) = lnods(3)
  lnodw(2,melem) = lnods(7)
  lnodw(3,melem) = lnods(6)
  lnodw(4,melem) = lnods(10)
  melem = melem+1
  lnodw(1,melem) = lnods(4)
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(9)
  lnodw(4,melem) = lnods(8)
  melem = melem+1
  lnodw(1,melem) = lnods(8)
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(9)
  lnodw(4,melem) = lnods(7)
  melem = melem+1
  lnodw(1,melem) = lnods(8)
  lnodw(2,melem) = lnods(9)
  lnodw(3,melem) = lnods(5)
  lnodw(4,melem) = lnods(7)
  melem = melem+1
  lnodw(1,melem) = lnods(9)
  lnodw(2,melem) = lnods(10)
  lnodw(3,melem) = lnods(6)
  lnodw(4,melem) = lnods(7)
  melem = melem+1
  lnodw(1,melem) = lnods(9)
  lnodw(2,melem) = lnods(6)
  lnodw(3,melem) = lnods(5)
  lnodw(4,melem) = lnods(7)

end subroutine top10b
