subroutine elsest_octfin2(ndime,ithre,xcoor,current)
  !
  ! Returns pointer CURRENT(ITHRE)%O to box containing XCOOR
  !
  use def_elsest, only : ip,rp,lg,poiarr
  use def_elsest, only : tree_root
  use def_elsest, only : comin,comax
  !use mod_elsest
  implicit none
  integer(ip),  intent(in)             :: ndime,ithre
  real(rp),     intent(in)             :: xcoor(ndime)
  type(poiarr), intent(inout), pointer :: current(:)
  integer(ip)                          :: ichild,generation,ipoin
  logical(lg)                          :: lfoun

  current(ithre) % o => tree_root
  generation = 1

  do while( current(ithre) % o % whoiam == 0 )

     if( ndime == 3 ) then
        !
        ! 3D case
        !
        lfoun = .false.
        childloop8: do ichild = 1,8

           if(    xcoor(1) >= current(ithre) % o % children(ichild) % minc(1) .and. &
                & xcoor(1) <= current(ithre) % o % children(ichild) % maxc(1) .and. &
                & xcoor(2) >= current(ithre) % o % children(ichild) % minc(2) .and. &
                & xcoor(2) <= current(ithre) % o % children(ichild) % maxc(2) .and. &
                & xcoor(3) >= current(ithre) % o % children(ichild) % minc(3) .and. &
                & xcoor(3) <= current(ithre) % o % children(ichild) % maxc(3) ) then
              current(ithre) % o => current(ithre) % o % children(ichild)
              lfoun      = .true.
              generation = generation + 1
              exit childloop8
           end if

        end do childloop8

        if( .not. lfoun ) then
           print*,'Generation  = ',generation
           print*,'Test point  = ',xcoor(1:3)
           write(100,'(a)') 'MESH elsest_octfin dimension 3 Elemtype Hexahedra Nnode 8'
           write(100,'(a)') 'coordinates'
           ipoin = 0
           do ichild = 1,8              
              ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,current(ithre) % o % children(ichild) % minc(1:3)
              ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,current(ithre) % o % children(ichild) % minc(1),current(ithre) % o % children(ichild) % maxc(2),current(ithre) % o % children(ichild) % minc(3) ! 2
              ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,current(ithre) % o % children(ichild) % minc(1),current(ithre) % o % children(ichild) % maxc(2),current(ithre) % o % children(ichild) % maxc(3) ! 3
              ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,current(ithre) % o % children(ichild) % minc(1),current(ithre) % o % children(ichild) % minc(2),current(ithre) % o % children(ichild) % maxc(3) ! 4
              ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,current(ithre) % o % children(ichild) % maxc(1),current(ithre) % o % children(ichild) % minc(2),current(ithre) % o % children(ichild) % minc(3) ! 5
              ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,current(ithre) % o % children(ichild) % maxc(1),current(ithre) % o % children(ichild) % maxc(2),current(ithre) % o % children(ichild) % minc(3) ! 6
              ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,current(ithre) % o % children(ichild) % maxc(1),current(ithre) % o % children(ichild) % maxc(2),current(ithre) % o % children(ichild) % maxc(3) ! 7
              ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,current(ithre) % o % children(ichild) % maxc(1),current(ithre) % o % children(ichild) % minc(2),current(ithre) % o % children(ichild) % maxc(3) ! 8
           end do
           ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,comin(1),comin(2),comin(3)
           ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,comin(1),comax(2),comin(3)
           ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,comin(1),comax(2),comax(3)
           ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,comin(1),comin(2),comax(3)
           ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,comax(1),comin(2),comin(3)
           ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,comax(1),comax(2),comin(3)
           ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,comax(1),comax(2),comax(3)
           ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,comax(1),comin(2),comax(3)
           ipoin = ipoin + 1 ; write(100,'(i3,1x,3(1x,e12.6))') ipoin,xcoor(1:3)
           write(100,'(a)') 'end coordinates'
           ipoin = 1
           write(100,'(a)') 'elements'
           ipoin = ipoin + 1 ; write(100,'(10(1x,i3))') ipoin, 1, 2, 3, 4, 5, 6, 7, 8      
           ipoin = ipoin + 1 ; write(100,'(10(1x,i3))') ipoin, 9,10,11,12,13,14,15,16
           ipoin = ipoin + 1 ; write(100,'(10(1x,i3))') ipoin,17,18,19,20,21,22,23,24
           ipoin = ipoin + 1 ; write(100,'(10(1x,i3))') ipoin,25,26,27,28,29,30,31,32
           ipoin = ipoin + 1 ; write(100,'(10(1x,i3))') ipoin,33,34,35,36,37,38,39,40
           ipoin = ipoin + 1 ; write(100,'(10(1x,i3))') ipoin,41,42,43,44,45,46,47,48
           ipoin = ipoin + 1 ; write(100,'(10(1x,i3))') ipoin,49,50,51,52,53,54,55,56
           ipoin = ipoin + 1 ; write(100,'(10(1x,i3))') ipoin,57,58,59,60,61,62,63,64
           write(100,'(a)') 'end elements'
           write(100,'(a)') 'MESH comin_comax dimension 3 Elemtype Hexahedra Nnode 8'
           write(100,'(a)') 'elements'
           ipoin = ipoin + 1 ; write(100,'(10(1x,i3))') ipoin,65,66,67,68,69,70,71,72
           write(100,'(a)') 'end elements'
           call runend('FATHER WITHOUT SON!')
        end if

     else if( ndime == 2 ) then
        !
        ! 2D case
        !
        childloop4: do ichild=1,4
           if(    xcoor(1) >= current(ithre) % o % children(ichild) % minc(1) .and. &
                & xcoor(1) <= current(ithre) % o % children(ichild) % maxc(1) .and. &
                & xcoor(2) >= current(ithre) % o % children(ichild) % minc(2) .and. &
                & xcoor(2) <= current(ithre) % o % children(ichild) % maxc(2)  ) then
              current(ithre) % o => current(ithre) % o % children(ichild)
              exit childloop4
           end if
        end do childloop4

     end if
  end do

end subroutine elsest_octfin2
