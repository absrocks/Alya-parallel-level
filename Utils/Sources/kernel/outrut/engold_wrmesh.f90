subroutine engold_wrmesh
  !-----------------------------------------------------------------------
  !****f* outrut/engold_wrmesh
  ! NAME 
  !    engold_wrmesh
  ! DESCRIPTION
  !    This routine writes the mesh in Ensight format
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use      def_kintyp
  use      def_parame
  use      def_master
  use      def_domain
  use      def_elmtyp
  use      mod_iofile
  use      def_postpr
  implicit none
  integer(ip)             :: ipart,ipoin,ielem,inode,idime,ielty
  integer(ip)             :: iesta,iesto
  real(rp)                :: xauxi
  character(40)           :: chens
  character(10)           :: cengo


  parts_pos(1)%npoin = npoin

  chens= adjustl(trim(namda))
  write(lun_outpu_dom,15) 'Problem name:  ',adjustl(trim(chens))
  write(lun_outpu_dom,10) 'Geometry file  '
  write(lun_outpu_dom,10) 'node id given'
  write(lun_outpu_dom,10) 'element id given'
  do ipart=1,npart_pos
     write(lun_outpu_dom,10) 'part'
     write(lun_outpu_dom,20) ipart
     write(lun_outpu_dom,10) parts_pos(ipart)%name
     write(lun_outpu_dom,10) 'coordinates'
     write(lun_outpu_dom,20) parts_pos(ipart)%npoin
     !
     ! Coordinates
     !
     do ipoin=1, parts_pos(ipart)%npoin
        write(lun_outpu_dom,20) ipoin
     end do
     do idime=1, ndime
        do ipoin=1, parts_pos(ipart)%npoin
           write(lun_outpu_dom,30) coord(idime,ipoin)
        end do
     end do
     if (ndime.eq.2) then
        xauxi= 0.0_rp
        do ipoin=1, parts_pos(ipart)%npoin
           write(lun_outpu_dom,30) xauxi
        end do
     end if
     !
     ! Volume elements
     !
     if(ndime==2) then
        iesta=10
        iesto=29
     else if(ndime==3) then
        iesta=30
        iesto=50
     end if
     do ielty=iesta,iesto
        if(lexis(ielty)>0) then
           cengo = cenal(ielty)
           if (cenal(ielty)=='tri3')  cengo = 'tria3'
           if (cenal(ielty)=='pyra5') cengo = 'pyramid5'
           write(lun_outpu_dom,10) trim(cengo)
           write(lun_outpu_dom,20) lnuty(ielty)
           do ielem=1,nelem
              if (ltype(ielem) == ielty) then
                 write(lun_outpu_dom,20) ielem                    
              end if
           end do
           do ielem=1,nelem
              if (ltype(ielem) == ielty) then
                 write(lun_outpu_dom,25) (lnods(inode,ielem) , inode=1,nnode(ielty))
              end if
           end do
        end if
     end do

  end do

10 format(a)
15 format(2a)
20 format(i10)
25 format(20i10)
30 format(e12.5)

end subroutine engold_wrmesh
