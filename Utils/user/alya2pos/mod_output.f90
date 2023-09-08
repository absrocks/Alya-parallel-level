module mod_output

  use def_kintyp, only :  ip,rp,lg,cenam,nnode,cetop
  use mod_maths
  implicit none

contains
  
  subroutine output_gid_mesh(&
       mnode,npoin,nelem,lunit,lexis,ltype,lnods,coord,&
       gisca,lelch,leinv,ieles,ipois,ndime,ndimb,kfl_markm,flag_coh,title)
    

  implicit none

  integer(ip),    intent(in)           :: mnode
  integer(ip),    intent(in)           :: npoin
  integer(ip),    intent(in)           :: nelem
  integer(ip),    intent(in)           :: lunit
  integer(ip),    intent(in)           :: kfl_markm
  integer(ip),    intent(in)           :: ndime
  integer(ip),    intent(in)           :: ndimb
  integer(ip),    intent(in)           :: lnods(mnode,*)
  integer(ip),    intent(in)           :: lexis(*)
  integer(ip),    intent(in)           :: ltype(*)
  integer(ip),    intent(in)           :: leinv(*)
  integer(ip),    intent(in)           :: gisca(*)
  integer(ip),    intent(in)           :: lelch(*)
  integer(ip),    intent(inout)        :: ieles
  integer(ip),    intent(in)           :: ipois
  real(rp),       intent(in)           :: coord(ndime,*)
  character(*),   intent(in)           :: title
  logical(lg),    intent(in)           :: flag_coh
  integer(ip)                          :: idime,ipoin,inode,ielem,ielty
  integer(ip)                          :: iesto,iesta
  integer(ip),    save                 :: ifirs 
  character(150)                       :: dumml
  !
  ! Element range
  !
  if(      ndimb == 1 ) then
     iesta =  2
     iesto =  9
  else if( ndimb == 2 ) then
     iesta = 10
     iesto = 29
  else if( ndimb == 3 ) then
     iesta = 30
     iesto = 50
  end if
  !
  ! Mesh
  !
  do ielty = iesta,iesto

     if( lexis(ielty) /= 0 ) then

        dumml = adjustl(trim(title)//'_'//cenam(ielty))
        !
        ! Header
        !
        write(lunit,1)&
             adjustl(trim(dumml)),ndime,&
             adjustl(trim(cetop(ielty))),nnode(ielty)
        !
        ! Coordinates
        !
        if( ifirs == 0 .and. npoin > 0 ) then
           ifirs = 1
           write(lunit,2) 'coordinates'
           if( ndime == 1 ) then
              do ipoin = 1,npoin
                 write(lunit,3) ipoin+ipois,coord(1,ipoin),0.0_rp
              end do
           else
              do ipoin=1,npoin
                write(lunit,3) ipoin+ipois,(coord(idime,ipoin),idime=1,ndime)
              end do
           end if
           write(lunit,2) 'end coordinates'
        end if
        !
        ! Connectivity
        !
        write(lunit,2) 'elements'
        if( kfl_markm == 3 ) then
           do ielem = 1,nelem
              if( abs(ltype(ielem)) == ielty ) then
                 if( ltype(ielem) < 0 ) then
                    write(lunit,4) leinv(ielem)+ieles,&
                         (lnods(inode,ielem)+ipois,inode=1,nnode(ielty)),1000+gisca(ielem)
                 else
                    write(lunit,4) leinv(ielem)+ieles,&
                         (lnods(inode,ielem)+ipois,inode=1,nnode(ielty)),gisca(ielem)
                 end if
              end if
           end do
        else
           do ielem = 1,nelem
              if( abs(ltype(ielem)) == ielty .and. lelch(ielem) /= 17_ip ) then
                 write(lunit,4) leinv(ielem)+ieles,&
                      (lnods(inode,ielem)+ipois,inode=1,nnode(ielty)),gisca(ielem)
              end if
           end do
        end if
        write(lunit,2) 'end elements'
        write(lunit,2) ''

     end if

  end do
  !
  ! Mesh for cohesive elements
  !
  if ( flag_coh ) then
     ifirs = 0
     do ielty = iesta,iesto
        if( lexis(ielty) /= 0 ) then
           !
           ! Header
           !
           write(lunit,1)&
                'coh',max(2_ip,ndimb),&
                'Hexahedra',8_ip
           !
           ! Coordinates
           !
           if( ifirs == 0 .and. npoin > 0 ) then
              ifirs = 1
              write(lunit,2) 'coordinates'
              write(lunit,2) 'end coordinates'
           end if
           !
           ! Connectivity
           !
           write(lunit,2) 'elements'
           do ielem = 1,nelem
              if( abs(ltype(ielem)) == ielty .and. lelch(ielem) == 17_ip) then
                 write(lunit,4) leinv(ielem)+ieles,&
                      (lnods(inode,ielem)+ipois,inode=1,nnode(ielty)),gisca(ielem)
              end if
           end do
           write(lunit,2) 'end elements'
           write(lunit,2) ''
        end if
     end do
  end if

  ieles = ieles + nelem

1 format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2 format(a)
3 format(i11, 3(1x,e16.8e3))
4 format(i11,50(1x,i9))

end subroutine output_gid_mesh

end module mod_output
