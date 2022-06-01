subroutine ensmsh(&
     kfl_bound,mnode,mnodb,npoin,nelem,nboun,lun_asc,&
     lun_bin,lexis,ltype,leinv,lnods,lbxis,ltypb,lnodb,lelch,coord,&
     ndime,title,kfl_markm,npart_par,lsubd,nelem_par)

  use def_kintyp, only          :  ip,rp,cenam,nnode,cetop,cepos
  use def_kintyp, only          :  nmax_ensi
  use def_kintyp, only          :  nelty,cenal,lnuty,varna_pos,varnu_pos
  use def_kintyp, only          :  tipoe_ens,nppti_ens,ncoun_pos
  use def_elmtyp

  implicit none
  integer(ip),    intent(in)    :: kfl_bound
  integer(ip),    intent(in)    :: mnode,npoin,nelem,lun_asc,lun_bin,ndime
  integer(ip),    intent(in)    :: mnodb,nboun
  integer(ip),    intent(in)    :: lnods(mnode,*)
  integer(ip),    intent(in)    :: lexis(*)
  integer(ip),    intent(in)    :: ltype(*)
  integer(ip),    intent(in)    :: leinv(*)
  integer(ip),    intent(in)    :: lnodb(mnodb,*)
  integer(ip),    intent(in)    :: lbxis(*)
  integer(ip),    intent(in)    :: ltypb(*)
  integer(ip),    intent(in)    :: lelch(*)
  real(rp),       intent(in)    :: coord(ndime,*)
  character(150), intent(in)    :: title
  integer(ip),    intent(in)    :: kfl_markm
  integer(ip),    intent(in)    :: npart_par
  integer(ip),    intent(in)    :: lsubd(*)
  integer(ip),    intent(in)    :: nelem_par(*)
  integer(ip)                   :: idime,ipoin,inode,ielem,pnode,ielty,ipart
  integer(ip)                   :: ifirs,ipoty,jboun,iblty,iboun,inodb,ii
  integer(ip)                   :: pnodb,jelty,ieles,istpp
  integer(ip)                   :: iesta,iesto,itise
  real(rp)                      :: xauxi,tiaux
  character(40)                 :: chens
  character(10)                 :: cengo
  
  !
  ! Initialize
  !
  ncoun_pos = 0
  nppti_ens = 0
  varnu_pos = 0
  do istpp=1,nmax_ensi
     varna_pos(1,istpp) = 'NULL'
     varna_pos(2,istpp) = 'NULL'
  end do
  do istpp=1,10000
     tipoe_ens(istpp) = -1.0_rp
  end do
  !
  ! Write a prelimiary case file just to postprocess the mesh (no variable)
  !
  itise = 1
  tiaux = 0.0_rp
  write(100,'(a)') '#'
  write(100,'(a)') '# Alya generated post-process files'
  write(100,'(a)') '# Ensight Gold Format'
  write(100,'(a)') '#'
  write(100,  50 ) '# Problem name:   ', adjustl(trim(title))
  write(100,'(a)') '#'
  write(100,'(a)') 'FORMAT'
  write(100,'(a)') 'type:    ensight gold'
  write(100,'(a)') 'GEOMETRY'
  write(100,  60 ) 'model:   ', itise, adjustl(trim(title))//'.ensi.geo'
  !write(100,  60 ) '#model:   ', itise, adjustl(trim(title))//'-filter.ensi.geo'
  if (kfl_markm == 4) then
     write(100,'(a)') 'VARIABLE'
     write(100,  50 ) 'scalar per element:   ', adjustl(trim(title))//'.ensi.LELCH'
  end if
  write(100,'(a)') 'TIME'
  write(100, 80  ) 'time set:               ',itise
  write(100, 80  ) 'number of steps:        ',itise
  write(100,'(a)') 'filename start number:        1 '
  write(100,'(a)') 'filename increment:           1 '
  write(100,'(a)') 'time values: '
! arnau
!  write(100,'(10(1x,f0.5))') tiaux
  write(100,'(10(1x,f5.0))') tiaux
  flush(100_4)
  !
  ! Write geometry
  !
  chens= adjustl(trim(title))
  write(101,15) 'Problem name:  ',adjustl(trim(chens))
  write(101,10) 'Geometry file  '
  write(101,10) 'node id given'
  write(101,10) 'element id given'
  do ipart=1,1
     write(101,10) 'part'
     write(101,20) ipart
     write(101,10) 'Volume Mesh'
     write(101,10) 'coordinates'
     write(101,20) npoin
     !
     ! Coordinates
     !
     do ipoin=1, npoin
        write(101,20) ipoin
     end do
     do idime=1, ndime
        do ipoin=1, npoin
           write(101,30) coord(idime,ipoin)
        end do
     end do
     if (ndime.eq.2) then
        xauxi= 0.0_rp
        do ipoin=1, npoin
           write(101,30) xauxi
        end do
     end if
     !
     ! Elements
     !
     if(ndime==2) then
        iesta=10
        iesto=29
     else if(ndime==3) then
        iesta=30
        iesto=50
     end if
     do ielty=1,nelty
        lnuty(ielty)=0
     end do
     do ielem = 1,nelem
        ielty = abs(ltype(ielem))
        lnuty(ielty)=lnuty(ielty)+1
     end do
     !
     ! Loop element types
     !
     do ielty=iesta,iesto
        if(lexis(ielty)>0) then
           cengo = cenal(ielty)
           if (cenal(ielty)=='tri3')  cengo = 'tria3'
           if (cenal(ielty)=='pyra5') cengo = 'pyramid5'
           ! Write element type
           write(101,10) trim(cengo)
           write(101,20) lnuty(ielty)
           ! Write element list
           do ielem=1,nelem
              if (abs(ltype(ielem)) == ielty) then
                 write(101,20) leinv(ielem)
              end if
           end do
           ! Write element node numbering
           do ielem=1,nelem
              if (abs(ltype(ielem)) == ielty) then
                 write(101,25) (lnods(inode,ielem) , inode=1,nnode(ielty))
              end if
           end do
        end if
     end do

     if (kfl_markm == 4) then
        !
        ! Mark type: ELCH
        !
        write(202,10) 'Alya Ensight Gold --- Scalar per-element variables file'
        write(202,10) 'part'
        write(202,20) 1_ip
        if(ndime==2) then
           iesta=10
           iesto=29
        else if(ndime==3) then
           iesta=30
           iesto=50
        end if
        !
        ! Loop element library
        !
        do ielty=iesta,iesto
           if(lexis(ielty)>0) then
              cengo = cenal(ielty)
              if (cenal(ielty)=='tri3')  cengo = 'tria3'
              if (cenal(ielty)=='pyra5') cengo = 'pyramid5'
              ! Write element type
              write(202,10) 'elements  ', trim(cengo)
              ! Write element characteristic
              do ielem=1,nelem
                 if (abs(ltype(ielem)) == ielty) then
                    write(202,20) lelch(leinv(ielem))
                 end if
              end do
           end if
        end do

     end if

  end do

10 format(a)
15 format(2a)
20 format(i10)
25 format(20i10)
30 format(e15.8)
!!! 30 format(e12.5)    ! this is the correct format from the ensight manual. try it when problems appear.
50 format(2a)
60 format(a,3x,i4,3x,a)
70 format(a,4x,i4,4x,a,4x,a)
80 format(a,3x,i4)

end subroutine ensmsh
