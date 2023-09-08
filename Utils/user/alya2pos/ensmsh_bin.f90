subroutine ensmsh_bin(&
     kfl_bound,mnode,mnodb,npoin,nelem,nboun,lun_asc,&
     lun_bin,lexis,ltype,lnods,lbxis,ltypb,lnodb,lelch,coord,&
     ndime,title,kfl_markm,npart_par,lsubd,nelem_par) 

  use,intrinsic :: ISO_C_BINDING
  use def_kintyp, only          :  ip,rp,cenam,nnode,cetop,cepos
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
  character(150)                :: fil_bin
  character(150)                :: creal,celem
  character(80)                 :: chaux
  integer(ip)                   :: iposi,istat,jelem,pelty
  !
  ! Initialize
  !
  ncoun_pos = 0
  nppti_ens = 0
  varnu_pos = 0
  do istpp=1,100
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
  write(100,'(a)') 'TIME' 
  write(100, 80  ) 'time set:               ',itise 
  write(100, 80  ) 'number of steps:        ',itise
  write(100,'(a)') 'filename start number:        1 '
  write(100,'(a)') 'filename increment:           1 '
  write(100,'(a)') 'time values: '
! arnau
!  write(100,'(10(1x,f0.5))') tiaux     
 write(100,'(10(1x,f5.0))') tiaux
  flush(100)
  !
  ! Write geometry
  !
  chens= adjustl(trim(title))
  chaux = 'C Binary'
  write(101) chaux  
  chaux = 'Problem name:  '//adjustl(trim(chens))
  write(101) chaux    
  chaux = 'Geometry file  '
  write(101) chaux
  chaux =  'node id given'
  write(101) chaux
  chaux =  'element id given'
  write(101) chaux

  do ipart=1,1
     chaux =  'part'
     write(101) chaux
     write(101) ipart
     chaux =  'Volume Mesh'
     write(101) chaux
     chaux =  'coordinates'
     write(101) chaux
     write(101) npoin
     !
     ! Coordinates
     !
     write(101) (ipoin, ipoin=1,npoin)
     do idime=1, ndime
        write(101) ( real(coord(idime,ipoin),4), ipoin=1, npoin)
     end do
     if (ndime.eq.2) then
        xauxi= 0.0_rp
        write(101) ( real(xauxi,4), ipoin=1, npoin)
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
     do ielty=1,nelty
        lnuty(ielty)=0
     end do
     do ielem = 1,nelem
        ielty = abs(ltype(ielem))
        lnuty(ielty)=lnuty(ielty)+1
     end do
     do ielty=iesta,iesto
        if(lexis(ielty)>0) then
           cengo = cenal(ielty)
           if (cenal(ielty)=='tri3')  cengo = 'tria3'
           if (cenal(ielty)=='pyra5') cengo = 'pyramid5'
           chaux =  trim(cengo)
           write(101) chaux
           write(101) lnuty(ielty)
           do ielem=1,nelem
              if (abs(ltype(ielem)) == ielty) then
                 write(101) ielem                    
              end if
           end do
           do ielem=1,nelem
              if (abs(ltype(ielem)) == ielty) then
                 write(101) (lnods(inode,ielem) , inode=1,nnode(ielty))
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
50 format(2a)
60 format(a,3x,i4,3x,a)
70 format(a,4x,i4,4x,a,4x,a)
80 format(a,3x,i4)

end subroutine ensmsh_bin
