subroutine vu_msh(&
     kfl_bound,mnode,mnodb,npoin,nelem,nboun,lun_asc,&
     lun_bin,lexis,ltype,lnods,lbxis,ltypb,lnodb,lelch,coord,&
     ndime,title,kfl_markm,npart_par,lsubd,nelem_par) 

  use def_kintyp, only          :  ip,rp,cenam,nnode,cetop,cepos,nelty
  use def_kintyp, only          :  ltyp2,lllll,nllll,llll2,intost
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
  integer(ip)                   :: pnodb,jelty,ieles
  character(150)                :: fil_bin
  character(150)                :: creal,celem
  integer(ip)                   :: iposi,istat,jelem,pelty
  integer(ip), pointer          :: lnods_tmp(:,:),label_tmp(:)
  integer(ip), pointer          :: lnodb_tmp(:,:),labbo_tmp(:)

  ifirs = 0
  ipoty = 0
  fil_bin = trim(title)//'.msh.vubin'

  !----------------------------------------------------------------------
  ! 
  ! Geometry file header
  !
  !----------------------------------------------------------------------

  if( rp == 4 ) then
     creal = 'float'
  else
     creal = 'double'
  end if
  iposi = 0_ip

  !----------------------------------------------------------------------
  !
  ! Zone definition
  !
  !----------------------------------------------------------------------

  if( kfl_markm == 4 ) then

     nllll =  1000
     allocate( ltyp2(nelem) )
     allocate( lllll(nllll) )
     allocate( llll2(nllll) )
     do ielty = 1,nllll
        lllll(ielty) = 0
        llll2(ielty) = 0
     end do 
     do ielem = 1,nelem
        ielty = 10 * ltype(ielem) + lelch(ielem)  
        ltyp2(ielem) = ielty     
        lllll(ielty) = lllll(ielty) + 1
        llll2(ielty) = ltype(ielem)    
     end do

  else if( kfl_markm == 3 ) then

     nllll = npart_par
     allocate( ltyp2(nelem) )
     allocate( lllll(nllll) )
     allocate( llll2(nllll) )
     do ielty = 1,nllll
        lllll(ielty) = 0
        llll2(ielty) = 0
     end do 
     ieles = 0
     do ipart = 1,npart_par
        if( lsubd(ipart) == 1 ) then
           do ielem = 1,nelem_par(ipart)
              ieles        = ieles + 1
              ielty        = ipart
              ltyp2(ieles) = ielty
              lllll(ielty) = lllll(ielty) + 1
              llll2(ielty) = ltype(ieles)    
           end do
        end if
     end do

  end if
 
  !----------------------------------------------------------------------
  !
  ! Coordinates
  !
  !----------------------------------------------------------------------

  iposi=iposi+4_ip
  write(lun_asc,10) trim(creal),trim(fil_bin),npoin*ndime,iposi 
  iposi=iposi+4_ip+npoin*ndime*rp
  write(lun_bin) ((coord(idime,ipoin),idime=1,ndime),ipoin=1,npoin)

  !----------------------------------------------------------------------
  !
  ! Allocate memory
  !
  !----------------------------------------------------------------------

  allocate(lnods_tmp(mnode,nelem),stat=istat)
  allocate(label_tmp(nelem),stat=istat)
  if(kfl_bound==1) then
     allocate(lnodb_tmp(mnodb,nboun),stat=istat)
     allocate(labbo_tmp(nboun),stat=istat)
  end if

  !----------------------------------------------------------------------
  !
  ! Reordering
  !
  !----------------------------------------------------------------------
  !
  ! Reorder elements: Zone is element type
  !
  jelem=0

  if( kfl_markm >= 3 ) then
     do ielty=1,nllll
        if( lllll(ielty) /= 0 ) then
           do ielem=1,nelem
              if(ltyp2(ielem)==ielty) then
                 pnode=nnode(ltype(ielem))
                 jelem=jelem+1
                 label_tmp(jelem)=ielem
                 do inode=1,pnode 
                    lnods_tmp(inode,jelem)=lnods(inode,ielem)
                 end do
              end if
           end do
        end if
     end do
  else
     do ielty=1,nelty
        if( lexis(ielty) /= 0 ) then
           pnode=nnode(ielty)
           do ielem=1,nelem
              if(ltype(ielem)==ielty) then
                 jelem=jelem+1
                 label_tmp(jelem)=ielem
                 do inode=1,pnode 
                    lnods_tmp(inode,jelem)=lnods(inode,ielem)
                 end do
              end if
           end do
        end if
     end do
  end if
  !
  ! Reorder boundaries
  ! 
  if(kfl_bound==1) then
     jboun=0
     do iblty=1,nelty
        if( lbxis(iblty) /= 0 ) then
           pnodb=nnode(iblty)
           do iboun=1,nboun
              if(ltypb(iboun)==iblty) then
                 jboun=jboun+1
                 labbo_tmp(jboun)=iboun
                 do inodb=1,pnodb
                    lnodb_tmp(inodb,jboun)=lnodb(inodb,iboun)
                 end do
              end if
           end do
        end if
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Connectivity in *.msh.vubin file
  !
  !----------------------------------------------------------------------

  !
  ! Element: Zone is element type
  !
  jelem=1
  if( kfl_markm >= 3 ) then
     do ielty=1,nllll
        if(lllll(ielty)/=0) then
           pnode=nnode(llll2(ielty))
           iposi=iposi+4_ip
           celem=intost(ielty)
           write(lun_asc,20) 'int',trim(celem),trim(fil_bin),lllll(ielty)*pnode,iposi
           iposi=iposi+4_ip+lllll(ielty)*pnode*ip
           write(lun_bin) ((lnods_tmp(inode,ielem),inode=1,pnode),ielem=jelem,jelem+lllll(ielty)-1)
           jelem=jelem+lllll(ielty)
        end if
     end do
  else
     do ielty=1,nelty
        if(lexis(ielty)/=0) then
           pnode=nnode(ielty)
           iposi=iposi+4_ip
           celem=intost(ielty)
           write(lun_asc,20) 'int',trim(celem),trim(fil_bin),lexis(ielty)*pnode,iposi
           iposi=iposi+4_ip+lexis(ielty)*pnode*ip
           write(lun_bin) ((lnods_tmp(inode,ielem),inode=1,pnode),ielem=jelem,jelem+lexis(ielty)-1)
           jelem=jelem+lexis(ielty)
        end if
     end do
  end if
  !
  ! Boundary 
  !
  if(kfl_bound==1) then
     jboun=1
     do iblty=1,nelty
        if(lbxis(iblty)/=0) then
           pnodb=nnode(iblty)
           iposi=iposi+4_ip
           celem=intost(iblty)
           write(lun_asc,21) 'int',trim(celem),trim(fil_bin),lbxis(iblty)*pnodb,iposi
           iposi=iposi+4_ip+lbxis(iblty)*pnodb*ip
           write(lun_bin) ((lnodb_tmp(inodb,iboun),inodb=1,pnodb),iboun=jboun,jboun+lbxis(iblty)-1)
           jboun=jboun+lbxis(iblty)
        end if
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Labeling 
  !
  !----------------------------------------------------------------------
  !
  ! Labeling of the elements: Zone is element type
  !
  jelem = 1
  if( kfl_markm >= 3 ) then
     do ielty = 1,nllll
        if( lllll(ielty) /= 0 ) then
           iposi = iposi + 4_ip
           celem = intost(ielty)
           write(lun_asc,25) 'int',trim(celem),trim(fil_bin),lllll(ielty),iposi
           iposi = iposi + 4_ip + lllll(ielty)*ip
           write(lun_bin) (label_tmp(ielem),ielem=jelem,jelem+lllll(ielty)-1) 
           jelem = jelem + lllll(ielty)
        end if
     end do
  else
     do ielty = 1,nelty
        if( lexis(ielty) /= 0 ) then
           iposi = iposi + 4_ip
           celem = intost(ielty)
           write(lun_asc,25) 'int',trim(celem),trim(fil_bin),lexis(ielty),iposi
           iposi = iposi + 4_ip + lexis(ielty)*ip
           write(lun_bin) (label_tmp(ielem),ielem=jelem,jelem+lexis(ielty)-1) 
           jelem = jelem + lexis(ielty)
        end if
     end do
  end if
  !
  ! Labeling of the boundaries
  !
  if( kfl_bound == 1 ) then
     jboun = 1
     do iblty = 1,nelty
        if( lbxis(iblty) /= 0 ) then
           iposi = iposi + 4_ip
           celem = intost(iblty)
           write(lun_asc,26) 'int',trim(celem),trim(fil_bin),lbxis(iblty),iposi
           iposi = iposi + 4_ip + lbxis(iblty)*ip
           write(lun_bin) (labbo_tmp(iboun),iboun=jboun,jboun+lbxis(iblty)-1)
           jboun = jboun + lbxis(iblty)
        end if
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Header for *.msh.vu file, definition of zones
  !
  !----------------------------------------------------------------------
  !
  ! Mesh zones 
  !
  write(lun_asc,30) 'MESH_'//trim(title) 
  write(lun_asc,40)
  !
  ! Zone is element type
  !
  if( kfl_markm >= 3 ) then
     do ielty=1,nllll
        if(lllll(ielty)/=0) then
           jelty=llll2(ielty)
           celem=intost(ielty)
           write(lun_asc,50) trim(cepos(jelty))//'_'//trim(celem),trim(cepos(jelty)),ndime,trim(celem),trim(celem)
        end if
     end do
  else
     do ielty=1,nelty
        if(lexis(ielty)/=0) then
           celem=intost(ielty)
           write(lun_asc,50) trim(cepos(ielty))//'_'//trim(celem),trim(cepos(ielty)),ndime,trim(celem),trim(celem)
        end if
     end do
  end if
  if(kfl_bound==1) then
     do iblty=1,nelty
        if(lbxis(iblty)/=0) then
           celem=intost(iblty)
           write(lun_asc,51) trim(cepos(iblty))//'_'//trim(celem),trim(cepos(iblty)),ndime,trim(celem),trim(celem)
        end if
     end do
  end if

  write(lun_asc,60)

  !----------------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !----------------------------------------------------------------------

  if( kfl_bound == 1 ) then
     deallocate(lnodb_tmp,stat=istat)
     deallocate(labbo_tmp,stat=istat)
  end if
  deallocate(lnods_tmp,stat=istat)
  deallocate(label_tmp,stat=istat)
  !
  ! Format
  !
10 format('FIELD<',a,'> Coord("',a,'",',i12,',',i1,');')
20 format('FIELD<',a,'> Connec',a,'("',a,'",',i12,',',i12,');')
21 format('FIELD<',a,'> Connecb',a,'("',a,'",',i12,',',i12,');')
25 format('FIELD<',a,'> Numer',a,'("',a,'",',i12,',',i12,');')
26 format('FIELD<',a,'> Numerb',a,'("',a,'",',i12,',',i12,');')
30 format('MESH ',a,'( ) =')
40 format('{')
50 format(' ZONE Zone',a,'( ',a,',Coord%',i1,',Connec',a,',,Numer',a,' );')
51 format(' ZONE Zone',a,'( ',a,',Coord%',i1,',Connecb',a,',,Numerb',a,' );')
60 format('};')

end subroutine vu_msh


