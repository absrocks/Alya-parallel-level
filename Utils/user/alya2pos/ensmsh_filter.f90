subroutine ensmsh_filter(&
     kfl_bound,mnode,mnodb,npoin,nelem,nboun,lun_asc,&
     lun_bin,lexis,ltype,lnods,lbxis,ltypb,lnodb,lelch,coord,&
     ndime,title,kfl_markm,npart_par,lsubd,nelem_par,&
     gesc3,gisc3,geve3) 

  use def_kintyp, only          :  ip,rp,cenam,nnode,cetop,cepos
  use def_kintyp, only          :  nelty,cenal,lnuty,varna_pos,varnu_pos
  use def_kintyp, only          :  tipoe_ens,nppti_ens,ncoun_pos
  use def_kintyp, only          :  pelpo,lelpo
  use def_elmtyp
  use def_kintyp, only          :  nmax_ensi

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
  integer(ip),    intent(in)    :: gisc3(*)
  real(rp),       intent(in)    :: gesc3(*)
  real(rp),       intent(in)    :: geve3(ndime,*)
!-------------------------------------------------------------------------------
  integer(ip)                   :: idime,ipoin,inode,ielem,pnode,ielty,ipart
  integer(ip)                   :: ifirs,ipoty,jboun,iblty,iboun,inodb,ii
  integer(ip)                   :: pnodb,jelty,ieles,istpp
  integer(ip)                   :: iesta,iesto,itise
  real(rp)                      :: xauxi,tiaux
  character(40)                 :: chens
  character(10)                 :: cengo
  character(150)                :: fil_bin
  character(150)                :: creal,celem
  integer(ip)                   :: iposi,istat,jelem,pelty
!---------------------------------------------------------------------------------
  integer(ip), pointer          :: perm1(:),perm2(:),perm3(:),perm4(:)
  integer(ip)                   :: melem,ielpo,mpoin
  integer(ip), pointer          :: lnods_tmp(:,:),label_tmp(:)
  integer(ip), pointer          :: lnodb_tmp(:,:),labbo_tmp(:)
!----------------------------------------------------------------------------------
  !write(*,*)'npoin=',npoin
  !write(*,*)'nelem=',nelem
  
              
                       
  allocate(perm1(nelem))
  do ielem = 1,nelem
     perm1(ielem) = 0_ip
  end do

  call elmtyp()
  call connpo(ielty,ltype,nelem,npoin,mnode,lnods)

  !do ipoin=1,npoin
  !   do ielpo=pelpo(ipoin),pelpo(ipoin+1)-1
  !      jelem=lelpo(ielpo)
  !   write(*,*) ipoin , jelem
  !   end do
  !end do 

  do ipoin=1,npoin
     if (gisc3(ipoin) /= 0) then
        !write(*,*) "ipoin=",ipoin
        !do ielpo=pelpo(gisc3(ipoin)),pelpo(gisc3(ipoin+1))-1
        do ielpo=pelpo(ipoin),pelpo(ipoin+1)-1   
           jelem=lelpo(ielpo)
           !write(*,*) "jelem=",jelem
           perm1(jelem)=1_ip
        end do
     end if
  end do

  melem=0

  do ielem=1,nelem
     if (perm1(ielem) /= 0) then
        melem=melem+1
     end if
  end do

  allocate(perm2(melem))
  do ielem = 1,melem
     perm2(ielem) = 0
  end do

  melem=0
  do ielem=1,nelem
     if (perm1(ielem) /= 0 ) then
        melem=melem+1
        perm2(melem)=ielem
     end if
  end do

  deallocate(perm1)  

  allocate(perm3(npoin))
  allocate(perm4(npoin))
  do ipoin = 1,npoin
     perm3(ipoin) = 0
     perm4(ipoin) = 0
  end do


  mpoin=0

  do jelem = 1,melem                            
     ielem = perm2(jelem)
     do inode = 1,nnode(ltype(ielem))
        ipoin=lnods(inode,ielem)
        if(perm3(ipoin)==0)then 
           mpoin=mpoin+1
           perm3(ipoin)= mpoin
           perm4(mpoin)= ipoin
        end if
     end do
  end do

  !write(*,*) "mpoin=",mpoin
  !write(*,*) "melem=",melem
  !write(*,*) "!-----------------------------"
  
  !write(*,*) 'size(perm2)',size(perm2)
  !write(*,*) 'size(perm3)',size(perm3)
  !write(*,*) 'size(perm4)',size(perm4)
  !write(*,*) 'ittim=',ittim
  !write(*,*) '1- kttim=',kttim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
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
  write(113,'(a)') '#' 
  write(113,'(a)') '# Alya generated post-process files' 
  write(113,'(a)') '# Ensight Gold Format' 
  write(113,'(a)') '#' 
  write(113,  50 ) '# Problem name:   ', adjustl(trim(title))
  write(113,'(a)') '#' 
  write(113,'(a)') 'FORMAT' 
  write(113,'(a)') 'type:    ensight gold' 
  write(113,'(a)') 'GEOMETRY' 
  write(113,  60 ) 'model:   ', itise, adjustl(trim(title))//'-filter.ensi.geo'
  if (kfl_markm == 4) then
     write(113,'(a)') 'VARIABLE' 
     write(113,  50 ) 'scalar per element:   ', adjustl(trim(title))//'.ensi.LELCH'     
  end if
  write(113,'(a)') 'TIME' 
  write(113, 80  ) 'time set:               ',itise 
  write(113, 80  ) 'number of steps:        ',itise
  write(113,'(a)') 'filename start number:        1 '
  write(113,'(a)') 'filename increment:           1 '
  write(113,'(a)') 'time values: '
! arnau
!  write(113,'(10(1x,f0.5))') tiaux     
  write(113,'(10(1x,f5.0))') tiaux
  flush(113)
  !
  ! Write geometry
  !
  chens= adjustl(trim(title))
  write(114,15) 'Problem name:  ',adjustl(trim(chens))
  write(114,10) 'Geometry file  '
  write(114,10) 'node id given'
  write(114,10) 'element id given'
  do ipart=1,1
     write(114,10) 'part'
     write(114,20) ipart
     write(114,10) 'Volume Mesh'
     write(114,10) 'coordinates'
     write(114,20) mpoin
     !
     ! Coordinates
     !
     do ipoin=1, mpoin
        write(114,20) ipoin
     end do
     do idime=1, ndime
        do ipoin=1, mpoin
           write(114,30) coord(idime,perm4(ipoin))
        end do
     end do
     if (ndime.eq.2) then
        xauxi= 0.0_rp
        do ipoin=1, mpoin
           write(114,30) xauxi
        end do
     end if
     !
     ! Volume elements
     !
!----------------------------------------------------------------------
!
! Allocate memory
!
!----------------------------------------------------------------------
           
     allocate(lnods_tmp(mnode,melem),stat=istat)
     allocate(label_tmp(melem),stat=istat)
     if(kfl_bound==1) then
        allocate(lnodb_tmp(mnodb,nboun),stat=istat)
        allocate(labbo_tmp(nboun),stat=istat)
     end if

!----------------------------------------------------------------------
!
!  Reordering
!
!----------------------------------------------------------------------
!
! Reorder elements: Zone is element type
!
     jelem=0
        do ielty=1,nelty 
        if( lexis(ielty) /= 0 ) then
           pnode=nnode(ielty)
           do ielem=1,melem
              if(ltype(perm2(ielem))==ielty) then
                 jelem=jelem+1
                 label_tmp(jelem)=ielem
                 do inode=1,pnode 
                    lnods_tmp(inode,jelem)=lnods(inode,perm2(ielem))
                 end do
              end if
           end do
        end if
     end do
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
!-----------------------------------------------------------------------------------

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
     do ielem = 1,melem
        ielty = abs(ltype(perm2(ielem)))
        lnuty(ielty)=lnuty(ielty)+1
     end do
     do ielty=iesta,iesto
        if(lexis(ielty)>0) then
           cengo = cenal(ielty)
           if (cenal(ielty)=='tri3')  cengo = 'tria3'
           if (cenal(ielty)=='pyra5') cengo = 'pyramid5'
           write(114,10) trim(cengo)
           write(114,20) lnuty(ielty)
           do ielem=1,melem
!!!              if (ltype(ielem) == ielty) then          antes estaba asi pero fallan los contactos... no se por que era asi
              if (abs(ltype(perm2(ielem))) == ielty) then
                 write(114,20) ielem                    
              end if
           end do
           do ielem=1,melem
              if (abs(ltype(perm2(ielem))) == ielty) then
                 write(114,25) (perm3(lnods(inode,perm2(ielem))) , inode=1,nnode(ielty))
                 !write(*,*) (perm3(lnods(inode,perm2(ielem))) , inode=1,nnode(ielty))
              end if
           end do
        end if
     end do

     if (kfl_markm == 4) then
        
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
        do ielty=iesta,iesto        
           if(lexis(ielty)>0) then
              cengo = cenal(ielty)
              if (cenal(ielty)=='tri3')  cengo = 'tria3'
              if (cenal(ielty)=='pyra5') cengo = 'pyramid5'
              write(114,10) 'elements  ', trim(cengo)
              do ielem=1,melem
                 if (abs(ltype(perm2(ielem))) == ielty) then
                    write(114,20) lelch(ielem)                    
                 end if
              end do
           end if
        end do
        

     end if

  end do
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
  !--------------------------------------------------------------
  deallocate(pelpo)
  deallocate(lelpo)
  deallocate(perm2)
  deallocate(perm3)  
  deallocate(perm4)
!---------------------------------------------------------------
10 format(a)
15 format(2a)
20 format(i10)
25 format(20i10)
30 format(e12.5)
50 format(2a)
60 format(a,3x,i4,3x,a)
70 format(a,4x,i4,4x,a,4x,a)
80 format(a,3x,i4)

end subroutine ensmsh_filter
