subroutine vu_filter(&
     kfl_bound,mnode,mnodb,npoin,nelem,nboun,lun_asc,&
     lun_bin,lexis,ltype,lnods,lbxis,ltypb,lnodb,lelch,coord,&
     ndime,namda,kfl_markm,npart_par,lsubd,nelem_par,&
     ittim,wopos,scavec,intrea,wbytes,&
     rttim,kfl_multi,gesc3,gisc3,geve3)


  use def_kintyp, only          :  ip,rp,cenam,nnode,cetop,cepos,nelty
  use def_kintyp, only          :  ltyp2,lllll,nllll,llll2,kttim,ncoun_pos
  use def_kintyp, only          :  pelpo,lelpo,intost
  use def_elmtyp

  implicit none
  integer(ip),    intent(in)    :: kfl_bound
  integer(ip),    intent(in)    :: mnode
  integer(ip),    intent(in)    :: npoin
  integer(ip),    intent(in)    :: nelem
  integer(ip),    intent(in)    :: lun_asc
  integer(ip),    intent(in)    :: lun_bin
  integer(ip),    intent(in)    :: ndime
  integer(ip),    intent(in)    :: mnodb,nboun
  integer(ip),    intent(in)    :: lnods(mnode,*)
  integer(ip),    intent(in)    :: lexis(*)
  integer(ip),    intent(in)    :: ltype(*)
  integer(ip),    intent(in)    :: lnodb(mnodb,*)
  integer(ip),    intent(in)    :: lbxis(*)
  integer(ip),    intent(in)    :: ltypb(*)
  integer(ip),    intent(in)    :: lelch(*)
  real(rp),       intent(in)    :: coord(ndime,*)
  character(150), intent(in)    :: namda
  integer(ip),    intent(in)    :: kfl_markm
  integer(ip),    intent(in)    :: npart_par
  integer(ip),    intent(in)    :: lsubd(*)
  integer(ip),    intent(in)    :: nelem_par(*)
  integer(ip),    intent(in)    :: ittim
  real(rp),       intent(in)    :: rttim
  character(5),   intent(in)    :: wopos
  character(5),   intent(in)    :: scavec
  character(5),   intent(in)    :: intrea
  character(5),   intent(in)    :: wbytes
  integer(ip),    intent(in)    :: kfl_multi
  integer(ip),    intent(in)    :: gisc3(*)
  real(rp),       intent(in)    :: gesc3(*)
  real(rp),       intent(in)    :: geve3(ndime,*)



  integer(ip)                   :: idime,ipoin,inode,ielem,pnode,ielty,ipart
  integer(ip)                   :: ifirs,ipoty,jboun,iblty,iboun,inodb,ii
  integer(ip)                   :: pnodb,jelty,ieles,nbyte,mpoin,pdime
  character(150)                :: fil_bin,fil_asc
  character(150)                :: creal,celem
  integer(ip)                   :: iposi,istat,jelem,pelty,test
  integer(ip), pointer          :: lnods_tmp(:,:),label_tmp(:)
  integer(ip), pointer          :: lnodb_tmp(:,:),labbo_tmp(:)
  character(8)                  :: chtim
  character(20)                 :: wopo2(3)
  integer(ip), pointer          :: perm1(:),perm2(:),perm3(:),perm4(:)
  integer(ip)                   :: melem,ielpo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  !write(*,*) "nelem=",nelem
  !write(*,*) "npoin=",npoin
  !write(*,*) "kfl_bound=",kfl_bound

  !do ipoin=1,npoin
  !   write(*,*) "ipoin=",ipoin
  !  write(*,*) "gisc3(ipoin)=",gisc3(ipoin)
  !end do   


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


  !write(*,*) "melem=",melem
  !write(*,*) "mpoin=",mpoin
  !write(*,*) 'size(perm2)',size(perm2)
  !write(*,*) 'size(perm3)',size(perm3)
  !write(*,*) 'size(perm4)',size(perm4)
  !write(*,*) 'ittim=',ittim
  !write(*,*) '1- kttim=',kttim
  !write(*,*)'nboun=',nboun
  !write(*,*)'mnodb=',mnodb
  !write(*,*)''
  !do iboun=1,nboun
  !   write(*,*)'lnodb(inodb,iboun)',(lnodb(inodb,iboun),inodb=1,mnodb) 
  !end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  ifirs = 0
  ipoty = 0

  pdime = ndime

  if( wbytes == '2BYTE' ) then
     nbyte = 2
  else if( wbytes == '4BYTE' ) then
     nbyte = 4
  else if( wbytes == '8BYTE' ) then
     nbyte = 8
  end if

  !
  ! File name
  !
  if(ittim<10) then
     write(chtim,'(a,i1)') '0000000',ittim
  else if(ittim<100) then
     write(chtim,'(a,i2)') '000000',ittim
  else if(ittim<1000) then
     write(chtim,'(a,i3)') '00000',ittim
  else if(ittim<10000) then
     write(chtim,'(a,i4)') '0000',ittim
  else if(ittim<100000) then
     write(chtim,'(a,i5)') '000',ittim
  else if(ittim<1000000) then
     write(chtim,'(a,i6)') '00',ittim
  else if(ittim<10000000) then
     write(chtim,'(a,i7)') '0',ittim
  end if
  if( kfl_multi == 0 ) then
     fil_asc = trim(namda)//'-'//trim(chtim)//'.filter.vu'
     fil_bin = trim(namda)//'-'//trim(chtim)//'.filter.vubin'
     if( ittim /= kttim ) then
        ncoun_pos = 0
        if( kttim /= -1 ) then
           close(lun_asc)
           close(lun_bin)
        end if
        open(unit=lun_asc,file=trim(fil_asc),form='formatted') 
        open(unit=lun_bin,file=trim(fil_bin),form='unformatted') 
     end if
  else
     ncoun_pos = 0
     fil_asc = trim(namda)//'-'//trim(wopos)//'-'//trim(chtim)//'.filter.vu'
     fil_bin = trim(namda)//'-'//trim(wopos)//'-'//trim(chtim)//'.filter.vubin'
     open(unit=lun_asc,file=trim(fil_asc),form='formatted') 
     open(unit=lun_bin,file=trim(fil_bin),form='unformatted') 
  end if

  

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
  if( ittim /= kttim ) then

     ncoun_pos=ncoun_pos+4_ip
     write(lun_asc,10) trim(creal),trim(fil_bin),mpoin*ndime,ncoun_pos 
     ncoun_pos=ncoun_pos+4_ip+mpoin*ndime*rp
     write(lun_bin) ((coord(idime,perm4(ipoin)),idime=1,ndime),ipoin=1,mpoin)


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
              do ielem=1,melem
                 if(ltyp2(perm2(ielem))==ielty) then
                    pnode=nnode(ltype(perm2(ielem)))
                    jelem=jelem+1
                    label_tmp(jelem)=ielem
                    do inode=1,pnode 
                       lnods_tmp(inode,jelem)=lnods(inode,perm2(ielem))
                    end do
                 end if
              end do
           end if
        end do
     else
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
     end if
     !
     ! Reorder boundaries
     ! 
     if(kfl_bound==1) then
        !
        !write(*,*)'nboun=',nboun
        !pnodb=3
        !do iboun=1,nboun
        !   do inodb=1,pnodb
        !      write(*,*)'lnodb(inodb,iboun)',lnodb(inodb,iboun)
        !   end do
        !end do
        !
        jboun=0
        do iblty=1,nelty
           if( lbxis(iblty) /= 0 ) then
              pnodb=nnode(iblty)
              !write(*,*)'nboun=',nboun
              !write(*,*)'pnodb',pnodb
              do iboun=1,nboun
                 if(ltypb(iboun)==iblty) then
                    jboun=jboun+1
                    labbo_tmp(jboun)=iboun
                    do inodb=1,pnodb
                       lnodb_tmp(inodb,jboun)=lnodb(inodb,iboun)
                       !write(*,*)'lnodb_tmp(inodb,jboun)',lnodb_tmp(inodb,jboun)
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
              ncoun_pos=ncoun_pos+4_ip
              celem=intost(ielty)
              !write(lun_asc,20) 'int',trim(celem),trim(fil_bin),lllll(ielty)*pnode,ncoun_pos
              write(lun_asc,20) 'int',trim(celem),trim(fil_bin),melem*pnode,ncoun_pos
              !ncoun_pos=ncoun_pos+4_ip+lllll(ielty)*pnode*ip
              ncoun_pos=ncoun_pos+4_ip+melem*pnode*ip
              !write(lun_bin) ((lnods_tmp(inode,ielem),inode=1,pnode),ielem=jelem,jelem+lllll(ielty)-1)
              write(lun_bin) ((lnods_tmp(inode,ielem),inode=1,pnode),ielem=jelem,jelem+melem-1)
              !jelem=jelem+lllll(ielty)
              jelem=jelem+melem
           end if
        end do
     else
        do ielty=1,nelty
           if(lexis(ielty)/=0) then
              pnode=nnode(ielty)
              ncoun_pos=ncoun_pos+4_ip
              celem=intost(ielty)
              !write(lun_asc,20) 'int',trim(celem),trim(fil_bin),lexis(ielty)*pnode,ncoun_pos
              write(lun_asc,20) 'int',trim(celem),trim(fil_bin),melem*pnode,ncoun_pos
              !ncoun_pos=ncoun_pos+4_ip+lexis(ielty)*pnode*ip
              ncoun_pos=ncoun_pos+4_ip+melem*pnode*ip
              !write(lun_bin) ((lnods_tmp(inode,ielem),inode=1,pnode),ielem=jelem,jelem+lexis(ielty)-1)
              write(lun_bin) ((perm3(lnods_tmp(inode,ielem)),inode=1,pnode),ielem=jelem,jelem+melem-1)
              !jelem=jelem+lexis(ielty)
              jelem=jelem+melem
           end if
        end do
     end if
     !
     ! Boundary 
     !
     if(kfl_bound==1) then
        write(*,*)'!----------'
        write(*,*)'debug'
        write(*,*)''
 
        !write(*,*)'nelty',nelty
         write(*,*)''
        jboun=1
        do iblty=1,nelty
           if(lbxis(iblty)/=0) then
              !write(*,*)'iblty',iblty
              !write(*,*)'lbxis(iblty)',lbxis(iblty)
              pnodb=nnode(iblty)
              ncoun_pos=ncoun_pos+4_ip
              celem=intost(iblty)
              !write(lun_asc,21) 'int',trim(celem),trim(fil_bin),lbxis(iblty)*pnodb,ncoun_pos
              write(lun_asc,21) 'int',trim(celem),trim(fil_bin),melem*pnodb,ncoun_pos
              !ncoun_pos=ncoun_pos+4_ip+lbxis(iblty)*pnodb*ip
              ncoun_pos=ncoun_pos+4_ip+melem*pnodb*ip
              !write(lun_bin) ((lnodb_tmp(inodb,iboun),inodb=1,pnodb),iboun=jboun,jboun+lbxis(iblty)-1)
!
!
              
              !write(*,*) ((lnodb_tmp(inodb,iboun),inodb=1,pnodb),iboun=jboun,jboun+melem-1)
              write(lun_bin) ((lnodb_tmp(inodb,iboun),inodb=1,pnodb),iboun=jboun,jboun+melem-1)
              !jboun=jboun+lbxis(iblty)
              jboun=jboun+melem
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
              ncoun_pos = ncoun_pos + 4_ip
              celem = intost(ielty)
              !write(lun_asc,25) 'int',trim(celem),trim(fil_bin),lllll(ielty),ncoun_pos
              !ncoun_pos = ncoun_pos + 4_ip + lllll(ielty)*ip
              !write(lun_bin) (label_tmp(ielem),ielem=jelem,jelem+lllll(ielty)-1) 
              !jelem = jelem + lllll(ielty)
              write(lun_asc,25) 'int',trim(celem),trim(fil_bin),melem,ncoun_pos
              ncoun_pos = ncoun_pos + 4_ip + melem*ip
              write(lun_bin) (label_tmp(ielem),ielem=jelem,jelem+melem-1) 
              jelem = jelem + melem
           end if
        end do
     else
        do ielty = 1,nelty
           if( lexis(ielty) /= 0 ) then
              ncoun_pos = ncoun_pos + 4_ip
              celem = intost(ielty)
              !write(lun_asc,25) 'int',trim(celem),trim(fil_bin),lexis(ielty),ncoun_pos
              !ncoun_pos = ncoun_pos + 4_ip + lexis(ielty)*ip
              !write(lun_bin) (label_tmp(ielem),ielem=jelem,jelem+lexis(ielty)-1)
              !jelem = jelem + lexis(ielty)
              write(lun_asc,25) 'int',trim(celem),trim(fil_bin),melem,ncoun_pos
              ncoun_pos = ncoun_pos + 4_ip + melem*ip
              write(lun_bin) (label_tmp(ielem),ielem=jelem,jelem+melem-1)
              jelem = jelem + melem
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
              ncoun_pos = ncoun_pos + 4_ip
              celem = intost(iblty)
              !write(lun_asc,26) 'int',trim(celem),trim(fil_bin),lbxis(iblty),ncoun_pos
              !ncoun_pos = ncoun_pos + 4_ip + lbxis(iblty)*ip
              !write(lun_bin) (labbo_tmp(iboun),iboun=jboun,jboun+lbxis(iblty)-1)
              !jboun = jboun + lbxis(iblty)
              write(lun_asc,26) 'int',trim(celem),trim(fil_bin),melem,ncoun_pos
              ncoun_pos = ncoun_pos + 4_ip + melem*ip
              write(lun_bin) (labbo_tmp(iboun),iboun=jboun,jboun+melem-1)
              jboun = jboun + melem
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
     write(lun_asc,30) 'MESH_'//trim(namda) 
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

     !----------------------------------------------------------------------
     ! 
     ! result part
     !
     !----------------------------------------------------------------------
  end if
  
  kttim=ittim
  write(lun_asc,205) rttim
  
  
  if( scavec == 'SCALA' ) then
     !
     ! SCALAR
     !
     ncoun_pos=ncoun_pos+4_ip 
     write(lun_asc,200) trim(wopos),trim(fil_bin),mpoin,ncoun_pos
     ncoun_pos=ncoun_pos+4_ip+nbyte*mpoin
     write(lun_asc,210) 
     if( kfl_markm >= 3 ) then
        do ielty=1,nllll
           if(lllll(ielty)/=0) then
              celem=intost(ielty)
              jelty=llll2(ielty)
              write(lun_asc,220) trim(wopos),&
                   trim(cepos(jelty)),trim(wopos),'Connec'//trim(celem),&
                   'Zone'//trim(cepos(jelty))//'_'//trim(celem)
           end if
        end do
     else
        do ielty=1,nelty
           if(lexis(ielty)/=0) then
              celem=intost(ielty)
              write(lun_asc,220) trim(wopos),&
                   trim(cepos(ielty)),trim(wopos),'Connec'//trim(celem),&
                   'Zone'//trim(cepos(ielty))//'_'//trim(celem)
           end if
        end do
     end if
     if( kfl_bound ==  1) then
        do ielty=1,nelty
           if(lbxis(ielty)/=0) then
              celem=intost(ielty)
              write(lun_asc,220) trim(wopos),&
                   trim(cepos(ielty)),trim(wopos),'Connecb'//trim(celem),&
                   'Zone'//trim(cepos(ielty))//'_'//trim(celem)
           end if
        end do
     end if

     write(lun_asc,230)
        if( intrea == 'INTEG' ) then
           write(lun_bin) (gesc3(perm4(ipoin)),ipoin=1,mpoin)
        else
           write(lun_bin) (gesc3(perm4(ipoin)),ipoin=1,mpoin)
        end if
    

  else if( scavec == 'VECTO' ) then
     !
     ! VECTOR
     !
     wopo2(1) = trim(wopos) // '_X'
     wopo2(2) = trim(wopos) // '_Y'
     wopo2(3) = trim(wopos) // '_Z'
     do idime = 1,pdime
        ncoun_pos=ncoun_pos+4_ip 
        write(lun_asc,200) trim(wopo2(idime)),trim(fil_bin),mpoin,ncoun_pos
        ncoun_pos=ncoun_pos+4_ip+nbyte*mpoin
        write(lun_asc,210) 
        if( kfl_markm >= 3 ) then
           do ielty=1,nllll
              if(lllll(ielty)/=0) then
                 celem=intost(ielty)
                 jelty=llll2(ielty)
                 write(lun_asc,220) trim(wopo2(idime)),&
                      trim(cepos(jelty)),trim(wopo2(idime)),'Connec'//trim(celem),&
                      'Zone'//trim(cepos(jelty))//'_'//trim(celem)
              end if
           end do
        else
           do ielty=1,nelty
              if(lexis(ielty)/=0) then
                 celem=intost(ielty)
                 write(lun_asc,220) trim(wopo2(idime)),&
                      trim(cepos(ielty)),trim(wopo2(idime)),'Connec'//trim(celem),&
                      'Zone'//trim(cepos(ielty))//'_'//trim(celem)
              end if
           end do
        end if
        if( kfl_bound ==  1) then
           do ielty=1,nelty
              if(lbxis(ielty)/=0) then
                 celem=intost(ielty)
                 write(lun_asc,220) trim(wopo2(idime)),&
                      trim(cepos(ielty)),trim(wopo2(idime)),'Connecb'//trim(celem),&
                      'Zone'//trim(cepos(ielty))//'_'//trim(celem)
              end if
           end do
        end if
        write(lun_asc,230)
        if( intrea == 'INTEG' ) then
           write(lun_bin) (gesc3((ipoin-1)*pdime+idime),ipoin=1,mpoin)
        else
           write(lun_bin) (gesc3((ipoin-1)*pdime+idime),ipoin=1,mpoin)
        end if
     end do
  end if

  deallocate(pelpo)
  deallocate(lelpo)
  deallocate(perm2)
  deallocate(perm3)  
  deallocate(perm4)

  if( kfl_multi == 1 ) then
     close(unit=lun_asc)
     close(unit=lun_bin)
  end if

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

200 format('FIELD<double> ',a,'("',a,'",',i7,',',i12,');')
201 format('FIELD<float>  ',a,'("',a,'",',i7,',',i12,');')
205 format('TEXTE Time(" ',e13.6,'");')
210 format('SOLUTION Solution( ) =',/,'{')
  !220 format('   VARIABLE ',a,'( ',a,',',a,',','Connec',i1,',Zone',i1,');')
220 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,');')
221 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,',',a,');')

230 format('};')



end subroutine vu_filter

