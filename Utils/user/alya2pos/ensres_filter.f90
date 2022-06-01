subroutine ensres_filter(&
     ittim,npoin,nelem,parr1,pari1,lexis,&
     lbxis,ltype,pdime,namda,wopos,scavec,intrea,wbytes,kfl_filt,&
     rttim,kfl_markm,kfl_multi,gisc3,mnode,lnods) 


  use def_kintyp, only          :  ip,rp,varna_pos,varnu_pos,nunam_pos
  use def_kintyp, only          :  tipoe_ens,nppti_ens,ncoun_pos,nunam_pos
  use def_kintyp, only          :  pelpo,lelpo,nnode
  use def_elmtyp
  use def_kintyp, only          :  nmax_ensi

  implicit none
  integer(ip),    intent(in)    :: ittim,npoin,nelem,pdime
  integer(ip),    intent(in)    :: lexis(*)
  integer(ip),    intent(in)    :: lbxis(*)
  integer(ip),    intent(in)    :: ltype(*)
  real(rp),       intent(in)    :: rttim
  real(rp),       intent(in)    :: parr1(*)
  integer(ip),    intent(in)    :: pari1(*)
  character(150), intent(in)    :: namda
  character(5),   intent(in)    :: wopos(*)
  character(5),   intent(in)    :: scavec
  character(5),   intent(in)    :: kfl_filt
  character(5),   intent(in)    :: intrea
  character(5),   intent(in)    :: wbytes
  integer(ip),    intent(in)    :: kfl_markm
  integer(ip),    intent(in)    :: kfl_multi
  character(8)                  :: chtim
  character(20)                 :: wopo2(3)
  integer(ip)                   :: ipoin,ielem,ielty,idime,i,kpoin
  integer(ip)                   :: mpoin,nbyte,jelty,istpp,itise
  real(rp)                      :: dummr
  character(150)                :: filva
!-------------------------------------------------------
  integer(ip),    intent(in)    :: gisc3(*),mnode,lnods(mnode,*)
  integer(ip),    pointer       :: perm1(:),perm2(:),perm3(:),perm4(:)
  integer(ip)                   :: melem,ielpo,inode,jelem,kfl_pod
!-------------------------------------------------------

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
     perm2(ielem) = 0_ip
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
     perm3(ipoin) = 0_ip
     perm4(ipoin) = 0_ip
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


  
  !write(*,*) "npoin=",npoin
  !write(*,*) "nelem=",nelem
  !write(*,*) "mpoin=",mpoin
  !write(*,*) "melem=",melem
  !write(*,*) "!-----------------------------"
  
  !write(*,*) 'size(perm2)',size(perm2)
  !write(*,*) 'size(perm3)',size(perm3)
  !write(*,*) 'size(perm4)',size(perm4)
  !write(*,*) 'ittim=',ittim
  !write(*,*) '1- kttim=',kttim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  if( ittim >= 0 ) then

     istpp=1
     do while(istpp<=nmax_ensi)        
        if(wopos(1)==varna_pos(1,istpp)) then
           istpp=nmax_ensi
        else if(varna_pos(1,istpp)=='NULL') then
           varna_pos(1,istpp) = wopos(1)
           varna_pos(2,istpp) = wopos(2)
           varna_pos(3,istpp) = kfl_filt
           varnu_pos          = varnu_pos+1
           istpp=nmax_ensi
        end if
        istpp=istpp+1
     end do


!write(*,*)'varna_pos',varna_pos
!write(*,*)'istpp=',istpp
!write(*,*)'rttim=',rttim
!write(*,*)'tipoe_ens(istpp)=',tipoe_ens(istpp)
!write(*,*)'size(tipoe_ens)=',size(tipoe_ens)
!write(*,*)''

     istpp=1
     do while(istpp<=10000)        
        if(rttim==tipoe_ens(istpp)) then
           istpp=10000
        else if(tipoe_ens(istpp)==-1.0_rp) then
           tipoe_ens(istpp)= rttim
           !write(*,*)'opt2 nppti_ens=',nppti_ens
           nppti_ens = nppti_ens + 1
           !write(*,*)'opt2 nppti_ens=',nppti_ens
           istpp     = 10000
        end if
        istpp=istpp+1
     end do

     write(nunam_pos,'(i6)') nppti_ens      !   <<-- to write an integer to a character              
     if(nppti_ens<10) then
        write(nunam_pos,'(a,i1)') '00000',nppti_ens
     else if(nppti_ens<100) then
        write(nunam_pos,'(a,i2)') '0000',nppti_ens
     else if(nppti_ens<1000) then
        write(nunam_pos,'(a,i3)') '000',nppti_ens
     else if(nppti_ens<10000) then
        write(nunam_pos,'(a,i4)') '00',nppti_ens
     else if(nppti_ens<100000) then
        write(nunam_pos,'(a,i5)') '0',nppti_ens
     end if
     !
     ! Write result
     !
     ncoun_pos = ncoun_pos + 1
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
     !
     ! POD flag kfl_pod=1 active
     ! POD flag kfl_pod=0 active
     !
     kfl_pod=0
     if (kfl_pod==1) then
        open(unit=51,file='mat',form='unformatted',action='write',position='append')
     endif
     !
     filva = trim(namda)//'.ensi.'//trim(wopos(1))//'-'//trim(nunam_pos)
     open(unit=102,file=trim(filva),status='unknown',form='formatted')
     if( scavec == 'SCALA' ) then
        write(102,100) 'Alya Ensight Gold --- Scalar per-node variables file'
        write(102,100) 'part'
        write(102,110) 1_ip
        write(102,100) 'coordinates'
        do ipoin = 1,mpoin
              write(102,120) parr1(perm4(ipoin))
        end do
     else if( scavec == 'VECTO' ) then
        write(102,100) 'Alya Ensight Gold --- Vector per-node variables file'
        write(102,100) 'part'
        write(102,110) 1_ip
        write(102,100) 'coordinates'
        do idime = 1,pdime
           do ipoin = 1,mpoin
              !raul save my life ! one more time !!!!!!!(hadrien)
              write(102,120) parr1((perm4(ipoin)-1)*pdime+idime)
           end do
        end do
        if (pdime == 2) then             
           dummr = 0.0_rp
           do ipoin = 1,mpoin
              write(102,120) dummr
           end do
        end if
        !
        !
        !
        if (kfl_pod==1) then
           do idime = 1,pdime
              do ipoin = 1,mpoin
                 write(51) parr1((perm4(ipoin)-1)*pdime+idime)
              enddo
           enddo
        end if
        !
        !
        !
     end if
     close(unit=102)
     !
     if (kfl_pod==1) then
        close(51)
     endif
     !
  !else if( ittim == -1 ) then
     !
     ! Rewrite geometry file
     !
     !rewind(113)
     itise                = 1
     !tipoe_ens(nppti_ens) = cutim
     rewind(113)
     write(113,'(a)') '#' 
     write(113,'(a)') '# Alya generated post-process files' 
     write(113,'(a)') '# Ensight Gold Format' 
     write(113,'(a)') '#' 
     write(113,  50 ) '# Problem name:   ', adjustl(trim(namda))  
     write(113,'(a)') '#' 
     write(113,'(a)') 'FORMAT' 
     write(113,'(a)') 'type:    ensight gold' 
     write(113,'(a)') 'GEOMETRY' 
     write(113,  60 ) 'model:   ', itise,  adjustl(trim(namda))//'-filter.ensi.geo'           
     write(113,'(a)') 'VARIABLE' 
!    if (kfl_markm == 4) then
!       write(113,  50 ) 'scalar per element:   ', adjustl(trim(title))//'.ensi.LELCH'     
!    end if


     flush(113_4)
     ! 
     ! Write case end of file
     !
     itise = 1
     if( ncoun_pos /= 0 ) then
        do istpp=1,varnu_pos
           
           !
           ! write the generic file name
           !
           
           filva=adjustl(trim(namda))//'.ensi.'//trim(varna_pos(1,istpp))//'-******'
           if(varna_pos(2,istpp)=='VECTO'.and.varna_pos(3,istpp)=='FILTE') then
              write(113,70) 'vector per node:',itise,varna_pos(1,istpp),adjustl(trim(filva))          
           else if(varna_pos(2,istpp)=='SCALA'.and.varna_pos(3,istpp)=='FILTE') then
              write(113,70) 'scalar per node:',itise,varna_pos(1,istpp),adjustl(trim(filva))       
           end if
        end do
        write(113,'(a)') 'TIME' 
        write(113, 80  ) 'time set:               ',itise
        write(113, 80  ) 'number of steps:        ',nppti_ens
        write(113,'(a)') 'filename start number:        1 '
        write(113,'(a)') 'filename increment:           1 '
        write(113,'(a)') 'time values: '
        write(113,'(10(1x,e12.5))') (tipoe_ens(i),i=1,nppti_ens)           
        flush(113_4)

        !nppti_ens     = nppti_ens + 1
        !nppva_ens     = 0
        !kfl_statu_ens = 2
        !ncoun_pos     = 0
     end if
  end if
!----------------------------------------------------------------------
!
! Deallocate memory
!
!----------------------------------------------------------------------
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
80 format(a,3x,i8)
100 format(a)
110 format(i10)
120 format(e16.8E3)

end subroutine ensres_filter
