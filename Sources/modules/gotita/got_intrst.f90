subroutine got_intrst()
  !------------------------------------------------------------------------
  !****f* Gotita/got_intrst
  ! NAME 
  !    got_intrst
  ! DESCRIPTION
  !    This routine interpolates initial boundary condition from a 
  !    formatted restart file
  ! USES
  !    runend
  !    elsest
  ! USED BY
  !    got_restar 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_gotita 
  use mod_memchk
  use mod_iofile
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip), parameter   :: mnode_b=27
  integer(ip), allocatable :: lnods_b(:,:)
  real(rp),    allocatable :: coord_b(:,:),vdrop_b(:,:),cdrop_b(:)
  integer(ip)              :: npoin_b,nelem_b,nnode_b,ltopo_b,iread,ielem
  integer(ip)              :: ipoin,idime,inode,ielem_b,inode_b,ipoin_b
  integer(ip)              :: ifoun,nelty_b
  integer(4)               :: istat
  real(rp)                 :: shape_b(mnode_b),deriv_b(3*mnode_b)
  real(rp)                 :: xjaci_b(9),xjacm_b(9),vdrop_f(3),cdrop_f
  real(rp)                 :: deltx_b(3),delts_b(3),xieta_b(3)
  real(rp)                 :: ltole_b,xbbox_b(9)
  !
  ! Read and write files
  !
  !lispa = 0
  !lisda = lun_rstar_got    ! Reading file
  !lisre = lun_outpu_got    ! Writing file
  call runend('GOT_INTRST: ERROR')
  !
  ! Read Dimensions
  !
  words(1)=''
  do while(words(1)/='DIMEN')
     call ecoute('got_intrst')
  end do
  do while(words(1)/='ENDDI') 
     call ecoute('got_intrst')
     if(words(1)=='NODAL') then
        npoin_b=getint('NODAL',0_ip,'*NUMBER OF NODAL POINTS')
     else if(words(1)=='ELEME') then
        nelem_b=getint('ELEME',0_ip,'*NUMBER OF ELEMENTS')
     else if(words(1)=='NODES') then
        nnode_b=getint('NODES',0_ip,'*NUMBER OF NODES PER ELEMENT')
     end if
  end do
  ! Allocate memory
  allocate(lnods_b(nnode_b,nelem_b),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'LNODS_B','got_intrst',lnods_b)
  allocate(coord_b(ndime,npoin_b),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'COORD_B','got_intrst',coord_b)
  allocate(vdrop_b(ndime,npoin_b),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'VDROP_B','got_intrst',vdrop_b)
  allocate(cdrop_b(npoin_b),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'CDROP_B','got_intrst',cdrop_b)
  ! Calculate ltopo_b
  ltopo_b = 1                                             ! Elements are simplices
  if(ndime==2) then
     if((nnode_b== 4).or.(nnode_b== 9).or.(nnode_b==16)) ltopo_b = 0
  else if(ndime==3) then
     if((nnode_b== 8).or.(nnode_b==27).or.(nnode_b==64)) ltopo_b = 0
     if((nnode_b== 6))                                   ltopo_b = 2
  end if
  iread=1
  !
  ! Read Geometry
  !
  do while(words(1)/='GEOME')
     call ecoute('got_intrst')
  end do
  if(iread==0) call runend('NSI_INTRST: WRONG INTERPOLATION FIELD.'&
       //'BACKGROUND DOMAIN DIMENSIONS WERE NOT DEFINED.')
  do while(words(1)/='ENDGE')
     call ecoute('got_intrst')
     if(words(1)=='ELEME') then
        call ecoute('got_intrst')
        do while(words(1)/='ENDEL')
           ielem=int(param(1))
           if(ielem<0.or.ielem>nelem_b) &
                call runend('NSI_INTRST: WRONG NUMBER OF ELEMENTS') 
           do inode=1,nnode_b
              lnods_b(inode,ielem)=int(param(inode+1))
           end do
           call ecoute('got_intrst')
        end do
     else if(words(1)=='COORD') then 
        call ecoute('got_intrst')
        do while(words(1)/='ENDCO')
           ipoin=int(param(1))
           if(ipoin<0.or.ipoin>npoin_b) &
                call runend('NSI_INTRST: WRONG NUMBER OF NODES') 
           do idime=1,ndime     
              coord_b(idime,ipoin)=param(idime+1)
           end do
           call ecoute('got_intrst')              
        end do
     end if
  end do
  !
  ! Read Results
  !
  do while(words(1)/='RESUL')
     call ecoute('got_intrst')
  end do
  do while(words(1)/='ENDRE')
     call ecoute('got_intrst')
     if(words(1)=='VDROP') then 
        call ecoute('got_intrst')
        do while(words(1)/='ENDVE')
           ipoin=int(param(1))
           if(ipoin<0.or.ipoin>npoin_b) &
                call runend('NSI_INTRST: WRONG NUMBER OF NODE') 
           do idime=1,ndime     
              vdrop_b(idime,ipoin)=param(idime+1)
           end do
           call ecoute('got_intrst')              
        end do
     else if(words(1)=='CDROP') then 
        call ecoute('got_intrst')
        do while(words(1)/='ENDPR')
           ipoin=int(param(1))
           if(ipoin<0.or.ipoin>npoin_b) &
                call runend('NSI_INTRST: WRONG NUMBER OF NODE') 
           cdrop_b(ipoin)=param(2)
           call ecoute('got_intrst')              
        end do
     end if
  end do
  !
  ! Change ltopo_b to nelty_b in Elsest format
  !
  if(ltopo_b==0) then
     nelty_b=1        ! Quadrilateral
  else
     nelty_b=2        ! Triangle
  end if
  !
  ! Loop over node with fixity flag=8
  !
  do ipoin=1,npoin
     ifoun=0
     ltole_b=1.0_rp
     !do while(ifoun==0.and.ltole_b<=100.0_rp)
     !   call Elsest(&
     !        one,    one,   zero,    two,&
     !        ielem_b,  ifoun,npoin_b,nnode_b,&
     !        nelem_b,nelty_b,nnode_b,  ndime,&     
     !        ip,     rp,   zero,   zero,&
     !        zero,   zero,   zero,   zero,&
     !        namod(modul),coord_b,lnods_b,coord(1,ipoin),&
     !        shape_b,deriv_b,xjacm_b,xjaci_b,&
     !        deltx_b,delts_b,xieta_b,ltole_b,&
     !        xbbox_b)
     !   ltole_b=ltole_b+10.0_rp
     !end do
     if(ifoun==0) then
        call runend('NSI_INTRST: COULD NOT FIND HOST ELEMENT IN BACKGROUND MESH'&
             //' FOR NODE '//adjustl(trim(intost(ipoin))))
     end if
     vdrop_f=0.0_rp 
     cdrop_f=0.0_rp
     do inode_b=1,nnode_b
        ipoin_b=lnods_b(inode_b,ielem_b)
        cdrop_f=cdrop_f+shape_b(inode_b)*cdrop_b(ipoin_b)
        do idime=1,ndime
           vdrop_f(idime)=vdrop_f(idime)+shape_b(inode_b)*vdrop_b(idime,ipoin_b)
        end do
     end do
     do idime=1,ndime
        vdrop(idime,ipoin,ncomp_got)=vdrop_f(idime)
     end do
     cdrop(ipoin,1)=cdrop_f
  end do
  !
  ! Deallocate Elsest memory.
  !
  !call Elsest(&
  !     two,    one,   zero,    two,&
  !     ielem_b,  ifoun,npoin_b,nnode_b,&
  !     nelem_b,nelty_b,nnode_b,  ndime,&     
  !     ip,     rp,   zero,   zero,&
  !     zero,   zero,   zero,   zero,&
  !     namod(modul),coord_b,lnods_b,coord(1,1),&
  !     shape_b,deriv_b,xjacm_b,xjaci_b,&
  !     deltx_b,delts_b,xieta_b,ltole_b,&
  !     xbbox_b)
  !
  ! Deallocate volatile memory
  !
  call memchk(two,istat,mem_modul(1:2,modul),'vdrop_b','got_intrst',vdrop_b)
  deallocate(vdrop_b,stat=istat)
  if(istat/=0)  call memerr(two,'VDROP_B','got_intrst',0_ip)   
  call memchk(two,istat,mem_modul(1:2,modul),'cdropb','got_intrst',cdrop_b)
  deallocate(cdrop_b,stat=istat)
  if(istat/=0)  call memerr(two,'CDROP_B','got_intrst',0_ip)   
  call memchk(two,istat,mem_modul(1:2,modul),'coord_b','got_intrst',coord_b)
  deallocate(coord_b,stat=istat)
  if(istat/=0)  call memerr(two,'COORD_B','got_intrst',0_ip)     
  call memchk(two,istat,mem_modul(1:2,modul),'lnods_b','got_intrst',lnods_b)
  deallocate(lnods_b,stat=istat)
  if(istat/=0)  call memerr(two,'LNODS_B','got_intrst',0_ip)

  lisda = momod(modul)%lun_pdata                                ! Reading file

end subroutine got_intrst
