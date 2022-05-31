subroutine mescek(itask)
  !-----------------------------------------------------------------------
  !****f* domain/mescek
  ! NAME
  !    mescek
  ! DESCRIPTION
  !    This routine checks if the mesh is correct
  ! OUTPUT
  !    VODOM ... Total domain volume
  !    VOAVE ... Averaged domain volume
  !    VOMIN ... Minimum element volume
  !    VOMAX ... Maximum element volume
  !    ELMIN ... Element with minimum volume
  !    ELMAX ... Element with maximum volume
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_domain
  use def_master
  use def_kermod
  use mod_memchk
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_MAX
  use mod_outfor,         only : outfor
  use mod_messages,       only : messages_live
  use mod_iofile,         only : iofile_flush_unit
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  
  implicit none
  
  integer(ip), intent(in)  :: itask
  integer(ip)              :: ielem,inode,pgaus,pnode,pelty,igaus
  integer(ip)              :: iboun,inodb,pgaub,pnodb,pblty,igaub
  integer(ip)              :: kstar,kzero,keror,ieror,ipoin,npoit,idime
  integer(ip)              :: neror(7),jpoin,iblty,nelet,izdom,pmate
  integer(ip)              :: iqual,kelem,istack,nstack,imesh,kpoin
  integer(ip)              :: kfl_paral_min,kfl_paral_max
  real(rp)                 :: detjm,volum,dummr,surfa
  real(rp)                 :: quali,qmini,qmaxi
  real(rp)                 :: cartd(ndime,mnode)
  real(rp)                 :: elcod(ndime,mnode),baloc(ndime,ndime)
  real(rp)                 :: xjaci(ndime,ndime),xjacm(ndime,ndime)
  character(10)            :: mess1,mess2
  logical(lg), pointer     :: touch(:)

  neror = 0
  nullify(touch)
  
  select case ( itask )

  case ( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Element Jacobian
     !
     !-------------------------------------------------------------------

     call messages_live('CHECK ELEMENT ORDERING')
     if( associated(vomat) ) then
        call memory_deallo(memor_dom,'VOMAT','mescek' ,vomat)
     end if
     call memory_alloca(memor_dom,'VOMAT','mescek' ,vomat,nmate)

     if( INOTMASTER ) then
        !
        ! Check element ordering by computing Jacobian sign
        !
        vomin    = huge(1.0_rp)
        vomax    = 0.0_rp
        vodom    = 0.0_rp

        do ielem = 1,nelem
           pelty = ltype(ielem)
           if( pelty > 0 ) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              pmate = lmate(ielem)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do
              volum = 0.0_rp
              if( pelty == SHELL ) then
                 call bouder(&
                      pnode,ndime,ndimb,elmar(pelty) % dercg,&
                      elcod,baloc,detjm)
                 volum = volum + elmar(pelty) % weicg * detjm
              else if( pelty == BAR3D ) then
                 !detjm = 1.0_rp
                 !call bouder(&
                 !     pnode,ndime,1_ip,elmar(pelty) % dercg,&
                 !     elcod,baloc,detjm)
                 volum = volum + 0.0_rp
              else
                 do igaus = 1,pgaus
                    call jacobi(&
                         ndime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                         xjacm,xjaci,cartd,detjm)
                    volum = volum + elmar(pelty)%weigp(igaus)*detjm
                 end do
                 !call jacobi(&
                 !     ndime,pnode,elcod,elmar(pelty) % dercg,&
                 !     xjacm,xjaci,cartd,detjm)
              end if

              if( detjm <= 0.0_rp .and. lelch(ielem) /= ELINT ) then
                 neror(7) = neror(7)+1
                 mess1    = intost(1_ip)
                 mess2    = intost(ielem)
                 call outfor(-2_ip,lun_livei,&
                      'JACOBIAN AT GAUSS POINT '//trim(mess1)&
                      //' OF ELEMENT '//trim(mess2)//' IS ZERO OR NEGATIVE')
              end if

              vomat(pmate) = vomat(pmate) + volum   ! Volume per material        
              vodom        = vodom + volum          ! Total volume
              if( volum >= vomax ) then
                 vomax = volum                      ! Maximum volume
                 elmax = ielem
              end if
              if( volum <= vomin) then
                 vomin = volum                      ! Minimum volume
                 elmin = ielem
              end if
           end if
        end do
     end if

     if( INOTMASTER ) nelet = nelem
     call PAR_SUM(neror(7))
     call PAR_SUM(nelet)
     call PAR_SUM(vodom)
     call PAR_SUM(vomat)
     call PAR_MIN(vomin,'IN MY CODE','EXCLUDE MASTER',kfl_paral_min)
     call PAR_MAX(vomax,'IN MY CODE','EXCLUDE MASTER',kfl_paral_max)
     voave = vodom/real(nelet,rp)               ! Averaged volume

     if( INOTMASTER ) then
        if( kfl_paral_min /= kfl_paral ) then
           elmin = 0
        else
           elmin = leinv_loc(elmin)
        end if
        if( kfl_paral_max /= kfl_paral ) then
           elmax = 0
        else
           elmax = leinv_loc(elmax)
        end if
     end if
     call PAR_MAX(elmin)
     call PAR_MAX(elmax)

     if( neror(7) /= 0 ) call runend('MESCEK: ELEMENT(S) WITH NEGATIVE JACOBIAN')

     call messages_live('CHECK BOUNDARY ORDERING')
     
     !-------------------------------------------------------------------
     !
     ! Boundary Jacobian
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then
        !
        ! Check boundary ordering by computing Jacobian sign
        !
        surfa = 0.0_rp
        do iboun = 1,nboun
           pblty = ltypb(iboun)
           ielem = lelbo(iboun)
           if( pblty > 0 ) then
              pnodb = nnode(pblty)
              pgaub = ngaus(pblty)
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 do idime = 1,ndime
                    elcod(idime,inodb) = coord(idime,ipoin)
                 end do
              end do
              do igaub = 1,pgaub
                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
                      elcod,baloc,detjm)                                   ! and Jacobian
                 !call bouder(&
                 !     pnodb,ndime,ndimb,elmar(pblty)%dercg,&               ! Cartesian derivative
                 !     elcod,baloc,detjm)                                   ! and Jacobian
                 surfa = surfa + elmar(pblty)%weigp(igaub) * detjm

                 if( detjm <= 0.0_rp .and. lelch(ielem) /= ELINT ) then
                    neror(7) = neror(7) + 1
                    mess1    = intost(1_ip)
                    mess2    = intost(iboun)
                    call outfor(-2_ip,lun_livei,&
                         'JACOBIAN AT GAUSS POINT '//trim(mess1)&
                         //' OF BOUNDARY '//trim(mess2)//' IS ZERO OR NEGATIVE')
                 end if
              end do
           end if
        end do

     end if

     call PAR_SUM(neror(7))
     if( neror(7) /= 0 ) call runend('MESCEK: BOUNDARY(IES) WITH NEGATIVE JACOBIAN')

     if( INOTMASTER ) then

        if( kfl_chege == 1 ) then
           !
           ! Checks against two identical nonzero nodal coordinates
           !
           npoit=npoin-1
           npoit=0   ! >>> Not activated. Too much work
           do ipoin=1,npoit
              do jpoin=ipoin+1,npoin
                 keror=0
                 do idime=1,ndime
                    if(coord(idime,ipoin)==coord(idime,jpoin)) keror=keror+1
                 end do
                 if (keror==ndime) then
                    neror(1)=neror(1)+1
                    mess1=intost(ipoin)
                    mess2=intost(jpoin)
                    call outfor(1_ip,lun_outpu,&
                         ' IDENTICAL COORDINATES HAVE BEEN FOUND FOR POINTS NUMBER '&
                         //trim(mess1)//','//trim(mess2))
                 end if
              end do
           end do
           !
           ! Check if all the nodes belong to an element
           !
           call messages_live('CHECK IF ALL NODES BELONG TO AN ELEMENT')
           call memory_alloca(memor_dom,'TOUCH','mescek' ,touch,npoin)
           touch=.false.
           do ielem=1,nelem
              do inode = 1,nnode(ltype(ielem))
                 ipoin=lnods(inode,ielem)
                 touch(ipoin)=.true.
              end do
           end do
           do ipoin=1,npoin
              if(.not.touch(ipoin)) then
                 mess1=intost(ipoin)
                 call outfor(1_ip,lun_outpu,'NODE '//trim(mess1)&
                      //' DOES NOT BELONG TO ANY ELEMENT')
                 neror(6)=neror(6)+1
              end if
           end do
           call memory_deallo(memor_dom,'TOUCH','mescek' ,touch)
           !
           ! Checks for any repetition of a node number within an element
           !
           npoit=npoin
           npoit=0   ! >>> Not activated. Too much work
           do ipoin=1,npoit
              !
              ! Seek first,last and intermediate appearances of node ipoin
              ! & calculate increase or decrease in frontwidth at each element stage
              !
              kstar=0
              do ielem=1,nelem
                 kzero=0
                 do inode=1,lnnod(ielem)
                    if(lnods(inode,ielem)==ipoin) then
                       kzero=kzero+1
                       if(kzero>1) then
                          neror(3)=neror(3)+1
                          mess1=intost(ipoin)
                          mess2=intost(ielem)
                          call outfor(1_ip,lun_outpu,&
                               'NODE '//trim(mess1)&
                               //' APPEARS MORE THAN ONCE IN THE LIST OF'&
                               // ' NODAL CONNECTIONS OF ELEMENT NUMBER '//trim(mess2))
                       end if
                       if(kstar==0) kstar=ielem
                    end if
                 end do
              end do
              if(kstar==0) then
                 !
                 ! Checks if this is an unused node & if it has non-zero coordinates
                 !
                 mess1=intost(ipoin)
                 call outfor(1_ip,lun_outpu,&
                      'NODE '//trim(mess1)//' NEVER APPEARS IN THE ELEMENT CONNECTIVITY')
                 neror(4)=neror(4)+1
              end if
           end do
        end if
     end if

     !-------------------------------------------------------------------
     !
     ! Compute and output element quality
     !
     !-------------------------------------------------------------------

     if( kfl_quali /= 0 ) then

        if( INOTMASTER ) then
           call memgen(0_ip,nelem,0_ip)
           call qualit(gesca,qmaxi,qmini)
        else
           call qualit(dummr,qmaxi,qmini)
        end if
        call openfi(13_ip)

        call memgen(1_ip,100_ip,0_ip)
        dummr = 99.0_rp/(qmaxi-qmini)

        if( INOTMASTER) then
           do ielem = 1,nelem
              if( lelch(ielem) /= ELHOL ) then
                 iqual = int((gesca(ielem)-qmini)*dummr,ip) + 1
                 gisca(iqual) = gisca(iqual) + 1
              end if
           end do
           kelem = 0
           ielem = 0
           do while( ielem < nelem )
              ielem = ielem + 1
              if( lelch(ielem) /= ELHOL ) then
                 if( gesca(ielem) == qmaxi ) then
                    kelem = ielem
                    ielem = nelem
                 end if
              end if
           end do
           call memgen(2_ip,nelem,0_ip)
        end if

        !
        ! The commented part is not scalable: do with allgtaher
        !

        !        if( IMASTER ) then
        !
        !           do kfl_desti_par = 1,npart
        !              call parari('RCV',0_ip,1_ip,ielem)
        !              if( ielem /= 0 ) kelem = kfl_desti_par
        !           end do
        !           do kfl_desti_par = 1,npart
        !              if( kfl_desti_par == kelem ) then
        !                 ielem = kfl_desti_par
        !              else
        !                 ielem = 0
        !              end if
        !              call parari('SND',0_ip,1_ip,ielem)
        !           end do
        !           kfl_desti_par = ielem
        !           call parari('RCV',0_ip,       1_ip,pnode)
        !           call pararr('RCV',0_ip,ndime*pnode,xcoor)
        !
        !        else
        !
        !           kfl_desti_par = 0
        !           call parari('SND',0_ip,1_ip,kelem)
        !           kfl_desti_par = 0
        !           call parari('RCV',0_ip,1_ip,ielem)
        !           if( ielem /= 0 ) then
        !              xcoor = 0.0_rp
        !              pelty = ltype(kelem)
        !              pnode = nnode(pelty)
        !              do inode = 1,pnode
        !                 ipoin = lnods(inode,kelem)
        !                 do idime = 1,ndime
        !                    xcoor(idime,inode) = coord(idime,ipoin)
        !                 end do
        !              end do
        !              call parari('SND',0_ip,       1_ip,pnode)
        !              call pararr('SND',0_ip,ndime*pnode,xcoor)
        !           end if
        !        end if

        call parari('SUM',0_ip,100_ip,gisca)
        if( INOTSLAVE ) then
           write(lun_quali,*) '# Min. quality=',qmini
           write(lun_quali,*) '# Max. quality=',qmaxi
           !write(lun_quali,*) '# Worse element:'
           !do inode = 1,pnode
           !   write(lun_quali,*) '# Node ',inode,'= ',xcoor(1:ndime,inode)
           !end do
           dummr = 1.0_rp / dummr
           do iqual = 1,100
              quali = real(iqual-1,rp) * dummr+qmini
              write(lun_quali,*) quali,gisca(iqual)
           end do
           call iofile_flush_unit(lun_quali)
        end if
        call memgen(3_ip,100_ip,0_ip)
        call openfi(-13_ip)

        !
        ! Compute min, max and average volume
        !
        !print*,'AQUI 1=',kfl_paral
        !call Parall(13_ip)
        !print*,'AQUI 2=',kfl_paral
     end if

  case(2)
     !
     ! Check element types
     !
     call messages_live('CHECK ELEMENT  TYPES')
     neror = 0
     do ielem = 1,nelem
        if( lexis(ltype(ielem)) == 0 ) then
           mess1 = intost(ielem)
           call runend('ERROR WHEN READING ELEMENT '//trim(mess1)//'. CHECK FORMAT')
        else if( ltype(ielem) < iesta_dom .or. ltype(ielem) > iesto_dom ) then
           neror(1) = neror(1)+1
           mess1    = intost(ltype(ielem))
           mess2    = intost(ielem)
           call runend('IMPOSSIBLE TYPE '//trim(mess1)//' OF ELEMENT '//trim(mess2))
        end if
     end do

  case(3)
     !
     ! Checks for impossible node numbers
     !
     call messages_live('CHECK ELEMENT  CONNECTIVITY')    
     neror=0
     do ielem=1,min(nelem,size(ltype,KIND=ip))
        pelty=ltype(ielem)
        pnode=nnode(pelty)
        do inode=1,nnode(pelty)
           if((lnods(inode,ielem)<0).or.(lnods(inode,ielem)>npoin)) then
              neror(1)=neror(1)+1
              mess1=intost(lnods(inode,ielem))
              mess2=intost(ielem)
              call runend('IMPOSSIBLE NODE NUMBER ('//trim(mess1)&
                   //') IN CONNECTIVITY OF ELEMENT '//trim(mess2))
           end if
        end do
     end do
     !
     ! Check nodes not belonging to any element
     !
     call memory_alloca(memor_dom,'TOUCH','mescek' ,touch,npoin)
     ieror = 0
     kpoin = 0
     do ielem = 1,min(nelem,size(ltype,KIND=ip))
        pelty = abs(ltype(ielem))
        do inode = 1,nnode(pelty)
           ipoin = lnods(inode,ielem)
           touch(ipoin) = .true.
        end do
     end do
     do ipoin = 1,npoin
        if( .not. touch(ipoin) ) then
           ieror = ieror + 1
           kpoin = ipoin
        end if
     end do
     if( ieror > 0 ) then
        call runend('NUMBER OF NODES, THAT DO NOT BELONG TO ANY ELEMENTS= '//trim(intost(ieror))//'. LAST NODE FOUND IS '//trim(intost(lninv_loc(kpoin))))
     end if
     call memory_deallo(memor_dom,'TOUCH','mescek' ,touch)
     
  case(4)
     !
     ! Checks for boundary types
     !
     call messages_live('CHECK BOUNDARY TYPES')    
     neror=0
     do iboun=1,nboun
        if(lexis(ltypb(iboun))==0) then
           mess1=intost(iboun)
           call runend('ERROR WHEN READING BOUNDARY '//trim(mess1)//'. CHECK FORMAT')
        else if(ltypb(iboun)<ibsta_dom.or.ltypb(iboun)>ibsto_dom) then
           neror(1)=neror(1)+1 
           mess1=intost(ltypb(iboun))
           mess2=intost(iboun)
           call runend('IMPOSSIBLE TYPE '//trim(mess1)//'OF BOUNDARY '//trim(mess2))
        end if
     end do

  case(5)
     !
     ! Boundary/element connectivity
     !
     call messages_live('CHECK CONNECTIVITY BOUNDARY/ELEMENT')    
     neror=0
     do iboun=1,nboun
        iblty=ltypb(iboun)
        ielem=lelbo(iboun)
        if(ielem<1.or.ielem>nelem) then
           neror(1)=neror(1)+1
           mess1=intost(iboun)
           mess2=intost(ielem)
           call runend('IMPOSSIBLE ELEMENT '//trim(mess2)//' CONNECTIVITY OF BOUNDARY '//trim(mess1))
        end if
     end do

  case(6)
     !
     ! Checks for boundary connectivity
     !    
     call messages_live('CHECK BOUNDARY CONNECTIVITY')    
     do iboun=1,nboun
        pblty=ltypb(iboun)
        pnodb=nnode(pblty)
        do inodb=1,nnode(pblty)
           if((lnodb(inodb,iboun)<0).or.(lnodb(inodb,iboun)>npoin)) then
              neror(2)=neror(2)+1
              mess1=intost(lnodb(inodb,iboun))
              mess2=intost(iboun)
              call runend('IMPOSSIBLE NODE NUMBER ('//trim(mess1)&
                   //') IN CONNECTIVITY OF BOUNDARY '//trim(mess2))
           end if
        end do
     end do

  case(7)
     !
     ! Element/Boundary local numbering
     !
     do iboun=1,nboun
        do inodb=1,nnode(ltypb(iboun))
           if(lboel(inodb,iboun)==0) then
              neror(1)=neror(1)+1
              ielem=lelbo(iboun)
              mess1=intost(iboun)
              mess2=intost(ielem)
              call runend('LOCAL NUMBERING NOT FOUND FOR BOUNDARY '//trim(mess1)//' IN ELEMENT '//trim(mess2))
           end if
        end do
     end do

  end select
  !
  ! Verifies if any mess errors have been detected
  !
  keror=0
  do ieror=1,7
     keror=keror+neror(ieror)
  end do
  if(keror/=0) then
     mess1=intost(keror)
     if(keror==1) then
        call runend('1 ERROR IN MESH DEFINITION HAS BEEN DETECTED')
     else
        call runend(trim(mess1)//' ERRORS IN MESH DEFINITION HAVE BEEN DETECTED')
     end if
  end if

end subroutine mescek
