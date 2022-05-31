subroutine divtet()  
  !-----------------------------------------------------------------------
  !****f* Domain/divtet
  ! NAME
  !    divtet
  ! DESCRIPTION
  !    Divide a mesh into tetrahedra
  ! OUTPUT
  !    LNODS, LTYPE, LNODB, LTYPB, LBOEL, LESET, LBSET, KFL_CODBO
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_domain
  use mod_memchk
  implicit none
  integer(ip)          :: ielem,inode,kelem,pelty,knode,pnode
  integer(ip)          :: iboun,inodb,kboun,pblty,knodb,pnodb
  integer(ip)          :: ipoin,jpoin,ifoun,lboun
  integer(ip)          :: qelem,tetra(4,6),ldia1(2),ldia2(2)
  integer(ip)          :: nelem_tmp,nboun_tmp
  integer(ip), pointer :: lnume(:) => null()
  integer(ip), pointer :: lnumb(:) => null()
  integer(ip), pointer :: lnumk(:) => null()
  integer(ip), pointer :: ltype_tmp(:) => null()
  integer(ip), pointer :: lnods_tmp(:,:) => null()
  integer(ip), pointer :: ltypb_tmp(:) => null()
  integer(ip), pointer :: lnodb_tmp(:,:) => null()
  integer(ip), pointer :: lboel_tmp(:,:) => null()
  integer(ip), pointer :: lelbo_tmp(:) => null()
  integer(ip), pointer :: leset_tmp(:) => null()
  integer(ip), pointer :: lbset_tmp(:) => null()
  integer(ip), pointer :: kfl_codbo_tmp(:) => null()
  integer(4)           :: istat

  if( kfl_divid == 0 ) return 
  !
  ! Calculate only the boundary/element nodal connectivity
  !
  if( kfl_bouel == 1 ) then
     do iboun = 1,nboun
        knodb = nnode(ltypb(iboun))
        ielem = lelbo(iboun)
        knode = nnode(ltype(ielem))
        do inodb = 1,knodb
           ipoin = lnodb(inodb,iboun)
           nodes2: do inode = 1,knode
              jpoin = lnods(inode,ielem)
              if( ipoin == jpoin ) then
                 lboel(inodb,iboun) = inode
                 exit nodes2
              end if
           end do nodes2
        end do
     end do
  else
     call runend('DIVTET: ELEMENT/BOUNDARY CONNECTIVITY MUST BE GIVEN')
  end if
  !
  ! Renumbering
  !
  allocate(lnume(nelem),stat=istat)
  call memchk(zero,istat,memor_dom,'LNUME','divtet',lnume)
  allocate(lnumb(nboun),stat=istat)
  call memchk(zero,istat,memor_dom,'LNUMB','divtet',lnumb)
  !
  ! Save old values
  !  
  nelem_tmp = nelem
  nboun_tmp = nboun
  
  !----------------------------------------------------------------------
  !
  ! Element: LNODS, LTYPE
  !
  !----------------------------------------------------------------------

  allocate(ltype_tmp(nelem),stat=istat)
  call memchk(zero,istat,memor_dom,'LTYPE_TMP','divtet',ltype_tmp)
  allocate(lnods_tmp(mnode,nelem),stat=istat)
  call memchk(zero,istat,memor_dom,'LNODS_TMP','divtet',lnods_tmp)

  kelem = 0
  do ielem = 1,nelem
     pelty = ltype(ielem)
     ltype_tmp(ielem) = pelty
     if( pelty == HEX08 ) then
        kelem = kelem + 6
     else
        kelem = kelem + 1
     end if
     do inode = 1,nnode(pelty)
        lnods_tmp(inode,ielem) = lnods(inode,ielem)       
     end do
  end do
  !
  ! Reallocate
  !
  deallocate(ltype,stat=istat)
  if(istat/=0) call memerr(two,'LTYPE','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LTYPE','divtet',ltype)  
  deallocate(lnods,stat=istat)
  if(istat/=0) call memerr(two,'LNODS','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LNODS','divtet',lnods)

  nelem = kelem
  mnode = nnode(TET04)
  lexis(HEX08) = 0

  allocate(ltype(nelem),stat=istat)
  call memchk(zero,istat,memor_dom,'LTYPE','divtet',ltype)
  allocate(lnods(mnode,nelem),stat=istat)
  call memchk(zero,istat,memor_dom,'LNODS','divtet',lnods)
  allocate(lnumk(nelem),stat=istat)
  call memchk(zero,istat,memor_dom,'LNUME','divtet',lnume)

  kelem = 0
  do ielem = 1,nelem_tmp

     pelty = ltype_tmp(ielem)
     pnode = nnode(pelty)

     if( pelty == HEX08 ) then

        tetra(1,1) = lnods_tmp(1,ielem)
        tetra(2,1) = lnods_tmp(3,ielem)
        tetra(3,1) = lnods_tmp(4,ielem)
        tetra(4,1) = lnods_tmp(5,ielem)

        tetra(1,2) = lnods_tmp(3,ielem)
        tetra(2,2) = lnods_tmp(4,ielem)
        tetra(3,2) = lnods_tmp(5,ielem)
        tetra(4,2) = lnods_tmp(8,ielem)

        tetra(1,3) = lnods_tmp(1,ielem)
        tetra(2,3) = lnods_tmp(2,ielem)
        tetra(3,3) = lnods_tmp(3,ielem)
        tetra(4,3) = lnods_tmp(5,ielem)
                              
        tetra(1,4) = lnods_tmp(2,ielem)
        tetra(2,4) = lnods_tmp(3,ielem)
        tetra(3,4) = lnods_tmp(5,ielem)
        tetra(4,4) = lnods_tmp(6,ielem)

        tetra(1,5) = lnods_tmp(3,ielem)
        tetra(2,5) = lnods_tmp(5,ielem)
        tetra(3,5) = lnods_tmp(6,ielem)
        tetra(4,5) = lnods_tmp(7,ielem)
                              
        tetra(1,6) = lnods_tmp(3,ielem)
        tetra(2,6) = lnods_tmp(5,ielem)
        tetra(3,6) = lnods_tmp(7,ielem)
        tetra(4,6) = lnods_tmp(8,ielem)

        lnume(ielem) = -(kelem+1)
        
        do qelem = 1,6
           kelem        = kelem  + 1
           lnumk(kelem) = ielem
           ltype(kelem) = TET04
           do knode = 1,4
              lnods(knode,kelem) = tetra(knode,qelem)
           end do
        end do

     else if( pelty == TET04 ) then

        kelem        = kelem + 1
        lnumk(kelem) = ielem
        lnume(ielem) = kelem
        ltype(kelem) = ltype_tmp(ielem)
        do inode = 1,pnode
           lnods(inode,kelem) = lnods_tmp(inode,ielem)
        end do

     else

        call runend('DIVTET: CAN ONLY DIVIDE HEX08')

     end if

  end do

  deallocate(lnods_tmp,stat=istat)
  if(istat/=0) call memerr(two,'LNODS_TMP','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LNODS_TMP','divtet',lnods_tmp)
  
  !----------------------------------------------------------------------
  !
  ! Boundary: LNODB, LTYPB, LBOEL
  !
  !----------------------------------------------------------------------

  allocate(ltypb_tmp(nboun),stat=istat)
  call memchk(zero,istat,memor_dom,'LTYPB_TMP','divtet',ltypb_tmp)
  allocate(lnodb_tmp(mnodb,nboun),stat=istat)
  call memchk(zero,istat,memor_dom,'LNODB_TMP','divtet',lnodb_tmp)
  allocate(lboel_tmp(mnodb,nboun),stat=istat)
  call memchk(zero,istat,memor_dom,'LBOEL_TMP','divtet',lboel_tmp)
  allocate(lelbo_tmp(nboun),stat=istat)
  call memchk(zero,istat,memor_dom,'LBOEL_TMP','divtet',lelbo_tmp)

  kboun = 0
  do iboun = 1,nboun
     pblty = ltypb(iboun)
     pnodb = nnode(pblty)
     ltypb_tmp(iboun) = pblty
     if( pblty == QUA04 ) then
        kboun = kboun + 2
     else if( pblty == TRI03 ) then
        kboun = kboun + 1
     else
        call runend('DIVTET: CAN ONLY DIVIDE QUA04')
     end if
     do inodb = 1,pnodb
        lnodb_tmp(inodb,iboun) = lnodb(inodb,iboun)
        lboel_tmp(inodb,iboun) = lboel(inodb,iboun)
     end do
     lelbo_tmp(iboun) = lelbo(iboun)
  end do
  !
  ! Reallocate
  !
  deallocate(ltypb,stat=istat)
  if(istat/=0) call memerr(two,'LTYPB','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LTYPB','divtet',ltypb)  
  deallocate(lnodb,stat=istat)
  if(istat/=0) call memerr(two,'LNODB','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LNODB','divtet',lnodb)
  deallocate(lboel,stat=istat)
  if(istat/=0) call memerr(two,'LBOEL','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LBOEL','divtet',lboel)
  deallocate(lelbo,stat=istat)
  if(istat/=0) call memerr(two,'LELBO','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LELBO','divtet',lelbo)

  nboun = kboun
  mnodb = nnode(TRI03)
  lexis(QUA04) = 0

  allocate(ltypb(nboun),stat=istat)
  call memchk(zero,istat,memor_dom,'LTYPB','divtet',ltypb)
  allocate(lnodb(mnodb,nboun),stat=istat)
  call memchk(zero,istat,memor_dom,'LNODB','divtet',lnodb)
  allocate(lboel(mnodb,nboun),stat=istat)
  call memchk(zero,istat,memor_dom,'LBOEL','divtet',lboel)
  allocate(lelbo(nboun),stat=istat)
  call memchk(zero,istat,memor_dom,'LELBO','divtet',lelbo)

  kboun = 0

  do iboun = 1,nboun_tmp

     pblty = ltypb_tmp(iboun)
     
     if( pblty == QUA04 ) then

        lnumb(iboun) = -QUA04
        ielem        = lelbo_tmp(iboun)
        kelem        = lnume(ielem)
        if( kelem > 0 ) call runend('DIVTET: QUA04 DOES NOT BELONG TO SPLIT ELEMENT')

        ldia1(1)   = lboel_tmp(1,iboun)
        ldia1(2)   = lboel_tmp(3,iboun)
        if( ldia1(2) < ldia1(1) ) then
           ldia1(2)   = lboel_tmp(1,iboun)
           ldia1(1)   = lboel_tmp(3,iboun)
        end if
        ldia2(1)   = lboel_tmp(2,iboun)
        ldia2(2)   = lboel_tmp(4,iboun)
        if( ldia2(2) < ldia2(1) ) then
           ldia2(1)   = lboel_tmp(4,iboun)
           ldia2(2)   = lboel_tmp(2,iboun)
        end if
        pnodb = 3

        if(  ( ldia1(1) == 4 .and. ldia1(2) == 5 ) .or. &
             ( ldia1(1) == 1 .and. ldia1(2) == 3 ) .or. &            
             ( ldia1(1) == 3 .and. ldia1(2) == 8 ) .or. &        
             ( ldia1(1) == 2 .and. ldia1(2) == 5 ) .or. &        
             ( ldia1(1) == 3 .and. ldia1(2) == 6 ) .or. &        
             ( ldia1(1) == 5 .and. ldia1(2) == 7 ) ) then

           kboun          = kboun + 1
           ltypb(kboun)   = TRI03
           lnodb(1,kboun) = lnodb_tmp(1,iboun)
           lnodb(2,kboun) = lnodb_tmp(2,iboun)
           lnodb(3,kboun) = lnodb_tmp(3,iboun)

           kboun          = kboun + 1
           ltypb(kboun)   = TRI03
           lnodb(1,kboun) = lnodb_tmp(1,iboun)
           lnodb(2,kboun) = lnodb_tmp(3,iboun)
           lnodb(3,kboun) = lnodb_tmp(4,iboun)

        else

           kboun          = kboun + 1
           ltypb(kboun)   = TRI03
           lnodb(1,kboun) = lnodb_tmp(1,iboun)
           lnodb(2,kboun) = lnodb_tmp(2,iboun)
           lnodb(3,kboun) = lnodb_tmp(4,iboun)

           kboun          = kboun + 1
           ltypb(kboun)   = TRI03
           lnodb(1,kboun) = lnodb_tmp(2,iboun)
           lnodb(2,kboun) = lnodb_tmp(3,iboun)
           lnodb(3,kboun) = lnodb_tmp(4,iboun)

        end if
        !
        ! Fill in lboel for the 2 boundaries generated (kboun-1 and kboun)
        !
        do lboun = kboun-1,kboun

           kelem = abs(lnume(ielem)) - 1
           qelem = 0
           do while( qelem < 6 )
              qelem = qelem + 1
              kelem = kelem + 1           
              knode = nnode(ltype(kelem))
              inodb = 0
              do while( inodb < 3 )
                 inodb = inodb + 1
                 ipoin = lnodb(inodb,lboun)
                 ifoun = 0
                 
                 nodes1: do inode = 1,knode
                    jpoin = lnods(inode,kelem)
                    if( ipoin == jpoin ) then
                       ifoun = 1
                       exit nodes1
                    end if
                 end do nodes1
                 if( ifoun == 0 ) inodb = 4
                 
              end do
              
              if( inodb /= 4 ) then
                 lelbo(lboun) = kelem
                 qelem = 7
              end if
           end do
           if( qelem /= 7 ) call runend('DIVTET: COULD NOT FIND ELEMENT')

        end do

     else

        kboun        = kboun + 1
        lnumb(iboun) = kboun
        ltypb(kboun) = ltypb_tmp(iboun)
        pnodb        = pnodb
        do inodb = 1,nnode(ltypb_tmp(iboun))
           lnodb(inodb,kboun) = lnodb_tmp(inodb,iboun)
           lboel(inodb,kboun) = lboel_tmp(inodb,iboun)
        end do
        lelbo(kboun) = lnume(ielem)

     end if

  end do

  deallocate(lnodb_tmp,stat=istat)
  if(istat/=0) call memerr(two,'LNODB_TMP','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LNODB_TMP','divtet',lnodb_tmp)
  deallocate(lboel_tmp,stat=istat)
  if(istat/=0) call memerr(two,'LBOEL_TMP','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LBOEL_TMP','divtet',lboel_tmp)
  deallocate(lelbo_tmp,stat=istat)
  if(istat/=0) call memerr(two,'LELBO_TMP','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LELBO_TMP','divtet',lelbo_tmp)
  
  !----------------------------------------------------------------------
  !
  ! Sets and boundary conditions: LESET, LBSET, KFL_CODBO
  !
  !----------------------------------------------------------------------

  if( neset > 0 ) then
     !
     ! LESET
     !
     allocate(leset_tmp(nelem_tmp),stat=istat)
     call memchk(zero,istat,memor_dom,'LESET_TMP','divtet',leset_tmp)

     do ielem = 1,nelem_tmp
        leset_tmp(ielem) = leset(ielem)
     end do

     deallocate(leset,stat=istat)
     if(istat/=0) call memerr(two,'LESET','divtet',0_ip)
     call memchk(two,istat,memor_dom,'LESET','divtet',leset)  
     allocate(leset(nelem),stat=istat)
     call memchk(zero,istat,memor_dom,'LESET','divtet',leset)

     kelem = 0
     do ielem = 1,nelem_tmp
        if( ltype_tmp(ielem) == TET04 ) then
           kelem = kelem + 1
           leset(kelem) = leset_tmp(ielem)
        else if ( ltype_tmp(ielem) == HEX08 ) then
           do knode = 1,6
              kelem = kelem + 1
              leset(kelem) = leset_tmp(ielem)
           end do
        end if
     end do

     deallocate(leset_tmp,stat=istat)
     if(istat/=0) call memerr(two,'LESET_TMP','divtet',0_ip)
     call memchk(two,istat,memor_dom,'LESET_TMP','divtet',leset_tmp)  
     
  end if

  if( nbset > 0 ) then
     !
     ! NBSET
     !
     allocate(lbset_tmp(nboun_tmp),stat=istat)
     call memchk(zero,istat,memor_dom,'LBSET_TMP','divtet',lbset_tmp)
     do iboun = 1,nboun_tmp
        lbset_tmp(iboun) = lbset(iboun)
     end do

     deallocate(lbset,stat=istat)
     if(istat/=0) call memerr(two,'LBSET','divtet',0_ip)
     call memchk(two,istat,memor_dom,'LBSET','divtet',lbset)
     allocate(lbset(nboun),stat=istat)
     call memchk(zero,istat,memor_dom,'LBSET','divtet',lbset)

     kboun = 0
     do iboun = 1,nboun_tmp
        if( ltypb_tmp(iboun) == TRI03 ) then
           kboun = kboun + 1
           lbset(kboun) = lbset_tmp(iboun)
        else if ( ltypb_tmp(iboun) == QUA04 ) then
           do knode = 1,2
              kboun = kboun + 1
              lbset(kboun) = lbset_tmp(iboun)
           end do
        end if
     end do
 
     deallocate(lbset_tmp,stat=istat)
     if(istat/=0) call memerr(two,'LBSET_TMP','divtet',0_ip)
     call memchk(two,istat,memor_dom,'LBSET_TMP','divtet',lbset_tmp)
    
  end if

  if( kfl_icodb > 0 ) then
     !
     ! KFL_CODBO
     !
     allocate(kfl_codbo_tmp(nboun_tmp),stat=istat)
     call memchk(zero,istat,memor_dom,'KFL_CODBO_TMP','divtet',kfl_codbo_tmp)
     do iboun = 1,nboun_tmp
        kfl_codbo_tmp(iboun) = kfl_codbo(iboun)
     end do

     deallocate(kfl_codbo,stat=istat)
     if(istat/=0) call memerr(two,'KFL_CODBO','divtet',0_ip)
     call memchk(two,istat,memor_dom,'KFL_CODBO','divtet',kfl_codbo)
     allocate(kfl_codbo(nboun),stat=istat)
     call memchk(zero,istat,memor_dom,'KFL_CODBO','divtet',kfl_codbo)

     kboun = 0
     do iboun = 1,nboun_tmp
        if( ltypb_tmp(iboun) == TRI03 ) then
           kboun = kboun + 1
           kfl_codbo(kboun) = kfl_codbo_tmp(iboun)
        else if ( ltypb_tmp(iboun) == QUA04 ) then
           do knode = 1,2
              kboun = kboun + 1
              kfl_codbo(kboun) = kfl_codbo_tmp(iboun)
           end do
        end if
     end do
 
     deallocate(kfl_codbo_tmp,stat=istat)
     if(istat/=0) call memerr(two,'KFL_CODBO_TMP','divtet',0_ip)
     call memchk(two,istat,memor_dom,'KFL_CODBO_TMP','divtet',kfl_codbo_tmp)

  end if

  !----------------------------------------------------------------------
  !
  ! Deallocate
  !
  !----------------------------------------------------------------------

  deallocate(ltypb_tmp,stat=istat)
  if(istat/=0) call memerr(two,'LTYPB_TMP','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LTYPB_TMP','divtet',ltypb_tmp)  

  deallocate(ltype_tmp,stat=istat)
  if(istat/=0) call memerr(two,'LTYPE_TMP','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LTYPE_TMP','divtet',ltype_tmp)

  deallocate(lnumb,stat=istat)
  if(istat/=0) call memerr(two,'LNUMB','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LNUMB','divtet',lnumb)  

  deallocate(lnume,stat=istat)
  if(istat/=0) call memerr(two,'LNUME','divtet',0_ip)
  call memchk(two,istat,memor_dom,'LNUME','divtet',lnume)  

end subroutine divtet
