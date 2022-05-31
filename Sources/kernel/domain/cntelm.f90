!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    cntelm.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Create contact related stuff
!> @details Create contact related stuff
!>    NECNT ............ Number of contact elements
!>
!>    Contact faces LFCNT
!>    LFCNT(1,IECNT) ... 1: Fluid face of contact element
!>    LFCNT(2,IECNT) ... 2: Solid face of contact element
!>    LFCNT(3,IECNT) ... IELEM: Corresponding element
!>    LFCNT(4,IECNE) ... FFACE: fluid face
!>    LFCNT(5,IECNE) ... FELEM: fluid element
!>
!>    LNCNT(1,INCNT) ... IPOIN: Fluid/ALE node contrained to JPOIN
!>    LNCNT(2,INCNT) ... JPOIN: Solid node
!>    LNOCH(  IPOIN) ... When ipoin is a contact: NODE_CONTACT_FLUID or NODE_CONTACT_SOLID
!>
!>        Solid     Contact     Fluid
!>    +------------+-------+-------------+
!>    |            | IECNT |   FELEM     |
!>    |            | IELEM |             |
!>    |            | 2   1 | FFACE       |
!>    |            |<=   =>|<=           |
!>    +------------+-------+-------------+ 
!>
!>       LNCNT(2,INCNT)  LNCNT(1,INCNT)
!>    -1          -1       0             1
!>    +------------X-------X-------------+
!>    |            |       |             |
!>    |            |       |             |
!>    |            |       |             |
!>    |            |       |             |
!>    +------------+-------+-------------+ 
!>    -1          -1       0             1
!>
!>    Note:
!>    LMATN(IPOIN) ... 1=fluid, -1=solid, 0=interface 
!>
!> @} 
subroutine cntelm()
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use mod_memory
  use mod_outfor, only : outfor
  use mod_messages, only : livinf
  implicit none
  integer(ip)          :: ielem,iface,pelty,pnodf,ipoin,inode,inodf,ierro
  integer(ip)          :: fnodb,snodb,fface,sface,jnodf,jpoin,jnode,iecnt,incnt
  integer(ip)          :: kfaci(mnodb),kfacj(mnodb),pnodi,pnodj,ii,jj,jelty,jelem
  integer(ip)          :: jface,ielty,ielel,kelem,kface
  integer(ip), pointer :: lzsol(:) => null()
  character(20)        :: mess1
  !
  ! Create list of faces LFCNT(3,NECNT) and nodes LNCNT(2,NECNT)
  !
  if( mecnt > 0 ) then
     !
     ! List of fluid nodes
     !
     !if( kfl_modul(ID_SOLIDZ) == 0 ) then
     !   call runend('SOLIDZ MODULE MUST BE SOLVED')
     !end if

     if( INOTMASTER ) then
        call memgen(1_ip,npoin,0_ip)   ! we need a gisca of size NPOIN
        call memory_alloca(memor_dom,'LZSOL','cntelm',lzsol,npoin)
        do ielem = 1,nelem
           if( lelch(ielem) /= ELCNT ) then
              !do kelem = 1,nelez(lzone(ID_NASTIN))
              !   ielem = lelez(lzone(ID_NASTIN)) % l(kelem)
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 lzsol(ipoin) = 1                               ! lzsol(ipoin)=1 means that ipoin is on the solid zone
              end do
           end if
        end do
     end if

     call livinf(0_ip,'DETECT FACES OF CONTACT ELEMENTS',0_ip)
     ierro = 0

     if( INOTMASTER ) then
        !
        ! Slaves interchange and normalization
        !
        call parari('SLX',NPOIN_TYPE,npoin,lzsol)
        do ipoin= 1,npoin
           if (lzsol(ipoin) > 1) lzsol(ipoin)=1
        end do
        !
        ! LFCNT: Contact faces
        !
        necnt = 0
        do ielem = 1,nelem
           if( lelch(ielem) == ELCNT ) then
              necnt = necnt + 1
              pelty = abs(ltype(ielem))
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 gisca(ipoin) = 1                              ! label in gisca the nodes in the contact elements
              end do
              do iface = 1,nface(pelty)
                 pnodf = nnodf(pelty) % l(iface)
                 fnodb = 0
                 snodb = 0
                 do inodf = 1,pnodf
                    inode = lface(pelty) % l(inodf,iface) 
                    ipoin = lnods(inode,ielem)
                    if( lzsol(ipoin) == 0 ) then
                       fnodb = fnodb + 1                       ! the node is on the fluids zone (it was not counted on lzsol)
                    else
                       snodb = snodb + 1                       ! the node is on the solidz zone
                    end if
                 end do                                        
                 if( fnodb == pnodf ) then
                    lfcnt(1,necnt) = iface                     ! this is the fluid face of the contact element
                    lfcnt(3,necnt) = ielem
                 else if( snodb == pnodf ) then
                    lfcnt(2,necnt) = iface                     ! this is the solid face of the contact element
                    lfcnt(3,necnt) = ielem
                 end if
              end do
           end if
        end do
        do iecnt = 1,necnt                                     ! loop over the contact elements to see if they have all
           ielem = lfcnt(3,iecnt)                              ! a face on the solid zone and the other one on the fluid zone...
           if( lfcnt(1,iecnt) == 0 ) then
              ierro = ierro + 1                                ! ... if not, stop the run
              mess1 = intost(ielem)
              call outfor(-2_ip,lun_livei,&
                   'CNTELM: CANNOT FIND FLUID FACE FOR CONTACT ELEMENT '//trim(mess1))
           else if( lfcnt(2,iecnt) == 0 ) then
              ierro = ierro + 1
              mess1 = intost(ielem)
              call outfor(-2_ip,lun_livei,&
                   'CNTELM: CANNOT FIND SOLID FACE FOR CONTACT ELEMENT '//trim(mess1))
           end if
        end do
     end if
     !
     ! Check errors
     !
     call parari('SUM',0_ip,1_ip,ierro)
     if( ierro > 0 ) call runend('CNTELM: WRONG CONTACT ELEMENTS')
     ierro = 0

     if( INOTMASTER ) then
        !
        ! Look for fluid face: 
        !      LFCNT(4,IECNE) = fluid face 
        !      LFCNT(5,IECNE) = fluid element
        ! The fluid face is the neighboring face to LFCNT(1,IECNE)
        !
        if( necnt > 0 ) then

           do iecnt = 1,necnt
              ielem = lfcnt(3,iecnt)
              iface = lfcnt(1,iecnt) 
              ielty = abs(ltype(ielem))
              pnodi = nnodf(ielty) % l(iface)
              do ii = 1,pnodi
                 kfaci(ii) = lnods(lface(ielty) % l(ii,iface),ielem)
              end do
              do ii = 1,pnodi-1
                 do jj = ii+1,pnodi
                    if( kfaci(ii) > kfaci(jj) ) then
                       kface     = kfaci(ii)
                       kfaci(ii) = kfaci(jj) 
                       kfaci(jj) = kface
                    end if
                 end do
              end do
              ielel = pelel_2(ielem)-1
              do while( ielel < pelel_2(ielem+1)-1 )
                 ielel = ielel + 1
                 jelem = lelel_2(ielel)
                 jelty = abs(ltype(jelem))
                 jface = 1
                 do while( jface <= nface(jelty) )
                    pnodj = nnodf(jelty) % l(jface)
                    if( pnodj == pnodi ) then
                       do ii = 1,pnodj
                          kfacj(ii) = lnods(lface(jelty) % l(ii,jface),jelem)
                       end do
                       do ii = 1,pnodj-1
                          do jj = ii+1,pnodj
                             if( kfacj(ii) > kfacj(jj) ) then
                                kface     = kfacj(ii)
                                kfacj(ii) = kfacj(jj) 
                                kfacj(jj) = kface
                             end if
                          end do
                       end do
                       inode = 0
                       do ii = 1,pnodi
                          if( kfaci(ii) == kfacj(ii) ) inode = inode + 1
                       end do
                       if( inode == pnodi ) then
                          lfcnt(4,iecnt) = jface
                          lfcnt(5,iecnt) = jelem
                          jface          = 100
                          jelem          = pelel_2(ielem+1)
                       end if
                    end if
                    jface = jface + 1
                 end do
              end do
           end do
        end if
        !
        ! LNCNT: Contact nodes
        ! 
        !  They are ordered by pairs: lncnt(1:2,1:nncnt) 
        !
        incnt = 0
        do iecnt = 1,necnt
           ielem = lfcnt(3,iecnt)      ! for element ielem...
           fface = lfcnt(1,iecnt)      ! ... fface is its fluid face
           sface = lfcnt(2,iecnt)       ! ... and sface is its solid face
           do inodf = 1,nnodf(pelty) % l(fface)             ! loop over the nodes of the fluid face
              inode = lface(pelty) % l(inodf,fface) 
              ipoin = lnods(inode,ielem)                    
              if( gisca(ipoin) == 1 ) then                  ! was the node labelled before (so it belongs to any of the contacts)?
                 jnodf = 0
                 do while( jnodf < nnodf(pelty) % l(sface) )
                    jnodf = jnodf + 1                       
                    jnode = lface(pelty) % l(jnodf,sface) 
                    jpoin = lnods(jnode,ielem)
                    if(    abs(coord(    1,ipoin)-coord(    1,jpoin)) < 1.0e-12_rp .and. &
                         & abs(coord(    2,ipoin)-coord(    2,jpoin)) < 1.0e-12_rp .and. &
                         & abs(coord(ndime,ipoin)-coord(ndime,jpoin)) < 1.0e-12_rp ) then
                       ! the correspondace is found
                       incnt = incnt + 1
                       lncnt(1,incnt) = ipoin 
                       lncnt(2,incnt) = jpoin
                       jnodf = 2*nnodf(pelty) % l(sface)
                    end if
                 end do
                 if( jnodf /= 2*nnodf(pelty) % l(sface) ) then
                    ! no correspondence, so not all the contacts where well defined
                    ierro = ierro + 1                             
                 else
                    ! correspondence found, so label both nodes with a 0
                    gisca(ipoin) = 0
                    gisca(jpoin) = 0
                 end if
              end if
           end do
        end do
        !
        ! Error check
        !
        do iecnt = 1,necnt
           if(    lfcnt(1,iecnt) == 0 .or. &
                & lfcnt(4,iecnt) == 0 .or. &
                & lfcnt(5,iecnt) == 0 ) ierro = ierro + 1
        end do
     end if
     !
     ! Check errors
     !
     call parari('SUM',0_ip,1_ip,ierro)
     if( ierro > 0 ) call runend('CNTELM: COULD NOT FIND FACES FOR SOME CONTACT ELEMENTS')
     !
     ! LNOCH: node characteristics, i.e, what zone each node belongs to
     !
     if( INOTMASTER ) then
        do incnt = 1,nncnt
           ipoin = lncnt(1,incnt)
           jpoin = lncnt(2,incnt)
           if( ipoin /= 0 ) gisca(ipoin) =  1                  ! fluid node
           if( jpoin /= 0 ) gisca(jpoin) = -1                  ! solid node
           if( gisca(ipoin) == gisca(jpoin) ) ierro = ierro + 1
        end do
        call parari('SLX',NPOIN_TYPE,npoin,gisca)
        do ipoin = 1,npoin
           if( gisca(ipoin) > 0 ) then
              lnoch(ipoin) =  NODE_CONTACT_FLUID
           else if( gisca(ipoin) < 0 ) then
              lnoch(ipoin) =  NODE_CONTACT_SOLID
           end if
        end do
     end if
     !
     ! LZSOL: Deallocate
     !
     if( INOTMASTER ) then
        call memgen(3_ip,npoin,0_ip)  ! deallocate gisca
        call memory_deallo(memor_dom,'LZSOL','cntelm',lzsol)
     end if

  end if

end subroutine cntelm
