subroutine witnes()
  !------------------------------------------------------------------------
  !****f* domain/witnes
  ! NAME 
  !    witnes
  ! DESCRIPTION
  !    Get witness point information
  !    LEWIT ... Host element
  !    SHWIT ... Shape function in host element
  ! USES
  ! USED BY
  !    witnes
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_communications, only : PAR_ALL_TO_ALL_ARRAY_OPERATION
  use mod_elsest,         only : elsest_host_element
  use mod_messages,       only : messages_live
  use mod_outfor,         only : outfor
  use mod_messages,       only : livinf
  implicit none
  integer(ip)          :: iwitn,kwitn,idime,inode,ierro
  integer(ip)          :: ielem,ipoin,pnode
  real(rp)             :: coloc(3),dNds(ndime,mnode),elcod(ndime,mnode)
  real(rp)             :: xjaci(9),xjacm(9),gpdet,dista
  integer(ip), pointer :: lnwit_tmp(:)

  if( nwitn > 0 ) then

     nullify( lnwit_tmp )
     call messages_live('COMPUTE WITNESS POINT INFORMATION')
     !
     ! Find host element: LEWIT
     !
     if( INOTMASTER ) then

        call memose(6_ip) 

        do iwitn = 1,nwitn
           call elsest_host_element(&
                ielse,relse,1_ip,meshe(ndivi),cowit(:,iwitn),lewit(iwitn),&
                shwit(:,iwitn),dNds,coloc,dista)
           if( lewit(iwitn) < 1 ) then
              lewit(iwitn) = 0 
           else
              ielem = lewit(iwitn)
              pnode = lnnod(ielem)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 elcod(1:ndime,inode) = coord(1:ndime,ipoin)
              end do
              call jacobi(&
                   ndime,pnode,elcod,dNds,&
                   xjacm,xjaci,dewit(1,1,iwitn),gpdet) 
           end if

           !if( lewit(iwitn) /= 0 ) &
           !     print*,dot_product( shwit(1:pnode,iwitn),coord(1,lnods(1:pnode,ielem))),&
           !     &      dot_product( shwit(1:pnode,iwitn),coord(2,lnods(1:pnode,ielem))),&
           !     &      dot_product( shwit(1:pnode,iwitn),coord(3,lnods(1:pnode,ielem)))

        end do

     else

        allocate( lewit(nwitn) )

     end if
     !
     ! LNWIT: Allocate memory
     !
     call memose(11_ip)
     !
     ! NWITN: Count my own witness points
     !    
     call PAR_ALL_TO_ALL_ARRAY_OPERATION(nwitn,lewit,'CHOOSE ONLY ONE')
     if( INOTMASTER ) then  
        kwitn = 0
        do iwitn = 1,nwitn
           if( lewit(iwitn) /= 0 ) kwitn = kwitn + 1
        end do
     end if
     !
     ! Renumber my witness nodes
     !
     if( ISLAVE ) then
        kwitn = 0
        do iwitn = 1,nwitn
           if( lewit(iwitn) /= 0 ) then
              kwitn        = kwitn + 1
              lewit(kwitn) = lewit(iwitn)
              lnwit(kwitn) = iwitn
              do inode = 1,mnode
                 shwit(inode,kwitn) = shwit(inode,iwitn)
                 do idime = 1,ndime
                    dewit(idime,inode,kwitn) = dewit(idime,inode,iwitn)
                 end do
              end do
              do idime = 1,ndime
                 cowit(idime,kwitn) = cowit(idime,iwitn)
              end do
           end if
        end do
        nwitn = kwitn
     end if
     !
     ! Prepare all gather to replace call Parall(45_ip)
     !
     !if( IPARALL ) then
     !   allocate(nwitn_par(0:npart))
     !   if( IMASTER ) then
     !      nwitn_par = 0
     !      do iwitn = 1,nwitn
     !         ipart = lewit(iwitn)
     !         if( ipart /= 0 ) then
     !            nwitn_par(ipart) = nwitn_par(ipart) + 1
     !         end if
     !      end do
     !   else
     !      nwitn_par = 0
     !      allocate( lnwit_tmp(nwitn) )
     !      lnwit_tmp(1:nwitn) = lnwit(1:nwitn)
     !   end if
     !   call PAR_GATHERV(lnwit_tmp,lnwit,nwitn_par,'IN MY CODE')
     !   if( ISLAVE ) then
     !      deallocate( lnwit_tmp )
     !      deallocate( nwitn_par )
     !   end if
     !end if
     !
     ! Master gets numbering LNWIT
     !
     call Parall(45_ip)
     !
     ! Detect a possible problem
     !
     ierro = 0
     if( IMASTER ) then
        if( mwitn /= nwitn ) ierro = mwitn - nwitn
        if (ierro < 0) ierro= -ierro
     else if( ISEQUEN ) then
        if( kwitn /= nwitn ) ierro = kwitn - nwitn
        if (ierro < 0) ierro= -ierro
     end if
     call outfor(50_ip,lun_outpu,'')
     ioutp(2) = 0
     if( ierro > 0 ) then
        call livinf(0_ip,&
             'WARNING: '//trim(intost(ierro))//' WITNESS POINTS ARE LOST: ZERO RESULTS WILL APPEAR',&
             0_ip)        
        if( IMASTER ) then
           call memgen(1_ip,nwitn,0_ip)
           do iwitn = 1,nwitn
              if( lnwit(iwitn) /= 0 ) gisca(lnwit(iwitn)) = 1
           end do
           do iwitn = 1,nwitn
              if( gisca(iwitn) == 0 )then
                 ioutp(1) = iwitn
                 call outfor(51_ip,lun_outpu,'')
                 ioutp(2) = ioutp(2) + 1
              end if
           end do
           call memgen(3_ip,nwitn,0_ip)
        else if( ISEQUEN ) then        
           do iwitn = 1,nwitn
              if( lewit(iwitn) == 0 )then
                 ioutp(1) = iwitn
                 call outfor(51_ip,lun_outpu,'')
                 ioutp(2) = ioutp(2) + 1
              end if
           end do
        end if
     else
        call outfor(52_ip,lun_outpu,'')
     end if

  end if

end subroutine witnes
