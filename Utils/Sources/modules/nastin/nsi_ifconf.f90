subroutine nsi_ifconf(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_ifconf
  ! NAME 
  !    nsi_ifconf
  ! DESCRIPTION
  !    This routine checks if the flow is confined or not and look for a node
  !    to impose pressure 
  ! USED BY
  !    nsi_solite
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_nastin
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idime,ielem,kbopo,iperi
  integer(ip)             :: inode

  return

  select case(itask)

  case(2_ip)

     if( INOTMASTER .and. kfl_confi_nsi == 1 .and. nodpr_nsi /= 0 .and. nperi > 0 ) then
        do iperi = 1,nperi
           if( lperi(2,iperi) == nodpr_nsi ) then
              nodpr_nsi = lperi(1,iperi)
           end if
        end do

     end if

  end select

  return

  select case(itask)

  case(1)
     !
     ! Look for first node on boundary
     !
     if(kfl_paral==0) then
        if(nodpr_nsi==0.and.kfl_confi_nsi==1) then
           if(nmate<=1) then
              ipoin=0
              do while(ipoin<npoin) 
                 ipoin=ipoin+1
                 if(lpoty(ipoin)/=0) then
                    nodpr_nsi=ipoin
                    ipoin=npoin
                 end if
              end do
           else
              elements2: do ielem=1,nelem
                 if(lmate(ielem)==1) then
                    inode=0
                    do while(inode<nnode(ltype(ielem)))
                       inode=inode+1
                       ipoin=lnods(inode,ielem)
                       if(lpoty(ipoin)/=0) then
                          nodpr_nsi=lnods(inode,ielem)
                          exit elements2
                       end if
                    end do
                 end if
              end do elements2
           end if
        end if
     end if

  case(2)
     !
     ! Check if flow is confined
     !
     if(kfl_paral/=0) then
        if(kfl_confi_nsi/=1.and.kfl_confi_nsi/=-1.and.kfl_regim_nsi/=1) then
           !
           ! Automatic check
           !
           kfl_confi_nsi=0
           do ipoin=1,npoin
              dimensions1: do idime=1,ndime
                 if(kfl_fixno_nsi(idime,ipoin)==1) then
                    kfl_confi_nsi=kfl_confi_nsi+1
                    exit dimensions1
                 end if
              end do dimensions1
           end do
           kbopo=nbopo
           if(kfl_confi_nsi==kbopo) then
              kfl_confi_nsi=1
           else
              kfl_confi_nsi=0
           end if
        end if
        !
        ! If pressure node is not prescribed, find one
        !  
        if(nodpr_nsi==0.and.kfl_confi_nsi==1) then
           ipoin=0
           do while(ipoin<npoin) 
              ipoin=ipoin+1
              if(lpoty(ipoin)/=0) then
                 nodpr_nsi=ipoin
                 ipoin=npoin
              end if
           end do
        end if
        !
        ! Pressure is prescribed on NODPR_NSI
        !
        if(kfl_confi_nsi==1.and.nodpr_nsi/=0.and.nperi>0) then
           do iperi = 1,nperi
              if( lperi(2,iperi) == nodpr_nsi ) then
                 nodpr_nsi = lperi(1,iperi)
              end if
           end do
           !call runend('NSI_IFCONF: TO BE CODED?')
           !if(ipass==0) then
           !   ipass=1
           !   loop_iperi: do iperi=1,nperi
           !      kperi=size(lperi(iperi)%l,KIND=ip)
           !      ipres=0
           !      loop_ii: do ii=1,kperi
           !         ipoin=lperi(iperi)%l(ii)
           !         if(ipoin==nodpr_nsi) then
           !            ipres=1
           !            exit loop_ii
           !         end if
           !      end do loop_ii
           !      if(ipres==1) then
           !         kfl_perip_nsi=1
           !         do ii=1,kperi
           !            lperp_nsi(ii)=lperi(iperi)%l(ii)
           !         end do
           !         exit loop_iperi
           !      end if
           !   end do loop_iperi
           !end if
        end if
     end if

  end select

end subroutine nsi_ifconf
