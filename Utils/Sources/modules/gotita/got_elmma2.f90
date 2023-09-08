subroutine got_elmma2(&
     pgaus,pnode,pevat,ndofr,gpvol,resim,tesmv,gprhs,&
     elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmma2
  ! NAME 
  !    got_elmma2
  ! DESCRIPTION
  !    Assemble the elemental matrix from Gauss point contributions
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  use def_gotita, only     :  kfl_artif_got,artif_got,&
       &                      kfl_staty_got,kfl_diffu_got,kfl_probl_got,&
       &                      kfl_linea_got
  implicit none
  integer(ip), intent(in)  :: pgaus,pnode,pevat,ndofr
  real(rp),    intent(in)  :: gpvol(pgaus),resim(ndofr,pnode,pgaus)
  real(rp),    intent(in)  :: tesmv(pnode,pgaus),gprhs(ndime,pgaus)
  real(rp),    intent(out) :: elmat(pevat,pevat),elrhs(pevat)
  integer(ip)              :: igaus,idime,inode,jnode
  integer(ip)              :: idofn,jdofn,ievat,jevat,jdime,idofr
  real(rp)                 :: fact1,fact2
  !
  ! Initialization
  !
  do ievat=1,pevat
     elrhs(ievat)=0.0_rp
     do jevat=1,pevat
        elmat(jevat,ievat)=0.0_rp
     end do
  end do

  if(kfl_linea_got==1) then
     !
     ! Picard linearization
     !
     do igaus=1,pgaus
        do inode=1,pnode
           idofn=(inode-1)*ndime
           fact1=gpvol(igaus)*tesmv(inode,igaus)
           do idime=1,ndime
              idofn=idofn+1
              do jnode=1,pnode
                 jdofn=(jnode-1)*ndime+idime
                 elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                      +fact1*resim(1,jnode,igaus)
              end do
              elrhs(idofn)=elrhs(idofn)+fact1*gprhs(idime,igaus)
           end do
        end do
     end do

     !do igaus=1,pgaus
     !   do inode=1,pnode
     !      fact1=gpvol(igaus)*tesmv(inode,igaus)           
     !      do jnode=1,pnode
     !         fact2=fact1*resim(1,jnode,igaus)
     !         jdofn=(jnode-1)*ndime
     !         do idime=1,ndime
     !            idofn=(inode-1)*ndime+idime
     !            jdofn=jdofn+1
     !            elmat(idofn,jdofn)=elmat(idofn,jdofn)+fact2
     !         end do
     !      end do
     !   end do
     !end do

  else 
     !
     ! Newton-Raphson linearization
     !
     do igaus=1,pgaus
        do inode=1,pnode
           idofn=(inode-1)*ndime
           fact1=gpvol(igaus)*tesmv(inode,igaus)
           do idime=1,ndime
              idofn=idofn+1
              do jnode=1,pnode
                 jdofn=(jnode-1)*ndime                 
                 do jdime=1,ndime
                    idofr=(idime-1)*ndime+jdime
                    jdofn=jdofn+1
                    elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                         +fact1*resim(idofr,jnode,igaus)
                 end do
              end do
              elrhs(idofn)=elrhs(idofn)+fact1*gprhs(idime,igaus)
           end do
        end do
     end do
  end if

end subroutine got_elmma2
