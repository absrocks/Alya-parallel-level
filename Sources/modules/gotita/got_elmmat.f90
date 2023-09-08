subroutine got_elmmat(&
     itask,pgaus,pnode,pevat,ndofn,gpvol,resim,resic,&
     tesmv,tesmb,tescv,tescb,gprhs,gpcar,gpgvd,gpgcd,&
     gpdif,elvdr,elcdr,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmmat
  ! NAME 
  !    got_elmmat
  ! DESCRIPTION
  !    Assemble the elemental matrix from Gauss point contributions
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  use def_gotita, only     :  ndofn_got,kfl_artif_got,artif_got,&
       &                      kfl_staty_got,kfl_diffu_got,kfl_coupl_got
  implicit none
  integer(ip), intent(in)  :: itask,pgaus,pnode,pevat,ndofn
  real(rp),    intent(in)  :: gpvol(pgaus)
  real(rp),    intent(in)  :: resim(ndofn_got(3),pnode,pgaus)
  real(rp),    intent(in)  :: resic(ndofn_got(3),pnode,pgaus)
  real(rp),    intent(in)  :: tesmv(pnode,pgaus)
  real(rp),    intent(in)  :: tesmb(ndime,pnode,pgaus)
  real(rp),    intent(in)  :: tescv(ndime,pnode,pgaus)
  real(rp),    intent(in)  :: tescb(pnode,pgaus)
  real(rp),    intent(in)  :: gprhs(ndofn_got(3),pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gpgvd(ndime,ndime,pgaus)
  real(rp),    intent(in)  :: gpgcd(ndime,pgaus),gpdif(pgaus)
  real(rp),    intent(in)  :: elvdr(ndime,pnode),elcdr(pnode)
  real(rp),    intent(out) :: elmat(pevat,pevat)
  real(rp),    intent(out) :: elrhs(pevat)
  integer(ip)              :: igaus,kdime,idime,inode,jnode
  integer(ip)              :: idofn,jdofn,ievat,jevat
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

  !----------------------------------------------------------------------
  !
  ! COUPLED EQUATIONS
  !
  !----------------------------------------------------------------------

  if(kfl_coupl_got==1) then

     if(kfl_staty_got==0.or.kfl_staty_got==1) then
        !
        ! NO STABILIZATION AND SUPG: MOMENTUM MATRIX
        !
        if(itask==1) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofn=(inode-1)*ndofn
                 do idime=1,ndime
                    idofn=idofn+1
                    fact1=gpvol(igaus)*tesmv(inode,igaus)
                    kdime=1+idime
                    do jnode=1,pnode
                       jdofn=(jnode-1)*ndofn+idime
                       elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                            +fact1*resim(1,jnode,igaus)
                       jdofn=jnode*ndofn
                       elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                            +fact1*resim(kdime,jnode,igaus)
                    end do
                 end do
              end do
           end do
        else if(itask==2) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofn=(inode-1)*ndofn
                 do idime=1,ndime
                    idofn=idofn+1
                    fact1=gpvol(igaus)*tesmv(inode,igaus)
                    kdime=1+idime
                    do jnode=1,pnode
                       jdofn=(jnode-1)*ndofn+idime
                       elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                            +fact1*resim(1,jnode,igaus)
                       elrhs(idofn)=elrhs(idofn)&
                            -fact1*resim(kdime,jnode,igaus)&
                            *elcdr(jnode)
                    end do
                 end do
              end do
           end do
        end if
        !
        ! NO STABILIZATION AND SUPG: CONTINUITY MATRIX
        !
        if(itask==1) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofn=inode*ndofn
                 fact1=gpvol(igaus)*tescb(inode,igaus)
                 do jnode=1,pnode
                    jdofn=(jnode-1)*ndofn
                    do idime=1,ndime
                       jdofn=jdofn+1
                       elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                            +fact1*resic(idime,jnode,igaus)
                    end do
                    jdofn=jdofn+1
                    elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                         +fact1*resic(ndofn_got(3),jnode,igaus)
                 end do
              end do
           end do
        else if(itask==3) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofn=inode*ndofn
                 fact1=gpvol(igaus)*tescb(inode,igaus)
                 do jnode=1,pnode
                    do idime=1,ndime
                       elrhs(idofn)=elrhs(idofn)&
                            -fact1*resic(idime,jnode,igaus)&
                            *elvdr(idime,jnode)
                    end do
                    jdofn=jnode
                    elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                         +fact1*resic(ndofn_got(3),jnode,igaus)
                 end do
              end do
           end do
        end if
        !
        ! NO STABILIZATION AND SUPG: MOMENTUM RHS
        !
        if(itask==1.or.itask==2) then
           do igaus=1,pgaus
              do inode=1,pnode     
                 idofn=(inode-1)*ndofn
                 fact1=gpvol(igaus)*tesmv(inode,igaus)
                 do idime=1,ndime
                    idofn=idofn+1
                    elrhs(idofn)=elrhs(idofn)&
                         +fact1*gprhs(idime,igaus)
                 end do
              end do
           end do
        end if
        !
        ! NO STABILIZATION AND SUPG: CONTINUITY RHS
        !
        if(itask==1.or.itask==3) then
           do igaus=1,pgaus
              do inode=1,pnode     
                 idofn=(inode-1)*ndofn
                 fact1=gpvol(igaus)*tesmv(inode,igaus)
                 idofn=inode*ndofn
                 elrhs(idofn)=elrhs(idofn)+gpvol(igaus)&
                      *gprhs(ndofn_got(3),igaus)*tescb(inode,igaus)
              end do
           end do
        end if

     else if(kfl_staty_got==2.or.kfl_staty_got==3) then
        !
        !           +---+---+---+
        !           | x |   |   |
        !           +---+---+---+
        ! MOMENTUM: |   | x |   | 
        !           +---+---+---+
        !           |   |   |   |
        !           +---+---+---+
        !
        if(itask==1.or.itask==2) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofn=(inode-1)*ndofn
                 do idime=1,ndime
                    idofn=idofn+1
                    fact1=gpvol(igaus)*tesmv(inode,igaus)
                    fact2=gpvol(igaus)*tescv(idime,inode,igaus)
                    do jnode=1,pnode
                       jdofn=(jnode-1)*ndofn+idime
                       elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                            +fact1*resim(1,    jnode,igaus)&
                            +fact2*resic(idime,jnode,igaus)
                    end do
                 end do
              end do
           end do
        end if
        !
        !           +---+---+---+
        !           |   |   | x |
        !           +---+---+---+
        ! MOMENTUM: |   |   | x |
        !           +---+---+---+
        !           |   |   |   |
        !           +---+---+---+
        !
        if(itask==1) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofn=(inode-1)*ndofn
                 fact1=gpvol(igaus)*tesmv(inode,igaus)
                 do idime=1,ndime
                    idofn=idofn+1
                    kdime=1+idime
                    fact2=gpvol(igaus)*tescv(idime,inode,igaus)
                    do jnode=1,pnode
                       jdofn=jnode*ndofn
                       elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                            +fact1*resim(kdime,jnode,igaus)&
                            +fact2*resic(ndofn_got(3),jnode,igaus)
                    end do
                 end do
              end do
           end do
        else if(itask==2) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofn=(inode-1)*ndofn
                 fact1=gpvol(igaus)*tesmv(inode,igaus)
                 do idime=1,ndime
                    idofn=idofn+1
                    kdime=1+idime
                    fact2=gpvol(igaus)*tescv(idime,inode,igaus)
                    do jnode=1,pnode
                       elrhs(idofn)=elrhs(idofn)&
                            &  -( fact1*resim(kdime,jnode,igaus)&
                            &    +fact2*resic(ndofn_got(3),jnode,igaus))*elcdr(jnode)
                    end do
                 end do
              end do
           end do
        end if
        !
        !             +---+---+---+
        !             |   |   |   |
        !             +---+---+---+
        ! CONTINUITY: |   |   |   |
        !             +---+---+---+
        !             | x | x |   |
        !             +---+---+---+
        !
        if(itask==1) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofn=inode*ndofn
                 fact2=gpvol(igaus)*tescb(inode,igaus)
                 do idime=1,ndime
                    fact1=gpvol(igaus)*tesmb(idime,inode,igaus)
                    do jnode=1,pnode
                       jdofn=(jnode-1)*ndofn+idime
                       elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                            +fact1*resim(    1,jnode,igaus)&
                            +fact2*resic(idime,jnode,igaus)
                    end do
                 end do
              end do
           end do
        else if(itask==3) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofn=inode*ndofn
                 fact2=gpvol(igaus)*tescb(inode,igaus)
                 do idime=1,ndime
                    fact1=gpvol(igaus)*tesmb(idime,inode,igaus)
                    do jnode=1,pnode
                       jdofn=(jnode-1)*ndofn+idime
                       elrhs(idofn)=elrhs(idofn)&
                            & -( fact1*resim(    1,jnode,igaus)&
                            &   +fact2*resic(idime,jnode,igaus))&
                            *elvdr(idime,jnode)
                    end do
                 end do
              end do
           end do
        end if
        !
        !             +---+---+---+
        !             |   |   |   |
        !             +---+---+---+
        ! CONTINUITY: |   |   |   |
        !             +---+---+---+
        !             |   |   | x |
        !             +---+---+---+
        !
        do igaus=1,pgaus 
           do inode=1,pnode
              idofn=inode*ndofn
              fact1=gpvol(igaus)*tescb(inode,igaus)
              do jnode=1,pnode
                 jdofn=jnode*ndofn
                 elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                      +fact1*resic(ndofn_got(3),jnode,igaus)
                 do idime=1,ndime
                    elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                         +gpvol(igaus)*tesmb(idime,inode,igaus)&
                         *resim(1+idime,jnode,igaus)
                 end do
              end do
           end do
        end do
        !
        ! RHS: MOMENTUM
        !
        if(itask==1.or.itask==2) then
           do igaus=1,pgaus
              do inode=1,pnode     
                 idofn=(inode-1)*ndofn
                 jdofn=inode*ndofn
                 do idime=1,ndime
                    idofn=idofn+1
                    elrhs(idofn)=elrhs(idofn)+gpvol(igaus)&
                         *(gprhs(idime,igaus)*tesmv(inode,igaus)&
                         + gprhs(ndofn_got(3),igaus)*tescv(idime,inode,igaus))
                 end do
              end do
           end do
        end if
        !
        ! RHS: CONTINUITY
        !
        if(itask==1.or.itask==3) then
           do igaus=1,pgaus
              do inode=1,pnode     
                 idofn=(inode-1)*ndofn
                 jdofn=inode*ndofn
                 do idime=1,ndime
                    elrhs(jdofn)=elrhs(jdofn)+gpvol(igaus)&
                         *gprhs(idime,igaus)*tesmb(idime,inode,igaus)
                 end do
                 elrhs(jdofn)=elrhs(jdofn)+gpvol(igaus)&
                      *gprhs(ndofn_got(3),igaus)*tescb(inode,igaus)
              end do
           end do
        end if
     end if

  !----------------------------------------------------------------------
  !
  ! UNCOUPLED EQUATIONS
  !
  !----------------------------------------------------------------------

  else if(kfl_coupl_got==0) then
     !
     !           +---+---+---+
     !           | x |   |   |
     !           +---+---+---+
     ! MOMENTUM: |   | x |   | 
     !           +---+---+---+
     !           |   |   |   |
     !           +---+---+---+
     !
     if(itask==1.or.itask==2) then
        do igaus=1,pgaus
           do inode=1,pnode
              idofn=(inode-1)*ndofn
              do idime=1,ndime
                 idofn=idofn+1
                 fact1=gpvol(igaus)*tesmv(inode,igaus)
                 do jnode=1,pnode
                    jdofn=(jnode-1)*ndofn+idime
                    elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                         +fact1*resim(1,jnode,igaus)
                 end do
              end do
           end do
        end do
     end if
     !
     !             +---+---+---+
     !             |   |   |   |
     !             +---+---+---+
     ! CONTINUITY: |   |   |   |
     !             +---+---+---+
     !             |   |   | x |
     !             +---+---+---+
     !
     do igaus=1,pgaus 
        do inode=1,pnode
           idofn=inode*ndofn
           fact1=gpvol(igaus)*tescb(inode,igaus)
           do jnode=1,pnode
              jdofn=jnode*ndofn
              elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                   +fact1*resic(ndofn_got(3),jnode,igaus)
           end do
        end do
     end do
     !
     ! RHS: MOMENTUM
     !
     if(itask==1.or.itask==2) then
        do igaus=1,pgaus
           do inode=1,pnode     
              idofn=(inode-1)*ndofn
              do idime=1,ndime
                 idofn=idofn+1
                 elrhs(idofn)=elrhs(idofn)+gpvol(igaus)&
                      *gprhs(idime,igaus)*tesmv(inode,igaus)
              end do
           end do
        end do
     end if
     !
     ! RHS: CONTINUITY
     !
     if(itask==1.or.itask==3) then
        do igaus=1,pgaus
           do inode=1,pnode     
              idofn=inode*ndofn
              elrhs(idofn)=elrhs(idofn)+gpvol(igaus)&
                   *gprhs(ndofn_got(3),igaus)*tescb(inode,igaus)
           end do
        end do
     end if
  end if
  !
  ! MOMENTUM: Diffusion
  !
  if((itask==1.or.itask==2).and.kfl_diffu_got/=0) then
     do igaus=1,pgaus
        fact1=gpdif(igaus)*gpvol(igaus)
        do inode=1,pnode
           idofn=(inode-1)*ndofn
           do idime=1,ndime
              idofn=idofn+1
              do jnode=1,pnode
                 jdofn=(jnode-1)*ndofn+idime
                 do kdime=1,ndime
                    elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                         +fact1*gpcar(kdime,inode,igaus)&
                         *gpcar(kdime,jnode,igaus)
                 end do
              end do
           end do
        end do
     end do
  end if
  !
  ! MOMENTUM: Artificial viscosity
  !
  if(itask==1.or.itask==2) then
     if(kfl_artif_got/=0.and.artif_got(1)/=0.0_rp) then
        do igaus=1,pgaus
           fact1=artif_got(1)*gpvol(igaus)
           do inode=1,pnode
              idofn=(inode-1)*ndofn
              do idime=1,ndime
                 idofn=idofn+1
                 do jnode=1,pnode
                    jdofn=(jnode-1)*ndofn+idime
                    do kdime=1,ndime
                       elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                            +fact1*gpcar(kdime,inode,igaus)&
                            *gpcar(kdime,jnode,igaus)
                    end do
                 end do
              end do
           end do
        end do
     end if
  end if
  !
  ! CONTINUITY: Artificial viscosity
  !
  if(itask==1.or.itask==3) then
     if(kfl_artif_got/=0.and.artif_got(2)/=0.0_rp) then
        do igaus=1,pgaus
           fact1=artif_got(2)*gpvol(igaus)
           do inode=1,pnode
              idofn=inode*ndofn
              do jnode=1,pnode
                 jdofn=jnode*ndofn
                 do kdime=1,ndime
                    elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                         +fact1*gpcar(kdime,inode,igaus)&
                         *gpcar(kdime,jnode,igaus)
                 end do
              end do
           end do
        end do
     end if
  end if
  !
  ! RHS iterative artificial viscosity
  !
  if(kfl_artif_got==2) then

     if(itask==1.or.itask==2) then
        do igaus=1,pgaus
           fact1=artif_got(1)*gpvol(igaus)
           do inode=1,pnode     
              idofn=(inode-1)*ndofn
              do idime=1,ndime
                 idofn=idofn+1
                 do kdime=1,ndime
                    elrhs(idofn)=elrhs(idofn)&
                         +fact1*gpgvd(kdime,idime,igaus)&
                         *gpcar(kdime,inode,igaus)
                 end do
              end do
           end do
        end do
     end if

     if(itask==1.or.itask==3) then
        do igaus=1,pgaus
           fact2=artif_got(2)*gpvol(igaus)
           do inode=1,pnode     
              idofn=inode*ndofn
              do kdime=1,ndime
                 elrhs(idofn)=elrhs(idofn)&
                      +fact2*gpgcd(kdime,igaus)&
                      *gpcar(kdime,inode,igaus)
              end do
           end do
        end do
     end if

  end if

end subroutine got_elmmat
