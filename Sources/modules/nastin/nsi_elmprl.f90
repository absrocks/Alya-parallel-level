subroutine nsi_elmprl(&
     itask,pnode,pgaus,igaui,igauf,elden,elvis,gpden,gpvis,&
     gpgvi,gpsha,gpcar,lnods,hleng)
  !****f* Nastin/nsi_elmprl
  ! NAME  
  !    nsi_elmprl
  ! DESCRIPTION
  !    Compute properties for bifluid:
  !            1 ... GPDEN=Density
  !            2 ... GPVIS=Viscosity
  !
  ! USED BY
  !    nsi_elmope
  !***
  use def_kintyp, only     :  ip,rp  
  use def_domain, only     :  mnode,ndime
  use def_nastin, only     :  kfl_grvis_nsi
  implicit none
  integer(ip), intent(in)  :: itask,pnode,pgaus,igaui,igauf
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: hleng(ndime)
  real(rp),    intent(in)  :: elden(pnode,*),elvis(pnode,*)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(out) :: gpden(pgaus,*),gpvis(pgaus)
  real(rp),    intent(out) :: gpgvi(ndime,pgaus)
  integer(ip)              :: igaus,inode,idime,jdime,ipoin

  do igaus=igaui,igauf
     gpden(igaus,1)=0.0_rp
  end do
  do inode=1,pnode
     do igaus=igaui,igauf
        gpden(igaus,1)=gpden(igaus,1)&
             +gpsha(inode,igaus)*elden(inode,1)
     end do
  end do
  do igaus=igaui,igauf
     gpvis(igaus)=0.0_rp
  end do
  do inode=1,pnode
     do igaus=igaui,igauf
        gpvis(igaus)=gpvis(igaus)&
             +gpsha(inode,igaus)*elvis(inode,1)
     end do
  end do

  if(itask==1_ip) then
     kfl_grvis_nsi=1
     do igaus=igaui,igauf
        do idime=1,ndime
           gpgvi(idime,igaus)=0.0_rp
        end do
     end do
     do inode=1,pnode
        do igaus=igaui,igauf
           do idime=1,ndime
              gpgvi(idime,igaus)=gpgvi(idime,igaus)&
                   +gpcar(idime,inode,igaus)*elvis(inode,1)
           end do
        end do
     end do
  end if

end subroutine nsi_elmprl
