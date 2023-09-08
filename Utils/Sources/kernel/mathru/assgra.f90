subroutine assgra(&
     pnode,pgaus,lnods,gpden,gpvis,gpvol,gpsha,gpcar,&
     elunk,elrhs,rhsid)
  !-----------------------------------------------------------------------
  !****f* mathru/assgra
  ! NAME 
  !    assgra
  ! DESCRIPTION
  !    Assembly of properties
  ! USES
  ! USED BY
  !    *_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only        :  ip,rp
  use def_domain, only        :  mnode,ntens,ndime
  implicit none
  integer(ip),  intent(in)    :: pnode,pgaus
  integer(ip),  intent(in)    :: lnods(pnode)
  real(rp),     intent(in)    :: gpden(pgaus)
  real(rp),     intent(in)    :: gpvis(pgaus)
  real(rp),     intent(in)    :: gpvol(pgaus)
  real(rp),     intent(in)    :: gpsha(pnode,pgaus)
  real(rp),     intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)    :: elunk(ndime,pnode)
  real(rp),     intent(out)   :: elrhs(ntens,pnode)
  real(rp),     intent(inout) :: rhsid(*)
  integer(ip)                 :: inode,igaus,itens,idime,jdime
  real(rp)                    :: fact1,elgra(3,3)

  do inode = 1,pnode
     do itens = 1,ntens
        elrhs(itens,inode) = 0.0_rp
     end do
  end do

  do igaus = 1,pgaus

     elgra(1:ndime,1:ndime) = 0.0_rp

     do idime = 1,ndime
        do jdime = 1,ndime
           do inode = 1,pnode
              elgra(idime,jdime) = elgra(idime,jdime) &
                   + gpcar(idime,inode,igaus) * elunk(jdime,inode)
           end do
        end do
     end do

     do inode = 1,pnode
        fact1          = gpvol(igaus) * gpvis(igaus) * gpsha(inode,igaus)
        elrhs(1,inode) = elrhs(1,inode) + fact1 * 2.0_rp * elgra(1,1)
        elrhs(2,inode) = elrhs(2,inode) + fact1 * 2.0_rp * elgra(2,2)
        elrhs(3,inode) = elrhs(3,inode) + fact1 * ( elgra(1,2) + elgra(2,1) )
        if( ndime == 3 ) then
           elrhs(4,inode) = elrhs(4,inode) + fact1 * 2.0_rp * elgra(3,3)
           elrhs(5,inode) = elrhs(5,inode) + fact1 * ( elgra(1,3) + elgra(3,1) )
           elrhs(6,inode) = elrhs(6,inode) + fact1 * ( elgra(2,3) + elgra(3,2) )
        end if
     end do

  end do

  call assrhs(ntens,pnode,lnods,elrhs,rhsid)

end subroutine assgra
