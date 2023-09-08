subroutine got_elmdir(&
     itask,pnode,pevat,ndofn,lnods,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmdir
  ! NAME 
  !    got_elmdir
  ! DESCRIPTION
  !    This routines prescribes the boundary conditions on u and alpha
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  cdrop
  use def_gotita, only       :  kfl_fixno_got,bvess_got,almax_got,&
       &                        cutof_got,ndofn_got,kfl_probl_got
  implicit none
  integer(ip), intent(in)    :: itask,pnode,pevat,ndofn
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pevat,pevat),elrhs(pevat)
  real(rp)                   :: adiag
  integer(ip)                :: inode,ipoin,idofn,ievat,jevat,jnode
  !
  ! Inflow boundary conditions: u, alpha or (u,alpha)
  !
  if(itask==1.or.itask==2) then
     do inode=1,pnode
        ipoin=lnods(inode)        
        if(kfl_fixno_got(ipoin)==1) then
           ievat=(inode-1)*ndofn
           do idofn=1,ndofn
              ievat=ievat+1
              adiag=elmat(ievat,ievat)
              if(adiag==0.0_rp) adiag=1.0_rp
              do jevat=1,pevat
                 elmat(ievat,jevat)=0.0_rp
              end do
              elmat(ievat,ievat)=adiag
              elrhs(ievat)=adiag*bvess_got(idofn,ipoin)
           end do
        end if
     end do
  end if

  if(kfl_probl_got/=2) then
     !
     ! Inflow boundary conditions: alpha 
     !
     if(itask==3) then
        do inode=1,pnode
           ipoin=lnods(inode)        
           if(kfl_fixno_got(ipoin)==1) then
              adiag=elmat(inode,inode)
              if(adiag==0.0_rp) adiag=1.0_rp
              do jnode=1,pnode
                 elmat(inode,jnode)=0.0_rp
              end do
              elmat(inode,inode)=adiag
              elrhs(inode)=adiag*bvess_got(ndofn_got(3),ipoin)
           end if
        end do
     end if
     !
     ! Iterative inlet on alpha: imposed when u.n<0 in got_updfix
     !
     if(itask==1.or.itask==3) then
        do inode=1,pnode
           ipoin=lnods(inode)
           if(kfl_fixno_got(ipoin)==2) then
              ievat=inode*ndofn
              adiag=elmat(ievat,ievat)
              if(adiag==0.0_rp) adiag=1.0_rp
              do jevat=1,pevat
                 elmat(ievat,jevat)=0.0_rp
              end do
              elmat(ievat,ievat)=adiag
              elrhs(ievat)=adiag*almax_got*cutof_got           
           end if
        end do
     end if
  end if

end subroutine got_elmdir
