subroutine qua_iniunk()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_iniunk
  ! NAME 
  !    qua_iniunk
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    Quanty
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_quanty
  implicit none
  integer(ip)  :: ipoin,idime,inode,ielem,izmat,itask
  integer(ip)  :: pgaus,pnode,pelty
  real(rp)     :: gpvol(mgaus),dummr
  real(rp)     :: gpcar(ndime,mnode,mgaus)
  real(rp)     :: elcod(ndime,mnode),elmat(mnode,mnode),elrhs(mnode)
  !
  ! Compute B matrix
  !
  if( INOTMASTER ) then

      if( eigen_qua(1)%kfl_massm == 0 ) then
         !
         ! Diagonal mass matrix
         !
         do ipoin = 1,npoin
            bmatr(ipoin) = vmasc(ipoin)
         end do

      else if( eigen_qua(1)%kfl_massm == 1 ) then       
        !
        ! Consistent mass matrix
        !
        do izmat = 1,solve_qua(1)%nzmat
           bmatr(izmat) = 0.0_rp
        end do
        elements: do ielem=1,nelem
           pelty = ltype(ielem) 
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           call elmcar(&
                pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                elmar(pelty)%deriv,dummr,elcod,gpvol,gpcar,dummr,ielem)
           call elmmas(pnode,pgaus,gpvol,elmar(pelty)%shape,elmat)
		   itask=0_ip
           call qua_elmdir(itask,pnode,lnods(1,ielem),elmat,elrhs)
           call assmat(&
                solve_qua(1)%ndofn,pnode,pnode,npoin,solve_qua(1)%kfl_algso,&
                ielem,lnods(1,ielem),elmat,bmatr)
        end do elements

      end if

  end if

end subroutine qua_iniunk
