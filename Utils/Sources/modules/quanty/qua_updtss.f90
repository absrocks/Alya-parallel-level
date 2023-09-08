subroutine qua_updtss()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_updtss
  ! NAME 
  !    qua_updtss
  ! DESCRIPTION
  !    This routine computes the time step size 
  ! USES  
  !   qua_updunk
  ! USED BY
  !    qua_timste
  !    qua_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  use mod_communications, only : PAR_MIN
  implicit none 
  integer(ip)      :: ielem,idime,inode,ipoin
  integer(ip)      :: pnode,pmate,pelty,porde
  real(rp)         :: dtcri,dtmin
  real(rp)         :: chale(2),hleng(3),tragl(9),chave(6)
  real(rp)         :: cartd(ndime,mnode) 
  real(rp)         :: elcod(ndime,mnode),elvel(ndime,mnode)
  real(rp)         :: eledd(mnode),eltem(mnode)
  real(rp), target :: dtpar(1)

  if( kfl_timei_qua /= 0 ) then

     dtmin = 1e6

     if( INOTMASTER ) then

        do ielem = 1,nelem
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           porde = lorde(pelty)
           pmate = 1
           if(nmate>1) pmate=lmate(ielem)
           !
           ! Gather
           !
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
              end do
           end do

           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              eltem(inode)=rhoon(ipoin,1) 
           end do

        end do
     end if
     !
     ! Look for minimum over whole mesh
     !
     if( IPARALL ) call PAR_MIN(dtmin)

     !dtcri_qua = dtmin
     !if(dtcri_qua/=0.0_rp) dtinv_qua = 1.0_rp/(dtcri_qua*safet_qua)
     !if(kfl_timco==1) dtinv=max(dtinv,dtinv_qua)

  end if

end subroutine qua_updtss
