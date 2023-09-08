subroutine tem_boumat(&
     pnode,pnodb, igaub,iboun,lboel,lelbo,xmmat,xmrhs,gbsha,&
     gbsur,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* temper/tem_boumat
  ! NAME 
  !    tem_boumat
  ! DESCRIPTION
  !    Assemble boundary contribution
  ! USES
  ! USED BY
  !    tem_bouope
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ltype, nnode, lnods
  use def_kermod, only       :  kfl_waexl_ker, kfl_waexl_imp_ker, lexlo_ker, shape_waexl, temel_ker
  use def_master, only       :  tempe

  implicit none
  integer(ip), intent(in)    :: pnode,pnodb, igaub, iboun
  integer(ip), intent(in)    :: lboel(pnodb)
  integer(ip), intent(in)    :: lelbo
  real(rp),    intent(in)    :: xmmat,xmrhs,gbsha(pnodb),gbsur
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inodb,jnodb,inode,jnode, ielem, pelty, index
  real(rp)                   :: xmuit, wetem

  if ( kfl_waexl_ker == 0.or.kfl_waexl_imp_ker==0) then 
     do inodb=1,pnodb  
        inode=lboel(inodb)
        elrhs(inode)=elrhs(inode)+gbsha(inodb)*xmrhs*gbsur
        xmuit=xmmat*gbsha(inodb)*gbsur
        do jnodb=1,pnodb
           jnode=lboel(jnodb)
           elmat(inode,jnode)=elmat(inode,jnode)&
                +xmuit*gbsha(jnodb)
        end do
     end do
  else      ! WALL EXCHANGHE LOCATION implicit 
     ielem = lelbo
     pelty = ltype(ielem)
     index = lexlo_ker(igaub,iboun)
     do jnode =1, pnode                
        xmuit=xmmat* shape_waexl(jnode, index)*gbsur
        wetem = wetem +  shape_waexl(jnode, index)*tempe(lnods(jnode, ielem), 1)
        do inodb = 1,pnodb
           inode=lboel(inodb) 
           elmat(inode,jnode)=elmat(inode,jnode)&
                +xmuit*gbsha(inodb)           
        end do
     end do
  end if

end subroutine tem_boumat
