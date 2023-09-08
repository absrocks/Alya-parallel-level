!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_elenod.f90
!> @author  Mariano Vazquez
!> @date    04/05/2016   
!> @brief   Projection DONE ONLY IF INOTMASTER 
!> @details Projection from element to node DONE ONLY IF INOTMASTER 
!> @} 
!-----------------------------------------------------------------------
subroutine sld_elenod(ktask,pnode,pgaus,ielem,gpvol,xshap)
  use def_kintyp, only : ip,rp
  use def_master, only : dtime,rhsid
  use def_domain, only : nelem,npoin, lnods, vmass
  use def_solidz, only : vdiag_sld, dttau_sld, accel_sld, kfl_pseud_sld,ndofn_sld,vmass_sld
  implicit none
  integer(ip) :: ktask  !> ktask=1,2
  integer(ip) :: pnode  !> pnode
  integer(ip) :: pgaus  !> pgaus
  integer(ip) :: ielem  !> ielem
  real(rp)    :: gpvol(pgaus)   !> dvolu * weight
  real(rp)    :: xshap(pnode,pgaus)   !> shape
  integer(ip) :: inode,igaus,ipoin,idofn,itott
  real(rp)    :: dt2, dtps2


  ! vmass: geometrical lumped mass matrix
  ! vmass_sld: geometrical lumped mass matrix times density


  if (kfl_pseud_sld == 0) return

  if (ktask==1) then
     !
     ! Project, done on each element and projecting to the nodes through asssembly
     !
call runend('MARIANO: CHEQUEA ESTO, Y SACA VMASS A FUERA EN SLD_MATRIX')
     dt2 = dtime*dtime     
     dtps2 = dttau_sld(ielem) * dttau_sld(ielem)
     do igaus = 1,pgaus
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
#ifdef NO_COLORING
           !$OMP ATOMIC
#endif
           vdiag_sld(ipoin) = vdiag_sld(ipoin) + (dt2/dtps2) * xshap(inode,igaus) * gpvol(igaus) / vmass(ipoin)
        end do
     end do

  else if (ktask==3) then
     !
     ! Compute, compute and correct the rhsid for the pseudotime step explicit scheme
     !
     ! a^i+1 = (rhsid + P a^i) / (P + 1)
     !
     
     do ipoin=1,npoin
        do idofn=1,ndofn_sld
           itott= (ipoin-1)*ndofn_sld + idofn
#ifdef NO_COLORING
           !$OMP ATOMIC
#endif
           rhsid(itott)= rhsid(itott) + vmass_sld(ipoin) * vdiag_sld(ipoin) * accel_sld(idofn,ipoin,1) 
        end do
     end do


  end if

end subroutine sld_elenod
