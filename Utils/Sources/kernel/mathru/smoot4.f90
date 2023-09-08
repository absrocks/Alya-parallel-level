!------------------------------------------------------------------------
!> @addtogroup Mathru 
!> @{
!> @file    smoot3.f90
!> @date    17/11/2014
!> @author  Mariano Vazquez
!> @brief   Idem Smooth but computing local velem from vsour and iterating
!> @details Idem Smooth but computing local velem from vsour and iterating
!> @} 
!------------------------------------------------------------------------
  subroutine smoot4(ndofn_smoo,vsour,vmass_smoo,vtarg)
    use def_kintyp, only              :  ip,rp,r1p
    use def_master, only              :  INOTMASTER,mem_modul,modul
    use def_domain, only              :  ndime,npoin,nelem,nnode,mnode,ntens
    use def_domain, only              :  lnods,ltype,coord,elmar
    use def_domain, only              :  lexis,ngaus,kfl_naxis,mgaus,lelch
    use def_domain, only              :  lmate,lmatn,nmate
    use def_elmtyp  
    use mod_memchk
    implicit none
    integer(ip),intent(in)            :: ndofn_smoo    !> local size of the vector
    real(rp)   ,intent(in)            :: vsour(*)      !> source vector, the vector to be smoothed
    real(rp)   ,intent(in)            :: vmass_smoo(*) !> lumped mass matrix
    real(rp)   ,intent(out)           :: vtarg(*) !> target vector, the smoothed source vector
    integer(ip)                       :: ipoin,jpoin,idime,inode,jnode,ielem,igaus,itott,jtott
    integer(ip)                       :: pnode,pelty,pgaus
    integer(ip)                       :: kdofn
    real(rp)                          :: detjm,gpvol,gpcar(ndime,mnode)
    real(rp)                          :: elcod(ndime,mnode)
    real(rp)                          :: xjaci(9),xjacm(9),shapi,shapj

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do itott = 1,(npoin*ndofn_smoo)
          vtarg(itott)     = 0.0_rp
       end do
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          pelty = ltype(ielem)
          if( pelty > 0 ) then
             pnode = nnode(pelty)
             pgaus = ngaus(pelty)
             !
             ! Gather vectors
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                do idime = 1,ndime
                   elcod(idime,inode) = coord(idime,ipoin)
                end do
             end do
             !
             ! Loop over Gauss points 
             !
             gauss_points: do igaus = 1,pgaus
                call elmder(pnode,ndime,elmar(pelty)%deriv(1,1,igaus),elcod,gpcar,detjm,xjacm,xjaci)
                gpvol = elmar(pelty) % weigp(igaus) * detjm
                !
                ! Assemble
                !
                do kdofn = 1,ndofn_smoo
                   do inode = 1,pnode
                      ipoin        = lnods(inode,ielem)
                      itott        = (ipoin-1)*ndofn_smoo+kdofn
                      shapi        = elmar(pelty) % shape(inode,igaus)
                      do jnode = 1,pnode
                         jpoin        = lnods(jnode,ielem)
                         jtott        = (jpoin-1)*ndofn_smoo+kdofn
                         shapj        = elmar(pelty) % shape(jnode,igaus)
                         vtarg(itott) = vtarg(itott) + vsour(jtott) * shapj * shapi * gpvol / vmass_smoo(ipoin)
                      end do
                   end do
                end do

             end do gauss_points
          end if
       end do elements
       !
       ! Parallelization
       !
       call rhsmod(ndofn_smoo,vtarg)

    end if

  end subroutine smoot4
