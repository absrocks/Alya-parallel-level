!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_elmset.f90
!> @author  Alfonso Santiago and Gerard Guillamet
!> @date    July, 2017
!>          - Subroutine creation
!> @author  Alfonso Santiago
!> @date    July, 2017
!>          - Sets invariants and Green components
!> @author  Gerard Guillamet and Guido Giuntoli
!> @date    March, 2019
!>          - Subroutine refactoring
!> @brief   This routine computes variables on an element set.
!> @details
!>          The variable are:
!>          1. SETVO:     Volume of the set.
!>          2. SETI1:     Set of the 1st invariant Green Lagrange
!>          3. SETI2:     Set of the 2nd invariant Green Lagrange
!>          4. SETI3:     Set of the 3rd invariant Green Lagrange
!>
!>         14. SETEP_AVE: Set average log strain in polar cylindrical (eps_rr + eps_tt)
!>         15. SETVE:     Set mean velocity
!>         16. SETMA:     Set total mass
!>         17. SETKI:     Set Kinetic energy
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_elmset(iesec,ieset)

  use def_kintyp,   only : ip, rp
  use def_master,   only : ITER_K
  use def_master,   only : postp, veset
  use def_domain,   only : ndime, mnode, nnode, mgaus, ngaus, nelem
  use def_domain,   only : leset, ltype, lnods, lmate
  use def_domain,   only : elmar, coord
  use def_domain,   only : lobas
  use def_solidz,   only : densi_sld, veloc_sld
  use def_solidz,   only : grlst_sld, lepse_sld, epsel_sld
  use mod_sld_csys, only : sld_csys_tensor_cylindrical_to_cartesian

  implicit none

  integer(ip), intent(in)  :: iesec                     !< Element number
  integer(ip), intent(in)  :: ieset                     !< Code element set
  real(rp),    pointer     :: setvo(:),seti1(:),seti2(:),seti3(:)
  real(rp),    pointer     :: setgrl_11(:), setgrl_12(:), setgrl_13(:)
  real(rp),    pointer     :: setgrl_21(:), setgrl_22(:), setgrl_23(:)
  real(rp),    pointer     :: setgrl_31(:), setgrl_32(:), setgrl_33(:)
  real(rp),    pointer     :: setep_ave(:)              !< Sum of logarithmic strain rr + tt components
  real(rp),    pointer     :: setve(:)                  !< Mean velocity
  real(rp),    pointer     :: setma(:)                  !< Total Mass
  real(rp),    pointer     :: setki(:)                  !< All Kinetic energy
  real(rp)                 :: elcog(ndime)              !< Element center of gravity
  real(rp)                 :: elcod(ndime,mnode)        !< Nodal coordinates at element level
  real(rp)                 :: elvel(ndime,mnode)        !< Nodal velocity at element level
  real(rp)                 :: csysp(3*ndime)            !< Parameters coordinate system
  real(rp)                 :: lepsi_car(ndime,ndime)    !< Logarithmic strain tensor in cartesian
  real(rp)                 :: lepsi_cyl(ndime,ndime)    !< Logarithmic strain tensor in cylindrical
  integer(ip)              :: pnode,pgaus,pelty,pmate,nvabi
  integer(ip)              :: ielem,igaus,inode,ipoin
  real(rp)                 :: gpcar(ndime,mnode,mgaus)
  real(rp)                 :: xjaci(ndime,ndime),xjacm(ndime,ndime)
  real(rp)                 :: gpvol, gpdet, gpvel(ndime), gpmas, veloc

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  nvabi     =  postp(1) % nvaes + 1
  setvo     => veset( nvabi:nvabi, ieset )
  seti1     => postp(1) % veset(2:2,ieset)
  seti2     => postp(1) % veset(3:3,ieset)
  seti3     => postp(1) % veset(4:4,ieset)
  setgrl_11 => postp(1) % veset(5:5,ieset)
  setgrl_12 => postp(1) % veset(6:6,ieset)
  setgrl_13 => postp(1) % veset(7:7,ieset)
  setgrl_21 => postp(1) % veset(8:8,ieset)
  setgrl_22 => postp(1) % veset(9:9,ieset)
  setgrl_23 => postp(1) % veset(10:10,ieset)
  setgrl_31 => postp(1) % veset(11:11,ieset)
  setgrl_32 => postp(1) % veset(12:12,ieset)
  setgrl_33 => postp(1) % veset(13:13,ieset)
  setep_ave => postp(1) % veset(14:14,ieset)
  setve     => postp(1) % veset(15:15,ieset)
  setma     => postp(1) % veset(16:16,ieset)
  setki     => postp(1) % veset(17:17,ieset)
  setvo     = 0.0_rp  ! Set volume
  seti1     = 0.0_rp  ! Set first invariant of green-lagrange
  seti2     = 0.0_rp  ! Set 2nd invariant of green-lagrange
  seti3     = 0.0_rp  ! Set 3rd invariant of green-lagrange
  setgrl_11 = 0.0_rp  !
  setgrl_12 = 0.0_rp  !
  setgrl_13 = 0.0_rp  !
  setgrl_21 = 0.0_rp  !
  setgrl_22 = 0.0_rp  !
  setgrl_23 = 0.0_rp  !
  setgrl_31 = 0.0_rp  !
  setgrl_32 = 0.0_rp  !
  setgrl_33 = 0.0_rp  !
  setep_ave = 0.0_rp  ! Set average log strain in polar cylindrical (eps_rr + eps_tt)
  setve     = 0.0_rp  ! Set mean velocity
  setma     = 0.0_rp  ! Set all mass
  setki     = 0.0_rp  ! Set all kinetic energy

  elements: do ielem = 1,nelem

     if( leset(ielem) == iesec ) then

        !----------------------------------------------------------------
        !
        ! Element properties, dimensions and gather
        !
        !----------------------------------------------------------------

        pelty = ltype(ielem)
        pmate = lmate(ielem)

        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           !
           ! Gather operations
           !
           elcog(:) = 0.0_rp
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
              elcog(1:ndime)       = elcog(1:ndime) + elmar(pelty)%shacg(inode) * elcod(1:ndime,inode)
              elvel(1:ndime,inode) = veloc_sld(1:ndime,ipoin,ITER_K)
           end do

           !----------------------------------------------------------------
           !
           ! Loop over Gauss points
           !
           !----------------------------------------------------------------

           gauss_points: do igaus = 1, pgaus

              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian

              gpvol = elmar(pelty)%weigp(igaus)*gpdet                 ! dV:=|J|*wg
              setvo = setvo + gpvol

              if( postp(1) % npp_setse(14) /= 0 ) then
                 !
                 ! LEPRT: eps_RR + eps_TT
                 !
                 !lepsi_car(1:ndime,1:ndime) = lepse_sld(ielem)%a(1:ndime,1:ndime,igaus)*gpvol
                 lepsi_car(1:ndime,1:ndime) = epsel_sld(ielem)%a(1:ndime,1:ndime,igaus)*gpvol   ! green
                 !
                 ! Get coordinate system parameters
                 csysp = lobas(iesec) % param(1:3*ndime)
                 ! Logarithmic strain tensor from Cartesian to Polar-Cylindrical
                 lepsi_cyl(:,:) = 0.0_rp
                 call sld_csys_tensor_cylindrical_to_cartesian(2_ip, ndime, &
                      elcog(:), csysp(:), lepsi_car(:,:), lepsi_cyl(:,:))
                 setep_ave = setep_ave + lepsi_cyl(1,1) + lepsi_cyl(2,2)

              end if

              if( postp(1) % npp_setse(15) /= 0 ) then
                 !
                 ! VELOC: Mean velocity
                 !
                 gpvel(:) = 0.0_rp
                 do inode = 1,pnode
                    gpvel(1:ndime) = gpvel(1:ndime) + elmar(pelty)%shape(inode,igaus)*elvel(1:ndime,inode)
                 end do
                 veloc = sqrt(sum(gpvel(1:ndime)**2))
                 setve = setve + gpvol*veloc

              end if

              if( postp(1) % npp_setse(16) /= 0 ) then
                 !
                 ! TMASS: Total mass
                 !
                 gpmas = densi_sld(1,pmate)*gpvol
                 setma = setma + gpmas

              end if

              if( postp(1) % npp_setse(17) /= 0 ) then
                 !
                 ! ALLKE: All kinetic energy
                 !
                 gpvel(:) = 0.0_rp
                 do inode = 1,pnode
                    gpvel(1:ndime) = gpvel(1:ndime) + elmar(pelty)%shape(inode,igaus)*elvel(1:ndime,inode)
                 end do
                 veloc = sqrt(sum(gpvel(1:ndime)**2))
                 gpmas = densi_sld(1,pmate)*gpvol
                 setki = setki + 0.5_rp * gpmas * veloc * veloc

              end if

           end do gauss_points

           if( postp(1) % npp_setse(2) /= 0 ) then
              !
              ! 1st inv GL strain tensor: SET1E
              !
              seti1 = seti1 &
                   +  grlst_sld(ielem,1) &
                   +  grlst_sld(ielem,5)
              if(ndime==3_ip) seti1 = seti1 +  grlst_sld(ielem,9)
           end if
           if( postp(1) % npp_setse(3) /= 0 ) then
              !
              ! 2nd inv GL strain tensor: SET2E
              !
              if(ndime==3_ip) then
                 seti2 = seti2 + grlst_sld(ielem,1)*grlst_sld(ielem,5) &
                      + grlst_sld(ielem,5)*grlst_sld(ielem,9) &
                      + grlst_sld(ielem,9)*grlst_sld(ielem,1) &
                      - grlst_sld(ielem,2)*grlst_sld(ielem,4) &
                      - grlst_sld(ielem,6)*grlst_sld(ielem,8) &
                      - grlst_sld(ielem,3)*grlst_sld(ielem,7)

              elseif(ndime==2_ip) then
                 seti2 = seti2 + grlst_sld(ielem,1)*grlst_sld(ielem,4) &
                      - grlst_sld(ielem,2)*grlst_sld(ielem,3)
              endif

           end if
           if( postp(1) % npp_setse(4) /= 0 ) then
              !
              ! 3rd inv GL strain tensor: SET3E
              !
              if(ndime==2_ip) seti3 = 0.0_rp !2D

              if(ndime==3_ip) then
                 seti3 = seti3 &
                      + grlst_sld(ielem,1)*grlst_sld(ielem,5)*grlst_sld(ielem,9) &
                      + grlst_sld(ielem,2)*grlst_sld(ielem,6)*grlst_sld(ielem,7) &
                      + grlst_sld(ielem,3)*grlst_sld(ielem,4)*grlst_sld(ielem,8) &
                      - grlst_sld(ielem,3)*grlst_sld(ielem,5)*grlst_sld(ielem,7) &
                      - grlst_sld(ielem,2)*grlst_sld(ielem,4)*grlst_sld(ielem,9) &
                      - grlst_sld(ielem,1)*grlst_sld(ielem,6)*grlst_sld(ielem,8)
              endif
           end if
           if( postp(1) % npp_setse(5) /= 0 ) then
              !
              ! Component 11 of green lagrange stress tensor
              !
              setgrl_11 = setgrl_11 + grlst_sld(ielem,1)
           end if
           if( postp(1) % npp_setse(6) /= 0 ) then
              !
              ! Component 12 of green lagrange stress tensor
              !
              setgrl_12 = setgrl_12 + grlst_sld(ielem,2)
           end if
           if( postp(1) % npp_setse(7) /= 0 ) then
              !
              ! Component 13 of green lagrange stress tensor
              !
              setgrl_13 = setgrl_13 + grlst_sld(ielem,3)
           end if
           if( postp(1) % npp_setse(8) /= 0 ) then
              !
              ! Component 21 of green lagrange stress tensor
              !
              setgrl_21 = setgrl_21 + grlst_sld(ielem,4)
           end if
           if( postp(1) % npp_setse(9) /= 0 ) then
              !
              ! Component 22 of green lagrange stress tensor
              !
              setgrl_22 = setgrl_22 + grlst_sld(ielem,5)
           end if
           if( postp(1) % npp_setse(10) /= 0 ) then
              !
              ! Component 23 of green lagrange stress tensor
              !
              setgrl_23 = setgrl_23 + grlst_sld(ielem,6)
           end if
           if( postp(1) % npp_setse(11) /= 0 ) then
              !
              ! Component 31 of green lagrange stress tensor
              !
              setgrl_31 = setgrl_31 + grlst_sld(ielem,7)
           end if
           if( postp(1) % npp_setse(12) /= 0 ) then
              !
              ! Component 32 of green lagrange stress tensor
              !
              setgrl_32 = setgrl_32 + grlst_sld(ielem,8)
           end if
           if( postp(1) % npp_setse(13) /= 0 ) then
              !
              ! Component 33 of green lagrange stress tensor
              !
              setgrl_33 = setgrl_33 + grlst_sld(ielem,9)
           end if
        end if
     end if

  end do elements

end subroutine sld_elmset
