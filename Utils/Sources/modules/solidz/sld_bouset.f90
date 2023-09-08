!------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_bouset.f90
!> @author  Mariano Vazquez
!> @date    April, 2010
!>          - Subroutine creation
!> @author  Gerard Guillamet
!>          August, 2018
!>          - Subroutine refactoring
!> @brief   This routine computes variables on a boundary set W.
!> @details
!>
!>          This routine computes variables on a boundary set W.
!>          The variable are:
!>          0        SETSU: surface                         =  int_W  = meas(W)
!>          1  -> 3  SETFO: Force                           =
!>          4        SETNF: Normal force                    =
!>          5  -> 7  SETDI: Averaged disp.                  =  Sum_{ui in W}/nodb (u)|_i
!>          8        SETDM: Averaged disp. (Magnitude)      =  sqrt(ux^2 + uy^2 + uz^2)
!>          9  -> 11 SETFR: Sum. Reaction force             =  Sum_{i in W} (Fext - Fint - M·a)|_i
!>          12       SETFM: Sum. reaction force (Magnitude) =  sqrt(fr_x^2 + fr_y^2 + fr_z^2)
!>          13 -> 15 SETCF: Sum. contact force              =  Sum(fc)
!>          16       SETFC: Sum. contact force (Magnitude)  =  sqrt(fc_x^2 + fc_y^2 + fc_z^2)
!>
!>          Values of SETDI and SETDM are averaged further on in sld_outset
!>
!>@warning  <GGU> SETFO and SETNF requires a revision
!>
!> @}
!------------------------------------------------------------------------------

subroutine sld_bouset(ibsec,ibset)

  use def_kintyp,         only : ip, rp
  use def_master,         only : ITER_K
  use def_master,         only : postp, displ, gisca, vbset
  use def_domain,         only : ndime, npoin, nboun, ndimb
  use def_domain,         only : elmar, coord
  use def_domain,         only : mnodb, mnode, mgaus, mgaub
  use def_domain,         only : lnodb, lboel, ltype, lnods, lnnob, ltypb, lbset, lelbo
  use def_domain,         only : nnode, ngaus
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use def_solidz,         only : gppio_sld
  use def_solidz,         only : kfl_fixno_sld, frxid_sld
  use def_solidz,         only : kfl_conta_sld, fcont_sld

  implicit none

  integer(ip), intent(in)  :: ibsec         !< Boundary number
  integer(ip), intent(in)  :: ibset         !< Set code
  real(rp),    pointer     :: setsu(:),setfo(:),setnf(:),setfr(:),setfm(:),setdi(:),setdm(:)
  real(rp),    pointer     :: setfc(:),setcf(:)
  integer(ip)              :: ielem,inode,ipoin,igaus,idime,jdime,nn
  integer(ip)              :: pnode,pgaus,iboun,igaub,inodb
  integer(ip)              :: pelty,pblty,pnodb,pgaub
  real(rp)                 :: baloc(ndime,ndime)
  real(rp)                 :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)                 :: bodis(ndime,mnodb),eldis(ndime,mnode)
  real(rp)                 :: gppio(ndime,ndime,mgaus)              ! Piola tensor P
  real(rp)                 :: gbsur,eucta
  real(rp)                 :: gbdis(ndime,mgaub)
  real(rp)                 :: dummr(3),tract(3),tmatr
  real(rp)                 :: frxid_aux(ndime,npoin),fcont_aux(ndime,npoin)

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  nn    =  postp(1) % nvabs + 1
  setsu => vbset( nn:nn, ibset ) ! Surface
  setfo => vbset(  1:3 , ibset ) ! Force
  setnf => vbset(  4:4 , ibset ) ! Normal Force
  setdi => vbset(  5:7 , ibset ) ! Averaged displacement
  setdm => vbset(  8:8 , ibset ) ! Averaged displacement (Magnitude)
  setfr => vbset(  9:11, ibset ) ! Sum of reaction forces
  setfm => vbset( 12:12, ibset ) ! Sum of reaction forces (Magnitude)
  setcf => vbset( 13:15, ibset ) ! Sum of contact forces
  setfc => vbset( 16:16, ibset ) ! Sum of contact forces (Magnitude)
  setsu =  0.0_rp
  setfo =  0.0_rp
  setnf =  0.0_rp
  setdi =  0.0_rp
  setdm =  0.0_rp
  setfr =  0.0_rp
  setfm =  0.0_rp
  setcf =  0.0_rp

  boundaries: do iboun = 1,nboun

     if ( lbset(iboun) == ibsec ) then

        !----------------------------------------------------------------
        !
        ! Element properties, dimensions and gather
        !
        !----------------------------------------------------------------

        pblty = ltypb(iboun)
        pnodb = nnode(pblty)
        pgaub = ngaus(pblty)

        do inodb = 1,pnodb
           ipoin = lnodb(inodb,iboun)
           bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
           bodis(1:ndime,inodb) = displ(1:ndime,ipoin,ITER_K)
        end do

        ielem = lelbo(iboun)
        pelty = ltype(ielem)
        if ( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
              eldis(1:ndime,inode) = displ(1:ndime,ipoin,ITER_K)
           end do

           !----------------------------------------------------------------
           !
           ! Values at Gauss points:
           !
           !----------------------------------------------------------------

           gppio(:,:,:) = 0.0_rp
           gbdis(:,:) = 0.0_rp
           do igaub = 1,pgaub
              do inodb = 1,pnodb
                 do idime = 1,ndime
                    gbdis(idime,igaub) = gbdis(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * bodis(idime,inodb)
                 end do
              end do
           end do

           !----------------------------------------------------------------
           !
           ! Loop over Gauss points
           !
           !----------------------------------------------------------------

           gauss_points: do igaub = 1,pgaub

              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&  ! Cartesian derivative
                   bocod,baloc,eucta)                                 ! and Jacobian
              gbsur = elmar(pblty)%weigp(igaub)*eucta
              setsu = setsu + gbsur
              call chenor(pnode,baloc,bocod,elcod)                    ! Check normal

              !-------------------------------------------------------------
              !
              ! Compute traction
              ! t = P . N  (1st piola . normal in the reference system)
              !
              !-------------------------------------------------------------

              if ( postp(1) % npp_setsb(1) /= 0 .or. postp(1) % npp_setsb(4) /= 0 ) then
                 !
                 ! Extrapolate from Gauss points igauss to node inode and
                 ! then interpolation at boundary Gauss points igaub
                 !
                 do igaus = 1,pgaus
                    do inodb = 1,pnodb
                       tmatr =  elmar(pelty)%shaga(igaus,lboel(inodb,iboun))&
                            & * elmar(pblty)%shape(inodb,igaub)
                       do jdime = 1,ndime
                          do idime = 1,ndime
                             gppio(idime,jdime,igaub) = gppio(idime,jdime,igaub) &
                                  + gppio_sld(ielem) % a(idime,jdime,igaus)*tmatr
                          end do
                       end do
                    end do
                 end do

                 do idime = 1,ndime
                    tract(idime)= 0.0_rp
                    do jdime = 1,ndime
                       tract(idime) = tract(idime) &
                            + gppio(jdime,idime,igaub) * baloc(jdime,ndime)
                    end do
                 end do
                 !
                 ! FORCE, F_y, F_z
                 !
                 if ( postp(1) % npp_setsb(1) /= 0  ) then
                    ! OPTION 1 : (x y z) component of the force
                    setfo(1) = setfo(1) + gbsur * tract(1)
                    setfo(2) = setfo(2) + gbsur * tract(2)
                    setfo(3) = setfo(3) + gbsur * tract(3)
                 end if
                 !
                 ! NORMA (Normal force)
                 !
                 if ( postp(1) % npp_setsb(4) /= 0  ) then
                    dummr(1) = 0.0_rp
                    do idime = 1,ndime
                       dummr(1) = dummr(1) + tract(idime) * baloc(idime,ndime)
                    end do
                    setnf(1) = setnf(1) + dummr(1) * gbsur
                 end if

              end if

              !-------------------------------------------------------------
              !
              ! Averaged displacement
              !
              !-------------------------------------------------------------

              if( any(postp(1) % npp_setsb(5:8) /= 0_ip) ) then
                 !
                 ! DIBOX, DIBOY, DIBOZ
                 !
                 setdi(1:ndime) = setdi(1:ndime) + gbsur * gbdis(1:ndime,igaub)
                 !
                 ! DIBOU (Magnitude)
                 !
                 if ( postp(1) % npp_setsb(8) /= 0 ) setdm = sqrt(sum(setdi(1:ndime)**2))

              end if

           end do gauss_points

        end if

     end if

  end do boundaries

  !----------------------------------------------------------
  !
  ! Sum of the Reaction and Contact forces
  !
  !----------------------------------------------------------

  if( any(postp(1) % npp_setsb(9:16) /= 0_ip) ) then

     call memgen(1_ip,npoin,0_ip)
     do iboun = 1,nboun
        if ( lbset(iboun) == ibset ) then
           do inodb = 1,lnnob(iboun)
              ipoin = lnodb(inodb,iboun)
              gisca(ipoin) = 1
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM','IN MY CODE')

     !
     ! FRBOX, FRBOY, FRBOZ (and FRBOU)
     !
     if( any(postp(1) % npp_setsb(9:12) /= 0_ip) ) then

        frxid_aux(1:ndime,1:npoin) = reshape(frxid_sld,shape=(/ndime,npoin/))
        do ipoin = 1,npoin
           if ( gisca(ipoin) /= 0_ip ) then
              do idime = 1,ndime
                 if ( kfl_fixno_sld(idime,ipoin) == 1_ip ) then
                    setfr(idime) = setfr(idime) - frxid_aux(idime,ipoin) / real(gisca(ipoin),rp)
                 end if
              end do
           end if
        end do
        ! FRBOU
        if ( postp(1) % npp_setsb(12) /= 0_ip ) setfm = sqrt(sum(setfr(1:ndime)**2))

     end if
     !
     ! FCONX, FCONY, FCONZ (and FCONT)
     !
     if( any(postp(1) % npp_setsb(13:16) /= 0_ip) ) then

        if ( kfl_conta_sld /= 0_ip ) then
           fcont_aux(1:ndime,1:npoin) = reshape(fcont_sld,shape=(/ndime,npoin/))
           do ipoin = 1,npoin
              if ( gisca(ipoin) /= 0_ip .and. kfl_fixno_sld(1,ipoin) == 3_ip ) then
                 setcf(1:ndime) = setcf(1:ndime) + fcont_aux(1:ndime,ipoin) / real(gisca(ipoin),rp)
              end if
           end do
           ! FCONT
           if ( postp(1) % npp_setsb(16) /= 0_ip ) setfc = sqrt(sum(setcf(1:ndime)**2))
        end if

     end if

     call memgen(3_ip,npoin,0_ip)

  end if

end subroutine sld_bouset
