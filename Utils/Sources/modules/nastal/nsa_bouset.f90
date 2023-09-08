subroutine nsa_bouset(ibsec,ibset)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_bouset
  ! NAME 
  !    nsa_bouset
  ! DESCRIPTION
  !    This routine computes variables on a boundary set W.
  !    The variable are: 
  !     1       SETMA: mass               =  int_W  rho*u.n                 
  !     2       SETRO: density            =  int_W  rho.n          
  !     3       SETBU: bulk velocity      =  int_W  Ub = int_W  rho*u.n / int_W  rho.n                 
  !
  !    Force and torque are exerted by the solid on the fluid.
  !    Values of SETMP and SETYP are averaged further on in nsa_outset
  !
  ! USES 
  !    bouder
  !    chenor
  ! USED BY
  !    nsa_outset
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastal
  implicit none

  integer(ip), intent(in)  :: ibsec,ibset
  real(rp),    pointer     :: setsu(:),setde(:),setma(:),setbu(:)
  integer(ip)              :: ielem,inode,ipoin,idime,nn
  integer(ip)              :: pnode,pgaus,iboun,igaub,inodb,kboun
  integer(ip)              :: pelty,pblty,pnodb,pgaub,pmate
  real(rp)                 :: baloc(ndime,ndime),bopre(mnodb),boden(mnodb),botem(mnodb),bovis(mnodb)
  real(rp)                 :: bovel(ndime,mnodb),elvel(ndime,mnode)
  real(rp)                 :: bovfi(ndime,mnodb),tragl(9)     
  real(rp)                 :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)                 :: gbpre(mgaus),gbtem(mgaus),hleng(3)
  real(rp)                 :: gbsur,eucta,gbden(mgaus),gbvis(mgaus)
  real(rp)                 :: gbvel(ndime,mgaus),xfact
  real(rp)                 :: gbvdt(ndime,mgaus)            ! tangent component of velocity - prescribed velocity.

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  nn    =  postp(1) % nvabs + 1
  setsu => vbset( nn:nn , ibset ) ! Surface           
  setma => vbset(  1: 1 , ibset ) ! Mass      
  setde => vbset(  2: 2 , ibset ) ! Density
  setbu => vbset(  3: 3 , ibset ) ! Bulk Velocity
  setsu = 0.0_rp
  setma = 0.0_rp
  setde = 0.0_rp
  setbu = 0.0_rp

  boundaries: do iboun = 1,nboun

     if( lbset(iboun) == ibsec ) then

        !----------------------------------------------------------------
        !
        ! Element properties and dimensions and gather
        !
        !----------------------------------------------------------------

        pblty = ltypb(iboun) 
        pnodb = nnode(pblty)
        pgaub = ngaus(pblty)
        pmate = 1

        do inodb = 1,pnodb
           ipoin = lnodb(inodb,iboun)
           bopre(inodb) = press(ipoin,1)
           boden(inodb) = densi(ipoin,1)
           botem(inodb) = tempe(ipoin,1)
           bovis(inodb) = visco(ipoin,1)
           do idime = 1,ndime
              bocod(idime,inodb) = coord(idime,ipoin)
              bovel(idime,inodb) = veloc(idime,ipoin,1)
           end do
           if( kfl_coupl(ID_NASTAL,ID_ALEFOR) /= 0 ) then 
              do idime = 1,ndime
                 bovfi(idime,inodb) = velom(idime,ipoin)    ! see comments in nsi_bouope
              end do
           else
              do idime = 1,ndime
                 bovfi(idime,inodb) = 0.0_rp
              end do
           end if
        end do    

        ielem = lelbo(iboun)
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elvel(idime,inode) = veloc(idime,ipoin,1)
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           !
           ! Element length HLENG
           !
           call elmlen(&
                ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                hnatu(pelty),hleng)

           !----------------------------------------------------------------
           !
           ! Values at Gauss points: GBPRE, GBVEL, GBTEM, GBVIS
           !
           !----------------------------------------------------------------
           do igaub = 1,pgaub
              gbpre(igaub) = 0.0_rp
              gbden(igaub) = 0.0_rp              
              gbtem(igaub) = 0.0_rp
              gbvis(igaub) = 0.0_rp

              do idime = 1,ndime
                 gbvel(idime,igaub) = 0.0_rp
                 gbvdt(idime,igaub) = 0.0_rp
              end do

              do inodb = 1,pnodb
                 gbpre(igaub) = gbpre(igaub) + elmar(pblty)%shape(inodb,igaub) * bopre(inodb)
                 gbden(igaub) = gbden(igaub) + elmar(pblty)%shape(inodb,igaub) * boden(inodb)                 
                 gbtem(igaub) = gbtem(igaub) + elmar(pblty)%shape(inodb,igaub) * botem(inodb)
                 gbvis(igaub) = gbvis(igaub) + elmar(pblty)%shape(inodb,igaub) * bovis(inodb)

                 do idime = 1,ndime
                    gbvel(idime,igaub) = gbvel(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * bovel(idime,inodb)
                    gbvdt(idime,igaub) = gbvdt(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * ( bovel(idime,inodb)  &
                                    &    - bovfi(idime,inodb) )  ! for the moment it includes normal component (substracted later)
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
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
                   bocod,baloc,eucta)                                   ! and Jacobian
              gbsur = elmar(pblty)%weigp(igaub)*eucta 
              setsu = setsu + gbsur
              call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
              
              !-------------------------------------------------------------
              !
              ! Mass
              !
              !-------------------------------------------------------------

              if( postp(1) % npp_setsb(1) /= 0 ) then
                 xfact = gbden(igaub) * gbsur
                 do idime = 1,ndime
                    setma = setma + xfact * gbvel(idime,igaub) * baloc(idime,ndime)
                 end do
              end if

              !-------------------------------------------------------------
              !
              ! Density
              !
              !-------------------------------------------------------------

              if( postp(1) % npp_setsb(2) /= 0 ) then
                 setde = setde + gbden(igaub) * gbsur
              end if

           end do gauss_points

        end if

     end if

  end do boundaries

!!DMM  call runend('NSA_BOUSET: TO BE DONE')

end subroutine nsa_bouset
