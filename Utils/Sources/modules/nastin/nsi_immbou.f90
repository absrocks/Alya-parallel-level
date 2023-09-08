subroutine nsi_immbou()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_immbou
  ! NAME
  !    nsi_immbou
  ! DESCRIPTION
  !    Compute force on particles
  ! USES
  !    bouder
  !    chenor
  ! USED BY
  !    nsi_outset
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_gradie
  use mod_ker_proper
  use mod_messages, only : livinf
  implicit none
  real(rp)             :: elvel(ndime,mnode)
  real(rp)             :: elpre(mnode)
  real(rp)             :: elcod(ndime,mnode)
  real(rp)             :: elmut(mnode)
  real(rp)             :: elfle(mnode)
  real(rp)             :: eltem(mnode)
  real(rp)             :: elgra(ntens,mnode)

  real(rp)             :: bocod(ndime,mnodb)
  real(rp)             :: bovel(ndime,mnodb)
  real(rp)             :: bovfi(ndime,mnodb)     
  real(rp)             :: bopre(mnodb)
  real(rp)             :: baloc(ndime,ndime)

  real(rp)             :: gbmut
  real(rp)             :: gblev(mgaus)
  real(rp)             :: gbgve(9)
  real(rp)             :: gbrvi(3)
  real(rp)             :: gbgvi
  real(rp)             :: gbden(mgaus)
  real(rp)             :: gbvis(mgaus)
  real(rp)             :: gbpor(mgaus)
  real(rp)             :: gbpre(mgaus)
  real(rp)             :: gbtem(mgaus)
  real(rp)             :: gbvel(ndime,mgaus)
  real(rp)             :: gbvdt(ndime,mgaus)   ! tangent component of velocity - prescribed velocity.
  real(rp)             :: gbcoo(3)
  real(rp)             :: grave(3,3)
  real(rp)             :: cartb(ndime,mnode)
  real(rp)             :: gpgvi(ndime,mgaus)
  real(rp)             :: gpcar(ndime,mnode,mgaus)

  integer(ip)          :: ielem,inode,ipoin,idime,iimbo,imeth,kgaus
  integer(ip)          :: pnode,iboun,igaub,inodb,dummi,itens,jelem
  integer(ip)          :: pelty,pblty,pnodb,pgaub,pmate,igaus,pgaus
  integer(ip)          :: jdime,kboun,izdom,jzdom,jpoin,kauxi
  real(rp)             :: xjaci(9),xjacm(9),hleng(3),velno,tauwa,ustar
  real(rp)             :: gbsur,eucta,gbdet,dummr
  real(rp)             :: tragl(9),F1v(3),F1p(3),T1p(3),T1v(3),x(3)  
  real(rp)             :: gbstr(ntens),bouno(3),gbpos(ndime)
  real(rp)             :: nx,ny,nz
  real(rp)             :: deriv(3,64),coloc(3),shaib(mnode)
  real(rp)             :: veaux(3),fauxi(3),fauxn
  real(rp),    pointer :: Fp(:),Fv(:),Tp(:),Tv(:)
  !
  ! New
  !
  integer(ip)          :: iperi,islav,nslav,pinfo
  integer(ip)          :: lslav(npoin_2)

  if( nimbo == 0 ) return

  call livinf(70_ip,' ',0_ip)
  !
  ! Initialization
  !
  imeth = 1
  do iimbo = 1,nimbo
     Fv => imbou(iimbo) % vforce
     Fp => imbou(iimbo) % pforce
     Tv => imbou(iimbo) % vtorqu
     Tp => imbou(iimbo) % ptorqu
     do idime = 1,3
        Fp(idime)  = 0.0_rp
        Fv(idime)  = 0.0_rp
        Tp(idime)  = 0.0_rp
        Tv(idime)  = 0.0_rp

        F1v(idime) = 0.0_rp
        F1p(idime) = 0.0_rp

        T1p(idime) = 0.0_rp
        T1v(idime) = 0.0_rp

        x(idime)   = 0.0_rp
     end do
  end do
  if (ittim <= 1 ) return
  !
  ! Allocate memory: If gradients are smoothed
  !
  if( INOTMASTER ) then
     !
     ! Loop over boundaries
     !
     do iimbo = 1,nimbo

        Fv     => imbou(iimbo) % vforce
        Fp     => imbou(iimbo) % pforce
        Tv     => imbou(iimbo) % vtorqu
        Tp     => imbou(iimbo) % ptorqu

        !-------------------------------------------------------------
        !
        ! Embedded bodies
        !
        !-------------------------------------------------------------

        do ipoin = 1,npoin
           if( lntib(ipoin) == -iimbo) then
              do idime = 1,ndime
                 F1v(idime) = 0.0_rp
                 F1p(idime) = 0.0_rp 
              end do
              do idime = 1,ndime
                 F1v(idime) = F1v(idime) + intfo_nsi(ipoin) % bu(idime)
                 jzdom = 0
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    jpoin = c_dom(izdom)
                    jzdom = jzdom + 1
                    do jdime = 1,ndime
                       F1v(idime) = F1v(idime) - intfo_nsi(ipoin) % Auu(jdime,idime,jzdom) * veloc(jdime,jpoin,1) 
                    end do
                    F1p(idime) = F1p(idime) - intfo_nsi(ipoin) % Aup(idime,jzdom) * press(jpoin,1)                           
                 end do
              end do

              x(3) = 0.0_rp   ! so that it has the correct value in the 2d case
              do idime = 1,ndime
                 x(idime) = coord(idime,ipoin) - imbou(iimbo) % posil(idime,1) 
              end do
              call vecpro(x,F1p,T1p,3_ip)                 ! T1 = (X-Xg) x F1  (pressure) 
              call vecpro(x,F1v,T1v,3_ip)                 ! T1 = (X-Xg) x F1  (viscous) 
              !
              ! F, T: Actualize force and torque 
              !
              do idime = 1,3
                 Fv(idime) = Fv(idime) + F1v(idime)
                 Tv(idime) = Tv(idime) + T1v(idime)
                 Fp(idime) = Fp(idime) + F1p(idime)
                 Tp(idime) = Tp(idime) + T1p(idime)
              end do
           end if
        end do
     end do
     
  end if
  !
  ! Reduce sum in Parallel
  !
  call ibmdef(6_ip) 

end subroutine nsi_immbou

