subroutine rad_bouset(ibset,setsu,setmt,sethf)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_bouset
  ! NAME 
  !    rad_bouset
  ! DESCRIPTION
  !    This routine computes variables on a boundary set W.
  !    The variable are:
  !    1. setsu: set surface          =  meas(W)=int_W
  !    2. setmt: set mean temprature =  int_W T/meas(W)
  !    3. sethf: set heat flux        =  int_W k*grad(T).n
  !    The heat flux represent the flo from the solid to the fluid.
  !    Positive value means incomming positive heat.
  ! USES
  !    bouder
  !    chenor
  ! USED BY
  !    rad_outset
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  implicit none

  integer(ip), intent(in)  :: ibset
  real(rp),    intent(out) :: setsu,setmt,sethf
  real(rp)                 :: baloc(ndime,ndime)
  real(rp)                 :: bocod(ndime,mnodb)
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: gbcar(ndime,mnode)
  real(rp)                 :: xjaci(ndime,ndime) 
  real(rp)                 :: xjacm(ndime,ndime) 
  integer(ip)              :: ielem,inode,ipoin
  integer(ip)              :: igaus,idime,igaub,iboun,inodb,pblty
  integer(ip)              :: pnodb,pmate,pnode,pelty,pgaus,pgaub
  real(rp)                 :: eucta,tmatr,dsurf,detjm
  real(rp)                 :: gbsph,gbden,gbdif
  real(rp)                 :: gpsph(mgaus),gpden(mgaus)
  real(rp)                 :: gpdif(mgaus),gpcon(mgaus)
  real(rp)                 :: gpabs(mgaus),gprad(mgaus)
  real(rp)                 :: gbrad,gbgrt(ndime)
  real(rp)                 :: eledd(mnode),gpcar(ndime,mnode,mgaus)
  real(rp)                 :: dummr(ndime,mnode)
  real(rp)                 :: gpbbr(mgaus),gptem(mgaus),gpgrd(ndime,mgaus) 

  gpcar=0.0_rp
  setsu=0.0_rp
  setmt=0.0_rp
  sethf=0.0_rp
  !
  ! Loop over elements
  !
  boundaries: do iboun=1,nboun

     if(lbset(iboun)==ibset) then

        pblty=ltypb(iboun)
        pnodb=nnode(pblty)
        pgaub=ngaus(pblty)
        ielem=lelbo(iboun)
        pelty=ltype(ielem)
        pnode=nnode(pelty)
        pgaus=ngaus(pelty)
        !
        ! Gather operations
        !
        do inodb=1,pnodb
           ipoin=lnodb(inodb,iboun)
           do idime=1,ndime
              bocod(idime,inodb)=coord(idime,ipoin)
           end do
        end do
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
        end do

        !
        ! Cartesian derivatives
        !
        do igaus=1,pgaus
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),& 
                elcod,gpcar(1,1,igaus),detjm,xjacm,xjaci) 
        end do
        !
        ! GPRAD: Radiation at Gauss point
        !
        call gather(&
             1_ip,pgaus,pnode,1_ip,lnods(1,ielem),&
             elmar(pelty)%shape,radav_rad,gprad)
        if(kfl_sgsno_rad==1)&
             call rad_sgsope(&
             2_ip,ielem,pgaus,rasgs_rad(ielem)%a,gprad)
        !
        ! Properties
        !
        call rad_elmpro(&
             ielem,pmate,pnode,pgaus,1_ip,pgaus,&
             elmar(pelty)%shape,gpcar,gpdif,gpabs,gpbbr,gptem,gpgrd)

        gauss_points: do igaub=1,pgaub
           !
           ! Properties
           !
           gbden=0.0_rp
           gbsph=0.0_rp
           gbdif=0.0_rp
           do igaus=1,pgaus
              do inodb=1,pnodb                  
                 tmatr=elmar(pelty)%shaga(igaus,lboel(inodb,iboun))&
                      *elmar(pblty)%shape(inodb,igaub)
                 gbden=gbden+gpden(igaus)*tmatr
                 gbsph=gbsph+gpsph(igaus)*tmatr
                 gbdif=gbdif+gpdif(igaus)*tmatr
              end do
           end do
           
           ! Jacobian
           call bouder(&
                pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&
                bocod,baloc,eucta)
           call chenor(pnode,baloc,bocod,elcod)
           dsurf=elmar(pblty)%weigp(igaub)*eucta 
           setsu=setsu+dsurf

           if(postp(1)%npp_setsb(1)/=0) then
              !
              ! Mean radiation intensity
              !
              gbrad=0.0_rp 
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun)
                 gbrad=gbrad+elmar(pblty)%shape(inodb,igaub)*radav_rad(ipoin,1)
              end do
              setmt=setmt+gbrad*dsurf
           end if

           if(postp(1)%npp_setsb(2)/=0) then
              !
              ! Heat flux
              !
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 elcod(1:ndime,inode)=coord(1:ndime,ipoin)
              end do
              !
              ! GBCAR: Cartesian derivatives at boundary Gauss point
              !
              call gpcabo(&
                   pnode,pgaus,pnodb,lboel(1,iboun),elmar(pelty)%shaga,&
                   gpcar,elmar(pblty)%shape(1,igaub),gbcar)
              gbgrt=0.0_rp
              do inode=1,pnode
                 ipoin=lnods(inode,ielem) 
                 do idime=1,ndime
                    gbgrt(idime)=gbgrt(idime)+gbcar(idime,inode)*radav_rad(ipoin,1)
                 end do
              end do
              sethf=sethf+gbdif*dsurf&
                   *dot_product(baloc(1:ndime,ndime),gbgrt(1:ndime))

           end if

        end do gauss_points

     end if

  end do boundaries

 
end subroutine rad_bouset
