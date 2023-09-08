subroutine rad_bouope()
  !------------------------------------------------------------------------
  !****f* Radiat/rad_bouope
  ! NAME 
  !    rad_bouope
  ! DESCRIPTION
  !    Assemble boundary elements
  ! USES
  ! USED BY
  !    rad_matrix 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  implicit none
  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  real(rp)    :: baloc(ndime,ndime)
  real(rp)    :: elvel(ndime,mnode),gprad(mgaus)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: bocod(ndime,mnodb)
  real(rp)    :: botem(mnodb)
  integer(ip) :: ielem,inode,ipoin,kfl_gobou
  integer(ip) :: igaus,igaub,iboun,inodb,pblty,idime
  integer(ip) :: pnodb,pmate,pnode,pelty,pgaus
  real(rp)    :: eucta,tmatr,gbsur,gpdet,adotn
  real(rp)    :: gbsph,gbden,gbcon,gbtem,gbvel(3)
  real(rp)    :: gpsph(mgaus),gpden(mgaus),gpdif(mgaus)
  real(rp)    :: gpabs(mgaus),gpbbr(mgaus),gptem(mgaus)
  real(rp)    :: arobi,trobi,qrobi
  real(rp)    :: emisi
  real(rp)    :: xmrhs,xmmat
  real(rp)    :: eledd(mnode),gpcar(ndime,mnode,mgaus),gpgrd(ndime,mgaus) 
  real(rp)    :: gpcon(mgaus),dummr(ndime*mnode),gpcod
  real(rp)    :: xjaci(9),xjacm(9)
  !
  ! Loop over elements  
  !
  boundaries: do iboun=1,nboun

     kfl_gobou = 1

     pblty = ltypb(iboun)
     pnodb = nnode(pblty)
     ielem = lelbo(iboun)
     pelty = ltype(ielem)

     if( pelty > 0 ) then

        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        !
        ! Inititalize
        !
        elmat = 0.0_rp
        elrhs = 0.0_rp
        !
        ! Gather at boundary nodes
        !
        bocod(1:ndime,1:pnodb) = coord(1:ndime,lnodb(1:pnodb,iboun))
        elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
        botem(1:pnodb)         = tempe_rad(lnodb(1:pnodb,iboun))
        !
        ! GPTEM+GPSGS: Radiation at Gauss point
        !
        call gather(&
             1_ip,pgaus,pnode,1_ip,lnods(1,ielem),&
             elmar(pelty)%shape,radav_rad,gprad)
        if(kfl_sgsno_rad==1)&
             call rad_sgsope(&
             2_ip,ielem,pgaus,rasgs_rad(ielem)%a,gprad)
        !
        ! Cartesian derivatives
        !
        do igaus=1,pgaus
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&      ! Cartesian derivative
                elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)        ! and Jacobian
        end do
        !
        ! Properties 
        !
!!$        call rad_elmpro(&
!!$             ielem,pmate,pnode,pgaus,1_ip,pgaus,&
!!$             elmar(pelty)%shape,gpcar,gpdif,gpabs,gpbbr,gptem,&  
!!$             gpgrd)          
        !
        ! Emissivity epsilon = 1 - reflectivity (for gray walls)
        !
        emisi = 0.75_rp !bvnat_rad(1,iboun,1)   !!F  BUG HERE BVNAT IS ALWAYS ZERO
        if (emisi <= 1.0_rp) then ! Transform emisi into correct prefactor for boundary
           emisi = 2.0_rp * emisi/ (2.0_rp - emisi)
        else
           call runend('RADIAT: Emmisivity should be less than 1')
        endif
        !
        ! Loop over Gauss points
        !
        gauss_points: do igaub=1,ngaus(pblty)
           !
           ! Jacobian EUCTA
           !
           call bouder(&
                pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&
                bocod,baloc,eucta)
           gbsur=elmar(pblty)%weigp(igaub)*eucta 
           call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
           !
           ! Cylindrical coordinates
           !
           if(kfl_naxis==1) then
              gpcod=0.0_rp
              do inodb=1,pnodb
                 gpcod=gpcod+bocod(1,inodb)*elmar(pblty)%shape(inodb,igaub)
              end do
              gbsur=gbsur*gpcod*twopi
           end if

           ! Convert temperature^4 at boundary nodes
           gbtem = 0.0_rp
           do inodb=1,pnodb           
              gbtem=gbtem+botem(inodb) &
                   *elmar(pblty)%shape(inodb,igaub)
           end do
           qrobi = -emisi * steph_rad * gbtem**4   ! qr is (2 epsilon/(2-epsilon) sigma T^4
           arobi = emisi/4.0_rp                          ! ar
           xmrhs = -qrobi                          
           xmmat = arobi                         

           call rad_boumat(&  
                pnode,pnodb,lboel(1,iboun),xmmat,xmrhs,&
                elmar(pblty)%shape(1,igaub),gbsur,elmat,elrhs)
        end do gauss_points

        !
        ! Prescribe Dirichlet boundary conditions
        !
!!F        call rad_elmdir(&
!!F             pnode,lnods(1,ielem),elmat,elrhs)
        !
        ! Assembly
        !
        call assrhs(solve(1)%ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
        call assmat(&
             solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
             ielem,lnods(1,ielem),elmat,amatr)

     end if

  end do boundaries

end subroutine rad_bouope
