subroutine hlm_bouope()
  !------------------------------------------------------------------------
  !****f* Temper/hlm_bouope
  ! NAME 
  !    hlm_bouope
  ! DESCRIPTION
  !    ORDER=1:
  !      Temperature equation, boundary operations
  ! USES
  ! USED BY
  !    hlm_matrix 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_helmoz
  implicit none
  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  real(rp)    :: baloc(ndime,ndime)
  real(rp)    :: elvel(ndime,mnode),gptem(mgaus)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: bocod(ndime,mnodb),botem(mnodb)
  real(rp)    :: bovel(ndime,mnodb)
  integer(ip) :: ielem,inode,ipoin,kfl_gobou
  integer(ip) :: igaus,igaub,iboun,inodb,pblty,idime
  integer(ip) :: pnodb,pmate,pnode,pelty,pgaus
  real(rp)    :: eucta,tmatr,gbsur,gpdet,adotn
  real(rp)    :: gbsph,gbden,gbcon,gbtem,gbvel(3)
  real(rp)    :: gpsph(mgaus),gpden(mgaus),gpdif(mgaus)
  real(rp)    :: gprea(mgaus)
  real(rp)    :: arobi,trobi,qrobi,twall,acvis
  real(rp)    :: para1,para2,para3,para4
  real(rp)    :: xmrhs,xmmat
  real(rp)    :: eledd(mnode),gpcar(ndime,mnode,mgaus)
  real(rp)    :: gpcon(mgaus),dummr(ndime*mnode),gpcod
  real(rp)    :: xjaci(9),xjacm(9)
  !
  ! Loop over elements  
  !
  boundaries: do iboun=1,nboun

        pblty = ltypb(iboun)
        pnodb = nnode(pblty)
        ielem = lelbo(iboun)
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        pmate = 1
        if( nmate > 1 ) pmate = lmate(ielem)
        !
        ! Inititalize
        !
        elmat = 0.0_rp
        elrhs = 0.0_rp
        !
        ! Gather operations
        !
        bocod(1:ndime,1:pnodb) = coord(1:ndime,lnodb(1:pnodb,iboun))
        elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
        !
        ! GPTEM+GPSGS: Temperature at Gauss point
        !
        !call gather(&
        !     1_ip,pgaus,pnode,1_ip,lnods(1,ielem),&
        !     elmar(pelty)%shape,tempe,gptem)
        !
        ! Cartesian derivatives
        !
        do igaus = 1,pgaus
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&      ! Cartesian derivative
                elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)        ! and Jacobian
        end do
        !
        ! Loop over Gauss points
        !
        gauss_points: do igaub = 1,ngaus(pblty)
           !
           ! Jacobian EUCTA
           !
           call bouder(&
                pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&
                bocod,baloc,eucta)
           gbsur=elmar(pblty)%weigp(igaub)*eucta 
           call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
 
           !call hlm_boumat(&
           !     pnode,pnodb,lboel(1,iboun),xmmat,xmrhs,&
           !     elmar(pblty)%shape(1,igaub),gbsur,elmat,elrhs)

        end do gauss_points
        !
        ! Assembly
        !
        !call hlm_assemb(1_ip,pnode,lnods(1,ielem),elmat,elrhs)


  end do boundaries

end subroutine hlm_bouope
