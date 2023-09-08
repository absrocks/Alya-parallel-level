subroutine tem_elmset(ieset,setsu,setmc)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_elmset
  ! NAME 
  !    tem_elmset
  ! DESCRIPTION
  !    This routine computes variables on an element set W.
  !    The variable are:
  !    1. setsu: set surface          =  meas(W)=int_W
  !    2. setmt: set mean convective temperature 
  !       = int_W rho*cp*(u.t)*T / int_W rho*cp*(u.t)
  ! USES
  !    elmder
  !    tem_elmpro
  ! USED BY
  !    tem_outset
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  use mod_ker_proper 
  implicit none

  integer(ip), intent(in)  :: ieset
  real(rp),    intent(out) :: setsu,setmc

  real(rp)    :: cartd(ndime,mnode)
  real(rp)    :: xjaci(ndime,ndime) 
  real(rp)    :: xjacm(ndime,ndime) 
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: denom,udott
  integer(ip) :: ielem,inode,ipoin                     ! Indices and dimensions
  integer(ip) :: igaus,idime,dummi
  integer(ip) :: pelty,pmate,pnode,pgaus
  real(rp)    :: gpvol,gpdet                           ! Values at Gauss points
  real(rp)    :: gpcon(mgaus)                          ! k
  real(rp)    :: gpdif(mgaus)                          ! k+kt, grad(k+kt)
  real(rp)    :: gptur(mgaus)
  real(rp)    :: gpsph(mgaus),gpden(mgaus)
  real(rp)    :: gprea(mgaus),gptem(mgaus)
  real(rp)    :: eledd(mnode)
  real(rp)    :: gpcar(ndime,mnode,mgaus)
  real(rp)    :: gpvel(ndime),dummr(ndime*mnode)

  setsu = 0.0_rp 
  setmc = 0.0_rp
  denom = 0.0_rp
  !
  ! Loop over elements
  !
  elements: do ielem=1,nelem

     if(leset(ielem)==ieset) then

        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        pmate = 1
        if(nmate>1) pmate=lmate(ielem)
        !
        ! Gather operations
        !
        elcod(1:ndime,1:pnode)=coord(1:ndime,lnods(1:pnode,ielem))
        call gather(&
             1_ip,pgaus,pnode,1_ip,lnods(1,ielem),&
             elmar(pelty)%shape,tempe,gptem)
        !
        ! GPTEM+GPSGS: Temperature needed when k=k(T), Cp=Cp(T)
        !
        call gather(&
             1_ip,pgaus,pnode,1_ip,lnods(1,ielem),&
             elmar(pelty)%shape,tempe,gptem)
        !
        ! Properties
        !
        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
        call ker_proper('CONDU','PGAUS',dummi,ielem,gpcon,pnode,pgaus,elmar(pelty)%shape,gpcar)
        call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)
        call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty) % shape,gpcar)
        ! 
        ! Coupling with turbul
        !
        call tem_turbul(&
             ielem,pnode,pgaus,1_ip,pgaus,elmar(pelty)%shape,gpcon,gpsph,gpdif,dummr,gpden,gptur)
              
        !
        ! Reaction term
        !
        call tem_elmrea( &
             1_ip,pnode,pgaus,1_ip,pgaus,dummr,gpden,gpcar,&
             gprea)
       

        gauss_points: do igaus=1,pgaus 
           !
           ! Cartesian derivatives and Jacobian
           !
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                elcod,cartd,gpdet,xjacm,xjaci)
           gpvol=elmar(pelty)%weigp(igaus)*gpdet
           setsu=setsu+gpvol
           !
           ! Set calculations
           !
           if(postp(1)%npp_setse(1)/=0) then
              if(kfl_advec_tem==1) then
                 gpvel=0.0_rp
                 do inode=1,pnode
                    ipoin=lnods(inode,ielem)
                    do idime=1,ndime
                       gpvel(idime)=gpvel(idime)+veloc(idime,ipoin,1)&
                            *elmar(pelty)%shape(inode,igaus)
                    end do
                 end do
                 udott=dot_product(gpvel(1:ndime),postp(1)%paese(1:ndime,1))
                 setmc=setmc+gpden(igaus)*gpsph(igaus)*udott*gptem(igaus)*gpvol
                 denom=denom+gpden(igaus)*gpsph(igaus)*udott*gpvol
              end if
           end if

        end do gauss_points

     end if

  end do elements

end subroutine tem_elmset
