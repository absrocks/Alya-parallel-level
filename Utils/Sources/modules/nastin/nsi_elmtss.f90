subroutine nsi_elmtss(&
     pelty,pmate,pnode,lnods,ielem,elcod,elvel,&
     gpcar,chale,hleng,dtcri)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmtss
  ! NAME
  !    nsi_elmtss
  ! DESCRIPTION
  !    This routine computes the element time step
  ! USED BY
  !    nsi_updtss
  !***
  !-----------------------------------------------------------------------
  use def_kintyp,     only   :  ip,rp
  use def_master,     only   :  densi,visco,wmean,kfl_coupl,&
       &                        ID_NASTIN,ID_CHEMIC,kfl_paral
  use def_domain,     only   :  mnode,ndime,elmar
  use def_nastin,     only   :  ncoef_nsi,kfl_taust_nsi,&
       &                        staco_nsi,gamth_nsi,&
       &                        kfl_cotur_nsi,kfl_advec_nsi,&
       &                        turbu_nsi,corio_nsi,kfl_regim_nsi
  use def_kermod,     only   :  kfl_prope
  use mod_ker_proper, only   :  ker_proper,kfl_kemod_ker
  use mod_tauadr, only       :  tauadr

  implicit none
  integer(ip), intent(in)    :: pelty,pmate,pnode
  integer(ip), intent(in)    :: lnods(pnode),ielem
  real(rp),    intent(in)    :: elcod(ndime,mnode)
  real(rp),    intent(in)    :: elvel(ndime,mnode)
  real(rp),    intent(in)    :: chale(ndime),hleng(ndime)
  real(rp),    intent(inout) :: dtcri
  real(rp),    intent(out)   :: gpcar(ndime,mnode)
  integer(ip)                :: idime,jdime,inode, taust
  real(rp)                   :: xjaci(9),xjacm(9)
  real(rp)                   :: gpdet,gpvis(1),gpden(1),gppor(1),gpgvi(3)
  real(rp)                   :: gpvel(3),gplev,gpvno,adv,dif,rea
  real(rp)                   :: gpgve(ndime,ndime),gppre,gptem,rnode
  real(rp)                   :: elmut(mnode)
  real(rp)                   :: gpmut(1),grvis(ndime,ndime)
  real(rp)                   :: fake_dtinv
  !
  ! Velocity
  !
  gplev = 0.0_rp
  gppre = 0.0_rp
  gpvel = 0.0_rp
  gpgve = 0.0_rp
  gptem = 0.0_rp
  gpmut = 0.0_rp
  gpvno = 0.0_rp
  dtcri = 0.0_rp
  rnode = 1.0_rp/real(pnode,rp)
  !
  ! GPCAR: Cartesian derivative
  !
  call elmder(&
       pnode,ndime,elmar(pelty)%dercg,elcod,&
       gpcar,gpdet,xjacm,xjaci)
  !
  ! Values at center of gravity
  !
  if( kfl_advec_nsi /= 0 ) then                        ! GPVEL: Velocity
     do inode = 1,pnode
        gpvel(1:ndime) = gpvel(1:ndime) + elvel(1:ndime,inode)
     end do
     gpvel = rnode*gpvel
     gpvno = sqrt(dot_product(gpvel(1:ndime),gpvel(1:ndime)))
  end if
  !
  ! GPVIS, GPDEN, GPPOR: Compute properties
  !
  call ker_proper('DENSI','COG  ',1_ip,ielem,gpden)
  call ker_proper('VISCO','COG  ',1_ip,ielem,gpvis)
  call ker_proper('POROS','COG  ',1_ip,ielem,gppor)
  call ker_proper('TURBU','COG  ',1_ip,ielem,gpmut)

  !gpden = 1.0_rp
  !gpvis = 1.0_rp
  !gppor = 0.0_rp
  !gpmut = 0.0_rp
  !return
  
  if( kfl_cotur_nsi /= 0 ) then
     gpgvi = 0.0_rp
     grvis = 0.0_rp
     call nsi_turbul(&
          1_ip,0_ip,pnode,1_ip,1_ip,1_ip,kfl_cotur_nsi,&
          elmar(pelty)%shacg,gpcar,elvel,gpden(1),gpvis(1),gpmut,&
          gpgvi,grvis,gpgve,ielem,kfl_kemod_ker)
  end if

  !
  ! DTCRI: Critical time step
  !
  adv        = gpden(1) * gpvno                      ! Convective term
  dif        = gpvis(1)                              ! Viscous term
  rea        = gpden(1) * corio_nsi + abs(gppor(1))  ! Porosity + Coriolis term

  ! When tau depends on time step. Let us use codina
  if (kfl_taust_nsi.gt.4) then
     taust=1
  else
     taust = kfl_taust_nsi
  endif
     
  
  call tauadr(&
       taust,staco_nsi,adv,dif,rea,&
       chale(1),chale(2),dtcri) 

  if( adv == 0.0_rp .and. dif == 0.0_rp .and. rea == 0.0_rp ) dtcri = huge(1.0_rp)
  
  dtcri = gpden(1)*dtcri

end subroutine nsi_elmtss
