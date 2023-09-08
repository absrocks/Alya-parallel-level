subroutine chm_heapre(&
     ielem,pnode,pgaus,lnods,elcon,eltem,gpsha,gpcar, gplap, &
                gp_enthalpy,gp_Cp,gp_grad_conce,gp_grad_temper, grad_enthalpy, gp_lap_con)
 
  !-----------------------------------------------------------------------
  !****f* chemic/chm_heapre
  ! NAME 
  !    chm_elmpre
  ! DESCRIPTION
  !    Compute quantities at Gauss points
  ! USES
  ! USED BY
  !    chm_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_domain, only      : ndime,mnode,mgaus
  use def_chemic, only      : nspec_chm,entha_chm
  use def_master,     only      : sphek
  implicit none
  integer(ip),  intent(in)  :: ielem,pnode,pgaus
  integer(ip),  intent(in)  :: lnods(pnode)
  real(rp),     intent(in)  :: elcon(pnode,nspec_chm)
  real(rp),     intent(in)  :: eltem(pnode)
  real(rp),     intent(in)  :: gpsha(pnode,pgaus)
  real(rp),     intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)  :: gplap(mgaus,mnode)                    ! Laplacian
  real(rp),     intent(out) :: gp_enthalpy(pgaus,nspec_chm)    ! Enthalpy per species
  real(rp),     intent(out) :: gp_Cp(pgaus,nspec_chm)          ! Specific heat
  real(rp),     intent(out) :: gp_grad_temper(ndime,pgaus), gp_grad_conce(ndime,pgaus,nspec_chm) 
  real(rp),     intent(out) :: grad_enthalpy(ndime,pgaus,nspec_chm), gp_lap_con(pgaus,nspec_chm) 
  integer(ip)               :: igaus,ispec,inode,idime,ipoin

  ! Convert to Gauss points and compute concentration and temperature gradients

  gp_enthalpy=0.0_rp
  gp_cp=0.0_rp
  gp_grad_conce=0.0_rp
  gp_grad_temper=0.0_rp
  grad_enthalpy=0.0_rp
  gp_lap_con=0.0_rp
  do ispec = 1,nspec_chm
     do igaus = 1,pgaus
        gp_enthalpy(igaus,ispec) = 0.0_rp
        gp_Cp(igaus,ispec) = 0.0_rp
        do inode = 1,pnode
           ipoin = lnods(inode) 
           gp_enthalpy(igaus,ispec) = gp_enthalpy(igaus,ispec) + gpsha(inode,igaus) * entha_chm(ipoin,ispec) 
           gp_Cp(igaus,ispec) = gp_Cp(igaus,ispec) + gpsha(inode,igaus) * sphek(ipoin,ispec) 
           do idime=1,ndime
              gp_grad_conce(idime,igaus,ispec) = gp_grad_conce(idime,igaus,ispec) + &
                   &  gpcar(idime,inode,igaus) * elcon(inode,ispec)
              grad_enthalpy(idime,igaus,ispec) = grad_enthalpy(idime,igaus,ispec) + &
                   gpcar(idime,inode,igaus) * entha_chm(ipoin,ispec) 
           enddo
           gp_lap_con (igaus,ispec) = gp_lap_con (igaus,ispec) + gplap(igaus,inode) * elcon(inode,ispec)         
        enddo
     enddo
  enddo

  do igaus = 1,pgaus
     do inode = 1,pnode
        do idime = 1,ndime
           gp_grad_temper(idime,igaus) = gp_grad_temper(idime,igaus) + gpcar(idime,inode,igaus) * eltem(inode)
        enddo
     enddo
  enddo
  
end subroutine chm_heapre
