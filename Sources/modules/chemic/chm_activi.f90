subroutine chm_activi(&
     pnode,pgaus,gp_nabla,gp_laplacian,el_temperature,el_massfrac,gp_gra_activ,gp_lap_activ)
  !-----------------------------------------------------------------------
  !****f* partis/chm_activi
  !    Compute chemical activities from coupling to LEVEL SET module
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mgaus,mnode
  use def_chemic, only      :  nspec_chm,interaction_chm
  use def_master, only      :  speci,wmean
  use mod_ker_proper, only  :  gasco

  implicit none
  integer(ip),  intent(in)  :: pnode,pgaus
  real(rp),     intent(in)  :: gp_nabla(ndime,pnode,pgaus)
  real(rp),     intent(in)  :: gp_laplacian(pnode,pgaus)                  
  real(rp),intent(in)       :: el_temperature(pnode)
  real(rp),intent(in)       :: el_massfrac(pnode,nspec_chm)
  real(rp),intent(out)      :: gp_gra_activ(ndime,pgaus,nspec_chm)          ! Gradient(activity) / activity
  real(rp),intent(out)      :: gp_lap_activ(pgaus,nspec_chm)                ! Laplacian(activity)/activity
  integer(ip)               :: igaus,ispec,inode,idime,jspec
  real(rp)                  :: zeta2, L(nspec_chm), sum_r(pnode), sum_q(pnode), sum_tau(pnode,nspec_chm),sum_L(pnode)
  real(rp)                  :: tau_interaction(pnode,nspec_chm,nspec_chm)
  real(rp)                  :: logActivity(pnode,nspec_chm),  el_molefrac(pnode,nspec_chm)

  zeta2 = 5.0_rp !   Z/2, but Z is always 10
  ! We need mole fractions everywhere
  do ispec=1,nspec_chm
     do inode = 1,pnode
        el_molefrac(inode,ispec) = el_massfrac(inode,ispec) * wmean(inode,1) / speci(ispec)%weigh
     enddo
  enddo

  sum_r = 0.0_rp
  sum_q = 0.0_rp
  sum_L = 0.0_rp
  tau_interaction=0.0_rp
  do ispec=1,nspec_chm
     L(ispec) = zeta2 * (speci(ispec)%activ(1)-speci(ispec)%activ(2))+1.0_rp-speci(ispec)%activ(1)
     do inode = 1,pnode
        do jspec=1, nspec_chm
           tau_interaction(inode,ispec,jspec) = exp(-interaction_chm(ispec,jspec)/(gasco*el_temperature(inode)))
        enddo
        sum_r(inode) = sum_r(inode) + speci(ispec)%activ(1) * el_molefrac(inode,ispec)
        sum_q(inode) = sum_q(inode) + speci(ispec)%activ(2) * el_molefrac(inode,ispec)
        sum_L(inode) = sum_L(inode) + L(ispec) * el_molefrac(inode,ispec)
     enddo
  enddo
  sum_tau = 0.0_rp
  do jspec=1,nspec_chm
     do ispec=1,nspec_chm
        do inode = 1,pnode
           sum_tau(inode,jspec) = sum_tau(inode,jspec) + speci(ispec)%activ(2) * el_molefrac(inode,ispec) &
                * tau_interaction(inode,ispec,jspec)/sum_q(inode)
        enddo
     enddo
  enddo
logActivity=0.0_rp
  do ispec=1,nspec_chm
     do inode = 1,pnode
        logActivity(inode,ispec) = logActivity(inode,ispec) + log( speci(ispec)%activ(1) / sum_r(inode) ) &
             + zeta2 *  speci(ispec)%activ(2) * log (speci(ispec)%activ(2)*sum_r(inode)/(speci(ispec)%activ(1)*sum_q(inode)) ) &
             + L(ispec) - (speci(ispec)%activ(1) / sum_r(inode))*sum_L(inode) &
             + speci(ispec)%activ(2) * (1.0_rp - log(sum_tau(inode,ispec))) 
        do jspec = 1,nspec_chm
           logActivity(inode,ispec) = logActivity(inode,ispec) -  speci(ispec)%activ(2) &
               *  (speci(jspec)%activ(2) * el_molefrac(inode,jspec) / sum_q(inode) ) &
               *  tau_interaction(inode,ispec,jspec) * sum_tau(inode,jspec)
        enddo
     enddo
  enddo

  gp_gra_activ = 0.0_rp
  gp_lap_activ = 0.0_rp
  do ispec=1,nspec_chm
     do igaus=1,pgaus
        do inode = 1,pnode
           do idime = 1,ndime
              gp_gra_activ(idime,igaus,ispec)= gp_gra_activ(idime,igaus,ispec) + gp_nabla(idime,inode,igaus) * logActivity(inode,ispec)
           enddo
           gp_lap_activ(igaus,ispec) = gp_lap_activ(igaus,ispec) + gp_laplacian(igaus,inode) * logActivity(inode,ispec)
        enddo
        do idime=1,ndime
           gp_lap_activ(igaus,ispec) = gp_lap_activ(igaus,ispec) + gp_gra_activ(idime,igaus,ispec)*gp_gra_activ(idime,igaus,ispec)
        enddo
     enddo
  enddo

!  gp_gra_activ = 0.0_rp
!  gp_lap_activ = 0.0_rp

end subroutine chm_activi
