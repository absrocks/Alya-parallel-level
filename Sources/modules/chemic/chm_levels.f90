subroutine chm_levels(&
     pnode,pgaus,gplap,ellev,gplev,gpgle,gpgac,gplac)

  !-----------------------------------------------------------------------
  !****f* partis/chm_activi
  !    Compute chemical activities from coupling to LEVEL SET module
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mgaus,mnode
  use def_chemic, only      :  nspec_chm
  use def_master, only      :  speci

  implicit none
  integer(ip),  intent(in)  :: pnode,pgaus
  real(rp),     intent(in)  :: gplap(mgaus,mnode)                  
  real(rp),intent(in) :: ellev(pnode)
  real(rp),intent(in) :: gplev(pgaus)
  real(rp),intent(in) :: gpgle(ndime,pgaus)
  real(rp),intent(out)    :: gpgac(ndime,pgaus,nspec_chm)          ! Gradient(activity) / activity
  real(rp),intent(out)    :: gplac(pgaus,nspec_chm)                ! Laplacian(activity)/activity
  integer(ip)               :: igaus,ispec,inode
  real(rp) :: Activity_1, Activity_2,delta_A, gp_Activity

  do ispec=1,nspec_chm
     Activity_1 = speci(ispec)%activ(1) ! Phase 1
     Activity_2 = speci(ispec)%activ(2) ! Phase 2
     delta_A =  Activity_2 - Activity_1 
     do igaus = 1,pgaus
        gp_Activity = ( gplev(igaus) * Activity_1 + (1.0_rp - gplev(igaus)) * delta_A)
        gplac (igaus,ispec) = 0.0_rp
        do inode = 1,pnode
           gplac (igaus,ispec) = gplac (igaus,ispec) + gplap(igaus,inode) * ( ellev(inode) * Activity_1 + (1.0_rp - ellev(inode)) * delta_A)
        enddo
        gplac (igaus,ispec) = ( gplac (igaus,ispec) / gp_Activity )
        gpgac(:,igaus,ispec) = ( gpgle(:,igaus) * Activity_1 + (1.0_rp - gpgle(:,igaus)) * delta_A) / gp_Activity
     enddo
  enddo

end subroutine chm_levels
