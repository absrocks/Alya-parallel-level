subroutine chm_post_scalar_dissipation_rate(ivari) 
   !------------------------------------------------------------------------
   ! NAME 
   !    chm_post_scalar_dissipation_rate
   ! DESCRIPTION
   !    Elemental operations for post-processing scalar dissipation rates 
   !    in flamelet combustion model
   ! USES
   ! USED BY
   !    chm_outvar
   !***
   !------------------------------------------------------------------------
   use def_kintyp, only      : ip,rp
   use def_domain, only      : ndime,mnode,mgaus,nelem,ltype,   &
                               nnode,ngaus,llapl,lorde,ltopo,   &
                               lnods,elmar,hnatu,ntens

   use def_chemic, only      : kfl_advec_chm,kfl_ellen_chm,     &
                               xYr_chm,xZr_chm,xYs_chm,xZs_chm, &
                               ADR_chm,nclas_chm,ncomp_chm
   use mod_ker_proper, only  : ker_proper
   use mod_solver, only      : solver_lumped_mass_system 

   implicit none
   integer(ip),  intent(in) :: ivari

   integer(ip) :: ielem
   integer(ip) :: pelty,pnode
   integer(ip) :: pgaus,plapl,porde,ptopo
   integer(ip) :: dummi

   real(rp)    :: elcon(mnode,nclas_chm,ADR_chm(1) % ntime)
   real(rp)    :: elcod(ndime,mnode)
   real(rp)    :: elvel(ndime,mnode)
   real(rp)    :: gpvol(mgaus)
   real(rp)    :: gphco(mgaus)
   real(rp)    :: gpsph(mgaus)
   real(rp)    :: gpcon(mgaus,nclas_chm,ncomp_chm)
   real(rp)    :: gpden(mgaus)
   real(rp)    :: gptur(mgaus)
   real(rp)    :: gpcar(ndime,mnode,mgaus)                 
   real(rp)    :: gphes(ntens,mnode,mgaus)
   real(rp)    :: gp_res_Chi_Yc(mgaus)
   real(rp)    :: gp_sgs_Chi_Yc(mgaus)
   real(rp)    :: gp_res_Chi_z(mgaus)
   real(rp)    :: gp_sgs_Chi_z(mgaus)
   real(rp)    :: gp_Chi_out(mgaus)
 

   real(rp)    :: dummr(mgaus*ndime)
   real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)

   !
   ! Initialization
   !
   select case (ivari)

     case (46_ip) 
       xYr_chm = 0.0_rp

     case (47_ip) 
       xZr_chm = 0.0_rp  

     case (48_ip) 
       xYs_chm = 0.0_rp

     case (49_ip) 
       xZs_chm = 0.0_rp

   end select

   !
   ! Loop over elements
   !

   elements: do ielem = 1,nelem

      !
      ! Element dimensions
      !
      pelty = ltype(ielem)
      pnode = nnode(pelty)
      pgaus = ngaus(pelty)
      plapl = llapl(pelty) 
      porde = lorde(pelty)
      ptopo = ltopo(pelty)

      !
      ! Initialization variables
      !
      gptur = 0.0_rp
      gpcon = 0.0_rp 

      gp_res_Chi_Yc = 0.0_rp
      gp_sgs_Chi_Yc = 0.0_rp
      gp_res_Chi_z  = 0.0_rp
      gp_sgs_Chi_z  = 0.0_rp

      !
      ! Gather all
      !
      call chm_post_gather(&
           pnode,lnods(1:pnode,ielem),elcon(1:pnode,:,:),elcod,elvel)

      !
      ! CHALE, HLENG and TRAGL 
      !
      call elmlen(&
           ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
           hleng)
      call elmchl(&
           tragl,hleng,elcod,dummr,chave,chale,pnode,&
           porde,hnatu(pelty),kfl_advec_chm,kfl_ellen_chm)

      !
      ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, GPVOL
      !
      call elmcar(&
           pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,elmar(pelty)%deriv, &
           elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)

      call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
      call ker_proper('CONDU','PGAUS',dummi,ielem,gphco,pnode,pgaus,elmar(pelty)%shape,gpcar)
      call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)
      call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty)%shape,gpcar)

      !
      ! Compute scalar dissipation rate at gauss points
      !
      call chm_calc_scalar_dissip(&
           pnode,pgaus,elcon(1:pnode,:,:),elvel,elmar(pelty)%shape,gpcar,hleng,gptur,&
           gphco,gpsph,gpden,gp_res_Chi_Yc,gp_sgs_Chi_Yc,gp_res_Chi_z,gp_sgs_Chi_z)           

      select case (ivari)

        case (46_ip) 
           gp_Chi_out = gp_res_Chi_Yc

        case (47_ip) 
           gp_Chi_out = gp_res_Chi_Z

        case (48_ip) 
           gp_Chi_out = gp_sgs_Chi_Yc

        case (49_ip) 
           gp_Chi_out = gp_sgs_Chi_Z

      end select

      !
      ! Projection 
      !
      call chm_post_projection(ivari,ielem,pnode,pgaus,elmar(pelty)%shape,gp_Chi_out,gpvol)

   end do elements

   !
   ! Lumped mass matrix 
   ! 
   select case (ivari)

     case (46_ip) 
        call solver_lumped_mass_system(1_ip,xYr_chm)

     case (47_ip) 
        call solver_lumped_mass_system(1_ip,xZr_chm)

     case (48_ip) 
        call solver_lumped_mass_system(1_ip,xYs_chm)

     case (49_ip) 
        call solver_lumped_mass_system(1_ip,xZs_chm)

   end select

end subroutine chm_post_scalar_dissipation_rate
