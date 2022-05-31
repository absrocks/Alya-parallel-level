subroutine chm_gp_reatab()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_reatab
  ! NAME 
  !    chm_reatab
  ! DESCRIPTION
  !    Read table properties for Flamelet combustion model
  ! USES
  ! USED BY
  !    chm_iniunk: initialize flow field
  !    chm_endite: update table properties to start doiter in temper with updated coefficients 
  !    chm_endste: properties available for all modules at the end of time-step
  !                (there is no update of properties during chemic iteration)   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only          : ip,rp
  use def_master, only          : wmean,conce,inotmaster, &
       sphec,therm,kfl_paral
  use def_master, only          : condu_gp,sphec_gp,visco_gp,wmean_gp,&
       sphec_gp_ht,sphec_gp_lt 
  use def_chemic, only          : mass_gp,rrt_gp 

  use def_domain, only          : npoin,ltypb,nnode,lnodb,nboun,nelem,ltype,&
       ngaus,llapl,lorde,ltopo,elmar,hnatu,lnods,&
       mgaus,ndime,mnode,ntens
  use def_chemic, only          : kfl_radia_chm,kfl_wallc_chm, &
       kfl_premix_chm,kfl_varYc_chm,kfl_varZ_chm, &
       kfl_ufpv_chm, xZr_chm, xZs_chm, kfl_lookg_chm, &
       kfl_advec_chm,kfl_ellen_chm,ADR_chm,nclas_chm,&
       kfl_heat_loss_chm
  use def_chemic,      only     : table_fw
  use mod_ker_proper,  only     : ker_proper
  use mod_interp_tab,  only     : fw_lookup 

  implicit none
  integer(ip)               :: ipoin,iclas,ivalu
  integer(ip)               :: iboun,pblty,inodb,pnodb
  real(rp)                  :: retva(table_fw % main_table % nvar)              ! Values read in from flamelet table: 
  !                                                                               (1) S_c (2) Wmean (3) Lambda (4) Mu (5-10) cp_low (11-16) cp_high (17) S_c*c
  real(rp)                  :: tab_scale_control(table_fw % main_table % ndim)  ! Concentration for each variable at each node
  real(rp)                  :: control(table_fw % main_table % ndim)            ! input of table lookup function 

  integer(ip) :: ielem,igaus,inode
  integer(ip) :: dummi
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo
  real(rp)    :: gphco(mgaus)
  real(rp)    :: gpsph(mgaus)
  real(rp)    :: gp_res_Chi_Yc(mgaus)
  real(rp)    :: gp_sgs_Chi_Yc(mgaus)
  real(rp)    :: gp_res_Chi_z(mgaus)
  real(rp)    :: gp_sgs_Chi_z(mgaus)
  real(rp)    :: gpcar(ndime,mnode,mgaus)                 
  real(rp)    :: gphes(ntens,mnode,mgaus)
  real(rp)    :: gpcon(mgaus,nclas_chm)
  real(rp)    :: gpthe(mgaus)
  real(rp)    :: gpden(mgaus)
  real(rp)    :: gptur(mgaus)
  real(rp)    :: gpvol(mgaus)
  real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)
  real(rp)    :: elcon(mnode,nclas_chm,ADR_chm(1) % ntime)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: dummr(ndime,mnode)
  integer(ip) :: lnods_loc(mnode)
  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem
     !
     ! Element dimensions
     !
     pelty = ltype(ielem)
     if( pelty > 0 ) then
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = llapl(pelty) 
        porde = lorde(pelty)
        ptopo = ltopo(pelty)

        !
        ! Gather all
        !
        lnods_loc(1:pnode) = lnods(1:pnode,ielem)
        call chm_post_gather(&
             pnode,lnods_loc,elcon(1:pnode,:,:),elcod,dummr)


        if (kfl_ufpv_chm > 0) then
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
                pnode,pgaus,elcon(1:pnode,:,:),dummr,elmar(pelty)%shape,gpcar,hleng,gptur,&
                gphco,gpsph,gpden,gp_res_Chi_Yc,gp_sgs_Chi_Yc,gp_res_Chi_z,gp_sgs_Chi_z)      
        endif

        !
        ! Initialization variables
        !
        gpcon = 0.0_rp 
        gpthe = 0.0_rp 
        condu_gp(ielem) % a      = 0.0_rp
        sphec_gp(ielem) % a      = 0.0_rp
        visco_gp(ielem) % a      = 0.0_rp
        do igaus = 1,pgaus
           !
           ! Only put current value to 0
           ! We need to save old steps
           !
           wmean_gp(ielem) % a (igaus,1,1)     = 0.0_rp
        enddo
        sphec_gp_ht(ielem) % a  = 0.0_rp
        sphec_gp_lt(ielem) % a  = 0.0_rp
        mass_gp(ielem) % a      = 0.0_rp
        rrt_gp(ielem) % a       = 0.0_rp

        !
        ! Concentration
        !
        do iclas = 1,nclas_chm
           do igaus = 1,pgaus
              do inode = 1,pnode
                 gpcon(igaus,iclas) = gpcon(igaus,iclas)&
                      + elmar(pelty)%shape(inode,igaus) * elcon(inode,iclas,1)
              end do
           end do
        end do

        if (kfl_heat_loss_chm /= 0) then
           do igaus = 1,pgaus
              do inode = 1,pnode
                 ipoin=lnods_loc(inode)
                 gpthe(igaus) = gpthe(igaus)&
                      +  elmar(pelty)%shape(inode,igaus) * therm(ipoin,1)
              end do
           end do
        endif


        !
        ! Lookup
        !
        do igaus = 1,pgaus
           control = 0.0_rp
           do iclas = 1, table_fw % main_table % ndim
              select case (table_fw % main_table % coords(iclas) % name)
              case ('CMEAN','C    ')
                  control(iclas) = gpcon(igaus,1)
              case ('CVAR ')
                  control(iclas) = gpcon(igaus,2)
              case ('CHIST')
                  control(iclas) = gp_res_Chi_z(igaus) + gp_sgs_Chi_z(igaus)
              case ('ZMEAN','Z    ')
                  control(iclas) = gpcon(igaus,3)
              case ('ZVAR ')
                  control(iclas) = gpcon(igaus,4)
              case ('IMEAN','I    ')
                  control(iclas) = gpthe(igaus)
              end select
           enddo
           call fw_lookup( control, tab_scale_control, table_fw, retva )

           !
           ! update properties 
           !
           select case(kfl_ufpv_chm )
           case (0)
              mass_gp(ielem) % a(igaus,1,1)  = retva(1)
              rrt_gp(ielem) % a(igaus,2,1)   = retva(17)  ! filter{w_c * c}
           case (1)
              mass_gp(ielem) % a(igaus,1,1)  = 0.0_rp
              rrt_gp(ielem) % a(igaus,2,1)   = 0.0_rp
           case (2)
              mass_gp(ielem) % a(igaus,1,1)  = retva(1)   ! omega_Yc
              rrt_gp(ielem) % a(igaus,2,1)   = 0.0_rp
           case (3)
              mass_gp(ielem) % a(igaus,1,1)  = retva(17)  ! rho * dYc/dt
              rrt_gp(ielem) % a(igaus,2,1)   = 0.0_rp
           end select

           wmean_gp(ielem) % a(igaus,1,1)   = retva(2)
           condu_gp(ielem) % a(igaus,1,1)   = retva(3)
           visco_gp(ielem) % a(igaus,1,1)   = retva(4)

           do ivalu = 1,6 
              !
              ! First is low temperature, second is high temperature in table
              !
              sphec_gp_lt(ielem) % a(igaus,ivalu,1) = retva(4+ivalu)
              sphec_gp_ht(ielem) % a(igaus,ivalu,1) = retva(10+ivalu)
           end do

           if (kfl_radia_chm > 0) then
              call runend('Chemic gp_reatab: Optically thin radiation model is not implemented.')
           end if
        enddo   !igaus loop

     endif
  end do elements
  !
  ! LOOKUP Cp coefficients on boundary for temperature BC:
  !
  boundaries: do iboun = 1,nboun
     pblty = ltypb(iboun)
     pnodb = nnode(pblty)
     do inodb = 1,pnodb
        ipoin   = lnodb(inodb,iboun)
        control = 0.0_rp

        do iclas = 1, table_fw % main_table % ndim
           select case (table_fw % main_table % coords(iclas) % name)
           case ('CMEAN','C    ')
               control(iclas) = conce(ipoin,1,1)
           case ('CVAR ')
               control(iclas) = conce(ipoin,2,1)
           case ('CHIST')
               control(iclas) = xZr_chm(ipoin) + xZs_chm(ipoin)
           case ('ZMEAN','Z    ')
               control(iclas) = conce(ipoin,3,1)
           case ('ZVAR ')
               control(iclas) = conce(ipoin,4,1)
           case ('IMEAN','I    ')
               if (kfl_heat_loss_chm /= 0) control(iclas) = therm(ipoin,1)
           end select
        enddo

        call fw_lookup( control, tab_scale_control, table_fw, retva )

        do ivalu = 1,6 
           sphec(ipoin,ivalu,1) = retva(4+ivalu)
           sphec(ipoin,ivalu,2) = retva(10+ivalu)
        end do

     end do
  end do boundaries

end subroutine chm_gp_reatab
